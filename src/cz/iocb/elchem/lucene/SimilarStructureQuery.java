package cz.iocb.elchem.lucene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.TimeoutException;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.index.BinaryDocValues;
import org.apache.lucene.index.DocValues;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.Term;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.BooleanQuery.Builder;
import org.apache.lucene.search.ConstantScoreQuery;
import org.apache.lucene.search.DisjunctionMaxQuery;
import org.apache.lucene.search.DocIdSetIterator;
import org.apache.lucene.search.Explanation;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreMode;
import org.apache.lucene.search.Scorer;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.Weight;
import org.apache.lucene.util.BytesRef;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.elchem.fingerprint.IOCBFingerprint;
import cz.iocb.elchem.molecule.AromaticityMode;
import cz.iocb.elchem.molecule.BinaryMolecule;
import cz.iocb.elchem.molecule.BinaryMoleculeBuilder;
import cz.iocb.elchem.molecule.MoleculeCreator;
import cz.iocb.elchem.molecule.TautomerMode;



public class SimilarStructureQuery extends Query
{
    private final String field;
    private final String query;

    private final AromaticityMode aromaticityMode;
    private final TautomerMode tautomerMode;
    private final float threshold;
    private final Query subquery;


    public SimilarStructureQuery(String field, String query, float threshold, AromaticityMode aromaticityMode,
            TautomerMode tautomerMode) throws CDKException, IOException, TimeoutException
    {
        this.field = field;
        this.query = query;
        this.threshold = threshold;
        this.aromaticityMode = aromaticityMode;
        this.tautomerMode = tautomerMode;

        List<IAtomContainer> queryMolecules = MoleculeCreator.translateQuery(query, aromaticityMode, tautomerMode);
        ArrayList<Query> subqueries = new ArrayList<Query>(queryMolecules.size());

        for(IAtomContainer molecule : queryMolecules)
            subqueries.add(new SingleSimilarityQuery(molecule));

        this.subquery = new DisjunctionMaxQuery(subqueries, 0);
    }


    @Override
    public Weight createWeight(IndexSearcher searcher, ScoreMode scoreMode, float boost) throws IOException
    {
        return subquery.createWeight(searcher, scoreMode, boost);
    }


    @Override
    public boolean equals(Object other)
    {
        return sameClassAs(other) && equalsTo(getClass().cast(other));
    }


    private boolean equalsTo(SimilarStructureQuery other)
    {
        return field.equals(other.field) && query.equals(other.query) && aromaticityMode.equals(other.aromaticityMode)
                && tautomerMode.equals(other.tautomerMode);
    }


    @Override
    public int hashCode()
    {
        int result = classHash();
        result = 31 * result + field.hashCode();
        result = 31 * result + query.hashCode();
        result = 3 * result + aromaticityMode.hashCode();
        result = 3 * result + tautomerMode.hashCode();
        return result;
    }


    @Override
    public String toString(String field)
    {
        //TODO:
        return "SimilaritySearchQuery(...)";
    }


    class SingleSimilarityQuery extends Query
    {
        private final Query parentQuery;
        private final IAtomContainer tautomer;

        private final Set<Integer> fp;


        SingleSimilarityQuery(IAtomContainer tautomer) throws CDKException, IOException
        {
            this.parentQuery = SimilarStructureQuery.this;
            this.tautomer = tautomer;

            byte[] moleculeData = (new BinaryMoleculeBuilder(tautomer)).asBytes(false);

            BinaryMolecule molecule = new BinaryMolecule(moleculeData, null, false, false, false, false, false, false);

            this.fp = IOCBFingerprint.getSimilarityFingerprint(molecule);
        }


        @Override
        public Weight createWeight(IndexSearcher searcher, ScoreMode scoreMode, float boost) throws IOException
        {
            return new SingleSimilarityWeight(searcher, scoreMode, boost);
        }


        @Override
        public boolean equals(Object other)
        {
            return sameClassAs(other) && equalsTo(getClass().cast(other));
        }


        private boolean equalsTo(SingleSimilarityQuery other)
        {
            return parentQuery.equals(other.parentQuery) && tautomer.equals(other.tautomer);
        }


        @Override
        public int hashCode()
        {
            int result = classHash();
            result = 31 * result + parentQuery.hashCode();
            result = 31 * result + tautomer.hashCode();
            return result;
        }


        @Override
        public String toString(String field)
        {
            //TODO:
            return "SimilaritySearchQuery.SingleSimilarityQuery(...)";
        }


        class SingleSimilarityWeight extends Weight
        {
            private final Weight innerWeight;


            public SingleSimilarityWeight(IndexSearcher searcher, ScoreMode scoreMode, float boost) throws IOException
            {
                super(SimilarStructureQuery.this);

                Builder builder = new BooleanQuery.Builder();
                FingerprintBitMapping mapping = new FingerprintBitMapping();

                int min = (int) Math.floor(fp.size() * threshold);
                int max = (int) Math.ceil(fp.size() / threshold);

                builder.add(IntPoint.newRangeQuery(field, min, max), BooleanClause.Occur.MUST);

                for(int bit : selectFingerprintBits(searcher))
                    builder.add(new TermQuery(new Term(field, mapping.bitAsString(bit))), BooleanClause.Occur.SHOULD);

                builder.setMinimumNumberShouldMatch(1);

                this.innerWeight = new ConstantScoreQuery(builder.build()).createWeight(searcher, scoreMode, boost);
            }


            @Override
            public Scorer scorer(LeafReaderContext context) throws IOException
            {
                Scorer scorer = innerWeight.scorer(context);

                if(scorer == null)
                    return null;

                return new SingleSimilarityScorer(context, scorer);
            }


            @Override
            public boolean isCacheable(LeafReaderContext context)
            {
                return false;
            }


            @Override
            public Explanation explain(LeafReaderContext context, int doc) throws IOException
            {
                Scorer scorer = scorer(context);

                if(scorer != null && doc != scorer.iterator().advance(doc))
                    return Explanation.match(scorer.score(), "match");

                return Explanation.noMatch("no match");
            }


            @Deprecated
            @Override
            public void extractTerms(Set<Term> set)
            {
                innerWeight.extractTerms(set);
            }


            private Set<Integer> selectFingerprintBits(IndexSearcher searcher) throws IOException
            {
                int limit = 1 + (int) Math.ceil(fp.size() * (1 - threshold));

                Map<Integer, Integer> ordered = new TreeMap<Integer, Integer>();
                FingerprintBitMapping mapping = new FingerprintBitMapping();

                for(int i : fp)
                    ordered.put(searcher.getIndexReader().docFreq(new Term(field, mapping.bitAsString(i))), i);

                Set<Integer> selected = new HashSet<Integer>();

                int count = 0;

                for(Integer bit : ordered.values())
                {
                    selected.add(bit);

                    if(++count == limit)
                        break;
                }

                return selected;
            }


            class SingleSimilarityScorer extends Scorer
            {
                private int docID = -1;
                private float score = 0;
                private final Scorer innerScorer;
                private final BinaryDocValues molDocValue;


                protected SingleSimilarityScorer(LeafReaderContext context, Scorer scorer) throws IOException
                {
                    super(SingleSimilarityWeight.this);
                    this.innerScorer = scorer;
                    this.molDocValue = DocValues.getBinary(context.reader(), field);
                }


                @Override
                public int docID()
                {
                    return docID;
                }


                @Override
                public float getMaxScore(int upTo) throws IOException
                {
                    return 1.0f;
                }


                @Override
                public float score() throws IOException
                {
                    return score;
                }


                boolean isValid() throws IOException
                {
                    molDocValue.advanceExact(docID);
                    BytesRef data = molDocValue.binaryValue();

                    int shared = 0;

                    for(int i = 0; i < data.length / 4; i++)
                    {
                        int value = 0;

                        for(int j = 0; j < 4; j++)
                            value |= Byte.toUnsignedInt(data.bytes[i * 4 + j]) << (j * 8);

                        if(fp.contains(value))
                            shared++;
                    }

                    float similarity = shared / (float) (fp.size() + data.length / 4 - shared);

                    if(similarity < threshold)
                        return false;

                    score = similarity;
                    return true;
                }


                @Override
                public DocIdSetIterator iterator()
                {
                    DocIdSetIterator innerDocIdSetIterator = innerScorer.iterator();

                    return new DocIdSetIterator()
                    {
                        @Override
                        public int advance(int target) throws IOException
                        {
                            docID = innerDocIdSetIterator.advance(target);

                            if(docID != NO_MORE_DOCS && !isValid())
                                nextDoc();

                            return docID;
                        }


                        @Override
                        public int nextDoc() throws IOException
                        {
                            while(true)
                            {
                                docID = innerDocIdSetIterator.nextDoc();

                                if(docID == NO_MORE_DOCS || isValid())
                                    return docID;
                            }
                        }


                        @Override
                        public int docID()
                        {
                            return docID;
                        }


                        @Override
                        public long cost()
                        {
                            return innerDocIdSetIterator.cost();
                        }
                    };
                }
            }
        }
    }
}
