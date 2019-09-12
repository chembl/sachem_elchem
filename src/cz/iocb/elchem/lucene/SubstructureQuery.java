package cz.iocb.elchem.lucene;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.TimeoutException;
import org.apache.lucene.index.BinaryDocValues;
import org.apache.lucene.index.DocValues;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.index.Term;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.BooleanQuery.Builder;
import org.apache.lucene.search.ConstantScoreQuery;
import org.apache.lucene.search.DisjunctionMaxQuery;
import org.apache.lucene.search.DocIdSetIterator;
import org.apache.lucene.search.DocValuesFieldExistsQuery;
import org.apache.lucene.search.Explanation;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreMode;
import org.apache.lucene.search.Scorer;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.Weight;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.elchem.fingerprint.IOCBFingerprint;
import cz.iocb.elchem.molecule.AromaticityMode;
import cz.iocb.elchem.molecule.BinaryMolecule;
import cz.iocb.elchem.molecule.BinaryMoleculeBuilder;
import cz.iocb.elchem.molecule.ChargeMode;
import cz.iocb.elchem.molecule.Isomorphism;
import cz.iocb.elchem.molecule.IsotopeMode;
import cz.iocb.elchem.molecule.MoleculeCreator;
import cz.iocb.elchem.molecule.SearchMode;
import cz.iocb.elchem.molecule.StereoMode;
import cz.iocb.elchem.molecule.TautomerMode;



public class SubstructureQuery extends Query
{
    private final String field;
    private final String query;
    private final SearchMode searchMode;
    private final ChargeMode chargeMode;
    private final IsotopeMode isotopeMode;
    private final StereoMode stereoMode;
    private final AromaticityMode aromaticityMode;
    private final TautomerMode tautomerMode;

    private final Query subquery;


    public SubstructureQuery(String field, String query, SearchMode graphMode, ChargeMode chargeMode,
            IsotopeMode isotopeMode, StereoMode stereoMode, AromaticityMode aromaticityMode, TautomerMode tautomerMode)
            throws CDKException, IOException, TimeoutException
    {
        this.field = field;
        this.query = query;
        this.searchMode = graphMode;
        this.chargeMode = chargeMode;
        this.isotopeMode = isotopeMode;
        this.stereoMode = stereoMode;
        this.aromaticityMode = aromaticityMode;
        this.tautomerMode = tautomerMode;

        List<IAtomContainer> queryMolecules = MoleculeCreator.translateQuery(query, aromaticityMode, tautomerMode);
        ArrayList<Query> subqueries = new ArrayList<Query>(queryMolecules.size());

        for(IAtomContainer molecule : queryMolecules)
            subqueries.add(new SingleSubstructureQuery(molecule));

        this.subquery = new DisjunctionMaxQuery(subqueries, 0);
    }


    @Override
    public Query rewrite(IndexReader reader)
    {
        return subquery;
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


    private boolean equalsTo(SubstructureQuery other)
    {
        return field.equals(other.field) && query.equals(other.query) && searchMode.equals(other.searchMode)
                && chargeMode.equals(other.chargeMode) && isotopeMode.equals(other.isotopeMode)
                && stereoMode.equals(other.stereoMode) && aromaticityMode.equals(other.aromaticityMode)
                && tautomerMode.equals(other.tautomerMode);
    }


    @Override
    public int hashCode()
    {
        int result = classHash();
        result = 31 * result + field.hashCode();
        result = 31 * result + query.hashCode();
        result = 3 * result + searchMode.hashCode();
        result = 3 * result + chargeMode.hashCode();
        result = 3 * result + isotopeMode.hashCode();
        result = 3 * result + stereoMode.hashCode();
        result = 3 * result + aromaticityMode.hashCode();
        result = 3 * result + tautomerMode.hashCode();
        return result;
    }


    @Override
    public String toString(String field)
    {
        //TODO:
        return "SubstructureSearchQuery(...)";
    }


    class SingleSubstructureQuery extends Query
    {
        private final Query parentQuery;
        private final IAtomContainer tautomer;

        boolean extended;
        private final Set<Integer> fp;
        private final Map<Integer, Set<Integer>> info;
        private final BinaryMolecule molecule;
        private final byte[] moleculeData;
        private final int originalMoleculeSize;
        private final boolean[] restH;
        private final Isomorphism isomorphism;


        SingleSubstructureQuery(IAtomContainer tautomer) throws CDKException, IOException
        {
            this.parentQuery = SubstructureQuery.this;
            this.tautomer = tautomer;

            this.moleculeData = (new BinaryMoleculeBuilder(tautomer)).asBytes(searchMode == SearchMode.EXACT);
            this.restH = new boolean[tautomer.getAtomCount()];

            for(int i = 0; i < tautomer.getAtomCount(); i++)
                restH[i] = Boolean.TRUE.equals(tautomer.getAtom(i).getProperty(CDKConstants.REST_H));


            this.extended = BinaryMolecule.isExtended(moleculeData, chargeMode != ChargeMode.IGNORE,
                    isotopeMode != IsotopeMode.IGNORE);

            this.molecule = new BinaryMolecule(moleculeData, restH, extended, chargeMode != ChargeMode.IGNORE,
                    isotopeMode != IsotopeMode.IGNORE, stereoMode != StereoMode.IGNORE, false, false);

            this.originalMoleculeSize = molecule.getOriginalAtomCount() + molecule.getOriginalBondCount();
            this.isomorphism = new Isomorphism(molecule, searchMode, chargeMode, isotopeMode, stereoMode);

            this.info = new HashMap<Integer, Set<Integer>>();
            this.fp = IOCBFingerprint.getSubstructureFingerprint(molecule, info);
        }


        @Override
        public Weight createWeight(IndexSearcher searcher, ScoreMode scoreMode, float boost) throws IOException
        {
            return new SingleSubstructureWeight(searcher, scoreMode, boost);
        }


        @Override
        public boolean equals(Object other)
        {
            return sameClassAs(other) && equalsTo(getClass().cast(other));
        }


        private boolean equalsTo(SingleSubstructureQuery other)
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
            return "SubstructureSearchQuery.SingleSubstructureQuery(...)";
        }


        class SingleSubstructureWeight extends Weight
        {
            private final Weight innerWeight;


            public SingleSubstructureWeight(IndexSearcher searcher, ScoreMode scoreMode, float boost) throws IOException
            {
                super(SubstructureQuery.this);

                if(!fp.isEmpty())
                {
                    Builder builder = new BooleanQuery.Builder();
                    FingerprintBitMapping mapping = new FingerprintBitMapping();

                    for(int bit : selectFingerprintBits(searcher))
                        builder.add(new TermQuery(new Term(field, mapping.bitAsString(bit))), BooleanClause.Occur.MUST);

                    this.innerWeight = new ConstantScoreQuery(builder.build()).createWeight(searcher,
                            ScoreMode.COMPLETE_NO_SCORES, boost);
                }
                else
                {
                    this.innerWeight = new DocValuesFieldExistsQuery(field).createWeight(searcher,
                            ScoreMode.COMPLETE_NO_SCORES, boost);
                }
            }


            @Override
            public Scorer scorer(LeafReaderContext context) throws IOException
            {
                Scorer scorer = innerWeight.scorer(context);

                if(scorer == null)
                    return null;

                return new SingleSubstructureScorer(context, scorer);
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


            private List<Integer> selectFingerprintBits(IndexSearcher searcher) throws IOException
            {
                final int maxSize = 32;
                final int atomCoverage = 2;

                Map<Integer, Integer> ordered = new TreeMap<Integer, Integer>();
                FingerprintBitMapping mapping = new FingerprintBitMapping();

                for(int i : fp)
                    ordered.put(searcher.getIndexReader().docFreq(new Term(field, mapping.bitAsString(i))), i);

                List<Integer> selected = new ArrayList<Integer>(maxSize);
                int[] coverage = new int[molecule.getAtomCount()];
                int uncovered = molecule.getAtomCount();

                for(Integer i : ordered.values())
                {
                    if(uncovered <= 0 || selected.size() >= maxSize)
                        break;

                    boolean found = false;

                    for(int a : info.get(i))
                    {
                        if(coverage[a] < atomCoverage)
                        {
                            found = true;
                            coverage[a]++;

                            if(coverage[a] == atomCoverage)
                                uncovered--;
                        }
                    }

                    if(found)
                        selected.add(i);
                }

                return selected;
            }


            class SingleSubstructureScorer extends Scorer
            {
                private int docID = -1;
                private float score = 0;
                private final Scorer innerScorer;
                private final BinaryDocValues molDocValue;


                protected SingleSubstructureScorer(LeafReaderContext context, Scorer scorer) throws IOException
                {
                    super(SingleSubstructureWeight.this);
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
                    byte[] data = molDocValue.binaryValue().bytes;

                    BinaryMolecule target;
                    boolean match;

                    if(!extended && (BinaryMolecule.hasPseudoAtom(data) || BinaryMolecule.hasMultivalentHydrogen(data)))
                    {
                        BinaryMolecule molecule = new BinaryMolecule(moleculeData, restH, true,
                                chargeMode != ChargeMode.IGNORE, isotopeMode != IsotopeMode.IGNORE,
                                stereoMode != StereoMode.IGNORE, false, false);

                        Isomorphism isomorphism = new Isomorphism(molecule, searchMode, chargeMode, isotopeMode,
                                stereoMode);


                        target = new BinaryMolecule(data, null, true, chargeMode != ChargeMode.IGNORE,
                                isotopeMode != IsotopeMode.IGNORE, stereoMode != StereoMode.IGNORE, false, false);

                        match = isomorphism.match(target);
                    }
                    else
                    {
                        target = new BinaryMolecule(data, null, extended, chargeMode != ChargeMode.IGNORE,
                                isotopeMode != IsotopeMode.IGNORE, stereoMode != StereoMode.IGNORE,
                                chargeMode == ChargeMode.DEFAULT_AS_UNCHARGED,
                                isotopeMode == IsotopeMode.DEFAULT_AS_STANDARD);

                        match = isomorphism.match(target);
                    }

                    if(match)
                    {
                        int size = (target.getOriginalAtomCount() + target.getOriginalBondCount());
                        score = size == 0 ? 0.0f : originalMoleculeSize / (float) size;
                    }

                    return match;
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
