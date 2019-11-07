package cz.iocb.elchem.elasticsearch;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.lucene.document.BinaryDocValuesField;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.IndexOptions;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.search.Query;
import org.apache.lucene.util.BytesRef;
import org.elasticsearch.common.settings.Settings;
import org.elasticsearch.index.mapper.FieldMapper;
import org.elasticsearch.index.mapper.MappedFieldType;
import org.elasticsearch.index.mapper.Mapper;
import org.elasticsearch.index.mapper.MapperParsingException;
import org.elasticsearch.index.mapper.ParseContext;
import org.elasticsearch.index.mapper.TypeParsers;
import org.elasticsearch.index.query.QueryShardContext;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import cz.iocb.elchem.fingerprint.IOCBFingerprint;
import cz.iocb.elchem.lucene.FingerprintTokenStream;
import cz.iocb.elchem.lucene.SimilarStructureQuery;
import cz.iocb.elchem.molecule.AromaticityMode;
import cz.iocb.elchem.molecule.BinaryMolecule;
import cz.iocb.elchem.molecule.BinaryMoleculeBuilder;
import cz.iocb.elchem.molecule.MoleculeCreator;



public class SimilarityFingerprintFieldMapper extends FieldMapper
{
    public static final String CONTENT_TYPE = "similarity_fingerprint";
    public static final int maximumDepth = 3;


    public static class Defaults
    {
        public static final MappedFieldType FIELD_TYPE = new FieldType();

        static
        {
            FIELD_TYPE.freeze();
        }
    }


    public static class Builder extends FieldMapper.Builder<Builder, SimilarityFingerprintFieldMapper>
    {
        public Builder(String name)
        {
            super(name, Defaults.FIELD_TYPE, Defaults.FIELD_TYPE);
            builder = this;
        }


        @Override
        public SimilarityFingerprintFieldMapper build(BuilderContext context)
        {
            setupFieldType(context);
            return new SimilarityFingerprintFieldMapper(name, fieldType, defaultFieldType, context.indexSettings(),
                    multiFieldsBuilder.build(this, context), copyTo);
        }


        @Override
        protected void setupFieldType(BuilderContext context)
        {
            super.setupFieldType(context);

            fieldType.setIndexOptions(IndexOptions.NONE);
            fieldType.setHasDocValues(false);
            fieldType.setStored(false);

            defaultFieldType.setIndexOptions(IndexOptions.NONE);
            defaultFieldType.setHasDocValues(false);
            defaultFieldType.setStored(false);
        }
    }


    public static class TypeParser implements Mapper.TypeParser
    {
        @Override
        public Mapper.Builder<?, ?> parse(String name, Map<String, Object> node, ParserContext parserContext)
                throws MapperParsingException
        {
            Builder builder = new Builder(name);

            TypeParsers.parseField(builder, name, node, parserContext);

            return builder;
        }
    }


    public static class FieldType extends MappedFieldType
    {
        public FieldType()
        {
        }


        protected FieldType(FieldType ref)
        {
            super(ref);
        }


        @Override
        public String typeName()
        {
            return CONTENT_TYPE;
        }


        @Override
        public FieldType clone()
        {
            return new FieldType(this);
        }


        @Override
        public Query existsQuery(QueryShardContext arg0)
        {
            throw new UnsupportedOperationException("Cannot run exists query on [" + CONTENT_TYPE + "]");
        }


        @Override
        public Query termQuery(Object arg0, QueryShardContext arg1)
        {
            throw new IllegalArgumentException("Queries on [" + CONTENT_TYPE + "] fields are not supported");
        }
    }


    protected SimilarityFingerprintFieldMapper(String simpleName, MappedFieldType fieldType,
            MappedFieldType defaultFieldType, Settings indexSettings, MultiFields multiFields, CopyTo copyTo)
    {
        super(simpleName, fieldType, defaultFieldType, indexSettings, multiFields, copyTo);
    }


    @Override
    protected String contentType()
    {
        return CONTENT_TYPE;
    }


    @Override
    protected void parseCreateField(ParseContext context, List<IndexableField> fields) throws IOException
    {
        String sdf = context.externalValueSet() ? context.externalValue().toString() : context.parser().textOrNull();

        if(sdf == null)
            return;

        try
        {
            IAtomContainer container = MoleculeCreator.getMoleculeFromMolfile(sdf, AromaticityMode.AUTO);
            BinaryMoleculeBuilder builder = new BinaryMoleculeBuilder(container);
            byte[] binary = builder.asBytes(true);

            String name = fieldType().name();
            BinaryMolecule molecule = new BinaryMolecule(binary);
            List<List<Integer>> fp = IOCBFingerprint.getSimilarityFingerprint(molecule, maximumDepth);


            byte[] array = new byte[fp.stream().map(i -> i.size() + 1).reduce(0, Integer::sum) * Integer.BYTES];

            for(int pos = 0, i = 0; i < fp.size(); i++)
            {
                for(int b = 0; b < Integer.BYTES; b++)
                    array[pos++] = (byte) (fp.get(i).size() >> (8 * b));

                for(int bit : fp.get(i))
                    for(int b = 0; b < Integer.BYTES; b++)
                        array[pos++] = (byte) (bit >> (8 * b));
            }


            Set<Integer> bits = new HashSet<Integer>();

            for(List<Integer> seg : fp)
                bits.addAll(seg);


            fields.add(new TextField(name, new FingerprintTokenStream(bits)));
            fields.add(new StoredField(name, array));
            fields.add(new BinaryDocValuesField(name, new BytesRef(array)));

            for(int size = 0, i = 0; i < fp.size(); i++)
            {
                size += fp.get(i).size();
                fields.add(new IntPoint(name, size));
                size += SimilarStructureQuery.iterationSizeOffset;
            }
        }
        catch(CDKException e)
        {
            throw new IOException(e);
        }
    }
}
