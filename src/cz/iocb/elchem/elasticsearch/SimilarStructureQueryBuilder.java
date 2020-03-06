package cz.iocb.elchem.elasticsearch;

import java.io.IOException;
import java.util.Objects;
import java.util.concurrent.TimeoutException;
import org.apache.lucene.search.Query;
import org.elasticsearch.common.ParseField;
import org.elasticsearch.common.ParsingException;
import org.elasticsearch.common.io.stream.StreamInput;
import org.elasticsearch.common.io.stream.StreamOutput;
import org.elasticsearch.common.xcontent.XContentBuilder;
import org.elasticsearch.common.xcontent.XContentParser;
import org.elasticsearch.index.query.AbstractQueryBuilder;
import org.elasticsearch.index.query.ExistsQueryBuilder;
import org.elasticsearch.index.query.QueryShardContext;
import org.openscience.cdk.exception.CDKException;
import cz.iocb.elchem.lucene.SimilarStructureQuery;
import cz.iocb.elchem.molecule.AromaticityMode;
import cz.iocb.elchem.molecule.QueryFormat;
import cz.iocb.elchem.molecule.TautomerMode;



public class SimilarStructureQueryBuilder extends AbstractQueryBuilder<SimilarStructureQueryBuilder>
{
    public static final String NAME = "match_similar_structures";

    public static final ParseField FIELD_FIELD = new ParseField("field");
    public static final ParseField MOLECULE_FIELD = new ParseField("molecule");
    public static final ParseField FORMAT_FIELD = new ParseField("format");
    public static final ParseField THRESHOLD_FIELD = new ParseField("threshold");
    public static final ParseField MAXIMUM_DEPTH_FIELD = new ParseField("maximum_depth");
    public static final ParseField AROMATICITY_MODE_FIELD = new ParseField("aromaticity_mode");
    public static final ParseField TAUTOMER_MODE_FIELD = new ParseField("tautomer_mode");


    private String fieldName;
    private String molecule;
    private QueryFormat queryFormat;
    private float threshold = 0.8f;
    private int maximumDepth = 1;
    private AromaticityMode aromaticityMode = AromaticityMode.AUTO;
    private TautomerMode tautomerMode = TautomerMode.IGNORE;


    public SimilarStructureQueryBuilder()
    {
    }


    public SimilarStructureQueryBuilder(StreamInput in) throws IOException
    {
        super(in);

        fieldName = in.readString();
        molecule = in.readString();
        queryFormat = in.readEnum(QueryFormat.class);
        threshold = in.readFloat();
        maximumDepth = in.readInt();
        aromaticityMode = in.readEnum(AromaticityMode.class);
        tautomerMode = in.readEnum(TautomerMode.class);
    }


    @Override
    protected void doWriteTo(StreamOutput out) throws IOException
    {
        out.writeString(fieldName);
        out.writeString(molecule);
        out.writeEnum(queryFormat);
        out.writeFloat(threshold);
        out.writeInt(maximumDepth);
        out.writeEnum(aromaticityMode);
        out.writeEnum(tautomerMode);
    }


    public static SimilarStructureQueryBuilder fromXContent(XContentParser parser) throws IOException
    {
        String fieldPattern = null;
        String moleculePattern = null;
        QueryFormat queryFormatPattern = QueryFormat.UNSPECIFIED;
        float thresholdPattern = 0.8f;
        int maximumDepthPattern = 1;
        AromaticityMode aromaticityModePattern = AromaticityMode.AUTO;
        TautomerMode tautomerModePattern = TautomerMode.IGNORE;

        String queryName = null;
        float boost = AbstractQueryBuilder.DEFAULT_BOOST;

        XContentParser.Token token;
        String currentFieldName = null;

        while((token = parser.nextToken()) != XContentParser.Token.END_OBJECT)
        {
            if(token == XContentParser.Token.FIELD_NAME)
            {
                currentFieldName = parser.currentName();
            }
            else if(token.isValue())
            {
                if(FIELD_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    fieldPattern = parser.text();
                }
                else if(MOLECULE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    moleculePattern = parser.text();
                }
                else if(FORMAT_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    queryFormatPattern = QueryFormat.valueOf(value.toUpperCase());

                    if(queryFormatPattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown query format [{}]", value);
                }
                else if(THRESHOLD_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    thresholdPattern = parser.floatValue();

                    if(thresholdPattern < 0.5f || thresholdPattern > 1.0f)
                        throw new ParsingException(parser.getTokenLocation(), "wrong threshold value [{}]",
                                thresholdPattern);
                }
                else if(MAXIMUM_DEPTH_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    maximumDepthPattern = parser.intValue();

                    if(maximumDepthPattern < 0 || maximumDepthPattern > SimilarityFingerprintFieldMapper.maximumDepth)
                        throw new ParsingException(parser.getTokenLocation(), "wrong maximum depth value [{}]",
                                maximumDepthPattern);
                }
                else if(AROMATICITY_MODE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    aromaticityModePattern = AromaticityMode.valueOf(value.toUpperCase());

                    if(aromaticityModePattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown aromaticity mode [{}]", value);
                }
                else if(TAUTOMER_MODE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    tautomerModePattern = TautomerMode.valueOf(value.toUpperCase());

                    if(tautomerModePattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown tautomer mode [{}]", value);
                }
                else if(AbstractQueryBuilder.NAME_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    queryName = parser.text();
                }
                else if(AbstractQueryBuilder.BOOST_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    boost = parser.floatValue();
                }
                else
                {
                    throw new ParsingException(parser.getTokenLocation(),
                            "[" + ExistsQueryBuilder.NAME + "] query does not support [" + currentFieldName + "]");
                }
            }
            else
            {
                throw new ParsingException(parser.getTokenLocation(), "[" + ExistsQueryBuilder.NAME
                        + "] unknown token [" + token + "] after [" + currentFieldName + "]");
            }
        }

        if(fieldPattern == null)
        {
            throw new ParsingException(parser.getTokenLocation(),
                    "[" + ExistsQueryBuilder.NAME + "] must be provided with a [field]");
        }

        if(moleculePattern == null)
        {
            throw new ParsingException(parser.getTokenLocation(),
                    "[" + ExistsQueryBuilder.NAME + "] must be provided with a [molecule]");
        }


        SimilarStructureQueryBuilder builder = new SimilarStructureQueryBuilder();
        builder.fieldName = fieldPattern;
        builder.molecule = moleculePattern;
        builder.queryFormat = queryFormatPattern;
        builder.threshold = thresholdPattern;
        builder.maximumDepth = maximumDepthPattern;
        builder.aromaticityMode = aromaticityModePattern;
        builder.tautomerMode = tautomerModePattern;
        builder.queryName(queryName);
        builder.boost(boost);
        return builder;
    }


    @Override
    protected void doXContent(XContentBuilder builder, Params params) throws IOException
    {
        builder.startObject(NAME);
        builder.field(FIELD_FIELD.getPreferredName(), fieldName);
        builder.field(MOLECULE_FIELD.getPreferredName(), molecule);
        builder.field(FORMAT_FIELD.getPreferredName(), queryFormat.name().toLowerCase());
        builder.field(THRESHOLD_FIELD.getPreferredName(), threshold);
        builder.field(MAXIMUM_DEPTH_FIELD.getPreferredName(), maximumDepth);
        builder.field(AROMATICITY_MODE_FIELD.getPreferredName(), aromaticityMode.name().toLowerCase());
        builder.field(TAUTOMER_MODE_FIELD.getPreferredName(), tautomerMode.name().toLowerCase());
        printBoostAndQueryName(builder);
        builder.endObject();
    }


    @Override
    protected Query doToQuery(QueryShardContext arg0) throws IOException
    {
        try
        {
            return new SimilarStructureQuery(fieldName, molecule, queryFormat, threshold, maximumDepth, aromaticityMode,
                    tautomerMode);
        }
        catch(CDKException | TimeoutException e)
        {
            throw new IOException(e);
        }
    }


    @Override
    public String getWriteableName()
    {
        return NAME;
    }


    @Override
    protected boolean doEquals(SimilarStructureQueryBuilder other)
    {
        return Objects.equals(fieldName, other.fieldName) && Objects.equals(molecule, other.molecule);
    }


    @Override
    protected int doHashCode()
    {
        return Objects.hash(molecule);
    }
}
