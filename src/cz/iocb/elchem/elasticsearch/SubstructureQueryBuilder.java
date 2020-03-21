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
import cz.iocb.elchem.lucene.SubstructureQuery;
import cz.iocb.elchem.molecule.AromaticityMode;
import cz.iocb.elchem.molecule.ChargeMode;
import cz.iocb.elchem.molecule.IsotopeMode;
import cz.iocb.elchem.molecule.QueryFormat;
import cz.iocb.elchem.molecule.SearchMode;
import cz.iocb.elchem.molecule.StereoMode;
import cz.iocb.elchem.molecule.TautomerMode;



public class SubstructureQueryBuilder extends AbstractQueryBuilder<SubstructureQueryBuilder>
{
    public static final String NAME = "match_substructures";

    public static final ParseField FIELD_FIELD = new ParseField("field");
    public static final ParseField MOLECULE_FIELD = new ParseField("molecule");
    public static final ParseField FORMAT_FIELD = new ParseField("format");
    public static final ParseField SEARCH_MODE_FIELD = new ParseField("search_mode");
    public static final ParseField CHARGE_MODE_FIELD = new ParseField("charge_mode");
    public static final ParseField ISOTOPE_MODE_FIELD = new ParseField("isotope_mode");
    public static final ParseField STEREO_MODE_FIELD = new ParseField("stereo_mode");
    public static final ParseField AROMATICITY_MODE_FIELD = new ParseField("aromaticity_mode");
    public static final ParseField TAUTOMER_MODE_FIELD = new ParseField("tautomer_mode");
    public static final ParseField MATCHING_LIMIT_FIELD = new ParseField("matching_limit");


    private String fieldName;
    private String molecule;
    private QueryFormat queryFormat;
    private SearchMode searchMode = SearchMode.SUBSTRUCTURE;
    private ChargeMode chargeMode = ChargeMode.DEFAULT_AS_ANY;
    private IsotopeMode isotopeMode = IsotopeMode.IGNORE;
    private StereoMode stereoMode = StereoMode.IGNORE;
    private AromaticityMode aromaticityMode = AromaticityMode.AUTO;
    private TautomerMode tautomerMode = TautomerMode.IGNORE;
    private long matchingLimit = 0;


    public SubstructureQueryBuilder()
    {
    }


    public SubstructureQueryBuilder(StreamInput in) throws IOException
    {
        super(in);

        fieldName = in.readString();
        molecule = in.readString();
        queryFormat = in.readEnum(QueryFormat.class);
        searchMode = in.readEnum(SearchMode.class);
        chargeMode = in.readEnum(ChargeMode.class);
        isotopeMode = in.readEnum(IsotopeMode.class);
        stereoMode = in.readEnum(StereoMode.class);
        aromaticityMode = in.readEnum(AromaticityMode.class);
        tautomerMode = in.readEnum(TautomerMode.class);
        matchingLimit = in.readLong();
    }


    @Override
    protected void doWriteTo(StreamOutput out) throws IOException
    {
        out.writeString(fieldName);
        out.writeString(molecule);
        out.writeEnum(queryFormat);
        out.writeEnum(searchMode);
        out.writeEnum(chargeMode);
        out.writeEnum(isotopeMode);
        out.writeEnum(stereoMode);
        out.writeEnum(aromaticityMode);
        out.writeEnum(tautomerMode);
        out.writeLong(matchingLimit);
    }


    public static SubstructureQueryBuilder fromXContent(XContentParser parser) throws IOException
    {
        String fieldPattern = null;
        String moleculePattern = null;
        QueryFormat queryFormatPattern = QueryFormat.UNSPECIFIED;
        SearchMode searchModePattern = SearchMode.SUBSTRUCTURE;
        ChargeMode chargeModePattern = ChargeMode.DEFAULT_AS_ANY;
        IsotopeMode isotopeModePattern = IsotopeMode.IGNORE;
        StereoMode stereoModePattern = StereoMode.IGNORE;
        AromaticityMode aromaticityModePattern = AromaticityMode.AUTO;
        TautomerMode tautomerModePattern = TautomerMode.IGNORE;
        long matchingLimitPattern = 0;

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
                else if(SEARCH_MODE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    searchModePattern = SearchMode.valueOf(value.toUpperCase());

                    if(searchModePattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown search mode [{}]", value);
                }
                else if(CHARGE_MODE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    chargeModePattern = ChargeMode.valueOf(value.toUpperCase());

                    if(chargeModePattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown charge mode [{}]", value);
                }
                else if(ISOTOPE_MODE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    isotopeModePattern = IsotopeMode.valueOf(value.toUpperCase());

                    if(isotopeModePattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown isotope mode [{}]", value);
                }
                else if(STEREO_MODE_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    String value = parser.text();
                    stereoModePattern = StereoMode.valueOf(value.toUpperCase());

                    if(stereoModePattern == null)
                        throw new ParsingException(parser.getTokenLocation(), "unknown stereo mode [{}]", value);
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
                else if(MATCHING_LIMIT_FIELD.match(currentFieldName, parser.getDeprecationHandler()))
                {
                    matchingLimitPattern = parser.longValue();

                    if(matchingLimitPattern < 0)
                        throw new ParsingException(parser.getTokenLocation(), "wrong matching limit value [{}]",
                                matchingLimitPattern);
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


        SubstructureQueryBuilder builder = new SubstructureQueryBuilder();
        builder.fieldName = fieldPattern;
        builder.molecule = moleculePattern;
        builder.queryFormat = queryFormatPattern;
        builder.searchMode = searchModePattern;
        builder.chargeMode = chargeModePattern;
        builder.isotopeMode = isotopeModePattern;
        builder.stereoMode = stereoModePattern;
        builder.aromaticityMode = aromaticityModePattern;
        builder.tautomerMode = tautomerModePattern;
        builder.matchingLimit = matchingLimitPattern;
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
        builder.field(SEARCH_MODE_FIELD.getPreferredName(), searchMode.name().toLowerCase());
        builder.field(CHARGE_MODE_FIELD.getPreferredName(), chargeMode.name().toLowerCase());
        builder.field(ISOTOPE_MODE_FIELD.getPreferredName(), isotopeMode.name().toLowerCase());
        builder.field(STEREO_MODE_FIELD.getPreferredName(), stereoMode.name().toLowerCase());
        builder.field(AROMATICITY_MODE_FIELD.getPreferredName(), aromaticityMode.name().toLowerCase());
        builder.field(TAUTOMER_MODE_FIELD.getPreferredName(), tautomerMode.name().toLowerCase());
        builder.field(MATCHING_LIMIT_FIELD.getPreferredName(), matchingLimit);
        printBoostAndQueryName(builder);
        builder.endObject();
    }


    @Override
    protected Query doToQuery(QueryShardContext arg0) throws IOException
    {
        try
        {
            return new SubstructureQuery(fieldName, molecule, queryFormat, searchMode, chargeMode, isotopeMode,
                    stereoMode, aromaticityMode, tautomerMode, matchingLimit);
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
    protected boolean doEquals(SubstructureQueryBuilder other)
    {
        return Objects.equals(fieldName, other.fieldName) && Objects.equals(molecule, other.molecule)
                && Objects.equals(queryFormat, other.queryFormat) && Objects.equals(searchMode, other.searchMode)
                && Objects.equals(chargeMode, other.chargeMode) && Objects.equals(isotopeMode, other.isotopeMode)
                && Objects.equals(stereoMode, other.stereoMode)
                && Objects.equals(aromaticityMode, other.aromaticityMode)
                && Objects.equals(tautomerMode, other.tautomerMode) && matchingLimit == other.matchingLimit;
    }


    @Override
    protected int doHashCode()
    {
        return Objects.hash(molecule);
    }
}
