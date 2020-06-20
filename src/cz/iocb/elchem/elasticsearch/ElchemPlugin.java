package cz.iocb.elchem.elasticsearch;

import java.security.AccessController;
import java.security.PrivilegedActionException;
import java.security.PrivilegedExceptionAction;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.elasticsearch.index.mapper.Mapper;
import org.elasticsearch.plugins.MapperPlugin;
import org.elasticsearch.plugins.Plugin;
import org.elasticsearch.plugins.SearchPlugin;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import com.google.common.collect.ImmutableMap;
import net.sf.jnati.deploy.NativeLibraryLoader;



public class ElchemPlugin extends Plugin implements MapperPlugin, SearchPlugin
{
    static
    {
        try
        {
            AccessController.doPrivileged((PrivilegedExceptionAction<Void>) () -> {
                SilentChemObjectBuilder.getInstance();
                NativeLibraryLoader.loadLibrary("elchem", "2.5.0");
                return null;
            });
        }
        catch(PrivilegedActionException e)
        {
            throw new RuntimeException(e);
        }
    }


    @Override
    public Map<String, Mapper.TypeParser> getMappers()
    {
        return ImmutableMap.<String, Mapper.TypeParser> builder()
                .put(StructureFingerprintFieldMapper.CONTENT_TYPE, new StructureFingerprintFieldMapper.TypeParser())
                .put(SimilarityFingerprintFieldMapper.CONTENT_TYPE, new SimilarityFingerprintFieldMapper.TypeParser())
                .build();
    }


    @Override
    public List<QuerySpec<?>> getQueries()
    {
        return Arrays.asList(
                new QuerySpec<SubstructureQueryBuilder>(SubstructureQueryBuilder.NAME, SubstructureQueryBuilder::new,
                        SubstructureQueryBuilder::fromXContent),
                new QuerySpec<SimilarStructureQueryBuilder>(SimilarStructureQueryBuilder.NAME,
                        SimilarStructureQueryBuilder::new, SimilarStructureQueryBuilder::fromXContent));
    }
}
