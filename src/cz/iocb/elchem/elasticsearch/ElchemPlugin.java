package cz.iocb.elchem.elasticsearch;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.elasticsearch.index.mapper.Mapper;
import org.elasticsearch.plugins.MapperPlugin;
import org.elasticsearch.plugins.Plugin;
import org.elasticsearch.plugins.SearchPlugin;
import com.google.common.collect.ImmutableMap;



public class ElchemPlugin extends Plugin implements MapperPlugin, SearchPlugin
{
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
