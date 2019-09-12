package cz.iocb.elchem.lucene;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;



public class FingerprintTokenStream extends TokenStream
{
    private final CharTermAttribute charTermAttribute = addAttribute(CharTermAttribute.class);
    private final FingerprintBitMapping mapping = new FingerprintBitMapping();
    private Iterator<Integer> iterator;


    public FingerprintTokenStream(Collection<Integer> fp)
    {
        iterator = fp.iterator();
    }


    @Override
    public boolean incrementToken() throws IOException
    {
        charTermAttribute.setEmpty();

        if(!iterator.hasNext())
            return false;

        charTermAttribute.append(mapping.bitAsString(iterator.next()));

        return true;
    }
}
