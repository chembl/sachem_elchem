package cz.iocb.elchem.lucene;



public final class FingerprintBitMapping
{
    private static final char[] b64str = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
            'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4',
            '5', '6', '7', '8', '9', '+', '/' };

    private final char[] buffer = new char[6];


    public final String bitAsString(int bit)
    {
        for(int i = 0; i < 6; i++)
            buffer[i] = b64str[bit >>> 6 * i & 0x3f];

        return new String(buffer);
    }
}
