package cz.iocb.elchem.molecule;

import java.nio.ByteBuffer;
import java.security.AccessController;
import java.security.PrivilegedActionException;
import java.security.PrivilegedExceptionAction;
import net.sf.jnati.deploy.NativeLibraryLoader;



public class NativeIsomorphism
{
    @SuppressWarnings("serial")
    public static class IterationLimitExceededException extends Exception
    {
    }


    private final ByteBuffer implementation;


    public NativeIsomorphism(byte[] query, boolean[] restH, SearchMode searchMode, ChargeMode chargeMode,
            IsotopeMode isotopeMode, RadicalMode radicalMode, StereoMode stereoMode)
    {
        implementation = create(query, restH, searchMode.ordinal(), chargeMode.ordinal(), isotopeMode.ordinal(),
                radicalMode.ordinal(), stereoMode.ordinal());
    }


    public float match(byte[] query, long limit) throws IterationLimitExceededException
    {
        return match(implementation, query, limit);
    }


    private static native float match(ByteBuffer implementation, byte[] query, long limit)
            throws IterationLimitExceededException;


    private static native ByteBuffer create(byte[] query, boolean[] restH, int searchMode, int chargeMode,
            int isotopeMode, int radicalMode, int stereoMode);

    private static native void init();


    static
    {
        try
        {
            AccessController.doPrivileged((PrivilegedExceptionAction<Void>) () -> {
                NativeLibraryLoader.loadLibrary("isomorphism", "2.0");
                return null;
            });

            init();
        }
        catch(PrivilegedActionException e)
        {
            throw new RuntimeException(e);
        }
    }
}
