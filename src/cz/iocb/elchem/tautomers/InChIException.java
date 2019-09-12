package cz.iocb.elchem.tautomers;

import org.openscience.cdk.exception.CDKException;



@SuppressWarnings("serial")
public class InChIException extends CDKException
{
    public InChIException(String message)
    {
        super(message);
    }
}
