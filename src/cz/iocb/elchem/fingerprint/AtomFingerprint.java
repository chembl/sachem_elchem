package cz.iocb.elchem.fingerprint;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import cz.iocb.elchem.molecule.Molecule;



public class AtomFingerprint extends Fingerprint
{
    public static final Map<Integer, Integer> getFingerprint(Molecule molecule, Map<Integer, Set<Integer>> info)
    {
        Map<Integer, Integer> fp = new HashMap<Integer, Integer>();

        for(int i = 0; i < molecule.getAtomCount(); i++)
        {
            int hsh = molecule.getAtomNumber(i);

            if(hsh <= Molecule.AtomType.H)
                continue;

            setFp(fp, hsh);

            if(info != null)
                setInfo(info, hsh, i);
        }

        return fp;
    }
}
