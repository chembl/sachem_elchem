package cz.iocb.elchem.fingerprint;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import cz.iocb.elchem.molecule.BinaryMolecule;
import cz.iocb.elchem.molecule.Isomorphism;
import cz.iocb.elchem.molecule.Molecule;



public class CRNGFingerprint extends Fingerprint
{
    private static final Isomorphism[] patterns;


    static
    {
        try(ObjectInputStream in = new ObjectInputStream(CRNGFingerprint.class.getResourceAsStream("/patterns.bin")))
        {
            byte[][] molecules = (byte[][]) in.readObject();
            patterns = new Isomorphism[molecules.length];

            for(int i = 0; i < molecules.length; i++)
                patterns[i] = new Isomorphism(new BinaryMolecule(molecules[i]));
        }
        catch(IOException | ClassNotFoundException e)
        {
            throw new RuntimeException(e);
        }
    }


    public static final Map<Integer, Integer> getFingerprint(Molecule molecule, Map<Integer, Set<Integer>> info)
    {
        Map<Integer, Integer> fp = new HashMap<Integer, Integer>();

        for(int i = 0; i < patterns.length; i++)
        {
            List<int[]> matches = patterns[i].matchAll(molecule, 256);

            if(matches.size() == 0)
                continue;

            setFp(fp, i, matches.size());

            if(info != null)
            {
                for(int[] match : matches)
                    for(int atom : match)
                        setInfo(info, i, atom);
            }
        }

        return fp;
    }
}
