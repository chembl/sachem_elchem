package cz.iocb.elchem.fingerprint;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import cz.iocb.elchem.molecule.BinaryMolecule;
import cz.iocb.elchem.molecule.Molecule;



public class IOCBFingerprint extends Fingerprint
{
    private static final void processElements(Map<Integer, Integer> var, Map<Integer, Set<Integer>> vari, int n,
            Set<Integer> fp, int maxFeatLogCount, boolean forQuery, Map<Integer, Set<Integer>> info)
    {
        for(Entry<Integer, Integer> i : var.entrySet())
        {
            int h = i.getKey();
            int cnt = i.getValue();

            if(forQuery)
            {
                int lc = 0;

                for(int c = 0; cnt != 0 && c < maxFeatLogCount; c++, cnt /= 2)
                    lc = c;

                int hsh = hash(n, lc, h);
                fp.add(hsh);

                if(info != null)
                    setInfo(info, hsh, vari.get(h));
            }
            else
            {
                for(int c = 0; cnt != 0 && c < maxFeatLogCount; c++, cnt /= 2)
                {
                    int hsh = hash(n, c, h);

                    fp.add(hsh);

                    if(info != null)
                        setInfo(info, hsh, vari.get(h));
                }
            }
        }
    }


    public static final Set<Integer> getSubstructureFingerprint(Molecule molecule, int graphSize, int maxFeatLogCount,
            boolean forQuery, Map<Integer, Set<Integer>> info)
    {
        Set<Integer> fp = new HashSet<Integer>();

        Map<Integer, Set<Integer>> sgi = info != null ? new HashMap<Integer, Set<Integer>>() : null;
        Map<Integer, Integer> sg = SGFingerprint.getFingerprint(molecule, 0, graphSize, forQuery, sgi);
        processElements(sg, sgi, 1, fp, maxFeatLogCount, forQuery, info);

        Map<Integer, Set<Integer>> crngi = info != null ? new HashMap<Integer, Set<Integer>>() : null;
        Map<Integer, Integer> crng = CRNGFingerprint.getFingerprint(molecule, crngi);
        processElements(crng, crngi, 2, fp, maxFeatLogCount, forQuery, info);

        Map<Integer, Set<Integer>> atomi = info != null ? new HashMap<Integer, Set<Integer>>() : null;
        Map<Integer, Integer> atom = AtomFingerprint.getFingerprint(molecule, atomi);
        processElements(atom, atomi, 3, fp, maxFeatLogCount, forQuery, info);

        return fp;
    }


    public static Set<Integer> getSubstructureFingerprint(BinaryMolecule molecule)
    {
        return getSubstructureFingerprint(molecule, 7, 5, false, null);
    }


    public static Set<Integer> getSubstructureFingerprint(BinaryMolecule molecule, Map<Integer, Set<Integer>> info)
    {
        return getSubstructureFingerprint(molecule, 7, 5, true, info);
    }


    public static final List<List<Integer>> getSimilarityFingerprint(Molecule molecule, int circSize)
    {
        return RCFingerprint.getFingerprint(molecule, 0, circSize);
    }
}
