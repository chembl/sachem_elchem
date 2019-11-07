package cz.iocb.elchem.fingerprint;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import cz.iocb.elchem.molecule.Molecule;



public class RCFingerprint extends Fingerprint
{
    private final static class AtomDesc
    {
        private final int hash;
        private final Set<Integer> cover;

        private AtomDesc(int hash)
        {
            this.hash = hash;
            this.cover = new HashSet<Integer>();
        }

        private AtomDesc(int hash, Set<Integer> cover)
        {
            this.hash = hash;
            this.cover = cover;
        }
    }


    public static final List<List<Integer>> getFingerprint(Molecule molecule, int minRadius, int maxRadius)
    {
        Map<Integer, AtomDesc> desc = new HashMap<Integer, AtomDesc>();
        List<List<Integer>> fp = new ArrayList<List<Integer>>(maxRadius - minRadius + 1);


        if(minRadius == 0)
            fp.add(new ArrayList<Integer>(molecule.getAtomCount()));

        for(int a = 0; a < molecule.getAtomCount(); a++)
        {
            if(molecule.getAtomNumber(a) == Molecule.AtomType.H /*|| molecule.isAtomPseudo(a)*/)
                continue;

            AtomDesc atomDesc = new AtomDesc(hashAtom(molecule, a));
            desc.put(a, atomDesc);

            if(minRadius == 0)
                fp.get(fp.size() - 1).add(atomDesc.hash);
        }


        for(int radius = 1; radius <= maxRadius; radius++)
        {
            Map<Integer, AtomDesc> newdesc = new HashMap<Integer, AtomDesc>();

            if(radius >= minRadius)
                fp.add(new ArrayList<Integer>(molecule.getAtomCount()));


            for(int i = 0; i < molecule.getAtomCount(); i++)
            {
                if(!desc.containsKey(i))
                    continue;

                List<Integer> hs = new ArrayList<Integer>(molecule.getAtomCount());
                Set<Integer> newcover = new HashSet<Integer>();

                for(int a : molecule.getBondedAtoms(i))
                {
                    int b = molecule.getBond(i, a);

                    if(!desc.containsKey(a) /*|| molecule.isQueryBond(b)*/)
                        continue;

                    newcover.add(b);
                    newcover.addAll(desc.get(a).cover);
                    hs.add(hash(hashBond(molecule, b), desc.get(a).hash));
                }


                if(!newcover.equals(desc.get(i).cover))
                {
                    AtomDesc atomDesc = new AtomDesc(hash(desc.get(i).hash, newcover.size(), hs), newcover);

                    newdesc.put(i, atomDesc);

                    if(radius >= minRadius)
                        fp.get(fp.size() - 1).add(atomDesc.hash);
                }
            }

            desc = newdesc;
        }


        for(List<Integer> list : fp)
            Collections.sort(list);


        return fp;
    }
}
