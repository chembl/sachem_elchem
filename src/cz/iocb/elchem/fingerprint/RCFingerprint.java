package cz.iocb.elchem.fingerprint;

import java.util.ArrayList;
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
        private final int atom;
        private final Set<Integer> cover;

        private AtomDesc(int hash, int atom)
        {
            this.hash = hash;
            this.atom = atom;
            this.cover = new HashSet<Integer>();
        }

        private AtomDesc(int hash, int atom, Set<Integer> cover)
        {
            this.hash = hash;
            this.atom = atom;
            this.cover = cover;
        }
    }


    public static final Map<Integer, Integer> getFingerprint(Molecule molecule, int minRadius, int maxRadius,
            Map<Integer, Set<Integer>> info)
    {
        Map<Integer, AtomDesc> desc = new HashMap<Integer, AtomDesc>();
        List<AtomDesc> result = new ArrayList<AtomDesc>(molecule.getAtomCount());

        for(int a = 0; a < molecule.getAtomCount(); a++)
        {
            if(molecule.getAtomNumber(a) == Molecule.AtomType.H || isDummyAtom(molecule, a))
                continue;

            AtomDesc atomDesc = new AtomDesc(hashAtom(molecule, a), a);
            desc.put(a, atomDesc);

            if(minRadius == 0)
                result.add(atomDesc);
        }


        for(int radius = 1; radius <= maxRadius; radius++)
        {
            Map<Integer, AtomDesc> newdesc = new HashMap<Integer, AtomDesc>();
            int updates = 0;


            for(int i = 0; i < molecule.getAtomCount(); i++)
            {
                if(molecule.getAtomNumber(i) == Molecule.AtomType.H || !desc.containsKey(i))
                    continue;

                List<Integer> hs = new ArrayList<Integer>(molecule.getAtomCount());
                Set<Integer> newcover = new HashSet<Integer>();
                boolean acceptable = true;

                for(int a : molecule.getBondedAtoms(i))
                {
                    int b = molecule.getBond(i, a);

                    if(molecule.getAtomNumber(a) == Molecule.AtomType.H)
                        continue;

                    if(isDummyAtom(molecule, a) || isDummyBond(molecule, b) || !desc.containsKey(a))
                    {
                        acceptable = false;
                        break;
                    }

                    newcover.add(b);
                    newcover.addAll(desc.get(a).cover);
                    hs.add(hash(hashBond(molecule, b), desc.get(a).hash));
                }


                if(!acceptable || newcover.isEmpty())
                    continue;


                if(newcover.equals(desc.get(i).cover))
                {
                    newdesc.put(i, desc.get(i));
                }
                else
                {
                    AtomDesc atomDesc = new AtomDesc(hash(desc.get(i).hash, newcover.size(), hs), i, newcover);

                    newdesc.put(i, atomDesc);

                    if(radius >= minRadius)
                        result.add(atomDesc);

                    updates++;
                }
            }

            desc = newdesc;

            if(updates == 0)
                break;
        }


        Map<Integer, Integer> fp = new HashMap<Integer, Integer>();

        for(AtomDesc i : result)
        {
            setFp(fp, i.hash);

            if(info != null)
            {
                setInfo(info, i.hash, i.atom);

                for(int b : i.cover)
                {
                    setInfo(info, i.hash, molecule.getBondAtom(b, 0));
                    setInfo(info, i.hash, molecule.getBondAtom(b, 1));
                }
            }
        }

        return fp;
    }


    private static final boolean isDummyAtom(Molecule molecule, int atom)
    {
        return molecule.isAtomPseudo(atom);
    }


    private static final boolean isDummyBond(Molecule molecule, int bond)
    {
        byte type = molecule.getBondType(bond);

        return type == Molecule.BondType.NONE || type > Molecule.BondType.AROMATIC;
    }
}
