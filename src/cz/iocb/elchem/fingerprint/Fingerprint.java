package cz.iocb.elchem.fingerprint;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import cz.iocb.elchem.molecule.Molecule;



public abstract class Fingerprint
{
    protected static final void setFp(Map<Integer, Integer> fp, int hash)
    {
        Integer value = fp.get(hash);

        if(value == null)
            fp.put(hash, 1);
        else
            fp.put(hash, value + 1);
    }


    protected static final void setFp(Map<Integer, Integer> fp, int hash, int size)
    {
        Integer value = fp.get(hash);

        if(value == null)
            fp.put(hash, size);
        else
            fp.put(hash, value + size);
    }


    protected static final void setInfo(Map<Integer, Set<Integer>> info, int hash, int atom)
    {
        Set<Integer> set = info.get(hash);

        if(set == null)
        {
            set = new HashSet<Integer>();
            info.put(hash, set);
        }

        set.add(atom);
    }


    protected static final void setInfo(Map<Integer, Set<Integer>> info, int hash, Set<Integer> atoms)
    {
        Set<Integer> set = info.get(hash);

        if(set == null)
        {
            set = new HashSet<Integer>();
            info.put(hash, set);
        }

        set.addAll(atoms);
    }


    protected static final int hashAtom(Molecule molecule, int atom)
    {
        return molecule.getAtomNumber(atom);
    }


    protected static final int hashBond(Molecule molecule, int bond)
    {
        return molecule.getBondType(bond);
    }


    protected static final int hash(int a, int b)
    {
        int seed = 0;
        seed = updateSeed(a, seed);
        seed = updateSeed(b, seed);
        return seed;
    }


    protected static final int hash(int a, int b, int c)
    {
        int seed = 0;
        seed = updateSeed(a, seed);
        seed = updateSeed(b, seed);
        seed = updateSeed(c, seed);
        return seed;
    }


    protected static final int hash(int a, int b, List<Integer> l)
    {
        Collections.sort(l, Integer::compareUnsigned);

        int seed = 0;
        seed = updateSeed(a, seed);
        seed = updateSeed(b, seed);

        for(int i : l)
            seed = updateSeed(i, seed);

        return seed;
    }


    protected static final int updateSeed(int x, int seed)
    {
        long xl = Integer.toUnsignedLong(x);
        long seedl = Integer.toUnsignedLong(seed);

        seed ^= (int) (xl * 2654435761l + 2654435769l + (seedl << 6) + (seedl >> 2));
        return seed;
    }
}
