/*
 * Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
 * Copyright (C) 2016-2017 Miroslav Kratochvil
 * Copyright (C) 2016-2019 Jakub Galgonek
 *
 * Subgraph functions of this file are part of the RDKit.
 *
 * The contents are covered by the terms of the BSD license which is included in the file license.txt, found at the root
 * of the RDKit source tree.
 */
package cz.iocb.elchem.fingerprint;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import cz.iocb.elchem.molecule.Molecule;



public class SGFingerprint extends Fingerprint
{
    private static final class AtomDesc
    {
        private final int hash;
        private final List<Integer> cover;

        private AtomDesc(int hash, ArrayList<Integer> cover)
        {
            this.hash = hash;
            this.cover = cover;
        }
    }


    private static final Integer hashSubmolecule(Collection<Integer> bondIds, Molecule molecule)
    {
        Map<Integer, List<Integer>> preatoms = new HashMap<Integer, List<Integer>>();

        int max = 1;

        for(int bid : bondIds)
        {
            for(int i = 0; i < 2; i++)
            {
                List<Integer> list = preatoms.get(molecule.getBondAtom(bid, i));

                if(list == null)
                {
                    list = new ArrayList<>();
                    preatoms.put(molecule.getBondAtom(bid, i), list);
                }

                list.add(bid);

                if(list.size() > max)
                    max = list.size();
            }
        }


        List<Map<Integer, AtomDesc>> atoms = new ArrayList<Map<Integer, AtomDesc>>(max);
        Map<Integer, Map<Integer, Integer>> bonds = new HashMap<Integer, Map<Integer, Integer>>();


        for(int i = 0; i <= max; i++)
            atoms.add(i, new HashMap<Integer, AtomDesc>());


        for(Entry<Integer, List<Integer>> a : preatoms.entrySet())
        {
            int aid = a.getKey();

            Map<Integer, Integer> map = new HashMap<Integer, Integer>();
            bonds.put(aid, map);

            for(int bid : a.getValue())
                map.put(molecule.getOtherBondAtom(bid, aid), hashBond(molecule, bid));

            atoms.get(a.getValue().size()).put(aid, new AtomDesc(hashAtom(molecule, aid), new ArrayList<Integer>()));
        }


        // purge the leaves until there is nothing left
        while(!atoms.get(1).isEmpty())
        {
            Map<Integer, AtomDesc> ats = new HashMap<Integer, AtomDesc>();
            ats.putAll(atoms.get(1));

            for(Entry<Integer, AtomDesc> a : ats.entrySet())
            {
                int aid = a.getKey();

                // there is just one thing in bonds.get(aid)
                int addToID = bonds.get(aid).keySet().iterator().next();
                int bondHash = bonds.get(aid).values().iterator().next();

                int resultHash = hash(ats.get(aid).hash, bondHash, ats.get(aid).cover);

                bonds.remove(aid);
                atoms.get(1).remove(aid);

                int addToDeg = bonds.get(addToID).size();

                if(ats.containsKey(addToID))
                {
                    // final doublet handling!
                    AtomDesc other = atoms.get(addToDeg).get(addToID); // clone?

                    int otherAtomHash = hash(other.hash, bondHash, other.cover);
                    atoms.get(addToDeg).remove(addToID);
                    bonds.get(addToID).remove(aid);

                    if(atoms.get(0).get(addToID) != null)
                        throw new RuntimeException();

                    ArrayList<Integer> list = new ArrayList<Integer>();
                    list.add(resultHash);
                    list.add(otherAtomHash);
                    atoms.get(0).put(addToID, new AtomDesc(666, list));

                    break;
                }

                // "normal" leaf
                AtomDesc resAtom = atoms.get(addToDeg).get(addToID); // clone?
                atoms.get(addToDeg).remove(addToID);
                bonds.get(addToID).remove(aid);
                resAtom.cover.add(resultHash);
                atoms.get(addToDeg - 1).put(addToID, resAtom);
            }
        }

        if(atoms.get(0).isEmpty())
        {
            // there must be a single cycle!
            for(int i = 0; i < atoms.size(); i++)
                if(i != 2 && atoms.get(i).size() != 0)
                    return null;

            // (note that the graph is connected)

            Map<Integer, AtomDesc> ats = atoms.get(2);

            if(ats.isEmpty())
                return null; // this would be just weird.

            int curId = ats.keySet().iterator().next();
            int lastId = -1;
            int startId = curId;
            boolean foundNext = true;
            List<Integer> cycle = new ArrayList<Integer>(2 * ats.size());

            AtomDesc firstAtom = ats.values().iterator().next();
            cycle.add(hash(firstAtom.hash, 0, firstAtom.cover));

            while(foundNext)
            {
                foundNext = false;

                for(Entry<Integer, Integer> i : bonds.get(curId).entrySet())
                {
                    if(i.getKey() != lastId)
                    {
                        cycle.add(i.getValue());

                        if(i.getKey() == startId)
                            break;

                        AtomDesc theAtom = ats.get(i.getKey());
                        cycle.add(hash(theAtom.hash, 0, theAtom.cover));
                        lastId = curId;
                        curId = i.getKey();
                        foundNext = true;
                        break;
                    }
                }
            }


            int minrot = 0;
            int mindir = -1;
            int n = cycle.size();

            for(int rot = 0; rot < n; rot++)
            {
                for(int dir = -1; dir <= 1; dir += 2)
                {
                    for(int i = 0; i < n; i++)
                    {
                        int cmp = Integer.compareUnsigned(cycle.get((n + minrot + i * mindir) % n),
                                cycle.get((n + rot + i * dir) % n));

                        if(cmp < 0)
                            break; // must be ok

                        if(cmp == 0)
                            continue; // ok

                        // found better!
                        minrot = rot;
                        mindir = dir;
                        break;
                    }
                }
            }


            int seed = 0;

            for(int i = 0; i < n; i++)
                seed = updateSeed(cycle.get((n + minrot + i * mindir) % n), seed);

            return seed;
        }
        else
        {
            AtomDesc a = atoms.get(0).values().iterator().next();
            return hash(a.hash, 0, a.cover);
        }
    }


    private static final List<List<Integer>> getNeighborLists(Molecule molecule)
    {
        List<List<Integer>> nbrs = new ArrayList<List<Integer>>(molecule.getBondCount());

        for(int i = 0; i < molecule.getBondCount(); i++)
            nbrs.add(new ArrayList<Integer>());


        // create a list of neighbors for each bond
        for(int i = 0; i < molecule.getAtomCount(); i++)
        {
            if(molecule.getAtomNumber(i) <= Molecule.AtomType.H)
                continue;

            int[] bondedAtoms = molecule.getBondedAtoms(i);

            for(int l : bondedAtoms)
            {
                if(molecule.getAtomNumber(l) <= Molecule.AtomType.H)
                    continue;

                int bid1 = molecule.getBond(i, l);

                if(molecule.getBondType(bid1) > Molecule.BondType.AROMATIC)
                    continue;

                for(int k : bondedAtoms)
                {
                    if(molecule.getAtomNumber(k) <= Molecule.AtomType.H)
                        continue;

                    int bid2 = molecule.getBond(i, k);

                    if(molecule.getBondType(bid2) > Molecule.BondType.AROMATIC)
                        continue;

                    if(bid1 != bid2)
                        nbrs.get(bid1).add(bid2);
                }
            }
        }

        return nbrs;
    }


    private static final void recurseWalkRange(List<List<Integer>> nbrs, List<Integer> spath, ArrayDeque<Integer> cands,
            int lowerLen, int upperLen, boolean[] forbidden, List<List<Integer>> res)
    {
        int nsize = spath.size();

        if(nsize >= lowerLen && nsize <= upperLen)
            res.add(spath);

        // end case for recursion
        if(nsize >= upperLen)
            return;


        // we  have the candidates that can be used to add to the existing path try extending the subgraphs
        while(cands.size() != 0)
        {
            int next = cands.removeLast(); // start with the last one in the candidate list

            if(!forbidden[next])
            {
                // this bond should not appear in the later subgraphs
                forbidden[next] = true;

                // update a local stack before the next recursive call
                ArrayDeque<Integer> tstack = new ArrayDeque<Integer>(cands);

                for(int bid : nbrs.get(next))
                    if(!forbidden[bid])
                        tstack.add(bid);

                ArrayList<Integer> tpath = new ArrayList<Integer>(spath);
                tpath.add(next);

                recurseWalkRange(nbrs, tpath, tstack, lowerLen, upperLen, forbidden.clone(), res);
            }
        }
    }


    private static final List<List<Integer>> findSubgraphs(Molecule molecule, int lowerLen, int upperLen,
            int rootedAtAtom)
    {
        boolean[] forbidden = new boolean[molecule.getBondCount()];

        List<List<Integer>> nbrs = getNeighborLists(molecule);

        // start path at each bond
        List<List<Integer>> res = new ArrayList<List<Integer>>();

        // start paths at each bond:
        for(int i = 0; i < molecule.getBondCount(); i++)
        {
            if(molecule.getBondType(i) > Molecule.BondType.AROMATIC)
                continue;

            if(molecule.getAtomNumber(molecule.getBondAtom(i, 0)) <= Molecule.AtomType.H)
                continue;

            if(molecule.getAtomNumber(molecule.getBondAtom(i, 1)) <= Molecule.AtomType.H)
                continue;


            // if we are only returning paths rooted at a particular atom, check now that this bond involves that atom:
            if(rootedAtAtom >= 0 && molecule.getBondAtom(i, 0) != rootedAtAtom
                    && molecule.getBondAtom(i, 1) != rootedAtAtom)
                continue;

            // do not come back to this bond in the later subgraphs
            if(forbidden[i])
                continue;

            forbidden[i] = true;

            // start the recursive path building with the current bond
            List<Integer> spath = new ArrayList<Integer>();
            spath.add(i);

            // neighbors of this bond are the next candidates
            ArrayDeque<Integer> cands = new ArrayDeque<Integer>(nbrs.get(i));

            // now call the recursive function little bit different from the python version
            // the result list of paths is passed as a reference, instead of on the fly appending
            recurseWalkRange(nbrs, spath, cands, lowerLen, upperLen, forbidden.clone(), res);
        }

        return res;
    }


    private static final void addMoleculeFingerprint(Molecule molecule, Map<Integer, Integer> fp, int minLen,
            int maxLen, Map<Integer, Set<Integer>> info)
    {
        List<List<Integer>> allSGs = findSubgraphs(molecule, minLen, maxLen, -1);


        for(List<Integer> i : allSGs)
        {
            Integer hash = hashSubmolecule(i, molecule);

            if(hash != null)
            {
                setFp(fp, hash);

                if(info != null)
                {
                    for(int b : i)
                    {
                        setInfo(info, hash, molecule.getBondAtom(b, 0));
                        setInfo(info, hash, molecule.getBondAtom(b, 1));
                    }
                }
            }
        }
    }


    private static final int fragmentWalk(Molecule molecule, int atom, boolean[] visitedAtoms, boolean[] visitedBonds,
            int visited)
    {
        visitedAtoms[atom] = true;

        if(molecule.getAtomNumber(atom) <= Molecule.AtomType.H)
            return visited;

        for(int other : molecule.getBondedAtoms(atom))
        {
            if(molecule.getAtomNumber(other) <= Molecule.AtomType.H)
                continue;

            int bond = molecule.getBond(atom, other);

            if(molecule.getBondType(bond) > Molecule.BondType.AROMATIC)
                continue;

            if(!visitedBonds[bond])
            {
                visitedBonds[bond] = true;
                visited++;
            }

            if(!visitedAtoms[other])
                visited += fragmentWalk(molecule, other, visitedAtoms, visitedBonds, visited);
        }

        return visited;
    }


    public static final Map<Integer, Integer> getFingerprint(Molecule molecule, int minLen, int maxLen,
            boolean forQuery, Map<Integer, Set<Integer>> info)
    {
        if(forQuery)
        {
            int minQueryLen = maxLen;

            boolean[] visitedAtoms = new boolean[molecule.getAtomCount()];
            boolean[] visitedBonds = new boolean[molecule.getBondCount()];

            for(int a = 0; a < molecule.getAtomCount(); a++)
            {
                if(!visitedAtoms[a])
                {
                    int visited = fragmentWalk(molecule, a, visitedAtoms, visitedBonds, 0);

                    if(visited > 0 && visited < minQueryLen && visited >= minLen)
                        minQueryLen = visited;
                }
            }

            minLen = minQueryLen;
        }


        Map<Integer, Integer> fp = new HashMap<Integer, Integer>();
        addMoleculeFingerprint(molecule, fp, minLen, maxLen, info);

        return fp;
    }
}
