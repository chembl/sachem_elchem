/*
 * Copyright (C) 2015-2020 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2011-2011 Mark Rijnbeek    markr@ebi.ac.uk
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 */
package cz.iocb.elchem.molecule;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import cz.iocb.elchem.molecule.InChITools.TautomericGroup;



public class InChITautomerGenerator
{
    @SuppressWarnings("serial")
    public static class InChITautomerException extends CDKException
    {
        public InChITautomerException(String message)
        {
            super(message);
        }
    }


    private final IAtomContainer skeleton;
    private final List<TautomericGroup> groups;
    private final LinkedHashSet<Integer> selectedAtoms;
    private final LinkedHashSet<Integer> selectedBonds;
    private final LinkedHashMap<Integer, List<Integer>> atomAlternBonds;

    private final int[] bondAtoms;
    private final int[] groupMobiles;
    private final int[] groupCharges;
    private final int[] losts;
    private final int[] opens;
    private final int[] mobiles;
    private final int[] charges;
    private final int[] orders;
    private final boolean[] aromatics;

    private int combination;

    private int queueIdx;
    private int determinedBondIdx;
    private final int[] determined;
    private final int[] queue;
    private final boolean[] queuing;

    private final List<IAtomContainer> tautomers = new ArrayList<IAtomContainer>();
    private final List<int[]> filters = new ArrayList<int[]>();


    public static List<IAtomContainer> generate(IAtomContainer molecule, List<Integer> alternBonds,
            List<TautomericGroup> groups) throws CDKException
    {
        InChITautomerGenerator generator = new InChITautomerGenerator(molecule, alternBonds, groups);

        return generator.tautomers;
    }


    private InChITautomerGenerator(IAtomContainer molecule, List<Integer> alternBonds, List<TautomericGroup> groups)
            throws CDKException
    {
        this.groups = groups;

        skeleton = clone(molecule);
        bondAtoms = new int[2 * molecule.getBondCount()];
        losts = new int[molecule.getAtomCount()];
        opens = new int[molecule.getAtomCount()];
        mobiles = new int[molecule.getAtomCount()];
        charges = new int[molecule.getAtomCount()];
        orders = new int[molecule.getBondCount()];
        aromatics = new boolean[molecule.getBondCount()];
        groupMobiles = new int[groups.size()];
        groupCharges = new int[groups.size()];

        atomAlternBonds = new LinkedHashMap<Integer, List<Integer>>(molecule.getAtomCount());
        int[] singleBondCounts = new int[molecule.getAtomCount()];

        for(int bond : alternBonds)
        {
            IBond moleculeBond = molecule.getBond(bond);
            aromatics[bond] = moleculeBond.isAromatic();

            for(int i = 0; i < 2; i++)
            {
                int atom = molecule.indexOf(moleculeBond.getAtom(i));
                List<Integer> bondList = atomAlternBonds.get(atom);

                if(bondList == null)
                    atomAlternBonds.put(atom, bondList = new ArrayList<Integer>(4));

                bondList.add(bond);
                bondAtoms[2 * bond + i] = atom;

                if(moleculeBond.getOrder() == Order.SINGLE)
                    singleBondCounts[atom]++;
            }
        }


        selectedAtoms = new LinkedHashSet<Integer>(molecule.getAtomCount());
        selectedBonds = new LinkedHashSet<Integer>(molecule.getBondCount());
        LinkedList<Integer> visitQueue = new LinkedList<Integer>();

        for(TautomericGroup group : groups)
            visitQueue.addAll(group.atoms);

        while(!visitQueue.isEmpty())
        {
            Integer atom = visitQueue.remove();

            if(!selectedAtoms.contains(atom))
            {
                selectedAtoms.add(atom);

                for(Integer bond : atomAlternBonds.get(atom))
                {
                    if(!selectedBonds.contains(bond))
                    {
                        selectedBonds.add(bond);
                        visitQueue.add(getOtherBondAtom(bond, atom));
                    }
                }
            }
        }


        for(int atom : selectedAtoms)
        {
            opens[atom] = atomAlternBonds.get(atom).size();
            losts[atom] = atomAlternBonds.get(atom).size() - singleBondCounts[atom];
        }

        for(int group = 0; group < groups.size(); group++)
        {
            for(int atom : groups.get(group).atoms)
            {
                IAtom skeletonAtom = skeleton.getAtom(atom);
                int charge = skeletonAtom.getFormalCharge();
                int hydrogen = skeletonAtom.getImplicitHydrogenCount();

                if(charge < 0)
                {
                    if(skeletonAtom.getValency() != null)
                        skeletonAtom.setValency(skeletonAtom.getValency() - charge);

                    skeletonAtom.setFormalCharge(0);
                    hydrogen -= charge;
                    groupCharges[group] -= charge;
                }

                int mobile = Math.min(hydrogen, singleBondCounts[atom]);

                skeletonAtom.setImplicitHydrogenCount(hydrogen - mobile);
                groupMobiles[group] += mobile;
                losts[atom] += mobile;
            }
        }


        determined = new int[selectedBonds.size()];
        queue = new int[selectedAtoms.size()];
        queuing = new boolean[molecule.getAtomCount()];

        generateHydrogenCombinations(0, 0, 0);
    }


    private void generateHydrogenCombinations(int group, int groupAtom, int setted) throws CDKException
    {
        if(group < groups.size())
        {
            List<Integer> groupAtoms = groups.get(group).atoms;
            int atom = groupAtoms.get(groupAtom);
            int limit = Math.min(losts[atom], groupMobiles[group] - setted);

            for(int i = 0; i <= limit; i++)
            {
                mobiles[atom] += i > 0 ? 1 : 0;
                losts[atom] -= i > 0 ? 1 : 0;

                if(groupAtom < groupAtoms.size() - 1)
                    generateHydrogenCombinations(group, groupAtom + 1, setted + i);
                else if(setted + i == groupMobiles[group])
                    generateHydrogenCombinations(group + 1, 0, 0);
            }

            mobiles[atom] -= limit;
            losts[atom] += limit;
        }
        else
        {
            combination++;

            if(combination == 1024 * 1)
                throw new InChITautomerException("too many hydrogen combinations");

            for(int atom : selectedAtoms)
                if(losts[atom] == 0 && opens[atom] > 0 || losts[atom] == opens[atom])
                    addToQueue(atom);

            filters.clear();
            findBondOrders();
        }
    }


    private void findBondOrders() throws CDKException
    {
        int bondIdxOrig = determinedBondIdx;

        try
        {
            while(queueIdx > 0)
            {
                int atom = queue[--queueIdx];
                queuing[atom] = false;

                int order = losts[atom] == 0 ? 1 : 2;

                for(int bond : atomAlternBonds.get(atom))
                {
                    int other = getOtherBondAtom(bond, atom);

                    if(orders[bond] == 0)
                    {
                        orders[bond] = order;

                        opens[atom]--;
                        opens[other]--;

                        losts[atom] -= order - 1;
                        losts[other] -= order - 1;

                        determined[determinedBondIdx++] = bond;

                        if(losts[other] < 0 || losts[other] > opens[other])
                            return;

                        if(losts[other] == 0 && opens[other] > 0 || losts[other] == opens[other])
                            addToQueue(other);
                    }
                }
            }

            if(determinedBondIdx == selectedBonds.size())
            {
                loop:
                for(int[] other : filters)
                {
                    for(int bond : selectedBonds)
                        if(!aromatics[bond] && orders[bond] != other[bond])
                            continue loop;

                    return;
                }

                filters.add(orders.clone());
                generateChargeCombinations(0, 0, 0);
            }
            else
            {
                int atom = 0;

                while(opens[atom] == 0)
                    atom++;

                int bond = -1;

                for(int b : atomAlternBonds.get(atom))
                {
                    if(orders[b] == 0)
                    {
                        bond = b;
                        break;
                    }
                }

                int other = getOtherBondAtom(bond, atom);

                // set the bond as single
                determined[determinedBondIdx++] = bond;
                orders[bond] = 1;

                if(losts[atom] == --opens[atom])
                    addToQueue(atom);

                if(losts[other] == --opens[other])
                    addToQueue(other);

                findBondOrders();

                // set the bond as double
                orders[bond] = 2;

                if(--losts[atom] == 0 && opens[atom] > 0)
                    addToQueue(atom);

                if(--losts[other] == 0 && opens[other] > 0)
                    addToQueue(other);

                findBondOrders();
            }
        }
        finally
        {
            while(queueIdx > 0)
                queuing[queue[--queueIdx]] = false;

            while(determinedBondIdx != bondIdxOrig)
            {
                int bond = determined[--determinedBondIdx];
                int atom1 = bondAtoms[2 * bond + 0];
                int atom2 = bondAtoms[2 * bond + 1];

                losts[atom1] += orders[bond] - 1;
                losts[atom2] += orders[bond] - 1;

                opens[atom1]++;
                opens[atom2]++;

                orders[bond] = 0;
            }
        }
    }


    private void generateChargeCombinations(int group, int groupAtom, int setted) throws CDKException
    {
        if(group < groups.size())
        {
            List<Integer> groupAtoms = groups.get(group).atoms;
            int atom = groupAtoms.get(groupAtom);
            int limit = Math.min(mobiles[atom], groupCharges[group] - setted);

            for(int i = 0; i <= limit; i++)
            {
                charges[atom] += i > 0 ? 1 : 0;

                if(groupAtom < groupAtoms.size() - 1)
                    generateChargeCombinations(group, groupAtom + 1, setted + i);
                else if(setted + i == groupCharges[group])
                    generateChargeCombinations(group + 1, 0, 0);
            }

            charges[atom] -= limit;
        }
        else
        {
            createTautomerMolecule();
        }
    }


    int count = 0;

    private void createTautomerMolecule() throws CDKException
    {
        IAtomContainer tautomer = clone(skeleton);

        for(int bond : selectedBonds)
        {
            tautomer.getBond(bond).setOrder(orders[bond] == 1 ? Order.SINGLE : Order.DOUBLE);
        }

        for(int atom : selectedAtoms)
        {
            IAtom tautomerAtom = tautomer.getAtom(atom);

            int hydrogens = tautomerAtom.getImplicitHydrogenCount() + mobiles[atom] - charges[atom];
            int charge = tautomerAtom.getFormalCharge() - charges[atom];

            tautomerAtom.setImplicitHydrogenCount(hydrogens);
            tautomerAtom.setFormalCharge(charge);

            if(tautomerAtom.getValency() != null)
                tautomerAtom.setValency(tautomerAtom.getValency() - charges[atom]);
        }

        //if(tautomers.size() == 1024 * 1)
        //    throw new InChITautomerException("too many tautomers");

        tautomers.add(tautomer);
    }


    private IAtomContainer clone(IAtomContainer molecule) throws CDKException
    {
        try
        {
            return molecule.clone();
        }
        catch(CloneNotSupportedException e)
        {
            throw new CDKException("molecule cloning is not supported", e);
        }
    }


    private int getOtherBondAtom(int bond, int atom)
    {
        return bondAtoms[2 * bond] == atom ? bondAtoms[2 * bond + 1] : bondAtoms[2 * bond];
    }

    private void addToQueue(int atom)
    {
        if(!queuing[atom])
        {
            queue[queueIdx++] = atom;
            queuing[atom] = true;
        }
    }
}
