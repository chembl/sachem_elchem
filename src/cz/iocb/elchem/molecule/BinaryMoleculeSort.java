/*
 * Copyright (C) 2015-2019 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2009-2009 Mark Rijnbeek    markr@ebi.ac.uk
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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;



public class BinaryMoleculeSort
{
    public static IAtom[] atomsByFrequency(IAtomContainer iac)
    {
        // create a map with (key,value) being (atom symbol, overall count)
        Map<Integer, Integer> elementCounts = new TreeMap<Integer, Integer>();

        for(IAtom atom : iac.atoms())
        {
            Integer count = elementCounts.get(atom.getAtomicNumber());

            if(count == null)
                elementCounts.put(atom.getAtomicNumber(), new Integer(1));
            else
                elementCounts.put(atom.getAtomicNumber(), ++count);
        }

        // create a map with (key,value) being (IAtom, number of bond IAtom occurs in)
        Map<IAtom, Integer> bondParticipationCount = new HashMap<IAtom, Integer>();

        for(IBond bond : iac.bonds())
        {
            if(bond.getAtom(0).getAtomicNumber() == Molecule.AtomType.H
                    || bond.getAtom(1).getAtomicNumber() == Molecule.AtomType.H)
                continue;

            for(IAtom atomInBond : bond.atoms())
            {
                if(!bondParticipationCount.containsKey(atomInBond))
                    bondParticipationCount.put(atomInBond, 1);
                else
                    bondParticipationCount.put(atomInBond, bondParticipationCount.get(atomInBond) + 1);
            }
        }

        // mind the atoms not in any bond
        for(IAtom atom : iac.atoms())
        {
            if(!bondParticipationCount.containsKey(atom))
                bondParticipationCount.put(atom, 0);
        }

        // we now have to maps that will be used to sort the incoming atom container
        List<AtomForIsomorphismSort> atomList = new ArrayList<AtomForIsomorphismSort>();

        for(IAtom atom : iac.atoms())
        {
            AtomForIsomorphismSort afis = new AtomForIsomorphismSort(atom, elementCounts.get(atom.getAtomicNumber()),
                    bondParticipationCount.get(atom));
            atomList.add(afis);
        }

        Collections.sort(atomList, new AtomForIsomorphismSortComparator());

        // create an output atom array based on the sorted list
        IAtom[] iAtomArray = new IAtom[iac.getAtomCount()];
        int iacSortedIdx = 0;

        for(AtomForIsomorphismSort afis : atomList)
        {
            iAtomArray[iacSortedIdx] = afis.iatom;
            iacSortedIdx++;
        }

        return iAtomArray;
    }


    private static class AtomForIsomorphismSort
    {
        IAtom iatom;
        Integer overallElementCount;
        Integer atomBondParticipationCount;

        public AtomForIsomorphismSort(IAtom _iatom, int _overallElementCount, int _atomBondParticipationCount)
        {
            iatom = _iatom;
            overallElementCount = _overallElementCount;
            atomBondParticipationCount = _atomBondParticipationCount;
        }
    }


    private static class AtomForIsomorphismSortComparator implements Comparator<AtomForIsomorphismSort>
    {
        @Override
        public int compare(AtomForIsomorphismSort e1, AtomForIsomorphismSort e2)
        {
            if(e1.iatom.getAtomicNumber() == Molecule.AtomType.H && e2.iatom.getAtomicNumber() != Molecule.AtomType.H)
                return 1;

            if(e1.iatom.getAtomicNumber() != Molecule.AtomType.H && e2.iatom.getAtomicNumber() == Molecule.AtomType.H)
                return -1;


            if(e1.iatom.getAtomicNumber() == Molecule.AtomType.C && e2.iatom.getAtomicNumber() != Molecule.AtomType.C)
                return 1;

            if(e1.iatom.getAtomicNumber() != Molecule.AtomType.C && e2.iatom.getAtomicNumber() == Molecule.AtomType.C)
                return -1;


            if(e1.overallElementCount.compareTo(e2.overallElementCount) == 0)
                return e2.atomBondParticipationCount.compareTo(e1.atomBondParticipationCount);
            else
                return e1.overallElementCount.compareTo(e2.overallElementCount);
        }
    }
}
