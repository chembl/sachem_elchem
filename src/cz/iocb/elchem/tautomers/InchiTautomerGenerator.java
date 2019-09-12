/*
 * Copyright (C) 2015-2019 Jakub Galgonek   galgonek@uochb.cas.cz
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
package cz.iocb.elchem.tautomers;


import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.TimeoutException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.AtomContainer;
import cz.iocb.elchem.molecule.BinaryMoleculeBuilder;
import cz.iocb.elchem.molecule.Molecule;
import cz.iocb.elchem.tautomers.InChI.Fragment;



/**
 * Creates tautomers for a given input molecule, based on the mobile H atoms listed in the InChI. Algorithm described in
 * {@cdk.cite Thalheim2010}.
 */
public class InchiTautomerGenerator
{
    private static class MobileHydrogenConfig
    {
        int hydrogenCount;
        int negativeChargeCount;
        List<Integer> positions;
    }


    private final long time = System.currentTimeMillis();
    private final long tautomerCombinationLimit;
    private final BigInteger hydrogenCombinationLimit;
    private final long timeLimit;


    public InchiTautomerGenerator(long tautomerCombinationLimit, long hydrogenCombinationLimit, long timeLimit)
    {
        this.tautomerCombinationLimit = tautomerCombinationLimit;
        this.hydrogenCombinationLimit = BigInteger.valueOf(((Long) hydrogenCombinationLimit).longValue());
        this.timeLimit = timeLimit;
    }


    public InchiTautomerGenerator()
    {
        this(16 * 1024, 16 * 1024, 5 * 60000);
    }


    /**
     * Produces a list of tautomers generated based on the input molecule.
     *
     * @param kekuleMDL molecule in MDL format, use kekule notation for any sort of success
     * @return tautomers generated for the input molecule
     * @throws IOException
     * @throws CDKException
     * @throws CombinationCountException
     * @throws CloneNotSupportedException
     * @throws TimeoutException
     * @throws InChIException
     */
    public List<IAtomContainer> getTautomers(IAtomContainer molecule) throws CDKException, IOException,
            TimeoutException, CloneNotSupportedException, CombinationCountException, InChIException
    {
        List<IAtomContainer> tautomers = new ArrayList<IAtomContainer>();
        molecule = molecule.clone();

        InChI inchi = new InChI(molecule);
        assignAtomLabel(molecule, inchi.getAuxInfo());

        List<IAtom> remainingAtoms = new LinkedList<IAtom>();

        for(IAtom a : molecule.atoms())
            if(a.getID() == null)
                remainingAtoms.add(a);

        for(IBond b : molecule.bonds())
        {
            IAtom a0 = b.getAtom(0);
            IAtom a1 = b.getAtom(1);

            if(a0.getID() != null && a1.getAtomicNumber() == Molecule.AtomType.H)
                remainingAtoms.remove(a1);

            if(a1.getID() != null && a0.getAtomicNumber() == Molecule.AtomType.H)
                remainingAtoms.remove(a0);
        }

        /*
        for(IAtom atom : remainingAtoms)
        {
            for(IBond bond : molecule.bonds())
                for(IAtom a : bond.atoms())
                    if(a == atom)
                        throw new CDKException("remaining atom has a bond");
        }
        */

        List<IBond> crossBonds = getCrossBondLabels(molecule);

        tautomers.add(new AtomContainer());

        int idx = 0;

        for(Fragment localInchi : inchi.decompose())
        {
            List<IAtomContainer> localTautomers = getTautomers(molecule, localInchi, Integer.toString(++idx));

            List<IAtomContainer> tmps = tautomers;
            tautomers = new ArrayList<IAtomContainer>();

            if(tmps.size() * localTautomers.size() > tautomerCombinationLimit)
                throw new CombinationCountException("too many tautomer combinations");

            for(IAtomContainer tmp : tmps)
            {
                for(IAtomContainer t : localTautomers)
                {
                    IAtomContainer container = tmp.clone();
                    container.add(t.clone());
                    tautomers.add(container);
                }
            }
        }

        for(IBond bond : crossBonds)
        {
            String label0 = bond.getAtom(0).getID();
            String label1 = bond.getAtom(1).getID();

            for(IAtomContainer tautomer : tautomers)
            {
                IBond newBond = tautomer.getBuilder().newInstance(IBond.class, getAtomById(tautomer, label0),
                        getAtomById(tautomer, label1), bond.getOrder());
                tautomer.addBond(newBond);
            }
        }

        for(IAtom a : remainingAtoms)
            for(IAtomContainer tautomer : tautomers)
                tautomer.addAtom(a.clone());


        return tautomers;
    }


    public static IAtom getAtomById(IAtomContainer molecule, String id)
    {
        for(IAtom atom : molecule.atoms())
        {
            if(id.equals(atom.getID()))
                return atom;
        }

        return null;
    }


    private void assignAtomLabel(IAtomContainer molecule, String auxInfo)
    {
        if(auxInfo.matches("AuxInfo=1/[0-9]*/N:.*"))
        {
            String positions = auxInfo.substring(14).replaceFirst("/.*", "");

            int major = 0;

            for(String majorPositions : positions.split(";"))
            {
                major++;

                int minor = 0;

                for(String minorPositions : majorPositions.split(","))
                {
                    minor++;

                    IAtom atom = molecule.getAtom(Integer.valueOf(minorPositions) - 1);
                    atom.setID(major + ":" + minor);
                }
            }
        }
    }


    private List<IBond> getCrossBondLabels(IAtomContainer molecule)
    {
        ArrayList<IBond> selectedBonds = new ArrayList<IBond>();

        for(IBond bond : molecule.bonds())
        {
            String id0 = bond.getAtom(0).getID();
            String id1 = bond.getAtom(1).getID();

            if(id0 == null || id1 == null)
                continue;

            if(!id0.split(":")[0].equals(id1.split(":")[0]))
                selectedBonds.add(bond);
        }

        return selectedBonds;
    }


    /**
     * Get tautomers for an input molecule with the InChI already provided as input argument.
     *
     * @param inputMolecule and input molecule for which to generate tautomers
     * @param fragment InChI for the inpute molecule
     * @param aux
     * @return a list of tautomers
     * @throws TimeoutException
     * @throws CDKException
     * @throws CloneNotSupportedException
     * @throws CombinationCountException
     */
    private List<IAtomContainer> getTautomers(IAtomContainer inputMolecule, Fragment fragment, String prefix)
            throws TimeoutException, CDKException, CloneNotSupportedException, CombinationCountException
    {
        List<IAtomContainer> tautomers = new ArrayList<IAtomContainer>();

        inputMolecule = getMoleculeFragment(inputMolecule, fragment, prefix);
        setUnsetBondOrders(inputMolecule);

        List<MobileHydrogenConfig> mobileHydrogens = parseMobileHydrogens(fragment.hydrogens);
        int[] valency = getValency(inputMolecule);
        resetMobileHydrogens(inputMolecule, mobileHydrogens, valency);

        tautomers = constructTautomers(inputMolecule, mobileHydrogens, valency);

        // remove duplicates
        return removeDuplicatesSimple(tautomers);
    }


    /**
     * Gets the part of the input molecule that corresponds to the given InChI fragment.
     *
     * @param inputMolecule user input molecule
     * @throws CDKException
     */
    private static IAtomContainer getMoleculeFragment(IAtomContainer inputMolecule, Fragment fragment, String prefix)
            throws CDKException, CloneNotSupportedException
    {
        inputMolecule = inputMolecule.clone();

        for(IAtom atom : inputMolecule.atoms())
            atom.setID(null);

        ArrayList<IAtom> atoms = new ArrayList<IAtom>();

        StringTokenizer tokenizer = new StringTokenizer(fragment.aux, ",");


        int id = 0;

        while(tokenizer.hasMoreTokens())
        {
            id++;
            Integer position = Integer.valueOf(tokenizer.nextToken());

            IAtom atom = inputMolecule.getAtom(position - 1);
            atom.setID(prefix + ":" + id);
            atoms.add(atom);
        }

        for(IBond bond : inputMolecule.bonds())
        {
            IAtom atom0 = bond.getAtom(0);
            IAtom atom1 = bond.getAtom(1);

            if(atom0.getAtomicNumber() == Molecule.AtomType.H && atom1.getID() != null && atom0.getID() == null)
            {
                atoms.add(atom0);
                atom0.setID(prefix + ":" + ++id);
            }

            if(atom1.getAtomicNumber() == Molecule.AtomType.H && atom0.getID() != null && atom1.getID() == null)
            {
                atoms.add(atom1);
                atom1.setID(prefix + ":" + ++id);
            }
        }


        List<IAtom> remove = new ArrayList<IAtom>();

        for(IAtom a : inputMolecule.atoms())
            if(a.getID() == null)
                remove.add(a);


        for(IAtom a : remove)
            inputMolecule.removeAtom(a);


        List<IBond> removeBonds = new ArrayList<IBond>();

        for(IBond bond : inputMolecule.bonds())
        {
            boolean rem = false;

            for(IAtom a : bond.atoms())
                if(remove.contains(a))
                    rem = true;

            if(rem)
                removeBonds.add(bond);
        }

        for(IBond b : removeBonds)
        {
            IAtom atom0 = b.getAtom(0);
            IAtom atom1 = b.getAtom(1);

            atom0.setValency(atom0.getValency() - 1);
            atom1.setValency(atom1.getValency() - 1);

            inputMolecule.removeBond(b);
        }

        inputMolecule.setAtoms(atoms.toArray(new IAtom[0]));


        ArrayList<IAtom> chiralAtoms = new ArrayList<IAtom>();

        if(fragment.tetrahedralStereo != null)
        {
            for(String definition : fragment.tetrahedralStereo.split(","))
            {
                String index = definition.substring(0, definition.length() - 1);
                IAtom atom = inputMolecule.getAtom(Integer.parseInt(index) - 1);
                chiralAtoms.add(atom);
            }
        }

        for(IAtom atom : inputMolecule.atoms())
            if(!chiralAtoms.contains(atom))
                atom.setProperty(BinaryMoleculeBuilder.STEREO_PROPERTY, BinaryMoleculeBuilder.IGNORE_STEREO);


        ArrayList<IBond> stereoBonds = new ArrayList<IBond>();

        if(fragment.doubleBondStereo != null)
        {
            for(String definition : fragment.doubleBondStereo.split(","))
            {
                String[] indexes = definition.substring(0, definition.length() - 1).split("-");

                IAtom a0 = inputMolecule.getAtom(Integer.parseInt(indexes[0]) - 1);
                IAtom a1 = inputMolecule.getAtom(Integer.parseInt(indexes[1]) - 1);

                IBond bond = inputMolecule.getBond(a0, a1);
                stereoBonds.add(bond);
            }
        }

        for(IBond bond : inputMolecule.bonds())
            if(!stereoBonds.contains(bond))
                bond.setProperty(BinaryMoleculeBuilder.STEREO_PROPERTY, BinaryMoleculeBuilder.IGNORE_STEREO);


        return inputMolecule;
    }


    private static void setUnsetBondOrders(IAtomContainer inputMolecule) throws CDKException
    {
        for(IBond b : inputMolecule.bonds())
        {
            if(b.getOrder() == IBond.Order.UNSET && b.isAromatic())
            {
                for(IBond bond : inputMolecule.bonds())
                    if(bond.isAromatic())
                        for(IAtom atom : bond.atoms())
                            atom.setFlag(CDKConstants.ISAROMATIC, true);

                Kekulization.kekulize(inputMolecule);
                break;
            }
        }
    }


    /**
     * Parses mobile H group(s) in an InChI String.
     *
     * Multiple InChI sequences of mobile hydrogens are joined into a single sequence (list), see step 1 of algorithm in
     * paper. Mobile H group has syntax (H[n][-[m]],a1,a2[,a3[,a4...]]) Brackets [ ] surround optional terms.
     * <ul>
     * <li>Term H[n] stands for 1 or, if the number n (n>1) is present, n mobile hydrogen atoms.</li>
     * <li>Term [-[m]], if present, stands for 1 or, if the number m (m>1) is present, m mobile negative charges.</li>
     * <li>a1,a2[,a3[,a4...]] are canonical numbers of atoms in the mobile H group.</li>
     * <li>no two mobile H groups may have an atom (a canonical number) in common.</li>
     * </ul>
     *
     * @param mobHydrAttachPositions list of positions where mobile H can attach
     * @param inputInchi InChI input
     * @return overall count of hydrogens to be dispersed over the positions
     */
    private static List<MobileHydrogenConfig> parseMobileHydrogens(String hydrogens)
    {
        List<MobileHydrogenConfig> mobileHydrogenList = new ArrayList<MobileHydrogenConfig>();

        if(hydrogens == null || hydrogens.indexOf("(H") == -1)
            return mobileHydrogenList;


        String mobileHydrogens = hydrogens.substring(hydrogens.indexOf('('));
        Pattern mobileHydrPattern = Pattern.compile("\\((.)*?\\)");
        Matcher match = mobileHydrPattern.matcher(mobileHydrogens);

        while(match.find())
        {
            String mobileHGroup = match.group();
            String head = mobileHGroup.substring(0, mobileHGroup.indexOf(',') + 1);

            if(head.contains("H,"))
                head = head.replace("H,", "H1,");
            if(head.contains("-,"))
                head = head.replace("-,", "-1,");

            head = head.substring(2);

            /*
             Pragmatically, also add any delocalised neg charge to the mobile H count.
             Based on examples like:
                C[N+](C)(C)CCCCC\C=C(/NC(=O)C1CC1(Cl)Cl)\C(=O)O  -> ...(H-,18,20,21,22)
                COc1cc(N)c(Cl)cc1C(=O)NC2C[N+]3(CCl)CCC2CC3      -> ...(H2-,19,20,22)
            */

            String[] counts = head.substring(0, head.length() - 1).split("-");

            MobileHydrogenConfig conf = new MobileHydrogenConfig();

            if(counts.length == 2)
                conf.negativeChargeCount = Integer.valueOf(counts[1]);

            if(!counts[0].isEmpty())
                conf.hydrogenCount = Integer.valueOf(counts[0]);


            mobileHGroup = mobileHGroup.substring(mobileHGroup.indexOf(',') + 1).replace(")", "");
            StringTokenizer tokenizer = new StringTokenizer(mobileHGroup, ",");

            List<Integer> positions = new ArrayList<Integer>();

            while(tokenizer.hasMoreTokens())
            {
                Integer position = Integer.valueOf(tokenizer.nextToken());
                positions.add(position);
            }

            conf.positions = positions;
            mobileHydrogenList.add(conf);
        }

        return mobileHydrogenList;
    }


    private static int[] getValency(IAtomContainer inputMolecule)
    {
        int[] valency = new int[inputMolecule.getAtomCount()];
        int idx = 0;

        for(IAtom atom : inputMolecule.atoms())
        {
            int val = atom.getImplicitHydrogenCount() + getConnectivity(atom, inputMolecule);
            valency[idx++] = val;
        }

        return valency;
    }


    private static void resetMobileHydrogens(IAtomContainer inputMolecule, List<MobileHydrogenConfig> mobileHydrogens,
            int[] valency)
    {
        for(MobileHydrogenConfig mobile : mobileHydrogens)
        {
            int mobileHydrogenCount = 0;
            int mobileChargeCount = 0;

            for(int l : mobile.positions)
            {
                IAtom atom = inputMolecule.getAtom(l - 1);

                mobileHydrogenCount += atom.getImplicitHydrogenCount();
                atom.setImplicitHydrogenCount(0);
            }

            for(int l : mobile.positions)
            {
                IAtom atom = inputMolecule.getAtom(l - 1);

                if(atom.getFormalCharge() < 0)
                {
                    mobileChargeCount -= atom.getFormalCharge();
                    atom.setValency(atom.getValency() - atom.getFormalCharge());
                    valency[l - 1] -= atom.getFormalCharge();
                    atom.setFormalCharge(0);
                }
            }

            mobile.hydrogenCount = mobileHydrogenCount;
            mobile.negativeChargeCount = mobileChargeCount;
        }
    }


    /**
     * Constructs tautomers following (most) steps of the algorithm in {@cdk.cite Thalheim2010}.
     *
     * @param inputMolecule input molecule
     * @param valency
     * @param mobHydrAttachPositions mobile H positions
     * @param totalMobHydrCount count of mobile hydrogens in molecule
     * @return tautomers
     * @throws TimeoutException
     * @throws CloneNotSupportedException
     * @throws CombinationCountException
     */
    private List<IAtomContainer> constructTautomers(IAtomContainer inputMolecule,
            List<MobileHydrogenConfig> mobileHydrogens, int[] valency)
            throws CloneNotSupportedException, TimeoutException, CombinationCountException
    {
        List<IAtomContainer> tautomers = new ArrayList<IAtomContainer>();

        //Tautomeric skeleton generation
        IAtomContainer skeleton = inputMolecule.clone();


        int doubleBondCount = 0;
        for(IBond bond : skeleton.bonds())
        {
            if(bond.getOrder().equals(IBond.Order.DOUBLE))
            {
                doubleBondCount++;
                bond.setOrder(IBond.Order.SINGLE);
            }
        }


        int[] freeConnectivityBase = new int[skeleton.getAtomCount()];

        for(int i = 0; i < skeleton.getAtomCount(); i++)
            freeConnectivityBase[i] = valency[i] - getConnectivity(skeleton.getAtom(i), skeleton);


        // Make combinations for mobile Hydrogen attachments
        BigInteger combinationCount = countHydrogenPositions(skeleton, mobileHydrogens, freeConnectivityBase);

        if(combinationCount.compareTo(hydrogenCombinationLimit) == 1)
            throw new CombinationCountException("too many mobile hydrogen combinations: " + combinationCount);

        List<List<Integer>> combinations = combineHydrogenPositions(skeleton, mobileHydrogens, freeConnectivityBase);

        for(List<Integer> hPositions : combinations)
        {
            int[] freeConnectivity = freeConnectivityBase.clone();


            IAtomContainer tautomerSkeleton = skeleton.clone();
            for(Integer hPos : hPositions)
            {
                if(hPos > 0)
                {
                    IAtom atom = tautomerSkeleton.getAtom(hPos - 1);
                    atom.setImplicitHydrogenCount(atom.getImplicitHydrogenCount() + 1);
                }
                else
                {
                    IAtom atom = tautomerSkeleton.getAtom(-hPos - 1);
                    atom.setFormalCharge(atom.getFormalCharge() - 1);
                    atom.setValency(atom.getValency() - 1);
                    freeConnectivity[-hPos - 1]--;
                }
            }


            generateDoubleBondCombinations(tautomers, 0, tautomerSkeleton, 0, 0, doubleBondCount, freeConnectivity);
        }


        return tautomers;
    }


    private static BigInteger countHydrogenPositions(IAtomContainer skeleton,
            List<MobileHydrogenConfig> mobileHydrogens, int[] freeConnectivity)
    {
        BigInteger count = new BigInteger("1");

        for(MobileHydrogenConfig mobile : mobileHydrogens)
        {
            if(mobile.hydrogenCount > 0)
                count = count.multiply(countHydrogenPositions(0, 0, skeleton, mobile.hydrogenCount, mobile.positions,
                        freeConnectivity));

            if(mobile.negativeChargeCount > 0)
                count = count.multiply(countHydrogenPositions(0, 0, skeleton, -mobile.negativeChargeCount,
                        mobile.positions, freeConnectivity));
        }

        return count;
    }


    private static BigInteger countHydrogenPositions(int offset, int taken, IAtomContainer skeleton,
            int totalMobHydrCount, List<Integer> mobHydrAttachPositions, int[] freeConnectivity)
    {
        BigInteger count = new BigInteger("0");
        BigInteger one = new BigInteger("1");

        if(offset >= mobHydrAttachPositions.size())
            return count;

        int pos = mobHydrAttachPositions.get(offset);
        int maxCnt = freeConnectivity[pos - 1];


        count = count.add(countHydrogenPositions(offset + 1, taken, skeleton, totalMobHydrCount, mobHydrAttachPositions,
                freeConnectivity));

        for(int i = 0; i < maxCnt; i++)
        {
            taken++;

            if(taken == Math.abs(totalMobHydrCount))
            {
                count = count.add(one);
                break;
            }

            count = count.add(countHydrogenPositions(offset + 1, taken, skeleton, totalMobHydrCount,
                    mobHydrAttachPositions, freeConnectivity));
        }

        return count;
    }


    private List<List<Integer>> combineHydrogenPositions(IAtomContainer skeleton,
            List<MobileHydrogenConfig> mobileHydrogens, int[] freeConnectivity) throws TimeoutException
    {
        List<List<Integer>> combinations = new ArrayList<List<Integer>>();
        combinations.add(new ArrayList<Integer>());

        for(MobileHydrogenConfig mobile : mobileHydrogens)
        {
            if(mobile.hydrogenCount > 0)
            {
                List<List<Integer>> localCombinations = combineHydrogenPositions(skeleton, mobile.hydrogenCount,
                        mobile.positions, freeConnectivity);

                combinations = combineHydrogenPositions(combinations, localCombinations);
            }

            if(mobile.negativeChargeCount > 0)
            {
                List<List<Integer>> localCombinations = combineHydrogenPositions(skeleton, -mobile.negativeChargeCount,
                        mobile.positions, freeConnectivity);

                combinations = combineHydrogenPositions(combinations, localCombinations);
            }
        }

        return combinations;
    }


    private List<List<Integer>> combineHydrogenPositions(IAtomContainer skeleton, int totalMobHydrCount,
            List<Integer> mobHydrAttachPositions, int[] freeConnectivity) throws TimeoutException
    {
        List<List<Integer>> combinations = new ArrayList<List<Integer>>();

        combineHydrogenPositions(0, new ArrayList<Integer>(), combinations, skeleton, totalMobHydrCount,
                mobHydrAttachPositions, freeConnectivity);

        return combinations;
    }


    /**
     * Makes combinations recursively of all possible mobile Hydrogen positions.
     *
     * @param taken positions taken by hydrogen
     * @param combinations combinations made so far
     * @param skeleton container to work on
     * @param totalMobHydrCount
     * @param mobHydrAttachPositions
     * @throws TimeoutException
     */
    private void combineHydrogenPositions(int offset, List<Integer> taken, List<List<Integer>> combinations,
            IAtomContainer skeleton, int totalMobHydrCount, List<Integer> mobHydrAttachPositions,
            int[] freeConnectivity) throws TimeoutException
    {
        checkTime();


        if(offset >= mobHydrAttachPositions.size())
            return;

        int pos = mobHydrAttachPositions.get(offset);
        int maxCnt = freeConnectivity[pos - 1];


        combineHydrogenPositions(offset + 1, taken, combinations, skeleton, totalMobHydrCount, mobHydrAttachPositions,
                freeConnectivity);

        int added = 0;

        for(int i = 0; i < maxCnt; i++)
        {
            taken.add(Integer.signum(totalMobHydrCount) * pos);
            added++;

            if(taken.size() == Math.abs(totalMobHydrCount))
            {
                List<Integer> addList = new ArrayList<Integer>(taken.size());
                addList.addAll(taken);
                combinations.add(addList);
                break;
            }

            combineHydrogenPositions(offset + 1, taken, combinations, skeleton, totalMobHydrCount,
                    mobHydrAttachPositions, freeConnectivity);
        }

        for(int i = 0; i < added; i++)
            taken.remove(taken.size() - 1);
    }


    private static List<List<Integer>> combineHydrogenPositions(List<List<Integer>> combinations,
            List<List<Integer>> localCombinations)
    {
        List<List<Integer>> oldCombinations = combinations;
        combinations = new ArrayList<List<Integer>>();

        for(List<Integer> newCombination : localCombinations)
        {
            for(List<Integer> oldCombination : oldCombinations)
            {
                ArrayList<Integer> l = new ArrayList<Integer>();
                l.addAll(oldCombination);
                l.addAll(newCombination);

                combinations.add(l);
            }
        }

        return combinations;
    }


    /**
     * Tries double bond combinations for a certain input container of which the double bonds have been stripped around
     * the mobile hydrogen positions. Recursive.
     *
     * @param container
     * @param dblBondsAdded counts double bonds added so far
     * @param bondOffSet offset for next double bond position to consider
     * @param doubleBondMax maximum number of double bonds to add
     * @param atomsInNeedOfFix atoms that require more bonds
     * @throws CloneNotSupportedException
     * @throws TimeoutException
     * @throws CombinationCountException
     */
    private void generateDoubleBondCombinations(List<IAtomContainer> result, int iter, IAtomContainer container,
            int dblBondsAdded, int bondOffSet, int doubleBondMax, int[] freeConnectivity)
            throws CloneNotSupportedException, TimeoutException, CombinationCountException
    {
        checkTime();

        HashSet<IAtom> atomsInNeedOfFix = new HashSet<IAtom>();
        List<IBond> clears = new ArrayList<IBond>();

        if(!clearDoubleBonds(bondOffSet - 1, container, atomsInNeedOfFix, clears, freeConnectivity))
        {
            // this configuration cannot be finished

            // revert container ...
            for(IBond b : clears)
            {
                b.setOrder(IBond.Order.SINGLE);
                freeConnectivity[container.indexOf(b.getAtom(0))]++;
                freeConnectivity[container.indexOf(b.getAtom(1))]++;
            }

            return;
        }
        else
        {
            dblBondsAdded += clears.size();
        }


        if(dblBondsAdded == doubleBondMax)
        {
            boolean validDoubleBondConfig = true;
            int idx = 0;

            for(IAtom atom : container.atoms())
            {
                if(freeConnectivity[idx] != atom.getImplicitHydrogenCount())
                {
                    validDoubleBondConfig = false;
                    break;
                }

                idx++;
            }

            if(!validDoubleBondConfig)
                return;

            IAtomContainer clone = container.clone();
            result.add(clone);

            if(result.size() > tautomerCombinationLimit)
                throw new CombinationCountException("too many tautomers");
        }


        if(doubleBondMax - dblBondsAdded > container.getBondCount() - bondOffSet)
        {
            // this configuration cannot be finished

            // revert container ...
            for(IBond b : clears)
            {
                b.setOrder(IBond.Order.SINGLE);
                freeConnectivity[container.indexOf(b.getAtom(0))]++;
                freeConnectivity[container.indexOf(b.getAtom(1))]++;
            }

            return;
        }


        for(int offSet = bondOffSet; offSet < container.getBondCount(); offSet++)
        {
            IBond bond = container.getBond(offSet);

            if(bond.getOrder() == IBond.Order.SINGLE && atomsInNeedOfFix.contains(bond.getAtom(0))
                    && atomsInNeedOfFix.contains(bond.getAtom(1)))
            {
                IAtom a0 = bond.getAtom(0);
                IAtom a1 = bond.getAtom(1);

                int a0idx = container.indexOf(a0);
                int a1idx = container.indexOf(a1);


                if(freeConnectivity[a0idx] > a0.getImplicitHydrogenCount()
                        && freeConnectivity[a1idx] > a1.getImplicitHydrogenCount())
                {
                    bond.setOrder(IBond.Order.DOUBLE);
                    dblBondsAdded = dblBondsAdded + 1;

                    freeConnectivity[a0idx]--;
                    freeConnectivity[a1idx]--;

                    generateDoubleBondCombinations(result, iter + 1, container, dblBondsAdded, offSet + 1,
                            doubleBondMax, freeConnectivity);

                    bond.setOrder(IBond.Order.SINGLE);
                    dblBondsAdded = dblBondsAdded - 1;

                    freeConnectivity[a0idx]++;
                    freeConnectivity[a1idx]++;
                }
            }
        }


        for(IBond b : clears)
        {
            b.setOrder(IBond.Order.SINGLE);
            freeConnectivity[container.indexOf(b.getAtom(0))]++;
            freeConnectivity[container.indexOf(b.getAtom(1))]++;
        }
    }


    private static boolean clearDoubleBonds(int offset, IAtomContainer tautomerSkeleton,
            HashSet<IAtom> atomsInNeedOfFix, List<IBond> bonds, int[] freeConnectivity)
    {
        atomsInNeedOfFix.clear();

        for(int idx = 0; idx < tautomerSkeleton.getAtomCount(); idx++)
        {
            IAtom atom = tautomerSkeleton.getAtom(idx);

            if(freeConnectivity[idx] != atom.getImplicitHydrogenCount())
                atomsInNeedOfFix.add(atom);
        }


        while(true)
        {
            boolean changed = false;

            for(int idx = 0; idx < tautomerSkeleton.getAtomCount(); idx++)
            {
                IAtom atom = tautomerSkeleton.getAtom(idx);
                int diff = freeConnectivity[idx] - atom.getImplicitHydrogenCount();


                if(diff < 0)
                    return false;

                if(diff > 0)
                {
                    List<IBond> possibilities = new LinkedList<IBond>();

                    for(int i = offset + 1; i < tautomerSkeleton.getBondCount(); i++)
                    {
                        IBond bond = tautomerSkeleton.getBond(i);

                        IAtom a0 = bond.getAtom(0);
                        IAtom a1 = bond.getAtom(1);

                        if(bond.getOrder() == IBond.Order.SINGLE && (a0 == atom && atomsInNeedOfFix.contains(a1)
                                || a1 == atom && atomsInNeedOfFix.contains(a0)))
                        {
                            possibilities.add(bond);
                        }
                    }


                    if(diff > possibilities.size())
                        return false;

                    if(diff == possibilities.size())
                    {
                        for(IBond bond : possibilities)
                        {
                            bond.setOrder(IBond.Order.DOUBLE);
                            bonds.add(bond);
                            changed = true;

                            freeConnectivity[tautomerSkeleton.indexOf(bond.getAtom(0))]--;
                            freeConnectivity[tautomerSkeleton.indexOf(bond.getAtom(1))]--;
                        }
                    }
                }
            }


            if(!changed)
                break;


            atomsInNeedOfFix.clear();

            for(int idx = 0; idx < tautomerSkeleton.getAtomCount(); idx++)
            {
                IAtom atom = tautomerSkeleton.getAtom(idx);

                if(freeConnectivity[idx] != atom.getImplicitHydrogenCount())
                    atomsInNeedOfFix.add(atom);
            }
        }

        return true;
    }


    private List<IAtomContainer> removeDuplicatesSimple(List<IAtomContainer> tautomers) throws TimeoutException
    {
        List<IAtomContainer> unique = new ArrayList<IAtomContainer>();

        if(tautomers.size() == 0)
            return unique;


        BitSet removed = new BitSet(tautomers.size());

        int atomCount = tautomers.get(0).getAtomCount();
        int bondCount = tautomers.get(0).getBondCount();

        byte[][] hydrogens = new byte[tautomers.size()][atomCount];
        byte[][] charges = new byte[tautomers.size()][atomCount];
        byte[][] bonds = new byte[tautomers.size()][bondCount];

        for(int idx = 0; idx < tautomers.size(); idx++)
        {
            IAtomContainer tautomer = tautomers.get(idx);

            for(int i = 0; i < atomCount; i++)
            {
                checkTime();

                IAtom atom = tautomer.getAtom(i);

                hydrogens[idx][i] = (byte) (int) atom.getImplicitHydrogenCount();
                charges[idx][i] = (byte) (int) atom.getFormalCharge();
            }

            for(int i = 0; i < bondCount; i++)
            {
                IBond bond = tautomer.getBond(i);

                if(bond.isAromatic())
                    bonds[idx][i] = -1;
                else
                    bonds[idx][i] = (byte) bond.getOrder().ordinal();
            }
        }

        for(int idx = 0; idx < tautomers.size(); idx++)
        {
            if(removed.get(idx))
                continue;

            byte[] hydrogens1 = hydrogens[idx];
            byte[] charges1 = charges[idx];
            byte[] bonds1 = bonds[idx];

            for(int idx2 = idx + 1; idx2 < tautomers.size(); idx2++)
            {
                if(removed.get(idx2))
                    continue;

                byte[] hydrogens2 = hydrogens[idx2];
                byte[] charges2 = charges[idx2];
                byte[] bonds2 = bonds[idx2];

                if(Arrays.equals(bonds1, bonds2) && Arrays.equals(hydrogens1, hydrogens2)
                        && Arrays.equals(charges1, charges2))
                {
                    removed.set(idx2);
                }
            }

            unique.add(tautomers.get(idx));
        }

        return unique;
    }


    /**
     * Sums the number of bonds (counting order) an atom is hooked up with.
     *
     * @param atom
     * @param container
     * @return
     */
    private static int getConnectivity(IAtom atom, IAtomContainer container)
    {
        int connectivity = 0;

        for(IBond bond : container.bonds())
        {
            if(bond.contains(atom))
            {
                switch(bond.getOrder())
                {
                    case SINGLE:
                        connectivity++;
                        break;
                    case DOUBLE:
                        connectivity += 2;
                        break;
                    case TRIPLE:
                        connectivity += 3;
                        break;
                    case QUADRUPLE:
                        connectivity += 4;
                        break;
                    case QUINTUPLE:
                        connectivity += 5;
                        break;
                    case SEXTUPLE:
                        connectivity += 6;
                        break;
                    case UNSET:
                        connectivity += 10;
                }
            }
        }

        return connectivity;
    }


    private void checkTime() throws TimeoutException
    {
        long duration = System.currentTimeMillis() - time;

        if(duration > timeLimit)
            throw new TimeoutException("exceeded time limit: " + duration);
    }
}
