/*
 * Copyright (C) 2015-2020 Jakub Galgonek   galgonek@uochb.cas.cz
 * Copyright (C) 2008-2009 Mark Rijnbeek    markr@ebi.ac.uk
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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeoutException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.RGroupQueryReader;
import org.openscience.cdk.io.formats.RGroupQueryFormat;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.isomorphism.matchers.RGroupQuery;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import cz.iocb.elchem.molecule.InChITautomerGenerator.InChITautomerException;
import cz.iocb.elchem.molecule.InChITools.InChIException;



public class MoleculeCreator
{
    public static class QueryMolecule
    {
        final public String name;
        final public List<IAtomContainer> tautomers;

        QueryMolecule(String name, List<IAtomContainer> tautomers)
        {
            this.name = name;
            this.tautomers = tautomers;
        }
    }


    public static final String STEREO_FLAG = "STEREO_FLAG";
    private static final int REQ_MODE_MIN_SB_RING_SHFT = 16;
    private static final int INCHI_QUERY_MODE = 8 << REQ_MODE_MIN_SB_RING_SHFT;
    private static final int INCHI_DB_MODE = 0;


    private static final ThreadLocal<Aromaticity> aromaticity = new ThreadLocal<Aromaticity>()
    {
        @Override
        protected Aromaticity initialValue()
        {
            return new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
        }
    };


    public static QueryMolecule translateQuery(String query, QueryFormat format, ChargeMode chargeMode,
            IsotopeMode isotopeMode, RadicalMode radicalMode, StereoMode stereoMode, AromaticityMode aromaticityMode,
            TautomerMode tautomerMode) throws CDKException, IOException, TimeoutException
    {
        List<IAtomContainer> molecules = null;

        switch(format)
        {
            case MOLFILE:
                molecules = Arrays.asList(MoleculeCreator.getMoleculeFromMolfile(query));
                break;

            case SMILES:
                molecules = Arrays.asList(MoleculeCreator.getMoleculeFromSmiles(query));
                break;

            case RGROUP:
                molecules = getMoleculesFromRGroupQuery(query);
                break;

            case UNSPECIFIED:
                List<String> lines = Arrays.asList(query.split("\\n"));

                if(((RGroupQueryFormat) RGroupQueryFormat.getInstance()).matches(lines).matched())
                    molecules = getMoleculesFromRGroupQuery(query);
                else if(lines.size() > 1)
                    molecules = Arrays.asList(MoleculeCreator.getMoleculeFromMolfile(query));
                else
                    molecules = Arrays.asList(MoleculeCreator.getMoleculeFromSmiles(query));

                break;
        }


        for(IAtomContainer molecule : molecules)
            configureAromaticity(molecule, aromaticityMode);

        String name = molecules.size() > 0 ? molecules.get(0).getTitle() : null;

        if(name == null)
            name = "";


        List<IAtomContainer> queries = new ArrayList<IAtomContainer>();

        if(tautomerMode == TautomerMode.IGNORE)
        {
            for(IAtomContainer molecule : molecules)
            {
                if(chargeMode == ChargeMode.IGNORE)
                    for(IAtom atom : molecule.atoms())
                        atom.setFormalCharge(0);

                if(radicalMode == RadicalMode.IGNORE)
                    for(IAtom atom : molecule.atoms())
                        atom.removeProperty(CDKConstants.SPIN_MULTIPLICITY);

                if(isotopeMode == IsotopeMode.IGNORE)
                    for(IAtom atom : molecule.atoms())
                        atom.setMassNumber(null);

                if(stereoMode == StereoMode.STRICT)
                    setStereo(molecule, new InChITools(molecule, INCHI_QUERY_MODE, false));

                queries.add(molecule);
            }
        }
        else
        {
            for(IAtomContainer molecule : molecules)
            {
                boolean resetStereo = false;

                if(isotopeMode == IsotopeMode.IGNORE)
                    for(IAtom atom : molecule.atoms())
                        atom.setMassNumber(null);

                if(stereoMode == StereoMode.STRICT && chargeMode == ChargeMode.IGNORE)
                    for(IAtom atom : molecule.atoms())
                        if(atom.getFormalCharge() != 0)
                            resetStereo = true;

                if(stereoMode == StereoMode.STRICT && radicalMode == RadicalMode.IGNORE)
                    for(IAtom atom : molecule.atoms())
                        if(atom.getProperty(CDKConstants.SPIN_MULTIPLICITY) != null)
                            resetStereo = true;


                try
                {
                    InChITools inchi = new InChITools(molecule, INCHI_QUERY_MODE, true);
                    resetStereo &= !inchi.getStereoAtoms().isEmpty() || !inchi.getStereoBonds().isEmpty();

                    if(stereoMode == StereoMode.STRICT)
                        setStereo(molecule, inchi);

                    for(IAtomContainer tautomer : InChITautomerGenerator.generate(molecule, inchi.getTautomericBonds(),
                            inchi.getTautomericGroups()))
                    {
                        if(chargeMode == ChargeMode.IGNORE)
                            for(IAtom atom : tautomer.atoms())
                                atom.setFormalCharge(0);

                        if(radicalMode == RadicalMode.IGNORE)
                            for(IAtom atom : tautomer.atoms())
                                atom.removeProperty(CDKConstants.SPIN_MULTIPLICITY);

                        if(resetStereo)
                        {
                            /* fix stereo */
                            InChITools tautomerInchi = new InChITools(tautomer, INCHI_QUERY_MODE, false);

                            HashSet<IChemObject> stereo = new HashSet<IChemObject>();
                            stereo.addAll(tautomerInchi.getStereoAtoms());
                            stereo.addAll(tautomerInchi.getStereoBonds());

                            for(IAtom atom : inchi.getStereoAtoms())
                                if(!stereo.contains(atom))
                                    atom.removeProperty(STEREO_FLAG);

                            for(IBond bond : inchi.getStereoBonds())
                                if(!stereo.contains(bond))
                                    bond.removeProperty(STEREO_FLAG);
                        }

                        queries.add(tautomer);
                    }
                }
                catch(InChIException e)
                {
                    throw new InChITautomerException(e.getMessage());
                }
            }
        }

        return new QueryMolecule(name, queries);
    }


    public static IAtomContainer translateMolecule(String mol, boolean inchiStereo) throws CDKException, IOException
    {
        IAtomContainer molecule = mol.contains("\n") ? getMoleculeFromMolfile(mol) : getMoleculeFromSmiles(mol);

        List<IBond> otherBonds = new LinkedList<IBond>();

        for(IBond bond : molecule.bonds())
            if(bond instanceof QueryBond && ((QueryBond) bond).getExpression().type() == Expr.Type.TRUE)
                otherBonds.add(bond);

        for(IBond bond : otherBonds)
            molecule.removeBond(bond);

        configureAromaticity(molecule, AromaticityMode.AUTO);

        if(inchiStereo)
            setStereo(molecule, new InChITools(molecule, INCHI_DB_MODE, false));

        for(IBond bond : otherBonds)
            molecule.addBond(bond);

        return molecule;
    }


    private static IAtomContainer getMoleculeFromMolfile(String mol) throws CDKException, IOException
    {
        DefaultChemObjectReader mdlReader;

        if(mol.contains("M  V30 BEGIN CTAB"))
        {
            mdlReader = new MDLV3000Reader();
        }
        else
        {
            mdlReader = new MDLV2000Reader();
            mdlReader.getSetting("AddStereoElements").setSetting("false");
        }

        mdlReader.setReader(new StringReader(mol));
        IAtomContainer molecule;


        try
        {
            molecule = mdlReader.read(new AtomContainer());
        }
        catch(Exception e)
        {
            throw new CDKException("cannot parse molfile: " + e.getMessage());
        }
        finally
        {
            mdlReader.close();
        }

        sanitizeMolecule(molecule);

        return molecule;
    }


    private static IAtomContainer getMoleculeFromSmiles(String smiles) throws CDKException
    {
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule;

        try
        {
            molecule = sp.parseSmiles(smiles);
        }
        catch(InvalidSmilesException e)
        {
            sp.kekulise(false);
            molecule = sp.parseSmiles(smiles);
        }

        sanitizeMolecule(molecule);

        molecule.setTitle(smiles);

        return molecule;
    }


    private static List<IAtomContainer> getMoleculesFromRGroupQuery(String rgroup) throws CDKException, IOException
    {
        try(RGroupQueryReader reader = new RGroupQueryReader(new ByteArrayInputStream(rgroup.getBytes())))
        {
            RGroupQuery rGroupQuery = reader.read(new RGroupQuery(SilentChemObjectBuilder.getInstance()));
            List<IAtomContainer> molecules = rGroupQuery.getAllConfigurations();

            for(IAtomContainer molecule : molecules)
                sanitizeMolecule(molecule);

            return molecules;
        }
    }


    private static void sanitizeMolecule(IAtomContainer molecule) throws CDKException
    {
        IChemObjectBuilder builder = molecule.getBuilder();

        for(int i = 0; i < molecule.getAtomCount(); i++)
        {
            IAtom atom = molecule.getAtom(i);

            if(atom instanceof IPseudoAtom)
            {
                atom.setAtomicNumber(0);

                String symbol = ((IPseudoAtom) atom).getLabel();

                if(symbol.equals("D") || symbol.equals("T"))
                {
                    IIsotope isotope = builder.newInstance(IIsotope.class, "H", symbol.equals("D") ? 2 : 3);
                    IAtom hydrogen = builder.newInstance(IAtom.class, isotope);
                    hydrogen.setFormalCharge(atom.getFormalCharge());
                    hydrogen.setPoint3d(atom.getPoint3d());
                    hydrogen.setPoint2d(atom.getPoint2d());

                    molecule.setAtom(i, hydrogen);
                }
                else if(symbol.equals("Mu-"))
                {
                    ((IPseudoAtom) atom).setLabel("Mu");
                    atom.setCharge(-1.0);
                }
            }
        }


        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(molecule.getBuilder());

        for(IAtom atom : molecule.atoms())
        {
            if(atom.getImplicitHydrogenCount() == null)
            {
                if(molecule.getConnectedBondsList(atom).stream().noneMatch(b -> b instanceof QueryBond))
                {
                    IAtomType type = matcher.findMatchingAtomType(molecule, atom);
                    AtomTypeManipulator.configure(atom, type);
                    CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
                    adder.addImplicitHydrogens(molecule, atom);
                }
                else
                {
                    atom.setImplicitHydrogenCount(0);
                }
            }
        }

        kekulize(molecule);

        molecule.setAtoms(BinaryMoleculeSort.atomsByFrequency(molecule));
    }


    private static void configureAromaticity(IAtomContainer molecule, AromaticityMode aromaticityMode)
            throws CDKException
    {
        for(IBond bond : molecule.bonds())
            if(bond instanceof IQueryBond || bond.getOrder() == Order.UNSET)
                return;

        switch(aromaticityMode)
        {
            case PRESERVE:
                return;

            case AUTO:
                if(molecule.getFlag(CDKConstants.ISAROMATIC))
                    return;

                for(IBond bond : molecule.bonds())
                    if(bond.isAromatic())
                        return;

            case DETECT:
                aromaticity.get().apply(molecule);
        }
    }


    private static void kekulize(IAtomContainer molecule)
    {
        List<IBond> unset = new LinkedList<IBond>();

        for(IBond bond : molecule.bonds())
            if(bond.getOrder() == Order.UNSET && bond.isAromatic())
                unset.add(bond);

        if(!unset.isEmpty())
        {
            try
            {
                Kekulization.kekulize(molecule);
            }
            catch(CDKException e)
            {
                for(IBond bond : unset)
                    bond.setOrder(Order.UNSET);
            }
        }
    }


    private static void setStereo(IAtomContainer molecule, InChITools inchi)
    {
        molecule.setStereoElements(inchi.getStereoElements());

        for(IAtom atom : inchi.getStereoAtoms())
            atom.setProperty(STEREO_FLAG, Boolean.TRUE);

        for(IBond bond : inchi.getStereoBonds())
            if(!bond.isAromatic())
                bond.setProperty(STEREO_FLAG, Boolean.TRUE);
    }
}
