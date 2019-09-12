/*
 * Copyright (C) 2015-2019 Jakub Galgonek   galgonek@uochb.cas.cz
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
import java.security.AccessController;
import java.security.PrivilegedActionException;
import java.security.PrivilegedExceptionAction;
import java.util.ArrayList;
import java.util.Arrays;
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
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.DefaultChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.RGroupQueryReader;
import org.openscience.cdk.io.formats.RGroupQueryFormat;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.isomorphism.matchers.RGroupQuery;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import cz.iocb.elchem.tautomers.InchiTautomerGenerator;



public class MoleculeCreator
{
    private static final ThreadLocal<Aromaticity> aromaticity = new ThreadLocal<Aromaticity>()
    {
        @Override
        protected Aromaticity initialValue()
        {
            return new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), Cycles.cdkAromaticSet());
        }
    };


    public static List<IAtomContainer> translateQuery(String query, AromaticityMode aromaticityMode,
            TautomerMode tautomerMode) throws CDKException, IOException, TimeoutException
    {
        try
        {
            List<IAtomContainer> queries = null;
            List<String> lines = Arrays.asList(query.split("\\n"));


            if(((RGroupQueryFormat) RGroupQueryFormat.getInstance()).matches(lines).matched())
                queries = getMoleculesFromRGroupQuery(query, aromaticityMode);
            else if(lines.size() > 1)
                queries = Arrays.asList(MoleculeCreator.getMoleculeFromMolfile(query, aromaticityMode));
            else
                queries = Arrays.asList(MoleculeCreator.getMoleculeFromSmiles(query, aromaticityMode));


            if(tautomerMode == TautomerMode.INCHI)
            {
                InchiTautomerGenerator generator = new InchiTautomerGenerator();
                List<IAtomContainer> tautomers = new ArrayList<IAtomContainer>();

                for(IAtomContainer molecule : queries)
                    tautomers.addAll(generator.getTautomers(molecule));

                for(IAtomContainer tautomer : tautomers)
                    tautomer.setAtoms(BinaryMoleculeSort.atomsByFrequency(tautomer));

                queries = tautomers;
            }


            return queries;
        }
        catch(CloneNotSupportedException e)
        {
            return null;
        }
    }


    public static IAtomContainer getMoleculeFromMolfile(String mol, AromaticityMode aromaticityMode)
            throws CDKException, IOException
    {
        DefaultChemObjectReader mdlReader;

        if(mol.contains("M  V30 BEGIN CTAB"))
            mdlReader = new MDLV3000Reader();
        else
            mdlReader = new MDLV2000Reader();

        mdlReader.setReader(new StringReader(mol));
        IAtomContainer molecule;


        try
        {
            molecule = AccessController.doPrivileged((PrivilegedExceptionAction<IAtomContainer>) () -> {
                return mdlReader.read(new AtomContainer());
            });
        }
        catch(PrivilegedActionException e)
        {
            throw(CDKException) e.getException();
        }
        catch(Exception e)
        {
            throw new CDKException("cannot parse molfile: " + e.getMessage());
        }
        finally
        {
            mdlReader.close();
        }

        IChemObjectBuilder builder = molecule.getBuilder();

        for(int i = 0; i < molecule.getAtomCount(); i++)
        {
            IAtom atom = molecule.getAtom(i);

            if(atom instanceof IPseudoAtom)
            {
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


        MoleculeCreator.configureMolecule(molecule, aromaticityMode);
        molecule.setAtoms(BinaryMoleculeSort.atomsByFrequency(molecule));

        return molecule;
    }


    public static IAtomContainer getMoleculeFromSmiles(String smiles, AromaticityMode aromaticityMode)
            throws CDKException
    {
        try
        {
            IAtomContainer molecule = AccessController.doPrivileged((PrivilegedExceptionAction<IAtomContainer>) () -> {
                SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

                try
                {
                    return sp.parseSmiles(smiles);
                }
                catch(InvalidSmilesException e)
                {
                    sp.kekulise(false);
                    return sp.parseSmiles(smiles);
                }
            });


            MoleculeCreator.configureMolecule(molecule, aromaticityMode);
            molecule.setAtoms(BinaryMoleculeSort.atomsByFrequency(molecule));

            return molecule;
        }
        catch(PrivilegedActionException e)
        {
            throw(CDKException) e.getException();
        }
    }


    public static List<IAtomContainer> getMoleculesFromRGroupQuery(String rgroup, AromaticityMode aromaticityMode)
            throws CDKException, IOException
    {
        try(RGroupQueryReader reader = new RGroupQueryReader(new ByteArrayInputStream(rgroup.getBytes())))
        {
            RGroupQuery rGroupQuery = reader.read(new RGroupQuery(SilentChemObjectBuilder.getInstance()));
            List<IAtomContainer> molecules = rGroupQuery.getAllConfigurations();

            for(IAtomContainer molecule : molecules)
                configureMolecule(molecule, aromaticityMode);

            return molecules;
        }
    }


    private static void configureMolecule(IAtomContainer molecule, AromaticityMode aromaticityMode) throws CDKException
    {
        kekulize(molecule);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        configureAromaticity(molecule, aromaticityMode);
    }


    private static void configureAromaticity(IAtomContainer molecule, AromaticityMode aromaticityMode)
            throws CDKException
    {
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
}
