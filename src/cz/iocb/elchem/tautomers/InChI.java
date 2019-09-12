/*
 * Copyright (C) 2015-2019 Jakub Galgonek   galgonek@uochb.cas.cz
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

import java.security.AccessController;
import java.security.PrivilegedActionException;
import java.security.PrivilegedExceptionAction;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import net.sf.jniinchi.INCHI_OPTION;



public class InChI
{
    public static class Fragment
    {
        public String formula;
        public String connections;
        public String hydrogens;
        public String formalCharge;
        public String doubleBondStereo;
        public String tetrahedralStereo;
        String aux;


        void setPart(char name, String part)
        {
            switch(name)
            {
                case 'c':
                    connections = part;
                    break;
                case 'h':
                    hydrogens = part;
                    break;
                case 'q':
                    formalCharge = part;
                    break;
                case 'b':
                    doubleBondStereo = part;
                    break;
                case 't':
                    tetrahedralStereo = part;
                    break;
            }
        }
    }


    static private final List<INCHI_OPTION> options;

    private String aux = null;
    private String value = null;


    static
    {
        options = new ArrayList<INCHI_OPTION>();
        options.add(INCHI_OPTION.SUU);
        options.add(INCHI_OPTION.RecMet);
    }


    public InChI(IAtomContainer molecule) throws CDKException
    {
        try
        {
            InChIGenerator generator = AccessController.doPrivileged((PrivilegedExceptionAction<InChIGenerator>) () -> {
                return InChIGeneratorFactory.getInstance().getInChIGenerator(molecule, options);
            });

            value = generator.getInchi();
            aux = generator.getAuxInfo();

            // workaround
            if(value != null && value.contains("/r"))
            {
                value = "InChI=1/" + value.substring(value.indexOf("/r") + 2);
                aux = "AuxInfo=1/" + aux.substring(aux.indexOf("/R:") + 4);
            }

            if(value == null || aux == null)
            {
                if(generator.getMessage() != null)
                    throw new InChIException("cannot generate InChI: " + generator.getMessage());

                throw new InChIException("cannot generate InChI");
            }
        }
        catch(PrivilegedActionException e)
        {
            throw(CDKException) e.getException();
        }
        catch(IllegalArgumentException e)
        {
            throw new InChIException("cannot generate InChI");
        }
    }


    public List<Fragment> decompose()
    {
        LinkedList<Fragment> fragments = new LinkedList<Fragment>();

        if(!aux.contains("/N:"))
            return fragments;

        String[] components = value.split("/");

        String formulaComponent = components[1];
        String[] formulas = formulaComponent.split("\\.");

        for(String formula : formulas)
        {
            Pattern formulaPattern = Pattern.compile("^[0-9]+");
            Matcher match = formulaPattern.matcher(formula);

            if(!match.find())
            {
                Fragment fragment = new Fragment();
                fragment.formula = formula;
                fragments.add(fragment);
            }
            else
            {
                String countString = match.group();
                String value = formula.substring(countString.length());
                int count = Integer.valueOf(countString);

                for(int i = 0; i < count; i++)
                {
                    Fragment fragment = new Fragment();
                    fragment.formula = value;
                    fragments.add(fragment);
                }
            }
        }


        boolean hasIsomericPart = false;

        for(int i = 2; i < components.length; i++)
        {
            char name = components[i].charAt(0);

            if(name == 'p' || name == 'h' && hasIsomericPart)
                continue;

            if(name == 'i')
                hasIsomericPart = true;

            String[] parts = components[i].substring(1).split(";");
            int index = 0;

            for(String part : parts)
            {
                Pattern formulaPattern = Pattern.compile("^[0-9]+\\*");
                Matcher match = formulaPattern.matcher(part);

                if(!match.find())
                {
                    if(!part.isEmpty())
                        fragments.get(index).setPart(name, part);

                    index++;
                }
                else
                {
                    String countString = match.group();
                    String value = part.substring(countString.length());
                    int count = Integer.valueOf(countString.substring(0, countString.length() - 1));

                    for(int j = 0; j < count; j++)
                    {
                        if(!value.isEmpty())
                            fragments.get(index).setPart(name, value);

                        index++;
                    }
                }
            }
        }


        String[] auxComponents = aux.split("/");

        if(auxComponents[2].startsWith("N:"))
        {
            String[] values = auxComponents[2].substring(2).split(";");

            int index = 0;
            for(String value : values)
                fragments.get(index++).aux = value;
        }

        return fragments;
    }


    public String getAuxInfo()
    {
        return aux;
    }


    public String getValue()
    {
        return value;
    }
}
