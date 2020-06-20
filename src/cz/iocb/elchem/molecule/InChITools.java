package cz.iocb.elchem.molecule;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry.Conformation;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.interfaces.ITetrahedralChirality.Stereo;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.ExtendedCisTrans;
import org.openscience.cdk.stereo.ExtendedTetrahedral;
import org.openscience.cdk.stereo.TetrahedralChirality;
import cz.iocb.elchem.molecule.Molecule.AtomType;



public class InChITools
{
    @SuppressWarnings("serial")
    public static class InChIException extends CDKException
    {
        public InChIException(String message)
        {
            super(message);
        }
    }


    public static class TautomericGroup
    {
        public int mobiles;
        public int charges;
        public List<Integer> atoms;
    }


    private static abstract class Field
    {
        private static final int elname = 0;
        private static final int el_number = 6;
        private static final int neighbor = 8;
        private static final int orig_at_number = 48;
        private static final int bond_stereo = 52;
        private static final int bond_type = 72;
        private static final int valence = 92;
        private static final int chem_bonds_valence = 93;
        private static final int num_H = 94;
        private static final int iso_atw_diff = 98;
        private static final int charge = 99;
        private static final int radical = 100;
        private static final int x = 112;
        private static final int y = 120;
        private static final int z = 128;
    }


    private static final int ZERO_ATW_DIFF = 127;
    private static final int MAXVAL = 20;
    private static final int RECORD_SIZE = 176;
    private static final int STEREO_RECORD_SIZE = 12;
    private static final short NO_ATOM = -1;
    private static Isotopes isotopes;

    private final IAtomContainer molecule;

    @SuppressWarnings("rawtypes") private final List<IStereoElement> stereo = new ArrayList<IStereoElement>();
    private final Set<IAtom> stereoAtoms = new HashSet<IAtom>();
    private final Set<IBond> stereoBonds = new HashSet<IBond>();

    private final List<TautomericGroup> tautomericGroups = new ArrayList<TautomericGroup>();
    private final List<Integer> tautomericBonds = new ArrayList<Integer>();


    static
    {
        try
        {
            isotopes = Isotopes.getInstance();
        }
        catch(IOException e)
        {
        }
    }


    public InChITools(IAtomContainer molecule, int mode, boolean tautomers) throws CDKException
    {
        this.molecule = molecule;

        boolean is0D = true;

        for(IAtom atom : molecule.atoms())
        {
            if(atom.getPoint3d() != null || atom.getPoint2d() != null)
            {
                is0D = false;
                break;
            }
        }

        if(process(moleculeAsBytes(molecule), is0D ? stereoAsBytes(molecule) : null, mode, tautomers) < 0)
            throw new InChIException("inchi generation failed");
    }


    private native int process(byte[] molecule, byte[] stereo, int mode, boolean tautomerism);


    void setStereoAtoms(short[] data)
    {
        for(int i = 0; i < data.length;)
        {
            final IAtom atom = molecule.getAtom(data[i++] - 1);
            final int parity = data[i++];

            if(parity == -1 || parity == -2)
            {
                List<IBond> bonds = molecule.getConnectedBondsList(atom);
                IAtom[] peripherals = new IAtom[4];
                int position = 0;

                for(int b = 0; b < 2; b++)
                {
                    IAtom terminal = getTerminalAtom(atom, bonds.get(b));
                    List<IAtom> atoms = new ArrayList<IAtom>(2);

                    for(IBond bond : molecule.getConnectedBondsList(terminal))
                        if(bond.getOrder() != Order.DOUBLE)
                            atoms.add(bond.getOther(terminal));

                    if(atoms.size() == 1)
                        peripherals[position++] = terminal;

                    for(int diff = 0; diff < 3; diff++)
                        for(IAtom a : atoms)
                            if(a.getAtomicNumber() == AtomType.H && getMassDiff(a) == diff)
                                peripherals[position++] = a;

                    for(IAtom a : atoms)
                        if(!isInchiH(a))
                            peripherals[position++] = a;
                }

                stereo.add(new ExtendedTetrahedral(atom, peripherals, parity == -1 ?
                        ITetrahedralChirality.Stereo.CLOCKWISE : ITetrahedralChirality.Stereo.ANTI_CLOCKWISE));
            }
            else if(parity == 1 || parity == 2)
            {
                IAtom[] ligands = new IAtom[4];

                List<IAtom> neighbours = molecule.getConnectedAtomsList(atom);

                int position = 0;

                if(neighbours.size() == 3)
                    ligands[position++] = atom;

                for(int diff = 0; diff < 3; diff++)
                    for(IAtom a : neighbours)
                        if(a.getAtomicNumber() == AtomType.H && getMassDiff(a) == diff)
                            ligands[position++] = a;

                for(IAtom a : neighbours)
                    if(!isInchiH(a))
                        ligands[position++] = a;

                stereo.add(new TetrahedralChirality(atom, ligands,
                        parity == 1 ? Stereo.ANTI_CLOCKWISE : Stereo.CLOCKWISE));
            }

            stereoAtoms.add(atom);
        }
    }


    void setStereoBonds(short[] data)
    {
        for(int i = 0; i < data.length; i++)
        {
            final IBond bond = molecule.getBond(molecule.getAtom(data[i++] - 1), molecule.getAtom(data[i++] - 1));
            final int parity = data[i];

            if(parity == -1 || parity == -2)
            {
                IBond bonds[] = new IBond[2];

                for(int j = 0; j < 2; j++)
                {
                    IAtom terminal = getTerminalAtom(bond.getAtom(j), bond);
                    List<IBond> bondsList = molecule.getConnectedBondsList(terminal);

                    for(Iterator<IBond> it = bondsList.iterator(); it.hasNext();)
                        if(it.next().getOrder() == Order.DOUBLE)
                            it.remove();

                    IAtom other0 = bondsList.get(0).getOther(terminal);
                    IAtom other1 = bondsList.size() > 1 ? bondsList.get(1).getOther(terminal) : null;

                    if(other1 == null)
                        bonds[j] = bondsList.get(0);
                    else if(!isInchiH(other1))
                        bonds[j] = bondsList.get(1);
                    else if(!isInchiH(other0))
                        bonds[j] = bondsList.get(0);
                    else if(getMassDiff(other1) > getMassDiff(other0))
                        bonds[j] = bondsList.get(1);
                    else
                        bonds[j] = bondsList.get(0);
                }

                stereo.add(new ExtendedCisTrans(bond, bonds,
                        parity == -1 ? IStereoElement.TOGETHER : IStereoElement.OPPOSITE));
            }
            else if(parity == 1 || parity == 2)
            {
                IBond bonds[] = new IBond[2];

                for(int j = 0; j < 2; j++)
                {
                    IAtom atom = bond.getAtom(j);
                    List<IBond> bondsList = molecule.getConnectedBondsList(atom);
                    bondsList.remove(bond);

                    IAtom other0 = bondsList.get(0).getOther(atom);
                    IAtom other1 = bondsList.size() > 1 ? bondsList.get(1).getOther(atom) : null;

                    if(other1 == null)
                        bonds[j] = bondsList.get(0);
                    else if(!isInchiH(other1))
                        bonds[j] = bondsList.get(1);
                    else if(!isInchiH(other0))
                        bonds[j] = bondsList.get(0);
                    else if(getMassDiff(other1) > getMassDiff(other0))
                        bonds[j] = bondsList.get(1);
                    else
                        bonds[j] = bondsList.get(0);
                }

                stereo.add(new DoubleBondStereochemistry(bond, bonds,
                        parity == 2 ? Conformation.OPPOSITE : Conformation.TOGETHER));
            }

            stereoBonds.add(bond);
        }
    }


    void setAlternatingBonds(short[] data)
    {
        for(int i = 0; i < data.length; i++)
            tautomericBonds.add(
                    molecule.indexOf(molecule.getBond(molecule.getAtom(data[i++] - 1), molecule.getAtom(data[i] - 1))));
    }


    void setTautomericGroup(short[] data)
    {
        TautomericGroup group = new TautomericGroup();
        tautomericGroups.add(group);

        group.mobiles = data[0];
        group.charges = data[1];
        group.atoms = new ArrayList<Integer>(data.length - 2);


        for(int i = 2; i < data.length; i++)
            group.atoms.add(data[i] - 1);
    }


    private final IAtom getTerminalAtom(IAtom atom, IBond bond)
    {
        while(true)
        {
            List<IBond> bonds = molecule.getConnectedBondsList(atom);
            bonds.remove(bond);

            if(bonds.size() != 1)
                return atom;

            bond = bonds.get(0);

            if(bond.getOrder() != Order.DOUBLE)
                return atom;

            atom = bond.getOther(atom);
        }
    }


    private static final byte[] stereoAsBytes(IAtomContainer molecule) throws CDKException
    {
        ArrayList<IStereoElement<?, ?>> elements = new ArrayList<IStereoElement<?, ?>>();

        for(IStereoElement<?, ?> e : molecule.stereoElements())
            elements.add(e);

        byte array[] = new byte[STEREO_RECORD_SIZE * elements.size()];
        ByteBuffer buffer = ByteBuffer.wrap(array);
        buffer.order(ByteOrder.nativeOrder());

        for(IStereoElement<?, ?> element : elements)
        {
            if(element instanceof TetrahedralChirality)
            {
                TetrahedralChirality e = (TetrahedralChirality) element;

                int implicitH = -1;

                for(int i = 0; i < 4; i++)
                    if(e.getLigands()[i] == e.getFocus())
                        implicitH = i;

                if(implicitH != -1)
                    buffer.putShort((short) molecule.indexOf(e.getFocus()));

                for(int i = 0; i < 4; i++)
                    if(e.getLigands()[i] != e.getFocus())
                        buffer.putShort((short) molecule.indexOf(e.getLigands()[i]));

                buffer.putShort((short) molecule.indexOf(e.getFocus()));

                buffer.put(/* INCHI_StereoType_Tetrahedral */ (byte) 2);

                if(implicitH == -1 || implicitH % 2 == 0)
                    buffer.put(e.getStereo() == Stereo.CLOCKWISE ? (byte) 2 : (byte) 1);
                else
                    buffer.put(e.getStereo() == Stereo.CLOCKWISE ? (byte) 1 : (byte) 2);
            }
            else if(element instanceof DoubleBondStereochemistry)
            {
                DoubleBondStereochemistry e = (DoubleBondStereochemistry) element;

                IBond bond = e.getStereoBond();
                IBond[] bonds = e.getBonds();

                if(!bonds[0].contains(bond.getAtom(0)))
                    bonds = new IBond[] { bonds[1], bonds[0] };

                buffer.putShort((short) molecule.indexOf(bonds[0].getOther(bond.getAtom(0))));
                buffer.putShort((short) molecule.indexOf(bond.getAtom(0)));
                buffer.putShort((short) molecule.indexOf(bond.getAtom(1)));
                buffer.putShort((short) molecule.indexOf(bonds[1].getOther(bond.getAtom(1))));

                buffer.putShort(NO_ATOM);

                buffer.put(/* INCHI_StereoType_DoubleBond */ (byte) 1);
                buffer.put(e.getStereo() == Conformation.OPPOSITE ? (byte) 2 : (byte) 1);
            }
            else if(element instanceof ExtendedTetrahedral)
            {
                ExtendedTetrahedral e = (ExtendedTetrahedral) element;

                IAtom[] terminals = e.findTerminalAtoms(molecule);
                IAtom[] ligands = e.getCarriers().toArray(new IAtom[0]);

                if((terminals[0] != ligands[0] && molecule.getBond(terminals[0], ligands[0]) == null)
                        || (terminals[0] != ligands[1] && molecule.getBond(terminals[0], ligands[1]) == null))
                    terminals = new IAtom[] { terminals[1], terminals[0] };

                buffer.putShort((short) molecule.indexOf(terminals[0] != ligands[0] ? ligands[0] : ligands[1]));
                buffer.putShort((short) molecule.indexOf(terminals[0]));
                buffer.putShort((short) molecule.indexOf(terminals[1]));
                buffer.putShort((short) molecule.indexOf(terminals[1] != ligands[2] ? ligands[2] : ligands[3]));

                buffer.putShort((short) molecule.indexOf(e.focus()));

                buffer.put(/* INCHI_StereoType_Allene */ (byte) 3);

                if(terminals[0] == ligands[0] ^ terminals[1] == ligands[2])
                    buffer.put(e.winding() == Stereo.CLOCKWISE ? (byte) 2 : (byte) 1);
                else
                    buffer.put(e.winding() == Stereo.CLOCKWISE ? (byte) 1 : (byte) 2);
            }
            else if(element instanceof ExtendedCisTrans)
            {
                ExtendedCisTrans e = (ExtendedCisTrans) element;

                IAtom[] terminals = ExtendedCisTrans.findTerminalAtoms(molecule, e.getFocus());
                IBond[] bonds = e.getCarriers().toArray(new IBond[0]);

                if(!bonds[0].contains(terminals[0]))
                    bonds = new IBond[] { bonds[1], bonds[0] };

                buffer.putShort((short) molecule.indexOf(bonds[0].getOther(terminals[0])));
                buffer.putShort((short) molecule.indexOf(terminals[0]));
                buffer.putShort((short) molecule.indexOf(terminals[1]));
                buffer.putShort((short) molecule.indexOf(bonds[1].getOther(terminals[1])));

                buffer.putShort(NO_ATOM);

                buffer.put(/* INCHI_StereoType_DoubleBond */ (byte) 1);
                buffer.put(e.getConfigOrder() == IStereoElement.OPPOSITE ? (byte) 2 : (byte) 1);
            }
        }

        return array;
    }


    private static final byte[] moleculeAsBytes(IAtomContainer molecule) throws CDKException
    {
        if(molecule.getAtomCount() >= Short.MAX_VALUE)
            throw new InChIException("too many atoms");

        int pseudoAtomId = Byte.MIN_VALUE;

        byte array[] = new byte[RECORD_SIZE * molecule.getAtomCount()];
        ByteBuffer buffer = ByteBuffer.wrap(array);
        buffer.order(ByteOrder.nativeOrder());

        for(int index = 0; index < molecule.getAtomCount(); index++)
        {
            int offset = index * RECORD_SIZE;
            IAtom atom = molecule.getAtom(index);
            List<IBond> bonds = molecule.getConnectedBondsList(atom);
            int multiplicity = molecule.getConnectedSingleElectronsCount(atom);

            int valence = 0;
            boolean hasKnownBonds = true;

            for(Iterator<IBond> it = bonds.iterator(); it.hasNext();)
            {
                IBond b = it.next();
                valence += getBondType(b);

                if(b.getOrder() != Order.SINGLE && b.getOrder() != Order.DOUBLE && b.getOrder() != Order.TRIPLE)
                {
                    it.remove();
                    hasKnownBonds = false;
                }
            }


            if(atom.getImplicitHydrogenCount() + bonds.size() > MAXVAL)
                throw new InChIException("too large number of bonds");


            buffer.putShort(offset + Field.orig_at_number, (short) (index + 1));
            buffer.put(offset + Field.valence, (byte) bonds.size());
            buffer.put(offset + Field.chem_bonds_valence, (byte) valence);
            buffer.put(offset + Field.num_H, (byte) (int) atom.getImplicitHydrogenCount());
            buffer.put(offset + Field.charge, (byte) (int) atom.getFormalCharge());
            buffer.put(offset + Field.radical, (byte) (multiplicity > 0 ? multiplicity + 1 : 0));

            if(!(atom instanceof IPseudoAtom) && hasKnownBonds)
            {
                byte[] name = atom.getSymbol().getBytes();

                for(int i = 0; i < name.length; i++)
                    buffer.put(offset + Field.elname + i, name[i]);

                int iso = getMassDiff(atom);
                buffer.put(offset + Field.iso_atw_diff, (byte) (iso == ZERO_ATW_DIFF ? 1 : iso > 0 ? iso + 1 : iso));
                buffer.put(offset + Field.el_number, (byte) (int) atom.getAtomicNumber());
            }
            else if(pseudoAtomId < Byte.MAX_VALUE)
            {
                buffer.put(offset + Field.iso_atw_diff, (byte) (pseudoAtomId++));
            }
            else
            {
                throw new InChIException("too many pseudo atoms or query bonds");
            }

            if(atom.getPoint3d() != null)
            {
                buffer.putDouble(offset + Field.x, atom.getPoint3d().x);
                buffer.putDouble(offset + Field.y, atom.getPoint3d().y);
                buffer.putDouble(offset + Field.z, atom.getPoint3d().z);
            }
            else if(atom.getPoint2d() != null)
            {
                buffer.putDouble(offset + Field.x, atom.getPoint2d().x);
                buffer.putDouble(offset + Field.y, atom.getPoint2d().y);
            }

            for(int i = 0; i < bonds.size(); i++)
            {
                IBond bond = bonds.get(i);

                buffer.putShort(offset + Field.neighbor + 2 * i, (short) molecule.indexOf(bond.getOther(atom)));
                buffer.put(offset + Field.bond_type + i, (byte) getBondType(bond));
                buffer.put(offset + Field.bond_stereo + i, (byte) getStereo(atom, bond));
            }
        }

        return array;
    }


    private static final int getMassDiff(IAtom atom)
    {
        if(atom.getMassNumber() == null || atom.getMassNumber() == -1)
            return 0;

        if(atom.getAtomicNumber() == 0)
            return 0;

        IIsotope majorIsotope = isotopes.getMajorIsotope(atom.getAtomicNumber());

        if(majorIsotope == null)
            return 0;

        return atom.getMassNumber() - majorIsotope.getMassNumber();
    }


    private static final int getBondType(IBond bond) throws CDKException
    {
        if(bond.getOrder() == Order.SINGLE)
            return 1;
        else if(bond.getOrder() == Order.DOUBLE)
            return 2;
        else if(bond.getOrder() == Order.TRIPLE)
            return 3;
        else
            return 0;
    }


    private static final int getStereo(IAtom atom, IBond bond)
    {
        int multiplier = bond.getAtom(0) == atom ? 1 : -1;

        switch(bond.getStereo())
        {
            case UP:
                return multiplier * 1;

            case UP_INVERTED:
                return multiplier * -1;

            case UP_OR_DOWN:
                return multiplier * 4;

            case UP_OR_DOWN_INVERTED:
                return multiplier * -4;

            case DOWN:
                return multiplier * 6;

            case DOWN_INVERTED:
                return multiplier * -6;

            case E_OR_Z:
                return 3;

            default:
                return 0;
        }
    }


    private static final boolean isInchiH(IAtom atom)
    {
        return atom.getAtomicNumber() == AtomType.H && getMassDiff(atom) < 3 && getMassDiff(atom) >= 0;
    }


    @SuppressWarnings("rawtypes")
    public List<IStereoElement> getStereoElements()
    {
        return stereo;
    }


    public Set<IAtom> getStereoAtoms()
    {
        return stereoAtoms;
    }


    public Set<IBond> getStereoBonds()
    {
        return stereoBonds;
    }


    public List<TautomericGroup> getTautomericGroups()
    {
        return tautomericGroups;
    }


    public List<Integer> getTautomericBonds()
    {
        return tautomericBonds;
    }


    private static native void init();


    static
    {
        init();
    }
}
