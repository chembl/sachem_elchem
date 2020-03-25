package cz.iocb.elchem.molecule;



public abstract class Molecule
{
    public static class SGroup
    {
        byte type;
        byte subtype;
        byte connectivity;

        int atoms[];
        int bonds[][];
    }


    public static class AtomLabel
    {
        int atom;
        byte[] label;
    }


    public static abstract class AtomType
    {
        public static final byte H = 1;
        public static final byte C = 6;
        public static final byte N = 7;
        public static final byte O = 8;
        public static final byte F = 9;
        public static final byte P = 15;
        public static final byte S = 16;
        public static final byte Cl = 17;
        public static final byte Br = 35;
        public static final byte I = 53;

        public static final byte R = -'R';
        public static final byte Q = -'Q';
        public static final byte M = -'M';
        public static final byte X = -'X';

        public static final byte A = -'A';
        public static final byte G = -'G';

        public static final byte POSITRONIUM = -'p';
        public static final byte ELECTRON = -'e';
        public static final byte PHOTON = -'h';
        public static final byte MUONIUM = -'m';

        public static final byte ENZYME = -'z';
        public static final byte ACP = -'a';

        public static final byte EMPTY = -' ';
        public static final byte UNKNOWN = -'?';
    }


    public static abstract class BondType
    {
        public static final byte NONE = 0;
        public static final byte SINGLE = 1;
        public static final byte DOUBLE = 2;
        public static final byte TRIPLE = 3;
        public static final byte QUADRUPLE = 4;
        public static final byte QUINTUPLE = 5;
        public static final byte SEXTUPLE = 6;
        public static final byte AROMATIC = 11;
        public static final byte SINGLE_OR_DOUBLE = 12;
        public static final byte SINGLE_OR_AROMATIC = 13;
        public static final byte DOUBLE_OR_AROMATIC = 14;
        public static final byte ANY = 15;
    }


    public static abstract class TetrahedralStereo
    {
        public static final byte NONE = 0;
        public static final byte CLOCKWISE = 1;
        public static final byte ANTI_CLOCKWISE = 2;
        public static final byte UNDEFINED = 3;
    }


    public static abstract class BondStereo
    {
        public static final byte NONE = 0;
        public static final byte OPPOSITE = 1;
        public static final byte TOGETHER = 2;
        public static final byte UNDEFINED = 3;
    }


    public static abstract class SgroupType
    {
        public static final byte NONE = 0;
        public static final byte SRU = 1;
        public static final byte MOD = 2;
        public static final byte MON = 3;
        public static final byte COP = 4;
        public static final byte GEN = 5;
        public static final byte ANY = 6;
        public static final byte CRO = 7;
        public static final byte MER = 8;
        public static final byte GRA = 9;
        public static final byte COM = 10;
        public static final byte FOR = 11;
        public static final byte MIX = 12;
    }


    public static abstract class SgroupSubtype
    {
        public static final byte NONE = 0;
        public static final byte ALT = 1;
        public static final byte RAN = 2;
        public static final byte BLO = 3;
        public static final byte UNKNOWN = -1;
    }


    public static abstract class SgroupConnectivity
    {
        public static final byte NONE = 0;
        public static final byte HH = 1;
        public static final byte HT = 2;
        public static final byte EU = 3;
        public static final byte UNKNOWN = -1;
    }


    public static final int MAX_ATOM_IDX = Short.MAX_VALUE;


    public abstract int getOriginalAtomCount();

    public abstract int getOriginalBondCount();

    public abstract int getAtomCount();

    public abstract int getBondCount();

    public abstract boolean hasPseudoAtom();

    public abstract boolean hasRestHydrogenFlags();

    public abstract byte getAtomNumber(int atom);

    public abstract AtomLabel getAtomLabel(int atom);

    public abstract byte getAtomHydrogenCount(int atom);

    public abstract byte getAtomFormalCharge(int atom);

    public abstract byte getAtomMass(int atom);

    public abstract byte getAtomRadicalType(int atom);

    public abstract byte getAtomStereo(int atom);

    public abstract boolean getAtomRestHydrogenFlag(int atom);

    public abstract int getBond(int atom0, int atom1);

    public abstract byte getBondType(int bond);

    public abstract byte getBondStereo(int atom);

    public abstract int getBondAtom(int bond, int atom);

    public abstract boolean isAtomInBond(int bond, int atom);

    public abstract int[] getBondedAtoms(int atom);

    public abstract SGroup[] getSGroups();


    public final boolean isAtomPseudo(int atom)
    {
        return getAtomNumber(atom) < 0;
    }


    public final boolean isAtomMetal(int atom)
    {
        byte number = getAtomNumber(atom);

        return (number > 2 && number < 5) || (number > 10 && number < 14) || (number > 18 && number < 32)
                || (number > 36 && number < 51) || (number > 54 && number < 85) || number > 86;
    }


    public final boolean isAtomHalogen(int atom)
    {
        byte number = getAtomNumber(atom);

        return number == 9 || number == 17 || number == 35 || number == 53 || number == 85;
    }


    public final boolean isQueryBond(int bond)
    {
        return getBondType(bond) >= BondType.SINGLE_OR_DOUBLE;
    }


    public final int getOtherBondAtom(int bond, int atom)
    {
        if(getBondAtom(bond, 0) == atom)
            return getBondAtom(bond, 1);
        else if(getBondAtom(bond, 1) == atom)
            return getBondAtom(bond, 0);
        else
            return -1;
    }


    public final int getOppositeAtom(int centre, int atom)
    {
        int[] bonded = getBondedAtoms(centre);

        if(bonded.length != 2)
            return -1;

        if(bonded[0] == atom)
            return bonded[1];
        else if(bonded[1] == atom)
            return bonded[0];
        else
            return -1;
    }


    public final int getLastStereoBondLigand(int atom, int other, int ligand)
    {
        int[] list = this.getBondedAtoms(atom);

        if(list.length == 2)
        {
            return MAX_ATOM_IDX;
        }
        else if(list.length == 3)
        {
            for(int i = 0; i < 3; i++)
            {
                int a = list[i];

                if(a != other && a != ligand)
                    return a;
            }
        }

        return MAX_ATOM_IDX; //unreachable
    }


    public final int getLastChiralLigand(int centre, int[] ligands)
    {
        int[] list = this.getBondedAtoms(centre);

        if(list.length == 3)
        {
            return MAX_ATOM_IDX;
        }
        else if(list.length == 4)
        {
            for(int i = 0; i < 4; i++)
            {
                int ligand = list[i];
                boolean contains = false;

                for(int j = 0; j < 3; j++)
                    if(ligands[j] == ligand)
                        contains = true;

                if(!contains)
                    return ligand;
            }
        }

        return MAX_ATOM_IDX; //unreachable
    }


    public final boolean isExtendedTetrahedralCentre(int centre)
    {
        int[] list = getBondedAtoms(centre);


        if(list.length != 2)
            return false;

        for(int ligand : list)
        {
            int bond = getBond(centre, ligand);

            if(getBondType(bond) != Molecule.BondType.DOUBLE)
                return false;
        }

        return true;
    }


    public final boolean isExtendedCisTrans(int centre)
    {
        for(int i = 0; i < 2; i++)
        {
            int end = getBondAtom(centre, i);
            int[] list = getBondedAtoms(end);

            if(list.length != 2)
                return false;

            for(int ligand : list)
            {
                int bond = getBond(end, ligand);

                if(getBondType(bond) != Molecule.BondType.DOUBLE)
                    return false;
            }
        }

        return true;
    }


    public static byte normalizeAtomStereo(int[] indexes, byte stereo)
    {
        int order = 0;

        for(int i = 0; i < 4; i++)
        {
            int value = 0;

            for(int j = 0; j < 4; j++)
                if(indexes[j] <= indexes[i])
                    value++;

            order = (order << 4) + value;
        }

        boolean reverse = true;

        int[] validReorder = { 0x1234, 0x1423, 0x1342, 0x2314, 0x2431, 0x2143, 0x3124, 0x3412, 0x3241, 0x4213, 0x4321,
                0x4132 };

        for(int i = 0; i < 12; i++)
            if(validReorder[i] == order)
                reverse = false;

        if(reverse)
            return (byte) (~stereo & 0x03);

        return stereo;
    }


    public static byte normalizeBondStereo(int[] indexes, byte conformation)
    {
        boolean reverse = false;

        if(indexes[0] > indexes[1])
            reverse = !reverse;

        if(indexes[2] > indexes[3])
            reverse = !reverse;

        if(reverse)
            return (byte) (~conformation & 0x03);

        return conformation;
    }
}
