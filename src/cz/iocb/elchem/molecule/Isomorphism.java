/*
 * Copyright (C) 2001
 *   Dipartimento di Informatica e Sistemistica,
 *   Universita degli studi di Napoli ``Federico II'
 *   <http://amalfi.dis.unina.it>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 *  1. The above copyright notice and this permission notice shall be included in all copies or substantial
 *      portions of the Software, together with the associated disclaimers.
 *  2. Any modification to the standard distribution of the Software shall be mentioned in a prominent notice
 *      in the documentation provided with the modified distribution, stating clearly how, when and by
 *      whom the Software has been modified.
 *  3. Either the modified distribution shall contain the entire sourcecode of the standard distribution of the
 *      Software, or the documentation shall provide instructions on where the source code of the standard
 *      distribution of the Software can be obtained.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */
package cz.iocb.elchem.molecule;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import cz.iocb.elchem.molecule.Molecule.AtomType;
import cz.iocb.elchem.molecule.Molecule.BondType;



public class Isomorphism
{
    private static final int UNDEFINED_CORE = -1;

    private final SearchMode searchMode;
    private final ChargeMode chargeMode;
    private final IsotopeMode isotopeMode;
    private final StereoMode stereoMode;

    private final Molecule query;
    private final int queryAtomCount;
    private final int[] queryOrder;
    private final int[] queryParents;


    public Isomorphism(Molecule query)
    {
        this(query, SearchMode.SUBSTRUCTURE, ChargeMode.IGNORE, IsotopeMode.IGNORE, StereoMode.IGNORE);
    }


    public Isomorphism(Molecule query, SearchMode graphMode, ChargeMode chargeMode, IsotopeMode isotopeMode,
            StereoMode stereoMode)
    {
        int queryAtomCount = query.getAtomCount();

        this.searchMode = graphMode;
        this.chargeMode = chargeMode;
        this.isotopeMode = isotopeMode;
        this.stereoMode = stereoMode;
        this.query = query;
        this.queryAtomCount = queryAtomCount;
        this.queryOrder = new int[queryAtomCount];
        this.queryParents = new int[queryAtomCount];


        for(int i = 0; i < queryAtomCount; i++)
            this.queryParents[i] = -1;

        byte[] queryFlags = new byte[queryAtomCount];

        for(int idx = 0; idx < queryAtomCount; idx++)
        {
            int selected = -1;
            int fallback = -1;

            for(int i = 0; i < queryAtomCount; i++)
            {
                if(selected == -1 && queryFlags[i] == 1)
                {
                    selected = i;
                    break;
                }

                if(fallback == -1 && queryFlags[i] == 0)
                    fallback = i;
            }

            if(selected == -1)
                selected = fallback;


            int[] queryBondedAtoms = query.getBondedAtoms(selected);

            queryFlags[selected] = 2;

            for(int i = 0; i < queryBondedAtoms.length; i++)
            {
                int idx2 = queryBondedAtoms[i];

                if(queryFlags[idx2] == 0)
                {
                    queryFlags[idx2] = 1;
                    this.queryParents[idx2] = selected;
                }
            }

            this.queryOrder[idx] = selected;
        }
    }


    public final boolean match(Molecule target)
    {
        return new VF2State().match(target, null, 0);
    }


    public List<int[]> matchAll(Molecule target, int limit)
    {
        ArrayList<int[]> result = new ArrayList<int[]>(limit);
        new VF2State().match(target, result, limit);
        return result;
    }


    private final class VF2State
    {
        private int[] queryCore;
        private int[] undosTargetSelector;
        private int[] undosTargetIdx;
        private Molecule target;
        private int targetAtomCount;
        private int[] targetCore;
        private int targetSelector;
        private int targetIdx;
        private int coreLength;
        private int queryIdx;


        private VF2State()
        {
            this.coreLength = 0;
            this.queryCore = new int[queryAtomCount];
            this.undosTargetSelector = new int[queryAtomCount];
            this.undosTargetIdx = new int[queryAtomCount];
        }


        private final boolean nextQuery()
        {
            if(coreLength >= queryAtomCount)
                return false;

            queryIdx = queryOrder[coreLength];
            targetIdx = -1;
            targetSelector = -1;

            return true;
        }


        private final boolean nextTarget()
        {
            int queryParent = queryParents[queryIdx];

            if(queryParent >= 0)
            {
                int targetParent = queryCore[queryParent];
                int[] targetBondedAtoms = target.getBondedAtoms(targetParent);

                for(targetSelector++; targetSelector < targetBondedAtoms.length; targetSelector++)
                {
                    int newTargetIdx = targetBondedAtoms[targetSelector];

                    if(!isCoreDefined(targetCore[newTargetIdx]))
                    {
                        targetIdx = newTargetIdx;
                        return true;
                    }
                }
            }
            else
            {
                for(targetIdx++; targetIdx < targetAtomCount; targetIdx++)
                {
                    if(!isCoreDefined(targetCore[targetIdx]))
                        return true;
                }
            }

            return false;
        }


        private final boolean atomMatches(int queryAtom, int targetAtom)
        {
            byte queryNumber = query.getAtomNumber(queryAtom);
            byte targetNumber = target.getAtomNumber(targetAtom);

            if(queryNumber == AtomType.UNKNOWN || targetNumber == AtomType.UNKNOWN)
                return false;
            else if(queryNumber == targetNumber || queryNumber == AtomType.R)
                return true;
            else if(target.isAtomPseudo(targetAtom))
                return queryNumber == AtomType.Q && (targetNumber == AtomType.M || targetNumber == AtomType.X);
            else if(queryNumber == AtomType.Q)
                return targetNumber != AtomType.C && targetNumber != AtomType.H;
            else if(queryNumber == AtomType.M)
                return target.isAtomMetal(targetAtom);
            else if(queryNumber == AtomType.X)
                return target.isAtomHalogen(targetAtom);
            else
                return false;
        }


        private final boolean bondMatches(int qIdx1, int qIdx2, int tIdx1, int tIdx2)
        {
            int queryBond = query.getBond(qIdx1, qIdx2);
            int targetbond = target.getBond(tIdx1, tIdx2);

            if(queryBond < 0 || targetbond < 0)
                return false;

            byte queryBondType = query.getBondType(queryBond);
            byte targetbondType = target.getBondType(targetbond);


            if(queryBondType == targetbondType || queryBondType == BondType.ANY)
                return true;
            else if(queryBondType == BondType.SINGLE_OR_DOUBLE)
                return targetbondType == BondType.SINGLE || targetbondType == BondType.DOUBLE;
            else if(queryBondType == BondType.SINGLE_OR_AROMATIC)
                return targetbondType == BondType.SINGLE || targetbondType == BondType.AROMATIC;
            else if(queryBondType == BondType.DOUBLE_OR_AROMATIC)
                return targetbondType == BondType.DOUBLE || targetbondType == BondType.AROMATIC;
            else
                return false;
        }


        private final boolean isFeasiblePair()
        {
            if(!atomMatches(queryIdx, targetIdx))
                return false;


            if(chargeMode != ChargeMode.IGNORE)
            {
                byte queryCharge = query.getAtomFormalCharge(queryIdx);
                byte targetCharge = target.getAtomFormalCharge(targetIdx);

                if(queryCharge != targetCharge && (queryCharge != 0 || chargeMode == ChargeMode.DEFAULT_AS_UNCHARGED))
                    return false;
            }


            if(isotopeMode != IsotopeMode.IGNORE)
            {
                byte queryMass = query.getAtomMass(queryIdx);
                byte targetMass = target.getAtomMass(targetIdx);

                if(queryMass != targetMass && (queryMass != 0 || isotopeMode == IsotopeMode.DEFAULT_AS_STANDARD))
                    return false;
            }


            if(!query.hasPseudoAtom() && !target.hasPseudoAtom())
            {
                if(searchMode != SearchMode.EXACT)
                {
                    if(query.getAtomHydrogenCount(queryIdx) > target.getAtomHydrogenCount(targetIdx))
                        return false;
                }
                else
                {
                    if(query.getAtomHydrogenCount(queryIdx) != target.getAtomHydrogenCount(targetIdx))
                        return false;
                }
            }


            int newQuery = 0;
            int newTarget = 0;

            int[] queryBondedAtoms = query.getBondedAtoms(queryIdx);

            for(int other1 : queryBondedAtoms)
            {
                if(isCoreDefined(queryCore[other1]))
                {
                    int other2 = queryCore[other1];

                    if(!bondMatches(queryIdx, other1, targetIdx, other2))
                        return false;
                }
                else
                {
                    newQuery++;
                }
            }


            int[] targetBondedAtoms = target.getBondedAtoms(targetIdx);

            for(int other2 : targetBondedAtoms)
            {
                if(isCoreDefined(targetCore[other2]))
                {
                    if(searchMode == SearchMode.EXACT)
                    {
                        int other1 = targetCore[other2];

                        if(!bondMatches(queryIdx, other1, targetIdx, other2))
                            return false;
                    }
                }
                else
                {
                    newTarget++;
                }
            }

            if(searchMode == SearchMode.EXACT)
                return newQuery == newTarget;
            else
                return newQuery <= newTarget;
        }


        private final void undoAddPair()
        {
            coreLength--;
            int undoTargetSelector = undosTargetSelector[coreLength];
            int undoTargetIdx = undosTargetIdx[coreLength];

            queryIdx = queryOrder[coreLength];

            queryCore[queryIdx] = UNDEFINED_CORE;
            targetCore[undoTargetIdx] = UNDEFINED_CORE;

            targetSelector = undoTargetSelector;
            targetIdx = undoTargetIdx;
        }


        private final void addPair()
        {
            undosTargetSelector[coreLength] = targetSelector;
            undosTargetIdx[coreLength] = targetIdx;

            coreLength++;
            queryCore[queryIdx] = targetIdx;
            targetCore[targetIdx] = queryIdx;
        }


        private final boolean isStereoValid()
        {
            int queryAtomCount = query.getAtomCount();
            int queryBondCount = query.getBondCount();


            for(int queryAtomIdx = 0; queryAtomIdx < queryAtomCount; queryAtomIdx++)
            {
                byte queryStereo = query.getAtomStereo(queryAtomIdx);

                int targetAtomIdx = queryCore[queryAtomIdx];
                byte targetStereo = target.getAtomStereo(targetAtomIdx);

                if(queryStereo != Molecule.TetrahedralStereo.NONE
                        && queryStereo != Molecule.TetrahedralStereo.UNDEFINED)
                {
                    if(targetStereo == Molecule.TetrahedralStereo.NONE
                            || targetStereo == Molecule.TetrahedralStereo.UNDEFINED)
                        continue;


                    if(query.isExtendedTetrahedralCentre(queryAtomIdx))
                    {
                        int[] queryTerminalAtoms = new int[2];
                        int[] queryPreTerminalAtoms = new int[2];
                        int[] queryAtoms = new int[4];
                        int listSize = 0;

                        for(int i = 0; i < 2; i++)
                        {
                            int atom = queryAtomIdx;
                            int bonded = query.getBondedAtoms(queryAtomIdx)[i];

                            while(true)
                            {
                                int[] newList = query.getBondedAtoms(bonded);

                                if(newList.length == 3)
                                {
                                    queryTerminalAtoms[i] = bonded;
                                    queryPreTerminalAtoms[i] = atom;

                                    for(int j = 0; j < 3; j++)
                                    {
                                        int o = newList[j];

                                        if(o == atom)
                                            continue;

                                        queryAtoms[listSize++] = o;
                                    }

                                    break;
                                }
                                else if(newList.length == 2)
                                {
                                    int next = query.getOppositeAtom(bonded, atom);

                                    if(query.getBondType(query.getBond(bonded, next)) != BondType.DOUBLE)
                                    {
                                        queryTerminalAtoms[i] = bonded;
                                        queryPreTerminalAtoms[i] = atom;
                                        queryAtoms[listSize++] = next;
                                        queryAtoms[listSize++] = Molecule.MAX_ATOM_IDX;
                                        break;
                                    }

                                    atom = bonded;
                                    bonded = next;
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }

                        if(listSize < 4)
                            continue;

                        sortBondAtoms(queryAtoms);


                        int targetTerminalAtom0 = queryCore[queryTerminalAtoms[0]];
                        int targetTerminalAtom1 = queryCore[queryTerminalAtoms[1]];
                        int targetPreTerminalAtom0 = queryCore[queryPreTerminalAtoms[0]];
                        int targetPreTerminalAtom1 = queryCore[queryPreTerminalAtoms[1]];

                        int[] targetAtoms = new int[] { -1, -1, -1, -1 };

                        for(int i = 0; i < 4; i++)
                            if(queryAtoms[i] != Molecule.MAX_ATOM_IDX)
                                targetAtoms[i] = queryCore[queryAtoms[i]];

                        if(queryAtoms[1] == Molecule.MAX_ATOM_IDX)
                            targetAtoms[1] = target.getLastStereoBondLigand(targetTerminalAtom0, targetPreTerminalAtom0,
                                    targetAtoms[0]);

                        if(queryAtoms[3] == Molecule.MAX_ATOM_IDX)
                            targetAtoms[3] = target.getLastStereoBondLigand(targetTerminalAtom1, targetPreTerminalAtom1,
                                    targetAtoms[2]);

                        if(Molecule.normalizeBondStereo(targetAtoms, targetStereo) != queryStereo)
                            return false;
                    }
                    else
                    {
                        int[] bondedAtoms = query.getBondedAtoms(queryAtomIdx);

                        if(bondedAtoms.length < 3)
                            continue;

                        int[] queryAtoms = new int[4];

                        for(int i = 0; i < bondedAtoms.length; i++)
                            queryAtoms[i] = bondedAtoms[i];

                        if(bondedAtoms.length == 3)
                            queryAtoms[3] = Molecule.MAX_ATOM_IDX;

                        sortStereoAtoms(queryAtoms);


                        int[] targetAtoms = new int[] { -1, -1, -1, -1 };

                        for(int i = 0; i < bondedAtoms.length; i++)
                            targetAtoms[i] = queryCore[queryAtoms[i]];

                        if(bondedAtoms.length == 3)
                            targetAtoms[3] = target.getLastChiralLigand(queryCore[queryAtomIdx], targetAtoms);

                        if(Molecule.normalizeAtomStereo(targetAtoms, targetStereo) != queryStereo)
                            return false;
                    }
                }
            }


            for(int queryBondIdx = 0; queryBondIdx < queryBondCount; queryBondIdx++)
            {
                byte queryStereo = query.getBondStereo(queryBondIdx);

                int queryBondAtom0 = query.getBondAtom(queryBondIdx, 0);
                int queryBondAtom1 = query.getBondAtom(queryBondIdx, 1);

                int targetBondAtom0 = queryCore[queryBondAtom0];
                int targetBondAtom1 = queryCore[queryBondAtom1];
                int targetBondIdx = target.getBond(targetBondAtom0, targetBondAtom1);
                byte targetStereo = target.getBondStereo(targetBondIdx);

                if(queryStereo != Molecule.BondStereo.NONE && queryStereo != Molecule.BondStereo.UNDEFINED)
                {
                    if(targetStereo == Molecule.BondStereo.NONE || targetStereo == Molecule.BondStereo.UNDEFINED)
                        continue;

                    if(query.isExtendedCisTrans(queryBondIdx))
                    {
                        int[] queryTerminalAtoms = new int[2];
                        int[] queryPreTerminalAtoms = new int[2];
                        int[] queryAtoms = new int[4];
                        int listSize = 0;

                        for(int i = 0; i < 2; i++)
                        {
                            int atom = query.getBondAtom(queryBondIdx, 0);
                            int bonded = query.getOtherBondAtom(queryBondIdx, atom);

                            while(true)
                            {
                                int[] newList = query.getBondedAtoms(bonded);

                                if(newList.length == 3)
                                {
                                    queryTerminalAtoms[i] = bonded;
                                    queryPreTerminalAtoms[i] = atom;

                                    for(int j = 0; j < 3; j++)
                                        if(newList[j] != atom)
                                            queryAtoms[listSize++] = newList[j];

                                    break;
                                }
                                else if(newList.length == 2)
                                {
                                    int next = query.getOppositeAtom(bonded, atom);

                                    if(query.getBondType(query.getBond(bonded, next)) != BondType.DOUBLE)
                                    {
                                        queryTerminalAtoms[i] = bonded;
                                        queryPreTerminalAtoms[i] = atom;
                                        queryAtoms[listSize++] = next;
                                        queryAtoms[listSize++] = Molecule.MAX_ATOM_IDX;
                                        break;
                                    }

                                    atom = bonded;
                                    bonded = next;
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }

                        if(listSize < 4)
                            continue;

                        sortBondAtoms(queryAtoms);


                        int targetTerminalAtom0 = queryCore[queryTerminalAtoms[0]];
                        int targetTerminalAtom1 = queryCore[queryTerminalAtoms[1]];
                        int targetPreTerminalAtom0 = queryCore[queryPreTerminalAtoms[0]];
                        int targetPreTerminalAtom1 = queryCore[queryPreTerminalAtoms[1]];

                        int[] targetAtoms = new int[] { -1, -1, -1, -1 };

                        for(int i = 0; i < 4; i++)
                            if(queryAtoms[i] != Molecule.MAX_ATOM_IDX)
                                targetAtoms[i] = queryCore[queryAtoms[i]];

                        if(queryAtoms[1] == Molecule.MAX_ATOM_IDX)
                            targetAtoms[1] = target.getLastStereoBondLigand(targetTerminalAtom0, targetPreTerminalAtom0,
                                    targetAtoms[0]);

                        if(queryAtoms[3] == Molecule.MAX_ATOM_IDX)
                            targetAtoms[3] = target.getLastStereoBondLigand(targetTerminalAtom1, targetPreTerminalAtom1,
                                    targetAtoms[2]);

                        if(Molecule.normalizeBondStereo(targetAtoms, targetStereo) != queryStereo)
                            return false;
                    }
                    else
                    {
                        int[] bondedAtoms0 = query.getBondedAtoms(queryBondAtom0);
                        int[] bondedAtoms1 = query.getBondedAtoms(queryBondAtom1);

                        if(bondedAtoms0.length < 2 || bondedAtoms1.length < 2)
                            continue;


                        int[] queryAtoms = new int[4];

                        int idx = 0;

                        for(int i = 0; i < bondedAtoms0.length; i++)
                            if(bondedAtoms0[i] != queryBondAtom1)
                                queryAtoms[idx++] = bondedAtoms0[i];

                        if(bondedAtoms0.length == 2)
                            queryAtoms[idx++] = Molecule.MAX_ATOM_IDX;


                        for(int i = 0; i < bondedAtoms1.length; i++)
                            if(bondedAtoms1[i] != queryBondAtom0)
                                queryAtoms[idx++] = bondedAtoms1[i];

                        if(bondedAtoms1.length == 2)
                            queryAtoms[idx] = Molecule.MAX_ATOM_IDX;

                        sortBondAtoms(queryAtoms);


                        int[] targetAtoms = new int[] { -1, -1, -1, -1 };

                        for(int i = 0; i < 4; i++)
                            if(queryAtoms[i] != Molecule.MAX_ATOM_IDX)
                                targetAtoms[i] = queryCore[queryAtoms[i]];

                        if(queryAtoms[1] == Molecule.MAX_ATOM_IDX)
                            targetAtoms[1] = target.getLastStereoBondLigand(targetBondAtom0, targetBondAtom1,
                                    targetAtoms[0]);

                        if(queryAtoms[3] == Molecule.MAX_ATOM_IDX)
                            targetAtoms[3] = target.getLastStereoBondLigand(targetBondAtom1, targetBondAtom0,
                                    targetAtoms[2]);

                        if(Molecule.normalizeBondStereo(targetAtoms, targetStereo) != queryStereo)
                            return false;
                    }
                }
            }


            return true;
        }


        private final boolean isMatchValid()
        {
            int queryAtomCount = query.getAtomCount();

            int queryBondCount = query.getBondCount();
            int targetBondCount = target.getBondCount();

            /*
            The goal has been reached, so the query has been mapped to target,
            therfore query is a substructure of the target.

            However, if this was an R-group query, the result could still be
            rejected.If the RestH property is true for some atom with an R-group
            linked, then the R-group may only be substituted with a member
            of the Rgroup or with H..
            This can be verified:
            - find any atom in the query with RestH flagged
            - find the atom mapped to it in the target container
            - see if the target has more (non hydrogen) bonds than the query.
              if so,discard it.
            */
            if(query.hasRestHydrogenFlags())
            {
                for(int i = 0; i < queryAtomCount; i++)
                {
                    if(query.getAtomRestHydrogenFlag(i) == true)
                    {
                        int targetAtomIdx = queryCore[i];

                        int queryConnectivityCount = 0;
                        int targetConnectivityCount = 0;

                        for(int b = 0; b < queryBondCount; b++)
                            if(query.isAtomInBond(i, b)
                                    && !(query.getAtomNumber(query.getOtherBondAtom(b, i)) == AtomType.H))
                                queryConnectivityCount++;

                        for(int b = 0; b < targetBondCount; b++)
                            if(target.isAtomInBond(targetAtomIdx, b)
                                    && !(target.getAtomNumber(target.getOtherBondAtom(b, targetAtomIdx)) == AtomType.H))
                                targetConnectivityCount++;

                        if(targetConnectivityCount > queryConnectivityCount)
                            return false;
                    }
                }
            }


            if(stereoMode == StereoMode.STRICT && !isStereoValid())
                return false;


            return true;
        }


        private final boolean matchCore(List<int[]> result, int limit)
        {
            recursion:
            while(true)
            {
                if(coreLength == query.getAtomCount() && isMatchValid())
                {
                    if(result == null)
                        return true;

                    int[] match = queryCore.clone();

                    Arrays.sort(match);

                    boolean included = false;

                    for(int[] r : result)
                    {
                        boolean same = true;

                        for(int i = 0; i < match.length; i++)
                        {
                            if(match[i] != r[i])
                            {
                                same = false;
                                break;
                            }
                        }

                        if(same)
                        {
                            included = true;
                            break;
                        }
                    }

                    if(!included)
                        result.add(match);

                    if(limit == result.size())
                        return false;
                }

                if(coreLength == query.getAtomCount() || !nextQuery())
                {
                    if(coreLength == 0)
                        return false;

                    undoAddPair();
                }

                while(true)
                {
                    while(nextTarget())
                    {
                        if(Thread.currentThread().isInterrupted())
                            return false;

                        if(isFeasiblePair())
                        {
                            addPair();
                            continue recursion;
                        }
                    }

                    if(coreLength == 0)
                        return false;

                    undoAddPair();
                }
            }
        }


        private final boolean match(Molecule target, List<int[]> result, int limit)
        {
            if(searchMode != SearchMode.EXACT)
            {
                if(query.getOriginalAtomCount() > target.getOriginalAtomCount())
                    return false;

                if(query.getOriginalBondCount() > target.getOriginalBondCount())
                    return false;

                if(queryAtomCount > target.getAtomCount() || query.getBondCount() > target.getBondCount())
                    return false;
            }
            else
            {
                if(query.getOriginalAtomCount() != target.getOriginalAtomCount())
                    return false;

                if(query.getOriginalBondCount() != target.getOriginalBondCount())
                    return false;

                if(queryAtomCount != target.getAtomCount() || query.getBondCount() != target.getBondCount())
                    return false;
            }

            int targetAtomCount = target.getAtomCount();

            this.target = target;
            this.targetAtomCount = targetAtomCount;
            coreLength = 0;

            targetCore = new int[targetAtomCount];

            for(int i = 0; i < targetAtomCount; i++)
                targetCore[i] = UNDEFINED_CORE;

            for(int i = 0; i < queryAtomCount; i++)
                queryCore[i] = UNDEFINED_CORE;


            return matchCore(result, limit);
        }
    }


    private static final boolean isCoreDefined(int value)
    {
        return value >= 0;
    }


    private static final void swapIndexes(int[] array, int a, int b)
    {
        int tmp = array[a];
        array[a] = array[b];
        array[b] = tmp;
    }


    private static final void sortStereoAtoms(int[] array)
    {
        for(int n = 4; n > 0; n--)
            for(int i = 1; i < n; i++)
                if(array[i - 1] > array[i])
                    swapIndexes(array, i, i - 1);
    }


    private static final void sortBondAtoms(int[] array)
    {
        if(array[0] > array[1])
            swapIndexes(array, 0, 1);

        if(array[2] > array[3])
            swapIndexes(array, 2, 3);
    }
}
