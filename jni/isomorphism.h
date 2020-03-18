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

#ifndef ISOMORPHISM_H__
#define ISOMORPHISM_H__

#include <stdbool.h>
#include "memory.h"
#include "compiler.h"
#include "molecule.h"


#define UNDEFINED_CORE              -1
#define is_core_defined(value)      ((value) >= 0)


typedef enum
{
    GRAPH_SUBSTRUCTURE = 0,
    GRAPH_EXACT = 1
}
GraphMode;


typedef enum
{
    CHARGE_IGNORE = 0,
    CHARGE_DEFAULT_AS_UNCHARGED = 1,
    CHARGE_DEFAULT_AS_ANY = 2
}
ChargeMode;


typedef enum
{
    ISOTOPE_IGNORE = 0,
    ISOTOPE_DEFAULT_AS_STANDARD = 1,
    ISOTOPE_DEFAULT_AS_ANY = 2
}
IsotopeMode;


typedef enum
{
    STEREO_IGNORE = 0,
    STEREO_STRICT = 1
}
StereoMode;


typedef struct
{
    AtomIdx targetSelector;
    AtomIdx targetIdx;
}
VF2Undo;


typedef struct
{
    uint counter;

    GraphMode graphMode;
    ChargeMode chargeMode;
    IsotopeMode isotopeMode;
    StereoMode stereoMode;

    const Molecule *restrict query;
    const Molecule *restrict target;

    int queryAtomCount;
    int targetAtomCount;

    AtomIdx *restrict queryCore;
    AtomIdx *restrict targetCore;
    int coreLength;

    AtomIdx *restrict queryOrder;
    AtomIdx *restrict queryParents;
    AtomIdx queryIdx;

    AtomIdx targetSelector;
    AtomIdx targetIdx;

    VF2Undo *restrict undos;
}
VF2State;


static inline void swap_idx(AtomIdx *restrict a, AtomIdx *restrict b)
{
    AtomIdx t = *a;
    *a = *b;
    *b = t;
}


static inline void sort_stereo_atoms(AtomIdx array[4])
{
    for(int n = 4; n > 0; n--)
        for(int i = 1; i < n; i++)
            if(array[i-1] > array[i])
                swap_idx(array + i, array + i - 1);
}


static inline void sort_bond_atoms(AtomIdx array[4])
{
    if(array[0] > array[1])
        swap_idx(array + 0, array + 1);

    if(array[2] > array[3])
        swap_idx(array + 2, array + 3);
}


static inline bool is_extended_tetrahedral_centre(const Molecule *restrict molecule, AtomIdx centre)
{
    MolSize listSize = molecule_get_bonded_atom_list_size(molecule, centre);

    if(listSize != 2)
        return false;

    for(int i = 0; i < listSize; i++)
    {
        AtomIdx ligand = molecule_get_bonded_atom_list(molecule, centre)[i];
        BondIdx bond = molecule_get_bond(molecule, centre, ligand);

        if(molecule_get_bond_type(molecule, bond) != BOND_DOUBLE)
            return false;
    }

    return true;
}


static inline bool is_extended_cis_trans(const Molecule *restrict molecule, BondIdx centre)
{
    for(int i = 0; i < 2; i++)
    {
        AtomIdx end = molecule_bond_atoms(molecule, centre)[i];
        MolSize listSize = molecule_get_bonded_atom_list_size(molecule, end);

        if(listSize != 2)
            return false;

        for(int i = 0; i < listSize; i++)
        {
            AtomIdx ligand = molecule_get_bonded_atom_list(molecule, end)[i];
            BondIdx bond = molecule_get_bond(molecule, end, ligand);

            if(molecule_get_bond_type(molecule, bond) != BOND_DOUBLE)
                return false;
        }
    }

    return true;
}


static inline size_t vf2state_mem_size(const uint8_t *restrict data, bool extended)
{
    int atomCount = (*(data + 0) << 8 | *(data + 1)) + (*(data + 2) << 8 | *(data + 3));

    if(extended)
        atomCount += *(data + 4) << 8 | *(data + 5);

    return align_size(sizeof(VF2State)) + 3 * align_size(atomCount * sizeof(AtomIdx)) + align_size(atomCount * sizeof(VF2Undo));
}


static inline size_t vf2state_extended_mem_size(const Molecule *restrict molecule)
{
    int atomCount = molecule->heavyAtomCount + molecule->hydrogenAtomCount;
    return align_size(sizeof(VF2State)) + 3 * align_size(atomCount * sizeof(AtomIdx)) + align_size(atomCount * sizeof(VF2Undo));
}


static inline size_t vf2state_match_mem_size(const uint8_t *restrict data, bool extended)
{
    int atomCount = (*(data + 0) << 8 | *(data + 1)) + (*(data + 2) << 8 | *(data + 3));

    if(extended)
        atomCount += *(data + 4) << 8 | *(data + 5);

    return align_size(atomCount * sizeof(AtomIdx));
}


static inline VF2State *vf2state_create(void *memory, const Molecule *restrict query, GraphMode graphMode,
        ChargeMode chargeMode, IsotopeMode isotopeMode, StereoMode stereoMode)
{
    VF2State *restrict vf2state = (VF2State *) alloc_memory(&memory, sizeof(VF2State));

    int queryAtomCount = query->atomCount;

    AtomIdx *restrict queryCore = (AtomIdx *) alloc_memory(&memory, queryAtomCount * sizeof(AtomIdx));
    AtomIdx *restrict queryOrder = (AtomIdx *) alloc_memory(&memory, queryAtomCount * sizeof(AtomIdx));
    AtomIdx *restrict queryParents = (AtomIdx *) alloc_memory(&memory, queryAtomCount * sizeof(AtomIdx));
    VF2Undo *restrict undos = (VF2Undo *) alloc_memory(&memory, queryAtomCount * sizeof(VF2Undo));


    for(int i = 0; i < queryAtomCount; i++)
        queryParents[i] = -1;

    uint8_t queryFlags[queryAtomCount];

    for(int i = 0; i < queryAtomCount; i++)
        queryFlags[i] = 0;

    for(int idx = 0; idx < queryAtomCount; idx++)
    {
        AtomIdx selected = -1;
        AtomIdx fallback = -1;

        for(AtomIdx i = 0; i < queryAtomCount; i++)
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


        AtomIdx *restrict queryBondedAtomList = molecule_get_bonded_atom_list(query, selected);
        MolSize queryBondedAtomListSize = molecule_get_bonded_atom_list_size(query, selected);

        queryFlags[selected] = 2;

        for(int i = 0; i < queryBondedAtomListSize; i++)
        {
            AtomIdx idx = queryBondedAtomList[i];

            if(queryFlags[idx] == 0)
            {
                queryFlags[idx] = 1;
                queryParents[idx] = selected;
            }
        }

        queryOrder[idx] = selected;
    }


    vf2state->graphMode = graphMode;
    vf2state->chargeMode = chargeMode;
    vf2state->isotopeMode = isotopeMode;
    vf2state->stereoMode = stereoMode;
    vf2state->query = query;
    vf2state->queryAtomCount = queryAtomCount;
    vf2state->coreLength = 0;
    vf2state->queryCore = queryCore;
    vf2state->queryOrder = queryOrder;
    vf2state->queryParents = queryParents;
    vf2state->undos = undos;

    return vf2state;
}


static inline bool vf2state_next_query(VF2State *restrict vf2state)
{
    if(unlikely(vf2state->coreLength >= vf2state->queryAtomCount))
        return false;

    vf2state->queryIdx = vf2state->queryOrder[vf2state->coreLength];
    vf2state->targetIdx = -1;
    vf2state->targetSelector = -1;

    return true;
}


static inline bool vf2state_next_target(VF2State *restrict vf2state)
{
    AtomIdx query_parent = vf2state->queryParents[vf2state->queryIdx];

    if(likely(query_parent >= 0))
    {
        AtomIdx target_parent = vf2state->queryCore[query_parent];
        AtomIdx *restrict targetBondedAtomList = molecule_get_bonded_atom_list(vf2state->target, target_parent);
        MolSize targetBondedAtomListSize = molecule_get_bonded_atom_list_size(vf2state->target, target_parent);

        for(vf2state->targetSelector++; vf2state->targetSelector < targetBondedAtomListSize; vf2state->targetSelector++)
        {
            AtomIdx targetIdx = targetBondedAtomList[vf2state->targetSelector];

            if(!is_core_defined(vf2state->targetCore[targetIdx]))
            {
                vf2state->targetIdx = targetIdx;
                return true;
            }
        }
    }
    else
    {
        for(vf2state->targetIdx++; vf2state->targetIdx < vf2state->targetAtomCount; vf2state->targetIdx++)
        {
            if(!is_core_defined(vf2state->targetCore[vf2state->targetIdx]))
                return true;
        }
    }

    return false;
}


static inline bool vf2state_atom_matches(const VF2State *restrict vf2state, AtomIdx queryAtom, AtomIdx targetAtom)
{
    int8_t queryAtomNumber = molecule_get_atom_number(vf2state->query, queryAtom);
    int8_t targetAtomNumber = molecule_get_atom_number(vf2state->target, targetAtom);

    if(vf2state->graphMode == GRAPH_EXACT)
        return queryAtomNumber == targetAtomNumber;
    else if(queryAtomNumber == UNKNOWN_ATOM_NUMBER || targetAtomNumber == UNKNOWN_ATOM_NUMBER)
        return false;
    else if(queryAtomNumber == targetAtomNumber || queryAtomNumber == R_ATOM_NUMBER)
        return true;
    else if(molecule_is_pseudo_atom(vf2state->target, targetAtom))
        return queryAtomNumber == Q_ATOM_NUMBER && (targetAtomNumber == M_ATOM_NUMBER || targetAtomNumber == X_ATOM_NUMBER);
    else if(queryAtomNumber == Q_ATOM_NUMBER)
        return targetAtomNumber != C_ATOM_NUMBER && targetAtomNumber != H_ATOM_NUMBER;
    else if(queryAtomNumber == M_ATOM_NUMBER)
        return molecule_is_metal(vf2state->target, targetAtom);
    else if(queryAtomNumber == X_ATOM_NUMBER)
        return molecule_is_halogen(vf2state->target, targetAtom);
    else
        return false;
}


static inline bool vf2state_bond_matches(const VF2State *restrict vf2state, AtomIdx qIdx1, AtomIdx qIdx2, AtomIdx tIdx1, AtomIdx tIdx2)
{
    BondIdx queryBond = molecule_get_bond(vf2state->query, qIdx1, qIdx2);
    BondIdx targetbond = molecule_get_bond(vf2state->target, tIdx1, tIdx2);

    if(likely(queryBond < 0 || targetbond < 0))
        return false;

    uint8_t queryBondType = molecule_get_bond_type(vf2state->query, queryBond);
    uint8_t targetbondType = molecule_get_bond_type(vf2state->target, targetbond);

    if(vf2state->graphMode == GRAPH_EXACT)
        return queryBondType == targetbondType;
    else if(queryBondType == targetbondType || queryBondType == BOND_ANY)
        return true;
    else if(queryBondType == BOND_SINGLE_OR_DOUBLE)
        return targetbondType == BOND_SINGLE || targetbondType == BOND_DOUBLE;
    else if(queryBondType == BOND_SINGLE_OR_AROMATIC)
        return targetbondType == BOND_SINGLE || targetbondType == BOND_AROMATIC;
    else if(queryBondType == BOND_DOUBLE_OR_AROMATIC)
        return targetbondType == BOND_DOUBLE || targetbondType == BOND_AROMATIC;
    else
        return false;
}


static inline bool vf2state_is_feasible_pair(const VF2State *restrict vf2state)
{
    if(likely(!vf2state_atom_matches(vf2state, vf2state->queryIdx, vf2state->targetIdx)))
        return false;


    if(vf2state->chargeMode != CHARGE_IGNORE)
    {
        int8_t queryCharge = molecule_get_formal_charge(vf2state->query, vf2state->queryIdx);
        int8_t targetCharge = molecule_get_formal_charge(vf2state->target, vf2state->targetIdx);

        if(queryCharge != targetCharge && (queryCharge != 0 || vf2state->chargeMode == CHARGE_DEFAULT_AS_UNCHARGED))
            return false;
    }


    if(vf2state->isotopeMode != ISOTOPE_IGNORE)
    {
        int8_t queryMass = molecule_get_atom_mass(vf2state->query, vf2state->queryIdx);
        int8_t targetMass = molecule_get_atom_mass(vf2state->target, vf2state->targetIdx);

        if(queryMass != targetMass && (queryMass != 0 || vf2state->isotopeMode == ISOTOPE_DEFAULT_AS_STANDARD))
            return false;
    }


    if(likely(vf2state->graphMode == GRAPH_EXACT))
    {
        if(unlikely(molecule_get_hydrogen_count(vf2state->query, vf2state->queryIdx) !=
                molecule_get_hydrogen_count(vf2state->target, vf2state->targetIdx)))
            return false;
    }
    else if(!vf2state->query->hasPseudoAtom && !vf2state->target->hasPseudoAtom)
    {
        if(unlikely(molecule_get_hydrogen_count(vf2state->query, vf2state->queryIdx) >
                molecule_get_hydrogen_count(vf2state->target, vf2state->targetIdx)))
            return false;
    }


    int newQuery = 0;
    int newTarget = 0;

    AtomIdx *restrict queryBondedAtomList = molecule_get_bonded_atom_list(vf2state->query, vf2state->queryIdx);
    MolSize queryBondedAtomListSize = molecule_get_bonded_atom_list_size(vf2state->query, vf2state->queryIdx);

    for(int i = 0; i < queryBondedAtomListSize; i++)
    {
        AtomIdx other1 = queryBondedAtomList[i];

        if(is_core_defined(vf2state->queryCore[other1]))
        {
            AtomIdx other2 = vf2state->queryCore[other1];

            if(!vf2state_bond_matches(vf2state, vf2state->queryIdx, other1, vf2state->targetIdx, other2))
                return false;
        }
        else
        {
            newQuery++;
        }
    }


    AtomIdx *restrict targetBondedAtomList = molecule_get_bonded_atom_list(vf2state->target, vf2state->targetIdx);
    MolSize targetBondedAtomListSize = molecule_get_bonded_atom_list_size(vf2state->target, vf2state->targetIdx);

    for(int i = 0; i < targetBondedAtomListSize; i++)
    {
        AtomIdx other2 = targetBondedAtomList[i];

        if(is_core_defined(vf2state->targetCore[other2]))
        {
            if(unlikely(vf2state->graphMode == GRAPH_EXACT))
            {
                AtomIdx other1 = vf2state->targetCore[other2];

                if(!vf2state_bond_matches(vf2state, vf2state->queryIdx, other1, vf2state->targetIdx, other2))
                    return false;
            }
        }
        else
        {
            newTarget++;
        }
    }

    if(unlikely(vf2state->graphMode == GRAPH_EXACT))
        return newQuery == newTarget;
    else
        return newQuery <= newTarget;
}


static inline void vf2state_undo_add_pair(VF2State *restrict vf2state)
{
    VF2Undo *restrict undo = &vf2state->undos[--vf2state->coreLength];

    vf2state->queryIdx = vf2state->queryOrder[vf2state->coreLength];

    vf2state->queryCore[vf2state->queryIdx] = UNDEFINED_CORE;
    vf2state->targetCore[undo->targetIdx] = UNDEFINED_CORE;

    vf2state->targetSelector = undo->targetSelector;
    vf2state->targetIdx = undo->targetIdx;
}


static inline void vf2state_add_pair(VF2State *restrict vf2state)
{
    VF2Undo *restrict undo = &vf2state->undos[vf2state->coreLength];

    vf2state->coreLength++;
    vf2state->queryCore[vf2state->queryIdx] = vf2state->targetIdx;
    vf2state->targetCore[vf2state->targetIdx] = vf2state->queryIdx;

    undo->targetSelector = vf2state->targetSelector;
    undo->targetIdx = vf2state->targetIdx;
}


static inline bool vf2state_is_stereo_valid(const VF2State *restrict vf2state)
{
    const Molecule *restrict query = vf2state->query;
    const Molecule *restrict target = vf2state->target;

    int queryAtomCount = query->atomCount;
    int queryBondCount = query->bondCount;


    for(AtomIdx queryAtomIdx = 0; queryAtomIdx < queryAtomCount; queryAtomIdx++)
    {
        uint8_t queryStereo = molecule_get_atom_stereo(query, queryAtomIdx);

        AtomIdx targetAtomIdx = vf2state->queryCore[queryAtomIdx];
        uint8_t targetStereo = molecule_get_atom_stereo(target, targetAtomIdx);

        if(queryStereo == TETRAHEDRAL_STEREO_UNDEFINED)
        {
            if(vf2state->graphMode == GRAPH_EXACT && targetStereo != TETRAHEDRAL_STEREO_UNDEFINED)
                return false;
        }
        else if(queryStereo != TETRAHEDRAL_STEREO_NONE)
        {
            if(vf2state->graphMode == GRAPH_EXACT && targetStereo == TETRAHEDRAL_STEREO_UNDEFINED)
                return false;

            if(targetStereo == TETRAHEDRAL_STEREO_NONE || targetStereo == TETRAHEDRAL_STEREO_UNDEFINED)
                continue;


            if(is_extended_tetrahedral_centre(query, queryAtomIdx))
            {
                AtomIdx queryTerminalAtoms[2];
                AtomIdx queryPreTerminalAtoms[2];
                AtomIdx queryAtoms[4];
                MolSize listSize = 0;

                for(int i = 0; i < 2; i++)
                {
                    AtomIdx atom = queryAtomIdx;
                    AtomIdx bonded = molecule_get_bonded_atom_list(query, queryAtomIdx)[i];

                    while(true)
                    {
                        MolSize newListSize = molecule_get_bonded_atom_list_size(query, bonded);

                        if(newListSize == 3)
                        {
                            queryTerminalAtoms[i] = bonded;
                            queryPreTerminalAtoms[i] = atom;

                            for(int j = 0; j < 3; j++)
                            {
                                AtomIdx o = molecule_get_bonded_atom_list(query, bonded)[j];

                                if(o == atom)
                                    continue;

                                queryAtoms[listSize++] = o;
                            }

                            break;
                        }
                        else if(newListSize == 2)
                        {
                            AtomIdx next = molecule_get_opposite_atom(query, bonded, atom);

                            if(molecule_get_bond_type(query, molecule_get_bond(query, bonded, next)) != BOND_DOUBLE)
                            {
                                queryTerminalAtoms[i] = bonded;
                                queryPreTerminalAtoms[i] = atom;
                                queryAtoms[listSize++] = next;
                                queryAtoms[listSize++] = MAX_ATOM_IDX;
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

                sort_bond_atoms(queryAtoms);


                AtomIdx targetTerminalAtom0 = vf2state->queryCore[queryTerminalAtoms[0]];
                AtomIdx targetTerminalAtom1 = vf2state->queryCore[queryTerminalAtoms[1]];
                AtomIdx targetPreTerminalAtom0 = vf2state->queryCore[queryPreTerminalAtoms[0]];
                AtomIdx targetPreTerminalAtom1 = vf2state->queryCore[queryPreTerminalAtoms[1]];

                AtomIdx targetAtoms[4] = { -1, -1, -1, -1 };

                for(int i = 0; i < 4; i++)
                    if(queryAtoms[i] != MAX_ATOM_IDX)
                        targetAtoms[i] = vf2state->queryCore[queryAtoms[i]];

                if(queryAtoms[1] == MAX_ATOM_IDX)
                    targetAtoms[1] = molecule_get_last_stereo_bond_ligand(target, targetTerminalAtom0, targetPreTerminalAtom0, targetAtoms[0]);

                if(queryAtoms[3] == MAX_ATOM_IDX)
                    targetAtoms[3] = molecule_get_last_stereo_bond_ligand(target, targetTerminalAtom1, targetPreTerminalAtom1, targetAtoms[2]);

                if(normalize_bond_stereo(targetAtoms, targetStereo) != queryStereo)
                    return false;
            }
            else
            {
                AtomIdx *bondedAtomList = molecule_get_bonded_atom_list(query, queryAtomIdx);
                MolSize listSize = molecule_get_bonded_atom_list_size(query, queryAtomIdx);

                if(listSize < 3)
                    continue;

                AtomIdx queryAtoms[4];

                for(int i = 0; i < listSize; i++)
                    queryAtoms[i] = bondedAtomList[i];

                if(listSize == 3)
                    queryAtoms[3] = MAX_ATOM_IDX;

                sort_stereo_atoms(queryAtoms);


                AtomIdx targetAtoms[4] = { -1, -1, -1, -1 };

                for(int i = 0; i < listSize; i++)
                    targetAtoms[i] = vf2state->queryCore[queryAtoms[i]];

                if(listSize == 3)
                    targetAtoms[3] = molecule_get_last_chiral_ligand(target, vf2state->queryCore[queryAtomIdx], targetAtoms);

                if(normalize_atom_stereo(targetAtoms, targetStereo) != queryStereo)
                    return false;
            }
        }
    }


    for(BondIdx queryBondIdx = 0; queryBondIdx < queryBondCount; queryBondIdx++)
    {
        uint8_t queryStereo = molecule_get_bond_stereo(query, queryBondIdx);

        AtomIdx *queryBondAtoms = molecule_bond_atoms(query, queryBondIdx);
        AtomIdx targetBondAtom0 = vf2state->queryCore[queryBondAtoms[0]];
        AtomIdx targetBondAtom1 = vf2state->queryCore[queryBondAtoms[1]];
        BondIdx targetBondIdx = molecule_get_bond(target, targetBondAtom0, targetBondAtom1);
        uint8_t targetStereo = molecule_get_bond_stereo(target, targetBondIdx);

        if(queryStereo == BOND_STEREO_UNDEFINED)
        {
            if(vf2state->graphMode == GRAPH_EXACT && targetStereo != BOND_STEREO_UNDEFINED)
                return false;
        }
        else if(queryStereo != BOND_STEREO_NONE)
        {
            if(vf2state->graphMode == GRAPH_EXACT && targetStereo == BOND_STEREO_UNDEFINED)
                return false;

            if(targetStereo == BOND_STEREO_NONE || targetStereo == BOND_STEREO_UNDEFINED)
                continue;


            if(is_extended_cis_trans(query, queryBondIdx))
            {
                AtomIdx queryTerminalAtoms[2];
                AtomIdx queryPreTerminalAtoms[2];
                AtomIdx queryAtoms[4];
                int listSize = 0;

                for(int i = 0; i < 2; i++)
                {
                    AtomIdx atom = molecule_bond_atoms(query, queryBondIdx)[i];
                    AtomIdx bonded = molecule_get_other_bond_atom(query, queryBondIdx, atom);


                    while(true)
                    {
                        AtomIdx *newList = molecule_get_bonded_atom_list(query, bonded);
                        MolSize newListSize = molecule_get_bonded_atom_list_size(query, bonded);

                        if(newListSize == 3)
                        {
                            queryTerminalAtoms[i] = bonded;
                            queryPreTerminalAtoms[i] = atom;

                            for(int j = 0; j < 3; j++)
                                if(newList[j] != atom)
                                    queryAtoms[listSize++] = newList[j];

                            break;
                        }
                        else if(newListSize == 2)
                        {
                            AtomIdx next = molecule_get_opposite_atom(query, bonded, atom);

                            if(molecule_get_bond_type(query, molecule_get_bond(query, bonded, next)) != BOND_DOUBLE)
                            {
                                queryTerminalAtoms[i] = bonded;
                                queryPreTerminalAtoms[i] = atom;
                                queryAtoms[listSize++] = next;
                                queryAtoms[listSize++] = MAX_ATOM_IDX;
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

                sort_bond_atoms(queryAtoms);


                AtomIdx targetTerminalAtom0 = vf2state->queryCore[queryTerminalAtoms[0]];
                AtomIdx targetTerminalAtom1 = vf2state->queryCore[queryTerminalAtoms[1]];
                AtomIdx targetPreTerminalAtom0 = vf2state->queryCore[queryPreTerminalAtoms[0]];
                AtomIdx targetPreTerminalAtom1 = vf2state->queryCore[queryPreTerminalAtoms[1]];

                AtomIdx targetAtoms[4] = { -1, -1, -1, -1 };

                for(int i = 0; i < 4; i++)
                    if(queryAtoms[i] != MAX_ATOM_IDX)
                        targetAtoms[i] = vf2state->queryCore[queryAtoms[i]];

                if(queryAtoms[1] == MAX_ATOM_IDX)
                    targetAtoms[1] = molecule_get_last_stereo_bond_ligand(target, targetTerminalAtom0, targetPreTerminalAtom0, targetAtoms[0]);

                if(queryAtoms[3] == MAX_ATOM_IDX)
                    targetAtoms[3] = molecule_get_last_stereo_bond_ligand(target, targetTerminalAtom1, targetPreTerminalAtom1, targetAtoms[2]);

                if(normalize_bond_stereo(targetAtoms, targetStereo) != queryStereo)
                    return false;
            }
            else
            {
                AtomIdx *bondedAtomList0 = molecule_get_bonded_atom_list(query, queryBondAtoms[0]);
                MolSize bondedAtomListSize0 = molecule_get_bonded_atom_list_size(query, queryBondAtoms[0]);

                AtomIdx *bondedAtomList1 = molecule_get_bonded_atom_list(query, queryBondAtoms[1]);
                MolSize bondedAtomListSize1 = molecule_get_bonded_atom_list_size(query, queryBondAtoms[1]);

                if(bondedAtomListSize0 < 2 || bondedAtomListSize1 < 2)
                    continue;


                AtomIdx queryAtoms[4];

                int idx = 0;

                for(int i = 0; i < bondedAtomListSize0; i++)
                    if(bondedAtomList0[i] != queryBondAtoms[1])
                        queryAtoms[idx++] = bondedAtomList0[i];

                if(bondedAtomListSize0 == 2)
                    queryAtoms[idx++] = MAX_ATOM_IDX;


                for(int i = 0; i < bondedAtomListSize1; i++)
                    if(bondedAtomList1[i] != queryBondAtoms[0])
                        queryAtoms[idx++] = bondedAtomList1[i];

                if(bondedAtomListSize1 == 2)
                    queryAtoms[idx] = MAX_ATOM_IDX;

                sort_bond_atoms(queryAtoms);


                AtomIdx targetAtoms[4] = { -1, -1, -1, -1 };

                for(int i = 0; i < 4; i++)
                    if(queryAtoms[i] != MAX_ATOM_IDX)
                        targetAtoms[i] = vf2state->queryCore[queryAtoms[i]];

                if(queryAtoms[1] == MAX_ATOM_IDX)
                    targetAtoms[1] = molecule_get_last_stereo_bond_ligand(target, targetBondAtom0, targetBondAtom1, targetAtoms[0]);

                if(queryAtoms[3] == MAX_ATOM_IDX)
                    targetAtoms[3] = molecule_get_last_stereo_bond_ligand(target, targetBondAtom1, targetBondAtom0, targetAtoms[2]);

                if(normalize_bond_stereo(targetAtoms, targetStereo) != queryStereo)
                    return false;
            }
        }
    }


    return true;
}


static inline bool vf2state_is_match_valid(const VF2State *restrict vf2state)
{
    const Molecule *restrict query = vf2state->query;
    const Molecule *restrict target = vf2state->target;

    int queryAtomCount = query->atomCount;

    int queryBondCount = query->bondCount;
    int targetBondCount = target->bondCount;

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

    if(molecule_has_restH_flags(query))
    {
        for(AtomIdx i = 0; i < queryAtomCount; i++)
        {
            if(molecule_get_atom_restH_flag(query, i) == true)
            {
                AtomIdx targetAtomIdx = vf2state->queryCore[i];

                int queryConnectivityCount = 0;
                int targetConnectivityCount = 0;

                for(BondIdx b = 0; b < queryBondCount; b++)
                    if(molecule_bond_contains(query, b, i) &&
                            !(molecule_get_atom_number(query, molecule_get_other_bond_atom(query, b, i)) == H_ATOM_NUMBER))
                        queryConnectivityCount++;

                for(BondIdx b = 0; b < targetBondCount; b++)
                    if(molecule_bond_contains(target, b, targetAtomIdx) &&
                            !(molecule_get_atom_number(target, molecule_get_other_bond_atom(target, b, targetAtomIdx)) == H_ATOM_NUMBER))
                        targetConnectivityCount++;

                if(targetConnectivityCount > queryConnectivityCount)
                    return false;
            }
        }
    }


    if(unlikely(vf2state->stereoMode == STEREO_STRICT) && !vf2state_is_stereo_valid(vf2state))
            return false;


    return true;
}


static inline bool vf2state_match_core(VF2State *restrict vf2state)
{
    while(true)
    {
        recursion_entry:

        if(unlikely(vf2state->coreLength == vf2state->query->atomCount))
        {
            if(vf2state_is_match_valid(vf2state))
                return true;

            goto recursion_return;
        }

        if(!vf2state_next_query(vf2state))
            goto recursion_return;


        while(vf2state_next_target(vf2state))
        {
            if(unlikely(--vf2state->counter == 0))
                return false;

            if(vf2state_is_feasible_pair(vf2state))
            {
                vf2state_add_pair(vf2state);
                goto recursion_entry;

                recursion_return:

                if(vf2state->coreLength == 0)
                    return false;

                vf2state_undo_add_pair(vf2state);
            }
        }

        goto recursion_return;
    }

    return false;
}


static inline bool vf2state_match(VF2State *restrict vf2state, const Molecule *restrict target, void *memory, int limit)
{
    vf2state->counter = limit > 0 ? limit : (uint) -1;

    if(likely(vf2state->graphMode != GRAPH_EXACT))
    {
        if(vf2state->query->heavyAtomCount + vf2state->query->hydrogenAtomCount > target->heavyAtomCount + target->hydrogenAtomCount)
            return false;

        if(vf2state->queryAtomCount > target->atomCount)
            return false;

        if(vf2state->query->bondCount > target->bondCount)
            return false;
    }
    else
    {
        if(vf2state->query->heavyAtomCount != target->heavyAtomCount)
            return false;

        if(vf2state->query->hydrogenAtomCount != target->hydrogenAtomCount)
            return false;

        if(vf2state->query->heavyBondCount != target->heavyBondCount)
            return false;

        if(vf2state->query->hydrogenBondCount != target->hydrogenBondCount)
            return false;
    }

    int targetAtomCount = target->atomCount;

    vf2state->target = target;
    vf2state->targetAtomCount = targetAtomCount;
    vf2state->coreLength = 0;

    vf2state->targetCore = (AtomIdx *) memory;

    for(int i = 0; i < targetAtomCount; i++)
        vf2state->targetCore[i] = UNDEFINED_CORE;

    for(int i = 0; i < vf2state->queryAtomCount; i++)
        vf2state->queryCore[i] = UNDEFINED_CORE;

    return vf2state_match_core(vf2state);
}

#endif /* ISOMORPHISM_H__ */
