#include <string.h>
#include <jni.h>
#include <stdbool.h>
#include "ichitime.h"
#include "ichister.h"
#include "ichi.h"
#include "strutil.h"
#include "inchi_api.h"
#include "readinch.h"


int inp2spATOM(inp_ATOM *inp_at, int num_inp_at, sp_ATOM *at);


int (*ConsoleQuit)(void) = NULL;
int (*UserAction)(void)  = NULL;

static jmethodID  setStereoAtomsMethod;
static jmethodID  setStereoBondsMethod;
static jmethodID  setAlternatingBondsMethod;
static jmethodID  setTautomericGroupMethod;


static int process_component(JNIEnv *env, jobject object, inp_ATOM *inp_at, int num_atoms, int mode, bool tautomer)
{
    BCN Bcn;
    ATOM_SIZES s;
    CANON_STAT CS;
    CANON_GLOBALS CG;
    T_GROUP_INFO group_info;
    T_GROUP_INFO *group_info_p = tautomer ? &group_info : NULL;

    memset(&Bcn, 0, sizeof(BCN));
    memset(&s, 0, sizeof(ATOM_SIZES));
    memset(&CS, 0, sizeof(CANON_STAT));
    memset(&CG, 0, sizeof(CANON_GLOBALS));
    memset(&group_info, 0, sizeof(T_GROUP_INFO));

    INCHI_CLOCK ic;
    memset(&ic, 0, sizeof(ic));

    sp_ATOM  *at;
    int ret=0;

    if(!(at = (sp_ATOM *) inchi_malloc(num_atoms * sizeof(*at))))
        ret = -1;

    if(ret)
        goto exit_function;

    /* remove hydrogens */
    group_info.tni.nNumRemovedExplicitH = num_atoms;
    num_atoms = remove_terminal_HDT(num_atoms, inp_at, 1 /* FIX_TERM_H_CHRG_BUG */);
    group_info.tni.nNumRemovedExplicitH -= num_atoms;
    add_DT_to_num_H(num_atoms, inp_at);
    int num_at_tg =  num_atoms;

    /* mark rings */
    MarkRingSystemsInp(inp_at, num_atoms, 0);

    /* normalization */
    group_info.bTautFlags = tautomer ? TG_FLAG_TEST_TAUT__ATOMS | TG_FLAG_KETO_ENOL_TAUT | TG_FLAG_1_5_TAUT : 0;

    if((ret = mark_alt_bonds_and_taut_groups(&ic, &CG, inp_at, NULL, num_atoms, NULL /*&ulMaxTime*/, group_info_p, &group_info.bTautFlags, &group_info.bTautFlagsDone, 0, NULL)) < 0)
        goto exit_function; /* out of RAM or other normalization problem */

    /* create internal sp_ATOM at[] out of inp_at */
    inp2spATOM(inp_at, num_atoms, at);

    /* set stereo parities to at[] */
    int bPointedEdgeStereo = /* PES_BIT_POINT_EDGE_STEREO | */ STEREO_WEDGE_ONLY | PES_BIT_PHOSPHINE_STEREO | PES_BIT_ARSINE_STEREO | PES_BIT_FIX_SP3_BUG;

    ret = set_stereo_parity(&CG, inp_at, at, num_atoms, group_info.tni.nNumRemovedExplicitH, &s.nMaxNumStereoAtoms, &s.nMaxNumStereoBonds, mode, bPointedEdgeStereo, AB_PARITY_UNDF /*AB_PARITY_UNKN*/);

    if(RETURNED_ERROR(ret))
        goto exit_function; /* stereo bond error */

    s.bMayHaveStereo = s.nMaxNumStereoAtoms || s.nMaxNumStereoBonds;

    /* mark isotopic atoms and atoms that have non-tautomeric isotopic terminal hydrogen atoms 1H, 2H(D), 3H(T) */
    s.num_isotopic_atoms = set_atom_iso_sort_keys(num_atoms, at, group_info_p, &s.bHasIsotopicTautGroups);


    /* prepare tautomeric (if no tautomerism found then prepare non-tautomeric) structure for canonicalizaton */
    ret = tautomer ? s.nLenLinearCTTautomer = CountTautomerGroups(at, num_atoms, group_info_p) : 0;

    if(RETURNED_ERROR(ret))
        goto exit_function;

    if(s.nLenLinearCTTautomer > 0)
    {
        num_at_tg = num_atoms + group_info.num_t_groups;

        /*  mark isotopic tautomer groups: calculate t_group->iWeight */
        s.nLenLinearCTIsotopicTautomer = set_tautomer_iso_sort_keys(group_info_p);

        if(s.nLenLinearCTIsotopicTautomer < 0)
            s.nLenLinearCTIsotopicTautomer = 0; /* ??? -error cannot happen- error has happened; no breakpoint here */
    }

    GetCanonLengths(num_atoms, at, &s, group_info_p);

    int bHasIsotopicAtoms = s.num_isotopic_atoms > 0 || s.bHasIsotopicTautGroups > 0 || (s.nLenIsotopicEndpoints > 1 && (group_info.bTautFlagsDone & (TG_FLAG_FOUND_ISOTOPIC_H_DONE|TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE))) ||
            (s.nLenLinearCTTautomer > 0 && (group_info.tni.nNumRemovedProtonsIsotopic[0] || group_info.tni.nNumRemovedProtonsIsotopic[1] || group_info.tni.nNumRemovedProtonsIsotopic[2]));

    if(!bHasIsotopicAtoms)
    {
        if(group_info.nIsotopicEndpointAtomNumber)
            group_info.nIsotopicEndpointAtomNumber[0] = inchi_min(1, group_info.nIsotopicEndpointAtomNumber[0]);

        memset(group_info.num_iso_H, 0, sizeof(group_info.num_iso_H));
        memset(group_info.tni.nNumRemovedProtonsIsotopic, 0, sizeof(group_info.tni.nNumRemovedProtonsIsotopic));

        group_info.bTautFlagsDone &= ~(TG_FLAG_FOUND_ISOTOPIC_H_DONE | TG_FLAG_FOUND_ISOTOPIC_ATOM_DONE);
        group_info.bIgnoreIsotopic = 1;

        s.bIgnoreIsotopic = 1;
        s.nLenIsotopic = 0;
        s.nLenIsotopicEndpoints = 0;
        s.nLenLinearCTIsotopicTautomer = 0;
    }

    INCHI_MODE nMode = ((s.nLenLinearCTTautomer == 0) ? CANON_MODE_CT : CANON_MODE_TAUT) /*| REQ_MODE_DIFF_UU_STEREO | CMODE_NO_ALT_SBONDS*/;
    nMode |= bHasIsotopicAtoms ? CANON_MODE_ISO : 0;
    nMode |= bHasIsotopicAtoms && s.bMayHaveStereo ? CANON_MODE_ISO_STEREO : 0;
    nMode |= !bHasIsotopicAtoms && s.bMayHaveStereo ? CANON_MODE_STEREO : 0;

    /* obtain all non-stereo canonical numberings */
    if((ret = GetBaseCanonRanking(&ic, num_atoms, num_at_tg, (sp_ATOM* []) { tautomer ? NULL : at, tautomer ? at : NULL }, group_info_p, &s - (tautomer != 0), &Bcn, NULL /*&ulMaxTime*/, &CG, 1 /* FIX_ISO_FIXEDH_BUG */, num_atoms <= MAX_ATOMS)) < 0)
        goto exit_function; /*  program error */

    if(ret = AllocateCS(&CS, num_atoms, num_at_tg, s.nLenCT, s.nLenCTAtOnly, s.nLenLinearCTStereoDble, s.nMaxNumStereoBonds, s.nLenLinearCTStereoCarb, s.nMaxNumStereoAtoms, s.nLenLinearCTTautomer, s.nLenLinearCTIsotopicTautomer, s.nLenIsotopic, nMode, &Bcn))
        goto exit_function;

    /*  settings */
    CS.lNumDecreasedCT      = -1;
    CS.bDoubleBondSquare    = DOUBLE_BOND_NEIGH_LIST ? 2 : 0;  /*  2 => special mode */
    CS.bIgnoreIsotopic      = s.bIgnoreIsotopic;
    CS.ulTimeOutTime        = Bcn.ulTimeOutTime;
    CS.t_group_info         = group_info_p;

    /* the last canonicalization step */
    CS.NeighList  = NULL;
    CS.pBCN       = &Bcn;

    CANON_STAT CS2 = CS;
    ret = Canon_INChI(&ic, num_atoms, num_at_tg, at, &CS, &CG, nMode, tautomer);

    /* export atom stereo */
    AT_RANK *nCanonOrdStereo = bHasIsotopicAtoms ? CS.nCanonOrdIsotopicStereo : CS.nCanonOrdStereo;
    AT_STEREO_CARB *LinearCTStereoCarb = bHasIsotopicAtoms ? CS.LinearCTIsotopicStereoCarb : CS.LinearCTStereoCarb;
    AT_STEREO_DBLE *LinearCTStereoDble = bHasIsotopicAtoms ? CS.LinearCTIsotopicStereoDble : CS.LinearCTStereoDble;
    int nLenLinearCTStereoCarb = bHasIsotopicAtoms ? CS.nLenLinearCTIsotopicStereoCarb : CS.nLenLinearCTStereoCarb;
    int nLenLinearCTStereoDble = bHasIsotopicAtoms ? CS.nLenLinearCTIsotopicStereoDble : CS.nLenLinearCTStereoDble;

    int stereo_carb_len = 0;
    int stereo_dble_len = 0;
    int i, j;


    int bond_count = 0;

    for(i = 0; i < num_atoms; i++)
        bond_count += tautomer ? at[i].valence : 0;

    int max_num_endpoints = 0;

    for(i = 0; i < group_info.num_t_groups; i++)
        if(max_num_endpoints < group_info.t_group[i].nNumEndpoints)
            max_num_endpoints = group_info.t_group[i].nNumEndpoints;

    jshort *stereo_carb_buffer = (jshort *) inchi_malloc((nLenLinearCTStereoCarb + nLenLinearCTStereoDble) * 2 * sizeof(jshort));
    jshort *stereo_dble_buffer = (jshort *) inchi_malloc(3 * nLenLinearCTStereoDble * sizeof(jshort));
    jshort *altern_bond_buffer = bond_count ? (jshort *) inchi_malloc(2 * bond_count * sizeof(jshort)) : NULL;
    jshort *t_group_buffer = group_info.num_t_groups ? (jshort *) inchi_malloc((2 + max_num_endpoints) * sizeof(jshort)) : NULL;

    if(!stereo_carb_buffer || !stereo_dble_buffer || bond_count && !altern_bond_buffer || group_info.num_t_groups && !t_group_buffer)
    {
        ret = -1;
        goto java_exception;
    }

    for(i = 0; i < nLenLinearCTStereoCarb; i++)
    {
        int at_num = nCanonOrdStereo[LinearCTStereoCarb[i].at_num - 1];
        stereo_carb_buffer[stereo_carb_len++] = at[at_num].orig_at_number;
        stereo_carb_buffer[stereo_carb_len++] = bHasIsotopicAtoms ? at[at_num].parity2 : at[at_num].parity;
    }

    for(i = 0; i < nLenLinearCTStereoDble; i++)
    {
        int at_num1 = nCanonOrdStereo[(int) LinearCTStereoDble[i].at_num1 - 1];
        int at_num2 = nCanonOrdStereo[(int) LinearCTStereoDble[i].at_num2 - 1];

        int at_num1_parity = bHasIsotopicAtoms ? at[at_num1].parity2 : at[at_num1].parity;
        int at_num2_parity = bHasIsotopicAtoms ? at[at_num2].parity2 : at[at_num2].parity;


        int n1 = 0, n2 = 0;

        for(; (bHasIsotopicAtoms ? at[at_num1].stereo_bond_neighbor2[n1] : at[at_num1].stereo_bond_neighbor[n1]) != at_num2 + 1; n1++);
        for(; (bHasIsotopicAtoms ? at[at_num2].stereo_bond_neighbor2[n2] : at[at_num2].stereo_bond_neighbor[n2]) != at_num1 + 1; n2++);

        int at_num1_bond_parity = PARITY_VAL((bHasIsotopicAtoms ? at[at_num1].stereo_bond_parity2 : at[at_num1].stereo_bond_parity)[n1]);
        int at_num2_bond_parity = PARITY_VAL((bHasIsotopicAtoms ? at[at_num2].stereo_bond_parity2 : at[at_num2].stereo_bond_parity)[n2]);

        int parity = 0;

        if(at_num1_parity == AB_PARITY_UNKN || at_num2_parity == AB_PARITY_UNKN || at_num1_bond_parity == AB_PARITY_UNKN || at_num2_bond_parity == AB_PARITY_UNKN)
        {
            parity = AB_PARITY_UNKN;
        }
        else if(at_num1_parity == AB_PARITY_UNDF || at_num2_parity == AB_PARITY_UNDF || at_num1_bond_parity == AB_PARITY_UNDF || at_num2_bond_parity == AB_PARITY_UNDF)
        {
            parity = AB_PARITY_UNDF;
        }
        else if(PARITY_WELL_DEF(at_num1_parity) && PARITY_WELL_DEF(at_num2_parity))
        {
            parity += at_num1_parity + at_num2_parity + at[at_num1].valence + at[at_num2].valence;
            parity += bHasIsotopicAtoms ? (int)at[at_num1].stereo_bond_ord2[n1] + (int)at[at_num2].stereo_bond_ord2[n2] : (int) at[at_num1].stereo_bond_ord[n1] + (int) at[at_num2].stereo_bond_ord[n2];

            int z_prod = bHasIsotopicAtoms ? at[at_num1].stereo_bond_z_prod2[n1] : at[at_num1].stereo_bond_z_prod[n1];

            if(MIN_DOT_PROD > abs(z_prod))
                parity = AB_PARITY_UNDF;
            else
                parity = (z_prod > 0)? 2 - parity % 2 : 1 + parity % 2;
        }


        int cumulene_len = BOND_CHAIN_LEN((bHasIsotopicAtoms ? at[at_num1].stereo_bond_parity2 : at[at_num1].stereo_bond_parity)[0]);

        if(cumulene_len)
        {
            int next_at_num = at[at_num1].neighbor[(bHasIsotopicAtoms ? at[at_num1].stereo_bond_ord2 : at[at_num1].stereo_bond_ord)[0]];

            int c;
            for(c = cumulene_len/2; c > 0 && 2==at[next_at_num].valence; c--)
            {
                int next_neighbor = (at_num1 == at[next_at_num].neighbor[0]);
                at_num1 = next_at_num;
                next_at_num = at[next_at_num].neighbor[next_neighbor];
            }

            if(cumulene_len % 2)
            {
                stereo_carb_buffer[stereo_carb_len++] = at[next_at_num].orig_at_number;
                stereo_carb_buffer[stereo_carb_len++] = -parity;
            }
            else
            {
                stereo_dble_buffer[stereo_dble_len++] = at[at_num1].orig_at_number;
                stereo_dble_buffer[stereo_dble_len++] = at[next_at_num].orig_at_number;
                stereo_dble_buffer[stereo_dble_len++] = -parity;
            }
        }
        else
        {
            stereo_dble_buffer[stereo_dble_len++] = at[at_num1].orig_at_number;
            stereo_dble_buffer[stereo_dble_len++] = at[at_num2].orig_at_number;
            stereo_dble_buffer[stereo_dble_len++] = parity;
        }
    }

    jshortArray array;

    if((array = (*env)->NewShortArray(env, stereo_carb_len)) == NULL)
    {
        ret = -1;
        goto java_exception;
    }

    (*env)->SetShortArrayRegion(env, array, 0, stereo_carb_len, stereo_carb_buffer);
    (*env)->CallVoidMethod(env, object, setStereoAtomsMethod, array);

    if((array = (*env)->NewShortArray(env, stereo_dble_len)) == NULL)
    {
        ret = -1;
        goto java_exception;
    }

    (*env)->SetShortArrayRegion(env, array, 0, stereo_dble_len, stereo_dble_buffer);
    (*env)->CallVoidMethod(env, object, setStereoBondsMethod, array);

    if(tautomer)
    {
        int altern_bond_len = 0;

        for(i = 0; i < num_atoms; i++)
        {
            for(j = 0; j < at[i].valence; j++)
            {
                if(i < at[i].neighbor[j] && ((at[i].bond_type[j] & BOND_TYPE_MASK) == BOND_TAUTOM || (at[i].bond_type[j] & BOND_TYPE_MASK) == BOND_ALT12NS))
                {
                    altern_bond_buffer[altern_bond_len++] = at[i].orig_at_number;
                    altern_bond_buffer[altern_bond_len++] = at[at[i].neighbor[j]].orig_at_number;
                }
            }
        }


        if((array = (*env)->NewShortArray(env, altern_bond_len)) == NULL)
        {
            ret = -1;
            goto java_exception;
        }

        (*env)->SetShortArrayRegion(env, array, 0, altern_bond_len, altern_bond_buffer);
        (*env)->CallVoidMethod(env, object, setAlternatingBondsMethod, array);

        for(i = 0; i < group_info.num_t_groups; i++)
        {
            int group_size = 0;

            t_group_buffer[group_size++] = group_info.t_group[i].num[0];
            t_group_buffer[group_size++] = group_info.t_group[i].num[1];

            for(j = 0; j < group_info.t_group[i].nNumEndpoints; j++)
            {
                int at_num = group_info.nEndpointAtomNumber[group_info.t_group[i].nFirstEndpointAtNoPos + j];
                t_group_buffer[group_size++] = at[at_num].orig_at_number;
            }

            if((array = (*env)->NewShortArray(env, group_size)) == NULL)
            {
                ret = -1;
                goto java_exception;
            }

            (*env)->SetShortArrayRegion(env, array, 0, group_size, t_group_buffer);
            (*env)->CallVoidMethod(env, object, setTautomericGroupMethod, array);
        }
    }

java_exception:
    FreeNeighList(CS.NeighList);
    DeAllocateCS(&CS2);

    if(stereo_carb_buffer)
        inchi_free(stereo_carb_buffer);

    if(stereo_dble_buffer)
        inchi_free(stereo_dble_buffer);

    if(altern_bond_buffer)
        inchi_free(altern_bond_buffer);

    if(t_group_buffer)
        inchi_free(t_group_buffer);

exit_function:

    DeAllocBCN(&Bcn);

    if(at)
        inchi_free(at);

    free_t_group_info(&group_info);

    return ret;
}


static int process_molecule(JNIEnv *env, jobject object, inp_ATOM *inp_at, int num_atoms, int mode, bool tautomer)
{
    int num_components = 0;

    ORIG_ATOM_DATA prep_inp_data;
    memset(&prep_inp_data, 0, sizeof(ORIG_ATOM_DATA));

    prep_inp_data.at = inp_at;
    prep_inp_data.num_inp_atoms = num_atoms;


    if((num_components = MarkDisconnectedComponents(&prep_inp_data, 0)) < 0)
        return -1;

    inp_ATOM *component_at = NULL;

    if(!(component_at = (inp_ATOM *) inchi_malloc(num_atoms * sizeof(inp_ATOM))))
        return -1;

    int i;
    int ret = 0;
    for(i = 1; i <= num_components; i++)
    {
        int component_num_atoms = ExtractConnectedComponent(inp_at, num_atoms, i, component_at);

        if(component_num_atoms <= 0)
        {
            ret = -i;
            goto exit_function;
        }

        if((ret = process_component(env, object, component_at, component_num_atoms, mode, tautomer)) < 0)
            goto exit_function;
    }

    exit_function:
        inchi_free(component_at);

    return ret;
}


JNIEXPORT jint JNICALL Java_cz_iocb_elchem_molecule_InChITools_process(JNIEnv *env, jobject object, jbyteArray array, jbyteArray stereo, jint mode, jboolean tautomer)
{
    inp_ATOM *inp_at = (inp_ATOM *) (*env)->GetByteArrayElements(env, array, NULL);
    int num_atoms = (*env)->GetArrayLength(env, array) / sizeof(inp_ATOM);

    if(inp_at == NULL)
        return -1;

    if(stereo != NULL)
    {
        inchi_Stereo0D *stereo0D = (inchi_Stereo0D *) (*env)->GetByteArrayElements(env, stereo, NULL);
        int num_stereos = (*env)->GetArrayLength(env, stereo) / sizeof(inchi_Stereo0D);

        if(stereo0D == NULL)
        {
            (*env)->ReleaseByteArrayElements(env, array, (jbyte *) inp_at, JNI_ABORT);
            return -1;
        }

        int err = 0;

        if(Extract0DParities(inp_at, num_atoms, stereo0D, num_stereos, NULL, &err, AB_PARITY_UNDF /*AB_PARITY_UNKN*/))
        {
            (*env)->ReleaseByteArrayElements(env, stereo, (jbyte *) stereo0D, JNI_ABORT);
            (*env)->ReleaseByteArrayElements(env, array, (jbyte *) inp_at, JNI_ABORT);
            return -1;
        }

        (*env)->ReleaseByteArrayElements(env, stereo, (jbyte *) stereo0D, JNI_ABORT);
    }

    int ret = process_molecule(env, object, inp_at, num_atoms, mode, tautomer);

    (*env)->ReleaseByteArrayElements(env, array, (jbyte *) inp_at, JNI_ABORT);

    return ret;
}


JNIEXPORT void JNICALL Java_cz_iocb_elchem_molecule_InChITools_init(JNIEnv *env, jclass iclass)
{
    jclass inchiToolsClass = (*env)->FindClass(env, "cz/iocb/elchem/molecule/InChITools");
    setStereoAtomsMethod = (*env)->GetMethodID(env, inchiToolsClass, "setStereoAtoms", "([S)V");
    setStereoBondsMethod = (*env)->GetMethodID(env, inchiToolsClass, "setStereoBonds", "([S)V");
    setAlternatingBondsMethod = (*env)->GetMethodID(env, inchiToolsClass, "setAlternatingBonds", "([S)V");
    setTautomericGroupMethod = (*env)->GetMethodID(env, inchiToolsClass, "setTautomericGroup", "([S)V");
}
