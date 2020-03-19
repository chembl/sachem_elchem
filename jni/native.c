#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <jni.h>
#include "molecule.h"
#include "isomorphism.h"


static jclass byteBufferClass;
static jclass outOfMemoryErrorClass;
static jclass iterationLimitExceededExceptionClass;
static jmethodID allocateDirectMethod;
static jmethodID outOfMemoryErrorConstructor;
static jmethodID iterationLimitExceededExceptionConstructor;


JNIEXPORT jobject JNICALL Java_cz_iocb_elchem_molecule_NativeIsomorphism_create(JNIEnv *env, jclass clazz,
        jbyteArray queryArray, jbooleanArray restHArray, jint searchMode, jint chargeMode, jint isotopeMode, jint stereoMode)
{
    uint8_t *query = (uint8_t *) (*env)->GetByteArrayElements(env, queryArray, NULL);
    uint8_t *restH = (uint8_t *) (*env)->GetBooleanArrayElements(env, restHArray, NULL);

    if(unlikely(query == NULL || restH == NULL))
    {
        if(query != NULL)
            (*env)->ReleaseByteArrayElements(env, queryArray, (jbyte *) query, JNI_ABORT);

        return NULL;
    }

    bool extended = molecule_is_extended_search_needed(query, searchMode != SEARCH_EXACT, chargeMode != CHARGE_IGNORE, isotopeMode != ISOTOPE_IGNORE);

    size_t isosize = vf2state_mem_size(query, extended);
    size_t molsize = molecule_mem_size(query, restH, extended, chargeMode != CHARGE_IGNORE, isotopeMode != ISOTOPE_IGNORE, stereoMode != STEREO_IGNORE, false, false);

    jobject buffer = (*env)->CallStaticObjectMethod(env, byteBufferClass, allocateDirectMethod, (jint) (isosize + molsize));

    if(likely(!(*env)->ExceptionCheck(env)))
    {
        void *memory = (*env)->GetDirectBufferAddress(env, buffer);
        Molecule *molecule = molecule_create(memory + isosize, query, restH, extended, chargeMode != CHARGE_IGNORE, isotopeMode != ISOTOPE_IGNORE, stereoMode != STEREO_IGNORE, false, false);
        vf2state_create(memory, molecule, searchMode, chargeMode, isotopeMode, stereoMode);
    }

    (*env)->ReleaseByteArrayElements(env, queryArray, (jbyte *) query, JNI_ABORT);
    (*env)->ReleaseBooleanArrayElements(env, restHArray, (jboolean *) restH, JNI_ABORT);

    return buffer;
}


JNIEXPORT jfloat JNICALL Java_cz_iocb_elchem_molecule_NativeIsomorphism_match(JNIEnv *env, jclass clazz, jobject buffer,
        jbyteArray targetArray, jint limit)
{
    VF2State *isomorphism = (VF2State *) (*env)->GetDirectBufferAddress(env, buffer);
    uint8_t *target = (uint8_t *) (*env)->GetByteArrayElements(env, targetArray, NULL);

    if(unlikely(target == NULL))
    {
        jobject error = (*env)->NewObject(env, outOfMemoryErrorClass, outOfMemoryErrorConstructor);

        if(!(*env)->ExceptionCheck(env))
            (*env)->Throw(env, error);

        return -INFINITY;
    }


    bool extend = !isomorphism->query->extended && isomorphism->query->hydrogenAtomCount &&
            (molecule_has_multivalent_hydrogen(target) ||
                    (isomorphism->searchMode == SEARCH_EXACT && (
                            (isomorphism->chargeMode == CHARGE_DEFAULT_AS_UNCHARGED && molecule_has_charged_hydrogen(target)) ||
                            (isomorphism->isotopeMode == ISOTOPE_DEFAULT_AS_STANDARD && molecule_has_hydrogen_isotope(target)))));

    size_t targetsize = molecule_mem_size(target, NULL, extend || isomorphism->query->extended,
            isomorphism->chargeMode != CHARGE_IGNORE, isomorphism->isotopeMode != ISOTOPE_IGNORE,
            isomorphism->stereoMode != STEREO_IGNORE, isomorphism->chargeMode == CHARGE_DEFAULT_AS_UNCHARGED,
            isomorphism->isotopeMode == ISOTOPE_DEFAULT_AS_STANDARD);

    size_t isosize = extend ? vf2state_extended_mem_size(isomorphism->query) : 0;
    size_t molsize = extend ? molecule_extended_mem_size(isomorphism->query) : 0;
    size_t matchsize = vf2state_match_mem_size(target, isomorphism->query->extended || extend);

    void *molmemory = malloc(targetsize + isosize + molsize + matchsize);

    if(unlikely(molmemory == NULL))
    {
        (*env)->ReleaseByteArrayElements(env, targetArray, (jbyte *) target, JNI_ABORT);

        jobject error = (*env)->NewObject(env, outOfMemoryErrorClass, outOfMemoryErrorConstructor);

        if(!(*env)->ExceptionCheck(env))
            (*env)->Throw(env, error);

        return -INFINITY;
    }

    if(extend)
    {
        Molecule *query = molecule_extend(molmemory + targetsize + isosize, isomorphism->query);
        isomorphism = vf2state_create(molmemory + targetsize, query, isomorphism->searchMode, isomorphism->chargeMode, isomorphism->isotopeMode, isomorphism->stereoMode);
    }

    Molecule *molecule = molecule_create(molmemory, target, NULL, isomorphism->query->extended,
            isomorphism->chargeMode != CHARGE_IGNORE, isomorphism->isotopeMode != ISOTOPE_IGNORE,
            isomorphism->stereoMode != STEREO_IGNORE, isomorphism->chargeMode == CHARGE_DEFAULT_AS_UNCHARGED,
            isomorphism->isotopeMode == ISOTOPE_DEFAULT_AS_STANDARD);

    (*env)->ReleaseByteArrayElements(env, targetArray, (jbyte *) target, JNI_ABORT);

    float score = NAN;

    if(vf2state_match(isomorphism, molecule, molmemory + targetsize + isosize + molsize, limit))
    {
        double heavyAtom = molecule->heavyAtomCount ? isomorphism->query->heavyAtomCount / (double) molecule->heavyAtomCount : 1.0;
        double hydrogenAtom = molecule->hydrogenAtomCount ? isomorphism->query->hydrogenAtomCount / (double) molecule->hydrogenAtomCount : 1.0;
        double heavyBond = molecule->heavyBondCount ? isomorphism->query->heavyBondCount / (double) molecule->heavyBondCount : 1.0;
        double hydrogenBond = molecule->hydrogenBondCount ? isomorphism->query->hydrogenBondCount / (double) molecule->hydrogenBondCount : 1.0;

        score = heavyAtom / 2 + heavyBond / 4 + hydrogenAtom / 8 + hydrogenBond / 16;
    }
    else if(unlikely(isomorphism->counter == 0))
    {
        jobject exception = (*env)->NewObject(env, iterationLimitExceededExceptionClass, iterationLimitExceededExceptionConstructor);

        if(!(*env)->ExceptionCheck(env))
            (*env)->Throw(env, exception);

        score = -INFINITY;
    }

    free(molmemory);
    return score;
}


JNIEXPORT void JNICALL Java_cz_iocb_elchem_molecule_NativeIsomorphism_init(JNIEnv *env, jclass clazz)
{
    outOfMemoryErrorClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "java/lang/OutOfMemoryError"));
    outOfMemoryErrorConstructor = (*env)->GetMethodID(env, outOfMemoryErrorClass, "<init>", "()V");

    iterationLimitExceededExceptionClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "cz/iocb/elchem/molecule/NativeIsomorphism$IterationLimitExceededException"));
    iterationLimitExceededExceptionConstructor = (*env)->GetMethodID(env, iterationLimitExceededExceptionClass, "<init>", "()V");

    byteBufferClass = (jclass) (*env)->NewGlobalRef(env, (*env)->FindClass(env, "java/nio/ByteBuffer"));
    allocateDirectMethod = (*env)->GetStaticMethodID(env, byteBufferClass, "allocateDirect", "(I)Ljava/nio/ByteBuffer;");
}
