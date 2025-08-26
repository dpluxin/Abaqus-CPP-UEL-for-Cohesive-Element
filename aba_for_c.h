//=============================================================================
// COPYRIGHT DASSAULT SYSTEMES 2001-2013
/**
* @CAA2Level L0
* @CAA2Usage U0
*/
//=============================================================================
/****************************************************************************
 *                                                                          *
 *            ABAQUS FORTRAN-to-C calling convensions macros                *
 *                                                                          *
 ****************************************************************************/

#ifndef FOR_C_H

#ifndef ABA_FOR_C_H
#define ABA_FOR_C_H

#if defined(ABQ_WIN86_64) || defined(ABQ_WIN86_32)
//#include <omi_types.h> copied as followed

// =====================================================================

    //=============================================================================
    // COPYRIGHT DASSAULT SYSTEMES 2001-2012
    /**
    * @CAA2Level L0
    * @CAA2Usage U0
    */
    //=============================================================================
    #ifndef omi_TYPES_h
    #define omi_TYPES_h

    #if defined(_WINDOWS_SOURCE) || defined (_WIN64)
    #include <stddef.h> /* intptr_t */
    #include <BaseTsd.h> /* SSIZE_T */
    #else
    #include <inttypes.h> /* intptr_t */
    #endif

    /* Begin local includes */

    /* Convenience typedef */
    #if defined(_WINDOWS_SOURCE) || defined (_WIN64)
    typedef SSIZE_T ssize_t;
    #endif

    /* The int64 and uint64 data types should be used only for special
       purposes, such as addressing large files. Discuss any proposed use
       of this type with Product Architecture. */

    /* uint64: */
    #if defined _WINDOWS_SOURCE
    typedef unsigned __int64 uint64;
    #elif defined _LINUX_SOURCE
    typedef unsigned long long uint64;
    #endif

    /* int64: */
    #if defined _WINDOWS_SOURCE
    typedef __int64 int64;
    #elif defined _LINUX_SOURCE
    typedef long long int64;
    #endif

    /* Convenience typedef to save typing
       Only use our typedef if it has not already been set by the system */
    #if (defined(_LINUX_SOURCE) && defined(__USE_MISC) && defined(_SYS_TYPES_H))
    #else
    /* uint */
    typedef unsigned int uint;
    #endif

    /* Convenience typedef to save typing
       Only use our typedef if it has not already been set by the system */
    #if (defined(_LINUX_SOURCE) && defined(__USE_MISC) && defined(_SYS_TYPES_H))
    #else
    /* ushort */
    typedef unsigned short int ushort;
    #endif

    /* Convenience typedef to explicitely indicate that char is used as int */
    typedef signed char schar;
    typedef unsigned char uchar;

    /* Special typedef to allow us to deal with MS special types where we need a
       32 bit type that is disguised as a long. This type should not show up in
       many places. If you want to use it get the blessing from architecture. */
    typedef int omiInt32_t;
    typedef unsigned int omiUint32_t;

    /* Special typedef to deal with long in External interfaces that we do not
       control or in interfaces we do not want to change. This type should not
       show up in many places. If you want to use it get the blessing from
       architecture.
       Note that these types are not portable:
       4 byte types on 32 bit systems;
       8 byte types on 64 bit UNIX and Linux;
       4 byte types on 64 bit windows;
       Never ever use these for a pointer in non platform specific code. */
    #if defined _WINDOWS_SOURCE
    typedef int omi_interfLong_t;
    typedef unsigned int omi_interfULong_t;
    #elif defined _LINUX_SOURCE
    typedef  long long omi_interfLong_t;
    typedef unsigned long long omi_interfULong_t;
    #endif

    #endif

// end copy of omi_types.h
// =====================================================================

#endif

// #include <omi_for_types.h> copied as followed
// =====================================================================
    //=============================================================================
    // COPYRIGHT DASSAULT SYSTEMES 2001-2014
    /**
    * @CAA2Level L0
    * @CAA2Usage U0
    */
    //=============================================================================
    #ifndef omi_for_types_h
    #define omi_for_types_h


    /* HACK around IBM compiler bug to avoid hugh build */
    /* Begin local includes */
    //#include <omi_types.h>

    /* FORTRAN Data type mapping */

    typedef int Integer;
    typedef char          Char;
    typedef unsigned int  UInt;
    typedef int           Int;

    typedef double        Double;
    typedef double        Real64;
    typedef float         Single;

    #if defined(DPX)
    typedef double        Real;
    #else
    typedef float         Real;
    #endif

    /* NOTE if we have a system with truly large memory space, i.e. > LONG_MAX
       we potentially will run into trouble and we will have to change Address
       to be uintptr_t */
    typedef intptr_t Address;

    #endif
// end copy of omi_for_types.h
// =====================================================================


#if defined FOR_TRAIL

#define FOR_NAME(lower_case_name,UPPER_CASE_NAME)  lower_case_name ## _
#define CALL_NAME(lower_case_name,UPPER_CASE_NAME) lower_case_name ## _

#elif defined (ABQ_WIN86_32) || defined (ABQ_WIN86_64)

#define FOR_NAME(lower_case_name,UPPER_CASE_NAME)  _cdecl lower_case_name
#define CALL_NAME(lower_case_nane,UPPER_CASE_NAME) lower_case_name

#else

#define FOR_NAME(lower_case_name,upper_case_name)  lower_case_name
#define CALL_NAME(lower_case_name,upper_case_name) lower_case_name

#endif

/******************************************************************************
 * The rules for passing type CHARACTER from FORTRAN to C are as follows:     *
 *   1) All type CHARACTER arguments must appear before any other arguments   *
 *      in the list. The parameter should be declared as type "Character".    *
 *   2) A Length parameter for each Character parameter must be appended to   *
 *	the C parameter list. The name of the Length parameter should be      *
 *	constructed by adding "_len" to the Character parameter name.	      *
 *   3) The CHARACTER data must be accessed via the PTR(cname) macro, which   *
 *	defines a pointer to cname.					      *
 *   4) The length of the CHARACTER string may be obtained via the	      *
 *	LEN(cname,cnname_len) macro, which defines the length of cname.	      *
 * This allows capture of the ancilliary length arguments generated by UNIX   *
 * and C/370 systems. The VMS implementation ignores these extra parameters.  *
 ******************************************************************************/


#if defined(ABQ_WIN86_64)
typedef char *Character;
typedef int64 Length;
#define PTR(id) id
#define LEN(id,len) len
#define GETLEN(id) id ## _len
#define CHNAME(id) Character id
#define CHNAME_C(id) id
#define CHLEN(id) , const Length id ## _len
#define CHLEN_C(id) , id ## _len
#elif defined(ABQ_WIN86_32)
typedef char *Character;
typedef unsigned int Length;
#define PTR(id) id
#define LEN(id,len) len
#define GETLEN(id) id ## _len
#define CHNAME(id) Character id
#define CHNAME_C(id) id
#define CHLEN(id) , const Length id ## _len
#define CHLEN_C(id) , id ## _len
#elif defined(_WINDOWS_SOURCE)
typedef char *Character;
typedef unsigned int Length;
#define PTR(id) id
#define LEN(id,len) len
#define GETLEN(id) id ## _len
#define CHNAME(id) Character id, const Length id ## _len
#define CHNAME_C(id) id, id ## _len
#define CHLEN(id)
#define CHLEN_C(id)
#else
typedef char *Character;
typedef int Length;
#define PTR(id) id
#define LEN(id,len) len
#define GETLEN(id) id ## _len
#define CHNAME(id) Character id
#define CHNAME_C(id) id
#define CHLEN(id) , const Length id ## _len
#define CHLEN_C(id) , id ## _len
#endif


#endif // ABA_FOR_C_H

#endif // FOR_C_H
