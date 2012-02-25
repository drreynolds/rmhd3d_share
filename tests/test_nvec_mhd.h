 /******************************************************************
 *                                                                 *
 * File          : nvec_mhd_test.h                                 *
 * Programmers   : Daniel R. Reynolds @ SMU                        *
 * Version of    : 20 April 2004                                   *
 *-----------------------------------------------------------------*
 *  This is the header file for the nvec_mhd_test.c routine, which *
 *  is used to test both the internal nvector_mhd NVector library, *
 *  as well as its corresponding Fortran interface.                *
 *                                                                 *
 *******************************************************************/

#ifndef _NVEC_MHD_TEST_H
#define _NVEC_MHD_TEST_H

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include <sundials_config.h>
#endif

#include "nvector.h"       /* definition of type N_Vector  */
#include "sundialstypes.h" /* definitions of type realtype */

#if defined(F77_FUNC)

#define F_NVECMHDTEST  F77_FUNC(fnvecmhdtest,   FNVECMHDTEST)
#define PRINT_FNVECMHD F77_FUNC(print_fnvecmhd, PRINT_FNVECMHD)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define F_NVECMHDTEST    fnvecmhdtest
#define PRINT_FNVECMHD   print_fnvecmhd

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define F_NVECMHDTEST    FNVECMHDTEST
#define PRINT_FNVECMHD   PRINT_FNVECMHD

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define F_NVECMHDTEST    fnvecmhdtest_
#define PRINT_FNVECMHD   print_fnvecmhd_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define F_NVECMHDTEST    FNVECMHDTEST_
#define PRINT_FNVECMHD   PRINT_FNVECMHD_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define F_NVECMHDTEST    fnvecmhdtest__
#define PRINT_FNVECMHD   print_fnvecmhd__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define F_NVECMHDTEST    FNVECMHDTEST__
#define PRINT_FNVECMHD   PRINT_FNVECMHD__

#endif



/* prototypes for internal functions */
void FortPrint(N_Vector);
void FilePrint(N_Vector, FILE *, int, int, int, int, int, int, int, int);
void F_NVECMHDTEST(realtype *, realtype *, int *, int *, int *, 
		   int *, int *, int *, int *, int *, int *, int *);
void PRINT_FNVECMHD(realtype *);


#endif
