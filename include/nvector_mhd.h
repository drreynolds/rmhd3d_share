/*******************************************************************
 * Daniel R. Reynolds                                              *
 * UC San Diego, Mathematics                                       *
 *-----------------------------------------------------------------*
 * This is the header file for an implementation of an NVECTOR     *
 * package, specifically suited to interface with a supplied       *
 * data structure.  In these vectors, local data consists of an    *
 * arbitrary array of a given length.                              *
 *                                                                 *
 * Part I of this file contains definitions and prototypes of the  *
 * Fortran interface functions for this implementation.            *
 *                                                                 *
 * Part II of this file contains declarations which are specific   *
 * to the particular vector specification in which this version    *
 * of the NVECTOR module is to be used. This includes the          *
 * typedef for the 'content' fields of the Vector.                 *
 *                                                                 *
 * Part III of this file defines accessor macros that allow the    *
 * user to use efficiently the type N_Vector without making        *
 * explicit references to its underlying representation.           *
 *                                                                 *
 * Part IV of this file contains the prototypes for                *
 * initialization and printing routines specific to this           *
 * implementation.                                                 *
 *                                                                 *
 * Part V of this file contains prototypes for the vector kernels  *
 * which operate on the N_Vector.  These prototypes are unique to  *
 * this particular implementation of the vector package, and are   *
 * attached to the generic N_Vector structure for use within the   *
 * various SUNDIALS solvers.                                       *
 *                                                                 *
 * NOTES:                                                          *
 *                                                                 *
 * The definitions of the generic N_Vector structure is in the     *
 * SUNDIALS header file nvector.h.                                 *
 *                                                                 *
 * The definition of the type realtype is in the SUNDIALS header   *
 * file sundialstypes.h and may be changed according to the user's *
 * needs. The SUNDIALS file sundialstypes.h also contains the      *
 * definition for the type booleantype.                            *
 *                                                                 *
 *******************************************************************/


#ifndef _NVECTOR_MHD_H
#define _NVECTOR_MHD_H

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_fnvector.h>



/****************************************************************
 * PART I:                                                      *
 * Fortran callable wrappers for the interface routines, as     *
 * well as prototypes of the supplied Fortran vector operations *
 * to be called by C routines.                                  *
 ****************************************************************/

/* Define Fortran/C interface wrappers */
#if defined(F77_FUNC)

#define FNVMHD_INIT         F77_FUNC(fnvinitmhd,         FNVINITMHD)
#define FNVMHD_PRINT        F77_FUNC(fnvmhdprint,        FNVMHDPRINT)
#define FNVMHD_LINSUM       F77_FUNC(fnvmhdlinsum,       FNVMHDLINSUM)
#define FNVMHD_CONST        F77_FUNC(fnvmhdconst,        FNVMHDCONST)
#define FNVMHD_PROD         F77_FUNC(fnvmhdprod,         FNVMHDPROD)
#define FNVMHD_DIV          F77_FUNC(fnvmhddiv,          FNVMHDDIV)
#define FNVMHD_SCALE        F77_FUNC(fnvmhdscale,        FNVMHDSCALE)
#define FNVMHD_ABS          F77_FUNC(fnvmhdabs,          FNVMHDABS)
#define FNVMHD_INV          F77_FUNC(fnvmhdinv,          FNVMHDINV)
#define FNVMHD_ADDCONST     F77_FUNC(fnvmhdaddconst,     FNVMHDADDCONST)
#define FNVMHD_DOTPROD      F77_FUNC(fnvmhddotprod,      FNVMHDDOTPROD)
#define FNVMHD_MAXNORM      F77_FUNC(fnvmhdmaxnorm,      FNVMHDMAXNORM)
#define FNVMHD_WRMSNORM     F77_FUNC(fnvmhdwrmsnorm,     FNVMHDWRMSNORM)
#define FNVMHD_WRMSNORMMASK F77_FUNC(fnvmhdwrmsnormmask, FNVMHDWRMSNORMMASK)
#define FNVMHD_MIN          F77_FUNC(fnvmhdmin,          FNVMHDMIN)
#define FNVMHD_WL2NORM      F77_FUNC(fnvmhdwl2norm,      FNVMHDWL2NORM)
#define FNVMHD_L1NORM       F77_FUNC(fnvmhdl1norm,       FNVMHDL1NORM)
#define FNVMHD_COMPARE      F77_FUNC(fnvmhdcompare,      FNVMHDCOMPARE)
#define FNVMHD_INVTEST      F77_FUNC(fnvmhdinvtest,      FNVMHDINVTEST)
#define FNVMHD_CONSTRMASK   F77_FUNC(fnvmhdconstrmask,   FNVMHDCONSTRMASK)
#define FNVMHD_MINQUOTIENT  F77_FUNC(fnvmhdminquotient,  FNVMHDMINQUOTIENT)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FNVMHD_INIT         fnvinitmhd
#define FNVMHD_PRINT        fnvmhdprint
#define FNVMHD_LINSUM       fnvmhdlinsum
#define FNVMHD_CONST        fnvmhdconst
#define FNVMHD_PROD         fnvmhdprod
#define FNVMHD_DIV          fnvmhddiv
#define FNVMHD_SCALE        fnvmhdscale
#define FNVMHD_ABS          fnvmhdabs
#define FNVMHD_INV          fnvmhdinv
#define FNVMHD_ADDCONST     fnvmhdaddconst
#define FNVMHD_DOTPROD      fnvmhddotprod
#define FNVMHD_MAXNORM      fnvmhdmaxnorm
#define FNVMHD_WRMSNORM     fnvmhdwrmsnorm
#define FNVMHD_WRMSNORMMASK fnvmhdwrmsnormmask
#define FNVMHD_MIN          fnvmhdmin
#define FNVMHD_WL2NORM      fnvmhdwl2norm
#define FNVMHD_L1NORM       fnvmhdl1norm
#define FNVMHD_COMPARE      fnvmhdcompare
#define FNVMHD_INVTEST      fnvmhdinvtest
#define FNVMHD_CONSTRMASK   fnvmhdconstrmask
#define FNVMHD_MINQUOTIENT  fnvmhdminquotient

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FNVMHD_INIT         FNVINITMHD
#define FNVMHD_PRINT        FNVMHDPRINT
#define FNVMHD_LINSUM       FNVMHDLINSUM
#define FNVMHD_CONST        FNVMHDCONST
#define FNVMHD_PROD         FNVMHDPROD
#define FNVMHD_DIV          FNVMHDDIV
#define FNVMHD_SCALE        FNVMHDSCALE
#define FNVMHD_ABS          FNVMHDABS
#define FNVMHD_INV          FNVMHDINV
#define FNVMHD_ADDCONST     FNVMHDADDCONST
#define FNVMHD_DOTPROD      FNVMHDDOTPROD
#define FNVMHD_MAXNORM      FNVMHDMAXNORM
#define FNVMHD_WRMSNORM     FNVMHDWRMSNORM
#define FNVMHD_WRMSNORMMASK FNVMHDWRMSNORMMASK
#define FNVMHD_MIN          FNVMHDMIN
#define FNVMHD_WL2NORM      FNVMHDWL2NORM
#define FNVMHD_L1NORM       FNVMHDL1NORM
#define FNVMHD_COMPARE      FNVMHDCOMPARE
#define FNVMHD_INVTEST      FNVMHDINVTEST
#define FNVMHD_CONSTRMASK   FNVMHDCONSTRMASK
#define FNVMHD_MINQUOTIENT  FNVMHDMINQUOTIENT

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FNVMHD_INIT         fnvinitmhd_
#define FNVMHD_PRINT        fnvmhdprint_
#define FNVMHD_LINSUM       fnvmhdlinsum_
#define FNVMHD_CONST        fnvmhdconst_
#define FNVMHD_PROD         fnvmhdprod_
#define FNVMHD_DIV          fnvmhddiv_
#define FNVMHD_SCALE        fnvmhdscale_
#define FNVMHD_ABS          fnvmhdabs_
#define FNVMHD_INV          fnvmhdinv_
#define FNVMHD_ADDCONST     fnvmhdaddconst_
#define FNVMHD_DOTPROD      fnvmhddotprod_
#define FNVMHD_MAXNORM      fnvmhdmaxnorm_
#define FNVMHD_WRMSNORM     fnvmhdwrmsnorm_
#define FNVMHD_WRMSNORMMASK fnvmhdwrmsnormmask_
#define FNVMHD_MIN          fnvmhdmin_
#define FNVMHD_WL2NORM      fnvmhdwl2norm_
#define FNVMHD_L1NORM       fnvmhdl1norm_
#define FNVMHD_COMPARE      fnvmhdcompare_
#define FNVMHD_INVTEST      fnvmhdinvtest_
#define FNVMHD_CONSTRMASK   fnvmhdconstrmask_
#define FNVMHD_MINQUOTIENT  fnvmhdminquotient_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FNVMHD_INIT         FNVINITMHD_
#define FNVMHD_PRINT        FNVMHDPRINT_
#define FNVMHD_LINSUM       FNVMHDLINSUM_
#define FNVMHD_CONST        FNVMHDCONST_
#define FNVMHD_PROD         FNVMHDPROD_
#define FNVMHD_DIV          FNVMHDDIV_
#define FNVMHD_SCALE        FNVMHDSCALE_
#define FNVMHD_ABS          FNVMHDABS_
#define FNVMHD_INV          FNVMHDINV_
#define FNVMHD_ADDCONST     FNVMHDADDCONST_
#define FNVMHD_DOTPROD      FNVMHDDOTPROD_
#define FNVMHD_MAXNORM      FNVMHDMAXNORM_
#define FNVMHD_WRMSNORM     FNVMHDWRMSNORM_
#define FNVMHD_WRMSNORMMASK FNVMHDWRMSNORMMASK_
#define FNVMHD_MIN          FNVMHDMIN_
#define FNVMHD_WL2NORM      FNVMHDWL2NORM_
#define FNVMHD_L1NORM       FNVMHDL1NORM_
#define FNVMHD_COMPARE      FNVMHDCOMPARE_
#define FNVMHD_INVTEST      FNVMHDINVTEST_
#define FNVMHD_CONSTRMASK   FNVMHDCONSTRMASK_
#define FNVMHD_MINQUOTIENT  FNVMHDMINQUOTIENT_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FNVMHD_INIT         fnvinitmhd__
#define FNVMHD_PRINT        fnvmhdprint__
#define FNVMHD_LINSUM       fnvmhdlinsum__
#define FNVMHD_CONST        fnvmhdconst__
#define FNVMHD_PROD         fnvmhdprod__
#define FNVMHD_DIV          fnvmhddiv__
#define FNVMHD_SCALE        fnvmhdscale__
#define FNVMHD_ABS          fnvmhdabs__
#define FNVMHD_INV          fnvmhdinv__
#define FNVMHD_ADDCONST     fnvmhdaddconst__
#define FNVMHD_DOTPROD      fnvmhddotprod__
#define FNVMHD_MAXNORM      fnvmhdmaxnorm__
#define FNVMHD_WRMSNORM     fnvmhdwrmsnorm__
#define FNVMHD_WRMSNORMMASK fnvmhdwrmsnormmask__
#define FNVMHD_MIN          fnvmhdmin__
#define FNVMHD_WL2NORM      fnvmhdwl2norm__
#define FNVMHD_L1NORM       fnvmhdl1norm__
#define FNVMHD_COMPARE      fnvmhdcompare__
#define FNVMHD_INVTEST      fnvmhdinvtest__
#define FNVMHD_CONSTRMASK   fnvmhdconstrmask__
#define FNVMHD_MINQUOTIENT  fnvmhdminquotient__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FNVMHD_INIT         FNVINITMHD__
#define FNVMHD_PRINT        FNVMHDPRINT__
#define FNVMHD_LINSUM       FNVMHDLINSUM__
#define FNVMHD_CONST        FNVMHDCONST__
#define FNVMHD_PROD         FNVMHDPROD__
#define FNVMHD_DIV          FNVMHDDIV__
#define FNVMHD_SCALE        FNVMHDSCALE__
#define FNVMHD_ABS          FNVMHDABS__
#define FNVMHD_INV          FNVMHDINV__
#define FNVMHD_ADDCONST     FNVMHDADDCONST__
#define FNVMHD_DOTPROD      FNVMHDDOTPROD__
#define FNVMHD_MAXNORM      FNVMHDMAXNORM__
#define FNVMHD_WRMSNORM     FNVMHDWRMSNORM__
#define FNVMHD_WRMSNORMMASK FNVMHDWRMSNORMMASK__
#define FNVMHD_MIN          FNVMHDMIN__
#define FNVMHD_WL2NORM      FNVMHDWL2NORM__
#define FNVMHD_L1NORM       FNVMHDL1NORM__
#define FNVMHD_COMPARE      FNVMHDCOMPARE__
#define FNVMHD_INVTEST      FNVMHDINVTEST__
#define FNVMHD_CONSTRMASK   FNVMHDCONSTRMASK__
#define FNVMHD_MINQUOTIENT  FNVMHDMINQUOTIENT__

#endif

  /* Declarations of global variables */
  extern N_Vector F2C_CVODE_vec;
  extern N_Vector F2C_KINSOL_vec;

  /* Prototype of initialization function */
  void FNVMHD_INIT(int *code, long int *length, int *ier);
  
  /* Prototypes of the Fortran routines -- these are provided by 
     the user, and define the variables in the supplied vector 
     depending on the physical problem under consideration */
  void FNVMHD_PRINT(realtype*);
  void FNVMHD_LINSUM(realtype*, realtype*, realtype*, realtype*, realtype*);
  void FNVMHD_CONST(realtype*, realtype*);
  void FNVMHD_PROD(realtype*, realtype*, realtype*);
  void FNVMHD_DIV(realtype*, realtype*, realtype*);
  void FNVMHD_SCALE(realtype*, realtype*, realtype*);
  void FNVMHD_ABS(realtype*, realtype*);
  void FNVMHD_INV(realtype*, realtype*);
  void FNVMHD_ADDCONST(realtype*, realtype*, realtype*);
  void FNVMHD_DOTPROD(realtype*, realtype*, realtype*);
  void FNVMHD_MAXNORM(realtype*, realtype*);
  void FNVMHD_WRMSNORM(realtype*, realtype*, realtype*);
  void FNVMHD_WRMSNORMMASK(realtype*, realtype*, realtype*, realtype*);
  void FNVMHD_MIN(realtype*, realtype*);
  void FNVMHD_WL2NORM(realtype*, realtype*, realtype*);
  void FNVMHD_L1NORM(realtype*, realtype*);
  void FNVMHD_COMPARE(realtype*, realtype*, realtype*);
  void FNVMHD_INVTEST(realtype*, realtype*, int*);
  void FNVMHD_CONSTRMASK(realtype*, realtype*, realtype*, int*);
  void FNVMHD_MINQUOTIENT(realtype*, realtype*, realtype*);



/****************************************************************
 * PART III:                                                    *
 * implementation of N_Vector                                   *
 ****************************************************************/

/* The FNVMHD implementation of the N_Vector 'content' structure 
   contains the local data length of the vector, and a pointer 
   to an array of real constants */
struct _N_VectorContent_MHD {
  long int length;         /* local data length          */
  realtype *data;          /* local data array           */
  booleantype own_data;    /* flag for ownership of data */
};

typedef struct _N_VectorContent_MHD *N_VectorContent_MHD;


/****************************************************************
 *                                                              *
 * PART IV: Macros                                              *
 *    NV_CONTENT_MHD, NV_DATA_MHD, NV_LENGTH_MHD,               *
 *    NV_OWN_DATA_MHD, NV_Ith_MHD                               *
 *--------------------------------------------------------------*
 * In the descriptions below, the following user                *
 * declarations are assumed:                                    *
 *                                                              *
 * N_Vector     v;                                              *
 * long int     v_len, s_len, i;                                *
 *                                                              *
 * (1) NV_CONTENT_MHD                                           *
 *                                                              *
 *     This routine gives access to the contents of the         *
 *     N_Vector.                                                *
 *                                                              *
 *     The assignment v_cont = NV_CONTENT_MHD(v) sets           *
 *     v_cont to be a pointer to the MHD N_Vector               *
 *     content structure.                                       *
 *                                                              *
 * (2) NV_DATA_MHD, NV_LENGTH_MHD                               *
 *                                                              *
 *     These routines give individual access to the parts of    *
 *     the content of a MHD N_Vector.                           *
 *                                                              *
 *     The assignment v_data = NV_DATA_MHD(v) sets v_data to    *
 *     be a pointer to the first component of the local data    *
 *     for the vector v. The assignment                         *
 *     NV_DATA_MHD(v) = v_data  sets the component array of     *
 *     v to be v_data by storing the pointer v_data.            *  
 *                                                              *
 *     The assignment v_len = NV_LENGTH_MHD(v) sets             *
 *     v_len to be the length of the vector v.                  *
 *     The call NV_LENGTH_MHD(v) = len_v sets the               *
 *     length of v to be len_v.                                 *
 *                                                              *
 * (3) NV_OWN_DATA_MHD                                          *
 *                                                              *
 *     This routine gives access to the individual part of the  *
 *     N_Vector that determines ownership of the data.          *
 *                                                              *
 * (4) NV_Ith_MHD                                               *
 *                                                              *
 *     In the following description, the components of the      *
 *     local part of an N_Vector are numbered 0..n-1, where n   *
 *     is the local length of (the local part of) v.            *
 *                                                              *
 *     The assignment r = NV_Ith_MHD(v,i) sets r to be the      *
 *     value of the ith component of the local part of the      *
 *     vector v.  The assignment NV_Ith_MHD(v,i) = r sets       *
 *     the value of the ith local component of v to be r.       *
 *                                                              *
 * Notes..                                                      *
 *                                                              *
 * When looping over the components of an N_Vector v, it is     *
 * more efficient to first obtain the component array via       *
 * v_data = NV_DATA_MHD(v) and then access v_data[i] within     *
 * the loop than it is to use NV_Ith_MHD(v,i) within the        *
 * loop.                                                        *
 *                                                              *
 ****************************************************************/ 

#define NV_CONTENT_MHD(v) ( (N_VectorContent_MHD)(v->content) )

#define NV_DATA_MHD(v) ( NV_CONTENT_MHD(v)->data )

#define NV_LENGTH_MHD(v) ( NV_CONTENT_MHD(v)->length )

#define NV_Ith_MHD(v,i) ( NV_DATA_MHD(v)[i] )

#define NV_OWN_DATA_MHD(v) ( NV_CONTENT_MHD(v)->own_data )



/****************************************************************
 * PART V:                                                      *
 * Functions exported by nvector_mhd                            *
 *                                                              *
 * CONSTRUCTORS:                                                *
 *    N_VNew_MHD                                                *
 *    N_VNewEmpty_MHD                                           *
 *    N_VMake_MHD                                               *
 * EXTRA OPERATIONS:                                            *
 *    N_VPrint_MHD                                              *
 *--------------------------------------------------------------*/

/* N_VNew_MHD */
/* This function creates and allocates memory for an MHD vector */
N_Vector N_VNew_MHD(long int len);

/* N_VNewEmpty_MHD */ 
/* This function creates a new MHD vector with an empty (NULL) 
   data array (sets own_data to FALSE) */
N_Vector N_VNewEmpty_MHD(long int len);

/* N_VMake_MHD */
/* This function creates and allocates memory for a MHD vector 
   with a user-provided data array */
N_Vector N_VMake_MHD(long int len, realtype *data);

/* N_VPrint_MHD */
/* This function prints the content of a MHD vector to stdout */
void N_VPrint_MHD(N_Vector v);



/****************************************************************
 * PART VI:                                                     *
 * vector operations in nvector_mhd                             *
 *--------------------------------------------------------------*/

N_Vector    N_VCloneEmpty_MHD(N_Vector);
N_Vector    N_VClone_MHD(N_Vector);
void        N_VDestroy_MHD(N_Vector);
void        N_VSpace_MHD(N_Vector, long int *, long int *);
realtype   *N_VGetArrayPointer_MHD(N_Vector);
void        N_VSetArrayPointer_MHD(realtype *, N_Vector);
void        N_VLinearSum_MHD(realtype, N_Vector, realtype, 
			     N_Vector, N_Vector);
void        N_VConst_MHD(realtype, N_Vector);
void        N_VProd_MHD(N_Vector, N_Vector, N_Vector);
void        N_VDiv_MHD(N_Vector, N_Vector, N_Vector);
void        N_VScale_MHD(realtype, N_Vector, N_Vector);
void        N_VAbs_MHD(N_Vector, N_Vector);
void        N_VInv_MHD(N_Vector, N_Vector);
void        N_VAddConst_MHD(N_Vector, realtype, N_Vector);
realtype    N_VDotProd_MHD(N_Vector, N_Vector);
realtype    N_VMaxNorm_MHD(N_Vector);
realtype    N_VWrmsNorm_MHD(N_Vector, N_Vector);
realtype    N_VWrmsNormMask_MHD(N_Vector, N_Vector, N_Vector);
realtype    N_VMin_MHD(N_Vector);
realtype    N_VWL2Norm_MHD(N_Vector, N_Vector);
realtype    N_VL1Norm_MHD(N_Vector);
void        N_VCompare_MHD(realtype, N_Vector, N_Vector);
booleantype N_VInvTest_MHD(N_Vector, N_Vector);
booleantype N_VConstrMask_MHD(N_Vector, N_Vector, N_Vector);   
realtype    N_VMinQuotient_MHD(N_Vector, N_Vector);


#ifdef __cplusplus
}
#endif

#endif



/*********************** END OF FILE ***********************/
