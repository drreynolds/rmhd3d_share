/*******************************************************************
 * Daniel R. Reynolds                                              *
 * UC San Diego, Mathematics                                       *
 *-----------------------------------------------------------------*
 * This is the implementation file for an implementation           *
 * of the NVECTOR package, specifically suited to interface with a *
 * Fortran data structure used by a collaborator's MHD code from   *
 * PPPL, where local spatial data includes ghost cell values.      *
 *                                                                 *
 * This vector kernel can be used in either parallel with MPI or   *
 * in serial, depending on the prepocessor variable 'PARALLEL'.    *
 *                                                                 *
 * It contains the N_Vector kernels listed in nvector_mhd.h,       *
 * as well as implementations of the Fortran interface to the      *
 * create and destroy routines N_VNew_MHD and N_VDestroy_MHD.      *
 *                                                                 *
 * It uses MPI for message-passing in the case that a parallel     *
 * environment is desired.                                         *
 *                                                                 *
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "nvector_mhd.h"
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)



/********************* Fortran Interface Functions ***************/

/* define global N_Vector variables */
N_Vector F2C_CVODE_vec;
N_Vector F2C_KINSOL_vec;


/* FNVMHD_INIT is the Fortran interface routine to initialize the 
   template NVector */
void FNVMHD_INIT(int *code, long int *length, int *ier)
{
  /* initialize return value to success */
  *ier = 0;
  
  /* Create solver-specific global NVector */
  /*   code = FCMIX_CVODE  (1) -> use with CVODE */
  /*   code = FCMIX_KINSOL (5) -> use with KINSOL */
  /*   [others not implemented] */
  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = NULL;
    F2C_CVODE_vec = N_VNewEmpty_MHD(*length);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = NULL;
    F2C_KINSOL_vec = N_VNewEmpty_MHD(*length);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}



/********************* Exported Functions ************************/


/* Function to create a new, empty MHD vector */
N_Vector N_VNewEmpty_MHD(long int length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_MHD content;

  /* Initialize vector pointers to NULL */
  v = NULL;
  ops = NULL;
  content = NULL;
  
  /* Create vector */
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) {free(v); return(NULL); }

  /* Attach custom vector routines to N_Vector_Ops structure */
  ops->nvclone           = N_VClone_MHD;
  ops->nvcloneempty      = N_VCloneEmpty_MHD;
  ops->nvdestroy         = N_VDestroy_MHD;
  ops->nvspace           = N_VSpace_MHD;
  ops->nvgetarraypointer = N_VGetArrayPointer_MHD;
  ops->nvsetarraypointer = N_VSetArrayPointer_MHD;
  ops->nvlinearsum       = N_VLinearSum_MHD;
  ops->nvconst           = N_VConst_MHD;
  ops->nvprod            = N_VProd_MHD;
  ops->nvdiv             = N_VDiv_MHD;
  ops->nvscale           = N_VScale_MHD;
  ops->nvabs             = N_VAbs_MHD;
  ops->nvinv             = N_VInv_MHD;
  ops->nvaddconst        = N_VAddConst_MHD;
  ops->nvdotprod         = N_VDotProd_MHD;
  ops->nvmaxnorm         = N_VMaxNorm_MHD;
  ops->nvwrmsnorm        = N_VWrmsNorm_MHD;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_MHD;
  ops->nvmin             = N_VMin_MHD;
  ops->nvwl2norm         = N_VWL2Norm_MHD;
  ops->nvl1norm          = N_VL1Norm_MHD;
  ops->nvcompare         = N_VCompare_MHD;
  ops->nvinvtest         = N_VInvTest_MHD;
  ops->nvconstrmask      = N_VConstrMask_MHD;
  ops->nvminquotient     = N_VMinQuotient_MHD;

  /* Create content */
  content = 
    (N_VectorContent_MHD) malloc(sizeof(struct _N_VectorContent_MHD));
  if (content == NULL) {free(ops); free(v); return(NULL);}

  /* Attach lengths and communicator to N_Vector_Content_MHD structure */
  content->length   = length;
  content->data     = NULL;
  content->own_data = FALSE;

  /* Attach content and ops to generic N_Vector */
  v->content = content;
  v->ops     = ops;

  return(v);
}




/* N_VNew_MHD (or nvnew) returns a new N_Vector.  If there is not 
   enough memory for a new N_Vector, then nvnew returns NULL */
N_Vector N_VNew_MHD(long int length)
{
  N_Vector v;
  realtype *data;

  v = NULL;
  data = NULL;

  /* Create the new N_Vector */
  v = N_VNewEmpty_MHD(length);
  if (v == NULL) return(NULL);

  /* Create data */
  if ( length > 0 ) {

    /* Allocate data memory */
    data = (realtype *) calloc(length,sizeof(realtype));
    if(data == NULL) { N_VDestroy_MHD(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_MHD(v) = TRUE;
    NV_DATA_MHD(v)     = data;
  }

  return(v);
}




/* N_VMake_MHD (or nvmake) creates a MHD N_Vector with component 
   array allocated by the user. */
N_Vector N_VMake_MHD(long int length, realtype *v_data)
{
  N_Vector v;

  /* initialize v to NULL */
  v = NULL;

  /* Create vector */
  v = N_VNewEmpty_MHD(length);
  if (v == NULL) return(NULL);

  /* Attach data if it has nonzero size */
  if ( length > 0 ) {
    NV_OWN_DATA_MHD(v) = FALSE;
    NV_DATA_MHD(v)     = v_data;
  }

  return(v);
}



/* N_VPrint_MHD prints the N_Vector v to stdout.  This routine is 
   provided to aid in debugging code using this vector package. */
void N_VPrint_MHD(N_Vector v)
{
  realtype *vd;

  /* initialize vd to NULL */
  vd = NULL;

  /* extract data array from vector */
  vd = NV_DATA_MHD(v);

  /* call Fortran subroutine to output vector */
  FNVMHD_PRINT(vd);
}




/***************************************************************************/


/* BEGIN implementation of vector operations */


/* N_VCloneEmpty_MHD returns a new N_Vector of the same form as the 
   input N_Vector, but with empty data container */ 
N_Vector N_VCloneEmpty_MHD(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_MHD content;

  /* initialize pointers to NULL */
  v = NULL;
  ops = NULL;
  content = NULL;

  /* Check that w has been created */
  if (w == NULL) return(NULL);

  /* Create vector */
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  /* Attach operations */
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;
  ops->nvprod            = w->ops->nvprod;
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale;
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */
  content = 
    (N_VectorContent_MHD) malloc(sizeof(struct _N_VectorContent_MHD));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and data to content structure */
  content->length   = NV_LENGTH_MHD(w);
  content->own_data = FALSE;
  content->data     = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}




/* N_VClone_MHD returns a new N_Vector of the same form as the 
   input N_Vector. */ 
N_Vector N_VClone_MHD(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int length;

  /* initialize pointers to NULL */
  v = NULL;
  data = NULL;

  /* Create vector */
  v = N_VCloneEmpty_MHD(w);
  if (v == NULL) return(NULL);

  /* Get local data size, create local data */
  length = NV_LENGTH_MHD(w);
  if ( length > 0 ) {

    /* Allocate memory */
    data = (realtype *) calloc(length,sizeof(realtype));
    if(data == NULL) { N_VDestroy_MHD(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_MHD(v) = TRUE;
    NV_DATA_MHD(v)     = data;
  }

  return(v);
}




/* N_VDestroy_MHD frees the data storage for an [empty] N_Vector */ 
void N_VDestroy_MHD(N_Vector v)
{
  if ( (NV_OWN_DATA_MHD(v) == TRUE) && (NV_DATA_MHD(v) != NULL) ) {
    free(NV_DATA_MHD(v));
    NV_DATA_MHD(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}




/* N_VSpace_MHD returns the space requirements for one N_Vector.  The 
   amount of realtype data is given in lrw, and long int data in liw.  
   Note: this includes ghost cell data storage as well */
void N_VSpace_MHD(N_Vector v, long int *lrw, long int *liw)
{
  *lrw = NV_LENGTH_MHD(v);
  *liw = 2;

  return;
}




/* N_VGetArrayPointer_MHD (or nvgetarraypointer) extracts the 
   data component array from the N_Vector v.  */
realtype *N_VGetArrayPointer_MHD(N_Vector v)
{
  return((realtype *) NV_DATA_MHD(v));
}



/* N_VSetArrayPointer_MHD or (nvsetarraypointer) attaches the 
   data component array v_data to the N_Vector v.  */
void N_VSetArrayPointer_MHD(realtype *v_data, N_Vector v)
{
  if (NV_LENGTH_MHD(v) > 0) NV_DATA_MHD(v) = v_data;

  return;
}




/* N_VLinearSum_MHD (or nvlinearsum) calculates z = a*x + b*y */
void N_VLinearSum_MHD(realtype a, N_Vector x, realtype b, 
		      N_Vector y, N_Vector z)
{
  /* extract data array handles from vectors */
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;
  xd = NV_DATA_MHD(x);
  yd = NV_DATA_MHD(y);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do operation */
  FNVMHD_LINSUM(&a, xd, &b, yd, zd);

  return;
}





/* N_VConst_MHD (or nvconst) calculates z[i] = c for all i */
void N_VConst_MHD(realtype c, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *zd;
  zd = NULL;
  zd = NV_DATA_MHD(z);

  /* call fortran routine to set z to c */
  FNVMHD_CONST(&c, zd);

  return;
}





/* N_VProd_MHD (or nvprod) calculates z[i] = x[i]*y[i] */
void N_VProd_MHD(N_Vector x, N_Vector y, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;
  xd = NV_DATA_MHD(x);
  yd = NV_DATA_MHD(y);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do product */
  FNVMHD_PROD(xd, yd, zd);

  return;
}





/* N_VDiv_MHD (or nvdiv) calculates z[i] = x[i]/y[i] */
void N_VDiv_MHD(N_Vector x, N_Vector y, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *yd, *zd;
  xd = yd = zd = NULL;
  xd = NV_DATA_MHD(x);
  yd = NV_DATA_MHD(y);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do division */
  FNVMHD_DIV(xd, yd, zd);

  return;
}





/* N_VScale_MHD (or nvscale) calculates z = c*x */
void N_VScale_MHD(realtype c, N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_MHD(x);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do operation */
  FNVMHD_SCALE(&c, xd, zd);

  return;
}





/* N_VAbs_MHD or (nvabs) calculates z[i] = |x[i]| */
void N_VAbs_MHD(N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_MHD(x);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do operation */
  FNVMHD_ABS(xd, zd);

  return;
}






/* N_VInv_MHD (or nvinv) calculates z[i] = 1/x[i].  
   Note: it does not check for division by 0.  It should be called only 
   with an N_Vector x which is guaranteed to have all non-zero components. */
void N_VInv_MHD(N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_MHD(x);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do operation */
  FNVMHD_INV(xd, zd);
  
  return;
}






/* N_VAddConst_MHD (or nvaddconst) calculates z[i] = x[i] + b */
void N_VAddConst_MHD(N_Vector x, realtype b, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_MHD(x);
  zd = NV_DATA_MHD(z);
  
  /* call fortran routine to do operation */
  FNVMHD_ADDCONST(&b, xd, zd);

  return;
}






/* N_VDotProd_MHD (or nvdotprod) returns the value of the 
   ordinary dot product of x and y, i.e. sum (i=0 to N-1) {x[i] * y[i]} */
realtype N_VDotProd_MHD(N_Vector x, N_Vector y)
{
  /* extract data array handles from N_Vectors */
  realtype sum, *xd, *yd;
  sum = ZERO;
  xd = yd = NULL;
  xd = NV_DATA_MHD(x);
  yd = NV_DATA_MHD(y);

  /* call fortran routine to do operation */
  FNVMHD_DOTPROD(xd, yd, &sum);
  return(sum);
}





/* N_VMaxNorm_MHD (or nvmaxnorm) returns the maximum norm of x, 
   i.e. max(i=1 to N-1) |x[i]|  */
realtype N_VMaxNorm_MHD(N_Vector x)
{ 
  /* extract data array handles from N_Vectors */
  realtype maxval, *xd;
  maxval = ZERO;
  xd = NULL;
  xd = NV_DATA_MHD(x);

  /* call fortran routine to do operation */
  FNVMHD_MAXNORM(xd, &maxval);
  return(maxval);
}






/* N_VWrmsNorm_MHD (or nvwrmsnorm) returns the weighted root 
   mean square norm of x with weight factor w, 
   i.e. sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N] */
realtype N_VWrmsNorm_MHD(N_Vector x, N_Vector w)
{
  /* extract data array handles from N_Vectors */
  realtype wrmsval, *xd, *wd;
  wrmsval = ZERO;
  xd = wd = NULL;
  xd = NV_DATA_MHD(x);
  wd = NV_DATA_MHD(w);

  /* call fortran routine to do operation */
  FNVMHD_WRMSNORM(xd, wd, &wrmsval);
  return(wrmsval);
}





/* N_VWrmsNormMask_MHD or (nvwrmsnormmask) returns wrms norm over 
   indices indicated by id */
realtype N_VWrmsNormMask_MHD(N_Vector x, N_Vector w, N_Vector id)
{
  /* extract data array handles from N_Vectors */
  realtype wrmsval, *xd, *wd, *idd;
  wrmsval = ZERO;
  xd = wd = idd = NULL;
  xd  = NV_DATA_MHD(x);
  wd  = NV_DATA_MHD(w);
  idd = NV_DATA_MHD(id);

  /* call fortran routine to do operation */
  FNVMHD_WRMSNORMMASK(xd, wd, idd, &wrmsval);
  return(wrmsval);
}




/* N_VMin_MHD (or nvmin) returns the smallest element of x */
realtype N_VMin_MHD(N_Vector x)
{
  /* extract data array handles from N_Vectors */
  realtype minval, *xd;
  xd = NULL;
  xd = NV_DATA_MHD(x);

  /* call fortran routine to do operation */
  FNVMHD_MIN(xd, &minval);
  return(minval);
}





/* N_VWL2Norm_MHD (or nvwl2norm) returns the weighted 
   Euclidean L2 norm of x with weight factor w, 
   i.e. sqrt [(sum (i=0 to N-1) {(x[i]*w[i])^2}) ] */
realtype N_VWL2Norm_MHD(N_Vector x, N_Vector w)
{
  /* extract data array handles from N_Vectors */
  realtype wl2val, *xd, *wd;
  wl2val = ZERO;
  xd = wd = NULL;
  xd = NV_DATA_MHD(x);
  wd = NV_DATA_MHD(w);

  /* call fortran routine to do operation */
  FNVMHD_WL2NORM(xd, wd, &wl2val);
  return(wl2val);
}





/* N_VL1Norm_MHD (or nvl1norm) returns the L1 norm of x, 
   i.e. sum (i=0 to N-1) {|x[i]|} */
realtype N_VL1Norm_MHD(N_Vector x)
{
  /* extract data array handles from N_Vectors */
  realtype l1norm, *xd;
  l1norm = ZERO;
  xd = NULL;
  xd = NV_DATA_MHD(x);

  /* call fortran routine to do operation */
  FNVMHD_L1NORM(xd, &l1norm);
  return(l1norm);
}





/* N_VCompare_MHD (or nvcompare) calculates 
   z[i] = 1 if |x[i]| > c, z[i] = 0 otherwise */
void N_VCompare_MHD(realtype c, N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  realtype *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_MHD(x);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do operation */
  FNVMHD_COMPARE(&c, xd, zd);

  return;
}





/* N_VInvTest_MHD (or nvinvtest) computes z[i] = 1/x[i] 
   with a test for x[i] == 0 before inverting x[i].  This routine 
   returns TRUE if all components of x are nonzero (successful 
   inversion) and returns FALSE otherwise. */
booleantype N_VInvTest_MHD(N_Vector x, N_Vector z)
{
  /* extract data array handles from N_Vectors */
  int testval;
  realtype *xd, *zd;
  xd = zd = NULL;
  xd = NV_DATA_MHD(x);
  zd = NV_DATA_MHD(z);

  /* call fortran routine to do operation */
  FNVMHD_INVTEST(xd, zd, &testval);

  if (testval == ZERO)  return(TRUE);
  else  return(FALSE);
}





/* N_VConstrMask_MHD (or nvconstrmask) returns a boolean FALSE 
   if any element fails the constraint test, and TRUE if all passed.  The 
   constraint test is as follows: 
         if c[i] =  2.0, then x[i] must be >  0.0
         if c[i] =  1.0, then x[i] must be >= 0.0
         if c[i] = -1.0, then x[i] must be <= 0.0
         if c[i] = -2.0, then x[i] must be <  0.0
   It also sets a mask vector m, with elements equal to 1.0 where the 
   corresponding constraint test failed, and equal to 0.0 where the 
   constraint test passed.  This routine is specialized in that it is 
   used only for constraint checking. */
booleantype N_VConstrMask_MHD(N_Vector c, N_Vector x, N_Vector m)
{
  /* extract data array handles from N_Vectors */
  int testval;
  realtype *cd, *xd, *md;
  cd = xd = md = NULL;
  cd = NV_DATA_MHD(c);
  xd = NV_DATA_MHD(x);
  md = NV_DATA_MHD(m);

  /* call fortran routine to do operation */
  FNVMHD_CONSTRMASK(cd, xd, md, &testval);

  if (testval == ZERO)  return(TRUE);
  else  return(FALSE);
}





/* N_VMinQuotient_MHD (or nvminquotient) returns 
   min(num[i]/denom[i]) over all i such that denom[i] != 0. */
realtype N_VMinQuotient_MHD(N_Vector num, N_Vector denom)
{
  /* extract data array handles from N_Vectors */
  realtype *nd, *dd, minquot;
  nd = dd = NULL;
  nd = NV_DATA_MHD(num);
  dd = NV_DATA_MHD(denom);
 
  /* call fortran routine to do operation */
  FNVMHD_MINQUOTIENT(nd, dd, &minquot);
  return(minquot);
}

 
/*********************** END OF FILE ***********************/
