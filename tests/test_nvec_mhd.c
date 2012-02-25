/*************************************************************************
 * File        : nvec_mhd_test.c                                         *
 * Programmers : Daniel R. Reynolds @ SMU                                *
 * Version of  : 2 October 2003                                          *
 *-----------------------------------------------------------------------*
 * This testing routine is to be called by the Fortran testing routine   *
 * fnvec_mhd_test.F with an array of initial data.  On that data, we     *
 * test the various vector kernel routines from the nvector_mhd_parallel *
 * and output the results.  From there, the Fortran calling routine will *
 * run the same tests and output the results as well.  These results     *
 * must be checked by hand to ensure that, first, the vector kernel      *
 * nvector_mhd_parallel is error-free, and second, that the associated   *
 * Fortran interface works as well.                                      *
 *                                                                       *
 * Note: this routine is called from the Fortran routine test_fnvec_mhd, *
 * which already makes certain that the NVector kernel has been compiled *
 * to use MPI.                                                           *
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "nvector_mhd.h"
#include "sundialsmath.h"
#include "test_nvec_mhd.h"
#include "mpi.h"


extern N_Vector F2C_KINSOL_vec;

/* C-testing routine */
void F_NVECMHDTEST( realtype *u, realtype *urp, int *mype, int *outproc, int *nx, 
		    int *ny, int *nz, int *xghost, int *yghost, int *zghost, 
		    int *ns, int *ier ) 
{
  /* Initialize local variables */
  N_Vector Uvec, Xvec, Yvec, Zvec;
  long int nfloat, nint;
  realtype a, b;
  FILE *fhandle;
  booleantype booltest;
  MPI_Comm comm;

  /* Make NVector out of input data */
  Uvec = N_VCloneEmpty_MHD(F2C_KINSOL_vec);
  N_VSetArrayPointer(u,Uvec);
  NV_RADIUS_MHD(Uvec) = *urp;
    
  /* get MPI communicator from MHD NVector */
  comm = MPI_COMM_WORLD;

  /* Make auxiliary test vectors */
  Xvec = N_VClone(Uvec);
  Yvec = N_VClone(Uvec);
  Zvec = N_VClone(Uvec);


  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 0                                         */
  /*   Print vector data to screen with (FortPrint) */
  if (*mype == *outproc) {
    printf("   C-Test 0: printing the input vector using FortPrint\n");
    printf("      also printing results to file C_test00\n");
    FortPrint(Uvec);
    fhandle = fopen("C_test00","w");  
    FilePrint( Uvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 1     */
  /*   N_VSpace */
  N_VSpace(Uvec,&nfloat,&nint);
  if (*mype == *outproc) {
    printf("   C-Test 1: the N_Vector requires %i reals and %i ints\n",
	   nfloat,nint);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 2     */
  /*   N_VScale: first copy U into X, then put 3*U into Y and output */
  a = 1.0; 
  b = 3.0;
  N_VScale( a, Uvec, Xvec );
  N_VScale( b, Xvec, Yvec );
  if (*mype == *outproc) {
    printf("   C-Test 2: printing results to file C_test02\n");
    fhandle = fopen("C_test02","w");  
    FilePrint( Yvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 3                                                      */
  /*   N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=0.5, b=-0.5 */
  a = 0.5; 
  b = -0.5;
  N_VLinearSum( a, Xvec, b, Yvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 3: printing results to file C_test03\n");
    fhandle = fopen("C_test03","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);

  
  /* Test 4                                                      */
  /*   N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=-1.0, b=1.0 */
  a = -1.0; 
  b = 1.0;
  N_VLinearSum( a, Xvec, b, Yvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 4: printing results to file C_test04\n");
    fhandle = fopen("C_test04","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 5                                                     */
  /*   N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=1.0, b=0.5 */
  a = 1.0; 
  b = 0.5;
  N_VLinearSum( a, Xvec, b, Yvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 5: printing results to file C_test05\n");
    fhandle = fopen("C_test05","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 6                                                             */
  /*   N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=5.0, b=1.0, Z=Yvec */
  a = 5.0; 
  b = 1.0;
  N_VLinearSum( a, Xvec, b, Yvec, Yvec );
  if (*mype == *outproc) {
    printf("   C-Test 6: printing results to file C_test06\n");
    fhandle = fopen("C_test06","w");  
    FilePrint( Yvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 7            */
  /*   N_VConst: a=8.2 */
  a = 8.2; 
  N_VConst( a, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 7: printing results to file C_test07\n");
    fhandle = fopen("C_test07","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 8                                  */
  /*   N_VProd: X=Xvec==Uvec, Y=Yvec==3*Uvec */
  b = 3.0;
  N_VScale( b, Uvec, Yvec );
  N_VProd( Xvec, Yvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 8: printing results to file C_test08\n");
    fhandle = fopen("C_test08","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 9                                */
  /*   N_VDiv: X=Xvec=Uvec, Y=Yvec==3*Uvec */
  N_VDiv( Xvec, Yvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 9: printing results to file C_test09\n");
    fhandle = fopen("C_test09","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 10                                       */
  /*   N_VAddConst: X=Xvec==Uvec, b=-20000, Z=Yvec */
  b = -20000.0;
  N_VAddConst( Xvec, b, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 10: printing results to file C_test10\n");
    fhandle = fopen("C_test10","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 11          */
  /*   N_VAbs: X=Yvec */
  N_VAbs( Zvec, Yvec );
  if (*mype == *outproc) {
    printf("   C-Test 11: printing results to file C_test11\n");
    fhandle = fopen("C_test11","w");  
    FilePrint( Yvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 12          */
  /*   N_VInv: X=Yvec */
  N_VInv( Yvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 12: printing results to file C_test12\n");
    fhandle = fopen("C_test12","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 13                                    */
  /*   N_VDotProd: X=Xvec==Uvec, Y=Yvec==3*Uvec */
  b = 3.0;
  N_VScale( b, Uvec, Yvec );
  a = N_VDotProd( Xvec, Yvec );
  if (*mype == *outproc)
    printf("   C-Test 13: DotProd(X,Y) = %g\n",a);

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 14                    */
  /*   N_VMaxNorm: X=Xvec==Uvec */
  a = N_VMaxNorm( Xvec );
  if (*mype == *outproc)
    printf("   C-Test 14: MaxNorm(X) = %g\n",a);

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 15                                  */
  /*   N_VWrmsNorm: X=Xvec==Uvec, W=Yvec==1.5 */
  b = 1.5;
  N_VConst( b, Yvec );
  a = N_VWrmsNorm( Xvec, Yvec );
  if (*mype == *outproc)
    printf("   C-Test 15: WrmsNorm(X,W) = %g\n",a);

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 16                                                  */
  /*   N_VWrmsNormMask: X=Xvec==Uvec, W=Yvec==1.5, ID=Zvec==1 */
  b = 1.0;
  N_VConst( b, Zvec );
  a = N_VWrmsNormMask( Xvec, Yvec, Zvec );
  if (*mype == *outproc)
    printf("   C-Test 16: WrmsNormMask(X,W,Id) = %g\n",a);

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 17                */
  /*   N_VMin: X=Xvec==Uvec */
  a = N_VMin( Xvec );
  if (*mype == *outproc)
    printf("   C-Test 17: Min(X) = %g\n",a);

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 18                                    */
  /*   N_VWL2Norm: X=Xvec==Uvec, W=Yvec==2*Uvec */
  b = 2.0;
  N_VScale( b, Uvec, Yvec );
  a = N_VWL2Norm( Xvec, Yvec );
  if (*mype == *outproc)
    printf("   C-Test 18: WL2Norm(X,W) = %g\n",a);


  /* Test 19                   */
  /*   N_VL1Norm: X=Xvec==Uvec */
  a = N_VL1Norm( Xvec );
  if (*mype == *outproc)
    printf("   C-Test 19: L1Norm(X) = %g\n",a);

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 20                                  */
  /*   N_VCompare: X=Xvec==Uvec, c=a==20000.0 */
  a = 2000.0;
  N_VCompare( a, Xvec, Zvec );
  if (*mype == *outproc) {
    printf("   C-Test 20: printing results to file C_test20\n");
    fhandle = fopen("C_test20","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 21                    */
  /*   N_VInvTest: X=Xvec==Uvec */
  booltest = N_VInvTest( Xvec, Zvec );
  printf("   C-Test 21: proc%i, InvTest(X) = %i\n",*mype, booltest);
  if (*mype == *outproc) {
    printf("     printing results to file C_test21\n");
    fhandle = fopen("C_test21","w");  
    FilePrint( Zvec, fhandle, *mype, *nx, *ny, *nz, *xghost, *yghost, 
	       *zghost, *ns );
    fclose(fhandle);
  }

  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Test 22          */
  /*   N_VMinQuotient: num=Xvec==Uvec, denom=Yvec==4d-5*Uvec */
  b = 0.00004;
  N_VScale( b, Uvec, Yvec );
  a = N_VMinQuotient( Xvec, Yvec );
  if (*mype == *outproc)
    printf("   C-Test 22: MinQuotient(X) = %g\n\n",a);
  
  /* pause test execution for everyone to catch up */
  MPI_Barrier(comm);


  /* Free homemade test vectors, Dispose of input vector, and return */
  N_VDestroy(Uvec);
  N_VDestroy(Xvec);
  N_VDestroy(Yvec);
  N_VDestroy(Zvec);  
  *ier = 0;
  
} 



/*************************************************************************
 * C function FortPrint acts as an interface between the C-testing       *
 * routine and the Fortran user-supplied subroutine print_fnvecmhd.      *
 *************************************************************************/

void FortPrint(N_Vector uu) 
{
  realtype *udata, urp;
  udata = N_VGetArrayPointer(uu);
  urp = NV_RADIUS_MHD(uu);
  PRINT_FNVECMHD(udata,&urp);
}



/*************************************************************************
 * C function FilePrint outputs the vector to the desired file           *
 *   Note: fhandle must be a file handle to a previously opened file     *
 *************************************************************************/

void FilePrint( N_Vector uu, FILE *fhandle, int mype, int Nx, int Ny, 
	       int Nz, int NGx, int NGy, int NGz, int Ns ) 
{
  /* Internal variables */
  realtype *udata, urp;
  int i, j, k, l, Sblock, Zblock, Yblock, loc;

  /* Get local access to the vector data */
  udata = N_VGetArrayPointer(uu);
  urp = NV_RADIUS_MHD(uu);

  fprintf(fhandle,"  %i %i %i %i %i %g\n",mype,0,0,0,0,urp);
  for (i=NGx; i<Nx+NGx; i++) {
    for (j=NGy; j<Ny+NGy; j++) { 
      Yblock = j * (Nx + 2*NGx);
      for (k=NGz; k<Nz+NGz; k++) {
	Zblock = k * (Nx + 2*NGx) * (Ny + 2*NGy);
	for (l=0; l<Ns; l++) {
	  Sblock = l * (Nx + 2*NGx) * (Ny + 2*NGy) * (Nz + 2*NGz);
	  loc = Sblock + Zblock + Yblock + i;
	  fprintf(fhandle,"  %i %i %i %i %i %g\n",mype,i-NGx+1,j-NGy+1,
		  k-NGz+1,l+1,udata[loc]);
	}
      }
    }
  }
  
}


/* end of testing routine ************************************************/

