/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Daniel R. Reynolds
 * UC San Diego, Mathematics
 * -----------------------------------------------------------------
 * Description: 
 *   This example sets up (and solves) a set of 3 matrices, each of 
 *   which occur in the viscous portion of a resistive MHD system.  
 *   This is strictly a small, single-processor test code, through 
 *   which we may test our interface to the HYPRE solver libraries.
 * -------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "viscous_prec.h"

/*--------------------------------------------------------------------
 * MAIN PROGRAM
 *------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  void *Amem, *Bmem, *Cmem;
  double *uu, *vv, *ww;
  double dtfac, Mu, Lu, Eta, Re, Gamma, Kappa, Pr, RGas;
  double delta;
  int local_N, i, j, k, l, indx, ier;
  long int Nx, Ny, Nz, Ns, NGx, NGy, NGz;
  double xL, xR, yL, yR, zL, zR, dx, dy, dz;
  MPI_Comm comm;
  int my_pe, npes, xprocs, yprocs, zprocs;


  /* Initialize matrix/solver memory to NULL */
  Amem = NULL;
  Bmem = NULL;
  Cmem = NULL;

  /* Set grid information */
  xprocs = 2; /* number of x mesh points per subgrid */
  yprocs = 1; /* number of y mesh points per subgrid */
  zprocs = 1; /* number of z mesh points per subgrid */
  Nx = 16;    /* number of x mesh points per subgrid */
  Ny = 32;    /* number of y mesh points per subgrid */
  Nz = 32;    /* number of z mesh points per subgrid */
  NGx = 4;    /* number of x mesh ghost cells */
  NGy = 4;    /* number of y mesh ghost cells */
  NGz = 4;    /* number of z mesh ghost cells */
  Ns = 8;     /* number of MHD species */
  xL = 0.0;   /* x-lower boundary */
  xR = 2.0;   /* x-upper boundary */
  yL = -1.0;  /* y-lower boundary */
  yR = 1.0;   /* y-upper boundary */
  zL = -1.0;  /* z-lower boundary */
  zR = 1.0;   /* z-upper boundary */

  /* Initialize MPI environment and Cartesian communicator */
  int MPIdims[3] = {xprocs, yprocs, zprocs};
  int periods[3] = {0, 0, 0};
  int three = 3;
  int zero = 0;
  MPI_Comm comm3d;
  int xbc=0, ybc=0, zbc=0;
#ifdef XPERIODIC
  xbc = 1;
  periods[0] = 1;
#elif XREFLECT
  xbc = 2;
#endif
#ifdef YPERIODIC
  ybc = 1;
  periods[1] = 1;
#elif YREFLECT
  ybc = 2;
#endif
#ifdef ZPERIODIC
  zbc = 1;
  periods[2] = 1;
#elif ZREFLECT
  zbc = 2;
#endif
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  if (npes != xprocs*yprocs*zprocs) {
    printf("ERROR: %i procs provided, run requires %i\n",
	   npes,xprocs*yprocs*zprocs);
    return -1;
  }
  MPI_Comm_rank(comm, &my_pe);
  printf("hypre_test: starting proc %i of %i\n",my_pe,npes);
  MPI_Cart_create(MPI_COMM_WORLD, three, MPIdims, periods, zero, &comm3d);
  int MPdims[3], MPperiods[3], MPcoords[3];
  MPI_Cart_get(comm3d, three, MPdims, MPperiods, MPcoords);
  MPcoords[0] += 1; MPcoords[1] += 1; MPcoords[2] += 1;
  printf("  p%i: iprocs = (%i,%i,%i), periods = (%i,%i,%i), coords = (%i,%i,%i)\n", my_pe, MPdims[0], MPdims[1], MPdims[2], MPperiods[0], MPperiods[1], MPperiods[2], MPcoords[0], MPcoords[1], MPcoords[2]);
  long int NxGlobal = Nx*MPdims[0];
  long int NyGlobal = Ny*MPdims[1];
  long int NzGlobal = Nz*MPdims[2];
  long int iXL = Nx*(MPcoords[0]-1) + 1;
  long int iXR = iXL + Nx - 1;
  long int iYL = Ny*(MPcoords[1]-1) + 1;
  long int iYR = iYL + Ny - 1;
  long int iZL = Nz*(MPcoords[2]-1) + 1;
  long int iZR = iZL + Nz - 1;
  printf("  p%i: globdims = (%li,%li,%li), extents = (%li:%li,%li:%li,%li:%li)\n", 
	 my_pe, NxGlobal, NyGlobal, NzGlobal, iXL, iXR, iYL, iYR, iZL, iZR); 
  int NBlt, NBrt, NBbt, NBtp, NBbk, NBft;
  int source;
  MPI_Cart_shift(comm3d, 0, -1, &source, &NBlt);
  MPI_Cart_shift(comm3d, 0,  1, &source, &NBrt);
  MPI_Cart_shift(comm3d, 1, -1, &source, &NBbt);
  MPI_Cart_shift(comm3d, 1,  1, &source, &NBtp);
  MPI_Cart_shift(comm3d, 2, -1, &source, &NBbk);
  MPI_Cart_shift(comm3d, 2,  1, &source, &NBft);
  printf("  p%i: neighbors = (%i:%i,%i:%i,%i:%i)\n", my_pe, NBlt, NBrt, NBbt, NBtp, NBbk, NBft);
  int outproc = (my_pe == 0);
  fflush(NULL);
  MPI_Barrier(comm3d);

  /* initialize domain information */
  dx = (xR-xL)/NxGlobal;
  dy = (yR-yL)/NyGlobal;
  dz = (zR-zL)/NzGlobal;

  /* Allocate vector memory, and set initial values */ 
  local_N = Ns * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
  uu = (double *) malloc(local_N * sizeof(double));
  vv = (double *) malloc(local_N * sizeof(double));
  ww = (double *) malloc(local_N * sizeof(double));
  if ((uu==NULL) || (vv==NULL) || (ww==NULL)) {
    if (outproc) printf("hypre_test error: failed to allocate array data\n");
    MPI_Abort(comm, 1);
    return(1);
  }
  for (i=0; i<local_N; i++) {
    uu[i] = 1.0;
    vv[i] = 0.0;
    ww[i] = 1.0;
  }
  
  /* set MHD parameters */
  dtfac = 1.0e-1;    /* time step factor */
  Mu    = 1.0e-1;    /* fluid viscosity */
  Lu    = 1.0;       /* Lundquist number */
  Eta   = 1.0e-0;    /* magnetic resistivity */
  Re    = 1.0;       /* Reynolds number */
  Gamma = 1.667;     /* specific heat ratio */
  Kappa = 1.0e-2;    /* thermal conductivity */
  Pr    = 1.0;       /* Prandtl number */
  RGas  = 1.0;       /* gas constant */

  /* set Hypre solver tolerance */
  delta = 1.0e-5;


  /* Allocate matrix/solver memory */
  if (outproc) printf("\n  hypre_test: allocating Du matrix/solver memory\n");
  Amem = VPrecDuAlloc(Nx, Ny, Nz, dx, dy, dz, Ns, NGx, NGy, NGz, 
		      MPdims[0], MPdims[1], MPdims[2],
		      MPcoords[0], MPcoords[1], MPcoords[2],
		      NBlt, NBrt, NBtp, NBbt, NBft, NBbk, xbc, ybc, zbc);
  fflush(NULL);
  MPI_Barrier(comm3d);

  if (outproc) printf("\n  hypre_test: allocating Db matrix/solver memory\n");
  Bmem = VPrecDbAlloc(Nx, Ny, Nz, dx, dy, dz, Ns, NGx, NGy, NGz, 
		      MPdims[0], MPdims[1], MPdims[2],
		      MPcoords[0], MPcoords[1], MPcoords[2],
		      NBlt, NBrt, NBtp, NBbt, NBft, NBbk, xbc, ybc, zbc);
  fflush(NULL);
  MPI_Barrier(comm3d);

  if (outproc) printf("\n  hypre_test: allocating De matrix/solver memory\n");
  Cmem = VPrecDeAlloc(Nx, Ny, Nz, dx, dy, dz, Ns, NGx, NGy, NGz, 
		      MPdims[0], MPdims[1], MPdims[2],
		      MPcoords[0], MPcoords[1], MPcoords[2],
		      NBlt, NBrt, NBtp, NBbt, NBft, NBbk, xbc, ybc, zbc);
  fflush(NULL);
  MPI_Barrier(comm3d);


  /* loop over "time" steps */
  int TimeSteps = 5;
  int NewtonIts = 5;
  int tstep, newt, ix, iy, iz;
  double alternating = 1.0;
  for (tstep=1; tstep<=TimeSteps; tstep++) {

    if (outproc) printf("\n\n Beginning time step %i/%i:\n\n",tstep,TimeSteps);

    /* update state vector uu */
    alternating *= -1;
    for (l=0; l<Ns; l++)
      for (k=0; k<Nz+2*NGz; k++) 
	for (j=0; j<Ny+2*NGy; j++) 
	  for (i=0; i<Nx+2*NGx; i++) {
	    iz = iZL+k-NGz;
	    iy = iYL+j-NGy;
	    ix = iXL+i-NGx;
	    uu[((l*(Nz+2*NGz)+k)*(Ny+2*NGy)+j)*(Nx+2*NGx)+i] = 
	      2.0+alternating*tstep/TimeSteps*ix/NxGlobal*iy/NyGlobal;
	  }
    

    /* Set up matrices */
    if (outproc) printf("\n  hypre_test: setting up Du matrix\n");
    ier = VPrecDuSetup(uu, dtfac, Mu, Re, vv, ww, Amem);
    if (ier != 0) {
      if (outproc) printf("  hypre_test error: Du matrix setup returned flag %i\n",ier);
      MPI_Abort(comm, 1);
      return(ier);
    }
    fflush(NULL);
    MPI_Barrier(comm3d);
    
    if (outproc) printf("\n  hypre_test: setting up Db matrix\n");
    ier = VPrecDbSetup(uu, dtfac, Lu, Eta, vv, ww, Bmem);
    if (ier != 0) {
      if (outproc) printf("  hypre_test error: Db matrix setup returned flag %i\n",ier);
      MPI_Abort(comm, 1);
      return(ier);
    }
    fflush(NULL);
    MPI_Barrier(comm3d);
    
    if (outproc) printf("\n  hypre_test: setting up De matrix\n");
    ier = VPrecDeSetup(uu, dtfac, Gamma, Kappa, Re, Pr, RGas, vv, ww, Cmem);
    if (ier != 0) {
      if (outproc) printf("  hypre_test error: De matrix setup returned flag %i\n",ier);
      MPI_Abort(comm, 1);
      return(ier);
    }
    fflush(NULL);
    MPI_Barrier(comm3d);
    

    /* take 'Newton' steps (solve multiple systems with same matrix) */
    for (newt=1; newt<=NewtonIts; newt++) {

      if (outproc) printf("\n  Newton Iteration %i/%i:\n",newt,NewtonIts);
      alternating *= -1;

      /* Solve Du system */
      if (outproc) printf("\n  hypre_test: solving Du matrix system \n");
      for (l=0; l<Ns; l++) {
	for (k=0; k<Nz; k++) { 
	  for (j=0; j<Ny; j++) {
	    for (i=0; i<Nx; i++) {
	      indx = l*(Nz+2*NGz)*(Ny+2*NGy)*(Nx+2*NGx) 
		+ (k+NGz)*(Ny+2*NGy)*(Nx+2*NGx)
		+ (j+NGy)*(Nx+2*NGx)
		+ i + NGx;
	      ww[indx] = 6.4*i*(Nx-1-i)/Nx/Nx*j*(Ny-1-j)/Ny/Ny*k*(Nz-1-k)/Nz/Nz*alternating*tstep*newt;
	    }
	  }
	}
      }
      ier = VPrecDuSolve(vv, ww, uu, delta, Amem);
      if (ier != 0) {
	if (outproc) printf("  hypre_test error: Du matrix solve returned flag %i\n",ier);
	VPrecDuFree(Amem);
	VPrecDbFree(Bmem);
	VPrecDeFree(Cmem);
	free(uu);
	free(vv);
	free(ww);
	MPI_Finalize();
	return(ier);
      }
      fflush(NULL);
      MPI_Barrier(comm3d);
    
      /* Solve Db system */
      if (outproc) printf("\n  hypre_test: solving Db matrix system \n");
      for (l=0; l<Ns; l++) {
	for (k=0; k<Nz; k++) { 
	  for (j=0; j<Ny; j++) {
	    for (i=0; i<Nx; i++) {
	      indx = l*(Nz+2*NGz)*(Ny+2*NGy)*(Nx+2*NGx) 
		+ (k+NGz)*(Ny+2*NGy)*(Nx+2*NGx)
		+ (j+NGy)*(Nx+2*NGx)
		+ i + NGx;
	      ww[indx] = 6.4*i*(Nx-1-i)/Nx/Nx*j*(Ny-1-j)/Ny/Ny*k*(Nz-1-k)/Nz/Nz*alternating*tstep*newt;
	    }
	  }
	}
      }
      ier = VPrecDbSolve(vv, ww, uu, delta, Bmem);
      if (ier != 0) {
	if (outproc) printf("  hypre_test error: Db matrix solve returned flag %i\n",ier);
	VPrecDuFree(Amem);
	VPrecDbFree(Bmem);
	VPrecDeFree(Cmem);
	free(uu);
	free(vv);
	free(ww);
	MPI_Finalize();
	return(ier);
      }
      fflush(NULL);
      MPI_Barrier(comm3d);
    
      /* Solve De system */
      if (outproc) printf("\n  hypre_test: solving De matrix system \n");
      for (l=0; l<Ns; l++) {
	for (k=0; k<Nz; k++) { 
	  for (j=0; j<Ny; j++) {
	    for (i=0; i<Nx; i++) {
	      indx = l*(Nz+2*NGz)*(Ny+2*NGy)*(Nx+2*NGx) 
		+ (k+NGz)*(Ny+2*NGy)*(Nx+2*NGx)
		+ (j+NGy)*(Nx+2*NGx)
		+ i + NGx;
	      ww[indx] = 6.4*i*(Nx-1-i)/Nx/Nx*j*(Ny-1-j)/Ny/Ny*k*(Nz-1-k)/Nz/Nz*alternating*tstep*newt;
	    }
	  }
	}
      }
      ier = VPrecDeSolve(vv, ww, uu, delta, Cmem);
      if (ier != 0) {
	if (outproc) printf("  hypre_test error: De matrix solve returned flag %i\n",ier);
	VPrecDuFree(Amem);
	VPrecDbFree(Bmem);
	VPrecDeFree(Cmem);
	free(uu);
	free(vv);
	free(ww);
	MPI_Finalize();
	return(ier);
      }
      fflush(NULL);
      MPI_Barrier(comm3d);


    }   /* end "Newton" steps */

  }  /* end "time" steps */


  /* Freeing matrix/solver memory */
  if (outproc) printf("\n  hypre_test: freeing Du matrix/solver memory\n");
  VPrecDuFree(Amem);
  fflush(NULL);
  MPI_Barrier(comm3d);

  if (outproc) printf("\n  hypre_test: freeing Db matrix/solver memory\n");
  VPrecDbFree(Bmem);
  fflush(NULL);
  MPI_Barrier(comm3d);

  if (outproc) printf("\n  hypre_test: freeing De matrix/solver memory\n");
  VPrecDeFree(Cmem);
  fflush(NULL);
  MPI_Barrier(comm3d);


  /* Freeing vector memory */
  free(uu);
  free(vv);
  free(ww);


  /* finalize MPI */
  MPI_Finalize();

  return(0);
}

