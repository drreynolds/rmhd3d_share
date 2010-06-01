/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Daniel R. Reynolds
 * UC San Diego, Mathematics
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "viscous_prec.h"     /* precond. routine prototypes      */


static void *VPDuData;


/********************************************************************/
/* Fortran callable interface routines                              */
/********************************************************************/


/* Du preconditioner dataspace allocation wrapper routine */
void VISCPREC_DU_INIT(long int *Nx, long int *Ny, long int *Nz, 
		      double *dx, double *dy, double *dz,
		      long int *Ns, long int *NGx, long int *NGy, 
		      long int *NGz, int *NPx, int *NPy, 
		      int *NPz, int *iPx, int *iPy, int *iPz, 
		      int *NBlt, int *NBrt, int *NBtp, int *NBbt, 
		      int *NBft, int *NBbk, int *XBcond, 
		      int *YBcond, int *ZBcond, int *ier)
{
  /* allocate preconditioner data */
  VPDuData = VPrecDuAlloc(*Nx, *Ny, *Nz, *dx, *dy, *dz, *Ns, *NGx, 
			  *NGy, *NGz, *NPx, *NPy, *NPz, *iPx, 
			  *iPy, *iPz, *NBlt, *NBrt, *NBtp, *NBbt, 
			  *NBft, *NBbk, *XBcond, *YBcond, *ZBcond);
  if (VPDuData == NULL) *ier = -1; 
  else                  *ier = 0;

  return;
}



/* Du preconditioner dataspace deallocation wrapper routine */
void VISCPREC_DU_FREE()
{
  VPrecDuFree(VPDuData);
  return;
}



/* Du preconditioner setup wrapper routine */
void VISCPREC_DU_SETUP(double *uu, double *gamdt, double *Mu, 
		       double *Re, double *v1, double *v2, int *ier)
{
  /* call the C preconditioner setup routine */
  *ier = VPrecDuSetup(uu, *gamdt, *Mu, *Re, v1, v2, VPDuData);

  return;
}



/* momentum preconditioner solve wrapper routine */
void VISCPREC_DU_SOLVE(double *xx, double *bb, double *tmp, double *delta, int *ier)
{
  /* call the C preconditioner solve routine */
  *ier = VPrecDuSolve(xx, bb, tmp, *delta, VPDuData);

  return;
}



/* momentum preconditioner multiply wrapper routine */
void VISCPREC_DU_MULTIPLY(double *xx, double *bb, double *tmp, int *ier)
{
  /* call the C preconditioner multiply routine */
  *ier = VPrecDuMultiply(xx, bb, tmp, VPDuData);

  return;
}



/* Du preconditioner options routine */
void SET_SOL_DU_OPTS(int *iopt, int *ier)
{
  if (VPDuData == NULL) *ier = 1;
  else {
    /* cast VPDuData as the correct structure */
    ViscPrecDuData pdata;
    pdata = (ViscPrecDuData) VPDuData;
    /* set maximum number of iterations into pdata */
    if (iopt[0] != 0) pdata->sol_maxit = iopt[0];
    /* set relative change stopping criteria into pdata */
    if (iopt[1] != 0) pdata->sol_relch = iopt[1];
    /* set relaxation type into pdata */
    if (iopt[2] != 0) pdata->sol_rlxtype = iopt[2];
    /* set number of pre-relaxation sweeps into pdata */
    if (iopt[3] != 0) pdata->sol_npre = iopt[3];
    /* set number of post-relaxation sweeps into pdata */
    if (iopt[4] != 0) pdata->sol_npost = iopt[4];
    /* set print level into pdata */
    if (iopt[5] != 0) pdata->sol_printl = iopt[5];
    /* set amount of logging into pdata */
    if (iopt[6] != 0) pdata->sol_log = iopt[6];
    /* set choice of zero initial guess */
    if (iopt[7] != 0) pdata->sol_zeroguess = 1;
    *ier = 0;
  }
  return;
}



/* Du preconditioner diagnostic output routine */
void VISCPREC_DU_NUMITERS(int *Niters)
{
  *Niters = VPrecDuNumIters(VPDuData);
  return;
}




/********************************************************************/
/* Internal Preconditioner Routines                                 */
/********************************************************************/


/* -----------------------------------------------------------------
 * Function : VPrecDuAlloc
 * -----------------------------------------------------------------
 * VPrecDuAlloc is called at initialization to set aside space for 
 * any internal storage that will be required by VPrecDuSetup and
 * VPrecDuSolve.
 * -------------------------------------------------------------- */
void *VPrecDuAlloc(long int Nx, long int Ny, long int Nz, 
		   double dx, double dy, double dz, 
		   long int Ns, long int NGx, long int NGy, 
		   long int NGz, int NPx, int NPy, int NPz, 
		   int iPx, int iPy, int iPz, int NBlt, int NBrt, 
		   int NBtp, int NBbt, int NBft, int NBbk, 
		   int XBcond, int YBcond, int ZBcond) 
{
  /* define necessary local variables, output variable */
  ViscPrecDuData pdata;
  int xtag, ytag, ztag;
  int ilower[3], iupper[3];
  long int Nxl, Nyl, Nzl;

  /* allocate preconditioner data, cast as ViscPrecData */
  pdata = (ViscPrecDuData) malloc(sizeof *pdata);
  if (pdata == NULL) return(NULL);

  /* local domain information */
  pdata->ndim = 3;          /* 3D grid */
  pdata->Nx   = Nx;         /* num points in x-dir. */
  pdata->Ny   = Ny;         /* num points in y-dir. */
  pdata->Nz   = Nz;         /* num points in z-dir. */
  pdata->dx   = dx;         /* mesh size in x-dir.  */
  pdata->dy   = dy;         /* mesh size in y-dir.  */
  pdata->dz   = dz;         /* mesh size in z-dir.  */
  pdata->Ns   = Ns;         /* num species          */
  pdata->NGx  = NGx;        /* num x-ghost points   */
  pdata->NGy  = NGy;        /* num y-ghost points   */
  pdata->NGz  = NGz;        /* num z-ghost points   */
  pdata->xbc  = XBcond;     /* x-boundary condition */
  pdata->ybc  = YBcond;     /* y-boundary condition */
  pdata->zbc  = ZBcond;     /* z-boundary condition */

  /* processor layout information */
  pdata->Xprocs = NPx;      /* num procs in x-dir. */
  pdata->Yprocs = NPy;      /* num procs in y-dir. */
  pdata->Zprocs = NPz;      /* num procs in z-dir. */
  pdata->iprocx = iPx;      /* x-loc in proc. grid */
  pdata->iprocy = iPy;      /* y-loc in proc. grid */
  pdata->iprocz = iPz;      /* z-loc in proc. grid */
  if ((iPx==1) && (iPy==1) && (iPz==1)) {pdata->outproc=1;}
  else {pdata->outproc=0;}

  /* global domain information */
  Nxl  = 0;    /* size of left neighbor's x-grid */
  Nyl  = 0;    /* size of bottom neighbor's y-grid */
  Nzl  = 0;    /* size of back neighbor's y-grid */
  xtag = 101;  /* x-communication tag */
  ytag = 202;  /* y-communication tag */
  ztag = 303;  /* z-communication tag */
#ifdef PARALLEL
  pdata->comm = MPI_COMM_WORLD;
  MPI_Status status;
  if (iPx < NPx) 
    {MPI_Send(&Nx, 1, MPI_LONG, NBrt, xtag, pdata->comm);}
  if (iPx > 1) 
    {MPI_Recv(&Nxl, 1, MPI_LONG, NBlt, xtag, pdata->comm, &status);}
  if (iPy < NPy) 
    {MPI_Send(&Ny, 1, MPI_LONG, NBtp, ytag, pdata->comm);}
  if (iPy > 1) 
    {MPI_Recv(&Nyl, 1, MPI_LONG, NBbt, ytag, pdata->comm, &status);}
  if (iPz < NPz) 
    {MPI_Send(&Nz, 1, MPI_LONG, NBft, ztag, pdata->comm);}
  if (iPz > 1) 
    {MPI_Recv(&Nzl, 1, MPI_LONG, NBbk, ztag, pdata->comm, &status);}
#else
  pdata->comm = 0;
#endif
  pdata->iXL = Nxl*(iPx-1) + 1;
  pdata->iXR = Nxl*(iPx-1) + Nx;
  pdata->iYL = Nyl*(iPy-1) + 1;
  pdata->iYR = Nyl*(iPy-1) + Ny;
  pdata->iZL = Nzl*(iPz-1) + 1;
  pdata->iZR = Nzl*(iPz-1) + Nz;


  /* HYPRE-specific information */

  /*     set the matrix type: HYPRE_SSTRUCT or HYPRE_PARCSR */
  pdata->mattype = HYPRE_SSTRUCT;

  /*    initializing solver diagnostic information */
  pdata->totIters = 0;

  /*    set up the grid */
  /*       create grid object */
  HYPRE_SStructGridCreate(pdata->comm, pdata->ndim, 1, &(pdata->grid));

  /*       set my grid extents as if we have one part with multiple boxes.
	   Have each processor describe it's own global extents */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructGridSetExtents(pdata->grid, 0, ilower, iupper);

    /*       set grid variables for this part */
  int vartypes[3] = {HYPRE_SSTRUCT_VARIABLE_CELL,
		     HYPRE_SSTRUCT_VARIABLE_CELL,
		     HYPRE_SSTRUCT_VARIABLE_CELL};
  HYPRE_SStructGridSetVariables(pdata->grid, 0, 3, vartypes);

  /*       set grid periodicity */
  /*         [this currently does not work simply in the interface,  */
  /*          so we must manually set neighbor boxes into the grid]  */
  int periodicity[3] = {0, 0, 0};
  if (XBcond == 1) {
    long int xlo_idx, xhi_idx;
    MPI_Allreduce(&(pdata->iXL), &xlo_idx, 1, MPI_LONG, MPI_MIN, pdata->comm);
    MPI_Allreduce(&(pdata->iXR), &xhi_idx, 1, MPI_LONG, MPI_MAX, pdata->comm);
    periodicity[0] = xhi_idx-xlo_idx+1;
  }
  if (YBcond == 1) {
    long int ylo_idx, yhi_idx;
    MPI_Allreduce(&(pdata->iYL), &ylo_idx, 1, MPI_LONG, MPI_MIN, pdata->comm);
    MPI_Allreduce(&(pdata->iYR), &yhi_idx, 1, MPI_LONG, MPI_MAX, pdata->comm);
    periodicity[1] = yhi_idx-ylo_idx+1;
  }
  if (ZBcond == 1) {
    long int zlo_idx, zhi_idx;
    MPI_Allreduce(&(pdata->iZL), &zlo_idx, 1, MPI_LONG, MPI_MIN, pdata->comm);
    MPI_Allreduce(&(pdata->iZR), &zhi_idx, 1, MPI_LONG, MPI_MAX, pdata->comm);
    periodicity[2] = zhi_idx-zlo_idx+1;
  }
  HYPRE_SStructGridSetPeriodic(pdata->grid, 0, periodicity);

  /*       assemble the grid */
  HYPRE_SStructGridAssemble(pdata->grid);
  
  /*    set up the stencils */
  pdata->wStSize = 25;

  /*       create wx stencil */
  HYPRE_SStructStencilCreate(pdata->ndim, pdata->wStSize, &(pdata->wxStencil));
  /*       create wy stencil */
  HYPRE_SStructStencilCreate(pdata->ndim, pdata->wStSize, &(pdata->wyStencil));
  /*       create wz stencil */
  HYPRE_SStructStencilCreate(pdata->ndim, pdata->wStSize, &(pdata->wzStencil));

  /*       set stencil entries */
  int offset[3];
  /*         wx dependency on wx to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 0, offset, 0);
  /*         wx dependency on wx to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 1, offset, 0);
  /*         wx dependency on wx to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 2, offset, 0);
  /*         wx dependency on wx to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 3, offset, 0);
  /*         wx dependency on wx to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 4, offset, 0);
  /*         wx dependency on wx to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 5, offset, 0);
  /*         wx dependency on wx to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 6, offset, 0);

  /*         wx dependency on wy to left bottom */
  offset[0] = -1;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 7, offset, 1);
  /*         wx dependency on wy to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 8, offset, 1);
  /*         wx dependency on wy to right bottom*/
  offset[0] = 1;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 9, offset, 1);
  /*         wx dependency on wy to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 10, offset, 1);
  /*         wx dependency on wy to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 11, offset, 1);
  /*         wx dependency on wy to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 12, offset, 1);
  /*         wx dependency on wy to left top */
  offset[0] = -1;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 13, offset, 1);
  /*         wx dependency on wy to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 14, offset, 1);
  /*         wx dependency on wy to right top */
  offset[0] = 1;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 15, offset, 1);

  /*         wx dependency on wz to left back */
  offset[0] = -1;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 16, offset, 2);
  /*         wx dependency on wz to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 17, offset, 2);
  /*         wx dependency on wz to right back */
  offset[0] = 1;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 18, offset, 2);
  /*         wx dependency on wz to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 19, offset, 2);
  /*         wx dependency on wz to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 20, offset, 2);
  /*         wx dependency on wz to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 21, offset, 2);
  /*         wx dependency on wz to left front */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 22, offset, 2);
  /*         wx dependency on wz to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 23, offset, 2);
  /*         wx dependency on wz to right front */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wxStencil, 24, offset, 2);

  /*         wy dependency on wy to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 0, offset, 1);
  /*         wy dependency on wy to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 1, offset, 1);
  /*         wy dependency on wy to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 2, offset, 1);
  /*         wy dependency on wy to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 3, offset, 1);
  /*         wy dependency on wy to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 4, offset, 1);
  /*         wy dependency on wy to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 5, offset, 1);
  /*         wy dependency on wy to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 6, offset, 1);

  /*         wy dependency on wz to bottom back */
  offset[0] = 0;  offset[1] = -1;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 7, offset, 2);
  /*         wy dependency on wz to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 8, offset, 2);
  /*         wy dependency on wz to top back */
  offset[0] = 0;  offset[1] = 1;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 9, offset, 2);
  /*         wy dependency on wz to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 10, offset, 2);
  /*         wy dependency on wz to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 11, offset, 2);
  /*         wy dependency on wz to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 12, offset, 2);
  /*         wy dependency on wz to bottom front */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 13, offset, 2);
  /*         wy dependency on wz to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 14, offset, 2);
  /*         wy dependency on wz to top front */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 15, offset, 2);

  /*         wy dependency on wx to left bottom */
  offset[0] = -1;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 16, offset, 0);
  /*         wy dependency on wx to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 17, offset, 0);
  /*         wy dependency on wx to right bottom */
  offset[0] = 1;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 18, offset, 0);
  /*         wy dependency on wx to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 19, offset, 0);
  /*         wy dependency on wx to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 20, offset, 0);
  /*         wy dependency on wx to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 21, offset, 0);
  /*         wy dependency on wx to left top */
  offset[0] = -1;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 22, offset, 0);
  /*         wy dependency on wx to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 23, offset, 0);
  /*         wy dependency on wx to right top */
  offset[0] = 1;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wyStencil, 24, offset, 0);

  /*         wz dependency on wz to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 0, offset, 2);
  /*         wz dependency on wz to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 1, offset, 2);
  /*         wz dependency on wz to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 2, offset, 2);
  /*         wz dependency on wz to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 3, offset, 2);
  /*         wz dependency on wz to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 4, offset, 2);
  /*         wz dependency on wz to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 5, offset, 2);
  /*         wz dependency on wz to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 6, offset, 2);

  /*         wz dependency on wx to left back */
  offset[0] = -1;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 7, offset, 0);
  /*         wz dependency on wx to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 8, offset, 0);
  /*         wz dependency on wx to right back */
  offset[0] = 1;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 9, offset, 0);
  /*         wz dependency on wx to left */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 10, offset, 0);
  /*         wz dependency on wx to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 11, offset, 0);
  /*         wz dependency on wx to right */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 12, offset, 0);
  /*         wz dependency on wx to left front */
  offset[0] = -1;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 13, offset, 0);
  /*         wz dependency on wx to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 14, offset, 0);
  /*         wz dependency on wx to right front */
  offset[0] = 1;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 15, offset, 0);

  /*         wz dependency on wy to bottom back */
  offset[0] = 0;  offset[1] = -1;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 16, offset, 1);
  /*         wz dependency on wy to back */
  offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 17, offset, 1);
  /*         wz dependency on wy to top back */
  offset[0] = 0;  offset[1] = 1;  offset[2] = -1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 18, offset, 1);
  /*         wz dependency on wy to bottom */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 19, offset, 1);
  /*         wz dependency on wy to self */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 20, offset, 1);
  /*         wz dependency on wy to top */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 21, offset, 1);
  /*         wz dependency on wy to bottom front */
  offset[0] = 0;  offset[1] = -1;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 22, offset, 1);
  /*         wz dependency on wy to front */
  offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 23, offset, 1);
  /*         wz dependency on wy to top front */
  offset[0] = 0;  offset[1] = 1;  offset[2] = 1;
  HYPRE_SStructStencilSetEntry(pdata->wzStencil, 24, offset, 1);
  
  /*    set up the graph */
  /*       create the graph object */
  HYPRE_SStructGraphCreate(pdata->comm, pdata->grid, &(pdata->graph));
  
  /*       set graph type according to solver desired */
  HYPRE_SStructGraphSetObjectType(pdata->graph, pdata->mattype);
  
  /*       set stencils into graph */
  /*          set wx stencil */
  HYPRE_SStructGraphSetStencil(pdata->graph, 0, 0, pdata->wxStencil);
  /*          set wy stencil */
  HYPRE_SStructGraphSetStencil(pdata->graph, 0, 1, pdata->wyStencil);
  /*          set wz stencil */
  HYPRE_SStructGraphSetStencil(pdata->graph, 0, 2, pdata->wzStencil);

  /*       add additional non-stencil entries into graph */
  /*       (none that I can think of) */

  /*       assemble the graph */
  HYPRE_SStructGraphAssemble(pdata->graph);


  /********************************/
  /*  continue with general setup */

  /*    set Du, b and x init flags to 0 at first */
  pdata->DuInit = 0;

  /* set Solver default options into pdata */
  pdata->sol_maxit     = 50;
  pdata->sol_relch     = 0;
  pdata->sol_rlxtype   = 1;
  pdata->sol_npre      = 1;
  pdata->sol_npost     = 1;
  pdata->sol_printl    = 1;
  pdata->sol_log       = 1;
  pdata->sol_zeroguess = 0;

  /* return preconditioner data, recast as void *. */
  return((void *)pdata);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuFree
 * -----------------------------------------------------------------
 * VPrecDuFree frees storage allocated by VPrecDuAlloc
 * -------------------------------------------------------------- */
void VPrecDuFree(void *P_data)
{
  /* ensure that P_data is non-null, and free space as required */
  if ( P_data != NULL ) {

    /* cast P_data as the correct structure */
    ViscPrecDuData pdata;
    pdata = (ViscPrecDuData) P_data;

    /* free the various components */
    HYPRE_SStructMatrixDestroy(pdata->Du);
    HYPRE_SStructGraphDestroy(pdata->graph);
    HYPRE_SStructStencilDestroy(pdata->wzStencil);
    HYPRE_SStructStencilDestroy(pdata->wyStencil);
    HYPRE_SStructStencilDestroy(pdata->wxStencil);
    HYPRE_SStructGridDestroy(pdata->grid);

    /* finally, free the pdata structure */
    free(pdata);
  }
}



/* -----------------------------------------------------------------
 * Function : VPrecDuSetup
 * -----------------------------------------------------------------
 * VPrecDuSetup sets up the viscous preconditioning matrix for the 
 * momentum equations.
 *
 * The parameters of VPrecSetup used here are as follows:
 *
 * uu      is the current state of the system
 * gamdt   is the time scaling in the Newton matrix, M = I+gamdt*J
 * Mu      is the plasma viscosity coefficient
 * Re      is the plasma Reynolds number
 * v1, v2  give already allocated arrays which may be used as 
 *         temporary storage or work space
 * pdata   is the preconditioner data returned by VPrecAlloc.
 *
 * Return value:
 * The value returned by this VPrecDuSetup function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDuSetup(double *uu, double gamdt, double Mu, 
		 double Re, double *v1, double *v2, void *P_data)
{
  /* local variables */
  long int Nx, Ny, Nz, NGx, NGy, NGz;
  long int ix, iy, iz, idx, Zbl, Ybl;
  int ilower[3], iupper[3], selfentries[7], nextentries[9], lastentries[9];
  int xLface, xRface, yLface, yRface, zLface, zRface, IdLoc;
  double dx, dy, dz, dxi2, dyi2, dzi2, dxdyfac, dxdzfac, dydzfac, rhofact, IdVal;
  double *matvals, uxvals[25], uyvals[25], uzvals[25];


  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;

  /* destroy old matrix if necessary */
  if (pdata->DuInit == 1) {
    HYPRE_SStructMatrixDestroy(pdata->Du);
    pdata->DuInit = 0;
  }
    
  /* create the SStruct matrix, and set init flag */
  HYPRE_SStructMatrixCreate(pdata->comm, pdata->graph, &(pdata->Du));
  pdata->DuInit = 1;

  /* set matrix storage type */
  HYPRE_SStructMatrixSetObjectType(pdata->Du, pdata->mattype);

/*   /\* set matrix symmetry *\/ */
/*   HYPRE_SStructMatrixSetSymmetric(pdata->Du, 0, 0, 0, 1); */
/*   HYPRE_SStructMatrixSetSymmetric(pdata->Du, 0, 1, 1, 1); */
/*   HYPRE_SStructMatrixSetSymmetric(pdata->Du, 0, 2, 2, 1); */
    
  /* initialize matrix */
  HYPRE_SStructMatrixInitialize(pdata->Du);

  /* get grid information shortcuts */
  Nx  = pdata->Nx;
  Ny  = pdata->Ny;
  Nz  = pdata->Nz;
  NGx = pdata->NGx;
  NGy = pdata->NGy;
  NGz = pdata->NGz;
  dx  = pdata->dx;
  dy  = pdata->dy;
  dz  = pdata->dz;
  dxi2 = 1.0/dx/dx;
  dyi2 = 1.0/dy/dy;
  dzi2 = 1.0/dz/dz;
  dxdyfac = 1.0/dx/dy/12.0;
  dxdzfac = 1.0/dx/dz/12.0;
  dydzfac = 1.0/dy/dz/12.0;

  /* allocate temporary array containing matrix values */
  /* (since the temporary vectors are only 8*Nx*Ny*Nz long) */
  matvals = (double *) malloc(9*Nx*Ny*Nz*sizeof(double));


  /* uxvals holds [unscaled] template for ux stencil */
  uxvals[0]  = dzi2;                              /* ux, back */
  uxvals[1]  = dyi2;                              /* ux, bottom */
  uxvals[2]  = dxi2*4.0/3.0;                      /* ux, left */
  uxvals[3]  = -2.0*(4.0/3.0*dxi2 + dyi2 + dzi2); /* ux, self */
  uxvals[4]  = dxi2*4.0/3.0;                      /* ux, right */
  uxvals[5]  = dyi2;                              /* ux, top */
  uxvals[6]  = dzi2;                              /* ux, front */
  uxvals[7]  = dxdyfac;                           /* uy, left bottom */
  uxvals[8]  = 0.0;                               /* uy, bottom */
  uxvals[9]  = -dxdyfac;                          /* uy, right bottom */
  uxvals[10] = 0.0;                               /* uy, left */
  uxvals[11] = 0.0;                               /* uy, self */
  uxvals[12] = 0.0;                               /* uy, right */
  uxvals[13] = -dxdyfac;                          /* uy, left top */
  uxvals[14] = 0.0;                               /* uy, top */
  uxvals[15] = dxdyfac;                           /* uy, right top */
  uxvals[16] = dxdzfac;                           /* uz, left back */
  uxvals[17] = 0.0;                               /* uz, back */
  uxvals[18] = -dxdzfac;                          /* uz, right back */
  uxvals[19] = 0.0;                               /* uz, left */
  uxvals[20] = 0.0;                               /* uz, self */
  uxvals[21] = 0.0;                               /* uz, right */
  uxvals[22] = -dxdzfac;                          /* uz, left front */
  uxvals[23] = 0.0;                               /* uz, front */
  uxvals[24] = dxdzfac;                           /* uz, right front */

  /* uyvals holds [unscaled] template for uy stencil */
  uyvals[0]  = dzi2;                              /* uy, back */
  uyvals[1]  = dyi2*4.0/3.0;                      /* uy, bottom */
  uyvals[2]  = dxi2;                              /* uy, left */
  uyvals[3]  = -2.0*(dxi2 + 4.0/3.0*dyi2 + dzi2); /* uy, self */
  uyvals[4]  = dxi2;                              /* uy, right */
  uyvals[5]  = dyi2*4.0/3.0;                      /* uy, top */
  uyvals[6]  = dzi2;                              /* uy, front */
  uyvals[7]  = dydzfac;                           /* uz, bottom back */
  uyvals[8]  = 0.0;                               /* uz, back */
  uyvals[9]  = -dydzfac;                          /* uz, top back */
  uyvals[10] = 0.0;                               /* uz, bottom */
  uyvals[11] = 0.0;                               /* uz, self */
  uyvals[12] = 0.0;                               /* uz, top */
  uyvals[13] = -dydzfac;                          /* uz, bottom front */
  uyvals[14] = 0.0;                               /* uz, front */
  uyvals[15] = dydzfac;                           /* uz, top front */
  uyvals[16] = dxdyfac;                           /* ux, left bottom */
  uyvals[17] = 0.0;                               /* ux, bottom */
  uyvals[18] = -dxdyfac;                          /* ux, right bottom */
  uyvals[19] = 0.0;                               /* ux, left */
  uyvals[20] = 0.0;                               /* ux, self */
  uyvals[21] = 0.0;                               /* ux, right */
  uyvals[22] = -dxdyfac;                          /* ux, left top */
  uyvals[23] = 0.0;                               /* ux, top */
  uyvals[24] = dxdyfac;                           /* ux, right top */

  /* uzvals holds [unscaled] template for uz stencil */
  uzvals[0]  = dzi2*4.0/3.0;                      /* uz, back */
  uzvals[1]  = dyi2;                              /* uz, bottom */
  uzvals[2]  = dxi2;                              /* uz, left */
  uzvals[3]  = -2.0*(dxi2 + dyi2 + 4.0/3.0*dzi2); /* uz, self */
  uzvals[4]  = dxi2;                              /* uz, right */
  uzvals[5]  = dyi2;                              /* uz, top */
  uzvals[6]  = dzi2*4.0/3.0;                      /* uz, front */
  uzvals[7]  = dxdzfac;                           /* ux, left back */
  uzvals[8]  = 0.0;                               /* ux, back */
  uzvals[9]  = -dxdzfac;                          /* ux, right back */
  uzvals[10] = 0.0;                               /* ux, left */
  uzvals[11] = 0.0;                               /* ux, self */
  uzvals[12] = 0.0;                               /* ux, right */
  uzvals[13] = -dxdzfac;                          /* ux, left front */
  uzvals[14] = 0.0;                               /* ux, front */
  uzvals[15] = dxdzfac;                           /* ux, right front */
  uzvals[16] = dydzfac;                           /* uy, bottom back */
  uzvals[17] = 0.0;                               /* uy, back */
  uzvals[18] = -dydzfac;                          /* uy, top back */
  uzvals[19] = 0.0;                               /* uy, bottom */
  uzvals[20] = 0.0;                               /* uy, self */
  uzvals[21] = 0.0;                               /* uy, top */
  uzvals[22] = -dydzfac;                          /* uy, bottom front */
  uzvals[23] = 0.0;                               /* uy, front */
  uzvals[24] = dydzfac;                           /* uy, top front */

  /* selfentries holds stencil locations for variable couplings to self */
  selfentries[0] = 0;
  selfentries[1] = 1;
  selfentries[2] = 2;
  selfentries[3] = 3;
  selfentries[4] = 4;
  selfentries[5] = 5;
  selfentries[6] = 6;

  /* nextentries holds stencil locations for variable couplings to next var */
  nextentries[0] = 7;
  nextentries[1] = 8;
  nextentries[2] = 9;
  nextentries[3] = 10;
  nextentries[4] = 11;
  nextentries[5] = 12;
  nextentries[6] = 13;
  nextentries[7] = 14;
  nextentries[8] = 15;

  /* lastentries holds stencil locations for variable couplings to last var */
  lastentries[0] = 16;
  lastentries[1] = 17;
  lastentries[2] = 18;
  lastentries[3] = 19;
  lastentries[4] = 20;
  lastentries[5] = 21;
  lastentries[6] = 22;
  lastentries[7] = 23;
  lastentries[8] = 24;

  /* IdLoc, IdVal hold identity location, value */
  IdLoc = 3;
  IdVal = 1.0;

  /* set flags determining whether proc owns external faces */
  xLface = (pdata->iprocx == 1);
  xRface = (pdata->iprocx == pdata->Xprocs);
  yLface = (pdata->iprocy == 1);
  yRface = (pdata->iprocy == pdata->Yprocs);
  zLface = (pdata->iprocz == 1);
  zRface = (pdata->iprocz == pdata->Zprocs);


  /* set ux stencil couplings to self */
  /*       internal cells */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uxvals[0];
	matvals[idx++] = rhofact*uxvals[1];
	matvals[idx++] = rhofact*uxvals[2];
	matvals[idx++] = rhofact*uxvals[3];
	matvals[idx++] = rhofact*uxvals[4];
	matvals[idx++] = rhofact*uxvals[5];
	matvals[idx++] = rhofact*uxvals[6];
      }
    }
  }
  /*       ix=0 face adjustment */
  if (xLface) {
    ix=0;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    else if (pdata->xbc == 0) {  /* zero-gradient */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
  }
  /*       ix=Nx-1 face adjustment */
  if (xRface) {
    ix=Nx-1;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+4];  matvals[idx+4] = 0.0;
	}
      }
    }
    else if (pdata->xbc == 0) {  /* zero-gradient */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+4];  matvals[idx+4] = 0.0;
	}
      }
    }
  }
  if (pdata->ybc != 1) {       /* not periodic */
    /*       iy=0 face adjustment */
    if (yLface) {
      iy=0;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero-gradient the same at this face */
	  matvals[idx+3] += matvals[idx+1];  matvals[idx+1] = 0.0;
	}
      }
    }
    /*       iy=Ny-1 face adjustment */
    if (yRface) {
      iy=Ny-1;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero-gradient the same at this face */
	  matvals[idx+3] += matvals[idx+5];  matvals[idx+5] = 0.0;
	}
      }
    }
  }
  if (pdata->zbc != 1) {       /* not periodic */
    /*       iz=0 face adjustment */
    if (zLface) {
      iz=0;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero-gradient the same at this face */
	  matvals[idx+3] += matvals[idx];  matvals[idx] = 0.0;
	}
      }
    }
    /*       iz=Nz-1 face adjustment */
    if (zRface) {
      iz=Nz-1;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero-gradient the same at this face */
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 0, 7, selfentries, matvals);


  /* set ux stencil couplings to uy */
  /*       internal cells */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uxvals[7];
	matvals[idx++] = rhofact*uxvals[8];
	matvals[idx++] = rhofact*uxvals[9];
	matvals[idx++] = rhofact*uxvals[10];
	matvals[idx++] = rhofact*uxvals[11];
	matvals[idx++] = rhofact*uxvals[12];
	matvals[idx++] = rhofact*uxvals[13];
	matvals[idx++] = rhofact*uxvals[14];
	matvals[idx++] = rhofact*uxvals[15];
      }
    }
  }
  if (pdata->xbc != 1) {       /* not periodic */
    /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+1] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+1] += matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] += matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  /*       iy=0 face adjustment */
  if (yLface) {
    iy=0;
    if (pdata->ybc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] -= matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (yRface) {
    iy=Ny-1;
    if (pdata->ybc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] -= matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] -= matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] += matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 0, 9, nextentries, matvals);


  /* set ux stencil couplings to uz */
  /*       internal cells */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uxvals[16];
	matvals[idx++] = rhofact*uxvals[17];
	matvals[idx++] = rhofact*uxvals[18];
	matvals[idx++] = rhofact*uxvals[19];
	matvals[idx++] = rhofact*uxvals[20];
	matvals[idx++] = rhofact*uxvals[21];
	matvals[idx++] = rhofact*uxvals[22];
	matvals[idx++] = rhofact*uxvals[23];
	matvals[idx++] = rhofact*uxvals[24];
      }
    }
  }
  if (pdata->xbc != 1) {      /* not periodic */
    /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+1] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+1] += matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] += matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  /*       iz=0 face adjustment */
  if (zLface) {
    iz=0;
    if (pdata->zbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] -= matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    else if (pdata->zbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
  }
  /*       iz=Nz-1 face adjustment */
  if (zRface) {
    iz=Nz-1;
    if (pdata->zbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] -= matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] -= matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
    else if (pdata->zbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] += matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 0, 9, lastentries, matvals);


  /* add one to matrix diagonal for identity contribution */
  for (ix=0; ix<Nx*Ny*Nz; ix++)  matvals[ix] = IdVal;
  HYPRE_SStructMatrixAddToBoxValues(pdata->Du, 0, ilower, 
				    iupper, 0, 1, &IdLoc, matvals);



  /* set uy stencil couplings to uy */
  /*       internal cells */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uyvals[0];
	matvals[idx++] = rhofact*uyvals[1];
	matvals[idx++] = rhofact*uyvals[2];
	matvals[idx++] = rhofact*uyvals[3];
	matvals[idx++] = rhofact*uyvals[4];
	matvals[idx++] = rhofact*uyvals[5];
	matvals[idx++] = rhofact*uyvals[6];
      }
    }
  }
  if (pdata->xbc != 1) {      /* not periodic */
  /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+4];  matvals[idx+4] = 0.0;
	}
      }
    }
  }
  /*       iy=0 face adjustment */
  if (yLface) {
    iy=0;
    if (pdata->ybc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+1];  matvals[idx+1] = 0.0;
	}
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+1];  matvals[idx+1] = 0.0;
	}
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (yRface) {
    iy=Ny-1;
    if (pdata->ybc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+5];  matvals[idx+5] = 0.0;
	}
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+5];  matvals[idx+5] = 0.0;
	}
      }
    }
  }
  if (pdata->zbc != 1) {      /* not periodic */
    /*       iz=0 face adjustment */
    if (zLface) {
      iz=0;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx];  matvals[idx] = 0.0;
	}
      }
    }
    /*       iz=Nz-1 face adjustment */
    if (zRface) {
      iz=Nz-1;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 1, 7, selfentries, matvals);


  /* set uy stencil couplings to uz */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uyvals[7];
	matvals[idx++] = rhofact*uyvals[8];
	matvals[idx++] = rhofact*uyvals[9];
	matvals[idx++] = rhofact*uyvals[10];
	matvals[idx++] = rhofact*uyvals[11];
	matvals[idx++] = rhofact*uyvals[12];
	matvals[idx++] = rhofact*uyvals[13];
	matvals[idx++] = rhofact*uyvals[14];
	matvals[idx++] = rhofact*uyvals[15];
      }
    }
  }
  if (pdata->ybc != 1) {      /* not periodic */
    /*       iy=0 face adjustment */
    if (yLface) {
      iy=0;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+1] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    /*       iy=Ny-1 face adjustment */
    if (yRface) {
      iy=Ny-1;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+1] += matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] += matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  /*       iz=0 face adjustment */
  if (zLface) {
    iz=0;
    if (pdata->zbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] -= matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    else if (pdata->zbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
  }
  /*       iz=Nz-1 face adjustment */
  if (zRface) {
    iz=Nz-1;
    if (pdata->zbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] -= matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] -= matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
    else if (pdata->zbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] += matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 1, 9, nextentries, matvals);


  /* set uy stencil couplings to ux */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uyvals[16];
	matvals[idx++] = rhofact*uyvals[17];
	matvals[idx++] = rhofact*uyvals[18];
	matvals[idx++] = rhofact*uyvals[19];
	matvals[idx++] = rhofact*uyvals[20];
	matvals[idx++] = rhofact*uyvals[21];
	matvals[idx++] = rhofact*uyvals[22];
	matvals[idx++] = rhofact*uyvals[23];
	matvals[idx++] = rhofact*uyvals[24];
      }
    }
  }
  /*       ix=0 face adjustment */
  if (xLface) {
    ix=0;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] -= matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] -= matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    else if (pdata->xbc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
  }
  /*       ix=Nx-1 face adjustment */
  if (xRface) {
    ix=Nx-1;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] -= matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] -= matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
    else if (pdata->xbc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] += matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] += matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  if (pdata->ybc != 1) {       /* not periodic */
    /*       iy=0 face adjustment */
    if (yLface) {
      iy=0;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    /*       iy=Ny-1 face adjustment */
    if (yRface) {
      iy=Ny-1;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] += matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 1, 9, lastentries, matvals);


  /* add one to matrix diagonal for identity contribution */
  for (ix=0; ix<Nx*Ny*Nz; ix++)  matvals[ix] = IdVal;
  HYPRE_SStructMatrixAddToBoxValues(pdata->Du, 0, ilower, 
				    iupper, 1, 1, &IdLoc, matvals);



  /* set uz stencil couplings to uz */
  /*       internal cells */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uzvals[0];
	matvals[idx++] = rhofact*uzvals[1];
	matvals[idx++] = rhofact*uzvals[2];
	matvals[idx++] = rhofact*uzvals[3];
	matvals[idx++] = rhofact*uzvals[4];
	matvals[idx++] = rhofact*uzvals[5];
	matvals[idx++] = rhofact*uzvals[6];
      }
    }
  }
  if (pdata->xbc != 1) {       /* not periodic */
    /*       ix=0 face adjustment */
    if (xLface) {
      ix=0;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    /*       ix=Nx-1 face adjustment */
    if (xRface) {
      ix=Nx-1;
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+4];  matvals[idx+4] = 0.0;
	}
      }
    }
  }
  if (pdata->ybc != 1) {       /* not periodic */
    /*       iy=0 face adjustment */
    if (yLface) {
      iy=0;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+1];  matvals[idx+1] = 0.0;
	}
      }
    }
    /*       iy=Ny-1 face adjustment */
    if (yRface) {
      iy=Ny-1;
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+5];  matvals[idx+5] = 0.0;
	}
      }
    }
  }
  /*       iz=0 face adjustment */
  if (zLface) {
    iz=0;
    if (pdata->zbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+0];  matvals[idx+0] = 0.0;
	}
      }
    }
    else if (pdata->zbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+0];  matvals[idx+0] = 0.0;
	}
      }
    }
  }
  /*       iz=Nz-1 face adjustment */
  if (zRface) {
    iz=Nz-1;
    if (pdata->zbc == 2) {       /* reflecting */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    else if (pdata->zbc == 0) {  /* zero gradient */
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 7*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 2, 7, selfentries, matvals);


  /* set uz stencil couplings to ux */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uzvals[7];
	matvals[idx++] = rhofact*uzvals[8];
	matvals[idx++] = rhofact*uzvals[9];
	matvals[idx++] = rhofact*uzvals[10];
	matvals[idx++] = rhofact*uzvals[11];
	matvals[idx++] = rhofact*uzvals[12];
	matvals[idx++] = rhofact*uzvals[13];
	matvals[idx++] = rhofact*uzvals[14];
	matvals[idx++] = rhofact*uzvals[15];
      }
    }
  }
  /*       ix=0 face adjustment */
  if (xLface) {
    ix=0;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] -= matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] -= matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    else if (pdata->xbc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
  }
  /*       ix=Nx-1 face adjustment */
  if (xRface) {
    ix=Nx-1;
    if (pdata->xbc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] -= matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] -= matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
    else if (pdata->xbc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (iy=0; iy<Ny; iy++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] += matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] += matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  if (pdata->zbc != 1) {       /* not periodic */
    /*       iz=0 face adjustment */
    if (zLface) {
      iz=0;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero grandient are the same at this face */
	  matvals[idx+3] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    /*       iz=Nz-1 face adjustment */
    if (zRface) {
      iz=Nz-1;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero grandient are the same at this face */
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] += matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 2, 9, nextentries, matvals);


  /* set uz stencil couplings to uy */
  idx = 0;
  for (iz=0; iz<Nz; iz++) {
    Zbl = (iz+NGz) * (Ny+2*NGy) * (Nx+2*NGx);
    for (iy=0; iy<Ny; iy++) {
      Ybl = (iy+NGy) * (Nx+2*NGx);
      for (ix=0; ix<Nx; ix++) {
	rhofact = -gamdt*Mu*uu[Zbl + Ybl + ix + NGx];
	matvals[idx++] = rhofact*uzvals[16];
	matvals[idx++] = rhofact*uzvals[17];
	matvals[idx++] = rhofact*uzvals[18];
	matvals[idx++] = rhofact*uzvals[19];
	matvals[idx++] = rhofact*uzvals[20];
	matvals[idx++] = rhofact*uzvals[21];
	matvals[idx++] = rhofact*uzvals[22];
	matvals[idx++] = rhofact*uzvals[23];
	matvals[idx++] = rhofact*uzvals[24];
      }
    }
  }
  /*       iy=0 face adjustment */
  if (yLface) {
    iy=0;
    if (pdata->ybc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] -= matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] -= matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] -= matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+3];  matvals[idx+3] = 0.0;
	  matvals[idx+7] += matvals[idx+6];  matvals[idx+6] = 0.0;
	}
      }
    }
  }
  /*       iy=Ny-1 face adjustment */
  if (yRface) {
    iy=Ny-1;
    if (pdata->ybc == 2) {       /* reflecting */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] -= matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] -= matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] -= matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
    else if (pdata->ybc == 0) {  /* zero gradient */
      for (iz=0; iz<Nz; iz++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  matvals[idx+1] += matvals[idx+2];  matvals[idx+2] = 0.0;
	  matvals[idx+4] += matvals[idx+5];  matvals[idx+5] = 0.0;
	  matvals[idx+7] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  if (pdata->zbc != 1) {       /* not periodic */
    /*       iz=0 face adjustment */
    if (zLface) {
      iz=0;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx];    matvals[idx]   = 0.0;
	  matvals[idx+4] += matvals[idx+1];  matvals[idx+1] = 0.0;
	  matvals[idx+5] += matvals[idx+2];  matvals[idx+2] = 0.0;
	}
      }
    }
    /*       iz=Nz-1 face adjustment */
    if (zRface) {
      iz=Nz-1;
      for (iy=0; iy<Ny; iy++) {
	for (ix=0; ix<Nx; ix++) {
	  idx = 9*((iz*Ny + iy)*Nx + ix);
	  /* reflecting and zero gradient are the same at this face */
	  matvals[idx+3] += matvals[idx+6];  matvals[idx+6] = 0.0;
	  matvals[idx+4] += matvals[idx+7];  matvals[idx+7] = 0.0;
	  matvals[idx+5] += matvals[idx+8];  matvals[idx+8] = 0.0;
	}
      }
    }
  }
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  HYPRE_SStructMatrixSetBoxValues(pdata->Du, 0, ilower, 
				  iupper, 2, 9, lastentries, matvals);


  /* add one to matrix diagonal for identity contribution */
  for (ix=0; ix<Nx*Ny*Nz; ix++)  matvals[ix] = IdVal;
  HYPRE_SStructMatrixAddToBoxValues(pdata->Du, 0, ilower, 
				    iupper, 2, 1, &IdLoc, matvals);


  /* assemble matrix */
  HYPRE_SStructMatrixAssemble(pdata->Du);
    
  /* free temporary array */
  free(matvals);


/*   if ((pdata->outproc)==1) */
/*     {printf("      printing HYPRE Du matrix to file \n");} */
/*   char *fname = "Du_precmat"; */
/*   HYPRE_SStructMatrixPrint(fname, pdata->Du, 0); */

  /* return success */ 
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuSolve
 * -----------------------------------------------------------------
 * VPrecDuSolve solves a linear system P z = r, with the
 * preconditioner matrix P generated by VPrecSetup and solved
 * using the HYPRE library.
 *
 * The parameters of VPrecDuSolve used here are as follows:
 *
 * xx      is the sol vector on output
 * bb      is the rhs vector on input
 * tmp     is a temporary vector the size of xx
 * delta   is the desired linear solve tolerance (if used in 
 *         some iterative method)
 * pdata   is the pre-computed Du preconditioner data
 *
 * The value returned by this VPrecDuSolve function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDuSolve(double *xx, double *bb, double *tmp, double delta, void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;

  /* local variables */
  int its, ilower[3], iupper[3];
  long int iz, iy, ix, iv, idx, Nx, NGx, Ny, NGy, Nz, NGz;
  long int Vbl, Ybl, Zbl;
  double finalresid, resid, val, bnorm;
  HYPRE_SStructVector bvec, xvec;
  HYPRE_SStructSolver solver;
  int printl = ((pdata->outproc)==1) ? pdata->sol_printl : 0;

/*   if (printl) printf("     solving Du preconditioner system\n"); */

  /* check that Du matrix initialized */
  if (pdata->DuInit == 0) {
    printf("VPrecDuSolve error: Du matrix uninitialized!\n");
    return(1);
  }

  /* set local variables */
  Nx  = pdata->Nx;
  NGx = pdata->NGx;
  Ny  = pdata->Ny;
  NGy = pdata->NGy;
  Nz  = pdata->Nz;
  NGz = pdata->NGz;

  /* create the SStruct vectors and set init flags to 1 */
  HYPRE_SStructVectorCreate(pdata->comm, pdata->grid, &bvec);
  HYPRE_SStructVectorCreate(pdata->comm, pdata->grid, &xvec);

  /* set vector storage type */
  HYPRE_SStructVectorSetObjectType(bvec, pdata->mattype);
  HYPRE_SStructVectorSetObjectType(xvec, pdata->mattype);
    
  /* initialize vectors */
  HYPRE_SStructVectorInitialize(bvec);
  HYPRE_SStructVectorInitialize(xvec);
  

  /* convert rhs, solution vectors to HYPRE format:         */
  /*    insert rhs vector entries into HYPRE vectors bvec   */
  /*    and xvec (use the rhs as the sol initial guess)     */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  for (iv=1; iv<=3; iv++) {
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
    idx = 0;
    for (iz=NGz; iz<Nz+NGz; iz++) {
      Zbl = iz * (Nx+2*NGx) * (Ny+2*NGy);
      for (iy=NGy; iy<Ny+NGy; iy++) {
	Ybl = iy * (Nx+2*NGx);
	for (ix=NGx; ix<Nx+NGx; ix++)
	  tmp[idx++] = bb[Vbl+Zbl+Ybl+ix];
      }
    }
    HYPRE_SStructVectorSetBoxValues(bvec, 0, ilower, iupper, iv-1, tmp);
    HYPRE_SStructVectorSetBoxValues(xvec, 0, ilower, iupper, iv-1, tmp);
  }

  /*    assemble vectors */
  HYPRE_SStructVectorAssemble(xvec);
  HYPRE_SStructVectorAssemble(bvec);

  /* set up the solver [SysPFMG for now] */
  /*    create the solver */
/*   if (printl) printf("      creating SysPFMG solver \n"); */
  HYPRE_SStructSysPFMGCreate(pdata->comm, &solver);
 
  /*    set solver options */
  /*    [could the first 8 of these be done in the PrecAlloc routine?] */
/*   if (printl) printf("      setting SysPFMG options \n"); */
  HYPRE_SStructSysPFMGSetMaxIter(solver, pdata->sol_maxit);
  HYPRE_SStructSysPFMGSetRelChange(solver, pdata->sol_relch);
  HYPRE_SStructSysPFMGSetRelaxType(solver, pdata->sol_rlxtype);
  HYPRE_SStructSysPFMGSetNumPreRelax(solver, pdata->sol_npre);
  HYPRE_SStructSysPFMGSetNumPostRelax(solver, pdata->sol_npost);
  HYPRE_SStructSysPFMGSetPrintLevel(solver, pdata->sol_printl);
  HYPRE_SStructSysPFMGSetLogging(solver, pdata->sol_log);
  if (delta != 0.0)  HYPRE_SStructSysPFMGSetTol(solver, delta);
  if (pdata->sol_zeroguess) 
    HYPRE_SStructSysPFMGSetZeroGuess(solver);
/*   if (printl) printf("      calling SysPFMG setup \n"); */
  HYPRE_SStructSysPFMGSetup(solver, pdata->Du, bvec, xvec);

  /* solve the linear system */
/*   if (printl) printf("      calling SysPFMG solver \n"); */
  HYPRE_SStructSysPFMGSolve(solver, pdata->Du, bvec, xvec);

  /* extract solver statistics */
/*   if (printl) printf("      extracting SysPFMG statistics \n"); */
  HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_SStructSysPFMGGetNumIterations(solver, &its);
  pdata->totIters += its;

  /* gather the solution vector and extract values */
/*   if (printl) printf("      extracting solution vector \n"); */
  HYPRE_SStructVectorGather(xvec);
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  for (iv=1; iv<=3; iv++) {
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
    HYPRE_SStructVectorGetBoxValues(xvec, 0, ilower, iupper, iv-1, tmp);
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
    idx = 0;
    for (iz=NGz; iz<Nz+NGz; iz++) {
      Zbl = iz * (Nx+2*NGx) * (Ny+2*NGy);
      for (iy=NGy; iy<Ny+NGy; iy++) {
	Ybl = iy * (Nx+2*NGx);
	for (ix=NGx; ix<Nx+NGx; ix++) 
	  xx[Vbl+Zbl+Ybl+ix] = tmp[idx++];
      }
    }
  }

  /* destroy vector and solver structures */
  HYPRE_SStructSysPFMGDestroy(solver);
  HYPRE_SStructVectorDestroy(bvec);
  HYPRE_SStructVectorDestroy(xvec);
      
  /* return success */
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuMultiply
 * -----------------------------------------------------------------
 * VPrecDuMultiply performs the matrix-vector product P x = b, with 
 * the preconditioner matrix P generated by VPrecSetup and 
 * multiplied using the HYPRE library.
 *
 * The parameters of VPrecDuMultiply used here are as follows:
 *
 * xx      is the x vector on input
 * bb      is the b vector on output
 * tmp     is a temporary vector the size of xx
 * pdata   is the pre-computed Du preconditioner data
 *
 * The value returned by this VPrecDuMultiply function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried),
 *  -1  for a non-recoverable error.
 * -------------------------------------------------------------- */
int VPrecDuMultiply(double *xx, double *bb, double *tmp, void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;

  /* local variables */
  int ilower[3], iupper[3];
  long int ix, iy, iz, iv, idx, Nx, NGx, Ny, NGy, Nz, NGz;
  long int Vbl, Ybl, Zbl;
  HYPRE_SStructVector bvec, xvec;
  double val;

  /* check that Du matrix initialized */
  if (pdata->DuInit == 0) {
    printf("VPrecDuMultiply error: Du matrix uninitialized!\n");
    return(1);
  }

  /* set local variables */
  Nx  = pdata->Nx;
  NGx = pdata->NGx;
  Ny  = pdata->Ny;
  NGy = pdata->NGy;
  Nz  = pdata->Nz;
  NGz = pdata->NGz;

  /* create the SStruct vectors and set init flags to 1 */
  HYPRE_SStructVectorCreate(pdata->comm, pdata->grid, &bvec);
  HYPRE_SStructVectorCreate(pdata->comm, pdata->grid, &xvec);

  /* set vector storage type */
  HYPRE_SStructVectorSetObjectType(bvec, pdata->mattype);
  HYPRE_SStructVectorSetObjectType(xvec, pdata->mattype);
    
  /* initialize vectors */
  HYPRE_SStructVectorInitialize(bvec);
  HYPRE_SStructVectorInitialize(xvec);
  
  /* convert product, result vectors to HYPRE format        */
  /*    insert product vector entries into HYPRE vector x   */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  for (iv=1; iv<=3; iv++) {
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
    idx = 0;
    for (iz=NGz; iz<Nz+NGz; iz++) {
      Zbl = iz * (Nx+2*NGx) * (Ny+2*NGy);
      for (iy=NGy; iy<Ny+NGy; iy++) {
	Ybl = iy * (Nx+2*NGx);
	for (ix=NGx; ix<Nx+NGx; ix++)
	  tmp[idx++] = xx[Vbl+Zbl+Ybl+ix];
      }
    }
    HYPRE_SStructVectorSetBoxValues(xvec, 0, ilower, iupper, iv-1, tmp);
  }

  /*    assemble vectors */
  HYPRE_SStructVectorAssemble(xvec);
  HYPRE_SStructVectorAssemble(bvec);

  /* computing the matvec */
  hypre_SStructMatvec(1.0, pdata->Du, xvec, 1.0, bvec);

  /* gather the solution vector before extracting values */
  HYPRE_SStructVectorGather(bvec);

  /* extract solution vector into bb */
  ilower[0] = pdata->iXL;  ilower[1] = pdata->iYL;  ilower[2] = pdata->iZL;
  iupper[0] = pdata->iXR;  iupper[1] = pdata->iYR;  iupper[2] = pdata->iZR;
  for (iv=1; iv<=3; iv++) {
    Vbl = iv * (Nx+2*NGx) * (Ny+2*NGy) * (Nz+2*NGz);
    HYPRE_SStructVectorGetBoxValues(bvec, 0, ilower, iupper, iv-1, tmp);
    idx = 0;
    for (iz=NGz; iz<Nz+NGz; iz++) {
      Zbl = iz * (Nx+2*NGx) * (Ny+2*NGy);
      for (iy=NGy; iy<Ny+NGy; iy++) {
	Ybl = iy * (Nx+2*NGx);
	for (ix=NGx; ix<Nx+NGx; ix++) 
	  bb[Vbl+Zbl+Ybl+ix] = tmp[idx++];
      }
    }
  }

  /* destroy old vectors */
  HYPRE_SStructVectorDestroy(bvec);
  HYPRE_SStructVectorDestroy(xvec);
      
  /* return success.  */
  return(0);
}



/* -----------------------------------------------------------------
 * Function : VPrecDuNumIters
 * -----------------------------------------------------------------
 * VPrecDuNumIters returns the current cumulative number of 
 * Multigrid Iterations used by the HYPRE preconditioner for 
 * the Du solves.
 * -------------------------------------------------------------- */
int VPrecDuNumIters(void *P_data)
{
  /* recast P_data as the correct structure */
  ViscPrecDuData pdata;
  pdata = (ViscPrecDuData) P_data;
  
  /* return with the result from pdata */
  return(pdata->totIters);
}


/********************************************************************/
