/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Daniel R. Reynolds
 * UC San Diego, Mathematics
 * -----------------------------------------------------------------
 * This is the header file for viscous_prec.c.
 * -----------------------------------------------------------------
 */

/*
 * ===========================================================================
 *
 *         Description and Usage of the VISCOUS_PREC Interface Package
 *         [THESE NEED TO BE UPDATED TO DESCRIBE NEW INTERNAL FUNCTIONS]
 *
 * ===========================================================================
 */

#ifndef _VISCOUS_PREC_H
#define _VISCOUS_PREC_H

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* HYPRE header files  */
#include <stdio.h>
#include <math.h>
#include "HYPRE_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"

/* MPI header file */
#ifdef PARALLEL
#include "mpi.h"
#else
typedef int MPI_Comm;
#endif


/* Definitions of interface function names */

#if defined(F77_FUNC)

#define VISCPREC_DU_INIT       F77_FUNC(fvprecduinit,     FVPRECDUINIT)
#define VISCPREC_DU_SETUP      F77_FUNC(fvprecdusetup,    FVPRECDUSETUP)
#define VISCPREC_DU_SOLVE      F77_FUNC(fvprecdusolve,    FVPRECDUSOLVE)
#define VISCPREC_DU_MULTIPLY   F77_FUNC(fvprecdumultiply, FVPRECDUMULTIPLY)
#define VISCPREC_DU_FREE       F77_FUNC(fvprecdufree,     FVPRECDUFREE)
#define VISCPREC_DU_NUMITERS   F77_FUNC(fvprecdunumiters, FVPRECDUNUMITERS)
#define SET_SOL_DU_OPTS        F77_FUNC(fvsetsolduopts,   FVSETSOLDUOPTS)

#define VISCPREC_DB_INIT       F77_FUNC(fvprecdbinit,     FVPRECDBINIT)
#define VISCPREC_DB_SETUP      F77_FUNC(fvprecdbsetup,    FVPRECDBSETUP)
#define VISCPREC_DB_SOLVE      F77_FUNC(fvprecdbsolve,    FVPRECDBSOLVE)
#define VISCPREC_DB_MULTIPLY   F77_FUNC(fvprecdbmultiply, FVPRECDBMULTIPLY)
#define VISCPREC_DB_FREE       F77_FUNC(fvprecdbfree,     FVPRECDBFREE)
#define VISCPREC_DB_NUMITERS   F77_FUNC(fvprecdbnumiters, FVPRECDBNUMITERS)
#define SET_SOL_DB_OPTS        F77_FUNC(fvsetsoldbopts,   FVSETSOLDBOPTS)

#define VISCPREC_DE_INIT       F77_FUNC(fvprecdeinit,     FVPRECDEINIT)
#define VISCPREC_DE_SETUP      F77_FUNC(fvprecdesetup,    FVPRECDESETUP)
#define VISCPREC_DE_SOLVE      F77_FUNC(fvprecdesolve,    FVPRECDESOLVE)
#define VISCPREC_DE_MULTIPLY   F77_FUNC(fvprecdemultiply, FVPRECDEMULTIPLY)
#define VISCPREC_DE_FREE       F77_FUNC(fvprecdefree,     FVPRECDEFREE)
#define VISCPREC_DE_NUMITERS   F77_FUNC(fvprecdenumiters, FVPRECDENUMITERS)
#define SET_SOL_DE_OPTS        F77_FUNC(fvsetsoldeopts,   FVSETSOLDEOPTS)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define VISCPREC_DU_INIT       fvprecduinit
#define VISCPREC_DU_SETUP      fvprecdusetup
#define VISCPREC_DU_SOLVE      fvprecdusolve
#define VISCPREC_DU_MULTIPLY   fvprecdumultiply
#define VISCPREC_DU_FREE       fvprecdufree
#define VISCPREC_DU_NUMITERS   fvprecdunumiters
#define SET_SOL_DU_OPTS        fvsetsolduopts

#define VISCPREC_DB_INIT       fvprecdbinit
#define VISCPREC_DB_SETUP      fvprecdbsetup
#define VISCPREC_DB_SOLVE      fvprecdbsolve
#define VISCPREC_DB_MULTIPLY   fvprecdbmultiply
#define VISCPREC_DB_FREE       fvprecdbfree
#define VISCPREC_DB_NUMITERS   fvprecdbnumiters
#define SET_SOL_DB_OPTS        fvsetsoldbopts

#define VISCPREC_DE_INIT       fvprecdeinit
#define VISCPREC_DE_SETUP      fvprecdesetup
#define VISCPREC_DE_SOLVE      fvprecdesolve
#define VISCPREC_DE_MULTIPLY   fvprecdemultiply
#define VISCPREC_DE_FREE       fvprecdefree
#define VISCPREC_DE_NUMITERS   fvprecdenumiters
#define SET_SOL_DE_OPTS        fvsetsoldeopts

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define VISCPREC_DU_INIT       FVPRECDUINIT
#define VISCPREC_DU_SETUP      FVPRECDUSETUP
#define VISCPREC_DU_SOLVE      FVPRECDUSOLVE
#define VISCPREC_DU_MULTIPLY   FVPRECDUMULTIPLY
#define VISCPREC_DU_FREE       FVPRECDUFREE
#define VISCPREC_DU_NUMITERS   FVPRECDUNUMITERS
#define SET_SOL_DU_OPTS        FVSETSOLDUOPTS

#define VISCPREC_DB_INIT       FVPRECDBINIT
#define VISCPREC_DB_SETUP      FVPRECDBSETUP
#define VISCPREC_DB_SOLVE      FVPRECDBSOLVE
#define VISCPREC_DB_MULTIPLY   FVPRECDBMULTIPLY
#define VISCPREC_DB_FREE       FVPRECDBFREE
#define VISCPREC_DB_NUMITERS   FVPRECDBNUMITERS
#define SET_SOL_DB_OPTS        FVSETSOLDBOPTS

#define VISCPREC_DE_INIT       FVPRECDEINIT
#define VISCPREC_DE_SETUP      FVPRECDESETUP
#define VISCPREC_DE_SOLVE      FVPRECDESOLVE
#define VISCPREC_DE_MULTIPLY   FVPRECDEMULTIPLY
#define VISCPREC_DE_FREE       FVPRECDEFREE
#define VISCPREC_DE_NUMITERS   FVPRECDENUMITERS
#define SET_SOL_DE_OPTS        FVSETSOLDEOPTS

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define VISCPREC_DU_INIT       fvprecduinit_
#define VISCPREC_DU_SETUP      fvprecdusetup_
#define VISCPREC_DU_SOLVE      fvprecdusolve_
#define VISCPREC_DU_MULTIPLY   fvprecdumultiply_
#define VISCPREC_DU_FREE       fvprecdufree_
#define VISCPREC_DU_NUMITERS   fvprecdunumiters_
#define SET_SOL_DU_OPTS        fvsetsolduopts_

#define VISCPREC_DB_INIT       fvprecdbinit_
#define VISCPREC_DB_SETUP      fvprecdbsetup_
#define VISCPREC_DB_SOLVE      fvprecdbsolve_
#define VISCPREC_DB_MULTIPLY   fvprecdbmultiply_
#define VISCPREC_DB_FREE       fvprecdbfree_
#define VISCPREC_DB_NUMITERS   fvprecdbnumiters_
#define SET_SOL_DB_OPTS        fvsetsoldbopts_

#define VISCPREC_DE_INIT       fvprecdeinit_
#define VISCPREC_DE_SETUP      fvprecdesetup_
#define VISCPREC_DE_SOLVE      fvprecdesolve_
#define VISCPREC_DE_MULTIPLY   fvprecdemultiply_
#define VISCPREC_DE_FREE       fvprecdefree_
#define VISCPREC_DE_NUMITERS   fvprecdenumiters_
#define SET_SOL_DE_OPTS        fvsetsoldeopts_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define VISCPREC_DU_INIT       FVPRECDUINIT_
#define VISCPREC_DU_SETUP      FVPRECDUSETUP_
#define VISCPREC_DU_SOLVE      FVPRECDUSOLVE_
#define VISCPREC_DU_MULTIPLY   FVPRECDUMULTIPLY_
#define VISCPREC_DU_FREE       FVPRECDUFREE_
#define VISCPREC_DU_NUMITERS   FVPRECDUNUMITERS_
#define SET_SOL_DU_OPTS        FVSETSOLDUOPTS_

#define VISCPREC_DB_INIT       FVPRECDBINIT_
#define VISCPREC_DB_SETUP      FVPRECDBSETUP_
#define VISCPREC_DB_SOLVE      FVPRECDBSOLVE_
#define VISCPREC_DB_MULTIPLY   FVPRECDBMULTIPLY_
#define VISCPREC_DB_FREE       FVPRECDBFREE_
#define VISCPREC_DB_NUMITERS   FVPRECDBNUMITERS_
#define SET_SOL_DB_OPTS        FVSETSOLDBOPTS_

#define VISCPREC_DE_INIT       FVPRECDEINIT_
#define VISCPREC_DE_SETUP      FVPRECDESETUP_
#define VISCPREC_DE_SOLVE      FVPRECDESOLVE_
#define VISCPREC_DE_MULTIPLY   FVPRECDEMULTIPLY_
#define VISCPREC_DE_FREE       FVPRECDEFREE_
#define VISCPREC_DE_NUMITERS   FVPRECDENUMITERS_
#define SET_SOL_DE_OPTS        FVSETSOLDEOPTS_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define VISCPREC_DU_INIT       fvprecduinit__
#define VISCPREC_DU_SETUP      fvprecdusetup__
#define VISCPREC_DU_SOLVE      fvprecdusolve__
#define VISCPREC_DU_MULTIPLY   fvprecdumultiply__
#define VISCPREC_DU_FREE       fvprecdufree__
#define VISCPREC_DU_NUMITERS   fvprecdunumiters__
#define SET_SOL_DU_OPTS        fvsetsolduopts__

#define VISCPREC_DB_INIT       fvprecdbinit__
#define VISCPREC_DB_SETUP      fvprecdbsetup__
#define VISCPREC_DB_SOLVE      fvprecdbsolve__
#define VISCPREC_DB_MULTIPLY   fvprecdbmultiply__
#define VISCPREC_DB_FREE       fvprecdbfree__
#define VISCPREC_DB_NUMITERS   fvprecdbnumiters__
#define SET_SOL_DB_OPTS        fvsetsoldbopts__

#define VISCPREC_DE_INIT       fvprecdeinit__
#define VISCPREC_DE_SETUP      fvprecdesetup__
#define VISCPREC_DE_SOLVE      fvprecdesolve__
#define VISCPREC_DE_MULTIPLY   fvprecdemultiply__
#define VISCPREC_DE_FREE       fvprecdefree__
#define VISCPREC_DE_NUMITERS   fvprecdenumiters__
#define SET_SOL_DE_OPTS        fvsetsoldeopts__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define VISCPREC_DU_INIT       FVPRECDUINIT__
#define VISCPREC_DU_SETUP      FVPRECDUSETUP__
#define VISCPREC_DU_SOLVE      FVPRECDUSOLVE__
#define VISCPREC_DU_MULTIPLY   FVPRECDUMULTIPLY__
#define VISCPREC_DU_FREE       FVPRECDUFREE__
#define VISCPREC_DU_NUMITERS   FVPRECDUNUMITERS__
#define SET_SOL_DU_OPTS        FVSETSOLDUOPTS__

#define VISCPREC_DB_INIT       FVPRECDBINIT__
#define VISCPREC_DB_SETUP      FVPRECDBSETUP__
#define VISCPREC_DB_SOLVE      FVPRECDBSOLVE__
#define VISCPREC_DB_MULTIPLY   FVPRECDBMULTIPLY__
#define VISCPREC_DB_FREE       FVPRECDBFREE__
#define VISCPREC_DB_NUMITERS   FVPRECDBNUMITERS__
#define SET_SOL_DB_OPTS        FVSETSOLDBOPTS__

#define VISCPREC_DE_INIT       FVPRECDEINIT__
#define VISCPREC_DE_SETUP      FVPRECDESETUP__
#define VISCPREC_DE_SOLVE      FVPRECDESOLVE__
#define VISCPREC_DE_MULTIPLY   FVPRECDEMULTIPLY__
#define VISCPREC_DE_FREE       FVPRECDEFREE__
#define VISCPREC_DE_NUMITERS   FVPRECDENUMITERS__
#define SET_SOL_DE_OPTS        FVSETSOLDEOPTS__

#endif



/* Data structure for internal precondtioner data */
typedef struct {

  /* Domain-specific variables */
  int     ndim;    /* dimensional size of grid */
  long int  Nx;    /* mesh cells in x-direction */
  long int  Ny;    /* mesh cells in y-direction */
  long int  Nz;    /* mesh cells in z-direction */
  double    dx;    /* mesh size in x-direction */
  double    dy;    /* mesh size in y-direction */
  double    dz;    /* mesh size in z-direction */
  long int  Ns;    /* number of MHD species */
  long int  NGx;   /* ghost cells in x-direction */
  long int  NGy;   /* ghost cells in y-direction */
  long int  NGz;   /* ghost cells in z-direction */
  long int  iXL;   /* lower x-bound for this processor in global array */
  long int  iXR;   /* upper x-bound for this processor in global array */
  long int  iYL;   /* lower y-bound for this processor in global array */
  long int  iYR;   /* upper y-bound for this processor in global array */
  long int  iZL;   /* lower z-bound for this processor in global array */
  long int  iZR;   /* upper z-bound for this processor in global array */
  int       xbc;   /* x-boundary condition (0->0-grad; 1->perdc; 2->reflect) */
  int       ybc;   /* y-boundary condition (0->0-grad; 1->perdc; 2->reflect) */
  int       zbc;   /* z-boundary condition (0->0-grad; 1->perdc; 2->reflect) */

  /* MPI-specific variables */
  int Xprocs;     /* number of processors in x-direction */
  int Yprocs;     /* number of processors in y-direction */
  int Zprocs;     /* number of processors in z-direction */
  int iprocx;     /* current process location in x-direction */
  int iprocy;     /* current process location in y-direction */
  int iprocz;     /* current process location in z-direction */
  int outproc;    /* integer representing output processor */
  MPI_Comm comm;  /* MPI communicator */

  /* HYPRE SStruct-specific data */
  int  mattype;                     /* HYPRE matrix type for Du solve */
  int  wStSize;                     /* stencil size */
  HYPRE_SStructGrid    grid;        /* HYPRE grid object for Du solve */
  HYPRE_SStructStencil wxStencil;   /* wx stencil object */
  HYPRE_SStructStencil wyStencil;   /* wy stencil object */
  HYPRE_SStructStencil wzStencil;   /* wz stencil object */
  HYPRE_SStructGraph   graph;       /* HYPRE graph object for Du solve */

  /* HYPRE Solver-specific data */
  HYPRE_SStructMatrix Du;      /* Jac. for momentum vars. */
  int sol_zeroguess;           /* use zero initial guess */
  int sol_maxit;               /* max iterations */
  int sol_relch;               /* rel. change stopping criteria */
  int sol_rlxtype;             /* relaxation type */
  int sol_npre;                /* num. pre-relaxation sweeps */
  int sol_npost;               /* num. post-relaxation sweeps */
  int sol_printl;              /* print level */
  int sol_log;                 /* amount of logging */

  /* Extra variables for solver diagnostics */
  int DuInit;    /* flag denoting initialization of Du matrix */
  int totIters;  /* total MG iterations for Du solves */

} *ViscPrecDuData;


/* Data structure for internal precondtioner data */
typedef struct {

  /* Domain-specific variables */
  int     ndim;    /* dimensional size of grid */
  long int  Nx;    /* mesh cells in x-direction */
  long int  Ny;    /* mesh cells in y-direction */
  long int  Nz;    /* mesh cells in z-direction */
  double    dx;    /* mesh size in x-direction */
  double    dy;    /* mesh size in y-direction */
  double    dz;    /* mesh size in z-direction */
  long int  Ns;    /* number of MHD species */
  long int  NGx;   /* ghost cells in x-direction */
  long int  NGy;   /* ghost cells in y-direction */
  long int  NGz;   /* ghost cells in z-direction */
  long int  iXL;   /* lower x-bound for this processor in global array */
  long int  iXR;   /* upper x-bound for this processor in global array */
  long int  iYL;   /* lower y-bound for this processor in global array */
  long int  iYR;   /* upper y-bound for this processor in global array */
  long int  iZL;   /* lower z-bound for this processor in global array */
  long int  iZR;   /* upper z-bound for this processor in global array */
  int       xbc;   /* x-boundary condition (0->0-grad; 1->perdc; 2->reflect) */
  int       ybc;   /* y-boundary condition (0->0-grad; 1->perdc; 2->reflect) */
  int       zbc;   /* z-boundary condition (0->0-grad; 1->perdc; 2->reflect) */

  /* MPI-specific variables */
  int Xprocs;     /* number of processors in x-direction */
  int Yprocs;     /* number of processors in y-direction */
  int Zprocs;     /* number of processors in z-direction */
  int iprocx;     /* current process location in x-direction */
  int iprocy;     /* current process location in y-direction */
  int iprocz;     /* current process location in z-direction */
  int outproc;    /* integer representing output processor */
  MPI_Comm comm;  /* MPI communicator */

  /* HYPRE SStruct-specific data */
  int  mattype;                     /* HYPRE matrix type for Db solve */
  int  bStSize;                     /* bx stencil size */
  HYPRE_SStructGrid    grid;        /* HYPRE grid object for Db solve */
  HYPRE_SStructStencil bxStencil;   /* bx stencil object */
  HYPRE_SStructStencil byStencil;   /* by stencil object */
  HYPRE_SStructStencil bzStencil;   /* bz stencil object */
  HYPRE_SStructGraph   graph;       /* HYPRE graph object for Db solve */

  /* HYPRE Solver-specific data */
  HYPRE_SStructMatrix Db;       /* Jac. for mag. field vars. */
  int sol_zeroguess;           /* use zero initial guess */
  int sol_maxit;               /* max iterations */
  int sol_relch;               /* rel. change stopping criteria */
  int sol_rlxtype;             /* relaxation type */
  int sol_npre;                /* num. pre-relaxation sweeps */
  int sol_npost;               /* num. post-relaxation sweeps */
  int sol_printl;              /* print level */
  int sol_log;                 /* amount of logging */

  /* Extra variables for solver diagnostics */
  int DbInit;    /* flag denoting initialization of Db matrix */
  int totIters;  /* total MG iterations for Db solves */

} *ViscPrecDbData;


/* Data structure for internal precondtioner data */
typedef struct {

  /* Domain-specific variables */
  int     ndim;    /* dimensional size of grid */
  long int  Nx;    /* mesh cells in x-direction */
  long int  Ny;    /* mesh cells in y-direction */
  long int  Nz;    /* mesh cells in z-direction */
  double    dx;    /* mesh size in x-direction */
  double    dy;    /* mesh size in y-direction */
  double    dz;    /* mesh size in z-direction */
  long int  Ns;    /* number of MHD species */
  long int  NGx;   /* ghost cells in x-direction */
  long int  NGy;   /* ghost cells in y-direction */
  long int  NGz;   /* ghost cells in z-direction */
  long int  iXL;   /* lower x-bound for this processor in global array */
  long int  iXR;   /* upper x-bound for this processor in global array */
  long int  iYL;   /* lower y-bound for this processor in global array */
  long int  iYR;   /* upper y-bound for this processor in global array */
  long int  iZL;   /* lower z-bound for this processor in global array */
  long int  iZR;   /* upper z-bound for this processor in global array */
  int       xbc;   /* x-boundary condition (0->0-grad; 1->perdc; 2->reflect) */
  int       ybc;   /* y-boundary condition (0->0-grad; 1->perdc; 2->reflect) */
  int       zbc;   /* z-boundary condition (0->0-grad; 1->perdc; 2->reflect) */

  /* MPI-specific variables */
  int Xprocs;     /* number of processors in x-direction */
  int Yprocs;     /* number of processors in y-direction */
  int Zprocs;     /* number of processors in z-direction */
  int iprocx;     /* current process location in x-direction */
  int iprocy;     /* current process location in y-direction */
  int iprocz;     /* current process location in z-direction */
  int outproc;    /* integer representing output processor */
  MPI_Comm comm;  /* MPI communicator */

  /* HYPRE SStruct-specific data */
  int  mattype;                    /* HYPRE matrix type for De solve */
  int  eStSize;                    /* e stencil size */
  HYPRE_SStructGrid    grid;       /* HYPRE grid object for De solve */
  HYPRE_SStructStencil eStencil;   /* e stencil object */
  HYPRE_SStructGraph   graph;      /* HYPRE graph object for De solve */

  /* HYPRE Solver-specific data */
  HYPRE_SStructMatrix De;      /* Jac. for energy */
  int sol_zeroguess;           /* use zero initial guess */
  int sol_maxit;               /* max iterations */
  int sol_relch;               /* rel. change stopping criteria */
  int sol_rlxtype;             /* relaxation type */
  int sol_npre;                /* num. pre-relaxation sweeps */
  int sol_npost;               /* num. post-relaxation sweeps */
  int sol_printl;              /* print level */
  int sol_log;                 /* amount of logging */

  /* Extra variables for solver diagnostics */
  int DeInit;    /* flag denoting initialization of De matrix */
  int totIters;  /* total MG iterations for De solves */

} *ViscPrecDeData;




/* Prototypes of C Preconditioner Functions */
void *VPrecDuAlloc(long int Nx, long int Ny, long int Nz, double dx, 
		   double dy, double dz, long int Ns, long int NGx, 
		   long int NGy, long int NGz, int NPx, int NPy, 
		   int NPz, int iPx, int iPy, int iPz, int NBlt, 
		   int NBrt, int NBtp, int NBbt, int NBft, int NBbk,
		   int xbcond, int ybcond, int zbcond);
void *VPrecDbAlloc(long int Nx, long int Ny, long int Nz, double dx, 
		   double dy, double dz, long int Ns, long int NGx, 
		   long int NGy, long int NGz, int NPx, int NPy, 
		   int NPz, int iPx, int iPy, int iPz, int NBlt, 
		   int NBrt, int NBtp, int NBbt, int NBft, int NBbk,
		   int xbcond, int ybcond, int zbcond);
void *VPrecDeAlloc(long int Nx, long int Ny, long int Nz, 
		   double dx, double dy, double dz, 
		   long int Ns, long int NGx, long int NGy, 
		   long int NGz, int NPx, int NPy, int NPz, 
		   int iPx, int iPy, int iPz, int NBlt, int NBrt, 
		   int NBtp, int NBbt, int NBft, int NBbk,
		   int xbcond, int ybcond, int zbcond);


void VPrecDuFree(void *P_data);
void VPrecDbFree(void *P_data);
void VPrecDeFree(void *P_data);


int VPrecDuSetup(double *uu, double gamdt, double Mu, 
		 double Re, double *v1, double *v2, void *pdata);
int VPrecDbSetup(double *uu, double gamdt, double Lu, 
		 double Eta, double *v1, double *v2, void *pdata);
int VPrecDeSetup(double *uu, double gamdt, double Gamma, 
		 double Kappa, double Re, double Pr, 
		 double RGas, double *v1, double *v2, void *pdata);


int VPrecDuSolve(double *xx, double *bb, double *tmp, double delta, void *pdata);
int VPrecDbSolve(double *xx, double *bb, double *tmp, double delta, void *pdata);
int VPrecDeSolve(double *xx, double *bb, double *tmp, double delta, void *pdata);


int VPrecDuMultiply(double *xx, double *bb, double *tmp, void *pdata);
int VPrecDbMultiply(double *xx, double *bb, double *tmp, void *pdata);
int VPrecDeMultiply(double *xx, double *bb, double *tmp, void *pdata);


int VPrecDuNumIters(void *P_data);
int VPrecDbNumIters(void *P_data);
int VPrecDeNumIters(void *P_data);



/* Prototypes of Fortran-callable C interface functions */
void VISCPREC_DU_INIT(long int *Nx, long int *Ny, long int *Nz, 
		      double *dx, double *dy, double *dz,
		      long int *Ns, long int *NGx, long int *NGy, 
		      long int *NGz, int *NPx, int *NPy, int *NPz, 
		      int *iPx, int *iPy, int *iPz, int *NBlt, 
		      int *NBrt, int *NBtp, int *NBbt, int *NBft, 
		      int *NBbk, int *XBcond, int *YBcond, 
		      int *ZBcond, int *ier);
void VISCPREC_DB_INIT(long int *Nx, long int *Ny, long int *Nz, 
		      double *dx, double *dy, double *dz,
		      long int *Ns, long int *NGx, long int *NGy, 
		      long int *NGz, int *NPx, int *NPy, int *NPz, 
		      int *iPx, int *iPy, int *iPz, int *NBlt, 
		      int *NBrt, int *NBtp, int *NBbt, int *NBft, 
		      int *NBbk, int *XBcond, int *YBcond, 
		      int *ZBcond, int *ier);
void VISCPREC_DE_INIT(long int *Nx, long int *Ny, long int *Nz, 
		      double *dx, double *dy, double *dz,
		      long int *Ns, long int *NGx, long int *NGy, 
		      long int *NGz, int *NPx, int *NPy, int *NPz, 
		      int *iPx, int *iPy, int *iPz, int *NBlt, 
		      int *NBrt, int *NBtp, int *NBbt, int *NBft, 
		      int *NBbk, int *XBcond, int *YBcond, 
		      int *ZBcond, int *ier);


void VISCPREC_DU_FREE();
void VISCPREC_DB_FREE();
void VISCPREC_DE_FREE();


void VISCPREC_DU_SETUP(double *uu, double *gamdt, double *Mu, 
		       double *Re, double *v1, double *v2, 
		       int *ier);
void VISCPREC_DB_SETUP(double *uu, double *gamdt, double *Lu, 
		       double *Eta, double *v1, double *v2, 
		       int *ier);
void VISCPREC_DE_SETUP(double *uu, double *gamdt, 
		       double *Gamma, double *Kappa, 
		       double *Re, double *Pr, double *RGas, 
		       double *v1, double *v2, int *ier);


void VISCPREC_DU_SOLVE(double *xx, double *bb, double *tmp, 
		       double *delta, int *ier);
void VISCPREC_DB_SOLVE(double *xx, double *bb, double *tmp, 
		       double *delta, int *ier);
void VISCPREC_DE_SOLVE(double *xx, double *bb, double *tmp, 
		       double *delta, int *ier);


void VISCPREC_DU_MULTIPLY(double *xx, double *bb, double *tmp, int *ier);
void VISCPREC_DB_MULTIPLY(double *xx, double *bb, double *tmp, int *ier);
void VISCPREC_DE_MULTIPLY(double *xx, double *bb, double *tmp, int *ier);


void VISCPREC_DU_NUMITERS(int *Niters);
void VISCPREC_DB_NUMITERS(int *Niters);
void VISCPREC_DE_NUMITERS(int *Niters);


void SET_SOL_DU_OPTS(int *iopt, int *ier);
void SET_SOL_DB_OPTS(int *iopt, int *ier);
void SET_SOL_DE_OPTS(int *iopt, int *ier);



#endif
/********************************************************************/
