/* -----------------------------------------------------------------
 * Daniel R. Reynolds
 * SMU Mathematics
 * Copyright 2010; all rights reserved
 * -----------------------------------------------------------------
 * Fortran interface routine to output simulation data to SILO 
 * output files (in parallel). 
 *
 *    Usage:  
 *       call WRITE_SILO(vx, div, jcurrent, index, xc, yc, zc, nx, 
 *                       ny, nz, nvar, nxtot, nytot, nztot, xoff, 
 *                       yoff, zoff, dx, dy, dz, iprocx, iprocy, 
 *                       iprocz, nodal, myid, nprocs, info)
 *    Input Arguments:
 *          vx -- [double] array of primitive variable solution data
 *         div -- [double] array of div(B) values
 *    jcurrent -- [double] array of z-directional current (Jz)
 *       index -- [int] number for output
 *          *c -- [double] arrays of cell centers in each direction
 *          n* -- [int] active size of grid in each direction
 *        nvar -- [int] number of variables per spatial cell
 *       n*tot -- [int] total array size in each direction
 *        *off -- [int] offset to beginning of active data in each dir
 *          d* -- [double] spatial mesh size in each direction
 *      iproc* -- [int] location of this process in Cartesian grid
 *       nodal -- [int] flag to use nodal output (1) or FV output (0)
 *        myid -- [int] task ID for this process
 *      nprocs -- [int] total number of processors
 *    Output Arguments:
 *        info -- [int] flag denoting success (0) or failure (1)
 * -----------------------------------------------------------------
 */

#ifdef SILO

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* headers */
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <silo.h>
#include <sys/stat.h>
  
/* Definitions of interface function names */
#if defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)
#define SILOWRITE write_silo

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)
#define SILOWRITE WRITE_SILO

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)
#define SILOWRITE write_silo_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)
#define SILOWRITE WRITE_SILO_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)
#define SILOWRITE write_silo__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)
#define SILOWRITE WRITE_SILO__

#else
#define SILOWRITE write_silo_

#endif

/* prototypes */
int write_my_silo(double *vx, double *div, double *jcurrent, int index, 
		  double *xc, double *yc, double *zc, int nx, int ny, 
		  int nz, int nvar, int nxtot, int nytot, int nztot, 
		  int xoff, int yoff, int zoff, double dx, 
		  double dy, double dz, int iprocx, int iprocy, 
		  int iprocz, int nodal, int myid);
int write_multimesh(int index, int myid, int nprocs, int nvar);


/* Fortran-callable C interface function */
void SILOWRITE(double *vx, double *div, double *jcurrent, int *index, 
	       double *xc, double *yc, double *zc, int *nx, int *ny, 
	       int *nz, int *nv, int *nxtot, int *nytot, int *nztot, 
	       int *xoff, int *yoff, int *zoff, double *dx, double *dy, 
	       double *dz, int *iprocx, int *iprocy, int *iprocz, 
	       int *nodal, int *myid, int *nprocs, int *info) 
{
  int ierr;

  /* initialze return value to success */
  *info = 0;
  
  /* have each process write it's own silo output file */
  ierr = write_my_silo(vx, div, jcurrent, *index, xc, yc, zc, *nx, 
		       *ny, *nz, *nv, *nxtot, *nytot, *nztot, *xoff, 
		       *yoff, *zoff, *dx, *dy, *dz, *iprocx, *iprocy, 
		       *iprocz, *nodal, *myid);
  if (ierr != 0) {
    fprintf(stderr, "p%i: Error writing silo output file: %i\n",*myid,ierr);
    *info = 1;
    return;
  }

  /* root process then outputs multimesh and multivar file */
  if (*myid == 0) {
    ierr = write_multimesh(*index, *myid, *nprocs, *nv);
    if (ierr != 0) {
      fprintf(stderr, "Error writing multimesh file: %i\n",ierr);
      *info = 1;
      return;
    }
  }
  return;
}
  
/* C internals */
int write_my_silo(double *vx, double *div, double *jcurrent, int index, 
		  double *xc, double *yc, double *zc, int nx, int ny, 
		  int nz, int nvar, int nxtot, int nytot, int nztot, 
		  int xoff, int yoff, int zoff, double dx, 
		  double dy, double dz, int iprocx, int iprocy, 
		  int iprocz, int nodal, int myid) 
{
  int ierr;

  /* determine problem dimensionality */
  int ndims = 3;
  if (nztot == 1)  ndims -= 1;
  if (nytot == 1)  ndims -= 1;
  if (nxtot == 1) {
    fprintf(stderr, "Write_my_silo error: proc %i has no x-dir. extents!\n",myid);
    return 1;
  }

  /* File name for data */
  char outfile[100];
  sprintf(outfile, "data%06d/cpu.%06d.silo", index, myid);

  /* create output file */
  DBfile *dbfile = NULL;
  dbfile = DBCreate(outfile, DB_CLOBBER, DB_LOCAL, "Box Data", DB_HDF5);
  if (dbfile == NULL) {
    fprintf(stderr, "Proc %i could not create silo file %s\n",myid,outfile);
    return 1;
  }

  /* create an option list to save cycle value */
  DBoptlist *optlist = DBMakeOptlist(2);
  ierr = DBAddOption(optlist, DBOPT_CYCLE, &index);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBAddOption = %i\n",myid,ierr);
    return ierr;
  }

  /* local data for writing file output */
  double *coords[3];
  int dims[3], idx, i, j, k, Vbuf, Zbuf, Ybuf;
  double *tmp = (double *) malloc(nx*ny*nz*sizeof(double));
  double *xn, *yn, *zn;

  /***** finite-difference plots (interpolate solution between values) *****/
  if (nodal) {

    /* for a nodal mesh, extend lower subdomain face to connect with neighbor */
    if (iprocx > 1) {
      nx += 1;
      xoff -= 1;
    }
    if ((iprocy > 1) && (ndims > 1)) {
      ny += 1;
      yoff -= 1;
    }
    if ((iprocz > 1) && (ndims > 2)) {
      nz += 1;
      zoff -= 1;
    }

    /* create rectilinear mesh */
    xn = (double *) malloc(nx*sizeof(double));
    for (i=0; i<nx; i++)  xn[i] = xc[i+xoff];
    yn = (double *) malloc(ny*sizeof(double));
    for (i=0; i<ny; i++)  yn[i] = yc[i+yoff];
    zn = (double *) malloc(nz*sizeof(double));
    for (i=0; i<nz; i++)  zn[i] = zc[i+zoff];

    /* store pointers for rectilinear mesh */
    coords[0] = xn;
    coords[1] = yn;
    coords[2] = zn;

    /* set local mesh dims */
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;

    /* write rectilinear mesh */
    ierr = DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, 
			 ndims, DB_DOUBLE, DB_COLLINEAR, optlist);
    if (ierr != 0) {
      fprintf(stderr, "Proc %i: error in DBPutQuadmesh = %i\n",myid,ierr);
      return ierr;
    }

  /***** finite-volume plots (constant values in each cell) *****/
  } else {

    /* create nodal mesh (wireframe containing all simulation cells) */
    xn = (double *) malloc((nx+1)*sizeof(double));
    for (i=0; i<=nx; i++)  xn[i] = xc[i+xoff]-0.5*dx;
    yn = (double *) malloc((ny+1)*sizeof(double));
    for (i=0; i<=ny; i++)  yn[i] = yc[i+yoff]-0.5*dy;
    zn = (double *) malloc((nz+1)*sizeof(double));
    for (i=0; i<=nz; i++)  zn[i] = zc[i+zoff]-0.5*dz;

    /* store pointers for rectilinear mesh */
    coords[0] = xn;
    coords[1] = yn;
    coords[2] = zn;

    /* set local mesh dims */
    dims[0] = nx+1;
    dims[1] = ny+1;
    dims[2] = nz+1;

    /* write rectilinear mesh */
    ierr = DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, 
			 ndims, DB_DOUBLE, DB_COLLINEAR, optlist);
    if (ierr != 0) {
      fprintf(stderr, "Proc %i: error in DBPutQuadmesh = %i\n",myid,ierr);
      return ierr;
    }

    /* reset local mesh dims for cell-centered variables */
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;

  } /* if (nodal) */



  /* write each variable to silo file */

  /**** rho ****/
  /* store values in temporary array */
  Vbuf = 0*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "rho", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "rho", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }
  
  /**** u ****/
  /* store values in temporary array */
  Vbuf = 1*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "u", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "u", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** v ****/
  /* store values in temporary array */
  Vbuf = 2*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "v", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "v", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** w ****/
  /* store values in temporary array */
  Vbuf = 3*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "w", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "w", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** bx ****/
  /* store values in temporary array */
  Vbuf = 4*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "bx", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "bx", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** by ****/
  /* store values in temporary array */
  Vbuf = 5*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "by", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "by", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** bz ****/
  /* store values in temporary array */
  Vbuf = 6*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "bz", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "bz", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** p ****/
  /* store values in temporary array */
  Vbuf = 7*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = vx[Vbuf+Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "p", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "p", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** divB ****/
  /* store values in temporary array */
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = div[Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "divB", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "divB", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** Jz ****/
  /* store values in temporary array */
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) 
	tmp[idx++] = jcurrent[Zbuf+Ybuf+i];
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "Jz", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "Jz", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }

  /**** ptot ****/
  /* store values in temporary array */
  double press, bx, by, bz;
  Vbuf = 0*nxtot*nytot*nztot;
  idx = 0;
  for (k=zoff; k<nz+zoff; k++) {
    Zbuf = k*nxtot*nytot;
    for (j=yoff; j<ny+yoff; j++) {
      Ybuf = j*nxtot;
      for (i=xoff; i<nx+xoff; i++) {
	bx    = vx[4*nxtot*nytot*nztot+Zbuf+Ybuf+i];
	by    = vx[5*nxtot*nytot*nztot+Zbuf+Ybuf+i];
	bz    = vx[6*nxtot*nytot*nztot+Zbuf+Ybuf+i];
	press = vx[7*nxtot*nytot*nztot+Zbuf+Ybuf+i];
	tmp[idx++] = press+0.5*(bx*bx + by*by + bz*bz);
      }
    }
  }
  /* write values to silo file */
  if (nodal) {
    ierr = DBPutQuadvar1(dbfile, "ptot", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
  } else {
    ierr = DBPutQuadvar1(dbfile, "ptot", "quadmesh", (float*) tmp, dims, 
			 ndims, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
  }
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutQuadvar1 = %i\n",myid,ierr);
    return ierr;
  }


  /* free the option list */
  ierr = DBFreeOptlist(optlist);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBFreeOptlist = %i\n",myid,ierr);
    return ierr;
  }

  /* free rectilinear mesh arrays */
  free(xn);
  free(yn);
  free(zn);

  /* free the temporary data array */
  free(tmp);

  /* close silo file */
  ierr = DBClose(dbfile);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBClose = %i\n",myid,ierr);
    return ierr;
  }

  return 0;
};

  
int write_multimesh(int index, int myid, int nprocs, int nvar) 
{

  /* only root outputs the multimesh */
  if (myid != 0)  return 0;

  /* allocate arrays for meshnames, and data types */
  char **meshnames = (char **) malloc(nprocs*sizeof(char *));
  int *meshtypes = (int *) malloc(nprocs*sizeof(int));

  /* fill in these arrays */
  int id, ierr;
  char tmpchar[100];
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's output file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:quadmesh", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUAD_RECT;
  }

  /* multimesh file name */
  char outfile[100];
  sprintf(outfile, "multimesh%06d.root", index);

  /* create an option list to save cycle and time values */
  DBoptlist *optlist = DBMakeOptlist(2);
  ierr = DBAddOption(optlist, DBOPT_CYCLE, &index);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBAddOption = %i\n",myid,ierr);
    return ierr;
  }

  /* create multimesh file */
  DBfile *dbfile = NULL;
  dbfile = DBCreate(outfile, DB_CLOBBER, DB_LOCAL, "multimesh root", DB_HDF5);
  if (dbfile == NULL) {
    fprintf(stderr, "Could not create multimesh file %s\n",outfile);
    return 1;
  }

  /* write the multimesh object */
  ierr = DBPutMultimesh(dbfile, "quadmesh", nprocs, meshnames, meshtypes, optlist);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultimesh = %i\n",myid,ierr);
    return ierr;
  }

  /* free the option list */
  ierr = DBFreeOptlist(optlist);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBFreeOptlist = %i\n",myid,ierr);
    return ierr;
  }


  /* now add the variables */

  /**** rho ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:rho", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "rho", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** u ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:u", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "u", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** v ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:v", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "v", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** w ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:w", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "w", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** bx ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:bx", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "bx", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** by ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:by", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "by", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** bz ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:bz", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "bz", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** p ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:p", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "p", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** divB ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:divB", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "divB", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** Jz ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:Jz", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "Jz", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }

  /**** ptot ****/
  /* fill in descriptor arrays */
  for (id=0; id<nprocs; id++) {
    /* set name for this processor's variable & file */
    sprintf(tmpchar, "data%06d/cpu.%06d.silo:ptot", index, id);
    meshnames[id] = strdup(tmpchar);
    meshtypes[id] = DB_QUADVAR;
  }
  ierr = DBPutMultivar(dbfile, "ptot", nprocs, meshnames, meshtypes, NULL);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBPutMultivar = %i\n",myid,ierr);
    return ierr;
  }


  /* close silo file */
  ierr = DBClose(dbfile);
  if (ierr != 0) {
    fprintf(stderr, "Proc %i: error in DBClose = %i\n",myid,ierr);
    return ierr;
  }

  /* free the allocated arrays */
  free(meshnames);
  free(meshtypes);

  return 0;
};


#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif
#endif
/********************************************************************/
