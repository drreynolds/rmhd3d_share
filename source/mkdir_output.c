/* -----------------------------------------------------------------
 * Daniel R. Reynolds
 * SMU Mathematics
 * Copyright 2010; all rights reserved
 * -----------------------------------------------------------------
 * Fortran interface routine to create directory for simulation 
 * output (in parallel). 
 *
 *    Usage:  
 *       call MKDIR_OUTPUT(index, myid, info)
 *    Input Arguments:
 *       index -- [int] number for output
 *        myid -- [int] MPI task ID for this process
 *    Output Arguments:
 *        info -- [int] flag denoting success (0) or failure (1)
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* headers */
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
  
/* Definitions of interface function names */
#if defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)
#define MAKEOUTPUTDIR mkdir_output

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)
#define MAKEOUTPUTDIR MKDIR_OUTPUT

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)
#define MAKEOUTPUTDIR mkdir_output_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)
#define MAKEOUTPUTDIR MKDIR_OUTPUT_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)
#define MAKEOUTPUTDIR mkdir_output__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)
#define MAKEOUTPUTDIR MKDIR_OUTPUT__

#else
#define MAKEOUTPUTDIR mkdir_output_

#endif

/* prototype */
int make_output_directory(int index, int myid);


/* Fortran-callable C interface function */
void MAKEOUTPUTDIR(int *index, int *myid, int *info) 
{
  int ierr;

  /* call C routine to create directory */
  ierr = make_output_directory(*index, *myid);

  /* set return value to match ierr */
  *info = ierr;

  return;
}
  

/* Fortran-callable C interface function */
int make_output_directory(int index, int myid) 
{
  int ierr;

  /* root process creates a directory for this set of outputs */
  char dirname[100];
  sprintf(dirname, "data%06d", index);
  struct stat st;
  if (myid == 0) {
    if(stat(dirname,&st) != 0) {
      ierr = mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (ierr != 0) {
	fprintf(stderr, "Could not create output directory %s\n",dirname);
	return 1;
      }
    }
  /* other processes wait for this directory to exist before proceeding */
  } else {
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 100;  /* pause 100 nanoseconds between checks */
    while(stat(dirname,&st) != 0)  nanosleep(&ts, NULL);
  }

  return 0;
}
  

#ifdef __cplusplus /* wrapper to enable C++ usage */
}
#endif
/********************************************************************/
