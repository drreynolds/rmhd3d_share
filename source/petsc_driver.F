c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c=================================================================
#define PETSC_AVOID_MPIF_H

      subroutine PetscDriver
c-----------------------------------------------------------------
c     Description: main MHD driver, handles initialization of 
c        variables, timestepping and data I/O.  This MHD driver
c        routine solves an implicit formulation of 2D, 2.5D, 
c        3D ideal, nonlinearly-resistive, MHD equations.  All 
c        initial setup and physics calculations in this 
c        formulation rely on a previous explicit formulation of 
c        these equations by Ravi Samtaney, Princeton Plasma 
c        Physics Laboratory.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use boundary_conds
      use iounits
      use tvdcoeffs
      use sundials_fcmix
      use profiling
      use petsc_data
      use rmhd3_petsc
      use mpistuff
      
c======= Declarations =========
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscdm.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscts.h>
#include <finclude/petscvec.h90>

c     local driver variables
      TS ts
      Vec uvec
      double precision, pointer :: uarr(:)
      double precision :: dt, ttot, ftime
      integer          :: istart, maxiter, lastiter
      integer          :: new
      integer          :: tvdFlag
      integer          :: timestep

#ifdef DYNAMIC
      double precision, allocatable :: abstol(:,:,:,:)
#else
      double precision :: abstol(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif

c     declare internal namelists for problem specification
      namelist /inparams/ maxiter, new, output_flag, binary_flag
      namelist /dumps/ ndump, ndiag
      namelist /gridbounds/ xl, xr, yl, yr, zl, zr


c======= Internals ============

c     see if this is the output node
      outnode = .false.
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1)  outnode = .true.
      
c     if necessary, allocate local data arrays
#ifdef DYNAMIC
      allocate(abstol(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#endif

c     input general MHD simulation information using namelists
      open(16,file='mhd.inp',form='formatted')
      read(16,inparams)
      read(16,dumps)
      read(16,gridbounds) 
      close(16)

c     root node outputs problem parameters
      if(outnode) then
         write(6,'(A,A,i6,3(A,i2))') '\nInput parameters:',
     &        '\n        maxiter =',maxiter,
     &        '\n            new =',new,
     &        '\n    output_flag =',output_flag,
     &        '\n    binary_flag =',binary_flag
         write(6,'(A,2(A,i4))') '\nDump intervals:',
     &        '\n    ndump =',ndump,
     &        '\n    ndiag =',ndiag
         write(6,'(A,6(A,es12.5))') '\nGrid bounds:',
     &        '\n    xl =',xl,',  xr =',xr,
     &        '\n    yl =',yl,',  yr =',yr,
     &        '\n    zl =',zl,',  zr =',zr
         write(6,'(A)') '\nBoundary conditions:'
         if (xbc == BCperiodic)  then
            write(6,'(A)') '    x: periodic'
         else if (xbc == BCreflecting) then
            write(6,'(A)') '    x: reflecting'
         else 
            write(6,'(A)') '    x: zero-gradient'
         endif
         if (ybc == BCperiodic) then
            write(6,'(A)') '    y: periodic'
         else if (ybc == BCreflecting) then
            write(6,'(A)') '    y: reflecting'
         else 
            write(6,'(A)') '    y: zero-gradient'
         endif
         if (zbc == BCperiodic) then
            write(6,'(A)') '    z: periodic'
         else if (zbc == BCreflecting) then
            write(6,'(A)') '    z: reflecting'
         else 
            write(6,'(A)') '    z: zero-gradient'
         endif
         write(6,'(A,3(A,i4),2(A,i2))') '\nGrid indices:',
     &        '\n    nx =',nx,',  ny =',ny,',  nz =',nz,
     &        ',  nvar =',nvar,',  nghost =',nghost
         write(6,'(A,3(A,i4))') '\nParallelism information:',
     &        '\n    xprocs =',xprocs,',  yprocs =',yprocs,
     &        ',  zprocs =',zprocs
      endif

c     initialize timestepping and data output information
      lastiter = 0
      istart   = 1

c     set up domain, grid
      call SetupDomain
      if(outnode) then
         write(6,*) '\nAfter domain setup...'
         write(6,'(A,3(A,es12.5))') '\nMesh spacing:',
     &        '\n   dx =',dx,',  dy =',dy,',  dz =',dz
      end if

c     determine local domain location in general domain
      call SetupLocalDomain
      if(outnode) then
         write(6,*) '\nAfter local domain setup...'
      endif
      write(6,'(A,i4,A,3(2(es8.1,1x),A))') 
     &     '\nLocal subdomain, p',my_id,' = [',XLloc,XRloc,
     &     '] x [',YLloc,YRloc,'] x [',ZLloc,ZRloc,']'

c     initialize PETSc solver structures
      if(outnode) 
     &     write(6,*) 'Initializing PETSc solver...'
      call RMHD3PetscInitialize
      
c     set up TVD coefficients
c      tvdFlag=1
c      if(tvdFlag.eq.1) call SetTVDCoeffs

c     initialize state variables
      call Initialize(new,ulocal,phi,dt,ttot,istart,lastiter)
      if(outnode) 
     &     write(6,*) '\nAfter local dataspace initialization...'

c     communicate initial boundary data among processors
      if(outnode) 
     &     write(6,*) 'Communicating initial condition to neighbors...'
      call CommState(ulocal)

c     create vector to hold global state
      call VecCreate(comm3d,uvec,ierr)
      CHKERRQ(ierr)
      call VecSetSizes(uvec,nxlocal*nylocal*nzlocal*nvar,
     &     nx*ny*nz*nvar,ierr)
      CHKERRQ(ierr)
      call VecSetBlockSize(uvec,nvar,ierr)
      CHKERRQ(ierr)
      call VecSetUp(uvec,ierr)
      CHKERRQ(ierr)

c     set initial state
      call VecGetArrayF90(uvec,uarr,ierr)
      CHKERRQ(ierr)
      call LocalToGlobal(ulocal,uarr)
      call VecRestoreArrayF90(uvec,uarr,ierr)
      CHKERRQ(ierr)

c     create time stepping context
      call TSCreate(comm3d,ts,ierr)
      CHKERRQ(ierr)
      call TSSetType(ts,TSROSW,ierr)
      CHKERRQ(ierr)
      call TSSetRHSFunction(ts,PETSC_NULL_OBJECT,TSRHSFunction_RMHD3,
     &     PETSC_NULL_OBJECT,ierr)
      CHKERRQ(ierr)

c     just a default, actual final time determined by -ts_final_time
      ftime = 1.0
      call TSSetDuration(ts,maxiter,ftime,ierr)
      CHKERRQ(ierr)

      call TSSetSolution(ts,uvec,ierr)
      CHKERRQ(ierr)
      call TSMonitorSet(ts,TSMonitor_RMHD3,PETSC_NULL_OBJECT,
     &     PETSC_NULL_FUNCTION,ierr)
      CHKERRQ(ierr)

c     take options from options file or PETSc command line
      call TSSetFromOptions(ts,ierr)
      CHKERRQ(ierr)
c     view TS configuration before starting stepping
      call TSView(ts,PETSC_NULL_OBJECT,ierr)
      CHKERRQ(ierr)

      if(outnode)  write(6,*) 
     &     '\nAfter PETSc initialization, starting time steps'

c     start MPI timer
      call prof_start(1)

c     integrate from uvec, returns actual final time in ftime
      call TSSolve(ts,uvec,ftime,ierr)
      CHKERRQ(ierr)

c     stop MPI timer
      call prof_stop(1)

c     output final solution
      call prof_start(6)
      if(output_flag==1) then
         if(binary_flag==0) then
            call WriteGnuplotFile(ulocal,timeStep)
         else
            call WriteBinaryFileParallel(ulocal,timeStep)
         endif
      endif
      call prof_stop(6)

c     output final solver diagnostics
      if (outnode) write(6,'(A,A,11(A,es12.5),A)') 
     &     '\n  ---------------------------------------',
     &     '\n  Runtime statistics:\n',
     &     '\n    Total Solution Time     =',prof_time(1),
     &     '\n    Total RHS Time          =',prof_time(2),
     &     '\n    Total Comm Time         =',prof_time(3),
     &     '\n    Total Pcond Setup Time  =',prof_time(4),
     &     '\n    Total Pcond Solve Time  =',prof_time(5),
     &     '\n    Total Diags & I/O Time  =',prof_time(6),
     &     '\n    Avg RHS Time            =',prof_time(2)/prof_count(2),
     &     '\n    Avg Comm Time           =',prof_time(3)/prof_count(3),
     &     '\n    Avg Pcond Setup Time    =',prof_time(4)/prof_count(4),
     &     '\n    Avg Pcond Solve Time    =',prof_time(5)/prof_count(5),
     &     '\n    Avg Diags & I/O Time    =',prof_time(6)/prof_count(6),
     &     '\n  ---------------------------------------\n'

c     dump all remaining data, output diagnostics, close I/O files
      call WrapUp(ulocal,phi,dt,ttot,timestep)

c     free up allocated memory
      call TSDestroy(ts,ierr)
      CHKERRQ(ierr)
      call VecDestroy(uvec,ierr)
      CHKERRQ(ierr)
      call RMHD3PetscFinalize
#ifdef DYNAMIC
      deallocate(abstol)
#endif


      return
      end subroutine PetscDriver
c=================================================================




      subroutine MHDMain
c-----------------------------------------------------------------
c     Description: main MHD routine, initializes more local grid
c        variables and allocates directional meshpoints.
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common

c======= Declarations =========
      implicit none

c======= Internals ============

c     mesh variables for a dyanmic grid
#ifdef DYNAMIC
      ixlo=1-nghost
      iylo=1-nghost
      izlo=1-nghost
      ixhi=nxlocal+nghost
      iyhi=nylocal+nghost
      izhi=nzlocal+nghost
      
c     mesh variable adjustments for 2-D, dynamic grid
#ifdef TWO_D
      izlo=1; izhi=nzlocal
#endif

c     mesh variables for general dynamic grid
      inlo=min(ixlo,iylo,izlo)
      inhi=max(ixhi,iyhi,izhi)
      nxlsize=nxlocal
      nylsize=nylocal
      nzlsize=nzlocal
#endif

c     Allocate directional meshpoints
      allocate(xc(ixlo:ixhi))
      allocate(yc(iylo:iyhi))
      allocate(zc(izlo:izhi))

c     call the MHD driver routine
      call PetscDriver

c     free space used by directional meshpoints
      deallocate(xc,yc,zc)

      return
      end subroutine MHDMain
c=================================================================
