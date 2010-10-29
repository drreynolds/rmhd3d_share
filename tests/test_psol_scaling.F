c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------
c     $Log: impRMHD_driver.F,v $
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
      call TestDriver

c     free space used by directional meshpoints
      deallocate(xc,yc,zc)

      return
      end subroutine MHDMain
c=================================================================





      subroutine TestDriver
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
      use fastwaveprec_data
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      implicit none

#ifdef DYNAMIC
      double precision, dimension(:,:,:,:), allocatable :: 
     &                  ux, fu, vv, tmp, utrue
      double precision, allocatable :: phi(:,:,:)
#else
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) :: 
     &                  ux, fu, vv, tmp, utrue
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
#endif
      double precision :: ttot, dt, dtexp, dtfac, tstop
      integer          :: istart, maxiter, lastiter
      integer          :: output_flag, binary_flag
      integer          :: ndump, ndiag, new
      integer          :: timestep
      logical          :: outnode


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
      allocate(ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(phi(ixlo:ixhi,iylo:iyhi,izlo:izhi))
      allocate(fu(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(tmp(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(utrue(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
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
         write(6,'(A,3(i3,1x),A,3(2(es12.5,1x),A))') 
     &        '\nLocal subdomain for iproc (',
     &        iprocx,iprocy,iprocz,') is the box [',XLloc,XRloc,
     &        '] x [',YLloc,YRloc,'] x [',ZLloc,ZRloc,']'
      endif

c     initialize state variables, Root node outputs progress
      call Initialize(new,ux,phi,dtexp,ttot,istart,lastiter)
      if(outnode)  write(6,*) '\nFinished local initialization'
      
c     communicate initial boundary data among processors
      if(outnode) 
     &     write(6,*) 'Communicating initial condition to neighbors...'
c     do the communications
      call mesh_update_bdry_async(ux, nvar)
      call SetBoundaryValues(ux)

c     wait for all processes to catch up
#ifdef PARALLEL
      call mpi_barrier(comm3d, ierr)
#endif

c     set up preconditioner system and data

c        get timestep size as a scaled explicit time step
      call newdt(dtexp,ux)
      dt = dtfac*dtexp

c        set fastwave preconditioner options
      NWaves = 1
      ordering = 0
      ADIcorr = .false.

c        set up preconditioner data
      if(outnode) 
     &     write(6,*) '\n Setting up preconditioner data'
      call fastwave_pset(ux, fu, dt, vv, tmp, ierr)


c     call directionally-split preconditioner
      if(outnode) 
     &     write(6,*) '\nRunning the new preconditioner'
      call mpi_barrier(comm3d, ierr)
      call prof_start(2)
      call fastwave_psol(ux, fu, dtexp, vv, tmp, ierr)
      call prof_stop(2)
      if (outnode) write(6,'(4(A,es12.5))') 
     &     '   x-directional solve time =',prof_time(20),
     &     '\n   y-directional solve time =',prof_time(21),
     &     '\n   z-directional solve time =',prof_time(22),
     &     '\n   Total Solution Time =',prof_time(2)

#ifdef DYNAMIC
      deallocate(ux,phi,fu,vv,tmp,utrue)
#endif      


      return
      end subroutine TestDriver
c=================================================================
