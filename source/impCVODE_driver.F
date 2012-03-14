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

      subroutine ImpCVODEDriver
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
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      implicit none

c     local driver variables
#ifdef DYNAMIC
      double precision, allocatable :: ux(:,:,:,:)
      double precision, allocatable :: phi(:,:,:)
#else
      double precision :: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
#endif
      double precision :: ttot, told, dt, dthyp, dtvisc, dtexp, tdesired
      integer          :: istart, maxiter, lastiter
      integer          :: output_flag, binary_flag
      integer          :: ndump, ndiag, new
      integer          :: tvdFlag
      integer          :: timestep
      logical          :: outnode

c     cvode timestepper/solver information
      integer          :: taskchoice, cvierr
      double precision :: tstop, dtfac
      double precision :: ropt(40)
      integer(kind=sund_long_int) :: iopt(40)
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
      allocate(ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(phi(ixlo:ixhi,iylo:iyhi,izlo:izhi))
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
      
c     set up TVD coefficients
c      tvdFlag=1
c      if(tvdFlag.eq.1) call SetTVDCoeffs

c     initialize state variables
      call Initialize(new,ux,phi,dt,ttot,istart,lastiter)
      if(outnode) 
     &     write(6,*) '\nAfter local dataspace initialization...'

c     communicate initial boundary data among processors
      if(outnode) 
     &     write(6,*) 'Communicating initial condition to neighbors...'
      call CommState(ux)


c     initialize CVODE solver structures
      if(outnode) 
     &     write(6,*) 'Initializing CVODE solver...'
      call CvodeInit(ttot, ux, abstol, iopt, ropt, 
     &     dtfac, tstop, taskchoice)


c     output initial state and diagnostics (energies, mass, etc)
      call OutputData(ux, phi, 0, output_flag, binary_flag, 
     &     0.d0, ttot)
      call Diagnostics(ux,phi,ttot)

      if(outnode)  write(6,*) 
     &     '\nAfter CVODE initialization, starting time steps'

c     start MPI timer
      call prof_start(1)

c     Start time stepping
      do timeStep = istart,maxiter

c        set next desired output time
         call newdt(ux,dthyp,dtvisc)
         dtexp = min(dthyp,dtvisc)
         told = ttot
         dt = min(dtfac, tstop-ttot)
         tdesired = ttot + dt
         if(outnode) then
            write(6,*) 'Explicit time step =', dtexp
            write(6,*) 'Desired, actual time steps =', dtfac, dt
         endif

c        call CVODE to solve for new time step
         if(outnode)  write(6,*) '\n'
         call fcvode(tdesired, ttot, ux, taskchoice, cvierr)

c        wait for all processes to catch up
#ifdef PARALLEL
         call mpi_barrier(comm3d, ierr)
#endif

c        do a final update of boundary information [?is this necessary?]
         call CommState(ux)
         
c        in one-step modes, output diagnostics and data every
c        ndump steps; output BDF, dt info on other steps.  
         if (taskchoice==2 .or. taskchoice==4) then
            if (mod(timeStep,ndump) == 0) then
               call fcvdiags(iopt,ropt)
               call OutputData(ux, phi, timeStep, output_flag, 
     &              binary_flag, dt, ttot)
            else
               if (outnode) then
                  write(6,'(4x,A,i3)') 'Current BDF Order =',iopt(9)
                  write(6,'(4x,A,es12.5)') 'Current step size =',ropt(3)
               endif
            endif
c        in normal modes, output diagnostics 
         else
            call fcvdiags(iopt,ropt)
         endif

c        check error flag and exit if required
         call fcverr(cvierr)
         if ( cvierr == 1 ) then
            call OutputData(ux, phi, timeStep, output_flag, 
     &           binary_flag, dt, ttot)
            exit
         end if

c        periodically dump data
         if(mod(timeStep,ndump)==0) then
            call OutputData(ux, phi, timeStep, output_flag, 
     &           binary_flag, dt, ttot)
         endif

c        root outputs timing information
         if(outnode) 
     &        write(6,'(A,i7,2(A,es12.5))') 'timeStep =',timeStep,
     &           ',   dt =',ttot-told,',   ttot =',ttot


c        periodically output diagnostic information
         if(mod(timeStep,ndiag)==0) then
            call Diagnostics(ux,phi,ttot)
c            call WaveSolutionError(ux,ttot)
         endif

c        exit time stepping loop if we have exceeded tstop
         if ( ttot >= tstop ) then
            if (outnode) 
     &           write(6,'(A,es9.2,A)') '\n Physical stop time, t =',
     &              tstop,' reached.  Stopping simulation.'
            call Diagnostics(ux,phi,ttot)
c            call WaveSolutionError(ux,ttot)
            exit
         endif

      enddo

c     stop MPI timer
      call prof_stop(1)

c     output final solution
       call OutputData(ux, phi, timeStep, output_flag, binary_flag, 
     &     dt, ttot)

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
      call WrapUp(ux,phi,dt,ttot,timestep)

c     free up allocated memory
      call CvodeFree
#ifdef DYNAMIC
      deallocate(ux,phi,abstol)
#endif


      return
      end subroutine ImpCVODEDriver
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
      call ImpCVODEDriver

c     free space used by directional meshpoints
      deallocate(xc,yc,zc)

      return
      end subroutine MHDMain
c=================================================================
