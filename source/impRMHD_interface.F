c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c=======================================================================

      module timelevel_data
c-----------------------------------------------------------------
c     Description: holds necessary timelevel data for nonlinear
c        function evaluation
c
c     Contains the internal variables:
c         uold - (dbl(*)) state vector for previous full timestep
c        dFold - (dbl(*)) divFlux for previous full timestep
c        theta - (dbl) theta-scheme parameter
c           dt - (dbl) current timestep increment
c           vv - primitive variables
c         *inv - hyperbolic fluxes
c         *vis - parabolic fluxes
c         dudx - velocity gradient 
c         dbdx - magnetic field gradient
c     
c     Note:  the choice of theta will determine the two-step 
c        time-updating method.  theta = 0 corresponds to forward
c        Euler, theta = 1 corresponds to backward Euler, and 
c        theta = 0.5 corresponds to Crank-Nicholson.  The rest of 
c        these variables are merely used for internal data 
c        storage. 
c-----------------------------------------------------------------
      use mesh_parms

      implicit none
      save
      
      integer :: freezePrec=0
      double precision :: theta, dt
#ifdef DYNAMIC
      double precision, dimension(:,:,:,:), allocatable :: uold, dFold, 
     &     dFnew, vv, finv, ginv, hinv, fvis, gvis, hvis
      double precision, dimension(:,:,:,:,:), allocatable :: dudx, dbdx
#else
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &     :: uold, dFold, dFnew, vv, finv
#ifndef NOVISCFLUX
      double precision :: fvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif
#ifndef ONE_D
      double precision :: ginv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#ifndef NOVISCFLUX
      double precision :: gvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif
#ifndef TWO_D
      double precision :: hinv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#ifndef NOVISCFLUX
      double precision :: hvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif
#endif
#endif

#ifndef NOVISCFLUX
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3)
     &     :: dudx, dbdx
#endif
#endif
      end module timelevel_data
c=================================================================




      subroutine fkfun(uu, fval)
c-----------------------------------------------------------------------
c     Description: fkfun provides the nonlinear system function 
c        f(u) for the implicit time integration of one time step.
c        It follows a theta-scheme for the nonlinear evaluation, in
c        which the nonlinear residual fval is calculated as
c            fval = uu - uold - dt*(th*f(uu) + (1-th)*f(uold)),
c        and the rhs function f may be calculated as the divergence 
c        of the flux for the MHD system.  All of these variables are 
c        discussed below.
c
c     Arguments:
c          uu - (dbl(*), input) new fluid variables
c        fval - (dbl(*), output) nonlinear residual function
c
c     Note: all of the data required to evaluate the system function, 
c        is handled through Fortran modules, since the data cannot be 
c        sent through the C interface.
c
c     Includes modules mesh.h, [mesh_parms.h, mesh_uparms.h], 
c        properties.h, mpistuff.h and timelevel_data.h.  These bring 
c        the following relevant variables into the available use:
c
c        mesh.h: includes mesh_parms.h and mesh_uparms.h
c
c        mesh_parms.h: provides the local domain, discretization vars.
c           ixlo,ixhi,iylo,iyhi,izlo,izhi - (int) array bounds (w/ ghost)
c
c        mesh_uparms.h: provides the further local domain variables
c           nvar - (int) number of fluid species in MHD model
c
c        timelevel_data.h: provides the timelevel data
c           uold  - state vector for previous full timestep
c           dFold - div(flux) for previous full timestep
c           dFnew - div(flux) for the new timestep
c           theta - theta-scheme parameter
c           dt    - current timestep increment
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use timelevel_data
      use profiling

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &     :: uu, fval

c     local variables
      double precision :: onemth
      integer :: i, j, k, l

c======= Internals ============

c     exchange ghost values and set the boundary conditions
      call CommState(uu)

c     start the RHS profiler 
      call prof_start(2)

c     get the divergence of the flux at the new state uu
      call divFlux(uu,dFnew)

c     combine all parts into the return vector fval, only adjusting 
c     the interior volume elements (implements a standard theta-scheme)
      onemth = 1.0d0-theta
      fval = uu - uold - dt*(theta*dFnew + onemth*dFold)
      
c     stop the RHS profiler 
      call prof_stop(2)

      return
      end subroutine fkfun
c=======================================================================




      subroutine CommState(uu)
c-----------------------------------------------------------------------
c     Description: CommState handles the inter-process 
c        communications necessary for the computation of the 
c        nonlinear residual in the subroutine fkfun.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at new time
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      use profiling

c======= Declarations =========
      implicit none
      double precision :: uu(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)

c======= Internals ============

c     start the communication profiler
      call prof_start(3)

c     do the communications
      call mesh_update_bdry_async(uu, nvar)
      call SetBoundaryValues(uu)

c     wait for all processes to catch up
#ifdef PARALLEL
      call mpi_barrier(comm3d, ierr)
#endif

c     stop the communication profiler
      call prof_stop(3)

      return
      end subroutine CommState
c=================================================================





      subroutine KinsolInit(uscale, fscale, iopt, ropt, 
     &     dtfac, tstop, globalstrat)
c-----------------------------------------------------------------------
c     Description: KinsolInit handles the initialization of the KINSOL 
c     solver and vector data structures.
c
c     Arguments:
c           uscale - (dbl(*)) variable scalings
c           fscale - (dbl(*)) equation scalings
c           iopt   - (int(*)) integer KINSOL options array
c           ropt   - (dbl(*)) real KINSOL options array
c           dtfac  - (dbl) time step factor
c           tstop  - (dbl) stopping time for simulation
c           globalstrat  - (int) globalization strategy
c
c
c     Explanation of nonlinear solver parameters in kincontrol.inp:
c
c     specifing default vs custom routines (0=>default, 1=>custom)
c         useMyJac  = 1  =>  use customized fkjtimes routine
c         useMyPrec = 1  =>  use customized fkpset, fkpsol routines
c         freezePrec = 1 =>  only set up preconditioner on init.
c
c     general Newton-Krylov information
c         globalstrat:  globalization strategy, 
c                       0 => inexact newton, 1 => linesearch
c         fnormtol:  scaled nonlinear residual norm tolerance
c         scsteptol: scaled step norm tolerance
c         lsolver: choice of iterative linear solver
c                       0 => default (GMRES)
c                       1 => GMRES
c                       2 => Bi-CGStab
c                       3 => TFQMR
c         maxl:      max number of linear iterations
c         maxlrst:   max num GMRES restarts (GMRES only)
c         msbpre:    num Newton steps between precond setups
c         kinprint:  kinsol output print level with choices
c                       0 => none, 
c                       1 => newt iters, ||resid||, # f calls
c                       2 => 1 + globalstrat statistics
c                       3 => 2 + Krylov iteration stats
c         mxnewtit:  max num Newton iterations
c         initpcond: flag to suppress initial pcond setup call:
c                       0 => allow, 1 => suppress
c         inchoice:  choice of inexact forcing term:
c                       1 => eta = abs(||f{k+1}||-||fk+Jk*p||/||fk||)
c                       2 => eta = etagam*(||f_{k+1}||/||fk||)^etaalph
c                       3 => eta = etaconst
c         suppmintol: flag to suppress minimum tolerance:
c                       0 => allow, 1 => suppress
c         mxnewtstep: maximum allowable size of Newton step (0=>default)
c         frelerr:   relative error in computation of residual f(x)
c         etaconst, etagamma, etaalpha: Inexact Newton forcing params
c         u1sc, ... u8sc: scaling values for each species of state 
c                    vector, so that elements of uXsc*u(X,:) have 
c                    similar magnitude at solution
c         f1sc, ...,  f8sc: scaling values for each species of residual, 
c                    so that elements of fXsc*f(:,:,:,X) have similar 
c                    magnitude away from solution. 
c         u1con, ..., u8con: constraint values for each species of 
c                    state vector
c                       0 => no constraints
c                       1 => value must be >= 0
c                      -1 => value must be <= 0
c                       2 => value must be >  0
c                      -2 => value must be <  0
c
c     implicit time-stepping parameters
c         theta: implicit scheme:
c                       0.0d0 => Explicit Euler
c                       0.5d0 => Crank-Nicholson
c                       1.0d0 => Implicit Euler
c         dtfac: time step multiplication factor (dt = dtfac*dt_exp)
c         tstop: physical stopping time
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use boundary_conds
      use sundials_fcmix
      use timelevel_data
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uscale, fscale
      double precision :: dtfac, tstop, ropt(40)
      integer(kind=sund_long_int) :: iopt(40)
      integer          :: globalstrat

      integer(kind=sund_long_int) :: msbpre, nglobal
      integer(kind=sund_long_int) :: xlocal, ylocal, zlocal, nspecies 
      integer(kind=sund_long_int) :: xghost, yghost, zghost
      integer          :: useMyJac, useMyPrec
      integer          :: lsolver, maxl, maxlrst
      integer          :: kinprint, mxnewtit, initpcond
      integer          :: inchoice, suppmintol, kier, kiertot
      integer          :: ADIcorrection 
      logical          :: outnode, iofile 
      double precision :: fnormtol, scsteptol
      double precision :: mxnewtstep, frelerr, etaconst
      double precision :: etagamma, etaalpha, u1sc, u2sc
      double precision :: u3sc, u4sc, u5sc, u6sc, u7sc, u8sc
      double precision :: f1sc, f2sc, f3sc, f4sc, f5sc, f6sc
      double precision :: f7sc, f8sc, u1con, u2con, u3con
      double precision :: u4con, u5con, u6con, u7con, u8con
      namelist /kinsol_controls/ useMyJac, useMyPrec, freezePrec,
     &     globalstrat, lsolver, maxl, maxlrst, msbpre, kinprint, 
     &     mxnewtit, initpcond, inchoice, suppmintol, fnormtol, 
     &     scsteptol, mxnewtstep, frelerr, etaconst, etagamma,  
     &     etaalpha, u1sc, u2sc, u3sc, u4sc, u5sc, u6sc, u7sc, u8sc, 
     &     f1sc, f2sc, f3sc, f4sc, f5sc, f6sc, f7sc, f8sc, u1con, 
     &     u2con, u3con, u4con, u5con, u6con, u7con, u8con, 
     &     theta, dtfac, tstop

c======= Internals ============

c     if using dynamic meshing, allocate timelevel_data structures
#ifdef DYNAMIC
      allocate( uold(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      uold = 0.d0
      allocate(dFold(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      dFold = 0.d0
      allocate(dFnew(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      dFnew = 0.d0
      allocate(   vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      vv = 0.d0
      allocate( finv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      finv = 0.d0
#ifndef NOVISCFLUX
      allocate( fvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      fvis = 0.d0
      allocate( dudx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3))
      dudx = 0.d0
      allocate( dbdx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3))
      dbdx = 0.d0
#endif
#ifndef ONE_D
      allocate( ginv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      ginv = 0.d0
#ifndef NOVISCFLUX
      allocate( gvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      gvis = 0.d0
#endif
#ifndef TWO_D
      allocate( hinv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      hinv = 0.d0
#ifndef NOVISCFLUX
      allocate( hvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      hvis = 0.d0
#endif
#endif
#endif
#ifdef ENTROPY_FORMULATION
      allocate( esrc(nxlsize,nylsize,nzlsize))
      esrc = 0.d0
#endif
#endif

c     get KINSOL parameters from file
      open(unit=411,file='kincontrol.inp',status='old')
      read(411,nml=kinsol_controls)
      close(unit=411,status='keep')

c     set up kinsol nonlinear solver information
      outnode = .false.
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1)  outnode = .true.

      if (outnode) then
         write(6,'(A,A,es9.2)') '\nImplicit Scheme:',
     &        '\n    theta =',theta
         if (lsolver == 2) then
            write(6,'(A)') '\nBiCG-Stab iterative linear solver'
         else if (lsolver == 3) then
            write(6,'(A)') '\nTFQMR iterative linear solver'
         else 
            write(6,'(A)') '\nGMRES iterative linear solver'
         endif
         write(6,'(A,A,i2,A,i2)') '\nPreconditioning:',
     &        '\n    useMyPrec =',useMyPrec,
     &        '\n   freezePrec =',freezePrec
         write(6,'(A,6(A,i6),7(A,es9.2),3(A,8(es8.1),A),A,es9.2)')
     &        '\nKINSOL solver parameters:',
     &        '\n         maxl =',maxl,
     &        '\n      maxlrst =',maxlrst,
     &        '\n       msbpre =',msbpre,
     &        '\n     mxnewtit =',mxnewtit,
     &        '\n   suppmintol =',suppmintol,
     &        '\n     inchoice =',inchoice,
     &        '\n     etaconst =',etaconst,
     &        '\n     etagamma =',etagamma,
     &        '\n     etaalpha =',etaalpha,
     &        '\n     fnormtol =',fnormtol,
     &        '\n    scsteptol =',scsteptol,
     &        '\n   mxnewtstep =',mxnewtstep,
     &        '\n      frelerr =',frelerr,
     &        '\n       uscale = (',u1sc, u2sc, u3sc, u4sc, u5sc, 
     &                              u6sc, u7sc, u8sc,')',
     &        '\n       fscale = (',f1sc, f2sc, f3sc, f4sc, f5sc, 
     &                              f6sc, f7sc, f8sc,')',
     &        '\n       constr = (',u1con, u2con, u3con, u4con, 
     &                              u5con, u6con, u7con, u8con,')',
     &        '\n        dtfac =',dtfac
         write(6,*) '  '
      endif

c     initialize Kinsol data & operators
c        Nvector specs
      nglobal = (ixhi-ixlo+1)*(iyhi-iylo+1)*(izhi-izlo+1)*nvar
      call fnvinitmhd(3, nglobal, kier)
#ifdef PARALLEL
      call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (kiertot /= 0) then
         if(outnode) print *,' impRMHD_driver: fnvinitmhd failed'
         if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &        ': ier =',kier
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (kier /= 0) then
         print '(A,A,i2)',' impRMHD_driver: fnvinitmhd',
     &        ' failed, ier =',kier 
         return
      end if
#endif

c        kinsol dataspace
      if (outnode)  write(6,*) 'Allocating KINSOL dataspace'
      call fkinmalloc(iopt, ropt, kier)
#ifdef PARALLEL
      call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (kiertot /= 0) then
         if(outnode) print *,' impRMHD_driver: fkinmalloc failed'
         if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &        ': ier =',kier
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (kier /= 0) then
         print '(A,A,i2)',' impRMHD_driver: fkinmalloc',
     &        ' failed, ier =',kier 
         return
      end if
#endif

      
c        kinsol options
      if (outnode)  write(6,*) 'Setting solver options'
      if (msbpre /= 0) then
         call fkinsetiin('MAX_SETUPS', msbpre, kier)
      endif
      if (fnormtol > 0.d0) then
         call fkinsetrin('FNORM_TOL', fnormtol, kier)
      endif
      if (scsteptol > 0.d0) then
         call fkinsetrin('SSTEP_TOL', scsteptol, kier)
      endif
      if (outnode) then
         call fkinsetiin('PRNT_LEVEL', kinprint, kier)
      else
         call fkinsetiin('PRNT_LEVEL', 0, kier)
      endif
      if (mxnewtit /= 0) then
         call fkinsetiin('MAX_NITERS', mxnewtit, kier)
      endif
      if (initpcond /= 0) then
         call fkinsetiin('NO_INIT_SETUP', initpcond, kier)
      endif
      if (suppmintol /= 0) then
         call fkinsetiin('NO_MIN_EPS', suppmintol, kier)
      endif
      if (mxnewtstep > 0.d0)  then
         call fkinsetrin('MAX_STEP', mxnewtstep, kier)
      endif
      if (frelerr > 0.d0) then
         call fkinsetrin('RERR_FUNC', frelerr, kier)
      endif
      if (inchoice /= 0) then
         call fkinsetiin('ETA_FORM', inchoice, kier)
      endif
      if (etaconst > 0.d0) then
         call fkinsetrin('ETA_CONST', etaconst, kier)
      endif
      if ((etagamma > 0.d0) .and. (etaalpha > 0.d0)) then
         call fkinsetrin('ETA_PARAMS', (/etagamma, etaalpha/), kier)
      endif
         

c        constraints vector
      uscale(:,:,:,1) = u1con;  uscale(:,:,:,2) = u2con
      uscale(:,:,:,3) = u3con;  uscale(:,:,:,4) = u4con
      uscale(:,:,:,5) = u5con;  uscale(:,:,:,6) = u6con
      uscale(:,:,:,7) = u7con;  uscale(:,:,:,8) = u8con
c        [we must access uscale values before sending into 
c         fkinmalloc, or else the setting of these entries 
c         will be optimized elsewhere.]
      if (maxval(uscale) > 0.d0) then
         if (outnode)  write(6,*) 'Setting constraints vector'
      endif
      call fkinsetvin('CONSTR_VEC', uscale, kier)
#ifdef PARALLEL
      call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (kiertot /= 0) then
         if(outnode) print *,' impRMHD_driver: fkinsetvin failed'
         if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &        ': ier =',kier
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (kier /= 0) then
         print '(A,A,i2)',' impRMHD_driver: fkinsetvin',
     &        ' failed, ier =',kier 
         return
      end if
#endif


c        iterative solver
c           BiCG-Stab
      if (lsolver == 2) then
         call fkinspbcg(maxl, kier)
#ifdef PARALLEL
         call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (kiertot /= 0) then
            if(outnode) print *,' impRMHD_driver: fkinspbcg failed'
            if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &           ': ier =',kier
            return
         end if
         call mpi_barrier(comm3d,ierr)
#else
         if (kier /= 0) then
            print '(A,A,i2)',' impRMHD_driver: fkinspbcg',
     &           ' failed, ier =',kier 
            return
         end if
#endif
c           TFQMR
      else if (lsolver == 3) then
         call fkinsptfqmr(maxl, kier)
#ifdef PARALLEL
         call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (kiertot /= 0) then
            if(outnode) print *,' impRMHD_driver: fkinsptfqmr failed'
            if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &           ': ier =',kier
            return
         end if
         call mpi_barrier(comm3d,ierr)
#else
         if (kier /= 0) then
            print '(A,A,i2)',' impRMHD_driver: fkinsptfqmr',
     &           ' failed, ier =',kier 
            return
         end if
#endif
c           GMRES
      else
         call fkinspgmr(maxl, maxlrst, kier)
#ifdef PARALLEL
         call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (kiertot /= 0) then
            if(outnode) print *,' impRMHD_driver: fkinspgmr failed'
            if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &           ': ier =',kier
            return
         end if
         call mpi_barrier(comm3d,ierr)
#else
         if (kier /= 0) then
            print '(A,A,i2)',' impRMHD_driver: fkinspgmr',
     &           ' failed, ier =',kier 
            return
         end if
#endif
      endif         

c        J*v routine specification
      call fkinspilssetjac(useMyJac, kier)
#ifdef PARALLEL
      call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (kiertot /= 0) then
         if(outnode) print *,' impRMHD_driver: fkinspilssetjac failed'
         if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &        ': ier =',kier
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (kier /= 0) then
         print '(A,A,i2)',' impRMHD_driver: fkinspilssetjac',
     &        ' failed, ier =',kier 
         return
      end if
#endif

c        preconditioner setup/solve (0 => default, 1 => custom routine)
      call fkinspilssetprec(useMyPrec, kier)
#ifdef PARALLEL
      call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (kiertot /= 0) then
         if(outnode) print *,' impRMHD_driver: fkinspilssetprec failed'
         if (kier /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &        ': ier =',kier
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (kier /= 0) then
         print '(A,A,i2)',' impRMHD_driver: fkinspilssetprec',
     &        ' failed, ier =',kier 
         return
      end if
#endif


c        set uscale and fscale vectors
      uscale(:,:,:,1) = u1sc;  uscale(:,:,:,2) = u2sc
      uscale(:,:,:,3) = u3sc;  uscale(:,:,:,4) = u4sc
      uscale(:,:,:,5) = u5sc;  uscale(:,:,:,6) = u6sc
      uscale(:,:,:,7) = u7sc;  uscale(:,:,:,8) = u8sc
      fscale(:,:,:,1) = f1sc;  fscale(:,:,:,2) = f2sc
      fscale(:,:,:,3) = f3sc;  fscale(:,:,:,4) = f4sc
      fscale(:,:,:,5) = f5sc;  fscale(:,:,:,6) = f6sc
      fscale(:,:,:,7) = f7sc;  fscale(:,:,:,8) = f8sc
c        [we must access uscale, fscale values before sending
c         to fkinsol, or else the setting of these entries 
c         will be optimized elsewhere (or so it seems).]
      if (maxval(uscale) > 0.d0) then
         if (outnode)  write(6,*) 'Set uscale vector'
      endif
      if (maxval(fscale) > 0.d0) then
         if (outnode)  write(6,*) 'Set fscale vector'
      endif
      if(outnode)  write(6,*) '\nFinished KINSOL initialization'

      return
      end subroutine KinsolInit
c=================================================================






      subroutine KinsolFree()
c-----------------------------------------------------------------------
c     Description: KinsolFree handles the initialization of the KINSOL 
c     solver and vector data structures.
c
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use timelevel_data

c======= Declarations =========
      implicit none

c======= Internals ============


c     free kinsol internal data
      call fkinfree

c     if using dynamic meshing, free timelevel_data structures
#ifdef DYNAMIC
      deallocate(uold,dFold,dFnew,vv,finv)
#ifndef NOVISCFLUX
      deallocate(fvis,dudx,dbdx)
#endif
#ifndef ONE_D
      deallocate(ginv)
#ifndef NOVISCFLUX
      deallocate(gvis)
#endif
#ifndef TWO_D
      deallocate(hinv)
#ifndef NOVISCFLUX
      deallocate(hvis)
#endif
#endif
#endif
#endif


      return
      end subroutine KinsolFree
c=================================================================






      subroutine fkjtimes(ww, Jw, new_uu, uu, kier)
c-----------------------------------------------------------------------
c     Description: fkjtimes provides the Jacobian-vector product routine
c        for the linearized Newton system, i.e. it calculates Jw = J*w, 
c        where J is equal, or approximately equal, to dF/du.
c     
c     Note: here this is a dummy routine, provided only to allow for 
c        compilation.  It is never called.
c
c     Arguments:
c            ww - (dbl(*), input) vector to multiply
c            Jw - (dbl(*), output) result of Jac*ww operation
c        new_uu - (int, input) flag dictating if state uu has changed
c            uu - (dbl(*), input) current iterate for system
c           kier - (int, output) return flag: success(0), failure(1)
c
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use timelevel_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &        :: ww, Jw, uu
      integer :: new_uu, kier

c======= Internals ============

c     No matrix-vector product implemented.  Return failure if called
      kier = 1

      return
      end subroutine fkjtimes
c=======================================================================





      subroutine fkpset(uu, usc, fu, fsc, v1, v2, kier)
c-----------------------------------------------------------------------
c     Description: fkpset provides the preconditioner setup 
c        routine for the linearized Newton system.
c
c     Arguments:
c          uu - (dbl(*), input) current iterate for system
c         usc - (dbl(*), input) scaling coefficients for uu
c          fu - (dbl(*), input) current nonlinear function value
c         fsc - (dbl(*), input) scaling coefficients for fu
c          v1 - (dbl(*)) temporary work array (structure of uu)
c          v2 - (dbl(*)) temporary work array (structure of uu)
c        kier - (int, output) return flag: success(0), failure(1)
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use timelevel_data
      use properties
      use profiling

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &                 :: uu, usc, fu, fsc, v1, v2
      integer          :: kier
      
c     local variables
      double precision :: gamdt, Re, Pr, Lu
      integer          :: ier

c======= Internals ============

c     initialize return flag to success
      kier = 0

c     if freezePrec is set to 2, return with success
      if (freezePrec == 2) then
         return
      endif

c     start the preconditioner setup profiler
      call prof_start(4)

c     if freezePrec is set to 1, set to 2
      if (freezePrec == 1)  then
         freezePrec = 2
      endif

c     set gamdt (in the general Jacobian sense of Jtot = I - gamdt*Jf)
      gamdt = theta*dt

c     set up preconditioners

c     stop the preconditioner setup profiler
      call prof_stop(4)

      return
      end subroutine fkpset
c=======================================================================




      subroutine fkpsol(uu, usc, fu, fsc, zz, vt, kier)
c-----------------------------------------------------------------------
c     Description: fkpsol provides the preconditioner solve 
c        routine for the preconditioning of the linearized Newton
c        system.  (solves P*x = z)
c
c     Arguments:
c         uu - (dbl(*), input) current iterate for system
c        usc - (dbl(*), input) scaling coefficients for uu
c         fu - (dbl(*), input) current nonlinear function value
c        fsc - (dbl(*), input) scaling coefficients for fu 
c         zz - (dbl(*), I/O) work array, holds the RHS vector w 
c               on input, and the result x on output/return 
c         vt - (dbl(*)) work array of same size as uu
c       kier - (int, output) return flag: success(0), failure(1)
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use timelevel_data
      use profiling

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, usc, fu, fsc, zz, vt
      integer          :: kier

c     local variables
      double precision :: gamdt, delta, fnorm

c======= Internals ============

c     set return flag to success
      kier = 0

c     start the preconditioner solve profiler
      call prof_start(5)

c     set gamdt (in the general Jacobian sense of Jtot = I - gamdt*Jf)
      gamdt = theta*dt

c     call preconditioner solvers

c     stop the preconditioner solve profiler
      call prof_stop(5)

      return
      end subroutine fkpsol
c=======================================================================




      subroutine kinerr(kier,iopt,ropt)
c-----------------------------------------------------------------------
c     Description: function to check output from fkinsol for 
c     success/failure as well as interpreting the various output flags
c
c     Arguments:
c            kier - (in/out) conserved variables before/after filter
c                   output value of 0=>success, 1=>failure
c            iopt - (in) integer optional outputs from fkinsol
c            ropt - (in) real optional outputs from fkinsol
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
      use sundials_fcmix

c======= Declarations =========
      implicit none

      integer(kind=sund_long_int)  :: iopt(40)
      double precision, intent(in) :: ropt(40)
      integer                      :: kier, kierout, kiermin


c======= Internals ============


c     if parallel, gather all kier together
#ifdef PARALLEL
      call MPI_Allreduce( kier, kiermin, 1, MPI_INTEGER, MPI_MIN, 
     &     comm3d, ierr )
      if ( ierr /= 0 ) then
         print '(A,i4,A,i2)', 'kinerr: mpi_allreduce failed, my_id =',
     &        my_id,',  ierr =', ierr
         kiermin = -1
         return
      end if
#else
      kiermin = kier
#endif
      
c     error if failure anywhere
      if ( kiermin < 0 ) then
         kierout = 1
      
#ifdef PARALLEL
         if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            write(6,'(A,i4,A)')' kinerr: fkinsol failed, my_id =',
     &           my_id,':'
         end if
#else
         write(6,*) ' kinerr: fkinsol failed:'
#endif
       
c        output relevant error message
         select case (kier)
         case (-1) 
            write(6,*)'     [not enough memory]\n'
         case (-2) 
            write(6,*)'     [illegal input]\n'
         case (-3) 
            write(6,*)'     [no malloc available]\n'
         case (-4) 
            write(6,*)'     [insuff. memory for linear solver]\n'
         case (-5) 
            write(6,*)'     [nonconvergent linesearch]\n'
         case (-6) 
            write(6,*)'     [maxiter reached]\n'
         case (-7) 
            write(6,*)'     [mxnewt 5X exceeded]\n'
         case (-8) 
            write(6,*)'     [linesearch bc fail]\n'
         case (-9) 
            write(6,*)'     [kyrlov failure]\n'
         case (-10) 
            write(6,*)'     [linear solver init fail]\n'
         case (-11) 
            write(6,*)'     [preconditioner setup failure]\n'
         case (-12) 
            write(6,*)'     [precond./linear solve failure]\n'
         case (-13) 
            write(6,*)'     [pdata NULL]\n'
         end select
      else
         kierout = 0
      end if

c     root outputs solver statistics
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) ' KINSol output values:'
         write(6,'(4x,A,i5)')'Num. Newton Iterations    =',iopt(3)
         write(6,'(4x,A,i5)')'Num. Nonlin. Func. Evals  =',iopt(4)
         if (iopt(5) /= 0) then
            write(6,'(4x,A,i5)')'Num. linsearch beta fails =',iopt(5)
         endif
         if (iopt(6) /= 0) then
            write(6,'(4x,A,i5)')'Num. linesearch bktracks  =',iopt(6)
         endif
         write(6,'(4x,A,es12.5)')'Final scaled norm of f(u) =',ropt(1)
         write(6,'(4x,A,es12.5)')'Scaled final steplength   =',ropt(2)
         write(6,'(4x,A,i5)')'Number Krylov iterations  =',iopt(14)
         if (iopt(12) /= 0) then
            write(6,'(4x,A,i5)')'Num. Precond. evaluations =',iopt(12)
         endif
         if (iopt(13) /= 0) then
            write(6,'(4x,A,i5)')'Num. Precond. solves      =',iopt(13)
         endif
         if (iopt(15) /= 0) then
            write(6,'(4x,A,i5)')'Num. Lin. Conv. Failures  =',iopt(15)
         endif
         if (iopt(9) /= 0) then
            write(6,'(4x,A,i5)')'Final Lin Solver Ret Flag =',iopt(9)
         endif
         print *, '   '
      end if

c     set final output
      kier = kierout

      
      return
      end subroutine kinerr
c=================================================================




      subroutine fkfunDiags(ux, fsc)
c-----------------------------------------------------------------
c     Description: calls fkfun and outputs all relevant 
c     diagnostics on the nonlinear residual function
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none

      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &                 :: ux, fsc, fval
      double precision :: f2n, f2nall, fIn, fInall
      double precision :: rho2n, rho2nall, rhoIn, rhoInall
      double precision :: ux2n, ux2nall, uxIn, uxInall
      double precision :: uy2n, uy2nall, uyIn, uyInall
      double precision :: uz2n, uz2nall, uzIn, uzInall
      double precision :: bx2n, bx2nall, bxIn, bxInall
      double precision :: by2n, by2nall, byIn, byInall
      double precision :: bz2n, bz2nall, bzIn, bzInall
      double precision :: e2n, e2nall, eIn, eInall
      integer          :: maxproc, maxnormloc(4)
      integer          :: maxlocsend(5), maxlocrecv(5)

c======= Internals ============

c     call fkfun for the current nonlinear residual, scale appropriately
      fval = 0.d0
      call fkfun(ux,fval)
      fval = fval*fsc

c     compute local residual norms (total and for each equation)
c        ||f||
      call fnvmhdwrmsnorm(fval,fsc,f2n)
      call fnvmhdmaxnorm(fval,fIn)

c        ||rho||
      call fnvmhdwrmsnormvar(fval,fsc,1,rho2n)
      call fnvmhdmaxnormvar(fval,1,rhoIn)

c        ||ux||
      call fnvmhdwrmsnormvar(fval,fsc,2,ux2n)
      call fnvmhdmaxnormvar(fval,2,uxIn)

c        ||uy||
      call fnvmhdwrmsnormvar(fval,fsc,3,uy2n)
      call fnvmhdmaxnormvar(fval,3,uyIn)

c        ||uz||
      call fnvmhdwrmsnormvar(fval,fsc,4,uz2n)
      call fnvmhdmaxnormvar(fval,4,uzIn)

c        ||bx||
      call fnvmhdwrmsnormvar(fval,fsc,5,bx2n)
      call fnvmhdmaxnormvar(fval,5,bxIn)

c        ||by||
      call fnvmhdwrmsnormvar(fval,fsc,6,by2n)
      call fnvmhdmaxnormvar(fval,6,byIn)

c        ||bz||
      call fnvmhdwrmsnormvar(fval,fsc,7,bz2n)
      call fnvmhdmaxnormvar(fval,7,bzIn)

c        ||e||
      call fnvmhdwrmsnormvar(fval,fsc,8,e2n)
      call fnvmhdmaxnormvar(fval,8,eIn)


c     root node outputs nonlinear residual information
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(2(A,es11.4))')
     &        '  ||f(u)||  =',f2n,  ',   max|f(u)|  =',fIn
         write(6,'(A)') 
     &        '  -------------------------------------------------'
         write(6,'(2(A,es11.4))')
     &        '  ||f_rho|| =',rho2n,',   max|f_rho| =',rhoIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_ux||  =',ux2n, ',   max|f_ux|  =',uxIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_uy||  =',uy2n, ',   max|f_uy|  =',uyIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_uz||  =',uz2n, ',   max|f_uz|  =',uzIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_Bx||  =',bx2n, ',   max|f_Bx|  =',bxIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_By||  =',by2n, ',   max|f_By|  =',byIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_Bz||  =',bz2n, ',   max|f_Bz|  =',bzIn
         write(6,'(2(A,es11.4))')
     &        '  ||f_e||   =',e2n,  ',   max|f_e|   =',eIn
      end if

      return
      end subroutine fkfunDiags
c=================================================================





      subroutine divFlux(uu,dFlux)
c-----------------------------------------------------------------------
c     Description: divFlux provides the divergence of the flux for 
c        the conservation equation describing the nonlinearly-resistive,
c        ideal MHD system.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at new time
c        dFlux - (dbl(*), output) divergence of the flux
c
c     Includes modules mesh.h, [mesh_parms.h, mesh_uparms.h] and 
c        timelevel_data.h.  These bring the following relevant 
c        variables into the available use:
c
c        mesh.h: provides the general domain, discretization variables
c           dx,dy,dz - (dbl) mesh spacing in each direction
c
c        mesh_parms.h: provides the local domain, discretization vars.
c           ixlo,ixhi,iylo,iyhi,izlo,izhi - (int) array bounds (w/ ghost)
c           nxlocal,nylocal,nzlocal - (int) local mesh lengths (no ghost)
c
c        mesh_uparms.h: provides the further local domain variables
c           nvar - (int) number of fluid species in MHD model
c
c        timelevel_data.h: provides the timelevel data
c           dt    - current timestep increment
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use timelevel_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &     :: uu, dFlux

c     differencing factors (depending on dimension)
      double precision :: xfac, yfac, zfac

c     indexes for flux calculations
      integer :: indx, indy, indz, i, j, k, l, m
      


c======= Internals ============


c     set differencing factors
      xfac = 1.0d0/dx;  yfac = 1.0d0/dy;  zfac = 1.0d0/dz


c     construct the vector of primitive variables
      call ConservedToPrimitive(uu,vv)


c     zero the fluxes (necessary???)
      finv = 0.0d0
#ifndef ONE_D
      ginv = 0.0d0
#ifndef TWO_D
      hinv = 0.0d0
#endif
#endif
#ifndef NOVISCFLUX
      fvis = 0.0d0
#ifndef ONE_D
      gvis = 0.0d0
#ifndef TWO_D
      hvis = 0.0d0
#endif
#endif
#endif

c     get the viscous flux requirements dudx and dbdx
#ifndef NOVISCFLUX
      do m=1,3,1
         dudx(:,:,:,1,m) = vv(:,:,:,m+1)
         dbdx(:,:,:,1,m) = vv(:,:,:,m+4)
         dudx(:,:,:,2,m) = 0.d0
         dbdx(:,:,:,2,m) = 0.d0
         dudx(:,:,:,3,m) = 0.d0
         dbdx(:,:,:,3,m) = 0.d0
      end do
      call GradientVector(dudx)
      call GradientVector(dbdx)
#endif

c     construct the face-centered x-directional fluxes
c
c     Hyperbolic fluxes
      inlo = ixlo;  inhi = ixhi
      call InviscidFluxWrapper(finv,vv,1,dt)
c
c     Parabolic fluxes
      inlo = ixlo;  inhi = ixhi
      indx = 1;  indy = 0;  indz = 0
      call ViscousFluxWrapper(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)


c     construct the face-centered y-directional fluxes (2D,2.5D,3D)
#ifndef ONE_D
c
c     Hyperbolic fluxes
      inlo = iylo;  inhi = iyhi
      call InviscidFluxWrapper(ginv,vv,2,dt)
c
c     Parabolic fluxes
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)


c     construct the face-centered z-directional fluxes (3D)
#ifndef TWO_D
c
c     Hyperbolic fluxes
      inlo = izlo;  inhi = izhi
      call InviscidFluxWrapper(hinv,vv,3,dt)
c
c     Parabolic fluxes
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif


c     do the final construction of the new div(flux) term
c
c     Density
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,1) = -xfac*(finv(i,j,k,1) - finv(i-1,j,k,1))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,1) - ginv(i,j-1,k,1) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,1) - hinv(i,j,k-1,1) )
#endif
#endif
            end do
         end do
      end do
c
c     X-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,2) = -xfac*(finv(i,j,k,2) - finv(i-1,j,k,2))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,2) - fvis(i-1,j,k,2) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,4) - ginv(i,j-1,k,4) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,4) - gvis(i,j-1,k,4) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,3) - hinv(i,j,k-1,3) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,3) - hvis(i,j,k-1,3) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Y-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,3) = -xfac*(finv(i,j,k,3) - finv(i-1,j,k,3))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,3) - fvis(i-1,j,k,3) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,2) - ginv(i,j-1,k,2) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,2) - gvis(i,j-1,k,2) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,4) - hinv(i,j,k-1,4) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,4) - hvis(i,j,k-1,4) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Z-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,4) = -xfac*(finv(i,j,k,4) - finv(i-1,j,k,4))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,4) - fvis(i-1,j,k,4) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,3) - ginv(i,j-1,k,3) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,3) - gvis(i,j-1,k,3) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,2) - hinv(i,j,k-1,2) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,2) - hvis(i,j,k-1,2) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     X-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,5) = -xfac*(finv(i,j,k,5) - finv(i-1,j,k,5))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,5) - fvis(i-1,j,k,5) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,7) - ginv(i,j-1,k,7) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,7) - gvis(i,j-1,k,7) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,6) - hinv(i,j,k-1,6) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,6) - hvis(i,j,k-1,6) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Y-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,6) = -xfac*(finv(i,j,k,6) - finv(i-1,j,k,6))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,6) - fvis(i-1,j,k,6) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,5) - ginv(i,j-1,k,5) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,5) - gvis(i,j-1,k,5) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,7) - hinv(i,j,k-1,7) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,7) - hvis(i,j,k-1,7) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Z-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,7) = -xfac*(finv(i,j,k,7) - finv(i-1,j,k,7))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,7) - fvis(i-1,j,k,7) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,6) - ginv(i,j-1,k,6) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,6) - gvis(i,j-1,k,6) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,5) - hinv(i,j,k-1,5) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,5) - hvis(i,j,k-1,5) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Energy
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,8) = -xfac*(finv(i,j,k,8) - finv(i-1,j,k,8))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,8) - fvis(i-1,j,k,8) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,8) - ginv(i,j-1,k,8) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,8) - gvis(i,j-1,k,8) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,8) - hinv(i,j,k-1,8) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,8) - hvis(i,j,k-1,8) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c

      return
      end subroutine divFlux
c=======================================================================
