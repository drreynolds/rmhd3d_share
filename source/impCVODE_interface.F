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

      module cvode_data
c-----------------------------------------------------------------
c     Description: holds necessary timelevel data for nonlinear
c        function evaluation
c
c     Contains the internal variables:
c           vv - primitive variables
c         *inv - hyperbolic fluxes
c         *vis - parabolic fluxes
c         dudx - velocity gradient 
c         dbdx - magnetic field gradient 
c-----------------------------------------------------------------
      use mesh_parms

      implicit none
      save
      
      integer :: freezePrec=0
#ifdef DYNAMIC
      double precision, dimension(:,:,:,:), allocatable ::
     &     vv, finv, fvis, ginv, gvis, hinv, hvis
      double precision, dimension(:,:,:,:,:), allocatable ::
     &     dudx, dbdx

#else

      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &     :: vv, finv
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

      end module cvode_data
c=================================================================




      subroutine fcvfun(t, uu, fdot, ipar, rpar, cvier)
c-----------------------------------------------------------------------
c     Description: fcvfun provides the right hand side function for the 
c        ODE:     du/dt = fdot(u,t)
c        For the current formulation of our PDE in conservation form, 
c        we have  du/dt = div(F(u)),  where F is the (nonlinear) flux.
c
c     Arguments:
c           t - (dbl, input) current time
c          uu - (dbl(*), input) fluid variables at current time
c        fdot - (dbl(*), output) right-hand side function
c
c     Note: all of the data required to evaluate the system function, 
c        is handled through Fortran modules, since the data cannot be 
c        sent through the C interface.
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use profiling

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, fdot
      double precision :: t, rpar(*)
      integer          :: i, j, k, l, ipar(*), cvier


c======= Internals ============

c     set return value to success
      cvier = 0

c     exchange ghost values and set the boundary conditions
      call CommState(uu)

c     start the RHS profiler 
      call prof_start(2)

c     initialize fdot to zero
      fdot = 0.d0

c     get the divergence of the flux at the new state uu
      call divFlux(uu,fdot)

c     stop the RHS profiler 
      call prof_stop(2)
      
      return
      end subroutine fcvfun
c=======================================================================




      subroutine CommState(uu)
c-----------------------------------------------------------------------
c     Description: CommState handles the inter-process 
c        communication necessary for the computation of the 
c        rhs function in the subroutine fcvfun.
c
c     Arguments:
c           uu - (dbl(*), input) current fluid variables
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
      call mesh_update_bdry_async(uu,nvar)
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





      subroutine CvodeInit(ttot, uu, abstol, iopt, ropt, 
     &     dtfac, tstop, taskchoice)
c-----------------------------------------------------------------------
c     Description: CvodeInit initializes the CVODE solver.
c
c     Arguments:
c         ttot - (dbl) simulation start time
c           uu - (dbl(*)) current fluid variables
c       abstol - (dbl(*)) absolute tolerance array
c         iopt - (int(*)) integer solver parameters
c         ropt - (dbl(*)) real solver parameters
c        dtfac - (dbl) dump times
c        tstop - (dbl) simulation end time
c   taskchoice - (int) type of integration scheme
c
c------- Implicit Scheme / Nonlinear Solver Parameters -------------
c        (all parameters set in file cvcontrol.inp)
c
c     specifying default vs custom routines (0=>default, 1=>custom)
c         useMyJv   = 1  =>  use my fcvjtimes routine
c         useMyPrec = 1  =>  use my fcvpset, fcvpsol routines
c         freezePrec = 1 =>  only set up preconditioner on init.
c
c     required CVODE solver information
c         intmethod: integration method 
c                       1 => Adams (nonstiff), 2 => BDF (stiff)
c         itermethod: nonlinear iteration method
c                       1 => functional iteration, 2 => Newton
c         taskchoice: CVODE task indicator
c                       1 => normal mode (overshoot & interpolate)
c                       2 => return after each internal step taken
c                       3 => like 1, but must stop at tstop
c                       4 => like 2, but must stop at tstop
c                      NOTE: modes 3 and 4 imply that tstop/=0
c         lsolver: choice of iterative linear solver
c                       0 => default (GMRES)
c                       1 => GMRES
c                       2 => Bi-CGStab
c                       3 => TFQMR
c         gstype: Gram-Schmidt orthogonalization method (GMRES only)
c                       1 => modified GS, 2 => classical GS
c         rtol: relative tolerance for nonlinear time step, should 
c               be on the order of spatial accuracy
c         atol1->atol8: absolute tolerances for each variable, 
c               on the order of the noise in each variable
c         tstop: physical stopping time
c         dtfac: desired output time factor (dtout = dtcfl*dtfac)
c
c     optional CVODE solver information [0 => default]
c         prectype: preconditioning type (default=none)
c                       1 => left, 2 => right, 3 => both
c         maxl: maximum Krylov subspace dimension (default=5)
c         mxbdf: maximum order of BDF method (1-5, default=5)
c         mxsteps: max time steps to reach solution (default=500)
c         mxwarn: max warnings for tn + h = tn (default=10)
c         stlimit: use stability limit detection (0=>off, 1=>on)
c         mxerrfails: max error test fails per step (default=7)
c         mxnewtit: max nonlin solver iters per step (default=3)
c         mxconvfails: max convergence fails per step (default=10)
c         lintol: ratio between lin and nonlin tols (default=0.05)
c         dtinit: initial timestep size (default is estimated)
c         maxdt: maximum absolute step size (default=infinity)
c         mindt: minimum absolute step size (default=0.0)
c         nonlincoeff: coeff in nonlin conv test (default=0.1)
c-------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use boundary_conds
      use sundials_fcmix
      use cvode_data
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      double precision :: uu(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: ttot
      double precision :: abstol(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: tstop, dtfac, ropt(40)
      integer          :: taskchoice
      integer(kind=sund_long_int) :: iopt(40), idummy(1)

      integer          :: useMyJv, useMyPrec, intmethod, itermethod
      integer          :: prectype, lsolver, gstype, maxl, mxbdf
      integer          :: mxsteps, mxwarn, stlimit, mxerrfails
      integer          :: mxnewtit, mxconvfails, cvierr, cvierrtot
      logical          :: outnode, iofile
      double precision :: lintol, dtinit, maxdt, mindt, nonlincoeff
      double precision :: rtol, atol1, atol2, atol3, atol4, atol5
      double precision :: atol6, atol7, atol8, rdummy(1)
      integer(kind=sund_long_int) :: xlocal, ylocal, zlocal, nspecies
      integer(kind=sund_long_int) :: xghost, yghost, zghost, nglobal
      namelist /cvode_controls/ useMyJv, useMyPrec, freezePrec,
     &     intmethod, itermethod, taskchoice, prectype, lsolver, gstype, 
     &     maxl, mxbdf, mxsteps, mxwarn, stlimit, mxerrfails, mxnewtit,  
     &     mxconvfails, lintol, dtinit, maxdt, mindt, nonlincoeff, rtol, 
     &     atol1, atol2, atol3, atol4, atol5, atol6, atol7, atol8, 
     &     tstop, dtfac

c======= Internals ============

c     if using dynamic meshing, allocate cvode_data structures
#ifdef DYNAMIC
      allocate(   vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate( finv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#ifndef NOVISCFLUX
      allocate( fvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate( dudx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3))
      allocate( dbdx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3))
#endif
#ifndef ONE_D
      allocate( ginv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#ifndef NOVISCFLUX
      allocate( gvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#endif
#ifndef TWO_D
      allocate( hinv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#ifndef NOVISCFLUX
      allocate( hvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#endif
#endif
#endif
#endif

c     see if this is the output node
      outnode = .false.
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1)  outnode = .true.

c     input cvode solver parameters
      open(411,file='cvcontrol.inp',form='formatted')
      read(411,cvode_controls)
      close(411)

c     root node outputs problem parameters
      if(outnode) then
         write(6,'(A,A,es9.2)') '\nImplicit Scheme:',
     &        '\n    Adaptive 5th-order BDF'
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
         write(6,'(A,3(A,i6),11(A,es9.2))') 
     &        '\nCVODE Solver parameters:',
     &        '\n           maxl =',maxl,
     &        '\n          mxbdf =',mxbdf,
     &        '\n        mxsteps =',mxsteps,
     &        '\n         lintol =',lintol,
     &        '\n    nonlincoeff =',nonlincoeff,
     &        '\n           rtol =',rtol,
     &        '\n          atol1 =',atol1,
     &        '\n          atol2 =',atol2,
     &        '\n          atol3 =',atol3,
     &        '\n          atol4 =',atol4,
     &        '\n          atol5 =',atol5,
     &        '\n          atol6 =',atol6,
     &        '\n          atol7 =',atol7,
     &        '\n          atol8 =',atol8
      endif


c     if preconditioning is to be used, set it to right preconditioning
      prectype = 0
      if (useMyPrec == 1) prectype = 2
      

c     set up CVODE NVector interface information using mesh variables
      nglobal = nprocs*nxlsize*nylsize*nzlsize*nvar
      xghost = nghost;  yghost = 0;  zghost = 0
#ifndef ONE_D
      yghost = nghost
#ifndef TWO_D 
      zghost = nghost
#endif
#endif

c     ensure that we send the correct integer kind through the interface
      xlocal = nxlsize
      ylocal = nylsize
      zlocal = nzlsize
      nspecies = nvar


c     initialize CVODE data & operators

c        Nvector specs
      nglobal = (xlocal+2*xghost)*(ylocal+2*yghost)*(zlocal+2*zghost)
     &     *nspecies
      call fnvinitmhd(1, nglobal, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(outnode)  print '(A)',' impCVODE_driver: fnvinitmhd failed'
         if (cvierr /= 0) 
     &        print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVODE_driver: fnvinitmhd',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        CVODE dataspace
c        [we must access abstol values before sending
c         to fcvmalloc.]
      abstol(:,:,:,1) = atol1;  abstol(:,:,:,2) = atol2
      abstol(:,:,:,3) = atol3;  abstol(:,:,:,4) = atol4
      abstol(:,:,:,5) = atol5;  abstol(:,:,:,6) = atol6
      abstol(:,:,:,7) = atol7;  abstol(:,:,:,8) = atol8
      if (maxval(abstol) > 0.d0) then
         if (outnode)  write(6,*) 'Setting abstol vector'
      endif
      call fcvmalloc(ttot, uu, intmethod, itermethod, 2,
     &     rtol, abstol, iopt, ropt, idummy, rdummy, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(outnode)  print '(A,A)',' impCVODE_driver: fcvmalloc failed'
         if (cvierr /= 0) 
     &        print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVODE_driver: fcvmalloc',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        CVODE options
      if (outnode)  write(6,*) 'Setting CVODE options'
      if (mxbdf /= 0) then
         call fcvsetiin('MAX_ORD',mxbdf,cvierr)
      endif
      if (mxsteps /= 0) then
         call fcvsetiin('MAX_NSTEPS',mxsteps,cvierr)
      endif
      if (mxwarn /= 0) then
         call fcvsetiin('HNIL_WARNS',mxwarn,cvierr)
      endif
      if (stlimit /= 0) then
         call fcvsetiin('STAB_LIM',stlimit,cvierr)
      endif
      if (mxerrfails /= 0) then
         call fcvsetiin('MAX_ERRFAIL',mxerrfails,cvierr)
      endif
      if (mxnewtit /= 0) then
         call fcvsetiin('MAX_NITERS',mxnewtit,cvierr)
      endif
      if (mxconvfails /= 0) then
         call fcvsetiin('MAX_CONVFAIL',mxconvfails,cvierr)
      endif
      if (dtinit > 0.d0) then
         call fcvsetrin('INIT_STEP',dtinit,cvierr)
      endif
      if (maxdt > 0.d0) then
         call fcvsetrin('MAX_STEP',maxdt,cvierr)
      endif
      if (mindt > 0.d0) then
         call fcvsetrin('MIN_STEP',mindt,cvierr)
      endif
      if (tstop > ttot) then
         call fcvsetrin('STOP_TIME',tstop,cvierr)
      endif
      if (nonlincoeff > 0.d0) then
         call fcvsetrin('NLCONV_COEF',nonlincoeff,cvierr)
      endif


c        iterative solver
c           BiCG-Stab
      if (lsolver == 2) then
         call fcvspbcg(prectype, maxl, lintol, cvierr)
#ifdef PARALLEL
         call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (cvierrtot /= 0) then
            if(outnode)  print '(A,A)',
     &           ' impCVODE_driver: fcvspbcg failed'
            if (cvierr /= 0) 
     &           print '(A,i4,A,i2)','    proc ',my_id,'
     &           : ier =',cvierr
            return
         end if
         call mpi_barrier(comm3d,ierr)
#else
         if (cvierr /= 0) then
            print '(A,A,i2)',' impCVODE_driver: fcvspbcg',
     &           ' failed, ier =',cvierr 
            return
         end if
#endif
c           TFQMR 
      else if (lsolver == 3) then
         call fcvsptfqmr(prectype, maxl, lintol, cvierr)
#ifdef PARALLEL
         call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (cvierrtot /= 0) then
            if(outnode)  print '(A,A)',
     &           ' impCVODE_driver: fcvsptfqmr failed'
            if (cvierr /= 0) 
     &           print '(A,i4,A,i2)','    proc ',my_id,'
     &           : ier =',cvierr
            return
         end if
         call mpi_barrier(comm3d,ierr)
#else
         if (cvierr /= 0) then
            print '(A,A,i2)',' impCVODE_driver: fcvsptfqmr',
     &           ' failed, ier =',cvierr 
            return
         end if
#endif
c           GMRES
      else
         call fcvspgmr(prectype, gstype, maxl, lintol, cvierr)
#ifdef PARALLEL
         call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (cvierrtot /= 0) then
            if(outnode)  print '(A,A)',
     &           ' impCVODE_driver: fcvspgmr failed'
            if (cvierr /= 0) 
     &           print '(A,i4,A,i2)','    proc ',my_id,'
     &           : ier =',cvierr
            return
         end if
         call mpi_barrier(comm3d,ierr)
#else
         if (cvierr /= 0) then
            print '(A,A,i2)',' impCVODE_driver: fcvspgmr',
     &           ' failed, ier =',cvierr 
            return
         end if
#endif
      endif


c        J*v routine specification
      call fcvspilssetjac(useMyJv, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(outnode)
     &        print '(A,A)',' impCVODE_driver: fcvspilssetjac failed'
         if (cvierr /= 0) 
     &        print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVODE_driver: fcvspilssetjac',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        preconditioner solver specification
      call fcvspilssetprec(useMyPrec, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce(cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(outnode)
     &        print '(A,A)',' impCVODE_driver: fcvspilssetprec failed'
         if (cvierr /= 0) 
     &        print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVODE_driver: fcvspilssetprec',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif


c     call custom preconditioner initialization routines
      call combo_prec_init(cvierr)
#ifdef PARALLEL
      call MPI_Allreduce(cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(outnode) print *,' impCVODE_driver: combo_prec_init failed'
         if (cvierr /= 0) print '(A,i4,A,i2)','    proc ',my_id,
     &        ': ier =',cvierr
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVODE_driver: combo_prec_init failed, ',
     &        ' ier =',cvierr 
         return
      end if
#endif

      return
      end subroutine CvodeInit
c=================================================================




      subroutine CvodeFree()
c-----------------------------------------------------------------------
c     Description: CvodeFree deallocates all dynamic solver memory.
c-------------------------------------------------------------------
c======= Inclusions ===========
      use cvode_data

c======= Declarations =========
      implicit none

c======= Internals ============


c     free preconditioner data
      call combo_prec_free()

c     free cvode internal data
      call fcvfree

c     if using dynamic meshing, deallocate cvode_data structures
#ifdef DYNAMIC
      deallocate(vv,finv)
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
      end subroutine CvodeFree
c=================================================================




      subroutine fcvjtimes(vv, Jv, t, uu, fu, h, ipar, rpar, wk, cvier)
c-----------------------------------------------------------------------
c     Description: fcvjtimes provides the Jacobian-vector product
c        routine for the linearized Newton system, i.e. it calculates 
c        Jv = J*vv, where J is [approximately] equal to dF/du.
c     
c     Note: here this is a dummy routine, provided only to allow for 
c        compilation.  It is never called.
c
c     Arguments:
c            vv - (dbl(*), input) vector to multiply
c            Jv - (dbl(*), output) result of Jac*vv operation
c             t - (dbl, input) current time
c            uu - (dbl(*), input) current state for system
c            fu - (dbl(*), input) current rhs for ODE, i.e. f(u,t)
c             h - (dbl, input) step size for the last internal step
c            wk - (dbl(*)) free work space (size of uu)
c         cvier - (int, output) return flag: success(0), failure(1)
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: vv, Jv, uu, fu, wk
      double precision :: t, h, rpar(*)
      integer          :: ipar(*), cvier


c======= Internals ============

c     No matrix-vector product implemented.  Return failure if called
      cvier = 1

      return
      end subroutine fcvjtimes
c=======================================================================




      subroutine fcvpset(t, uu, fu, jok, jcur, gamdt, 
     &     h, ipar, rpar, v1, v2, v3, cvier)
c-----------------------------------------------------------------------
c     Description: fcvpset provides the preconditioner setup 
c        routine for the Newton system.
c
c     Arguments:
c             t - (dbl, input) current time
c            uu - (dbl(*), input) current state of system
c            fu - (dbl(*), input) rhs of ODE, i.e. f(u,t)
c           jok - (int, input) flag denoting whether to recompute 
c                  Jacobian-related data: 0=>recompute, 1=>unnecessary
c          jcur - (int, output) output flag to say if Jacobian data 
c                  was recomputed: 1=>was recomputed, 0=>was not
c         gamdt - (dbl, input) the scalar appearing in the Newton matrix
c                  M = I-gamdt*J
c             h - (dbl, input) time step size for last internal step
c            v1 - (dbl(*)) scratch space with same size as uu
c            v2 - (dbl(*)) scratch space with same size as uu
c            v3 - (dbl(*)) scratch space with same size as uu
c          cvier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c
c     Note: all of the data required to store the preconditioner setup 
c        information is handled through Fortran modules, since the data 
c        cannot be sent through the interface.
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use cvode_data
      use profiling

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &                 :: uu, fu, v1, v2, v3
      double precision :: t, gamdt, h, rpar(*)
      integer          :: jok, jcur, cvier, ipar(*)
      
c     internal variables
      double precision :: Re, Pr, Lu
      integer          :: ier

c======= Internals ============

c     initialize output flag to success
      cvier = 0

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

c     set jcur flag
      jcur = 0

c     if required (jok), call preconditioning setup routine
      if (jok == 0) then

c        set up preconditioners
         call combo_prec_setup(uu, fu, gamdt, v1, v2, cvier)

c        update jcur flag
         jcur = 1
      endif
      
c     stop the preconditioner setup profiler
      call prof_stop(4)

      return
      end subroutine fcvpset
c=======================================================================




      subroutine fcvpsol(t, uu, fu, rr, zz, gamdt, delta, lr, ipar, 
     &     rpar, vt, cvier)
c-----------------------------------------------------------------------
c     Description: fcvpsol provides the preconditioner solve 
c        routine for the preconditioning of the Newton system.
c             i.e. solves P*zz = rr
c
c     Arguments:
c             t - (dbl, input) current time
c            uu - (dbl(*), input) current state of system
c            fu - (dbl(*), input) rhs of ODE, i.e. f(u,t)
c            vt - (dbl(*), input) scratch space with same size as uu
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I-gamdt*J
c         delta - (dbl, input) desired tolerance for use in an
c                  iterative method.  In that case, solve until 
c                  Sqrt[Sum((r-Pz).*ewt)^2] < delta
c            rr - (dbl(*), input) rhs vector of prec. system 
c            lr - (int, input) flag indicating left or right 
c                  preconditioning: lr=1 => left,  lr=2 => right
c            zz - (dbl(*), output) solution vector of prec. system
c          cvier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c
c     Note: all of the data required to store the preconditioner 
c        information is handled through Fortran modules, since the 
c        data cannot be sent through the interface.
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use cvode_data
      use profiling

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, fu, vt, rr, zz
      double precision :: t, gamdt, delta, rpar(*)
      integer          :: lr, cvier, ipar(*)

c     local variables
      integer :: i, j, k, l

c======= Internals ============

c     initialize output flag to success
      cvier = 0

c     start the preconditioner solve profiler
      call prof_start(5)

c     reduce delta by a reasonable amount
      delta = delta*1.d-4

c     copy rhs into solution so that solver only uses one vector
      zz = rr

c     call preconditioner solvers
      call combo_prec_solve(uu, fu, gamdt, delta, zz, vt, cvier)

c     re-set vt to ameliorate cvode bug
      vt = 0.d0

c     stop the preconditioner solve profiler
      call prof_stop(5)

      return
      end subroutine fcvpsol
c=======================================================================




      subroutine fcverr(cvierr)
c-----------------------------------------------------------------------
c     Description: subroutine to check cvode error flags
c
c     Arguments:
c          cvierr - (in) conserved variables before/after filter
c          fcverr - (out) flag denoting success/failure of cvode solver
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      integer :: cvierr, cvierrmin, cvierrout

c======= Internals ============

c     collect error values for everyone to return
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrmin, 1, MPI_INTEGER, MPI_MIN, 
     &     comm3d, ierr)
      if ( ierr /= 0 ) then
         write(6,'(A,i4,A,i2)')'fcverr: mpi_allreduce failed, my_id =',
     &        my_id,',  ierr =', ierr
         cvierrmin = -1
         return
      end if
#else
      cvierrmin = cvierr
#endif

c     error if failure anywhere
      if ( cvierrmin < 0 ) then
         cvierrout = 1

#ifdef PARALLEL
         write(6,'(A,i4,A)')' CVODE ERROR, my_id =',my_id,':'
#else
         write(6,*)' CVODE ERROR:'
#endif

c        output relevant error message
         select case (cvierr)
         case(-1)
            write(6,*)'      [not enough memory]\n'
         case(-2)
            write(6,*)'      [problem in malloc routine]\n'
         case(-3)
            write(6,*)'      [illegal input, see docs]\n'
         case(-4)
            write(6,*)'      [too many internal steps]\n'
         case(-5)
            write(6,*)'      [could not satisfy accuracy req.]\n'
         case(-6)
            write(6,*)'      [exceed max error test failures]\n'
         case(-7)
            write(6,*)'      [exceed max conv. test failures]\n'
         case(-8)
            write(6,*)'      [linear solver setup failure]\n'
         case(-9)
            write(6,*)'      [linear solver solve failure]\n'
         end select
      else
         cvierrout = 0
      end if

c     set final output
      cvierr = cvierrout


      return
      end subroutine fcverr
c=================================================================





      subroutine fcvdiags(iopt,ropt)
c-----------------------------------------------------------------------
c     Description: subroutine to output cvode diagnostics
c
c     Arguments:
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

c======= Internals ============

c     root outputs CVODE statistics
      if (iprocx==1 .and. iprocy==1 .and. iprocz==1) then
c        general output values
         write(6,*) ' CVODE output values:'
         write(6,'(4x,A,i9)')'Total internal steps taken =',iopt(3)
         write(6,'(4x,A,i9)')'Total rhs function calls   =',iopt(4)
         write(6,'(4x,A,i9)')'Total nonlinear iterations =',iopt(7)
         write(6,'(4x,A,i9)')'Total linear iterations    =',iopt(20)

c        custom routine call information
         if (iopt(8)+iopt(18)+iopt(19) > 0) then
            write(6,'(4x,A,i9)')'Num lin solver setup calls =',iopt(8)
            write(6,'(4x,A,i9)')'Num PSet routine calls     =',iopt(18)
            write(6,'(4x,A,i9)')'Num PSolve routine calls   =',iopt(19)
         end if

c        error information
         if ( iopt(6) > 0 ) then
            write(6,'(4x,A,i9)')'Num nonlin conv failures   =',iopt(6)
         end if
         if ( iopt(5) > 0 ) then
            write(6,'(4x,A,i9)')'Num error test failures    =',iopt(5)
         end if
         if ( iopt(21) > 0 ) then
            write(6,'(4x,A,i9)')'Num linear conv failures   =',iopt(21)
         end if

c        general time-stepping information
         write(6,'(4x,A,i9)')'BDF order on last timestep =',iopt(9)
         write(6,'(4x,A,i9)')'BDF order on next timestep =',iopt(10)
         write(6,'(4x,A,i9)')'Num BDF order reductions   =',iopt(11)
         write(6,'(4x,A,es12.5)')'Suggested tol scale factor  =',ropt(5)
         write(6,'(4x,A,es12.5)')'Last internal timestep size =',ropt(2)
         write(6,'(4x,A,es12.5)')'Next internal timestep size =',ropt(3)
         write(6,'(4x,A,es12.5)')'Current internal time       =',ropt(4)
         write(6,*)'   '
      end if

      return
      end subroutine fcvdiags
c=================================================================





      subroutine divFlux(uu,dFlux)
c-----------------------------------------------------------------------
c     Description: divFlux provides a centered-difference spatial 
c        approximation to the divergence of the flux for the 
c        conservation equation describing the nonlinearly-resistive,
c        ideal MHD system.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at new time
c        dFlux - (dbl(*), output) divergence of the flux
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use cvode_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &     :: uu, dFlux

c     differencing factors
      double precision :: xfac, yfac, zfac, dt

c     indexes for flux calculations
      integer :: indx, indy, indz, i, j, k, l, m
      


c======= Internals ============


c     set differencing factors
      xfac = 1.0d0/dx;  yfac = 1.0d0/dy;  zfac = 1.0d0/dz

c     set dt value (unused, but passed to flux routines anyway)
      dt = 0.0d0

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


c------------ x-directional fluxes ------------
c     Hyperbolic fluxes
      inlo = ixlo;  inhi = ixhi
      call InviscidFluxWrapper(finv,vv,1,dt)

c     Parabolic fluxes
      inlo = ixlo;  inhi = ixhi
      indx = 1;  indy = 0;  indz = 0
      call ViscousFluxWrapper(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)


c------------ y-directional fluxes ------------
#ifndef ONE_D
c     Hyperbolic fluxes
      inlo = iylo;  inhi = iyhi
      call InviscidFluxWrapper(ginv,vv,2,dt)

c     Parabolic fluxes
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)


c------------ z-directional fluxes ------------
#ifndef TWO_D
c     Hyperbolic fluxes
      inlo = izlo;  inhi = izhi
      call InviscidFluxWrapper(hinv,vv,3,dt)

c     Parabolic fluxes
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif


c------------ construction of div(flux) ------------
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





      subroutine divFlux_viscous(uu,dFlux)
c-----------------------------------------------------------------------
c     Description: divFlux_viscous provides a centered-difference 
c        spatial approximation to the divergence of the flux for the 
c        viscous portion of the conservation equation describing the 
c        single-fluid resistive MHD system.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at new time
c        dFlux - (dbl(*), output) divergence of the flux
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use cvode_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &     :: uu, dFlux

c     differencing factors
      double precision :: xfac, yfac, zfac, dt

c     indexes for flux calculations
      integer :: indx, indy, indz, i, j, k, l, m
      


c======= Internals ============

#ifndef NOVISCFLUX

c     set differencing factors
      xfac = 1.0d0/dx;  yfac = 1.0d0/dy;  zfac = 1.0d0/dz

c     set dt value (unused, but passed to flux routines anyway)
      dt = 0.0d0

c     construct the vector of primitive variables
      call ConservedToPrimitive(uu,vv)

c     zero the fluxes (necessary???)
      fvis = 0.0d0
#ifndef ONE_D
      gvis = 0.0d0
#ifndef TWO_D
      hvis = 0.0d0
#endif
#endif

c     get the viscous flux requirements dudx and dbdx
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


c------------ x-directional fluxes ------------
c     Parabolic fluxes
      inlo = ixlo;  inhi = ixhi
      indx = 1;  indy = 0;  indz = 0
      call ViscousFluxWrapper(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)


c------------ y-directional fluxes ------------
#ifndef ONE_D
c     Parabolic fluxes
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)


c------------ z-directional fluxes ------------
#ifndef TWO_D
c     Parabolic fluxes
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif


c------------ construction of div(flux) ------------
c
c     Density
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,1) = 0.d0
            end do
         end do
      end do
c
c     X-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,2) = 
     &               xfac*( fvis(i,j,k,2) - fvis(i-1,j,k,2) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,4) - gvis(i,j-1,k,4) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,3) - hvis(i,j,k-1,3) )
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
               dFlux(i,j,k,3) =
     &               xfac*( fvis(i,j,k,3) - fvis(i-1,j,k,3) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,2) - gvis(i,j-1,k,2) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,4) - hvis(i,j,k-1,4) )
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
               dFlux(i,j,k,4) = 
     &               xfac*( fvis(i,j,k,4) - fvis(i-1,j,k,4) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,3) - gvis(i,j-1,k,3) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,2) - hvis(i,j,k-1,2) )
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
               dFlux(i,j,k,5) = 
     &               xfac*( fvis(i,j,k,5) - fvis(i-1,j,k,5) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,7) - gvis(i,j-1,k,7) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,6) - hvis(i,j,k-1,6) )
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
               dFlux(i,j,k,6) = 
     &               xfac*( fvis(i,j,k,6) - fvis(i-1,j,k,6) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,5) - gvis(i,j-1,k,5) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,7) - hvis(i,j,k-1,7) )
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
               dFlux(i,j,k,7) = 
     &               xfac*( fvis(i,j,k,7) - fvis(i-1,j,k,7) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,6) - gvis(i,j-1,k,6) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,5) - hvis(i,j,k-1,5) )
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
               dFlux(i,j,k,8) = 
     &               xfac*( fvis(i,j,k,8) - fvis(i-1,j,k,8) )
#ifndef ONE_D
     &              +yfac*( gvis(i,j,k,8) - gvis(i,j-1,k,8) )
#ifndef TWO_D
     &              +zfac*( hvis(i,j,k,8) - hvis(i,j,k-1,8) )
#endif
#endif
            end do
         end do
      end do

#endif

      return
      end subroutine divFlux_viscous
c=======================================================================

