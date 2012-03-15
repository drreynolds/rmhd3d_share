c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     SMU Mathematics
c
c     This file implements a solver for the viscous preconditioning 
c     components, which calls multigrid solver modules from the C 
c     routines defined in the file viscous_prec.h.  These routines
c     handle the block-diagonal solves in the viscous system.  The
c     lower-triangular block products are computed via a standard 
c     finite-difference Jacobian-vector product, albeit used in a 
c     non-standard way.  For discussion of the solver, see the 
c     paper.
c=======================================================================


      module combination_prec
c-----------------------------------------------------------------
c     Description: holds information on how the combination 
c        preconditioning should be applied
c
c     Contains the internal variables:
c         hprec - use the hyperbolic (fastwave) preconditioner
c         vprec - use the viscous (HYPRE) preconditioner
c        method - choice of how to combine preconditioners 
c         tmp1, tmp2 - temporary arrays used with combination 
c                preconditioning methods 3 and 4
c-----------------------------------------------------------------
      use mesh
      implicit none
      save

      integer :: hprec=0
      integer :: vprec=0
      integer :: method=0
#ifdef DYNAMIC
      real*8, dimension(:,:,:,:), allocatable :: tmp1, tmp2
#else
      real*8, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) :: tmp1,tmp2
#endif

      end module combination_prec
c=================================================================



      subroutine combo_prec_init(kier)
c-----------------------------------------------------------------------
c     Description: combo_prec_init inputs the relevant input files to 
c        determine which type of preconditining (if any) to use, and 
c        then calls the preconditioner-specific initialization routines.
c     Output: kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use fastwaveprec_data
      use sundials_fcmix
      use boundary_conds
      use combination_prec
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      integer :: kier, kiertot, ADIcorrection
      logical :: outnode, iofile
      integer(kind=sund_long_int) :: xlocal, ylocal, zlocal, nspecies
      integer(kind=sund_long_int) :: xghost, yghost, zghost
      namelist /prec/   hprec, vprec, method
      namelist /fwprec/ NWaves, ordering, ADIcorrection 


c======= Internals ============

      outnode = (iprocx==1 .and. iprocy==1 .and. iprocz==1)

c     input preconditioner options 
      inquire(file='prec.inp',exist=iofile) 
      if (iofile) then 
         open(unit=201,file='prec.inp',status='old') 
         read(201,nml=prec) 
         if (ADIcorrection == 0) then 
            ADIcorr = .false. 
         endif 
         if ((Nwaves<0) .or. (Nwaves>8)) then 
            Nwaves = 8 
         endif 
         if ((ordering<0) .or. (ordering>2)) then 
            ordering = 0 
         endif 

c        if using fastwave preconditioner, input options
         if (hprec == 1) then
            read(201,nml=fwprec) 
            if (ADIcorrection == 0) then 
               ADIcorr = .false. 
            endif 
            if ((Nwaves<0) .or. (Nwaves>8)) then 
               Nwaves = 8 
            endif 
            if ((ordering<0) .or. (ordering>2)) then 
               ordering = 0 
            endif 
         endif
         close(unit=201,status='keep') 
      endif 

c     ouptut preconditioner configuration to stdout
      if (outnode) then
         write(*,*) 'Overall preconditioner options:' 
         write(*,'(A,i2)') '    hprec =',hprec
         write(*,'(A,i2)') '    vprec =',vprec
         write(*,'(A,i2)') '   method =',method
         if (hprec == 1) then
            write(*,*) '   '
            write(*,*) 'FastWave Preconditioner options:' 
            write(*,'(A,i2)') '    NWaves =',NWaves 
            write(*,'(A,i2)') '    ordering =',ordering 
            write(*,*) '   ADIcorr =',ADIcorr 
         endif
      endif

c     initialize viscous preconditioner components (if used)
      if (vprec == 1) then
         xlocal = nxlsize
         ylocal = nylsize
         zlocal = nzlsize
         nspecies = nvar
         xghost = nghost;  yghost = 0;  zghost = 0
#ifndef ONE_D
         yghost = nghost
#ifndef TWO_D 
         zghost = nghost
#endif
#endif

#ifdef PARALLEL
c        viscous preconditioner for momentum equation
         call FVPrecDuInit(xlocal, ylocal, zlocal, dx, dy, dz, 
     &        nspecies, xghost, yghost, zghost, Xprocs, Yprocs, 
     &        Zprocs, iprocx, iprocy, iprocz, left, right, top, 
     &        bottom, forward, behind, xbc, ybc, zbc, kier)
         call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (kiertot /= 0) then
            if (outnode) 
     &           print '(A,A)',' combo_prec_init: FVPrecDuInit failed'
            if (kier /= 0) 
     &           print '(A,i4,A,i2)','    proc ',my_id,': ier =',kier
            return
         endif
         call mpi_barrier(comm3d,ierr)

c        resistive preconditioner for magnetic field equation
         call FVPrecDbInit(xlocal, ylocal, zlocal, dx, dy, dz, 
     &        nspecies, xghost, yghost, zghost, Xprocs, Yprocs, 
     &        Zprocs, iprocx, iprocy, iprocz, left, right, top, 
     &        bottom, forward, behind, xbc, ybc, zbc, kier)
         call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (kiertot /= 0) then
            if (outnode) 
     &           print '(A,A)',' combo_prec_init: FVPrecDbInit failed'
            if (kier /= 0) 
     &           print '(A,i4,A,i2)','    proc ',my_id,': ier =',kier
            return
         endif
         call mpi_barrier(comm3d,ierr)
         
c        heat conduction preconditioner for total energy equation
         call FVPrecDeInit(xlocal, ylocal, zlocal, dx, dy, dz, 
     &        nspecies, xghost, yghost, zghost, Xprocs, Yprocs, 
     &        Zprocs, iprocx, iprocy, iprocz, left, right, top, 
     &        bottom, forward, behind, xbc, ybc, zbc, kier)
         call MPI_Allreduce( kier, kiertot, 1, MPI_INTEGER, MPI_SUM, 
     &        comm3d, ierr)
         if (kiertot /= 0) then
            if (outnode) 
     &           print '(A,A)',' combo_prec_init: FVPrecDeInit failed'
            if (kier /= 0) 
     &           print '(A,i4,A,i2)','    proc ',my_id,': ier =',kier
            return
         endif
         call mpi_barrier(comm3d,ierr)
         
#else

c        viscous preconditioner for momentum equation
         call FVPrecDuInit(xlocal, ylocal, zlocal, dx, dy, dz, 
     &        nspecies, xghost, yghost, zghost, Xprocs, Yprocs, 
     &        Zprocs, 0, 0, 0, 0, 0, 0, 0, 0, 0, xbc, ybc, zbc, kier)
         if (kier /= 0) then
            print '(A,A,i2)',' combo_prec_init: FVPrecDuInit',
     &           ' failed, ier =',kier 
            return
         endif
         
c        resistive preconditioner for magnetic field equation
         call FVPrecDbInit(xlocal, ylocal, zlocal, dx, dy, dz, 
     &        nspecies, xghost, yghost, zghost, Xprocs, Yprocs, 
     &        Zprocs, 0, 0, 0, 0, 0, 0, 0, 0, 0, xbc, ybc, zbc, kier)
         if (kier /= 0) then
            print '(A,A,i2)',' combo_prec_init: FVPrecDbInit',
     &           ' failed, ier =',kier 
            return
         endif
         
c        heat conduction preconditioner for total energy equation
         call FVPrecDeInit(xlocal, ylocal, zlocal, dx, dy, dz, 
     &        nspecies, xghost, yghost, zghost, Xprocs, Yprocs, 
     &        Zprocs, 0, 0, 0, 0, 0, 0, 0, 0, 0, xbc, ybc, zbc, kier)
         if (kier /= 0) then
            print '(A,A,i2)',' combo_prec_init: FVPrecDeInit',
     &           ' failed, ier =',kier 
            return
         endif
#endif
      endif


c     initialize local temporary array storage (if needed)
#ifdef DYNAMIC
      if ((hprec*vprec == 1) .and. (method > 2)) then
         if (.not. allocated(tmp1))
     &        allocate(tmp1(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar))
         if (.not. allocated(tmp2))
     &        allocate(tmp2(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar))
      endif
#endif



      return
      end subroutine combo_prec_init
c=======================================================================




      subroutine combo_prec_free()
c-----------------------------------------------------------------------
c     Description: combo_prec_free deallocates all memory associated 
c        with the preconditioners.
c-----------------------------------------------------------------------
      use combination_prec

c======= Declarations =========
      implicit none

c======= Internals ============

c     free preconditioner data
      if (vprec == 1) then
         call FVPrecDuFree
         call FVPrecDbFree
         call FVPrecDeFree
      endif
      if (hprec == 1) then
         call fastwave_free
      endif

c     free local temporary array storage (if used)
#ifdef DYNAMIC
      if (allocated(tmp1))  deallocate(tmp1)
      if (allocated(tmp2))  deallocate(tmp2)
#endif

      return
      end subroutine combo_prec_free
c=======================================================================




      subroutine combo_prec_setup(uu, fu, gamdt, v1, v2, kier)
c-----------------------------------------------------------------------
c     Description: combo_prec_setup updates the relevant preconditioner
c        setup routines using the current solution.
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c            fu - (dbl(*), input) current residual or rhs f(uu)
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            v1 - (dbl(*)) scratch space with same size as uu
c            v2 - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use combination_prec

c======= Declarations =========
      implicit none

c     calling variables
      real*8, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &        :: uu, fu, v1, v2
      real*8  :: gamdt
      integer :: kier

c     local variables
      real*8 :: Re, Pr, Lu

c======= Internals ============


c     update hyperbolic preconditioner, if required
      if (hprec == 1) then
         call fastwave_pset(uu, fu, gamdt, v1, v2, kier)
      endif
      
c     update viscous preconditioner, if required
      if (vprec == 1) then

c        set non-dimensionalization parameters to 1 (unused)
         Re = 1.D0
         Pr = 1.D0
         Lu = 1.D0

         call fvprecdusetup(uu, gamdt, mu, Re, v1, v2, kier)
         if (kier /= 0) then
            write(0,*) '    fvprecdusetup error =', kier
            return
         endif
         
         call fvprecdbsetup(uu, gamdt, Lu, eta, v1, v2, kier)
         if (kier /= 0) then
            write(0,*) '    fvprecdbsetup error =', kier
            return
         endif
         
         call fvprecdesetup(uu, gamdt, gamma, kappa, Re, Pr, rgas, 
     &        v1, v2, kier)
         if (kier /= 0) then
            write(0,*) '    fvprecdesetup error =', kier
            return
         endif
      endif         

      return
      end subroutine combo_prec_setup
c=======================================================================




      subroutine combo_prec_solve(uu, fu, gamdt, delta, ww, tmp3, kier)
c-----------------------------------------------------------------------
c     Description: combo_prec_solve calls the preconditioner solution
c        routines to perform the combination preconditioning comprised 
c        of both hyperbolic and parabolic preconditioning components.  
c        This operates on the system  P z = r.
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c            fu - (dbl(*), input) current residual or rhs, i.e. f(uu)
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c         delta - (dbl, input) desired tolerance for precond. solve
c                  0.0 results in using HYPRE default
c            ww - (dbl(*), in/out) on input the rhs vector r, and on 
c                  output the solution vector z
c          tmp3 - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use combination_prec

c======= Declarations =========
      implicit none

c     calling variables
      real*8, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &        :: uu, fu, ww, tmp3
      real*8  :: gamdt, delta
      integer :: kier

c     local variables
      real*8  :: sigma

c======= Internals ============

c==== Mode 1: no/single preconditioner ====

c     call the hyperbolic preconditioner, if required
      if (hprec == 1 .and. vprec == 0) then
         call fastwave_psol(uu, fu, gamdt, ww, tmp3, kier)
      endif

c     call viscous preconditioner, if required
      if (vprec == 1 .and. hprec == 0) then
         call vprec_solver(uu, gamdt, delta, ww, tmp3, kier)
      endif

c==== Mode 2: combination preconditioner ====

      if (hprec*vprec == 1) then

c        choose approach based on "method" input parameter
         select case (method)

c        method 4: predictor-corrector preconditoner, v then h
         case (4)

c           store ww
            tmp1 = ww

c           first solve with viscous preconditioner
            call vprec_solver(uu, gamdt, delta, ww, tmp3, kier)

c           approximate tmp2 = tmp1 - J*w (the resulting linear residual), 
c           where J is full Jacobian, via FD
            sigma = 1.d-8
            tmp3 = uu + sigma*ww
            call divFlux(tmp3, tmp2)
            call divFlux(uu, tmp3)
            tmp2 = tmp1 + gamdt/sigma*(tmp3 - tmp2) - ww 

c           now solve for correction with hyperbolic preconditioner
            call fastwave_psol(uu, fu, gamdt, tmp2, tmp3, kier)
            ww = ww + tmp2

c        method 3: predictor-corrector preconditoner, h then v
         case (3)

c           store ww
            tmp1 = ww

c           first solve with hyperbolic preconditioner
            call fastwave_psol(uu, fu, gamdt, ww, tmp3, kier)

c           approximate tmp2 = tmp1 - J*w (the resulting linear residual), 
c           where J is full Jacobian, via FD
            sigma = 1.d-8
            tmp3 = uu + sigma*ww
            call divFlux(tmp3, tmp2)
            call divFlux(uu, tmp3)
            tmp2 = tmp1 + gamdt/sigma*(tmp3 - tmp2) - ww 

c           now solve for correction with viscous preconditioner
            call vprec_solver(uu, gamdt, delta, tmp2, tmp3, kier)
            ww = ww + tmp2

c        method 2: standard operator-split preconditioner v then h
         case (2)
            call vprec_solver(uu, gamdt, delta, ww, tmp3, kier)
            call fastwave_psol(uu, fu, gamdt, ww, tmp3, kier)

c        method 1 (default): standard op-split preconditioner h then v
         case default
            call fastwave_psol(uu, fu, gamdt, ww, tmp3, kier)
            call vprec_solver(uu, gamdt, delta, ww, tmp3, kier)

         end select

      endif
         
      return
      end subroutine combo_prec_solve
c=======================================================================
