c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
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
      use iounits
      use properties
      use sundials_fcmix
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      implicit none

c     local driver variables
      double precision :: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: x(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: b_mv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: b_fd(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: tmp1(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: tmp2(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      double precision :: delta
      double precision :: ttot, told, dt, dtfac, tdesired, tstop
      double precision :: Re, Lu, Pr, sigma, rmsdiff, infdiff, pi
      integer          :: istart, maxiter, lastiter
      integer          :: output_flag, binary_flag
      integer          :: ndump, ndiag, new
      integer          :: timestep, i, j, k, l

c     cvode timestepper/solver information
      double precision :: abstol(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: reltol, lintol, ropt(40), rdummy(1)
      integer          :: intmethod, itermethod, atoltype, taskindicator
      integer          :: prectype, gstype, maxl, cvierr
      integer          :: useMyJv, useMyPrec
      integer          :: cvierrtot, globalstrat
      integer          :: maxlrst, msbpre
      integer(kind=sund_long_int) :: nglobal, iopt(40), idummy(1)
      integer(kind=sund_long_int) :: xlocal, ylocal, zlocal, nspecies
      integer(kind=sund_long_int) :: xghost, yghost, zghost

c     variable output information
      double precision :: val(1:NXLsize-2,1:NYLsize-2)
      double precision :: dxval(1:NXLsize-3,1:NYLsize-2)
      double precision :: dyval(1:NXLsize-2,1:NYLsize-3)
      double precision :: dxxval(1:NXLsize-4,1:NYLsize-2)
      double precision :: dyyval(1:NXLsize-2,1:NYLsize-4)
      double precision :: dxyval(1:NXLsize-3,1:NYLsize-3)
      integer          :: gridsize, ix, iy, iz

c     declare internal namelists for problem specification
      namelist /inparams/ maxiter, new, output_flag, binary_flag
      namelist /dumps/ ndump, ndiag
      namelist /gridbounds/ xl, xr, yl, yr, zl, zr


c======= Implicit Scheme / Nonlinear Solver Parameters ============

c     set integration method
c     [1 => Adams method (nonstiff), 2 => BDF (stiff)]
      intmethod = 2

c     set nonlinear iteration method
c     [1 => functional iteration, 2 => Newton iteration]
      itermethod = 2

c     set type of absolute tolerance
c     [1 => scalar, 2 => array]
      atoltype = 2

c     set relative and absolute tolerances for nonlinear time step
c       * relative tolerance should be on the order of spatial accuracy
      reltol          = 1.0d-7    ! relative tolerance
c       * absolute tols should be floor values for each variable
c         (i.e. a measure of the allowable noise in solution)
      abstol(:,:,:,1) = 1.0d-11   ! rho absolute tolerance
      abstol(:,:,:,2) = 1.0d-11   ! ux absolute tolerance
      abstol(:,:,:,3) = 1.0d-11   ! uy absolute tolerance
      abstol(:,:,:,4) = 1.0d-11   ! uz absolute tolerance
      abstol(:,:,:,5) = 1.0d-11   ! Bx absolute tolerance
      abstol(:,:,:,6) = 1.0d-11   ! By absolute tolerance
      abstol(:,:,:,7) = 1.0d-11   ! Bz absolute tolerance
      abstol(:,:,:,8) = 1.0d-11   ! e absolute tolerance

c     set CVODE task indicator (NO DEFAULT)
c     [1 => normal mode (overshoot & interpolate),               ]
c     [2 => one-step mode (return after each internal step taken)]
c     [3 => normal tstop mode (like 1, but must stop at ropt(8)) ]
c     [4 => one-step tstop mode (like 2, but stops at ropt(8))   ]
c     [NOTE: modes 3 and 4 imply that ropt(8)/=0                 ]
      taskindicator = 1

c     set SPGMR options (0 => use default)
c        [preconditioning: 0=>none, 1=>left, 2=>right, 3=>both]
      prectype = 0
c        [1=>modified Gram-Schmidt, 2=>classical G-S; NO DEFAULT]
      gstype = 1
c        [maximum Krylov subspace dimension]
      maxl = 10
c        [ratio between linear & nonlinear tolerances]
      lintol = 0.0d0


c     set up time step multiplication factor (dt = dtfac)
      dtfac = 1.0d-3

c     set physical stopping time
      tstop = 2.5d3

c     use custom J*v, preconditioning routines (0=>default, 1=>custom)
      useMyJv   = 0   ! 1 => use my fcvjtimes routine
      useMyPrec = 0   ! 1 => use my fcvpset/fcvpsol routines

c     set up cvode solver options (0 => use default)
      iopt = 0;  ropt = 0.0d0
      iopt(1)  = 5       ! maximum order of BDF method (1-5, default: 5)
      iopt(2)  = 5000    ! max # internal steps to reach tdesired (def: 500)
      iopt(3)  = 0       ! max # warnings for  tn + h = tn  (default: 10)
      iopt(14) = 0       ! activate stability limit detection (0=off, 1=on)
      iopt(22) = 0       ! max # error test failures per step (default: 7)
      iopt(23) = 0       ! max # nonlin solver iters per step (default: 3)
      iopt(24) = 0       ! max # allowed conv fails per step (default: 10)
      ropt(1)  = 0.0d0   ! initial timestep size (default: estimated)
      ropt(2)  = 0.0d0   ! maximum absolute step size (default: infinity)
      ropt(3)  = 0.0d0   ! minimum absolute step size (default: 0.0)
      ropt(8)  = tstop   ! max time solution cannot pass (default: infinity)
      ropt(9)  = 0.0d0   ! coeff in nonlin conv test (default: 0.1)


c======= Internals ============


c     input general MHD simulation information using namelists
      open(16,file='mhd.inp',form='formatted')
      read(16,inparams)
      read(16,dumps)
      read(16,gridbounds) 
      close(16)

c     root node opens I/O streams and outputs problem parameters
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
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
         write(6,'(A,2(A,i4),3(A,i2))') '\nGrid indices:',
     &        '\n    nx =',nx,',  ny =',ny,',  nz =',nz,
     &        ',  nvar =',nvar,',  nghost =',nghost
         write(6,'(A,A,es9.2)') '\nImplicit Scheme:',
     &        '\n    Adaptive 5th-order BDF'
      endif

c     initialize timestepping and data output information
      lastiter = 0
      istart   = 1

c     set up domain, grid
      call SetupDomain
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter domain setup...'
         write(6,'(A,3(A,es12.5))') '\nMesh spacing:',
     &        '\n   dx =',dx,',  dy =',dy,',  dz =',dz
      end if

c     determine local domain location in general domain
      call SetupLocalDomain
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter local domain setup...'
         write(6,'(A,3(i2,1x),A,3(2(es9.2,1x),A))') 
     &        '\nLocal subdomain for iproc (',
     &        iprocx,iprocy,iprocz,') is the box\n      [',
     &        XLloc,XRloc,'] x [',YLloc,YRloc,'] x [',ZLloc,ZRloc,']'
      endif

c     initialize state variables
      call Initialize(new,ux,phi,dt,ttot,istart,lastiter)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,*) '\nAfter local dataspace initialization...'
      end if

c     communicate initial boundary data among processors
      call CommState(ux)

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


c     set up CVODE NVector interface information using mesh variables
      nglobal = (xlocal+2*xghost)*(ylocal+2*yghost)*(zlocal+2*zghost)
     &     *nvar

c     initialize CVODE data & operators

c        Nvector specs
      nglobal = (nxlsize+2*xghost)*(nylsize+2*yghost)*(nzlsize+2*zghost)
     &     *nvar
      call fnvinitmhd(1, nglobal, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A)',' test_driver: fnvinitmhd failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' test_driver: fnvinitmhd',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        cvode dataspace
c        [we must access abstol values before sending to fcvmalloc]
      if (maxval(abstol) > 0.d0)
     &     write(6,*) '\nSetting abstol vector'
      call fcvmalloc(ttot, ux, intmethod, itermethod, 2,
     &     reltol, abstol, iopt, ropt, idummy, rdummy, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' test_driver: fcvmalloc failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' test_driver: fcvmalloc',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        spgmr dataspace
      call fcvspgmr(prectype, gstype, maxl, lintol, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' test_driver: fcvspgmr failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' test_driver: fcvspgmr',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        J*v routine specification
      call fcvspgmrsetjac(useMyJv, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' test_driver: fcvspgmrsetjac failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' test_driver: fcvspgmrsetjac',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

c        preconditioner setup/solver specification
      call fcvspgmrsetprec(useMyPrec, cvierr)
#ifdef PARALLEL
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' test_driver: fcvspgmrsetprec failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      if (cvierr /= 0) then
         print '(A,A,i2)',' test_driver: fcvspgmrsetprec',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif


      if (my_id == 0) print *,' driver: calling FVPrecDuInit'
c        custom viscous momentum preconditioner initialization
#ifdef PARALLEL
      call FVPrecDuInit(nxlocal, nylocal, nzlocal, dx, dy, dz, 
     &     nvar, nghost, nghost, nghost, Xprocs, Yprocs, 
     &     Zprocs, iprocx, iprocy, iprocz, left, right, top, 
     &     bottom, forward, behind, cvierr)
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' impCVPRE_driver: FVPrecDuInit failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      call FVPrecDuInit(nxlocal, nylocal, nzlocal, dx, dy, dz, 
     &     nvar, nghost, nghost, nghost, Xprocs, Yprocs, 
     &     Zprocs, 0, 0, 0, 0, 0, 0, 0, 0, 0, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVPRE_driver: FVPrecDuInit',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

      if (my_id == 0) print *,' driver: calling FVPrecDbInit'
c        custom viscous magnetic field preconditioner initialization
#ifdef PARALLEL
      call FVPrecDbInit(nxlocal, nylocal, nzlocal, dx, dy, dz, 
     &     nvar, nghost, nghost, nghost, Xprocs, Yprocs, 
     &     Zprocs, iprocx, iprocy, iprocz, left, right, top, 
     &     bottom, forward, behind, cvierr)
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' impCVPRE_driver: FVPrecDbInit failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      call FVPrecDbInit(nxlocal, nylocal, nzlocal, dx, dy, dz, 
     &     nvar, nghost, nghost, nghost, Xprocs, Yprocs, 
     &     Zprocs, 0, 0, 0, 0, 0, 0, 0, 0, 0, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVPRE_driver: FVPrecDbInit',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif

      if (my_id == 0) print *,' driver: calling FVPrecDeInit'
c        custom viscous energy preconditioner initialization
#ifdef PARALLEL
      call FVPrecDeInit(nxlocal, nylocal, nzlocal, dx, dy, dz, 
     &     nvar, nghost, nghost, nghost, Xprocs, Yprocs, 
     &     Zprocs, iprocx, iprocy, iprocz, left, right, top, 
     &     bottom, forward, behind, cvierr)
      call MPI_Allreduce( cvierr, cvierrtot, 1, MPI_INTEGER, MPI_SUM, 
     &     comm3d, ierr)
      if (cvierrtot /= 0) then
         if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
            print '(A,A)',' impCVPRE_driver: FVPrecDeInit failed'
         end if
         if (cvierr /= 0) then
            print '(A,i4,A,i2)','    proc ',my_id,': ier =',cvierr
         endif
         return
      end if
      call mpi_barrier(comm3d,ierr)
#else
      call FVPrecDeInit(nxlocal, nylocal, nzlocal, dx, dy, dz, 
     &     nvar, nghost, nghost, nghost, Xprocs, Yprocs, 
     &     Zprocs, 0, 0, 0, 0, 0, 0, 0, 0, 0, cvierr)
      if (cvierr /= 0) then
         print '(A,A,i2)',' impCVPRE_driver: FVPrecDeInit',
     &        ' failed, ier =',cvierr 
         return
      end if
#endif



c     if parallel, wait until all processors are ready
#ifdef PARALLEL
      call mpi_barrier(comm3d, ierr)
#endif


cccc  FIRST, GET THE PRECONDITIONER MATRIX-VECTOR PRODUCT  cccc

c     initialize non-dimensionalization parameters to 1.d0
      Re = 1.d0
      Pr = 1.d0
      Lu = 1.d0

c     call viscous preconditioner matrix setup routines
      write(6,*) '    calling fvprecdusetup'
      call fvprecdusetup(ux, dtfac, mu, Re, tmp1, tmp2, cvierr)
      write(6,*) '    fvprecdusetup complete'
      if (cvierr /= 0) then
         write(0,*) '    fvprecdusetup error =', cvierr
      endif
      
      write(6,*) '    calling fvprecdbsetup'
      call fvprecdbsetup(ux, dtfac, Lu, eta, tmp1, tmp2, cvierr)
      write(6,*) '    fvprecdbsetup complete'
      if (cvierr /= 0) then
         write(0,*) '    fvprecdbsetup error =', cvierr
      endif
      
      write(6,*) '    calling fvprecdesetup'
      call fvprecdesetup(ux, dtfac, gamma, kappa, Re, Pr, rgas, 
     &     tmp1, tmp2, cvierr)
      write(6,*) '    fvprecdesetup complete'
      if (cvierr /= 0) then
         write(0,*) '    fvprecdesetup error =', cvierr
      endif

c     set up vector to multiply 
c$$$c        [set to ones for matrix row sums]
c$$$      x = 1.d0
c        [set to smoothly-varying field]
      pi = (4.d0*datan(1.d0))
      do k=1,nzlocal,1
         iz = ZLloc + k*dz*0.5d0
         do j=1,nylocal,1
            iy = YLloc + j*dy*0.5d0
            do i=1,nxlocal,1
               ix = XLloc + i*dx*0.5d0
                x(i,j,k,:) = (1.0d-1
     &              + 1.d0*(ix-xl)*(xr-ix)/(xr-xl)/(xr-xl)
     &              + 1.d0*(iy-yl)*(yr-iy)/(yr-yl)/(yr-yl)
     &              + 1.d0*(iz-zl)*(zr-iz)/(zr-zl)/(zr-zl))
     &              *(ux(i,j,k,:) + 1.d-3)
            enddo
         enddo
      enddo
      call CommState(x)

c     call vprec_mult to get Preconditioner matrix-vector product
      b_mv = x
      call vprec_mult_hypre(ux, dtfac, b_mv, tmp1, cvierr)
      if (cvierr /= 0) then
         write(0,*) '    vprec_mult error =',cvierr
      endif


cccc  SECOND, GET THE FINITE-DIFFERENCE MATRIX-VECTOR PRODUCT  cccc

c     set FD product differencing parameter
      sigma = 1.0d-8

c     get viscous flux divergence at ux -> tmp2
      tmp2 = 0.d0
      call CommState(ux)
      call divFlux_viscous(ux,tmp2)

c     set differencing state b_fd = ux + sigma*x
      b_fd = ux + sigma*x

c     get viscous flux divergence at ux + sigma*v -> tmp1
      tmp1 = 0.d0
      call CommState(b_fd)
      call divFlux_viscous(b_fd,tmp1)

c     set FD matrix-vector product
      b_fd = x - dtfac/sigma*(tmp1 - tmp2)


cccc  CHECK THAT THESE MATRIX-VECTOR PRODUCTS AGREE  cccc
      write(6,*) '  Results of viscous matrix-vector product tests:'
      write(6,*) '          dtfac = ',dtfac
      write(6,*) '          sigma = ',sigma
      
c     compute b_fd-b_mv -> tmp1
      tmp1 = abs(b_fd - b_mv)

c     compute (b_fd-b_mv)/b_fd -> tmp2
      tmp2 = abs(tmp1/b_fd)

c     output matvecs and difference to file
      write(6,*) '  Printing multiplier vec x to file mv_x.'
      open(1000,file='mv_x')
      do l=2,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(1000,'(3x,4(i4,1x),es20.13)') 
     &                 i, j, k, l, x(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      close(1000)

      write(6,*) '  Printing FD matvec to file mv_fd.'
      open(999,file='mv_fd')
      do l=2,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(999,'(3x,4(i4,1x),es20.13)') 
     &                 i, j, k, l, b_fd(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      close(999)

      write(6,*) '  Printing Preconditioner matvec to file mv_pc.'
      open(998,file='mv_pc') 
      do l=2,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(998,'(3x,4(i4,1x),es20.13)') 
     &                 i, j, k, l, b_mv(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      close(998)

      write(6,*) '  Printing absolute difference to file mvdiff_abs.'
      open(997,file='mvdiff_abs')
      do l=2,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(997,'(3x,4(i4,1x),es20.13)') 
     &                 i, j, k, l, tmp1(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      close(997)
c      pause

c     get absolute difference in rho values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,1)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,1))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_rho values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in ux values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,2)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,2))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_ux values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in uy values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,3)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,3))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_uy values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in uz values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,4)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,4))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_uz values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in Bx values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,5)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,5))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_Bx values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in By values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,6)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,6))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_By values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in Bz values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,7)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,7))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_Bz values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff

c     get absolute difference in e values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,8)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,8))
      write(6,*) '   '
      write(6,*) '  Absolute difference in Matvec_e values:'
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||b_fd - b_mv||_inf =',infdiff


c     get relative difference in rho values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,1)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,1))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_rho values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in ux values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,2)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,2))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_ux values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in uy values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,3)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,3))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_uy values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in uz values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,4)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,4))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_uz values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in Bx values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,5)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,5))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_Bx values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in By values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,6)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,6))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_By values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in Bz values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,7)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,7))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_Bz values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff

c     get relative difference in e values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,8)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,8))
      write(6,*) '   '
      write(6,*) '  Relative difference in Matvec_e values:'
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(b_fd - b_mv)/b_fd||_inf =',infdiff


cccc  Test that Pv^{-1}*Pv*x = x  cccc
      write(6,*) '   '
      write(6,*) '  Testing that Pv^{-1}*Pv*x = x:'
c      pause
     
c     set b_mv = Pv^{-1}*x
      b_mv = x

c     set b_mv = Pv*b_mv
      delta = 1.d-12
      call vprec_solver(ux, dtfac, delta, b_mv, tmp1, cvierr)

c     set b_mv = Pv*b_mv
      call vprec_mult(ux, dtfac, b_mv, tmp1, cvierr)

c     compute differences between b_mv and x
c     compute x-b_mv -> tmp1
      tmp1 = abs(x - b_mv)

c     compute (x-b_mv)/x -> tmp2
      tmp2 = abs(tmp1/x)

c     output matvecs and difference to file
      write(6,*) '  Printing Pv*Pv^{-1}*x file mv_ppx.'
      open(996,file='mv_ppx') 
      do l=2,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(996,'(3x,4(i4,1x),es20.13)') 
     &                 i, j, k, l, b_mv(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      close(996)

      write(6,*) '  Printing absolute difference to file ppdiff_abs.'
      open(995,file='ppdiff_abs')
      do l=2,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(995,'(3x,4(i4,1x),es20.13)') 
     &                 i, j, k, l, tmp1(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      close(995)

c     get absolute difference in rho values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,1)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,1))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_rho values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in ux values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,2)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,2))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_ux values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in uy values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,3)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,3))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_uy values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in uz values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,4)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,4))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_uz values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in Bx values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,5)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,5))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_Bx values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in By values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,6)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,6))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_By values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in Bz values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,7)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,7))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_Bz values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff

c     get absolute difference in e values
      rmsdiff = sqrt(sum(tmp1(1:nxlocal,1:nylocal,1:nzlocal,8)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp1(1:nxlocal,1:nylocal,1:nzlocal,8))
      write(6,*) '   '
      write(6,*) '  Absolute difference in P*Pi*x_e values:'
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||P*Pi*x - x||_inf =',infdiff


c     get relative difference in rho values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,1)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,1))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_rho values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in ux values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,2)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,2))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_ux values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in uy values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,3)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,3))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_uy values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in uz values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,4)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,4))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_uz values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in Bx values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,5)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,5))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_Bx values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in By values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,6)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,6))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_By values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in Bz values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,7)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,7))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_Bz values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff

c     get relative difference in e values
      rmsdiff = sqrt(sum(tmp2(1:nxlocal,1:nylocal,1:nzlocal,8)**2)
     &     / (nxlocal*nylocal*nzlocal))
      infdiff = maxval(tmp2(1:nxlocal,1:nylocal,1:nzlocal,8))
      write(6,*) '   '
      write(6,*) '  Relative difference in P*Pi*x_e values:'
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_rms =',rmsdiff
      write(6,'(A,es12.5)') '     ||(P*Pi*x - x)/x||_inf =',infdiff
      

      write(6,*) '   '
      write(6,*) '  Freeing solver memory'
c     free up solver dataspace memory
      call FVPrecDuFree
      call FVPrecDbFree
      call FVPrecDeFree
      call fcvfree

c     dump all remaining data, output diagnostics, close I/O files
      call WrapUp(ux,phi,dt,ttot,timestep)


      return
      end subroutine TestDriver
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
