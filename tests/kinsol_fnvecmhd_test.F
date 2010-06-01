      program kinsol_fnvecmhd_test
c=======================================================================c
c     File        : kinsol_fnvecmhd_test.F                              c
c     Programmers : Daniel R. Reynolds @ LLNL                           c
c     Version of  : 2 October 2003                                      c
c-----------------------------------------------------------------------c
c     This testing routine is used to ensure the interface of the       c
c     new N_Vector kernel with the Fortran Kinsol interface.            c
c                                                                       c
c     Uses the internals of a fortran kinsol example provided as part   c
c     of the SUNDIALS library.  Here, it performs a simple diagonal     c
c     test using user-supplied preconditioner setup and solve routines, c
c     solving the nonlinear system                                      c
c                       f(u) = u(i)^2 - i^2 = 0                         c
c     The solve does not incorporate scaling, and an approximate        c
c     diagonal preconditioner is used.                                  c
c                                                                       c
c     Execute line:  mpirun -np N kinsol_fnvecmhd_test.x                c
c         where N is the number of processors to use                    c
c=======================================================================c
c     all variables must be declared
      implicit none

c     modules 
      include "mpif.h"

c     mesh variables
      integer, parameter :: nx = 32
      integer, parameter :: ny = 1
      integer, parameter :: nz = 1
      integer, parameter :: ns = 1
      integer, parameter :: xghost = 0
      integer, parameter :: yghost = 0
      integer, parameter :: zghost = 0
      integer            :: ixlo = 1-xghost, ixhi = nx+xghost
      integer            :: iylo = 1-yghost, iyhi = ny+xghost
      integer            :: izlo = 1-zghost, izhi = nz+xghost
      integer            :: islo = 1,        ishi = ns

c     internal variables
      integer            :: i, ii, ier, size, globalstrat, rank
      integer            :: mype, nprocs, baseadd, nlocal, nglobal
      integer            :: iopt(40)
      double precision   :: fnormtol, scsteptol
      double precision   :: ropt(40)
      double precision   :: uu(nx,ny,nz,ns), scale(nx,ny,nz,ns)
      double precision   :: rp, rpscale, frpscale
      double precision   :: constr(nx,ny,nz,ns), rpconstr
      double precision   :: pp(nx,ny,nz,ns)

c     kinsol optimization parameters
      integer  ::  inopt   = 0
      integer  ::  maxl    = 10
      integer  ::  maxlrst = 2
      integer  ::  msbpre  = 5

      common /pcom/ pp, mype, nprocs, baseadd, ixlo, ixhi, iylo, iyhi, 
     &     izlo, izhi, islo, ishi

      globalstrat = 1
      fnormtol    = 1d-5
      scsteptol   = 1d-4
      nlocal  = nx*ny*nz*ns
      inopt   = 0
      maxl    = 10
      maxlrst = 2
      msbpre  = 5

c     make certain that everything's been compiled with the 'PARALLEL' 
c     compile-time variable
#ifndef PARALLEL
      print *, 'Kinsol_fnvecmhd_test Error: compile with -DPARALLEL'
#endif

c     The user MUST call mpi_init, Fortran binding, for the fkinsol package
c     to work. The communicator, MPI_COMM_WORLD, is the only one common 
c     between the Fortran and C bindings. So in the following, the communicator
c     MPI_COMM_WORLD is used in calls to mpi_comm_size and mpi_comm_rank
c     to determine the total number of processors and the rank (0 ... size-1) 
c     number of this process.
      call mpi_init(ier)
      if (ier /= 0) then
         print '(A,i2)', 'Kinsol_fnvecmhd_test: mpi_init failed, ier ='
     &        ,ier
         stop
      endif
      
c     get number of processors
      call mpi_comm_size(MPI_COMM_WORLD,size,ier)
      if (ier /= 0) then
         print '(A,i2)', 'Kinsol_fnvecmhd_test: ',
     &        'mpi_comm_size failed, ier =',ier
         stop
      endif      
      nprocs = size

c     get rank of local process
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ier)
      if (ier /= 0) then
         print '(A,i2)', 'Kinsol_fnvecmhd_test: ',
     &        'mpi_comm_rank failed, ier =',ier
         stop
      endif
      mype = rank
      baseadd = mype * nlocal 

c     initialize the new nvspec
      nglobal = nlocal*nprocs
      call fnvinitmhd(nx, ny, nz, xghost, yghost, zghost, ns, 
     &     nglobal, ier)
      if (ier /= 0) then
         print '(A,A,i2)','Kinsol_fnvecmhd_test: ',
     &        'fnvinitmhd failed, ier =', ier
         stop
      endif

c     initialize solution guess, scaling and constraint arrays
      do ii = 1,nlocal
         i = ii + baseadd
         uu(ii,1,1,1)     = 2.0d0*i
         scale(ii,1,1,1)  = 1.0d0
         constr(ii,1,1,1) = 0.0d0
      end do
      rp = 1.0d0
      rpscale = 1.0d0
      rpconstr = 0.0d0
      frpscale = 1.0d0
      
c     initialize kinsol data space
      call fkinmalloc(msbpre, fnormtol, scsteptol, constr, 
     &     rpconstr, inopt, iopt, ropt, ier)      
      if (ier /= 0) then
         print '(A,i2)', 'Kinsol_fnvecmhd_test: ',
     &        'fkinmalloc failed, ier =',ier
         stop
      endif

c     initialize sparse gmres dataspace, solver and preconditioner 
      call fkinspgmr(maxl, maxlrst, ier)
      call fkinspgmrsetjac(0, ier)
      call fkinspgmrsetpsol(1, ier)
      call fkinspgmrsetpset(1, ier)
      
c     output startup message
      if(mype == 0) then
         print *, '\nExample program kinsol_fnvecmhd_test'
         print *, '  This fkinsol/fnvec_mhd example code solves a'
         print *, '  n*32 equation diagonal algebraic system (here'
         print *, '  n is the number of processors).  Its purpose'
         print *, '  is to demonstrate the use of the fortran'
         print *, '  kinsol interface to the new nvec_mhd NVector'
         print *, '  kernel.\n'
      end if

c     use fkinsol to solve the nonlinear system
      call fkinsol(uu, rp, globalstrat, scale, 
     &     rpscale, scale, frpscale, ier)
      if (mype == 0) then
         print '(A,i5)', 'Kinsol_fnvecmhd_test: fkinsol return code is',
     &        ier
         print *, '\n  The resultant values on processor 0 are:\n'
         print '(4x,A,f10.6)', 'rp =',rp
         do i = 1,nlocal,4
            print '(4x,i4,4(1x,f10.6))', i+baseadd, uu(i,1,1,1), 
     &           uu(i+1,1,1,1), uu(i+2,1,1,1), uu(i+3,1,1,1)
         end do
         print *, '\n\n  The kinsol output values are:'
         print '(4x,6(A,i4))', 'nni =',iopt(4),',  nli =',iopt(11),
     &        ',  nfe =',iopt(5),',  npe =',iopt(12),',  nps =',
     &        iopt(13),',  ncfl =',iopt(14)
      end if

c     free up kinsol and nvspec data
      call fkinfree
      call fnvfreemhd
      
c     finalize mpi
      call mpi_finalize(ier)
      
      stop
      end program kinsol_fnvecmhd_test
c=======================================================================c
      




      subroutine fkfun(uu, rp, fval, frp)
c=======================================================================c
c     fkfun provides the nonlinear system function f(u).                c
c=======================================================================c
      implicit none

c     problem parameters
      integer, parameter :: nx = 32
      integer, parameter :: ny = 1
      integer, parameter :: nz = 1
      integer, parameter :: ns = 1
      integer, parameter :: xghost = 0
      integer, parameter :: yghost = 0
      integer, parameter :: zghost = 0
      
c     calling variables
      double precision, dimension( 1-xghost:nx+xghost, 
     &                             1-yghost:ny+yghost, 
     &                             1-zghost:nz+zghost, 1:ns) 
     &                   :: uu, fval
      double precision   :: rp, frp
      
c     internal, common variables
      double precision   :: pp(nx,ny,nz,ns)
      integer            :: i, j, k, l, mype, nprocs, baseadd
      integer            :: ixlo,ixhi,iylo,iyhi,izlo,izhi,islo,ishi
      common /pcom/ pp, mype, nprocs, baseadd, ixlo, ixhi, iylo, iyhi, 
     &     izlo, izhi, islo, ishi
      
      do l = islo,ishi
         do k = izlo,izhi
            do j = iylo,iyhi
               do i = ixlo,ixhi
                  fval(i,j,k,l) = uu(i,j,k,l)*uu(i,j,k,l) - 
     &                 (i+baseadd)*(i+baseadd)
               end do
            end do
         end do
      end do
      frp = (rp-3.0d0)**3 + 5.0d-1
      
      return
      end subroutine fkfun
c=======================================================================c
      
      
      


      subroutine fkpset(udata, urp, uscale, urpsc, fdata, frp, 
     &     fscale, frpsc, vtmp1, vtmp2, ier)
c=======================================================================c
c     fkpset provides the preconditioner setup routine.                 c
c=======================================================================c
      implicit none

c     problem parameters
      integer, parameter :: nx = 32
      integer, parameter :: ny = 1
      integer, parameter :: nz = 1
      integer, parameter :: ns = 1
      integer, parameter :: xghost = 0
      integer, parameter :: yghost = 0
      integer, parameter :: zghost = 0
      
c     calling variables
      double precision, dimension( 1-xghost:nx+xghost, 
     &                             1-yghost:ny+yghost, 
     &                             1-zghost:nz+zghost, 1:ns) 
     &                   :: udata, uscale, fdata, fscale
      integer            :: ier
      double precision   :: vtmp1(*), vtmp2(*)
      double precision   :: urp, urpsc, frp, frpsc
      
c     internal, common variables
      double precision   :: pp(nx,ny,nz,ns)
      integer            :: mype, nprocs, baseadd
      integer            :: ixlo,ixhi,iylo,iyhi,izlo,izhi,islo,ishi
      common /pcom/ pp, mype, nprocs, baseadd, ixlo, ixhi, iylo, iyhi, 
     &     izlo, izhi, islo, ishi

      pp = 5d-1 / (udata+5d0)
      ier = 0
      
      return
      end subroutine fkpset
c=======================================================================c
      
      
      

      
      subroutine fkpsol(udata, urp, uscale, urpsc, fdata, frp, 
     &     fscale, frpsc, vv, vrp, ftem, ier)
c=======================================================================c
c     fkpsol provides the preconditioner solve routine.                 c
c=======================================================================c
      implicit none

c     problem parameters
      integer, parameter :: nx = 32
      integer, parameter :: ny = 1
      integer, parameter :: nz = 1
      integer, parameter :: ns = 1
      integer, parameter :: xghost = 0
      integer, parameter :: yghost = 0
      integer, parameter :: zghost = 0
      
c     calling variables
      double precision, dimension( 1-xghost:nx+xghost, 
     &                             1-yghost:ny+yghost, 
     &                             1-zghost:nz+zghost, 1:ns) 
     &                   :: udata, uscale, fdata, fscale, vv
      integer            :: ier
      double precision   :: ftem(*), urp, urpsc, frp, frpsc, vrp
      
c     internal, common variables
      double precision   :: pp(nx,ny,nz,ns)
      integer            :: mype, nprocs, baseadd
      integer            :: ixlo,ixhi,iylo,iyhi,izlo,izhi,islo,ishi
      common /pcom/ pp, mype, nprocs, baseadd, ixlo, ixhi, iylo, iyhi, 
     &     izlo, izhi, islo, ishi
      
      vv = vv*pp
      ier = 0

      return
      end subroutine fkpsol
c=======================================================================c
      
      
      

      
      subroutine fkjtimes(ww, wrp, Jw, Jwrp, new_uu, uu, urp, ier)
c=======================================================================c
c     fkpsol provides the preconditioner solve routine.                 c
c=======================================================================c
      implicit none

c     problem parameters
      integer, parameter :: nx = 32
      integer, parameter :: ny = 1
      integer, parameter :: nz = 1
      integer, parameter :: ns = 1
      integer, parameter :: xghost = 0
      integer, parameter :: yghost = 0
      integer, parameter :: zghost = 0
      
c     calling variables
      double precision, dimension( 1-xghost:nx+xghost, 
     &                             1-yghost:ny+yghost, 
     &                             1-zghost:nz+zghost, 1:ns) 
     &                   :: ww, Jw, uu
      integer            :: new_uu, ier
      double precision   :: wrp, Jwrp, urp

c     dummy routine, so return as successful
      ier = 0

      return
      end subroutine fkjtimes
c=======================================================================c
