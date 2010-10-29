c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Daniel R. Reynolds
c     SMU Mathematics
c-----------------------------------------------------------------
c     $Log:  $
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
c     Description: test driver, handles initialization of 
c        variables, timestepping and data I/O.  This 
c        routine tests the accuracy of my new ptri_parallel
c        solver.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use iounits
      use mpistuff
      use ptri_solver
      use ptri_product
      
c======= Declarations =========
      implicit none

      real*8, dimension(nxlocal,nylocal,nzlocal) :: ux, utrue
      real*8  :: lerr, gerr
      integer :: dims, i, j, k, order(3)
      logical :: outnode

c     declare internal namelists for problem specification
      namelist /gridbounds/ xl, xr, yl, yr, zl, zr

c======= Internals ============


c     see if this is the output node
      outnode = .false.
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1)  outnode = .true.

c     input general MHD simulation information using namelists
      open(16,file='mhd.inp',form='formatted')
      read(16,gridbounds) 
      close(16)

c     root node outputs problem parameters
      if(outnode) then
         write(6,'(A,3(A,i4))') '\nGrid indices:',
     &        '\n    nx =',nx,',  ny =',ny,',  nz =',nz
         write(6,'(A,3(A,i4))') '\nParallelism information:',
     &        '\n    xprocs =',xprocs,',  yprocs =',yprocs,
     &        ',  zprocs =',zprocs
         write(6,*) '  '
      endif

c     set up domain, grid
      call SetupDomain

c     determine local domain location in general domain
      call SetupLocalDomain

c     wait for all processes to catch up
#ifdef PARALLEL
      call mpi_barrier(comm3d, ierr)
#endif

c     set up solver data
      if(outnode) 
     &     write(6,*) '\n Setting up ptri solver data'
         dims = 3
#ifdef ONE_D
         dims = 1
#endif
#ifdef TWO_D
         dims = 2
#endif
      call ptri_setup(comm3d,dims,nxlocal,nylocal,nzlocal,1,1,1,ierr)
      call ptri_mult_setup(comm3d,dims,nxlocal,nylocal,nzlocal,ierr)




c     call directionally-split solver for 1D test problems
c        create the true solution (put in utrue)
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
c$$$               utrue(i,j,k) = 1.d0*xl + dx*(i-1) + dy*(j-1) + dz*(k-1)
               utrue(i,j,k) = 1.d0*i
            end do
         end do      
      end do

c        x-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 1D tests (x)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 1, 0, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 1, 0, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr

c        y-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 1D tests (y)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 2, 0, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 2, 0, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        z-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 1D tests (z)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 3, 0, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 3, 0, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr






c     call directionally-split solver for 2D test problems

c        xy-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 2D tests (xy)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 2, 1, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 1, 2, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        xz-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 2D tests (xz)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 3, 1, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 1, 3, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        yz-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 2D tests (yz)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 3, 2, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 2, 3, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        yx-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 2D tests (yx)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 1, 2, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 2, 1, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        zx-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 2D tests (zx)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 1, 3, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 3, 1, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        zy-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 2D tests (zy)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 2, 3, 0 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 3, 2, 0 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr







c     call directionally-split solver for 3D test problems

c        xyz-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 3D tests (xyz)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 3, 2, 1 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 1, 2, 3 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        xzy-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 3D tests (xzy)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 2, 3, 1 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 1, 3, 2 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        yzx-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 3D tests (yzx)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 1, 3, 2 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 2, 3, 1 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        yxz-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 3D tests (yxz)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 3, 1, 2 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 2, 1, 3 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        zxy-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 3D tests (zxy)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 2, 1, 3 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 3, 1, 2 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr


c        zyx-direction
      if(outnode) 
     &     write(6,*) '\nRunning the new solver, 3D tests (zyx)'
c        compute the rhs (put in ux)
      ux = utrue
      order = (/ 1, 2, 3 /)
      call ptri_mult(ux, order, ierr)
c        call the solver
      order = (/ 3, 2, 1 /)
      call ptri_sol(ux, order, ierr)
c        compute the errors
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = lerr + (utrue(i,j,k)-ux(i,j,k))**2
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_SUM,comm3d,ierr)
      if(outnode)  write(6,*) '   l2 error =',sqrt(gerr)
      lerr = 0.d0
      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               lerr = max(lerr,abs(utrue(i,j,k)-ux(i,j,k)))
            end do
         end do
      end do
      call mpi_allreduce(lerr,gerr,1,MPI_REAL8,MPI_MAX,comm3d,ierr)
      if(outnode)  write(6,*) '  max error =',gerr



c     clean up
      call ptri_mult_free(ierr)
      call ptri_free(ierr)

      return
      end subroutine TestDriver
c=================================================================




      subroutine mult_strip_setup(idir, j, k, ww, L, D, U, s, ier)
c-----------------------------------------------------------------------
c     Description: mult_strip_setup fills in the periodic tridiagonal 
c     system 
c             L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1)
c     
c     Arguments:
c        idir - (int, input) direction for upcoming multiply
c           j - (int, input) fast index orthogonal to idir
c           k - (int, input) slow index orthogonal to idir
c          ww - (dbl(*), input) array to multiply
c           L - (dbl(*), output) lower diagonal for matrix
c           D - (dbl(*), output) diagonal for matrix
c           U - (dbl(*), output) upper diagonal for matrix
c           s - (dbl(*), output) vector to multiply
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common
      use boundary_conds
      use mpistuff

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in) :: idir, j, k
      integer, intent(out) :: ier
      real*8, dimension(*) :: L, D, U, s
      real*8 :: ww(nxlocal,nylocal,nzlocal)

c     Local variables
      integer :: i, iLeft
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     set up system based on direction
      if (idir == 1) then

c        compute number of cells to the left of this processor
         iLeft = xLloc/dx

c        extract vector to multiply
         do i=1,nxlocal
            s(i) = ww(i,j,k)
         end do

c        set matrix arrays and vector to multiply
         do i=1,nxlocal
c$$$            L(i) = 1.d0 + (i+iLeft)*dx
c$$$            D(i) = 1.d1 + 1.d1*(i+iLeft)*dx
c$$$            U(i) = 2.d0 - (i+iLeft)*dx
            L(i) = 1.d0
            D(i) = 10.d0
            U(i) = 2.d0
         enddo

c        undo couplings for boundary, non-periodic case
         if (xbc /= BCperiodic) then
            if (iprocx == 1)       L(1) = 0.d0
            if (iprocx == xprocs)  U(nxlocal) = 0.d0
         endif

      else if (idir == 2) then

c        compute number of cells to the left of this processor
         iLeft = yLloc/dy

c        extract vector to multiply
         do i=1,nylocal
            s(i) = ww(j,i,k)
         end do

c        set matrix arrays and vector to multiply
         do i=1,nylocal
c$$$            L(i) = 1.d0 + (i+iLeft)*dy
c$$$            D(i) = 1.d1 + 1.d1*(i+iLeft)*dy
c$$$            U(i) = 2.d0 - (i+iLeft)*dy
            L(i) = 1.d0
            D(i) = 10.d0
            U(i) = 2.d0
         enddo

c        undo couplings for boundary, non-periodic case
         if (ybc /= BCperiodic) then
            if (iprocy == 1)       L(1) = 0.d0
            if (iprocy == yprocs)  U(nylocal) = 0.d0
         endif

      else

c        compute number of cells to the left of this processor
         iLeft = zLloc/dz

c        extract vector to multiply
         do i=1,nzlocal
            s(i) = ww(j,k,i)
         end do

c        set matrix arrays and vector to multiply
         do i=1,nzlocal
c$$$            L(i) = 1.d0 + (i+iLeft)*dz
c$$$            D(i) = 1.d1 + 1.d1*(i+iLeft)*dz
c$$$            U(i) = 2.d0 - (i+iLeft)*dz
            L(i) = 1.d0
            D(i) = 10.d0
            U(i) = 2.d0
         enddo

c        undo couplings for boundary, non-periodic case
         if (zbc /= BCperiodic) then
            if (iprocz == 1)       L(1) = 0.d0
            if (iprocz == zprocs)  U(nzlocal) = 0.d0
         endif

      end if

      return

      end subroutine mult_strip_setup
c=======================================================================





      subroutine mult_strip_pack(idir, j, k, ww, res, ier)
c-----------------------------------------------------------------------
c     Description: mult_strip_pack puts the periodic tridiagonal system 
c     product into the multi-dimensional ww array.
c     
c     Arguments:
c        idir - (int, input) direction for multiply result
c           j - (int, input) fast index orthogonal to idir
c           k - (int, input) slow index orthogonal to idir
c          ww - (dbl(*), input) multi-dimensional product array
c         res - (dbl(*), input) product of 1D system
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in) :: idir, j, k
      integer, intent(out) :: ier
      real*8, dimension(*) :: res
      real*8 :: ww(nxlocal,nylocal,nzlocal)

c     Local variables
      integer :: i
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     pack solution based on direction
      if (idir == 1) then

c        put solution back into ww
         do i=1,nxlocal
            ww(i,j,k) = res(i)
         enddo
         
      else if (idir == 2) then

c        put solution back into ww
         do i=1,nylocal
            ww(j,i,k) = res(i)
         enddo

      else 

c        put solution back into ww
         do i=1,nzlocal
            ww(j,k,i) = res(i)
         enddo

      end if

      return

      end subroutine mult_strip_pack
c=======================================================================





      subroutine prepare_rhs(idir, ww, ier)
c-----------------------------------------------------------------------
c     Description: prepare_rhs is a dummy routine here
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c          ww - (dbl(*)) rhs array
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in)  :: idir
      integer, intent(out) :: ier
      double precision     :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

      return
      end subroutine prepare_rhs
c=======================================================================





      subroutine process_sol(idir, ww, ier)
c-----------------------------------------------------------------------
c     Description: process_sol is a dummy routine here
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c          ww - (dbl(*)) rhs array
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in)  :: idir
      integer, intent(out) :: ier
      double precision     :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

      return
      end subroutine process_sol
c=======================================================================





      subroutine setup_strip(idir, j, k, iw, ww, L, D, U, r, ier)
c-----------------------------------------------------------------------
c     Description: setup_strip fills in the periodic tridiagonal system 
c             r(i) = L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1)
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c           j - (int, input) fast index orthogonal to idir
c           k - (int, input) slow index orthogonal to idir
c          iw - (int, input) wave number to work on
c          ww - (dbl(*), input) rhs array
c           L - (dbl(*), output) lower diagonal for matrix
c           D - (dbl(*), output) diagonal for matrix
c           U - (dbl(*), output) upper diagonal for matrix
c           r - (dbl(*), output) rhs for system
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common
      use mpistuff
      use boundary_conds

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in) :: idir, j, k, iw
      integer, intent(out) :: ier
      real*8, dimension(*) :: D, U, L, r
      real*8 :: ww(nxlocal,nylocal,nzlocal)

c     Local variables
      integer :: i, iLeft
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     set up system based on direction
      if (idir == 1) then

c        compute number of cells to the left of this processor
         iLeft = xLloc/dx

c        set matrix arrays and vector to multiply
         do i=1,nxlocal
c$$$            L(i) = 1.d0 + (i+iLeft)*dx
c$$$            D(i) = 1.d1 + 1.d1*(i+iLeft)*dx
c$$$            U(i) = 2.d0 - (i+iLeft)*dx
            L(i) = 1.d0
            D(i) = 10.d0
            U(i) = 2.d0
         enddo

c        set r to be rhs on the spatial strip
         do i=1,nxlocal
            r(i) = ww(i,j,k)
         enddo

c        undo couplings for boundary, non-periodic case
         if (xbc /= BCperiodic) then
            if (iprocx == 1)       L(1) = 0.d0
            if (iprocx == xprocs)  U(nxlocal) = 0.d0
         endif

      else if (idir == 2) then

c        compute number of cells to the left of this processor
         iLeft = yLloc/dy

c        set matrix arrays and vector to multiply
         do i=1,nylocal
c$$$            L(i) = 1.d0 + (i+iLeft)*dy
c$$$            D(i) = 1.d1 + 1.d1*(i+iLeft)*dy
c$$$            U(i) = 2.d0 - (i+iLeft)*dy
            L(i) = 1.d0
            D(i) = 10.d0
            U(i) = 2.d0
         enddo

c        set r to be rhs on the spatial strip
         do i=1,nylocal
            r(i) = ww(j,i,k)
         enddo

c        undo couplings for boundary, non-periodic case
         if (ybc /= BCperiodic) then
            if (iprocy == 1)       L(1) = 0.d0
            if (iprocy == yprocs)  U(nylocal) = 0.d0
         endif

      else

c        compute number of cells to the left of this processor
         iLeft = zLloc/dz

c        set matrix arrays and vector to multiply
         do i=1,nzlocal
c$$$            L(i) = 1.d0 + (i+iLeft)*dz
c$$$            D(i) = 1.d1 + 1.d1*(i+iLeft)*dz
c$$$            U(i) = 2.d0 - (i+iLeft)*dz
            L(i) = 1.d0
            D(i) = 10.d0
            U(i) = 2.d0
         enddo

c        set r to be rhs on the spatial strip
         do i=1,nzlocal
            r(i) = ww(j,k,i)
         enddo

c        undo couplings for boundary, non-periodic case
         if (zbc /= BCperiodic) then
            if (iprocz == 1)       L(1) = 0.d0
            if (iprocz == zprocs)  U(nzlocal) = 0.d0
         endif

      end if

      return

      end subroutine setup_strip
c=======================================================================





      subroutine pack_strip(idir, j, k, iw, ww, sol, ier)
c-----------------------------------------------------------------------
c     Description: pack_strip puts the periodic tridiagonal system 
c     solution into the multi-dimensional ww array.
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c           j - (int, input) fast index orthogonal to idir
c           k - (int, input) slow index orthogonal to idir
c          iw - (int, input) wave number to work on
c          ww - (dbl(*), input) multi-dimensional sol array
c         sol - (dbl(*), input) sol for 1D system
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in) :: idir, j, k, iw
      integer, intent(out) :: ier
      real*8, dimension(*) :: sol
      real*8 :: ww(nxlocal,nylocal,nzlocal)

c     Local variables
      integer :: i
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     pack solution based on direction
      if (idir == 1) then

c        put solution back into ww
         do i=1,nxlocal
            ww(i,j,k) = sol(i)
         enddo
         
      else if (idir == 2) then

c        put solution back into ww
         do i=1,nylocal
            ww(j,i,k) = sol(i)
         enddo

      else 

c        put solution back into ww
         do i=1,nzlocal
            ww(j,k,i) = sol(i)
         enddo

      end if

      return

      end subroutine pack_strip
c=======================================================================
