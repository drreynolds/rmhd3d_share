c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------
c     $Log: fhypre_test.F,v $
c=================================================================




      program fhypre_test
c-----------------------------------------------------------------
#ifdef PARALLEL
      use mpi
#endif
      
c======= Declarations =========
      implicit none

c     parameters for problem domain, etc
      integer, parameter :: ndims=3
      integer, parameter :: xprocs=2
      integer, parameter :: yprocs=1
      integer, parameter :: zprocs=1
      integer, parameter :: NxL=16
      integer, parameter :: NyL=32
      integer, parameter :: NzL=32
      integer, parameter :: NGx=4
      integer, parameter :: NGy=4
      integer, parameter :: NGz=4
      integer, parameter :: Nvar=8
      integer, parameter :: ixlo=1-NGx
      integer, parameter :: ixhi=NxL+NGx
      integer, parameter :: iylo=1-NGy
      integer, parameter :: iyhi=NyL+NGy
      integer, parameter :: izlo=1-NGz
      integer, parameter :: izhi=NzL+NGz
      double precision   :: xL =  0.0d0
      double precision   :: xR =  2.0d0
      double precision   :: yL = -1.0d0
      double precision   :: yR =  1.0d0
      double precision   :: zL = -1.0d0
      double precision   :: zR =  1.0d0
      double precision   :: dx, dy, dz
      integer :: iprocx, iprocy, iprocz, npes
      integer :: left, right, top, bottom, forward, behind
      integer :: iXL, iXR, iYL, iYR, iZL, iZR
      integer :: Nx, Ny, Nz, xbc, ybc, zbc
      integer :: comm3d, NProcAxis(3), my_id, ierr
      logical :: periods(3), reorder

c     local driver variables
      double precision :: uu(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: mu, eta, kappa, gamma, rgas, Re, Lu, Pr
      double precision :: dtfac, delta, alternating
      integer          :: i, j, k, l
      integer          :: TimeSteps, tstep, NewtonIts, newt
      integer          :: ix, iy, iz
      integer          :: MPcoords(3)

c======= Internals ============


      
c     initialize MPI and set up cartesian communicator
      NprocAxis(1) = xprocs
      NprocAxis(2) = yprocs
      NprocAxis(3) = zprocs
      xbc = 0; ybc = 0; zbc = 0
      periods = .false.
#ifdef XPERIODIC
      xbc = 1
      periods(1) = .true.
#elif XREFLECT
      xbc = 2
#endif
#ifdef YPERIODIC
      ybc = 1
      periods(2) = .true.
#elif YREFLECT
      ybc = 2
#endif
#ifdef ZPERIODIC
      zbc = 1
      periods(3) = .true.
#elif ZREFLECT
      zbc = 2
#endif
      call MPI_init(ierr)
      call MPI_comm_size(MPI_COMM_WORLD, npes, ierr)
      if (npes /= xprocs*yprocs*zprocs) then
         write(*,'(2(A,i3))') 'ERROR: ',npes,
     &        ' procs provided, run requires ',xprocs*yprocs*zprocs
         pause
      endif
      call MPI_comm_rank(MPI_COMM_WORLD, my_id, ierr)
      write(*,'(2(A,i3))') 'fhypre_test: starting proc ',
     &     my_id,' of ',npes
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call MPI_cart_create(MPI_COMM_WORLD, ndims, NprocAxis, periods, 
     &     reorder, comm3D, ierr)
      call MPI_cart_get(comm3D, 3, NprocAxis, periods, MPcoords, ierr)
      iprocx = MPcoords(1)+1
      iprocy = MPcoords(2)+1
      iprocz = MPcoords(3)+1
      write(*,'(A,i3,A,3(i2,1x),A,3(i2,1x),A,3(i2,1x),A)') 
     &     '  p',my_id,': iprocs = (',
     &     NprocAxis(1),NprocAxis(2),NprocAxis(3),'), periods = (',
     &     periods(1),periods(2),periods(3),'), coords = (',
     &     MPcoords(1),MPcoords(2),MPcoords(3),')'
      call MPI_Barrier(comm3d, ierr)

      Nx = NxL*NprocAxis(1)
      Ny = NyL*NprocAxis(2)
      Nz = NzL*NprocAxis(3)
      iXL = NxL*(iprocx-1)+1
      iXR = iXL + NxL - 1
      iYL = NyL*(iprocy-1)+1
      iYR = iYL + NyL - 1
      iZL = NzL*(iprocz-1)+1
      iZR = iZL + NzL - 1
      write(*,'(A,i3,A,3(i4,1x),A,6(i4,A))') 
     &     '  p',my_id,': globdims = (',Nx,Ny,Nz,'), extents = (',
     &     iXL, ':', iXR, ',', iYL, ':', iYR, ',', iZL, ':', iZR,')'
      call MPI_Cart_shift(comm3D, 0, -1, i, left,    ierr)
      call MPI_Cart_shift(comm3D, 0,  1, i, right,   ierr)
      call MPI_Cart_shift(comm3D, 1, -1, i, bottom,  ierr)
      call MPI_Cart_shift(comm3D, 1,  1, i, top,     ierr)
      call MPI_Cart_shift(comm3D, 2, -1, i, behind,  ierr)
      call MPI_Cart_shift(comm3D, 2,  1, i, forward, ierr)
      write(*,'(A,i3,A,6(i2,A))') 
     &     '  p',my_id,': neighbors = (', left,':', right,',', 
     &     bottom,':', top,',', behind,':', forward,')'
      call MPI_Barrier(comm3D, ierr)

c     initialize domain information
      if (my_id == 0) print *,' driver: initializing domain info.'
      dx = (xR-xL)/Nx
      dy = (yR-yL)/Ny
      dz = (zR-zL)/Nz

c     initialize state variables
      uu(:,:,:,:) = 1.0d0
      vv(:,:,:,:) = 0.0d0
      ww(:,:,:,:) = 1.0d0

c     set MHD parameters
      dtfac = 1.0d-1
      mu    = 1.0d-1
      Lu    = 1.0d0
      eta   = 1.0d-0
      Re    = 1.0d0
      gamma = 1.667d0
      kappa = 1.0d-2
      Pr    = 1.0d0
      rgas  = 1.0d0

c     set HYPRE solver tolerance
      delta = 1.0d-5;


c     custom viscous momentum preconditioner initialization
      if (my_id == 0) print *,' driver: calling FVPrecDuInit'
      call FVPrecDuInit(NxL, NyL, NzL, dx, dy, dz, 
     &     nvar, NGx, NGy, NGz, NprocAxis(1), NprocAxis(2), 
     &     NprocAxis(3), iprocx, iprocy, iprocz, left, right, top, 
     &     bottom, forward, behind, xbc, ybc, zbc, ierr)
      call MPI_Barrier(comm3D, ierr)

c     custom viscous magnetic field preconditioner initialization
      if (my_id == 0) print *,' driver: calling FVPrecDbInit'
      call FVPrecDbInit(NxL, NyL, NzL, dx, dy, dz, 
     &     nvar, NGx, NGy, NGz, NprocAxis(1), NprocAxis(2), 
     &     NprocAxis(3), iprocx, iprocy, iprocz, left, right, top, 
     &     bottom, forward, behind, xbc, ybc, zbc, ierr)
      call MPI_Barrier(comm3D, ierr)


c     custom viscous energy preconditioner initialization
      if (my_id == 0) print *,' driver: calling FVPrecDeInit'
      call FVPrecDeInit(NxL, NyL, NzL, dx, dy, dz, 
     &     nvar, NGx, NGy, NGz, NprocAxis(1), NprocAxis(2), 
     &     NprocAxis(3), iprocx, iprocy, iprocz, left, right, top, 
     &     bottom, forward, behind, xbc, ybc, zbc, ierr)
      call MPI_Barrier(comm3D, ierr)


c     loop over 'time' steps
      TimeSteps = 5
      NewtonIts = 5
      alternating = 1.0
      do tstep=1,TimeSteps,1

         call MPI_Barrier(comm3D, ierr)
         if (my_id == 0) print '(A,2(i2,A))',
     &        '\n Beginning time step',tstep,'/',TimeSteps,':\n'

c        update state vector uu
         alternating = -alternating
         do l=1,nvar,1
            do k=izlo,izhi,1
               do j=iylo,iyhi,1
                  do i=ixlo,ixhi,1
                     iz = iZL + k - 1
                     iy = iYL + j - 1
                     ix = iXL + i - 1
                     uu(i,j,k,l) = 2.d0 
     &                    + alternating*tstep/TimeSteps*ix/Nx*iy/Ny
                  enddo
               enddo
            enddo
         enddo

c        call viscous preconditioner matrix setup routines
         if (my_id == 0) print *,'    calling fvprecdusetup'
         call fvprecdusetup(uu, dtfac, mu, Re, vv, ww, ierr)
         if (ierr /= 0) then
            write(0,*) '    fvprecdusetup error =', ierr
         endif
         call MPI_Barrier(comm3d, ierr)
         
         if (my_id == 0) print *,'    calling fvprecdbsetup'
         call fvprecdbsetup(uu, dtfac, Lu, eta, vv, ww, ierr)
         if (ierr /= 0) then
            write(0,*) '    fvprecdbsetup error =', ierr
         endif
         call MPI_Barrier(comm3d, ierr)
         
         if (my_id == 0) print *,'    calling fvprecdesetup'
         call fvprecdesetup(uu, dtfac, gamma, kappa, Re, Pr, rgas, 
     &        vv, ww, ierr)
         if (ierr /= 0) then
            write(0,*) '    fvprecdesetup error =', ierr
         endif
         call MPI_Barrier(comm3d, ierr)
         
         
c        take 'Newton' steps (solve multiple systems with same matrix)
         do newt=1,NewtonIts,1

            if (my_id == 0) print '(A,2(i2,A))',
     &           '\n  Newton Iteration',newt,'/',NewtonIts,':\n'
            alternating = -alternating

c           solve systems, get iteration counts
            if (my_id == 0) print *,
     &           '\n  fhypre_test: solving Du matrix system:\n'
            do l=1,nvar,1
               do k=1,NzL,1
                  do j=1,NyL,1
                     do i=1,NxL,1
                        ww(i,j,k,l)=6.4d0*(i-1)*(NxL-i)/NxL/NxL
     &                       *(j-1)*(NyL-j)/NyL/NyL
     &                       *(k-1)*(NzL-j)/NzL/NzL
     &                       *alternating*tstep*newt
                     enddo
                  enddo
               enddo
            enddo
            call FVPrecDuSolve(vv, ww, uu, delta, ierr)
            if (ierr /= 0) then
               if (my_id == 0) print *,
     &              '  fhypre_test error: Du solve returned flag ',
     &              ierr
               call FVPrecDuFree
               call FVPrecDbFree
               call FVPrecDeFree
               call MPI_Finalize(ierr)
            endif
            call MPI_Barrier(comm3d, ierr)

            if (my_id == 0) print *,
     &           '\n  fhypre_test: solving Db matrix system:\n'
            do l=1,nvar,1
               do k=1,NzL,1
                  do j=1,NyL,1
                     do i=1,NxL,1
                        ww(i,j,k,l)=6.4d0*(i-1)*(NxL-i)/NxL/NxL
     &                       *(j-1)*(NyL-j)/NyL/NyL
     &                       *(k-1)*(NzL-j)/NzL/NzL
     &                       *alternating*tstep*newt
                     enddo
                  enddo
               enddo
            enddo
            call FVPrecDbSolve(vv, ww, uu, delta, ierr)
            if (ierr /= 0) then
               if (my_id == 0) print *,
     &              '  fhypre_test error: Db solve returned flag ',
     &              ierr
               call FVPrecDuFree
               call FVPrecDbFree
               call FVPrecDeFree
               call MPI_Finalize(ierr)
            endif
            call MPI_Barrier(comm3d, ierr)
            
            if (my_id == 0) print *,
     &           '\n  fhypre_test: solving De matrix system:\n'
            do l=1,nvar,1
               do k=1,NzL,1
                  do j=1,NyL,1
                     do i=1,NxL,1
                        ww(i,j,k,l)=6.4d0*(i-1)*(NxL-i)/NxL/NxL
     &                       *(j-1)*(NyL-j)/NyL/NyL
     &                       *(k-1)*(NzL-j)/NzL/NzL
     &                       *alternating*tstep*newt
                     enddo
                  enddo
               enddo
            enddo
            call FVPrecDeSolve(vv, ww, uu, delta, ierr)
            if (ierr /= 0) then
               if (my_id == 0) print *,
     &              '  fhypre_test error: De matrix returned flag ',
     &              ierr
               call FVPrecDuFree
               call FVPrecDbFree
               call FVPrecDeFree
               call MPI_Finalize(ierr)
            endif
            call MPI_Barrier(comm3d, ierr)
            
         enddo
      enddo         

c     free up solver dataspace memory
      if (my_id == 0) print *,
     &     '\n  fhypre_test: freeing Du matrix/solver memory\n'
      call FVPrecDuFree
      call MPI_Barrier(comm3d, ierr)
      if (my_id == 0) print *,
     &     '\n  fhypre_test: freeing Db matrix/solver memory\n'
      call FVPrecDbFree
      call MPI_Barrier(comm3d, ierr)
      if (my_id == 0) print *,
     &     '\n  fhypre_test: freeing De matrix/solver memory\n'
      call FVPrecDeFree
      call MPI_Barrier(comm3d, ierr)

#ifdef PARALLEL
c     finalize MPI
      call MPI_Finalize(ierr)

#endif
      end program fhypre_test
c=================================================================
