      program fnvec_mhd_test
c=======================================================================c
c     File        : fnvec_mhd_test.F                                    c
c     Programmers : Daniel R. Reynolds @ SMU                            c
c     Version of  : 2 October 2003                                      c
c-----------------------------------------------------------------------c
c     This testing routine is used to ensure the validity of the        c
c     Fortran-C interface to the new N_Vector kernel.                   c
c                                                                       c
c     Execute line:  mpirun -np N fnvec_mhd_test.x                      c
c         where N is the number of processors to use                    c
c=======================================================================c
c------ Inclusions ------
      use mesh
      use mesh_common
      use mpistuff

c------ Declarations ------
      implicit none

c     mesh variables
c$$$      integer, parameter :: nx = 3
c$$$      integer, parameter :: ny = 2
c$$$      integer, parameter :: nz = 1
c$$$      integer, parameter :: ns = 4
      integer :: xghost, yghost, zghost

c     internal variables
      integer :: nglobal, ier, outputproc
      integer :: i, j, k, l
      integer :: NProcAxis(3)
      logical :: periods(3), reorder
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: u, xvec, yvec, zvec
      double precision :: urp, xrp, yrp, zrp
      double precision :: a, b, reduction


c     determine directional ghost cells
      xghost = 1-ixlo
      yghost = 1-iylo
      zghost = 1-izlo

c     determine the processor that will output local information
      outputproc = 0

c     this test must be compiled with to be used in parallel
#ifndef PARALLEL
      print *, 'Test_fnvec_mhd Error: this test must be run with MPI'
      stop
#endif

c     output progress
      print *, 'Test_fnvec_mhd: Local variables initialized:'
      print *, 'Test_fnvec_mhd: Initializing MPI...'

c     The user MUST call mpi_init, Fortran binding, for the interface 
c     package to work. The communicator, MPI_COMM_WORLD, is the only 
c     one common between the Fortran and C bindings. So in the following, 
c     the communicator MPI_COMM_WORLD is used in calls to mpi_comm_size 
c     and mpi_comm_rank to determine the total number of processors and 
c     the rank (0 ... size-1) number of this process.
      call MPI_Init(ier)
      NProcAxis = (/ Xprocs, Yprocs, Zprocs /)
      periods = (/ .false., .false., .false. /)
#ifdef XPERIODIC
      periods(1) = .true.
#endif
#ifdef YPERIODIC
      periods(2) = .true.
#endif
#ifdef ZPERIODIC
      periods(3) = .true.
#endif
      reorder = .false.
      call MPI_cart_create( MPI_COMM_WORLD, 3, NProcAxis, periods,
     1     reorder, comm3D, ierr)

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      print *, '   '
      
c     output progress
      print *, 'Test_fnvec_mhd: Obtaining number of processors...'

c     get number of processors, etc
      call mpi_comm_size(MPI_COMM_WORLD,nprocs,ier)
      if (ier .ne. 0) then
         print '(A,i2)','Test_fnvec_mhd: mpi_comm_size failed, ier = '
     &        ,ier
         stop
      endif

c     output progress
      print '(A,i1)', '     mpi comm size = ',nprocs
      print *, 'Test_fnvec_mhd: Obtaining processor rank...'

c     get local processor rank
      call mpi_comm_rank(MPI_COMM_WORLD, my_id, ier)
      if (ier .ne. 0) then
         print '(A,i2)','Test_fnvec_mhd: mpi_comm_rank failed, ier = '
     &        ,ier
         stop
      endif

c     output progress
      print '(A,i1)', '     mpi processor rank = ',my_id

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)

c     initialize mesh parameters
      nxlocal = nxlsize
      nylocal = nylsize
      nzlocal = nzlsize

c     initialize the fortran-C nvec_mhd_parallel specification
      nglobal = (nxlocal+2*xghost)*(nylocal+2*yghost)*
     &     (nzlocal+2*zghost)*nvar
      print '(A,i2,A)', 'Test_fnvec_mhd: proc',my_id
      call fnvinitmhd(5, nglobal, ier)
      if (ier .ne. 0) then
         print '(A,i2)','Test_fnvec_mhd: proc',my_id,
     &        ', fnvinitmhd failed, ier = ',ier
         stop
      endif

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      
c     output progress
      print '(A,i2,A)', 'Test_fnvec_mhd: proc',my_id,
     &     ', filling in local array...'
      
c     fill in local array, not including ghost cells
      u = 0d0
      urp = 0.d0
      do i = 1,nxlocal
         do j = 1,nylocal
            do k = 1,nzlocal
               do l = 1,nvar
                  u(i,j,k,l) = 1d4*my_id + 1d3*i + 1d2*j + 1d1*k + 1d0*l
               enddo
            enddo
         enddo
      enddo
      urp = 5.4321d0

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      
c     root node prints local array to screen for sanity check
      if (my_id == outputproc) then
         if (ier == 0) then
            print *, '  '
            call print_fnvecmhd(u,urp)
            print *, '  '
         end if 
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      
c     output progress
      if (my_id == outputproc) then
         print '(A,i2,A)', 'Test_fnvec_mhd: proc',my_id,
     &        ', running C test routine fnvecmhdtest...'
      end if

c     run C testing routine
      call fnvecmhdtest( u, urp, my_id, outputproc, nxlocal, 
     &     nylocal, nzlocal, xghost, yghost, zghost, nvar, ier)
      if (ier .ne. 0) then
         print '(A,i2)','Test_fnvec_mhd: proc',my_id,
     &        ', fnvecmhdtest failed, ier = ',ier
         stop
      endif

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      
c     output progress
      if (my_id == outputproc) then
         print '(A,i2,A)', 'Test_fnvec_mhd: proc',my_id,
     &        ', running local tests...'
      end if

c     run local tests to check that nvec_mhd library is
c     calculating things correctly
      
c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 0:
      if (my_id == outputproc) then
         print *, '  F-Test 0: printing initial vector to file F_test00'
         open(unit=100,action='write',file='F_test00')
         call print_fnvecmhd_file(u,urp,100)
         close(unit=100)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)

      
c     TEST 1: skipped
      if (my_id == outputproc) then
         print *, '  F-Test 1: skipped'
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)      


c     TEST 2:
      a = 1d0; b = 3d0
      xvec = a*u
      xrp = a*urp
      yvec = b*xvec
      yrp = b*xrp
      if (my_id == outputproc) then
         print *, '  F-Test 2: printing results to file F_test02'
         open(unit=101,action='write',file='F_test02')
         call print_fnvecmhd_file(yvec,yrp,101)
         close(unit=101)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 3:
      a = 5d-1; b = -5d-1
      zvec = a*xvec + b*yvec
      zrp = a*xrp + b*yrp
      if (my_id == outputproc) then
         print *, '  F-Test 3: printing results to file F_test03'
         open(unit=103,action='write',file='F_test03')
         call print_fnvecmhd_file(zvec,zrp,103)
         close(unit=103)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 4:
      a = -1d0; b = 1d0
      zvec = a*xvec + b*yvec
      zrp = a*xrp + b*yrp
      if (my_id == outputproc) then
         print *, '  F-Test 4: printing results to file F_test04'
         open(unit=104,action='write',file='F_test04')
         call print_fnvecmhd_file(zvec,zrp,104)
         close(unit=104)
      end if 

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 5:
      a = 1d0; b = 5d-1
      zvec = a*xvec + b*yvec
      zrp = a*xrp + b*yrp
      if (my_id == outputproc) then
         print *, '  F-Test 5: printing results to file F_test05'
         open(unit=105,action='write',file='F_test05')
         call print_fnvecmhd_file(zvec,zrp,105)
         close(unit=105)
      end if 

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 6:
      a = 5d0; b = 1d0
      yvec = a*xvec + b*yvec
      yrp = a*xrp + b*yrp
      if (my_id == outputproc) then
         print *, '  F-Test 6: printing results to file F_test06'
         open(unit=106,action='write',file='F_test06')
         call print_fnvecmhd_file(yvec,yrp,106)
         close(unit=106)
      end if 

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 7:
      a = 8.2d0
      zvec = a
      zrp = a
      if (my_id == outputproc) then
         print *, '  F-Test 7: printing results to file F_test07'
         open(unit=107,action='write',file='F_test07')
         call print_fnvecmhd_file(zvec,zrp,107)
         close(unit=107)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 8:
      b = 3d0
      yvec = b*u
      yrp = b*urp
      zvec = xvec*yvec
      zrp = xrp*yrp
      if (my_id == outputproc) then
         print *, '  F-Test 8: printing results to file F_test08'
         open(unit=108,action='write',file='F_test08')
         call print_fnvecmhd_file(zvec,zrp,108)
         close(unit=108)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 9:
      where ( yvec /= 0.d0 )
         zvec = xvec/yvec
      elsewhere
         zvec = 0.0d0
      end where
      if (yrp /= 0.d0 ) then
         zrp = xrp/yrp
      else
         zrp = 0.d0
      endif
      if (my_id == outputproc) then
         print *, '  F-Test 9: printing results to file F_test09'
         open(unit=109,action='write',file='F_test09')
         call print_fnvecmhd_file(zvec,zrp,109)
         close(unit=109)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 10:
      b = -2d4
      zvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar) = 
     &     xvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar) + b
      zrp = xrp + b
      if (my_id == outputproc) then
         print *, '  F-Test 10: printing results to file F_test10'
         open(unit=110,action='write',file='F_test10')
         call print_fnvecmhd_file(zvec,zrp,110)
         close(unit=110)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 11:
      yvec = abs(zvec)
      yrp = abs(zrp)
      if (my_id == outputproc) then
         print *, '  F-Test 11: printing results to file F_test11'
         open(unit=111,action='write',file='F_test11')
         call print_fnvecmhd_file(yvec,yrp,111)
         close(unit=111)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 12:
      where ( yvec /= 0.0d0 )
         zvec = 1.d0/yvec
      elsewhere
         zvec = 0.0d0
      end where   
      if (yrp /= 0.d0 ) then
         zrp = 1.d0/yrp
      else
         zrp = 0.d0
      endif
      if (my_id == outputproc) then
         print *, '  F-Test 12: printing results to file F_test12'
         open(unit=112,action='write',file='F_test12')
         call print_fnvecmhd_file(zvec,zrp,112)
         close(unit=112)
      end if 

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 13:
      b = 3d0
      yvec = b*u
      yrp = b*urp
      a = sum(xvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar)
     &     *yvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar))
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, MPI_COMM_WORLD, ier)
      a = a + xrp*yrp
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 13: dot_product(xvec,yvec) = ',
     &        reduction
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 14:
      a = maxval(xvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar))
      a = max(a, xrp)
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MAX, MPI_COMM_WORLD, ier)
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 14: max(xvec) = ',reduction
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 15:
      b = 1.5d0
      yvec = b
      yrp = b
      zvec = xvec*yvec
      zrp = xrp*yrp
      a = sum(zvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar)
     &     *zvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar))
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, MPI_COMM_WORLD, ier)
      a = a + zrp*zrp
      reduction = sqrt(reduction/(nx*ny*nz*nvar+1))
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 15: WrmsNorm(xvec,yvec) = ',
     &        reduction
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 16:
      if (my_id == outputproc) then
         print *, '  F-Test 16: skipped'
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 17:
      a = minval(xvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar))
      a = min(a,xrp)
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MIN, MPI_COMM_WORLD, ier)
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 17: min(xvec) = ',reduction
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 18:
      b = 2d0
      yvec = b*u
      yrp = b*urp
      zvec = xvec * yvec
      zrp = xrp * yrp
      a = sum(zvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar)
     &     *zvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar))
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, MPI_COMM_WORLD, ier)
      a = a + zrp*zrp
      reduction = sqrt(reduction)
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 18: WL2Norm(xvec,yvec) = ',
     &        reduction
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 19:
      a = sum( abs(xvec(1:nxlocal,1:nylocal,1:nzlocal,1:nvar)) )
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, MPI_COMM_WORLD, ier)
      a = a + abs(xrp)
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 19: L1Norm(xvec) = ',reduction
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 20:
      a = 2.d3
      zvec = 0.d0
      zrp = 0.d0
      do i = 1,nxlocal
         do j = 1,nylocal
            do k = 1,nzlocal
               do l = 1,nvar
                  if ( abs(xvec(i,j,k,l)) > a ) zvec(i,j,k,l) = 1d0
               end do
            end do
         end do
      end do
      if (abs(xrp) > a)  zrp = 1.d0
      if (my_id == outputproc) then
         print *, '  F-Test 20: printing results to file F_test20'
         open(unit=120,action='write',file='F_test20')
         call print_fnvecmhd_file(zvec,zrp,120)
         close(unit=120)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 21:
      a = 1.d0
      zvec = 0.d0
      zrp = 0.d0
      do i = 1,nxlocal
         do j = 1,nylocal
            do k = 1,nzlocal
               do l = 1,nvar
                  if ( xvec(i,j,k,l) == 0.d0 ) then
                     a = 0.d0
                  else  
                     zvec(i,j,k,l) = 1.d0/xvec(i,j,k,l)
                  end if
               end do
            end do
         end do
      end do      
      if (xrp == 0.d0) then
         a = 0.d0
      else
         zrp = 1.d0/xrp
      endif
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MIN, MPI_COMM_WORLD, ier)
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 21: InvTest(xvec) = ',reduction
         print *, '    printing results to file F_test21'
         open(unit=121,action='write',file='F_test21')
         call print_fnvecmhd_file(zvec,zrp,121)
         close(unit=121)
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     TEST 22:
      b = 4.d-5
      yvec = b*u
      yrp = b*urp
      b = 1.d0
      a = 0.d0
      do i = 1,nxlocal
         do j = 1,nylocal
            do k = 1,nzlocal
               do l = 1,nvar
                  if ( yvec(i,j,k,l) /= 0.d0 ) then
                     if ( b == 1.d0 ) then
                        a = xvec(i,j,k,l) / yvec(i,j,k,l)
                        b = 0.d0
                     else
                        a = min( a, xvec(i,j,k,l)/yvec(i,j,k,l) )
                     end if
                  end if
               end do
            end do
         end do
      end do      
      if (yrp /= 0.d0) then
         if (b == 1.d0) then
            a = xrp / yrp
            b = 0.d0
         else
            a = min(a, xrp/yrp)
         endif
      endif
      call MPI_Allreduce( a, reduction, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MIN, MPI_COMM_WORLD, ier)
      if (my_id == outputproc) then
         print '(A,ES12.5)', '   F-Test 22: MinQuotient(xvec,yvec) = ',
     &        reduction
         print *,'  '
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      

c     output progress
      if (my_id == outputproc) then
         print '(A,i2,A)', 'Test_fnvec_mhd: proc',my_id,
     &        ', finalizing MPI...\n'
      end if

c     pause test execution for everyone to catch up
      call mpi_barrier(MPI_COMM_WORLD,ier)
      
c     An explicit call to mpi_finalize (Fortran binding) is required by 
c     the constructs used in the interface
      call mpi_finalize(ier)
      
c     output progress
      if (my_id == outputproc) then
         print '(A,i2,A)', 'Test_fnvec_mhd: Tests Complete.\n'
      end if

      end program fnvec_mhd_test
c=======================================================================c




c     The fortran vector output function (to screen)
      subroutine print_fnvecmhd(uu,urp)

c     modules
      use mesh
      use mesh_common
      use mpistuff

c     calling variable
      double precision :: uu(*), urp

c     local variables
      integer :: i, j, k, l, jbase, kbase, sbase, loc

c     internals
      print '(A6,i1,A3,ES16.5)',
     &     '   urp',my_id,' = ',urp
      do i = 1,nxlocal
         do j = 1,nylocal
            jbase = (j-iylo)*(ixhi-ixlo+1)
            do k = 1,nzlocal
               kbase = (k-izlo)*(ixhi-ixlo+1)*(iyhi-iylo+1)
               do l = 1,nvar
                  sbase = (l-1)*(ixhi-ixlo+1)
     &                 *(iyhi-iylo+1)*(izhi-izlo+1)
                  loc = sbase + kbase + jbase + (i-ixlo+1)
                  print '(A6,i1,4(A2,i6),A4,ES16.5)', 
     &                 '     u',my_id,'(',i,', ',j,', ',k,', ',l,') = ',
     &                 uu(loc)
               enddo
            enddo
         enddo
      enddo
      print *, '  '

      end subroutine print_fnvecmhd
c=======================================================================c




c     The fortran vector output function (to file)
      subroutine print_fnvecmhd_file(uu,urp,fhandle)

c     modules
      use mesh
      use mesh_common
      use mpistuff

c     calling variables
      double precision :: uu(*), urp
      integer          :: fhandle

c     local variables
      integer :: i, j, k, l, jbase, kbase, sbase, loc

c     internals
      write(fhandle,'(A3,5(i6,1x),ES16.5)') 
     &     '   ',my_id, 0, 0, 0, 0, urp
      do i = 1,nxlocal
         do j = 1,nylocal
            jbase = (j-iylo)*(ixhi-ixlo+1)
            do k = 1,nzlocal
               kbase = (k-izlo)*(ixhi-ixlo+1)*(iyhi-iylo+1)
               do l = 1,nvar
                  sbase = (l-1)*(ixhi-ixlo+1)
     &                 *(iyhi-iylo+1)*(izhi-izlo+1)
                  loc = sbase + kbase + jbase + (i-ixlo+1)
                  write(fhandle,'(A3,5(i6,1x),ES16.5)') 
     &                 '   ',my_id, i, j, k, l, uu(loc)
               enddo
            enddo
         enddo
      enddo

      end subroutine print_fnvecmhd_file
c=======================================================================c
