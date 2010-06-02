c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     SMU Mathematics
c     5/1/2009
c
c     This file implements a parallelized solver for ADI-like solves 
c            T1 * T2 * T3 * sol = rhs
c     where each of the T matrices are periodic tridiagonal linear 
c     systems of the form
c             r(i) = a(i)*s(i-1) + b(i)*s(i) + c(i)*s(i+1).
c     The implicit assumption is that each processor will own a large
c     number of these systems, within either a 1D, 2D or 3D processor 
c     topology.  Under this assumption, parallelization is obtained 
c     through a pipelining process, whereby each periodic tridiagonal 
c     solve is broken into three components: 
c            (i) a forward sweep, 
c           (ii) a reverse sweep, and 
c          (iii) a forward update.  
c     The pipeline is achieved through performing (i) for all systems
c     in a given coordinate direction, followed by (ii), and followed 
c     by (iii).  
c
c     Once these solves are completed in one direction, the solver 
c     proceeds to the next coordinate direction and repeats the process.
c
c     We have a number of options:
c
c     - Dimension - This solver may be used for 1, 2 or 3 dimensional 
c       systems.
c
c     - Multiple directional systems - If a user wishes to solve multiple
c       systems in a given coordinate direction before moving along to
c       the next coordinate direction, they may specify that number and 
c       the code will attack these one a time, reusing the data structures
c       for the multiple solves.
c
c     - Customized solver order - A user may specify the order of the 
c       three (or less) solves, so that they proceed in a desired order, 
c       e.g. (x,z,y), (z,y,x), (x), (y), or even (x,x).
c
c     - User-provided pre-processing, post-processing and system setup 
c       routines - The user is expected to provide the following routines 
c       that are called by this solver:
c
c       prepare_rhs(idir, vec, ier) -- Given the direction idir {0,1,2} 
c          and the handle to the user-supplied right-hand-side array vec, 
c          any necessary pre-processing of the user's data to prepare for 
c          solving idir-directional systems may be performed.
c
c       setup_strip(idir, j, k, iw, vec, b, c, a, r, ier) -- Given the 
c          description of which 1D tridiagonal system to set up (j,k,iw 
c          correspond to the 2 orthogonal spatial indices and the given 
c          'wave' in that direction), and the user-supplied right-hand-side 
c          array vec, this routine should return the arrays b, c, a and r 
c          for the tridiagonal system 
c             r(i) = a(i)*s(i-1) + b(i)*s(i) + c(i)*s(i+1)
c          on this processor.
c
c       pack_strip(idir, j, k, iw, vec, s, ier) -- Given the description 
c          of which 1D tridiagonal system has just been solved, the 
c          computed solution s, and the user-supplied solution array vec 
c          (overwrites the user-supplied rhs array), this routine should 
c          pack the 1D solution values from s into vec.
c
c       process_sol(idir, vec, ier) -- Given teh direction idir and the 
c          user-supplied solution array vec, any necessary post-processing 
c          of the solution data to convert back from solving the
c          directional systems may be performed here.
c
c       All of these routines should also return the success/failure 
c       flag ier, where:
c              0 -> success
c              1 -> recoverable failure
c             -1 -> fatal failure
c
c-----------------------------------------------------------------------
c     $aog: ptri_parallel.F,v $
c=======================================================================



      module ptri_solver
c-----------------------------------------------------------------------
c     Description: holds necessary temporary data structures and 1D 
c     communicators.
c-----------------------------------------------------------------------
      use mpi
      implicit none
      save

c     problem dimension
      integer, private :: dims
c     grid dimensions
      integer, private :: nx, ny, nz
c     # of solves in each direction (for multiple wave problems)
      integer, private :: nwave_x, nwave_y, nwave_z
c     1D communicators
      integer, target, private :: x_comm, y_comm, z_comm
c     parallelism information
      integer, private :: pdims(3), pcoords(3), myid_x, myid_y, myid_z
      integer, private :: pxl, pxr, pyl, pyr, pzl, pzr
      logical, private :: periods(3)
c     setup flag (pre-initialized)
      logical, private :: comms_set = .false.

c     temporary arrays
      real*8, dimension(:), allocatable, private :: b, c, a, r, s
      real*8, dimension(:), allocatable, private :: tmp1, tmp2, tmp3
c-----------------------------------------------------------------------



      contains

      
      
      subroutine ptri_setup(comm3d, dims_, nx_, ny_, nz_, 
     &                      nw_x, nw_y, nw_z, ier)
c-----------------------------------------------------------------------
c     Description: sets up the solvers
c
c     Arguments:
c        comm3d - (int, input) 3D cartesian MPI communicator
c         dims_ - (int, input) # of problem dimensions
c           nx_ - (int, input) local x-dimension length
c           ny_ - (int, input) local y-dimension length (ignored if dims<2)
c           nz_ - (int, input) local z-dimension length (ignored if dims<3)
c          nw_x - (int, input) number of systems to solve in x-direction for 
c                 each (y,z) pencil (in case of multiple waves)
c          nw_y - (int, input) number of systems to solve in y-direction for 
c                 each (x,z) pencil (ignored if dims<2)
c          nw_z - (int, input) number of systems to solve in z-direction for 
c                 each (x,y) pencil (ignored if dims<3)
c           ier - (int, output) return flag: 0=>success, otherwise error
c-----------------------------------------------------------------------
c======= Declarations =========
      implicit none

c     calling variables
      integer, intent(in)  :: comm3d, dims_
      integer, intent(in)  :: nx_, ny_, nz_, nw_x, nw_y, nw_z
      integer, intent(out) :: ier

c     local variables
      logical :: remain_dims(3)
      integer :: len, mpier, myid, nprocs, loc(1)

c======= Internals ============

c     initialize return value to success
      ier = 0

c     if we're already set up, just return
      if (comms_set)  return

c     check for legal inputs
      if ((dims_ < 1) .or. (dims_ > 3)) then
         write(0,*) 'ptri_setup error: illegal dims =',dims_
         ier = -1
         return
      end if
      if (nx_ < 1) then
         write(0,*) 'ptri_setup error: illegal nx =',nx_
         ier = -1
         return
      end if
      if ((ny_ < 1) .and. (dims_ > 1)) then
         write(0,*) 'ptri_setup error: illegal ny =',ny_
         ier = -1
         return
      end if
      if ((nz_ < 1) .and. (dims_ > 2))then
         write(0,*) 'ptri_setup error: illegal nz =',nz_
         ier = -1
         return
      end if
      if (nw_x < 1) then
         write(0,*) 'ptri_setup error: illegal nw_x =',nw_x
         ier = -1
         return
      end if
      if ((nw_y < 1) .and. (dims_ > 1)) then
         write(0,*) 'ptri_setup error: illegal nw_y =',nw_y
         ier = -1
         return
      end if
      if ((nw_z < 1) .and. (dims_ > 2))then
         write(0,*) 'ptri_setup error: illegal nw_z =',nw_z
         ier = -1
         return
      end if

c     set inputs into module variables
      dims = dims_
      nx = nx_
      ny = ny_
      nz = nz_
      nwave_x = nw_x
      nwave_y = nw_y
      nwave_z = nw_z
      if (dims < 2) then 
         ny = 1
         nwave_y = 0
      end if
      if (dims < 3) then
         nz = 1
         nwave_z = 0
      end if

c     get information from 3D communicator
      call mpi_comm_size(comm3d, nprocs, mpier)
      call mpi_comm_rank(comm3d, myid, mpier)
      call mpi_cart_get(comm3d, 3, pdims, periods, pcoords, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup: mpi_cart_get error =',mpier
         ier = -1
         return
      end if
      
c     construct 1D communicators for this process

c        x-directional communicator
      remain_dims = (/ .true., .false., .false. /)
      call mpi_cart_sub(comm3d, remain_dims, x_comm, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_setup mpi_sub x-dir err =',mpier
         ier = -1
         return
      endif

c        y-directional communicator
      if (dims > 1) then
         remain_dims = (/ .false., .true., .false. /)
         call mpi_cart_sub(comm3d, remain_dims, y_comm, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_setup mpi_sub y-dir err =',mpier
            ier = -1
            return
         endif
      end if

c        z-directional communicator
      if (dims > 2) then
         remain_dims = (/ .false., .false., .true. /)
         call mpi_cart_sub(comm3d, remain_dims, z_comm, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,'ptri_setup mpi_sub z-dir err =',mpier
            ier = -1
            return
         endif
      end if



c     get my locations in each directional communicator

c        x-directional rank
      call mpi_comm_rank(x_comm, myid_x, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_setup mpi_rank x-dir err =',mpier
         ier = -1
         return
      endif

c        y-directional rank
      if (dims > 1) then
         call mpi_comm_rank(y_comm, myid_y, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_setup mpi_rank y-dir err =',mpier
            ier = -1
            return
         endif  
      else
         myid_y = 0
      end if

c        z-directional rank
      if (dims > 2) then
         call mpi_comm_rank(z_comm, myid_z, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_setup mpi_rank z-dir err =',mpier
            ier = -1
            return
         endif  
      else
         myid_z = 0
      end if


c     get neighbor processor information
      pxl = MPI_PROC_NULL
      pxr = MPI_PROC_NULL
      pyl = MPI_PROC_NULL
      pyr = MPI_PROC_NULL
      pzl = MPI_PROC_NULL
      pzr = MPI_PROC_NULL

c        x-neighbors
      if (myid_x /= 0) then
         loc(1) = myid_x-1
         call mpi_cart_rank(x_comm, loc, pxl, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_setup mpi_rank xl err =',mpier
            ier = -1
            return
         endif  
      end if
      if (myid_x /= pdims(1)-1) then
         loc(1) = myid_x+1
         call mpi_cart_rank(x_comm, loc, pxr, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_setup mpi_rank xr err =',mpier
            ier = -1
            return
         endif  
      end if

c        y-neighbors
      if (dims > 1) then
         if (myid_y /= 0) then
            loc(1) = myid_y-1
            call mpi_cart_rank(y_comm, loc, pyl, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_setup mpi_rank yl err =',mpier
               ier = -1
               return
            endif  
         end if
         if (myid_y /= pdims(2)-1) then
            loc(1) = myid_y+1
            call mpi_cart_rank(y_comm, loc, pyr, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_setup mpi_rank yr err =',mpier
               ier = -1
               return
            endif  
         end if
      end if

c        z-neighbors
      if (dims > 2) then
         if (myid_z /= 0) then
            loc(1) = myid_z-1
            call mpi_cart_rank(z_comm, loc, pzl, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_setup mpi_rank zl err =',mpier
               ier = -1
               return
            endif  
         end if
         if (myid_z /= pdims(3)-1) then
            loc(1) = myid_z+1
            call mpi_cart_rank(z_comm, loc, pzr, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_setup mpi_rank zr err =',mpier
               ier = -1
               return
            endif  
         end if
      end if


c     allocate temporary arrays
      len = max(nx,ny)
      len = max(len,nz)
      allocate(b(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup b allocation error =',mpier
         ier = -1
         return
      endif  
      allocate(c(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup c allocation error =',mpier
         ier = -1
         return
      endif  
      allocate(a(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup a allocation error =',mpier
         ier = -1
         return
      endif  
      allocate(r(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup r allocation error =',mpier
         ier = -1
         return
      endif  
      allocate(s(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup s allocation error =',mpier
         ier = -1
         return
      endif  
      allocate(tmp1(nx*ny*nz), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup tmp1 alloc error =',mpier
         ier = -1
         return
      endif  
      allocate(tmp2(nx*ny*nz), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup tmp2 alloc error =',mpier
         ier = -1
         return
      endif  
      allocate(tmp3(nx*ny*nz), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_setup tmp3 alloc error =',mpier
         ier = -1
         return
      endif  
      
c     set initialization flag to true
      comms_set = .true.
      
      return
      end subroutine ptri_setup
c-----------------------------------------------------------------------




      subroutine ptri_free(ier)
c-----------------------------------------------------------------------
c     Description: deallocates all internal data
c-----------------------------------------------------------------------
      implicit none
      integer :: ierr, ier
      
      ier = 0
      if (allocated(b))     deallocate(b, stat=ierr)
      ier = ier + ierr
      if (allocated(c))     deallocate(c, stat=ierr)
      ier = ier + ierr
      if (allocated(a))     deallocate(a, stat=ierr)
      ier = ier + ierr
      if (allocated(r))     deallocate(r, stat=ierr)
      ier = ier + ierr
      if (allocated(s))     deallocate(s, stat=ierr)
      ier = ier + ierr
      if (allocated(tmp1))  deallocate(tmp1, stat=ierr)
      ier = ier + ierr
      if (allocated(tmp2))  deallocate(tmp2, stat=ierr)
      ier = ier + ierr
      if (allocated(tmp3))  deallocate(tmp3, stat=ierr)
      ier = ier + ierr
      if (ier /= 0)  write(0,*) 'ptri_free deallocation error =',ier
      return

      end subroutine ptri_free
c-----------------------------------------------------------------------




      subroutine ptri_sol(vec, order, ier)
c-----------------------------------------------------------------------
c     Description: provides the driver for the consecutive periodic 
c     tridiagonal solves, 
c            T1 * T2 * T3 * sol = rhs
c
c     Arguments:
c           vec - (dbl(*)) system rhs on input, solution on output; 
c                 this is just passed to the user constructor routine
c                 for each sub-system
c         order - (int(*), input) order for solving the three 
c                 directional systems, with 0 denoting to stop:
c                     (1,2,3)  =>  Tz * Ty * Tx * sol = rhs
c                     (3,1,2)  =>  Ty * Tx * Tz * sol = rhs
c                     (1,3,2)  =>  Ty * Tz * Tx * sol = rhs
c                     (1,0,0)  =>  Tx * sol = rhs
c                     (2,0,0)  =>  Ty * sol = rhs
c                     (0,0,0)  =>  sol = rhs
c                 Note: 2 and 3 are disabled if dims<2 or dims<3, 
c                 respectively.
c           ier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Declarations =========
      use profiling
      implicit none

c     calling variables
      real*8  :: vec(*)
      integer :: order(3)
      integer, intent(out) :: ier

c     internal variables
      integer :: i

c======= Internals ============

c     initialize return flag to success
      ier = 0

c     check that order has allowable arguments
      do i=1,3
         if ((order(i) < 0) .or. (order(i) > dims)) then
            write(0,'(3(A,i2))') 'ptri_sol warning: order(',i,') =',
     &           order(i),' is not between 0 and ',dims,', setting to 0'
            order(i) = 0
         endif
      enddo

c     solve for first direction
      if (order(1) == 0)  return
      call prof_start(20)
      call ptri_solve1d(order(1), vec, ier)
      if (ier /= 0) then
         write(0,'(A,i4)') 'ptri_sol first solve error =',ier
         return
      end if      
      call prof_stop(20)

c     solve for second direction
      if (order(2) == 0)  return
      call prof_start(21)
      call ptri_solve1d(order(2), vec, ier)
      if (ier /= 0) then
         write(0,'(A,i4)') 'ptri_sol second solve error =',ier
         return
      end if      
      call prof_stop(21)

c     solve for third direction
      if (order(3) == 0)  return
      call prof_start(22)
      call ptri_solve1d(order(3), vec, ier)
      if (ier /= 0) then
         write(0,'(A,i4)') 'ptri_sol third solve error =',ier
         return
      end if      
      call prof_stop(22)

      return
      end subroutine ptri_sol
c-----------------------------------------------------------------------




      subroutine ptri_solve1d(idir, vec, ier)
c-----------------------------------------------------------------------
c     Description: sets up and solves each set of directional systems.
c     
c     Arguments:
c        idir - (int, input) solution direction {1, 2 or 3}
c         vec - (dbl(*)) on input, system rhs; on output, system sol
c         ier - (int, output) return flag: 0=>success, 
c                1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Declarations =========
      implicit none

c     calling variables
      double precision     :: vec(*)
      integer, intent(in)  :: idir
      integer, intent(out) :: ier
      integer, pointer     :: comm
      
c     local variables
      integer :: ihi, jhi, khi, nw, i, j, k, iw, off, is, ie
      integer :: msg_id, myid, pl, pr
      integer :: status(MPI_STATUS_SIZE)
      real*8  :: rbuf(4), sbuf(4)
      real*8  :: c_l, a_l, a_r, a_n, r_l, r_r, r_n, s1

c======= Internals ============

c     set return flag to success
      ier = 0

c     set solver parameters for this direction
      if (idir == 1) then
         ihi = nx
         jhi = ny
         khi = nz
         nw = nwave_x
         comm => x_comm
         myid = myid_x
         pl = pxl
         pr = pxr
         call mpi_barrier(x_comm, ier)
         if (ier /= MPI_SUCCESS) then
            write(0,*) 'ptri_solve1d: mpi_barrier error =',ier
            return
         end if
      else if (idir == 2) then
         ihi = ny
         jhi = nx
         khi = nz
         nw = nwave_y
         comm => y_comm
         myid = myid_y
         pl = pyl
         pr = pyr
         call mpi_barrier(y_comm, ier)
         if (ier /= MPI_SUCCESS) then
            write(0,*) 'ptri_solve1d: mpi_barrier error =',ier
            return
         end if
      else
         ihi = nz
         jhi = nx
         khi = ny
         nw = nwave_z
         comm => z_comm
         myid = myid_z
         pl = pzl
         pr = pzr
         call mpi_barrier(z_comm, ier)
         if (ier /= MPI_SUCCESS) then
            write(0,*) 'ptri_solve1d: mpi_barrier error =',ier
            return
         end if
      end if
      
c     call user routine to prepare vec for this direction
      call prepare_rhs(idir, vec, ier)
      if (ier /= 0) then
         write(0,*) 'ptri_solve1d: prepare_rhs error =',ier
         return
      end if

c     iterate over the desired number of waves
      do iw=1,nw

c        perform initial forward phase of solver on strips
         do k=1,khi
            do j=1,jhi

c              call user routine to set up system on this strip for this wave
c                      r(i) = a(i)*s(i-1) + b(i)*s(i) + c(i)*s(i+1)
               call setup_strip(idir, j, k, iw, vec, a, b, c, r, ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_solve1d: setup_strip error =',ier
                  return
               end if

c              set unique ID for this pencil
               msg_id = k*jhi+j

c                 if not the first proc, wait for input from left proc
               if (pl /= MPI_PROC_NULL) then
                  call mpi_recv(rbuf, 3, MPI_REAL8, pl, msg_id+myid, 
     &                          comm, status, ier)
                  if (ier /= MPI_SUCCESS) then
                     write(*,*) 'ptri_solve1d: mpi_recv error =',ier
                     return
                  endif
                  a_l = rbuf(1)
                  c_l = rbuf(2)
                  r_l = rbuf(3)

c                 eliminate subdiagonal in first col, storing new last column in a
                  b(1) = b(1)-a(1)*c_l
                  r(1) = (r(1)-a(1)*r_l)/b(1)
                  a(1) = -a(1)*a_l/b(1)
                  c(1) = c(1)/b(1)

               else

c                 normalize the first row
                  a(1) = a(1)/b(1)
                  c(1) = c(1)/b(1)
                  r(1) = r(1)/b(1)

               endif

c                 eliminate subdiagonal, storing new last column in a
               ie = ihi
               if (pr == MPI_PROC_NULL)  ie = ihi-2
               do i = 2, ie
                  b(i) = b(i)-a(i)*c(i-1)
                  r(i) = (r(i)-a(i)*r(i-1))/b(i)
                  a(i) = -a(i)*a(i-1)/b(i)
                  c(i) = c(i)/b(i)
               end do

c                 send right-most b,c,r values to right neighbor
               if (pr /= MPI_PROC_NULL) then
                  sbuf(1) = a(ihi)
                  sbuf(2) = c(ihi)
                  sbuf(3) = r(ihi)
                  call mpi_send(sbuf, 3, MPI_REAL8, pr, msg_id+pr, 
     &                          comm, ier)
                  if (ier /= MPI_SUCCESS) then
                     write(*,*) 'ptri_solve1d: mpi_send error =',ier
                     return
                  end if

               else

c                 eliminate subdiagonal in last two rows
                  b(ihi-1) = b(ihi-1) - a(ihi-1)*c(ihi-2)
                  r(ihi-1) = (r(ihi-1) - a(ihi-1)*r(ihi-2))/b(ihi-1)
                  a(ihi-1) = (c(ihi-1) - a(ihi-1)*a(ihi-2))/b(ihi-1)
                  b(ihi) = b(ihi) - a(ihi-1)*a(ihi)
                  r(ihi) = (r(ihi) - r(ihi-1)*a(ihi))/b(ihi)
                  a(ihi) = c(ihi)/b(ihi)
                  
c                 eliminate superdiagonal in next to last row
                  r(ihi-1) = r(ihi-1) - a(ihi-1)*r(ihi)
                  a(ihi-1) = -a(ihi-1)*a(ihi)

               end if

c              store lower and upper diagonals in tmp1, tmp2; store rhs in tmp3
               off = ((k-1)*jhi + j-1)*ihi
               tmp1(off+1:off+ihi) = a(1:ihi)
               tmp2(off+1:off+ihi) = c(1:ihi)
               tmp3(off+1:off+ihi) = r(1:ihi)

c           end loop over directional strips
            end do
         end do



c        perform reverse phase of solver on strips
         do k=1,khi
            do j=1,jhi

c              unpack lower, upper diagonals, rhs
               off = ((k-1)*jhi + j-1)*ihi
               a(1:ihi) = tmp1(off+1:off+ihi)
               c(1:ihi) = tmp2(off+1:off+ihi)
               r(1:ihi) = tmp3(off+1:off+ihi)

c              set unique ID for this pencil
               msg_id = (khi+k)*jhi+j

c                 if not the last proc, wait for input from right proc
               if (pr /= MPI_PROC_NULL) then
                  call mpi_recv(rbuf, 4, MPI_REAL8, pr, 
     &                          msg_id+pr, comm, status, ier)
                  if (ier /= MPI_SUCCESS) then
                     write(*,*) 'ptri_solve1d: mpi_recv error =',ier
                     return
                  end if
                  a_r = rbuf(1)
                  a_n = rbuf(2)
                  r_r = rbuf(3)
                  r_n = rbuf(4)
                  
c                 eliminate last 2 superdiagonals and last col, storing new 1st col in a
                  r(ihi) = r(ihi) - c(ihi)*r_r - a(ihi)*r_n;
                  a(ihi) = -c(ihi)*a_r - a(ihi)*a_n;
                  r(ihi-1) = r(ihi-1) - c(ihi-1)*r(ihi) - a(ihi-1)*r_n;
                  a(ihi-1) = -c(ihi-1)*a(ihi) - a(ihi-1)*a_n;

               else

c                 otherwise just set these variables appropriately
                  a_n = a(ihi)
                  r_n = r(ihi)

               end if
      
c                 eliminate superdiagonal and last col, storing new 1st col in a
               do i = ihi-2,2,-1
                  r(i) = r(i) - c(i)*r(i+1) - a(i)*r_n
                  a(i) = -c(i)*a(i+1) - a(i)*a_n
               end do
               if (pl /= MPI_PROC_NULL) then
                  r(1) = r(1) - c(1)*r(2) - a(1)*r_n
                  a(1) = -c(1)*a(2) - a(1)*a_n
               endif

c                 send left-most a, r values to left neighbor
               if (pl /= MPI_PROC_NULL) then
                  sbuf(1) = a(1)
                  sbuf(2) = a_n
                  sbuf(3) = r(1)
                  sbuf(4) = r_n
                  call mpi_send(sbuf, 4, MPI_REAL8, pl, msg_id+myid, 
     &                          comm, ier)
                  if (ier /= MPI_SUCCESS) then
                     write(*,*) 'ptri_solve1d: mpi_send error =',ier
                     return
                  end if
               else
         
c                 first proc computes s at left end of overall domain
                  s(1) = (r(1) - c(1)*r(2) - a(1)*r_n)
     &                 / (1.d0 - c(1)*a(2) - a(1)*a_n)

               endif

c              store lower diagonal in tmp1, sol in tmp2, rhs in tmp3
               off = ((k-1)*jhi + j-1)*ihi
               tmp1(off+1:off+ihi) = a(1:ihi)
               tmp2(off+1:off+ihi) = s(1:ihi)
               tmp3(off+1:off+ihi) = r(1:ihi)

c           end loop over directional strips
            end do
         end do



c        perform final forward phase over strips
         do k=1,khi
            do j=1,jhi

c              unpack lower diagonal, sol, rhs
               off = ((k-1)*jhi + j-1)*ihi
               a(1:ihi) = tmp1(off+1:off+ihi)
               s(1:ihi) = tmp2(off+1:off+ihi)
               r(1:ihi) = tmp3(off+1:off+ihi)

c              set unique ID for this pencil
               msg_id = (2*khi+k)*jhi+j

c                 receive s left end of overall domain
               if (pl /= MPI_PROC_NULL) then
                  call mpi_recv(s1, 1, MPI_REAL8, pl, msg_id+myid, 
     &                          comm, status, ier)
                  if (ier /= MPI_SUCCESS) then
                     write(*,*) 'ptri_solve1d: mpi_recv error =',ier
                     return
                  endif
                  s(1) = r(1) - a(1)*s1
               else
                  s1 = s(1)
               endif

c                 send s at left end of overall domain
               if (pr /= MPI_PROC_NULL) then
                  call mpi_send(s1, 1, MPI_REAL8, pr, msg_id+pr, 
     &                          comm, ier)
                  if (ier /= MPI_SUCCESS) then
                     write(*,*) 'ptri_solve1d: mpi_send error =',ier
                     return
                  endif
               endif
               
c                 perform local forward sweep, computing solution
               do i = 2,ihi
                  s(i) = r(i) - a(i)*s1
               end do

c              call user routine to store solution on this strip for this wave
               call pack_strip(idir, j, k, iw, vec, s, ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_solve1d: pack_strip error =',ier
                  return
               end if

c           end loop over x-directional strips
            end do
         end do

c     end loop over waves
      end do


c     call user routine to process solution for this direction
      call process_sol(idir, vec, ier)
      if (ier /= 0) then
         write(0,*) 'ptri_solve1d: process_sol error =',ier
         return
      end if

c     nullify comm pointer
      nullify(comm)

      return

      end subroutine ptri_solve1d
c-----------------------------------------------------------------------




      end module ptri_solver
c=======================================================================
