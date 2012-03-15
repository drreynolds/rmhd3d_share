c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     SMU Mathematics
c     6/3/2009
c
c     This file implements a parallelized routine for ADI-like products
c             r = T1 * T2 * T3 * x
c     where each of the T matrices are periodic tridiagonal linear 
c     systems of the form
c             r(i) = L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1).
c     The implicit assumption is that each processor will own a large
c     number of these systems, within either a 1D, 2D or 3D processor 
c     topology.
c
c     Once these products are completed in one direction, the routine
c     proceeds to the next coordinate direction and repeats the process.
c
c     We have a number of options:
c
c     - Dimension - This solver may be used for 1, 2 or 3 dimensional 
c       systems.
c
c     - Customized product order - A user may specify the order of the 
c       three (or less) products, so that they proceed in a desired order, 
c       e.g. (x,z,y), (z,y,x), (x), (y), or even (x,x).
c
c     - User-provided pre-processing, post-processing and system setup 
c       routines - The user is expected to provide the following routines 
c       that are called by this solver:
c
c       mult_strip_setup(idir, j, k, vec, L, D, U, s, ier) -- Given the 
c          description of which 1D tridiagonal system to set up (j,k
c          correspond to the 2 orthogonal spatial indices), and the 
c          user-supplied multiplying array vec, this routine should 
c          return the arrays D, U, L and r for the tridiagonal system 
c             r(i) = L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1)
c          on this processor.
c
c       mult_strip_pack(idir, j, k, vec, r, ier) -- Given the description 
c          of which 1D tridiagonal system has just been multiplied, the 
c          computed product r, and the user-supplied solution array vec 
c          (overwrites the user-supplied array), this routine should 
c          pack the 1D product values from r into vec.
c
c       All of these routines should also return the success/failure 
c       flag ier, where:
c              0 -> success
c              1 -> recoverable failure
c             -1 -> fatal failure
c
c-----------------------------------------------------------------------
c     $Log: ptri_product.F,v $
c=======================================================================



      module ptri_product
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
c     1D communicators
      integer, target, private :: x_comm, y_comm, z_comm
c     parallelism information
      integer, private :: pdims(3), pcoords(3), myid_x, myid_y, myid_z
      integer, private :: pxl, pxr, pyl, pyr, pzl, pzr
      logical, private :: periods(3)
c     setup flag (pre-initialized)
      logical, private :: comms_set = .false.

c     temporary arrays
      real*8, dimension(:), allocatable, private :: D, U, L, r, s
c-----------------------------------------------------------------------



      contains

      
      
      subroutine ptri_mult_setup(comm3d, dims_, nx_, ny_, nz_, ier)
c-----------------------------------------------------------------------
c     Description: sets up the solvers
c
c     Arguments:
c        comm3d - (int, input) 3D cartesian MPI communicator
c         dims_ - (int, input) # of problem dimensions
c           nx_ - (int, input) local x-dimension length
c           ny_ - (int, input) local y-dimension length (ignored if dims<2)
c           nz_ - (int, input) local z-dimension length (ignored if dims<3)
c           ier - (int, output) return flag: 0=>success, otherwise error
c-----------------------------------------------------------------------
c======= Declarations =========
      implicit none

c     calling variables
      integer, intent(in)  :: comm3d, dims_
      integer, intent(in)  :: nx_, ny_, nz_
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
         write(0,*) 'ptri_mult_setup error: illegal dims =',dims_
         ier = -1
         return
      end if
      if (nx_ < 1) then
         write(0,*) 'ptri_mult_setup error: illegal nx =',nx_
         ier = -1
         return
      end if
      if ((ny_ < 1) .and. (dims_ > 1)) then
         write(0,*) 'ptri_mult_setup error: illegal ny =',ny_
         ier = -1
         return
      end if
      if ((nz_ < 1) .and. (dims_ > 2))then
         write(0,*) 'ptri_mult_setup error: illegal nz =',nz_
         ier = -1
         return
      end if

c     set inputs into module variables
      dims = dims_
      nx = nx_
      ny = ny_
      nz = nz_
      if (dims < 2)  ny = 1
      if (dims < 3)  nz = 1

c     get information from 3D communicator
      call mpi_comm_size(comm3d, nprocs, mpier)
      call mpi_comm_rank(comm3d, myid, mpier)

      call mpi_cart_get(comm3d, dims, pdims, periods, pcoords, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_mult_setup: cart_get error =',mpier
         ier = -1
         return
      end if
      
c     construct 1D communicators for this process

c        x-directional communicator
      remain_dims = (/ .true., .false., .false. /)
      call mpi_cart_sub(comm3d, remain_dims, x_comm, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_mult_setup cart_sub x err =',mpier
         ier = -1
         return
      endif

c        y-directional communicator
      if (dims > 1) then
      remain_dims = (/ .false., .true., .false. /)
      call mpi_cart_sub(comm3d, remain_dims, y_comm, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_mult_setup cart_sub y err =',mpier
         ier = -1
         return
      endif
      endif

c        z-directional communicator
      if (dims > 2) then
      remain_dims = (/ .false., .false., .true. /)
      call mpi_cart_sub(comm3d, remain_dims, z_comm, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_mult_setup cart_sub z err =',mpier
         ier = -1
         return
      endif
      endif



c     get my locations in each directional communicator

c        x-directional rank
      call mpi_comm_rank(x_comm, myid_x, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_mult_setup comm_rank x err =',mpier
         ier = -1
         return
      endif

c        y-directional rank
      if (dims > 1) then
      call mpi_comm_rank(y_comm, myid_y, mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) 'p',myid,' ptri_mult_setup comm_rank y err =',mpier
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
         write(0,*) 'p',myid,' ptri_mult_setup comm_rank z err =',mpier
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
      if ((myid_x /= 0) .or. periods(1)) then
         loc(1) = myid_x-1
         call mpi_cart_rank(x_comm, loc, pxl, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_mult_setup cart_rank xl err =',
     &           mpier
            ier = -1
            return
         endif  
      end if
      if ((myid_x /= pdims(1)-1) .or. periods(1)) then
         loc(1) = myid_x+1
         call mpi_cart_rank(x_comm, loc, pxr, mpier)
         if (mpier /= MPI_SUCCESS) then
            write(0,*) 'p',myid,' ptri_mult_setup cart_rank xr err =',
     &           mpier
            ier = -1
            return
         endif  
      end if

c        y-neighbors
      if (dims > 1) then
         if ((myid_y /= 0) .or. periods(2)) then
            loc(1) = myid_y-1
            call mpi_cart_rank(y_comm, loc, pyl, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_mult_setup cart_rank yl err ='
     &              ,mpier
               ier = -1
               return
            endif  
         end if
         if ((myid_y /= pdims(2)-1) .or. periods(2)) then
            loc(1) = myid_y+1
            call mpi_cart_rank(y_comm, loc, pyr, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_mult_setup cart_rank yr err ='
     &              ,mpier
               ier = -1
               return
            endif  
         end if
      end if

c        z-neighbors
      if (dims > 2) then
         if ((myid_z /= 0) .or. periods(3)) then
            loc(1) = myid_z-1
            call mpi_cart_rank(z_comm, loc, pzl, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_mult_setup cart_rank zl err ='
     &              ,mpier
               ier = -1
               return
            endif  
         end if
         if ((myid_z /= pdims(3)-1) .or. periods(3)) then
            loc(1) = myid_z+1
            call mpi_cart_rank(z_comm, loc, pzr, mpier)
            if (mpier /= MPI_SUCCESS) then
               write(0,*) 'p',myid,' ptri_mult_setup cart_rank zr err ='
     &              ,mpier
               ier = -1
               return
            endif  
         end if
      end if


c     allocate temporary arrays
      len = max(nx,ny)
      len = max(len,nz)
      allocate(D(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_mult_setup D alloc error =',mpier
         ier = -1
         return
      endif  
      allocate(U(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_mult_setup U alloc error =',mpier
         ier = -1
         return
      endif  
      allocate(L(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_mult_setup L alloc error =',mpier
         ier = -1
         return
      endif  
      allocate(r(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_mult_setup r alloc error =',mpier
         ier = -1
         return
      endif  
      allocate(s(len), stat=mpier)
      if (mpier /= MPI_SUCCESS) then
         write(0,*) ' p',myid,' ptri_mult_setup s alloc error =',mpier
         ier = -1
         return
      endif  
      
c     set initialization flag to true
      comms_set = .true.

      return
      end subroutine ptri_mult_setup
c-----------------------------------------------------------------------




      subroutine ptri_mult_free(ier)
c-----------------------------------------------------------------------
c     Description: deallocates all internal data
c-----------------------------------------------------------------------
      implicit none
      integer :: mpier, ier
      
      deallocate(D, U, L, r, s, stat=mpier)
      if (mpier /= 0) then
         write(0,*) 'ptri_mult_free deallocation error =',mpier
         ier = 1
      endif  
      return

      end subroutine ptri_mult_free
c-----------------------------------------------------------------------




      subroutine ptri_mult(vec, order, ier)
c-----------------------------------------------------------------------
c     Description: provides the driver for the consecutive periodic 
c     tridiagonal matrix-vector product
c            rvec = T1 * T2 * T3 * xvec
c
c     Arguments:
c           vec - (dbl(*)) xvec on input, rvec on output
c         order - (int(*), input) order for multiplying the three 
c                 directional systems, with 0 denoting to stop:
c                     (1,2,3)  =>  rvec = Tz * Ty * Tx * xvec
c                     (3,1,2)  =>  rvec = Ty * Tx * Tz * xvec
c                     (1,3,2)  =>  rvec = Ty * Tz * Tx * xvec
c                     (1,0,0)  =>  rvec = Tx * xvec
c                     (2,0,0)  =>  rvec = Ty * xvec
c                     (0,0,0)  =>  rvec = xvec
c                 Note: 2 and 3 are disabled if dims<2 or dims<3, 
c                 respectively.
c           ier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Declarations =========
      implicit none

c     calling variables
      real*8 :: vec(*)
      integer, intent(in)  :: order(3)
      integer, intent(out) :: ier

c     internal variables
      integer :: mpier

c======= Internals ============

c     initialize return flag to success
      ier = 0

c     check that order has allowable arguments
      if ((order(1) < 0) .or. (order(1) > dims)) then
         write(0,'(A,i)') 'ptri_mult error: illegal order(1) =',order(1)
         ier = -1
         return
      end if
      if (order(1) /= 0) then
         if ((order(2)<0) .or. (order(2)>dims)) then
            write(0,'(A,i)') 'ptri_mult error: illegal order(2) =',
     &           order(2)
            ier = -1
            return
         end if
         if (order(2) /= 0) then
            if ((order(3)<0) .or. (order(3)>dims)) then
               write(0,'(A,i)') 'ptri_mult error: illegal order(3) =',
     &              order(2)
               ier = -1
               return
            end if
         end if
      end if

c     multiply for first direction
      if (order(1) == 0)  return
      call ptri_mult1d(order(1), vec, mpier)
      if (ier /= 0) then
         write(0,'(A,i4)') 'ptri_mult first multiply error =',mpier
         ier = mpier
         return
      end if      

c     multiply for second direction
      if (order(2) == 0)  return
      call ptri_mult1d(order(2), vec, mpier)
      if (ier /= 0) then
         write(0,'(A,i4)') 'ptri_mult second multiply error =',mpier
         ier = mpier
         return
      end if      

c     multiply for third direction
      if (order(3) == 0)  return
      call ptri_mult1d(order(3), vec, mpier)
      if (ier /= 0) then
         write(0,'(A,i4)') 'ptri_mult third multiply error =',mpier
         ier = mpier
         return
      end if      

      return
      end subroutine ptri_mult
c-----------------------------------------------------------------------




      subroutine ptri_mult1d(idir, vec, ier)
c-----------------------------------------------------------------------
c     Description: sets up and multiplies each set of directional 
c       systems,
c            rvec = T * xvec
c
c     Arguments:
c        idir - (int, input) solution direction {1, 2 or 3}
c         vec - (dbl(*)) xvec on input, rvec on output
c         ier - (int, output) return flag: 0=>success, 
c                1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Declarations =========
      implicit none

c     calling variables
      real*8               :: vec(*)
      integer, intent(in)  :: idir
      integer, intent(out) :: ier
      integer, pointer     :: comm
      
c     local variables
      integer :: ihi, jhi, khi, i, j, k, iw, off, is, ie
      integer :: msg_id, iprocs, myid, pl, pr
      integer :: id_recvL, id_recvR, id_sendL, id_sendR
      integer :: status(MPI_STATUS_SIZE)
      real*8  :: sbuf, sL, sR

c======= Internals ============

c     set return flag to success
      ier = 0

c     set multiply parameters for this direction
      if (idir == 1) then
         ihi = nx
         jhi = ny
         khi = nz
         comm => x_comm
         myid = myid_x
         pl = pxl
         pr = pxr
         call mpi_comm_size(comm, iprocs, ier)
         if (ier /= MPI_SUCCESS) then
            write(0,*) 'ptri_mult1d: mpi_comm_size error =',ier
            return
         end if
      else if (idir == 2) then
         ihi = ny
         jhi = nx
         khi = nz
         comm => y_comm
         myid = myid_y
         pl = pyl
         pr = pyr
         call mpi_comm_size(comm, iprocs, ier)
         if (ier /= MPI_SUCCESS) then
            write(0,*) 'ptri_mult1d: mpi_comm_size error =',ier
            return
         end if
      else
         ihi = nz
         jhi = nx
         khi = ny
         comm => z_comm
         myid = myid_z
         pl = pzl
         pr = pzr
         call mpi_comm_size(comm, iprocs, ier)
         if (ier /= MPI_SUCCESS) then
            write(0,*) 'ptri_mult1d: mpi_comm_size error =',ier
            return
         end if
      end if
      
c     perform multiplication on strips
      do k=1,khi
         do j=1,jhi

c           call user routine to set up system on this strip 
c              r(i) = L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1)
            call mult_strip_setup(idir, j, k, vec, L, D, U, s, ier)
            if (ier /= 0) then
               write(0,*) 'ptri_mult1d: mult_strip_setup error =',ier
               return
            end if

c           perform exchange to get sL and sR
            sL = 0.d0
            sR = 0.d0
            if (pl /= MPI_PROC_NULL) then
               call mpi_irecv(sL,1,MPI_REAL8,pl,100,comm,id_recvL,ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_mult1d: mpi_irecv (L) error =',ier
                  return
               end if
               call mpi_isend(s(1),1,MPI_REAL8,pl,101,comm,id_sendL,ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_mult1d: mpi_isend (L) error =',ier
                  return
               end if
            end if
            if (pr /= MPI_PROC_NULL) then
               call mpi_irecv(sR,1,MPI_REAL8,pr,101,comm,id_recvR,ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_mult1d: mpi_irecv (R) error =',ier
                  return
               end if
               call mpi_isend(s(ihi),1,MPI_REAL8,pr,100,comm,id_sendR,
     &                        ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_mult1d: mpi_isend (R) error =',ier
                  return
               end if
            end if

c           wait for arrival
            if (pl /= MPI_PROC_NULL) then
               call mpi_wait(id_recvL, status, ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_mult1d: mpi_wait (L) error =',ier
                  return
               end if
            end if
            if (pr /= MPI_PROC_NULL) then
               call mpi_wait(id_recvR, status, ier)
               if (ier /= 0) then
                  write(0,*) 'ptri_mult1d: mpi_wait (R) error =',ier
                  return
               end if
            end if


c           perform local matrix-vector multiply
            r(1) = L(1)*sL + D(1)*s(1) + U(1)*s(2)
            do i=2,ihi-1
               r(i) = L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1)
            end do
            r(ihi) = L(ihi)*s(ihi-1) + D(ihi)*s(ihi) + U(ihi)*sR


c           call user routine to store solution on this strip
            call mult_strip_pack(idir, j, k, vec, r, ier)
            if (ier /= 0) then
               write(0,*) 'ptri_mult1d: mult_strip_pack error =',ier
               return
            end if

c        end loop over strips
         end do
      end do

      return

      end subroutine ptri_mult1d
c-----------------------------------------------------------------------





      end module ptri_product
c=======================================================================
