c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     SMU Mathematics
c
c     This file implements a preconditioner in which we only treat the 
c     fastest linear waves in the hyperbolic portion of the resistive
c     MHD system.  The idea is as follows: we first decompose the 
c     hyperbolic portion of the RMHD system into x, y and z-directional
c     components, based on the pieces of the divergence terms.  We then 
c     precondition the overall system by using these individual pieces
c     in an operator-split fashion.
c
c-----------------------------------------------------------------------
c     $Log: fastwave_prec.F,v $
c=======================================================================




      module fastwaveprec_data
c-----------------------------------------------------------------------
c     Description: holds necessary data describing the preconditioning 
c     method and system, local reusable data for the preconditioner 
c     solve routine, and data stored by the preconditioner setup 
c     routine.
c
c     Contains the general preconditioning variables:
c
c        NWaves - number of fast waves the precond. should handle
c        wvindx - array of wave numbers in order of decreasing 
c                 magnitude
c
c     Contains the reusable psolve variables 
c     (*_x indicate x-directional variables, similar for *_y, *_z):
c           alamda - array of eigenvalues on spatial strip
c              evl - array of left eigenvectors on spatial strip
c              evr - array of right eigenvectors on spatial strip
c-----------------------------------------------------------------------
      use mesh
      use ptri_solver

      implicit none
      save
      
c     customization parameters (with defaults)
      integer :: NWaves  = 8
      integer :: ordering = 0
c$$$      logical :: ADICorr = .false.

c     internal module vars, set via fastwave_pset routine
      integer, dimension(8) :: wvindx = (/5, 6, 3, 4, 1, 2, 7, 8/)
      real*8  :: dxi2, dyi2, dzi2, gamdt
      logical :: ptri_set = .false.

#ifndef ONE_D
c     permutation matrices for transitioning between local/global coords
      real*8, dimension(nvar,nvar) :: PermY, PermYt
#ifndef TWO_D
      real*8, dimension(nvar,nvar) :: PermZ, PermZt
#endif
#endif

c     directional eigen-components
#ifndef DYNAMIC    
      real*8 :: alamda_x(nvar,0:nxlsize+1,1:nylsize,1:nzlsize)
      real*8 :: evl_x(nvar,nvar,0:nxlsize+1,1:nylsize,1:nzlsize), 
     &     evr_x(nvar,nvar,0:nxlsize+1,1:nylsize,1:nzlsize)
#ifndef ONE_D
      real*8 :: alamda_y(nvar,0:nylsize+1,1:nxlsize,1:nzlsize)
      real*8 :: evl_y(nvar,nvar,0:nylsize+1,1:nxlsize,1:nzlsize), 
     &     evr_y(nvar,nvar,0:nylsize+1,1:nxlsize,1:nzlsize)
#ifndef TWO_D
      real*8 :: alamda_z(nvar,0:nzlsize+1,1:nxlsize,1:nylsize)
      real*8 :: evl_z(nvar,nvar,0:nzlsize+1,1:nxlsize,1:nylsize), 
     &     evr_z(nvar,nvar,0:nzlsize+1,1:nxlsize,1:nylsize)
#endif
#endif

c     storage for block correction matrix and solves
      real*8 :: Mcorr(nvar,nvar,1:nxlsize,1:nylsize,1:nzlsize)

#else
      real*8, allocatable :: alamda_x(:,:,:,:)
      real*8, dimension(:,:,:,:,:), allocatable :: evl_x, evr_x
#ifndef ONE_D
      real*8, allocatable :: alamda_y(:,:,:,:)
      real*8, dimension(:,:,:,:,:), allocatable :: evl_y, evr_y
#ifndef TWO_D
      real*8, allocatable :: alamda_z(:,:,:,:)
      real*8, dimension(:,:,:,:,:), allocatable :: evl_z, evr_z
#endif
#endif

c     storage for block correction matrix and solves
      real*8, allocatable :: Mcorr(:,:,:,:,:)

#endif
cc endif DYNAMIC

      end module fastwaveprec_data
c=======================================================================




      subroutine fastwave_pset(uu, fu, gam, vv, tmp, ier)
c-----------------------------------------------------------------------
c     Description: fastwave_pset provides the fast wave preconditioner 
c        setup routine 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c            fu - (dbl(*), input) rhs of ODE, i.e. f(u,t)
c           gam - (dbl, input) the scalar appearing in the Newton matrix
c                  M = I+gam*J
c            vv - (dbl(*)) scratch space with same size as uu
c           tmp - (dbl(*)) scratch space with same size as uu
c           ier - (int, output) return flag: 0=>success, otherwise error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use fastwaveprec_data
      use mpistuff

c======= Declarations =========
      implicit none

c     calling variables
      real*8, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &        :: uu, fu, vv, tmp
      real*8  :: gam
      integer :: ier

c     local variables
      integer :: MPierr, dims

c======= Internals ============

c     initialize return value to success
      ier = 0

c     if data structures unallocated, allocate them now
#ifdef DYNAMIC
      if (.not. allocated(alamda_x) ) then
c$$$         if (ADIcorr) 
c$$$     &        allocate(Mcorr(nvar,nvar,1:nxlsize,1:nylsize,1:nzlsize))
         allocate(alamda_x(nvar,0:nxlsize+1,1:nylsize,1:nzlsize))
         allocate(evl_x(nvar,nvar,0:nxlsize+1,1:nylsize,1:nzlsize))
         allocate(evr_x(nvar,nvar,0:nxlsize+1,1:nylsize,1:nzlsize))
#ifndef ONE_D
         allocate(alamda_y(nvar,0:nylsize+1,1:nxlsize,1:nzlsize))
         allocate(evl_y(nvar,nvar,0:nylsize+1,1:nxlsize,1:nzlsize))
         allocate(evr_y(nvar,nvar,0:nylsize+1,1:nxlsize,1:nzlsize))
#ifndef TWO_D
         allocate(alamda_z(nvar,0:nzlsize+1,1:nxlsize,1:nylsize))
         allocate(evl_z(nvar,nvar,0:nzlsize+1,1:nxlsize,1:nylsize))
         allocate(evr_z(nvar,nvar,0:nzlsize+1,1:nxlsize,1:nylsize))
#endif
#endif
      endif
#endif

c     set differencing increments and time step factor into module
      dxi2 = 0.5d0/dx
      dyi2 = 0.5d0/dy
      dzi2 = 0.5d0/dz
      gamdt = gam


#ifndef ONE_D
c     construct the permutation matrices for y and z directions
      PermY = 0.d0;       PermY(1,1) = 1.d0;  PermY(3,2) = 1.d0
      PermY(4,3) = 1.d0;  PermY(2,4) = 1.d0;  PermY(6,5) = 1.d0
      PermY(7,6) = 1.d0;  PermY(5,7) = 1.d0;  PermY(8,8) = 1.d0
      PermYt = 0.d0;      PermYt(1,1) = 1.d0; PermYt(2,3) = 1.d0
      PermYt(3,4) = 1.d0; PermYt(4,2) = 1.d0; PermYt(5,6) = 1.d0
      PermYt(6,7) = 1.d0; PermYt(7,5) = 1.d0; PermYt(8,8) = 1.d0
#ifndef TWO_D
      PermZ = 0.d0;       PermZ(1,1) = 1.d0;  PermZ(4,2) = 1.d0
      PermZ(2,3) = 1.d0;  PermZ(3,4) = 1.d0;  PermZ(7,5) = 1.d0
      PermZ(5,6) = 1.d0;  PermZ(6,7) = 1.d0;  PermZ(8,8) = 1.d0
      PermZt = 0.d0;      PermZt(1,1) = 1.d0; PermZt(2,4) = 1.d0
      PermZt(3,2) = 1.d0; PermZt(4,3) = 1.d0; PermZt(5,7) = 1.d0
      PermZt(6,5) = 1.d0; PermZt(7,6) = 1.d0; PermZt(8,8) = 1.d0
#endif
#endif


c     construct/save the directional eigen-decompositions
      call SetADIEigensystem(uu,vv)


c$$$c     construct the non-constant coefficient correction matrix 
c$$$c     (sans I + gamma*(), so that we can remain valid with the 1D solves)
c$$$      if (ADICorr) then
c$$$         if (ordering == 0) then
c$$$            call construct_FWcorrection0(ier)
c$$$         elseif (ordering == 1) then
c$$$            call construct_FWcorrection1(ier)
c$$$         else
c$$$            call construct_FWcorrection2(ier)
c$$$         endif
c$$$      endif


c     construct/save the parallelism information and 1D process 
c     communicators for this process
      if (ptri_set == .false.) then

         dims = 3
#ifdef ONE_D
         dims = 1
#endif
#ifdef TWO_D
         dims = 2
#endif

c        set up the parallel solvers
         call ptri_setup(comm3d, dims, nxlocal, nylocal, nzlocal, 
     &                   NWaves, NWaves, NWaves, ier)
         if (ier /= 0) then
            write(0,*) 'fastwave_pset: ptri_setup error = ',ier
            return
         endif

         ptri_set = .true.
      endif
      
      return
      end subroutine fastwave_pset
c=======================================================================




      subroutine fastwave_free()
c-----------------------------------------------------------------------
c     Description: fastwave_free deallocates all internal fastwave data
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use fastwaveprec_data

c======= Declarations =========
      implicit none
      integer :: ier
      
c======= Internals ============

      call ptri_free(ier)
      if (ier /= 0)  write(0,*) 'fastwave_free: ptri_free error =',ier

#ifdef DYNAMIC
      if (allocated(Mcorr))     deallocate(Mcorr)
      if (allocated(alamda_x))  deallocate(alamda_x)
      if (allocated(evl_x))     deallocate(evl_x)
      if (allocated(evr_x))     deallocate(evr_x)
#ifndef ONE_D
      if (allocated(alamda_y))  deallocate(alamda_y)
      if (allocated(evl_y))     deallocate(evl_y)
      if (allocated(evr_y))     deallocate(evr_y)
#ifndef TWO_D
      if (allocated(alamda_z))  deallocate(alamda_z)
      if (allocated(evl_z))     deallocate(evl_z)
      if (allocated(evr_z))     deallocate(evr_z)
#endif
#endif
#endif
      
      return
      end subroutine fastwave_free
c=======================================================================




      subroutine fastwave_psol(uu, fu, gam, ww, vv, kier)
c-----------------------------------------------------------------------
c     Description: fastwave_psol provides the fast wave preconditioner 
c        solve routine for the preconditioning system:  P z = r
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c            fu - (dbl(*), input) rhs of ODE, i.e. f(u,t)
c           gam - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gam*J
c            ww - (dbl(*), in/out) on input the rhs vector r, and on 
c                  output the solution vector z
c           vv - (dbl(*)) scratch space with same size as uu
c                 vv holds primitive variables
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use fastwaveprec_data

c======= Declarations =========
      implicit none

c     calling variables
      real*8, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &        :: uu, fu, ww, vv
      real*8  :: gam
      integer :: kier

c     internal variables
      integer :: order(3)


c======= Internals ============

c     initialize return flag to success
      kier = 0

c     set up solver order
      order = (/ 1, 2, 3 /)
#ifdef TWO_D
      order = (/ 1, 2, 0 /)
#endif
#ifdef ONE_D
      order = (/ 1, 0, 0 /)
#endif

c     call periodic tridiagonal solvers
      call ptri_sol(ww, order, kier)
      if (kier /= 0) then
         write(0,*) 'fastwave_psol: ptri_sol error =',kier
         return
      endif

c$$$c     Perform correction solve to account for the spatial variation 
c$$$c     in the left eigenvectors
c$$$      if (ADICorr) then
c$$$         call fastwave_correction(ww, kier)
c$$$         if (kier /= 0) 
c$$$     &        write(0,'(A,i4)') 'fastwave_psol: correction error =',kier
c$$$      endif

      return
      end subroutine fastwave_psol
c=======================================================================





      subroutine prepare_rhs(idir, ww, ier)
c-----------------------------------------------------------------------
c     Description: prepare_rhs converts the rhs to characteristic 
c     variables
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c          ww - (dbl(*)) rhs array
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common
      use fastwaveprec_data

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in)  :: idir
      integer, intent(out) :: ier
      real*8 :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      
c     Local variables
      integer :: i, j, k, l
      real*8  :: wloc(nvar)
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     project system 
      if (idir == 1) then

c        convert ww to characteristic space over domain
         do k=1,nzlocal
            do j=1,nylocal
               do i=0,nxlocal+1
                  do l=1,nvar
                     wloc(l) = ww(i,j,k,l)
                  enddo
                  wloc = matmul(evl_x(:,:,i,j,k), wloc)
                  do l=1,nvar
                     ww(i,j,k,l) = wloc(l)
                  enddo
               enddo
            enddo
         enddo            

      else if (idir == 2) then

#ifndef ONE_D
c        remap rhs to local coordinate system
         call RemapToLocalCoords3D(ww,2)
         
c        convert ww to characteristic space over domain
         do k=1,nzlocal
            do i=1,nxlocal
               do j=0,nylocal+1
                  do l=1,nvar
                     wloc(l) = ww(i,j,k,l)
                  enddo
                  wloc = matmul(evl_y(:,:,j,i,k), wloc)
                  do l=1,nvar
                     ww(i,j,k,l) = wloc(l)
                  enddo
               enddo
            enddo
         enddo            
#endif

      else

#ifndef ONE_D
#ifndef TWO_D
c        remap rhs to local coordinate system
         call RemapToLocalCoords3D(ww,3)

c        convert ww to characteristic space over domain
         do j=1,nylocal
            do i=1,nxlocal
               do k=0,nzlocal+1
                  do l=1,nvar
                     wloc(l) = ww(i,j,k,l)
                  end do
                  wloc = matmul(evl_z(:,:,k,i,j), wloc)
                  do l=1,nvar,1
                     ww(i,j,k,l) = wloc(l)
                  end do
               end do
            end do
         end do
#endif
#endif

      end if
     
      return

      end subroutine prepare_rhs
c=======================================================================





      subroutine process_sol(idir, ww, ier)
c-----------------------------------------------------------------------
c     Description: process_sol converts the solution to conserved 
c     variables.
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c          ww - (dbl(*)) rhs array
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common
      use fastwaveprec_data

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in)  :: idir
      integer, intent(out) :: ier
      real*8  :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      
c     Local variables
      integer :: i, j, k, l
      real*8  :: wloc(nvar)
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     project system 
      if (idir == 1) then

c        convert ww back to conserved vars over domain
         do k=1,nzlocal
            do j=1,nylocal
               do i=0,nxlocal+1
                  do l=1,nvar
                     wloc(l) = ww(i,j,k,l)
                  enddo
                  wloc = matmul(evr_x(:,:,i,j,k), wloc)
                  do l=1,nvar
                     ww(i,j,k,l) = wloc(l)
                  enddo
               enddo
            enddo
         enddo            

      else if (idir == 2) then

#ifndef ONE_D
c        convert ww back to conserved vars over domain
         do k=1,nzlocal
            do i=1,nxlocal
               do j=0,nylocal+1
                  do l=1,nvar
                     wloc(l) = ww(i,j,k,l)
                  enddo
                  wloc = matmul(evr_y(:,:,j,i,k), wloc)
                  do l=1,nvar
                     ww(i,j,k,l) = wloc(l)
                  enddo
               enddo
            enddo
         enddo            

c        remap to the global coordinate system
         call RemapToGlobalCoords3D(ww,2)
#endif

      else

#ifndef ONE_D
#ifndef TWO_D
c        convert ww back to conserved vars over domain
         do j=1,nylocal
            do i=1,nxlocal
               do k=0,nzlocal+1
                  do l=1,nvar
                     wloc(l) = ww(i,j,k,l)
                  enddo
                  wloc = matmul(evr_z(:,:,k,i,j), wloc)
                  do l=1,nvar
                     ww(i,j,k,l) = wloc(l)
                  enddo
               enddo
            enddo
         enddo            

c        remap to the global coordinate system
         call RemapToGlobalCoords3D(ww,3)
#endif
#endif

      end if

      return

      end subroutine process_sol
c=======================================================================





      subroutine setup_strip(idir, j, k, nw, ww, L, D, U, r, ier)
c-----------------------------------------------------------------------
c     Description: setup_strip fills in the periodic tridiagonal system 
c             r(i) = L(i)*s(i-1) + D(i)*s(i) + U(i)*s(i+1)
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c           j - (int, input) fast index orthogonal to idir
c           k - (int, input) slow index orthogonal to idir
c          nw - (int, input) wave number to work on
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
      use fastwaveprec_data
      use boundary_conds

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in) :: idir, j, k, nw
      integer, intent(out) :: ier
      real*8, dimension(*) :: D, U, L, r
      real*8 :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)

c     Local variables
      integer :: i, iw
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     get index of wave to work on
      iw = wvindx(nw)

c     set up system based on direction
      if (idir == 1) then

c        set r to be this wave rhs on the spatial strip
         do i=1,nxlocal
            r(i) = ww(i,j,k,iw)
         enddo

c        reset matrix arrays
         D(1:nxlocal) = 1.d0
         L(1:nxlocal) = 0.d0
         U(1:nxlocal) = 0.d0

c        set arrays for the spatial couplings on the block matrix corresponding 
c        to this characteristic variable and spatial strip
         if (ordering < 2) then
            do i=1,nxlocal
               L(i) = -gamdt*alamda_x(iw,i,j,k)*dxi2
               U(i) =  gamdt*alamda_x(iw,i,j,k)*dxi2
            enddo
         else
            do i=1,nxlocal
               L(i) = -gamdt*alamda_x(iw,i-1,j,k)*dxi2
               U(i) =  gamdt*alamda_x(iw,i+1,j,k)*dxi2
            enddo
         endif

c        undo couplings for boundary, non-periodic case
         if (xbc /= BCPeriodic) then
            if (iprocx == 1) then
               L(1) = 0.d0
               U(1) = 0.d0
            endif
            if (iprocx == xprocs) then
               L(nxlocal) = 0.d0
               U(nxlocal) = 0.d0
            endif
         endif

      else if (idir == 2) then

#ifndef ONE_D
c        set r to be this wave rhs on the spatial strip
         do i=1,nylocal
            r(i) = ww(j,i,k,iw)
         enddo

c        reset matrix arrays
         D(1:nylocal) = 1.d0
         L(1:nylocal) = 0.d0
         U(1:nylocal) = 0.d0

c        set arrays for the spatial couplings on the block matrix corresponding 
c        to this characteristic variable and spatial strip
         if (ordering < 2) then
            do i=1,nylocal
               L(i) = -gamdt*alamda_y(iw,i,j,k)*dyi2
               U(i) =  gamdt*alamda_y(iw,i,j,k)*dyi2
            enddo
         else
            do i=1,nylocal
               L(i) = -gamdt*alamda_y(iw,i-1,j,k)*dyi2
               U(i) =  gamdt*alamda_y(iw,i+1,j,k)*dyi2
            enddo
         endif
         
c        undo couplings for boundary, non-periodic case
         if (ybc /= BCPeriodic) then
            if (iprocy == 1) then
               L(1) = 0.d0
               U(1) = 0.d0
            endif
            if (iprocy == yprocs) then
               L(nylocal) = 0.d0
               U(nylocal) = 0.d0
            endif
         endif
#endif

      else

#ifndef ONE_D
#ifndef TWO_D
c        set r to be this wave rhs on the spatial strip
         do i=1,nzlocal
            r(i) = ww(j,k,i,iw)
         enddo

c        reset matrix arrays
         D(1:nzlocal) = 1.d0
         L(1:nzlocal) = 0.d0
         U(1:nzlocal) = 0.d0

c        set arrays for the spatial couplings on the block matrix corresponding 
c        to this characteristic variable and spatial strip
         if (ordering < 2) then
            do i=1,nzlocal
               L(i) = -gamdt*alamda_z(iw,i,j,k)*dzi2
               U(i) =  gamdt*alamda_z(iw,i,j,k)*dzi2
            enddo
         else
            do i=1,nzlocal
               L(i) = -gamdt*alamda_z(iw,i-1,j,k)*dzi2
               U(i) =  gamdt*alamda_z(iw,i+1,j,k)*dzi2
            enddo
         endif

c        undo couplings for boundary, non-periodic case
         if (zbc /= BCPeriodic) then
            if (iprocz == 1) then
               L(1) = 0.d0
               U(1) = 0.d0
            endif
            if (iprocz == zprocs) then
               L(nzlocal) = 0.d0
               U(nzlocal) = 0.d0
            endif
         endif
#endif
#endif

      end if

      return

      end subroutine setup_strip
c=======================================================================





      subroutine pack_strip(idir, j, k, nw, ww, sol, ier)
c-----------------------------------------------------------------------
c     Description: pack_strip puts the periodic tridiagonal system 
c     solution into the multi-dimensional ww array.
c     
c     Arguments:
c        idir - (int, input) direction for upcoming solve
c           j - (int, input) fast index orthogonal to idir
c           k - (int, input) slow index orthogonal to idir
c          nw - (int, input) wave number to work on
c          ww - (dbl(*), input) multi-dimensional sol array
c         sol - (dbl(*), input) sol for 1D system
c         ier - (int, output) success(0) or failure(1)
c-----------------------------------------------------------------------
c=======Inclusions ===========
      use mesh
      use mesh_common
      use fastwaveprec_data

c=======Declarations =========
      implicit none

c     calling variables
      integer, intent(in) :: idir, j, k, nw
      integer, intent(out) :: ier
      real*8, dimension(*) :: sol
      real*8 :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)

c     Local variables
      integer :: i, iw
      
c=======Internals ============

c     initialize output flag to success
      ier = 0

c     extract wave index to work on
      iw = wvindx(nw)

c     pack solution based on direction
      if (idir == 1) then

c        put solution back into ww
         do i=1,nxlocal
            ww(i,j,k,iw) = sol(i)
         enddo
         
      else if (idir == 2) then

c        put solution back into ww
         do i=1,nylocal
            ww(j,i,k,iw) = sol(i)
         enddo

      else 

c        put solution back into ww
         do i=1,nzlocal
            ww(j,k,i,iw) = sol(i)
         enddo

      end if

      return

      end subroutine pack_strip
c=======================================================================





c$$$      subroutine fastwave_correction(ww, ier)
c$$$c-----------------------------------------------------------------------
c$$$c     Description: fastwave_correction provides the preconditioner 
c$$$c     correction routine.
c$$$c     
c$$$c     Arguments:
c$$$c          ww - (dbl(*)) on input, rhs of correction system;
c$$$c                on output, sol of precond system
c$$$c         ier - (int, output) return flag: 0=>success, 
c$$$c                1=>recoverable error, -1=>non-recoverable error
c$$$c-----------------------------------------------------------------------
c$$$c=======Inclusions ===========
c$$$      use mesh
c$$$      use fastwaveprec_data
c$$$
c$$$c=======Declarations =========
c$$$      implicit none
c$$$
c$$$c     calling variables
c$$$      real*8  :: ww(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
c$$$      integer :: ier
c$$$
c$$$c     Local variables
c$$$      integer :: i, j, k, l, m, mier
c$$$      integer :: ipiv(nvar)
c$$$      real*8  :: Mtmp(nvar,nvar), wloc(nvar)
c$$$      
c$$$c=======Internals ============
c$$$
c$$$c     initialize return to success
c$$$      ier = 0
c$$$
c$$$c     loop over spatial indices
c$$$      do k=1,nzlocal,1
c$$$         do j=1,nylocal,1
c$$$            do i=1,nxlocal,1
c$$$
c$$$c              clear out ipiv
c$$$               ipiv = 0
c$$$
c$$$c              set up local block matrix, M = I+gamdt*Mcorr
c$$$               do m=1,nvar,1
c$$$                  do l=1,nvar,1
c$$$                     Mtmp(l,m) = gamdt*Mcorr(l,m,i,j,k)
c$$$                  enddo
c$$$               enddo
c$$$               do m=1,nvar,1
c$$$                  Mtmp(m,m) = Mtmp(m,m) + 1.d0
c$$$               enddo
c$$$
c$$$c              set up block rhs vector
c$$$               do l=1,nvar
c$$$                  wloc(l) = ww(i,j,k,l)
c$$$               enddo
c$$$
c$$$c              obtain LU factorization of M via LAPACK
c$$$               call dgetrf(nvar,nvar,Mtmp,nvar,ipiv,mier)
c$$$               if (mier < 0) then
c$$$                  write(0,*) 'dgetrf: illegal argument #',-mier
c$$$                  ier = -1
c$$$                  return
c$$$               elseif (mier > 0) then
c$$$                  write(0,*) 'dgetrf: singular row #',mier
c$$$                  ier = 1
c$$$                  return
c$$$               endif
c$$$
c$$$c              solve correction system via LAPACK
c$$$               call dgetrs('N',nvar,1,Mtmp,nvar,ipiv,wloc,nvar,mier)
c$$$               if (mier < 0) then
c$$$                  write(0,*) 'dgetrs: illegal argument #',-mier
c$$$                  ier = -1
c$$$                  return
c$$$               endif        
c$$$
c$$$c              put solution back into ww
c$$$               do l=1,nvar,1
c$$$                  ww(i,j,k,l) = wloc(l)
c$$$               enddo
c$$$
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$      return
c$$$         
c$$$      end subroutine fastwave_correction
c$$$c=======================================================================






      subroutine construct_FWcorrection0(ier)
c-----------------------------------------------------------------------
c     Description: construct_FWcorrection0 provides the non-constant 
c        coefficient preconditioner correction matrix (sans the 
c        I-gamma part), based on the derivation 0,
c            Mcorr = -R_i*Lambda_i*d[L_i]/dx_i  (sum over i)
c     
c     NOTE: All operations are performed on data within the 
c        fastwaveprec_data module
c
c     Arguments: ier - return flag (0=>success, 1=>failure)
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use fastwaveprec_data

c======= Declarations =========
      implicit none

c     calling variables
      integer :: ier

c     local variables 
      real*8, dimension(nvar,nvar) :: Mt0, Mt1, Mt2
      integer :: i, j, k, l, m

c======= Internals ============

c     initialize ier to success, since all computations either succeed 
c     or cause a segmentation fault (no errors to catch).
      ier = 0

c     loop over spatial indices
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1

c              set up x piece of block matrix, Mx = -R*Lambda*d[L]/dx
c                 Mt0 = R*Lambda
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = evr_x(l,m,i,j,k)*alamda_x(m,i,j,k)
                  enddo
               enddo

c                 Mt1 = R*Lambda*L_{i-1}
c                 Mt2 = R*Lambda*L_{i+1}
               Mt1(:,:) = matmul(Mt0(:,:),evl_x(:,:,i-1,j,k))
               Mt2(:,:) = matmul(Mt0(:,:),evl_x(:,:,i+1,j,k))

c                 Mc = -R*Lambda*d[L]/dx
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = dxi2*(Mt1(l,m) - Mt2(l,m))
                  enddo
               enddo

#ifndef ONE_D
c              add in y piece, My = -PermYt*R*Lambda*d[L]/dy*PermY

c                 Mt0 = R*Lambda
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = evr_y(l,m,j,i,k)*alamda_y(m,j,i,k)
                  enddo
               enddo      

c                 Mt1 = R*Lambda*L_{j-1}
c                 Mt2 = R*Lambda*L_{j+1}
               Mt1(:,:) = matmul(Mt0(:,:),evl_y(:,:,j-1,i,k))
               Mt2(:,:) = matmul(Mt0(:,:),evl_y(:,:,j+1,i,k))
               
c                 Mt0 = -R*Lambda*d[L]/dy
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = dyi2*(Mt1(l,m) - Mt2(l,m))
                  enddo
               enddo

c                 Mt1 = PermYt*Mt0
               Mt1 = matmul(PermYt, Mt0)

c                 Mt0 = Mt1*PermY
               Mt0 = matmul(Mt1, PermY)

c                 Mc = Mc + Mt0
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = Mcorr(l,m,i,j,k) + Mt0(l,m)
                  enddo
               enddo
               
#ifndef TWO_D

c              add in z piece, Mz = -PermZt*R*Lambda*d[L]/dz*PermZ

c                 Mt0 = R*Lambda
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = evr_z(l,m,k,i,j)*alamda_z(m,k,i,j)
                  enddo
               enddo
               
c                 Mt1 = R*Lambda*L_{k-1}
c                 Mt2 = R*Lambda*L_{k+1}
               Mt1(:,:) = matmul(Mt0(:,:),evl_z(:,:,k-1,i,j))
               Mt2(:,:) = matmul(Mt0(:,:),evl_z(:,:,k+1,i,j))
               
c                 Mt0 = -R*Lambda*d[L]/dz
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = dzi2*(Mt1(l,m) - Mt2(l,m))
                  enddo
               enddo

c                 Mt1 = PermZt*Mt0
               Mt1 = matmul(PermZt, Mt0)

c                 Mt0 = Mt1*PermZ
               Mt0 = matmul(Mt1, PermZ)
               
c                 Mc = Mc + Mt0
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = Mcorr(l,m,i,j,k) + Mt0(l,m)
                  enddo
               enddo

#endif
#endif
            enddo
         enddo
      enddo

      return
         
      end subroutine construct_FWcorrection0
c=======================================================================






      subroutine construct_FWcorrection1(ier)
c-----------------------------------------------------------------------
c     Description: construct_FWcorrection1 provides the non-constant 
c        coefficient preconditioner correction matrix (sans the 
c        I-gamma part), based on the derivation 1,
c            Mcorr = d[R_i*Lambda_i]/dx_i*L_i  (sum over i)
c     
c     NOTE: All operations are performed on data within the 
c        fastwaveprec_data module
c
c     Arguments: ier - return flag (0=>success, 1=>failure)
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use fastwaveprec_data

c======= Declarations =========
      implicit none

c     calling variables
      integer :: ier

c     local variables 
      real*8, dimension(nvar,nvar) :: Mt0, Mt1, Mt2
      integer :: i, j, k, l, m

c======= Internals ============

c     initialize ier to success, since all computations either succeed 
c     or cause a segmentation fault (no errors to catch).
      ier = 0

c     loop over spatial indices
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1

c              set up x piece of block matrix, Mx = d[R*Lambda]/dx*L

c                 Mt0 = (R*Lambda)_{i-1}
c                 Mt1 = (R*Lambda)_{i+1}
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = evr_x(l,m,i-1,j,k)
     &                    *alamda_x(m,i-1,j,k)
                     Mt1(l,m) = evr_x(l,m,i+1,j,k)
     &                    *alamda_x(m,i+1,j,k)
                  enddo
               enddo

c                 Mt2 = d[R*Lambda]/dx
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt2(l,m) = dxi2*(Mt1(l,m) - Mt0(l,m))
                  enddo
               enddo

c                 Mc =  d[R*Lambda]/dx*L_i
               Mcorr(:,:,i,j,k) = matmul(Mt2(:,:),evl_x(:,:,i,j,k))

#ifndef ONE_D
c              add in y piece, My = PermYt*d[R*Lambda]/dy*L*PermY
c                 Mt0 = (R*Lambda)_{j-1}
c                 Mt1 = (R*Lambda)_{j+1}
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = evr_y(l,m,j-1,i,k)
     &                    *alamda_y(m,j-1,i,k)
                     Mt1(l,m) = evr_y(l,m,j+1,i,k)
     &                    *alamda_y(m,j+1,i,k)
                  enddo
               enddo

c                 Mt2 = d[R*Lambda]/dy
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt2(l,m) = dyi2*(Mt1(l,m) - Mt0(l,m))
                  enddo
               enddo

c                 Mt0 =  d[R*Lambda]/dy*L
               Mt0(:,:) = matmul(Mt2(:,:),evl_y(:,:,j,i,k))

c                 Mt1 = PermYt*Mt0
               Mt1 = matmul(PermYt, Mt0)

c                 Mt0 = Mt1*PermY
               Mt0 = matmul(Mt1, PermY)

c                 Mc = Mc + Mt0
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = Mcorr(l,m,i,j,k) + Mt0(l,m)
                  enddo
               enddo
               
#ifndef TWO_D

c              add in z piece, Mz = PermZt*d[R*Lambda]/dz*L*PermZ
c                 Mt0 = (R*Lambda)_{k-1}
c                 Mt1 = (R*Lambda)_{k+1}
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = evr_z(l,m,k-1,i,j)
     &                    *alamda_z(m,k-1,i,j)
                     Mt1(l,m) = evr_z(l,m,k+1,i,j)
     &                    *alamda_z(m,k+1,i,j)
                  enddo
               enddo

c                 Mt2 = d[R*Lambda]/dz
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt2(l,m) = dzi2*(Mt1(l,m) - Mt0(l,m))
                  enddo
               enddo

c                 Mt0 =  d[R*Lambda]/dz*L
               Mt0(:,:) = matmul(Mt2(:,:),evl_z(:,:,k,i,j))

c                 Mt1 = PermZt*Mt0
               Mt1 = matmul(PermZt, Mt0)

c                 Mt0 = Mt1*PermZ
               Mt0 = matmul(Mt1, PermZ)

c                 Mc = Mc + Mt0
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = Mcorr(l,m,i,j,k) + Mt0(l,m)
                  enddo
               enddo

#endif
#endif
            enddo
         enddo
      enddo

      return
         
      end subroutine construct_FWcorrection1
c=======================================================================






      subroutine construct_FWcorrection2(ier)
c-----------------------------------------------------------------------
c     Description: construct_FWcorrection2 provides the non-constant 
c        coefficient preconditioner correction matrix (sans the 
c        I-gamma part), based on the derivation 2,
c           Mcorr = -R_i*d[L_i]/dx_i*R_i*Lambda_i*L_i  (sum over i)
c     
c     NOTE: All operations are performed on data within the 
c        fastwaveprec_data module
c
c     Arguments: ier - return flag (0=>success, 1=>failure)
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use fastwaveprec_data

c======= Declarations =========
      implicit none

c     calling variables
      integer :: ier

c     local variables 
      real*8, dimension(nvar,nvar) :: Mt0, Mt1, Mt2
      integer :: i, j, k, l, m

c======= Internals ============

c     initialize ier to success, since all computations either succeed 
c     or cause a segmentation fault (no errors to catch).
      ier = 0

c     loop over spatial indices
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1

c              set up x piece of block matrix, Mx = -R*d[L]/dx*R*Lambda*L
c                 Mt0 = -d[L]/dx
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = dxi2*(evl_x(l,m,i-1,j,k)
     &                               -evl_x(l,m,i+1,j,k))
                  enddo
               enddo

c                 Mt1 = -R*d[L]/dx*R
               Mt2 = matmul(Mt0,evr_x(:,:,i,j,k))
               Mt1 = matmul(evr_x(:,:,i,j,k),Mt2)

c                 Mt2 = Lambda*L
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt2(l,m) = alamda_x(l,i,j,k)*evl_x(l,m,i,j,k)
                  enddo
               enddo

c                 Mc = -R*d[L]/dx*R*Lambda*L
               Mcorr(:,:,i,j,k) = matmul(Mt1,Mt2)
                  
#ifndef ONE_D
c              add in y piece, My = -PermYt*R*d[L]/dy*R*Lambda*L*PermY
c                 Mt0 = -d[L]/dy
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = dyi2*(evl_y(l,m,j-1,i,k)
     &                               -evl_y(l,m,j+1,i,k))
                  enddo
               enddo

c                 Mt1 = -R*d[L]/dy*R
               Mt2 = matmul(Mt0,evr_y(:,:,j,i,k))
               Mt1 = matmul(evr_y(:,:,j,i,k),Mt2)

c                 Mt2 = Lambda*L
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt2(l,m) = alamda_y(l,j,i,k)*evl_y(l,m,j,i,k)
                  enddo
               enddo

c                 Mt0 = -R*d[L]/dy*R*Lambda*L
               Mt0 = matmul(Mt1,Mt2)

c                 Mt1 = PermYt*Mt0
               Mt1 = matmul(PermYt, Mt0)

c                 Mt0 = Mt1*PermY
               Mt0 = matmul(Mt1, PermY)

c                 Mc = Mc + Mt0
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = Mcorr(l,m,i,j,k) + Mt0(l,m)
                  enddo
               enddo
               
#ifndef TWO_D
c              add in z piece, Mz = -PermZt*R*d[L]/dz*R*Lambda*L*PermZ
c                 Mt0 = -d[L]/dz
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt0(l,m) = dzi2*(evl_z(l,m,k-1,i,j)
     &                               -evl_z(l,m,k+1,i,j))
                  enddo
               enddo

c                 Mt1 = -R*d[L]/dz*R
               Mt2 = matmul(Mt0,evr_z(:,:,k,i,j))
               Mt1 = matmul(evr_z(:,:,k,i,j),Mt2)

c                 Mt2 = Lambda*L
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mt2(l,m) = alamda_z(l,k,i,j)*evl_z(l,m,k,i,j)
                  enddo
               enddo

c                 Mt0 = -R*d[L]/dz*R*Lambda*L
               Mt0 = matmul(Mt1,Mt2)

c                 Mt1 = PermZt*Mt0
               Mt1 = matmul(PermZt, Mt0)

c                 Mt0 = Mt1*PermZ
               Mt0 = matmul(Mt1, PermZ)

c                 Mc = Mc + Mt0
               do m=1,nvar,1
                  do l=1,nvar,1
                     Mcorr(l,m,i,j,k) = Mcorr(l,m,i,j,k) + Mt0(l,m)
                  enddo
               enddo

#endif
#endif
            enddo
         enddo
      enddo

      return
         
      end subroutine construct_FWcorrection2
c=======================================================================






      subroutine SetADIEigenSystem(uu, vv)
c-----------------------------------------------------------------------
c     Description: SetADIEigenSystem computes the x, y, and z-
c        directional eigensystems over the entire domain, for use 
c        within the fast wave preconditioning solution approach.
c
c     Arguments:
c          uu - (dbl(*), input) conserved state variables
c          vv - (dbl(*)) scratch space with same size as uu
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use fastwaveprec_data
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      real*8, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &        :: uu, vv

c     local variables 
      integer :: i, j, k, l, m
      integer :: xlbdry, xhbdry, ylbdry, yhbdry, zlbdry, zhbdry
      real*8  :: dtmp, kenergy, gamm, gammInv, Benergy
      real*8  :: rhoInv, rhoSqrt, rhoSqrtInv, rhoCSqrt
      real*8  :: TwoSqrtInv, rhoCSqrtInv, Asq, csndsq
      real*8  :: csnd, csndinv
      real*8  :: csndsqInv, rho, p, betaPerpMagInv
      real*8  :: u_x, v_x, w_x, bx_x, by_x, bz_x
      real*8  :: u_y, v_y, w_y, bx_y, by_y, bz_y
      real*8  :: u_z, v_z, w_z, bx_z, by_z, bz_z
      real*8  :: axsq, alphaS, alphaF, Sbx, cfast, cslow
      real*8  :: betaY_z, betaY_y, betaY_x
      real*8  :: betaZ_x, betaZ_y, betaZ_z
      logical :: x_int, y_int, z_int

c======= Internals ============

c     first convert conserved variables to primitive variables
      call ConservedToPrimitive(uu, vv)

c     General constants
      gamm = gamma-1.D0
      gammInv = 1.d0/gamm
      TwoSqrtInv = dsqrt(0.5D0)

c     set directional extents for performing eigen-decomposition
c     (include one ghost zone for all nontrivial dimensions)
      ylbdry = 1; yhbdry = nylocal; zlbdry = 1; zhbdry = nzlocal
      xlbdry = 0
      xhbdry = nxlocal+1
#ifndef ONE_D
      ylbdry = 0
      yhbdry = nylocal+1
#ifndef TWO_D
      zlbdry = 0
      zhbdry = nzlocal+1
#endif
#endif      

c     loop over space
      do k=zlbdry,zhbdry,1

c        determine whether we are in the interior of z
         z_int = .false.
         if ((k>=1) .and. (k<=nzlocal)) z_int = .true.

         do j=ylbdry,yhbdry,1

c           determine whether we are in the interior of y
            y_int = .false.
            if ((j>=1) .and. (j<=nylocal)) y_int = .true.

            do i=xlbdry,xhbdry,1

c              determine whether we are in the interior of x
               x_int = .false.
               if ((i>=1) .and. (i<=nxlocal)) x_int = .true.
               
c              set internal primitive vars
               rho  = vv(i,j,k,1)
               p    = vv(i,j,k,8)
               u_x  = vv(i,j,k,2)
               v_x  = vv(i,j,k,3)
               w_x  = vv(i,j,k,4)
               bx_x = vv(i,j,k,5)
               by_x = vv(i,j,k,6)
               bz_x = vv(i,j,k,7)
#ifndef ONE_D
               u_y  = vv(i,j,k,3)
               v_y  = vv(i,j,k,4)
               w_y  = vv(i,j,k,2)
               bx_y = vv(i,j,k,6)
               by_y = vv(i,j,k,7)
               bz_y = vv(i,j,k,5)
#ifndef TWO_D
               u_z  = vv(i,j,k,4)
               v_z  = vv(i,j,k,2)
               w_z  = vv(i,j,k,3)
               bx_z = vv(i,j,k,7)
               by_z = vv(i,j,k,5)
               bz_z = vv(i,j,k,6)
#endif
#endif
               
c              Plasma beta
               Benergy = dsqrt(by_x*by_x + bz_x*bz_x)
               if (Benergy == 0.d0) then
                  betaY_x = TwoSqrtInv
                  betaZ_x = TwoSqrtInv
               else
                  betaPerpMagInv = 1.D0/Benergy
                  betaY_x = by_x*betaPerpMagInv
                  betaZ_x = bz_x*betaPerpMagInv
               endif
#ifndef ONE_D
               Benergy = dsqrt(by_y*by_y + bz_y*bz_y)
               if (Benergy == 0.d0) then
                  betaY_y = TwoSqrtInv
                  betaZ_y = TwoSqrtInv
               else
                  betaPerpMagInv = 1.D0/Benergy
                  betaY_y = by_y*betaPerpMagInv
                  betaZ_y = bz_y*betaPerpMagInv
               endif

#ifndef TWO_D
               Benergy = dsqrt(by_z*by_z + bz_z*bz_z)
               if (Benergy == 0.d0) then
                  betaY_z = TwoSqrtInv
                  betaZ_z = TwoSqrtInv
               else
                  betaPerpMagInv = 1.D0/Benergy
                  betaY_z = by_z*betaPerpMagInv
                  betaZ_z = bz_z*betaPerpMagInv
               endif
#endif
#endif
               
c              Set some constants
               rhoInv = 1.D0/rho
               csndsq = gamma*p*rhoInv
               Asq = csndsq + (bx_x**2 + by_x**2 + bz_x**2)*rhoInv
               csndsqInv = 1.D0/csndsq
               rhoSqrtInv = dsqrt(rhoInv)
               rhoSqrt = dsqrt(rho)
               rhoCSqrt = dsqrt(rho*csndsq)
               rhoCSqrtInv = 1.D0/rhoCSqrt
               kenergy = u_x**2 + v_x**2 + w_x**2



               
c              x-directional eigen-decomposition
               if (y_int .and. z_int) then

c                 x-direction specific constants
                  axsq = rhoInv*bx_x**2
                  dtmp = dsqrt(dmax1((Asq**2-4.D0*csndsq*axsq),0.D0))
                  cfast = 0.5D0*(Asq+dtmp)
                  cslow = 0.5D0*(Asq-dtmp)
                  dtmp = 1.d0/(cfast-cslow)
                  alphaS = dsqrt((dmax1(cfast-csndsq,0.D0))*dtmp)
                  alphaF = dsqrt((dmax1(csndsq-cslow,0.D0))*dtmp)
                  Sbx = dsign(1.D0,bx_x)
                  
c                 Rescale some constants for waves
                  cfast = dsqrt(cfast)
                  cslow = dsqrt(dmax1(cslow,0.D0))
                  axsq  = dsqrt(axsq)
                  
c                 x-directional Eigenvalues
                  alamda_x(1,i,j,k) = u_x
                  alamda_x(2,i,j,k) = u_x
                  alamda_x(3,i,j,k) = u_x + axsq
                  alamda_x(4,i,j,k) = u_x - axsq
                  alamda_x(5,i,j,k) = u_x + cfast
                  alamda_x(6,i,j,k) = u_x - cfast
                  alamda_x(7,i,j,k) = u_x + cslow
                  alamda_x(8,i,j,k) = u_x - cslow
     
c                 Right eigenvectors
c                 (convert from characteristic to primitive variables)
                  evr_x(:,:,i,j,k) = 0.D0

c                 Right eigenvector corresponding to u
                  evr_x(1,1,i,j,k) = 1.D0

c                 Right eigenvector corresponding to u
                  evr_x(5,2,i,j,k) = 1.D0

c                 Right eigenvector corresponding to u+ax
                  evr_x(3,3,i,j,k) = -betaZ_x*TwoSqrtInv
                  evr_x(4,3,i,j,k) =  betaY_x*TwoSqrtInv
                  evr_x(6,3,i,j,k) =  rhoSqrt*TwoSqrtInv*betaZ_x*Sbx
                  evr_x(7,3,i,j,k) = -rhoSqrt*TwoSqrtInv*betaY_x*Sbx

c                 Right eigenvector corresponding to u-ax
                  evr_x(3,4,i,j,k) =  evr_x(3,3,i,j,k)
                  evr_x(4,4,i,j,k) =  evr_x(4,3,i,j,k)
                  evr_x(6,4,i,j,k) = -evr_x(6,3,i,j,k)
                  evr_x(7,4,i,j,k) = -evr_x(7,3,i,j,k)

c                 Right eigenvector corresponding to u+cfast
                  evr_x(1,5,i,j,k) =  rho*alphaF
                  evr_x(2,5,i,j,k) =  alphaF*cfast
                  evr_x(3,5,i,j,k) = -alphaS*cslow*betaY_x*Sbx
                  evr_x(4,5,i,j,k) = -alphaS*cslow*betaZ_x*Sbx
                  evr_x(6,5,i,j,k) =  alphaS*rhoCSqrt*betaY_x
                  evr_x(7,5,i,j,k) =  alphaS*rhoCSqrt*betaZ_x
                  evr_x(8,5,i,j,k) =  rho*alphaF*csndsq
                  
c                 Right eigenvector corresponding to u-cfast
                  evr_x(1,6,i,j,k) =  evr_x(1,5,i,j,k)
                  evr_x(2,6,i,j,k) = -evr_x(2,5,i,j,k)
                  evr_x(3,6,i,j,k) = -evr_x(3,5,i,j,k)
                  evr_x(4,6,i,j,k) = -evr_x(4,5,i,j,k)
                  evr_x(6,6,i,j,k) =  evr_x(6,5,i,j,k)
                  evr_x(7,6,i,j,k) =  evr_x(7,5,i,j,k)
                  evr_x(8,6,i,j,k) =  evr_x(8,5,i,j,k)
                  
c                 Right eigenvector corresponding to u+cslow
                  evr_x(1,7,i,j,k) =  rho*alphaS
                  evr_x(2,7,i,j,k) =  alphaS*cslow
                  evr_x(3,7,i,j,k) =  alphaF*cfast*betaY_x*Sbx
                  evr_x(4,7,i,j,k) =  alphaF*cfast*betaZ_x*Sbx
                  evr_x(6,7,i,j,k) = -alphaF*rhoCSqrt*betaY_x
                  evr_x(7,7,i,j,k) = -alphaF*rhoCSqrt*betaZ_x
                  evr_x(8,7,i,j,k) =  rho*alphaS*csndsq
                  
c                 Right eigenvector corresponding to u-cslow
                  evr_x(1,8,i,j,k) =  evr_x(1,7,i,j,k)
                  evr_x(2,8,i,j,k) = -evr_x(2,7,i,j,k)
                  evr_x(3,8,i,j,k) = -evr_x(3,7,i,j,k)
                  evr_x(4,8,i,j,k) = -evr_x(4,7,i,j,k)
                  evr_x(6,8,i,j,k) =  evr_x(6,7,i,j,k)
                  evr_x(7,8,i,j,k) =  evr_x(7,7,i,j,k)
                  evr_x(8,8,i,j,k) =  evr_x(8,7,i,j,k)

c                 Convert primitive to conservative eigenvectors
                  do l=1,nvar,1
                     evr_x(8,l,i,j,k) = 0.5D0*kenergy*evr_x(1,l,i,j,k)
     &                    + rho*(u_x*evr_x(2,l,i,j,k)
     &                         + v_x*evr_x(3,l,i,j,k)
     &                         + w_x*evr_x(4,l,i,j,k))
     &                    + bx_x*evr_x(5,l,i,j,k)
     &                    + by_x*evr_x(6,l,i,j,k)
     &                    + bz_x*evr_x(7,l,i,j,k)
     &                    + evr_x(8,l,i,j,k)*gammInv
                     evr_x(2,l,i,j,k) = u_x*evr_x(1,l,i,j,k) 
     &                                + rho*evr_x(2,l,i,j,k)
                     evr_x(3,l,i,j,k) = v_x*evr_x(1,l,i,j,k) 
     &                                + rho*evr_x(3,l,i,j,k)
                     evr_x(4,l,i,j,k) = w_x*evr_x(1,l,i,j,k) 
     &                                + rho*evr_x(4,l,i,j,k)
                  enddo
                  
c                 Left eigenvectors
                  evl_x(:,:,i,j,k) = 0.D0

c                 Left eigenvector corresponding to u
                  evl_x(1,1,i,j,k) =  1.D0
                  evl_x(1,8,i,j,k) = -csndsqInv
                  
c                 Left eigenvector corresponding to u
                  evl_x(2,5,i,j,k) = 1.D0
                  
c                 Left eigenvector corresponding to u+ax
                  evl_x(3,3,i,j,k) =-betaZ_x*TwoSqrtInv
                  evl_x(3,4,i,j,k) = betaY_x*TwoSqrtInv
                  evl_x(3,6,i,j,k) = rhoSqrtInv*TwoSqrtInv*betaZ_x*Sbx
                  evl_x(3,7,i,j,k) =-rhoSqrtInv*TwoSqrtInv*betaY_x*Sbx
                  
c                 Left eigenvector corresponding to u-ax
                  evl_x(4,3,i,j,k) =  evl_x(3,3,i,j,k)
                  evl_x(4,4,i,j,k) =  evl_x(3,4,i,j,k)
                  evl_x(4,6,i,j,k) = -evl_x(3,6,i,j,k)
                  evl_x(4,7,i,j,k) = -evl_x(3,7,i,j,k)
                  
c                 Left eigenvector corresponding to u+cfast
                  evl_x(5,2,i,j,k) =  0.5D0*alphaF*cfast*csndsqInv
                  evl_x(5,3,i,j,k) = -0.5D0*alphaS*cslow*betaY_x 
     &                 * Sbx*csndsqInv
                  evl_x(5,4,i,j,k) = -0.5D0*alphaS*cslow*betaZ_x 
     &                 * Sbx*csndsqInv
                  evl_x(5,6,i,j,k) =  0.5D0*alphaS*rhoCSqrtInv*betaY_x
                  evl_x(5,7,i,j,k) =  0.5D0*alphaS*rhoCSqrtInv*betaZ_x
                  evl_x(5,8,i,j,k) =  0.5D0*alphaF*rhoInv*csndsqInv
                  
c                 Left eigenvector corresponding to u-cfast
                  evl_x(6,2,i,j,k) = -evl_x(5,2,i,j,k)
                  evl_x(6,3,i,j,k) = -evl_x(5,3,i,j,k)
                  evl_x(6,4,i,j,k) = -evl_x(5,4,i,j,k)
                  evl_x(6,6,i,j,k) =  evl_x(5,6,i,j,k)
                  evl_x(6,7,i,j,k) =  evl_x(5,7,i,j,k)
                  evl_x(6,8,i,j,k) =  evl_x(5,8,i,j,k)
                  
c                 Left eigenvector corresponding to u+cslow 
                  evl_x(7,2,i,j,k) =  0.5D0*alphaS*cslow*csndsqInv
                  evl_x(7,3,i,j,k) =  0.5D0*alphaF*cfast*betaY_x 
     &                 * Sbx*csndsqInv
                  evl_x(7,4,i,j,k) =  0.5D0*alphaF*cfast*betaZ_x 
     &                 * Sbx*csndsqInv
                  evl_x(7,6,i,j,k) = -0.5D0*alphaF*rhoCSqrtInv*betaY_x
                  evl_x(7,7,i,j,k) = -0.5D0*alphaF*rhoCSqrtInv*betaZ_x
                  evl_x(7,8,i,j,k) =  0.5D0*alphaS*rhoInv*csndsqInv
                  
c                 Left eigenvector corresponding to u-cslow
                  evl_x(8,2,i,j,k) = -evl_x(7,2,i,j,k)
                  evl_x(8,3,i,j,k) = -evl_x(7,3,i,j,k)
                  evl_x(8,4,i,j,k) = -evl_x(7,4,i,j,k)
                  evl_x(8,6,i,j,k) =  evl_x(7,6,i,j,k)
                  evl_x(8,7,i,j,k) =  evl_x(7,7,i,j,k)
                  evl_x(8,8,i,j,k) =  evl_x(7,8,i,j,k)
                  
c                 Convert primitive to conservative eigenvectors
                  do l=1,nvar,1
                     evl_x(l,1,i,j,k) = evl_x(l,1,i,j,k) 
     &                    - rhoInv*(u_x*evl_x(l,2,i,j,k)
     &                            + v_x*evl_x(l,3,i,j,k)
     &                            + w_x*evl_x(l,4,i,j,k))
     &                    +0.5D0*gamm*kenergy*evl_x(l,8,i,j,k)
                     evl_x(l,2,i,j,k) = evl_x(l,2,i,j,k)*rhoInv
     &                       - gamm*u_x*evl_x(l,8,i,j,k)
                     evl_x(l,3,i,j,k) = evl_x(l,3,i,j,k)*rhoInv
     &                       - gamm*v_x*evl_x(l,8,i,j,k)
                     evl_x(l,4,i,j,k) = evl_x(l,4,i,j,k)*rhoInv
     &                       - gamm*w_x*evl_x(l,8,i,j,k)
                     evl_x(l,5,i,j,k) = evl_x(l,5,i,j,k) 
     &                      - gamm*bx_x*evl_x(l,8,i,j,k)
                     evl_x(l,6,i,j,k) = evl_x(l,6,i,j,k) 
     &                      - gamm*by_x*evl_x(l,8,i,j,k)
                     evl_x(l,7,i,j,k) = evl_x(l,7,i,j,k) 
     &                      - gamm*bz_x*evl_x(l,8,i,j,k)
                     evl_x(l,8,i,j,k) = gamm*evl_x(l,8,i,j,k)
                  enddo
               endif



               
c              y-directional eigen-decomposition
#ifndef ONE_D
               if (x_int .and. z_int) then

c                 y-direction specific constants
                  axsq = rhoInv*bx_y**2
                  dtmp = dsqrt(dmax1((Asq**2-4.D0*csndsq*axsq),0.D0))
                  cfast = 0.5D0*(Asq+dtmp)
                  cslow = 0.5D0*(Asq-dtmp)
                  dtmp = 1.d0/(cfast-cslow)
                  alphaS = dsqrt((dmax1(cfast-csndsq,0.D0))*dtmp)
                  alphaF = dsqrt((dmax1(csndsq-cslow,0.D0))*dtmp)
                  Sbx = dsign(1.D0,bx_y)

c                 Rescale some constants for waves
                  cfast = dsqrt(cfast)
                  cslow = dsqrt(dmax1(cslow,0.D0))
                  axsq  = dsqrt(axsq)
                  
c                 y-directional Eigenvalues
                  alamda_y(1,j,i,k) = u_y
                  alamda_y(2,j,i,k) = u_y
                  alamda_y(3,j,i,k) = u_y + axsq
                  alamda_y(4,j,i,k) = u_y - axsq
                  alamda_y(5,j,i,k) = u_y + cfast
                  alamda_y(6,j,i,k) = u_y - cfast
                  alamda_y(7,j,i,k) = u_y + cslow
                  alamda_y(8,j,i,k) = u_y - cslow
                  
c                 Right eigenvectors
c                 (convert from characteristic to primitive variables)
                  evr_y(:,:,j,i,k) = 0.D0

c                 Right eigenvector corresponding to u
                  evr_y(1,1,j,i,k) = 1.D0

c                 Right eigenvector corresponding to u
                  evr_y(5,2,j,i,k) = 1.D0

c                 Right eigenvector corresponding to u+ax
                  evr_y(3,3,j,i,k) = -betaZ_y*TwoSqrtInv
                  evr_y(4,3,j,i,k) =  betaY_y*TwoSqrtInv
                  evr_y(6,3,j,i,k) =  rhoSqrt*TwoSqrtInv*betaZ_y*Sbx
                  evr_y(7,3,j,i,k) = -rhoSqrt*TwoSqrtInv*betaY_y*Sbx
                  
c                 Right eigenvector corresponding to u-ax
                  evr_y(3,4,j,i,k) =  evr_y(3,3,j,i,k)
                  evr_y(4,4,j,i,k) =  evr_y(4,3,j,i,k)
                  evr_y(6,4,j,i,k) = -evr_y(6,3,j,i,k)
                  evr_y(7,4,j,i,k) = -evr_y(7,3,j,i,k)
                  
c                 Right eigenvector corresponding to u+cfast
                  evr_y(1,5,j,i,k) =  rho*alphaF
                  evr_y(2,5,j,i,k) =  alphaF*cfast
                  evr_y(3,5,j,i,k) = -alphaS*cslow*betaY_y*Sbx
                  evr_y(4,5,j,i,k) = -alphaS*cslow*betaZ_y*Sbx
                  evr_y(6,5,j,i,k) =  alphaS*rhoCSqrt*betaY_y
                  evr_y(7,5,j,i,k) =  alphaS*rhoCSqrt*betaZ_y
                  evr_y(8,5,j,i,k) =  rho*alphaF*csndsq
                  
c                 Right eigenvector corresponding to u-cfast
                  evr_y(1,6,j,i,k) =  evr_y(1,5,j,i,k)
                  evr_y(2,6,j,i,k) = -evr_y(2,5,j,i,k)
                  evr_y(3,6,j,i,k) = -evr_y(3,5,j,i,k)
                  evr_y(4,6,j,i,k) = -evr_y(4,5,j,i,k)
                  evr_y(6,6,j,i,k) =  evr_y(6,5,j,i,k)
                  evr_y(7,6,j,i,k) =  evr_y(7,5,j,i,k)
                  evr_y(8,6,j,i,k) =  evr_y(8,5,j,i,k)
                  
c                 Right eigenvector corresponding to u+cslow
                  evr_y(1,7,j,i,k) =  rho*alphaS
                  evr_y(2,7,j,i,k) =  alphaS*cslow
                  evr_y(3,7,j,i,k) =  alphaF*cfast*betaY_y*Sbx
                  evr_y(4,7,j,i,k) =  alphaF*cfast*betaZ_y*Sbx
                  evr_y(6,7,j,i,k) = -alphaF*rhoCSqrt*betaY_y
                  evr_y(7,7,j,i,k) = -alphaF*rhoCSqrt*betaZ_y
                  evr_y(8,7,j,i,k) =  rho*alphaS*csndsq
                  
c                 Right eigenvector corresponding to u-cslow
                  evr_y(1,8,j,i,k) =  evr_y(1,7,j,i,k)
                  evr_y(2,8,j,i,k) = -evr_y(2,7,j,i,k)
                  evr_y(3,8,j,i,k) = -evr_y(3,7,j,i,k)
                  evr_y(4,8,j,i,k) = -evr_y(4,7,j,i,k)
                  evr_y(6,8,j,i,k) =  evr_y(6,7,j,i,k)
                  evr_y(7,8,j,i,k) =  evr_y(7,7,j,i,k)
                  evr_y(8,8,j,i,k) =  evr_y(8,7,j,i,k)
                  
c                 Convert primitive to conservative eigenvectors
                  do l=1,nvar,1
                     evr_y(8,l,j,i,k) = 0.5D0*kenergy*evr_y(1,l,j,i,k)
     &                    + rho*(u_y*evr_y(2,l,j,i,k)
     &                         + v_y*evr_y(3,l,j,i,k)
     &                         + w_y*evr_y(4,l,j,i,k))
     &                    + bx_y*evr_y(5,l,j,i,k)
     &                    + by_y*evr_y(6,l,j,i,k)
     &                    + bz_y*evr_y(7,l,j,i,k)
     &                    + evr_y(8,l,j,i,k)*gammInv
                     evr_y(2,l,j,i,k) = u_y*evr_y(1,l,j,i,k) 
     &                                + rho*evr_y(2,l,j,i,k)
                     evr_y(3,l,j,i,k) = v_y*evr_y(1,l,j,i,k) 
     &                                + rho*evr_y(3,l,j,i,k)
                     evr_y(4,l,j,i,k) = w_y*evr_y(1,l,j,i,k) 
     &                                + rho*evr_y(4,l,j,i,k)
                  enddo
                  
c                 Left eigenvectors
                  evl_y(:,:,j,i,k) = 0.D0
                  
c                 Left eigenvector corresponding to u
                  evl_y(1,1,j,i,k) =  1.D0
                  evl_y(1,8,j,i,k) = -csndsqInv
                  
c                 Left eigenvector corresponding to u
                  evl_y(2,5,j,i,k) = 1.D0
                  
c                 Left eigenvector corresponding to u+ax
                  evl_y(3,3,j,i,k) =-betaZ_y*TwoSqrtInv
                  evl_y(3,4,j,i,k) = betaY_y*TwoSqrtInv
                  evl_y(3,6,j,i,k) = rhoSqrtInv*TwoSqrtInv*betaZ_y*Sbx
                  evl_y(3,7,j,i,k) =-rhoSqrtInv*TwoSqrtInv*betaY_y*Sbx
                  
c                 Left eigenvector corresponding to u-ax
                  evl_y(4,3,j,i,k) =  evl_y(3,3,j,i,k)
                  evl_y(4,4,j,i,k) =  evl_y(3,4,j,i,k)
                  evl_y(4,6,j,i,k) = -evl_y(3,6,j,i,k)
                  evl_y(4,7,j,i,k) = -evl_y(3,7,j,i,k)
                  
c                 Left eigenvector corresponding to u+cfast
                  evl_y(5,2,j,i,k) =  0.5D0*alphaF*cfast*csndsqInv
                  evl_y(5,3,j,i,k) = -0.5D0*alphaS*cslow*betaY_y 
     &                 * Sbx*csndsqInv
                  evl_y(5,4,j,i,k) = -0.5D0*alphaS*cslow*betaZ_y 
     &                 * Sbx*csndsqInv
                  evl_y(5,6,j,i,k) =  0.5D0*alphaS*rhoCSqrtInv*betaY_y
                  evl_y(5,7,j,i,k) =  0.5D0*alphaS*rhoCSqrtInv*betaZ_y
                  evl_y(5,8,j,i,k) =  0.5D0*alphaF*rhoInv*csndsqInv
                  
c                 Left eigenvector corresponding to u-cfast
                  evl_y(6,2,j,i,k) = -evl_y(5,2,j,i,k)
                  evl_y(6,3,j,i,k) = -evl_y(5,3,j,i,k)
                  evl_y(6,4,j,i,k) = -evl_y(5,4,j,i,k)
                  evl_y(6,6,j,i,k) =  evl_y(5,6,j,i,k)
                  evl_y(6,7,j,i,k) =  evl_y(5,7,j,i,k)
                  evl_y(6,8,j,i,k) =  evl_y(5,8,j,i,k)
                  
c                 Left eigenvector corresponding to u+cslow 
                  evl_y(7,2,j,i,k) =  0.5D0*alphaS*cslow*csndsqInv
                  evl_y(7,3,j,i,k) =  0.5D0*alphaF*cfast*betaY_y 
     &                 * Sbx*csndsqInv
                  evl_y(7,4,j,i,k) =  0.5D0*alphaF*cfast*betaZ_y 
     &                 * Sbx*csndsqInv
                  evl_y(7,6,j,i,k) = -0.5D0*alphaF*rhoCSqrtInv*betaY_y
                  evl_y(7,7,j,i,k) = -0.5D0*alphaF*rhoCSqrtInv*betaZ_y
                  evl_y(7,8,j,i,k) =  0.5D0*alphaS*rhoInv*csndsqInv
                  
c                 Left eigenvector corresponding to u-cslow
                  evl_y(8,2,j,i,k) = -evl_y(7,2,j,i,k)
                  evl_y(8,3,j,i,k) = -evl_y(7,3,j,i,k)
                  evl_y(8,4,j,i,k) = -evl_y(7,4,j,i,k)
                  evl_y(8,6,j,i,k) =  evl_y(7,6,j,i,k)
                  evl_y(8,7,j,i,k) =  evl_y(7,7,j,i,k)
                  evl_y(8,8,j,i,k) =  evl_y(7,8,j,i,k)
                  
c                 Convert primitive to conservative eigenvectors
                  do l=1,nvar,1
                     evl_y(l,1,j,i,k) = evl_y(l,1,j,i,k) 
     &                    - rhoInv*(u_y*evl_y(l,2,j,i,k)
     &                            + v_y*evl_y(l,3,j,i,k)
     &                            + w_y*evl_y(l,4,j,i,k))
     &                    + 0.5D0*gamm*kenergy*evl_y(l,8,j,i,k)
                     evl_y(l,2,j,i,k) = evl_y(l,2,j,i,k)*rhoInv
     &                       - gamm*u_y*evl_y(l,8,j,i,k)
                     evl_y(l,3,j,i,k) = evl_y(l,3,j,i,k)*rhoInv
     &                       - gamm*v_y*evl_y(l,8,j,i,k)
                     evl_y(l,4,j,i,k) = evl_y(l,4,j,i,k)*rhoInv
     &                       - gamm*w_y*evl_y(l,8,j,i,k)
                     evl_y(l,5,j,i,k) = evl_y(l,5,j,i,k) 
     &                      - gamm*bx_y*evl_y(l,8,j,i,k)
                     evl_y(l,6,j,i,k) = evl_y(l,6,j,i,k) 
     &                      - gamm*by_y*evl_y(l,8,j,i,k)
                     evl_y(l,7,j,i,k) = evl_y(l,7,j,i,k) 
     &                      - gamm*bz_y*evl_y(l,8,j,i,k)
                     evl_y(l,8,j,i,k) = gamm*evl_y(l,8,j,i,k)
                  enddo
               endif



               
c              z-directional eigen-decomposition
#ifndef TWO_D
               if (x_int .and. y_int) then

c                 z-direction specific constants
                  axsq = rhoInv*bx_z**2
                  dtmp = dsqrt(dmax1((Asq**2-4.D0*csndsq*axsq),0.D0))
                  cfast = 0.5D0*(Asq+dtmp)
                  cslow = 0.5D0*(Asq-dtmp)
                  dtmp = 1.d0/(cfast-cslow)
                  alphaS = dsqrt((dmax1(cfast-csndsq,0.D0))*dtmp)
                  alphaF = dsqrt((dmax1(csndsq-cslow,0.D0))*dtmp)
                  Sbx = dsign(1.D0,bx_z)

c                 Rescale some constants for waves
                  cfast = dsqrt(cfast)
                  cslow = dsqrt(dmax1(cslow,0.D0))
                  axsq  = dsqrt(axsq)
                  
c                 z-directional Eigenvalues
                  alamda_z(1,k,i,j) = u_z
                  alamda_z(2,k,i,j) = u_z
                  alamda_z(3,k,i,j) = u_z + axsq
                  alamda_z(4,k,i,j) = u_z - axsq
                  alamda_z(5,k,i,j) = u_z + cfast
                  alamda_z(6,k,i,j) = u_z - cfast
                  alamda_z(7,k,i,j) = u_z + cslow
                  alamda_z(8,k,i,j) = u_z - cslow
                  
c                 Right eigenvectors
c                 (convert from characteristic to primitive variables)
                  evr_z(:,:,k,i,j) = 0.D0

c                 Right eigenvector corresponding to u
                  evr_z(1,1,k,i,j) = 1.D0

c                 Right eigenvector corresponding to u
                  evr_z(5,2,k,i,j) = 1.D0

c                 Right eigenvector corresponding to u+ax
                  evr_z(3,3,k,i,j) = -betaZ_z*TwoSqrtInv
                  evr_z(4,3,k,i,j) =  betaY_z*TwoSqrtInv
                  evr_z(6,3,k,i,j) =  rhoSqrt*TwoSqrtInv*betaZ_z*Sbx
                  evr_z(7,3,k,i,j) = -rhoSqrt*TwoSqrtInv*betaY_z*Sbx
                  
c                 Right eigenvector corresponding to u-ax
                  evr_z(3,4,k,i,j) =  evr_z(3,3,k,i,j)
                  evr_z(4,4,k,i,j) =  evr_z(4,3,k,i,j)
                  evr_z(6,4,k,i,j) = -evr_z(6,3,k,i,j)
                  evr_z(7,4,k,i,j) = -evr_z(7,3,k,i,j)
                  
c                 Right eigenvector corresponding to u+cfast
                  evr_z(1,5,k,i,j) =  rho*alphaF
                  evr_z(2,5,k,i,j) =  alphaF*cfast
                  evr_z(3,5,k,i,j) = -alphaS*cslow*betaY_z*Sbx
                  evr_z(4,5,k,i,j) = -alphaS*cslow*betaZ_z*Sbx
                  evr_z(6,5,k,i,j) =  alphaS*rhoCSqrt*betaY_z
                  evr_z(7,5,k,i,j) =  alphaS*rhoCSqrt*betaZ_z
                  evr_z(8,5,k,i,j) =  rho*alphaF*csndsq
                  
c                 Right eigenvector corresponding to u-cfast
                  evr_z(1,6,k,i,j) =  evr_z(1,5,k,i,j)
                  evr_z(2,6,k,i,j) = -evr_z(2,5,k,i,j)
                  evr_z(3,6,k,i,j) = -evr_z(3,5,k,i,j)
                  evr_z(4,6,k,i,j) = -evr_z(4,5,k,i,j)
                  evr_z(6,6,k,i,j) =  evr_z(6,5,k,i,j)
                  evr_z(7,6,k,i,j) =  evr_z(7,5,k,i,j)
                  evr_z(8,6,k,i,j) =  evr_z(8,5,k,i,j)
                  
c                 Right eigenvector corresponding to u+cslow
                  evr_z(1,7,k,i,j) =  rho*alphaS
                  evr_z(2,7,k,i,j) =  alphaS*cslow
                  evr_z(3,7,k,i,j) =  alphaF*cfast*betaY_z*Sbx
                  evr_z(4,7,k,i,j) =  alphaF*cfast*betaZ_z*Sbx
                  evr_z(6,7,k,i,j) = -alphaF*rhoCSqrt*betaY_z
                  evr_z(7,7,k,i,j) = -alphaF*rhoCSqrt*betaZ_z
                  evr_z(8,7,k,i,j) =  rho*alphaS*csndsq
                  
c                 Right eigenvector corresponding to u-cslow
                  evr_z(1,8,k,i,j) =  evr_z(1,7,k,i,j)
                  evr_z(2,8,k,i,j) = -evr_z(2,7,k,i,j)
                  evr_z(3,8,k,i,j) = -evr_z(3,7,k,i,j)
                  evr_z(4,8,k,i,j) = -evr_z(4,7,k,i,j)
                  evr_z(6,8,k,i,j) =  evr_z(6,7,k,i,j)
                  evr_z(7,8,k,i,j) =  evr_z(7,7,k,i,j)
                  evr_z(8,8,k,i,j) =  evr_z(8,7,k,i,j)
                  
c                 Convert primitive to conservative eigenvectors
                  do l=1,nvar,1
                     evr_z(8,l,k,i,j) = 0.5D0*kenergy*evr_z(1,l,k,i,j)
     &                    + rho*(u_z*evr_z(2,l,k,i,j)
     &                         + v_z*evr_z(3,l,k,i,j)
     &                         + w_z*evr_z(4,l,k,i,j))
     &                    + bx_z*evr_z(5,l,k,i,j)
     &                    + by_z*evr_z(6,l,k,i,j)
     &                    + bz_z*evr_z(7,l,k,i,j)
     &                    + evr_z(8,l,k,i,j)*gammInv
                     evr_z(2,l,k,i,j) = u_z*evr_z(1,l,k,i,j) 
     &                                + rho*evr_z(2,l,k,i,j)
                     evr_z(3,l,k,i,j) = v_z*evr_z(1,l,k,i,j) 
     &                                + rho*evr_z(3,l,k,i,j)
                     evr_z(4,l,k,i,j) = w_z*evr_z(1,l,k,i,j) 
     &                                + rho*evr_z(4,l,k,i,j)
                  enddo
                  
c                 Left eigenvectors
                  evl_z(:,:,k,i,j) = 0.D0
                  
c                 Left eigenvector corresponding to u
                  evl_z(1,1,k,i,j) =  1.D0
                  evl_z(1,8,k,i,j) = -csndsqInv
                  
c                 Left eigenvector corresponding to u
                  evl_z(2,5,k,i,j) = 1.D0
                  
c                 Left eigenvector corresponding to u+ax
                  evl_z(3,3,k,i,j) =-betaZ_z*TwoSqrtInv
                  evl_z(3,4,k,i,j) = betaY_z*TwoSqrtInv
                  evl_z(3,6,k,i,j) = rhoSqrtInv*TwoSqrtInv*betaZ_z*Sbx
                  evl_z(3,7,k,i,j) =-rhoSqrtInv*TwoSqrtInv*betaY_z*Sbx
                  
c                 Left eigenvector corresponding to u-ax
                  evl_z(4,3,k,i,j) =  evl_z(3,3,k,i,j)
                  evl_z(4,4,k,i,j) =  evl_z(3,4,k,i,j)
                  evl_z(4,6,k,i,j) = -evl_z(3,6,k,i,j)
                  evl_z(4,7,k,i,j) = -evl_z(3,7,k,i,j)
                  
c                 Left eigenvector corresponding to u+cfast
                  evl_z(5,2,k,i,j) =  0.5D0*alphaF*cfast*csndsqInv
                  evl_z(5,3,k,i,j) = -0.5D0*alphaS*cslow*betaY_z 
     &                 * Sbx*csndsqInv
                  evl_z(5,4,k,i,j) = -0.5D0*alphaS*cslow*betaZ_z 
     &                 * Sbx*csndsqInv
                  evl_z(5,6,k,i,j) =  0.5D0*alphaS*rhoCSqrtInv*betaY_z
                  evl_z(5,7,k,i,j) =  0.5D0*alphaS*rhoCSqrtInv*betaZ_z
                  evl_z(5,8,k,i,j) =  0.5D0*alphaF*rhoInv*csndsqInv
                  
c                 Left eigenvector corresponding to u-cfast
                  evl_z(6,2,k,i,j) = -evl_z(5,2,k,i,j)
                  evl_z(6,3,k,i,j) = -evl_z(5,3,k,i,j)
                  evl_z(6,4,k,i,j) = -evl_z(5,4,k,i,j)
                  evl_z(6,6,k,i,j) =  evl_z(5,6,k,i,j)
                  evl_z(6,7,k,i,j) =  evl_z(5,7,k,i,j)
                  evl_z(6,8,k,i,j) =  evl_z(5,8,k,i,j)
                  
c                 Left eigenvector corresponding to u+cslow 
                  evl_z(7,2,k,i,j) =  0.5D0*alphaS*cslow*csndsqInv
                  evl_z(7,3,k,i,j) =  0.5D0*alphaF*cfast*betaY_z 
     &                 * Sbx*csndsqInv
                  evl_z(7,4,k,i,j) =  0.5D0*alphaF*cfast*betaZ_z 
     &                 * Sbx*csndsqInv
                  evl_z(7,6,k,i,j) = -0.5D0*alphaF*rhoCSqrtInv*betaY_z
                  evl_z(7,7,k,i,j) = -0.5D0*alphaF*rhoCSqrtInv*betaZ_z
                  evl_z(7,8,k,i,j) =  0.5D0*alphaS*rhoInv*csndsqInv
                  
c                 Left eigenvector corresponding to u-cslow
                  evl_z(8,2,k,i,j) = -evl_z(7,2,k,i,j)
                  evl_z(8,3,k,i,j) = -evl_z(7,3,k,i,j)
                  evl_z(8,4,k,i,j) = -evl_z(7,4,k,i,j)
                  evl_z(8,6,k,i,j) =  evl_z(7,6,k,i,j)
                  evl_z(8,7,k,i,j) =  evl_z(7,7,k,i,j)
                  evl_z(8,8,k,i,j) =  evl_z(7,8,k,i,j)
                  
c                 Convert primitive to conservative eigenvectors
                  do l=1,nvar,1
                     evl_z(l,1,k,i,j) = evl_z(l,1,k,i,j) 
     &                    - rhoInv*(u_z*evl_z(l,2,k,i,j)
     &                            + v_z*evl_z(l,3,k,i,j)
     &                            + w_z*evl_z(l,4,k,i,j))
     &                    + 0.5D0*gamm*kenergy*evl_z(l,8,k,i,j)
                     evl_z(l,2,k,i,j) = evl_z(l,2,k,i,j)*rhoInv
     &                       - gamm*u_z*evl_z(l,8,k,i,j)
                     evl_z(l,3,k,i,j) = evl_z(l,3,k,i,j)*rhoInv
     &                       - gamm*v_z*evl_z(l,8,k,i,j)
                     evl_z(l,4,k,i,j) = evl_z(l,4,k,i,j)*rhoInv
     &                       - gamm*w_z*evl_z(l,8,k,i,j)
                     evl_z(l,5,k,i,j) = evl_z(l,5,k,i,j) 
     &                      - gamm*bx_z*evl_z(l,8,k,i,j)
                     evl_z(l,6,k,i,j) = evl_z(l,6,k,i,j) 
     &                      - gamm*by_z*evl_z(l,8,k,i,j)
                     evl_z(l,7,k,i,j) = evl_z(l,7,k,i,j) 
     &                      - gamm*bz_z*evl_z(l,8,k,i,j)
                     evl_z(l,8,k,i,j) = gamm*evl_z(l,8,k,i,j)
                  enddo
               endif
#endif
#endif
            enddo
         enddo
      enddo

      return

      end subroutine SetADIEigenSystem
c=======================================================================
