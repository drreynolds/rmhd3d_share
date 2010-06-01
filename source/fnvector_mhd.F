c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------------
c     $Log: fnvector_mhd.F,v $
c=======================================================================



      subroutine fnvmhdprint(xvec)
c-----------------------------------------------------------------------
c     Print routine for MHD vector
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in) :: xvec
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdprint'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  write(*,*) '  ', i, j, k, l, xvec(i,j,k,l)
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdprint
c=======================================================================



      subroutine fnvmhdlinsum(aval, xvec, bval, yvec, zvec)
c-----------------------------------------------------------------------
c     z = a*x + b*y
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      double precision :: aval, bval
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdlinsum, a =',aval,', b =',bval

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = aval*xvec(i,j,k,l)+bval*yvec(i,j,k,l)
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdlinsum
c=======================================================================



      subroutine fnvmhdconst(cval, zvec)
c-----------------------------------------------------------------------
c     z = c 
c     (also perform on ghost cells, since often follows initialization)
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      double precision :: cval
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdconst'

      do l=1,nvar,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  zvec(i,j,k,l) = cval
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdconst
c=======================================================================



      subroutine fnvmhdprod(xvec, yvec, zvec)
c-----------------------------------------------------------------------
c     z = x.*y
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdprod'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = xvec(i,j,k,l)*yvec(i,j,k,l)
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdprod
c=======================================================================



      subroutine fnvmhddiv(xvec, yvec, zvec)
c-----------------------------------------------------------------------
c     z = x./y (no check for zero denominator)
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhddiv'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = xvec(i,j,k,l)/yvec(i,j,k,l)
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhddiv
c=======================================================================



      subroutine fnvmhdscale(cval, xvec, zvec)
c-----------------------------------------------------------------------
c     z = c*x
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      double precision :: cval
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdscale'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = cval*xvec(i,j,k,l)
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdscale
c=======================================================================



      subroutine fnvmhdabs(xvec, zvec)
c-----------------------------------------------------------------------
c     z = |x|
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdabs'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = dabs(xvec(i,j,k,l))
               end do
            end do
         end do
      end do

      return
      end subroutine fnvmhdabs
c=======================================================================



      subroutine fnvmhdinv(xvec, zvec)
c-----------------------------------------------------------------------
c     z = 1./x (no check for division by zero)
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdinv'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = 1.D0/xvec(i,j,k,l)
               end do
            end do
         end do
      end do

      return
      end subroutine fnvmhdinv
c=======================================================================



      subroutine fnvmhdaddconst(cval, xvec, zvec)
c-----------------------------------------------------------------------
c     z = x+c
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      double precision :: cval
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdaddconst'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  zvec(i,j,k,l) = xvec(i,j,k,l) + cval
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdaddconst
c=======================================================================



      subroutine fnvmhddotprod(xvec, yvec, cval)
c-----------------------------------------------------------------------
c     c = <x,y>
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhddotprod'

      cvalloc = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cvalloc = cvalloc + xvec(i,j,k,l)*yvec(i,j,k,l)
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhddotprod MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhddotprod
c=======================================================================



      subroutine fnvmhddotprod_rho(xvec, yvec, cval)
c-----------------------------------------------------------------------
c     c = <x_rho,y_rho>
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision :: cval, cvalloc
      integer :: i, j, k, mperr

c=======Internals ============

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cvalloc = cvalloc + xvec(i,j,k,1)*yvec(i,j,k,1)
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhddotprod_rho MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhddotprod_rho
c=======================================================================



      subroutine fnvmhddotprod_u(xvec, yvec, cval)
c-----------------------------------------------------------------------
c     c = <x_u,y_u>
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============

      cvalloc = 0.d0
      do l=2,4,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cvalloc = cvalloc + xvec(i,j,k,l)*yvec(i,j,k,l)
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhddotprod_u MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhddotprod_u
c=======================================================================



      subroutine fnvmhddotprod_B(xvec, yvec, cval)
c-----------------------------------------------------------------------
c     c = <x_B,y_B>
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============

      cvalloc = 0.d0
      do l=5,7,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cvalloc = cvalloc + xvec(i,j,k,l)*yvec(i,j,k,l)
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhddotprod_B MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhddotprod_B
c=======================================================================



      subroutine fnvmhddotprod_e(xvec, yvec, cval)
c-----------------------------------------------------------------------
c     c = <x_e,y_e>
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision :: cval, cvalloc
      integer :: i, j, k, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhddotprod'

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cvalloc = cvalloc + xvec(i,j,k,8)*yvec(i,j,k,8)
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhddotprod_e MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhddotprod_e
c=======================================================================



      subroutine fnvmhdmaxnorm(xvec, cval)
c-----------------------------------------------------------------------
c     c = max(|x|)
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdmaxnorm'

      cvalloc = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cval = dabs(xvec(i,j,k,l))
                  if (cval > cvalloc) cvalloc = cval
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      cval = 0.d0
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MAX, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdmaxnorm MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdmaxnorm
c=======================================================================



      subroutine fnvmhdmaxnormvar(xvec, var, cval)
c-----------------------------------------------------------------------
c     c = max(|x_var|)
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      integer, intent(in) :: var
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdmaxnormvar'

      if ((var < 0) .or. (var > nvar)) then
         write(0,*) 'fnvmaxnormvar ERROR, illegal var = ',var
         cval = 0.d0
         return
      endif

      cvalloc = xvec(1,1,1,var)
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cval = dabs(xvec(i,j,k,var))
               if (cval > cvalloc) cvalloc = cval
            end do
         end do
      end do

#ifdef PARALLEL
      cval = 0.d0
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MAX, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdmaxnormvar MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdmaxnormvar
c=======================================================================



      subroutine fnvmhdwrmsnorm(xvec, wvec, cval)
c-----------------------------------------------------------------------
c     cval = ||x.*w||_2 / Ntotal
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, wvec
      double precision :: cval, cvalloc, dtmp
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdwrmsnorm'

      cvalloc = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cval = xvec(i,j,k,l)*wvec(i,j,k,l)
                  cvalloc = cvalloc + cval*cval
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, dtmp, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdwrmsnorm MPI_AllReduce error = ',mperr
      endif
#else
      dtmp = cvalloc
#endif
      cval = dsqrt(dtmp/nx/ny/nz/8)
      
      return
      end subroutine fnvmhdwrmsnorm
c=======================================================================



      subroutine fnvmhdwrmsnormvar(xvec, wvec, var, cval)
c-----------------------------------------------------------------------
c     cval = ||x_var.*w_var||_2 / Ntotal
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, wvec
      integer, intent(in) :: var
      double precision :: cval, cvalloc, dtmp
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdwrmsnormvar'

      if ((var < 0) .or. (var > nvar)) then
         write(0,*) 'fnvwrmsnormvar ERROR, illegal var = ',var
         cval = 0.d0
         return
      endif

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cval = xvec(i,j,k,var)*wvec(i,j,k,var)
               cvalloc = cvalloc + cval*cval
            end do
         end do
      end do

#ifdef PARALLEL
      if (my_id == 0)  cvalloc = cvalloc
      call MPI_AllReduce(cvalloc, dtmp, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdwrmsnormvar MPI_AllReduce error = ',mperr
      endif
#else
      cvalloc = cvalloc
      dtmp = cvalloc
#endif
      cval = dsqrt(dtmp/(nx*ny*nz))
      
      return
      end subroutine fnvmhdwrmsnormvar
c=======================================================================



      subroutine fnvmhdwrmsnormmask(xvec, wvec, locvec, cval)
c-----------------------------------------------------------------------
c     cval = ||x.*w||_2 / Ntotal (but only where locvec > 0)
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, wvec, locvec
      double precision :: cval, ctmp(2), dtmp(2)
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdwrmsnormmask'

      ctmp = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  if (locvec(i,j,k,l) > 0.d0) then
                     cval = xvec(i,j,k,l)*wvec(i,j,k,l)
                     ctmp(1) = ctmp(1) + cval*cval
                     ctmp(2) = ctmp(2) + 1.d0
                  endif
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(ctmp, dtmp, 2, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdwrmsnormmask MPI_AllReduce error = ',mperr
      endif
#else
      dtmp = ctmp
#endif
      cval = dsqrt(dtmp(1)/dtmp(2))
      
      return
      end subroutine fnvmhdwrmsnormmask
c=======================================================================



      subroutine fnvmhdmin(xvec, cval)
c-----------------------------------------------------------------------
c     cval = min(xvec)
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdmin'

      cvalloc = xvec(1,1,1,1)
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cval = xvec(i,j,k,l)
                  if (cval < cvalloc) cvalloc = cval
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      cval = 0.d0
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MIN, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdmin MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdmin
c=======================================================================



      subroutine fnvmhdminvar(xvec, var, cval)
c-----------------------------------------------------------------------
c     cval = min(xvec_var)
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      integer, intent(in) :: var
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdminvar'

      if ((var < 0) .or. (var > nvar)) then
         write(0,*) 'fnvmhdminvar ERROR, illegal var = ',var
         cval = 0.d0
         return
      endif

      cvalloc = xvec(1,1,1,var)
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cval = xvec(i,j,k,var)
               if (cval < cvalloc) cvalloc = cval
            end do
         end do
      end do

#ifdef PARALLEL
      cval = 0.d0
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MIN, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdminvar MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdminvar
c=======================================================================



      subroutine fnvmhdwl2norm(xvec, wvec, cval)
c-----------------------------------------------------------------------
c     c = ||x.*w||_2
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, wvec
      double precision :: cval, cvalloc, dtmp
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdwl2norm'

      cvalloc = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cval = xvec(i,j,k,l)*wvec(i,j,k,l)
                  cvalloc = cvalloc + cval*cval
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, dtmp, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdwl2norm MPI_AllReduce error = ',mperr
      endif
#else
      dtmp = cvalloc
#endif
      cval = dsqrt(dtmp)
      
      return
      end subroutine fnvmhdwl2norm
c=======================================================================



      subroutine fnvmhdwl2normvar(xvec, wvec, var, cval)
c-----------------------------------------------------------------------
c     c = ||x_var.*w_var||_2
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, wvec
      integer, intent(in) :: var
      double precision :: cval, cvalloc, dtmp
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdwl2normvar'

      if ((var < 0) .or. (var > nvar)) then
         write(0,*) 'fnvmhdwl2normvar ERROR, illegal var = ',var
         cval = 0.d0
         return
      endif

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cval = xvec(i,j,k,var)*wvec(i,j,k,var)
               cvalloc = cvalloc + cval*cval
            end do
         end do
      end do

#ifdef PARALLEL
      if (my_id == 0)  cvalloc = cvalloc
      call MPI_AllReduce(cvalloc, dtmp, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdwl2normvar MPI_AllReduce error = ',mperr
      endif
#else
      cvalloc = cvalloc
      dtmp = cvalloc
#endif
      cval = dsqrt(dtmp)
      
      return
      end subroutine fnvmhdwl2normvar
c=======================================================================



      subroutine fnvmhdl1norm(xvec, cval)
c-----------------------------------------------------------------------
c     c = ||x||_1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdl1norm'

      cvalloc = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cvalloc = cvalloc + dabs(xvec(i,j,k,l))
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdl1norm MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdl1norm
c=======================================================================



      subroutine fnvmhdl1normvar(xvec, var, cval)
c-----------------------------------------------------------------------
c     c = ||x_var||_1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      integer, intent(in) :: var
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdl1normvar'

      if ((var < 0) .or. (var > nvar)) then
         write(0,*) 'fnvmhdl1normvar ERROR, illegal var = ',var
         cval = 0.d0
         return
      endif

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cvalloc = cvalloc + dabs(xvec(i,j,k,var))
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdl1normvar MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdl1normvar
c=======================================================================



      subroutine fnvmhdl1norm_rho(xvec, cval)
c-----------------------------------------------------------------------
c     c = ||x_rho||_1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, mperr

c=======Internals ============

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cvalloc = cvalloc + dabs(xvec(i,j,k,1))
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdl1norm_rho MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdl1norm_rho
c=======================================================================



      subroutine fnvmhdl1norm_u(xvec, cval)
c-----------------------------------------------------------------------
c     c = ||x_u||_1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============

      cvalloc = 0.d0
      do l=2,4,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cvalloc = cvalloc + dabs(xvec(i,j,k,l))
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdl1norm_u MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdl1norm_u
c=======================================================================



      subroutine fnvmhdl1norm_B(xvec, cval)
c-----------------------------------------------------------------------
c     c = ||x_B||_1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============

      cvalloc = 0.d0
      do l=5,7,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  cvalloc = cvalloc + dabs(xvec(i,j,k,l))
               end do
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdl1norm_B MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdl1norm_B
c=======================================================================



      subroutine fnvmhdl1norm_e(xvec, cval)
c-----------------------------------------------------------------------
c     c = ||x_e||_1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision :: cval, cvalloc
      integer :: i, j, k, mperr

c=======Internals ============

      cvalloc = 0.d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               cvalloc = cvalloc + dabs(xvec(i,j,k,8))
            end do
         end do
      end do

#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_SUM, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdl1norm_e MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif
      
      return
      end subroutine fnvmhdl1norm_e
c=======================================================================



      subroutine fnvmhdcompare(cval, xvec, zvec)
c-----------------------------------------------------------------------
c     z(i) = 1 if x(i) > c, otherwise z(i) = 0
c-----------------------------------------------------------------------
      use mesh
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      double precision :: cval
      integer :: i, j, k, l

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdcompare'

      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  if (xvec(i,j,k,l) > cval) then
                     zvec(i,j,k,l) = 1.d0
                  else
                     zvec(i,j,k,l) = 0.d0
                  endif
               end do
            end do
         end do
      end do
      
      return
      end subroutine fnvmhdcompare
c=======================================================================



      subroutine fnvmhdinvtest(xvec, zvec, cval)
c-----------------------------------------------------------------------
c     z = 1./x, if all x nonzero, cval = 0, else cval = 1
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: zvec
      integer :: cval, cvalloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdinvtest'

      cvalloc = 0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  if (xvec(i,j,k,l) /= 0.d0) then
                     zvec(i,j,k,l) = 1.d0/xvec(i,j,k,l)
                  else
                     zvec(i,j,k,l) = 0.d0
                     cvalloc = 1
                  endif
               end do
            end do
         end do
      end do
      
#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_INTEGER, 
     &     MPI_MAX, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdinvtest MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhdinvtest
c=======================================================================



      subroutine fnvmhdconstrmask(cvec, xvec, mvec, testval)
c-----------------------------------------------------------------------
c     Returns testval = 1 if any element fails constraint test,
c     otherwise testval = 0.  Test is as follows:
c     if c(i) =  2, then x(i) >  0
c     if c(i) =  1, then x(i) >= 0
c     if c(i) = -1, then x(i) <= 0
c     if c(i) = -2, then x(i) <  0
c     fills m(i) = 1 where test fails, m(i) = 0 elsewhere.
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: cvec, xvec
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(out) :: mvec
      integer :: testval, testloc
      integer :: i, j, k, l, mperr

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdconstrmask'

      testloc = 0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  mvec(i,j,k,l) = 0.d0
                  if (dabs(cvec(i,j,k,l)) > 1.5d0) then
                     if (xvec(i,j,k,l)*cvec(i,j,k,l) <= 0.d0) then
                        testloc = 1
                        mvec(i,j,k,l) = 1.d0
                     endif
                  else if (dabs(cvec(i,j,k,l)) > 0.5d0) then
                     if (xvec(i,j,k,l)*cvec(i,j,k,l) < 0.d0) then
                        testloc = 1
                        mvec(i,j,k,l) = 1.d0
                     endif                     
                  endif
               end do
            end do
         end do
      end do
      
#ifdef PARALLEL
      call MPI_AllReduce(testloc, testval, 1, MPI_INTEGER, 
     &     MPI_MAX, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdconstrmask MPI_AllReduce error = ',mperr
      endif
#else
      testval = testloc
#endif

      return
      end subroutine fnvmhdconstrmask
c=======================================================================



      subroutine fnvmhdminquotient(xvec, yvec, cval)
c-----------------------------------------------------------------------
c     c = min(x./y), over all y /= 0
c-----------------------------------------------------------------------
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar),
     &     intent(in)  :: xvec, yvec
      double precision :: cval, cvalloc
      integer :: i, j, k, l, mperr
      logical :: notyet

c=======Internals ============
c$$$  write(*,*) 'entering fnvmhdminquotient'

      notyet = .true.
      cvalloc = 0.d0
      do l=1,nvar,1
         do k=1,nzlocal,1
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  if (yvec(i,j,k,l) /= 0.d0) then
                     cval = xvec(i,j,k,l)/yvec(i,j,k,l)
                     if (notyet) then
                        cvalloc = cval
                        notyet = .false.
                     else
                        cvalloc = min(cvalloc, cval)
                     endif
                  endif
               end do
            end do
         end do
      end do
      
#ifdef PARALLEL
      call MPI_AllReduce(cvalloc, cval, 1, MPI_DOUBLE_PRECISION, 
     &     MPI_MIN, comm3d, mperr)
      if (mperr /= 0 ) then
         write(0,*) 'fnvmhdminquotient MPI_AllReduce error = ',mperr
      endif
#else
      cval = cvalloc
#endif

      return
      end subroutine fnvmhdminquotient
c=======================================================================
