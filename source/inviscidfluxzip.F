c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------------
c     $Log: inviscidfluxzip.F,v $
c=======================================================================
c     Uses zip differencing for all except B



      subroutine InviscidFluxZip(finv,vx,direction,dt)
c-----------------------------------------------------------------------
c     Description:  Calculation of Hyperbolic fluxes 
c     Input: vx - primitive variables
c            direction - 1,2,3 (=x,y,z)
c            dt - time step (not used)
c     Output: finv - Hyperbolic fluxes 
c-----------------------------------------------------------------------
c======= Inclusions ==========
      use mesh

c======= Declarations ========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &                 :: finv, vx
      double precision :: dt
      integer          :: direction

c     local variables
      double precision, dimension(inlo:inhi,nvar) 
     &                 :: fx, vd, vl, vr
      integer          :: ilo, ihi, i, j, k


c======= Internals ===========

      select case (direction)

      case(1)
         do k = 1,nzlocal
            do j = 1,nylocal
               ilo = ixlo;  ihi = ixhi
               vd(ilo:ihi,:) = vx(ilo:ihi,j,k,:)
               do i=ilo+1,ihi-nghost,1
                  vl(i,:) = vd(i,:)
                  vr(i,:) = vd(i+1,:)
               end do
               ilo = 0;  ihi = nxlocal
               call InviscidFluxZipCal(fx,vl,vr,ilo,ihi)
               finv(ilo:ihi,j,k,:) = fx(ilo:ihi,:)
            enddo
         enddo

      case(2)
         do k = 1,nzlocal
            do i = 1,nxlocal
               ilo = iylo;  ihi = iyhi
               vd(ilo:ihi,:) = vx(i,ilo:ihi,k,:)
               call RemapToLocalCoords(vd,ilo,ihi,2)
               do j=ilo+1,ihi-nghost,1
                  vl(j,:) = vd(j,:)
                  vr(j,:) = vd(j+1,:)
               end do
               ilo = 0;  ihi = nylocal
               call InviscidFluxZipCal(fx,vl,vr,ilo,ihi)
               finv(i,ilo:ihi,k,:) = fx(ilo:ihi,:)
            enddo
         enddo

      case(3)
         do j = 1,nylocal
            do i = 1,nxlocal
               ilo = izlo;  ihi = izhi
               vd(ilo:ihi,:) = vx(i,j,ilo:ihi,:)
               call RemapToLocalCoords(vd,ilo,ihi,3)
               do k=ilo+1,ihi-nghost,1
                  vl(k,:) = vd(k,:)
                  vr(k,:) = vd(k+1,:)
               end do
               ilo = 0;  ihi = nzlocal
               call InviscidFluxZipCal(fx,vl,vr,ilo,ihi)
               finv(i,j,ilo:ihi,:) = fx(ilo:ihi,:)
            enddo
         enddo
     
      end select
     
      return
      end subroutine InviscidFluxZip
c=======================================================================





      subroutine InviscidFluxZipCal(finv,vvl,vvr,ilo,ihi)
c-----------------------------------------------------------------------
c     Description: Calculate inviscid flux based on centered-difference
c                  approximations to the fluid variables.
c     Input: vl, vr: left/right states
c     Output: finv: inviscid flux at cell faces
c-----------------------------------------------------------------------
c======= Inclusions ==========
c      use mesh_common
      use mesh
      use properties

c======= Declarations ========
      implicit none

c     calling variables
      double precision, dimension(inlo:inhi,nvar) :: finv, vvl, vvr
      integer :: ilo, ihi

c     local variables
      double precision, dimension(nvar) :: fl, fr
      double precision :: rhol, ul, vl, wl, bxl, byl, bzl, eintl
      double precision :: rhor, ur, vr, wr, bxr, byr, bzr, eintr
      double precision :: etotl,etotr
      double precision :: pressl,pressr
      double precision:: d1
      integer :: i,l


c     Second order standard central difference
      d1=0.5D0
      
c======= Internals ===========

c     loop over indices
      do i = ilo,ihi,1

c        flux from left side of face
         rhol  = vvl(i,1)
         ul    = vvl(i,2)  
         vl    = vvl(i,3)  
         wl    = vvl(i,4)
         bxl   = vvl(i,5)  
         byl   = vvl(i,6)  
         bzl   = vvl(i,7)
         pressl = vvl(i,8)
         eintl = vvl(i,8)/(gamma-1.D0)
         etotl = eintl+0.5D0*(rhol*(ul*ul+vl*vl+wl*wl)+
     &                 bxl*bxl+byl*byl+bzl*bzl)
c
         rhor  = vvr(i,1)
         ur    = vvr(i,2) 
         vr    = vvr(i,3) 
         wr    = vvr(i,4)
         bxr   = vvr(i,5) 
         byr   = vvr(i,6) 
         bzr   = vvr(i,7)
         eintr = vvr(i,8)/(gamma-1.D0)
         pressr = vvr(i,8)
         etotr = eintr+0.5D0*(rhor*(ur*ur+vr*vr+wr*wr)+
     &                 bxr*bxr+byr*byr+bzr*bzr)
c
         fl(1) = 0.5D0*(rhol*ur+rhor*ul)
c
         fl(2) = 0.5D0*(rhol*ul*ur+rhor*ur*ul)
     &        + 0.5D0*(pressl+pressr)
     &        + 0.5D0*(byl*byr+bzl*bzr-bxl*bxr)
c
         fl(3) = 0.25D0*(rhol*ul*vr+rhol*ur*vl+rhor*ul*vr+rhor*ur*vl)
     &        - 0.5D0*(bxl*byr+bxr*byl)
c
         fl(4) = 0.25D0*(rhol*ul*wr+rhol*ur*wl+rhor*ul*wr+rhor*ur*wl)
     &        - 0.5D0*(bxl*bzr+bxr*bzl)
c
         fl(5) = 0.0D0
         fl(6) = 0.5D0*(ul*byl - vl*bxl+ur*byr - vr*bxr)
         fl(7) = 0.5D0*(ul*bzl - wl*bxl+ur*bzr - wr*bxr)
c
         fl(8) = 0.5D0*((etotl+pressl)*ur+
     &                  (etotr+pressr)*ul)
     &          +0.25D0*((bxl*bxr+byl*byr+bzl*bzr)*(ul+ur))
     &          -0.5D0*(bxl*ul*ur+bxr*ur*ul)
     &          -0.25D0*(byl*vl*ur+byl*vr*ul+byr*vl*ur+byr*vr*ul)
     &          -0.25D0*(bxl*wl*ur+bzl*wr*ul+bzr*wl*ur+bzr*wr*ul)
c
         do l=1,nvar,1
            finv(i,l) = fl(l)
         enddo
      enddo

      return
      end subroutine InviscidFluxZipCal
c=======================================================================
