c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------------
c     $Log: inviscidfluxcd4.F,v $
c=======================================================================




      subroutine InviscidFluxCD4(finv,vx,direction,dt)
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
     &                 :: fx, vd, vl, vr, vll,vrr
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
                  vll(i,:) = vd(i-1,:)
                  vrr(i,:) = vd(i+2,:)
               end do
               ilo = 0;  ihi = nxlocal
               call InviscidFluxCD4Cal(fx,vl,vr,vll,vrr,ilo,ihi)
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
                  vll(j,:) = vd(j-1,:)
                  vrr(j,:) = vd(j+2,:)
               end do
               ilo = 0;  ihi = nylocal
               call InviscidFluxCD4Cal(fx,vl,vr,vll,vrr,ilo,ihi)
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
                  vll(k,:) = vd(k-1,:)
                  vrr(k,:) = vd(k+2,:)
               end do
               ilo = 0;  ihi = nzlocal
               call InviscidFluxCD4Cal(fx,vl,vr,vll,vrr,ilo,ihi)
               finv(i,j,ilo:ihi,:) = fx(ilo:ihi,:)
            enddo
         enddo
     
      end select
     
      return
      end subroutine InviscidFluxCD4
c=======================================================================





      subroutine InviscidFluxCD4Cal(finv,vl,vr,vll,vrr,ilo,ihi)
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
      double precision, dimension(inlo:inhi,nvar) :: finv, vl, vr,
     &     vll,vrr
      integer :: ilo, ihi

c     local variables
      double precision, dimension(nvar) :: fl, fr
      double precision, dimension(nvar) :: fll, frr
      double precision :: rho, u, v, w, bx, by, bz, eint
      double precision :: entropy, press
      double precision :: a2,a1,a0,am1
      double precision :: d2,d1
      integer :: i,l


c     Fourth order standard central difference
      d1=2.D0/3.D0
      d2=-1.D0/12.D0

c     compute the corresponding weights
      a2=d2
      a1=d1+d2
      a0=a1
      am1=a2
      
c======= Internals ===========

c     loop over indices
      do i = ilo,ihi,1

c        flux from left side of face minus one
         rho  = vll(i,1)
         u    = vll(i,2)  
         v    = vll(i,3)  
         w    = vll(i,4)
         bx   = vll(i,5)  
         by   = vll(i,6)  
         bz   = vll(i,7)
         eint = vll(i,8)/(gamma-1.D0)
         press = vll(i,8)
         fll(1) = rho*u
         fll(2) = rho*u*u + press + 0.5D0*(by*by+bz*bz-bx*bx)
         fll(3) = rho*u*v - bx*by
         fll(4) = rho*u*w - bx*bz
         fll(5) = 0.0D0
         fll(6) = u*by - v*bx
         fll(7) = u*bz - w*bx
         fll(8) = (0.5D0*rho*(u*u+v*v+w*w) + eint 
     &        + vll(i,8) + (bx*bx+by*by+bz*bz))*vll(i,2)
     &        - bx*(u*bx+v*by+w*bz)

c        flux from left side of face
         rho  = vl(i,1)
         u    = vl(i,2)  
         v    = vl(i,3)  
         w    = vl(i,4)
         bx   = vl(i,5)  
         by   = vl(i,6)  
         bz   = vl(i,7)
         eint = vl(i,8)/(gamma-1.D0)
         press = vl(i,8)
         fl(1) = rho*u
         fl(2) = rho*u*u + press + 0.5D0*(by*by+bz*bz-bx*bx)
         fl(3) = rho*u*v - bx*by
         fl(4) = rho*u*w - bx*bz
         fl(5) = 0.0D0
         fl(6) = u*by - v*bx
         fl(7) = u*bz - w*bx
         fl(8) = (0.5D0*rho*(u*u+v*v+w*w) + eint 
     &        + vl(i,8) + (bx*bx+by*by+bz*bz))*vl(i,2)
     &        - bx*(u*bx+v*by+w*bz)

         
c        flux from right side of face
         rho  = vr(i,1)
         u    = vr(i,2) 
         v    = vr(i,3) 
         w    = vr(i,4)
         bx   = vr(i,5) 
         by   = vr(i,6) 
         bz   = vr(i,7)
         eint = vr(i,8)/(gamma-1.D0)
         press = vr(i,8)
         fr(1) = rho*u
         fr(2) = rho*u*u + press + 0.5D0*(by*by+bz*bz-bx*bx)
         fr(3) = rho*u*v - bx*by
         fr(4) = rho*u*w - bx*bz
         fr(5) = 0.0D0
         fr(6) = u*by - v*bx
         fr(7) = u*bz - w*bx
         fr(8) = (0.5D0*rho*(u*u+v*v+w*w) + eint
     &        + vr(i,8) + (bx*bx+by*by+bz*bz))*vr(i,2)
     &        - bx*(u*bx+v*by+w*bz)

c        flux from right side of face plus one
         rho  = vrr(i,1)
         u    = vrr(i,2) 
         v    = vrr(i,3) 
         w    = vrr(i,4)
         bx   = vrr(i,5) 
         by   = vrr(i,6) 
         bz   = vrr(i,7)
         eint = vrr(i,8)/(gamma-1.D0)
         press = vrr(i,8)
         frr(1) = rho*u
         frr(2) = rho*u*u + press + 0.5D0*(by*by+bz*bz-bx*bx)
         frr(3) = rho*u*v - bx*by
         frr(4) = rho*u*w - bx*bz
         frr(5) = 0.0D0
         frr(6) = u*by - v*bx
         frr(7) = u*bz - w*bx
         frr(8) = (0.5D0*rho*(u*u+v*v+w*w) + eint
     &        + vrr(i,8) + (bx*bx+by*by+bz*bz))*vrr(i,2)
     &        - bx*(u*bx+v*by+w*bz)

c        central difference at face comes from arithmetic average
         do l=1,nvar,1
            finv(i,l) = a2*frr(l)+a1*fr(l)+a0*fl(l)+am1*fll(l)
         enddo
      enddo

      return
      end subroutine InviscidFluxCD4Cal
c=======================================================================
