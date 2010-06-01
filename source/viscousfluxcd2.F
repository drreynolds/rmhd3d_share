C -*- Mode: Fortran; -*-
c     Ravi Samtaney
c     Copyright 2004
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c     
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c     
c-----------------------------------------------------------------------
c     $Log: viscousfluxcd2.F,v $
c-----------------------------------------------------------------------
c
      subroutine ViscousFluxCD2(fvis,vx,dudx,dbdx,
     &     direction,dt,indx,indy,indz)
c-----------------------------------------------------------------------
c     Description:  Calculation of parabolic or diffusive fluxes 
c     Input: vx - primitive variables
c            direction - 1,2,3 (=x,y,z)
c            dt - time step (not used)
c            dudx - grad(velocity)
c            dbdx - grad(magnetic field)
c            dudx and dbdx are precomputed.
c            indx,indy,indz - specify which direction (sort of 
c            redundant to specify both direction and these)
c     Output: fvis - Parabolic fluxes 
c     Issues: eta is computed using Bz - this makes it nonlinear
c-----------------------------------------------------------------------
      use mesh
      use properties

      implicit none
c     
c     This subroutine determines the vector of inviscid fluxes.
      integer:: indx,indy,indz

      double precision:: fvis(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)

      double precision:: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: dudx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,3,3)
      double precision:: dbdx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,3,3)
      double precision:: dt
      integer:: direction,ilo,ihi
      integer:: i,j,k
c     
      double precision:: vd(INLO:INHI,nvar)
      double precision:: vdjp(INLO:INHI,nvar)
      double precision:: vdjm(INLO:INHI,nvar)
      double precision:: vdkp(INLO:INHI,nvar)
      double precision:: vdkm(INLO:INHI,nvar)
      double precision:: etaPrime(INLO:INHI)
      double precision:: fx(INLO:INHI,nvar)
      double precision:: dxi,dyi,dzi
c     This subroutine determines the vector of diffusive fluxes.
c     v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c     rho=density, u=x-velocity, v=y-velocity
c     w=z-velocity
c     Bx,By,Bz=x,y,z components of the magnetic field
c     p=thermal pressure
c
      dxi=1.D0/dx
      dyi=1.D0/dy
      dzi=1.D0/dz
      fx=0.D0
c     
      select case  (direction)
      case(1)
         do k=1,nzlocal,1
            do j=1,nylocal,1
               ilo=ixlo; ihi=ixhi;
               vd(ilo:ihi,:)=vx(ilo:ihi,j,k,:)
               vdjp(ilo:ihi,:)=vx(ilo:ihi,j+1,k,:)
               vdjm(ilo:ihi,:)=vx(ilo:ihi,j-1,k,:)
#if (defined(ONE_D) || defined(TWO_D))
               vdkp(ilo:ihi,:)=vx(ilo:ihi,j,k,:)
               vdkm(ilo:ihi,:)=vx(ilo:ihi,j,k,:)
#else
               vdkp(ilo:ihi,:)=vx(ilo:ihi,j,k+1,:)
               vdkm(ilo:ihi,:)=vx(ilo:ihi,j,k-1,:)
#endif
               ilo=0; ihi=nxlocal;
               call ComputeResistivity(etaPrime,vd(inlo,7),ilo,ihi)
               call ViscousFlux1DCD2(vd,vdjp,vdjm,vdkp,
     &              vdkm,etaPrime,fx,ilo,ihi,dxi,dyi,dzi)
               fvis(ilo:ihi,j,k,:)=fx(ilo:ihi,:)
            enddo
         enddo
      case(2)
         do k=1,nzlocal,1
            do i=1,nxlocal,1
               ilo=iylo; ihi=iyhi;
               vd(ilo:ihi,:)=vx(i,ilo:ihi,k,:)
#if (defined(ONE_D) || defined(TWO_D))
               vdjp(ilo:ihi,:)=vx(i,ilo:ihi,k,:)
               vdjm(ilo:ihi,:)=vx(i,ilo:ihi,k,:)
#else
               vdjp(ilo:ihi,:)=vx(i,ilo:ihi,k+1,:)
               vdjm(ilo:ihi,:)=vx(i,ilo:ihi,k-1,:)
#endif
               vdkp(ilo:ihi,:)=vx(i+1,ilo:ihi,k,:)
               vdkm(ilo:ihi,:)=vx(i-1,ilo:ihi,k,:)
               call RemapToLocalCoords(vd,ilo,ihi,2)
               call RemapToLocalCoords(vdjp,ilo,ihi,2)
               call RemapToLocalCoords(vdkp,ilo,ihi,2)
               call RemapToLocalCoords(vdjm,ilo,ihi,2)
               call RemapToLocalCoords(vdkm,ilo,ihi,2)
               ilo=0; ihi=nylocal;
               call ComputeResistivity(etaPrime,vd(inlo,6),ilo,ihi)
               call ViscousFlux1DCD2(vd,vdjp,vdjm,vdkp,
     &              vdkm,etaPrime,fx,ilo,ihi,dyi,dzi,dxi)
               fvis(i,ilo:ihi,k,:)=fx(ilo:ihi,:)
            enddo
         enddo
      case(3)
         do j=1,nylocal,1
            do i=1,nxlocal,1
               ilo=izlo; ihi=izhi;
               vd(ilo:ihi,:)=vx(i,j,ilo:ihi,:)
               vdjp(ilo:ihi,:)=vx(i+1,j,ilo:ihi,:)
               vdjm(ilo:ihi,:)=vx(i-1,j,ilo:ihi,:)
               vdkp(ilo:ihi,:)=vx(i,j+1,ilo:ihi,:)
               vdkm(ilo:ihi,:)=vx(i,j-1,ilo:ihi,:)
               call RemapToLocalCoords(vd,ilo,ihi,3)
               call RemapToLocalCoords(vdjp,ilo,ihi,3)
               call RemapToLocalCoords(vdkp,ilo,ihi,3)
               call RemapToLocalCoords(vdjm,ilo,ihi,3)
               call RemapToLocalCoords(vdkm,ilo,ihi,3)
               ilo=0; ihi=nzlocal;
               call ComputeResistivity(etaPrime,vd(inlo,5),ilo,ihi)
               call ViscousFlux1DCD2(vd,vdjp,vdjm,vdkp,
     &              vdkm,etaPrime,fx,ilo,ihi,dzi,dxi,dyi)
               fvis(i,j,ilo:ihi,:)=fx(ilo:ihi,:)
            enddo
         enddo
c     
      end select
c     
      return
      end subroutine ViscousFluxCD2
c
c-----------------------------------------------------------------------
      subroutine ViscousFlux1DCD2(vd,vdjp,vdjm,vdkp,
     &	       vdkm,etaPrime,fvis,ilo,ihi,adx,ady,adz)
c-----------------------------------------------------------------------
c     Description:  Calculation of parabolic or diffusive fluxes 
c     Input: vd - primitive variables
c            vdjp - vd + 1 in j dir
c            vdjm - vd - 1 in j dir
c            vdkp - vd + 1 in k dir
c            vdkm - vd - 1 in k dir
c            etaPrime - resistivity
c     Output: fvis - Parabolic fluxes 
c-----------------------------------------------------------------------
      use mesh
      use properties
      double precision:: fvis(inlo:inhi,nvar)
      double precision:: vd(inlo:inhi,nvar)
      double precision:: vdjp(inlo:inhi,nvar)
      double precision:: vdjm(inlo:inhi,nvar)
      double precision:: vdkp(inlo:inhi,nvar)
      double precision:: vdkm(inlo:inhi,nvar)
      double precision:: etaPrime(inlo:inhi)
      double precision:: tauxx,tauxy,tauxz
      double precision:: dbxdx,dbxdy,dbxdz
      double precision:: dbydx, dbzdx,dbmagdx,bi,bj,bk
      double precision:: dudx,dudy,dudz,div
      double precision:: dvdx,dwdx,dvdy,dwdz
      double precision:: dtdx
      double precision:: adx,ady,adz
      integer:: ilo,ihi
      integer:: i,j,k,l
c	 This subroutine determines the vector of inviscid fluxes.
c
c	 v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c	    rho=density, u=x-velocity, v=y-velocity
c           w=z-velocity
c           Bx,By,Bz= Magnetic Field
c	    p=pressure
c
         do i=ilo,ihi,1
c
            dudx=(vd(i+1,2)-vd(i,2))*adx
            dvdx=(vd(i+1,3)-vd(i,3))*adx
            dwdx=(vd(i+1,4)-vd(i,4))*adx
c     Note: these are multiplied by 0.25 since we are averaging the
c     upper and lower states (accounts for 0.5), and since they are
c     2*dy (or 2*dz, respectively) apart (the other 0.5)
            dudy=0.25D0*(vdjp(i,2)+vdjp(i+1,2)-vdjm(i,2)-vdjm(i+1,2))
     &           *ady
            dvdy=0.25D0*(vdjp(i,3)+vdjp(i+1,3)-vdjm(i,3)-vdjm(i+1,3))
     &           *ady
c            dwdy=0.25D0*(vdjp(i,4)+vdjp(i+1,4)-vdjm(i,4)-vdjm(i+1,4))
c
            dudz=0.25D0*(vdkp(i,2)+vdkp(i+1,2)-vdkm(i,2)-vdkm(i+1,2))
     &           *adz
c            dvdz=0.25D0*(vdkp(i,3)+vdkp(i+1,3)-vdkm(i,3)-vdkm(i+1,3))
            dwdz=0.25D0*(vdkp(i,4)+vdkp(i+1,4)-vdkm(i,4)-vdkm(i+1,4))
     &           *adz
c
            div=(dudx+dvdy+dwdz)/3.D0
c
            tauxx=2.d0*mu*(dudx-div)
            tauxy=mu*(dvdx+dudy)
            tauxz=mu*(dwdx+dudz)
c
            dbxdx=(vd(i+1,5)-vd(i,5))*adx
            dbydx=(vd(i+1,6)-vd(i,6))*adx
            dbzdx=(vd(i+1,7)-vd(i,7))*adx
c
            dbxdy=0.25D0*(vdjp(i,5)+vdjp(i+1,5)-vdjm(i,5)-vdjm(i+1,5))
     &           *ady
            dbxdz=0.25D0*(vdkp(i,5)+vdkp(i+1,5)-vdkm(i,5)-vdkm(i+1,5))
     &           *adz
c            dbmagdx=(sum(vd(i+1,5:7)*vd(i+1,5:7))-
c     &           (sum(vd(i,5:7)*vd(i,5:7))))*0.5D0*adx
            dbmagdx=(vd(i+1,5)**2+vd(i+1,6)**2+vd(i+1,7)**2-
     &           vd(i,5)**2-vd(i,6)**2-vd(i,7)**2)*0.5D0*adx
c
            bi=0.5D0*(vd(i,5)+vd(i+1,5))
            bj=0.5D0*(vd(i,6)+vd(i+1,6))
            bk=0.5D0*(vd(i,7)+vd(i+1,7))
c
            dtdx=(vd(i+1,8)/vd(i+1,1)-vd(i,8)/vd(i,1))*adx/rgas
c
            fvis(i,2)=tauxx
            fvis(i,3)=tauxy
            fvis(i,4)=tauxz
ccccccccccc  CHANGED FOLLOWING TO FIX THE DIV(B)=0 CONDITION, 
ccccccccccc  AS PER RAVI'S REQUEST ON 2/23/2004
c            fvis(i,6)=etaPrime(i)*(dbydx-dbxdy)
c            fvis(i,7)=etaPrime(i)*(dbzdx-dbxdz)
            fvis(i,5) = etaPrime(i)*(dbxdx)
            fvis(i,6) = etaPrime(i)*(dbydx)
            fvis(i,7) = etaPrime(i)*(dbzdx)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            fvis(i,8)=etaPrime(i)*(dbmagdx-bi*dbxdx-bj*dbxdy-bk*dbxdz)+
     &           0.5D0*((vd(i,2)+vd(i+1,2))*tauxx+
     &                  (vd(i,3)+vd(i+1,3))*tauxy+
     &                  (vd(i,4)+vd(i+1,4))*tauxz)+
     &           gamma/(gamma-1.D0)/0.7D0*kappa*dtdx
         enddo
c
c
	return
        end subroutine ViscousFlux1DCD2
c
c
c-----------------------------------------------------------------------
      subroutine computeResistivity(etaPrime,vd,ilo,ihi)
c-----------------------------------------------------------------------
c     Description: Compute the nonlinear resistivity
c     Input: vd (primitive variables at cell centers)
c     Output: etaPrime at cell face
c-----------------------------------------------------------------------
      use mesh
      use properties
      double precision:: etaPrime(inlo:inhi)
      double precision:: vd(inlo:inhi)
      integer:: ilo,ihi
      integer:: i,j,k,l
c     	 This subroutine determines the mag field dependent
c     resistivity
c
      etaPrime=eta
      return
c         do i=ilo,(ihi+ilo)/2+1,1
         do i=ilo,ihi,1
c            if(vd(i)+vd(i+1).le.-2.D-04) then
            if(vd(i)+vd(i+1).le.0.D0) then
               etaPrime(i)=eta+etaFactor*eta*(1.D0-
     &              dexp(-((vd(i)+vd(i+1))/2.D0/0.075D0/7.8D0)**2))
            else
              etaPrime(i)=eta
            endif
         enddo
c         do i=(ihi+ilo)/2+1,ihi,1
c            etaPrime(i)=etaPrime(ihi-i+ilo)
c         enddo
         return
         end subroutine computeResistivity
c-----------------------------------------------------------------------

