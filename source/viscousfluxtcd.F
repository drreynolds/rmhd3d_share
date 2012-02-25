c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c-----------------------------------------------------------------------

      subroutine ViscousFluxTCD(fvis,vx,dudx,dbdx,
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
      integer:: i,j,k,l,m,ll,mm
c     
      double precision:: vd(INLO:INHI,nvar)
      double precision:: dvdx(INLO:INHI,3,3)
      double precision:: dmagdx(INLO:INHI,3,3)
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
               dvdx(ilo:ihi,:,:)=dudx(ilo:ihi,j,k,:,:)
               dmagdx(ilo:ihi,:,:)=dbdx(ilo:ihi,j,k,:,:)
               ilo=0; ihi=nxlocal;
               call ComputeResistivity(etaPrime,vd(inlo,7),ilo,ihi)
               call ViscousFlux1DTCD(vd,dvdx,dmagdx,
     &	       etaPrime,fx,ilo,ihi,dxi,dyi,dzi)
               fvis(ilo:ihi,j,k,:)=fx(ilo:ihi,:)
            enddo
         enddo
      case(2)
         do k=1,nzlocal,1
            do i=1,nxlocal,1
               ilo=iylo; ihi=iyhi;
               vd(ilo:ihi,:)=vx(i,ilo:ihi,k,:)
               do l=1,3,1
                  do m=1,3,1
                     ll=mod(l,3)+1
                     mm=mod(m,3)+1
                     dvdx(ilo:ihi,l,m)=
     &              dudx(i,ilo:ihi,k,ll,mm)
                     dmagdx(ilo:ihi,l,m)=
     &              dbdx(i,ilo:ihi,k,ll,mm)
                  enddo
               enddo
               call RemapToLocalCoords(vd,ilo,ihi,2)
               ilo=0; ihi=nylocal;
               call ComputeResistivity(etaPrime,vd(inlo,6),ilo,ihi)
               call ViscousFlux1DTCD(vd,dvdx,dmagdx,
     &	       etaPrime,fx,ilo,ihi,dyi,dzi,dxi)
               fvis(i,ilo:ihi,k,:)=fx(ilo:ihi,:)
            enddo
         enddo
      case(3)
         do j=1,nylocal,1
            do i=1,nxlocal,1
               ilo=izlo; ihi=izhi;
               vd(ilo:ihi,:)=vx(i,j,ilo:ihi,:)
               do l=1,3,1
                  do m=1,3,1
                     ll=mod(l+1,3)+1
                     mm=mod(m+1,3)+1
                     dvdx(ilo:ihi,l,m)=
     &              dudx(i,j,ilo:ihi,ll,mm)
                     dmagdx(ilo:ihi,l,m)=
     &              dbdx(i,j,ilo:ihi,ll,mm)
                  enddo
               enddo
               call RemapToLocalCoords(vd,ilo,ihi,3)
               ilo=0; ihi=nzlocal;
               call ComputeResistivity(etaPrime,vd(inlo,5),ilo,ihi)
               call ViscousFlux1DTCD(vd,dvdx,dmagdx,
     &	       etaPrime,fx,ilo,ihi,dzi,dxi,dyi)
               fvis(i,j,ilo:ihi,:)=fx(ilo:ihi,:)
            enddo
         enddo
c     
      end select
c     
      return
      end subroutine ViscousFluxTCD
c
c-----------------------------------------------------------------------
      subroutine ViscousFlux1DTCD(vd,dpdx,dbdx,
     &	      etaPrime,fvis,ilo,ihi,adx,ady,adz)
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
      double precision:: fc(inlo:inhi,nvar)
      double precision:: vd(inlo:inhi,nvar)
      double precision:: dpdx(inlo:inhi,3,3)
      double precision:: dbdx(inlo:inhi,3,3)
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
      double precision:: a2,a1,a0,am1
      double precision:: d2,d1
      double precision:: dxi,dyi,dzi


c     Second order tuned central differnce 
      d1=0.894D0
      d2=-0.197D0
c     Second order standard central difference
c      d1=0.5D0
c      d2=0.D0
c     Fourth order standard central difference
c      d1=2.D0/3.D0
c      d2=-1.D0/12.D0
      a2=d2
      a1=d1+d2
      a0=a1
      am1=a2
c	 This subroutine determines the vector of inviscid fluxes.
c
c	 v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c	    rho=density, u=x-velocity, v=y-velocity
c           w=z-velocity
c           Bx,By,Bz= Magnetic Field
c	    p=pressure
c
         do i=ilo-1,ihi+2,1
c
            dudx=dpdx(i,1,1)
            dvdx=dpdx(i,2,1)
            dwdx=dpdx(i,3,1)
c
            dudy=dpdx(i,1,2)
            dvdy=dpdx(i,2,2)
            dudz=dpdx(i,1,3)
            dwdz=dpdx(i,3,3)
c
            div=(dudx+dvdy+dwdz)/3.D0
c            div=0.D0
c
cccc        D.R.R. - changed the following back to the correct
cccc        viscosity, 9/15/05
            tauxx=2.d0*mu*(dudx-div)
            tauxy=mu*(dvdx+dudy)
            tauxz=mu*(dwdx+dudz)
c$$$            tauxx=mu*dudx
c$$$            tauxy=mu*dvdx
c$$$            tauxz=mu*dwdx
c
            dbxdx=dbdx(i,1,1)
            dbydx=dbdx(i,2,1)
            dbzdx=dbdx(i,3,1)
c
            dbxdy=dbdx(i,1,2)
            dbxdz=dbdx(i,1,3)
            dbmagdx=d1*(vd(i+1,5)**2+vd(i+1,6)**2+vd(i+1,7)**2-
     &           vd(i-1,5)**2-vd(i-1,6)**2-vd(i-1,7)**2)+
     &           d2*(vd(i+2,5)**2+vd(i+2,6)**2+vd(i+2,7)**2-
     &           vd(i-2,5)**2-vd(i-2,6)**2-vd(i-2,7)**2)
            dbmagdx=dbmagdx*0.5D0*adx
c
            bi=vd(i,5)
            bj=vd(i,6)
            bk=vd(i,7)
c
            dtdx=d1*(vd(i+1,8)/vd(i+1,1)-
     &           vd(i-1,8)/vd(i-1,1))+
     &           d2*(vd(i+2,8)/vd(i+2,1)-
     &           vd(i-2,8)/vd(i-2,1))
            dtdx=dtdx*adx/rgas
c
            fc(i,2)=tauxx
            fc(i,3)=tauxy
            fc(i,4)=tauxz
            fc(i,5)=etaPrime(i)*(dbxdx)
            fc(i,6)=etaPrime(i)*(dbydx)
            fc(i,7)=etaPrime(i)*(dbzdx)
c            fc(i,6)=etaPrime(i)*(dbydx-dbxdy)
c            fc(i,7)=etaPrime(i)*(dbzdx-dbxdz)
            fc(i,8)=etaPrime(i)*(dbmagdx-bi*dbxdx-bj*dbxdy-bk*dbxdz)+
     &           vd(i,2)*tauxx+
     &                  vd(i,3)*tauxy+
     &                  vd(i,4)*tauxz+
     &           gamma/(gamma-1.D0)/0.7D0*kappa*dtdx
         enddo
c
         do i=ilo,ihi,1
            do l=2,nvar,1
               fvis(i,l) = a2*fc(i+2,l)+a1*fc(i+1,l)+
     &              a0*fc(i,l)+am1*fc(i-1,l)
            enddo
         enddo
c     
	return
        end subroutine ViscousFlux1DTCD

c-----------------------------------------------------------------------
