C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2001
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------------
c     $Log: inviscidfluxroe.F,v $
c-----------------------------------------------------------------------
      subroutine InviscidFluxRoe(finv,vx,direction,dt)
      use mesh
      use mesh_common
      use properties
      use LRStates

c     
c     This subroutine determines the vector of inviscid fluxes.
      double precision:: finv(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: dt
      integer:: direction,ilo,ihi
      integer:: i,j,k,l
c     
      double precision:: vd(INLO:INHI,nvar)
      double precision:: ud(INLO:INHI,nvar)
      double precision:: vl(INLO:INHI,nvar)
      double precision:: vr(INLO:INHI,nvar)
      double precision:: fx(INLO:INHI,nvar)
c     This subroutine determines the vector of inviscid fluxes.
c     v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c     rho=density, u=x-velocity, v=y-velocity
c     w=z-velocity
c     Bx,By,Bz=x,y,z components of the magnetic field
c     p=thermal pressure
c     
c     vl(i) = Left state after linear representation in cell, i .
c     at the i-th interface
c     vr(i) = Right state after linear representation in cell, i+1 .
c     at the i-th interface
c     
c     
      select case  (direction)
      case(1)
         do k=1,nzlocal,1
            do j=1,nylocal,1
               ilo=ixlo; ihi=ixhi;
               vd(ilo:ihi,:)=vx(ilo:ihi,j,k,:)
c               ud(ilo:ihi,:)=ux(ilo:ihi,j,k,:)
               call ConstructLRStates(vd,vl,vr,ilo,ihi)
               ilo=0; ihi=nxlocal;
c               write(6,*) 'Flux at j=',j
               call InviscidFluxRoeCal(fx,vl,vr,vd,ilo,ihi,direction)
               finv(ilo:ihi,j,k,:)=fx(ilo:ihi,:)
            enddo
         enddo
c         stop
      case(2)
         do k=1,nzlocal,1
            do i=1,nxlocal,1
               ilo=iylo; ihi=iyhi;
               vd(ilo:ihi,:)=vx(i,ilo:ihi,k,:)
c               ud(ilo:ihi,:)=ux(i,ilo:ihi,k,:)
               call RemapToLocalCoords(vd,ilo,ihi,2)
c               call RemapToLocalCoords(ud,ilo,ihi,2)
               call ConstructLRStates(vd,vl,vr,ilo,ihi)
               ilo=0; ihi=nylocal;
               call InviscidFluxRoeCal(fx,vl,vr,vd,ilo,ihi,direction)
               finv(i,ilo:ihi,k,:)=fx(ilo:ihi,:)
            enddo
         enddo
      case(3)
         do j=1,nylocal,1
            do i=1,nxlocal,1
               ilo=izlo; ihi=izhi;
               vd(ilo:ihi,:)=vx(i,j,ilo:ihi,:)
c               ud(ilo:ihi,:)=ux(i,j,ilo:ihi,:)
               call RemapToLocalCoords(vd,ilo,ihi,3)
c               call RemapToLocalCoords(ud,ilo,ihi,3)
               call ConstructLRStates(vd,vl,vr,ilo,ihi)
               ilo=0; ihi=nzlocal;
               call InviscidFluxRoeCal(fx,vl,vr,vd,ilo,ihi,direction)
               finv(i,j,ilo:ihi,:)=fx(ilo:ihi,:)
            enddo
         enddo
c     
      end select
c     
      return
      end
c
c-----------------------------------------------------------------------
      subroutine InviscidFluxRoeCal(finv,vl,vr,vd,ilo,ihi,direction)
      use mesh
      use mesh_common
      use properties
      use slamdas
      double precision:: finv(inlo:inhi,nvar)
      double precision:: vl(inlo:inhi,nvar)
      double precision:: vr(inlo:inhi,nvar)
      double precision:: vd(inlo:inhi,nvar)
      double precision:: sls(inlo:inhi,nvar)
      double precision:: finvl(nvar), finvr(nvar)
      double precision:: finvll(nvar), finvrr(nvar)
      double precision:: rho,u,v,w,bx,by,bz
      double precision:: eint, lam
      double precision:: a2,a1,a0,am1,d2,d1
      integer:: direction,ilo,ihi
      integer:: i,j,k,l
c	 This subroutine determines the vector of inviscid fluxes.
c
c	 v: Vector of variables - rho, u,v,w,Bx,By,Bz,p
c	    rho=density, u=x-velocity, v=y-velocity
c           w=z-velocity
c           Bx,By,Bz= Magnetic Field
c	    p=pressure
c
c	 vl(i) = Left state after linear representation in cell, i .
c		 at the i-th interface
c	 vr(i) = Right state after linear representation in cell, i+1 .
c		 at the i-th interface
c
c
      double precision:: utilde(inlo:inhi,nvar)
c	 utilde = vector of Roe averaged variables.
c	 utilde(1,2,3,4,5,6,7,8) = {rhotilde, utilde, vtilde, wtilde, 
c                           Bxtilde,Bytilde,Bztilde,ptilde}
c	 sls = vector of S |Lamda| S^-1
c	 S: matrix of left eigenvectors
c	 |Lamda|: absolute value of eigenvalues
c	 S^-1: Inverse of S.
c	 
c	 Fluxes finvl=f(u_l)
c		finvr=f(u_r)
c
c
c	 Set the Roe averaged variables depending  on vl,vr
c      write(6,*) 'In Inviscidflux Roe', indx,indy,indz
	 call SetRoeVariables(utilde,vl,vr,ilo,ihi)
c
	 call SetSLS(utilde,sls,vl,vr,ilo,ihi,direction)
c
#ifdef HYBRID
c        Set the centered differences parameters for B field
#ifdef TCD
c        Second order tuned central difference
         d1=0.894d0
         d2=-0.197d0
#elif CD4
c        Fourth order standard central difference
         d1=2.0d0/3.0d0
         d2=-1.0d0/12.0d0
#else
c        Second order standard central difference
         d1=0.5d0
         d2=0.0d0
#endif

c        Compute the corresponding CD weights
         a2=d2
         a1=d1+d2
         a0=a1
         am1=a2
#endif
c
c
         do i=ilo,ihi,1
            rho=vl(i,1)
            u=vl(i,2); v=vl(i,3); w=vl(i,4)
            bx=vl(i,5); by=vl(i,6); bz=vl(i,7)
            finvl(1)=rho*u
            finvl(2)=rho*u*u +vl(i,8)+0.5D0*(by*by+bz*bz-bx*bx)
            finvl(3)=rho*u*v-bx*by
            finvl(4)=rho*u*w-bx*bz
            eint=vl(i,8)/(gamma-1.D0)
            finvl(5)=0.D0
            finvl(6)=u*by-v*bx
            finvl(7)=u*bz-w*bx
            finvl(8)=(0.5D0*rho*(u*u+v*v+w*w)+eint+
     &           vl(i,8)+(bx*bx+by*by+bz*bz))*vl(i,2)-
     &           bx*(u*bx+v*by+w*bz)
c
            rho=vr(i,1)
            u=vr(i,2); v=vr(i,3); w=vr(i,4)
            bx=vr(i,5); by=vr(i,6); bz=vr(i,7)
            finvr(1)=rho*u
            finvr(2)=rho*u*u +vr(i,8)+0.5D0*(by*by+bz*bz-bx*bx)
            finvr(3)=rho*u*v-bx*by
            finvr(4)=rho*u*w-bx*bz
            eint=vr(i,8)/(gamma-1.D0)
            finvr(5)=0.D0
            finvr(6)=u*by-v*bx
            finvr(7)=u*bz-w*bx
            finvr(8)=(0.5D0*rho*(u*u+v*v+w*w)+eint+
     &           vr(i,8)+(bx*bx+by*by+bz*bz))*vr(i,2)-
     &           bx*(u*bx+v*by+w*bz)
c
#ifdef HYBRID
c     Central difference B field
c           flux from left side of face minus one
            rho=vd(i-1,1)
            u=vd(i-1,2); v=vd(i-1,3); w=vd(i-1,4)
            bx=vd(i-1,5); by=vd(i-1,6); bz=vd(i-1,7)
            finvll(5)=0.0d0
            finvll(6)=u*by-v*bx
            finvll(7)=u*bz-w*bx
         
c           flux from left side of face
            rho=vd(i,1)
            u=vd(i,2); v=vd(i,3); w=vd(i,4)
            bx=vd(i,5); by=vd(i,6); bz=vd(i,7)
            finvl(6)=u*by-v*bx
            finvl(7)=u*bz-w*bx
         
c           flux from right side of face
            rho=vd(i+1,1)
            u=vd(i+1,2); v=vd(i+1,3); w=vd(i+1,4)
            bx=vd(i+1,5); by=vd(i+1,6); bz=vd(i+1,7)
            finvr(6)=u*by-v*bx
            finvr(7)=u*bz-w*bx

c           flux from right side of face plus one
            rho=vd(i+2,1)
            u=vd(i+2,2); v=vd(i+2,3); w=vd(i+2,4)
            bx=vd(i+2,5); by=vd(i+2,6); bz=vd(i+2,7)
            finvrr(5)=0.0d0
            finvrr(6)=u*by-v*bx
            finvrr(7)=u*bz-w*bx
#endif
c
c
c     Compute Fluxes
            finv(i,1:4)=0.5D0*(finvl(1:4)+finvr(1:4))
     &           -0.5D0*sls(i,1:4)
            finv(i,5)=0.D0
#ifdef HYBRID
c     Centered-Differences B field
            finv(i,6:7) = a2*finvrr(6:7) + a1*finvr(6:7) 
     &           + a0*finvl(6:7) + am1*finvll(6:7)
#else
c     Upwind B field
            finv(i,6:7)=0.5D0*(finvl(6:7)+finvr(6:7))
     &           -0.5D0*sls(i,6:7)
#endif
            finv(i,nvar)=0.5D0*(finvl(nvar)+finvr(nvar))
     &           -0.5D0*sls(i,nvar)
         enddo
c
c
	return
        end

c$$$c-----------------------------------------------------------------------
c$$$c-----------------------------------------------------------------------
c$$$      subroutine RemapToLocalCoords(vx,ilo,ihi,direction)
c$$$c-----------------------------------------------------------------------
c$$$c     Description: Map to local directions
c$$$c                  x: u,v,w ---> u,v,w
c$$$c                  y: u,v,w ---> w,u,v
c$$$c                  z: u,v,w ---> v,w,u
c$$$c     Input: vx
c$$$c     Output: vx after mapping
c$$$c-----------------------------------------------------------------------
c$$$      use mesh
c$$$      use mesh_common
c$$$      use properties
c$$$c     
c$$$c     Routine to map velocities to local coordinates.
c$$$      double precision:: vx(INLO:INHI,nvar)
c$$$      double precision:: tmp(INLO:INHI,3)
c$$$      integer:: direction
c$$$      integer:: i,j,k,l
c$$$      integer:: ilo, ihi
c$$$c     
c$$$      select case (direction)
c$$$      case (2)
c$$$         tmp(ilo:ihi,1)=vx(ilo:ihi,3)
c$$$         tmp(ilo:ihi,2)=vx(ilo:ihi,4)
c$$$         tmp(ilo:ihi,3)=vx(ilo:ihi,2)
c$$$         vx(ilo:ihi,2:4)=tmp(ilo:ihi,:)
c$$$         tmp(ilo:ihi,1)=vx(ilo:ihi,6)
c$$$         tmp(ilo:ihi,2)=vx(ilo:ihi,7)
c$$$         tmp(ilo:ihi,3)=vx(ilo:ihi,5)
c$$$         vx(ilo:ihi,5:7)=tmp(ilo:ihi,:)
c$$$      case (3)
c$$$         tmp(ilo:ihi,1)=vx(ilo:ihi,4)
c$$$         tmp(ilo:ihi,2)=vx(ilo:ihi,2)
c$$$         tmp(ilo:ihi,3)=vx(ilo:ihi,3)
c$$$         vx(ilo:ihi,2:4)=tmp(ilo:ihi,:)
c$$$         tmp(ilo:ihi,1)=vx(ilo:ihi,7)
c$$$         tmp(ilo:ihi,2)=vx(ilo:ihi,5)
c$$$         tmp(ilo:ihi,3)=vx(ilo:ihi,6)
c$$$         vx(ilo:ihi,5:7)=tmp(ilo:ihi,:)
c$$$      end select
c$$$c     
c$$$      return
c$$$      end
c$$$c
