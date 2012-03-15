c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c Original code:
c     Ravi Samtaney
c     Copyright 2003
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c Revised by:
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------------
c     $Log: inviscidfluxrp.F,v $
c-----------------------------------------------------------------------
c      subroutine InviscidFluxRP(finv,vrp,ux,vx,direction,dt)
      subroutine InviscidFluxRP(finv,vx,direction,dt)
c-----------------------------------------------------------------------
c     Description:  Calculation of Hyperbolic fluxes 
c     Input: ux - conservative variables  (not used)
c            vx - primitive variables
c            direction - 1,2,3 (=x,y,z)
c            dt - time step (not used)
c     Output: finv - Hyperbolic fluxes 
c             vrp - Primitive variables at cell faces in direction
c-----------------------------------------------------------------------
c     04/14/04: Commented out vrp - Riemann problem solution not
c                an argument to inviscidFlux (makes argument list
c                consistent with inviscidFluxCD.
c                Upwinding on rho, u_i, e 
c                Central difference on B.
c                RP Solve changed to return average of normal B component
c                - makes it consistent with 7 wave formulation.
c-----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use LRStates

c     
      double precision :: finv(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c      double precision :: vrp(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: dt
      integer          :: direction
c
c     Local variables
c     ud: conserved var (not used)
c     vd: primitive var 
c     vl, vr: Left right states at cell interface 
c     fx, vm: Flux and solution of RP at cell interface
      integer :: ilo,ihi
      integer :: i,j,k,l
c     
      double precision :: vd(INLO:INHI,nvar)
      double precision :: ud(INLO:INHI,nvar)
      double precision :: vl(INLO:INHI,nvar)
      double precision :: vr(INLO:INHI,nvar)
      double precision :: fx(INLO:INHI,nvar)
      double precision :: vm(INLO:INHI,nvar)
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
               call InviscidFluxRPCal(fx,vm,vl,vr,vd,ilo,ihi,direction)
               finv(ilo:ihi,j,k,:)=fx(ilo:ihi,:)
c               vrp(ilo:ihi,j,k,:)=vm(ilo:ihi,:)
            enddo
         enddo
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
               call InviscidFluxRPCal(fx,vm,vl,vr,vd,ilo,ihi,direction)
               finv(i,ilo:ihi,k,:)=fx(ilo:ihi,:)
c               vrp(i,ilo:ihi,k,:)=vm(ilo:ihi,:)
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
               call InviscidFluxRPCal(fx,vm,vl,vr,vd,ilo,ihi,direction)
               finv(i,j,ilo:ihi,:)=fx(ilo:ihi,:)
c               vrp(i,j,ilo:ihi,:)=vm(ilo:ihi,:)
            enddo
         enddo
c     
      end select
c     
      return
      end
c
c-----------------------------------------------------------------------
      subroutine InviscidFluxRPCal(finv,vm,vl,vr,vd,ilo,ihi,direction)
c-----------------------------------------------------------------------
c     Description: Calculate flux based on linearized solution
c                  to the Riemann problem.
c     Issues: We are computing F=F(vm) where vm = RP(vl,vr) is the
c             solution of the Riemann problem using vl,vr - the
c             left and right states. F=F(vm) exactly iff we solve
c             the nonlinear RP exactly (e.g. as in Godunov for
c             gas dynamics). If we solve a linearized RP (as in our
c             MHD case) F!=F(vm) but we think to leading order this
c             is still OK (refer to discussion with P. Colella)
c     Input: vl, vr: left/right states
c     Output: finv, vm: flux and cell face primitive variables
c-----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use alphas
#ifdef PRIMITIVE
      use evectors_primitive
#else
      use evectors_conservative
#endif
      double precision :: finv(inlo:inhi,nvar)
      double precision :: vl(inlo:inhi,nvar)
      double precision :: vr(inlo:inhi,nvar)
      double precision :: vd(inlo:inhi,nvar)
      double precision :: vm(inlo:inhi,nvar)
c
c     Local variables
c     alamda: eigenvalues based on utilde
c     evl, evr: left/right eigenvectors
c     
      double precision :: rho,u,v,w,bx,by,bz
      double precision :: eint
      double precision :: alamda(INLO:INHI,nvar)
      double precision :: evl(INLO:INHI,nvar,nvar)
      double precision :: evr(INLO:INHI,nvar,nvar)
      double precision :: durl(INLO:INHI,nvar)
      double precision :: alpha(INLO:INHI,nvar)
c
      double precision :: finvl(nvar), finvll(nvar)
      double precision :: finvr(nvar), finvrr(nvar)
      double precision :: a2,a1,a0,am1,d2,d1
c
      integer :: i,j,k,l,m,ilo,ihi,direction
c
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
      double precision :: utilde(inlo:inhi,nvar)
c	 utilde = vector of Roe-averaged variables.
c	 utilde(1,2,3,4,5,6,7,8) = {rhotilde, utilde, vtilde, wtilde, 
c                           Bxtilde,Bytilde,Bztilde,ptilde}
c     Actually only a simple arithmetic average of vl and vr.
      call SetRoeVariables(utilde,vl,vr,ilo,ihi)
c
      alamda = 0.0D0
      durl   = 0.0D0
c
      call SetEigenSystem(utilde,evl,evr,alamda,ilo,ihi,direction)
      durl(ilo:ihi,:)=vr(ilo:ihi,:)-vl(ilo:ihi,:)
      call SetAlphas(utilde,alpha,durl,evl,evr,ilo,ihi)
      call RiemannProblem(evl,evr,alamda,alpha,vl,vr,vm,ilo,ihi)
c
#ifdef HYBRID
c     Set the central difference parameters for B field
#ifdef TCD
c     Second order tuned central difference
      d1=0.894d0
      d2=-0.197d0
#elif CD4
c     Fourth order standard central difference
      d1=2.0d0/3.0d0
      d2=-1.0d0/12.0d0
#else
c     Second order standard central difference
      d1=0.5d0
      d2=0.0d0
#endif
      
c     Compute the corresponding CD weights
      a2=d2
      a1=d1+d2
      a0=a1
      am1=a2
#endif
c     
c
      do i=ilo,ihi,1
#ifdef HYBRID
c        flux from left side of face minus one
         rho=vd(i-1,1)
         u=vd(i-1,2); v=vd(i-1,3); w=vd(i-1,4)
         bx=vd(i-1,5); by=vd(i-1,6); bz=vd(i-1,7)
         finvll(5)=0.0d0
         finvll(6)=u*by-v*bx
         finvll(7)=u*bz-w*bx

c        flux from left side of face
         rho=vd(i,1)
         u=vd(i,2); v=vd(i,3); w=vd(i,4)
         bx=vd(i,5); by=vd(i,6); bz=vd(i,7)
         finvl(5)=0.D0
         finvl(6)=u*by-v*bx
         finvl(7)=u*bz-w*bx
         
c        flux from right side of face
         rho=vd(i+1,1)
         u=vd(i+1,2); v=vd(i+1,3); w=vd(i+1,4)
         bx=vd(i+1,5); by=vd(i+1,6); bz=vd(i+1,7)
         finvr(5)=0.D0
         finvr(6)=u*by-v*bx
         finvr(7)=u*bz-w*bx

c        flux from right side of face plus one
         rho=vd(i+2,1)
         u=vd(i+2,2); v=vd(i+2,3); w=vd(i+2,4)
         bx=vd(i+2,5); by=vd(i+2,6); bz=vd(i+2,7)
         finvrr(5)=0.0d0
         finvrr(6)=u*by-v*bx
         finvrr(7)=u*bz-w*bx

c        central difference B field
         finv(i,5) = 0.D0
         finv(i,6:7) = a2*finvrr(6:7) + a1*finvr(6:7) 
     &        + a0*finvl(6:7) + am1*finvll(6:7)
#endif

c        upwind density, momentum and energy fields
         rho=vm(i,1)
         u=vm(i,2); v=vm(i,3); w=vm(i,4)
         bx=vm(i,5); by=vm(i,6); bz=vm(i,7)
         eint=vm(i,8)/(gamma-1.D0)
         finv(i,1)=rho*u
         finv(i,2)=rho*u*u +vm(i,8)+0.5D0*(by*by+bz*bz-bx*bx)
         finv(i,3)=rho*u*v-bx*by
         finv(i,4)=rho*u*w-bx*bz
#ifndef UPWIND
         finv(i,5)=0.D0
         finv(i,6)=u*by-v*bx
         finv(i,7)=u*bz-w*bx
#endif
         finv(i,8)=(0.5D0*rho*(u*u+v*v+w*w)+eint+
     &        vm(i,8)+(bx*bx+by*by+bz*bz))*vm(i,2)-
     &        bx*(u*bx+v*by+w*bz)

      enddo
c     
c     
      return
      end

c-----------------------------------------------------------------------
      subroutine RiemannProblem(evl,evr,alamda,alpha,vl,vr,vm,ilo,ihi)
c-----------------------------------------------------------------------
c     Description: Linearized Riemann Problem solver
c     Input: evl, evr - left/right eigenvectors (evl is not used)
c            alamda - eigenvalues
c            alpha - evl.durl 
c            vl,vr - left/right states
c     Output: vm 
c-----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      double precision:: vl(inlo:inhi,nvar)
      double precision:: vr(inlo:inhi,nvar)
      double precision:: vm(inlo:inhi,nvar)
      double precision:: alamda(INLO:INHI,nvar)
      double precision:: evl(INLO:INHI,nvar,nvar)
      double precision:: evr(INLO:INHI,nvar,nvar)
      double precision:: alpha(INLO:INHI,nvar)
c
      integer:: direction,ilo,ihi
c     
c     Local variables
      double precision:: vml(inlo:inhi,nvar)
      double precision:: vmr(inlo:inhi,nvar)
      integer:: i,j,k,l
      double precision:: alamdamin=1.0d-1
c
      vmr(ilo:ihi,:)=vr(ilo:ihi,:)
      vml(ilo:ihi,:)=vl(ilo:ihi,:)
      do i=ilo,ihi,1
         do k=1,nvar,1
c            if(alamda(i,k).gt.0.D0) then
            if(alamda(i,k).gt.alamdamin) then
c              upwind direction is to the left
               vmr(i,:)=vmr(i,:)-alpha(i,k)*evr(i,:,k)
c            else if(alamda(i,k).lt.0.D0) then
            else if(alamda(i,k).lt.-alamdamin) then
c              upwind direction is to the right
               vml(i,:)=vml(i,:)+alpha(i,k)*evr(i,:,k)
            else 
c              there is no upwind direction, take an average
               vmr(i,:)=vmr(i,:)-alpha(i,k)*evr(i,:,k)
               vml(i,:)=vml(i,:)+alpha(i,k)*evr(i,:,k)
            endif
         enddo
      enddo
c
c     I like this form because it preserves symmetry
      vm(ilo:ihi,:)=0.5D0*(vml(ilo:ihi,:)+vmr(ilo:ihi,:))
      vm(ilo:ihi,5)=0.5D0*(vml(ilo:ihi,5)+vmr(ilo:ihi,5))
c
      return
      end
c-----------------------------------------------------------------------
