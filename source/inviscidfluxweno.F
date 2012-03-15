C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c
c     Ravi Samtaney
c     Copyright 2001
c     California Institute of Technology
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------------
c     $Log: inviscidfluxweno.F,v $
c-----------------------------------------------------------------------
      subroutine InviscidFlux(finv,ux,vx,csnd,direction,dt)
      use mesh
      use mesh_common
      use properties

c     
c     This subroutine determines the vector of inviscid fluxes.
      double precision:: finv(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: csnd(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision:: vd(INLO:INHI,nvar)
      double precision:: ud(INLO:INHI,nvar)
      double precision:: lami(INLO:INHI,3)
      double precision:: lamj(INLO:INHI,3)
      double precision:: lamk(INLO:INHI,3)
      double precision:: lamjp(INLO:INHI,3)
      double precision:: lamjm(INLO:INHI,3)
      double precision:: lamkp(INLO:INHI,3)
      double precision:: lamkm(INLO:INHI,3)
      double precision:: eta(INLO:INHI)
      double precision:: dt
      integer:: direction
      integer:: i,j,k,l
      integer:: ilo, ihi
c     
c     vx: Vector of primitive variables - rho,u,v,w,p,psi
c     ux: Vector of conserved variables - rho,rhou,rhov,rhow,E,rhopsi
c     finv: Flux of conserved variables
c     csnd: sound speed
c     lami, lamj, lamk: eigenvalues using u,v,w
c     lamjp, lamjm, lamkp, lamkm: eigenvalues used for carbuncle fix.
c     eta: carbuncle fix variable.
c
c
      select case  (direction)
      case(1)
            do k=1,nzlocal,1
               do j=1,nylocal,1
                  lami=0.0; lamj=0.0; lamk=0.0; lamjp=0.0; lamjm=0.0;
                  lamkp=0.0; lamkm=0.0;
                  ilo=ixlo; ihi=ixhi;
c                  write(6,*) 'Before FluxWENO',ilo,ihi,direction,dt
                  vd(ilo:ihi,:)=vx(ilo:ihi,j,k,:)
                  ud(ilo:ihi,:)=ux(ilo:ihi,j,k,:)
c     
                  lami(ilo:ihi,1)=vx(ilo:ihi,j,k,2)
                  lami(ilo:ihi,2)=vx(ilo:ihi,j,k,2)+csnd(ilo:ihi,j,k)
                  lami(ilo:ihi,3)=vx(ilo:ihi,j,k,2)-csnd(ilo:ihi,j,k)
c
                  lamj(ilo:ihi,1)=vx(ilo:ihi,j,k,3)
                  lamj(ilo:ihi,2)=vx(ilo:ihi,j,k,3)+csnd(ilo:ihi,j,k)
                  lamj(ilo:ihi,3)=vx(ilo:ihi,j,k,3)-csnd(ilo:ihi,j,k)
c
                  lamjp(ilo:ihi,1)=vx(ilo:ihi,j+1,k,3)
                  lamjp(ilo:ihi,2)=vx(ilo:ihi,j+1,k,3)
     &                 +csnd(ilo:ihi,j+1,k)
                  lamjp(ilo:ihi,3)=vx(ilo:ihi,j+1,k,3)
     &                 -csnd(ilo:ihi,j+1,k)
c
                  lamjm(ilo:ihi,1)=vx(ilo:ihi,j-1,k,3)
                  lamjm(ilo:ihi,2)=vx(ilo:ihi,j-1,k,3)
     &                 +csnd(ilo:ihi,j-1,k)
                  lamjm(ilo:ihi,3)=vx(ilo:ihi,j-1,k,3)
     &                 -csnd(ilo:ihi,j-1,k)
c
#ifndef TWO_D
                  lamk(ilo:ihi,1)=vx(ilo:ihi,j,k,4)
                  lamk(ilo:ihi,2)=vx(ilo:ihi,j,k,4)+csnd(ilo:ihi,j,k)
                  lamk(ilo:ihi,3)=vx(ilo:ihi,j,k,4)-csnd(ilo:ihi,j,k)
c
                  lamkp(ilo:ihi,1)=vx(ilo:ihi,j,k+1,4)
                  lamkp(ilo:ihi,2)=vx(ilo:ihi,j,k+1,4)
     &                 +csnd(ilo:ihi,j,k+1)
                  lamkp(ilo:ihi,3)=vx(ilo:ihi,j,k+1,4)
     &                 -csnd(ilo:ihi,j,k+1)
c
                  lamkm(ilo:ihi,1)=vx(ilo:ihi,j,k-1,4)
                  lamkm(ilo:ihi,2)=vx(ilo:ihi,j,k-1,4)
     &                 +csnd(ilo:ihi,j,k-1)
                  lamkm(ilo:ihi,3)=vx(ilo:ihi,j,k-1,4)
     &                 -csnd(ilo:ihi,j,k-1)
#endif
c
                  ilo=0; ihi=nxlocal;
                  Call CalculateETA(lami,lamj,lamk,
     &                 lamjp,lamjm,lamkp,lamkm,
     &                 eta,ilo,ihi)
                  Call FluxWENO(vd,ud,lami,eta,ilo,ihi)
c                  write(6,*) 'After FluxWENO',ilo,ihi,direction,dt
                  finv(ilo:ihi,j,k,:)=vd(ilo:ihi,:)
               enddo
            enddo
c
      case(2)
            do k=1,nzlocal,1
               do i=1,nxlocal,1
            	ilo=iylo; ihi=iyhi;
                  lami=0.0; lamj=0.0; lamk=0.0; lamjp=0.0; lamjm=0.0;
                  lamkp=0.0; lamkm=0.0;
                  vd(ilo:ihi,:)=vx(i,ilo:ihi,k,:)
                  ud(ilo:ihi,:)=ux(i,ilo:ihi,k,:)
                  call RemapToLocalCoords(vd,ilo,ihi,2)
                  call RemapToLocalCoords(ud,ilo,ihi,2)
c
c
                  lami(ilo:ihi,1)=vx(i,ilo:ihi,k,3)
                  lami(ilo:ihi,2)=vx(i,ilo:ihi,k,3)+csnd(i,ilo:ihi,k)
                  lami(ilo:ihi,3)=vx(i,ilo:ihi,k,3)-csnd(i,ilo:ihi,k)
c
                  lamk(ilo:ihi,1)=vx(i,ilo:ihi,k,2)
                  lamk(ilo:ihi,2)=vx(i,ilo:ihi,k,2)+csnd(i,ilo:ihi,k)
                  lamk(ilo:ihi,3)=vx(i,ilo:ihi,k,2)-csnd(i,ilo:ihi,k)
c
#ifndef TWO_D
c
                  lamj(ilo:ihi,1)=vx(i,ilo:ihi,k,4)
                  lamj(ilo:ihi,2)=vx(i,ilo:ihi,k,4)+csnd(i,ilo:ihi,k)
                  lamj(ilo:ihi,3)=vx(i,ilo:ihi,k,4)-csnd(i,ilo:ihi,k)
                  lamjp(ilo:ihi,1)=vx(i,ilo:ihi,k+1,4)
                  lamjp(ilo:ihi,2)=vx(i,ilo:ihi,k+1,4)
     &                 +csnd(i,ilo:ihi,k+1)
                  lamjp(ilo:ihi,3)=vx(i,ilo:ihi,k+1,4)
     &                 -csnd(i,ilo:ihi,k+1)
c
                  lamjm(ilo:ihi,1)=vx(i,ilo:ihi,k-1,4)
                  lamjm(ilo:ihi,2)=vx(i,ilo:ihi,k-1,4)
     &                 +csnd(i,ilo:ihi,k-1)
                  lamjm(ilo:ihi,3)=vx(i,ilo:ihi,k-1,4)
     &                 -csnd(i,ilo:ihi,k-1)
#endif
c
                  lamkp(ilo:ihi,1)=vx(i+1,ilo:ihi,k,2)
                  lamkp(ilo:ihi,2)=vx(i+1,ilo:ihi,k,2)
     &                 +csnd(i+1,ilo:ihi,k)
                  lamkp(ilo:ihi,3)=vx(i+1,ilo:ihi,k,2)
     &                 -csnd(i+1,ilo:ihi,k)
c
                  lamkm(ilo:ihi,1)=vx(i-1,ilo:ihi,k,2)
                  lamkm(ilo:ihi,2)=vx(i-1,ilo:ihi,k,2)
     &                 +csnd(i-1,ilo:ihi,k)
                  lamkm(ilo:ihi,3)=vx(i-1,ilo:ihi,k,2)
     &                 -csnd(i-1,ilo:ihi,k)
c
                  ilo=0; ihi=nylocal;
                  Call CalculateETA(lami,lamj,lamk,
     &                 lamjp,lamjm,lamkp,lamkm,
     &                 eta,ilo,ihi)
                  Call FluxWENO(vd,ud,lami,eta,ilo,ihi)
                  finv(i,ilo:ihi,k,:)=vd(ilo:ihi,:)
               enddo
            enddo
c
      case(3)
            do j=1,nylocal,1
               do i=1,nxlocal,1
                  ilo=izlo; ihi=izhi;
                  vd(ilo:ihi,:)=vx(i,j,ilo:ihi,:)
                  ud(ilo:ihi,:)=ux(i,j,ilo:ihi,:)
                  call RemapToLocalCoords(vd,ilo,ihi,3)
                  call RemapToLocalCoords(ud,ilo,ihi,3)
c
                  lami(ilo:ihi,1)=vx(i,j,ilo:ihi,4)
                  lami(ilo:ihi,2)=vx(i,j,ilo:ihi,4)+csnd(i,j,ilo:ihi)
                  lami(ilo:ihi,3)=vx(i,j,ilo:ihi,4)-csnd(i,j,ilo:ihi)
c
                  lamj(ilo:ihi,1)=vx(i,j,ilo:ihi,2)
                  lamj(ilo:ihi,2)=vx(i,j,ilo:ihi,2)+csnd(i,j,ilo:ihi)
                  lamj(ilo:ihi,3)=vx(i,j,ilo:ihi,2)-csnd(i,j,ilo:ihi)
c
                  lamk(ilo:ihi,1)=vx(i,j,ilo:ihi,3)
                  lamk(ilo:ihi,2)=vx(i,j,ilo:ihi,3)+csnd(i,j,ilo:ihi)
                  lamk(ilo:ihi,3)=vx(i,j,ilo:ihi,3)-csnd(i,j,ilo:ihi)
c
                  lamjp(ilo:ihi,1)=vx(i+1,j,ilo:ihi,2)
                  lamjp(ilo:ihi,2)=vx(i+1,j,ilo:ihi,2)
     &                 +csnd(i+1,j,ilo:ihi)
                  lamjp(ilo:ihi,3)=vx(i+1,j,ilo:ihi,2)
     &                 -csnd(i+1,j,ilo:ihi)
c
                  lamjm(ilo:ihi,1)=vx(i-1,j,ilo:ihi,2)
                  lamjm(ilo:ihi,2)=vx(i-1,j,ilo:ihi,2)
     &                 +csnd(i-1,j,ilo:ihi)
                  lamjm(ilo:ihi,3)=vx(i-1,j,ilo:ihi,2)
     &                 -csnd(i-1,j,ilo:ihi)
c
                  lamkp(ilo:ihi,1)=vx(i,j+1,ilo:ihi,3)
                  lamkp(ilo:ihi,2)=vx(i,j+1,ilo:ihi,3)
     &                 +csnd(i,j+1,ilo:ihi)
                  lamkp(ilo:ihi,3)=vx(i,j+1,ilo:ihi,3)
     &                 -csnd(i,j+1,ilo:ihi)
c
                  lamkm(ilo:ihi,1)=vx(i,j-1,ilo:ihi,3)
                  lamkm(ilo:ihi,2)=vx(i,j-1,ilo:ihi,3)
     &                 +csnd(i,j-1,ilo:ihi)
                  lamkm(ilo:ihi,3)=vx(i,j-1,ilo:ihi,3)
     &                 -csnd(i,j-1,ilo:ihi)
c
                  ilo=0; ihi=nzlocal;
                  Call CalculateETA(lami,lamj,lamk,
     &                 lamjp,lamjm,lamkp,lamkm,
     &                 eta,ilo,ihi)
                  Call FluxWENO(vd,ud,lami,eta,ilo,ihi)
                  finv(i,j,ilo:ihi,:)=vd(ilo:ihi,:)
               enddo
            enddo
      end select
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine FluxWENO(vx,ux,lami,eta,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: vx(inlo:inhi,nvar)
      double precision:: ux(inlo:inhi,nvar)
      double precision:: eta(inlo:inhi)
      double precision:: alpha(inlo:inhi,nvar)
      double precision:: lami(inlo:inhi,3)
      double precision:: vroe(INLO:INHI,nvar+3)
      double precision:: fx(inlo:inhi,nvar)
      double precision:: fxp(inlo:inhi,nvar)
      double precision:: fxm(inlo:inhi,nvar)
      double precision:: lc(inlo:inhi,nvar,nvar)
      double precision:: lamda(inlo:inhi,nvar)
      double precision:: rc(inlo:inhi,nvar,nvar)
      double precision:: fs(-enoOrder+1:enoOrder,inlo:inhi,nvar)
      double precision:: us(-enoOrder+1:enoOrder,inlo:inhi,nvar)
      double precision:: fps(-enoOrder+1:enoOrder-1,inlo:inhi,nvar)
      double precision:: fms(-enoOrder+1:enoOrder-1,inlo:inhi,nvar)
      double precision:: isk(0:enoOrder-1,inlo:inhi,nvar)
      double precision:: omega(0:enoOrder-1,inlo:inhi,nvar)
      double precision:: qk(0:enoOrder-1,inlo:inhi,nvar)
      integer:: i,k,l
c
      fxp=0.D0; fxm=0.D0
c
      call Flux(vx,fx,ilo-nghost+1,ihi+nghost)
c
c
      call RoeVariables(vx,vroe,ilo,ihi)
c
      call EigenSystem(vroe,lc,rc,lamda,ilo,ihi)
c
c      call CheckProjection(fx,lc,rc,lamda,ilo,ihi)
c
      call CalculateAlpha(alpha,lamda,lami,eta,ilo,ihi)
      call ProjectFluxStencils(fs,fx,us,ux,lc,ilo,ihi)
c
      call SignedFluxes(fps,fms,fs,us,alpha,ilo,ihi)
c
      call SmoothnessMeasure(fms,isk,ilo,ihi)
c
      call CalculateWeights(omega,isk,ilo,ihi)
c
      call PolynomialFit(qk,fms,ilo,ihi)
c
      do i=ilo,ihi,1
         do l=1,nvar,1
            fxm(i,l)=sum(qk(0:enoOrder-1,i,l)
     &           *omega(0:enoOrder-1,i,l))
         enddo
      enddo
c
      qk=0.D0; omega=0.D0; isk=0.D0;
      call SmoothnessMeasure(fps,isk,ilo,ihi)
c
      call CalculateWeights(omega,isk,ilo,ihi)
c
      call PolynomialFit(qk,fps,ilo,ihi)
c

      do i=ilo,ihi,1
         do l=1,nvar,1
           fxp(i,l)=sum(qk(0:enoOrder-1,i,l)*
     &           omega(0:enoOrder-1,i,l))
         enddo
      enddo
c
c
c
      vx=fxp+fxm
      call ProjectFluxPrimitive(vx,rc,ilo,ihi)
c
c
      return
      end

c-----------------------------------------------------------------------
      subroutine  SignedFluxes(fps,fms,fs,us,alpha,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: alpha(inlo:inhi,nvar)
      double precision:: fs(-enoOrder+1:enoOrder,inlo:inhi,nvar)
      double precision:: us(-enoOrder+1:enoOrder,inlo:inhi,nvar)
      double precision:: fps(-enoOrder+1:enoOrder-1,inlo:inhi,nvar)
      double precision:: fms(-enoOrder+1:enoOrder-1,inlo:inhi,nvar)
c
      integer:: i,j,l
c
      do j=1,nvar,1
         do i=ilo,ihi,1
      do l=-enoOrder+1, enoOrder-1
c         fps(l,:,:) = 0.5*(fs(l,:,:) + us(l,:,:))
c         fps(l,:,:) = 0.5*(fs(l,:,:) + alpha(:,:)*us(l,:,:))
         fps(l,i,j) = 0.5*(fs(l,i,j) + alpha(i,j)*us(l,i,j))
      enddo
c
      do l=-enoOrder+1, enoOrder-1
         fms(l,i,j) = 0.5*(fs(-l+1,i,j) - alpha(i,j)*us(-l+1,i,j))
c         fms(l,:,:) = 0.5*(fs(-l+1,:,:) - alpha(:,:)*us(-l+1,:,:))
c         fms(l,:,:) = 0.5*(fs(-l+1,:,:) - us(-l+1,:,:))
      enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ProjectFluxStencils(fs,fx,us,vx,lc,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: fx(inlo:inhi,nvar)
      double precision:: fs(-enoOrder+1:enoOrder,inlo:inhi,nvar)
      double precision:: vx(inlo:inhi,nvar)
      double precision:: us(-enoOrder+1:enoOrder,inlo:inhi,nvar)
      double precision:: lc(inlo:inhi,nvar,nvar)
      double precision:: lamda(inlo:inhi,nvar)
      double precision:: alpha(inlo:inhi,nvar)
      integer:: i,k,l

      do i=ilo,ihi,1
         do l=1,nvar,1
            do k=-enoOrder+1,enoOrder,1
               fs(k,i,l)=lc(i,l,1)*fx(i+k,1)+
     &              lc(i,l,2)*fx(i+k,2)+
     &              lc(i,l,3)*fx(i+k,3)+
     &              lc(i,l,4)*fx(i+k,4)+
     &              lc(i,l,5)*fx(i+k,5)+
     &              lc(i,l,6)*fx(i+k,6)
c
               us(k,i,l)=lc(i,l,1)*vx(i+k,1)+
     &              lc(i,l,2)*vx(i+k,2)+
     &              lc(i,l,3)*vx(i+k,3)+
     &              lc(i,l,4)*vx(i+k,4)+
     &              lc(i,l,5)*vx(i+k,5)+
     &              lc(i,l,6)*vx(i+k,6)
c
c               us(k,i,l)=lc(i,l,1)*vx(i+k,1)*alpha(i,1)+
c     &              lc(i,l,2)*vx(i+k,2)*alpha(i,2)+
c     &              lc(i,l,3)*vx(i+k,3)*alpha(i,3)+
c     &              lc(i,l,4)*vx(i+k,4)*alpha(i,4)+
c     &              lc(i,l,5)*vx(i+k,5)*alpha(i,5)+
c     &              lc(i,l,6)*vx(i+k,6)*alpha(i,6)
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine SmoothnessMeasure(fs,isk,ilo,ihi)
      use mesh
      use mesh_common
      use wenoCoeffs
      integer:: ilo, ihi
      double precision,intent(in)::
     &     fs(-enoOrder+1:enoOrder-1,inlo:inhi,nvar)
      double precision:: isk(0:enoOrder-1,inlo:inhi,nvar)
      double precision:: tmpsum(0:enoOrder-1,enoOrder-1)
      integer:: i,j,k,l,m
c
      isk=0.D0
c
      do i=ilo,ihi,1
         do j=1,nvar,1
c
c     Compute inner sum for m=1,enoOrder-1
            do m=1,enoOrder-1,1
               do k=0,enoOrder-1,1
                  tmpsum(k,m)=(sum(dw(k,m,0:enoOrder-1)*
     &                 fs(-enoOrder+1+k:k,i,j)))**2
               enddo
            enddo
c
            do k=0,enoOrder-1,1
               isk(k,i,j)=sum(tmpsum(k,1:enoOrder-1))
            enddo
c
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine CalculateWeights(omega,isk,ilo,ihi)
      use mesh
      use mesh_common
      use wenoCoeffs
      integer:: ilo, ihi
      double precision:: isk(0:enoOrder-1,inlo:inhi,nvar)
      double precision:: omega(0:enoOrder-1,inlo:inhi,nvar)
      double precision:: alpha(0:enoOrder-1)
      integer:: i,j,k,l
      double precision, parameter:: eps=1.D-08
      integer, parameter:: pwr=2
c

      do i=ilo,ihi,1
         do j=1,nvar,1
c            isk(enoOrder,i,j)=maxval(isk(0:enoOrder,i,j))
            alpha(0:enoOrder-1) = cw(0:enoOrder-1)/
     &           (eps+isk(0:enoOrder-1,i,j))**pwr
c     
            omega(0:enoOrder-1,i,j)=
     &           alpha(0:enoOrder-1)/sum(alpha(0:enoOrder-1))
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine PolynomialFit(qk,fs,ilo,ihi)
      use mesh
      use mesh_common
      use wenoCoeffs
      integer:: ilo, ihi
      double precision,intent(in)::
     &     fs(-enoOrder+1:enoOrder-1,inlo:inhi,nvar)
      double precision:: qk(0:enoOrder-1,inlo:inhi,nvar)
      integer:: i,j,k,l
c

      do i=ilo,ihi,1
         do j=1,nvar,1
            do k=0,enoOrder-1
               qk(k,i,j)=sum(
     &           aw(k,0:enoOrder-1)*fs(-enoOrder+1+k:k,i,j))
            enddo
         enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine ProjectFluxPrimitive(fx,rc,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: fx(inlo:inhi,nvar)
      double precision:: rc(inlo:inhi,nvar,nvar)
      double precision:: tmp(nvar)
      integer:: i,k,l

      do i=ilo,ihi,1
         do l=1,nvar,1
            tmp(l)=rc(i,l,1)*fx(i,1)+
     &           rc(i,l,2)*fx(i,2)+
     &           rc(i,l,3)*fx(i,3)+
     &           rc(i,l,4)*fx(i,4)+
     &           rc(i,l,5)*fx(i,5)+
     &           rc(i,l,6)*fx(i,6)
         enddo
         fx(i,1:nvar)=tmp
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine EigenSystem(vx,lc,rc,lamda,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: lamda(inlo:inhi,nvar)
      double precision:: vx(inlo:inhi,nvar+3)
      double precision:: lc(inlo:inhi,nvar,nvar)
      double precision:: rc(inlo:inhi,nvar,nvar)
      double precision:: prod(nvar,nvar)
      integer:: i,j,k,l
      double precision:: rho, u, v, w, gamm, h, a,aisq,qsq
c
c     vx: rho,u,v,w,p,rho1,a,h,gamma 
c     vx: 1   2 3 4 5  6   7 8  9
c      csndsq(ilo:ihi)=(gamma(ilo:ihi)*vx(ilo:ihi,3)/vx(ilo:ihi,1))
c      csnd(ilo:ihi)=sqrt(csndsq(ilo:ihi))
c      enthalpy(ilo:ihi)=gamma*(gamma-1.0)*vx(ilo:ihi,3)/vx(ilo:ihi,1)
c
      lamda(ilo:ihi,1)=vx(ilo:ihi,2)
      lamda(ilo:ihi,2)=vx(ilo:ihi,2)
      lamda(ilo:ihi,3)=vx(ilo:ihi,2)
      lamda(ilo:ihi,4)=vx(ilo:ihi,2)+vx(ilo:ihi,7)
      lamda(ilo:ihi,5)=vx(ilo:ihi,2)-vx(ilo:ihi,7)
      lamda(ilo:ihi,6)=vx(ilo:ihi,2)
c
      do i=ilo,ihi,1
         rho=vx(i,1)
         u=vx(i,2); v=vx(i,3); w=vx(i,4)
c         a=vx(i,7)
c         h=vx(i,8)
         gamm=vx(i,9)-1.D0
         qsq=(u*u+v*v+w*w)
         a=sqrt(vx(i,9)*vx(i,5)/vx(i,1))
         h=vx(i,9)/gamm*vx(i,5)/vx(i,1)+0.5*qsq
         vx(i,8)=h; vx(i,7)=a
         aisq=1.D0/a/a
c
         rc(i,1,1)=0.D0
         rc(i,2,1)=0.D0
         rc(i,3,1)=0.D0
         rc(i,4,1)=1.D0
         rc(i,5,1)=w
         rc(i,6,1)=0.D0
c
         rc(i,1,2)=0.D0
         rc(i,2,2)=0.D0
         rc(i,3,2)=1.D0
         rc(i,4,2)=0.D0
         rc(i,5,2)=v
         rc(i,6,2)=0.D0
c
         rc(i,1,3)=1.D0
         rc(i,2,3)=u
         rc(i,3,3)=v
         rc(i,4,3)=w
         rc(i,5,3)=0.5D0*qsq
         rc(i,6,3)=0.D0
c
         rc(i,1,4)=1.D0
         rc(i,2,4)=u+a
         rc(i,3,4)=v
         rc(i,4,4)=w
         rc(i,5,4)=h+u*a
         rc(i,6,4)=0.D0
c
         rc(i,1,5)=1.D0
         rc(i,2,5)=u-a
         rc(i,3,5)=v
         rc(i,4,5)=w
         rc(i,5,5)=h-u*a
         rc(i,6,5)=0.D0
c
         rc(i,1,6)=0.D0
         rc(i,2,6)=0.D0
         rc(i,3,6)=0.D0
         rc(i,4,6)=0.D0
         rc(i,5,6)=0.D0
         rc(i,6,6)=1.D0
c
c
         lc(i,1,1)=-w
         lc(i,1,2)=0.D0
         lc(i,1,3)=0.D0
         lc(i,1,4)=1.D0
         lc(i,1,5)=0.D0
         lc(i,1,6)=0.D0
c
         lc(i,2,1)=-v
         lc(i,2,2)=0.D0
         lc(i,2,3)=1.D0
         lc(i,2,4)=0.D0
         lc(i,2,5)=0.D0
         lc(i,2,6)=0.D0
c
         lc(i,3,1)=1.D0-0.5D0*gamm*qsq*aisq
         lc(i,3,2)=gamm*u*aisq
         lc(i,3,3)=gamm*v*aisq
         lc(i,3,4)=gamm*w*aisq
         lc(i,3,5)=-gamm*aisq
         lc(i,3,6)=0.D0
c
         lc(i,4,1)=0.25D0*gamm*qsq*aisq-0.5D0*u/a
         lc(i,4,2)=0.5D0/a-0.5D0*gamm*u*aisq
         lc(i,4,3)=-0.5D0*gamm*v*aisq
         lc(i,4,4)=-0.5D0*gamm*w*aisq
         lc(i,4,5)=0.5D0*gamm*aisq
         lc(i,4,6)=0.D0
c
         lc(i,5,1)=0.25D0*gamm*qsq*aisq+0.5D0*u/a
         lc(i,5,2)=-0.5D0/a-0.5D0*gamm*u*aisq
         lc(i,5,3)=-0.5D0*gamm*v*aisq
         lc(i,5,4)=-0.5D0*gamm*w*aisq
         lc(i,5,5)=0.5D0*gamm*aisq
         lc(i,5,6)=0.D0
c
         lc(i,6,1)=0.D0
         lc(i,6,2)=0.D0
         lc(i,6,3)=0.D0
         lc(i,6,4)=0.D0
         lc(i,6,5)=0.D0
         lc(i,6,6)=1.D0
c
      enddo
c
      return

      end

c-----------------------------------------------------------------------
      subroutine RemapToLocalCoords(vx,ilo,ihi,direction)
      use mesh
      use mesh_common
      use properties
c     
c     Routine to map velocities to local coordinates.
      double precision:: vx(INLO:INHI,nvar)
      double precision:: tmp(INLO:INHI,3)
      integer:: direction
      integer:: i,j,k,l
      integer:: ilo, ihi
c
      select case (direction)
      case (2)
         tmp(ilo:ihi,1)=vx(ilo:ihi,3)
         tmp(ilo:ihi,2)=vx(ilo:ihi,4)
         tmp(ilo:ihi,3)=vx(ilo:ihi,2)
         vx(ilo:ihi,2:4)=tmp(ilo:ihi,:)
      case (3)
         tmp(ilo:ihi,1)=vx(ilo:ihi,4)
         tmp(ilo:ihi,2)=vx(ilo:ihi,2)
         tmp(ilo:ihi,3)=vx(ilo:ihi,3)
         vx(ilo:ihi,2:4)=tmp(ilo:ihi,:)
      end select
c
      return
      end

c-----------------------------------------------------------------------
      subroutine ConservedToPrimitive(u,v,csnd)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: u(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: v(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: csnd(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      integer:: i,j,k,l
      double precision:: gamma, rho_inv, vi,vj,vk
      double precision, external:: gammafun


#ifdef ONE_D
      do k=1,1
         do j=1,1
#else
#ifdef TWO_D
      do k=1,1
#else
      do k=1-NGHOST,NZlocal+NGHOST,1
#endif
         do j=1-NGHOST,NYlocal+NGHOST,1
#endif
            do i=1-NGHOST, NXlocal+NGHOST,1
               rho_inv=1.0D0/u(i,j,k,1)
               v(i,j,k,1)=u(i,j,k,1)
               v(i,j,k,2)=u(i,j,k,2)*rho_inv
               v(i,j,k,3)=u(i,j,k,3)*rho_inv
               v(i,j,k,4)=u(i,j,k,4)*rho_inv
               v(i,j,k,6)=u(i,j,k,6)*rho_inv
               vi=v(i,j,k,2)
               vj=v(i,j,k,3)
               vk=v(i,j,k,4)
               gamma=gammafun(v(i,j,k,6))
               v(i,j,k,5)=(gamma-1)*(u(i,j,k,5)-0.5D0*v(i,j,k,1)
     &              *(vi*vi+vj*vj+vk*vk))
               csnd(i,j,k)=sqrt(gamma*v(i,j,k,5)/v(i,j,k,1))
            enddo
         enddo
      enddo
c
      return
      end

c
c-----------------------------------------------------------------------
      subroutine CalculateETA(lami,lamj,lamk,lamjp,lamjm,lamkp,lamkm,
     &                 eta,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: lami(INLO:INHI,3)
      double precision:: lamj(INLO:INHI,3)
      double precision:: lamk(INLO:INHI,3)
      double precision:: lamjp(INLO:INHI,3)
      double precision:: lamjm(INLO:INHI,3)
      double precision:: lamkp(INLO:INHI,3)
      double precision:: lamkm(INLO:INHI,3)
      double precision:: eta(INLO:INHI)
      double precision:: etaH(9)
      integer:: i
c
      eta=0.D0
      do i=ilo,ihi,1
         etaH(1)=maxval(abs(lami(i+1,:)-lami(i,:)))
         etaH(2)=maxval(abs(lamjp(i,:)-lamj(i,:)))
         etaH(3)=maxval(abs(lamj(i,:)-lamjm(i,:)))
         etaH(4)=maxval(abs(lamjp(i+1,:)-lamj(i+1,:)))
         etaH(5)=maxval(abs(lamj(i+1,:)-lamjm(i+1,:)))
         etaH(6)=maxval(abs(lamkp(i,:)-lamk(i,:)))
         etaH(7)=maxval(abs(lamk(i,:)-lamkm(i,:)))
         etaH(8)=maxval(abs(lamkp(i+1,:)-lamk(i+1,:)))
         etaH(9)=maxval(abs(lamk(i+1,:)-lamkm(i+1,:)))
c
         eta(i)=0.5*maxval(etaH)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine CalculateAlpha(alpha,lamda,lami,eta,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: eta(inlo:inhi)
      double precision:: alpha(inlo:inhi,nvar)
      double precision:: lami(inlo:inhi,3)
      double precision:: lamw(inlo:inhi,nvar)
      double precision:: lamda(inlo:inhi,nvar)
c
      lamw(:,1)=lami(:,1)
      lamw(:,2)=lami(:,1)
      lamw(:,3)=lami(:,1)
      lamw(:,4)=lami(:,2)
      lamw(:,5)=lami(:,3)
      lamw(:,6)=lami(:,1)
      lamw=abs(lamw)
c
      do j=1,nvar,1
         do i=ilo,ihi
            alpha(i,j)=max(dabs(lamda(i,j)),
     &           maxval(lamw(i-enoOrder+1:i+enoOrder,j)),eta(i))
         enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine RoeVariables(vx,vroe,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision,intent(in):: vx(inlo:inhi,nvar)
      double precision:: vroe(inlo:inhi,nvar+3)
      double precision:: gammal, gammar
      double precision:: csndl, csndr
      double precision:: hl, hr
      double precision:: rhol,rhor,qsq
      double precision,external:: gammafun
      integer:: i,j
c
      do i=ilo,ihi
         rhol=sqrt(vx(i,1)); rhor=sqrt(vx(i+1,1))
c
         gammal=gammafun(vx(i,6))
         gammar=gammafun(vx(i+1,6))
c
         csndl=sqrt(gammal*vx(i,5)/vx(i,1))
         csndr=sqrt(gammar*vx(i+1,5)/vx(i+1,1))
c
         qsq=(vx(i,2)**2+vx(i,3)**2+vx(i,4)**2)
         hl=gammal/(gammal-1.0)*vx(i,5)/vx(i,1)+0.5*qsq
c
         qsq=(vx(i+1,2)**2+vx(i+1,3)**2+vx(i+1,4)**2)
         hr=gammar/(gammar-1.0)*vx(i+1,5)/vx(i+1,1)+0.5*qsq
c
         vroe(i,1:nvar)=(vx(i,1:nvar)*rhol+vx(i+1,1:nvar)*rhor)
     &        /(rhol+rhor)
         vroe(i,nvar+1)=(csndl*rhol+csndr*rhor)
     &        /(rhol+rhor)
         vroe(i,nvar+2)=(hl*rhol+hr*rhor)
     &        /(rhol+rhor)
         vroe(i,nvar+3)=(gammal*rhol+gammar*rhor)
     &        /(rhol+rhor)
c
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine Flux(vx,fx,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: vx(inlo:inhi,nvar)
      double precision:: fx(inlo:inhi,nvar)
      double precision:: gamma(inlo:inhi)
      double precision,external:: gammafun
      integer:: i,j
c
      do i=ilo,ihi
         gamma(i)=gammafun(vx(i,6))
      enddo
c
      fx(ilo:ihi,1)=vx(ilo:ihi,1)*vx(ilo:ihi,2)

      fx(ilo:ihi,2)=vx(ilo:ihi,1)*vx(ilo:ihi,2)**2+vx(ilo:ihi,5)

      fx(ilo:ihi,3)=fx(ilo:ihi,1)*vx(ilo:ihi,3)

      fx(ilo:ihi,4)=fx(ilo:ihi,1)*vx(ilo:ihi,4)

      fx(ilo:ihi,5)=(gamma(ilo:ihi)/(gamma(ilo:ihi)-1.0)*
     &     vx(ilo:ihi,5)+0.5*vx(ilo:ihi,1)*(vx(ilo:ihi,2)**2+
     &     vx(ilo:ihi,3)**2+vx(ilo:ihi,4)**2))*vx(ilo:ihi,2)

      fx(ilo:ihi,6)=fx(ilo:ihi,1)*vx(ilo:ihi,6)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine CheckProjection(fx,lc,rc,lamda,ilo,ihi)
      use mesh
      use mesh_common
      use properties
      integer:: ilo, ihi
      double precision:: fx(inlo:inhi,nvar)
      double precision:: lc(inlo:inhi,nvar,nvar)
      double precision:: rc(inlo:inhi,nvar,nvar)
      double precision:: lamda(inlo:inhi,nvar)
      double precision:: tmp(nvar), tmp2(nvar)
      integer:: i,k,l
c
         do i=ilo,ihi,1
            tmp=0.D0
            tmp2=0.D0
            do k=1,nvar,1
               do l=1,nvar,1
                  tmp(k)=tmp(k)+lc(i,k,l)*fx(i,l)
               enddo
            enddo
            do k=1,nvar,1
               do l=1,nvar,1
                  tmp2(k)=tmp2(k)+rc(i,k,l)*tmp(l)
               enddo
            enddo
               write(51,*) 'FX',i,fx(i,:)
               write(51,*) 'LF',i,tmp(:)
               write(51,*) 'RLF',i,tmp2(:)
         enddo
      
      return
      end
