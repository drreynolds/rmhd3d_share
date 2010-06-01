C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2004
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------
c     $Log: initwave3.F,v $
c-----------------------------------------------------------------
c Subroutine to set up the initial conditions
c---------------------------------------------------------------------------
      module waveparams
      save
      integer  pdir,kratio,waveNumber
      real*8 pertAmplitude 
      real*8 alpha
      real*8 perturbation
      end module waveparams
      
c
      subroutine InitialConditions(ux,phi)
c----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use iounits
      use waveparams
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
c
c=======Declarations=========
	integer ilo,ihi,jlo,jhi,klo,khi
	integer iulo,iuhi,julo,juhi,kulo,kuhi
c
	integer  i,j,k,l,m
        real*8 dxi,dyi
c
        real*8 x(3), kwave(3)
        real*8 v(8), v0(8), w0(8), w(8)
        real*8 evl(8,8), evr(8,8),lambda(8)
        real*8 theta, beta
        real*8 un,ut,bn,bt, uz,bz
        real*8 bx0,by0
        real*8 pi
c
c
	namelist/molwts/wmol
	namelist/gammas/gamma
	namelist/diffusiveproperties/mu,eta,etaFactor,kappa
        namelist/waveinit/pdir,alpha,waveNumber,kratio,pertAmplitude
c
        UX=0.d0
C
        open(16,file='prop.inp',form='formatted')
	read(16,molwts)
	read(16,gammas) 
	read(16,diffusiveproperties) 
        close(16)
        if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
           write(6,'(A,6(A,es12.5))') '\nPlasma properties:',
     &          '\n       wmol =',wmol,
     &          '\n      gamma =',gamma,
     &          '\n         mu =',mu,
     &          '\n        eta =',eta,
     &          '\n  etaFactor =',etaFactor,
     &          '\n      kappa =',kappa
        end if
c
        open(17,file='waveinit.inp',form='formatted')
        read(17,waveinit)
        close(17)
        if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
           write(6,'(A,A,i3,A,es10.2,A,i3,A,i3,A,es10.2)') 
     &          '\nWave Properties:',
     &          '\n           pdir =',pdir,
     &          '\n          alpha =',alpha,
     &          '\n         kratio =',kratio,
     &          '\n     waveNumber =',waveNumber,
     &          '\n  pertAmplitude =',pertAmplitude
        end if
c
        v(1)=1.D0
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
        v(5)=1.D0
        v(6)=0.D0
        v(7)=0.D0
        v(8)=0.01D0
c
        pi=4.D0*datan(1.D0)
        kwave(1)=pi*waveNumber
        kwave(2)=kratio*kwave(1)
c        kwave(2)=pi*waveNumber
c        kwave(1)=kratio*kwave(2)
        kwave(3)=0.D0
c        kwave(2)=kwave(1)
c        kwave(3)=kwave(1)
c        theta=datan(2.D0)
c     Wave vector is at an angle theta to x
        theta=datan(kwave(2)/kwave(1))
c
c     B0 is aligned at an angle alpha to x
c        alpha=-pi/4.D0
c        alpha=pi/180.D0*(44.D0)
c        alpha=0.D0
        alpha=pi/180.D0*alpha
        bx0=dcos(alpha)
        by0=dsin(alpha)
c
        beta=alpha+theta
c
c        kwave(1)=kwave(1)*dcos(beta)
c        kwave(2)=kwave(1)*dsin(beta)
c
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
        v(5)=bx0*dcos(theta)+by0*dsin(theta)
        v(6)=-bx0*dsin(theta)+by0*dcos(theta)
        v(7)=0.D0
c
        write(6,*) 'Wave no=',kwave
        write(6,*) 'Angles', theta*180.D0/pi, alpha*180.D0/pi,
     &       beta*180.D0/pi
c
        call ProjectToCharacteristicSpace(v,w0,evl,evr,lambda)
c
        write(6,*) 'Primitive Vars',v
        write(6,*) 'Characteristic Vars',w0
        write(6,*) 'Eigenvalues',lambda, dsqrt(gamma*v(8)/v(1))
        do l=1,8,1
           do k=1,8,1
              write(6,*) 'Eigenvectors',l,k,evl(l,k),evr(l,k)
           enddo
        enddo
c
c     
c        write(6,*) 'ORTHO CHECK', matmul(evl,evr)
c     
c        v0=v
        do k=1,nzlocal,1
c           x(3) = zl+(k - 0.5D0)*dz;
           x(3) = ZLloc+(k - 0.5D0)*dz;
           do j=1,nylocal,1
c              x(2) = yl+(j - 0.5D0)*dy;
              x(2) = YLloc+(j - 0.5D0)*dy;
              do i=1,nxlocal,1
c                 x(1) = xl+(i - 0.5D0)*dx;
                 x(1) = XLloc+(i - 0.5D0)*dx;
                 perturbation=pertAmplitude*
     &                dcos((kwave(1)*x(1)+kwave(2)*x(2)+kwave(3)*x(3)))
                 do l=1,8,1
                    w(l)=w0(l)
                    if(l.eq.pdir) then
                       w(l)=w0(l)+perturbation
                    endif
                 enddo
c
                 do l=1,8,1
                    v(l)=0.D0
                    do m=1,8,1
                       v(l)=v(l)+evr(l,m)*w(m)
                    enddo
                 enddo
c
c                 v(4)=0.D0
c                 v(7)=0.D0
                 un=v(2)
                 ut=v(3)
                 bn=v(5)
                 bt=v(6)
                 if(bn.eq.0.D0) then
                    write(6,*) 'ZERO',v
                 endif
                 v(2)=un*dcos(theta)-ut*dsin(theta)
                 v(3)=un*dsin(theta)+ut*dcos(theta)
                 v(5)=bn*dcos(theta)-bt*dsin(theta)
                 v(6)=bn*dsin(theta)+bt*dcos(theta)
c     
                 ux(i,j,k,1)=v(1)
                 ux(i,j,k,2)=v(1)*v(2)
                 ux(i,j,k,3)=v(1)*v(3)
                 ux(i,j,k,4)=v(1)*v(4)
                 ux(i,j,k,5)=v(5)
                 ux(i,j,k,6)=v(6)
                 ux(i,j,k,7)=v(7)
                 ux(i,j,k,8)=v(8)/(gamma-1.D0)+
     &                0.5D0*v(1)*(v(2)*v(2)+v(3)*v(3)+v(4)*v(4))
     &                +0.5D0*(v(5)*v(5)+v(6)*v(6)+v(7)*v(7))
              enddo
           enddo
        enddo
                 
c$$$        i=1
c$$$        do k=1,nzlocal,1
c$$$           do j=1,nylocal,1
c$$$              x(2) = (j - 0.5D0)*dy;
c$$$              x(1) = (i - 0.5D0)*dx;
c$$$              write(21,*) x(2),
c$$$     &             ux(i,j,k,1),ux(i,j,k,2),
c$$$     &             ux(i,j,k,3),
c$$$     &             ux(i,j,k,4),
c$$$     &             ux(i,j,k,5),
c$$$     &             ux(i,j,k,6),
c$$$     &             ux(i,j,k,7),
c$$$     &             ux(i,j,k,8)
c$$$           enddo
c$$$        enddo
c     
c$$$        j=1
c$$$        do k=1,nzlocal,1
c$$$           do i=1,nxlocal,1
c$$$              x(2) = (j - 0.5D0)*dy;
c$$$              x(1) = (i - 0.5D0)*dx;
c$$$              write(22,*) x(1),
c$$$     &             ux(i,j,k,1),ux(i,j,k,2),
c$$$     &             ux(i,j,k,3),
c$$$     &             ux(i,j,k,4),
c$$$     &             ux(i,j,k,5),
c$$$     &             ux(i,j,k,6),
c$$$     &             ux(i,j,k,7),
c$$$     &             ux(i,j,k,8)
c$$$           enddo
c$$$        enddo
c$$$        return
c     
      return
      end

c-----------------------------------------------------------------------
      subroutine ProjectToCharacteristicSpace(v,w,evl,evr,lambda)
      use properties
c     
      real*8 v(8)
      real*8 w(8)
      real*8 lambda(8)
      real*8 vi, vj, vk
      real*8 bi, bj, bk
      real*8 evl(8,8), evr(8,8), rhs(8)
      real*8 rhoInv
      real*8 Asq
      real*8 axsq
      real*8 csndsq
      real*8 cfast
      real*8 cslow
      real*8 tmp
      real*8 alphaS,alphaF
      real*8 betaY,betaZ
      real*8 betaPerpMagInv
      real*8 Sbx
      real*8 evec
      real*8 gamm
      real*8 dtodx
      real*8 pl(8), pr(8), up, um
      real*8 betaPerpMag
      integer l,k
c     
      real*8 TwoSqrtInv
      TwoSqrtInv=0.70710678D0
c     
         do k=1,8,1
            do l=1,8,1
               evr(k,l)=0.D0
               evl(k,l)=0.D0
            enddo
         enddo
         rhoInv=1.D0/v(1)
c     
         axsq=v(5)**2*rhoInv
         Asq=(v(5)**2+v(6)**2+v(7)**2)*rhoInv
         csndsq=gamma*v(8)*rhoInv
         tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.D0*csndsq*axsq),0.D0))
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
         alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
         Sbx=dsign(1.D0,v(5))
c     
         betaPerpMag=dsqrt(v(6)**2+v(7)**2)
c         if(v(6).eq.0.D0.and.v(7).eq.0.D0) then
         if(betaPerpMag.eq.0.D0) then
            betaPerpMagInv=1.D0
            betaY=TwoSqrtInv
            betaZ=TwoSqrtInv
         else
c            betaPerpMagInv=1.D0/dsqrt(v(6)**2+v(7)**2)
            betaPerpMagInv=1.D0/betaPerpMag
            betaY=v(6)*betaPerpMagInv
            betaZ=v(7)*betaPerpMagInv
         endif
c     
c         if(betaPerpMagInv.gt.1.D08) then
c            betaPerpMagInv=1.D0
c            betaY=TwoSqrtInv
c            betaZ=TwoSqrtInv
c         endif
c
         cfast=dsqrt(cfast)
         cslow=dsqrt(dmax1(cslow,0.D0))
c
         axsq=dsqrt(dmax1(axsq,0.D0))
c
         lambda(1)=v(2)
         lambda(2)=v(2)
         lambda(3)=v(2)+axsq
         lambda(4)=v(2)-axsq
         lambda(5)=v(2)+cfast
         lambda(6)=v(2)-cfast
         lambda(7)=v(2)+cslow
         lambda(8)=v(2)-cslow
c
c         write(6,*) 'LAMDAs',lambda
c     Eigenvector corresponding to u
         evl(1,1)=1.D0
         evl(1,8)=-1.D0/csndsq
c     
c     Eigenvector corresponding to u
         evl(2,5)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evl(3,3)=-betaZ*TwoSqrtInv
         evl(3,4)=betaY*TwoSqrtInv
         evl(3,6)=1.D0/dsqrt(v(1))*TwoSqrtInv*betaZ*SBx
         evl(3,7)=-1.D0/dsqrt(v(1))*TwoSqrtInv*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evl(4,3)=-betaZ*TwoSqrtInv
         evl(4,4)=betaY*TwoSqrtInv
         evl(4,6)=-evl(3,6)
         evl(4,7)=-evl(3,7)
c     
c     Eigenvector corresponding to u+cfast
         evl(5,1)=0.D0
         evl(5,2)=0.5D0*alphaF*cfast/csndsq
         evl(5,3)=-0.5D0*alphaS*cslow*betaY*Sbx/csndsq
         evl(5,4)=-0.5D0*alphaS*cslow*betaZ*Sbx/csndsq
         evl(5,6)=0.5D0*alphaS/dsqrt(v(1)*csndsq)*betaY
         evl(5,7)=0.5D0*alphaS/dsqrt(v(1)*csndsq)*betaZ
         evl(5,8)=0.5D0*alphaF/v(1)/csndsq
c     
c     Eigenvector corresponding to u-cfast
         evl(6,1)=evl(5,1)
         evl(6,2)=-evl(5,2)
         evl(6,3)=-evl(5,3)
         evl(6,4)=-evl(5,4)
         evl(6,6)=evl(5,6)
         evl(6,7)=evl(5,7)
         evl(6,8)=evl(5,8)
c     
c     Eigenvector corresponding to u+cslow
         evl(7,1)=0.D0
         evl(7,2)=0.5D0*alphaS*cslow/csndsq
         evl(7,3)=0.5D0*alphaF*cfast*betaY*Sbx/csndsq
         evl(7,4)=0.5D0*alphaF*cfast*betaZ*Sbx/csndsq
         evl(7,6)=-0.5D0*alphaF/dsqrt(v(1)*csndsq)*betaY
         evl(7,7)=-0.5D0*alphaF/dsqrt(v(1)*csndsq)*betaZ
         evl(7,8)=0.5D0*alphaS/v(1)/csndsq
c     
c     
c     Eigenvector corresponding to u-cslow
         evl(8,1)=evl(7,1)
         evl(8,2)=-evl(7,2)
         evl(8,3)=-evl(7,3)
         evl(8,4)=-evl(7,4)
         evl(8,6)=evl(7,6)
         evl(8,7)=evl(7,7)
         evl(8,8)=evl(7,8)
c
c     Right eigenvectors
c     Eigenvector corresponding to u
         evr(1,1)=1.D0
c     
c     Eigenvector corresponding to u
         evr(5,2)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evr(3,3)=-betaZ*TwoSqrtInv
         evr(4,3)=betaY*TwoSqrtInv
         evr(6,3)=dsqrt(v(1))*TwoSqrtInv*betaZ*SBx
         evr(7,3)=-dsqrt(v(1))*TwoSqrtInv*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evr(3,4)=-betaZ*TwoSqrtInv
         evr(4,4)=betaY*TwoSqrtInv
         evr(6,4)=-evr(6,3)
         evr(7,4)=-evr(7,3)
c     
c     Eigenvector corresponding to u+cfast
         evr(1,5)=v(1)*alphaF
         evr(2,5)=alphaF*cfast
         evr(3,5)=-alphaS*cslow*betaY*Sbx
         evr(4,5)=-alphaS*cslow*betaZ*Sbx
         evr(6,5)=alphaS*dsqrt(v(1)*csndsq)*betaY
         evr(7,5)=alphaS*dsqrt(v(1)*csndsq)*betaZ
         evr(8,5)=v(1)*alphaF*csndsq
c     
c     Eigenvector corresponding to u-cfast
         evr(1,6)=evr(1,5)
         evr(2,6)=-evr(2,5)
         evr(3,6)=-evr(3,5)
         evr(4,6)=-evr(4,5)
         evr(6,6)=evr(6,5)
         evr(7,6)=evr(7,5)
         evr(8,6)=evr(8,5)
c     
c     Eigenvector corresponding to u+cslow
         evr(1,7)=v(1)*alphaS
         evr(2,7)=alphaS*cslow
         evr(3,7)=alphaF*cfast*betaY*Sbx
         evr(4,7)=alphaF*cfast*betaZ*Sbx
         evr(6,7)=-alphaF*dsqrt(v(1)*csndsq)*betaY
         evr(7,7)=-alphaF*dsqrt(v(1)*csndsq)*betaZ
         evr(8,7)=v(1)*alphaS*csndsq
c     
c     
c     Eigenvector corresponding to u-cslow
         evr(1,8)=evr(1,7)
         evr(2,8)=-evr(2,7)
         evr(3,8)=-evr(3,7)
         evr(4,8)=-evr(4,7)
         evr(6,8)=evr(6,7)
         evr(7,8)=evr(7,7)
         evr(8,8)=evr(8,7)
c     
         do l=1,8,1
            w(l)=0.D0
            do k=1,8,1
               w(l)=w(l)+evl(l,k)*v(k)
            enddo
         enddo
      return
      end subroutine ProjectToCharacteristicSpace
c-----------------------------------------------------------------------
      subroutine WaveSolutionError(ux,ttot)
c----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use iounits
      use waveparams
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: uxa(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: ttot
c
c
        real*8 x(3), kwave(3)
        real*8 v(8), v0(8),w0(8), w(8)
        real*8 evl(8,8), evr(8,8),lambda(8)
        real*8 theta, beta
        real*8 un,ut,bn,bt, uz,bz
        real*8 bx0,by0
        real*8 pi
        real*8 lambdax,lambday
        real*8 linferr(8), l2err(8)
        real*8 l2sum(8)
c
        write(6,*) 'WAVE SOLUTION ERROR [L-2, L-inf]'
c
        v(1)=1.D0
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
        v(5)=1.D0
        v(6)=0.D0
        v(7)=0.D0
        v(8)=0.01D0
c
        pi=4.D0*datan(1.D0)
        kwave(1)=pi*waveNumber
c        kwave(1)=4.D0*datan(1.D0)*waveNumber
c        kwave(2)=0.D0
c        kwave(2)=kwave(1)*dtan(theta)
        kwave(2)=kratio*kwave(1)
        kwave(3)=0.D0
c        kwave(2)=kwave(1)
c        kwave(3)=kwave(1)
c        theta=datan(2.D0)
c     Wave vector is at an angle theta to x
        theta=datan(kwave(2)/kwave(1))
c
c     B0 is aligned at an angle alpha to x
c        alpha=-pi/4.D0
c        alpha=pi/180.D0*(44.D0)
c        alpha=0.D0
        bx0=dcos(alpha)
        by0=dsin(alpha)
c
        beta=alpha+theta
c
c        kwave(1)=kwave(1)*dcos(beta)
c        kwave(2)=kwave(1)*dsin(beta)
c
        v(2)=0.D0
        v(3)=0.D0
        v(4)=0.D0
        v(5)=bx0*dcos(theta)+by0*dsin(theta)
        v(6)=-bx0*dsin(theta)+by0*dcos(theta)
        v(7)=0.D0
c
        call ProjectToCharacteristicSpace(v,w0,evl,evr,lambda)
c
        call ConservedToPrimitive(ux,vx)
c
        l2sum=0.D0
        linferr=-1.D0
        lambdax=lambda(pdir)*cos(theta)
        lambday=lambda(pdir)*sin(theta)
        do k=izlo,izhi,1
c$$$           x(3) = zl+(k - 0.5D0)*dz;
           x(3) = ZLloc+(k - 0.5D0)*dz;
           do j=iylo,iyhi,1
c$$$              x(2) = yl+(j - 0.5D0)*dy;
              x(2) = YLloc+(j - 0.5D0)*dy;
              do i=ixlo,ixhi,1
c$$$                 x(1) = xl+(i - 0.5D0)*dx;
                 x(1) = XLloc+(i - 0.5D0)*dx;
                 perturbation=pertAmplitude*
     &                dcos((kwave(1)*x(1)+kwave(2)*x(2)+kwave(3)*x(3))
     &                -kwave(1)*lambdax*ttot
     &                -kwave(2)*lambday*ttot)
                 do l=1,8,1
                    w(l)=w0(l)
                    if(l.eq.pdir) then
                       w(l)=w0(l)+perturbation
                    endif
                 enddo
c
                 do l=1,8,1
                    v(l)=0.D0
                    do m=1,8,1
                       v(l)=v(l)+evr(l,m)*w(m)
                    enddo
                 enddo
c
                 un=v(2)
                 ut=v(3)
                 bn=v(5)
                 bt=v(6)
                 if(bn.eq.0.D0) then
                    write(6,*) '  ZERO',v
                 endif
                 v(2)=un*dcos(theta)-ut*dsin(theta)
                 v(3)=un*dsin(theta)+ut*dcos(theta)
                 v(5)=bn*dcos(theta)-bt*dsin(theta)
                 v(6)=bn*dsin(theta)+bt*dcos(theta)
c     
                 uxa(i,j,k,1)=v(1)
                 uxa(i,j,k,2)=v(2)*v(1)
                 uxa(i,j,k,3)=v(3)*v(1)
                 uxa(i,j,k,4)=v(4)*v(1)
                 uxa(i,j,k,5)=v(5)
                 uxa(i,j,k,6)=v(6)
                 uxa(i,j,k,7)=v(7)
                 uxa(i,j,k,8)=v(8)/(gamma-1.D0)+
     &                0.5D0*v(1)*(v(2)*v(2)+v(3)*v(3)+v(4)*v(4))
     &                +0.5D0*(v(5)*v(5)+v(6)*v(6)+v(7)*v(7))
c$$$                 l2sum(1)=l2sum(1)+(ux(i,j,k,1)-v(1))**2
c$$$                 linferr(1)=max(linferr(1),dabs(ux(i,j,k,1)-v(1)))
c$$$                 l2sum(2)=l2sum(2)+(ux(i,j,k,2)-v(2))**2
c$$$                 linferr(2)=max(linferr(2),dabs(ux(i,j,k,2)-v(2)))
c$$$                 l2sum(3)=l2sum(3)+(ux(i,j,k,3)-v(3))**2
c$$$                 linferr(3)=max(linferr(3),dabs(ux(i,j,k,3)-v(3)))
c$$$                 l2sum(4)=l2sum(4)+(ux(i,j,k,4)-v(4))**2
c$$$                 linferr(4)=max(linferr(4),dabs(ux(i,j,k,4)-v(4)))
c$$$                 l2sum(5)=l2sum(5)+(ux(i,j,k,5)-v(5))**2
c$$$                 linferr(5)=max(linferr(5),dabs(ux(i,j,k,5)-v(5)))
c$$$                 l2sum(6)=l2sum(6)+(ux(i,j,k,6)-v(6))**2
c$$$                 linferr(6)=max(linferr(6),dabs(ux(i,j,k,6)-v(6)))
c$$$                 l2sum(7)=l2sum(7)+(ux(i,j,k,7)-v(7))**2
c$$$                 linferr(7)=max(linferr(7),dabs(ux(i,j,k,7)-v(7)))
c$$$                 l2sum(8)=l2sum(8)+(ux(i,j,k,8)-v(8))**2
c$$$                 linferr(8)=max(linferr(8),dabs(ux(i,j,k,8)-v(8)))
                 l2sum(1)=l2sum(1)+(vx(i,j,k,1)-v(1))**2
                 linferr(1)=max(linferr(1),dabs(vx(i,j,k,1)-v(1)))
                 l2sum(2)=l2sum(2)+(vx(i,j,k,2)-v(2))**2
                 linferr(2)=max(linferr(2),dabs(vx(i,j,k,2)-v(2)))
                 l2sum(3)=l2sum(3)+(vx(i,j,k,3)-v(3))**2
                 linferr(3)=max(linferr(3),dabs(vx(i,j,k,3)-v(3)))
                 l2sum(4)=l2sum(4)+(vx(i,j,k,4)-v(4))**2
                 linferr(4)=max(linferr(4),dabs(vx(i,j,k,4)-v(4)))
                 l2sum(5)=l2sum(5)+(vx(i,j,k,5)-v(5))**2
                 linferr(5)=max(linferr(5),dabs(vx(i,j,k,5)-v(5)))
                 l2sum(6)=l2sum(6)+(vx(i,j,k,6)-v(6))**2
                 linferr(6)=max(linferr(6),dabs(vx(i,j,k,6)-v(6)))
                 l2sum(7)=l2sum(7)+(vx(i,j,k,7)-v(7))**2
                 linferr(7)=max(linferr(7),dabs(vx(i,j,k,7)-v(7)))
                 l2sum(8)=l2sum(8)+(vx(i,j,k,8)-v(8))**2
                 linferr(8)=max(linferr(8),dabs(vx(i,j,k,8)-v(8)))
              enddo
           enddo
        enddo
        l2err=dsqrt(l2sum/nx/ny)
        write(6,*) '  LERROR: rho =',l2err(1), linferr(1)
        write(6,*) '  LERROR:  ux =',l2err(2), linferr(2)
        write(6,*) '  LERROR:  uy =',l2err(3), linferr(3)
        write(6,*) '  LERROR:  uz =',l2err(4), linferr(4)
        write(6,*) '  LERROR:  Bx =',l2err(5), linferr(5)
        write(6,*) '  LERROR:  By =',l2err(6), linferr(6)
        write(6,*) '  LERROR:  Bz =',l2err(7), linferr(7)
        write(6,*) '  LERROR:   e =',l2err(8), linferr(8)
c                 
c$$$        i=1
c$$$        do k=1,nzlocal,1
c$$$           do j=1,nylocal,1
c$$$              x(2) = (j - 0.5D0)*dy;
c$$$              x(1) = (i - 0.5D0)*dx;
c$$$              write(21,*) x(2),
c$$$     &             ux(i,j,k,1),ux(i,j,k,2),
c$$$     &             ux(i,j,k,3),
c$$$     &             ux(i,j,k,4),
c$$$     &             ux(i,j,k,5),
c$$$     &             ux(i,j,k,6),
c$$$     &             ux(i,j,k,7),
c$$$     &             ux(i,j,k,8)
c$$$           enddo
c$$$        enddo

c$$$        j=iylo
c$$$        do k=izlo,izhi,1
c$$$           do i=ixlo,ixhi,1
c$$$              x(2) = (j - 0.5D0)*dy;
c$$$              x(1) = (i - 0.5D0)*dx;
c$$$              write(22,*) x(1),
c$$$     &             ux(i,j,k,1),ux(i,j,k,2),
c$$$     &             ux(i,j,k,3),
c$$$     &             ux(i,j,k,4),
c$$$     &             ux(i,j,k,5),
c$$$     &             ux(i,j,k,6),
c$$$     &             ux(i,j,k,7),
c$$$     &             ux(i,j,k,8)
c$$$           enddo
c$$$        enddo
c     
      return
      end

