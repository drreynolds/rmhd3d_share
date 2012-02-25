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

      module evectors_conservative
      contains
      subroutine SetEigenSystem(utilde,evl,evr,alamda,ilo,ihi,direction)
      use mesh
      use mesh_common
      use properties
c     
      double precision:: utilde(INLO:INHI,nvar)
      double precision:: evl(INLO:INHI,nvar,nvar)
      double precision:: evr(INLO:INHI,nvar,nvar)
      double precision:: alamda(INLO:INHI,nvar)
      integer:: i,j,k,l,m,ilo,ihi
      double precision:: uabs(INLO:INHI)
      double precision:: rhoInv(INLO:INHI)
      double precision:: Asq(INLO:INHI)
      double precision:: axsq(INLO:INHI)
      double precision:: csndsq(INLO:INHI)
      double precision:: cfast(INLO:INHI)
      double precision:: cslow(INLO:INHI)
      double precision:: tmp(INLO:INHI)
      logical:: Mask(INLO:INHI)
      double precision:: rho,u,v,w,p,c,eint
      double precision:: bx,by,bz
      double precision:: alphaS(INLO:INHI),alphaF(INLO:INHI)
      double precision:: betaY(INLO:INHI),betaZ(INLO:INHI)
      double precision:: betaPerpMagInv(INLO:INHI)
      double precision:: Sbx(INLO:INHI)
      double precision:: evec(nvar,nvar)
      double precision:: gamm
c     
      double precision:: TwoSqrtInv
      double precision,parameter:: epsilon=1.D-02
      integer:: ii,direction
c     utilde(1)=rho, 
c     utilde(2)=u
c     utilde(3)=v
c     utilde(4)=w
c     utilde(5)=Bx
c     utilde(6)=By
c     utilde(7)=Bz
c     utilde(8)=p
c     
c     
      TwoSqrtInv = 1.D0/dsqrt(2.D0)
      
      Mask=.false.
      Mask(ilo:ihi)=.true.
      where(mask)
c     
         uabs=abs(utilde(:,2))
         rhoInv(:)=1.D0/utilde(:,1)
c     
         axsq=utilde(:,5)**2*rhoInv(:)
         Asq=(utilde(:,5)**2+utilde(:,6)**2+utilde(:,7)**2)*rhoInv
         csndsq=gamma*utilde(:,8)*rhoInv
         tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.D0*csndsq*axsq),0.D0))
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaS=dsqrt((dmax1(cfast-csndsq,0.D0))/(cfast-cslow))
         alphaF=dsqrt((dmax1(csndsq-cslow,0.D0))/(cfast-cslow))
c     
         Sbx=dsign(1.D0,utilde(:,5))
c     
      endwhere
      betaY=TwoSqrtInv
      betaZ=TwoSqrtInv
      where(utilde(:,6).eq.0.D0.and.utilde(:,7).eq.0.D0)
         betaPerpMagInv=1.D0
         betaY=TwoSqrtInv
         betaZ=TwoSqrtInv
      elsewhere
         betaPerpMagInv=1.D0/dsqrt(utilde(:,6)**2+utilde(:,7)**2)
         betaY=utilde(:,6)*betaPerpMagInv
         betaZ=utilde(:,7)*betaPerpMagInv
      endwhere
c     
      evl=0.D0
      evr=0.D0
      where(mask)
         cfast=dsqrt(cfast)
	 cslow=dsqrt(dmax1(cslow,0.D0))
         axsq=dsqrt(axsq)
c     
         alamda(:,1)=uabs
         alamda(:,2)=uabs
         alamda(:,3)=abs(utilde(:,2)+axsq)
         alamda(:,4)=abs(utilde(:,2)-axsq)
         alamda(:,5)=abs(utilde(:,2)+cfast)
         alamda(:,6)=abs(utilde(:,2)-cfast)
         alamda(:,7)=abs(utilde(:,2)+cslow)
         alamda(:,8)=abs(utilde(:,2)-cslow)
c
c$$$         alamda(:,1)=utilde(:,2)
c$$$         alamda(:,2)=utilde(:,2)
c$$$         alamda(:,3)=(utilde(:,2)+axsq)
c$$$         alamda(:,4)=(utilde(:,2)-axsq)
c$$$         alamda(:,5)=(utilde(:,2)+cfast)
c$$$         alamda(:,6)=(utilde(:,2)-cfast)
c$$$         alamda(:,7)=(utilde(:,2)+cslow)
c$$$         alamda(:,8)=(utilde(:,2)-cslow)
c     
c     Eigenvector corresponding to u
         evr(:,1,1)=1.D0
c     
c     Eigenvector corresponding to u
         evr(:,5,2)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evr(:,3,3)=-betaZ*TwoSqrtInv
         evr(:,4,3)=betaY*TwoSqrtInv
         evr(:,6,3)=dsqrt(utilde(:,1))*TwoSqrtInv*betaZ*SBx
         evr(:,7,3)=-dsqrt(utilde(:,1))*TwoSqrtInv*betaY*SBx
c
c     Eigenvector corresponding to u-ax
        evr(:,3,4)=-betaZ*TwoSqrtInv
        evr(:,4,4)=betaY*TwoSqrtInv
        evr(:,6,4)=-evr(:,6,3)
        evr(:,7,4)=-evr(:,7,3)
c
c     Eigenvector corresponding to u+cfast
        evr(:,1,5)=utilde(:,1)*alphaF
        evr(:,2,5)=alphaF*cfast
        evr(:,3,5)=-alphaS*cslow*betaY*Sbx
        evr(:,4,5)=-alphaS*cslow*betaZ*Sbx
        evr(:,6,5)=alphaS*dsqrt(utilde(:,1)*csndsq)*betaY
        evr(:,7,5)=alphaS*dsqrt(utilde(:,1)*csndsq)*betaZ
        evr(:,8,5)=utilde(:,1)*alphaF*csndsq
c
c     Eigenvector corresponding to u-cfast
        evr(:,1,6)=evr(:,1,5)
        evr(:,2,6)=-evr(:,2,5)
        evr(:,3,6)=-evr(:,3,5)
        evr(:,4,6)=-evr(:,4,5)
        evr(:,6,6)=evr(:,6,5)
        evr(:,7,6)=evr(:,7,5)
        evr(:,8,6)=evr(:,8,5)
c
c     Eigenvector corresponding to u+cslow
        evr(:,1,7)=utilde(:,1)*alphaS
        evr(:,2,7)=alphaS*cslow
        evr(:,3,7)=alphaF*cfast*betaY*Sbx
        evr(:,4,7)=alphaF*cfast*betaZ*Sbx
        evr(:,6,7)=-alphaF*dsqrt(utilde(:,1)*csndsq)*betaY
        evr(:,7,7)=-alphaF*dsqrt(utilde(:,1)*csndsq)*betaZ
        evr(:,8,7)=utilde(:,1)*alphaS*csndsq
c
c
c     Eigenvector corresponding to u-cslow
        evr(:,1,8)=evr(:,1,7)
        evr(:,2,8)=-evr(:,2,7)
        evr(:,3,8)=-evr(:,3,7)
        evr(:,4,8)=-evr(:,4,7)
        evr(:,6,8)=evr(:,6,7)
        evr(:,7,8)=evr(:,7,7)
        evr(:,8,8)=evr(:,8,7)

c
c     Eigenvector corresponding to u
        evl(:,1,1)=1.D0
        evl(:,1,8)=-1.D0/csndsq
c
c     Eigenvector corresponding to u
        evl(:,2,5)=1.D0
c
c     Eigenvector corresponding to u+ax
        evl(:,3,3)=-betaZ*TwoSqrtInv
        evl(:,3,4)=betaY*TwoSqrtInv
        evl(:,3,6)=1.D0/dsqrt(utilde(:,1))*TwoSqrtInv*betaZ*SBx
        evl(:,3,7)=-1.D0/dsqrt(utilde(:,1))*TwoSqrtInv*betaY*SBx
c
c     Eigenvector corresponding to u-ax
        evl(:,4,3)=-betaZ*TwoSqrtInv
        evl(:,4,4)=betaY*TwoSqrtInv
        evl(:,4,6)=-evl(:,3,6)
        evl(:,4,7)=-evl(:,3,7)
c
c     Eigenvector corresponding to u+cfast
        evl(:,5,1)=0.D0
        evl(:,5,2)=0.5D0*alphaF*cfast/csndsq
        evl(:,5,3)=-0.5D0*alphaS*cslow*betaY*Sbx/csndsq
        evl(:,5,4)=-0.5D0*alphaS*cslow*betaZ*Sbx/csndsq
        evl(:,5,6)=0.5D0*alphaS/dsqrt(utilde(:,1)*csndsq)*betaY
        evl(:,5,7)=0.5D0*alphaS/dsqrt(utilde(:,1)*csndsq)*betaZ
        evl(:,5,8)=0.5D0*alphaF/utilde(:,1)/csndsq
c
c     Eigenvector corresponding to u-cfast
        evl(:,6,1)=evl(:,5,1)
        evl(:,6,2)=-evl(:,5,2)
        evl(:,6,3)=-evl(:,5,3)
        evl(:,6,4)=-evl(:,5,4)
        evl(:,6,6)=evl(:,5,6)
        evl(:,6,7)=evl(:,5,7)
        evl(:,6,8)=evl(:,5,8)
c
c     Eigenvector corresponding to u+cslow
        evl(:,7,1)=0.D0
        evl(:,7,2)=0.5D0*alphaS*cslow/csndsq
        evl(:,7,3)=0.5D0*alphaF*cfast*betaY*Sbx/csndsq
        evl(:,7,4)=0.5D0*alphaF*cfast*betaZ*Sbx/csndsq
        evl(:,7,6)=-0.5D0*alphaF/dsqrt(utilde(:,1)*csndsq)*betaY
        evl(:,7,7)=-0.5D0*alphaF/dsqrt(utilde(:,1)*csndsq)*betaZ
        evl(:,7,8)=0.5D0*alphaS/utilde(:,1)/csndsq
c
c
c     Eigenvector corresponding to u-cslow
        evl(:,8,1)=evl(:,7,1)
        evl(:,8,2)=-evl(:,7,2)
        evl(:,8,3)=-evl(:,7,3)
        evl(:,8,4)=-evl(:,7,4)
        evl(:,8,6)=evl(:,7,6)
        evl(:,8,7)=evl(:,7,7)
        evl(:,8,8)=evl(:,7,8)

c
c
      endwhere
c
c$$$c     Check l*lamda*r
c$$$        i=(ilo+ihi)/2
c$$$c
c$$$        rho=utilde(i,1)
c$$$        u=utilde(i,2)
c$$$        v=utilde(i,3)
c$$$        w=utilde(i,4)
c$$$        bx=utilde(i,5)
c$$$        by=utilde(i,6)
c$$$        bz=utilde(i,7)
c$$$        p=utilde(i,8)
c$$$        evec=0.D0
c$$$        evec(1,1)=u
c$$$        evec(1,2)=rho
c$$$        evec(2,2)=u
c$$$        evec(2,6)=by/rho
c$$$        evec(2,7)=bz/rho
c$$$        evec(2,8)=1.D0/rho
c$$$        evec(3,3)=u
c$$$        evec(3,6)=-bx/rho
c$$$        evec(4,4)=u
c$$$        evec(4,7)=-bx/rho
c$$$        evec(5,5)=u
c$$$        evec(6,2)=by
c$$$        evec(6,3)=-bx
c$$$        evec(6,6)=u
c$$$        evec(7,2)=bz
c$$$        evec(7,4)=-bx
c$$$        evec(7,7)=u
c$$$        evec(8,2)=gamma*p
c$$$        evec(8,8)=u
c$$$c
c$$$        write(6,*)
c$$$        write(6,*) 'U=',u
c$$$        evec=matmul(evl(i,:,:),matmul(evec,evr(i,:,:)))
c$$$        do k=1,nvar,1
c$$$           write(6,1001) evec(k,k),
c$$$     &          alamda(i,k)
c$$$        enddo
c$$$ 1001   format(1h ,2(e11.4,x))
c$$$        write(6,*)
c$$$        do k=1,nvar,1
c$$$           write(6,1000) evec(k,:)
c$$$        enddo
c$$$        write(6,*) 
c$$$c        write(6,*) 'Al=lamda l'
c$$$c        do k=1,nvar,1
c$$$c           write(6,1000) matmul(evl(i,k,:),evec)
c$$$c           write(6,1000) alamda(i,k)*evl(i,k,:)
c$$$c           write(6,1000) matmul(evec,evr(i,:,k))
c$$$c           write(6,1000) alamda(i,k)*evr(i,:,k)
c$$$c        enddo
c$$$        write(6,*)
c$$$        write(6,*) 'Eigenvector orthogonality'
c$$$        evec=matmul(evl(i,:,:),evr(i,:,:))
c$$$        do k=1,nvar,1
c$$$           write(6,1000) evec(k,:)
c$$$        enddo
c
c      Conservative eigenvectors
c
      gamm=gamma-1.D0
      do i=ilo,ihi,1
         do k=1,nvar,1
c     
            evec(k,1)=evl(i,k,1)
     &           -1.D0/utilde(i,1)*(utilde(i,2)*evl(i,k,2)
     &           +utilde(i,3)*evl(i,k,3)+utilde(i,4)*evl(i,k,4))
     &           +0.5D0*gamm*(utilde(i,2)**2+utilde(i,3)**2
     &           +utilde(i,4)**2)*evl(i,k,8)
c     
            evec(k,2)=evl(i,k,2)/utilde(i,1)
     &           -gamm*utilde(i,2)*evl(i,k,8)
c     
            evec(k,3)=evl(i,k,3)/utilde(i,1)
     &           -gamm*utilde(i,3)*evl(i,k,8)
c     
            evec(k,4)=evl(i,k,4)/utilde(i,1)
     &           -gamm*utilde(i,4)*evl(i,k,8)
c     
            evec(k,5:7)=evl(i,k,5:7)
     &           -gamm*utilde(i,5:7)*evl(i,k,8)
c     
            evec(k,8)=gamm*evl(i,k,8)
c     
         enddo
         evl(i,:,:)=evec
c     
         do k=1,nvar,1
c     
            evec(1,k)=evr(i,1,k)
c     
            evec(2:4,k)=utilde(i,2:4)*evr(i,1,k)
     &           +utilde(i,1)*evr(i,2:4,k)
c     
            evec(5:7,k)=evr(i,5:7,k)
c     
            evec(8,k)=0.5D0*(utilde(i,2)**2+utilde(i,3)**2
     &           +utilde(i,4)**2)*evr(i,1,k)
     &           +utilde(i,1)*(utilde(i,2)*evr(i,2,k)
     &           +utilde(i,3)*evr(i,3,k)+utilde(i,4)*evr(i,4,k))
     &           +utilde(i,5)*evr(i,5,k)+utilde(i,6)*evr(i,6,k)
     &           +utilde(i,7)*evr(i,7,k)+evr(i,8,k)/gamm
c     
         enddo
         evr(i,:,:)=evec
      enddo
c     
      return
c     
#ifdef HARTEN
c     Harten's entropy fix
c     HARTEN'S correction for transonic case.
      do l=1,nvar,1
         where(alamda(:,l).lt.epsilon)
            alamda(:,l)=0.5D0*(alamda(:,l)**2
     &           +epsilon**2)/epsilon
         endwhere
      enddo
#endif
c     

c$$$c     Check Orthogonality of eigenvectors
c$$$        i=(ilo+ihi)/2
c$$$        evec=matmul(evl(i,:,:),evr(i,:,:))
c$$$        write(6,*)
c$$$        write(6,*) 'Conservative Ev Ortho Check'
c$$$        do k=1,nvar,1
c$$$           write(6,1000) evec(k,:)
c$$$        enddo
c$$$ 1000   format(1h ,8(e11.4,x))

      return
      end subroutine SetEigenSystem
c-----------------------------------------------------------------------






c-----------------------------------------------------------------------
c     SetEigenSystem2: optimized version of SetEigenSystem
c-----------------------------------------------------------------------
      subroutine SetEigenSystem2(vv,evl,evr,alamda,ilo,ihi,idir)
      use mesh
      use mesh_common
      use properties
c     
      implicit none
      double precision, dimension(nvar,inlo:inhi) :: vv, alamda
      double precision, dimension(nvar,nvar,inlo:inhi) :: evl, evr
      integer:: ilo, ihi, idir
c
      double precision :: rhoInv, rhoSqrt, rhoSqrtInv, rhoCSqrt
      double precision :: TwoSqrtInv, rhoCSqrtInv
      double precision :: Asq, axsq, csndsq, csndsqInv, cfast, cslow
      double precision :: alphaS, alphaF, betaY, betaZ, betaPerpMagInv
      double precision :: dtmp, kenergy, Benergy, gamm, gammInv, Sbx
      integer :: i, k
c
c     General constants
      gamm = gamma-1.D0
      gammInv = 1.d0/gamm
      TwoSqrtInv = dsqrt(0.5D0)

c     loop over space
      do i=ilo,ihi,1
c      
c     Plasma beta
         Benergy = dsqrt(vv(6,i)*vv(6,i) + vv(7,i)*vv(7,i))
         if (Benergy == 0.d0) then
            betaY = TwoSqrtInv
            betaZ = TwoSqrtInv
         else
            betaPerpMagInv = 1.D0/Benergy
            betaY = vv(6,i)*betaPerpMagInv
            betaZ = vv(7,i)*betaPerpMagInv
         endif
c     
c     Set some constants
         rhoInv = 1.D0/vv(1,i)
         axsq = vv(5,i)*vv(5,i)*rhoInv
         csndsq = gamma*vv(8,i)*rhoInv
         Asq = csndsq + 
     &        (vv(5,i)*vv(5,i)+vv(6,i)*vv(6,i)+vv(7,i)*vv(7,i))*rhoInv
         dtmp = dsqrt(dmax1((Asq*Asq-4.D0*csndsq*axsq),0.D0))
         cfast = 0.5D0*(Asq+dtmp)
         cslow = 0.5D0*(Asq-dtmp)
         dtmp = 1.d0/(cfast-cslow)
         alphaS = dsqrt((dmax1(cfast-csndsq,0.D0))*dtmp)
         alphaF = dsqrt((dmax1(csndsq-cslow,0.D0))*dtmp)
         csndsqInv = 1.D0/csndsq
         rhoSqrtInv = dsqrt(rhoInv)
         rhoSqrt = dsqrt(vv(1,i))
         rhoCSqrt = dsqrt(vv(1,i)*csndsq)
         rhoCSqrtInv = 1.D0/rhoCSqrt
         kenergy = vv(2,i)*vv(2,i)+vv(3,i)*vv(3,i)+vv(4,i)*vv(4,i) 
c     
c     Sign of Bx magnetic field
         Sbx = dsign(1.D0,vv(5,i))
c     
c     Rescale some constants for waves
         cfast = dsqrt(cfast)
	 cslow = dsqrt(dmax1(cslow,0.D0))
         axsq  = dsqrt(axsq)
c     
c     Eigenvalues
         alamda(1,i) = vv(2,i)
         alamda(2,i) = vv(2,i)
         alamda(3,i) = vv(2,i) + axsq
         alamda(4,i) = vv(2,i) - axsq
         alamda(5,i) = vv(2,i) + cfast
         alamda(6,i) = vv(2,i) - cfast
         alamda(7,i) = vv(2,i) + cslow
         alamda(8,i) = vv(2,i) - cslow
c     
c     Right eigenvector corresponding to u
         evr(1,1,i) = 1.D0
         evr(2,1,i) = 0.D0
         evr(3,1,i) = 0.D0
         evr(4,1,i) = 0.D0
         evr(5,1,i) = 0.D0
         evr(6,1,i) = 0.D0
         evr(7,1,i) = 0.D0
         evr(8,1,i) = 0.D0
c     
c     Right eigenvector corresponding to u
         evr(1,2,i) = 0.D0
         evr(2,2,i) = 0.D0
         evr(3,2,i) = 0.D0
         evr(4,2,i) = 0.D0
         evr(5,2,i) = 1.D0
         evr(6,2,i) = 0.D0
         evr(7,2,i) = 0.D0
         evr(8,2,i) = 0.D0
c     
c     Right eigenvector corresponding to u+ax
         evr(1,3,i) =  0.D0
         evr(2,3,i) =  0.D0
         evr(3,3,i) = -betaZ*TwoSqrtInv
         evr(4,3,i) =  betaY*TwoSqrtInv
         evr(5,3,i) =  0.D0
         evr(6,3,i) =  rhoSqrt*TwoSqrtInv*betaZ*SBx
         evr(7,3,i) = -rhoSqrt*TwoSqrtInv*betaY*SBx
         evr(8,3,i) =  0.D0
c
c     Right eigenvector corresponding to u-ax
         evr(1,4,i) =  0.D0
         evr(2,4,i) =  0.D0
         evr(3,4,i) =  evr(3,3,i)
         evr(4,4,i) =  evr(4,3,i)
         evr(5,4,i) =  0.D0
         evr(6,4,i) = -evr(6,3,i)
         evr(7,4,i) = -evr(7,3,i)
         evr(8,4,i) =  0.D0
c
c     Right eigenvector corresponding to u+cfast
         evr(1,5,i) =  vv(1,i)*alphaF
         evr(2,5,i) =  alphaF*cfast
         evr(3,5,i) = -alphaS*cslow*betaY*Sbx
         evr(4,5,i) = -alphaS*cslow*betaZ*Sbx
         evr(5,5,i) =  0.D0
         evr(6,5,i) =  alphaS*rhoCSqrt*betaY
         evr(7,5,i) =  alphaS*rhoCSqrt*betaZ
         evr(8,5,i) =  vv(1,i)*alphaF*csndsq
c     
c     Right eigenvector corresponding to u-cfast
         evr(1,6,i) =  evr(1,5,i)
         evr(2,6,i) = -evr(2,5,i)
         evr(3,6,i) = -evr(3,5,i)
         evr(4,6,i) = -evr(4,5,i)
         evr(5,6,i) =  0.D0
         evr(6,6,i) =  evr(6,5,i)
         evr(7,6,i) =  evr(7,5,i)
         evr(8,6,i) =  evr(8,5,i)
c     
c     Right eigenvector corresponding to u+cslow
         evr(1,7,i) =  vv(1,i)*alphaS
         evr(2,7,i) =  alphaS*cslow
         evr(3,7,i) =  alphaF*cfast*betaY*Sbx
         evr(4,7,i) =  alphaF*cfast*betaZ*Sbx
         evr(5,7,i) =  0.D0
         evr(6,7,i) = -alphaF*rhoCSqrt*betaY
         evr(7,7,i) = -alphaF*rhoCSqrt*betaZ
         evr(8,7,i) =  vv(1,i)*alphaS*csndsq
c     
c     
c     Right eigenvector corresponding to u-cslow
         evr(1,8,i) =  evr(1,7,i)
         evr(2,8,i) = -evr(2,7,i)
         evr(3,8,i) = -evr(3,7,i)
         evr(4,8,i) = -evr(4,7,i)
         evr(5,8,i) =  0.D0
         evr(6,8,i) =  evr(6,7,i)
         evr(7,8,i) =  evr(7,7,i)
         evr(8,8,i) =  evr(8,7,i)
c
c     Convert primitive to conservative eigenvectors
         do k=1,nvar,1
c     
            evr(8,k,i) = 0.5D0*kenergy*evr(1,k,i)
     &           + vv(1,i)*(vv(2,i)*evr(2,k,i)
     &           + vv(3,i)*evr(3,k,i)+vv(4,i)*evr(4,k,i))
     &           + vv(5,i)*evr(5,k,i)+vv(6,i)*evr(6,k,i)
     &           + vv(7,i)*evr(7,k,i)+evr(8,k,i)*gammInv
c     
            evr(2:4,k,i) = vv(2:4,i)*evr(1,k,i)
     &           +vv(1,i)*evr(2:4,k,i)
         enddo
c
c     
c     Left eigenvector corresponding to u
         evl(1,1,i) =  1.D0
         evl(1,2,i) =  0.D0
         evl(1,3,i) =  0.D0
         evl(1,4,i) =  0.D0
         evl(1,5,i) =  0.D0
         evl(1,6,i) =  0.D0
         evl(1,7,i) =  0.D0
         evl(1,8,i) = -csndsqInv
c     
c     Left eigenvector corresponding to u
         evl(2,1,i) = 0.D0
         evl(2,2,i) = 0.D0
         evl(2,3,i) = 0.D0
         evl(2,4,i) = 0.D0
         evl(2,5,i) = 1.D0
         evl(2,6,i) = 0.D0
         evl(2,7,i) = 0.D0
         evl(2,8,i) = 0.D0
c
c     Left eigenvector corresponding to u+ax
         evl(3,1,i) =  0.D0
         evl(3,2,i) =  0.D0
         evl(3,3,i) = -betaZ*TwoSqrtInv
         evl(3,4,i) =  betaY*TwoSqrtInv
         evl(3,5,i) =  0.D0
         evl(3,6,i) =  rhoSqrtInv*TwoSqrtInv*betaZ*SBx
         evl(3,7,i) = -rhoSqrtInv*TwoSqrtInv*betaY*SBx
         evl(3,8,i) =  0.D0
c     
c     Left eigenvector corresponding to u-ax
         evl(4,1,i) =  0.D0
         evl(4,2,i) =  0.D0
         evl(4,3,i) =  evl(3,3,i)
         evl(4,4,i) =  evl(3,4,i)
         evl(4,5,i) =  0.D0
         evl(4,6,i) = -evl(3,6,i)
         evl(4,7,i) = -evl(3,7,i)
         evl(4,8,i) =  0.D0
c     
c     Left eigenvector corresponding to u+cfast
         evl(5,1,i) =  0.D0
         evl(5,2,i) =  0.5D0*alphaF*cfast*csndsqInv
         evl(5,3,i) = -0.5D0*alphaS*cslow*betaY*Sbx*csndsqInv
         evl(5,4,i) = -0.5D0*alphaS*cslow*betaZ*Sbx*csndsqInv
         evl(5,5,i) =  0.D0
         evl(5,6,i) =  0.5D0*alphaS*rhoCSqrtInv*betaY
         evl(5,7,i) =  0.5D0*alphaS*rhoCSqrtInv*betaZ
         evl(5,8,i) =  0.5D0*alphaF*rhoInv*csndsqInv
c     
c     Left eigenvector corresponding to u-cfast
         evl(6,1,i) =  0.D0
         evl(6,2,i) = -evl(5,2,i)
         evl(6,3,i) = -evl(5,3,i)
         evl(6,4,i) = -evl(5,4,i)
         evl(6,5,i) =  0.D0
         evl(6,6,i) =  evl(5,6,i)
         evl(6,7,i) =  evl(5,7,i)
         evl(6,8,i) =  evl(5,8,i)
c     
c     Left eigenvector corresponding to u+cslow
         evl(7,1,i) =  0.D0
         evl(7,2,i) =  0.5D0*alphaS*cslow*csndsqInv
         evl(7,3,i) =  0.5D0*alphaF*cfast*betaY*Sbx*csndsqInv
         evl(7,4,i) =  0.5D0*alphaF*cfast*betaZ*Sbx*csndsqInv
         evl(7,5,i) =  0.D0
         evl(7,6,i) = -0.5D0*alphaF*rhoCSqrtInv*betaY
         evl(7,7,i) = -0.5D0*alphaF*rhoCSqrtInv*betaZ
         evl(7,8,i) =  0.5D0*alphaS*rhoInv*csndsqInv
c     
c     Left eigenvector corresponding to u-cslow
         evl(8,1,i) =  0.D0
         evl(8,2,i) = -evl(7,2,i)
         evl(8,3,i) = -evl(7,3,i)
         evl(8,4,i) = -evl(7,4,i)
         evl(8,5,i) =  0.D0
         evl(8,6,i) =  evl(7,6,i)
         evl(8,7,i) =  evl(7,7,i)
         evl(8,8,i) =  evl(7,8,i)
c
c
c     Convert primitive to conservative eigenvectors
         do k=1,nvar,1
c
            evl(k,1,i) = evl(k,1,i) 
     &           - rhoInv*(vv(2,i)*evl(k,2,i)
     &           +vv(3,i)*evl(k,3,i)+vv(4,i)*evl(k,4,i))
     &           +0.5D0*gamm*kenergy*evl(k,8,i)
c
            evl(k,2:4,i) = evl(k,2:4,i)*rhoInv
     &           -gamm*vv(2:4,i)*evl(k,8,i)
c
            evl(k,5:7,i) = evl(k,5:7,i) 
     &           -gamm*vv(5:7,i)*evl(k,8,i)
c
            evl(k,8,i) = gamm*evl(k,8,i)
c
         enddo
c
      enddo
c
      return
      end subroutine SetEigenSystem2
c-----------------------------------------------------------------------
      end module evectors_conservative
c-----------------------------------------------------------------------
     
