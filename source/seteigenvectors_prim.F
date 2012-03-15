c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2003
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------
c     $Log: seteigenvectors_prim.F,v $
c-----------------------------------------------------------------
      module evectors_primitive
      contains
      subroutine SetEigenSystem(utilde,evl,evr,alamda,ilo,ihi,direction)
c-----------------------------------------------------------------------
c     Description: Calculate of eigenvalues and eigenvectors in the
c                  8-wave formulation.
c     Input: utilde: some average primitive variables at cell centers
c            direction: (not used)
c     Output: evl,evr: left/right eigenvectors
c             alamda: eigenvalues
c     Issues: Do we use the entropy fix? Spoils symmetry. For now
c             we don't use it.
c             Do we use eigenvalue degeneracy fixes (a la Crockett)?
c             These are commented out.
c
c-----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
c     
      double precision :: utilde(INLO:INHI,nvar)
      double precision :: evl(INLO:INHI,nvar,nvar)
      double precision :: evr(INLO:INHI,nvar,nvar)
      double precision :: alamda(INLO:INHI,nvar)
      integer          :: ilo,ihi
c     
c     Local variables
      integer          :: i,j,k,l,m
      double precision :: uabs(INLO:INHI)
      double precision :: rhoInv(INLO:INHI)
      double precision :: Asq(INLO:INHI)
      double precision :: axsq(INLO:INHI)
      double precision :: csndsq(INLO:INHI)
      double precision :: cfast(INLO:INHI)
      double precision :: cslow(INLO:INHI)
      double precision :: tmp(INLO:INHI)
      logical          :: Mask(INLO:INHI)
      double precision :: rho,u,v,w,p,c,eint
      double precision :: bx,by,bz
      double precision :: alphaS(INLO:INHI),alphaF(INLO:INHI)
      double precision :: betaY(INLO:INHI),betaZ(INLO:INHI)
      double precision :: betaPerpMagInv(INLO:INHI)
      double precision :: betaPerpMag(INLO:INHI)
      double precision :: Sbx(INLO:INHI)
      double precision :: evec(nvar,nvar)
      double precision :: gamm
c     
      double precision :: TwoSqrtInv
      double precision, parameter :: epsilon=1.0D-8
      double precision, parameter :: tolerance=1.0D-8
      integer                     :: ii,direction
c     utilde(1)=rho, 
c     utilde(2)=u
c     utilde(3)=v
c     utilde(4)=w
c     utilde(5)=Bx
c     utilde(6)=By
c     utilde(7)=Bz
c     utilde(8)=p
c     

      TwoSqrtInv = 1.0D0/dsqrt(2.0D0)
      betaPerpMag = 0.0D0
      Mask = .false.
      Mask(ilo:ihi) = .true.
      where(mask)
c     
         uabs = abs(utilde(:,2))
         rhoInv(:) = 1.0D0/utilde(:,1)
c     
         axsq = utilde(:,5)**2*rhoInv(:)
         Asq = (utilde(:,5)**2+utilde(:,6)**2+utilde(:,7)**2)*rhoInv
         csndsq = gamma*utilde(:,8)*rhoInv
         tmp = dsqrt(dmax1(((csndsq+Asq)**2-4.0D0*csndsq*axsq),0.0D0))
         cfast = 0.5D0*(csndsq+Asq+tmp)
         cslow = 0.5D0*(csndsq+Asq-tmp)
c     
         alphaS = dsqrt((dmax1(cfast-csndsq,0.0D0))/(cfast-cslow))
         alphaF = dsqrt((dmax1(csndsq-cslow,0.0D0))/(cfast-cslow))
c     
         Sbx = dsign(1.0D0,utilde(:,5))
c     
         betaPerpMag = dsqrt(utilde(:,6)**2+utilde(:,7)**2)
c     
         betaY = TwoSqrtInv
         betaZ = TwoSqrtInv
      endwhere
      where(betaPerpMag.lt.tolerance)
         betaPerpMagInv = 1.0D0
         betaY = TwoSqrtInv
         betaZ = TwoSqrtInv
      elsewhere
         betaPerpMagInv = 1.0D0/betaPerpMag
         betaY = utilde(:,6)*betaPerpMagInv
         betaZ = utilde(:,7)*betaPerpMagInv
      endwhere
c     
      evl = 0.0D0
      evr = 0.0D0
      where(mask)
         cfast = dsqrt(cfast)
	 cslow = dsqrt(dmax1(cslow,0.0D0))
         axsq = sqrt(axsq)
c     
         alamda(:,1) = utilde(:,2)
         alamda(:,2) = utilde(:,2)
         alamda(:,3) = (utilde(:,2)+axsq)
         alamda(:,4) = (utilde(:,2)-axsq)
         alamda(:,5) = (utilde(:,2)+cfast)
         alamda(:,6) = (utilde(:,2)-cfast)
         alamda(:,7) = (utilde(:,2)+cslow)
         alamda(:,8) = (utilde(:,2)-cslow)
c
c        Eigenvalue degeneracies:
c$$$         if(axsq(i).lt.epsilon*cfast(i)) then
c$$$            alamda(i,3) = alamda(i,2)
c$$$            alamda(i,4) = alamda(i,2)
c$$$            alamda(i,7) = alamda(i,2)
c$$$            alamda(i,8) = alamda(i,2)
c$$$         endif
c$$$         if(cslow(i).lt.epsilon*cfast(i)) then
c$$$            alamda(i,7) = alamda(i,2)
c$$$            alamda(i,8) = alamda(i,2)
c$$$         endif
c
c
c     Eigenvector corresponding to u
         evr(:,1,1) = 1.0D0
c     
c     Eigenvector corresponding to u
         evr(:,5,2) = 1.0D0
c     
c     Eigenvector corresponding to u+ax
         evr(:,3,3) = -betaZ*TwoSqrtInv
         evr(:,4,3) =  betaY*TwoSqrtInv
         evr(:,6,3) =  dsqrt(utilde(:,1))*TwoSqrtInv*betaZ*Sbx
         evr(:,7,3) = -dsqrt(utilde(:,1))*TwoSqrtInv*betaY*Sbx
c
c     Eigenvector corresponding to u-ax
        evr(:,3,4) = -betaZ*TwoSqrtInv
        evr(:,4,4) =  betaY*TwoSqrtInv
        evr(:,6,4) = -evr(:,6,3)
        evr(:,7,4) = -evr(:,7,3)
c
c     Eigenvector corresponding to u+cfast
        evr(:,1,5) =  utilde(:,1)*alphaF
        evr(:,2,5) =  alphaF*cfast
        evr(:,3,5) = -alphaS*cslow*betaY*Sbx
        evr(:,4,5) = -alphaS*cslow*betaZ*Sbx
        evr(:,6,5) =  alphaS*dsqrt(utilde(:,1)*csndsq)*betaY
        evr(:,7,5) =  alphaS*dsqrt(utilde(:,1)*csndsq)*betaZ
        evr(:,8,5) =  utilde(:,1)*alphaF*csndsq
c
c     Eigenvector corresponding to u-cfast
        evr(:,1,6) =  evr(:,1,5)
        evr(:,2,6) = -evr(:,2,5)
        evr(:,3,6) = -evr(:,3,5)
        evr(:,4,6) = -evr(:,4,5)
        evr(:,6,6) =  evr(:,6,5)
        evr(:,7,6) =  evr(:,7,5)
        evr(:,8,6) =  evr(:,8,5)
c
c     Eigenvector corresponding to u+cslow
        evr(:,1,7) =  utilde(:,1)*alphaS
        evr(:,2,7) =  alphaS*cslow
        evr(:,3,7) =  alphaF*cfast*betaY*Sbx
        evr(:,4,7) =  alphaF*cfast*betaZ*Sbx
        evr(:,6,7) = -alphaF*dsqrt(utilde(:,1)*csndsq)*betaY
        evr(:,7,7) = -alphaF*dsqrt(utilde(:,1)*csndsq)*betaZ
        evr(:,8,7) =  utilde(:,1)*alphaS*csndsq
c
c
c     Eigenvector corresponding to u-cslow
        evr(:,1,8) =  evr(:,1,7)
        evr(:,2,8) = -evr(:,2,7)
        evr(:,3,8) = -evr(:,3,7)
        evr(:,4,8) = -evr(:,4,7)
        evr(:,6,8) =  evr(:,6,7)
        evr(:,7,8) =  evr(:,7,7)
        evr(:,8,8) =  evr(:,8,7)

c
c     Eigenvector corresponding to u
        evl(:,1,1) =  1.0D0
        evl(:,1,8) = -1.0D0/csndsq
c
c     Eigenvector corresponding to u
        evl(:,2,5) = 1.0D0
c
c     Eigenvector corresponding to u+ax
        evl(:,3,3) = -betaZ*TwoSqrtInv
        evl(:,3,4) =  betaY*TwoSqrtInv
        evl(:,3,6) =  1.0D0/dsqrt(utilde(:,1))*TwoSqrtInv*betaZ*SBx
        evl(:,3,7) = -1.0D0/dsqrt(utilde(:,1))*TwoSqrtInv*betaY*SBx
c
c     Eigenvector corresponding to u-ax
        evl(:,4,3) = -betaZ*TwoSqrtInv
        evl(:,4,4) =  betaY*TwoSqrtInv
        evl(:,4,6) = -evl(:,3,6)
        evl(:,4,7) = -evl(:,3,7)
c
c     Eigenvector corresponding to u+cfast
        evl(:,5,1) =  0.0D0
        evl(:,5,2) =  0.5D0*alphaF*cfast/csndsq
        evl(:,5,3) = -0.5D0*alphaS*cslow*betaY*Sbx/csndsq
        evl(:,5,4) = -0.5D0*alphaS*cslow*betaZ*Sbx/csndsq
        evl(:,5,6) =  0.5D0*alphaS/dsqrt(utilde(:,1)*csndsq)*betaY
        evl(:,5,7) =  0.5D0*alphaS/dsqrt(utilde(:,1)*csndsq)*betaZ
        evl(:,5,8) =  0.5D0*alphaF/utilde(:,1)/csndsq
c
c     Eigenvector corresponding to u-cfast
        evl(:,6,1) =  evl(:,5,1)
        evl(:,6,2) = -evl(:,5,2)
        evl(:,6,3) = -evl(:,5,3)
        evl(:,6,4) = -evl(:,5,4)
        evl(:,6,6) =  evl(:,5,6)
        evl(:,6,7) =  evl(:,5,7)
        evl(:,6,8) =  evl(:,5,8)
c
c     Eigenvector corresponding to u+cslow
        evl(:,7,1) =  0.0D0
        evl(:,7,2) =  0.5D0*alphaS*cslow/csndsq
        evl(:,7,3) =  0.5D0*alphaF*cfast*betaY*Sbx/csndsq
        evl(:,7,4) =  0.5D0*alphaF*cfast*betaZ*Sbx/csndsq
        evl(:,7,6) = -0.5D0*alphaF/dsqrt(utilde(:,1)*csndsq)*betaY
        evl(:,7,7) = -0.5D0*alphaF/dsqrt(utilde(:,1)*csndsq)*betaZ
        evl(:,7,8) =  0.5D0*alphaS/utilde(:,1)/csndsq
c
c
c     Eigenvector corresponding to u-cslow
        evl(:,8,1) =  evl(:,7,1)
        evl(:,8,2) = -evl(:,7,2)
        evl(:,8,3) = -evl(:,7,3)
        evl(:,8,4) = -evl(:,7,4)
        evl(:,8,6) =  evl(:,7,6)
        evl(:,8,7) =  evl(:,7,7)
        evl(:,8,8) =  evl(:,7,8)

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
c$$$           write(6,'(1h,2(e11.4,x))') evec(k,k),
c$$$     &          alamda(i,k)
c$$$        enddo
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
#ifdef HARTEN
c     Harten's entropy fix
c     HARTEN'S correction for transonic case.
      do l=1,nvar
         where(alamda(:,l).lt.epsilon)
            alamda(:,l) = sign(alamda(:,l),1.D0)*
     &           0.5D0*(alamda(:,l)**2
     &           +epsilon**2)/epsilon
         endwhere
      enddo
#endif
c     

c$$$c     Check Orthogonality of eigenvectors
c$$$      i = (ilo+ihi)/2
c$$$      evec = matmul(evl(i,:,:),evr(i,:,:))
c$$$      write(6,*)
c$$$      write(6,*) 'Conservative Ev Ortho Check'
c$$$      do k=1,nvar,1
c$$$         write(6,'(1h,8(e11.4,x))') evec(k,:)
c$$$      enddo
      
      return
      end subroutine SetEigenSystem
      end module evectors_primitive
      
