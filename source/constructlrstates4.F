C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2001
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------
c     $Log: constructlrstates4.F,v $
c-----------------------------------------------------------------
      module LRStates
      contains
      subroutine ConstructLRStates(v,vl,vr,ilo,ihi)
      use mesh_common
      use mesh
      use properties
c     
c      double precision:: u(inlo:inhi,nvar)
      double precision:: v(inlo:inhi,nvar)
      double precision:: vl(inlo:inhi,nvar)
      double precision:: vr(inlo:inhi,nvar)
      double precision:: dvdx(inlo:inhi,nvar)
      double precision:: vc(nvar,3)
      double precision:: vi, vj, vk
      double precision:: bi, bj, bk
      integer:: i,j,k,l,ilo,ihi
      double precision:: evl(8,8), evr(inlo:inhi,8,8), rhs(8)
      double precision:: rhoInv
      double precision:: Asq
      double precision:: axsq
      double precision:: csndsq
      double precision:: cfast
      double precision:: cslow
      double precision:: tmp
      double precision:: alphaS,alphaF
      double precision:: betaY,betaZ
      double precision:: betaPerpMagInv
      double precision:: Sbx
      double precision:: evec
      double precision:: gamm
c     
      double precision,parameter:: TwoSqrtInv=0.70710678D0
c     
c     Evaluating v should be done only once. Because the
c     calculations are the same for all the directions.
c     We need to only switch the indices for the
c     velocity when doing the Y and Z sweeps.
c     
c$$$      do i=ilo,ihi,1
c$$$         rhoInv=1.0D0/u(i,1)
c$$$         v(i,1)=u(i,1)
c$$$         v(i,2)=u(i,2)*rhoInv
c$$$         v(i,3)=u(i,3)*rhoInv
c$$$         v(i,4)=u(i,4)*rhoInv
c$$$         v(i,5)=u(i,5)
c$$$         v(i,6)=u(i,6)
c$$$         v(i,7)=u(i,7)
c$$$         vi=v(i,2)
c$$$         vj=v(i,3)
c$$$         vk=v(i,4)
c$$$         bi=v(i,5)
c$$$         bj=v(i,6)
c$$$         bk=v(i,7)
c$$$         v(i,8)=(gamma-1.D0)*(u(i,8)-0.5D0*v(i,1)
c$$$     &        *(vi*vi+vj*vj+vk*vk)-0.5*(bi*bi+bj*bj+bk*bk))
c$$$      enddo
c     
c     
c     Set the characteristic form
c     
      vl=0.D0
      vr=0.D0
      do i=ilo+1,ihi-nghost,1
         vl(i,:)=v(i,:)
         vr(i,:)=v(i+1,:)
      enddo	
      return
c
      do i=ilo+1,ihi-1,1
         evl=0.D0
         rhoInv=1.D0/v(i,1)
c     
         axsq=v(i,5)**2*rhoInv
         Asq=(v(i,5)**2+v(i,6)**2+v(i,7)**2)*rhoInv
         csndsq=gamma*v(i,8)*rhoInv
         tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*axsq),0.D0))
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
         alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
         Sbx=dsign(1.D0,v(i,5))
c     
         if(v(i,6).eq.0.D0.and.v(i,7).eq.0.D0) then
            betaPerpMagInv=1.D0
            betaY=TwoSqrtInv
            betaZ=TwoSqrtInv
         else
            betaPerpMagInv=1.0/dsqrt(v(i,6)**2+v(i,7)**2)
            betaY=v(i,6)*betaPerpMagInv
            betaZ=v(i,7)*betaPerpMagInv
         endif
c     
c         if(betaPerpMagInv.gt.1.D-08) then
c         else
c         endif
c     
         cfast=dsqrt(cfast)
         cslow=dsqrt(dmax1(cslow,0.D0))
c     
c     

c     Eigenvector corresponding to u
         evl(1,1)=1.D0
c     
c     Eigenvector corresponding to u
         evl(2,5)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evl(3,3)=0.5D0*betaZ
         evl(3,4)=-0.5D0*betaY
         evl(3,6)=-0.5D0/dsqrt(v(i,1))*betaZ*SBx
         evl(3,7)=0.5D0/dsqrt(v(i,1))*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evl(4,3)=-0.5D0*betaZ
         evl(4,4)=0.5D0*betaY
         evl(4,6)=evl(3,6)
         evl(4,7)=evl(3,7)
c     
c     Eigenvector corresponding to u+cfast
         evl(5,1)=0.D0
         evl(5,2)=0.5D0*alphaF*cfast/csndsq
         evl(5,3)=-0.5D0*alphaS*cslow*betaY*Sbx/csndsq
         evl(5,4)=-0.5D0*alphaS*cslow*betaZ*Sbx/csndsq
         evl(5,6)=0.5D0*alphaS/dsqrt(v(i,1)*csndsq)*betaY
         evl(5,7)=0.5D0*alphaS/dsqrt(v(i,1)*csndsq)*betaZ
         evl(5,8)=0.5D0*alphaF/v(i,1)/csndsq
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
         evl(7,3)=-0.5D0*alphaF*cfast*betaY*Sbx/csndsq
         evl(7,4)=-0.5D0*alphaF*cfast*betaZ*Sbx/csndsq
         evl(7,6)=0.5D0*alphaF/dsqrt(v(i,1)*csndsq)*betaY
         evl(7,7)=0.5D0*alphaF/dsqrt(v(i,1)*csndsq)*betaZ
         evl(7,8)=0.5D0*alphaS/v(i,1)/csndsq
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
         rhs=v(i,:)
         vc(:,2)=0.D0
c         do l=1,nvar,1
c            vc(l,2)=vc(l,2)+sum(evl(l,:)*rhs(:))
c         enddo
         vc(:,2)=matmul(evl,rhs)
c     
         rhs=v(i+1,:)
         vc(:,3)=matmul(evl,rhs)
c         vc(:,3)=0.D0
c         do l=1,nvar,1
c            vc(l,3)=vc(l,3)+sum(evl(l,:)*rhs(:))
c         enddo
c     
         rhs=v(i-1,:)
         vc(:,1)=matmul(evl,rhs)
c         vc(:,1)=0.D0
c         do l=1,nvar,1
c            vc(l,1)=vc(l,1)+sum(evl(l,:)*rhs(:))
c         enddo
c     
c     
c     Compute the slopes and use a slope limiter
c     which is more lenient than the minmod slope limiter.
         do l=1,nvar,1
            dvdx(i,l)=0.5D0*dmin1( dabs(vc(l,3)-vc(l,2))/(1
     &           .0D0),dabs(vc(l,2)-vc(l,1)  )/(1.0D0))*dsign(1.D0
     &           ,vc(l,3)-vc(l,1))*dmax1(dsign(1.D0,(vc(l,3)-vc(l
     &           ,2))*(vc(l,2)-vc(l,1))),0.0D0)
c            dvdx(i,l)=0.D0
         enddo
      enddo
c     
c     Get the left and right states	
c     The left states are vl(i)=v(i)+dvdx(i)*evr(i)
c     The right states are vr(i)=v(i+1)-dvdx(i+1)*evr(i+1)
c     
c     Set up the eigenvectors for all points
      do i=ilo+1,ihi-1,1
c     
        rhoInv=1.D0/v(i,1)
c     
         axsq=v(i,5)**2*rhoInv
         Asq=(v(i,5)**2+v(i,6)**2+v(i,7)**2)*rhoInv
         csndsq=gamma*v(i,8)*rhoInv
         tmp=dsqrt(dmax1(((csndsq+Asq)**2-4.0*csndsq*axsq),0.D0))
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
         alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
         Sbx=dsign(1.D0,v(i,5))
c     
         if(v(i,6).eq.0.D0.and.v(i,7).eq.0.D0) then
            betaPerpMagInv=1.D0
            betaY=TwoSqrtInv
            betaZ=TwoSqrtInv
         else
            betaPerpMagInv=1.0/dsqrt(v(i,6)**2+v(i,7)**2)
            betaY=v(i,6)*betaPerpMagInv
            betaZ=v(i,7)*betaPerpMagInv
         endif
c     
         cfast=dsqrt(cfast)
         cslow=dsqrt(dmax1(cslow,0.D0))
c
         evr(i,:,:)=0.D0
c     Eigenvector corresponding to u
         evr(i,1,1)=1.D0
         evr(i,8,1)=1.D0/csndsq
c     
c     Eigenvector corresponding to u
         evr(i,5,2)=1.D0
c     
c     Eigenvector corresponding to u+ax
         evr(i,3,3)=betaZ
         evr(i,4,3)=-betaY
         evr(i,6,3)=-dsqrt(v(i,1))*betaZ*SBx
         evr(i,7,3)=dsqrt(v(i,1))*betaY*SBx
c     
c     Eigenvector corresponding to u-ax
         evr(i,3,4)=-betaZ
         evr(i,4,4)=betaY
         evr(i,6,4)=evr(i,6,3)
         evr(i,7,4)=evr(i,7,3)
c     
c     Eigenvector corresponding to u+cfast
         evr(i,1,5)=v(i,1)*alphaF
         evr(i,2,5)=alphaF*cfast
         evr(i,3,5)=-alphaS*cslow*betaY*Sbx
         evr(i,4,5)=-alphaS*cslow*betaZ*Sbx
         evr(i,6,5)=alphaS*dsqrt(v(i,1)*csndsq)*betaY
         evr(i,7,5)=alphaS*dsqrt(v(i,1)*csndsq)*betaZ
         evr(i,8,5)=v(i,1)*alphaF*csndsq
c     
c     Eigenvector corresponding to u-cfast
         evr(i,1,6)=evr(i,1,5)
         evr(i,2,6)=-evr(i,2,5)
         evr(i,3,6)=-evr(i,3,5)
         evr(i,4,6)=-evr(i,4,5)
         evr(i,6,6)=evr(i,6,5)
         evr(i,7,6)=evr(i,7,5)
         evr(i,8,6)=evr(i,8,5)
c     
c     Eigenvector corresponding to u+cslow
         evr(i,1,7)=v(i,1)*alphaS
         evr(i,2,7)=alphaS*cslow
         evr(i,3,7)=-alphaF*cfast*betaY*Sbx
         evr(i,4,7)=-alphaF*cfast*betaZ*Sbx
         evr(i,6,7)=alphaF*dsqrt(v(i,1)*csndsq)*betaY
         evr(i,7,7)=alphaF*dsqrt(v(i,1)*csndsq)*betaZ
         evr(i,8,7)=v(i,1)*alphaS*csndsq
c     
c     
c     Eigenvector corresponding to u-cslow
         evr(i,1,8)=evr(i,1,7)
         evr(i,2,8)=-evr(i,2,7)
         evr(i,3,8)=-evr(i,3,7)
         evr(i,4,8)=-evr(i,4,7)
         evr(i,6,8)=evr(i,6,7)
         evr(i,7,8)=evr(i,7,7)
         evr(i,8,8)=evr(i,8,7)
c     
      enddo
c
      do i=ilo+1,ihi-nghost,1
c
         rhs=dvdx(i,:)
         vl(i,:)=v(i,:)+matmul(evr(i,:,:),rhs)
c         do l=1,nvar,1
c            vl(i,l)=v(i,l)+sum(evr(i,l,:)*rhs(:))
c         enddo
         rhs=dvdx(i+1,:)
         vr(i,:)=v(i+1,:)-matmul(evr(i+1,:,:),rhs)
c         do l=1,nvar,1
c            vr(i,l)=v(i+1,l)-sum(evr(i+1,l,:)*rhs(:))
c         enddo
c     Additional Dissipation
         if(vr(i,8).le.0.D0) vr(i,:)=v(i+1,:)
         if(vl(i,8).le.0.D0) vl(i,:)=v(i,:)
      enddo
c$$$         i=0
c$$$      	 write(6,*) 'V',i,v(i,:)
c$$$         i=ihi-2
c$$$      	 write(6,*) 'V',i,v(i,:)
c$$$         i=1
c$$$      	 write(6,*) 'V',i,v(i,:)
c$$$         i=ihi-1
c$$$      	 write(6,*) 'V',i,v(i,:)
c$$$         i=2
c$$$      	 write(6,*) 'V',i,v(i,:)
c$$$         i=ihi
c$$$      	 write(6,*) 'V',i,v(i,:)
c$$$         i=1
c$$$      	 write(6,*) 'dvdx',i,dvdx(i,:)
c$$$         i=ihi-1
c$$$      	 write(6,*) 'dvdx',i,dvdx(i,:)
c$$$c
c$$$	 i=1
c$$$         do k=1,8,1
c$$$            write(6,1000) evr(i,k,:)
c$$$         enddo
c$$$         write(6,*)
c$$$	 i=ihi-1
c$$$         do k=1,8,1
c$$$            write(6,1000) evr(i,k,:)
c$$$         enddo
c$$$         write(6,*)
c$$$         i=0
c$$$      	 write(6,*) 'VR',i,VR(i,:)
c$$$         i=ihi-2
c$$$      	 write(6,*) 'VR',i,VR(i,:)
c$$$ 1000   format(1h ,8(e11.4,x))
c     
c      i=(ilo+ihi)/2
c      i=ihi-1
c      call CheckOrthogonality(v(i,1),v(i,2),v(i,3),v(i,4),
c     &	   v(i,5),v(i,6),v(i,7),v(i,8))
      return
      end subroutine ConstructLRStates
      end module LRStates
c-----------------------------------------------------------------------
	Subroutine CheckOrthogonality(v1,v2,v3,v4,v5,v6,v7,v8)
        use properties
	double precision:: v1,v2,v3,v4,v5,v6,v7,v8
	double precision:: evec(8,8),evl(8,8), evr(8,8)
        double precision:: rhoInv		       
        double precision:: Asq
        double precision:: axsq
        double precision:: csndsq
        double precision:: cfast
        double precision:: cslow
        double precision:: tmp
        double precision:: alphaS,alphaF
        double precision:: betaY,betaZ
        double precision:: betaPerpMagInv
        double precision:: Sbx
        integer:: k,l
c     
        double precision,parameter:: TwoSqrtInv=0.70710678
c
         evl=0.D0
         rhoInv=1.D0/v1
c     
         axsq=v5**2*rhoInv
         Asq=(v5**2+v6**2+v7**2)*rhoInv
         csndsq=gamma*v8*rhoInv
         tmp=dsqrt((csndsq+Asq)**2-4.0*csndsq*axsq)
         cfast=0.5D0*(csndsq+Asq+tmp)
         cslow=0.5D0*(csndsq+Asq-tmp)
c     
         alphaF=dsqrt(dmax1(csndsq-cslow,0.D0)/(cfast-cslow))
         alphaS=dsqrt(dmax1(cfast-csndsq,0.D0)/(cfast-cslow))
c     
         Sbx=dsign(1.D0,v5)
c     
         if(v6.eq.0.D0.and.v7.eq.0.D0) then
            betaPerpMagInv=1.0
            betaY=TwoSqrtInv
            betaZ=TwoSqrtInv
         else
            betaPerpMagInv=1.0/dsqrt(v6**2+v7**2)
            betaY=v6*betaPerpMagInv
            betaZ=v7*betaPerpMagInv
         endif
c     
c         if(betaPerpMagInv.gt.1.D-08) then
c         else
c         endif
c     
         cfast=dsqrt(cfast)
         cslow=dsqrt(cslow)
c     
c     

c     Eigenvector corresponding to u
         evl(1,1)=1.0
         evl(1,8)=-1.0/csndsq
c     
c     Eigenvector corresponding to u
         evl(2,5)=1.0
c     
c     Eigenvector corresponding to u+ax
         evl(3,3)=-betaZ*TwoSqrtInv
         evl(3,4)=betaY*TwoSqrtInv
         evl(3,6)=1.D0/dsqrt(v1)*TwoSqrtInv*betaZ
         evl(3,7)=-1.D0/dsqrt(v1)*TwoSqrtInv*betaY
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
         evl(5,6)=0.5D0*alphaS/dsqrt(v1*csndsq)*betaY
         evl(5,7)=0.5D0*alphaS/dsqrt(v1*csndsq)*betaZ
         evl(5,8)=0.5D0*alphaF/v1/csndsq
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
         evl(7,6)=-0.5D0*alphaF/dsqrt(v1*csndsq)*betaY
         evl(7,7)=-0.5D0*alphaF/dsqrt(v1*csndsq)*betaZ
         evl(7,8)=0.5D0*alphaS/v1/csndsq
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

         evr=0.0
c     Eigenvector corresponding to u
         evr(1,1)=1.0
c     
c     Eigenvector corresponding to u
         evr(5,2)=1.0
c     
c     Eigenvector corresponding to u+ax
         evr(3,3)=-betaZ*TwoSqrtInv
         evr(4,3)=betaY*TwoSqrtInv
         evr(6,3)=dsqrt(v1)*TwoSqrtInv*betaZ
         evr(7,3)=-dsqrt(v1)*TwoSqrtInv*betaY
c     
c     Eigenvector corresponding to u-ax
         evr(3,4)=-betaZ*TwoSqrtInv
         evr(4,4)=betaY*TwoSqrtInv
         evr(6,4)=-evr(6,3)
         evr(7,4)=-evr(7,3)
c     
c     Eigenvector corresponding to u+cfast
         evr(1,5)=v1*alphaF
         evr(2,5)=alphaF*cfast
         evr(3,5)=-alphaS*cslow*betaY*Sbx
         evr(4,5)=-alphaS*cslow*betaZ*Sbx
         evr(6,5)=alphaS*dsqrt(v1*csndsq)*betaY
         evr(7,5)=alphaS*dsqrt(v1*csndsq)*betaZ
         evr(8,5)=v1*alphaF*csndsq
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
         evr(1,7)=v1*alphaS
         evr(2,7)=alphaS*cslow
         evr(3,7)=alphaF*cfast*betaY*Sbx
         evr(4,7)=alphaF*cfast*betaZ*Sbx
         evr(6,7)=-alphaF*dsqrt(v1*csndsq)*betaY
         evr(7,7)=-alphaF*dsqrt(v1*csndsq)*betaZ
         evr(8,7)=v1*alphaS*csndsq
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
c
        evec=matmul(evl,evr)
        write(6,*)
        write(6,*) 'Eigenvector Orthogonality Check'
        do k=1,8,1
           do l=1,8,1
              if(abs(evec(k,l)).lt.1.D-14) evec(l,k)=0.D0
           enddo
        enddo
        do k=1,8,1
           write(6,1000) evec(k,:)
        enddo
 1000   format(1h ,8(e11.4,x))
        return
        end  subroutine CheckOrthogonality
