C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
      module slamdas
      contains
      subroutine SetSLS(utilde,sls,vl,vr,ilo,ihi,direction)
      use mesh
      use mesh_common
      use udiff
      use alphas
#ifdef PRIMITIVE
      use evectors_primitive
#else
      use evectors_conservative
#endif
      
      double precision:: vl(INLO:INHI,nvar)
      double precision:: vr(INLO:INHI,nvar)
      double precision:: sls(INLO:INHI,nvar)
      double precision:: utilde(INLO:INHI,nvar)
      double precision:: alamda(INLO:INHI,nvar)
      double precision:: evl(INLO:INHI,nvar,nvar)
      double precision:: evr(INLO:INHI,nvar,nvar)
      double precision:: alpha(INLO:INHI,nvar)
      double precision:: durl(INLO:INHI,nvar)
      integer:: i,j,k,l,m,ilo,ihi,direction,ii
c     
      sls=0.D0
      alamda=0.D0
      alpha=0.D0
      durl=0.D0
c     
c     write(6,*) 'SETSLS', indx,indy,indz
      call SetEigenSystem(utilde,evl,evr,alamda,ilo,ihi,direction)
      call SetDurl(durl,vl,vr,ilo,ihi)
      call SetAlphas(utilde,alpha,durl,evl,evr,ilo,ihi)
c     
c     Note that evr are stored as column vectors
c     evl are stored as row vector
c     
      do i=ilo,ihi,1
         do m=1,nvar,1
            sls(i,m)=0.D0
            do l=1,nvar,1
               sls(i,m)=sls(i,m)+alpha(i,l)*
     &              alamda(i,l)*evr(i,m,l)
            enddo
         enddo
      enddo
c     do i=ilo,ihi,1
c     write(6,*) 'SLS',i,sls(i,5)
c     sls(i,5)=0.D0
c     enddo
c$$$      if(direction.eq.1) then
c$$$         do i=ilo,ihi,1
c$$$            write(81,*) i,sls(i,:)
c$$$            write(82,*) i,alpha(i,:)
c$$$            write(83,*) i,alamda(i,:)
c$$$         enddo
c$$$         write(81,*)
c$$$         write(82,*)
c$$$         write(83,*)
c$$$      endif
      return
      end subroutine SetSLS
      end module slamdas
      




