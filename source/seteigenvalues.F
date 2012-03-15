c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney 
c     Copyright 2005
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------
c     $Log: seteigenvalues.F,v $
c-----------------------------------------------------------------
      module evaluesMax
      contains
      subroutine SetEigenValueMax(utilde,alamda,ilo,ihi,direction)
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
      implicit none
      double precision :: utilde(INLO:INHI,nvar)
      double precision :: alamda(INLO:INHI)
      integer          :: ilo,ihi,direction
c     
c     Local variables
      integer          :: i
      double precision :: rhoInv
      double precision :: Asq
      double precision :: axsq
      double precision :: csndsq
      double precision :: cfast
      double precision :: rho,u,v,w,p,c,eint
      double precision :: bx,by,bz
      double precision :: tmp, gamm
c     
c     utilde(1)=rho, 
c     utilde(2)=u
c     utilde(3)=v
c     utilde(4)=w
c     utilde(5)=Bx
c     utilde(6)=By
c     utilde(7)=Bz
c     utilde(8)=p
c     

      do i=ilo,ihi,1
c     
         rhoInv = 1.0D0/utilde(i,1)
c     
         axsq = utilde(i,5)**2*rhoInv
         Asq = (utilde(i,5)**2+utilde(i,6)**2+utilde(i,7)**2)*rhoInv
         csndsq = gamma*utilde(i,8)*rhoInv
         tmp = dsqrt(dmax1(((csndsq+Asq)**2-4.0D0*csndsq*axsq),0.0D0))
         cfast = dsqrt(0.5D0*(csndsq+Asq+tmp))
c     
         axsq = sqrt(axsq)
c     
         alamda(i) = max(dabs(utilde(i,2)+cfast), 
     &        dabs(utilde(i,2)-cfast))
      enddo
      return
      end subroutine SetEigenValueMax
      end module evaluesMax
      
