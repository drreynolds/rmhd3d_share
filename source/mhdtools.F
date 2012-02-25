c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c=======================================================================

      subroutine ConservedToPrimitive(u,v)
c-----------------------------------------------------------------------
c     Description: Convert conservative to primitive variables at cell 
c        centers
c
c     Arguments:
c        u - (dbl(*), input) conservative variables
c        v - (dbl(*), output) primitive variables
c
c     Note:  Code taken directly from Ravi's routine mhdsolverk2.F, 
c        with a few updates for commenting and readability.
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties


c======= Declarations =========
      implicit none
      double precision :: u(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: v(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer          :: i, j, k, l
      double precision :: rho_inv, vi, vj, vk, bi, bj, bk, pfloor, pcorr
      integer          :: pbreach


c======= Internals ============

c     set pressure floor
      pfloor = 1.0d-5
      pcorr  = 1.0d-5
      pbreach = 0

c     iterate through spatial cells
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
            do i=1-NGHOST,NXlocal+NGHOST,1
               rho_inv=1.0D0/u(i,j,k,1)
               v(i,j,k,1)=u(i,j,k,1)
               v(i,j,k,2)=u(i,j,k,2)*rho_inv
               v(i,j,k,3)=u(i,j,k,3)*rho_inv
               v(i,j,k,4)=u(i,j,k,4)*rho_inv
               v(i,j,k,5:7)=u(i,j,k,5:7)
               vi=v(i,j,k,2)
               vj=v(i,j,k,3)
               vk=v(i,j,k,4)
               bi=v(i,j,k,5)
               bj=v(i,j,k,6)
               bk=v(i,j,k,7)
c
               v(i,j,k,8)=(gamma-1.0D0)*(u(i,j,k,8)-0.5D0*v(i,j,k,1)
     &              *(vi*vi+vj*vj+vk*vk)-0.5d0*(bi*bi+bj*bj+bk*bk))
c               if(v(i,j,k,8).lt.0.0D0) then
c                  write(6,*) 'WARNING: Negative Pressure',i,j,k,
c     &                 v(i,j,k,8),0.5D0*v(i,j,k,1)
c     &              *(vi*vi+vj*vj+vk*vk),0.5d0*(bi*bi+bj*bj+bk*bk),
c     &                 u(i,j,k,8)
c               endif
c     Arbirtrary pressure floor - if breached print warning 
c     and recommend killing the process.
               if(v(i,j,k,8) .lt. pfloor) then
                  pbreach = 1
c$$$                  write(6,'(A,3(1x,i4),4(1x,A,1x,es12.5))') 
c$$$     &                 'WARNING: Pressure floor breached at (',
c$$$     &                 i, j, k, ')\n     p =', v(i,j,k,8), ', <v,v> =', 
c$$$     &                 0.5D0*v(i,j,k,1)*(vi*vi+vj*vj+vk*vk), 
c$$$     &                 ', <B,B> =', 0.5D0*(bi*bi+bj*bj+bk*bk), 
c$$$     &                 ', e =', u(i,j,k,8)
                  v(i,j,k,8)=pcorr
                  u(i,j,k,8)=v(i,j,k,8)/(gamma-1.0D0)
     &                 + 0.5D0*v(i,j,k,1)*(vi*vi+vj*vj+vk*vk)
     &                 + 0.5D0*(bi*bi+bj*bj+bk*bk)
               endif
            enddo
         enddo
      enddo
      if (pbreach == 1) then
         write(6,*)  'WARNING: Pressure floor breached'
      endif

c
      return
      end subroutine ConservedToPrimitive
c=======================================================================




      subroutine BfaceDivergenceF(divb,vface,idir)
c-----------------------------------------------------------------------
c     Description: Compute portion of divB using normal
c                  magnetic field components on faces in direction idir
c     Input: vface - face-centered primitive vars in idir-direction
c     Ouput: divb 
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties

c======= Declarations =========
      implicit none

      double precision :: divb(NXLSIZE,NYLSIZE,NZLSIZE)
      double precision :: vface(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer          :: idir
      integer          :: i,j,k
      integer          :: ioff,joff,koff


c======= Internals ============

      ioff=0
      joff=0
      koff=0
      if(idir.eq.1) ioff=1
      if(idir.eq.2) joff=1
      if(idir.eq.3) koff=1

c     BCs at faces (need to revisit - must have all BCs specified separately
c     from the solver code
      if(idir.eq.1) then
         if(iprocx.eq.1) then
            vface(0,:,:,2)=0.D0
            vface(0,:,:,5)=0.D0
         endif
         if(iprocx.eq.XPROCS) then
            vface(nxlocal,:,:,2)=0.D0
            vface(nxlocal,:,:,5)=0.D0
         endif
      endif

      if(idir.eq.2) then
         if(iprocy.eq.1) then
            vface(:,0,:,5)=0.D0
            vface(:,0,:,2)=0.D0
         endif
         if(iprocy.eq.XPROCS) then
            vface(:,nylocal,:,2)=0.D0
            vface(:,nylocal,:,5)=0.D0
         endif
      endif

      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               divb(i,j,k)=divb(i,j,k)+
     &              (vface(i,j,k,5)
     &              -vface(i-ioff,j-joff,k-koff,5))/dx
            enddo
         enddo
      enddo

      return
      end subroutine BfaceDivergenceF
c=======================================================================





      subroutine BFaceDivergenceF2(divb,ux)
c-----------------------------------------------------------------------
c     Description: Compute portion of divB using normal
c                  magnetic field components on faces which are
c                  obtained by averaging the cell centered values
c     Input: ux - cell-centered conservative vars
c     Ouput: divb 
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      
c======= Declarations =========
      double precision :: divb(NXLSIZE,NYLSIZE,NZLSIZE)
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer          :: i,j,k
c     
      divb=0.D0
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               divb(i,j,k)=
     &              (ux(i+1,j,k,5)
     &              -ux(i-1,j,k,5)+
     &              ux(i,j+1,k,6)
     &              -ux(i,j-1,k,6))/2.D0/dx
c               write(6,*) 'DIVB=',i,j,divb(i,j,k)
            enddo
         enddo
      enddo
      return
      end subroutine BFaceDivergenceF2
c=======================================================================





      subroutine DivSource(ux,ux0,dt)
c-----------------------------------------------------------------------
c     Description: Non-conservative source which uses divB based on 
c                  cell centered B
c     Arguments:
c              ux - (out) new fluid variables over entire mesh
c              dt - (in) current time step
c             ux0 - (in) previous fluid variables
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: ux0(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: dt
     
c     Local variables
c       divergB - divergence of B
      double precision :: divergB(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: rho_inv, vi,vj,vk,bi,bj,bk
      double precision :: dxi,dyi,dzi
      integer          :: i,j,k,l


c======= Internals ============

      dxi=0.5D0/dx
      dyi=0.5D0/dy
      dzi=0.5D0/dz
      divergB=0.D0

#ifdef TWO_D
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               divergB(i,j,k)=dxi*(ux0(i+1,j,k,5)-ux0(i-1,j,k,5))+
     &              dyi*(ux0(i,j+1,k,6)-ux0(i,j-1,k,6))
            enddo
         enddo
      enddo
c      write(6,*) 'DIV', maxval(divergB), minval(divergB)
#else
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               divergB(i,j,k)=dxi*(ux0(i+1,j,k,5)-ux0(i-1,j,k,5))+
     &              dyi*(ux0(i,j+1,k,6)-ux0(i,j-1,k,6))+
     &              dzi*(ux0(i,j,k+1,7)-ux0(i,j,k-1,7))
            enddo
         enddo
      enddo
#endif


      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
c               bi=ux(i,j,k,5); bj=ux(i,j,k,6); bk=ux(i,j,k,7)
c               vi=ux(i,j,k,2)/ux(i,j,k,1)
c               vj=ux(i,j,k,3)/ux(i,j,k,1)
c               vk=ux(i,j,k,4)/ux(i,j,k,1)
               bi=ux0(i,j,k,5); bj=ux0(i,j,k,6); bk=ux0(i,j,k,7)
               vi=ux0(i,j,k,2)/ux0(i,j,k,1)
               vj=ux0(i,j,k,3)/ux0(i,j,k,1)
               vk=ux0(i,j,k,4)/ux0(i,j,k,1)
               ux(i,j,k,2)=ux(i,j,k,2)-dt*divergB(i,j,k)*bi
               ux(i,j,k,3)=ux(i,j,k,3)-dt*divergB(i,j,k)*bj
#ifdef TWOHALF_D
               ux(i,j,k,4)=ux(i,j,k,4)-dt*divergB(i,j,k)*bk
#endif
               ux(i,j,k,5)=ux(i,j,k,5)-dt*divergB(i,j,k)*vi
               ux(i,j,k,6)=ux(i,j,k,6)-dt*divergB(i,j,k)*vj
#ifdef TWOHALF_D
               ux(i,j,k,7)=ux(i,j,k,7)-dt*divergB(i,j,k)*vk
#endif
               ux(i,j,k,8)=ux(i,j,k,8)
     &              -dt*divergB(i,j,k)*(vi*bi+vj*bj+vk*bk)
            enddo
         enddo
      enddo

      return
      end subroutine DivSource
c=================================================================




      subroutine DivBDiffusionFilter(ux)
c-----------------------------------------------------------------------
c     Description: diffuse away divergence; this has a large stabilizing
c                  effect on the code. Diffusing away divB implies 
c                  a certain filter on B noted below.
c
c     Arguments:
c              ux - (inout) conserved variables before/after filter
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties

c======= Declarations =========
      implicit none

      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: bnew(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,3)
      double precision :: rho_inv, vi, vj, vk, bi, bj, bk
      double precision :: dxi, dyi, dzi
      double precision :: dxi2, dyi2, dzi2, dxidyi
      double precision :: lamda
      integer          :: i, j, k, l

c======= Internals ============

      lamda = dx*dx/8.D0

      dxi  = 0.5D0/dx
      dyi  = 0.5D0/dy
      dzi  = 0.5D0/dz
      dxi2 = 1.D0/dx/dx
      dyi2 = 1.D0/dy/dy
      dzi2 = 1.D0/dz/dz
      dxidyi = 0.25D0/dx/dy

      bnew(:,:,:,1) = ux(:,:,:,5)
      bnew(:,:,:,2) = ux(:,:,:,6)
      bnew(:,:,:,3) = ux(:,:,:,7)
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               bnew(i,j,k,1)=ux(i,j,k,5)+lamda*(
     &              dxi2*(ux(i+1,j,k,5)+ux(i-1,j,k,5)-2.D0*ux(i,j,k,5))
     &              +dxidyi*(ux(i+1,j+1,k,6)+ux(i-1,j-1,k,6)-
     &              ux(i-1,j+1,k,6)-ux(i+1,j-1,k,6)))
               bnew(i,j,k,2)=ux(i,j,k,6)+lamda*(
     &              dyi2*(ux(i,j+1,k,6)+ux(i,j-1,k,6)-2.D0*ux(i,j,k,6))
     &              +dxidyi*(ux(i+1,j+1,k,5)+ux(i-1,j-1,k,5)-
     &              ux(i-1,j+1,k,5)-ux(i+1,j-1,k,5)))
            enddo
         enddo
      enddo

c     Is it advisable to readjust total energy? Don't think so. 
c     A place holder for this is below -just in case we need it later.
c$$$      call ConservedEnergyCorrection(ux(ixlo,iylo,izlo,8),
c$$$     &       bnew,ux(ixlo,iylo,izlo,5))
      ux(1:nxlocal,1:nylocal,1:nzlocal,5)
     &     =bnew(1:nxlocal,1:nylocal,1:nzlocal,1)
      ux(1:nxlocal,1:nylocal,1:nzlocal,6)
     &     =bnew(1:nxlocal,1:nylocal,1:nzlocal,2)

      return
      end subroutine DivBDiffusionFilter
c=================================================================




      subroutine RemapToLocalCoords(vx,ilo,ihi,direction)
c-----------------------------------------------------------------------
c     Description: Map to local directions
c                  x: u,v,w ---> u,v,w
c                  y: u,v,w ---> w,u,v
c                  z: u,v,w ---> v,w,u
c     Input: vx
c     Output: vx after mapping
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none
      double precision:: vx(INLO:INHI,nvar)
      double precision:: tmp(INLO:INHI,3)
      integer:: ilo, ihi, direction

c======= Internals ============

      select case (direction)
      case (2)
         tmp(ilo:ihi,1) = vx(ilo:ihi,3)
         tmp(ilo:ihi,2) = vx(ilo:ihi,4)
         tmp(ilo:ihi,3) = vx(ilo:ihi,2)
         vx(ilo:ihi,2:4) = tmp(ilo:ihi,1:3)

         tmp(ilo:ihi,1) = vx(ilo:ihi,6)
         tmp(ilo:ihi,2) = vx(ilo:ihi,7)
         tmp(ilo:ihi,3) = vx(ilo:ihi,5)
         vx(ilo:ihi,5:7) = tmp(ilo:ihi,1:3)
      case (3)
         tmp(ilo:ihi,1) = vx(ilo:ihi,4)
         tmp(ilo:ihi,2) = vx(ilo:ihi,2)
         tmp(ilo:ihi,3) = vx(ilo:ihi,3)
         vx(ilo:ihi,2:4) = tmp(ilo:ihi,1:3)

         tmp(ilo:ihi,1) = vx(ilo:ihi,7)
         tmp(ilo:ihi,2) = vx(ilo:ihi,5)
         tmp(ilo:ihi,3) = vx(ilo:ihi,6)
         vx(ilo:ihi,5:7) = tmp(ilo:ihi,1:3)         
      end select
c     
      return
      end subroutine RemapToLocalCoords
c=================================================================




      subroutine RemapToGlobalCoords(vx,ilo,ihi,direction)
c-----------------------------------------------------------------------
c     Description: Map to global directions
c                  x: u,v,w ---> u,v,w
c                  y: w,u,v ---> u,v,w
c                  z: v,w,u ---> u,v,w
c     Input: vx
c     Output: vx after mapping
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none
      double precision :: vx(INLO:INHI,nvar)
      double precision :: tmp(INLO:INHI,3)
      integer:: ilo, ihi, direction, i

c======= Internals ============

      select case (direction)
      case (2) 
         tmp(ilo:ihi,1) = vx(ilo:ihi,4)
         tmp(ilo:ihi,2) = vx(ilo:ihi,2)
         tmp(ilo:ihi,3) = vx(ilo:ihi,3)
         vx(ilo:ihi,2:4) = tmp(ilo:ihi,1:3)

         tmp(ilo:ihi,1) = vx(ilo:ihi,7)
         tmp(ilo:ihi,2) = vx(ilo:ihi,5)
         tmp(ilo:ihi,3) = vx(ilo:ihi,6)
         vx(ilo:ihi,5:7) = tmp(ilo:ihi,1:3)
      case (3)
         tmp(ilo:ihi,1) = vx(ilo:ihi,3)
         tmp(ilo:ihi,2) = vx(ilo:ihi,4)
         tmp(ilo:ihi,3) = vx(ilo:ihi,2)
         vx(ilo:ihi,2:4) = tmp(ilo:ihi,1:3)

         tmp(ilo:ihi,1) = vx(ilo:ihi,6)
         tmp(ilo:ihi,2) = vx(ilo:ihi,7)
         tmp(ilo:ihi,3) = vx(ilo:ihi,5)
         vx(ilo:ihi,5:7) = tmp(ilo:ihi,1:3)         
      end select
c     
      return
      end subroutine RemapToGlobalCoords
c=================================================================




      subroutine RemapToLocalCoords2(vx,ilo,ihi,direction)
c-----------------------------------------------------------------------
c     Description: Map to local directions
c                  x: u,v,w ---> u,v,w
c                  y: u,v,w ---> w,u,v
c                  z: u,v,w ---> v,w,u
c     Input: vx
c     Output: vx after mapping
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none
      double precision:: vx(nvar,INLO:INHI)
      double precision:: tmp(6)
      integer:: ilo, ihi, direction, i

c======= Internals ============

      select case (direction)
      case (2)
         do i=ilo,ihi,1
            tmp(1) = vx(3,i)
            tmp(2) = vx(4,i)
            tmp(3) = vx(2,i)
            tmp(4) = vx(6,i)
            tmp(5) = vx(7,i)
            tmp(6) = vx(5,i)
            vx(2:7,i) = tmp
         enddo
      case (3)
         do i=ilo,ihi,1
            tmp(1) = vx(4,i)
            tmp(2) = vx(2,i)
            tmp(3) = vx(3,i)
            tmp(4) = vx(7,i)
            tmp(5) = vx(5,i)
            tmp(6) = vx(6,i)
            vx(2:7,i) = tmp
         enddo
      end select

      return
      end subroutine RemapToLocalCoords2
c=================================================================




      subroutine RemapToGlobalCoords2(vx,ilo,ihi,direction)
c-----------------------------------------------------------------------
c     Description: Map to global directions
c                  x: u,v,w ---> u,v,w
c                  y: w,u,v ---> u,v,w
c                  z: v,w,u ---> u,v,w
c     Input: vx
c     Output: vx after mapping
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none
      double precision :: vx(nvar,INLO:INHI)
      double precision :: tmp(6)
      integer:: ilo, ihi, direction, i

c======= Internals ============

      select case (direction)
      case (2) 
         do i=ilo,ihi,1
            tmp(1) = vx(4,i)
            tmp(2) = vx(2,i)
            tmp(3) = vx(3,i)
            tmp(4) = vx(7,i)
            tmp(5) = vx(5,i)
            tmp(6) = vx(6,i)
            vx(2:7,i) = tmp
         enddo
      case (3)
         do i=ilo,ihi,1
            tmp(1) = vx(3,i)
            tmp(2) = vx(4,i)
            tmp(3) = vx(2,i)
            tmp(4) = vx(6,i)
            tmp(5) = vx(7,i)
            tmp(6) = vx(5,i)
            vx(2:7,i) = tmp
         enddo
      end select

      return
      end subroutine RemapToGlobalCoords2
c=================================================================




      subroutine RemapToLocalCoords3D(vx,direction)
c-----------------------------------------------------------------------
c     Description: Map to local directions
c                  x: u,v,w ---> u,v,w
c                  y: u,v,w ---> w,u,v
c                  z: u,v,w ---> v,w,u
c     Input: vx
c     Output: vx after mapping
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none
      double precision:: vx(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision:: tmp(6)
      integer:: direction, i, j, k

c======= Internals ============

      select case (direction)
      case (2)
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  tmp(1) = vx(i,j,k,3)
                  tmp(2) = vx(i,j,k,4)
                  tmp(3) = vx(i,j,k,2)
                  tmp(4) = vx(i,j,k,6)
                  tmp(5) = vx(i,j,k,7)
                  tmp(6) = vx(i,j,k,5)
                  vx(i,j,k,2) = tmp(1)
                  vx(i,j,k,3) = tmp(2)
                  vx(i,j,k,4) = tmp(3)
                  vx(i,j,k,5) = tmp(4)
                  vx(i,j,k,6) = tmp(5)
                  vx(i,j,k,7) = tmp(6)
               enddo
            enddo
         enddo
      case (3)
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  tmp(1) = vx(i,j,k,4)
                  tmp(2) = vx(i,j,k,2)
                  tmp(3) = vx(i,j,k,3)
                  tmp(4) = vx(i,j,k,7)
                  tmp(5) = vx(i,j,k,5)
                  tmp(6) = vx(i,j,k,6)
                  vx(i,j,k,2) = tmp(1)
                  vx(i,j,k,3) = tmp(2)
                  vx(i,j,k,4) = tmp(3)
                  vx(i,j,k,5) = tmp(4)
                  vx(i,j,k,6) = tmp(5)
                  vx(i,j,k,7) = tmp(6)
               enddo
            enddo
         enddo
      end select

      return
      end subroutine RemapToLocalCoords3D
c=================================================================




      subroutine RemapToGlobalCoords3D(vx,direction)
c-----------------------------------------------------------------------
c     Description: Map to global directions
c                  x: u,v,w ---> u,v,w
c                  y: w,u,v ---> u,v,w
c                  z: v,w,u ---> u,v,w
c     Input: vx
c     Output: vx after mapping
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none
      double precision :: vx(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: tmp(6)
      integer:: direction, i, j, k

c======= Internals ============

      select case (direction)
      case (2) 
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  tmp(1) = vx(i,j,k,4)
                  tmp(2) = vx(i,j,k,2)
                  tmp(3) = vx(i,j,k,3)
                  tmp(4) = vx(i,j,k,7)
                  tmp(5) = vx(i,j,k,5)
                  tmp(6) = vx(i,j,k,6)
                  vx(i,j,k,2) = tmp(1)
                  vx(i,j,k,3) = tmp(2)
                  vx(i,j,k,4) = tmp(3)
                  vx(i,j,k,5) = tmp(4)
                  vx(i,j,k,6) = tmp(5)
                  vx(i,j,k,7) = tmp(6)
               enddo
            enddo
         enddo
      case (3)
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  tmp(1) = vx(i,j,k,3)
                  tmp(2) = vx(i,j,k,4)
                  tmp(3) = vx(i,j,k,2)
                  tmp(4) = vx(i,j,k,6)
                  tmp(5) = vx(i,j,k,7)
                  tmp(6) = vx(i,j,k,5)
                  vx(i,j,k,2) = tmp(1)
                  vx(i,j,k,3) = tmp(2)
                  vx(i,j,k,4) = tmp(3)
                  vx(i,j,k,5) = tmp(4)
                  vx(i,j,k,6) = tmp(5)
                  vx(i,j,k,7) = tmp(6)
               enddo
            enddo
         enddo
      end select

      return
      end subroutine RemapToGlobalCoords3D
c=================================================================
