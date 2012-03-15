C -*- Mode: Fortran; -*-
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

c$$$      module TVDCOEFFS
c$$$      save
c$$$      double precision:: tvdCoeff1(3)
c$$$      double precision:: tvdCoeff2(3)
c$$$      double precision:: tvdCoeff3(3)
c$$$      end module TVDCOEFFS
c$$$c-----------------------------------------------------------------------
c$$$        subroutine SetTVDCoeffs
c$$$        use TVDCOEFFS
c$$$        tvdCoeff1(1)=1.D0
c$$$        tvdCoeff1(2)=0.D0
c$$$        tvdCoeff1(3)=1.D0
c$$$        tvdCoeff2(1)=0.75D0
c$$$        tvdCoeff2(2)=0.25D0
c$$$        tvdCoeff2(3)=0.25D0
c$$$        tvdCoeff3(1)=1.D0/3.D0
c$$$        tvdCoeff3(2)=2.D0/3.D0
c$$$        tvdCoeff3(3)=2.D0/3.D0
c$$$        return
c$$$        end
c-----------------------------------------------------------------
      subroutine MHDSolve(ux,phi,dt,ttot)
c-----------------------------------------------------------------
c     Description: Fourth order Runge-Kutta time integration 
c     of the semi-discrete form of the MHD equations
c     Input:  dt -  time step
c     ux - conserved variables at time ttot
c     Output: ux - conserved variables at time ttot+dt
c     phi - auxillary variable used for projection
c-----------------------------------------------------------------
      use mesh
      use properties
#ifdef PARALLEL
      use mpistuff
#endif
c     
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision:: dt, ttot
c     
c     Local variables
c     ux0 = storage for ux at t=ttot
      double precision:: ux0(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     vx - cell centered primitive variables
      double precision:: vx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     vrp{x,y,z}: Face centered primitive variables for faces {x,y,z}
      double precision:: vrpx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: vrpy(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: vrpz(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     {f,g,h}inv: Hyperbolic fluxes in the {x,y,z} directions
      double precision:: finv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: ginv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: hinv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     
c     {f,g,h}vis: Parabolic or diffusive fluxes in the {x,y,z} directions
      double precision:: fvis(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: gvis(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: hvis(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     
c     dudx, dbdx: Grad of U and B
      double precision:: dudx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,3,3)
      double precision:: dbdx(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,3,3)
c     RHS for RK4
      double precision:: rhs1(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: rhs2(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: rhs3(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: rhs4(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
c     divb: Cell centered divergence computed from face centered B
      double precision:: divb(NXLSIZE,NYLSIZE,NZLSIZE)
c     
      integer:: indx,indy,indz
      integer:: i,j,k,l
c     c
c     Flag for non-conservative source prop to divB: 0=no 1=yes
      integer:: doDivSource
c     Flag for divB Filter: 0=no 1=yes
      integer:: doDivFilter

c     Exchange ghost values and set BCs
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
      doDivSource=0
      doDivFilter=0
c     
      ux0=ux
      ginv=0.D0
      hinv=0.D0
      finv=0.D0
      gvis=0.D0
      hvis=0.D0
      fvis=0.D0
      dudx=0.D0
      dbdx=0.D0
      rhs1=0.D0
      rhs2=0.D0
      rhs3=0.D0
      rhs4=0.D0
c
c     
c     First stage in RK4
      call ConservedToPrimitive(ux,vx)
#ifndef NOVISCFLUX
      dudx(:,:,:,1,1:3)=vx(:,:,:,2:4)
      call GradientVector(dudx)
      dbdx(:,:,:,1,1:3)=vx(:,:,:,5:7)
      call GradientVector(dbdx)
#endif
c     
      inlo=ixlo; inhi=ixhi
      call InviscidFluxWrapper(finv,vx,1,dt)
      inlo=ixlo; inhi=ixhi
      indx=1; indy=0;indz=0;
      call ViscousFluxWrapper(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
c     
c     
#ifndef ONE_D
      inlo=iylo; inhi=iyhi
      call InviscidFluxWrapper(ginv,vx,2,dt)
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
c     
#ifndef TWO_D
      inlo=izlo; inhi=izhi
      call InviscidFluxWrapper(hinv,vx,3,dt)
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
c     
c     
      call ComputeRHS(rhs1,ux,finv,ginv,hinv,fvis,gvis,hvis,doDivSource)
      ux=ux0+0.5D0*dt*rhs1
c     
c     Exchange ghost values and set BCs
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
c     End of first stage of RK4
c     Second stage of the RK4 scheme
      call ConservedToPrimitive(ux,vx)
c     
#ifndef NOVISCFLUX
      dudx(:,:,:,1,1:3)=vx(:,:,:,2:4)
      call GradientVector(dudx)
      dbdx(:,:,:,1,1:3)=vx(:,:,:,5:7)
      call GradientVector(dbdx)
#endif
      divb=0.D0
      inlo=ixlo; inhi=ixhi
      call InviscidFluxWrapper(finv,vx,1,dt)
      inlo=ixlo; inhi=ixhi
      indx=1; indy=0;indz=0;
      call ViscousFluxWrapper(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
c     
#ifndef ONE_D
      inlo=iylo; inhi=iyhi
      call InviscidFluxWrapper(ginv,vx,2,dt)
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
c     
#ifndef TWO_D
      inlo=izlo; inhi=izhi
      call InviscidFluxWrapper(hinv,vx,3,dt)
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
c     
      call ComputeRHS(rhs2,ux,finv,ginv,hinv,fvis,gvis,hvis,doDivSource)
      ux=ux0+0.5D0*dt*rhs2
c     
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
c     
c     End of second stage of RK4
c     Third stage of the RK4 scheme
      call ConservedToPrimitive(ux,vx)
c     
#ifndef NOVISCFLUX
      dudx(:,:,:,1,1:3)=vx(:,:,:,2:4)
      call GradientVector(dudx)
      dbdx(:,:,:,1,1:3)=vx(:,:,:,5:7)
      call GradientVector(dbdx)
#endif
      divb=0.D0
      inlo=ixlo; inhi=ixhi
      call InviscidFluxWrapper(finv,vx,1,dt)
      inlo=ixlo; inhi=ixhi
      indx=1; indy=0;indz=0;
      call ViscousFluxWrapper(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
c     
#ifndef ONE_D
      inlo=iylo; inhi=iyhi
      call InviscidFluxWrapper(ginv,vx,2,dt)
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
c     
#ifndef TWO_D
      inlo=izlo; inhi=izhi
      call InviscidFluxWrapper(hinv,vx,3,dt)
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
c     
      call ComputeRHS(rhs3,ux,finv,ginv,hinv,fvis,gvis,hvis,doDivSource)
      ux=ux0+dt*rhs3
c     
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
c     
c     End of third stage of RK4
c     Fourth stage of the RK4 scheme
      call ConservedToPrimitive(ux,vx)
c     
#ifndef NOVISCFLUX
      dudx(:,:,:,1,1:3)=vx(:,:,:,2:4)
      call GradientVector(dudx)
      dbdx(:,:,:,1,1:3)=vx(:,:,:,5:7)
      call GradientVector(dbdx)
#endif
      divb=0.D0
      inlo=ixlo; inhi=ixhi
      call InviscidFluxWrapper(finv,vx,1,dt)
      inlo=ixlo; inhi=ixhi
      indx=1; indy=0;indz=0;
      call ViscousFluxWrapper(fvis,vx,dudx,dbdx,1,dt,indx,indy,indz)
c     
#ifndef ONE_D
      inlo=iylo; inhi=iyhi
      call InviscidFluxWrapper(ginv,vx,2,dt)
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
      call ViscousFluxWrapper(gvis,vx,dudx,dbdx,2,dt,indx,indy,indz)
c     
#ifndef TWO_D
      inlo=izlo; inhi=izhi
      call InviscidFluxWrapper(hinv,vx,3,dt)
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
      call ViscousFluxWrapper(hvis,vx,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
c     
      call ComputeRHS(rhs4,ux,finv,ginv,hinv,fvis,gvis,hvis,doDivSource)
      ux=ux0+dt*(rhs1/6.D0+rhs2/3.D0+rhs3/3.D0+rhs4/6.D0)
c     
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
c     
c     
c     Diffuse divergence of B based on cell-centered B
      if(doDivFilter.eq.1) then
         call DivBDiffusionFilter(ux)
      endif
c     
c     Exchange ghost values and set BCs one last time
      call mesh_update_bdry_async(ux,nvar)
      call SetBoundaryValues(ux)
c     
      return
      end
c     
c----------------------------------------------------------------------
      subroutine ComputeRHS(rhs,ux,finv,ginv,hinv,fvis,gvis,hvis,
     &     doDivSource)
      use mesh
      use mesh_common
      integer:: doDivSource
      double precision:: rhs(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     ux(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     finv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     ginv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     hinv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     fvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     gvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision,intent(in):: 
     &     hvis(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: divergB(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      integer:: i,j,k,l
      double precision:: dxi, dyi, dzi
c     
C     $DOACROSS LOCAL(i,j,k)
c
      dxi=1.D0/dx
      dyi=1.D0/dy
      dzi=1.D0/dz
c
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1,NXlocal,1
c
c     Density
               rhs(i,j,k,1)=
     &              -dxi*(
     &              finv(i,j,k,1)
     &              -finv(i-1,j,k,1))
     &              -dyi*(ginv(i,j,k,1)  
     &              -ginv(i,j-1,k,1))
     &              -dzi*(hinv(i,j,k,1)  
     &              -hinv(i,j,k-1,1))
c     momentum
               rhs(i,j,k,2)=
     &              -dxi*(
     &              finv(i,j,k,2)  
     &              -finv(i-1,j,k,2))
     &              -dyi*(
     &              ginv(i,j,k,4)  
     &              -ginv(i,j-1,k,4))
     &              -dzi*(hinv(i,j,k,3)  
     &              -hinv(i,j,k-1,3))
     &              +dxi*(
     &              fvis(i,j,k,2)  
     &              -fvis(i-1,j,k,2))
     &              +dyi*(
     &              gvis(i,j,k,4)  
     &              -gvis(i,j-1,k,4))
     &              +dzi*(
     &              hvis(i,j,k,3)  
     &              -hvis(i,j,k-1,3))
c     
               rhs(i,j,k,3)=
     &              -dxi*(
     &              finv(i,j,k,3)  
     &              -finv(i-1,j,k,3))
     &              -dyi*(
     &              ginv(i,j,k,2)  
     &              -ginv(i,j-1,k,2))
     &              -dzi*(hinv(i,j,k,4)  
     &              -hinv(i,j,k-1,4))
     &              +dxi*(
     &              fvis(i,j,k,3)  
     &              -fvis(i-1,j,k,3))
     &              +dyi*(
     &              gvis(i,j,k,2)  
     &              -gvis(i,j-1,k,2))
     &              +dzi*(
     &              hvis(i,j,k,4)  
     &              -hvis(i,j,k-1,4))
c
               rhs(i,j,k,4)=
     &              -dxi*(
     &              finv(i,j,k,4)  
     &              -finv(i-1,j,k,4))
     &              -dyi*(
     &              ginv(i,j,k,3)  
     &              -ginv(i,j-1,k,3))
     &              -dzi*(hinv(i,j,k,2)  
     &              -hinv(i,j,k-1,2))
     &              +dxi*(
     &              fvis(i,j,k,4)  
     &              -fvis(i-1,j,k,4))
     &              +dyi*(
     &              gvis(i,j,k,3)  
     &              -gvis(i,j-1,k,3))
     &              +dzi*(
     &              hvis(i,j,k,2)  
     &              -hvis(i,j,k-1,2))
c     X Magnetic Field
               rhs(i,j,k,5)=
     &              -dxi*(
     &              finv(i,j,k,5)  
     &              -finv(i-1,j,k,5))
     &              -dyi*(
     &              ginv(i,j,k,7)  
     &              -ginv(i,j-1,k,7))
     &              -dzi*(hinv(i,j,k,6)  
     &              -hinv(i,j,k-1,6))
     &              +dxi*(
     &              fvis(i,j,k,5)  
     &              -fvis(i-1,j,k,5))
     &              +dyi*(
     &              gvis(i,j,k,7)  
     &              -gvis(i,j-1,k,7))
     &              +dzi*(
     &              hvis(i,j,k,6)  
     &              -hvis(i,j,k-1,6))
c     
c     Y Magnetic Field
               rhs(i,j,k,6)=
     &              -dxi*(
     &              finv(i,j,k,6)  
     &              -finv(i-1,j,k,6))
     &              -dyi*(
     &              ginv(i,j,k,5)  
     &              -ginv(i,j-1,k,5))
     &              -dzi*(hinv(i,j,k,7)  
     &              -hinv(i,j,k-1,7))
     &              +dxi*(
     &              fvis(i,j,k,6)  
     &              -fvis(i-1,j,k,6))
     &              +dyi*(
     &              gvis(i,j,k,5)  
     &              -gvis(i,j-1,k,5))
     &              +dzi*(
     &              hvis(i,j,k,7)  
     &              -hvis(i,j,k-1,7))
c     
c     Z Magnetic Field
               rhs(i,j,k,7)=
     &              -dxi*(
     &              finv(i,j,k,7)  
     &              -finv(i-1,j,k,7))
     &              -dyi*(
     &              ginv(i,j,k,6)  
     &              -ginv(i,j-1,k,6))
     &              -dzi*(hinv(i,j,k,5)  
     &              -hinv(i,j,k-1,5))
     &              +dxi*(
     &              fvis(i,j,k,7)  
     &              -fvis(i-1,j,k,7))
     &              +dyi*(
     &              gvis(i,j,k,6)  
     &              -gvis(i,j-1,k,6))
     &              +dzi*(
     &              hvis(i,j,k,5)  
     &              -hvis(i,j,k-1,5))
c     
c     Energy
               rhs(i,j,k,8)=
     &              -dxi*(
     &              finv(i,j,k,8)  
     &              -finv(i-1,j,k,8))
     &              -dyi*(
     &              ginv(i,j,k,8)  
     &              -ginv(i,j-1,k,8))
     &              -dzi*(hinv(i,j,k,8)  
     &              -hinv(i,j,k-1,8))
     &              +dxi*(
     &              fvis(i,j,k,8)  
     &              -fvis(i-1,j,k,8))
     &              +dyi*(
     &              gvis(i,j,k,8)  
     &              -gvis(i,j-1,k,8))
     &              +dzi*(
     &              hvis(i,j,k,8)  
     &              -hvis(i,j,k-1,8))
c
            enddo
         enddo
      enddo
c
c
c
      if(doDivSource.eq.1) then
#ifdef TWO_D
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               divergB(i,j,k)=dxi*(ux(i+1,j,k,5)-ux(i-1,j,k,5))+
     &              dyi*(ux(i,j+1,k,6)-ux(i,j-1,k,6))
            enddo
         enddo
      enddo
#else
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               divergB(i,j,k)=dxi*(ux(i+1,j,k,5)-ux(i-1,j,k,5))+
     &              dyi*(ux(i,j+1,k,6)-ux(i,j-1,k,6))+
     &              dzi*(ux(i,j,k+1,7)-ux(i,j,k-1,7))
            enddo
         enddo
      enddo
#endif
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal,1
               bi=ux(i,j,k,5); bj=ux(i,j,k,6); bk=ux(i,j,k,7)
               vi=ux(i,j,k,2)/ux(i,j,k,1)
               vj=ux(i,j,k,3)/ux(i,j,k,1)
               vk=ux(i,j,k,4)/ux(i,j,k,1)
               rhs(i,j,k,2)=rhs(i,j,k,2)-divergB(i,j,k)*bi
               rhs(i,j,k,3)=rhs(i,j,k,3)-divergB(i,j,k)*bj
               rhs(i,j,k,4)=rhs(i,j,k,4)-divergB(i,j,k)*bk
               rhs(i,j,k,5)=rhs(i,j,k,5)-divergB(i,j,k)*vi
               rhs(i,j,k,6)=rhs(i,j,k,6)-divergB(i,j,k)*vj
               rhs(i,j,k,7)=rhs(i,j,k,7)-divergB(i,j,k)*vk
               rhs(i,j,k,8)=rhs(i,j,k,8)
     &              -divergB(i,j,k)*(vi*bi+vj*bj+vk*bk)
            enddo
         enddo
      enddo
      endif
      return
      end
c     
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c$$$      subroutine ConservedToPrimitive(u,v)
c$$$c-----------------------------------------------------------------------
c$$$c     Description: Convert conservative to primitive variables at cell 
c$$$c                  centers
c$$$c     Input: u - conservative
c$$$c     Ouput: v - primitive
c$$$c-----------------------------------------------------------------------
c$$$      use mesh
c$$$      use mesh_common
c$$$      use properties
c$$$      double precision:: u(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c$$$      double precision:: v(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c$$$      integer:: i,j,k,l
c$$$      double precision::  rho_inv, vi,vj,vk,bi,bj,bk
c$$$
c$$$
c$$$#ifdef ONE_D
c$$$      do k=1,1
c$$$         do j=1,1
c$$$#else
c$$$#ifdef TWO_D
c$$$      do k=1,1
c$$$#else
c$$$      do k=1-NGHOST,NZlocal+NGHOST,1
c$$$#endif
c$$$         do j=1-NGHOST,NYlocal+NGHOST,1
c$$$#endif
c$$$            do i=1-NGHOST, NXlocal+NGHOST,1
c$$$               rho_inv=1.0D0/u(i,j,k,1)
c$$$               v(i,j,k,1)=u(i,j,k,1)
c$$$               v(i,j,k,2)=u(i,j,k,2)*rho_inv
c$$$               v(i,j,k,3)=u(i,j,k,3)*rho_inv
c$$$               v(i,j,k,4)=u(i,j,k,4)*rho_inv
c$$$               v(i,j,k,5:7)=u(i,j,k,5:7)
c$$$               vi=v(i,j,k,2)
c$$$               vj=v(i,j,k,3)
c$$$               vk=v(i,j,k,4)
c$$$               bi=v(i,j,k,5)
c$$$               bj=v(i,j,k,6)
c$$$               bk=v(i,j,k,7)
c$$$               v(i,j,k,8)=(gamma-1.D0)*(u(i,j,k,8)-0.5D0*v(i,j,k,1)
c$$$     &              *(vi*vi+vj*vj+vk*vk)-0.5*(bi*bi+bj*bj+bk*bk))
c$$$c               if(v(i,j,k,8).lt.0.D0) then
c$$$c                  write(6,*) 'WARNING: Negative Pressure',i,j,k,
c$$$c     &                 v(i,j,k,8),0.5D0*v(i,j,k,1)
c$$$c     &              *(vi*vi+vj*vj+vk*vk),0.5*(bi*bi+bj*bj+bk*bk),
c$$$c     &                 u(i,j,k,8)
c$$$c               endif
c$$$c     Arbirtrary pressure floor - if breached print warning 
c$$$c     and recommend killing the process.
c$$$               if(v(i,j,k,8).lt.0.00001D0) then
c$$$                  write(6,*) 'WARNING: Pressure floor breached',
c$$$     &                 i,j,k,
c$$$     &                 v(i,j,k,8),0.5D0*v(i,j,k,1)
c$$$     &              *(vi*vi+vj*vj+vk*vk),0.5*(bi*bi+bj*bj+bk*bk),
c$$$     &                 u(i,j,k,8)
c$$$                  v(i,j,k,8)=0.00001D0
c$$$                  u(i,j,k,8)=v(i,j,k,8)/(gamma-1.D0)
c$$$     &                 +0.5D0*v(i,j,k,1)
c$$$     &                 *(vi*vi+vj*vj+vk*vk)+0.5*(bi*bi+bj*bj+bk*bk)
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$c
c$$$      return
c$$$      end
c$$$c----------------------------------------------------------------------
c$$$      subroutine BFaceDivergenceF(divb,vface,idir)
c$$$c-----------------------------------------------------------------------
c$$$c     Description: Compute portion of divB using normal
c$$$c                  magnetic field components on faces in direction idir
c$$$c     Input: vface - face-centered primitive vars in idir-direction
c$$$c     Ouput: divb 
c$$$c-----------------------------------------------------------------------
c$$$      use mesh
c$$$      use mesh_common
c$$$      use properties
c$$$      integer:: idir
c$$$      double precision:: divb(NXLSIZE,NYLSIZE,NZLSIZE)
c$$$      double precision:: vface(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c$$$      integer:: i,j,k
c$$$      integer:: ioff,joff,koff
c$$$      ioff=0
c$$$      joff=0
c$$$      koff=0
c$$$      if(idir.eq.1) ioff=1
c$$$      if(idir.eq.2) joff=1
c$$$      if(idir.eq.3) koff=1
c$$$c
c$$$c     BCs at faces (need to revisit - must have all BCs specified separately
c$$$c     from the solver code
c$$$c$$$      if(idir.eq.1) then
c$$$c$$$         if(iprocx.eq.1) then
c$$$c$$$            vface(0,:,:,2)=0.D0
c$$$c$$$            vface(0,:,:,5)=0.D0
c$$$c$$$         endif
c$$$c$$$         if(iprocx.eq.XPROCS) then
c$$$c$$$            vface(nxlocal,:,:,2)=0.D0
c$$$c$$$            vface(nxlocal,:,:,5)=0.D0
c$$$c$$$         endif
c$$$c$$$      endif
c$$$c
c$$$      if(idir.eq.2) then
c$$$         if(iprocy.eq.1) then
c$$$            vface(:,0,:,5)=0.D0
c$$$            vface(:,0,:,2)=0.D0
c$$$         endif
c$$$         if(iprocy.eq.YPROCS) then
c$$$            vface(:,nylocal,:,2)=0.D0
c$$$            vface(:,nylocal,:,5)=0.D0
c$$$         endif
c$$$      endif
c$$$c
c$$$      do k=1,NZlocal,1
c$$$         do j=1,NYlocal,1
c$$$            do i=1, NXlocal,1
c$$$               divb(i,j,k)=divb(i,j,k)+
c$$$     &              (vface(i,j,k,5)
c$$$     &              -vface(i-ioff,j-joff,k-koff,5))/dx
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      return
c$$$      end
c$$$c
c$$$c----------------------------------------------------------------------
c$$$      subroutine BFaceDivergenceF2(divb,ux)
c$$$c-----------------------------------------------------------------------
c$$$c     Description: Compute portion of divB using normal
c$$$c                  magnetic field components on faces which are
c$$$c                  obtained by averaging the cell centered values
c$$$c     Input: ux - cell-centered conservative vars
c$$$c     Ouput: divb 
c$$$c-----------------------------------------------------------------------
c$$$      use mesh
c$$$      use mesh_common
c$$$      use properties
c$$$      integer:: idir
c$$$      double precision:: divb(NXLSIZE,NYLSIZE,NZLSIZE)
c$$$      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c$$$      integer:: i,j,k
c$$$c
c$$$      divb=0.D0
c$$$      do k=1,NZlocal,1
c$$$         do j=1,NYlocal,1
c$$$            do i=1, NXlocal,1
c$$$               divb(i,j,k)=
c$$$     &              (ux(i+1,j,k,5)
c$$$     &              -ux(i-1,j,k,5)+
c$$$     &              ux(i,j+1,k,6)
c$$$     &              -ux(i,j-1,k,6))/2.D0/dx
c$$$c               write(6,*) 'DIVB=',i,j,divb(i,j,k)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      return
c$$$      end
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$      subroutine DivBDiffusionFilter(ux)
c$$$c-----------------------------------------------------------------------
c$$$c     Description: diffuse away divergence; this has a large stabilizing
c$$$c                  effect on the code. Diffusing away divB implies 
c$$$c                  a certain filter on B noted below.
c$$$c     Input: ux - conservative variables before filter
c$$$c     Ouput: ux - conservative variables after filter
c$$$c-----------------------------------------------------------------------
c$$$      use mesh
c$$$      use mesh_common
c$$$      use properties
c$$$      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
c$$$c
c$$$      double precision:: bnew(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,3)
c$$$      integer:: i,j,k,l
c$$$      double precision::  rho_inv, vi,vj,vk,bi,bj,bk
c$$$      double precision:: dxi,dyi,dzi
c$$$      double precision:: dxi2,dyi2,dzi2,dxidyi
c$$$      double precision:: lamda
c$$$
c$$$      lamda=dx*dx/80.D0
c$$$c
c$$$      dxi=0.5D0/dx
c$$$      dyi=0.5D0/dy
c$$$      dzi=0.5D0/dz
c$$$      dxi2=1.D0/dx/dx
c$$$      dyi2=1.D0/dy/dy
c$$$      dzi2=1.D0/dz/dz
c$$$      dxidyi=0.25D0/dx/dy
c$$$c
c$$$      bnew(:,:,:,1)=ux(:,:,:,5)
c$$$      bnew(:,:,:,2)=ux(:,:,:,6)
c$$$      bnew(:,:,:,3)=ux(:,:,:,7)
c$$$      do k=1,NZlocal,1
c$$$         do j=1,NYlocal,1
c$$$            do i=1, NXlocal,1
c$$$               bnew(i,j,k,1)=ux(i,j,k,5)+lamda*(
c$$$     &              0.25D0*
c$$$     &              dxi2*(ux(i+2,j,k,5)+ux(i-2,j,k,5)-2.D0*ux(i,j,k,5))
c$$$     &              +dxidyi*(ux(i+1,j+1,k,6)+ux(i-1,j-1,k,6)-
c$$$     &              ux(i-1,j+1,k,6)-ux(i+1,j-1,k,6)))
c$$$               bnew(i,j,k,2)=ux(i,j,k,6)+lamda*(
c$$$     &              0.25D0*
c$$$     &              dyi2*(ux(i,j+2,k,6)+ux(i,j-2,k,6)-2.D0*ux(i,j,k,6))
c$$$     &              +dxidyi*(ux(i+1,j+1,k,5)+ux(i-1,j-1,k,5)-
c$$$     &              ux(i-1,j+1,k,5)-ux(i+1,j-1,k,5)))
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$c$$$      do k=1,NZlocal,1
c$$$c$$$         do j=1,NYlocal,1
c$$$c$$$            do i=1, NXlocal,1
c$$$c$$$               bnew(i,j,k,1)=ux(i,j,k,5)+lamda*(
c$$$c$$$     &              dxi2*(ux(i+1,j,k,5)+ux(i-1,j,k,5)-2.D0*ux(i,j,k,5))
c$$$c$$$     &              +dxidyi*(ux(i+1,j+1,k,6)+ux(i-1,j-1,k,6)-
c$$$c$$$     &              ux(i-1,j+1,k,6)-ux(i+1,j-1,k,6)))
c$$$c$$$               bnew(i,j,k,2)=ux(i,j,k,6)+lamda*(
c$$$c$$$     &              dyi2*(ux(i,j+1,k,6)+ux(i,j-1,k,6)-2.D0*ux(i,j,k,6))
c$$$c$$$     &              +dxidyi*(ux(i+1,j+1,k,5)+ux(i-1,j-1,k,5)-
c$$$c$$$     &              ux(i-1,j+1,k,5)-ux(i+1,j-1,k,5)))
c$$$c$$$            enddo
c$$$c$$$         enddo
c$$$c$$$      enddo
c$$$c
c$$$c     Is it advisable to readjust total energy? Don't think so. 
c$$$c     A place holder for this is below -just in case we need it later.
c$$$c$$$      call ConservedEnergyCorrection(ux(ixlo,iylo,izlo,8),
c$$$c$$$     &       bnew,ux(ixlo,iylo,izlo,5))
c$$$      ux(1:nxlocal,1:nylocal,1:nzlocal,5)
c$$$     &     =bnew(1:nxlocal,1:nylocal,1:nzlocal,1)
c$$$      ux(1:nxlocal,1:nylocal,1:nzlocal,6)
c$$$     &     =bnew(1:nxlocal,1:nylocal,1:nzlocal,2)
c$$$      return
c$$$      end
c
c-----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine fluxInterpolateCD2(finv,ginv,hinv)
      use mesh
      use mesh_common
      double precision:: 
     &     finv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: 
     &     ginv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: 
     &     hinv(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI,nvar)
      integer:: i,j,k,l 
      double precision:: 
     &     fmod(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI)
      double precision:: 
     &     gmod(IXLO:IXHI, IYLO:IYHI, IZLO:IZHI)
c
      fmod=finv(:,:,:,6)
      gmod=ginv(:,:,:,7)
c
      do k=1,NZlocal,1
         do j=1,NYlocal,1
            do i=1, NXlocal-1,1
               fmod(i,j,k)=0.125D0*(
     &              2.D0*finv(i,j,k,6)
     &              +finv(i+1,j,k,6)
     &              +finv(i-1,j,k,6)
     &              -ginv(i+1,j-1,k,7)
     &              -ginv(i+1,j,k,7)
     &              -ginv(i,j-1,k,7)
     &              -ginv(i,j,k,7))

            enddo
         enddo
      enddo
c
      do k=1,NZlocal,1
         do j=1,NYlocal-1,1
            do i=1, NXlocal,1
               gmod(i,j,k)=0.125D0*(
     &              2.D0*ginv(i,j,k,7)
     &              +ginv(i,j+1,k,7)
     &              +ginv(i,j-1,k,7)
     &              -finv(i-1,j,k,6)
     &              -finv(i-1,j+1,k,6)
     &              -finv(i,j,k,6)
     &              -finv(i,j+1,k,6))
            enddo
         enddo
      enddo
c
      finv(:,:,:,6)=fmod(:,:,:)
      ginv(:,:,:,7)=gmod(:,:,:)
      return
      end

