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
c=================================================================

      subroutine InitialConditions(ux,phi)
c-----------------------------------------------------------------
c     Description: subroutine to set up the initial conditions
c        for the reconnection simulation
c
c     Arguments:
c          ux - (in) current state of all variables over local mesh
c         phi - (in) projection information for div B
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits
      
c======= Declarations =========
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: psi(IXLO-1:IXHI+1,IYLO-1:IYHI+1,IZLO:IZHI)
      double precision :: rho, press, Bi, Bj, Bk, u, v, w
      double precision :: dxi, dyi
      double precision :: x, y, z
      double precision :: pi, kx, ky, lamda
      integer :: i, j, k, l
      character*24 fdate


c======= Internals ============

c     set up namelists for plasma parameter data
      namelist /molwts/ wmol
      namelist /gammas/ gamma
      namelist /diffusiveproperties/ mu, eta, etaFactor, kappa

c     input plasma parameter data, root node outputs info
      open(16,file='prop.inp',form='formatted')
      read(16,molwts)
      read(16,gammas) 
      read(16,diffusiveproperties) 
      close(16)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,6(A,es12.5))') '\nPlasma properties:',
     &        '\n       wmol =',wmol,
     &        '\n      gamma =',gamma,
     &        '\n         mu =',mu,
     &        '\n        eta =',eta,
     &        '\n  etaFactor =',etaFactor,
     &        '\n      kappa =',kappa
      end if

c     write(6,*) 'In InitFluid after reading parameters'

c     setup additional plasma parameter data
ccccc CHANGED PER RAVI'S REQUEST, 6/17/2004
c      rgas = 8317.D0/wmol
      rgas = 1.0d0/wmol
      pi = 4.D0*datan(1.D0)
ccccc CHANGED PER RAVI'S REQUEST, 6/17/2004
c      dxi = 0.5D0/dx
c      dyi = 0.5D0/dy
      dxi = 1.0d0/dx
      dyi = 1.0d0/dy


c     initialize projection information for div B
      phi = 0.D0

c     initialize further plasma parameters
      rho = 1.D0
ccccc CHANGED PER RAVI'S REQUEST, 6/8/2004
c      press = 0.2D0
      press = 0.5D0
      u = 0.D0
      v = 0.D0
      w = 0.D0
     
c     For reconnection at x=0
c      kx = 3.D0*pi/2.D0
c     ky = 2.D0*pi
c      ky = pi
ccccc NEW FORM FOR MORE BENIGN INITIAL CONDITION
      kx = pi/xr
      ky = pi/(2.0d0*yr)
      psi0 = 0.1d0
cccccccccccccccccccccccccccccccccccccccccccccccc
     
c     For reconnection at y=0
c     ky = 3.D0*pi/2.D0
c     kx = pi

c     set up internal calculation variable
      lamda=dsqrt(kx*kx+ky*ky)

c     set up initial density psi
      do k=izlo,izhi,1
         do j=iylo-1,iyhi+1,1
            y=(j)*dy-0.5D0*dy
c            y=yc(j)
            do i=ixlo-1,ixhi+1,1
               x=-1.D0+(i)*dx-0.5D0*dx
c               x=xc(i)
c              For reconnection at x=0
c               psi(i,j,k) = dcos(kx*x)*dsin(ky*y)
c              For reconnection at y=0
c               psi(i,j,k) = -cos(ky*y)*sin(kx*x)
cccccccccccccc NEW FORM FOR MORE BENIGN INITIAL CONDITION
               psi(i,j,k) = 0.5d0*log(cosh(2.0d0*y)) + 
     &              psi0*dcos(kx*x)*dcos(ky*y)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            enddo
         enddo
      enddo

c     Initialize everything to zero
      ux=0.D0
c     set up initial density, velocity, magnetic field, energy 
      do k = izlo,izhi,1
ccccc CHANGED PER RAVI'S REQUEST, 6/17/2004
c         do j = iylo,iyhi,1
         do j = 1,Nylocal,1
            y = yc(j)
ccccc CHANGED PER RAVI'S REQUEST, 6/17/2004
c            do i = ixlo,ixhi,1
            do i = 1,nxlocal,1
               x = xc(i)
               
c              For reconnection at y=0
c               Bi = -ky*dsin(ky*y)*dsin(kx*x)
c               Bj = -kx*dcos(ky*y)*dcos(kx*x)

c              For reconnection at x=0
ccccc INITIAL CONDITION #1 
c               Bi = -ky*dcos(kx*x)*dcos(ky*y)
c               Bj = -kx*dsin(kx*x)*dsin(ky*y)
ccccc INITIAL CONDITION #2 (zero initial divB)
c               Bi = -(psi(i,j+1,k)-psi(i,j-1,k))*dyi
c               Bj = (psi(i+1,j,k)-psi(i-1,j,k))*dxi
ccccc INITIAL CONDITION #3 (Ravi's fix, 3/12/04)
#ifndef CD4
               Bi = tanh(y/0.5d0)-psi0*ky*dcos(kx*x)*dsin(ky*y)
               Bj = psi0*kx*dsin(kx*x)*dcos(ky*y)
#else
cccccccccccccc ADDED FOR 4TH-ORDER DIV(B)=0, 6/17/2004
               Bi = (2.0d0/3.0d0*(psi(i,j+1,k)-psi(i,j-1,k))*dyi 
     &              - 1.0d0/12.0d0*(psi(i,j+2,k)-psi(i,j-2,k))*dyi) 
               Bj = -(2.0d0/3.0d0*(psi(i+1,j,k)-psi(i-1,j,k))*dxi
     &              - 1.0d0/12.0d0*(psi(i+2,j,k)-psi(i-2,j,k))*dxi)
#endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccc NEW FORM FOR MORE BENIGN INITIAL CONDITION
cccccccccccccc CHANGED RHO PER RAVI'S REQUEST, 6/8/2004
c               rho = 1.0d0 + 0.2d0*(1.0d0 - tanh(y/0.5d0)**2)
               rho = 0.2d0+(1.0d0-tanh(y/0.5d0)**2)
cccccccccccccc ADDED THE FOLLOWING PER RAVI'S REQUEST, 6/15/2004
               press = 0.5d0*rho
c
               Bk = 0.0d0
c              For reconnection at x=0
c               Bk = -lamda*dcos(kx*x)*dsin(ky*y)
c              For reconnection at y=0
c               Bk = -lamda*dcos(ky*y)*dsin(kx*x)
c               if(cos(kx*x)*sin(ky*y).gt.0.D0) then
c                  Bi = 0.1*Bi; Bj=0.1*Bj; Bk=0.1*Bk;
c                  Bi=0.D0; Bj=0.D0; Bk=0.D0
c               endif
c	        if(Bk.lt.0.D0) then
c                  rho = 0.1D0
c                  press = 0.002D0
c               else
c                  rho = 1.D0
c                  press = 0.2D0
c               endif

c              finally insert initial state
               ux(i,j,k,1) = rho
               ux(i,j,k,2) = rho*u
               ux(i,j,k,3) = rho*v
               ux(i,j,k,4) = rho*w
               ux(i,j,k,5) = Bi
               ux(i,j,k,6) = Bj
               ux(i,j,k,7) = Bk
               ux(i,j,k,8) = press/(gamma-1.D0)+
     &              0.5D0*rho*(u*u+v*v+w*w)+0.5D0*(Bi*Bi+Bj*Bj+Bk*Bk)
            enddo
         enddo
      enddo

ccccc ADDED THE FOLLOWING PER RAVI'S REQUEST, 6/17/2004
      call SetBoundaryValues(ux)

c     check initial Divergence
      Bj = 0.0d0
      Bk = 0.0d0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               Bi=2.D0/3.D0*((ux(i+1,j,k,5)-ux(i-1,j,k,5))/dx+
     &                (ux(i,j+1,k,6)-ux(i,j-1,k,6))/dy)
     &            -1.D0/12.D0*((ux(i+2,j,k,5)-ux(i-2,j,k,5))/dx+
     &                (ux(i,j+2,k,6)-ux(i,j-2,k,6))/dy)
            enddo
         enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc CHANGED PER RAVI'S REQUEST, 6/17/2004
c      phi = 0.0D0

      return


c     [Dan: I don't know what these are here for, but they aren't used]
      do k=izlo,izhi,1
         do j=iylo,iyhi/2,1
            do i=ixhi/2,ixhi,1
               ux(i,j,k,5)=ux(ixhi-i+ixlo,j,k,5)
               ux(i,j,k,6)=-ux(ixhi-i+ixlo,j,k,6)
               ux(i,j,k,7)=ux(ixhi-i+ixlo,j,k,7)
               ux(i,j,k,8)=ux(ixhi-i+ixlo,j,k,8)
               psi(i,j,k)=psi(ixhi-i+ixlo,j,k)
            enddo
         enddo
      enddo
      do k=izlo,izhi,1
         do j=iyhi/2,iyhi,1
            do i=ixlo,ixhi,1
               ux(i,j,k,5)=-ux(i,iyhi-j+iylo,k,5)
               ux(i,j,k,6)=ux(i,iyhi-j+iylo,k,6)
               ux(i,j,k,7)=ux(i,iyhi-j+iylo,k,7)
               ux(i,j,k,8)=ux(i,iyhi-j+iylo,k,8)
               psi(i,j,k)=psi(i,iyhi-j+iylo,k)
            enddo
         enddo
      enddo
     
      phi = -psi
      return
      end
c=================================================================
