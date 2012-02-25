C -*- Mode: Fortran; -*-
c----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c----------------------------------------------------------------------
c Subroutine to set up initial conditions for the Tilt-Mode problem
c----------------------------------------------------------------------

      subroutine InitialConditions(ux,phi)
c----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use iounits
      implicit none

c=======Declarations=========
c
c     input variables
      double precision :: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
c
c     local variables
      integer          :: i, j, k, bier
      double precision :: b0, beta, aradius, pi, v0, press
      double precision :: radius, sintheta, j0, j1, kbess, d1, d2
      double precision :: psi(ixlo:ixhi,iylo:iyhi,izlo:izhi)

      double precision, dimension(1000) :: rvec, Jvec, Jpvec, Jppvec, 
     &     Bessvec
c     
      namelist /molwts/ wmol
      namelist /gammas/ gamma
      namelist /diffusiveproperties/ mu, eta, etaFactor, kappa
      namelist /tminit/ b0, beta, aradius
c

c=======Internals=========
     
c     input general problem parameters
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

c     setup additional plasma parameter data (for viscous fluxes)
      rgas = 1.d0/wmol
     
c     input tilt mode parameters
      open(17,file='tiltmode.inp',form='formatted')
      read(17,tminit)
      close(17)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,3(A,es12.5))') 
     &        '\nTile Mode Properties:',
     &        '\n       b0 =',b0,
     &        '\n     beta =',beta,
     &        '\n  aradius =',aradius
      end if

c     compute initial pressure based on plasma beta
c     Q: DO WE SET PRESSURE OFF OF SCALING PARAMETER B0, OR THE MAXIMUM VALUE OF B???
c     [B is somewhat large in central region]
c$$$      press = beta*b0*b0/2.d0
      press = beta*(1.d1*b0)*(1.d1*b0)/2.d0

c     compute velocity perturbation magnitude (Alfven velocity, rho=1)
      v0 = 1.d-3*b0

c     set bessel function input scaling, pi
      kbess = 3.83170597020751d0/aradius
      call besj( kbess*aradius, 0.d0, 1, j0, bier)
      if (bier /= 0) then
         write(0,*) 'Tilt Mode init: error in besj =',bier
      endif
      write(0,*) ' besj0(k*aradius) =',j0
      pi = 4.D0*datan(1.D0)

c     set differencing parameters
#ifdef CD4
c        fourth order standard central difference
      d1 = 2.d0/3.d0
      d2 = -1.d0/12.d0
#elif  TCD
c        second order tuned central difference
      d1 = 0.894d0
      d2 = -0.197d0
#else
c        second order standard central difference
      d1 = 0.5d0
      d2 = 0.0d0
#endif
     
c     initialize constant fields
      ux(:,:,:,1) = 1.d0
      ux(:,:,:,2) = 0.d0
      ux(:,:,:,3) = 0.d0
      ux(:,:,:,4) = 0.d0
      ux(:,:,:,5) = 0.d0
      ux(:,:,:,6) = 0.d0
      ux(:,:,:,7) = b0
      ux(:,:,:,8) = press/(gamma-1.D0)
     
c     iterate over the domain to set up magnetic flux function
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1

c              set polar coordinates
c                [uses sin(theta) = y/r, so don't need theta directly]
               radius = dsqrt(xc(i)**2 + yc(j)**2)
               sintheta = yc(j)/radius

c              compute Bessel function at this location
               call besj( kbess*radius, 1.d0, 1, j1, bier)
               if (bier /= 0) then
                  write(0,*) 'Tilt Mode init: error in besj =',bier
               endif

c              set psi
               if (radius <= aradius) then
c$$$                  psi(i,j,k) = -2.d0*b0/j0*j1*sintheta
                  psi(i,j,k) = 2.d0*b0/kbess/j0*j1*sintheta
               else
                  psi(i,j,k) = b0*(radius-aradius**2/radius)*sintheta
               endif

            enddo
         enddo
      enddo
     
c     iterate over the domain to set spatially-varying fields
      do k=1,Nzlocal,1
         do j=0,Nylocal+1,1
            do i=0,Nxlocal+1,1

c              set radial coordinate
               radius = dsqrt(xc(i)**2 + yc(j)**2)

c              set momentum perturbation (y-direction only)
               ux(i,j,k,3) =-ux(i,j,k,1)*v0*dexp(-radius**2)*dsin(xc(i))

c              compute initial perpendicular magnetic field
               if (radius <= aradius) then
                  ux(i,j,k,7) = dsqrt(b0**2 + kbess**2 * psi(i,j,k)**2)
               else
                  ux(i,j,k,7) = b0
               endif

c              set planar magnetic field
               ux(i,j,k,5) = -(d1*(psi(i,j+1,k)-psi(i,j-1,k))
     &                        +d2*(psi(i,j+2,k)-psi(i,j-2,k)))/dy
               ux(i,j,k,6) = (d1*(psi(i+1,j,k)-psi(i-1,j,k))
     &                       +d2*(psi(i+2,j,k)-psi(i-2,j,k)))/dx

c              set total energy
               ux(i,j,k,8) = press/(gamma-1.D0) 
     &              + 0.5d0*(ux(i,j,k,5)*ux(i,j,k,5)
     &                     + ux(i,j,k,6)*ux(i,j,k,6)
     &                     + ux(i,j,k,7)*ux(i,j,k,7))
     &              + 0.5d0*(ux(i,j,k,2)**2
     &                     + ux(i,j,k,3)**2
     &                     + ux(i,j,k,4)**2)/ux(i,j,k,1)

            enddo
         enddo
      enddo
c
      return
c     
      end
c
c======================================================================
