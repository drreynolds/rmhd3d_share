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
c Subroutine to set up the initial conditions for 
c Kelvin-Helmholtz instability test problem
c----------------------------------------------------------------------

      subroutine InitialConditions(ux,phi)
c----------------------------------------------------------------------
      use mesh
      use mesh_common
      use properties
      use iounits
      implicit none

c=======Declarations=========

c     input variables
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)

c     local variables
      integer          :: i,j,k
      double precision :: xloc, yloc, zloc, pi, den
      double precision :: lambda, epsilon, p0, perturb
      double precision :: press, Bx, By, etot
      double precision :: pE, pW, pN, pS, xE, xW, yN, yS

      namelist /molwts/ wmol
      namelist /gammas/ gamma
      namelist /diffusiveproperties/ mu, eta, etaFactor, kappa
      namelist /island/ lambda, epsilon, p0, perturb

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

c     input island coalescence parameters
      open(17,file='island.inp',form='formatted')
      read(17,island)
      close(17)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,4(A,es12.5))') 
     &        '\nIsland Coalescence Properties:',
     &        '\n         lambda =',lambda,
     &        '\n        epsilon =',epsilon,
     &        '\n             p0 =',p0,
     &        '\n        perturb =',perturb
      end if


c     initialize state to constants
      ux  = 0.d0
      phi = 0.d0

c     set pi
      pi = 4.D0*datan(1.D0)

c     set density
      ux(:,:,:,1) = 1.d0

c     iterate over the domain to set up x-velocity
      do k=izlo,izhi,1
         zloc = ZLloc + (k-0.5d0)*dz

         do j=iylo,iyhi,1
            yloc = YLloc + (j-0.5D0)*dy
            yS = yloc - dy
            yN = yloc + dy

            do i=ixlo,ixhi,1
               xloc = XLloc + (i-0.5D0)*dx
               xW = xloc - dx
               xE = xloc + dx

c              set pressure 
               press = p0 + (1-epsilon**2)/2.d0/(cosh(xloc/lambda) 
     &                       + epsilon*cos(yloc/lambda))**2


ccccc Analytical approach for initial conditions (more accurate) ccccc
c$$$c              compute denominator used in magnetic fields
c$$$               den = cosh(xloc/lambda) + epsilon*cos(yloc/lambda)
c$$$
c$$$c              set x-magnetic field
c$$$               Bx = -epsilon*sin(yloc/lambda)/den + perturb*pi
c$$$     &              * sin(0.5d0*pi*xloc)*sin(pi*yloc)*cos(pi*zloc)
c$$$
c$$$c              set y-magnetic field
c$$$               By = -sinh(xloc/lambda)/den + perturb*0.5d0*pi
c$$$     &              * cos(0.5d0*pi*xloc)*cos(pi*yloc)*cos(pi*zloc)


ccccc Numerical approach for initial conditions (better numerical 
ccccc preservation of divergence constraint on B field) ccccc

c              compute psi in neighboring cells
               pE = -lambda*log(cosh(xE/lambda) 
     &              + epsilon*cos(yloc/lambda)) 
     &              + perturb*sin(0.5d0*pi*xE)*cos(pi*yloc)*cos(pi*zloc)
               pW = -lambda*log(cosh(xW/lambda) 
     &              + epsilon*cos(yloc/lambda)) 
     &              + perturb*sin(0.5d0*pi*xW)*cos(pi*yloc)*cos(pi*zloc)
               pN = -lambda*log(cosh(xloc/lambda) 
     &              + epsilon*cos(yN/lambda)) 
     &              + perturb*sin(0.5d0*pi*xloc)*cos(pi*yN)*cos(pi*zloc)
               pS = -lambda*log(cosh(xloc/lambda) 
     &              + epsilon*cos(yS/lambda)) 
     &              + perturb*sin(0.5d0*pi*xloc)*cos(pi*yS)*cos(pi*zloc)

c              use centered differences for magnetic field entries
               Bx = -(pN - pS)/(2.d0*dy)
               By =  (pE - pW)/(2.d0*dx)


c              set total energy
               etot = press/(gamma-1.D0) + 0.5d0*(Bx*Bx+By*By)

c              store values
               ux(i,j,k,5) = Bx
               ux(i,j,k,6) = By
               ux(i,j,k,8) = etot

            enddo
         enddo
      enddo

      return
      end
c======================================================================
