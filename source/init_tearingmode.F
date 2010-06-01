C -*- Mode: Fortran; -*-
c----------------------------------------------------------------------
c     Daniel Reynolds
c     Copyright 2006
c     University of California, San Diego
c     All Rights Reserved
c----------------------------------------------------------------------
c     Subroutine to set up the initial conditions for resistive 
c     tearing-mode instability test problem.  Initial conditions as 
c     outlined in L. Chacon's 3d XMHD poster at Sherwood 2004 meeting.
c     
c     It looks like I am missing the form of perturbation to the 
c     initial state (also missing from Chacon's paper/poster), if 
c     there is one for this problem (there is in other tearing mode 
c     formulations).
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
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
c
c     local variables
      integer          :: i,j,k
      double precision :: xloc, yloc, zloc
c
c     tearing-mode initialization parameters
      double precision, parameter :: lambda=0.2d0
      double precision, parameter :: B0=1.d0
c     
      namelist/molwts/wmol
      namelist/gammas/gamma
      namelist/diffusiveproperties/mu,eta,etaFactor,kappa
c

c=======Internals=========
c     
c     initialize state to constants
      ux(:,:,:,1) = 1.d0
      ux(:,:,:,2) = 0.d0
      ux(:,:,:,3) = 0.d0
      ux(:,:,:,4) = 0.d0
      ux(:,:,:,5) = 0.d0
      ux(:,:,:,6) = 0.d0
      ux(:,:,:,7) = 0.d0
      ux(:,:,:,8) = 1.d0
c     
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
c
c     setup additional plasma parameter data (for viscous fluxes)
      rgas = 1.d0/wmol
c     
c     
c     iterate over the domain to set up non-homogeneous fields
      do k=izlo,izhi,1
         zloc = ZLloc + (k-0.5d0)*dz

         do j=iylo,iyhi,1
            yloc = YLloc + (j-0.5D0)*dy;

            do i=ixlo,ixhi,1
               xloc = XLloc + (i-0.5D0)*dx;

c              initial x and z-directional magnetic fields
               ux(i,j,k,5) = B0*tanh(yloc/lambda)
               ux(i,j,k,7) = dsqrt(B0*B0-ux(i,j,k,5)**2)

c              perturbation to initial field????

c              resulting total energy
               ux(i,j,k,8) = ux(i,j,k,8) + 
     &              0.5d0*(ux(i,j,k,5)**2 + ux(i,j,k,7)**2)
               
            enddo
         enddo
      enddo
c     
      return
      end
c
c======================================================================
