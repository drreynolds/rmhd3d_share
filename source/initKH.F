C -*- Mode: Fortran; -*-
c----------------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c----------------------------------------------------------------------
c     $Log: initKH.F,v $
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
c
c     input variables
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
c
c     local variables
      integer          :: i,j,k,xoscill,yoscill,zoscill
      double precision :: xloc, yloc, zloc, pi
      double precision :: v0, lambda, xperturb, yperturb, zperturb
      double precision :: press, Bx, Bz
c     
      namelist /molwts/ wmol
      namelist /gammas/ gamma
      namelist /diffusiveproperties/ mu, eta, etaFactor, kappa
      namelist /khinit/ v0, lambda, xperturb, yperturb, zperturb, 
     &     xoscill, yoscill, zoscill, press, Bz, Bx
c

c=======Internals=========
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
c     input Kelvin Helmholtz parameters
      open(17,file='kh_init.inp',form='formatted')
      read(17,khinit)
      close(17)
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(6,'(A,5(A,es12.5),3(A,i3),3(A,es12.5))') 
     &        '\nKelvin-Helmholtz Properties:',
     &        '\n             v0 =',v0,
     &        '\n         lambda =',lambda,
     &        '\n       xperturb =',xperturb,
     &        '\n       yperturb =',yperturb,
     &        '\n       zperturb =',zperturb,
     &        '\n        xoscill =',xoscill,
     &        '\n        yoscill =',yoscill,
     &        '\n        zoscill =',zoscill,
     &        '\n          press =',press,
     &        '\n             Bx =',Bx,
     &        '\n             Bz =',Bz
      end if
c     
c     
c     initialize state to constants
      ux(:,:,:,1) = 1.d0
      ux(:,:,:,2) = 0.d0
      ux(:,:,:,3) = 0.d0
      ux(:,:,:,4) = 0.d0
      ux(:,:,:,5) = Bx
      ux(:,:,:,6) = 0.d0
      ux(:,:,:,7) = Bz
      ux(:,:,:,8) = press/(gamma-1.D0) + 0.5d0*(Bx*Bx+Bz*Bz)
c     
c     set pi
      pi = 4.D0*datan(1.D0)
c     
c     iterate over the domain to set up x-velocity
      do k=izlo,izhi,1
         zloc = ZLloc + (k-0.5d0)*dz

         do j=iylo,iyhi,1
            yloc = YLloc + (j-0.5D0)*dy

            do i=ixlo,ixhi,1
               xloc = XLloc + (i-0.5D0)*dx

c              set x-momentum
               ux(i,j,k,2) = ux(i,j,k,1)*(v0*dtanh(yloc/lambda) 
     &              + xperturb*dcos(xloc/(xR-xL)*2.d0*pi*xoscill)
     &              + yperturb*dsin(yloc/(yR-yL)*pi*(2*yoscill-1))
     &              + zperturb*dcos(zloc/(zR-zL)*2.d0*pi*zoscill))

c              update total energy
               ux(i,j,k,8) = ux(i,j,k,8) 
     &                     + 0.5d0*ux(i,j,k,2)**2/ux(i,j,k,1)
            enddo
         enddo
      enddo
      return
c     
      end
c
c======================================================================
