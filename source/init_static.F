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
c-----------------------------------------------------------------
c     Subroutine to set up the initial conditions - equilibrium
c     solving grad(p)=JxB in Cartesian geometry
c---------------------------------------------------------------------------

      subroutine InitialConditions(ux,phi)
      use mesh
      use mesh_common
      use properties
      use iounits
#ifdef PARALLEL
      use mpistuff
#endif
      double precision:: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision:: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision:: rho,press, Bi, Bj, Bk, u,v,w
c     
      double precision:: x,y,z,xface
      double precision::  a
      double precision::  asq
      double precision:: alpha, pi, lamda,aspect
      double precision:: fun(ixlo:ixhi)
      double precision:: dfdxi(ixlo:ixhi)
      double precision:: d2fdxi2(ixlo:ixhi)
      double precision:: fmin
      double precision:: psi, dpsidx,dpsidy, g0, pisq
      double precision::  d2psidx2,d2psidy2
c
      double precision:: xi,zeta
      double precision:: coszeta, sinzeta, psifac, pio2ka
c
c      double precision:: lref, tauRef, velRef
c      double precision:: Bref, pref, n0ref, rhoRef, Tref, TrefeV
c      double precision:: Mi, mu0, K2eV, RUniversal
      double precision:: pAvgDim, plasmaBeta
c     
c     
c     Generic constants
	namelist/molwts/wmol
	namelist/gammas/gamma
	namelist/diffusiveproperties/mu,eta,etaFactor,kappa
c
        open(16,file='prop.inp',form='formatted')
	read(16,molwts)
	read(16,gammas) 
	read(16,diffusiveproperties) 
        close(16)
        if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
           write(6,molwts)
           write(6,gammas)
           write(6,diffusiveproperties)
        endif
c
c     Init ux to zero
        ux = 0.0d0
c
c     Init Phi to zero
        
c
        pi=4.D0*datan(1.D0)
        pisq=pi*pi
c
c
c     Maximum toroidal field in Tesla
      BToroidalMax=0.2D0
c     Reference quantities
      wmol=1.D0
      RUniversal=8314.5D0
c      rgas=RUniversal/wmol
c     This is a non-dim gas constant - Eos is p=rho T
      rgas=1.D0
c     Permeability of free space in Henry/m or kg-m/Coulomb^2
      mu0=4.D0*pi * 1.D-07
c
c     Assume the following plasma beta
      plasmaBeta=0.1D0
c     Reference Magnetic field in Tesla (kg/s/coulomb)
      Bref=0.02D0
      Bref=BToroidalMax/10.D0
c      Bref=BToroidalMax
c
c     Alternative pref is Bref^2/mu0 (in N/m^2)
      pref=Bref**2/mu0
c     Reference no. of molecules /m^3
      n0ref=1.5D19
c     Mass of each proton in Kg
      Mi=1.6605 D-27
c     Reference density rhoRef=n0ref*Mi*wmol
c      rhoRef=1.6605*wmol*1.D-08
      rhoRef=Mi*wmol*n0Ref
c     Reference temperature in K
      TRef=pref/(rhoRef*RUniversal/wmol)
c     Convert K to eV
      K2eV=1.1604D04
c     Reference temperature in eV
      TRefeV=Tref/K2eV
c
c     Reference velocity is the alfven velocity defined 
c        as Bref/sqrt(rhoRef*mu0)
      velRef=Bref/dsqrt(rhoRef*mu0)
c     
c     Reference length scale is the 1/2 the width of the tokamak in m
      lref=0.2D0
c
c     Reference time scale tauRef
      tauRef=lref/velRef
c
c     Assume the following plasma beta
      plasmaBeta=0.1D0
c     Initial average pressure should be 0.1 atm
      pAvgDim=rhoRef*250.D0*K2eV*RUniversal/Wmol
c      pAvgDim=0.01D0/2.D0*BToroidalMax**2/mu0
      pAvgDim=plasmaBeta/2.D0*BToroidalMax**2/mu0
      pAvg=pAvgDim/pref
c
      alpha=0.5D0
c      psi0=1.D0
c      psi0=0.1D0
c      psi0=2.D0/pi*dsqrt(pAvg/alpha)
      psi0 = 2.D0/pi
c
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(6,*) 'Reference Quantities'
         write(6,*) 'Density', rhoRef,n0Ref
         write(6,*) 'Pressure', pRef
         write(6,*) 'Temperature', TRef, TRefeV
         write(6,*) 'Velocity', velRef
         write(6,*) 'Length', lref
         write(6,*) 'time', tauRef
c     
         write(6,*) 'pAvg', pAvg, pAvgDim, pAvgDim*2.D0/
     &        (BToroidalMax)**2*mu0
         write(6,*) 'T Average (eV)', pAvg*TRefeV
      endif
c
      aspect=1.D0
      a=1.D0
      pio2ka=0.5D0*pi/(aspect*a)
      pfac=0.5D0*alpha*pi*pi/(a*a)
c      write(6,*) '2pfac=',2*pfac,pi*pi/(4.D0*a*a)+pi*pi/(4.D0*a*a*aspect**2)
c
c      
c$$$      fun=0.D0
c$$$      dfdxi=0.D0
c$$$      d2fdxi2=0.D0
c
      do i=ixlo,ixhi,1
         x = xc(i)
         fun(i)=-dsin(0.5D0*pi/a*x)
         dfdxi(i)=-0.5D0*pi/a*dcos(0.5D0*pi/a*x)
         d2fdxi2(i)=(0.5D0*pi/a)**2*dsin(0.5D0*pi/a*x)
      enddo
c
      g0=BToroidalMax/Bref
c      write(6,*) 'Max value of |fun|=',fmin
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(6,*) 'g0=',g0
      endif
c
c     Constant density, zero velocity
      rho=1.D0
      u=0.D0
      v=0.D0
      w=0.D0
c
c     constant Toroidal field
      bk=g0
c
      do k=izlo,izhi,1
         z = zc(k)
         do j=iylo,iyhi,1
            y = yc(j)
            zeta=pio2ka*y
            coszeta=dcos(zeta)
            sinzeta=dsin(zeta)
            do i=ixlo,ixhi,1
c
               x = xc(i)
               xi=pio2ka*x
c
               psi=psi0*fun(i)*coszeta
               dpsidx=psi0*dfdxi(i)*coszeta
               d2psidx2=psi0*d2fdxi2(i)*coszeta
               dpsidy=-psi0*pio2ka*fun(i)*sinzeta
               d2psidy2=-psi0*pio2ka**2*fun(i)*coszeta
c
               press=pavg/2.D0+pfac*psi*psi
               bi=-dpsidy
               bj=dpsidx
c
c
               ux(i,j,k,1)=rho
               ux(i,j,k,2)=rho*u
               ux(i,j,k,3)=rho*v
               ux(i,j,k,4)=rho*w
               ux(i,j,k,5)=Bi
               ux(i,j,k,6)=Bj
               ux(i,j,k,7)=Bk
               ux(i,j,k,8)=press/(gamma-1.D0)+
     &              0.5*rho*(u*u+v*v+w*w)+0.5D0*(bi*bi+bj*bj+bk*bk)
c               write(6,*) 'INIT EQCHECK',i,j,k,d2psidx2+d2psidy2
c     &              +2.D0*pfac*psi, bi,bj
c$$$               write(6,*) 'INIT',i,j,ux(i,j,k,5),ux(i,j,k,6),
c$$$     &              ux(i,j,k,8)
            enddo
         enddo
      enddo
c
      plasmaBeta=0.D0
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               Bi=ux(i,j,k,5)
               Bj=ux(i,j,k,6)
               Bk=ux(i,j,k,7)
               press=(gamma-1.D0)*
     &              (ux(i,j,k,8)-0.5D0*(bi*bi+bj*bj+bk*bk))
               plasmaBeta=plasmaBeta+press/(0.5D0*(bi*bi+bj*bj+bk*bk))
            enddo
         enddo
      enddo
      plasmaBeta=plasmaBeta/(nxlocal*nylocal*nzlocal)
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(6,*) 'Plasma Beta=',plasmaBeta
      endif
      return
      end










