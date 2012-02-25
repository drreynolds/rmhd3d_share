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

      subroutine Diagnostics(ux,phi,ttot)
c-----------------------------------------------------------------
c     Description: Diagnostics at regular intervals: e.g. 
c        Total energy, mass, kinetic energy, magnetic energy, etc.
c
c     Arguments:
c          ux - (in) current state of all variables over local mesh
c         phi - (in) projection information for div B
c        ttot - (in) current time
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: ttot
      double precision :: toteng
      double precision :: totmass
      double precision :: tke, tme, tie, flux, flux2, dpsidt
      double precision :: divBlocmax, divBlocmin, divBlocavg
      double precision :: divBmax, divBmin, divBavg
      integer          :: i, j, k
  

c======= Internals ============

c     'toteng' contains the total energy in the plasma
      call TotalEnergy(ux(ixlo,iylo,izlo,8),toteng)

c     'totmass' contains the total mass of the plasma
      call TotalEnergy(ux(ixlo,iylo,izlo,1),totmass)

c     'tke' contains the total kinetic energy in the plasma
      call TotalKineticEnergy(ux,tke)

c     'tme' contains the total magnetic energy in the plasma
      call TotalMagneticEnergy(ux,tme)

c     'flux' contains the reconnected flux
      call ReconnectedFlux(ux,flux)

c     'flux2' 
c      flux2 = -sum(ux(nxlocal/2,1:nylocal/2,1:nzlocal,5))*dy

c     'dpsidt' contains the reconnection rate
c      call ReconnectionRate(ux,dpsidt)

c     'tie' contains the total internal energy in the plasma
      tie = toteng - tke - tme
      
c     'div' contains the magnetic field divergence
      call DivergenceCD(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &                  ux(ixlo,iylo,izlo,7),div)	
      divBlocmax = maxval(abs(div))
      divBlocmin = minval(abs(div))
      divBlocavg = sum(abs(div))

c     get global max and min magnetic field divergence
#ifdef PARALLEL
      call MPI_Allreduce(divBlocmax, divBmax, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_MAX, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
      call MPI_Allreduce(divBlocmin, divBmin, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_MIN, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
      call MPI_Allreduce(divBlocavg, divBavg, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, comm3D, ierr)
      call ErrorHandler(ierr,ERROR_ALLREDUCE)
      divBavg = divBavg/nx/ny/nz
#else   
      divBmax = divBlocmax
      divBmin = divBlocmin
      divBavg = divBlocavg/nx/ny/nz
#endif

c     root processor outputs the max and min divergence to stdout
c        iostream 6  => istdout
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(6,'(3(A,es12.5))')  '   max(divB) =', divBmax,
     &        ',   min(divB) =', divBmin,',  avg(divB) =', divBavg
      endif

c     all processors write their flux2 value to stdout
c      write(6,*) 'flux2=',flux,flux2

c     root note outputs the global energy, mass and divB results
c        iostream 91 => idiag
c        iostream 92 => idiagCons
c        iostream 94 => idiagEng
c        iostream 95 => idiagflux
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(91,'(3(es12.5,1x))') divBmax, divBmin, divBavg
         write(92,'(3(es16.9,1x))') ttot, toteng, totmass
         write(94,'(4(es16.9,1x))') ttot, tke, tme, tie
         write(95,'(2(es16.9,1x))') ttot, flux
      endif

      return
      end
c=================================================================



      subroutine DivergenceCD(Bx,By,Bz,div)
c-----------------------------------------------------------------
c     Description: Calculate div B on local domain
c
c     Arguments:
c          Bx - (in) x-component of magnetic field
c          By - (in) y-component of magnetic field
c          Bz - (in) z-component of magnetic field
c         div - (out) divergence of magnetic field
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      double precision :: Bx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: By(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: Bz(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: div(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI)
      double precision :: dxi, dyi, dzi, d1, d2
      integer          :: i, j, k


c======= Internals ============

c     initialize internal derivative fractions, div
      dxi = 0.5D0/dx;  dyi = 0.5D0/dy;  dzi = 0.5D0/dz
      d1  = 2.d0/3.d0; d2  = -1.d0/12.d0
      div = 0.D0

c     2-D calculation of the local divergence
#ifdef TWO_D
      do k = 1,nzlocal,1	
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
#ifdef CD4
c               4th-order CD calculation of Div(B)
               div(i,j,k) = d1*dxi*(Bx(i+1,j,k)-Bx(i-1,j,k)) +
     &              d1*dyi*(By(i,j+1,k)-By(i,j-1,k)) +
     &              d2*dxi*(Bx(i+2,j,k)-Bx(i-2,j,k)) +
     &              d2*dyi*(By(i,j+2,k)-By(i,j-2,k)) 
#else
c               2nd-order CD calculation of Div(B)
               div(i,j,k) = dxi*(Bx(i+1,j,k)-Bx(i-1,j,k)) +
     &                      dyi*(By(i,j+1,k)-By(i,j-1,k))
#endif
            enddo
         enddo
      enddo

c     3-D calculation of the local divergence
#else
      do k = 1,nzlocal,1	
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
#ifdef CD4
               div(i,j,k)=d1*dxi*(Bx(i+1,j,k)-Bx(i-1,j,k))+
     &              d1*dyi*(By(i,j+1,k)-By(i,j-1,k))+
     &              d1*dzi*(Bz(i,1,k+1)-Bz(i,j,k-1))+
     &              d2*dxi*(Bx(i+2,j,k)-Bx(i-2,j,k))+
     &              d2*dyi*(By(i,j+2,k)-By(i,j-2,k))+
     &              d2*dzi*(Bz(i,j,k+2)-Bz(i,j,k-2))
#else
               div(i,j,k)=dxi*(Bx(i+1,j,k)-Bx(i-1,j,k))+
     &                    dyi*(By(i,j+1,k)-By(i,j-1,k))+
     &                    dzi*(Bz(i,j,k+1)-Bz(i,j,k-1))
#endif
            enddo
         enddo
      enddo
#endif

      return
      end
c=================================================================



      subroutine ReconnectionRate(ux,dpsidt)
c-----------------------------------------------------------------
c     Description: Calculate the reconnection rate
c
c     Arguments:
c            ux - (in) vector of conserved values
c        dpsidt - (out) reconnection rate
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      implicit none
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: dpsidt,dpsidtLocal,fluxLocal,flux
      integer:: jloc, iloc, k

c======= Internals ============
      call jlocal(ny/2,iprocy,jloc)
      call ilocal(nx/2,iprocx,iloc)
c     
      k=1
      write(6,*) 'PSIDOT',iprocx,iprocy,iloc,jloc
c     
      if((iloc.lt.1.or.iloc.gt.nxlocal).and. 
     &     (jloc.lt.1.or.jloc.gt.nylocal)) then 
         dpsidtlocal=0.D0
      else
         dpsidtlocal=-eta*0.5D0*(
     &        ((ux(iloc+1,jloc,k,6)+ux(iloc+1,jloc+1,k,6))-
     &        (ux(iloc,jloc,k,6)+ux(iloc,jloc+1,k,6)))/dx-
     &        ((ux(iloc,jloc+1,k,5)+ux(iloc+1,jloc+1,k,5))-
     &        (ux(iloc,jloc,k,5)+ux(iloc+1,jloc,k,5)))/dy)
c     
      endif
c     
c     
      dpsidt=0.D0
#ifdef PARALLEL
      call MPI_Allreduce(fluxLocal, flux, 1, MPI_DOUBLE_PRECISION,
     &     MPI_SUM, comm3D, ierr)
      Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else   
      dpsidt=dpsidtLocal
#endif  
c     
      return	
      end
c=================================================================




      subroutine ReconnectedFLux(ux,fluxtot)
c-----------------------------------------------------------------
c     Description: Calculate total energy on entire domain
c     
c     Arguments:
c           eng - (in) vector energy information on local mesh
c        toteng - (out) total energy on entire domain
c-----------------------------------------------------------------
c======= Inclusions ===========
       use mesh
       use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
      double precision:: flux(1:nx), fluxtot
      double precision:: fluxtotLocal
      integer:: jg, jl
      integer:: ilo,ihi, ipy


c======= Internals ============
      jg=ny/2
c      write(6,*) 'FLUX0',iprocx,iprocy
      call jlocal(jg,ipy,jl)
      call iglobal(iprocx,1,ilo)
      call iglobal(iprocx,nxlocal,ihi)
c      write(6,*) 'FLUX',iprocx,iprocy,ipy,ilo,ihi,jl
      if(ipy.eq.iprocy) then
         flux((iprocx-1)*NXlocal + 1:iprocx*Nxlocal)=
     &        ux(1:nxlocal,jl,1,6)
      else
         flux(ilo:ihi)=0.D0
      endif
      fluxtotLocal=sum(flux(1:nx/2))
#ifdef PARALLEL
      call MPI_Allreduce(fluxtotLocal, fluxtot,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM, comm3D, ierr)
      Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else   
      fluxtot=fluxtotLocal
#endif  
c     
      fluxtot=-fluxtot*dx
      return	
      end
c=================================================================



      subroutine TotalEnergy(eng,toteng)
c-----------------------------------------------------------------
c     Description: Calculate total energy on entire domain
c
c     Arguments:
c           eng - (in) vector energy information on local mesh
c        toteng - (out) total energy on entire domain
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      double precision :: eng(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: toteng
      double precision :: totengLocal


c======= Internals ============

c     initialize total energy
      toteng = 0.D0

c     obtain total energy on local domain
      totengLocal = sum(eng(1:NXlocal,1:NYlocal,1:NZlocal))

c     if there are multiple domains, combine all energy into one
#ifdef PARALLEL
      call MPI_Allreduce(totengLocal, toteng, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, comm3D, ierr)
      Call ErrorHandler(ierr,ERROR_ALLREDUCE)

c     otherwise, set output variable
#else   
      toteng = totengLocal
#endif

      return	
      end
c=================================================================



      subroutine TotalKineticEnergy(ux,tke)
c-----------------------------------------------------------------
c     Description: Calculate total kinetic energy on entire domain
c
c     Arguments:
c         ux - (in) all state information on local mesh 
c        tke - (out) total kinetic energy on entire domain
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: tke
      double precision :: tkeLocal


c======= Internals ============

c     initialize total kinetic energy
      tke = 0.D0

c     calculate local kinetic energy
      tkeLocal = 0.5D0*sum((ux(1:NXlocal,1:NYlocal,1:NZlocal,2)**2 +
     &                      ux(1:NXlocal,1:NYlocal,1:NZlocal,3)**2 +
     &                      ux(1:NXlocal,1:NYlocal,1:NZlocal,4)**2)
     &                     / ux(1:NXlocal,1:Nylocal,1:Nzlocal,1))

c     if there are multiple domains, combine all energy into one
#ifdef PARALLEL
      call MPI_Allreduce(tkeLocal, tke, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, comm3D, ierr)
      Call ErrorHandler(ierr,ERROR_ALLREDUCE)

c     otherwise, set output variable
#else   
      tke = tkeLocal
#endif  

      return	
      end
c=================================================================



      subroutine TotalMagneticEnergy(ux,tme)
c-----------------------------------------------------------------
c     Description: Calculate total magnetic energy on entire domain
c
c     Arguments:
c         ux - (in) all state information on local mesh 
c        tme - (out) total magnetic energy on entire domain
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: tme
      double precision :: tmeLocal


c======= Internals ============     

c     initialize total magnetic energy
      tme = 0.D0

c     calculate local magnetic energy
      tmeLocal = 0.5D0*sum((ux(1:NXlocal,1:NYlocal,1:NZlocal,5)**2 +
     &                      ux(1:NXlocal,1:NYlocal,1:NZlocal,6)**2 +
     &                      ux(1:NXlocal,1:NYlocal,1:NZlocal,7)**2))

c     if there are multiple domains, combine all energy into one
#ifdef PARALLEL
      call MPI_Allreduce(tmeLocal, tme, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, comm3D, ierr)
      Call ErrorHandler(ierr,ERROR_ALLREDUCE)

c     otherwise, set output variable
#else   
      tme = tmeLocal
#endif  

      return	
      end
c=================================================================



      subroutine ToroidalCurrent(Bx,By,Bz,Jcur)
c-----------------------------------------------------------------
c     Description: Calculate toroidal current on local domain
c
c     Arguments:
c          Bx - (in) x-component of magnetic field 
c          By - (in) y-component of magnetic field 
c          Bz - (in) z-component of magnetic field 
c        Jcur - (out) local toroidal current
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      double precision :: Bx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: By(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: Bz(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: Jcur(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: dxi, dyi, dzi
      integer          :: i, j, k


c======= Internals ============

c     initialize internal derivative fractions, Jcur
      dxi = 0.5D0/dx;  dyi = 0.5D0/dy;  dzi = 0.5D0/dz
      Jcur = 0.D0

c     2-D calculation of the local toroidal current
#ifdef TWO_D
      do k = 1,nzlocal,1	
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               Jcur(i,j,k) = dyi*(Bx(i,j+1,k)-Bx(i,j-1,k)) -
     &                       dxi*(By(i+1,j,k)-By(i-1,j,k))
            enddo
         enddo
      enddo
      
c     3-D calculation of the local toroidal current
#else
      do k = 1,nzlocal,1	
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               Jcur(i,j,k) = dxi*(Bx(i+1,j,k)-Bx(i-1,j,k)) +
     &                       dyi*(By(i,j+1,k)-By(i,j-1,k)) +
     &                       dzi*(Bz(i,j,k+1)-Bz(i,j,k-1))
            enddo
         enddo
      enddo
#endif
      
      return
      end
c=================================================================
