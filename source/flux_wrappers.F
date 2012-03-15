c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2012
c     All Rights Reserved
c=======================================================================


      subroutine InviscidFluxWrapper(finv,vv,idir,dt)
c-----------------------------------------------------------------------
c     Description: InviscidFluxWrapper calls the requested inviscid
c        flux routine, as determined by preprocessor directives.
c
c     Arguments:
c         finv - (dbl(*), output) flux
c           vv - (dbl(*), input)  primitive variables over domain
c         idir - (integer, input) flux direction
c           dt - (dbl, input)     time step size
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none

c     calling variables
      real*8,  intent(in)  ::   vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      real*8,  intent(out) :: finv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      real*8,  intent(in)  :: dt
      integer, intent(in)  :: idir

c======= Internals ============


#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
      call InviscidFluxRP(finv,vv,idir,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(finv,vv,idir,dt)
#else
c     using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(finv,vv,idir,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
c     2nd-order tuned central differences
      call InviscidFluxTCD(finv,vv,idir,dt)
#elif CD4
c     4th-order standard central differences
      call InviscidFluxCD4(finv,vv,idir,dt)
#elif ZIP
c     2nd-order ZIP difference
      call InviscidFluxZip(finv,vv,idir,dt)
#else
c     2nd-order standard central differences
      call InviscidFluxCD2(finv,vv,idir,dt)
#endif
#endif

      return
      end subroutine InviscidFluxWrapper
c-----------------------------------------------------------------------



      subroutine ViscousFluxWrapper(fvis, vv, dudx, dbdx, idir, dt, 
     &                              indx, indy, indz)
c-----------------------------------------------------------------------
c     Description: InviscidFluxWrapper calls the requested inviscid
c        flux routine, as determined by preprocessor directives.
c
c     Arguments:
c         fvis - (dbl(*), output) flux
c           vv - (dbl(*), input)  primitive variables over domain
c         dudx - (dbl(*), input)  velocity gradient
c         dbdx - (dbl(*), input)  magnetic field gradient
c         idir - (integer, input) flux direction
c           dt - (dbl, input)     time step size
c         ind* - (integer, input) indices for flux calculations
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh

c======= Declarations =========
      implicit none

c     calling variables
      real*8,  intent(in)  ::   vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      real*8,  intent(in)  :: dudx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3)
      real*8,  intent(in)  :: dbdx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3)
      real*8,  intent(out) :: fvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      real*8,  intent(in)  :: dt
      integer, intent(in)  :: idir, indx, indy, indz

c======= Internals ============

#ifndef NOVISCFLUX
#ifdef TCD
c     2nd-order tuned central differences
      call ViscousFluxTCD(fvis,vv,dudx,dbdx,idir,dt,indx,indy,indz)
#elif CD4
c     4th-order standard central differences
      call ViscousFluxCD4(fvis,vv,dudx,dbdx,idir,dt,indx,indy,indz)
#else
c     2nd-order standard central differences
      call ViscousFluxCD2(fvis,vv,dudx,dbdx,idir,dt,indx,indy,indz)
#endif
#endif

      return
      end subroutine ViscousFluxWrapper
c=======================================================================
