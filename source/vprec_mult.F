c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c     This file uses the components from vprec_solver.F to perform a 
c     (viscous Jacobian)-vector product, which calls HYPRE modules 
c     from the C routines defined in the file viscous_prec.c.  
c-----------------------------------------------------------------------
c     $Log: vprec_mult.F,v $
c=======================================================================




      subroutine vprec_mult_hypre(uu, gamdt, ww, utmp, kier)
c-----------------------------------------------------------------------
c     Description: vprec_mult provides the driver for the block 
c     lower-triangular (viscous Jacobian)-vector product
c        system:  P*z = r
c     Calls HYPRE preconditioner module for matrix-vector product.
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector z, 
c                  and on output the vector r
c          utmp - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, utmp
      double precision :: gamdt
      integer          :: kier

c     internal variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: fu, ftmp
      double precision :: sigma, utw, wtw, wnorm
      integer          :: i, j, k, l, ier

c======= Internals ============

c     initialize return flag to success
      kier = 0

c     multiply energy block to update ww_e
c        put product into utmp
      call FVPrecDeMultiply(ww, ww, ftmp, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif

c     update ww_e with contributions from other components:
c     ww_e = De*ww_e + L_rho*ww_rho + L_u*ww_u + L_b*ww_b 
c     [this update, or at least the block products, are 
c     accomplished via finite-differencing the viscous divflux 
c     terms with an appropriate perturbation, and using the 
c     resulting energy difference value]


c        put the viscous flux divergence of uu into fu
      call CommState(uu)
      call divFlux_viscous(uu,fu)

c        set sigma, the differencing parameter
      call fnvmhddotprod(uu,ww,utw)
      call fnvmhddotprod(ww,ww,wtw)
      call fnvmhdl1norm(ww,wnorm)
      sigma = sign(1.d-8,utw)*max(utw,wnorm)/wtw
      if (wtw < 1.d-32)  sigma = 1.d-8

c        set the perturbed state into utmp
      do l=1,7,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = uu(i,j,k,l) + sigma*ww(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               utmp(i,j,k,8) = uu(i,j,k,8)
            enddo
         enddo
      enddo
         
c        get the viscous flux divergence at the perturbed state -> ftmp
      call CommState(utmp)
      call divFlux_viscous(utmp,ftmp)

c        update ww_e with the F-D approximation to 
c           -gamdt*(L_rho*ww_rho + L_u*ww_u + L_b*ww_b)
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               ww(i,j,k,8) = ww(i,j,k,8) 
     &              - gamdt/sigma*(ftmp(i,j,k,8)-fu(i,j,k,8))
            enddo
         enddo
      enddo


c     multiply magnetic field block to update ww_b
c        put product into utmp
      call FVPrecDbMultiply(ww, ww, ftmp, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif

c     multiply momentum block to update ww_u
c        put product into utmp
      call FVPrecDuMultiply(ww, ww, ftmp, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif
      

      return
      end subroutine vprec_mult_hypre
c=======================================================================




      subroutine vprec_mult(uu, gamdt, ww, utmp, kier)
c-----------------------------------------------------------------------
c     Description: vprec_mult provides the driver for the block 
c     lower-triangular (viscous Jacobian)-vector product
c        system:  P*z = r
c     Performs matrix-vector product via finite-differences.
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector z, 
c                  and on output the vector r
c          utmp - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, utmp
      double precision :: gamdt
      integer          :: kier

c     internal variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: fu, ftmp
      double precision :: sigma, utw, wtw, wnorm
      integer          :: i, j, k, l, ier

c======= Internals ============

c     initialize return flag to success
      kier = 0

c     put the viscous flux divergence of uu into fu
      call CommState(uu)
      call divFlux_viscous(uu,fu)

c     multiply energy block to update ww_e
      call vprec_De_mult(uu, gamdt, ww, fu, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif

c     update ww_e with contributions from other components:
c     ww_e = De*ww_e + L_rho*ww_rho + L_u*ww_u + L_b*ww_b 
c     [this update, or at least the block products, are 
c     accomplished via finite-differencing the viscous divflux 
c     terms with an appropriate perturbation, and using the 
c     resulting energy difference value]

c        set sigma, the differencing parameter 
      call fnvmhddotprod(uu,ww,utw)
      call fnvmhddotprod(ww,ww,wtw)
      call fnvmhdl1norm(ww,wnorm)
      sigma = sign(1.d-8,utw)*max(utw,wnorm)/wtw
      if (wtw < 1.d-32)  sigma = 1.d-8

c        set the perturbed state into utmp
      do l=1,7,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = uu(i,j,k,l) + sigma*ww(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               utmp(i,j,k,8) = uu(i,j,k,8)
            enddo
         enddo
      enddo
      
c        get the viscous flux divergence at the perturbed state -> ftmp
      call CommState(utmp)
      call divFlux_viscous(utmp,ftmp)

c        update ww_e with the F-D approximation to 
c           -gamdt*(L_rho*ww_rho + L_u*ww_u + L_b*ww_b)
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               ww(i,j,k,8) = ww(i,j,k,8) 
     &              - gamdt/sigma*(ftmp(i,j,k,8)-fu(i,j,k,8))
            enddo
         enddo
      enddo


c     multiply magnetic field block to update ww_b
      call vprec_Db_mult(uu, gamdt, ww, fu, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif

c     multiply momentum block to update ww_u
      call vprec_Du_mult(uu, gamdt, ww, fu, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif


      return
      end subroutine vprec_mult
c=======================================================================




      subroutine vprec_Du_mult_hypre(uu, gamdt, ww, utmp, kier)
c-----------------------------------------------------------------------
c     Description: vprec_Du_mult provides the driver for a 
c     finite-difference approximation to the single-block 
c     multiplication yu = Du*xu 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector xu
c                  and on output the vector yu
c          utmp - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, utmp
      double precision :: gamdt
      integer          :: kier

c     internal variables
      integer          :: ier

c======= Internals ============
c
c     initialize return flag to success
      kier = 0
c
c     multiply momentum block to update ww_u
      call FVPrecDuMultiply(ww, ww, utmp, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif
c
      return
      end subroutine vprec_Du_mult_hypre
c=======================================================================




      subroutine vprec_Du_mult(uu, gamdt, ww, fu, kier)
c-----------------------------------------------------------------------
c     Description: vprec_Du_mult provides the driver for a 
c     finite-difference approximation to the single-block 
c     multiplication yu = Du*xu 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector xu
c                  and on output the vector yu
c            fu - (dbl(*)) viscous fluxes evaluated at uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, fu
      double precision :: gamdt
      integer          :: kier

c     internal variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: utmp, ftmp
      integer          :: i, j, k, l, ier
      double precision :: sigma, utw, wtw, wnorm

c======= Internals ============

c     initialize return flag to success
      kier = 0

c     set sigma, the differencing parameter 
      call fnvmhddotprod_u(uu,ww,utw)
      call fnvmhddotprod_u(ww,ww,wtw)
      call fnvmhdl1norm_u(ww,wnorm)
      sigma = sign(1.d-8,utw)*max(utw,wnorm)/wtw
      if (wtw < 1.d-32)  sigma = 1.d-8

c     set the perturbed state into utmp
      do l=1,8,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = uu(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      do l=2,4,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = utmp(i,j,k,l) + sigma*ww(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

c     get the viscous flux divergence at the perturbed state -> ftmp
      call CommState(utmp)
      call divFlux_viscous(utmp,ftmp)

c     update ww_u with the F-D approximation to (- gamdt*Du*ww_u)
      do l=2,4,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  ww(i,j,k,l) = ww(i,j,k,l) 
     &                 - gamdt/sigma*(ftmp(i,j,k,l)-fu(i,j,k,l))
               enddo
            enddo
         enddo
      enddo


      return
      end subroutine vprec_Du_mult
c=======================================================================




      subroutine vprec_Db_mult_hypre(uu, gamdt, ww, utmp, kier)
c-----------------------------------------------------------------------
c     Description: vprec_Db_mult provides the driver for a 
c     finite-difference approximation to the single-block 
c     multiplication yb = Db*xb 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector xb
c                  and on output the vector yb
c          utmp - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, utmp
      double precision :: gamdt
      integer          :: kier

c     internal variables
      integer          :: ier

c======= Internals ============
c
c     initialize return flag to success
      kier = 0
c
c     multiply magnetic field block to update ww_b
      call FVPrecDbMultiply(ww, ww, utmp, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif
c
      return
      end subroutine vprec_Db_mult_hypre
c=======================================================================




      subroutine vprec_Db_mult(uu, gamdt, ww, fu, kier)
c-----------------------------------------------------------------------
c     Description: vprec_Db_mult provides the driver for a 
c     finite-difference approximation to the single-block 
c     multiplication yb = Db*xb 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector xb
c                  and on output the vector yb
c            fu - (dbl(*)) viscous fluxes evaluated at uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, fu
      double precision :: gamdt
      integer          :: kier

c     internal variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: utmp, ftmp
      integer          :: i, j, k, l, ier
      double precision :: sigma, utw, wtw, wnorm

c======= Internals ============

c     initialize return flag to success
      kier = 0

c     set sigma, the differencing parameter 
      call fnvmhddotprod_B(uu,ww,utw)
      call fnvmhddotprod_B(ww,ww,wtw)
      call fnvmhdl1norm_B(ww,wnorm)
      sigma = sign(1.d-8,utw)*max(utw,wnorm)/wtw
      if (wtw < 1.d-32)  sigma = 1.d-8

c     set the perturbed state into utmp
      do l=1,8,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = uu(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      do l=5,7,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = utmp(i,j,k,l) + sigma*ww(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

c     get the viscous flux divergence at the perturbed state -> ftmp
      call CommState(utmp)
      call divFlux_viscous(utmp,ftmp)

c     update ww_u with the F-D approximation to (- gamdt*Db*ww_b)
      do l=5,7,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  ww(i,j,k,l) = ww(i,j,k,l) 
     &              - gamdt/sigma*(ftmp(i,j,k,l)-fu(i,j,k,l))
               enddo
            enddo
         enddo
      enddo


      return
      end subroutine vprec_Db_mult
c=======================================================================




      subroutine vprec_De_mult_hypre(uu, gamdt, ww, utmp, kier)
c-----------------------------------------------------------------------
c     Description: vprec_De_mult provides the driver for a 
c     finite-difference approximation to the single-block 
c     multiplication ye = De*xe 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector xe
c                  and on output the vector ye
c          utmp - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, utmp
      double precision :: gamdt
      integer          :: kier

c     internal variables
      integer          :: ier

c======= Internals ============
c
c     initialize return flag to success
      kier = 0
c
c     multiply magnetic field block to update ww_b
      call FVPrecDeMultiply(ww, ww, utmp, ier)
      if (ier /= 0) then
         kier = ier
         return
      endif
c
      return
      end subroutine vprec_De_mult_hypre
c=======================================================================




      subroutine vprec_De_mult(uu, gamdt, ww, fu, kier)
c-----------------------------------------------------------------------
c     Description: vprec_De_mult provides the driver for a 
c     finite-difference approximation to the single-block 
c     multiplication ye = De*xe 
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c            ww - (dbl(*), in/out) on input the vector xe
c                  and on output the vector ye
c            fu - (dbl(*)) viscous fluxes evaluated at uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, fu
      double precision :: gamdt
      integer          :: kier

c     internal variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: utmp, ftmp
      integer          :: i, j, k, l, ier
      double precision :: sigma, utw, wtw, wnorm

c======= Internals ============

c     initialize return flag to success
      kier = 0

c     set sigma, the differencing parameter 
      call fnvmhddotprod_e(uu,ww,utw)
      call fnvmhddotprod_e(ww,ww,wtw)
      call fnvmhdl1norm_e(ww,wnorm)
      sigma = sign(1.d-8,utw)*max(utw,wnorm)/wtw
      if (wtw < 1.d-32)  sigma = 1.d-8

c     set the perturbed state into utmp
      do l=1,7,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  utmp(i,j,k,l) = uu(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               utmp(i,j,k,8) = uu(i,j,k,8) + sigma*ww(i,j,k,8)
            enddo
         enddo
      enddo

c     get the viscous flux divergence at the perturbed state -> ftmp
      call CommState(utmp)
      call divFlux_viscous(utmp,ftmp)

c     update ww_e with the F-D approximation to (- gamdt*De*ww_e)
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               ww(i,j,k,8) = ww(i,j,k,8) 
     &              - gamdt/sigma*(ftmp(i,j,k,8)-fu(i,j,k,8))
            enddo
         enddo
      enddo


      return
      end subroutine vprec_De_mult
c=======================================================================
