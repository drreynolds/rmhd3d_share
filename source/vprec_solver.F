c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c     This file implements a solver for the viscous preconditioning 
c     components, which calls multigrid solver modules from the C 
c     routines defined in the file viscous_prec.h.  These routines
c     handle the block-diagonal solves in the viscous system.  The
c     lower-triangular block products are computed via a standard 
c     finite-difference Jacobian-vector product, albeit used in a 
c     non-standard way.  For discussion of the solver, see the 
c     paper.
c
c-----------------------------------------------------------------------
c     $Log: vprec_solver.F,v $
c=======================================================================




      subroutine vprec_solver(uu, gamdt, delta, ww, vtmp, kier)
c-----------------------------------------------------------------------
c     Description: vprec_solver provides the overall solution driver 
c        for the block lower-triangular viscous preconditioning 
c        system:  P z = r
c
c     Arguments:
c            uu - (dbl(*), input) current state of system
c         gamdt - (dbl, input) scalar appearing in the Newton Matrix
c                  M = I+gamdt*J
c         delta - (dbl, input) desired tolerance for precond. solve
c                  0.0 results in using HYPRE default
c            ww - (dbl(*), in/out) on input the rhs vector r, and on 
c                  output the solution vector z
c          vtmp - (dbl(*)) scratch space with same size as uu
c          kier - (int, output) return flag: 0=>success, 
c                  1=>recoverable error, -1=>non-recoverable error
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: uu, ww, vtmp
      double precision :: gamdt, delta
      integer          :: kier

c     internal variables
#ifdef DYNAMIC
      double precision, allocatable, dimension(:,:,:,:)
     &                 :: fu, ftmp1, ftmp2, wwsave
#else
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &                 :: fu, ftmp1, ftmp2, wwsave
#endif
      double precision :: sigma, dtmp, utw, wtw, wnorm
      integer          :: i, j, k, l, ier
      logical          :: debug=.false.

c======= Internals ============

      if (debug)  write(6,'(2(A, es12.5))') 
     &     '\n  Entering vprec_solver routine, gamdt =',
     &     gamdt,', delta =', delta

#ifdef DYNAMIC
c     allocate temporary arrays
      if (.not. allocated(fu)) then
         allocate(fu(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      endif
      if (.not. allocated(ftmp2)) then
         allocate(ftmp2(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      endif
      if (debug) then
         if (.not. allocated(ftmp1)) then
            allocate(ftmp1(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
         endif
         if (.not. allocated(wwsave)) then
            allocate(wwsave(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
         endif
      endif
#endif

c     initialize return flag to success
      kier = 0

c     put the viscous flux divergence of uu into fu
      call CommState(uu)
      call divFlux_viscous(uu,fu)

c     no need to update ww_rho, since that is multiplied by I


c     solve momentum system to update ww_u
      if (debug) then
         if (iprocx==1 .and. iprocy==1 .and. iprocz==1) 
     &        write(6,*) '\n  vprec_solver: calling FVPrecDuSolve'
      endif
      
c        check initial residual of Du system (rhs in ftmp1, sol in ftmp2)
      if (debug) then
         wwsave = ww
         ftmp1 = wwsave
         ftmp2 = ww
c         call vprec_Du_mult(uu, gamdt, ftmp2, fu, ier)
         call vprec_Du_mult_hypre(uu, gamdt, ftmp2, vtmp, ier)
         vtmp = ftmp2-ftmp1
         dtmp = sqrt(sum(vtmp(1:nxlocal,1:nylocal,1:nzlocal,2:4)
     &        *vtmp(1:nxlocal,1:nylocal,1:nzlocal,2:4)))/
     &        sqrt(sum(ftmp1(1:nxlocal,1:nylocal,1:nzlocal,2:4)
     &        *ftmp1(1:nxlocal,1:nylocal,1:nzlocal,2:4)))
         write(0,'(A,es15.8)') 
     &        '    vprec_solver: Du initial residual error =',dtmp
      endif

c        solve the Du system
      call FVPrecDuSolve(ww, ww, ftmp2, delta, ier)
      if (ier /= 0) then
         kier = ier
         write(0,*) '    vprec_solver: FVPrecDuSolve error =',kier
         return
      endif

c        check linear residual of Du system (rhs in ftmp1, sol in ftmp2)
      if (debug) then
         ftmp1 = wwsave
         ftmp2 = ww
c         call vprec_Du_mult(uu, gamdt, ftmp2, fu, ier)
         call vprec_Du_mult_hypre(uu, gamdt, ftmp2, vtmp, ier)
         vtmp = ftmp2-ftmp1
         dtmp = sqrt(sum(vtmp(1:nxlocal,1:nylocal,1:nzlocal,2:4)
     &        * vtmp(1:nxlocal,1:nylocal,1:nzlocal,2:4)))/
     &        sqrt(sum(ftmp1(1:nxlocal,1:nylocal,1:nzlocal,2:4)
     &        * ftmp1(1:nxlocal,1:nylocal,1:nzlocal,2:4)))
         write(0,'(A,es15.8)') 
     &        '    vprec_solver: Du residual error =',dtmp
      endif



c     solve magnetic field system to update ww_b
      if (debug) then
         if (iprocx==1 .and. iprocy==1 .and. iprocz==1) 
     &        write(6,*) '\n  vprec_solver: calling FVPrecDbSolve'
      endif

c        check initial residual of Db system (rhs in ftmp1, sol in ftmp2)
      if (debug) then
         ftmp1 = wwsave
         ftmp2 = ww
c         call vprec_Db_mult(uu, gamdt, ftmp2, fu, ier)
         call vprec_Db_mult_hypre(uu, gamdt, ftmp2, vtmp, ier)
         vtmp = ftmp2-ftmp1
         dtmp = sqrt(sum(vtmp(1:nxlocal,1:nylocal,1:nzlocal,5:7)
     &        * vtmp(1:nxlocal,1:nylocal,1:nzlocal,5:7)))/
     &        sqrt(sum(ftmp1(1:nxlocal,1:nylocal,1:nzlocal,5:7)
     &        * ftmp1(1:nxlocal,1:nylocal,1:nzlocal,5:7)))
         write(0,'(A,es15.8)') 
     &        '    vprec_solver: Db initial residual error =',dtmp
      endif

c        solve the Db system
      call FVPrecDbSolve(ww, ww, ftmp2, delta, ier)
      if (ier /= 0) then
         kier = ier
         write(0,*) '    vprec_solver: FVPrecDbSolve error =',kier
         return
      endif

c        check linear residual of Db system (rhs in ftmp1, sol in ftmp2)
      if (debug) then
         ftmp1 = wwsave
         ftmp2 = ww
c         call vprec_Db_mult(uu, gamdt, ftmp2, fu, ier)
         call vprec_Db_mult_hypre(uu, gamdt, ftmp2, vtmp, ier)
         vtmp = ftmp2-ftmp1
         dtmp = sqrt(sum(vtmp(1:nxlocal,1:nylocal,1:nzlocal,5:7)
     &        * vtmp(1:nxlocal,1:nylocal,1:nzlocal,5:7)))/
     &        sqrt(sum(ftmp1(1:nxlocal,1:nylocal,1:nzlocal,5:7)
     &        * ftmp1(1:nxlocal,1:nylocal,1:nzlocal,5:7)))
         write(0,'(A,es15.8)') 
     &        '    vprec_solver: Db residual error =',dtmp
      endif


c     update ww_e with contributions from other components:
c     ww_e = ww_e - L_rho*ww_rho - L_u*ww_u - L_b*ww_b
c     [this update, or at least the block products, are 
c     accomplished via finite-differencing the viscous divflux 
c     terms with an appropriate perturbation, and using the 
c     resulting energy difference value]
      if (debug) then
         if (iprocx==1 .and. iprocy==1 .and. iprocz==1) 
     &        write(6,*) '\n  vprec_solver: updating ww_e'
      endif

c        set sigma, the differencing parameter 
      call fnvmhddotprod(uu,ww,utw)
      call fnvmhddotprod(ww,ww,wtw)
      call fnvmhdl1norm(ww,wnorm)
      sigma = sign(1.d-8,utw)*max(utw,wnorm)/wtw
      if (wtw < 1.d-32)  sigma = 1.d-8

c        set the perturbed state into vtmp
      do l=1,7,1
         do k=izlo,izhi,1
            do j=iylo,iyhi,1
               do i=ixlo,ixhi,1
                  vtmp(i,j,k,l) = uu(i,j,k,l) + sigma*ww(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               vtmp(i,j,k,8) = uu(i,j,k,8)
            enddo
         enddo
      enddo

c        get the viscous flux divergence at the perturbed state -> ftmp2
      call CommState(vtmp)
      call divFlux_viscous(vtmp,ftmp2)

c        update ww_e with the F-D approximation to 
c          -gamdt*(L_who*ww_rho + L_u*ww_u + L_b*ww_b)
      do k=izlo,izhi,1
         do j=iylo,iyhi,1
            do i=ixlo,ixhi,1
               ww(i,j,k,8) = ww(i,j,k,8) 
     &              + gamdt/sigma*(ftmp2(i,j,k,8)-fu(i,j,k,8))
            enddo
         enddo
      enddo


c     solve energy system to update ww_e
      if (debug) then
         if (iprocx==1 .and. iprocy==1 .and. iprocz==1) 
     &        write(6,*) '\n  vprec_solver: calling FVPrecDeSolve'
      endif

c        check initial residual of De system (rhs in ftmp1, sol in ftmp2)
      if (debug) then
         wwsave(:,:,:,8) = ww(:,:,:,8)
         ftmp1 = wwsave
         ftmp2 = ww
c         call vprec_De_mult(uu, gamdt, ftmp2, fu, ier)
         call vprec_De_mult_hypre(uu, gamdt, ftmp2, vtmp, ier)
         vtmp = ftmp2-ftmp1
         dtmp = sqrt(sum(vtmp(1:nxlocal,1:nylocal,1:nzlocal,8)
     &        * vtmp(1:nxlocal,1:nylocal,1:nzlocal,8)))/
     &        sqrt(sum(ftmp1(1:nxlocal,1:nylocal,1:nzlocal,8)
     &        * ftmp1(1:nxlocal,1:nylocal,1:nzlocal,8)))
         write(0,'(A,es15.8)') 
     &        '    vprec_solver: De initial residual error =',dtmp
      endif

c        solve the De system
      call FVPrecDeSolve(ww, ww, ftmp2, delta, ier)
      if (ier /= 0) then
         kier = ier
         write(0,*) '    vprec_solver: FVPrecDeSolve error =',kier
         return
      endif

c        check linear residual of De system (rhs in ftmp1, sol in ftmp2)
      if (debug) then
         ftmp1 = wwsave
         ftmp2 = ww
c         call vprec_De_mult(uu, gamdt, ftmp2, fu, ier)
         call vprec_De_mult_hypre(uu, gamdt, ftmp2, vtmp, ier)
         vtmp = ftmp2-ftmp1
         dtmp = sqrt(sum(vtmp(1:nxlocal,1:nylocal,1:nzlocal,8)
     &        * vtmp(1:nxlocal,1:nylocal,1:nzlocal,8)))/
     &        sqrt(sum(ftmp1(1:nxlocal,1:nylocal,1:nzlocal,8)
     &        * ftmp1(1:nxlocal,1:nylocal,1:nzlocal,8)))
         write(0,'(A,es15.8)') 
     &        '    vprec_solver: De residual error =',dtmp
         write(0,*) '  '
      endif

#ifdef DYNAMIC
c     free temporary arrays
      if (allocated(fu)) then
         deallocate(fu)
      endif
      if (allocated(ftmp2)) then
         deallocate(ftmp2)
      endif
      if (allocated(ftmp1)) then
         deallocate(ftmp1)
      endif
      if (allocated(wwsave)) then
         deallocate(wwsave)
      endif
#endif


      return
      end subroutine vprec_solver
c=======================================================================
