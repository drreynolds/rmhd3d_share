c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
c     Ravi Samtaney
c     KAUST, Mechanical Engineering
c
c     Daniel R. Reynolds
c     SMU, Mathematics
c
c     Copyright 2004
c     All Rights Reserved
c=======================================================================

      module petsc_data
c-----------------------------------------------------------------
c     Description: holds necessary timelevel data for nonlinear
c        function evaluation
c
c     Contains the internal variables:
c           vv - primitive variables
c         *inv - hyperbolic fluxes
c         *vis - parabolic fluxes
c         dudx - velocity gradient 
c         dbdx - magnetic field gradient 
c-----------------------------------------------------------------
      use mesh_parms

      implicit none
      save
      
      logical :: outnode
      integer :: binary_flag, ndump, output_flag, ndiag
#ifdef DYNAMIC
      double precision, dimension(:,:,:), allocatable ::
     &     phi
      double precision, dimension(:,:,:,:), allocatable ::
     &     ulocal, vv, finv, fvis, ginv, gvis, hinv, hvis
      double precision, dimension(:,:,:,:,:), allocatable ::
     &     dudx, dbdx

#else
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &     :: ulocal, vv, finv
#ifndef NOVISCFLUX
      double precision :: fvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif
#ifndef ONE_D
      double precision :: ginv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#ifndef NOVISCFLUX
      double precision :: gvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif
#ifndef TWO_D
      double precision :: hinv(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#ifndef NOVISCFLUX
      double precision :: hvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
#endif
#endif
#endif

#ifndef NOVISCFLUX
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3)
     &     :: dudx, dbdx
#endif
#endif

c     PETSc distributed array used for parallel diagnostic output
#include <finclude/petscdmdef.h>
      PetscInt ndiagvar
      DM da
c     communicator for DMDA, ranks reordered from DM
      MPI_Comm dacomm3d

      end module petsc_data

      module rmhd3_petsc
      interface TSMonitor_RMHD3
      subroutine TSMonitor_RMHD3(ts,step,time,uvec,dummy,ierr)
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscts.h>
      TS :: ts
      PetscReal :: step,time
      Vec :: uvec
      PetscObject :: dummy
      PetscErrorCode :: ierr
      end subroutine TSMonitor_RMHD3
      end interface TSMonitor_RMHD3
      interface TSRHSFunction_RMHD3
      subroutine TSRHSFunction_RMHD3(ts,time,uvec,fvec,dummy,ierr)
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscts.h>
      TS :: ts
      PetscReal :: time
      Vec :: uvec,fvec
      PetscObject :: dummy
      PetscErrorCode :: ierr
      end subroutine TSRHSFunction_RMHD3
      end interface TSRHSFunction_RMHD3
      interface TSIFunction_RMHD3
      subroutine TSIFunction_RMHD3(ts,time,uvec,udotvec,gvec,dummy,ierr)
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscts.h>
      TS :: ts
      PetscReal :: time
      Vec :: uvec,udotvec,gvec
      PetscObject :: dummy
      PetscErrorCode :: ierr
      end subroutine TSIFunction_RMHD3
      end interface TSIFunction_RMHD3
      end module rmhd3_petsc
c=================================================================




      subroutine TSRHSFunction_RMHD3(ts, t, uvec, fvec, user, ierr)
c-----------------------------------------------------------------------
c     Description: TSRHSFunction_RMHD3 provides the right hand side function for the 
c        ode/dae:  g(t,u,du/dt) = f(u,t)
c        for the current formulation of our pde in conservation form, 
c        we have  du/dt = div(F(u)),  where F is the (nonlinear) flux.
c
c     arguments:
c          ts - time stepping context
c           t - (dbl, input) current time
c        uvec - (vec, input) fluid variables at current time
c        fvec - (vec, output) right-hand side function
c        user - optional user context, unused because we use f90 module
c        ierr - return code
c
c     note: all of the data required to evaluate the system function, is
c        handled through fortran modules. it could also be passed
c        through the user argument.
c-----------------------------------------------------------------------
c======= inclusions ===========
      use mesh
      use profiling
      use petsc_data

c======= declarations =========
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscsnes.h>
#include <finclude/petscts.h>
#include <finclude/petscvec.h90>

c     calling variables
      TS ts
      Vec uvec,fvec
      PetscScalar,pointer :: uarr(:),farr(:)
      PetscObject user
      double precision :: t
      PetscInt         :: i, j, k, l
      PetscErrorCode ierr


c======= Internals ============

c     exchange ghost values and set the boundary conditions
      call VecGetArrayF90(uvec,uarr,ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(fvec,farr,ierr)
      CHKERRQ(ierr)
      call GlobalToLocal(uarr,ulocal)
      call CommState(ulocal)

c     start the RHS profiler
      call prof_start(2)

c     get the negative divergence of the flux at the new state uu
      call divFluxInviscid(ulocal,farr)

c     stop the RHS profiler 
      call prof_stop(2)

      call VecRestoreArrayF90(uvec,uarr,ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayF90(fvec,farr,ierr)
      CHKERRQ(ierr)
      
      return
      end subroutine TSRHSFunction_RMHD3
c=======================================================================

      subroutine TSIFunction_RMHD3(ts,t,uvec,udotvec,gvec,user,ierr)
c-----------------------------------------------------------------------
c     Description: TSIFunction_RMHD3 provides the left hand side function for the
c        ode/dae:  g(t,u,du/dt) = f(u,t)
c        for the current formulation of our PDE in conservation form,
c        we have  du/dt = div(F(u)),  where F is the (nonlinear) flux.
c
c     arguments:
c          ts - time stepping context
c           t - (dbl, input) current time
c        uvec - (vec, input) fluid variables at current time
c     udotvec - (vec, input) time derivative of fluid variables at current time
c        gvec - (vec, output) implicit part of ODE
c        user - optional user context, unused because we use f90 module
c        ierr - return code
c
c     note: all of the data required to evaluate the system function, is
c        handled through fortran modules. it could also be passed
c        through the user argument.
c-----------------------------------------------------------------------
c======= inclusions ===========
      use mesh
      use profiling
      use petsc_data

c======= declarations =========
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscsnes.h>
#include <finclude/petscts.h>
#include <finclude/petscvec.h90>

c     calling variables
      TS ts
      Vec uvec,udotvec,gvec
      PetscScalar,pointer :: uarr(:),garr(:)
      PetscObject user
      double precision :: t
      PetscInt         :: i, j, k, l
      PetscErrorCode ierr


c======= Internals ============

c     exchange ghost values and set the boundary conditions
      call VecGetArrayF90(uvec,uarr,ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(gvec,garr,ierr)
      CHKERRQ(ierr)
      call GlobalToLocal(uarr,ulocal)
      call CommState(ulocal)

c     start the RHS profiler
      call prof_start(2)

c     get the negative divergence of the flux at the new state uu
      call divFluxViscous(ulocal,garr)

c     stop the RHS profiler 
      call prof_stop(2)

      call VecRestoreArrayF90(uvec,uarr,ierr)
      CHKERRQ(ierr)
      call VecRestoreArrayF90(gvec,garr,ierr)
      CHKERRQ(ierr)

      call VecAYPX(gvec,-1.0d0,udotvec,ierr)
      CHKERRQ(ierr)

      return
      end subroutine TSIFunction_RMHD3
c=======================================================================

      subroutine TSMonitor_RMHD3(ts,timeStep,time,uvec,dummy,ierr)
      use petsc_data, only : ulocal, phi, binary_flag, ndiag, ndump,
     &     outnode, output_flag
      use iounits
      use profiling
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscdmda.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>
#include <finclude/petscsnes.h>
#include <finclude/petscts.h>
#include <finclude/petscvec.h90>
      TS ts
      PetscInt timeStep
      PetscObject dummy
      PetscReal time
      Vec      uvec
      PetscErrorCode ierr
      PetscScalar,pointer :: uarr(:)
      PetscReal dt,dthyp,dtvisc
      PetscBool :: output_vts

      call TSGetTimeStep(ts,dt,ierr)
      call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-output_vts',
     &     output_vts,PETSC_NULL_BOOL,ierr)
      CHKERRQ(ierr)

c     periodically dump data
      if(mod(timeStep,ndump)==0
     &     .and. (output_flag > 0 .or. output_vts)) then
c     Update ghost values. They will normally already be updated.
         call VecGetArrayF90(uvec,uarr,ierr)
         call GlobalToLocal(uarr,ulocal)
         call VecRestoreArrayF90(uvec,uarr,ierr)

c     compute explicit time step
         call newdt(ulocal,dthyp,dtvisc)
         if(outnode) then
            write(6,*) 'Hyperbolic/viscous explicit time step: ',
     &           dthyp, dtvisc
         endif

c     Any upstream output
         call OutputData(ulocal,phi,timeStep,
     &           output_flag,binary_flag,min(dthyp,dtvisc),time)
         output_vts = PETSC_FALSE

         call RMHD3PetscWriteVTKFile(ulocal,timeStep,time)
      endif

c     periodically output diagnostic information
      if(mod(timeStep,ndiag)==0) then
         call Diagnostics(ulocal,phi,time)
      endif
      end subroutine TSMonitor_RMHD3
c=======================================================================

      subroutine GlobalToLocal(uglobal,ulocal)
c     Map global state to local state (which includes ghost values)
      use mesh
      use mesh_parms
      use mpistuff
      implicit none
      PetscInt :: i,j,k,l
      double precision :: uglobal(nxlocal,nylocal,nzlocal,nvar)
      double precision :: ulocal(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
c     copy the global state into the local state
      do l=1,nvar
         do k=1,nzlocal
            do j=1,nylocal
               do i=1,nxlocal
                  ulocal(i,j,k,l) = uglobal(i,j,k,l)
               end do
            end do
         end do
      end do
      call CommState(ulocal)
      end subroutine GlobalToLocal

      subroutine LocalToGlobal(ulocal,uglobal)
c     Map local state (includes ghost values) to global state
      use mesh
      use mesh_parms
      use mpistuff
      implicit none
      PetscInt :: i,j,k,l
      double precision :: uglobal(nxlocal,nylocal,nzlocal,nvar)
      double precision :: ulocal(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
c     copy the global state into the local state
      do l=1,nvar
         do k=1,nzlocal
            do j=1,nylocal
               do i=1,nxlocal
                  uglobal(i,j,k,l) = ulocal(i,j,k,l)
               end do
            end do
         end do
      end do
      end subroutine LocalToGlobal

      subroutine CommState(uu)
c-----------------------------------------------------------------------
c     Description: CommState updates the ghost values in-place.
c
c     Argument:
c           uu - (dbl(*), input) current fluid variables
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
#ifdef PARALLEL
      use mpistuff
#endif
      use profiling

c======= Declarations =========
      implicit none
      PetscInt :: i,j,k,l
      double precision :: uu(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)


c======= Internals ============

c     start the communication profiler
      call prof_start(3)

c     do the communications
      call mesh_update_bdry_async(uu,nvar)
      call SetBoundaryValues(uu)

c     wait for all processes to catch up
#ifdef PARALLEL
      call mpi_barrier(comm3d, ierr)
#endif

c     stop the communication profiler
      call prof_stop(3)

      return
      end subroutine CommState
c=================================================================



      subroutine RMHD3PetscInitialize
c-----------------------------------------------------------------------
c     Description: RMHD3PetscInitialize allocates work vectors and initializes PETSc.

c-------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use boundary_conds
      use mpistuff, only: comm3d
      use petsc_data

c======= Declarations =========
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
      double precision :: uu(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: ttot
      double precision :: abstol(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      double precision :: tstop, dtfac, ropt(40)
      integer          :: taskchoice
      integer          :: i,j,k,newrank(0:Nprocs-1)
      PetscErrorCode   :: ierr
      MPI_Group        :: group,newgroup

c======= Internals ============

c     if using dynamic meshing, allocate work space for PETSc function evaluation
#ifdef DYNAMIC
      allocate(   phi(ixlo:ixhi,iylo:iyhi,izlo:izhi))
      allocate(ulocal(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(    vv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(  finv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#ifndef NOVISCFLUX
      allocate(  fvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
      allocate(  dudx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3))
      allocate(  dbdx(ixlo:ixhi,iylo:iyhi,izlo:izhi,3,3))
#endif
#ifndef ONE_D
      allocate(  ginv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#ifndef NOVISCFLUX
      allocate(  gvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#endif
#ifndef TWO_D
      allocate(  hinv(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#ifndef NOVISCFLUX
      allocate(  hvis(ixlo:ixhi,iylo:iyhi,izlo:izhi,1:nvar))
#endif
#endif
#endif
#endif

      call PetscInitialize('petsccontrol.inp',ierr)

      call MPI_Comm_group(comm3d,group,ierr)
      do i=0,Xprocs-1
         do j=0,Yprocs-1
            do k=0,Zprocs-1
               newrank((i*Yprocs+j)*Zprocs+k) = i+Xprocs*(j+Yprocs*k)
            end do
         end do
      end do
      call MPI_Group_incl(group,Nprocs,newrank,newgroup,ierr)
      call MPI_Comm_create(comm3d,newgroup,dacomm3d,ierr)
      call MPI_Group_free(group,ierr)
      call MPI_Group_free(newgroup,ierr)

      ndiagvar = 11
      call DMDACreate3d(dacomm3d,
     &     DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
     &     DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,nx,ny,nz,
     &     Xprocs,Yprocs,Zprocs,
     &     ndiagvar,0,
     &     PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,
     &     da,ierr)
      CHKERRQ(ierr)
      call DMDASetFieldName(da,0,'rho',ierr)
      call DMDASetFieldName(da,1,'Vx',ierr)
      call DMDASetFieldName(da,2,'Vy',ierr)
      call DMDASetFieldName(da,3,'Vz',ierr)
      call DMDASetFieldName(da,4,'Bx',ierr)
      call DMDASetFieldName(da,5,'By',ierr)
      call DMDASetFieldName(da,6,'Bz',ierr)
      call DMDASetFieldName(da,7,'Pressure',ierr)
      call DMDASetFieldName(da,8,'divB',ierr)
      call DMDASetFieldName(da,9,'Jcurrent',ierr)
      call DMDASetFieldName(da,10,'te',ierr)
      call DMDASetUniformCoordinates(da,xl,xr,yl,yr,zl,zr,ierr)
      CHKERRQ(ierr)

      return
      end subroutine RMHD3PetscInitialize
c=================================================================


      subroutine RMHD3PetscFinalize()
c-----------------------------------------------------------------------
c     Description: RMHD3PetscFinalize deallocates all dynamic solver memory.
c-------------------------------------------------------------------
c======= Inclusions ===========
      use petsc_data

c======= Declarations =========
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscdm.h>
      PetscErrorCode ierr

c======= Internals ============

      call DMDestroy(da,ierr)
      CHKERRQ(ierr)

      call MPI_Comm_free(dacomm3d,ierr)

c     free PETSc internal data
      call PetscFinalize(ierr)

c     if using dynamic meshing, deallocate petsc_data structures
#ifdef DYNAMIC
      deallocate(ulocal,vv,finv)
#ifndef NOVISCFLUX
      deallocate(fvis,dudx,dbdx)
#endif
#ifndef ONE_D
      deallocate(ginv)
#ifndef NOVISCFLUX
      deallocate(gvis)
#endif
#ifndef TWO_D
      deallocate(hinv)
#ifndef NOVISCFLUX
      deallocate(hvis)
#endif
#endif
#endif
#endif

      return
      end subroutine RMHD3PetscFinalize
c=================================================================


      subroutine divFluxBoth(uu,dFlux)
c-----------------------------------------------------------------------
c     Description: divFlux provides a centered-difference spatial 
c        approximation to the negative divergence of the flux for the 
c        conservation equation describing the nonlinearly-resistive,
c        ideal MHD system.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at current time
c        dFlux - (dbl(*), output) negative divergence of the flux
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use petsc_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar) 
     &     :: uu
      double precision, dimension(nxlocal,nylocal,nzlocal,nvar)
     &     :: dFlux

c     differencing factors
      double precision :: xfac, yfac, zfac, dt

c     indexes for flux calculations
      integer :: indx, indy, indz, i, j, k, l, m
      


c======= Internals ============


c     set differencing factors
      xfac = 1.0d0/dx;  yfac = 1.0d0/dy;  zfac = 1.0d0/dz

c     set dt value (unused, but passed to flux routines anyway)
      dt = 0.0d0

c     construct the vector of primitive variables
      call ConservedToPrimitive(uu,vv)


c     zero the fluxes (necessary???)
      finv = 0.0d0
#ifndef ONE_D
      ginv = 0.0d0
#ifndef TWO_D
      hinv = 0.0d0
#endif
#endif
#ifndef NOVISCFLUX
      fvis = 0.0d0
#ifndef ONE_D
      gvis = 0.0d0
#ifndef TWO_D
      hvis = 0.0d0
#endif
#endif
#endif


c     get the viscous flux requirements dudx and dbdx
#ifndef NOVISCFLUX
      do m=1,3,1
         dudx(:,:,:,1,m) = vv(:,:,:,m+1)
         dbdx(:,:,:,1,m) = vv(:,:,:,m+4)
         dudx(:,:,:,2,m) = 0.d0
         dbdx(:,:,:,2,m) = 0.d0
         dudx(:,:,:,3,m) = 0.d0
         dbdx(:,:,:,3,m) = 0.d0
      end do
      call GradientVector(dudx)
      call GradientVector(dbdx)
#endif


c------------ x-directional fluxes ------------
c     Hyperbolic fluxes
      inlo = ixlo;  inhi = ixhi
#ifdef UPWIND
#ifdef PRIMITIVE
c        using upwind with eigensystem from primitive variables
      call InviscidFluxRP(finv,vv,1,dt)
#elif  ROE
c        using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(finv,vv,1,dt)
#else
c        using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(finv,vv,1,dt)
#endif
#else
c        using centered-differences (of some form)
#ifdef TCD
c        2nd-order tuned central differences
      call InviscidFluxTCD(finv,vv,1,dt)
#elif CD4
c        4th-order standard central differences
      call InviscidFluxCD4(finv,vv,1,dt)
#else
c        2nd-order standard central differences
      call InviscidFluxCD2(finv,vv,1,dt)
#endif
#endif

c     Parabolic fluxes
      inlo = ixlo;  inhi = ixhi
      indx = 1;  indy = 0;  indz = 0
#ifndef NOVISCFLUX
#ifdef TCD
c        2nd-order tuned central differences
      call ViscousFluxTCD(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)
#elif CD4
c        4th-order standard central differences
      call ViscousFluxCD4(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)
#else
c        2nd-order standard central differences
      call ViscousFluxCD2(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)
#endif
#endif


c------------ y-directional fluxes ------------
#ifndef ONE_D
c     Hyperbolic fluxes
      inlo = iylo;  inhi = iyhi
#ifdef UPWIND
#ifdef PRIMITIVE
c        using upwind with eigensystem from primitive variables
      call InviscidFluxRP(ginv,vv,2,dt)
#elif  ROE
c        using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(ginv,vv,2,dt)
#else
c        using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(ginv,vv,2,dt)
#endif
#else
c        using centered-differences (of some form)
#ifdef TCD
c        2nd-order tuned central differences
      call InviscidFluxTCD(ginv,vv,2,dt)
#elif CD4
c        4th-order standard central differences
      call InviscidFluxCD4(ginv,vv,2,dt)
#else
c        2nd-order standard central differences
      call InviscidFluxCD2(ginv,vv,2,dt)
#endif
#endif

c     Parabolic fluxes
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
#ifndef NOVISCFLUX
#ifdef TCD
c        2nd-order tuned central differences
      call ViscousFluxTCD(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)
#elif CD4
c        4th-order standard central differences
      call ViscousFluxCD4(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)
#else
c        2nd-order standard central differences
      call ViscousFluxCD2(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)
#endif
#endif


c------------ z-directional fluxes ------------
#ifndef TWO_D
c     Hyperbolic fluxes
      inlo = izlo;  inhi = izhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
      call InviscidFluxRP(hinv,vv,3,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(hinv,vv,3,dt)
#else
c        using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(hinv,vv,3,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
c        2nd-order tuned central differences
      call InviscidFluxTCD(hinv,vv,3,dt)
#elif CD4
c        4th-order standard central differences
      call InviscidFluxCD4(hinv,vv,3,dt)
#else
c        2nd-order standard central differences
      call InviscidFluxCD2(hinv,vv,3,dt)
#endif
#endif

c     Parabolic fluxes
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
#ifndef NOVISCFLUX
#ifdef TCD
c        2nd-order tuned central differences
      call ViscousFluxTCD(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#elif CD4
c        4th-order standard central differences
      call ViscousFluxCD4(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#else
c        2nd-order standard central differences
      call ViscousFluxCD2(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
#endif
#endif


c------------ construction of div(flux) ------------
c
c     Density
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,1) = -xfac*(finv(i,j,k,1) - finv(i-1,j,k,1))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,1) - ginv(i,j-1,k,1) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,1) - hinv(i,j,k-1,1) )
#endif
#endif
            end do
         end do
      end do
c
c     X-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,2) = -xfac*(finv(i,j,k,2) - finv(i-1,j,k,2))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,2) - fvis(i-1,j,k,2) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,4) - ginv(i,j-1,k,4) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,4) - gvis(i,j-1,k,4) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,3) - hinv(i,j,k-1,3) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,3) - hvis(i,j,k-1,3) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Y-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,3) = -xfac*(finv(i,j,k,3) - finv(i-1,j,k,3))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,3) - fvis(i-1,j,k,3) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,2) - ginv(i,j-1,k,2) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,2) - gvis(i,j-1,k,2) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,4) - hinv(i,j,k-1,4) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,4) - hvis(i,j,k-1,4) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Z-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,4) = -xfac*(finv(i,j,k,4) - finv(i-1,j,k,4))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,4) - fvis(i-1,j,k,4) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,3) - ginv(i,j-1,k,3) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,3) - gvis(i,j-1,k,3) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,2) - hinv(i,j,k-1,2) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,2) - hvis(i,j,k-1,2) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     X-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,5) = -xfac*(finv(i,j,k,5) - finv(i-1,j,k,5))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,5) - fvis(i-1,j,k,5) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,7) - ginv(i,j-1,k,7) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,7) - gvis(i,j-1,k,7) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,6) - hinv(i,j,k-1,6) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,6) - hvis(i,j,k-1,6) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Y-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,6) = -xfac*(finv(i,j,k,6) - finv(i-1,j,k,6))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,6) - fvis(i-1,j,k,6) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,5) - ginv(i,j-1,k,5) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,5) - gvis(i,j-1,k,5) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,7) - hinv(i,j,k-1,7) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,7) - hvis(i,j,k-1,7) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Z-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,7) = -xfac*(finv(i,j,k,7) - finv(i-1,j,k,7))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,7) - fvis(i-1,j,k,7) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,6) - ginv(i,j-1,k,6) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,6) - gvis(i,j-1,k,6) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,5) - hinv(i,j,k-1,5) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,5) - hvis(i,j,k-1,5) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Energy
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,8) = -xfac*(finv(i,j,k,8) - finv(i-1,j,k,8))
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,8) - fvis(i-1,j,k,8) )
#endif
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,8) - ginv(i,j-1,k,8) )
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,8) - gvis(i,j-1,k,8) )
#endif
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,8) - hinv(i,j,k-1,8) )
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,8) - hvis(i,j,k-1,8) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c


      return
      end subroutine divFluxBoth
c=======================================================================

      subroutine divFluxInviscid(uu,dFlux)
c-----------------------------------------------------------------------
c     Description: divFlux provides a centered-difference spatial
c        approximation to the negative divergence of the flux for the
c        conservation equation describing the nonlinearly-resistive,
c        ideal MHD system.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at current time
c        dFlux - (dbl(*), output) negative divergence of the inviscid flux
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use petsc_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &     :: uu
      double precision, dimension(nxlocal,nylocal,nzlocal,nvar)
     &     :: dFlux

c     differencing factors
      double precision :: xfac, yfac, zfac, dt

c     indexes for flux calculations
      integer :: indx, indy, indz, i, j, k, l, m


c======= Internals ============


c     set differencing factors
      xfac = 1.0d0/dx;  yfac = 1.0d0/dy;  zfac = 1.0d0/dz

c     set dt value (unused, but passed to flux routines anyway)
      dt = 0.0d0

c     construct the vector of primitive variables
      call ConservedToPrimitive(uu,vv)


c     zero the fluxes (necessary???)
      finv = 0.0d0
#ifndef ONE_D
      ginv = 0.0d0
#ifndef TWO_D
      hinv = 0.0d0
#endif
#endif

c------------ x-directional fluxes ------------
c     Hyperbolic fluxes
      inlo = ixlo;  inhi = ixhi
#ifdef UPWIND
#ifdef PRIMITIVE
c        using upwind with eigensystem from primitive variables
      call InviscidFluxRP(finv,vv,1,dt)
#elif  ROE
c        using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(finv,vv,1,dt)
#else
c        using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(finv,vv,1,dt)
#endif
#else
c        using centered-differences (of some form)
#ifdef TCD
c        2nd-order tuned central differences
      call InviscidFluxTCD(finv,vv,1,dt)
#elif CD4
c        4th-order standard central differences
      call InviscidFluxCD4(finv,vv,1,dt)
#else
c        2nd-order standard central differences
      call InviscidFluxCD2(finv,vv,1,dt)
#endif
#endif

c------------ y-directional fluxes ------------
#ifndef ONE_D
c     Hyperbolic fluxes
      inlo = iylo;  inhi = iyhi
#ifdef UPWIND
#ifdef PRIMITIVE
c        using upwind with eigensystem from primitive variables
      call InviscidFluxRP(ginv,vv,2,dt)
#elif  ROE
c        using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(ginv,vv,2,dt)
#else
c        using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(ginv,vv,2,dt)
#endif
#else
c        using centered-differences (of some form)
#ifdef TCD
c        2nd-order tuned central differences
      call InviscidFluxTCD(ginv,vv,2,dt)
#elif CD4
c        4th-order standard central differences
      call InviscidFluxCD4(ginv,vv,2,dt)
#else
c        2nd-order standard central differences
      call InviscidFluxCD2(ginv,vv,2,dt)
#endif
#endif

c------------ z-directional fluxes ------------
#ifndef TWO_D
c     Hyperbolic fluxes
      inlo = izlo;  inhi = izhi
#ifdef UPWIND
#ifdef PRIMITIVE
c     using upwind with eigensystem from primitive variables
      call InviscidFluxRP(hinv,vv,3,dt)
#elif  ROE
c     using upwind with eigensystem from conservative variables
      call InviscidFluxRoe(hinv,vv,3,dt)
#else
c        using upwind with Lax-Friedrichs scheme
      call InviscidFluxLF(hinv,vv,3,dt)
#endif
#else
c     using centered-differences (of some form)
#ifdef TCD
c        2nd-order tuned central differences
      call InviscidFluxTCD(hinv,vv,3,dt)
#elif CD4
c        4th-order standard central differences
      call InviscidFluxCD4(hinv,vv,3,dt)
#else
c        2nd-order standard central differences
      call InviscidFluxCD2(hinv,vv,3,dt)
#endif
#endif

c------------ construction of div(flux) ------------
c
c     Density
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,1) = -xfac*(finv(i,j,k,1) - finv(i-1,j,k,1))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,1) - ginv(i,j-1,k,1) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,1) - hinv(i,j,k-1,1) )
#endif
#endif
            end do
         end do
      end do
c
c     X-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,2) = -xfac*(finv(i,j,k,2) - finv(i-1,j,k,2))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,4) - ginv(i,j-1,k,4) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,3) - hinv(i,j,k-1,3) )
#endif
#endif
            end do
         end do
      end do
c
c     Y-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,3) = -xfac*(finv(i,j,k,3) - finv(i-1,j,k,3))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,2) - ginv(i,j-1,k,2) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,4) - hinv(i,j,k-1,4) )
#endif
#endif
            end do
         end do
      end do
c
c     Z-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,4) = -xfac*(finv(i,j,k,4) - finv(i-1,j,k,4))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,3) - ginv(i,j-1,k,3) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,2) - hinv(i,j,k-1,2) )
#endif
#endif
            end do
         end do
      end do
c
c     X-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,5) = -xfac*(finv(i,j,k,5) - finv(i-1,j,k,5))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,7) - ginv(i,j-1,k,7) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,6) - hinv(i,j,k-1,6) )
#endif
#endif
            end do
         end do
      end do
c
c     Y-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,6) = -xfac*(finv(i,j,k,6) - finv(i-1,j,k,6))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,5) - ginv(i,j-1,k,5) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,7) - hinv(i,j,k-1,7) )
#endif
#endif
            end do
         end do
      end do
c
c     Z-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,7) = -xfac*(finv(i,j,k,7) - finv(i-1,j,k,7))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,6) - ginv(i,j-1,k,6) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,5) - hinv(i,j,k-1,5) )
#endif
#endif
            end do
         end do
      end do
c
c     Energy
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,8) = -xfac*(finv(i,j,k,8) - finv(i-1,j,k,8))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,8) - ginv(i,j-1,k,8) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,8) - hinv(i,j,k-1,8) )
#endif
#endif
            end do
         end do
      end do
c
c


      return
      end subroutine divFluxInviscid
c=======================================================================

      subroutine divFluxViscous(uu,dFlux)
c-----------------------------------------------------------------------
c     Description: divFlux provides a centered-difference spatial
c        approximation to the negative divergence of the flux for the
c        conservation equation describing the nonlinearly-resistive,
c        ideal MHD system.
c
c     Arguments:
c           uu - (dbl(*), input) fluid variables at current time
c        dFlux - (dbl(*), output) negative divergence of the viscous flux
c-----------------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use properties
      use petsc_data

c======= Declarations =========
      implicit none

c     calling variables
      double precision, dimension(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
     &     :: uu
      double precision, dimension(nxlocal,nylocal,nzlocal,nvar)
     &     :: dFlux

c     differencing factors
      double precision :: xfac, yfac, zfac, dt

c     indexes for flux calculations
      integer :: indx, indy, indz, i, j, k, l, m


c======= Internals ============


c     set differencing factors
      xfac = 1.0d0/dx;  yfac = 1.0d0/dy;  zfac = 1.0d0/dz

c     set dt value (unused, but passed to flux routines anyway)
      dt = 0.0d0

c     construct the vector of primitive variables
      call ConservedToPrimitive(uu,vv)


c     zero the fluxes (necessary???)
#ifndef NOVISCFLUX
      fvis = 0.0d0
#ifndef ONE_D
      gvis = 0.0d0
#ifndef TWO_D
      hvis = 0.0d0
#endif
#endif
#endif


c     get the viscous flux requirements dudx and dbdx
#ifndef NOVISCFLUX
      do m=1,3,1
         dudx(:,:,:,1,m) = vv(:,:,:,m+1)
         dbdx(:,:,:,1,m) = vv(:,:,:,m+4)
         dudx(:,:,:,2,m) = 0.d0
         dbdx(:,:,:,2,m) = 0.d0
         dudx(:,:,:,3,m) = 0.d0
         dbdx(:,:,:,3,m) = 0.d0
      end do
      call GradientVector(dudx)
      call GradientVector(dbdx)
#endif


c------------ x-directional fluxes ------------
c     Parabolic fluxes
      inlo = ixlo;  inhi = ixhi
      indx = 1;  indy = 0;  indz = 0
#ifndef NOVISCFLUX
#ifdef TCD
c        2nd-order tuned central differences
      call ViscousFluxTCD(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)
#elif CD4
c        4th-order standard central differences
      call ViscousFluxCD4(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)
#else
c        2nd-order standard central differences
      call ViscousFluxCD2(fvis,vv,dudx,dbdx,1,dt,indx,indy,indz)
#endif
#endif


c------------ y-directional fluxes ------------
c     Parabolic fluxes
      inlo = iylo;  inhi = iyhi
      indx = 0;  indy = 1;  indz = 0
#ifndef NOVISCFLUX
#ifdef TCD
c        2nd-order tuned central differences
      call ViscousFluxTCD(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)
#elif CD4
c        4th-order standard central differences
      call ViscousFluxCD4(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)
#else
c        2nd-order standard central differences
      call ViscousFluxCD2(gvis,vv,dudx,dbdx,2,dt,indx,indy,indz)
#endif
#endif


c------------ z-directional fluxes ------------
c     Parabolic fluxes
      inlo = izlo;  inhi = izhi
      indx = 0;  indy = 0;  indz = 1
#ifndef NOVISCFLUX
#ifdef TCD
c        2nd-order tuned central differences
      call ViscousFluxTCD(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#elif CD4
c        4th-order standard central differences
      call ViscousFluxCD4(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#else
c        2nd-order standard central differences
      call ViscousFluxCD2(hvis,vv,dudx,dbdx,3,dt,indx,indy,indz)
#endif
#endif
#endif
#endif


c------------ construction of div(flux) ------------
c
c     Density
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,1) = -xfac*(finv(i,j,k,1) - finv(i-1,j,k,1))
#ifndef ONE_D
     &              -yfac*( ginv(i,j,k,1) - ginv(i,j-1,k,1) )
#ifndef TWO_D
     &              -zfac*( hinv(i,j,k,1) - hinv(i,j,k-1,1) )
#endif
#endif
            end do
         end do
      end do
c
c     X-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,2) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,2) - fvis(i-1,j,k,2) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,4) - gvis(i,j-1,k,4) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,3) - hvis(i,j,k-1,3) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Y-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,3) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,3) - fvis(i-1,j,k,3) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,2) - gvis(i,j-1,k,2) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,4) - hvis(i,j,k-1,4) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Z-Momentum
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,4) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,4) - fvis(i-1,j,k,4) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,3) - gvis(i,j-1,k,3) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,2) - hvis(i,j,k-1,2) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     X-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,5) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,5) - fvis(i-1,j,k,5) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,7) - gvis(i,j-1,k,7) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,6) - hvis(i,j,k-1,6) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Y-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,6) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,6) - fvis(i-1,j,k,6) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,5) - gvis(i,j-1,k,5) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,7) - hvis(i,j,k-1,7) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Z-Magnetic field
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,7) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,7) - fvis(i-1,j,k,7) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,6) - gvis(i,j-1,k,6) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,5) - hvis(i,j,k-1,5) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c     Energy
      do k=1,nzlocal,1
         do j=1,nylocal,1
            do i=1,nxlocal,1
               dFlux(i,j,k,8) =
#ifndef NOVISCFLUX
     &              +xfac*( fvis(i,j,k,8) - fvis(i-1,j,k,8) )
#endif
#ifndef ONE_D
#ifndef NOVISCFLUX
     &              +yfac*( gvis(i,j,k,8) - gvis(i,j-1,k,8) )
#endif
#ifndef TWO_D
#ifndef NOVISCFLUX
     &              +zfac*( hvis(i,j,k,8) - hvis(i,j,k-1,8) )
#endif
#endif
#endif
            end do
         end do
      end do
c
c


      return
      end subroutine divFluxViscous
c=======================================================================

      subroutine RMHD3PetscWriteVTKFile(ulocal,timeStep,time)
      use mesh
      use petsc_data, only: da, dacomm3d, outnode
      implicit none
#include <finclude/petscsys.h>
#include <finclude/petscviewer.h>
#include <finclude/petscvec.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscvec.h90>
      double precision     :: ulocal(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      PetscInt             :: timeStep
      PetscReal            :: time
      PetscViewer          :: viewer
      Vec                  :: diagvec
      PetscScalar, pointer :: darr(:)
      PetscErrorCode       :: ierr
      character*256        :: outfile

      call DMGetGlobalVector(da,diagvec,ierr)
      CHKERRQ(ierr)
      call VecGetArrayF90(diagvec,darr,ierr)
      CHKERRQ(ierr)
      call RMHD3PetscComputeDiagonstics(ulocal,darr)
      call VecRestoreArrayF90(diagvec,darr,ierr)
      CHKERRQ(ierr)

      write(outfile,1000) timeStep
 1000 format(7Hisland-,i4.4,4H.vts)
      if (outnode) then
         write(6,*) 'Writing ', outfile
      end if

      call PetscViewerVTKOpen(dacomm3d,outfile,FILE_MODE_WRITE,
     &     viewer,ierr)
      CHKERRQ(ierr)
      call VecView(diagvec,viewer,ierr)
      CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr)
      CHKERRQ(ierr)

      call DMRestoreGlobalVector(da,diagvec,ierr)
      CHKERRQ(ierr)
      end subroutine RMHD3PetscWriteVTKFile

      subroutine RMHD3PetscComputeDiagonstics(ux,darr)
      use mesh_parms
      use petsc_data, only: ndiagvar
      use properties
      implicit none
      integer :: i,j,k
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      PetscScalar ::
     &     darr(0:ndiagvar-1,nxlocal,nylocal,nzlocal)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: jcurrent(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: rho,u,v,w,bi,bj,bk,press

c     get div B, toroidal current
      call DivergenceCD(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),div)
      call ToroidalCurrent(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),jcurrent)

      do k=1,nzlocal
         do j=1,nylocal
            do i=1,nxlocal
               rho = ux(i,j,k,1)
               u   = ux(i,j,k,2)/rho
               v   = ux(i,j,k,3)/rho
               w   = ux(i,j,k,4)/rho
               bi  = ux(i,j,k,5)
               bj  = ux(i,j,k,6)
               bk  = ux(i,j,k,7)
               press = (gamma-1.0D0)*(ux(i,j,k,8)
     &              -0.5D0*rho*(u*u+v*v+w*w)-0.5D0*(bi*bi+bj*bj+bk*bk))
               darr(0,i,j,k) = rho
               darr(1,i,j,k) = u
               darr(2,i,j,k) = v
               darr(3,i,j,k) = w
               darr(4,i,j,k) = bi
               darr(5,i,j,k) = bj
               darr(6,i,j,k) = bk
               darr(7,i,j,k) = press
               darr(8,i,j,k) = div(i,j,k)
               darr(9,i,j,k) = jcurrent(i,j,k)
               darr(10,i,j,k) = press+0.5D0*(bi*bi+bj*bj+bk*bk)
            end do
         end do
      end do
      end subroutine RMHD3PetscComputeDiagonstics
