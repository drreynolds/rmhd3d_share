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

      module mesh_uparms
c-----------------------------------------------------------------
c     Description: universal mesh and problem parameters
c
c     Contains:
c        Xprocs - (int) number of processors in x-direction
c        Yprocs - (int) number of processors in y-direction
c        Zprocs - (int) number of processors in z-direction
c            nx - (int) number of x-meshpoints
c            ny - (int) number of y-meshpoints
c            nz - (int) number of z-meshpoints
c        nghost - (int, param) number of ghost cells at boundaries
c          nvar - (int, param) number of fluid state variables
c-----------------------------------------------------------------
      save

      integer :: Xprocs, Yprocs, Zprocs
      integer :: debugprocs
      integer :: nx, ny, nz, nghost
      integer, parameter :: nvar=8

      end module mesh_uparms
c=================================================================
