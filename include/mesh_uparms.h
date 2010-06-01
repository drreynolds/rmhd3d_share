c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c Original code:
c     Ravi Samtaney
c     Copyright 2003
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c Revised for additional diagnositics, implicit time integration, 
c and clarity:
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------
c     $Log: mesh_uparms.h,v $
c=================================================================


      module mesh_uparms
c-----------------------------------------------------------------
c     Description: universal mesh and problem parameters
c
c     Contains:
c        Xprocs - (int) number of processors in x-direction
c        Yprocs - (int) number of processors in y-direction
c        Xprocs - (int) number of processors in z-direction
c            nx - (int) number of x-meshpoints
c            ny - (int) number of y-meshpoints
c            nz - (int) number of z-meshpoints
c        nghost - (int, param) number of ghost cells at boundaries
c          nvar - (int, param) number of fluid state variables
c-----------------------------------------------------------------
      save

      integer :: Xprocs, Yprocs, Zprocs
      integer :: nx, ny, nz, nghost
      integer, parameter :: nvar=8

      end module mesh_uparms
c=================================================================
