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
c=================================================================


      module mesh_uparms
c-----------------------------------------------------------------
c     Description: universal mesh and problem parameters for 
c        static, single-processor grid. 
c
c     Contains:
c        Xprocs - (int, param) number of processors in x-direction
c        Yprocs - (int, param) number of processors in y-direction
c        Xprocs - (int, param) number of processors in z-direction
c            nx - (int, param) number of x-meshpoints
c            ny - (int, param) number of y-meshpoints
c            nz - (int, param) number of z-meshpoints
c        nghost - (int, param) number of ghost cells at boundaries
c          nvar - (int, param) number of fluid state variables
c-----------------------------------------------------------------
      save

      integer, parameter :: XPROCS=1
      integer, parameter :: YPROCS=1
      integer, parameter :: ZPROCS=1
      integer, parameter :: nx=64
      integer, parameter :: ny=64
      integer, parameter :: nz=64
      integer, parameter :: nghost=4
      integer, parameter :: nvar=8
      
      end module mesh_uparms
c=================================================================
