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

      module mesh
c-----------------------------------------------------------------
c     Description: holds all of the general mesh information for 
c        entire domain
c
c     Contains:
c        xl - (dbl) x-direction lower endpoint
c        xr - (dbl) x-direction upper endpoint
c        yl - (dbl) y-direction lower endpoint
c        yr - (dbl) y-direction upper endpoint
c        zl - (dbl) z-direction lower endpoint
c        zr - (dbl) z-direction upper endpoint
c        dx - (dbl) x-direction mesh increment
c        dy - (dbl) y-direction mesh increment
c        dz - (dbl) z-direction mesh increment
c        xc - (dbl, allocatable vector) x-direction meshpoints
c        yc - (dbl, allocatable vector) y-direction meshpoints
c        zc - (dbl, allocatable vector) z-direction meshpoints
c-----------------------------------------------------------------
       use mesh_parms
       save

       double precision:: xl, xr, yl, yr, zl, zr
       double precision:: dx, dy,dz
       double precision, allocatable, dimension(:):: xc,yc,zc

       end module mesh
c=================================================================
