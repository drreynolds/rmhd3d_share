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
c     $Log: mesh_parms.h,v $
c=================================================================


       module mesh_parms
c-----------------------------------------------------------------
c     Description: holds all of the general mesh information for
c        a given processor 
c
c     Contains:
c         Nprocs - (int) number of processors
c        NXLsize - (int) x-size of local grid (no ghost cells)
c        NYLsize - (int) y-size of local grid (no ghost cells)
c        NZLsize - (int) z-size of local grid (no ghost cells)
c           IXLO - (int) x-lower bound for local grid (w/ ghosts)
c           IXHI - (int) x-upper bound for local grid (w/ ghosts)
c           IYLO - (int) y-lower bound for local grid (w/ ghosts)
c           IYHI - (int) y-upper bound for local grid (w/ ghosts)
c           IZLO - (int) z-lower bound for local grid (w/ ghosts)
c           IZHI - (int) z-upper bound for local grid (w/ ghosts)
c           INLO - (int) interface lower bound for local grid
c           INHI - (int) interface upper bound for local grid
c        NXlocal - (int) x-size of local grid
c        NYlocal - (int) y-size of local grid
c        NZlocal - (int) z-size of local grid
c        XLloc - (dbl) x-direction lower endpoint (local grid)
c        XRloc - (dbl) x-direction upper endpoint (local grid)
c        YLloc - (dbl) y-direction lower endpoint (local grid)
c        YRloc - (dbl) y-direction upper endpoint (local grid)
c        ZLloc - (dbl) z-direction lower endpoint (local grid)
c        ZRloc - (dbl) z-direction upper endpoint (local grid)
c        MSG_XCH_XLOW_TAG - (int, param) dir. of x-low neighbor
c        MSG_XCH_XHI_TAG  - (int, param) dir. of x-high neighbor
c        MSG_XCH_YLOW_TAG - (int, param) dir. of y-low neighbor
c        MSG_XCH_YHI_TAG  - (int, param) dir. of y-high neighbor
c        MSG_XCH_ZLOW_TAG - (int, param) dir. of z-low neighbor
c        MSG_XCH_ZHI_TAG  - (int, param) dir. of z-high neighbor
c        MSG_MAX_TAG      - (int, param) ???? (unused) ????
c-----------------------------------------------------------------
       use mesh_uparms
       save

       integer :: Nprocs
       integer :: NXlsize, NYlsize, NZlsize
       integer :: IXLO, IXHI, IYLO, IYHI, IZLO, IZHI, INLO, INHI
       integer :: NXlocal, NYlocal, NZlocal
       double precision   :: XLloc, XRloc, YLloc, YRloc, ZLloc, ZRloc
       integer, parameter :: MSG_XCH_XLOW_TAG=1, MSG_XCH_XHI_TAG=2
       integer, parameter :: MSG_XCH_YLOW_TAG=3, MSG_XCH_YHI_TAG=4
       integer, parameter :: MSG_XCH_ZLOW_TAG=5, MSG_XCH_ZHI_TAG=6
       integer, parameter :: MSG_MAX_TAG=10

       end module mesh_parms
c=================================================================
