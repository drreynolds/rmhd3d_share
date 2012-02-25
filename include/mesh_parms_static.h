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

      module mesh_parms
c-----------------------------------------------------------------
c     Description: holds all of the general mesh information for
c        a given processor.  In this case, they are set for a 
c        static grid and left unchanged.
c
c     Contains:
c         Nprocs - (int, param) number of processors
c        NXLsize - (int, param) x-size of local grid (no ghosts)
c        NYLsize - (int, param) y-size of local grid (no ghosts)
c        NZLsize - (int, param) z-size of local grid (no ghosts)
c           IXLO - (int, param) x-lower bound for local grid (w/ ghosts)
c           IXHI - (int, param) x-upper bound for local grid (w/ ghosts)
c           IYLO - (int, param) y-lower bound for local grid (w/ ghosts)
c           IYHI - (int, param) y-upper bound for local grid (w/ ghosts)
c           IZLO - (int, param) z-lower bound for local grid (w/ ghosts)
c           IZHI - (int, param) z-upper bound for local grid (w/ ghosts)
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

      integer, parameter :: Nprocs = Xprocs*Yprocs*Zprocs
      integer, parameter :: NXlsize = NX/XPROCS
      integer, parameter :: NYlsize = NY/YPROCS
      integer, parameter :: NZlsize = NZ/ZPROCS
      integer, parameter :: IXLO = 1-NGHOST
      integer, parameter :: IXHI = nxlsize+NGHOST
c$$$      integer, parameter :: IYLO = 1
c$$$      integer, parameter :: IYHI = nylsize
      integer, parameter :: IYLO = 1-NGHOST
      integer, parameter :: IYHI = nylsize+NGHOST
c$$$      integer, parameter :: IZLO = 1
c$$$      integer, parameter :: IZHI = nzlsize
      integer, parameter :: IZLO = 1-NGHOST
      integer, parameter :: IZHI = nzlsize+NGHOST
      integer :: INLO, INHI
      integer :: NXlocal, NYlocal, NZlocal
      double precision   :: XLloc, XRloc, YLloc, YRloc, ZLloc, ZRloc
      integer, parameter :: MSG_XCH_XLOW_TAG=1, MSG_XCH_XHI_TAG=2
      integer, parameter :: MSG_XCH_YLOW_TAG=3, MSG_XCH_YHI_TAG=4
      integer, parameter :: MSG_XCH_ZLOW_TAG=5, MSG_XCH_ZHI_TAG=6
      integer, parameter :: MSG_MAX_TAG=10

      end module mesh_parms
c=================================================================
