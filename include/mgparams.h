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

      module MGParams
c-----------------------------------------------------------------
c     Description: holds all of the general multigrid parameters
c
c     Contains:
c               lMax - (int) max number of multigrid levels
c          lMaxlocal - (int) max number of local multigrid levels
c            bcHFlag - (int) homogeneous BC at this level
c        xChangeFlag - (int) exchange ghost info at this level
c-----------------------------------------------------------------
      save

      integer :: lMax, lMaxLocal
      integer :: bcHFlag
      integer :: xChangeFlag

      end module MGParams
c=================================================================
