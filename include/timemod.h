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

c      module timemod
c-----------------------------------------------------------------
c     Description: holds time-stepping information (disabled).  
c        Now contains the global common block 'time'
c
c     Contains:
c          dt - (dbl) the time-step
c        ttot - (dbl) the current time
c-----------------------------------------------------------------

      double precision:: dt, ttot
      common /time/ dt, ttot

c      end module timemod
c=================================================================
