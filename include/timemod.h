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
c     $Log: timemod.h,v $
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
