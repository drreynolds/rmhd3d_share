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
c     $Log: mgparams.h,v $
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
