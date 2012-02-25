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

      module mesh_common
c-----------------------------------------------------------------
c     Description: holds all of the general mesh location 
c        information for a given processor 
c
c     Contains:
c           iprocx - (int) x-location in process grid
c           iprocy - (int) y-location in process grid
c           iprocz - (int) z-location in process grid
c        iproc_idx - (int) the global MPI process number
c-----------------------------------------------------------------
      save

      integer :: iprocx, iprocy, iprocz, iproc_idx

      end module mesh_common
c=================================================================
