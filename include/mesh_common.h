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
c     $Log: mesh_common.h,v $
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
