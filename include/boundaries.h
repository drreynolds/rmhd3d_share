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

      module boundary_conds
c-----------------------------------------------------------------
c     Description: holds all of the boundary condition choices:
c            0 => zero-gradient
c            1 => periodic
c            2 => reflecting
c
c     Contains:
c        xbc - (int) x-direction boundary condition choice
c        ybc - (int) x-direction boundary condition choice
c        zbc - (int) y-direction boundary condition choice
c-----------------------------------------------------------------
       save

       integer :: xbc
       integer :: ybc
       integer :: zbc
       integer, parameter :: BCperiodic=1, BCreflecting=2, BCzerograd=0

       end module boundary_conds
c=================================================================
