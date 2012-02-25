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

      module mpistuff
c-----------------------------------------------------------------
c     Description: holds all of the MPI internal information 
c        necessary for the communicator and error messages.
c
c     Contains:
c         comm3d - (int) mpi communicator
c         master - (int) master node information
c         status - (int vec(MPI_STATUS_SIZE)) MPI_status
c           ierr - (int) error flag
c          my_id - (int) local process identifier
c           left - (int) x-lower neighbor id
c          right - (int) x-upper neighbor id
c            top - (int) y-upper neighbor id
c         bottom - (int) y-lower neighbor id
c         behind - (int) z-lower neighbor id
c        forward - (int) z-upper neighbor id
c        ERROR_CARTCOORDS - (int, param) error flag for cartcoords
c         ERROR_CARTSHIFT - (int, param) error flag for cartshift
c              ERROR_WAIT - (int, param) error flag for wait
c              ERROR_SEND - (int, param) error flag for send
c              ERROR_RECV - (int, param) error flag for recv
c         ERROR_ALLREDUCE - (int, param) error flag for allreduce
c-----------------------------------------------------------------
      use mpi

      save
c      include "mpif.h"
      
      integer :: comm3d, master
      integer :: status(MPI_STATUS_SIZE), ierr
      integer :: my_id
      integer :: left, right, top, bottom, behind, forward
      integer, parameter :: ERROR_CARTCOORDS=1
      integer, parameter :: ERROR_CARTSHIFT=2
      integer, parameter :: ERROR_WAIT=3
      integer, parameter :: ERROR_SEND=4
      integer, parameter :: ERROR_RECV=5
      integer, parameter :: ERROR_ALLREDUCE=6
      
      end module mpistuff
c=================================================================
