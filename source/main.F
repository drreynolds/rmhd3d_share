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

      program MHD
c-----------------------------------------------------------------
c     Description: Code to solve 2D, 2.5D, 3D ideal, 
c        nonlinearly-resistive, MHD equations 
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh_parms
      use mesh_common
      use boundary_conds
      use iounits
#ifdef PARALLEL
      use mpistuff
#endif
      
c======= Declarations =========
      implicit none
      logical :: error
      integer :: num_nodes
      integer :: ipxnum
      integer :: ipynum
      integer :: ipznum
      integer :: ipx, ipy, ipz
      integer :: myid
      integer :: meshbuffer(10)

c     safety checks
      integer :: ndims
      logical :: periods(3), reorder
      integer :: NProcAxis(3)

c     namelists
      namelist /bdryconds/ xbc, ybc, zbc
#ifdef DYNAMIC
      namelist /meshsize/  nx, ny, nz, nghost
      namelist /meshprocs/ xprocs, yprocs, zprocs
#endif

c======= Internals ============

c     initialize the MPI communicator and get initial processor rank
      myid = 0
#ifdef PARALLEL
      call MPI_init(ierr)
      call MPI_comm_rank(MPI_COMM_WORLD, myid, ierr)
#endif

c     if mesh is dynamic, have root processor input initial parameters
#ifdef DYNAMIC
      if (myid == 0) then
         open(imesh,file='mesh.inp',form='formatted')
         read(imesh,meshsize)
         read(imesh,meshprocs)
         read(imesh,bdryconds)
         close(imesh)
      endif

c     have root processor send out mesh information to all other procs
#ifdef PARALLEL
      meshbuffer(1) = nx
      meshbuffer(2) = ny
      meshbuffer(3) = nz
      meshbuffer(4) = nghost
      meshbuffer(5) = xprocs
      meshbuffer(6) = yprocs
      meshbuffer(7) = zprocs
      meshbuffer(8) = xbc
      meshbuffer(9) = ybc
      meshbuffer(10) = zbc
      call MPI_Bcast(meshbuffer,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nx = meshbuffer(1)
      ny = meshbuffer(2)
      nz = meshbuffer(3)
      nghost = meshbuffer(4)
      xprocs = meshbuffer(5)
      yprocs = meshbuffer(6)
      zprocs = meshbuffer(7)
      xbc = meshbuffer(8)
      ybc = meshbuffer(9)
      zbc = meshbuffer(10)
#endif      
#else
      if (myid == 0) then
         open(imesh,file='mesh.inp',form='formatted')
         read(imesh,bdryconds)
         close(imesh)
      endif

c     have root processor send out mesh information to all other procs
#ifdef PARALLEL
      meshbuffer(1) = xbc
      meshbuffer(2) = ybc
      meshbuffer(3) = zbc
      call MPI_Bcast(meshbuffer,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      xbc = meshbuffer(1)
      ybc = meshbuffer(2)
      zbc = meshbuffer(3)
#endif      
#endif

c     set up dimension size
      ndims = 3

c     set up initial mesh parameters
#ifdef DYNAMIC
      Nprocs = Xprocs*Yprocs*Zprocs
      nxlocal = nx/xprocs
      nylocal = ny/yprocs
      nzlocal = nz/zprocs
      nxlsize = nxlocal
      nylsize = nylocal
      nzlsize = nzlocal
#else
      nxlocal = nxlsize
      nylocal = nylsize
      nzlocal = nzlsize
#endif


c     set up parallel, periodic sizes 
#ifdef PARALLEL
      NProcAxis(1) = Xprocs            
      NProcAxis(2) = Yprocs            
      NProcAxis(3) = Zprocs
      if (xbc == BCperiodic) then
         periods(1) = .true.
      else
         periods(1) = .false.
      endif
      if (ybc == BCperiodic) then
         periods(2) = .true.
      else
         periods(2) = .false.
      endif
      if (zbc == BCperiodic) then
         periods(3) = .true.
      else
         periods(3) = .false.
      endif
      reorder = .false.

c     set up the master MPI node
      master = 0
      
c     initialize the MPI communicator, call it comm3D, and get other 
c     process information
      call MPI_cart_create( MPI_COMM_WORLD, ndims, NProcAxis, periods,
     1     reorder, comm3D, ierr)
      call MPI_comm_rank( comm3D, my_id, ierr)
      call MPI_comm_size( comm3D, num_nodes, ierr)
      
c     output relevant processor information
      if (my_id == 0) then
         write(6,*) '  '
      endif
      write(6,'(A,i4,A,i4)') '  Started proc #', my_id, ' /', num_nodes
      
      if (num_nodes .lt. Nprocs) then
         write(istderr,*) 'SORRY: PROGRAM ONLY WORKS FOR ', Nprocs,
     -        ' PROCESSES.'
         call MPI_Finalize(ierr)
         stop
      endif
      
c     check for invalid combinations of dimensions and
c     process-grid dimensions
      error = .FALSE.
      if (NX - (Xprocs-1)*NXlsize .le. 0) then
         print *, 'invalid combination of NX, Xprocs'
         error = .TRUE.
      endif
      if (NY - (Yprocs-1)*NYlsize .le. 0) then
         print *, 'invalid combination of NY, Yprocs'
         error = .TRUE.
      endif
      if (NZ - (Zprocs-1)*NZlsize .le. 0) then
         print *, 'invalid combination of NZ, Zprocs'
         error = .TRUE.
      endif
      if (error) STOP
      
c     if a valid processor gridding, get process id and run simulation
      if (my_id .lt. Nprocs) then
         call pxpypz (my_id+1, ipx, ipy, ipz)
         call mesh_grid_proc(ipx,ipy,ipz)
      endif
      
c     close MPI channels
      call MPI_Finalize(ierr)
      

c     This is a single-processor simulation, set process id and run simulation
#else
      ipx = 1
      ipy = 1
      ipz = 1
      call mesh_grid_proc(ipx,ipy,ipz)
#endif
      
      end
c=================================================================



      subroutine mesh_grid_proc(a_iprocx, a_iprocy, a_iprocz)
c-----------------------------------------------------------------
c     Description: main program for grid process
c
c     Arguments:
c        a_iprocx - (int) x-location of subdomain in process grid
c        a_iprocy - (int) y-location of subdomain in process grid
c        a_iprocz - (int) z-location of subdomain in process grid
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh_parms
      use mesh_common
      
c======= Declarations =========
      implicit none
      integer           :: a_iprocx, a_iprocy, a_iprocz
      integer, external :: iprocnum


c======= Internals ============

c     set process grid information into shared variables for mesh
      iprocx = a_iprocx
      iprocy = a_iprocy
      iprocz = a_iprocz
      iproc_idx = iprocnum(iprocx, iprocy, iprocz)

c     set sizes of local grid into shared variables for mesh
      NXlocal = NXlsize
      if (iprocx .eq. Xprocs) NXlocal = NX-(Xprocs-1)*NXlsize
      NYlocal = NYlsize
      if (iprocy .eq. Yprocs) NYlocal = NY-(Yprocs-1)*NYlsize
      NZlocal = NZlsize
      if (iprocz .eq. Zprocs) NZlocal = NZ-(Zprocs-1)*NZlsize

c     call main MHD simulation routine
      call MHDMain
      
      return
      end
c=================================================================
