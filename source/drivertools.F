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

      subroutine Initialize(new,ux,phi,dt,ttot,istart,lastiter)
c-----------------------------------------------------------------
c     Description: initializes fluid variables, dt, opens I/O 
c        streams,  determines multigrid information 
c     
c     Arguments:
c             new - (out) init flag, 1 => do not input from file
c              ux - (out) fluid variables over entire mesh
c             phi - (out) magnetic field density
c              dt - (out) initial timestep size
c            ttot - (out) initial time (if input from file)
c          istart - (out) initial timelevel (if input from file)
c        lastiter - (out) previous timelevel (if input from file)
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh_parms
      use mesh_common
      use iounits
      use mgparams
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      double precision :: dt, ttot
      integer          :: istart, new, lastiter


c======= Internals ============

c     set initial experimental conditions, time, timestep
      call InitialConditions(ux,phi)
      ttot=0.D0
      call newdt(dt,ux)

c     input initial information from file if required
      if(new/=1) then
         open(itime,file='time.txt',form='formatted')
         read(itime,*) dt
         read(itime,*) ttot
         read(itime,*) lastiter
         close(itime)
         istart = istart + lastiter
         call DumpInFluidParallel(ux,phi)
      endif

c     root node opens output files
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         open(idiagCons,file='conservation.txt',form='formatted')
         open(idiagEng,file='energyHistory.txt',form='formatted')
         open(idiagflux,file='fluxHistory.txt',form='formatted')
         open(idiag,file='Diagnostic.history',form='formatted')
      endif

c     determine overall, local maximum levels for Multigrid
      call DetermineMaxLevel(lMax,nx,ny)
      call DetermineMaxLevel(lMaxLocal,nxlocal,nylocal)

      return
      end subroutine Initialize
c=================================================================




      subroutine WrapUp(ux,phi,dt,ttot,timestep)
c-----------------------------------------------------------------
c     Description: Outputs final fluid information, closes IO 
c        streams
c     
c     Arguments:
c              ux - (in) fluid variables over entire mesh
c             phi - (in) magnetic field density
c              dt - (in) current timestep size
c            ttot - (in) current time
c        timestep - (in) current timelevel
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh_parms
      use mesh_common
      use iounits
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(IXLO:IXHI,IYLO:IYHI)
      double precision :: dt, ttot
      integer          :: timestep


c======= Internals ============

c     root node outputs timelevel information
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         open(itime,file='time.txt',form='formatted')
         write(itime,*) dt
         write(itime,*) ttot
         write(itime,*) timestep-1
         close(itime)
      endif

c     output final fluid information
      call DumpOutFluidParallel(ux,phi)

c     root node closes files
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         close(idiagCons)
         close(idiagEng)
         close(idiagflux)
         close(idiag)
      endif

      return
      end subroutine WrapUp
c=================================================================




      subroutine Timing
c-----------------------------------------------------------------
c     Description: calls timing information to track performance
c
c     Note: this is virtually unchanged from Ravi's code, it has 
c        only been cleaned up for commenting and clarity.
c-----------------------------------------------------------------
c======= Declarations =========
      real :: extime(2)
      real :: etime, dtime
      character(len=24) :: fdate
      character(len=10) :: mytime, myzone
      integer :: mydate
      integer :: mydate_time(8)


c======= Internals ============

c     for each system, root node calls timing information
#ifdef INTEL
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(18,*) iprocx, iprocy, fdate()
         write(18,*) iprocx, iprocy, etime(extime), extime(1), extime(2)
         write(18,*) iprocx, iprocy, dtime(extime), extime(1), extime(2)
      endif
#endif
#ifdef CRAY
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         call date(mydate)
         write(18,*) mydate
         write(18,*) timef(), second(), tsecnd()
      endif  
#endif
#ifdef SGI
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(18,*)  etime(extime), extime(1), extime(2)
         write(18,*)  dtime(extime), extime(1), extime(2)
      endif
#endif
#ifdef PC
      if(iprocx==1 .and. iprocy==1 .and. iprocz==1) then
         write(18,*)  etime(extime), extime(1), extime(2)
         write(18,*)  dtime(extime), extime(1), extime(2)
      endif
#endif
      return
      end subroutine Timing
c=================================================================
