c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------
c     $Log: profiling.F,v $
c=================================================================


      module profiling
c-----------------------------------------------------------------
c     Description: contains information used in obtaining timing 
c                  profiling information.  
c
c     Contains:
c        nprof - (int) the number of total possible profiles
c        tprof - (dbl array) total timings for each profile
c        cprof - (int array) total calls of each profile
c       pstart - (dbl array) start times used by each profile
c       istart - (int array) start times used by each profile
c
c        prof_start (subroutine) begins the timing for a given 
c                   profile
c        prof_stop  (subroutine) ends the timing for a given 
c                   profile
c        prof_time  (function) returns the current profiler 
c                   timing information
c        prof_count (function) returns the current profiler
c                   count information
c
c-----------------------------------------------------------------
      save

      integer, parameter, private :: nprof = 100
      double precision,   private :: tprof(nprof)  = 0.d0
      integer,            private :: cprof(nprof)  = 0
      double precision,   private :: pstart(nprof) = -1.d0
      integer,            private :: istart(nprof) = -1
      
      contains
      

c     ------------------------------------------------------------
c     starts the clock for the given profiler
c     ------------------------------------------------------------
      subroutine prof_start(iprof)
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer :: iprof, count, rate, max

      if ((iprof < 1) .or. (iprof > nprof)) then
         write(0,*) 'prof_start error: iprof =',iprof,
     &        ' outside allowed range of 1 to',nprof
         return
      endif
#ifdef PARALLEL
      if (pstart(iprof) /= -1.d0) then
         write(0,*) 'prof_start error: cannot nest profiling calls',
     &        ', iprof =',iprof
         return
      else
         pstart(iprof) = mpi_wtime()
      endif
#else
      if (istart(iprof) /= -1) then
         write(0,*) 'prof_start error: cannot nest profiling calls',
     &        ', iprof =',iprof
         return
      else
         call system_clock(count, rate, max)
         istart(iprof) = count
      endif
#endif
      return

      end subroutine prof_start
c     ------------------------------------------------------------

      
c     ------------------------------------------------------------
c     stops given profiler clock and increments the counter
c     ------------------------------------------------------------
      subroutine prof_stop(iprof)
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer :: iprof, count, rate, max
      
      if ((iprof < 1) .or. (iprof > nprof)) then
         write(0,*) 'prof_stop error: iprof =',iprof,
     &        ' outside allowed range of 1 to',nprof
         return
      endif
#ifdef PARALLEL
      if (pstart(iprof) == -1.d0) then
         write(0,*) 'prof_stop error: timer was not started!',
     &        ', iprof =',iprof
         return
      else
         tprof(iprof) = tprof(iprof) + mpi_wtime() - pstart(iprof)
      endif
      pstart(iprof) = -1.d0
#else
      if (istart(iprof) == -1) then
         write(0,*) 'prof_stop error: timer was not started!',
     &        ', iprof =',iprof
         return
      else
         call system_clock(count, rate, max)
         if (count < istart(iprof)) then
            tprof(iprof) = (max-istart(iprof)+count)/rate
         else
            tprof(iprof) = (count-istart(iprof))/rate
         endif
      endif
      istart(iprof) = -1
#endif
      cprof(iprof) = cprof(iprof) + 1
      return

      end subroutine prof_stop
c     ------------------------------------------------------------

      
c     ------------------------------------------------------------
c     returns the current total time for the profiler
c     ------------------------------------------------------------
      function prof_time(iprof)
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer          :: iprof
      double precision :: prof_time

      prof_time = 0.d0
      if ((iprof < 1) .or. (iprof > nprof)) then
         write(0,*) 'prof_time error: iprof =',iprof,
     &        ' outside allowed range of 1 to',nprof
         return
      endif
      prof_time = tprof(iprof)
      return

      end function prof_time
c     ------------------------------------------------------------

      
c     ------------------------------------------------------------
c     returns the current total count for the profiler
c     ------------------------------------------------------------
      function prof_count(iprof)
#ifdef PARALLEL
      use mpistuff
#endif
      implicit none
      integer :: iprof, prof_count

      prof_count = 0
      if ((iprof < 1) .or. (iprof > nprof)) then
         write(0,*) 'prof_count error: iprof =',iprof,
     &        ' outside allowed range of 1 to',nprof
         return
      endif
      prof_count = cprof(iprof)
      return

      end function prof_count
c     ------------------------------------------------------------

      
      end module profiling
c=================================================================
