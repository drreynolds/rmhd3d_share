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

      subroutine OutputData(ux, phi, idump, output_flag, binary_flag, 
     &     dt, ttot)
c-----------------------------------------------------------------
c Description: wrapper to call desired data output routine.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_parms
      use mesh_common
      use mpistuff
      use profiling

c======= Declarations =========
      implicit none
      real*8 :: ux(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      real*8 :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      real*8 :: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      real*8 :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      real*8 :: jcurrent(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      real*8 :: dt, ttot
      integer :: idump, output_flag, binary_flag
      integer :: nxtot, nytot, nztot, xoff, yoff, zoff
      integer, parameter :: nodal=1

c======= Internals ============

c     start the I/O profiler
      call prof_start(6)

c     dump restart information only
      if (output_flag==2) then

         if(iprocx*iprocy*iprocz == 1) then
            open(444,file='time.txt',form='formatted')
            write(444,*) dt
            write(444,*) ttot
            write(444,*) idump-1
            close(444)
         endif
         call DumpOutFluidParallel(ux,phi)
      endif

c     normal output
      if (output_flag==1) then

c        create output directory for this dump
         call mkdir_output(idump, my_id, ierr)
         if (ierr /= 0) then
            write(0,*) 'p',my_id,': error in mkdir_output =',ierr
            call ErrorHandler(ierr, ERROR_OTHER)
         endif

c        if text output requested, write Matlab/Python-based output
         if (binary_flag==0) then

            call WriteGridFiles()
c            call WriteGnuplotFile(ux,idump)
            call WriteMatlabFile(ux,idump)

c        otherwise output AVS data or SILO data, depending on how 
c        the code was compiled
         else
#ifdef SILO
c           get primitive variables, div B, toroidal current
            call ConservedToPrimitive(ux,vx)
            call DivergenceCD(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &           ux(ixlo,iylo,izlo,7), div)	
            call ToroidalCurrent(ux(ixlo,iylo,izlo,5),
     &           ux(ixlo,iylo,izlo,6), ux(ixlo,iylo,izlo,7), jcurrent)	
      
            nxtot = ixhi-ixlo+1
            nytot = iyhi-iylo+1
            nztot = izhi-izlo+1
            xoff = 1-ixlo
            yoff = 1-iylo
            zoff = 1-izlo
            call Write_silo(vx, div, jcurrent, idump, xc, yc, zc, 
     &           nxlocal, nylocal, nzlocal, nvar, nxtot, nytot, nztot, 
     &           xoff, yoff, zoff, dx, dy, dz, iprocx, iprocy, iprocz, 
     &           nodal, my_id, nprocs, ierr)
            if (ierr /= 0) then
               write(0,*) 'p',my_id,': error in write_silo =',ierr
               call ErrorHandler(ierr, ERROR_OTHER)
            endif
#else
c            call WriteBinaryFileParallel(ux,idump)
            call WriteAVSFile(ux,idump)
#endif
         endif

      endif

c     stop the I/O profiler
      call prof_stop(6)

      return
      end subroutine OutputData
c=================================================================



      subroutine DumpOutFluidParallel(ux,phi)
c-----------------------------------------------------------------
c     Description: parallel file output of conserved variables, etc.
c
c     Arguments:
c          ux - (array input) conserved variables
c         phi - (array input) projection information for div B
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh_parms
      use mesh_common
      use properties
      use iounits

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      integer          :: FileNumber
      integer          :: i, j, k, l
      character*50     :: DumpFile
      

c======= Internals ============

c     set output filename
      FileNumber = iprocx*10000 + iprocy*100 + iprocz
      write(DumpFile,'(4Hdump,i6.6,4H.bin)')  FileNumber
      
c     output information to file
      open(idump, file=DumpFile, form='unformatted')
      write(6,*) 'p',iproc_idx-1,': writing to file',DumpFile
      write(idump) rgas, gamma, wmol
      write(idump) ux
      write(idump) phi
      close(idump)
      
      return
      end
c=================================================================



      subroutine DumpInFluidParallel(ux,phi)
c-----------------------------------------------------------------
c     Description: parallel file input of conserved variables, etc.
c
c     Arguments:
c          ux - (array output) conserved variables
c         phi - (array output) projection information for div B
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh_parms
      use mesh_common
      use properties
      use iounits

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: phi(ixlo:ixhi,iylo:iyhi,izlo:izhi)
      integer          :: FileNumber
      character*50     :: DumpFile
      
c======= Internals ============

c     set input filename
      FileNumber = iprocx*10000 + iprocy*100 + iprocz
      write(DumpFile,'(4Hdump,i6.6,4H.bin)')  FileNumber
      
c     input information from file
      open(idump, file=DumpFile, form='unformatted')
      write(6,*) 'p',iproc_idx-1,': reading from file',DumpFile
      read(idump) rgas, gamma, wmol
      read(idump) ux
      read(idump) phi
      close(idump)
      
      return
      end
c=================================================================



      subroutine WriteBinaryFileParallel(ux,index)
c-----------------------------------------------------------------
c     Description: binary parallel file output of conserved 
c        variables
c
c     Arguments:
c            ux - (array input) conserved variables
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use mpistuff
      use iounits

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      integer          :: index
      character*50     :: outfile


c======= Internals ============

c     set output filename
      write(outfile,'(4Hdata,i6.6,7H/output,i6.6,4H.bin)') 
     &     index, iproc_idx-1

c     output information to file
      open(ibinoutput,file=outfile,form='unformatted')
      write(ibinoutput) ux
      close(ibinoutput)

      end
c=================================================================



      subroutine WriteTextFileParallel(ux,index)
c-----------------------------------------------------------------
c     Description: text parallel file output of conserved 
c        variables
c
c     Arguments:
c            ux - (array input) conserved variables
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits
      
c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: press, rho, u, v, w, bi, bj, bk
      integer          :: i, j, k
      integer          :: index
      integer          :: itmp, ierr, nelements
      character*50     :: outfile


c======= Internals ============

c     set output filename, and open iostream
      write(outfile,'(4Hdata,i6.6,7H/output,i6.6,4H.txt)') 
     &     index, iproc_idx-1
      open(itecoutput,file=outfile,form='formatted')

c     get div B
      call DivergenceCD(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),div)	
c     output variables and grid location
c      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
      write(itecoutput,*) 'Variables=x,y,z,rho,u,v,w,bx,by,bz,p'
      write(itecoutput,*) 'Zone',' i=',nxlsize,
     &     ' j=',nylsize,' k=',nzlsize
c      endif
      
c     output variables to file
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               rho = ux(i,j,k,1)
               u   = ux(i,j,k,2)/ux(i,j,k,1)
               v   = ux(i,j,k,3)/ux(i,j,k,1)
               w   = ux(i,j,k,4)/ux(i,j,k,1)
               bi  = ux(i,j,k,5)
               bj  = ux(i,j,k,6)
               bk  = ux(i,j,k,7)
               press = (gamma-1.0D0)*(ux(i,j,k,8)
     &              -0.5D0*rho*(u*u+v*v+w*w)-0.5D0*(bi*bi+bj*bj+bk*bk))
               write(itecoutput,*) xc(i), yc(j), zc(k),
     &              rho, u, v, w, bi, bj, bk, press, div(i,j,k)
            enddo
         enddo
      enddo

c     close output file
      close(itecoutput)

      end
c=================================================================



      subroutine WriteGnuplotFile(ux,index)
c-----------------------------------------------------------------
c     Description: output conserved variables to file for use in 
c        gnuplot.
c
c     Arguments:
c            ux - (array input) conserved variables
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits
      
c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: jcurrent(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: press, rho, u, v, w, bi, bj, bk
      integer          :: i, j, k
      integer          :: index
      integer          :: nbytes, itmp, ierr, nelements
      character*50     :: outfile


c======= Internals ============

c     set output filename and open iostream
      write(outfile,'(4Hdata,i6.6,7H/output,i6.6,4H.txt)') 
     &     index, iproc_idx-1
      open(itecoutput,file=outfile,form='formatted')
      
c     get div B, toroidal current
      call DivergenceCD(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),div)	
      call ToroidalCurrent(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),jcurrent)	
      
c     output variables to file
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               rho = ux(i,j,k,1)
               u   = ux(i,j,k,2)/ux(i,j,k,1)
               v   = ux(i,j,k,3)/ux(i,j,k,1)
               w   = ux(i,j,k,4)/ux(i,j,k,1)
               bi  = ux(i,j,k,5)
               bj  = ux(i,j,k,6)
               bk  = ux(i,j,k,7)
               press = (gamma-1.0D0)*(ux(i,j,k,8)
     &              -0.5D0*rho*(u*u+v*v+w*w)-0.5D0*(bi*bi+bj*bj+bk*bk))
               write(itecoutput,2000) xc(i), yc(j), zc(k),
     &              rho, u, v, w, bi, bj, bk,
     &              press, div(i,j,k), jcurrent(i,j,k),
     &              press+0.5D0*(bi*bi+bj*bj+bk*bk)
            enddo
 2000       format(1h ,14(es21.12e3,1x))
            write(itecoutput,*)
         enddo
      enddo

c     close output file
      close(itecoutput)

      end
c=================================================================



      subroutine WriteMatlabFile(ux,index)
c-----------------------------------------------------------------
c     Description: output conserved variables to file for use in 
c        Matlab.
c
c     Arguments:
c            ux - (array input) conserved variables
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits
      
c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: jcurrent(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: press, rho, u, v, w, bi, bj, bk
      integer          :: i, j, k
      integer          :: index
      integer          :: nbytes, itmp, ierr, nelements
      character*50     :: outfile


c======= Internals ============

c     set output filename and open iostream
      write(outfile,'(4Hdata,i6.6,7H/output,i6.6,4H.txt)') 
     &     index, iproc_idx-1
      open(itecoutput,file=outfile,form='formatted')
      
c     get div B, toroidal current
      call DivergenceCD(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),div)	
      call ToroidalCurrent(ux(ixlo,iylo,izlo,5),ux(ixlo,iylo,izlo,6),
     &     ux(ixlo,iylo,izlo,7),jcurrent)	
      
c     output variables in one long column, one variable at a time
c        rho
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,1)
            enddo
         enddo
      enddo
c        u
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
              write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,2)/ux(i,j,k,1)
            enddo
         enddo
      enddo
c        v
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
              write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,3)/ux(i,j,k,1)
            enddo
         enddo
      enddo
c        w
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
              write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,4)/ux(i,j,k,1)
            enddo
         enddo
      enddo
c        bx
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,5)
            enddo
         enddo
      enddo
c        by
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,6)
            enddo
         enddo
      enddo
c        bz
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,'(1x,es21.12e3)') ux(i,j,k,7)
            enddo
         enddo
      enddo
c        press
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               rho = ux(i,j,k,1)
               u   = ux(i,j,k,2)/ux(i,j,k,1)
               v   = ux(i,j,k,3)/ux(i,j,k,1)
               w   = ux(i,j,k,4)/ux(i,j,k,1)
               bi  = ux(i,j,k,5)
               bj  = ux(i,j,k,6)
               bk  = ux(i,j,k,7)
               press = (gamma-1.0D0)*(ux(i,j,k,8)
     &              -0.5D0*rho*(u*u+v*v+w*w)-0.5D0*(bi*bi+bj*bj+bk*bk))
               write(itecoutput,'(1x,es21.12e3)') press
            enddo
         enddo
      enddo
c        div
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,'(1x,es21.12e3)') div(i,j,k)
            enddo
         enddo
      enddo
c        jcurrent
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,'(1x,es21.12e3)') jcurrent(i,j,k)
            enddo
         enddo
      enddo
c        total pressure
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               rho = ux(i,j,k,1)
               u   = ux(i,j,k,2)/ux(i,j,k,1)
               v   = ux(i,j,k,3)/ux(i,j,k,1)
               w   = ux(i,j,k,4)/ux(i,j,k,1)
               bi  = ux(i,j,k,5)
               bj  = ux(i,j,k,6)
               bk  = ux(i,j,k,7)
               press = (gamma-1.0D0)*(ux(i,j,k,8)
     &              -0.5D0*rho*(u*u+v*v+w*w)-0.5D0*(bi*bi+bj*bj+bk*bk))
               write(itecoutput,'(1x,es21.12e3)') 
     &              press+0.5D0*(bi*bi+bj*bj+bk*bk)
            enddo
         enddo
      enddo

c     close output file
      close(itecoutput)

      end
c=================================================================



      subroutine WriteAVSFile(ux,index)
c-----------------------------------------------------------------
c     Description: output conserved variables to file for use in 
c        AVS
c
c     Arguments:
c            ux - (array input) conserved variables
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision :: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: jcurrent(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: tmp(NXlsize,NYlsize,NZlsize)
      real             :: xn(NXlsize,NYlsize,NZlsize)
      real             :: yn(NXlsize,NYlsize,NZlsize)
      real             :: zn(NXlsize,NYlsize,NZlsize)
      double precision :: press, rho, u, v, w, bi, bj, bk
      integer          :: i, j, k
      integer          :: index
      integer          :: nbytes, itmp, ierr, nelements
      character*50     :: outfile, fldfile, meshfile
      integer, parameter :: ifld=91

     
c======= Internals ============

c     set output filename
      write(outfile,'(4Hdata,i6.6,7H/outavs,i6.6)') 
     &     index, iproc_idx-1

c     root node outputs his information
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(fldfile,2000) float(index)/1000000.0
 2000    format(6Houtavs,f7.6,'.fld')
         open(ifld,file=fldfile,form='formatted')
         write(ifld,1001)
 1001    format('# AVS field file',/""
     &        'ndim=3')
         write(ifld,1002) nx, ny, nz
 1002    format('dim1=',i4,/"",
     &        'dim2=',i4,/"",
     &        'dim3=',i4,/"",
     &        'nspace=3'/"",
     &        'veclen=8'/"",
     &        'data=double'/"",
     &        'field=uniform')

         do i=1,nvar,1
            write(ifld,1004) i,outfile,
     &           4+(i-1)*nxlocal*nylocal*nzlocal*8+8*(i-1)
 1004       format('variable',i2,' file=',a50,
     &           'filetype=binary skip=',i8
     &           ' offset=0')
         enddo
      endif

c     open iostream
      open(itecoutput,file=outfile,form='unformatted')
     
      call ConservedToPrimitive(ux,vx)

c     write primitive variables to output file 
      do i = 1,nvar,1
         tmp=vx(1:nxlocal,1:nylocal,1:nzlocal,i)
         write(itecoutput) tmp
      enddo

c     close output file
      close(itecoutput)

      end
c=================================================================



      subroutine WriteGnuplotFile2(phi,index)
c-----------------------------------------------------------------
c     Description: output div B information to file for use in 
c        gnuplot
c
c     Arguments:
c           phi - (array input) div B information
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits

c======= Declarations =========
      implicit none
      double precision :: phi(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      integer          :: i, j, k
      integer          :: index
      integer          :: nbytes, itmp, ierr, nelements
      character*50     :: outfile

     
c======= Internals ============

c     set output filename and open iostream
      write(outfile,'(4Hdata,i6.6,4H/phi,i6.6)') 
     &     index, iproc_idx-1
      open(itecoutput,file=outfile,form='formatted')

c     output data to file
      do k = 1,nzlocal,1
         do j = 1,nylocal,1
            do i = 1,nxlocal,1
               write(itecoutput,2001) xc(i), yc(j),
     &              phi(i,j,k)
            enddo
 2001       format(1h ,3(e11.4,1x))
            write(itecoutput,*)
         enddo
      enddo

c     close output file
      close(itecoutput)

      end
c=================================================================



c-----------------------------------------------------------------------
      subroutine WriteAVSFile2D(ux,index)
c-----------------------------------------------------------------
c     Description: output conserved variables to file for use in 
c        AVS
c
c     Arguments:
c            ux - (array input) conserved variables
c         index - (int input) process index
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: vx(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI,nvar)
      double precision :: div(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: jcurrent(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI)
      double precision :: tmp(NXlsize,NYlsize,NZlsize)
      double precision :: rhs(NXlsize,NYlsize,NZlsize,nvar)
      double precision :: press, rho, u, v, w, bi, bj, bk
      integer          :: i,j,k
      integer          :: index
      integer          :: nbytes, itmp, ierr, nelements
      character*50     :: outfile,fldfile, meshfile
      integer,parameter ::ifld=999

     
c======= Internals ============

c     File name for data 
      write(outfile,'(4Hdata,i6.6,7H/outavs,i6.6)') 
     &     index, iproc_idx-1

c     File name for coords
      write(meshfile,'(4Hdata,i6.6,4H/xyz,i6.6)') 
     &     index, iproc_idx-1

c     Write out fld file
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(fldfile,'(4Hdata,i6.6,7H/outavs,i6.6,4H.fld)') 
     &        index, iproc_idx-1
         open(ifld,file=fldfile,form='formatted')
         write(ifld,1001)
 1001    format('# AVS field file',/""
     &        'ndim=2')
         write(ifld,1002) nx, ny
 1002    format('dim1=',i4,/"",
     &        'dim2=',i4,/"",
     &        'nspace=2'/"",
     &        'veclen=8'/"",
     &        'data=double'/"",
     &        'field=uniform')
c     &       'field=uniform')
c     
         do i=1,nvar,1
            write(ifld,1004) i,outfile,
     &           4+(i-1)*nxlocal*nylocal*8+8*(i-1)
 1004       format('variable',i3,' file=',a50,
     &           'filetype=binary skip=',i8
     &           ' offset=0')
         enddo
      endif
      open(itecoutput,file=outfile,form='unformatted')
      
c     
      call ConservedToPrimitive(ux,vx)
c     
      do i=1,nvar,1
         tmp=vx(1:nxlocal,1:nylocal,1:nzlocal,i)
         write(itecoutput) tmp
      enddo
      close(itecoutput)
c
c     
      end
c=================================================================



      subroutine WriteGridFiles()
c-----------------------------------------------------------------
c     Description: output xc, yc, zc to files, for use with text
c     outputs.
c-----------------------------------------------------------------
c======= Inclusions ===========
      use mesh
      use mesh_common
      use properties
      use iounits
      use mpistuff

c======= Declarations =========
      implicit none
      integer :: i, coords(3)
      character*50 :: outfile

c======= Internals ============

c     write xc to file
      if (iprocy*iprocz == 1) then
         write(outfile,'(6Hx.grid,f4.3)')  float(iprocx)/1000.0
         open(1505, file=outfile, form='formatted')
         do i = 1,nxlocal
            write(1505,'(es21.12e3)') xc(i)
         enddo
         close(1505)
      endif

c     write yc to file
      if (iprocx*iprocz == 1) then
         write(outfile,'(6Hy.grid,f4.3)')  float(iprocy)/1000.0
         open(1506, file=outfile, form='formatted')
         do i = 1,nylocal
            write(1506,'(es21.12e3)') yc(i)
         enddo
         close(1506)
      endif
      
c     write zc to file
      if (iprocx*iprocy == 1) then
         write(outfile,'(6Hz.grid,f4.3)')  float(iprocz)/1000.0
         open(1507, file=outfile, form='formatted')
         do i = 1,nzlocal
            write(1507,'(es21.12e3)') zc(i)
         enddo
         close(1507)
      endif
      
c     write layout to file, along with process->location table
      if (iprocx*iprocy*iprocz == 1) then
         open(1508, file='proc_layout.txt', form='formatted')
         write(1508,'(4(i4,2x))') -1, xprocs, yprocs, zprocs
         do i=0,xprocs*yprocs*zprocs-1
            call MPI_Cart_coords(comm3d, i, 3, coords, ierr)
            call ErrorHandler(ierr, ERROR_CARTCOORDS)
            write(1508,'(4(i4,2x))') i, coords(1)+1, coords(2)+1, 
     &           coords(3)+1
         enddo
         close(1508)
      endif
      
      end subroutine WriteGridFiles
!=================================================================
