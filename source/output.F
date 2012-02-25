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
      integer          :: nbytes, offset,itmp
      integer          :: FileNumber
      integer          :: i, j, k, l
      character*50     :: infile, outfile
      character*50     :: DumpFile
      

c======= Internals ============

c     set output filename
      FileNumber = iprocx*10000 + iprocy*100 + iprocz
      write(DumpFile, 1000) float(FileNumber) / 1000000
 1000 format(4Hdump, f7.6)
      
c     output information to file
      open(idump, file=DumpFile, form='unformatted')
      write(6,*) 'writing to file',iproc_idx,DumpFile
      write(idump) rgas,gamma,wmol
      write(idump)ux
      write(idump)phi
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
      integer          :: nbytes, offset, itmp
      integer          :: FileNumber
      integer          :: i, j, k, l
      character*50     :: infile, outfile
      character*50     :: DumpFile
      
c======= Internals ============

c     set input filename
      FileNumber = iprocx*10000 + iprocy*100 + iprocz
      write(DumpFile, 1000) float(FileNumber) / 1000000
 1000 format(4Hdump, f7.6)
      
c     input information from file
      open(idump, file=DumpFile, form='unformatted')
      write(6,*) 'reading from file',iproc_idx,DumpFile
      read(idump) rgas,gamma,wmol
      read(idump)ux
      read(idump)phi
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
      use iounits

c======= Declarations =========
      implicit none
      double precision :: ux(IXLO:IXHI,IYLO:IYHI,IZLO:IZHI, nvar)
c      double precision :: rho(NXlsize,Nylsize,NZlsize)
      double precision :: u, v, w
      integer          :: i, j, k
      integer          :: index
      integer          :: nbytes, offset,itmp, ierr, nelements
      character*50     :: outfile


c======= Internals ============

c     set output filename
      write(outfile,1000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 1000 format(6Houtput,f4.3,f7.6)

c     collocate output data for binary dump
c$$$      do k=1,nzlocal,1
c$$$         do j=1,nylocal,1
c$$$            do i=1,nxlocal,1
c$$$               rho(i,j,k)=ux(i,j,k,1)
c$$$            enddo
c$$$         enddo
c$$$      enddo

c     output information to file
      open(ibinoutput,file=outfile,form='unformatted')
      write(ibinoutput)ux
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
      integer          :: nbytes, offset,itmp, ierr, nelements
      character*50     :: outfile


c======= Internals ============

c     set output filename, and open iostream
      write(outfile,1000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 1000 format(6Houtput,f4.3,f7.6)
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
               write(itecoutput,*) xc(i), yc(j),zc(k),
     &              rho, u,v,w,bi,bj,bk,press,div(i,j,k)
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
      integer          :: nbytes, offset,itmp, ierr, nelements
      character*50     :: outfile


c======= Internals ============

c     set output filename and open iostream
      write(outfile,1000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 1000 format(6Houtput,f4.3,f7.6)
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
     &              rho, u,v,w,bi,bj,bk,
     &              press,div(i,j,k),jcurrent(i,j,k),
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
c      double precision :: xn(NXlsize,NYlsize,NZlsize)
c      double precision :: yn(NXlsize,NYlsize,NZlsize)
c      double precision :: zn(NXlsize,NYlsize,NZlsize)
      real             :: xn(NXlsize,NYlsize,NZlsize)
      real             :: yn(NXlsize,NYlsize,NZlsize)
      real             :: zn(NXlsize,NYlsize,NZlsize)
      double precision :: press, rho, u, v, w, bi, bj, bk
      integer          :: i, j, k
      integer          :: index
      integer          :: nbytes, offset,itmp, ierr, nelements
      character*50     :: outfile, fldfile, meshfile
      integer, parameter :: ifld=91

     
c======= Internals ============

c     set output filename
      write(outfile,1000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 1000 format(6Houtavs,f4.3,f7.6)
c     File name for coords
c$$$      write(meshfile,3000) float(iproc_idx)/1000.0, 
c$$$     &     float(index)/1000000.0
c$$$ 3000 format(3Hxyz,f4.3,f7.6)

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
c
c$$$      do i=1,3,1
c$$$      write(ifld,1003) i,meshfile,
c$$$     &        4+(i-1)*nxlocal*nylocal*nzlocal*4+8*(i-1)
c$$$c     &        4+(i-1)*nxlocal*nylocal*nzlocal*8+8*(i-1)
c$$$ 1003 format('coord ',i2,' file=',a50,
c$$$     &        'filetype=binary skip=',i8
c$$$     &     ' offset=0')
c$$$      enddo
c
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
c     
      call ConservedToPrimitive(ux,vx)

c     write primitive variables to output file 
      do i = 1,nvar,1
         tmp=vx(1:nxlocal,1:nylocal,1:nzlocal,i)
         write(itecoutput) tmp
      enddo

c     close output file
      close(itecoutput)
c
c     
c$$$      do k=1,nzlocal,1
c$$$         do j=1,nylocal,1
c$$$            do i=1,nxlocal,1
c$$$               xn(i,j,k)=Rc(i,j)*dcos(phi(k))
c$$$               yn(i,j,k)=Rc(i,j)*dsin(phi(k))
c$$$c               xn(i,j,k)=Rc(i,j)
c$$$               zn(i,j,k)=Zc(i,j)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c
c$$$      open(itecoutput,file=meshfile,form='unformatted')
c$$$      write(itecoutput) xn
c$$$      write(itecoutput) yn
c$$$      write(itecoutput) zn
c$$$      close(itecoutput)
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
      integer          :: nbytes, offset,itmp, ierr, nelements
      character*50     :: outfile

     
c======= Internals ============

c     set output filename and open iostream
      write(outfile,1000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 1000 format(3Hphi,f4.3,f7.6)
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
      integer          :: nbytes, offset,itmp, ierr, nelements
      character*50     :: outfile,fldfile, meshfile
      integer,parameter ::ifld=999

     
c======= Internals ============

c     File name for data 
      write(outfile,1000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 1000 format(6Houtavs,f4.3,f7.6)

c     File name for coords
      write(meshfile,3000) float(iproc_idx)/1000.0, 
     &     float(index)/1000000.0
 3000 format(3Hxyz,f4.3,f7.6)

c     Write out fld file
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(fldfile,2000) float(index)/1000000.0
 2000    format(6Houtavs,f7.6,'.fld')
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
