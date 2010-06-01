c -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c-----------------------------------------------------------------
c     $Log: setuplocaldomain.F,v $
c=================================================================



      subroutine SetupLocalDomain
c-----------------------------------------------------------------
c     Description: wrapper for RectangularDomain routine
c-----------------------------------------------------------------     
c======= Inclusions ===========
      use mesh
      use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif

c======= Declarations =========
      integer :: Nxleft, NYleft, NZleft
      
c======= Internals ============

c     determine local domain location in general domain
#ifdef PARALLEL
      if (Xprocs > 1) then
         Nxleft = 0
         if (iprocx < Xprocs) then
            call MPI_Send(Nxlocal,1,MPI_INTEGER,right,101,comm3d,ierr)
         endif
         if (iprocx > 1) then
            call MPI_Recv(Nxleft, 1, MPI_INTEGER, left, 101, comm3d, 
     &           status, ierr)
         endif
         XLloc = xl + Nxleft*(iprocx-1)*dx
         XRloc = XLloc + Nxlocal*dx
      else
         XLloc = xl
         XRloc = xr
      endif

      if (Yprocs > 1) then
         Nyleft = 0
         if (iprocy < Yprocs) then
            call MPI_Send(Nylocal,1,MPI_INTEGER,top,202,comm3d,ierr)
         endif
         if (iprocy > 1) then
            call MPI_Recv(Nyleft, 1, MPI_INTEGER, bottom, 202, comm3d, 
     &           status, ierr)
         endif
         YLloc = yl + Nyleft*(iprocy-1)*dy
         YRloc = YLloc + Nylocal*dy
      else
         YLloc = yl
         YRloc = yr
      endif

      if (Zprocs > 1) then
         Nzleft = 0
         if (iprocz < Zprocs) then
            call MPI_Send(Nzlocal,1,MPI_INTEGER,forward,303,comm3d,ierr)
         endif
         if (iprocz > 1) then
            call MPI_Recv(Nzleft, 1, MPI_INTEGER, behind, 303, comm3d, 
     &           status, ierr)
         endif
         ZLloc = zl + Nzleft*(iprocz-1)*dz
         ZRloc = ZLloc + Nzlocal*dz
      else
         ZLloc = zl
         ZRloc = zr
      endif
#else
      XLloc = xl
      XRloc = xr
      YLloc = yl
      YRloc = yr
      ZLloc = zl
      ZRloc = zr
#endif


c     update xc, yc, zc arrays
      do k = IZLO,IZHI,1
         zc(k) = ZLloc + (k-0.5d0)*dz
      enddo
      do j = IYLO,IYHI,1
         yc(j) = YLloc + (j-0.5d0)*dy
      enddo
      do i = IXLO,IXHI,1
         xc(i) = XLloc + (i-0.5d0)*dx
      enddo
      

      return
      end subroutine SetupLocalDomain
c=================================================================



