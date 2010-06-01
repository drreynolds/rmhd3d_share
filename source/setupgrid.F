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
c     $Log: setupgrid.F,v $
c=================================================================



      subroutine SetupDomain
c-----------------------------------------------------------------
c     Description: wrapper for RectangularDomain routine
c-----------------------------------------------------------------     
c======= Inclusions ===========
      use mesh
      use mesh_common

c======= Internals ============
      call RectangularDomain
      
      return
      end
c=================================================================



      subroutine RectangularDomain
c-----------------------------------------------------------------
c     Description: sets up domain parameters and directional meshes
c-----------------------------------------------------------------
c======= Inclusions ===========    
      use mesh
      use mesh_common

c======= Declarations =========
      integer :: i,j, k
      integer :: ig, jg, kg
      double precision :: w, l


c======= Internals ============

c     set up directional mesh spacings
      dx = (xr-xl)/nx
      dy = (yr-yl)/ny
      dz = (zr-zl)/nz

c     root node outputs discretization parameters to stdout
      if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
         write(6,'(2(A,1x,es12.5,1x))') 'xl=',xl, 'xr=',xr
         write(6,'(2(A,1x,es12.5,1x))') 'yl=',yl, 'yr=',yr
         write(6,'(2(A,1x,es12.5,1x))') 'zl=',zl, 'zr=',zr
         write(6,'(A,1x,es12.5)') 'dx=', dx
         write(6,'(A,1x,es12.5)') 'dy=', dy
         write(6,'(A,1x,es12.5)') 'dz=', dz
      endif

c     get the global indices
      call iglobal(iprocx,1,ig)
      call jglobal(iprocy,1,jg)
      call kglobal(iprocz,1,kg)

c     set up the directional meshes for x,y,z
      do k = IZLO,IZHI,1
         zc(k) = zl+(kg+k-1)*dz-0.5d0*dz
      enddo
      do j = IYLO,IYHI,1
         yc(j) = yl+(jg+j-1)*dy-0.5d0*dy
      enddo
      do i = IXLO,IXHI,1
         xc(i) = xl+(ig+i-1)*dx-0.5d0*dx
      enddo

      return
      end
c=================================================================
