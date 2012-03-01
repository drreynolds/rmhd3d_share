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
c-----------------------------------------------------------------

      subroutine SetBoundaryValues(u)
c-----------------------------------------------------------------
c     Description: Boundary conditions enforcement routine.  Since
c       boundary conditions are enforced at the outer faces of the
c       domain, but values are stored at cell centers, all
c       conditions are enforced through manipulation of ghost zone
c       entries. 
c     
c     Types of boundary conditions:
c       Periodic -- this routine only enforces periodic BCs in 
c                   non-parallel directions.
c       Reflecting -- homogeneous Neumann BCs for the density, 
c                   tangential velocity, tangential magnetic field 
c                   and total energy;  homogeneous Dirichlet BCs 
c                   for the normal velocity and magnetic field.
c       Zero Gradient -- homogeneous Neumann for all fields
c
c     Input: u - conserved variables at cell centers
c     Output: Ghost values of u at physical boundaries are 
c             over-written
c-----------------------------------------------------------------
      use mesh_parms
      use mesh_common
      use boundary_conds
      double precision:: u(ixlo:ixhi,iylo:iyhi,izlo:izhi,nvar)
      integer i,j,k,l,n

c X-Boundary Conditions are performed for all problems

      if (xbc == BCperiodic) then
c     x-periodic: only set values in single-processor case
c	 xl Boundary
         if(xprocs.eq.1) then
            i=0
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do j=iylo,iyhi,1
                     do n=0,nghost-1,1
                        u(i-n,j,k,l) = u(nx-n,j,k,l)
                     enddo
                  enddo
               enddo
            enddo
c        xr Boundary
            i=nxlocal+1
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do j=iylo,iyhi,1
                     do n=0,nghost-1,1
                        u(i+n,j,k,l) = u(n+1,j,k,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
      else if (xbc == BCreflecting) then
c     x-reflecting
c        xl Boundary
         if (iprocx .eq. 1) then
            i=0
            do k=izlo,izhi,1
               do j=iylo,iyhi,1
                  do n=0,nghost-1,1
                     u(i-n,j,k,1) =  u(i+n+1,j,k,1)
                     u(i-n,j,k,2) = -u(i+n+1,j,k,2)
                     u(i-n,j,k,3) =  u(i+n+1,j,k,3)
                     u(i-n,j,k,4) =  u(i+n+1,j,k,4)
                     u(i-n,j,k,5) = -u(i+n+1,j,k,5)
                     u(i-n,j,k,6) =  u(i+n+1,j,k,6)
                     u(i-n,j,k,7) =  u(i+n+1,j,k,7)
                     u(i-n,j,k,8) =  u(i+n+1,j,k,8)
                  enddo
               enddo
            enddo
         endif
c     
c        xr Boundary
         if (iprocx .eq. xprocs) then
            i=nxlocal+1
            do k=izlo,izhi,1
               do j=iylo,iyhi,1
                  do n=0,nghost-1,1
                     u(i+n,j,k,1) =  u(i-n-1,j,k,1)
                     u(i+n,j,k,2) = -u(i-n-1,j,k,2)
                     u(i+n,j,k,3) =  u(i-n-1,j,k,3)
                     u(i+n,j,k,4) =  u(i-n-1,j,k,4)
                     u(i+n,j,k,5) = -u(i-n-1,j,k,5)
                     u(i+n,j,k,6) =  u(i-n-1,j,k,6)
                     u(i+n,j,k,7) =  u(i-n-1,j,k,7)
                     u(i+n,j,k,8) =  u(i-n-1,j,k,8)
                  enddo
               enddo
            enddo
         endif
      else
c     x zero gradient
c        xl Boundary
         if (iprocx .eq. 1) then
            i=0
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do j=iylo,iyhi,1
                     do n=0,nghost-1,1
                        u(i-n,j,k,l) = u(i+n+1,j,k,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
c     
c        xr Boundary
         if (iprocx .eq. xprocs) then
            i=nxlocal+1
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do j=iylo,iyhi,1
                     do n=0,nghost-1,1
                        u(i+n,j,k,l) = u(i-n-1,j,k,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
c     
c     
c
c Y-Boundary Conditions are performed for 2D, 2.5D and 3D simulations
#ifndef ONE_D
c
      if (ybc == BCperiodic) then
c     y-periodic: only set values in single-processor case
c        yl Boundary
         if(yprocs.eq.1) then
            j=0
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do n=0,nghost-1,1
                     do i=ixlo,ixhi,1
                        u(i,j-n,k,l) = u(i,ny-n,k,l)
                     enddo
                  enddo
               enddo
            enddo
c        yr Boundary
            j=nylocal+1
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do n=0,nghost-1,1
                     do i=ixlo,ixhi,1
                        u(i,j+n,k,l) = u(i,n+1,k,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
      else if (ybc == BCreflecting) then
c     y-reflecting
c        yl Boundary
         if (iprocy .eq. 1) then
            j=0
            do k=izlo,izhi,1
               do n=0,nghost-1,1
                  do i=ixlo,ixhi,1
                     u(i,j-n,k,1) =  u(i,j+n+1,k,1)
                     u(i,j-n,k,2) =  u(i,j+n+1,k,2)
                     u(i,j-n,k,3) = -u(i,j+n+1,k,3)
                     u(i,j-n,k,4) =  u(i,j+n+1,k,4)
                     u(i,j-n,k,5) =  u(i,j+n+1,k,5)
                     u(i,j-n,k,6) = -u(i,j+n+1,k,6)
                     u(i,j-n,k,7) =  u(i,j+n+1,k,7)
                     u(i,j-n,k,8) =  u(i,j+n+1,k,8)
                  enddo
               enddo
            enddo
         endif
c     
c        yr Boundary
         if (iprocy .eq. yprocs) then
            j=nylocal+1
            do k=izlo,izhi,1
               do n=0,nghost-1,1
                  do i=ixlo,ixhi,1
                     u(i,j+n,k,1) =  u(i,j-n-1,k,1)
                     u(i,j+n,k,2) =  u(i,j-n-1,k,2)
                     u(i,j+n,k,3) = -u(i,j-n-1,k,3)
                     u(i,j+n,k,4) =  u(i,j-n-1,k,4)
                     u(i,j+n,k,5) =  u(i,j-n-1,k,5)
                     u(i,j+n,k,6) = -u(i,j-n-1,k,6)
                     u(i,j+n,k,7) =  u(i,j-n-1,k,7)
                     u(i,j+n,k,8) =  u(i,j-n-1,k,8)
                  enddo
               enddo
            enddo
         endif
      else
c     y zero gradient
c        yl Boundary
         if (iprocy .eq. 1) then
            j=0
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do n=0,nghost-1,1
                     do i=ixlo,ixhi,1
                        u(i,j-n,k,l) = u(i,j+1+n,k,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
c     
c        yr Boundary
         if (iprocy .eq. yprocs) then
            j=nylocal+1
            do l=1,nvar,1
               do k=izlo,izhi,1
                  do n=0,nghost-1,1
                     do i=ixlo,ixhi,1
                        u(i,j+n,k,l) =  u(i,j-n-1,k,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
c       
c
c Z-Boundary Conditions are performed for 3D simulations
#ifndef TWO_D
c
      if (zbc == BCperiodic) then
c     z-periodic: only set values in single-processor case
c	 zl Boundary
         if(zprocs.eq.1) then
            k=0
            do l=1,nvar,1
               do n=0,nghost-1,1
                  do j=iylo,iyhi,1
                     do i=ixlo,ixhi,1
                        u(i,j,k-n,l) = u(i,j,nz-n,l)
                     enddo
                  enddo
               enddo
            enddo
c     
c	 zr Boundary
            k=nzlocal+1
            do l=1,nvar,1
               do n=0,nghost-1,1
                  do j=iylo,iyhi,1
                     do i=ixlo,ixhi,1
                        u(i,j,k+n,l) = u(i,j,n+1,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
      else if (zbc == BCreflecting) then
c     z-reflecting
c	 zl Boundary
         if (iprocz .eq. 1) then
            k=0
            do n=0,nghost-1,1
               do j=iylo,iyhi,1
                  do i=ixlo,ixhi,1
                     u(i,j,k-n,1) =  u(i,j,k+1+n,1)
                     u(i,j,k-n,2) =  u(i,j,k+1+n,2)
                     u(i,j,k-n,3) =  u(i,j,k+1+n,3)
                     u(i,j,k-n,4) = -u(i,j,k+1+n,4)
                     u(i,j,k-n,5) =  u(i,j,k+1+n,5)
                     u(i,j,k-n,6) =  u(i,j,k+1+n,6)
                     u(i,j,k-n,7) = -u(i,j,k+1+n,7)
                     u(i,j,k-n,8) =  u(i,j,k+1+n,8)
                  enddo
               enddo
            enddo
         endif
c     
c	 zr Boundary
         if (iprocz .eq. zprocs) then
            k=nzlocal+1
            do n=0,nghost-1,1
               do j=iylo,iyhi,1
                  do i=ixlo,ixhi,1
                     u(i,j,k+n,1) =  u(i,j,k-1-n,1)
                     u(i,j,k+n,2) =  u(i,j,k-1-n,2)
                     u(i,j,k+n,3) =  u(i,j,k-1-n,3)
                     u(i,j,k+n,4) = -u(i,j,k-1-n,4)
                     u(i,j,k+n,5) =  u(i,j,k-1-n,5)
                     u(i,j,k+n,6) =  u(i,j,k-1-n,6)
                     u(i,j,k+n,7) = -u(i,j,k-1-n,7)
                     u(i,j,k+n,8) =  u(i,j,k-1-n,8)
                  enddo
               enddo
            enddo
         endif
      else
c     z zero gradient 
c	 zl Boundary
         if (iprocz .eq. 1) then
            k=0
            do l=1,nvar,1
               do n=0,nghost-1,1
                  do j=iylo,iyhi,1
                     do i=ixlo,ixhi,1
                        u(i,j,k-n,l) = u(i,j,k+1+n,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
c       
c	 zr Boundary
         if (iprocz .eq. zprocs) then
            k=nzlocal+1
            do l=1,nvar,1
               do n=0,nghost-1,1
                  do j=iylo,iyhi,1
                     do i=ixlo,ixhi,1
                        u(i,j,k+n,l) = u(i,j,k-1-n,l)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
#endif
#endif
	return
	end
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



