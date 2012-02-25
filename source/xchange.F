C -*- Mode: Fortran; -*-
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

	subroutine ExchangePsi(psi,imax,jmax,xChangeFlag)
c-----------------------------------------------------------------
c       Description: Exchange values at overlapping processor
c                    boundaries for psi
c       Input: psi - field with ghost which must be exchanged.
c              imax,jmax,kmax: mesh extent
c              xChangeFlag: 0 (no exchange in which case this is 
c                              an empty routine)
c                         : 1 (do the exchange)
c-----------------------------------------------------------------
	use mesh_uparms
	use mesh_common
#ifdef PARALLEL
	use mpistuff
#endif
	implicit none

        integer,parameter:: MSG_XCH_XLOW_TAG=1,MSG_XCH_XHI_TAG=2
        integer,parameter:: MSG_XCH_YLOW_TAG=3,MSG_XCH_YHI_TAG=4
        integer,parameter:: MSG_XCH_ZLOW_TAG=5,MSG_XCH_ZHI_TAG=6
        integer,parameter:: MSG_MAX_TAG=10
c
	integer:: imax,jmax,xChangeFlag
	double precision:: psi(0:imax+1,0:jmax+1)
	integer, parameter:: ng=1

#ifdef PARALLEL
	double precision:: xbuffer_send(jmax)
	double precision:: xbuffer_recv(jmax)
	double precision:: ybuffer_send(imax)
	double precision:: ybuffer_recv(imax)

	integer:: XBUFFSIZE, YBUFFSIZE, ZBUFFSIZE
	integer:: ii,jj,mm,kk,l,idest
	integer:: iprocnum
	integer:: msg_id_send_x_low
	integer:: msg_id_send_x_hi
	integer:: msg_id_recv_x_low
	integer:: msg_id_recv_x_hi

	integer:: msg_id_send_y_low
	integer:: msg_id_send_y_hi
	integer:: msg_id_recv_y_low
	integer:: msg_id_recv_y_hi

	integer:: msg_id_send_z_low
	integer:: msg_id_send_z_hi
	integer:: msg_id_recv_z_low
	integer:: msg_id_recv_z_hi
	integer:: ybufsize, zbufsize
c       

c	write(6,*) 'PSI EXCHANGE',iprocx,iprocy,xChangeFlag
	if(xChangeFlag.eq.0) return

C       -------X DIRECTION COMMUNICATION
c       Update x-low boundaries

c	 XBUFFSIZE=(ng) * (jmax+2*ng)
c	 YBUFFSIZE=(ng) * (imax+2*ng)
	 XBUFFSIZE=(ng) * (jmax)
	 YBUFFSIZE=(ng) * (imax)
	if (iprocx .gt. 1) then
c	   isource = iprocnum (iprocx-1, iprocy, iprocz)-1	   
	   call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1	   left, MSG_XCH_XHI_TAG, comm3D, msg_id_recv_x_hi, ierr
     2	   )	
	   Call ErrorHandler(ierr,ERROR_RECV)
	endif
       
	if (iprocx .lt. XPROCS) then
	   xbuffer_send(:)=psi(imax,1:jmax)

	   call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1	        right, MSG_XCH_XHI_TAG, comm3D, msg_id_send_x_hi, ierr)
	   Call ErrorHandler(ierr,ERROR_SEND)
	endif
       
	if (iprocx .gt. 1) then

	   call MPI_Wait(msg_id_recv_x_hi, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)

	   psi(0,1:jmax)=xbuffer_recv(:)
	endif

	if (iprocx .lt. XPROCS) then
	   call MPI_Wait(msg_id_send_x_hi, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)
	endif

C	update x-high boundaries

	if (iprocx .lt. XPROCS) then
c	   isource = iprocnum (iprocx+1, iprocy, iprocz)-1	   
	   call MPI_Irecv(xbuffer_recv, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1	   right, MSG_XCH_XLOW_TAG, comm3D, msg_id_recv_x_low,
     2	   ierr) 
	   Call ErrorHandler(ierr,ERROR_RECV)
	endif

	if (iprocx .gt. 1) then
	   xbuffer_send(:)=psi(1,1:jmax)

	   call MPI_Isend(xbuffer_send, XBUFFSIZE, MPI_DOUBLE_PRECISION,
     1	        left, MSG_XCH_XLOW_TAG, comm3D, msg_id_send_x_low, ierr
     2	        )
	   Call ErrorHandler(ierr,ERROR_SEND)
	endif

	if (iprocx .lt. XPROCS) then

	   call MPI_Wait(msg_id_recv_x_low, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)

	   psi(imax+1,1:jmax)=xbuffer_recv(:)
	endif

	if (iprocx .gt. 1) then
	   call MPI_Wait(msg_id_send_x_low, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)
	endif

C       -------Y DIRECTION COMMUNICATION

C	update y-low boundaries
c	ybufsize=(ihi-ilo+1)*(IZHI-IZLO+1)*REAL_SIZE*nfields*NGHOST

	if (iprocy .gt. 1) then
	   call  MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION
     1	         , bottom, MSG_XCH_YHI_TAG, comm3D, msg_id_recv_y_hi,
     2	         ierr)
	   Call ErrorHandler(ierr,ERROR_RECV)
	endif
       
	if (iprocy .lt. YPROCS) then
	
	   ybuffer_send(1:imax)=psi(1:imax,jmax)

	   call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION
     1	        , top, MSG_XCH_YHI_TAG, comm3D, msg_id_send_y_hi, ierr)

	   Call ErrorHandler(ierr,ERROR_SEND)
	endif
       
	if (iprocy .gt. 1) then


	   call MPI_Wait(msg_id_recv_y_hi, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)
	   psi(1:imax,0)=ybuffer_recv(1:imax)
	endif
	
	if (iprocy .lt. YPROCS) then

	   call MPI_Wait(msg_id_send_y_hi, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)
	endif

C	update y-high boundaries

	if (iprocy .lt. YPROCS) then
	   call MPI_Irecv(ybuffer_recv, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     1	        top, MSG_XCH_YLOW_TAG, comm3D, msg_id_recv_y_low,
     2	        ierr)
	   Call ErrorHandler(ierr,ERROR_RECV)
	endif
       
	if (iprocy .gt. 1) then
	   ybuffer_send(1:imax)=psi(1:imax,1)

	   call MPI_Isend(ybuffer_send, YBUFFSIZE, MPI_DOUBLE_PRECISION,
     1		bottom, MSG_XCH_YLOW_TAG, comm3D, msg_id_send_y_low,
     2	        ierr)  
	   Call ErrorHandler(ierr,ERROR_SEND)
	endif

	if (iprocy .lt. YPROCS) then
	   call MPI_Wait(msg_id_recv_y_low, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)
	   psi(1:imax,jmax+1)=ybuffer_recv(1:imax)
	endif

	if (iprocy .gt. 1) then
	   call MPI_Wait(msg_id_send_y_low, status, ierr)
	   Call ErrorHandler(ierr,ERROR_WAIT)
	endif
#endif

	return
	end
c
c-----------------------------------------------------------------------
	subroutine GatherPsi(psi,imax,jmax,psig,imaxg,jmaxg)
c-----------------------------------------------------------------
c       Description: All gather of psi - collected when 
c                    lMaxLocal is reached in multigrid
c                    boundaries for psi
c       Input: psi - field with ghost which must be gathered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: psig - gathered field with global mesh extent
c-----------------------------------------------------------------
	use mesh_uparms
	use mesh_common
#ifdef PARALLEL
	use mpistuff
#endif

	implicit none

	integer:: imax, jmax, imaxg, jmaxg
	double precision:: psi(0:imax+1,0:jmax+1)
	double precision:: psig(0:imaxg+1,0:jmaxg+1)

#ifdef PARALLEL
	double precision:: buffer_send(imax,jmax)
	double precision:: buffer_recv(imax,jmax,XPROCS*YPROCS)

	integer:: bufferSize
	integer:: recvBufferSize
	integer:: ii,jj,mm,kk,l,idest
	integer:: iprocnum
	integer:: msg_id_send
	integer:: msg_id_recv
	integer:: ipx,ipy,ipz,coords(3)


	bufferSize=imax*jmax
	buffer_send=psi(1:imax,1:jmax)
	
	call MPI_AllGather(buffer_send,bufferSize,MPI_DOUBLE_PRECISION,
     &   buffer_recv,bufferSize,MPI_DOUBLE_PRECISION,
     &	       comm3D,
     &	        ierr)
c
	do l=1,XPROCS*YPROCS,1
 	call MPI_Cart_Coords(comm3D, l-1, 3, coords, ierr)
        Call ErrorHandler(ierr,ERROR_CARTCOORDS)
 	ipx = coords(1) + 1
 	ipy = coords(2) + 1
 	ipz = coords(3) + 1	   
	psig((ipx-1)*imax+1:ipx*imax,(ipy-1)*jmax+1:ipy*jmax)
     &         =buffer_recv(:,:,l)
	enddo
#else
	psig=psi
#endif
	return
	end
c-----------------------------------------------------------------------
	subroutine GatherRes(res,imax,jmax,resg,imaxg,jmaxg)
c-----------------------------------------------------------------
c       Description: All gather of psi - collected when 
c                    lMaxLocal is reached in multigrid
c                    boundaries for psi
c       Input: res - field with no ghost which must be gathered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: resg - gathered field with global mesh extent
c-----------------------------------------------------------------
	use mesh_uparms
	use mesh_common
#ifdef PARALLEL
	use mpistuff
#endif

	implicit none

	integer:: imax, jmax, imaxg, jmaxg
	double precision:: res(imax,jmax)
	double precision:: resg(imaxg,jmaxg)

#ifdef PARALLEL
	double precision:: buffer_send(imax,jmax)
	double precision:: buffer_recv(imax,jmax,XPROCS*YPROCS)

	integer:: bufferSize
	integer:: ii,jj,mm,kk,l,idest, i,j
	integer:: iprocnum
	integer:: msg_id_send
	integer:: msg_id_recv
	integer:: ipx,ipy,ipz,coords(3)


	bufferSize=imax*jmax
	buffer_send=res(1:imax,1:jmax)
c
	call MPI_AllGather(buffer_send,bufferSize,MPI_DOUBLE_PRECISION,
     &   buffer_recv,bufferSize,MPI_DOUBLE_PRECISION,
     &	       comm3D,
     &	        ierr)
c
	do l=1,XPROCS*YPROCS,1
 	call MPI_Cart_Coords(comm3D, l-1, 3, coords, ierr)
        Call ErrorHandler(ierr,ERROR_CARTCOORDS)
 	ipx = coords(1) + 1
 	ipy = coords(2) + 1
 	ipz = coords(3) + 1	   
c
	resg((ipx-1)*imax+1:ipx*imax,(ipy-1)*jmax+1:ipy*jmax)
     &         =buffer_recv(:,:,l)
c
	enddo
#else
	resg=res
#endif
	return
	end
c-----------------------------------------------------------------------
	subroutine ScatterPsi(psi,imax,jmax,psig,imaxg,jmaxg)
c-----------------------------------------------------------------
c       Description: All scatter of psi when going back to
c                    a level whence distributed multigrid is continued.
c       Input: psig - field with ghost which must be scattered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: psi - local field
c-----------------------------------------------------------------
	use mesh_common
	implicit none

	integer:: imax, jmax, imaxg, jmaxg
	double precision:: psi(0:imax+1,0:jmax+1)
	double precision:: psig(0:imaxg+1,0:jmaxg+1)

	psi(1:imax,1:jmax)=
     &   psig((iprocx-1)*imax+1:iprocx*imax,
     &        (iprocy-1)*jmax+1:iprocy*jmax)
	return
	end
c-----------------------------------------------------------------------
	subroutine ScatterRes(res,imax,jmax,resg,imaxg,jmaxg)
c-----------------------------------------------------------------
c       Description: All scatter of res when going back to
c                    a level whence distributed multigrid is continued.
c       Input: resg - field with no ghost which must be scattered
c              imax,jmax,kmax: local mesh extent 
c              imaxg,jmaxg,kmaxg: global mesh extent 
c       Output: res - local field
c-----------------------------------------------------------------
	use mesh_common
	implicit none

	integer:: imax, jmax, imaxg, jmaxg
	double precision:: res(imax,jmax)
	double precision:: resg(imaxg,jmaxg)

	res(1:imax,1:jmax)=
     &   resg((iprocx-1)*imax+1:iprocx*imax,
     &        (iprocy-1)*jmax+1:iprocy*jmax)
	return
	end
c--------------------------------------------------------------
      subroutine MaxNorm(res,imax,jmax,kmax,resMax)
c-----------------------------------------------------------------
c       Description: Compute max norm of res
c       Input: res - field with no ghost 
c              imax,jmax,kmax: local mesh extent 
c       Output: resMax - max norm
c-----------------------------------------------------------------
	use mesh_uparms
	use mesh_common
#ifdef PARALLEL
	use mpistuff
#endif
	implicit none

	integer:: imax, jmax, kmax
	double precision:: res(imax,jmax)
	double precision:: resMax
	double precision:: resMaxLocal
#ifdef PARALLEL
	resMaxLocal=maxval(dabs(res))
        call MPI_Allreduce(resMaxLocal, resMax, 1, MPI_DOUBLE_PRECISION,
     &       MPI_MAX, comm3D, ierr)
        Call ErrorHandler(ierr,ERROR_ALLREDUCE)
#else
	resMax=maxval(dabs(res))
#endif
	return
	end

