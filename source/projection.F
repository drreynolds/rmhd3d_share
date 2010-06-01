C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------
c     Ravi Samtaney
c     Copyright 2003
c     Princeton Plasma Physics Laboratory
c     All Rights Reserved
c
c     Daniel R. Reynolds
c     UC San Diego, Mathematics
c
c-----------------------------------------------------------------
c     Multigrid projection
c-----------------------------------------------------------------------
       subroutine Projection(phi,divb)
c-----------------------------------------------------------------
c     Input: divb - cell centered divergence calculated using face
c     centered values of the magnetic field (from a Riemann solve)
c     Output: phi -auxillary variable used for projection
c     Solves the Poisson equation L(phi)=divb where L is the
c     Laplacian. Boundary conditions on phi are Neumann.
c     Limitations: For now implemented in 2D.
c     
       use mesh
       use properties
#ifdef PARALLEL
      use mpistuff
#endif

c
c=======Declarations=========
c
	double precision:: divb(NXLSIZE,NYLSIZE,NZLSIZE)
	double precision:: res(NXLSIZE,NYLSIZE)
	double precision:: phi(IXLO:IXHI,IYLO:IYHI)
        integer:: i,j,k
c
        res(:,:)=divb(:,:,1)
c        res=0.D0
c        res(32,16)=1.D0
c        res(34,17)=-1.D0
        call Multigrid(phi,res)
c     
c$$$        do k=1,nzlocal,1
c$$$           do j=1,nylocal,1
c$$$              do i=1,nxlocal,1
c$$$                 write(6,*) 'PHI',i,j,phi(i,j),phi(i+1,j),phi(i,j+1)
c$$$                 write(6,*) 'RES',i,j,res(i,j)
c$$$              enddo
c$$$           enddo
c$$$        enddo
c
        return
        end
c-----------------------------------------------------------------------
       subroutine CorrectFlux(finv,phi,vface,idir)
c-----------------------------------------------------------------
c     
c     Input: phi (auxially variable used for projection)
c            vface (face centered primitive variable (after RP solve)
c            idir - direction(=1,2,3)
c     Output: finv - fluxes corrected post B projection
c     Issues: Do we also change the components of B which are not normal 
c             to the face in question? 
       use mesh
       use properties
#ifdef PARALLEL
      use mpistuff
#endif

c
c=======Declarations=========
c
	double precision:: finv(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
	double precision:: vface(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
	double precision:: phi(IXLO:IXHI,IYLO:IYHI)
c
        double precision:: bi,bj,bk
        double precision:: bio,bjo,bko
        double precision:: uu,vv,ww
c
        integer:: indx,indy,indz
c
        indx=0
        indy=0
        indz=0
        if(idir.eq.1) indx=1
        if(idir.eq.2) indy=1
        if(idir.eq.3) indz=1
c     Correct flux based on projected face Magnetic field
        do k=1-indz,nzlocal,1
           do j=1-indy,nylocal,1
              do i=1-indx,nxlocal,1
                 uu=vface(i,j,k,2)
                 vv=vface(i,j,k,3)
                 ww=vface(i,j,k,4)
                 bio=vface(i,j,k,5)
                 bjo=vface(i,j,k,6)
                 bko=vface(i,j,k,7)
                 bi=bio-(phi(i+indx,j+indy)-phi(i,j))/dx
c                 write(6,*) idir, i,j,'bi-bio',bi,bio,bi-bio
c                 write(6,*) 'phi',idir, i,j,phi(i+indx,j+indy),phi(i,j)
                 bj=bjo
                 bk=bko
                 finv(i,j,k,2)= finv(i,j,k,2)+
     &                0.5D0*(bj*bj+bk*bk-bjo*bjo-bko*bko
     &                +bio*bio-bi*bi)
                 finv(i,j,k,3)= finv(i,j,k,3)+
     &                bio*bjo-bi*bj
                 finv(i,j,k,4)= finv(i,j,k,4)+
     &                bio*bko-bi*bk
                 finv(i,j,k,6)= finv(i,j,k,6)+
     &                uu*(bj-bjo)+
     &                vv*bio-vv*bi
                 finv(i,j,k,7)= finv(i,j,k,7)+
     &                uu*(bk-bko)+
     &                ww*bio-ww*bi
                 finv(i,j,k,8) = finv(i,j,k,8)
     &                -uu*(bio*bio+bjo*bjo+bko*bko
     &                -bi*bi-bj*bj-bk*bk)
     &                +bio*(uu*bio+vv*bjo+ww*bko)
     &                -bi*(uu*bi+vv*bj+ww*bk)
              enddo
           enddo
        enddo
c        
        return
        end
c-----------------------------------------------------------------------
       subroutine CorrectMagneticField(ux,phi)
c-----------------------------------------------------------------
c     
c     Input: phi (auxially variable used for projection)
c     Output: ux (cell centered conserved variable) - only the 
c             magnetic field components are changed
       use mesh
       use properties
#ifdef PARALLEL
      use mpistuff
#endif

c
c=======Declarations=========
c
	double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
	double precision:: phi(IXLO:IXHI,IYLO:IYHI)
	double precision:: Bxf(IXLO:IXHI,IYLO:IYHI)
	double precision:: Byf(IXLO:IXHI,IYLO:IYHI)
c
        double precision:: bi,bj,bk
        double precision:: bio,bjo,bko
        double precision:: uu,vv,ww
c
c
        do k=1,nzlocal,1
           do j=1,nylocal,1
              do i=1,nxlocal,1
                 bio=ux(i,j,k,5)
                 bjo=ux(i,j,k,6)
c                 ux(i,j,k,5)=bio
c                 ux(i,j,k,6)=bjo
                 ux(i,j,k,5)=bio-0.5D0*(phi(i+1,j)-phi(i-1,j))/dx
                 ux(i,j,k,6)=bjo-0.5D0*(phi(i,j+1)-phi(i,j-1))/dx
              enddo
           enddo
        enddo
        return
        do k=1,nzlocal,1
           do j=1,nylocal,1
              do i=0,nxlocal,1
                 bio=(ux(i,j,k,5)+ux(i+1,j,k,5))/2.D0
                 Bxf(i,j)=bio
c-(phi(i+1,j)-phi(i,j))/dx
              enddo
           enddo
        enddo
        do k=1,nzlocal,1
           do j=0,nylocal,1
              do i=1,nxlocal,1
                 bjo=(ux(i,j,k,6)+ux(i,j+1,k,6))/2.D0
                 Byf(i,j)=bjo
c-(phi(i,j+1)-phi(i,j))/dx
              enddo
           enddo
        enddo
        do k=1,nzlocal,1
           do j=1,nylocal,1
              do i=1,nxlocal,1
                 ux(i,j,k,5)=0.5D0*(Bxf(i,j)+Bxf(i-1,j))
                 ux(i,j,k,6)=0.5D0*(Byf(i,j)+Byf(i,j-1))
              enddo
           enddo
        enddo
c        
        return
        end
c
c-----------------------------------------------------------------------
       subroutine Multigrid(phi,res)
c-----------------------------------------------------------------
c     Input: res - RHS of the equations L(phi)=res (=divb)
c     Implicit Input: lMax and lMaxLocal which are determined 
c                     outside of this subroutine, are the
c                     maximum no of global and local levels, resp.
c     Output: phi
c     Multigrid implementation in 2D 
c     Author: Ravi Samtaney
c     Assumed parameters: No of presmooth and postsmooth iterations 
c                         set to 10.
c                         Convergence tolerance set to 10^-8
c                         Max no. of V-cycles =10. 
c    
c     
c     
       use mesh
       use properties
       use mesh_common
#ifdef PARALLEL
      use mpistuff
#endif
      use MGParams

c
c=======Declarations=========
	double precision:: ux(IXLO:IXHI,IYLO:IYHI, IZLO:IZHI,nvar)
        double precision:: dt
c
	double precision:: phi(IXLO:IXHI,IYLO:IYHI)
c
	double precision,allocatable:: psi(:,:)
	double precision,allocatable:: psi0(:,:)
	double precision:: res(nxlocal,nylocal)
	double precision:: resx(nxlocal,nylocal)
c
	integer:: indx,indy,indz
	integer:: i,j,k,l
        integer:: maxOuterIter
        integer:: iter
c
        double precision:: resMax
        double precision:: tolerance
c
        allocate(psi(0:nxlocal+1,0:nylocal+1))
        allocate(psi0(0:nxlocal+1,0:nylocal+1))
        maxOuterIter=10
        maxGrid=5
cc
        psi=0.D0
        psi0=0.D0
        phi=0.D0
c        res=0.D0
c        res(10,10)=1.D0
c        res(11,11)=-1.D0
c
        bcHFlag=0
        xChangeFlag=1
        call BC(psi,nxlocal,nylocal,nzlocal,1)
        call BC(psi0,nxlocal,nylocal,nzlocal,1)
        resMax=1.D16
        tolerance=1.D-08
        iter=1
        do while(iter.le.maxOuterIter.and.resMax.gt.tolerance)
           call MGSolve(psi,psi0,res,nxlocal,nylocal,nzlocal,0,dx)
           psi0=psi
           resx=0.D0
           call BC(psi,nxlocal,nylocal,nzlocal,1)
           call BC(psi0,nxlocal,nylocal,nzlocal,1)
           call Residual(psi,res,resx,nxlocal,nylocal,nzlocal,dx)
           call MaxNorm(resx,nxlocal,nylocal,nzlocal,resMax)
           if(iprocx.eq.1.and.iprocy.eq.1.and.iprocz.eq.1) then
              write(6,*) 'Max Norm Residual=',iter,resMax
           endif
           iter=iter+1
        enddo
        phi(0:nxlocal+1,0:nylocal+1)=
     &       psi(0:nxlocal+1,0:nylocal+1)
c        stop
	return
	end
c     

c     
c-----------------------------------------------------------------------
      recursive subroutine MGSolve(psi,psiT,res,imax,jmax,kmax,level,dx)
      use MGParams
      use mesh_uparms
      use mesh_common
      integer:: imax,jmax,kmax,level
      double precision:: psi(0:imax+1,0:jmax+1)
      double precision:: psiT(0:imax+1,0:jmax+1)
      double precision:: res(imax,jmax)
      double precision:: resCorr(imax,jmax)
      double precision, allocatable:: psiCoarse(:,:)
      double precision, allocatable:: psiTCoarse(:,:)
      double precision, allocatable:: resCorrCoarse(:,:)
      double precision, allocatable:: resg(:,:)
      double precision, allocatable:: psig(:,:), psiTg(:,:)
      double precision:: dx
      integer:: maxPreSmoothIterations
      integer:: maxPostSmoothIterations
      integer:: gammaK
      integer:: iterGamma
      integer:: imaxc,jmaxc,kmaxc
      maxPreSmoothIterations=10
      maxPostSmoothIterations=10
      gammaK=1
      
      
c      write(6,*) 'MGSOLVE: level=',level,dx, imax,jmax,xChangeFlag
      if(level.eq.lMax) then
         if(level.eq.0) bcHFlag=0
         call BottomSolve(psiT,psi,res,imax,jmax,kmax,level,dx)
      else
         call Smooth(psiT,psi,res,maxPreSmoothIterations,
     &        imax,jmax,kmax,level,dx)
         call Residual(psi,res,resCorr,imax,jmax,kmax,dx)
         if(level.eq.lMaxLocal) then
            imaxg=(xprocs*imax)
            jmaxg=(yprocs*jmax)
            kmaxg=(zprocs*kmax)
            imaxc=imaxg/2
            jmaxc=jmaxg/2
            kmaxc=kmaxg/2
            allocate(psig(0:imaxg+1,0:jmaxg+1))
            allocate(psiTg(0:imaxg+1,0:jmaxg+1))
            allocate(resg(imaxg,jmaxg))
            allocate(resCorrCoarse(imaxc,jmaxc))
            allocate(psiCoarse(0:imaxc+1,0:jmaxc+1))
            allocate(psiTCoarse(0:imaxc+1,0:jmaxc+1))
            call GatherRes(resCorr,imax,jmax,resg,imaxg,jmaxg)
            call GatherPsi(psi,imax,jmax,psig,imaxg,jmaxg)
            call GatherPsi(psiT,imax,jmax,psiTg,imaxg,jmaxg)
            xChangeFlag=0
         else
            imaxg=(imax)
            jmaxg=(jmax)
            kmaxg=(kmax)
            imaxc=imaxg/2
            jmaxc=jmaxg/2
            kmaxc=kmaxg/2
            allocate(psig(0:imaxg+1,0:jmaxg+1))
            allocate(psiTg(0:imaxg+1,0:jmaxg+1))
            allocate(resg(imaxg,jmaxg))
            allocate(resCorrCoarse(imaxc,jmaxc))
            allocate(psiCoarse(0:imaxc+1,0:jmaxc+1))
            allocate(psiTCoarse(0:imaxc+1,0:jmaxc+1))
            resg=resCorr
            psig=psi
            psiTg=psiT
         endif
         call Restrict(resg,resCorrCoarse,
     &           imaxg,jmaxg,kmaxg)
         psiTCoarse=0.D0
         psiCoarse=0.D0
         bcHFlag=1
         do iterGamma=1,gammaK,1
               call MGSolve(psiCoarse,psiTCoarse,resCorrCoarse,
     &              imaxc,jmaxc,kmaxc,
     &              level+1,dx*2)
               psiTCoarse=psiCoarse
         enddo
         resCorr=0.D0
         call Prolong(resg,imaxg,jmaxg,kmaxg, psiCoarse)
         if(level.eq.lMaxLocal) then
            call ScatterRes(resCorr,imax,jmax,resg,imaxg,jmaxg)
            call ScatterPsi(psi,imax,jmax,psig,imaxg,jmaxg)
            xChangeFlag=1
         else
            resCorr=resg
            psi=psig
            psiT=psiTg
         endif
         psi(1:imax,1:jmax)=psi(1:imax,1:jmax)+
     &        resCorr(1:imax,1:jmax)
         bcHFlag=1
         psiT=psi
         if(level.eq.0) bcHFlag=0
         call Smooth(psi,psi,res,maxPostSmoothIterations,
     &        imax,jmax,kmax,level,dx)
      endif
      return
      end
c-----------------------------------------------------------------------
       subroutine DetermineMaxLevel(lMax,nx,ny)
c-----------------------------------------------------------------
c      Input: Domain size nx x ny
c      Output: Max no. of levels possible if each level is coarsened
c              by a factor of 2 from the previous finer level.
       integer:: lMax, nx,ny
       integer:: m,n
       integer:: i
c
       m=nx
       n=ny
       lMax=0
       do i=1,64,1
          if(mod(m,2).ne.0. or. mod(n,2).ne.0) then
             return
          else
             lmax=lmax+1
             m=m/2
             n=n/2
          endif
       enddo
       return
       end
c-----------------------------------------------------------------------
      subroutine Smooth(psiT,psi,res,maxIterations,
     &        imax,jmax,kmax,level,dx)
c-----------------------------------------------------------------------
c     Input: psiT - psi_tilde
c            res  - residual
c            maxIterations - max no. of smoothing iterations
c            imax,jmax,kmax : mesh size
c            level - multi-grid level (not used here)
c            dx - mesh spacing
c     Implicit input: xChangeFlag  =1,0. Do (not) Exchange ghost 
c                     information at this level if this is 1 (0).
c                     bchflag = 0,1. Use (do not use)  homogeneous 
c                     boundary conditions if this is  1 (0)
c     Output: psi - after smoothing of psiT
c     Description: Uses Gauss-Siedel red-black ordering
c-----------------------------------------------------------------------
      use mesh_uparms
      use mesh_common
      use MGParams
      double precision:: psi(0:imax+1,0:jmax+1)
      double precision:: psiT(0:imax+1,0:jmax+1)
      double precision:: res(imax,jmax)
      integer:: i,j,k
      double precision:: lamda,dxsq
      double precision:: dx
c
      dxsq=dx*dx
      lamda=1.D0/dxsq
         
c
      do iter=1,maxIterations,1
      if (bcHFlag.eq.0) then
         call BC(psiT,imax,jmax,kmax,xChangeFlag)
      else
         call BCHomogeneous(psiT,imax,jmax,kmax,xChangeFlag)
      endif
      Call ExchangePsi(psiT,imax,jmax,xChangeFlag)
c     Red-Black ordering
      do j=1,jmax,1
         do i=mod(j+1,2)+1,imax,2
            psi(i,j)=psiT(i,j)+0.25D0*(
     &           psiT(i+1,j)+psiT(i-1,j)+psiT(i,j-1)+psiT(i,j+1)
     &           -4.D0*psiT(i,j))-0.25D0*dxsq*res(i,j)
         enddo
      enddo
c
      do j=1,jmax,1
         do i=mod(j,2)+1,imax,2
            psi(i,j)=psiT(i,j)+0.25D0*(
     &           psiT(i+1,j)+psiT(i-1,j)+psiT(i,j-1)+psiT(i,j+1)
     &           -4.D0*psiT(i,j))-0.25D0*dxsq*res(i,j)
         enddo
      enddo
c
      do j=1,jmax,1
         do i=1,imax,1
            psiT(i,j)=psi(i,j)
         enddo
      enddo
      if (bcHFlag.eq.0) then
         call BC(psi,imax,jmax,kmax,xChangeFlag)
      else
         call BCHomogeneous(psi,imax,jmax,kmax,xChangeFlag)
      endif
      Call ExchangePsi(psi,imax,jmax,xChangeFlag)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine BottomSolve(psiT,psi,res,
     &        imax,jmax,kmax,level,dx)
c-----------------------------------------------------------------------
c     Description: Used when mesh cannot be coarsened anymore. 
c                  Still uses Gauss-Siedel Red-Black for
c                  100 iterations. 
c     Issues: 100 is hard -coded but seems to work for the reconnection
c             problem. A better bottom solver (e.g. biconjugate gradient)
c             would be desirable.
c     Input: psiT - psi_tilde
c            res - residual
c            imax,jmax,kmax: local mesh size
c            level - multi-grid level
c            dx - local mesh spacing
c     Output: psi
c-----------------------------------------------------------------------
      double precision:: psi(0:imax+1,0:jmax+1)
      double precision:: psiT(0:imax+1,0:jmax+1)
      double precision:: res(imax,jmax)
      double precision:: dx
      integer:: i,j,k
c
      call Smooth(psiT,psi,res,100,
     &        imax,jmax,kmax,level,dx)
      return
      end

c-----------------------------------------------------------------------
      Subroutine Residual(psi,res,resCorr,imax,jmax,kmax,dx)
c-----------------------------------------------------------------------
c     Description: Calculate the residual Residual=L(phi)- RHS
c     Input: psi 
c            res - rhs 
c            imax,jmax,kmax: local mesh size
c            dx - local mesh spacing
c     Output: resCorr  evualated as resCorr=L(psi)-res
c-----------------------------------------------------------------------
      integer:: imax,jmax,kmax
      double precision:: psi(0:imax+1,0:jmax+1)
      double precision:: resCorr(imax,jmax)
      double precision:: res(imax,jmax)
      double precision:: dx
c
      double precision:: lamda
      integer:: i,j
c
      lamda=1.D0/dx**2
c
      do j=1,jmax,1
         do i=1,imax,1
            resCorr(i,j)=res(i,j)-
     &           lamda*(psi(i+1,j)+psi(i-1,j)+psi(i,j-1)+psi(i,j+1)-
     &           4.D0*psi(i,j))
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      Subroutine Restrict(psi,psiCoarse,imax,jmax,kmax)
c-----------------------------------------------------------------------
c     Description: Restrict the solution onto the next coarser level
c     Input: psi - solution at current level
c            imax,jmax,kmax: mesh size at the current level
c     Output: phiCoarse: solution on the next coarser level
c-----------------------------------------------------------------------
      integer:: imax,jmax,kmax
      double precision:: psi(imax,jmax)
      double precision:: psiCoarse(imax/2,jmax/2)
c
      integer:: i,j
      integer:: ii,jj
c
      do j=1,jmax/2,1
         jj=(j-1)*2+1
         do i=1,imax/2,1
            ii=(i-1)*2+1
            psiCoarse(i,j)=0.25D0*(
     &           psi(ii,jj)+psi(ii+1,jj)+psi(ii+1,jj+1)+psi(ii,jj+1))
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine  Prolong(psi,imax,jmax,kmax, psiCoarse)
c-----------------------------------------------------------------------
c     Description: Prolong the solution onto the next finer level
c     Input: psiCoarse - solution at current level
c            imax,jmax,kmax: mesh size at the next finer level
c     Output: phi: solution prolonged on the next finer level
c     Issues:  Is this second order accurate? Do we care, because the 
c              entire MG solver is converging well?
c-----------------------------------------------------------------------
      integer:: imax,jmax,kmax
      double precision:: psi(imax,jmax)
      double precision:: psiCoarse(0:imax/2+1,0:jmax/2+1)
      integer:: i,j
      integer:: ii,jj
      do j=1,jmax/2,1
         jj=(j-1)*2+1
         do i=1,imax/2,1
            ii=(i-1)*2+1
            psi(ii,jj)=psiCoarse(i,j)
            psi(ii+1,jj)=psiCoarse(i,j)
            psi(ii+1,jj+1)=psiCoarse(i,j)
            psi(ii,jj+1)=psiCoarse(i,j)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine BCHomogeneous(psi,imax,jmax,kmax,xChangeFlag)
c-----------------------------------------------------------------------
c     Description: Homogeneous Neumann BCs for psi
c     Input: psi - solution
c            imax,jmax,kmax: mesh size at this level
c            xChangeFlag - whether we are doing this on 1 or many procs.
c     Output: psi: Modified at the ghost cells at physical boundaries
c     
c-----------------------------------------------------------------------
      use mesh_uparms
      use mesh_common
      integer:: imax,jmax,kmax
      integer:: xChangeFlag
      double precision:: psi(0:imax+1,0:jmax+1)
c
      psi(0,1:jmax)=psi(imax,1:jmax)
      psi(imax+1,1:jmax)=psi(1,1:jmax)
      psi(1:imax,0)=psi(1:imax,1)
      psi(1:imax,jmax+1)=psi(1:imax,jmax)
c$$$      psi(0,0:jmax+1)=psi(imax,0:jmax+1)
c$$$      psi(imax+1,0:jmax+1)=psi(1,0:jmax+1)
c$$$      psi(0:imax+1,0)=psi(0:imax+1,jmax)
c$$$      psi(0:imax+1,jmax+1)=psi(0:imax+1,1)
c      psi(0,0:jmax+1)=0.D0-psi(1,0:jmax+1)
c      psi(imax+1,0:jmax+1)=0.D0
c      psi(0:imax+1,0)=0.D0
c      psi(0:imax+1,jmax+1)=0.D0
c$$$      if(xChangeFlag.eq.0) then
c$$$         psi(0,0:jmax+1)=psi(1,0:jmax+1)
c$$$         psi(imax+1,0:jmax+1)=psi(imax,0:jmax+1)
c$$$         psi(0:imax+1,0)=psi(0:imax+1,1)
c$$$         psi(0:imax+1,jmax+1)=psi(0:imax+1,jmax)
c$$$      else
c$$$         if(iprocx .eq. 1) then
c$$$            psi(0,0:jmax+1)=psi(1,0:jmax+1)
c$$$         endif
c$$$         if(iprocx .eq. XPROCS) then
c$$$            psi(imax+1,0:jmax+1)=psi(imax,0:jmax+1)
c$$$         endif
c$$$         if(iprocy .eq. 1) then
c$$$            psi(0:imax+1,0)=psi(0:imax+1,1)
c$$$         endif
c$$$         if(iprocy .eq. YPROCS) then
c$$$            psi(0:imax+1,jmax+1)=psi(0:imax+1,jmax)
c$$$         endif
c$$$      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine BC(psi,imax,jmax,kmax,xChangeFlag)
c-----------------------------------------------------------------------
c     Description: Neumann BCs for psi - still homogeneous but it is 
c                  important to make the distinction with the 
c                  BCHomogeneous routine in case we have differenct BCs.
c     Input: psi - solution
c            imax,jmax,kmax: mesh size at this level
c            xChangeFlag - whether we are doing this on 1 or many procs.
c     Output: psi: Modified at the ghost cells at physical boundaries
c     
c-----------------------------------------------------------------------
      use mesh_uparms
      use mesh_common
      integer:: imax,jmax,kmax
      integer:: xChangeFlag
      double precision:: psi(0:imax+1,0:jmax+1)
c
      psi(0,1:jmax)=psi(imax,1:jmax)
      psi(imax+1,1:jmax)=psi(1,1:jmax)
      psi(1:imax,0)=psi(1:imax,1)
      psi(1:imax,jmax+1)=psi(1:imax,jmax)
c      psi(0,0:jmax+1)=2.D0-psi(1,0:jmax+1)
c      psi(imax+1,0:jmax+1)=0.D0
c      psi(0:imax+1,0)=0.D0
c      psi(0:imax+1,jmax+1)=0.D0
c$$$      if(xChangeFlag.eq.0) then
c$$$         psi(0,0:jmax+1)=psi(1,0:jmax+1)
c$$$         psi(imax+1,0:jmax+1)=psi(imax,0:jmax+1)
c$$$         psi(0:imax+1,0)=psi(0:imax+1,1)
c$$$         psi(0:imax+1,jmax+1)=psi(0:imax+1,jmax)
c$$$      else
c$$$         if(iprocx .eq. 1) then
c$$$            psi(0,0:jmax+1)=psi(1,0:jmax+1)
c$$$         endif
c$$$         if(iprocx .eq. XPROCS) then
c$$$            psi(imax+1,0:jmax+1)=psi(imax,0:jmax+1)
c$$$         endif
c$$$         if(iprocy .eq. 1) then
c$$$            psi(0:imax+1,0)=psi(0:imax+1,1)
c$$$         endif
c$$$         if(iprocy .eq. YPROCS) then
c$$$            psi(0:imax+1,jmax+1)=psi(0:imax+1,jmax)
c$$$         endif
c$$$      endif
      return
      end

         
