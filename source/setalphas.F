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
c     $Log: setalphas.F,v $
c-----------------------------------------------------------------
      module alphas
      contains
      subroutine SetAlphas(utilde,alpha,durl,evl,evr,ilo,ihi)
c-----------------------------------------------------------------
c     Description: Calculation of [evl]durl
c     Input: utilde - primitive variables at cell faces (not used)
c            durl - = vr-vl
c            evl, evr -left/right eigenvectors. evr is not used
c     Output: alpha = [evl]durl
c     
c-----------------------------------------------------------------
      use mesh
      use mesh_common

      integer:: ilo,ihi
      double precision:: utilde(INLO:INHI,nvar)
      double precision:: alpha(INLO:INHI,nvar)
      double precision:: durl(INLO:INHI,nvar)
      double precision:: evl(INLO:INHI,nvar,nvar)
      double precision:: evr(INLO:INHI,nvar,nvar)
c     
c     Local variables
      integer:: i,j,k,l,m
c     
      do i=ilo,ihi,1
         alpha(i,:)=matmul(evl(i,:,:),durl(i,:))
      enddo
c
      return
      end subroutine SetAlphas
      end module alphas




