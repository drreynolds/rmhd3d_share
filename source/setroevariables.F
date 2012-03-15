
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
c     $Log: setroevariables.F,v $
c-----------------------------------------------------------------
      subroutine SetRoeVariables(utilde,vl,vr,ilo,ihi)
c-----------------------------------------------------------------------
c     Description: compute average cell face primitive variables
c                  using simple arithmetic averages. This violates
c                  Roe's conditions but not all of Roe's conditions
c                  can be satisfied in MHD.
c     Input: vl,vr (left/right states)
c     Output: utilde - average of vl and vr.
c-----------------------------------------------------------------------
        use mesh
        double precision:: vl(INLO:INHI,nvar)
        double precision:: vr(INLO:INHI,nvar)
        double precision:: utilde(INLO:INHI,nvar)
        double precision:: rhol,rhor
        integer:: ilo,ihi
        integer:: i,j,k,l
        integer:: ierror

        utilde=0.D0
c     Simple arithmetic average
        do i=ilo,ihi,1
c           rhol=sqrt(vl(i,1))
c           rhor=sqrt(vr(i,1))
           do l=1,nvar,1
             utilde(i,l)=0.5D0*(vl(i,l)+ vr(i,l))
c              utilde(i,l)=(vl(i,l)*rhol+ vr(i,l)*rhor)/(rhol+rhor)
c              utilde(i,l)=vl(i,l)
           enddo
        enddo	
	return
	end
