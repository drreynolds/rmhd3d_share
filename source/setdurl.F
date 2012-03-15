C -*- Mode: Fortran; -*-
c-----------------------------------------------------------------------
        module udiff
        contains
        subroutine SetDURL(durl,vl,vr,ilo,ihi)
        use mesh
        use mesh_common
        use properties

        double precision:: durl(INLO:INHI,nvar)
	double precision:: vl(INLO:INHI,nvar)
	double precision:: vr(INLO:INHI,nvar)
        integer:: i,j,k,l,m,ilo,ihi
        logical:: Mask(INLO:INHI)
c
c
        Mask=.false.
        Mask(ilo:ihi)=.true.
c
        where(mask)
           durl(:,1)=vr(:,1)-vl(:,1)
           durl(:,2)=vr(:,1)*vr(:,2)-
     &          vl(:,1)*vl(:,2)
           durl(:,3)=vr(:,1)*vr(:,3)-
     &          vl(:,1)*vl(:,3)
           durl(:,4)=vr(:,1)*vr(:,4)-
     &          vl(:,1)*vl(:,4)
           durl(:,5)=vr(:,5)-vl(:,5)
           durl(:,6)=vr(:,6)-vl(:,6)
           durl(:,7)=vr(:,7)-vl(:,7)
        endwhere
c
        where(mask)
           durl(:,8)=vl(:,8)/(gamma-1.D0)+0.5D0*vl(:,1)
     &          *(vl(:,2)**2+vl(:,3)**2+
     &          vl(:,4)**2)+0.5D0*(vl(:,5)**2+
     &          vl(:,6)**2+vl(:,7)**2)
        endwhere
        where(mask)
           durl(:,8)=vr(:,8)/(gamma-1.D0)+0.5D0*vr(:,1)
     &          *(vr(:,2)**2+vr(:,3)**2+
     &          vr(:,4)**2)+0.5D0*(vr(:,5)**2+
     &          vr(:,6)**2+vr(:,7)**2)
     &          -durl(:,8)
        endwhere
	
        return
        end subroutine SetDURL
        end module udiff
c-----------------------------------------------------------------------
