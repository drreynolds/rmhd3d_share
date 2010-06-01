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
c     $Log: modules.F,v $
c-----------------------------------------------------------------
c     Description: chooses dynamic/static grid modules, inputs 
c        problem properties, etc.
c-----------------------------------------------------------------

c======= General Inclusions ===========
#ifdef DYNAMIC
	include 'mesh_uparms.h'
	include 'mesh_parms.h'
#else
	include 'mesh_uparms_static.h'
	include 'mesh_parms_static.h'
#endif
	include 'mesh_common.h'
	include 'mesh.h'
	include 'properties.h'
	include 'mgparams.h'
	include 'boundaries.h'

#ifdef LONGINT8
	module sundials_fcmix
	save
	integer, parameter :: sund_long_int=8
	end module sundials_fcmix
#else
	module sundials_fcmix
	save
	integer, parameter :: sund_long_int=4
	end module sundials_fcmix
#endif


#ifdef PARALLEL
	include 'mpistuff.h'
#endif


	module iounits
c-----------------------------------------------------------------
c     Description: maps output streams 
c        istderr, istdin, istdout, ihis, imesh, ifluid, idump,
c        ibinoutput, itecoutput, itime, idiag, idiagCons, idiagEng
c-----------------------------------------------------------------

        integer, parameter :: istderr    = 0
        integer, parameter :: istdin     = 5
        integer, parameter :: istdout    = 6
        integer, parameter :: ihis       = 8
        integer, parameter :: imesh      = 16
        integer, parameter :: ifluid     = 21
        integer, parameter :: idump      = 71
        integer, parameter :: ibinoutput = 41
        integer, parameter :: itecoutput = 61
        integer, parameter :: itime      = 81
        integer, parameter :: idiag      = 91
        integer, parameter :: idiagCons  = 92
        integer, parameter :: idiagEng   = 94
        integer, parameter :: idiagflux  = 95

	end module iounits
c=================================================================


	
	module TVDCOEFFS
c-----------------------------------------------------------------
c       Description: holds TVD coefficients 
c       tvdCoeff1(1:3), tvdCoeff2(1:3), tvdCoeff3(1:3)
c-----------------------------------------------------------------
	save
	double precision:: tvdCoeff1(3)
	double precision:: tvdCoeff2(3)
	double precision:: tvdCoeff3(3)
	end module TVDCOEFFS
c-----------------------------------------------------------------
        subroutine SetTVDCoeffs
        use TVDCOEFFS
        tvdCoeff1(1)=1.D0
        tvdCoeff1(2)=0.D0
        tvdCoeff1(3)=1.D0
        tvdCoeff2(1)=0.75D0
        tvdCoeff2(2)=0.25D0
        tvdCoeff2(3)=0.25D0
        tvdCoeff3(1)=1.D0/3.D0
        tvdCoeff3(2)=2.D0/3.D0
        tvdCoeff3(3)=2.D0/3.D0
        return
        end subroutine SetTVDCoeffs
c=================================================================
