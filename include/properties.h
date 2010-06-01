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
c     $Log: properties.h,v $
c=================================================================


      module properties
c-----------------------------------------------------------------
c     Description: contains physical constants for the plasma
c
c     Contains:
c             rgas - (dbl) universal gas constant
c            gamma - (dbl) ratio of specific heats
c             wmol - (dbl) molecular weight
c               mu - (dbl) viscosity coefficient
c              eta - (dbl) permittivity constant
c        etaFactor - (dbl) ??
c            kappa - (dbl) thermal conductivity
c-----------------------------------------------------------------
      save

      double precision :: rgas
      double precision :: gamma
      double precision :: wmol
      double precision :: mu
      double precision :: eta
      double precision :: etaFactor
      double precision :: kappa
c
      double precision :: lref, tauRef, velRef
      double precision :: Bref, pref, n0ref, rhoRef, Tref, TrefeV
      double precision :: Mi, mu0, K2eV, RUniversal
c     
      end module properties
c=================================================================
