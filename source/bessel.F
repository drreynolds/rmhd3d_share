
      subroutine besj ( x, alpha, n, y, nz )      
c*****************************************************************************80
c
cc BESJ computes a sequence of J Bessel functions of increasing order.
c
c  Discussion:
c
c    BESJ computes an N member sequence of J Bessel functions
c
c      J(ALPHA+K-1) (X)
c
c    for K=1,..,N for non-negative order ALPHA and X.
c
c    A combination of the power series, the asymptotic expansion for X
c    to infinity and the uniform asymptotic expansion for NU to infinity
c    are applied over subdivisions of the (NU,X) plane.  For values of
c    (NU,X) not covered by one of these formulae, the order is
c    incremented or decremented by integer values into a region where
c    one of the formulas apply.
c
c    Backward recursion is applied to reduce orders by integer values
c    except where the entire sequence lies in the oscillatory region.
c    In this case forward recursion is stable and values from the
c    asymptotic expansion for X to infinity start the recursion when it
c    is efficient to do so.
c
c    Leading terms of the series and uniform expansion are tested for
c    underflow.  If a sequence is requested and the last member would
c    underflow, the result is set to zero and the next lower order
c    tried, until a member comes on scale or all members are set
c    to zero.
c
c    Overflow cannot occur.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c    Donald Amos, SL Daniel, MK Weston,
c    CDC 6600 subroutines IBESS and JBESS for Bessel functions
c    I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0
c    ACM Transactions on Mathematical Software,
c    Volume 3, pages 76-92, 1977.
c
c    Frank Olver,
c    Tables of Bessel Functions of Moderate or Large Orders,
c    NPL Mathematical Tables, Volume 6,
c    Her Majesty's Stationery Office, London, 1962.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) X, the argument of the Bessel function.
c    X must be nonnegative.
c
c    Input, real ( kind = 8 ) ALPHA, the order of the first member of
c    the sequence.  ALPHA must be at least 0.0.
c
c    Input, integer N, the number of members in the sequence,
c    N must be at least 1.
c
c    Output, real ( kind = 8 ) Y(N), a vector whose first N components contain
c    values for J(ALPHA+K-1)(X), K=1,...,N
c
c    Output, integer NZ, the number of components of Y set to zero
c    due to underflow.
c
c    NZ=0, normal return, computation completed
c
c    NZ /= 0, Y(N-NZ+1) through Y(N) were set to 0.
c
      implicit none
      
      integer n
      
      real ( kind = 8 ) ak
      real ( kind = 8 ) akm
      real ( kind = 8 ) alngam
      real ( kind = 8 ) alpha
      real ( kind = 8 ) ans
      real ( kind = 8 ) ap
      real ( kind = 8 ) arg
      real ( kind = 8 ) coef
      real ( kind = 8 ) d1mach
      real ( kind = 8 ) dalpha
      real ( kind = 8 ) dfn
      real ( kind = 8 ) dtm
      real ( kind = 8 ) earg
      real ( kind = 8 ) elim1
      real ( kind = 8 ) etx
      real ( kind = 8 ) fidal
      real ( kind = 8 ) flgjy
      real ( kind = 8 ) fn
      real ( kind = 8 ) fnf
      real ( kind = 8 ) fni
      real ( kind = 8 ) fnp1
      real ( kind = 8 ) fnu
      real ( kind = 8 ), parameter, dimension ( 2 ) :: 
     &     fnulim = (/ .0D+00, 60.0D+00 /)
      real ( kind = 8 ) gln
      integer i
      integer i1
      integer i1mach
      integer i2
      integer ialp
      integer idalp
      integer iflw
      integer in
      integer, parameter :: inlim = 150
      integer is
      external jairy
      integer k
      integer kk
      integer km
      integer kt
      integer nn
      integer ns
      integer nz
      real ( kind = 8 ), parameter :: pdf = 0.785398163397448D+00
      real ( kind = 8 ), parameter :: pidt = 1.57079632679490D+00
      real ( kind = 8 ), parameter, dimension ( 4 ) :: pp = (/ 
     &     8.72909153935547D+00, 2.65693932265030D-01, 
     &     1.24578576865586D-01, 7.70133747430388D-04 /)
      real ( kind = 8 ) rden
      real ( kind = 8 ) relb
      real ( kind = 8 ), parameter :: rttp = 7.97884560802865D-01
      real ( kind = 8 ), parameter :: rtwo = 1.34839972492648D+00
      real ( kind = 8 ) rtx
      real ( kind = 8 ) rzden
      real ( kind = 8 ) s
      real ( kind = 8 ) sa
      real ( kind = 8 ) sb
      real ( kind = 8 ) sxo2
      real ( kind = 8 ) s1
      real ( kind = 8 ) s2
      real ( kind = 8 ) t
      real ( kind = 8 ) ta
      real ( kind = 8 ) tau
      real ( kind = 8 ) tb
      real ( kind = 8 ) temp(3)
      real ( kind = 8 ) tfn
      real ( kind = 8 ) tm
      real ( kind = 8 ) tol
      real ( kind = 8 ) tolln
      real ( kind = 8 ) trx
      real ( kind = 8 ) tx
      real ( kind = 8 ) t1
      real ( kind = 8 ) t2
      real ( kind = 8 ) wk(7)
      real ( kind = 8 ) x
      real ( kind = 8 ) xo2
      real ( kind = 8 ) xo2l
      real ( kind = 8 ) y(n)

      nz = 0
      kt = 1
c
c  I1MACH(14) replaces I1MACH(11) in a double precision code
c  I1MACH(15) replaces I1MACH(12) in a double precision code
c
      ta = epsilon ( ta )
      tol = max ( ta, 1.0D-15 )
      i1 = i1mach(14) + 1
      i2 = i1mach(15)
      tb = d1mach(5)
      elim1 = 2.303D+00 * ( real ( -i2, kind = 8 ) * tb - 3.0D+00 )
c
c  TOLLN = -ln(tol)
c
      tolln = 2.303D+00 * tb * real ( i1, kind = 8 )
      tolln = min ( tolln, 34.5388D+00 )

      if ( n < 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'BESJ - Fatal error!'
         write ( *, '(a)' ) '  N is less than 1.'
         return
      end if
      
      if ( n == 1 ) then
         kt = 2
      end if
      
      nn = n
      
      if ( x < 0.0D+00 ) then
         call xerror ( 'BESJ - X less than zero.', 2, 1 )
         return
      end if
      
      if ( x == 0.0D+00 ) then
         
         if ( alpha < 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'BESJ - Fatal error!'
            write ( *, '(a)' ) '  ALPHA less than zero.'
            return
         end if
         
         if ( alpha == 0.0D+00 ) then
            
            y(1) = 1.0D+00
            
            if ( n == 1 ) then
               return
            end if
            
            i1 = 2
            
         else
            
            i1 = 1
            
         end if
         
         y(i1:n) = 0.0D+00
         
         return
         
      end if
      
      if ( alpha < 0.0D+00 ) then
         call xerror ( 'BESJ - order, alpha, less than zero.', 2, 1)
         return
      end if
      
      ialp = int ( alpha )
      fni = real ( ialp + n - 1, kind = 8 )
      fnf = alpha - real ( ialp, kind = 8 )
      dfn = fni + fnf
      fnu = dfn
      xo2 = x * 0.5D+00
      sxo2 = xo2 * xo2
c     
c     Decision tree for region where series, asymptotic expansion for x
c     to infinity and asymptotic expansion for nu to infinity are applied.
c     
      if ( sxo2 <= ( fnu+1.0D+00 ) ) then
         go to 90
      end if
      
      ta = max ( 20.0D+00, fnu )
      
      if ( ta < x ) then
         go to 120
      end if
      
      if ( 12.0D+00 < x ) then
         go to 110
      end if
      
      xo2l = log ( xo2 )
      ns = int ( sxo2 - fnu ) + 1
      go to 100
      
 90   continue
      
      fn = fnu
      fnp1 = fn + 1.0D+00
      xo2l = log ( xo2 )
      is = kt
      
      if ( x <= 0.50D+00 ) then
         go to 330
      end if
      
      ns = 0
      
 100  continue
      
      fni = fni + real ( ns, kind = 8 )
      dfn = fni + fnf
      fn = dfn
      fnp1 = fn + 1.0D+00
      is = kt
      
      if ( 0 < n - 1 + ns ) then 
         is = 3
      end if
      
      go to 330
      
 110  continue
      
      ans = max ( 36.0D+00 - fnu, 0.0D+00 )
      ns = int ( ans )
      fni = fni + real ( ns, kind = 8 )
      dfn = fni + fnf
      fn = dfn
      is = kt
      
      if ( 0 < n - 1 + ns ) then
         is = 3
      end if
      
      go to 130
      
 120  continue
      
      rtx = sqrt ( x )
      tau = rtwo * rtx
      ta = tau + fnulim(kt)
      
      if ( fnu <= ta ) then
         go to 480
      end if
      
      fn = fnu
      is = kt
c     
c     Uniform asymptotic expansion for NU to infinity.
c     
 130  continue
      
      i1 = abs ( 3 - is )
      i1 = max ( i1, 1 )
      flgjy = 1.0D+00
      
      call asyjy ( jairy, x, fn, flgjy, i1, temp(is), wk, iflw )
      
      if ( iflw /= 0 ) then
         go to 380
      end if
      
      go to (320, 450, 620), is
      
 310  continue
      
      temp(1) = temp(3)
      kt = 1
      
 320  continue
      
      is = 2
      fni = fni - 1.0D+00
      dfn = fni + fnf
      fn = dfn
      
      if ( i1 == 2 ) then
         go to 450
      end if
      
      go to 130
c     
c     Series for (x/2)**2<=nu+1
c     
 330  continue
      
      gln = alngam ( fnp1 )
      arg = fn * xo2l - gln
      
      if ( arg < (-elim1) ) then
         go to 400
      end if
      
      earg = exp ( arg )
      
 340  continue
      
      s = 1.0D+00
      
      if ( x < tol ) then
         go to 360
      end if
      
      ak = 3.0D+00
      t2 = 1.0D+00
      t = 1.0D+00
      s1 = fn
      
      do k = 1, 17
         s2 = t2 + s1
         t = - t * sxo2 / s2
         s = s + t
         if ( abs ( t ) < tol ) then
            exit
         end if
         t2 = t2 + ak
         ak = ak + 2.0D+00
         s1 = s1 + fn
      end do
      
 360  continue
      
      temp(is) = s * earg
      go to (370, 450, 610), is
      
 370  continue
      earg = earg * fn / xo2
      fni = fni - 1.0D+00
      dfn = fni + fnf
      fn = dfn
      is = 2
      go to 340
c     
c     Set underflow value and update parameters
c     
 380  continue
      
      y(nn) = 0.0D+00
      nn = nn - 1
      fni = fni - 1.0D+00
      dfn = fni + fnf
      fn = dfn
      if ( nn-1 ) 440, 390, 130
      
 390  continue
      
      kt = 2
      is = 2
      go to 130
      
 400  continue
      
      y(nn) = 0.0D+00
      nn = nn - 1
      fnp1 = fn
      fni = fni - 1.0D+00
      dfn = fni + fnf
      fn = dfn
      if ( nn-1 ) 440, 410, 420
      
 410  continue
      kt = 2
      is = 2
      
 420  continue
      
      if ( sxo2 <= fnp1 ) then
         go to 430
      end if
      
      go to 130
      
 430  continue
      
      arg = arg - xo2l + log ( fnp1 )
      if ( arg < (-elim1) ) then
         go to 400
      end if
      
      go to 330
      
 440  nz = n - nn
      return
c     
c     Backward recursion section
c     
 450  continue
      nz = n - nn
      if ( kt == 2 ) then
         go to 470
      end if
c     
c     Backward recur from index ALPHA+NN-1 to ALPHA.
c     
      y(nn) = temp(1)
      y(nn-1) = temp(2)
      if ( nn == 2 ) then
         return
      end if
      
      trx = 2.0D+00 / x
      dtm = fni
      tm = ( dtm + fnf ) * trx
      k = nn + 1
      
      do i = 3, nn
         k = k - 1
         y(k-2) = tm * y(k-1) - y(k)
         dtm = dtm - 1.0D+00
         tm = ( dtm + fnf ) * trx
      end do
      
      return
      
 470  continue
      
      y(1) = temp(2)
      return
c     
c     Asymptotic expansion for X to infinity with forward recursion in
c     oscillatory region max ( 20, NU ) < X, provided the last member
c     of the sequence is also in the region.
c     
 480  continue
      
      in = int ( alpha - tau + 2.0D+00 )
      
      if ( in <= 0 ) then
         go to 490
      end if
      
      idalp = ialp - in - 1
      kt = 1
      go to 500
      
 490  continue
      
      idalp = ialp
      in = 0
      
 500  continue
      
      is = kt
      fidal = real ( idalp, kind = 8 )
      dalpha = fidal + fnf
      arg = x - pidt * dalpha - pdf
      sa = sin ( arg )
      sb = cos ( arg )
      coef = rttp / rtx
      etx = 8.0D+00 * x
      
 510  continue
      
      dtm = fidal + fidal
      dtm = dtm * dtm
      tm = 0.0D+00
      
      if ( fidal == 0.0D+00 .and. abs ( fnf ) < tol ) then
         go to 520
      end if
      
      tm = 4.0D+00 * fnf * ( fidal + fidal + fnf )
      
 520  continue
      
      trx = dtm - 1.0D+00
      t2 = ( trx + tm ) / etx
      s2 = t2
      relb = tol * abs ( t2 )
      t1 = etx
      s1 = 1.0D+00
      fn = 1.0D+00
      ak = 8.0D+00
      
      do k = 1, 13
         t1 = t1 + etx
         fn = fn + ak
         trx = dtm - fn
         ap = trx + tm
         t2 = -t2 * ap / t1
         s1 = s1 + t2
         t1 = t1 + etx
         ak = ak + 8.0D+00
         fn = fn + ak
         trx = dtm - fn
         ap = trx + tm
         t2 = t2 * ap / t1
         s2 = s2 + t2
         if ( abs ( t2 ) <= relb ) then
            exit
         end if
         ak = ak + 8.0D+00
      end do
      
 540  continue
      
      temp(is) = coef * ( s1 * sb - s2 * sa )
      
      if ( is == 2 ) then
         go to 560
      end if
      
 550  continue
      
      fidal = fidal + 1.0D+00
      dalpha = fidal + fnf
      is = 2
      tb = sa
      sa = -sb
      sb = tb
      go to 510
c     
c     Forward recursion section
c     
 560  continue
      
      if ( kt == 2 ) then
         go to 470
      end if
      
      s1 = temp(1)
      s2 = temp(2)
      tx = 2.0D+00 / x
      tm = dalpha * tx
      
      if ( in == 0 ) then
         go to 580
      end if
c     
c     Forward recursion to index alpha
c     
      do i = 1, in
         s = s2
         s2 = tm * s2 - s1
         tm = tm + tx
         s1 = s
      end do
      
      if ( nn == 1 ) then
         go to 600
      end if
      
      s = s2
      s2 = tm * s2 - s1
      tm = tm + tx
      s1 = s
      
 580  continue
c     
c     Forward recursion from index ALPHA to ALPHA+N-1.
c     
      y(1) = s1
      y(2) = s2
      
      do i = 3, nn
         y(i) = tm * y(i-1) - y(i-2)
         tm = tm + tx
      end do
      
      return
      
 600  continue
      
      y(1) = s2
      return
c     
c     Backward recursion with normalization by
c     asymptotic expansion for nu to infinity or power series.
c     
 610  continue
c     
c     Computation of last order for series normalization
c     
      akm = max ( 3.0D+00 - fn, 0.0D+00 )
      km = int ( akm )
      tfn = fn + real ( km, kind = 8 )
      ta = ( gln + tfn - 0.9189385332D+00 - 0.0833333333D+00 / tfn ) 
     &     / ( tfn + 0.5D+00 )
      ta = xo2l - ta
      tb = - ( 1.0D+00 -1.5D+00 / tfn ) / tfn
      akm = tolln / ( - ta + sqrt ( ta * ta - tolln * tb ) ) + 1.5D+00
      in = km + int ( akm )
      go to 660
      
 620  continue
c     
c     Computation of last order for asymptotic expansion normalization
c     
      gln = wk(3) + wk(2)
      if ( 30.0D+00 < wk(6) ) then
         go to 640
      end if
      
      rden = ( pp(4) * wk(6) + pp(3) ) * wk(6) + 1.0D+00
      rzden = pp(1) + pp(2) * wk(6)
      ta = rzden / rden
      
      if ( wk(1) < 0.10D+00 ) then
         go to 630
      end if
      
      tb = gln / wk(5)
      go to 650
      
 630  continue
      
      tb = ( 1.259921049D+00 + ( 0.1679894730D+00 + 0.0887944358D+00 
     &     * wk(1) ) * wk(1) ) / wk(7)
      go to 650
      
 640  continue
      ta = 0.5D+00 * tolln / wk(4)
      ta=( ( 0.0493827160D+00 * ta - 0.1111111111D+00 ) * ta 
     &     + 0.6666666667D+00 ) * ta * wk(6)
      
      if ( wk(1) < 0.10D+00 ) then
         go to 630
      end if
      
      tb = gln / wk(5)
      
 650  continue
      
      in = int ( ta / tb + 1.5D+00 )
      
      if ( inlim < in ) then
         go to 310
      end if
      
 660  continue
      dtm = fni + real ( in, kind = 8 )
      trx = 2.0D+00 / x
      tm = ( dtm + fnf ) * trx
      ta = 0.0D+00
      tb = tol
      kk = 1
      
 670  continue
c     
c     Backward recur unindexed
c     
      do i = 1, in
         s = tb
         tb = tm * tb - ta
         ta = s
         dtm = dtm - 1.0D+00
         tm = ( dtm + fnf ) * trx
      end do
c     
c     Normalization.
c     
      if ( kk == 1 ) then
         
         ta = ( ta / tb ) * temp(3)
         tb = temp(3)
         kk = 2
         in = ns
         
         if ( ns /= 0 ) then
            go to 670
         end if
         
      end if
      
      y(nn) = tb
      nz = n - nn
      
      if ( nn == 1 ) then
         return
      end if
      
      k = nn - 1
      y(k) = tm * tb - ta
      
      if ( nn == 2 ) then
         return
      end if
      
      dtm = dtm - 1.0D+00
      tm = ( dtm + fnf ) * trx
      km = k - 1
c     
c     Backward recur indexed
c     
      do i = 1, km
         y(k-1) = tm * y(k) - y(k+1)
         dtm = dtm - 1.0D+00
         tm = ( dtm + fnf ) * trx
         k = k - 1
      end do
      
      return
      end
c*****************************************************************************80




      function alngam ( x )
c*****************************************************************************80
c
cc ALNGAM computes the log of the absolute value of the Gamma function.
c
c  Definition:
c
c    The Gamma function is defined as
c
c      GAMMA(Z) = INTEGRAL ( 0 <= T < Infinity) T**(Z-1) EXP ( -T ) DT
c
c    If Z is a positive integer, GAMMA(Z) = (Z-1)!, the factorial.
c
c    There is a special value:
c
c      GAMMA(0.5) = SQRT ( PI ).
c
c  Modified:
c
c    31 May 2000
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) X, the argument of the gamma function.
c
c    Output, real ( kind = 8 ) ALNGAM, the logarithm of the absolute
c    value of GAMMA(X).
c
      implicit none
      
      real ( kind = 8 ) alngam
      real ( kind = 8 ) d9lgmc
      real ( kind = 8 ), save :: dxrel = 0.0D+00
      real ( kind = 8 ) gamma
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) sinpiy
      real ( kind = 8 ), parameter :: sq2pil = 0.91893853320467274D+00
      real ( kind = 8 ), parameter :: sqpi2l = 0.22579135264472743D+00
      real ( kind = 8 ) x
      real ( kind = 8 ), save :: xmax = 0.0D+00
      real ( kind = 8 ) y
      
      if ( xmax == 0.0D+00 ) then
         xmax = huge ( xmax ) / log ( huge ( xmax ) )
         dxrel = sqrt ( epsilon ( dxrel ) )
      end if
      
      y = abs ( x )
      
      if ( y <= 10.0D+00 ) then
         alngam = log ( abs ( gamma ( x ) ) )
         return
      end if
      
      if ( xmax < y ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'ALNGAM - Fatal error!'
         write ( *, '(a)' ) '  |X| is so big that ALNGAM will overflow.'
         stop
      end if
      
      if ( 0.0D+00 < x ) then
         alngam = sq2pil + ( x - 0.5D+00 ) * log ( x ) - x + d9lgmc (y)
         return
      end if
      
      sinpiy = abs ( sin ( pi * y ) )
      
      if ( sinpiy == 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'ALNGAM - Fatal error!'
         write ( *, '(a)' ) '  X is a negative integer.'
         stop
      end if
      
      if ( abs((x-real(int(x-0.5D0),kind=8))/x) < dxrel ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'ALNGAM - Warning:'
         write ( *, '(a)' ) '  Answer has reduced precision'
         write ( *, '(a)' ) '  X is very near a negative integer.'
      end if
      
      alngam = sqpi2l + ( x - 0.5D+00 ) * log ( y ) - x - log (sinpiy)
     &     - d9lgmc ( y )
      
      return
      end
c*****************************************************************************80




      subroutine asyjy ( funjy, x, fnu, flgjy, in, y, wk, iflw )
c*****************************************************************************80
c
cc ASYJY computes high order Bessel functions J and Y.
c
c  Description:
c
c    ASYJY implements the uniform asymptotic expansion of
c    the J and Y Bessel functions for 35 <= FNU and 0.0 < X.  
c
c    The forms are identical except for a change
c    in sign of some of the terms.  This change in sign is
c    accomplished by means of the flag FLGJY = 1 or -1.
c 
c    On FLGJY = 1 the Airy functions AI(X) and DAI(X) are
c    supplied by the external function JAIRY, and on
c    FLGJY = -1 the Airy functions BI(X) and DBI(X) are
c    supplied by the external funtion YAIRY.
c
c  Modified:
c
c    25 August 2001
c
c  Author:
c
c    Donald Amos
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, external FUNJY, is the function JAIRY or YAIRY.
c
c    Input, real ( kind = 8 ) X, the argument, which must be greater than 0.
c
c    Input, real ( kind = 8 ) FNU, the order of the first Bessel function.
c    FNU is generally at least 35.
c
c    Input, real ( kind = 8 ) FLGJY, a selection flag
c     1.0D+00 gives the J function
c    -1.0D+00 gives the Y function
c
c    Input, integer IN, the number of functions desired, which should be
c    1 or 2.
c
c    Output, real ( kind = 8 ) Y(IN), contains the desired function values.
c
c    Output, integer IFLW, a flag indicating underflow or overflow
c    return variables for BESJ only.
c
c    Output, real ( kind = 8 ) WK(7), contains the following values:
c
c      wk(1) = 1 - (x/fnu)**2 = w**2
c      wk(2) = sqrt ( abs ( wk(1) ) )
c      wk(3) = abs ( wk(2) - atan ( wk(2) ) )  or
c              abs ( ln((1 + wk(2) )/ ( x / fnu ) ) - wk(2))
c            = abs ( (2/3)*zeta**(3/2))
c      wk(4) = fnu*wk(3)
c      wk(5) = (1.5*wk(3) * fnu)**(1/3) = sqrt ( zeta ) * fnu**(1/3)
c      wk(6) = sign ( 1.0, w**2 ) * wk(5)**2 
c            = sign ( 1.0, w**2 ) * zeta * fnu**(2/3)
c      wk(7) = fnu**(1/3)
c
      implicit none
      
      real ( kind = 8 ) abw2
      real ( kind = 8 ) akm
      real ( kind = 8 ) alfa(26,4)
      real ( kind = 8 ) alfa1(26,2)
      real ( kind = 8 ) alfa2(26,2)
      real ( kind = 8 ) ap
      real ( kind = 8 ), parameter, dimension ( 8 ) :: ar = (/ 
     &     8.35503472222222D-02, 1.28226574556327D-01, 
     &     2.91849026464140D-01, 8.81627267443758D-01, 
     &     3.32140828186277D+00, 1.49957629868626D+01, 
     &     7.89230130115865D+01, 4.74451538868264D+02 /)
      real ( kind = 8 ) asum
      real ( kind = 8 ) az
      real ( kind = 8 ) beta(26,5)
      real ( kind = 8 ) beta1(26,2)
      real ( kind = 8 ) beta2(26,2)
      real ( kind = 8 ) beta3(26,1)
      real ( kind = 8 ) br(10)
      real ( kind = 8 ) bsum
      real ( kind = 8 ) c(65)
      real ( kind = 8 ), parameter :: con1 = 6.66666666666667D-01
      real ( kind = 8 ), parameter :: con2 = 3.33333333333333D-01
      real ( kind = 8 ), parameter :: con548 = 1.04166666666667D-01
      real ( kind = 8 ) cr(10)
      real ( kind = 8 ) crz32
      real ( kind = 8 ) d1mach
      real ( kind = 8 ) dfi
      real ( kind = 8 ) elim
      real ( kind = 8 ) dr(10)
      real ( kind = 8 ) fi
      real ( kind = 8 ) flgjy
      real ( kind = 8 ) fn
      real ( kind = 8 ) fnu
      real ( kind = 8 ) fn2
      external funjy
      real ( kind = 8 ) gama(26)
      integer i
      integer i1mach
      integer iflw
      integer in
      integer j
      integer jn
      integer jr
      integer ju
      integer k
      integer kb
      integer klast
      integer kmax(5)
      integer kp1
      integer ks
      integer ksp1
      integer kstemp
      integer l
      integer lr
      integer lrp1
      real ( kind = 8 ) phi
      real ( kind = 8 ) rcz
      real ( kind = 8 ) rden
      real ( kind = 8 ) relb
      real ( kind = 8 ) rfn2
      real ( kind = 8 ) rtz
      real ( kind = 8 ) rzden
      real ( kind = 8 ) sa
      real ( kind = 8 ) sb
      real ( kind = 8 ) suma
      real ( kind = 8 ) sumb
      real ( kind = 8 ) s1
      real ( kind = 8 ) ta
      real ( kind = 8 ) tau
      real ( kind = 8 ) tb
      real ( kind = 8 ) tfn
      real ( kind = 8 ) tol
      real ( kind = 8 ), save :: tols = -6.90775527898214D+00
      real ( kind = 8 ) t2
      real ( kind = 8 ) upol(10)
      real ( kind = 8 ) wk(*)
      real ( kind = 8 ) x
      real ( kind = 8 ) xx
      real ( kind = 8 ) y(*)
      real ( kind = 8 ) z
      real ( kind = 8 ) z32
      
      equivalence (alfa(1,1),alfa1(1,1))
      equivalence (alfa(1,3),alfa2(1,1))
      equivalence (beta(1,1),beta1(1,1))
      equivalence (beta(1,3),beta2(1,1))
      equivalence (beta(1,5),beta3(1,1))
      
      data  br(1),br(2),br(3),br(4),br(5),br(6),br(7),br(8),br(9),br(10) 
     &     /-1.45833333333333D-01, -9.87413194444444D-02,
     &     -1.43312053915895D-01, -3.17227202678414D-01,
     &     -9.42429147957120D-01, -3.51120304082635D+00,
     &     -1.57272636203680D+01,-8.22814390971859D+01,
     &     -4.92355370523671D+02,-3.31621856854797D+03/

      data c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10), 
     &     c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18), 
     &     c(19), c(20), c(21), c(22), c(23), c(24)
     &     /-2.08333333333333D-01,        1.25000000000000D-01, 
     &     3.34201388888889D-01,       -4.01041666666667D-01, 
     &     7.03125000000000D-02,       -1.02581259645062D+00, 
     &     1.84646267361111D+00,       -8.91210937500000D-01, 
     &     7.32421875000000D-02,        4.66958442342625D+00, 
     &     -1.12070026162230D+01,        8.78912353515625D+00, 
     &     -2.36408691406250D+00,        1.12152099609375D-01, 
     &     -2.82120725582002D+01,        8.46362176746007D+01, 
     &     -9.18182415432400D+01,        4.25349987453885D+01, 
     &     -7.36879435947963D+00,        2.27108001708984D-01, 
     &     2.12570130039217D+02,       -7.65252468141182D+02, 
     &     1.05999045252800D+03,       -6.99579627376133D+02/

      data c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32), 
     &     c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40), 
     &     c(41), c(42), c(43), c(44), c(45), c(46), c(47), c(48)
     &     / 2.18190511744212D+02,       -2.64914304869516D+01, 
     &     5.72501420974731D-01,       -1.91945766231841D+03, 
     &     8.06172218173731D+03,       -1.35865500064341D+04, 
     &     1.16553933368645D+04,       -5.30564697861340D+03, 
     &     1.20090291321635D+03,       -1.08090919788395D+02, 
     &     1.72772750258446D+00,        2.02042913309661D+04, 
     &     -9.69805983886375D+04,        1.92547001232532D+05, 
     &     -2.03400177280416D+05,        1.22200464983017D+05, 
     &     -4.11926549688976D+04,        7.10951430248936D+03, 
     &     -4.93915304773088D+02,        6.07404200127348D+00, 
     &     -2.42919187900551D+05,        1.31176361466298D+06, 
     &     -2.99801591853811D+06,        3.76327129765640D+06/

      data c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56), 
     &     c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64), 
     &     c(65)
     &     / -2.81356322658653D+06,        1.26836527332162D+06, 
     &     -3.31645172484564D+05,        4.52187689813627D+04, 
     &     -2.49983048181121D+03,        2.43805296995561D+01, 
     &     3.28446985307204D+06,       -1.97068191184322D+07, 
     &     5.09526024926646D+07,       -7.41051482115327D+07, 
     &     6.63445122747290D+07,       -3.75671766607634D+07, 
     &     1.32887671664218D+07,       -2.78561812808645D+06, 
     &     3.08186404612662D+05,       -1.38860897537170D+04, 
     &     1.10017140269247D+02/

      data alfa1(1,1), alfa1(2,1), alfa1(3,1), alfa1(4,1), alfa1(5,1), 
     &     alfa1(6,1), alfa1(7,1), alfa1(8,1), alfa1(9,1), alfa1(10,1), 
     &     alfa1(11,1),alfa1(12,1),alfa1(13,1),alfa1(14,1),alfa1(15,1), 
     &     alfa1(16,1),alfa1(17,1),alfa1(18,1),alfa1(19,1),alfa1(20,1), 
     &     alfa1(21,1),alfa1(22,1),alfa1(23,1),alfa1(24,1),alfa1(25,1), 
     &     alfa1(26,1) /
     &     -4.44444444444444D-03,-9.22077922077922D-04, 
     &     -8.84892884892885D-05, 1.65927687832450D-04, 
     &     2.46691372741793D-04, 2.65995589346255D-04, 
     &     2.61824297061501D-04, 2.48730437344656D-04, 
     &     2.32721040083232D-04, 2.16362485712365D-04, 
     &     2.00738858762752D-04, 1.86267636637545D-04, 
     &     1.73060775917876D-04, 1.61091705929016D-04, 
     &     1.50274774160908D-04, 1.40503497391270D-04, 
     &     1.31668816545923D-04, 1.23667445598253D-04, 
     &     1.16405271474738D-04, 1.09798298372713D-04, 
     &     1.03772410422993D-04, 9.82626078369363D-05, 
     &     9.32120517249503D-05, 8.85710852478712D-05, 
     &     8.42963105715700D-05, 8.03497548407791D-05/

      data alfa1(1,2), alfa1(2,2), alfa1(3,2), alfa1(4,2), alfa1(5,2), 
     &     alfa1(6,2), alfa1(7,2), alfa1(8,2), alfa1(9,2), alfa1(10,2), 
     &     alfa1(11,2),alfa1(12,2),alfa1(13,2),alfa1(14,2),alfa1(15,2), 
     &     alfa1(16,2),alfa1(17,2),alfa1(18,2),alfa1(19,2),alfa1(20,2), 
     &     alfa1(21,2),alfa1(22,2),alfa1(23,2),alfa1(24,2),alfa1(25,2), 
     &     alfa1(26,2)     / 6.93735541354589D-04, 2.32241745182922D-04, 
     &     -1.41986273556691D-05,-1.16444931672049D-04,
     &     -1.50803558053049D-04,-1.55121924918096D-04,
     &     -1.46809756646466D-04,-1.33815503867491D-04, 
     &     -1.19744975684254D-04,-1.06184319207974D-04,
     &     -9.37699549891194D-05,-8.26923045588193D-05,
     &     -7.29374348155221D-05,-6.44042357721016D-05, 
     &     -5.69611566009369D-05,-5.04731044303562D-05,
     &     -4.48134868008883D-05,-3.98688727717599D-05,
     &     -3.55400532972042D-05,-3.17414256609022D-05, 
     &     -2.83996793904175D-05,-2.54522720634871D-05,
     &     -2.28459297164725D-05,-2.05352753106481D-05,
     &     -1.84816217627666D-05,-1.66519330021394D-05/

      data alfa2(1,1), alfa2(2,1), alfa2(3,1), alfa2(4,1), alfa2(5,1), 
     &     alfa2(6,1), alfa2(7,1), alfa2(8,1), alfa2(9,1), alfa2(10,1), 
     &     alfa2(11,1),alfa2(12,1),alfa2(13,1),alfa2(14,1),alfa2(15,1), 
     &     alfa2(16,1),alfa2(17,1),alfa2(18,1),alfa2(19,1),alfa2(20,1), 
     &     alfa2(21,1),alfa2(22,1),alfa2(23,1),alfa2(24,1),alfa2(25,1), 
     &     alfa2(26,1) /-3.54211971457744D-04,-1.56161263945159D-04, 
     &     3.04465503594936D-05, 1.30198655773243D-04, 
     &     1.67471106699712D-04, 1.70222587683593D-04, 
     &     1.56501427608595D-04, 1.36339170977445D-04, 
     &     1.14886692029825D-04, 9.45869093034688D-05, 
     &     7.64498419250898D-05, 6.07570334965197D-05, 
     &     4.74394299290509D-05, 3.62757512005344D-05, 
     &     2.69939714979225D-05, 1.93210938247939D-05, 
     &     1.30056674793963D-05, 7.82620866744497D-06, 
     &     3.59257485819352D-06, 1.44040049814252D-07, 
     &     -2.65396769697939D-06,-4.91346867098486D-06,
     &     -6.72739296091248D-06, -8.17269379678658D-06,
     &     -9.31304715093561D-06,-1.02011418798016D-05/

      data alfa2(1,2), alfa2(2,2), alfa2(3,2), alfa2(4,2), alfa2(5,2), 
     &     alfa2(6,2), alfa2(7,2), alfa2(8,2), alfa2(9,2), alfa2(10,2), 
     &     alfa2(11,2),alfa2(12,2),alfa2(13,2),alfa2(14,2),alfa2(15,2), 
     &     alfa2(16,2),alfa2(17,2),alfa2(18,2),alfa2(19,2),alfa2(20,2), 
     &     alfa2(21,2),alfa2(22,2),alfa2(23,2),alfa2(24,2),alfa2(25,2), 
     &     alfa2(26,2) / 3.78194199201773D-04, 2.02471952761816D-04, 
     &     -6.37938506318862D-05,-2.38598230603006D-04,
     &     -3.10916256027362D-04, -3.13680115247576D-04,
     &     -2.78950273791323D-04,-2.28564082619141D-04, 
     &     -1.75245280340847D-04,-1.25544063060690D-04,
     &     -8.22982872820208D-05, -4.62860730588116D-05,
     &     -1.72334302366962D-05, 5.60690482304602D-06, 
     &     2.31395443148287D-05, 3.62642745856794D-05, 
     &     4.58006124490189D-05, 5.24595294959114D-05, 
     &     5.68396208545815D-05, 5.94349820393104D-05, 
     &     6.06478527578422D-05, 6.08023907788436D-05, 
     &     6.01577894539460D-05, 5.89199657344698D-05, 
     &     5.72515823777593D-05, 5.52804375585853D-05/

       data beta1(1,1), beta1(2,1), beta1(3,1), beta1(4,1), beta1(5,1), 
     &     beta1(6,1), beta1(7,1), beta1(8,1), beta1(9,1), beta1(10,1), 
     &     beta1(11,1),beta1(12,1),beta1(13,1),beta1(14,1),beta1(15,1), 
     &     beta1(16,1),beta1(17,1),beta1(18,1),beta1(19,1),beta1(20,1), 
     &     beta1(21,1),beta1(22,1),beta1(23,1),beta1(24,1),beta1(25,1), 
     &     beta1(26,1) / 1.79988721413553D-02, 5.59964911064388D-03, 
     &     2.88501402231133D-03, 1.80096606761054D-03, 
     &     1.24753110589199D-03, 9.22878876572938D-04, 
     &     7.14430421727287D-04, 5.71787281789705D-04, 
     &     4.69431007606482D-04, 3.93232835462917D-04, 
     &     3.34818889318298D-04, 2.88952148495752D-04, 
     &     2.52211615549573D-04, 2.22280580798883D-04, 
     &     1.97541838033063D-04, 1.76836855019718D-04, 
     &     1.59316899661821D-04, 1.44347930197334D-04, 
     &     1.31448068119965D-04, 1.20245444949303D-04, 
     &     1.10449144504599D-04, 1.01828770740567D-04, 
     &     9.41998224204238D-05, 8.74130545753834D-05, 
     &     8.13466262162801D-05, 7.59002269646219D-05/

       data beta1(1,2), beta1(2,2), beta1(3,2), beta1(4,2), beta1(5,2), 
     &      beta1(6,2), beta1(7,2), beta1(8,2), beta1(9,2), beta1(10,2), 
     &      beta1(11,2),beta1(12,2),beta1(13,2),beta1(14,2),beta1(15,2), 
     &      beta1(16,2),beta1(17,2),beta1(18,2),beta1(19,2),beta1(20,2), 
     &      beta1(21,2),beta1(22,2),beta1(23,2),beta1(24,2),beta1(25,2), 
     &      beta1(26,2) /-1.49282953213429D-03,-8.78204709546389D-04, 
     &      -5.02916549572035D-04,-2.94822138512746D-04,
     &      -1.75463996970783D-04, -1.04008550460816D-04,
     &      -5.96141953046458D-05,-3.12038929076098D-05, 
     &      -1.26089735980230D-05,-2.42892608575730D-07,
     &      8.05996165414274D-06,  1.36507009262147D-05,
     &      1.73964125472926D-05, 1.98672978842134D-05, 
     &      2.14463263790823D-05, 2.23954659232457D-05,
     &      2.28967783814713D-05,  2.30785389811178D-05,
     &      2.30321976080909D-05, 2.28236073720349D-05, 
     &      2.25005881105292D-05, 2.20981015361991D-05,
     &      2.16418427448104D-05,  2.11507649256221D-05,
     &      2.06388749782171D-05, 2.01165241997082D-05/

       data beta2(1,1), beta2(2,1), beta2(3,1), beta2(4,1), beta2(5,1), 
     &      beta2(6,1), beta2(7,1), beta2(8,1), beta2(9,1), beta2(10,1), 
     &      beta2(11,1),beta2(12,1),beta2(13,1),beta2(14,1),beta2(15,1), 
     &      beta2(16,1),beta2(17,1),beta2(18,1),beta2(19,1),beta2(20,1), 
     &      beta2(21,1),beta2(22,1),beta2(23,1),beta2(24,1),beta2(25,1), 
     &      beta2(26,1) / 5.52213076721293D-04, 4.47932581552385D-04, 
     &      2.79520653992021D-04, 1.52468156198447D-04,
     &      6.93271105657044D-05,  1.76258683069991D-05,
     &      -1.35744996343269D-05,-3.17972413350427D-05, 
     &      -4.18861861696693D-05,-4.69004889379141D-05,
     &      -4.87665447413787D-05, -4.87010031186735D-05,
     &      -4.74755620890087D-05,-4.55813058138628D-05, 
     &      -4.33309644511266D-05,-4.09230193157750D-05,
     &      -3.84822638603221D-05, -3.60857167535411D-05,
     &      -3.37793306123367D-05,-3.15888560772110D-05, 
     &      -2.95269561750807D-05,-2.75978914828336D-05,
     &      -2.58006174666884D-05, -2.41308356761280D-05,
     &      -2.25823509518346D-05,-2.11479656768913D-05/

       data beta2(1,2), beta2(2,2), beta2(3,2), beta2(4,2), beta2(5,2), 
     &      beta2(6,2), beta2(7,2), beta2(8,2), beta2(9,2), beta2(10,2), 
     &      beta2(11,2),beta2(12,2),beta2(13,2),beta2(14,2),beta2(15,2), 
     &      beta2(16,2),beta2(17,2),beta2(18,2),beta2(19,2),beta2(20,2), 
     &      beta2(21,2),beta2(22,2),beta2(23,2),beta2(24,2),beta2(25,2), 
     &      beta2(26,2) /-4.74617796559960D-04,-4.77864567147321D-04, 
     &      -3.20390228067038D-04,-1.61105016119962D-04,
     &      -4.25778101285435D-05,  3.44571294294968D-05,
     &      7.97092684075675D-05, 1.03138236708272D-04, 
     &      1.12466775262204D-04, 1.13103642108481D-04,
     &      1.08651634848774D-04,  1.01437951597662D-04,
     &      9.29298396593364D-05, 8.40293133016090D-05, 
     &      7.52727991349134D-05, 6.69632521975731D-05,
     &      5.92564547323195D-05,  5.22169308826976D-05,
     &      4.58539485165361D-05, 4.01445513891487D-05, 
     &      3.50481730031328D-05, 3.05157995034347D-05,
     &      2.64956119950516D-05,  2.29363633690998D-05,
     &      1.97893056664022D-05, 1.70091984636413D-05/

      data beta3(1,1), beta3(2,1), beta3(3,1), beta3(4,1), beta3(5,1), 
     &      beta3(6,1), beta3(7,1), beta3(8,1), beta3(9,1), beta3(10,1), 
     &      beta3(11,1),beta3(12,1),beta3(13,1),beta3(14,1),beta3(15,1), 
     &      beta3(16,1),beta3(17,1),beta3(18,1),beta3(19,1),beta3(20,1), 
     &      beta3(21,1),beta3(22,1),beta3(23,1),beta3(24,1),beta3(25,1), 
     &      beta3(26,1) / 7.36465810572578D-04, 8.72790805146194D-04, 
     &      6.22614862573135D-04, 2.85998154194304D-04,
     &      3.84737672879366D-06, -1.87906003636972D-04,
     &      -2.97603646594555D-04,-3.45998126832656D-04, 
     &      -3.53382470916038D-04,-3.35715635775049D-04,
     &      -3.04321124789040D-04, -2.66722723047613D-04,
     &      -2.27654214122820D-04,-1.89922611854562D-04, 
     &      -1.55058918599094D-04,-1.23778240761874D-04,
     &      -9.62926147717644D-05, -7.25178327714425D-05,
     &      -5.22070028895634D-05,-3.50347750511901D-05, 
     &      -2.06489761035552D-05,-8.70106096849767D-06,
     &      1.13698686675100D-06,  9.16426474122779D-06,
     &      1.56477785428873D-05, 2.08223629482467D-05/

      data gama(1),   gama(2),   gama(3),   gama(4),   gama(5), 
     &     gama(6),   gama(7),   gama(8),   gama(9),   gama(10), 
     &     gama(11),  gama(12),  gama(13),  gama(14),  gama(15), 
     &     gama(16),  gama(17),  gama(18),  gama(19),  gama(20), 
     &     gama(21),  gama(22),  gama(23),  gama(24),  gama(25), 
     &     gama(26) / 6.29960524947437D-01, 2.51984209978975D-01, 
     &     1.54790300415656D-01, 1.10713062416159D-01,
     &     8.57309395527395D-02,  6.97161316958684D-02,
     &     5.86085671893714D-02, 5.04698873536311D-02, 
     &     4.42600580689155D-02, 3.93720661543510D-02,
     &     3.54283195924455D-02,  3.21818857502098D-02,
     &     2.94646240791158D-02, 2.71581677112934D-02, 
     &     2.51768272973862D-02, 2.34570755306079D-02,
     &     2.19508390134907D-02,  2.06210828235646D-02,
     &     1.94388240897881D-02, 1.83810633800683D-02, 
     &     1.74293213231963D-02, 1.65685837786612D-02,
     &     1.57865285987918D-02,  1.50729501494096D-02,
     &     1.44193250839955D-02, 1.38184805735342D-02/
c     
c     I1MACH(14) replaces I1MACH(11) in a double precision code
c     I1MACH(15) replaces I1MACH(12) in a double precision code
c     
      ta = epsilon ( ta )
      tol = max ( ta, 1.0D-15 )
      tb = d1mach(5)
      ju = i1mach(15)
      
      if ( flgjy /= 1.0D+00 ) then
         jr = i1mach(14)
         elim = 2.303D+00 * tb * ( real ( - ju ) - real ( jr ) )
      else
         elim = 2.303D+00 * ( tb * real ( - ju ) - 3.0D+00 )
      end if
      
      fn = fnu
      iflw = 0
      
      do jn = 1, in
         
         xx = x / fn
         wk(1) = 1.0D+00 - xx * xx
         abw2 = abs ( wk(1) )
         wk(2) = sqrt ( abw2 )
         wk(7) = fn**con2
         
         if ( 0.27750D+00 < abw2 ) then
            go to 80
         end if
c     
c     Asymptotic expansion.
c     
c     Cases near x=fn, abs ( 1-(x/fn)**2 ) <= 0.2775
c     coefficients of asymptotic expansion by series
c     
c     ZETA and truncation for a(zeta) and b(zeta) series
c     
c     KMAX is the truncation index for a(zeta) and b(zeta) series = max ( 2, sa )
c     
         if ( abw2 == 0.0D+00 ) then
            sa = 0.0D+00
         else
            sa = tols / log ( abw2 )
         end if
         
         sb = sa
         
         do i = 1, 5
            akm = max ( sa, 2.0D+00 )
            kmax(i) = int ( akm )
            sa = sa + sb
         end do
         
         kb = kmax(5)
         klast = kb - 1
         sa = gama(kb)
         
         do k = 1, klast
            kb = kb - 1
            sa = sa * wk(1) + gama(kb)
         end do
         
         z = wk(1) * sa
         az = abs ( z )
         rtz = sqrt ( az )
         wk(3) = con1 * az * rtz
         wk(4) = wk(3) * fn
         wk(5) = rtz * wk(7)
         wk(6) = - wk(5) * wk(5)
         
         if ( 0.0D+00 < z ) then
            
            if ( elim < wk(4) ) then
               iflw = 1
               return
            end if
            
            wk(6) = -wk(6)
            
         end if
         
         phi = sqrt ( sqrt ( sa + sa + sa + sa ) )
c     
c     b(zeta) for s=0
c     
         kb = kmax(5)
         klast = kb - 1
         sb = beta(kb,1)
         
         do k = 1, klast
            kb = kb - 1
            sb = sb * wk(1) + beta(kb,1)
         end do
         
         ksp1 = 1
         fn2 = fn * fn
         rfn2 = 1.0D+00 / fn2
         rden = 1.0D+00
         asum = 1.0D+00
         relb = tol * abs ( sb )
         bsum = sb
         
         do ks = 1, 4
            
            ksp1 = ksp1 + 1
            rden = rden * rfn2
c     
c     a(zeta) and b(zeta) for s=1,2,3,4
c     
            kstemp = 5 - ks
            kb = kmax(kstemp)
            klast = kb - 1
            sa = alfa(kb,ks)
            sb = beta(kb,ksp1)
            
            do k = 1, klast
               kb = kb - 1
               sa = sa * wk(1) + alfa(kb,ks)
               sb = sb * wk(1) + beta(kb,ksp1)
            end do
            
            ta = sa * rden
            tb = sb * rden
            asum = asum + ta
            bsum = bsum + tb
            
            if ( abs ( ta ) <= tol .and. abs ( tb ) <= relb ) then
               exit
            end if
            
         end do
         
         bsum = bsum / ( fn * wk(7) )
         
         go to 160
         
 80      continue
         
         upol(1) = 1.0D+00
         tau = 1.0D+00 / wk(2)
         t2 = 1.0D+00 / wk(1)
c     
c     Cases for sqrt ( 1.2775 ) < (x/fn).
c     
         if ( wk(1) < 0.0D+00 ) then
            
            wk(3) = abs ( wk(2) - atan ( wk(2) ) )
            wk(4) = wk(3) * fn
            rcz = -con1 / wk(4)
            z32 = 1.5D+00 * wk(3)
            rtz = z32**con2
            wk(5) = rtz * wk(7)
            wk(6) = -wk(5) * wk(5)
c     
c     Cases for (x/fn) < sqrt ( 0.7225 )
c     
         else
            
            wk(3) = abs ( log ( ( 1.0D+00 + wk(2) ) / xx ) - wk(2) )
            wk(4) = wk(3) * fn
            rcz = con1 / wk(4)
            
            if ( elim < wk(4) ) then
               iflw = 1
               return
            end if
            
            z32 = 1.5D+00 * wk(3)
            rtz = z32**con2
            wk(7) = fn**con2
            wk(5) = rtz * wk(7)
            wk(6) = wk(5) * wk(5)
            
         end if
         
         phi = sqrt ( ( rtz + rtz ) * tau )
         tb = 1.0D+00
         asum = 1.0D+00
         tfn = tau / fn
         upol(2) = ( c(1) * t2 + c(2) ) * tfn
         crz32 = con548 * rcz
         bsum = upol(2) + crz32
         relb = tol * abs ( bsum )
         ap = tfn
         ks = 0
         kp1 = 2
         rzden = rcz
         l = 2
         
         do lr = 2, 8, 2
c     
c     Compute two U polynomials for next a(zeta) and b(zeta)
c     
            lrp1 = lr + 1
            
            do k = lr, lrp1
               ks = ks + 1
               kp1 = kp1 + 1
               l = l + 1
               s1 = c(l)
               do j = 2, kp1
                  l = l + 1
                  s1 = s1 * t2 + c(l)
               end do
               ap = ap * tfn
               upol(kp1) = ap * s1
               cr(ks) = br(ks) * rzden
               rzden = rzden * rcz
               dr(ks) = ar(ks) * rzden
            end do
            
            suma = upol(lrp1)
            sumb = upol(lr+2) + upol(lrp1) * crz32
            ju = lrp1
            
            do jr = 1, lr
               ju = ju - 1
               suma = suma + cr(jr) * upol(ju)
               sumb = sumb + dr(jr) * upol(ju)
            end do
            
            tb = -tb
            if ( 0.0D+00 < wk(1) ) then
               tb = abs ( tb )
            end if
            
            asum = asum + suma * tb
            bsum = bsum + sumb * tb
            
            if ( abs ( suma ) <= tol .and. abs ( sumb ) <= relb ) then
               exit
            end if
            
         end do
         
         tb = wk(5)
         if ( 0.0D+00 < wk(1) ) then
            tb = -tb
         end if
         
         bsum = bsum / tb
         
 160     continue
         
         call funjy ( wk(6), wk(5), wk(4), fi, dfi )
         y(jn) = flgjy * phi * ( fi * asum + dfi * bsum ) / wk(7)
         fn = fn - flgjy
         
      end do
      
      return
      end
c*****************************************************************************80




      function d1mach ( i )
c*****************************************************************************80
c
cc D1MACH returns double precision real machine constants.
c
c  Discussion:
c
c    Assuming that the internal representation of a double precision real
c    number is in base B, with T the number of base-B digits in the mantissa,
c    and EMIN the smallest possible exponent and EMAX the largest possible
c    exponent, then
c
c      D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c      D1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
c      D1MACH(3) = B**(-T), the smallest relative spacing.
c      D1MACH(4) = B**(1-T), the largest relative spacing.
c      D1MACH(5) = log10(B).
c
c  Modified:
c
c    24 April 2007
c
c  Author:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 5.
c
c    Output, real ( kind = 8 ) D1MACH, the value of the chosen parameter.
c
      implicit none
      
      real ( kind = 8 ) d1mach
      integer i
      
      if ( i < 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'D1MACH - Fatal error!'
         write ( *, '(a)' ) '  The input argument I is out of bounds.'
         write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
         write ( *, '(a,i12)' ) '  I = ', i
         d1mach = 0.0D+00
         stop
      else if ( i == 1 ) then
         d1mach = 4.450147717014403D-308
      else if ( i == 2 ) then
         d1mach = 8.988465674311579D+307
      else if ( i == 3 ) then
         d1mach = 1.110223024625157D-016
      else if ( i == 4 ) then
         d1mach = 2.220446049250313D-016
      else if ( i == 5 ) then
         d1mach = 0.301029995663981D+000
      else if ( 5 < i ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'D1MACH - Fatal error!'
         write ( *, '(a)' ) '  The input argument I is out of bounds.'
         write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
         write ( *, '(a,i12)' ) '  I = ', i
         d1mach = 0.0D+00
         stop
      end if
      
      return
      end
c*****************************************************************************80




      function i1mach ( i )
c*****************************************************************************80
c
cc I1MACH returns integer machine constants.
c
c  I/O unit numbers.
c
c    I1MACH(1) = the standard input unit.
c    I1MACH(2) = the standard output unit.
c    I1MACH(3) = the standard punch unit.
c    I1MACH(4) = the standard error message unit.
c
c  Words.
c
c    I1MACH(5) = the number of bits per integer storage unit.
c    I1MACH(6) = the number of characters per integer storage unit.
c
c  Integers.
c
c  Assume integers are represented in the S digit base A form:
c
c  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
c  where 0<=X(I)<A for I=0 to S-1.
c
c    I1MACH(7) = A, the base.
c    I1MACH(8) = S, the number of base A digits.
c    I1MACH(9) = A**S-1, the largest integer.
c
c  Floating point numbers
c
c  Assume floating point numbers are represented in the T digit base B form:
c
c    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
c
c  where 0<=X(I)<B for I = 1 to T, 0<X(1) and EMIN<=E<=EMAX
c
c    I1MACH(10) = B, the base.
c
c  Single precision
c
c    I1MACH(11) = T, the number of base B digits.
c    I1MACH(12) = EMIN, the smallest exponent E.
c    I1MACH(13) = EMAX, the largest exponent E.
c
c  Double precision
c
c    I1MACH(14) = T, the number of base B digits.
c    I1MACH(15) = EMIN, the smallest exponent E.
c    I1MACH(16) = EMAX, the largest exponent E.
c
c  To alter this function for a particular environment, the desired set of DATA
c  statements should be activated by removing the C from column 1.  On rare
c  machines, a STATIC statement may need to be added, but probably more systems
c  prohibit than require it.
c
c  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
c  consistency with the local operating system.  For FORTRAN 77, you may wish
c  to adjust the data statement so imach(6) is set to 1, and then to comment
c  out the executable test on I.EQ.6 below.
c
c  For IEEE-arithmetic machines (binary standard), the first set of constants
c  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
c
      implicit none
      
      integer i
      integer i1mach
      integer imach(16)
      integer output
      
      equivalence (imach(4),output)
c
c     IEEE arithmetic machines, such as the ATT 3B series, Motorola
c     68000 based machines such as the SUN 3 and ATT PC 7300, and
c     8087 based micros such asthe IBM PC and ATT 6300.
c     
      data imach( 1) /    5 /
      data imach( 2) /    6 /
      data imach( 3) /    7 /
      data imach( 4) /    6 /
      data imach( 5) /   32 /
      data imach( 6) /    4 /
      data imach( 7) /    2 /
      data imach( 8) /   31 /
      data imach( 9) / 2147483647 /
      data imach(10) /    2 /
      data imach(11) /   24 /
      data imach(12) / -125 /
      data imach(13) /  128 /
      data imach(14) /   53 /
      data imach(15) / -1021 /
      data imach(16) /  1024 /
c
c  ALLIANT FX/8 UNIX FORTRAN compiler.
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     6 /
c      data imach( 4) /     0 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    32 /
c      data imach( 9) /2147483647/
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -126 /
c      data imach(13) /   128 /
c      data imach(14) /    53 /
c      data imach(15) / -1022 /
c      data imach(16) /  1024 /
c
c  AMDAHL machines.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  32 /
c      data imach( 6) /   4 /
c      data imach( 7) /   2 /
c      data imach( 8) /  31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /  16 /
c      data imach(11) /   6 /
c      data imach(12) / -64 /
c      data imach(13) /  63 /
c      data imach(14) /  14 /
c      data imach(15) / -64 /
c      data imach(16) /  63 /
c
c  BURROUGHS 1700 system.
c
c      data imach( 1) /    7 /
c      data imach( 2) /    2 /
c      data imach( 3) /    2 /
c      data imach( 4) /    2 /
c      data imach( 5) /   36 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   33 /
c      data imach( 9) / Z1FFFFFFFF /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -256 /
c      data imach(13) /  255 /
c      data imach(14) /   60 /
c      data imach(15) / -256 /
c      data imach(16) /  255 /
c
c  BURROUGHS 5700 system.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  48 /
c      data imach( 6) /   6 /
c      data imach( 7) /   2 /
c      data imach( 8) /  39 /
c      data imach( 9) / O0007777777777777 /
c      data imach(10) /   8 /
c      data imach(11) /  13 /
c      data imach(12) / -50 /
c      data imach(13) /  76 /
c      data imach(14) /  26 /
c      data imach(15) / -50 /
c      data imach(16) /  76 /
c
c  BURROUGHS 6700/7700 systems.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  48 /
c      data imach( 6) /   6 /
c      data imach( 7) /   2 /
c      data imach( 8) /  39 /
c      data imach( 9) / O0007777777777777 /
c      data imach(10) /   8 /
c      data imach(11) /  13 /
c      data imach(12) / -50 /
c      data imach(13) /  76 /
c      data imach(14) /  26 /
c      data imach(15) / -32754 /
c      data imach(16) /  32780 /
c
c  CDC CYBER 170/180 series using NOS
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   60 /
c      data imach( 6) /   10 /
c      data imach( 7) /    2 /
c      data imach( 8) /   48 /
c      data imach( 9) / O"00007777777777777777" /
c      data imach(10) /    2 /
c      data imach(11) /   48 /
c      data imach(12) / -974 /
c      data imach(13) / 1070 /
c      data imach(14) /   96 /
c      data imach(15) / -927 /
c      data imach(16) / 1070 /
c
c  CDC CYBER 170/180 series using NOS/VE
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     7 /
c      data imach( 4) /     6 /
c      data imach( 5) /    64 /
c      data imach( 6) /     8 /
c      data imach( 7) /     2 /
c      data imach( 8) /    63 /
c      data imach( 9) / 9223372036854775807 /
c      data imach(10) /     2 /
c      data imach(11) /    47 /
c      data imach(12) / -4095 /
c      data imach(13) /  4094 /
c      data imach(14) /    94 /
c      data imach(15) / -4095 /
c      data imach(16) /  4094 /
c
c  CDC CYBER 200 series
c
c      data imach( 1) /      5 /
c      data imach( 2) /      6 /
c      data imach( 3) /      7 /
c      data imach( 4) /      6 /
c      data imach( 5) /     64 /
c      data imach( 6) /      8 /
c      data imach( 7) /      2 /
c      data imach( 8) /     47 /
c      data imach( 9) / X'00007FFFFFFFFFFF' /
c      data imach(10) /      2 /
c      data imach(11) /     47 /
c      data imach(12) / -28625 /
c      data imach(13) /  28718 /
c      data imach(14) /     94 /
c      data imach(15) / -28625 /
c      data imach(16) /  28718 /
c
c  CDC 6000/7000 series using FTN4.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   60 /
c      data imach( 6) /   10 /
c      data imach( 7) /    2 /
c      data imach( 8) /   48 /
c      data imach( 9) / 00007777777777777777B /
c      data imach(10) /    2 /
c      data imach(11) /   47 /
c      data imach(12) / -929 /
c      data imach(13) / 1070 /
c      data imach(14) /   94 /
c      data imach(15) / -929 /
c      data imach(16) / 1069 /
c
c  CDC 6000/7000 series using FTN5.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   60 /
c      data imach( 6) /   10 /
c      data imach( 7) /    2 /
c      data imach( 8) /   48 /
c      data imach( 9) / O"00007777777777777777" /
c      data imach(10) /    2 /
c      data imach(11) /   47 /
c      data imach(12) / -929 /
c      data imach(13) / 1070 /
c      data imach(14) /   94 /
c      data imach(15) / -929 /
c      data imach(16) / 1069 /
c
c  CONVEX C-1.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   53 /
c      data imach(15) /-1024 /
c      data imach(16) / 1023 /
c
c  CONVEX C-120 (native mode) without -R8 option
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    0 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   53 /
c      data imach(15) / -1023 /
c      data imach(16) /  1023 /
c
c  CONVEX C-120 (native mode) with -R8 option
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     0 /
c      data imach( 4) /     6 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    53 /
c      data imach(12) / -1023 /
c      data imach(13) /  1023 /
c      data imach(14) /    53 /
c      data imach(15) / -1023 /
c      data imach(16) /  1023 /
c
c  CONVEX C-120 (IEEE mode) without -R8 option
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    0 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -125 /
c      data imach(13) /  128 /
c      data imach(14) /   53 /
c      data imach(15) / -1021 /
c      data imach(16) /  1024 /
c
c  CONVEX C-120 (IEEE mode) with -R8 option
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     0 /
c      data imach( 4) /     6 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    53 /
c      data imach(12) / -1021 /
c      data imach(13) /  1024 /
c      data imach(14) /    53 /
c      data imach(15) / -1021 /
c      data imach(16) /  1024 /
c
c  CRAY 1, 2, XMP and YMP.
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /   102 /
c      data imach( 4) /     6 /
c      data imach( 5) /    64 /
c      data imach( 6) /     8 /
c      data imach( 7) /     2 /
c      data imach( 8) /    63 /
c      data imach( 9) /  777777777777777777777B /
c      data imach(10) /     2 /
c      data imach(11) /    47 /
c      data imach(12) / -8189 /
c      data imach(13) /  8190 /
c      data imach(14) /    94 /
c      data imach(15) / -8099 /
c      data imach(16) /  8190 /
c
c  DATA GENERAL ECLIPSE S/200.
c
c      data imach( 1) /   11 /
c      data imach( 2) /   12 /
c      data imach( 3) /    8 /
c      data imach( 4) /   10 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) /32767 /
c      data imach(10) /   16 /
c      data imach(11) /    6 /
c      data imach(12) /  -64 /
c      data imach(13) /   63 /
c      data imach(14) /   14 /
c      data imach(15) /  -64 /
c      data imach(16) /   63 /
c
c  ELXSI 6400
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     6 /
c      data imach( 4) /     6 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    32 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -126 /
c      data imach(13) /   127 /
c      data imach(14) /    53 /
c      data imach(15) / -1022 /
c      data imach(16) /  1023 /
c
c  HARRIS 220
c
c      data imach( 1) /       5 /
c      data imach( 2) /       6 /
c      data imach( 3) /       0 /
c      data imach( 4) /       6 /
c      data imach( 5) /      24 /
c      data imach( 6) /       3 /
c      data imach( 7) /       2 /
c      data imach( 8) /      23 /
c      data imach( 9) / 8388607 /
c      data imach(10) /       2 /
c      data imach(11) /      23 /
c      data imach(12) /    -127 /
c      data imach(13) /     127 /
c      data imach(14) /      38 /
c      data imach(15) /    -127 /
c      data imach(16) /     127 /
c
c  HARRIS SLASH 6 and SLASH 7.
c
c      data imach( 1) /       5 /
c      data imach( 2) /       6 /
c      data imach( 3) /       0 /
c      data imach( 4) /       6 /
c      data imach( 5) /      24 /
c      data imach( 6) /       3 /
c      data imach( 7) /       2 /
c      data imach( 8) /      23 /
c      data imach( 9) / 8388607 /
c      data imach(10) /       2 /
c      data imach(11) /      23 /
c      data imach(12) /    -127 /
c      data imach(13) /     127 /
c      data imach(14) /      38 /
c      data imach(15) /    -127 /
c      data imach(16) /     127 /
c
c  HONEYWELL DPS 8/70 and 600/6000 series.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /   43 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / O377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   63 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /
c
c  HP 2100, 3 word double precision option with FTN4
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    4 /
c      data imach( 4) /    1 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) / 32767 /
c      data imach(10) /    2 /
c      data imach(11) /   23 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   39 /
c      data imach(15) / -128 /
c      data imach(16) /  127 /
c
c  HP 2100, 4 word double precision option with FTN4
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    4 /
c      data imach( 4) /    1 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) / 32767 /
c      data imach(10) /    2 /
c      data imach(11) /   23 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   55 /
c      data imach(15) / -128 /
c      data imach(16) /  127 /
c
c  HP 9000
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     6 /
c      data imach( 4) /     7 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    32 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -126 /
c      data imach(13) /   127 /
c      data imach(14) /    53 /
c      data imach(15) / -1015 /
c      data imach(16) /  1017 /
c
c  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
c  and PERKIN ELMER (INTERDATA) 3230.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   7 /
c      data imach( 4) /   6 /
c      data imach( 5) /  32 /
c      data imach( 6) /   4 /
c      data imach( 7) /   2 /
c      data imach( 8) /  31 /
c      data imach( 9) / Z7FFFFFFF /
c      data imach(10) /  16 /
c      data imach(11) /   6 /
c      data imach(12) / -64 /
c      data imach(13) /  63 /
c      data imach(14) /  14 /
c      data imach(15) / -64 /
c      data imach(16) /  63 /
c
c  IBM PC - Microsoft FORTRAN
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     6 /
c      data imach( 4) /     0 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -126 /
c      data imach(13) /   127 /
c      data imach(14) /    53 /
c      data imach(15) / -1022 /
c      data imach(16) /  1023 /
c
c  IBM PC - Professional FORTRAN and Lahey FORTRAN
c
c      data imach( 1) /     4 /
c      data imach( 2) /     7 /
c      data imach( 3) /     7 /
c      data imach( 4) /     0 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -126 /
c      data imach(13) /   127 /
c      data imach(14) /    53 /
c      data imach(15) / -1022 /
c      data imach(16) /  1023 /
c
c  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
c  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
c  constants with Y's.
c
c      data imach( 1) /   5 /
c      data imach( 2) /   6 /
c      data imach( 3) /   6 /
c      data imach( 4) /   6 /
c      data imach( 5) /  32 /
c      data imach( 6) /   4 /
c      data imach( 7) /   2 /
c      data imach( 8) /  31 /
c      data imach( 9) / Z'7FFFFFFF' /
c      data imach(10) /  16 /
c      data imach(11) /   6 /
c      data imach(12) / -64 /
c      data imach(13) /  62 /
c      data imach(14) /  14 /
c      data imach(15) / -64 /
c      data imach(16) /  62 /
c
c  PDP-10 (KA processor).
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    5 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / "377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   54 /
c      data imach(15) / -101 /
c      data imach(16) /  127 /
c
c  PDP-10 (KI processor).
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    5 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / "377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   62 /
c      data imach(15) / -128 /
c      data imach(16) /  127 /
c
c  PDP-11 FORTRANS supporting 32-bit integer arithmetic.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /
c
c  PDP-11 FORTRANS supporting 16-bit integer arithmetic.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) / 32767 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /
c
c  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
c  supplied by Igor Bray.
c
c      data imach( 1) /            1 /
c      data imach( 2) /            1 /
c      data imach( 3) /            2 /
c      data imach( 4) /            1 /
c      data imach( 5) /           32 /
c      data imach( 6) /            4 /
c      data imach( 7) /            2 /
c      data imach( 8) /           31 /
c      data imach( 9) / :17777777777 /
c      data imach(10) /            2 /
c      data imach(11) /           23 /
c      data imach(12) /         -127 /
c      data imach(13) /         +127 /
c      data imach(14) /           47 /
c      data imach(15) /       -32895 /
c      data imach(16) /       +32637 /
c
c  SEQUENT BALANCE 8000.
c
c      data imach( 1) /     0 /
c      data imach( 2) /     0 /
c      data imach( 3) /     7 /
c      data imach( 4) /     0 /
c      data imach( 5) /    32 /
c      data imach( 6) /     1 /
c      data imach( 7) /     2 /
c      data imach( 8) /    31 /
c      data imach( 9) /  2147483647 /
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -125 /
c      data imach(13) /   128 /
c      data imach(14) /    53 /
c      data imach(15) / -1021 /
c      data imach(16) /  1024 /
c
c  SUN Microsystems UNIX F77 compiler.
c
c      data imach( 1) /     5 /
c      data imach( 2) /     6 /
c      data imach( 3) /     6 /
c      data imach( 4) /     0 /
c      data imach( 5) /    32 /
c      data imach( 6) /     4 /
c      data imach( 7) /     2 /
c      data imach( 8) /    32 /
c      data imach( 9) /2147483647/
c      data imach(10) /     2 /
c      data imach(11) /    24 /
c      data imach(12) /  -126 /
c      data imach(13) /   128 /
c      data imach(14) /    53 /
c      data imach(15) / -1022 /
c      data imach(16) /  1024 /
c
c  SUN 3 (68881 or FPA)
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    6 /
c      data imach( 4) /    0 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -125 /
c      data imach(13) /  128 /
c      data imach(14) /   53 /
c      data imach(15) / -1021 /
c      data imach(16) /  1024 /
c
c  UNIVAC 1100 series.
c  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
c  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
c  instead.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   36 /
c      data imach( 6) /    6 /
c      data imach( 7) /    2 /
c      data imach( 8) /   35 /
c      data imach( 9) / O377777777777 /
c      data imach(10) /    2 /
c      data imach(11) /   27 /
c      data imach(12) / -128 /
c      data imach(13) /  127 /
c      data imach(14) /   60 /
c      data imach(15) /-1024 /
c      data imach(16) / 1023 /
c
c  VAX.
c
c      data imach( 1) /    5 /
c      data imach( 2) /    6 /
c      data imach( 3) /    7 /
c      data imach( 4) /    6 /
c      data imach( 5) /   32 /
c      data imach( 6) /    4 /
c      data imach( 7) /    2 /
c      data imach( 8) /   31 /
c      data imach( 9) / 2147483647 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /
c
c  Z80 microprocessor.
c
c      data imach( 1) /    1 /
c      data imach( 2) /    1 /
c      data imach( 3) /    0 /
c      data imach( 4) /    1 /
c      data imach( 5) /   16 /
c      data imach( 6) /    2 /
c      data imach( 7) /    2 /
c      data imach( 8) /   15 /
c      data imach( 9) / 32767 /
c      data imach(10) /    2 /
c      data imach(11) /   24 /
c      data imach(12) / -127 /
c      data imach(13) /  127 /
c      data imach(14) /   56 /
c      data imach(15) / -127 /
c      data imach(16) /  127 /
c
      if ( i < 1 .or. 16 < i )then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'I1MACH - Fatal error!'
         write ( *, '(a,i6)' ) '  I is out of bounds:', i
         i1mach = 0
         stop
      else
         i1mach = imach(i)
      end if
      
      return
      end
c*****************************************************************************80




      subroutine jairy ( x, rx, c, ai, dai )
c*****************************************************************************80
c
cc JAIRY computes the Airy function and its derivative.
c
c  Author:
c
c    D E Amos,
c    S L Daniel,
c    M K Weston.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) X, the argument.
c
c    Input, real ( kind = 8 ) RX, sqrt ( abs ( X ) ).
c
c    Input, real ( kind = 8 ) C, = 2 * ( ABS ( X )**1.5 ) / 3, computed 
c    by ASYJY.
c
c    Output, real ( kind = 8 ) AI, the value of the Airy function.
c
c    Output, real ( kind = 8 ) DAI, the derivative of the Airy function.
c
      implicit none
      
      real ( kind = 8 ) a(15)
      real ( kind = 8 ) ai
      real ( kind = 8 ) ajn(19)
      real ( kind = 8 ) ajp(19)
      real ( kind = 8 ), parameter, dimension ( 14 ) :: ak1 = (/ 
     &     2.20423090987793D-01,-1.25290242787700D-01, 
     &     1.03881163359194D-02,     8.22844152006343D-04,
     &     -2.34614345891226D-04, 1.63824280172116D-05, 
     &     3.06902589573189D-07,-1.29621999359332D-07, 
     &     8.22908158823668D-09, 1.53963968623298D-11, 
     &     -3.39165465615682D-11, 2.03253257423626D-12, 
     &     -1.10679546097884D-14, -5.16169497785080D-15 /)
      real ( kind = 8 ), parameter, dimension ( 23 ) :: ak2 = (/ 
     &     2.74366150869598D-01, 5.39790969736903D-03, 
     &     -1.57339220621190D-03, 4.27427528248750D-04,
     &     -1.12124917399925D-04, 2.88763171318904D-05,
     &     -7.36804225370554D-06, 1.87290209741024D-06, 
     &     -4.75892793962291D-07, 1.21130416955909D-07,
     &     -3.09245374270614D-08, 7.92454705282654D-09,
     &     -2.03902447167914D-09, 5.26863056595742D-10, 
     &     -1.36704767639569D-10, 3.56141039013708D-11,
     &     -9.31388296548430D-12, 2.44464450473635D-12,
     &     -6.43840261990955D-13, 1.70106030559349D-13, 
     &     -4.50760104503281D-14, 1.19774799164811D-14,
     &     -3.19077040865066D-15 /)
      real ( kind = 8 ) ak3(14)
      real ( kind = 8 ) b(15)
      real ( kind = 8 ) c
      real ( kind = 8 ) ccv  
      real ( kind = 8 ), parameter :: con2 = 5.03154716196777D+00
      real ( kind = 8 ), parameter :: con3 = 3.80004589867293D-01
      real ( kind = 8 ), parameter :: con4 = 8.33333333333333D-01
      real ( kind = 8 ), parameter :: con5 = 8.66025403784439D-01
      real ( kind = 8 ) cv
      real ( kind = 8 ) da(15)
      real ( kind = 8 ) dai
      real ( kind = 8 ) dajn(19)
      real ( kind = 8 ) dajp(19)
      real ( kind = 8 ) dak1(14)
      real ( kind = 8 ) dak2(24)
      real ( kind = 8 ) dak3(14)
      real ( kind = 8 ) db(15)
      real ( kind = 8 ) ec
      real ( kind = 8 ) e1
      real ( kind = 8 ) e2
      real ( kind = 8 ), parameter :: fpi12 = 1.30899693899575D+00
      real ( kind = 8 ) f1
      real ( kind = 8 ) f2
      integer i
      integer j
      integer, parameter :: m1 = 12
      integer m1d
      integer, parameter :: m2 = 21
      integer m2d
      integer, parameter :: m3 = 17
      integer m3d
      integer, parameter :: m4 = 13
      integer m4d
      integer, parameter :: n1 = 14
      integer n1d
      integer, parameter :: n2 = 23
      integer n2d
      integer, parameter :: n3 = 19
      integer n3d
      integer, parameter :: n4 = 15
      integer n4d
      real ( kind = 8 ) rtrx
      real ( kind = 8 ) rx
      real ( kind = 8 ) scv
      real ( kind = 8 ) t
      real ( kind = 8 ) temp1
      real ( kind = 8 ) temp2
      real ( kind = 8 ) tt
      real ( kind = 8 ) x
      
      data ak3(1), ak3(2), ak3(3), ak3(4), ak3(5), ak3(6), ak3(7), 
     &     ak3(8), ak3(9), ak3(10),ak3(11),ak3(12),ak3(13), ak3(14) /
     &     2.80271447340791D-01,-1.78127042844379D-03, 
     &     4.03422579628999D-05,-1.63249965269003D-06, 
     &     9.21181482476768D-08,-6.52294330229155D-09, 
     &     5.47138404576546D-10,-5.24408251800260D-11, 
     &     5.60477904117209D-12,-6.56375244639313D-13, 
     &     8.31285761966247D-14,-1.12705134691063D-14, 
     &     1.62267976598129D-15,-2.46480324312426D-16/

      data ajp(1), ajp(2), ajp(3), ajp(4), ajp(5), ajp(6), ajp(7), 
     &     ajp(8), ajp(9), ajp(10),ajp(11),ajp(12),ajp(13),ajp(14), 
     &     ajp(15),ajp(16),ajp(17),ajp(18),ajp(19) /
     &     7.78952966437581D-02,-1.84356363456801D-01, 
     &     3.01412605216174D-02, 3.05342724277608D-02,
     &     -4.95424702513079D-03, -1.72749552563952D-03, 
     &     2.43137637839190D-04, 5.04564777517082D-05, 
     &     -6.16316582695208D-06,-9.03986745510768D-07, 
     &     9.70243778355884D-08, 1.09639453305205D-08,
     &     -1.04716330588766D-09,-9.60359441344646D-11, 
     &     8.25358789454134D-12, 6.36123439018768D-13,
     &     -4.96629614116015D-14,-3.29810288929615D-15, 
     &     2.35798252031104D-16/

      data ajn(1), ajn(2), ajn(3), ajn(4), ajn(5), ajn(6), ajn(7), 
     &     ajn(8), ajn(9), ajn(10),ajn(11),ajn(12),ajn(13),ajn(14), 
     &     ajn(15),ajn(16),ajn(17),ajn(18),ajn(19) / 
     &     3.80497887617242D-02,-2.45319541845546D-01, 
     &     1.65820623702696D-01, 7.49330045818789D-02,
     &     -2.63476288106641D-02,-5.92535597304981D-03, 
     &     1.44744409589804D-03, 2.18311831322215D-04, 
     &     -4.10662077680304D-05,-4.66874994171766D-06, 
     &     7.15218807277160D-07, 6.52964770854633D-08,
     &     -8.44284027565946D-09,-6.44186158976978D-10, 
     &     7.20802286505285D-11, 4.72465431717846D-12,
     &     -4.66022632547045D-13,-2.67762710389189D-14, 
     &     2.36161316570019D-15/

      data a(1), a(2), a(3),  a(4),  a(5),  a(6),  a(7), 
     &     a(8), a(9), a(10), a(11), a(12), a(13), a(14), a(15) / 
     &     4.90275424742791D-01, 1.57647277946204D-03, 
     &     -9.66195963140306D-05, 1.35916080268815D-07, 
     &     2.98157342654859D-07, -1.86824767559979D-08,
     &     -1.03685737667141D-09, 3.28660818434328D-10, 
     &     -2.57091410632780D-11,-2.32357655300677D-12,
     &     9.57523279048255D-13, -1.20340828049719D-13,
     &     -2.90907716770715D-15, 4.55656454580149D-15, 
     &     -9.99003874810259D-16/

      data b(1), b(2), b(3),  b(4),  b(5),  b(6),  b(7), 
     &     b(8), b(9), b(10), b(11), b(12), b(13), b(14), b(15) / 
     &     2.78593552803079D-01,-3.52915691882584D-03, 
     &     -2.31149677384994D-05, 4.71317842263560D-06,
     &     -1.12415907931333D-07, -2.00100301184339D-08, 
     &     2.60948075302193D-09,-3.55098136101216D-11, 
     &     -3.50849978423875D-11, 5.83007187954202D-12,
     &     -2.04644828753326D-13, -1.10529179476742D-13, 
     &     2.87724778038775D-14,-2.88205111009939D-15, 
     &     -3.32656311696166D-16/

      data n1d,n2d,n3d,n4d/14,24,19,15/
      data m1d,m2d,m3d,m4d/12,22,17,13/

      data dak1(1), dak1(2), dak1(3), dak1(4), dak1(5), dak1(6), 
     &     dak1(7), dak1(8), dak1(9), dak1(10), dak1(11), 
     &     dak1(12), dak1(13),dak1(14)/ 
     &     2.04567842307887D-01,-6.61322739905664D-02, 
     &     -8.49845800989287D-03, 3.12183491556289D-03,
     &     -2.70016489829432D-04,-6.35636298679387D-06, 
     &     3.02397712409509D-06,-2.18311195330088D-07, 
     &     -5.36194289332826D-10, 1.13098035622310D-09,
     &     -7.43023834629073D-11, 4.28804170826891D-13, 
     &     2.23810925754539D-13,-1.39140135641182D-14/

       data dak2(1), dak2(2), dak2(3), dak2(4), dak2(5), dak2(6), 
     &     dak2(7), dak2(8), dak2(9), dak2(10),dak2(11),dak2(12), 
     &     dak2(13),dak2(14),dak2(15),dak2(16),dak2(17),dak2(18), 
     &     dak2(19),dak2(20),dak2(21),dak2(22),dak2(23),dak2(24) / 
     &     2.93332343883230D-01,-8.06196784743112D-03, 
     &     2.42540172333140D-03,-6.82297548850235D-04, 
     &     1.85786427751181D-04,-4.97457447684059D-05, 
     &     1.32090681239497D-05,-3.49528240444943D-06, 
     &     9.24362451078835D-07,-2.44732671521867D-07, 
     &     6.49307837648910D-08,-1.72717621501538D-08, 
     &     4.60725763604656D-09,-1.23249055291550D-09, 
     &     3.30620409488102D-10,-8.89252099772401D-11, 
     &     2.39773319878298D-11,-6.48013921153450D-12, 
     &     1.75510132023731D-12,-4.76303829833637D-13, 
     &     1.29498241100810D-13,-3.52679622210430D-14, 
     &     9.62005151585923D-15,-2.62786914342292D-15/

       data dak3(1), dak3(2), dak3(3), dak3(4), dak3(5), dak3(6), 
     &      dak3(7), dak3(8), dak3(9), dak3(10),dak3(11),dak3(12), 
     &      dak3(13),dak3(14)/ 
     &      2.84675828811349D-01, 2.53073072619080D-03, 
     &      -4.83481130337976D-05, 1.84907283946343D-06,
     &      -1.01418491178576D-07, 7.05925634457153D-09,
     &      -5.85325291400382D-10, 5.56357688831339D-11, 
     &      -5.90889094779500D-12, 6.88574353784436D-13,
     &      -8.68588256452194D-14, 1.17374762617213D-14,
     &      -1.68523146510923D-15, 2.55374773097056D-16/

       data dajp(1), dajp(2), dajp(3), dajp(4), dajp(5), dajp(6), 
     &      dajp(7), dajp(8), dajp(9), dajp(10),dajp(11),dajp(12), 
     &      dajp(13),dajp(14),dajp(15),dajp(16),dajp(17),dajp(18), 
     &      dajp(19)        / 
     &      6.53219131311457D-02,-1.20262933688823D-01, 
     &      9.78010236263823D-03, 1.67948429230505D-02,
     &      -1.97146140182132D-03,-8.45560295098867D-04, 
     &      9.42889620701976D-05, 2.25827860945475D-05, 
     &      -2.29067870915987D-06,-3.76343991136919D-07, 
     &      3.45663933559565D-08, 4.29611332003007D-09,
     &      -3.58673691214989D-10,-3.57245881361895D-11, 
     &      2.72696091066336D-12, 2.26120653095771D-13,
     &      -1.58763205238303D-14,-1.12604374485125D-15, 
     &      7.31327529515367D-17/

       data dajn(1), dajn(2), dajn(3), dajn(4), dajn(5), dajn(6), 
     &      dajn(7), dajn(8), dajn(9), dajn(10),dajn(11),dajn(12), 
     &      dajn(13),dajn(14),dajn(15),dajn(16),dajn(17),dajn(18), 
     &      dajn(19)        / 
     &      1.08594539632967D-02, 8.53313194857091D-02, 
     &      -3.15277068113058D-01,-8.78420725294257D-02, 
     &      5.53251906976048D-02, 9.41674060503241D-03,
     &      -3.32187026018996D-03,-4.11157343156826D-04, 
     &      1.01297326891346D-04, 9.87633682208396D-06,
     &      -1.87312969812393D-06,-1.50798500131468D-07, 
     &      2.32687669525394D-08, 1.59599917419225D-09, 
     &      -2.07665922668385D-10,-1.24103350500302D-11, 
     &      1.39631765331043D-12, 7.39400971155740D-14,
     &      -7.32887475627500D-15/
       
       data da(1),  da(2),  da(3),  da(4),  da(5),  da(6),  da(7), 
     &      da(8),  da(9),  da(10), da(11), da(12), da(13), da(14), 
     &      da(15) / 
     &      4.91627321104601D-01, 3.11164930427489D-03, 
     &      8.23140762854081D-05,-4.61769776172142D-06,
     &      -6.13158880534626D-08, 2.87295804656520D-08,
     &      -1.81959715372117D-09,-1.44752826642035D-10, 
     &      4.53724043420422D-11,-3.99655065847223D-12,
     &      -3.24089119830323D-13, 1.62098952568741D-13,
     &      -2.40765247974057D-14, 1.69384811284491D-16, 
     &      8.17900786477396D-16/

       data db(1),  db(2),  db(3),  db(4),  db(5),  db(6),  db(7), 
     &      db(8),  db(9),  db(10), db(11), db(12), db(13), db(14), 
     &      db(15) /
     &      -2.77571356944231D-01, 4.44212833419920D-03, 
     &      -8.42328522190089D-05,-2.58040318418710D-06, 
     &      3.42389720217621D-07,-6.24286894709776D-09,
     &      -2.36377836844577D-09, 3.16991042656673D-10, 
     &      -4.40995691658191D-12,-5.18674221093575D-12, 
     &      9.64874015137022D-13,-4.90190576608710D-14,
     &      -1.77253430678112D-14, 5.55950610442662D-15, 
     &      -7.11793337579530D-16/

       if ( x < 0.0D+00 ) then
          go to 90
       end if
       
       if ( 5.0D+00 < c ) then
          go to 60
       end if
       
       if ( 1.20D+00 < x ) then
          go to 30
       end if
       
       t = ( x + x - 1.2D+00 ) * con4
       tt = t + t
       j = n1
       f1 = ak1(j)
       f2 = 0.0D+00
       
       do i = 1, m1
          j = j - 1
          temp1 = f1
          f1 = tt * f1 - f2 + ak1(j)
          f2 = temp1
       end do
       
       ai = t * f1 - f2 + ak1(1)
       
       j = n1d
       f1 = dak1(j)
       f2 = 0.0D+00
       do i = 1, m1d
          j = j - 1
          temp1 = f1
          f1 = tt * f1 - f2 + dak1(j)
          f2 = temp1
       end do
       
       dai = -( t * f1 - f2 + dak1(1) )
       
       return
       
 30    continue
       
       t = ( x + x - con2 ) * con3
       tt = t + t
       j = n2
       f1 = ak2(j)
       f2 = 0.0D+00
       
       do i = 1, m2
          j = j - 1
          temp1 = f1
          f1 = tt * f1 - f2 + ak2(j)
          f2 = temp1
       end do
       
       rtrx = sqrt ( rx )
       ec = exp ( -c )
       ai = ec * ( t * f1 - f2 + ak2(1) ) / rtrx
       j = n2d
       f1 = dak2(j)
       f2 = 0.0D+00
       
       do i = 1, m2d
          j = j - 1
          temp1 = f1
          f1 = tt * f1 - f2 + dak2(j)
          f2 = temp1
       end do
       
       dai = -ec * ( t * f1 - f2 + dak2(1) ) * rtrx
       return
       
 60    continue
       
       t = 10.0D+00 / c - 1.0D+00
       tt = t + t
       j = n1
       f1 = ak3(j)
       f2 = 0.0D+00
       
       do i = 1, m1
          j = j - 1
          temp1 = f1
          f1 = tt * f1 - f2 + ak3(j)
          f2 = temp1
       end do
       
       rtrx = sqrt ( rx )
       ec = exp ( -c )
       ai = ec * ( t * f1 - f2 + ak3(1) ) / rtrx
       j = n1d
       f1 = dak3(j)
       f2 = 0.0D+00
       
       do i = 1, m1d
          j = j - 1
          temp1 = f1
          f1 = tt * f1 - f2 + dak3(j)
          f2 = temp1
       end do
       
       dai = -rtrx * ec * ( t * f1 - f2 + dak3(1) )
       return
       
 90    continue
       
       if ( 5.0D+00 < c ) then
          go to 120
       end if
       
       t = 0.4D+00 * c - 1.0D+00
       tt = t + t
       j = n3
       f1 = ajp(j)
       e1 = ajn(j)
       f2 = 0.0D+00
       e2 = 0.0D+00
       
       do i = 1, m3
          j = j - 1
          temp1 = f1
          temp2 = e1
          f1 = tt * f1 - f2 + ajp(j)
          e1 = tt * e1 - e2 + ajn(j)
          f2 = temp1
          e2 = temp2
       end do
       
       ai = ( t * e1 - e2 + ajn(1) ) - x * ( t * f1 - f2 + ajp(1) )
       j = n3d
       f1 = dajp(j)
       e1 = dajn(j)
       f2 = 0.0D+00
       e2 = 0.0D+00
       
       do i = 1, m3d
          j = j - 1
          temp1 = f1
          temp2 = e1
          f1 = tt * f1 - f2 + dajp(j)
          e1 = tt * e1 - e2 + dajn(j)
          f2 = temp1
          e2 = temp2
       end do
       
       dai = x * x * ( t * f1 - f2 + dajp(1) ) + (t * e1 - e2 + dajn(1))
       
       return
       
 120   continue
       
       t = 10.0D+00 / c - 1.0D+00
       tt = t + t
       j = n4
       f1 = a(j)
       e1 = b(j)
       f2 = 0.0D+00
       e2 = 0.0D+00
       
       do i = 1, m4
          j = j - 1
          temp1 = f1
          temp2 = e1
          f1 = tt * f1 - f2 + a(j)
          e1 = tt * e1 - e2 + b(j)
          f2 = temp1
          e2 = temp2
       end do
       
       temp1 = t * f1 - f2 + a(1)
       temp2 = t * e1 - e2 + b(1)
       rtrx = sqrt ( rx )
       cv = c - fpi12
       ccv = cos ( cv )
       scv = sin ( cv )
       ai = ( temp1 * ccv - temp2 * scv ) / rtrx
       j = n4d
       f1 = da(j)
       e1 = db(j)
       f2 = 0.0D+00
       e2 = 0.0D+00
       
       do i = 1, m4d
          j = j - 1
          temp1 = f1
          temp2 = e1
          f1 = tt * f1 - f2 + da(j)
          e1 = tt * e1 - e2 + db(j)
          f2 = temp1
          e2 = temp2
       end do
       
       temp1 = t * f1 - f2 + da(1)
       temp2 = t * e1 - e2 + db(1)
       e1 = ccv * con5 + 0.5D+00 * scv
       e2 = scv * con5 - 0.5D+00 * ccv
       dai = ( temp1 * e1 - temp2 * e2 ) * rtrx
       
       return
       end
c*****************************************************************************80




      subroutine xerror ( messg, nerr, level )
c*****************************************************************************80
c
cc XERROR processes a diagnostic error message.
c
c  Discussion:
c
c    XERROR processes a diagnostic message, in a manner
c    determined by the value of level and the current value
c    of the library error control flag, kontrl.
c    See subroutine xsetf for details.
c
c  Example:
c
c    call xerror('smooth -- num was zero.',1,2)
c
c    call xerror('integ  -- less than full accuracy achieved.',2,1)
c
c    call xerror('rooter -- actual zero of f found before interval fully collapsed.',3,0)
c
c    call xerror('exp    -- underflows being set to zero.',1,-1)
c
c  Modified:
c
c    13 August 2005
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character ( len = * ) MESSG, the message to be processed, 
c    containing no more than 72 characters.
c
c    Input, integer NERR, the error number associated with this message.
c    NERR must not be zero.
c
c    Input, integer LEVEL, the error category.
c    2 means this is an unconditionally fatal error.
c    1 means this is a recoverable error.  (i.e., it is
c      non-fatal if XSETF has been appropriately called.)
c    0 means this is a warning message only.
c    -1 means this is a warning message which is to be printed at most once, 
c      regardless of how many times this call is executed.
c
      implicit none
      
      integer level
      character ( len = * ) messg
      integer nerr
      integer nmessg
      
      nmessg = len ( messg )
      
      call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.0D+00,0.0D+00)
      
      return
      end
c*****************************************************************************80




      subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
c*****************************************************************************80
c
cc XERRWV processes an error message that includes numeric information.
c
c  Discussion:
c
c    XERRWV processes a diagnostic message, in a manner
c    determined by the value of level and the current value
c    of the library error control flag, kontrl.
c    (see subroutine xsetf for details.)
c    in addition, up to two integer values and two real
c    values may be printed along with the message.
c
c  Example:
c
c    call xerrwv ( 'smooth -- num (=i1) was zero.',29,1,2,1,num,0,0,0.,0.)
c
c    call xerrwv ( 'quadxy -- requested error (r1) less than minimum(r2).',
c   &   54,77,1,0,0,0,2,errreq,errmin)
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character ( len = * ) MESSG, the message to be processed.
c
c    Input, integer NMESSG, the number of characters in MESSG.
c
c    Input, integer NERR, the error number associated with this message.
c    NERR must not be zero.
c
c    Input, integer LEVEL, the error category.
c    2 means this is an unconditionally fatal error.
c    1 means this is a recoverable error.  (i.e., it is
c      non-fatal if xsetf has been appropriately called.)
c    0 means this is a warning message only.
c    -1 means this is a warning message which is to be printed at most 
c      once, regardless of how many times this call is executed.
c
c    Input, integer NI, the number of integer values to be printed. (0 to 2)
c
c    Input, integer I1, I2, the first and second integer values.
c
c    Input, integer NR, the number of real values to be printed. (0 to 2)
c
c    Input, real ( kind = 8 ) R1, R2, the first and second real values.
c
      implicit none
      
      character ( len = 37 ) form
      integer i
      integer i1
      integer i1mach
      integer i2
      integer ifatal
      integer isizei
      integer isizef
      integer iunit
      integer j4save
      integer junk
      integer kdummy
      integer kount
      integer kunit
      integer lerr
      integer level
      character ( len = 20 ) lfirst
      integer lkntrl
      integer llevel
      integer lmessg
      integer lun(5)
      integer maxmes
      character ( len = * ) messg
      integer mkntrl
      integer nerr
      integer ni
      integer nmessg
      integer nr
      integer nunit
      real ( kind = 8 ) r1
      real ( kind = 8 ) r2
c     
c     Get flags
c     
      lkntrl = j4save ( 2, 0, .false. )
      maxmes = j4save ( 4, 0, .false. )
c     
c     Check for valid input
c     
      if ( 0<nmessg .and. nerr/=0 .and. -1<=level .and. level<=2 ) then
         go to 10
      end if
      
      if ( 0 < lkntrl ) then
         call xerprt('fatal error in...',17)
      end if
      
      call xerprt( 'XERROR -- invalid input', 23 )
      
      if ( 0 < lkntrl ) then
         call xerprt('job abort due to fatal error.',29)
      end if
      
      if ( 0 < lkntrl ) then
         call xersav ( ' ', 0, 0, 0, kdummy )
      end if
      
      call xerabt('XERROR -- invalid input',23)
      return
      
 10   continue
c     
c     Record the message.
c     
      junk = j4save(1,nerr,.true.)
      call xersav ( messg, nmessg, nerr, level, kount )
c     
c     Let user override
c     
      lfirst = messg
      lmessg = nmessg
      lerr = nerr
      llevel = level
      call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
c     
c     Reset to original values.
c     
      lmessg = nmessg
      lerr = nerr
      llevel = level
      lkntrl = max ( -2, min ( 2, lkntrl ) )
      mkntrl = abs ( lkntrl )
c     
c     Decide whether to print message
c     
      if ( llevel < 2 .and. lkntrl == 0 ) then
         go to 100
      end if
      
      if (((llevel == (-1)) .and. ( min ( 1, maxmes ) < kount ) ) 
     &     .or.((llevel == 0) .and. ( maxmes < kount )) 
     &     .or.((llevel == 1) .and. ( maxmes < kount ).and.(mkntrl==1) ) 
     &     .or.((llevel == 2) .and. ( max (1,maxmes) < kount ) ) ) then
      go to 100
      end if
      
      if ( 0 < lkntrl ) then
         
         call xerprt(' ',1)
         
         if ( llevel == -1 ) then
            
            call xerprt 
     &           ('this warning message will only be printed once.',57)
            
         end if
         
         if ( llevel == 0 ) then
            call xerprt ( 'warning in...', 13 ) 
         else if ( llevel == 1 ) then
            call xerprt ( 'recoverable error in...', 23 )
         else if ( llevel == 2 ) then
            call xerprt ( 'fatal error in...', 17 )
         end if
         
      end if
c     
c     Message
c     
      call xerprt(messg,lmessg)
      call xgetua(lun,nunit)
      isizei = log10 ( real (i1mach(9), kind=8 ) ) + 1.0D+00
      isizef = log10 ( real(i1mach(10), kind=8 )**i1mach(14) ) + 1.0D+00

      do kunit = 1, nunit
         
         iunit = lun(kunit)
         
         do i = 1, min ( ni, 2 )
            write (form,21) i,isizei
 21         format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
            if ( iunit == 0 ) then
               if (i == 1) write (*,form) i1
               if (i == 2) write (*,form) i2
            else
               if (i == 1) write (iunit,form) i1
               if (i == 2) write (iunit,form) i2
            end if
         end do
         
         do i = 1, min ( nr, 2 )
            write (form,23) i,isizef+10,isizef
 23         format('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')')
            if ( iunit == 0 ) then
               if ( i == 1 ) write (*,form) r1
               if ( i == 2 ) write (*,form) r2
            else
               if (i == 1) write (iunit,form) r1
               if (i == 2) write (iunit,form) r2
            end if
         end do
         
         if ( lkntrl <= 0 ) then
            go to 40
         end if
c     
c     error number
c     
         if ( iunit == 0 ) then
            write(*,30) lerr
         else
            write (iunit,30) lerr
         end if
 30      format (15h error number =,i10)
 40      continue

      end do
c     
c     Traceback
c     
 100  continue
      ifatal = 0
      if ((llevel == 2).or.((llevel==1) .and. (mkntrl==2))) then
         ifatal = 1
      end if
c     
c     quit here if message is not fatal
c     
      if ( ifatal <= 0 ) then
         return
      end if
      
      if ( lkntrl <= 0 .or. max ( 1, maxmes ) < kount ) then
         go to 120
      end if
c     
c     Print reason for abort
c     
      if ( llevel == 1 ) then
         call xerprt ('job abort due to unrecovered error.',35)
      end if
      
      if ( llevel == 2 ) then
         call xerprt('job abort due to fatal error.',29)
      end if
c     
c     Print error summary
c     
      call xersav ( ' ', -1, 0, 0, kdummy )
      
 120  continue
c     
c     Abort
c     
      if ( llevel == 2 .and. max ( 1, maxmes ) < kount ) then
         lmessg = 0
      end if
      
      call xerabt ( messg, lmessg )
      
      return
      end
c*****************************************************************************80




      subroutine xersav ( messg, nmessg, nerr, level, icount )
c*****************************************************************************80
c
cc XERSAV records that an error occurred.
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character ( len = * ) MESSG, as in XERROR.
c
c    Input, integer NMESSG, as in XERROR, except that, when NMESSG = 0,
c    the tables will be dumped and cleared; and when NMESSG < 0,
c    the tables will be dumped, but not cleared.
c
c    Input, integer NERR, as in XERROR.
c
c    Input, integer LEVEL, as in XERROR.
c
c    Output, integer ICOUNT, the number of times this message has
c    been seen, or zero if the table has overflowed and
c    does not contain this message specifically.
c    when nmessg=0, icount will not be altered.
c
      implicit none
      
      integer i
      integer i1mach
      integer icount
      integer ii
      integer iunit
      integer, save, dimension (10) :: kount=(/0,0,0,0,0,0,0,0,0,0/)
      integer, save :: kountx = 0
      integer kunit
      integer level
      integer, save, dimension ( 10 ) :: levtab
      integer lun(5)
      character ( len = 20 ) mes
      character ( len = * ) messg
      character ( len = 20 ), save, dimension ( 10 ) :: mestab
      integer nerr
      integer, save, dimension ( 10 ) :: nertab
      integer nmessg
      integer nunit
c     
c     Dump the table
c     
      if ( nmessg <= 0 ) then
         
         if ( kount(1) == 0 ) then
            return
         end if
c     
c     Print to each unit
c     
         call xgetua ( lun, nunit )
         
         do kunit = 1, nunit
            
            iunit = lun(kunit)
            
            if ( iunit == 0 ) then
               iunit = i1mach(4)
            end if
c     
c     Print table header
c     
            write ( iunit, 10 )
 10         format ('0          error message summary'/ 
     &           ' message start             nerr     level     count')
c     
c     print body of table
c     
            do i = 1, 10
               if ( kount(i) == 0 ) then
                  exit
               end if
               write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
 15            format (1x,a20,3i10)
            end do
c     
c     Print number of other errors
c     
            if ( kountx /= 0 ) then
               write (iunit,40) kountx
            end if
            
 40         format (41h0other errors not individually tabulated=,i10)
            write ( iunit, '(a)' ) ' '
         end do
c     
c     Clear the error tables
c     
         if ( nmessg == 0 ) then
            kount(1:10) = 0
            kountx = 0
         end if
         
         return
         
      end if
c     
c     process a message...
c     search for this message, or else an empty slot for this messg,
c     or else determine that the error table is full.
c     
      mes = messg
      
      do i = 1, 10
         
         ii = i
         
         if ( kount(i) == 0 ) then
            mestab(ii) = mes
            nertab(ii) = nerr
            levtab(ii) = level
            kount(ii)  = 1
            icount = 1
            return
         end if
         
         if ( mes /= mestab(i) ) then
            go to 90
         end if
         
         if (nerr /= nertab(i) ) then
            go to 90
         end if
         
         if (level /= levtab(i) ) then
            go to 90
         end if
         
         go to 100
         
 90      continue
         
      end do
c     
c     The table is full.
c     
      kountx = kountx + 1
      icount = 1
      return
c     
c     Message found in table
c     
 100  continue
      
      kount(ii) = kount(ii) + 1
      icount = kount(ii)
      
      return
      end
c*****************************************************************************80




      function d9lgmc ( x )
c*****************************************************************************80
c
cc D9LGMC computes the log gamma correction factor.
c
c  Discussion:
c
c    The routine computes the log gamma correction factor for 10 <= X 
c    so that
c
c      log ( gamma ( x ) ) = 
c        log ( sqrt ( 2 * pi ) ) + ( x - 0.5 ) * log ( x ) - x + d9lgmc ( x )
c
c  Parameters:
c
c    Input, real ( kind = 8 ) X, the argument of the log gamma function.
c    X must be at least 10.
c
c    Output, real ( kind = 8 ) D9LGMC, the correction.
c
      implicit none
      
      real ( kind = 8 ), save, dimension ( 6 ) :: algmcs = (/ 
     &     0.166638948045186D+00, -0.0000138494817606D+00,  
     &     0.0000000098108256D+00, -0.0000000000180912D+00, 
     &     0.0000000000000622D+00, -0.0000000000000003D+00 /)
      real ( kind = 8 ) arg
      real ( kind = 8 ) csevl
      real ( kind = 8 ) d9lgmc
      integer ( kind = 8 ) inits
      integer ( kind = 8 ), save :: nalgm = 0
      real ( kind = 8 ) x
      real ( kind = 8 ), save :: xbig = 0.0D+00
      real ( kind = 8 ), save :: xmax = 0.0D+00
      
      if ( nalgm == 0 ) then
         nalgm = inits ( algmcs, 6, epsilon ( algmcs ) )
         xbig = 1.0D+00 / sqrt ( epsilon ( xbig ) )
         xmax = exp ( min ( log ( huge ( xmax ) / 12.0D+00 ), 
     &        -log ( 12.0D+00 * tiny ( xmax ) ) ) )
      end if
      
      if ( x < 10.0D+00 ) then
         
         call xerror ( 'D9LGMC - 10 <= x required', 1, 2 )
         
      else if ( x < xbig ) then
         
         arg = 2.0D+00 * ( 10.0D+00 / x )**2 - 1.0D+00
         d9lgmc = csevl ( arg, algmcs, nalgm ) / x
         
      else if ( x < xmax ) then
         
         d9lgmc = 1.0D+00 / ( 12.0D+00 * x )
         
      else
         
         d9lgmc = 0.0D+00
         call xerror ( 'D9LGMC - X so big d9lgmc underflows', 2, 1)
         
      end if
      
      return
      end
c*****************************************************************************80




      function gamma ( x )
c*****************************************************************************80
c
cc GAMMA computes the gamma function.
c
c  Parameters:
c
c    Input, real X, the argument of the gamma function, which must not
c    be 0, -1, or any other negative integral value.
c
c    Output, real GAMMA, the value of the gamma function of X.
c
      implicit none
      
      real ( kind = 8 ) csevl
      real ( kind = 8 ) d9lgmc
      real ( kind = 8 ), save :: dxrel = 0.0D+00
      real ( kind = 8 ) gamma
      real ( kind = 8 ), parameter, dimension ( 23 ) :: gcs = (/ 
     &     0.008571195590989331D+00, 
     &     0.004415381324841007D+00, 
     &     0.05685043681599363D+00, 
     &     -0.004219835396418561D+00, 
     &     0.001326808181212460D+00, 
     &     -0.0001893024529798880D+00, 
     &     0.0000360692532744124D+00, 
     &     -0.0000060567619044608D+00, 
     &     0.0000010558295463022D+00, 
     &     -0.0000001811967365542D+00, 
     &     0.0000000311772496471D+00, 
     &     -0.0000000053542196390D+00, 
     &     0.0000000009193275519D+00, 
     &     -0.0000000001577941280D+00, 
     &     0.0000000000270798062D+00, 
     &     -0.0000000000046468186D+00, 
     &     0.0000000000007973350D+00, 
     &     -0.0000000000001368078D+00, 
     &     0.0000000000000234731D+00, 
     &     -0.0000000000000040274D+00, 
     &     0.0000000000000006910D+00, 
     &     -0.0000000000000001185D+00, 
     &     0.0000000000000000203D+00 /)
      integer i
      integer inits
      integer n
      integer, save :: ngcs = 0
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) sinpiy
      real ( kind = 8 ), parameter :: sq2pil = 0.91893853320467274D+00
      real ( kind = 8 ) x
      real ( kind = 8 ), save :: xmax = 0.0D+00
      real ( kind = 8 ), save :: xmin = 0.0D+00
      real ( kind = 8 ) y
c     
c     Initialize.  Find legal bounds for X, and determine the number of
c     terms in the series required to attain an accuracy ten times better
c     than machine precision.
c     
      if ( ngcs == 0 ) then
         ngcs = inits ( gcs, 23, 0.1D+00 * epsilon ( gcs ) )
         call gamlim ( xmin, xmax )
         dxrel = sqrt ( epsilon ( dxrel ) )
      end if
      
      y = abs ( x )
      if ( 10.0D+00 < y ) then
         go to 50
      end if
c     
c     Compute gamma(x) for abs ( x ) <= 10.0.  Reduce interval and
c     find gamma(1+y) for 0 <= y < 1. 
c     
      n = int ( x )
      
      if ( x < 0.0D+00 ) then
         n = n - 1
      end if
      
      y = x - real ( n, kind = 8 )
      n = n - 1
      gamma = 0.9375D+00 + csevl ( 2.0D+00 * y - 1.0D+00, gcs, ngcs )
      
      if ( n == 0 ) then
         
         return
         
      else if ( n < 0 ) then
         
         n = -n
         if ( x == 0.0D+00 ) then
            call xerror ( 'GAMMA - x is 0', 4, 2)
         end if
         
         if ( x<0.D0 .and. x + real(n-2,kind=8) == 0.0D+00 ) then
            call xerror (  'GAMMA - x is a negative integer', 4, 2 )
         end if
         
         if ( x < -0.5D+00 .and. 
     &        abs ( ( x - aint ( x - 0.5D+00) ) / x ) < dxrel ) then
            call xerror ( 
     &           'GAMMA < 1/2 precision since x nearly negative',
     &           1, 1 )
         end if
      
      do i = 1, n
         gamma = gamma / ( x + real ( i - 1, kind = 8 ) )
      end do
      
      return
      
      else
         
         do i = 1, n
            gamma = ( y + real ( i, kind = 8 ) ) * gamma
         end do
         
         return
         
      end if
c     
c     Compute gamma(x) for 10 < abs ( x ).  Recall y = abs ( x ).
c     
 50   continue
      
      if ( xmax < x ) then
         call xerror ( 'GAMMA - x so big gamma overflows', 3, 2)
      end if
      
      gamma = 0.0D+00
      
      if ( x < xmin ) then
         call xerror ( 'GAMMA - x so small gamma underflows', 2, 1)
         return
      end if
      
      gamma = exp( (y-0.5D0) * log(y) - y + sq2pil + d9lgmc(y) )

      if ( 0.0D+00 < x ) then
         return
      end if
      
      if ( abs ( ( x - aint ( x - 0.5D+00 ) ) / x ) < dxrel ) then
         call xerror ( 
     &        'GAMMA < half precision, x too near negative integer',
     &        1, 1)
      end if
      
      sinpiy = sin ( pi * y )
      
      if ( sinpiy == 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'GAMMA - Fatal error!'
         write ( *, '(a)' ) '  The argument X is a negative integer.'
         stop
      end if
      
      gamma = -pi / ( y * sinpiy * gamma )
      
      return
      end
c*****************************************************************************80




      function j4save ( iwhich, ivalue, iset )
c*****************************************************************************80
c
cc J4SAVE saves variables needed by the library error handling routines.
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer IWHICH, the index of the item desired.
c    1, the current error number.
c    2, the current error control flag.
c    3, the current unit number to which error messages are sent.  
c       (0 means use standard.)
c    4, the maximum times any message is printed (as set by xermax).
c    5, the number of units to which each error message is written.
c    6, the 2nd unit for error messages.
c    7, the 3rd unit for error messages.
c    8, the 4th unit for error messages.
c    9, the 5th unit for error messages.
c
c    Input, integer IVALUE, the value to be set for the IWHICH-th parameter,
c    if ISET is TRUE.
c
c    Input, logical ISET.
c    TRUE: the IWHICH-th parameter will be given the value, IVALUE.  
c
c    Output, integer J4SAVE, the old value of the IWHICH-th parameter.
c
      implicit none

      integer, save, dimension(9) :: iparam = (/ 0,2,0,10,1,0,0,0,0 /)
      logical iset
      integer ivalue
      integer iwhich
      integer j4save
      
      j4save = iparam(iwhich)
      
      if ( iset ) then
         iparam(iwhich) = ivalue
      end if
      
      return
      end
c*****************************************************************************80




      subroutine xerabt ( messg, nmessg )
c*****************************************************************************80
c
cc XERABT aborts program execution and prints an error message.
c
c  Discussion:
c
c    XERABT aborts the execution of the program.  The error message causing 
c    the abort is given in the calling sequence, in case one needs it for 
c    printing on a dayfile, for example.
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character ( len = * ) MESSG, the message to be processed, 
c    containing no more than 72 characters.
c
c    Input, integer NMESSG, the actual number of characters in MESSG.
c    If NMESSG is 0, no message is being supplied.
c
      implicit none
      
      character ( len = * ) messg
      integer nmessg
      
      if ( 0 < nmessg ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'XERABT - Termination after fatal error!'
         write ( *, '(a)' ) trim ( messg )
      end if
      
      stop
      end
c*****************************************************************************80




      subroutine xerprt ( messg, nmessg )
c*****************************************************************************80
c
cc XERPRT prints a message on each file indicated by xgetua.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character ( len = * ) MESSG, the message to be printed.
c
c    Input, integer NMESSG, the actual number of characters in MESSG.
c
      implicit none
      
      integer ichar
      integer iunit
      integer kunit
      integer last
      integer lenmes
      integer lun(5)
      character ( len = * ) messg
      integer nmessg
      integer nunit
c
c     Obtain unit numbers and write line to each unit
c
      call xgetua ( lun, nunit )
      
      lenmes = len ( messg )
      
      do kunit = 1, nunit
         
         iunit = lun(kunit)
         
         do ichar = 1, lenmes, 72
            last = min ( ichar+71 , lenmes )
            if ( iunit == 0 ) then
               write (*,'(1x,a)') messg(ichar:last)
            else
               write (iunit,'(1x,a)') messg(ichar:last)
            end if
         end do
         
      end do
      
      return
      end
c*****************************************************************************80




      subroutine xerctl ( messg1, nmessg, nerr, level, kontrl )
c*****************************************************************************80
c
cc XERCTL allows user control over handling of individual errors.
c
c  Discussion:
c
c    Allows user control over handling of individual errors.
c    Just after each message is recorded, but before it is
c    processed any further (i.e., before it is printed or
c    a decision to abort is made), a call is made to XERCTL.
c    If the user has provided his own version of XERCTL, he
c    can then override the value of KONTRL used in processing
c    this message by redefining its value.
c
c    KONTRL may be set to any value from -2 to 2.
c    The meanings for KONTRL are the same as in XSETF, except
c    that the value of KONTRL changes only for this message.
c    If KONTRL is set to a value outside the range from -2 to 2,
c    it will be moved back into that range.
c
c  Modified:
c
c    05 April 2007
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, character ( len = * ) MESSG1, the first word (only) of the error 
c    message.
c
c    Input, integer NMESSG, same as in the call to XERROR or XERRWV.
c
c    Input, integer NERR, same as in the call to XERROR or XERRWV.
c
c    Input, integer LEVEL, same as in the call to XERROR or XERRWV.
c
c    Input/output, integer KONTRL.  On input, the current value of the control 
c    flag as set by a call to XSETF.  On output, the new value of kontrl.  
c    If KONTRL is not defined, it will remain at its original value.
c    This changed value affects only the current occurrence of the current 
c    message.
c
      implicit none
      
      integer kontrl
      integer level
      character ( len = * ) messg1
      integer nerr
      integer nmessg
      
      return
      end
c*****************************************************************************80




      subroutine xgetua ( iunita, n )
c*****************************************************************************80
c
cc XGETUA returns the unit number(s) to which error messages are being sent.
c
c  Discussion:
c
c    XGETUA may be called to determine the unit number or numbers to which 
c    error messages are being sent.  These unit numbers may have been set 
c    by a call to XSETUN, or a call to XSETUA, or may be a default value.
c
c  Author:
c
c    Ron Jones
c
c  Reference:
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Technical Report SAND82-0800,
c    Sandia National Laboratories, 1982.
c
c    Ron Jones, David Kahaner,
c    XERROR, The SLATEC Error Handling Package,
c    Software: Practice and Experience,
c    Volume 13, Number 3, 1983, pages 251-257.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Output, integer IUNITA(N),  an array unit numbers, 
c    A value of zero refers to the default unit, as defined by the 
c    I1MACH machine constant routine.  Only IUNITA(1),..., IUNITA(N) are
c    defined by XGETUA.  The values of IUNITA(N+1),..., IUNITA(5) are 
c    not defined (for N < 5) or altered in any way by XGETUA.
c
c    Output, integer N, the number of units to which copies of the
c    error messages are being sent.  N will be in the range from 1 to 5.
c
      implicit none
      
      integer i
      integer index
      integer iunita(5)
      integer j4save
      integer n
      
      n = j4save ( 5, 0, .false. )
      
      do i = 1, n
         
         index = i+4
         if ( i == 1 ) then
            index = 3
         end if
         
         iunita(i) = j4save ( index, 0, .false. )
         
      end do
      
      return
      end
c*****************************************************************************80




      function csevl ( x, cs, n )
c*****************************************************************************80
c
cc CSEVL evaluates an N term Chebyshev series.
c
c  Modified:
c
c    15 April 2003
c
c  Reference:
c
c    R Broucke,
c    Algorithm 446, 
c    Communications of the ACM,
c    Volume 16, page 254, 1973.
c
c    Leslie Fox, Ian Parker,
c    Chebyshev Polynomials in Numerical Analysis,
c    Oxford Press, page 56.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) X, the argument at which the series is to be
c    evaluated.  X must satisfy -1.0 <= X <= 1.0.
c
c    Input, real ( kind = 8 ) CS(N), the array of N terms of a Chebyshev series.  
c    In evaluating CS, only half the first coefficient is summed.
c
c    Input, integer N, the number of terms in array CS.
c    N must be at least 1, and no more than 1000.
c
c    Output, real ( kind = 8 ) CSEVL, the value of the Chebyshev series.
c
      implicit none
      
      integer n
      
      real ( kind = 8 ) b0
      real ( kind = 8 ) b1
      real ( kind = 8 ) b2
      real ( kind = 8 ) cs(n)
      real ( kind = 8 ) csevl
      integer i
      real ( kind = 8 ) x
      
      if ( n < 1 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CSEVL - Fatal error!'
         write ( *, '(a)' ) '  Number of terms N is less than 1.'
         stop
      end if
      
      if ( 1000 < n ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'CSEVL - Fatal error!'
         write ( *, '(a)' ) '  The number of terms is more than 1000.'
         stop
      end if
      
      if ( x < -1.0D+00 .or. 1.0D+00 < x ) then
         write(*,'(a)') ' '
         write(*,'(a)') 'CSEVL - Fatal error!'
         write(*,'(a)') '  Input argument X outside interval [-1,1].'
         stop
      end if
      
      b1 = 0.0D+00
      b0 = 0.0D+00
      
      do i = n, 1, -1
         b2 = b1
         b1 = b0
         b0 = 2.0D+00 * x * b1 - b2 + cs(i)
      end do
      
      csevl = 0.5D+00 * ( b0 - b2 )
      
      return
      end
c*****************************************************************************80




      subroutine gamlim ( xmin, xmax )
c*****************************************************************************80
c
cc GAMLIM computes the minimum and maximum bounds for X in GAMMA(X).
c
c  Discussion:
c
c    GAMLIM calculates the minimum and maximum legal bounds for X in GAMMA(X).
c
c  Modified:
c
c    11 August 2001
c
c  Parameters:
c
c    Output, real ( kind = 8 ) XMIN, the minimum legal value of X in GAMMA(X).
c    Any smaller value might result in underflow.
c
c    Output, real ( kind = 8 ) XMAX, the maximum legal value of X in GAMMA(X).
c    Any larger value will cause overflow.
c
      implicit none
      
      real ( kind = 8 ) alnbig
      real ( kind = 8 ) alnsml
      logical converged
      integer i
      real ( kind = 8 ) xln
      real ( kind = 8 ) xmax
      real ( kind = 8 ) xmin
      real ( kind = 8 ) xold
      
      alnsml = log ( tiny ( alnsml ) )
      xmin = -alnsml
      
      converged = .false.
      
      do i = 1, 10
         
         xold = xmin
         xln = log ( xmin )
         
         xmin = xmin - xmin * 
     &        ( ( xmin + 0.5D+00 ) * xln - xmin - 0.2258D+00 + alnsml ) 
     &        / ( xmin * xln + 0.5D+00 )

         if ( abs ( xmin - xold ) < 0.005D+00 ) then
            converged = .true.
            exit
         end if
         
      end do
      
      if ( .not. converged ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'GAMLIM - Fatal error!'
         write ( *, '(a)' ) '  Unable to determine XMIN.'
         stop
      end if
      
      xmin = -xmin + 0.01D+00
      
      alnbig = log ( huge ( alnbig ) )
      xmax = alnbig
      
      converged = .false.
      
      do i = 1, 10
         
         xold = xmax
         xln = log ( xmax )
         xmax = xmax - xmax * 
     &        ( ( xmax - 0.5D+00 ) * xln - xmax + 0.9189D+00 - alnbig ) 
     &        / ( xmax * xln - 0.5D+00 )
         
         if ( abs ( xmax - xold ) < 0.005D+00 ) then
            converged = .true.
            exit
         end if
         
      end do
      
      if ( .not. converged ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'GAMLIM - Fatal error!'
         write ( *, '(a)' ) '  Unable to determine XMAX.'
         stop
      end if
      
      xmax = xmax - 0.01D+00
      xmin = max ( xmin, -xmax + 1.0D+00 )
      
      return
      end
c*****************************************************************************80




      function inits ( os, nos, eta )
c*****************************************************************************80
c
cc INITS estimates the order of an orthogonal series for a given accuracy.
c
c  Discussion:
c
c    Because this is a function, it is not possible to print out
c    warning error messages.  Therefore, if an error condition is
c    detected, a bogus value of INITS is returned.
c
c  Modified:
c
c    15 April 2003
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, real ( kind = 8 ) OS(NOS), the coefficients in the series.
c
c    Input, integer NOS, the number of coefficients.  NOS must be
c    at least 1, or an error condition arises.
c
c    Input, real ( kind = 8 ) ETA, the requested accuracy of the series.
c    Ordinarily, ETA will be chosen to be one-tenth machine precision.
c
c    Output, integer INITS, the order of the series guaranteeing the
c    given accuracy.  However, on error, INITS will be returned
c    as a negative huge number.
c
      implicit none
      
      integer nos
      
      real ( kind = 8 ) err
      real ( kind = 8 ) eta
      integer i
      integer ii
      integer inits
      real ( kind = 8 ) os(nos)
c     
c     Fatal error.  Number of coefficients less than 1.
c     
      if ( nos < 1 ) then
         inits = - huge ( i )
         return
      end if
      
      err = 0.0D+00
      
      i = 0
      
      do ii = 1, nos
         
         i = nos + 1 - ii
         err = err + abs ( os(i) )
         
         if ( eta < err ) then
            i = nos + 1 - ii
            exit
         end if
         
      end do
      
      if ( i == 0 ) then
         i = nos
         call xerror ( 'INITS - ETA may be too small', 1, 2)
      end if
      
      inits = i
      
      return
      end
c*****************************************************************************80
