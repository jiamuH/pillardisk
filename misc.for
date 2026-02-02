C*CGS ... physical constants in cgs units
C+
      REAL FUNCTION CGS( REQUEST )
*
* supplies value in cgs units of requested constant
C--
* 1988 Jul Keith Horne @ STScI
* 1989 Jul Keith Horne @ STScI - added MH, ME, E
* 1990 Jan KDH @ STScI - call double-precision version
	CHARACTER*(*) REQUEST
	REAL*8 DCGS
	CGS = DCGS( REQUEST )
	RETURN
	END

C*DCGS ... physical constants in double precision cgs units
C+
      REAL*8 FUNCTION DCGS( REQUEST )
*
* supplies value in cgs units of requested constant
C--
* Jul 1988 Keith Horne @ STScI
* Jul 1989 Keith Horne @ STScI - added MH, ME, E
* Jan 1990 KDH @ STScI - additional constants, double precision version
* 2001 Apr KDH @ Austin - MJY and ANGSTROM
* 2001 Jul KDH @ StAnd - g77 can't concatenate request
* 2013 Jul KDH @ StAnd - Thompson scattering cross-section
* 2014 Mar KDH @ StA - DAY
* 2014 Dec KDH @ StA - WIEN, ALPHA
* 2014 Dec KDH @ StA - STEFAN and THOMPSON aliases
	CHARACTER*(*)	REQUEST

* WARNING: long names must preceed otherwise identical short ones
* since otherwise wrong identification may occur

* pi
      IF( REQUEST .EQ. 'PI' ) THEN
        DCGS = 4.D0 * DATAN(1.D0)
* light speed c [cm/s]
      ELSE IF( REQUEST .EQ. 'C' ) THEN
        DCGS = 2.997925D10
* km/s ( KMS must precede K )
      ELSE IF( REQUEST .EQ. 'KMS' ) THEN
	DCGS = 1.E5

* Newton G [cm3/g/s2]
      ELSE IF( REQUEST .EQ. 'G' ) THEN
        DCGS = 6.673D-8
* Planck h [erg s]
      ELSE IF( REQUEST .EQ. 'H' ) THEN
        DCGS = 6.6262D-27
* Boltzman k [erg/K]
      ELSE IF( REQUEST .EQ. 'K' ) THEN
        DCGS = 1.3806D-16

* fine structure alpha = ( 2 pi e^2 ) / ( h c )
      ELSE IF( REQUEST .EQ. 'ALPHA' ) THEN
	DCGS = 1.d0 / 137.035965D0

* Thompson cross-section (8pi/3)(e4/m2c4) [cm2]
      ELSE IF( REQUEST .EQ. 'SIGMAT' .OR. REQUEST .EQ. 'THOMPSON' ) THEN
        DCGS = 0.66524E-24

* Stefan-Boltzman sigma = 2 pi^4 k^4 / ( 15 c^2 h^3 ) [erg/cm2/s/K4]
      ELSE IF( REQUEST .EQ. 'SIGMA' .or. REQUEST .EQ. 'STEFAN' ) THEN
        DCGS = 5.66956D-5

* Wien (lam_max of B_lambda) [cm Kelvin]
      ELSE IF( REQUEST .EQ. 'WIEN' ) THEN
	DCGS = 0.28979D0

* mJy [erg/cm2/s/Hz]
      ELSE IF( REQUEST .EQ. 'MJY' ) THEN
	DCGS = 1.D-26
* Angstroms [cm]
      ELSE IF( REQUEST .EQ. 'ANGSTROM' ) THEN
	DCGS = 1.D-8

* energies [erg = g cm2/s2]
      ELSE IF( REQUEST .EQ. 'RYD') THEN
        DCGS = 2.17992D-11
      ELSE IF( REQUEST .EQ. 'EV') THEN
        DCGS = 1.602192D-12

* Earth
      ELSE IF( REQUEST .EQ. 'MEARTH' ) THEN
        DCGS = 5.976D27
      ELSE IF( REQUEST .EQ. 'REARTH' ) THEN
        DCGS = 6.378164D8
* Sun
      ELSE IF( REQUEST .EQ. 'MSUN' ) THEN
        DCGS = 1.989D33
      ELSE IF( REQUEST .EQ. 'RSUN' ) THEN
        DCGS = 6.9599D10
      ELSE IF( REQUEST .EQ. 'LSUN' ) THEN
        DCGS = 3.826D33

* Jupiter
      ELSE IF( REQUEST .EQ. 'MJUP' ) THEN
	DCGS = 3.1783D2 * 5.976D27
      ELSE IF( REQUEST .EQ. 'RJUP' ) THEN
	DCGS = 7.13D9

* Eddington luminosity [erg/s/Msun]
*  Ledd = ( 4 pi G c m_p / sigma_T ) (M/Msun)
      ELSE IF( REQUEST .EQ. 'LEDSUN' ) THEN
        DCGS = 1.257D38

* Eddington accretion rate [g/s/Msun]
* ( Ledd / eta c^2 ) = ( 4 pi G m_p / eta c sigma_t ) [g/s/Msun]
* KDH: WHAT ETA IS ASSUMED ?
      ELSE IF( REQUEST .EQ. 'MDOTEDSUN' ) THEN
        DCGS = 1.399D17

* particle masses [g]
      ELSE IF( REQUEST .EQ. 'ME' ) THEN
        DCGS = 9.10956D-28
      ELSE IF( REQUEST .EQ. 'MP' ) THEN
        DCGS = 1.672661D-24
      ELSE IF( REQUEST .EQ. 'MH' ) THEN
        DCGS = 1.67352D-24
      ELSE IF( REQUEST .EQ. 'AMU' ) THEN
        DCGS = 1.660531D-24
* electron charge [ESU]
      ELSE IF( REQUEST .EQ. 'E') THEN
        DCGS = 4.80325D-10

* time
      ELSE IF( REQUEST .EQ. 'YR' .or. REQUEST .EQ. 'YEAR' ) THEN
        DCGS = 3.1556925D7
      ELSE IF( REQUEST .EQ. 'DY' .or. REQUEST .EQ. 'DAY' ) THEN
	DCGS = 24.D0 * 3600.D0
* distance
      ELSE IF( REQUEST .EQ. 'PC' ) THEN
        DCGS = 3.085678D18
      ELSE IF( REQUEST .EQ. 'AU' ) THEN
        DCGS = 1.49597D13

* help from the user
      ELSE
10      WRITE(*,'(A,A,A,$)')
     &	' Enter CGS value for "', REQUEST, '" : '
        READ(*,*,ERR=10 ) DCGS
      END IF

C	WRITE(*,*) 'CGS value for ", REQUEST, " is ', DCGS
      RETURN
      END		
	function extmag_mw( wave, ebmv )
	extmag_mw = extmag( wave, ebmv )
	return
	end

C*EXTMAG ... interstellar extinction function from Seaton(1979) and Nandy(1975)
C+
	FUNCTION EXTMAG( WAVE, EBMV )
*
* interstellar extinction law
*	1000 < lambda < 3704	Seaton(1979) MNRAS 187,73p.
*	3704 < lambda < 10,000	Nandy(1975) A+A 44, 195. (corrected to R=3.2)
*
* Input:
*	WAVE	R4 wavelength (Angstroms)
*	EBMV	R4 extinction parameter E(B-V) (mags)
*
* Output:
*	EXTMAG	= extinction (mags)
C--
* Mar 1986 Keith Horne @ STScI -- adapted from R.WADE routine SEATON
* May 1989 Keith Horne @ STScI -- improve extrapolation
* Jun 1991 KDH @ USM - minor changes
*
* COMMENTS FROM R.WADE SUBROUTINE:
c seaton's paper in m.n.r.a.s. vol 187, page 75p (1979).  
c the formulae are based on an adopted value of R = 3.20.
c
c note that seaton's representation of of the interstellar reddening law
c differs substantially from schild's representation (astron. j. 82, 339,
c table ii, 1977) in the region of overlap.  schild also adopted r = 3.20.
c
c for wavelengths > 3704 angstroms, the program interpolates 
c linearly in 1/lambda c in seaton's table 3.  
c for wavelengths < 3704 angstroms, the program uses the formulae 
c from seaton's table 2.  
c the formulae match at the endpoints of their respective intervals. 
c there is a mismatch of 0.009 mag/ebmv at nu=2.7 (lambda=3704 angstroms).
c seaton's tabulated value of 1.44 mags at 1/lambda = 1.1 may be in error;
c 1.64 seems more consistent with his other values.
c      
c wavelength range allowed is 0.1 to 1.0 microns.
c calls to the subroutine for wavelengths outside this range
c result in extrapolated values for the extinction being returned.
*
* sources:
*	lambda < 1000		same as lambda = 1000.
*	1000 < lambda < 3704	Seaton(1979) MNRAS 187,73p.
*	3704 < lambda < 10,000	Nandy(1975) A+A 44, 195. (corrected to R=3.2)
*	10000 < lambda		extrapolate linearly in 1/lam (can be improved)
* 2001 Jul KDH @ St.And - PARAMETER()

	PARAMETER ( NTABLE = 19 )
	REAL*4 XTABLE(NTABLE), ETABLE(NTABLE)

* tabulated inverse wavelengths (1/micron)
      DATA    XTABLE/ 0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
     &             1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
     &             2.2, 2.3, 2.4, 2.5, 2.6, 2.7	/

* tabulated extinctions at E(B-V)=1.
c      DATA    ETABLE/ 0., 1.36, 1.44, 1.84, 2.04, 2.24, 2.44,
      DATA    ETABLE/ 0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44,
     &          2.66, 2.88, 3.14, 3.36, 3.56, 3.77,
     &          3.96, 4.15, 4.26, 4.40, 4.52, 4.64/


	EXTMAG = 0.
	IF( WAVE.LE.0. ) RETURN
	IF( EBMV.EQ.0. ) RETURN
	X = 10000. / WAVE

* infrared - extend optical results linearly to 0 at 1/lam=0
      IF( X.LE.1.0 ) THEN
	EXTMAG = ETABLE(2) * X * X

* optical - interpolate linearly in magnitude vs 1/lam from Seaton's Table 3
      ELSE IF ( X.LT.2.7 ) THEN
   10	CALL LOOKUP( NTABLE, XTABLE, X, ILO, IHI, PART )
	EXTMAG = ETABLE(ILO)*(1.-PART) + ETABLE(IHI)*PART

* ultraviolet - use analytic formulae from Seaton's Table 2
      ELSE IF( X.LT.3.65 ) THEN
	DIFF = X - 4.6
	EXTMAG = 1.56 + 1.048*X + 1.01/( DIFF*DIFF + 0.280)

      ELSE IF( X.LT.7.14 ) THEN
	DIFF = X - 4.6
	EXTMAG = 2.29 + 0.848*X + 1.01/( DIFF*DIFF + 0.280)

      ELSE IF( X.LE.10. ) THEN
	EXTMAG = 16.17 + X*(-3.20 + 0.2975*X)

* far-uv - generally should not occur, Lyman edge is at 912 A.
      ELSE
	X = MIN( X, 50. )
	EXTMAG = 16.17 + X*(-3.20 + 0.2975*X)

      END IF

	EXTMAG = EBMV * EXTMAG
	RETURN
	END
	function extmag_smc( wave, ebmv )
* SMC Average dust extinction A(x)/A(V)
* from Table 7 in K.D.Gordon et al 2024 (ApJ 970, 51)
* 2025 Mar Keith Horne @ St Andrews
* In:	wave	r4 wavelength (Angstroms)
*	ebmv	r4 E(B-V) (mag)
* Out:	extmag_smc	r4 extinction in magnitudes
	parameter( nx = 30 )
	real*4 x( nx ), ax( nx ), sax( nx ), fx( nx )
* x = (micron/wavelength)
	data x / 0., 0.460,
     &	0.617, 0.812, 1.837, 2.284, 2.756,
     &	3.323, 3.462, 3.607, 3.759, 3.916,
     &	4.080, 4.251, 4.430, 4.615, 4.809,
     &	5.011, 5.221, 5.439, 5.668, 5.905,
     &	6.153, 6.411, 6.679, 6.959, 7.251,
     &	7.555, 7.872,        8.546       /
c     &	7.555, 7.872, 8.202, 8.546, 8.904/
* ax = A(x)/A(V)
	data ax / 0., 0.062,
     &	0.125, 0.324, 1.021, 1.349, 1.514,
     &	1.888, 1.968, 2.091, 2.174, 2.328,
     &	2.409, 2.571, 2.735, 2.894, 3.026,
     &	3.162, 3.272, 3.435, 3.615, 3.804,
     &	3.973, 4.167, 4.428, 4.654, 4.914,
     &	5.173, 5.511,        5.899       /
c     &	5.173, 5.511, 5.377, 5.899, 6.929/
* sax = 1-sigma uncertainty in ax
	data sax/ 0.001, 0.104,
     &	0.081, 0.069, 0.012, 0.012, 0.037,
     &	0.010, 0.009, 0.010, 0.009, 0.011,
     &	0.010, 0.015, 0.012, 0.014, 0.014,
     &	0.015, 0.016, 0.016, 0.018, 0.023,
     &	0.019, 0.024, 0.021, 0.021, 0.023,
     &	0.036, 0.031,        0.036       /
c     &	0.036, 0.031, 0.180, 0.036, 0.203/
* R = A(V)/E(B-V)
	data r, sigr / 3.02, 0.18 /
c	logical firstcall / .true. /
c      if( firstcall ) then
c	call splfit1( nx, x, ax, sax, fx, nx, ifail )
c	write(*,*) 'SMC spline fit IFAIL=', ifail
c	firstcall = .false.
c      end if


	xx = 1.e4 / wave
* interpolate
      if( xx .le. x( nx ) ) then
	call lookup( nx, x, xx, i1, i2, p )
	axx = terp1( nx, ax, i1, i2, p )
* extrapolate
      else
	nm = nx - 2
	dx = x( nx ) - x( nm )
	da = ax( nx ) - ax( nm )
	p = ( xx - x( nx ) ) / dx
	axx = ax( nx ) + da * p
      end if

c	call splcalc( nx, xx, axx, ifail )
c	write(*,*) 'lam', wave, ' X', xx, ' A', axx

	extmag_smc = axx * ebmv * r
	return
	end
*==================================
	FUNCTION RAN3(IDUM)
C			IMPLICIT REAL*4(M)
C			PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
	PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
	DIMENSION MA(55)
	DATA IFF /0/
	IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
		IFF=1
		MJ=MSEED-IABS(IDUM)
		MJ=MOD(MJ,MBIG)
		MA(55)=MJ
		MK=1
		DO I=1,54
			II=MOD(21*I,55)
			MA(II)=MK
			MK=MJ-MK
			IF(MK.LT.MZ)MK=MK+MBIG
			MJ=MA(II)
		ENDDO
		DO K=1,4
			DO I=1,55
				MA(I)=MA(I)-MA(1+MOD(I+30,55))
				IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
			ENDDO
		ENDDO
		INEXT=0
		INEXTP=31
		IDUM=1
	ENDIF
	INEXT=INEXT+1
	IF(INEXT.EQ.56)INEXT=1
	INEXTP=INEXTP+1
	IF(INEXTP.EQ.56)INEXTP=1
	MJ=MA(INEXT)-MA(INEXTP)
	IF(MJ.LT.MZ)MJ=MJ+MBIG
	MA(INEXT)=MJ
	RAN3=MJ*FAC
	RETURN
	END
************************************************
	function rang( avg, rms, iseed )
* Gaussian random number generator.
* Input:
*	avg	r4	mean value
*	rms	r4	standard deviation
*	iseed	i4	seed integer for ran( iseed ) 
* Output:
*	rang	r4	gaussian random number
*	iseed	i4	seed integer from ran( iseed )
*
* Box-muller transform of 2 independent random variables
* 1978 Peter Young @ Caltech -- original version GAUSS
* 1995 Keith Horne @ St-Andrews -- adapted from GAUSS.
* 2001 Sep KDH @ St-And - tidy up
* 2016 Jul KDH @ Iceland - ran -> rankdh (ran fails w gfortran)
	logical	newpair /.true./
	save y, u, newpair
	rang = avg
	if( rms .le. 0. ) return
      if ( newpair ) then
	r2 = 10.
      do while ( r2 .ge. 1. )
	x = 2. * rankdh( iseed ) - 1.
	y = 2. * rankdh( iseed ) - 1.
	r2 = x * x + y * y
      end do
	u = sqrt( -2. * alog( r2 ) / r2 )
	rang = x * u
      else
	rang = y * u
      end if
	newpair = .not. newpair
	rang = avg + rms * rang
	return
	end
*========================================
        subroutine rangtest( iseed )
* check avg and rms of gaussian random numbers
* 2016 Jul Keith Horne @ Iceland
        real*8 sum, sum2
        logical fail
* sample
        n = 10000
	sum = 0.d0
	sum2 = 0.d0
      do i = 1, n
        add = rang( 0., 1., iseed )
        sum = sum + add
        sum2 = sum2 + add * add
      end do
* mean
        avg = sum / n
	siga = 1. / sqrt( float( n ) )
	chia = avg / siga
* rms
	sum2 = ( sum2 - ( sum / n ) * sum ) / ( n - 1 )
        rms = dsqrt( sum2 )
	sigr = sqrt( 2. / ( n - 1 ) )
	chir = ( rms - 1. ) / sigr
* report
	write(*,*) 'avg', avg, ' =', chia, ' sigma'
	write(*,*) 'rms', rms, ' =', chir, ' sigma'
* test
        fail = max( abs( chia ), abs( chir ) ) .gt. 10.
      if( fail ) then
        write(*,*) '** RANG TEST FAILED  :('
	stop
      end if
        return
        end
*==================================
	real*4 function ranu( a, b, iseed )
* random number uniform on (a,b)
* In:	a	r4 lower limit
* 	b	r4 upper limit
*	iseed	i4 seed integer
* Out:  ranu	r4 random number uniform on (a,b)
*	iseed	i4 seed integer for use on next call
*
* 2020 May Keith Horne @ St.Andrews
	p = rankdh( iseed )
	ranu = a * ( 1. - p ) + p * b
	return
	end
*==================================
	real*4 function rankdh( iseed )
* random number uniform on (0-1)
* hi-level routines call this one, which can
* be used to select different lo-level generators.
* Input:
*	iseed	i4  repeat sequence if iseed < 0, otherwise ignore
*
* 2001 Jul Keith Horne @ St.Andrews
* 2006 Aug KDH @ St.A - use numerical recipes routine ran1
* 2008 Oct KDH @ St.A - try ran3
* 2016 Jan KDH @ Iceland - irand
	logical firstcall/.true./
	save firstcall
	data irand/0/
	is = iseed
      if( irand .ne. 0 ) then
* system-supplied
	ran = rand( is )
      else
* numerical recipes routine 
      if( firstcall ) then
	is = - iabs( is )
	firstcall = .false.
      end if
	ran = ran3( is )
      end if
	iseed = is
	rankdh = ran
	return
	end
************************************
	integer function len1( string )
* length up to last non-blank character
	character*(*) string
      do i = len( string ), 1, -1
      if( string( i:i ) .ne. ' ' ) then
	len1 = i
	return
      end if
      end do
	len1 = 0
	return
	end
C*WORD1 -- isolate first word of a text string
C+
	SUBROUTINE WORD1( STRING, L1, L2 )
*
* Isolates the first word of a character string 
*
* Input:
*	STRING	= character string
* Output:
*	L1,L2	= Inclusive index limits of first word
C--
* Jan 1986 Keith Horne @ STScI
* Jan 1989 KDH @ STScI - revise to allow BLANK or TAB separators
* Oct 2014 KDH @ StAnd - revise to allow comma separators
	CHARACTER*(*) STRING
	CHARACTER*1 TAB, BLANK, TEST, comma
	comma = ','
	TAB = CHAR( 9 )
	BLANK = ' '
	LMAX = LEN1( STRING )
      DO I1 = 1, LMAX
	TEST = STRING(I1:I1)
      IF(     TEST .NE. BLANK
     &	.AND. TEST .NE. TAB
     &	.and. test .ne. comma ) THEN
	L1 = I1
      DO I2 = L1 + 1, LMAX
	TEST = STRING(I2:I2)
      IF(    TEST .EQ. BLANK
     &	.OR. TEST .EQ. TAB
     &	.or. test .eq. comma ) THEN
	L2 = I2 - 1
	RETURN
      END IF
      END DO
	L2 = LMAX
	RETURN
      END IF
      END DO
	L1 = 0
	L2 = 0
	RETURN
	END
***********************************************
	SUBROUTINE UPCASE( TEXTIN, TEXTOUT )
* replacement for VMS library routine STR$UPCASE( TEXTOUT, TEXTIN )
	CHARACTER*(*) TEXTIN
	CHARACTER*(*) TEXTOUT
	TEXTOUT = TEXTIN
	JADD = ICHAR('A') - ICHAR('a')
      DO I=1,LEN1(TEXTOUT)
	JTEST = ICHAR(TEXTOUT(I:I) )
      IF( JTEST .GE. ICHAR('a') .AND. JTEST .LE. ICHAR('z') ) THEN
	TEXTOUT(I:I) = CHAR( JTEST + JADD )
      END IF
      END DO
	RETURN
	END
*************************
	subroutine append( str1, str2 )
* append two strings with a space in between
* 2003 Jun Keith Horne @ St.Andrews
	character*(*) str1, str2
	l1 = len1( str1 )
	l2 = len1( str2 )
      if( l1 .le. 0 ) then
	str1 = str2( : l2 )
      else
	str1 = str1( : l1 ) // ' ' // str2( : l2 )
      end if
	return
	end
****************************************
	subroutine append_i( str, i )
* append integer to string
* 2003 Jun Keith Horne @ St.Andrews
	character*(*) str
	character*40 number
	write( number, * ) i
	call word1( number, l1, l2 )
	call append( str, number( l1 : l2 ) )
	return
	end
****************************************
	subroutine append_f( str, x, fmt )
* append formatted real to string
* 2003 Jun Keith Horne @ St.Andrews
* 2008 Aug KDH @ St.A - add '*' format
	character*(*) str, fmt
	character*40 number
      if( fmt(1:1) .eq. '*' ) then
	write( number, * ) x
      else
	write( number, fmt ) x
      end if
	call word1( number, l1, l2 )
	call append( str, number( l1 : l2 ) )
	return
	end
****************************************
	subroutine append_d( str, x, fmt )
* append formatted real*8 to string
* 2004 Feb Keith Horne @ St.Andrews
	real*8 x
	character*(*) str, fmt
	character*40 number
	write( number, fmt ) x
	call word1( number, l1, l2 )
	call append( str, number( l1 : l2 ) )
	return
	end

****************************************
	subroutine append_n( title, val, idig )
* append val to title with idig significant digits
* 2008 Jun Keith Horne @ St.Andrews
* 2008 Jul KDH @ SAAO - idig
	character*(*) title
	character*20 fmt
	i = max( 2, idig )
	a= abs( val ) / 10**(i)
      if( a .eq. 0.0 ) then
	call append_f( title, val, '(f10.1)' )
      else if( a .lt. 9.9e-8) then
	fmt = '(1pe'
	call append_i( fmt, i + 7 )
	call append( fmt, '.' )
	call append_i( fmt, i - 1 )
	call append( fmt, ')' )
	call nospace( fmt )
	call append_f( title, val, fmt )
      else if( a .lt. 9.9e-7) then
	call append_f( title, val, '(f15.6)' )
      else if( a .lt. 9.9e-6) then
	call append_f( title, val, '(f15.5)' )
      else if( a .lt. 9.9e-5 ) then
	call append_f( title, val, '(f15.4)' )
      else if( a .lt. 9.9e-4 ) then
	call append_f( title, val, '(f15.3)' )
      else if( a .lt. 9.9e-3 ) then
	call append_f( title, val, '(f15.2)' )
      else if( a .lt. 0.099 ) then
	call append_f( title, val, '(f15.1)' )
      else
	call append_i( title, nint( val ) )
      end if
	return
	end

********************************************
        subroutine append_datsig( title, dat, sig, kfig )
* Input
*       title   c* string
*       dat     r4 data
*       sig     r4 sigma
*       kfig    i4 number of significant figures on sigma
* Output
*       title   c* string // data(sigma)
* 2004 Feb Keith Horne @ St-Andrews
* 2013 May KDH @ StA trap sig not positive
* 2017 Sep KDH @ StA revise for very small/large numbers
* 2022 Aug KDH @ StA use +/- symbol when advantageous
* 2022 Sep KDH @ StA fix bug (restrict final nospace)
* 2025 Jan KDH @ Crail revise final nospace
* 2025 Jan KDH @ SrA fix bug converting to X+/-Y
        character*(*) title
        real*8 hjd, big
	dmx = 1000
	dmn = 0.0001

	l0 = max( 0, len1( title ) )

      if( sig .le. 0. ) then
	call append_n( title, dat, max( 2, kfig ) )

      else if( abs(dat) .ge. dmn .and. abs( dat ) .le. dmx ) then
        hjd = dat
        call append_hjdsig( title, hjd, sig, kfig )

      else

* divide by a power of 10
c        mpow = nint( log10( abs( dat * 2. ) ) )
        mpow = nint( alog10( abs( dat * 2. ) ) )
        ipow = nint( alog10( sig * 2. ) )
        kpow = ipow - kfig

	kpow = int( alog10( max( abs(dat), sig ) ) ) - 1
        div = 10. ** kpow
	hjd = dat / div

	call append_hjdsig( title, hjd, sig/div, kfig )
* exponent
      if( kpow .ne. 0 ) then
	l = len1( title ) + 1
	title(l:l) = 'E'
	call append_i( title, kpow )
	call nospace( title( l : len1(title) ) )
      end if
      end if

* use +/- symbol where advantageous
*  find left and right parentheses
	l1 = l0 + index( title(l0+1:), '(' )
	l2 = l0 + index( title(l0+1:), ')' )
* these embrace the uncertainty
      if( l1 .gt. l0 .and. l2 .gt. l1 ) then
* look for decimal points
	k1 = l0 + index( title(l0+1:), '.' )
	k2 = l1 + index( title(l1+1:), '.' )
	write(*,*) '12345678901234567890'
	write(*,*) title(:len1(title))
	write(*,*) l1, k1, l2, k2
* if no decimal points
      if( ( k1 .le. l0 .and. k2 .le. l1 ) 
* or same number of digits follow both decimal points
     &	.or. ( k1 .gt. l0 .and. k2 .gt. l1
     &	.and. ( l1 - k1 .eq. l2 - k2 ) ) ) then
* replace "(" by "+/-" and omit ")".
	title(l2:l2) = ' '
	title = title(:l1-1) // '\(2233)' // title(l1+1:)
* embrace (X+/-Y) before exponent
* 2025-Jan-16 fix bug
c	l = index( title, 'E' )
	l = index( title(l0+1:), 'E' )
	if( l .gt. 0 ) title = title(:l0) // '(' // title(l0+1:l-1) // ')' // title(l:)
* restrict nospace to the new bits on the tail
c	call nospace( title(l0+1:) )
	call nospace( title(l0+2:) )

      end if
      end if

        return
        end

********************************************
        subroutine append_hjdsig( title, hjd, sig, kfig )
* Input
*       title   c* string
*       hjd     r8 data
*       sig     r4 sigma
*       kfig    i4 number of significant figures on sigma
* Output
*       title   c* string // data(sigma)
* 2004 Feb Keith Horne @ St-Andrews
* 2012 Dec KDH @ StA omit E0, skip if sig not positive
* 2013 May KDH @ StA trap sig not positive
        character*(*) title
        real*8 hjd
        character*50 fmt

* KDH : DEBUG REPORT
c	write(*,*) '** append_hjdsig(', hjd, sig, kfig, ' )'

      if( sig .le. 0. ) then
	call append_n( title, real(hjd), max( 2, kfig ) )
	return
      end if

        l1 = len1( title )
        mpow = nint( log10( abs( hjd * 2. ) ) )
        ipow = nint( alog10( sig * 2. ) )
        kpow = ipow - kfig
        div = 10. ** kpow

      if( kpow .lt. 0 ) then
        fmt = '(f20.'
        call append_i( fmt, -kpow )
        call append( fmt, ')' )
        call nospace( fmt )
        call append_d( title, hjd, fmt )

      if( sig .gt. 0. ) then
        sigdiv = sig / div
        call append( title, '(' )
      if( sig .lt. 1. ) then
        call append_i( title, nint( sigdiv ) )
      else
        call append_f( title, sig, fmt )
      end if
        call append( title, ')' )
      end if

      else if( kpow .gt. 3 ) then
        lpow = mpow
        if( iabs( ipow ) .gt. iabs( mpow ) ) lpow = ipow
        div = 10. ** lpow
        fmt = '(f20.'
        call append_i( fmt, kfig )
        call append( fmt, ')' )
        call nospace( fmt )
        call append_d( title, hjd / div, fmt )

      if( sig .gt. 0. ) then
        call append( title, '(' )
        call append_f( title, sig / div, fmt )
        call append( title, ')' )
      if( lpow .ne. 0 ) then
        call append( title, 'E' )
        call append_i( title, lpow )
      end if
      end if

      else
        idiv = 10 ** kpow
        call append_i( title, nint( hjd / idiv ) * idiv )
      if( sig .gt. 0. ) then
        call append( title, '(' )
        call append_i( title, nint( sig / idiv ) * idiv )
        call append( title, ')' )
      end if

      end if

* KDH : IS THIS THE CULPRIT?
c	write(*,*) title(:len1(title))
c	write(*,*) ' l1', l1, title(l1+2:len1(title))

        call nospace( title(l1+2:) )

c	write(*,*) title(:len1(title))
        return
        end

****************************************
	subroutine nospace( str )
* remove spaces from string
* 2003 Aug Keith Horne @ St.Andrews
	character*(*) str
	n = len1( str )
	if( n .lt. 1 ) return
	j = 0
      do i=1,n
      if( str( i : i ) .ne. ' ' ) then
	j = j + 1
	str( j : j ) = str( i : i )
      end if
      end do
	str = str( : j )
	return
	end
****************************************
	subroutine notail( str, chr )
* remove specified character from end of string
* 2003 Dec Keith Horne @ St.Andrews
* 2024-May-07 KDH @ St.A - avoid testing str(0:0)
	character*(*) str, chr
	n = len1( str )
	if( n .lt. 1 ) return
c     do while( n .gt. 0 .and. str(n:n) .eq. chr(1:1) )
      do while( n .gt. 0 )
	if( str( n : n ) .ne. chr( 1 : 1 ) ) return
	n = n - 1
	if( n .gt. 0 ) str = str(:n)
      end do
	return
	end
*------------------------------
	subroutine noeqbl( title )
* change '= ' to '='
* 2006 Feb Keith Horne @ St.Andrews
	character*(*) title
    1	l = index( title, '= ' )
	if( l .le. 0 ) return
	title = title(:l) // title(l+2:)
	goto 1
	end
*------------------------------
	subroutine nogreek( title )
* change pgplot greeks to text
* 2007 Aug Keith Horne @ LaSilla
	character*(*) title
	parameter ( nfix = 15 )
	character*3 old( nfix )
	character*5 new( nfix )
	data old/
     &		'\ga',	'\gb',	'\gg',	'\gg',	'\ge'
     &	,	'\gh',	'\gy',	'\gm',	'\gn',	'\gp'
     &	,	'\gr',	'\gs',	'\gx',	'\gD',	'\gW'
     &	/
	data new/
     &		'alp',	'bet',	'gam',	'del',	'eps'
     &	,	'theta', 'eta',	'mu',	'nu',	'pi'
     &	,	'rho',	'sig',	'chi',	'Del',	'Omega'
     &	/

      do k=1,nfix
	i = 1
      do while( i .gt. 0 )
	i = index( title, old(k) )
      if( i .gt. 0 ) then
	n = len1( new(k) )
	title = title(:i-1) // new(k)(:n) // title(i+3:)
      end if
      end do
      end do

	return
	end
*------------------------------
	subroutine noupdn( title )
* change pgplot sub and superscripts to text
* 2007 Aug Keith Horne @ LaSilla
	character*(*) title
	m = 1
      do while( m .gt. 0 )
	i = index( title, '\u' )
	j = index( title, '\d' )
      if( i .gt. 0 .and. j .gt. i ) then
	title = title(:i-1) // '^' // title(i+2:j-1) // title(j+2:)
      else if( j .gt. 0 .and. i .gt. j ) then
	title = title(:j-1) // '_' // title(j+2:i-1) // title(i+2:)
      else
	m = 0
      end if
      end do
	return
	end

*=============================================
        subroutine append_kmb( string, big, n )
* append big with units K=1e3, M=1e6, B=1e9
* In:   string  c* character string
*       big     r4 large number
*       n       i4 number of digits of precision
* Out:  label   c* appended string
* 2020 Apr Keith Horne @ St.Andrews
        character*(*) string
        x = big
        a = abs( x )
      if( a .lt. 1.e4 ) then
        call append_i( string, nint( x ) )
      else if( a .lt. 1.e6 ) then
        call append_n( string, x / 1.e3, n )
        call append( string, 'K' )
      else if( a .lt. 1.e9 ) then
        call append_n( string, x / 1.e6, n )
        call append( string, 'M' )
      else
        call append_n( string, x / 1.e9, n )
        call append( string, 'B' )
      end if
        return
  	end
	SUBROUTINE INQ( MESSAGE, DEFAULT, REPLY )
*
* get character string from user
*
* input:
*	MESSAGE	C* message to be sent to user
*	DEFAULT C* default reply
* output:
*	REPLY	C* user reply to the message
*
* 1988 Mar KDH @ STScI
* 2001 Jul KDH @ St.And - avoid char*(*) concatenation
	CHARACTER*(*) MESSAGE
	CHARACTER*(*) DEFAULT
	CHARACTER*(*) REPLY
	CHARACTER*256 LOCAL

      IF( DEFAULT.EQ.' ' ) THEN
	WRITE( *, '(1X,A,A,$)' ) MESSAGE( :LEN1(MESSAGE) ), ' '
      ELSE
	WRITE( *, '(1X,A,A,A,A,$)' ) MESSAGE( :LEN1(MESSAGE) ),
     #		' [' , DEFAULT( :LEN1(DEFAULT) ) , '] '
      END IF
	READ( *, '(A)' ) LOCAL
      IF( LOCAL.EQ.' ' ) THEN
	REPLY = DEFAULT
      ELSE
	REPLY = LOCAL
      END IF
	RETURN
	END
*===========================================================
	SUBROUTINE INQI( MESSAGE, IDEF, IVAL )
*
* get integer from the user
*
* input:
*	MESSAGE	C* message
*	IDEF 	I4 default value
* output:
*	IVAL	I4 returned value
*
* 1988 Mar KDH @ STScI
* 2002 Jul KDH @ St.And - add end=10
*
	CHARACTER*(*) MESSAGE
	CHARACTER*20 NUMBER
   10	WRITE( NUMBER, * ) IDEF
	CALL WORD1( NUMBER, L1, L2 )
	CALL INQ( MESSAGE, NUMBER(L1:L2), NUMBER )
	READ( NUMBER, *, ERR=10, END=10 ) NEW
	IVAL = NEW
	RETURN
	END
*===========================================================
	SUBROUTINE INQR( MESSAGE, DEF, VAL )
*
* get real value from the user
*
* input:
*	MESSAGE	C* message
*	DEF 	R4 default value
* output:
*	VAL	R4 returned value
*
* 1988 Mar KDH @ STScI
* 2002 Jul KDH @ St.And - add end=10
*
	CHARACTER*(*) MESSAGE
	CHARACTER*20 NUMBER
   10	WRITE( NUMBER, * ) DEF
	CALL INQ( MESSAGE, NUMBER, NUMBER )
	READ( NUMBER, *, ERR=10, END=10 ) WANT
	VAL = WANT
	RETURN
	END
*===========================================================
	SUBROUTINE INQD( MESSAGE, DEF, VAL )
*
* get real*8 value from the user
*
* input:
*	MESSAGE	C* message
*	DEF 	R8 default value
* output:
*	VAL	R8 returned value
*
* 1988 Mar KDH @ STScI
* 2002 Jul KDH @ St.And - add end=10
* 2004 Feb KDH @ St.And - adapt from INQR
* 2012 Nov KDH @ St.And - char*20->char*30 number
	REAL*8 DEF, VAL, WANT
	CHARACTER*(*) MESSAGE
	CHARACTER*30 NUMBER
   10	WRITE( NUMBER, * ) DEF
	CALL INQ( MESSAGE, NUMBER, NUMBER )
	READ( NUMBER, *, ERR=10, END=10 ) WANT
	VAL = WANT
	RETURN
	END
*===========================================================
	SUBROUTINE INQ2R( NAME, R1, R2 )
* get 2 reals from the user
* input:
*	NAME	C* name
*	R1 	R4 default value
*	R2 	R4 default value
* output:
*	R1	R4 returned value
*	R2	R4 returned value
* 2002 Jul KDH @ St.And
* 2008 Jul KDH @ SAAO - read v1, then v1,v2
	CHARACTER*(*) NAME
	CHARACTER*256 REPLY
   10	WRITE(*,*) 'Old ', NAME( :LEN1(NAME) ), ' : ', R1, R2
	WRITE(*,'(A,A,A,$)') ' New ', NAME( :LEN1(NAME) ), ' : '
	READ( *, '(A)' ) REPLY
      IF( REPLY .NE. ' ' ) THEN
	READ( REPLY ,*, ERR=10, END=20 ) V1
	R1 = V1
	READ( REPLY ,*, ERR=10, END=20 ) V1, V2
	R2 = V2
      END IF
   20	RETURN
	END
*===========================================================
	SUBROUTINE INQ2D( NAME, R1, R2 )
* get 2 real*8 values from the user
* In:	NAME	C* name
*	R1 	R8 default value
*	R2 	R8 default value
* Out:	R1	R8 returned value
*	R2	R8 returned value
* 2002 Jul KDH @ St.And
* 2008 Jul KDH @ SAAO - read v1, then v1,v2
* 2020 Jul KDH @ ST.And - adapt from INQ2R
	REAL*8 R1, R2, V1, V2
	CHARACTER*(*) NAME
	CHARACTER*256 REPLY
   10	WRITE(*,*) 'Old ', NAME( :LEN1(NAME) ), ' : ', R1, R2
	WRITE(*,'(A,A,A,$)') ' New ', NAME( :LEN1(NAME) ), ' : '
	READ( *, '(A)' ) REPLY
      IF( REPLY .NE. ' ' ) THEN
	READ( REPLY ,*, ERR=10, END=20 ) V1
	R1 = V1
	READ( REPLY ,*, ERR=10, END=20 ) V1, V2
	R2 = V2
      END IF
   20	RETURN
	END
*===========================================================
	SUBROUTINE INQ3R( NAME, R1, R2, R3 )
* get 3 reals from the user
* input:
*	NAME	C* name
*	R1 	R4 default value
*	R2 	R4 default value
*	R3 	R4 default value
* output:
*	R1	R4 returned value
*	R2	R4 returned value
*	R3	R4 returned value
* 2002 Jul KDH @ St.And
* 2008 Jul KDH @ SAAO - read v1, then v1,v2
* 2008 Jul KDH @ StA - adapt from INQ2R
	CHARACTER*(*) NAME
	CHARACTER*256 REPLY
   10	WRITE(*,*) 'Old ', NAME( :LEN1(NAME) ), ' : ', R1, R2, R3
	WRITE(*,'(A,A,A,$)') ' New ', NAME( :LEN1(NAME) ), ' : '
	READ( *, '(A)' ) REPLY
      IF( REPLY .NE. ' ' ) THEN
	READ( REPLY ,*, ERR=10, END=20 ) V1
	R1 = V1
	READ( REPLY ,*, ERR=10, END=20 ) V1, V2
	R2 = V2
	READ( REPLY ,*, ERR=10, END=20 ) V1, V2, V3
	R3 = V3
      END IF
   20	RETURN
	END
*===========================================================
	SUBROUTINE INQ2I( NAME, I1, I2 )
* get 2 integers from the user
* input:
*	NAME	C* name
*	I1 	I4 default value
*	I2 	I4 default value
* output:
*	I1	I4 returned value
*	I2	I4 returned value
* 2002 Jul KDH @ St.And
* 2008 Jul KDH @ SAAO - read j1, then j1,j2
	CHARACTER*(*) NAME
	CHARACTER*256 REPLY
   10	WRITE(*,*) 'Old ', NAME( :LEN1(NAME) ), ' : ', I1, I2
	WRITE(*,'(A,A,A,$)') ' New ', NAME( :LEN1(NAME) ), ' : '
	READ( *, '(A)' ) REPLY
      IF( REPLY .NE. ' ' ) THEN
	read( reply, *, err=10, end=20 ) j1
	i1 = j1
	READ( REPLY ,*, ERR=10, END=20 ) J1, J2
	I2 = J2
      END IF
   20	RETURN
	END

*******************************************
	subroutine mnmx( n, x, xmn, xmx )
* find min and max values
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	xmn	r4 minimum value
*	xmx	r4 maximum value
	real*4 x(*)
	xmn = x(1)
	xmx = xmn
      do i=1,n
	xmn = min( xmn, x(i) )
	xmx = max( xmx, x(i) )
      end do
	return
	end
*******************************************
	subroutine mnmx_index( n, x, imn, imx )
* find index of min and max values
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	imn	r4 index of first min value
*	imx	r4 index of last max value
* 2018 Oct Keith Horne @St.Andrews
	real*4 x(*)
	imn = 1
	imx = 1
      do i=1,n
	if( x(i) .lt. x(imn) ) imn = i
	if( x(i) .ge. x(imx) ) imx = i
      end do
	return
	end
*******************************************
	function minpos( n, x )
* find index of min value
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	minpos	r4 index of first min value
* 2018 Oct Keith Horne @St.Andrews
	real*4 x(*)
	m = 1
      do i=1,n
	if( x(i) .lt. x(m) ) m = i
      end do
	minpos = m
	return
	end
*******************************************
	function maxpos( n, x )
* find index of max value
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	maxpos	r4 index of last max value
* 2018 Oct Keith Horne @St.Andrews
	real*4 x(*)
	m = 1
      do i=1,n
	if( x(i) .gt. x(m) ) m = i
      end do
	maxpos = m
	return
	end
*******************************************
	subroutine mnmxpos( n, x, xmn, xmx )
* find min and max of positive values
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	xmn	r4 minimum value
*	xmx	r4 maximum value
* 2004 Apr Keith Horne @ St.And
* 2012 Sep KDH @ StA fix bug when x(1) < 0
	real*4 x(*), xmn, xmx
	npos = 0
      do i=1,n
      if( x(i) .gt. 0. ) then
	npos = npos + 1
      if( npos .eq. 1 ) then
	xmn = x(i)
	xmx = xmn
      else
	xmn = min( xmn, x(i) )
	xmx = max( xmx, x(i) )
      end if
      end if
      end do
      if( npos .lt. 1 ) then
	xmn = 0.
	xmx = 0.
      end if
	return
	end
*******************************************
	subroutine mnmxi( n, i, imn, imx )
* find min and max values for integer array 
* Input:
*	n	i4 number of data values
*	i(n)	i4 data values
* Output:
*	imn	r4 minimum value
*	imx	r4 maximum value
* 2004 Feb Keith Horne @ St-Andrews
	integer*4 i( * )
	imn = i( 1 )
	imx = imn
      do k=1,n
	imn = min( imn, i( k ) )
	imx = max( imx, i( k ) )
      end do
	return
	end

*******************************************
	subroutine dmnmx( n, x, xmn, xmx )
* find min and max of real*8 values
* Input:
*	n	i4 number of data values
*	x(n)	r8 data values
* Output:
*	xmn	r8 minimum value
*	xmx	r8 maximum value
	real*8 x(*), xmn, xmx
	xmn = x(1)
	xmx = xmn
      do i=1,n
	xmn = min( xmn, x(i) )
	xmx = max( xmx, x(i) )
      end do
	return
	end

*******************************************
	subroutine imnmx( n, x, imn, imx )
* find indices of min and max values
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	imn	r4 index of minimum value
*	imx	r4 index of maximum value
* 2012 Apr Keith Horne @ StA
	real*4 x( * )
	imn = 1
	imx = 1
      do i=1,n
	if( x( i ) .lt. x( imn ) ) imn = i
	if( x( i ) .gt. x( imx ) ) imx = i
      end do
	return
	end
*******************************************
	subroutine imnmxpos( n, x, imn, imx )
* find indices of min and max of positive values
* Input:
*	n	i4 number of data values
*	x(n)	r4 data values
* Output:
*	imn	r4 index of minimum positive value
*	imx	r4 index of maximum positive value
* 2012 Apr Keith Horne @ StA
	real*4 x( * )
	imn = 0
	imx = 0
	np = 0
      do i=1,n
      if( x( i ) .gt. 0. ) then
	np = np + 1
      if( np .eq. 1 ) then
	imn = i
	imx = i
      else
	if( x( i ) .lt. x( imn ) ) imn = i
	if( x( i ) .gt. x( imx ) ) imx = i
      end if
      end if
      end do
	return
	end
*******************************************
	subroutine imnmxi( n, i, imn, imx )
* find indeces of min and max values for integer array 
* Input:
*	n	i4 number of data values
*	i(n)	i4 data values
* Output:
*	imn	r4 index of minimum value
*	imx	r4 index of maximum value
* 2012 Apr Keith Horne @ StA
	integer*4 i( * )
	imn = 1
	imx = 1
      do k=1,n
	if( i( k ) .lt. i( imn ) ) imn = k
	if( i( k ) .gt. i( imx ) ) imx = k
      end do
	return
	end

*******************************************
	subroutine imnmx8( n, x, imn, imx )
* find indices of min and max of real*8 values
* Input:
*	n	i4 number of data values
*	x(n)	r8 data values
* Output:
*	imn	i4 index of minimum value
*	imx	i4 index of maximum value
* 2012 Apr Keith Horne @ StA
	real*8 x( * )
	imn = 1
	imx = 1
      do i=1,n
	if( x( i ) .lt. x( imn ) ) imn = i
	if( x( i ) .gt. x( imx ) ) imx = i
      end do
	return
	end
******************************************
	subroutine avgrms( n, x, bar, rms )
* unweighted mean and rms of n data points x(i)
* Input:
*	n	i4 number of data
*	x(n)	r4 data
* output:
*	bar	r4 mean value
*	rms	r4 root-mean-square
* Running update method (B.P.Welford 1962)
* Seems to be more robust to round-off errors.
* 2011 Feb Keith Horne @ St.Andrews
	real*4 x(*)
	real*8 sum, sum2, add
      if( n .le. 1 ) then
	bar = x(1)
	rms = 0.
	return
      end if
	sum1 = x(1)
	sum2 = 0.d0
      do i=2,n
	add = dble( x(i) ) - sum1
	sum1 = sum1 + add / i
	sum2 = sum2 + add * ( dble( x(i) ) - sum1 )
      end do
	bar = sum1
	var = sum2 / max( 1, n - 1 )	
	rms = sqrt( var )
	return
	end

******************************************
	subroutine avgrms1( n, x, bar, rms )
* unweighted mean and rms of n data points x(i)
* Input:
*	n	i4 number of data
*	x(n)	r4 data
* output:
*	bar	r4 mean value
*	rms	r4 root-mean-square
* 1999 Sep Keith Horne @ St.Andrews
* 2011 Feb KDH @ StA - x(i)**2 -> x(i)*x(i)
	real*4 x(*)
	real*8 sum, sum2
	bar = x(1)
	rms = 0.
	if( n.le.1 ) return
	sum = 0.d0
	sum2 = 0.d0
      do i=1,n
	sum = sum + x(i)
	sum2 = sum2 + x(i) * x(i)
      end do
	bar = sum / n
	var = ( sum2 - sum * sum / n ) / max( 1, n-1 )
	rms = sqrt( max( 0., var ) )
	return
	end
******************************************
	subroutine medmad_old( n, x, xmed, xmad )
* median and mean absolute deviation of n data points
* Input:
*	n	i4 number of data
*	x(n)	r4 data
* output:
*	xmed	r4 median value
*	xmad	r4 mean absolute deviation
* 2009 May Keith Horne @ St.Andrews
	real*4 x(*)
	real*8 sum
	xmed = xmedian( n, x )
	sum = 0.d0
      do i=1,n
	sum = sum + abs( x(i) - xmed )
      end do
	xmad = sum / max( 1, n )

* median absolute deviation = 1.4862 sigma
*   mean absolute deviation = 1.2533 sigma

	return
	end

******************************************
	subroutine medmad( n, d, dmed, dmad )
* median and MAD = median absolute deviation of d( n )
* or of a uniform sample from d( n ) if n > maxn
* In:	n	i4 number of data
*	d(n)	r4 data
* Out:	dmed	r4 median value
*	dmad	r4 median absolute deviation
* 2018 Feb Keith Horne @ St.Andrews
* 2019 Sep KDH @ St.A  maxn 1e5->1e6
	real*4 d( n )

	parameter( maxn = 1e6 )
	integer*4 ir( maxn )
	real*4 work( maxn )

	dmed = 0.
	dmad = 0.
	if( n .le. 0 ) return
	dmed = d( 1 )
	if( n .eq. 1 ) return

* short dataset
      if( n .le. maxn ) then 
	call medmad0( n, d, dmed, dmad, work, ir )
* sample long dataset
      else
	m = min( n, maxn )
      do i = 1, m
	p = ( i - 0.5 ) / m
	k = nint( 0.5 + p * n )
	work( i ) = d( k )
      end do
* medmad
	call medmad0( m, work, dmed, dmad, work, ir )
      end if

* median absolute deviation = 1.4862 sigma
*   mean absolute deviation = 1.2533 sigma

	return
	end

******************************************
	subroutine medmad0( n, d, dmed, dmad, a, ir )
* median and MAD = median absolute deviation of d( n )
* In:	n	i4 number of data
*	d( n )	r4 data
* Out:	dmed	r4 median value
*	dmad	r4 median absolute deviation
*	a( n )	r4 sorted absolute deviations
*	ir( n )	r4 dat(ir(i)) ascends
* 2018 Feb Keith Horne @ St.Andrews
	real*4 d( n ), a( n )
	integer*4 ir( n )
* median
	call median0( n, d, dmed, ir )
* mad
      do i = 1, n
	a( i ) = abs( d( i ) - dmed )
      end do
	call median0( n, a, dmad, ir )
	return
	end

******************************************
	subroutine median0( n, d, dmed, ir )
* median of d( n )
* In:	n	i4 number of data
*	d( n )	r4 data
* Out:	dmed	r4 median value
*	ir( n )	r4 d( ir( i ) ) ascends
* 2018 Feb Keith Horne @ St.Andrews
	real*4 d( n )
	integer*4 ir( n )

	dmed = 0.
	if( n .le. 0 ) return
	ir( 1 ) = 1
	if( n .eq. 1 ) return
* median
	call quicksort( n, d, ir )
	k = n / 2
	kp = k + 1
	dmed = d( ir( k ) )
	if( mod( n, 2 ) .eq. 0 ) dmed = ( d( ir( kp ) ) + dmed ) / 2.

	return
	end
***************************************************
	FUNCTION XMEDIAN( NDATA, DATA )
* find median of an array of real numbers
*  Input:
*	NDATA	= I4 number of data values
*	DATA	= R4 array of data values
*  Output:
*	MEDVAL	= R4 median value of data
*
* Mar 1988 Keith Horne @ STScI - add special cases N=1,2, trap invalid inputs
* Jan 1995 KDH @ St.And - fix bug when NPTS = 2, eliminate type.
* Mar 2000 KDH @ St.And - MAXDATA 10000->20000
* Mar 2000 KDH @ St.And - medval(x,n) -> xmedian(n,x)
* Sep 2000 KDH @ St.And - shellsort -> quicksort
* Sep 2003 KDH @ St.And - MAXDATA 20,000->100,000
* Jun 2012 KDH @ St.And - MAXDATA 100,000->200,000

	REAL*4 DATA(*)
	PARAMETER ( MAXDATA = 200000 )
	INTEGER*4 IRANK( MAXDATA )

* trap invalid dimensions

      IF( NDATA.LT.1 ) THEN
	WRITE(*,*) '** INVALID DIMENSION', NDATA 
	GOTO 999
      ELSE IF( NDATA.GT.MAXDATA ) THEN
	WRITE(*,*) '** WARNING: buffer overflow in XMEDIAN', NDATA, MAXDATA
      END IF
	NPTS = MIN( NDATA, MAXDATA )

      IF( NPTS.EQ.1 ) THEN
	XMEDIAN = DATA(1)
	GOTO 1000
      ELSE IF( NPTS.EQ.2 ) THEN
	XMEDIAN = 0.5 * ( DATA(1) + DATA(2) )
	GOTO 1000
      END IF

* compute rank of each data value using a sort

	CALL QUICKSORT( NPTS, DATA, IRANK )
	MEDRANK = ( NPTS + 1 ) / 2 
	XMEDIAN = DATA( IRANK( MEDRANK ) )

* average two values if the number of points is even 

      IF( MEDRANK + MEDRANK .EQ.NPTS ) THEN
	XMEDIAN = 0.5 * ( XMEDIAN + DATA( IRANK( MEDRANK+1 ) ) )
      END IF

* normal return

 1000	RETURN

* error return

  999	XMEDIAN = 0.
	WRITE(*,*) '** XMEDIAN aborted.'
	RETURN
	END
********************************************
	SUBROUTINE QUICKSORT( N, X, KEY )
*
* QUICKSORT subroutine (Hoare)
* Produces index such that X(KEY(I)) ascends with I
*
* Input:
*	N	= number of items to be sorted
*	X(N)	= REAL*4 values to be sorted
* Output:
*	KEY(N)	= sort key such that X(KEY(I)) ascends with I
*
* Adapted by T.R.Marsh from a Mike Irwin routine.
* X is left unchanged.  Compatible arguments with
* SHELLSORT and BUBSORT.  About 3 times faster than
* SHELLSORT for more than 5000 points, roughly the
* same below 2000.
*
* July 1985 KDH @ STScI -- typed in and tested.
* 1999 Oct KDH @ St.Andrews - eliminate TYPE
* 2001 Jul KDH @ St.And -- use .eqv. when comparing logicals
* 2018 Feb KDH @ St.And X(*),KEY(*)=>X(N),KEY(N)
	REAL*4 X( N )
	INTEGER*4 KEY( N )
	PARAMETER ( MAXSTACK=300 )
	INTEGER*4 ISTACK(2,MAXSTACK)
	LOGICAL*1 L1, L2, L3, CARRYON
      DO I = 1, N
	KEY(I) = I
      END DO
	CARRYON = .TRUE.
	M = 9
	IF( N.LE.M ) GOTO 500
	IL = 1
	IR = N
	ISP = 0			! Initialize stack pointer
      DO WHILE( CARRYON )
	I = IL
	J = IR + 1
	IM = ( IL + IR ) / 2	! Find median for partition
	L1 = X( KEY(IM) ) .LT. X( KEY(IR) )
	L2 = X( KEY(IR) ) .LT. X( KEY(IL) )
	L3 = X( KEY(IL) ) .LT. X( KEY(IM) )
      IF( L2 .EQV. L3 ) THEN
	IQ = IL
      ELSE IF( L1 .EQV. L3 ) THEN
	IQ = IM
      ELSE
	IQ = IR
      END IF
	XTEST = X( KEY(IQ) )
	ITEST = KEY(IQ)
	KEY(IQ) = KEY(IL)
	KEY(IL) = ITEST
  100	I = I + 1
      DO WHILE( X( KEY(I) ) .LT. XTEST )
	I = I + 1
      END DO
	J = J - 1
      DO WHILE( XTEST .LT. X( KEY(J) ) )
	J = J - 1
      END DO
      IF( J.LE.I ) THEN
	ITEST = KEY(J)
	KEY(J) = KEY(IL)
	KEY(IL) = ITEST
      ELSE
	ITEST = KEY(J)
	KEY(J) = KEY(I)
	KEY(I) = ITEST
	GOTO 100
      END IF

* Update stack

      IF( IR-J.GE.J-IL .AND. J-IL.GT.M ) THEN
	ISP = ISP + 1
	IF( ISP.GT.MAXSTACK ) write(*,*)
     &		CHAR(7),'*** STACK OVERFLOW IN QUICKSORT.'
	ISTACK(1,ISP) = J + 1
	ISTACK(2,ISP) = IR
	IR = J - 1
      ELSE IF( J-IL.GT.IR-J .AND. IR-J.GT.M ) THEN
	ISP = ISP + 1
	IF( ISP.GT.MAXSTACK ) write(*,*)
     &		CHAR(7),'*** STACK OVERFLOW IN QUICKSORT.'
	ISTACK(1,ISP) = IL
	ISTACK(2,ISP) = J - 1
	IL = J+1
      ELSE IF( IR-J.GT.M .AND. M.GE.J-IL ) THEN
	IL = J + 1
      ELSE IF( J-IL.GT.M .AND. M.GE.IR-J ) THEN
	IR = J-1
      ELSE IF( ISP.GT.0 ) THEN
	IL = ISTACK(1,ISP)
	IR = ISTACK(2,ISP)
	ISP = ISP - 1
      ELSE
	CARRYON = .FALSE.
      END IF
      END DO

* Straight insertion sort

  500	IF( M .EQ. 1 ) RETURN
	INT = 1
	IFIN = N - INT
      DO II = 1, IFIN
	I = II
	J = I + INT
      IF( X( KEY(I) ) .GT. X( KEY(J) ) ) THEN
	XTEST = X( KEY(J) )
	ITEST = KEY(J) 
	KEY(J) = KEY(I)
	J = I
	I = I - INT
      DO WHILE( I .GT. 0 .AND. X( KEY(MAX(1,I)) ) .GT. XTEST )
	KEY(J) = KEY(I)
	J = I
	I = I - INT
      END DO
	KEY(J) = ITEST
      END IF
      END DO
	RETURN
	END
C*LOOKUP -- inverse linear interpolation
C+
	SUBROUTINE LOOKUP( NPIX, DATA, TARGET, ILO, IHI, PART )
*
* Inverse linear interpolator. Finds ILO, IHI, PART such that
*   TARGET = DATA(ILO) * (1.-PART) + DATA(IHI) * PART
*
* This version uses safe but slow method of one-pixel steps
*
* Input:
*	NPIX	= Number of data points
*	DATA	= Data values ( must be monotonic )
*	TARGET	= Target value
* Output:
*	ILO	= Low pixel
*	IHI	= Hi pixel
*	PART	= fraction (0-1) for interpolation
C--
* Nov 1985 KDH @ STScI
* 2001 Aug KDH @ St.And - single point, IBEG
	REAL*4 DATA(*)
	DATA NLAST,IHILAST,ILOLAST/1,1,1/

* Null array
      IF( NPIX.LE.0 ) THEN
	ILO = 0
	IHI = 0
	PART = 0.
	RETURN
      END IF

* Single point
      IF( NPIX .EQ. 1 ) THEN
	ILO = 1
	IHI = 1
	PART = 0.
	RETURN
      END IF

* Use previous location
      IF( NPIX .EQ. NLAST ) THEN
      IF( DATA(IHILAST) .GE. TARGET .AND.
     *	  DATA(ILOLAST) .LE. TARGET ) THEN
	IHI = IHILAST
	ILO = ILOLAST
	GOTO 50
      END IF
	IBEG = ILOLAST
      ELSE
	IBEG = ( 1 + NPIX ) / 2
      END IF

* Determine uphill direction

      IF( DATA(1) .LE. DATA(NPIX) ) THEN
	ILOEND = 1
	IHIEND = NPIX
	IUP = 1
      ELSE
	ILOEND = NPIX
	IHIEND = 1
	IUP = -1
      END IF

* Locate an uphill point

      DO I = IBEG, IHIEND, IUP
	IF( DATA(I) .GE. TARGET ) GOTO 10
      END DO
	I = IHIEND
   10	IHI = I

* Locate nearest downhill point

      DO I = IHI, ILOEND, -IUP
	IF( DATA(I) .LE. TARGET ) GOTO 20
      END DO
	I = ILOEND
   20	ILO = I

* Locate nearest uphill point

      DO I = ILO, IHI, IUP
	IF( DATA(I) .GE. TARGET ) GOTO 30
      END DO
	I = IHI
   30	IHI = I

* Compute fractional part

   50	PART = DATA(IHI) - DATA(ILO)
	IF( PART.NE.0. ) PART = ( TARGET - DATA(ILO) ) / PART

	NLAST = NPIX
	ILOLAST = ILO
	IHILAST = IHI

	RETURN
	END
	function terp1( n, x, i1, i2, p )
* linear interpolation
* 2001 Aug Keith Horne @ St-Andrews
	real*4 x(n)
	terp1 = x(i1) * (1.-p) + x(i2) * p
	return
	end
***********************************
	subroutine peqpmc( n, x, y, c )
* x(i) = y(i) for i=1...n
	real*4 x(n), y(n)
	if( c .eq. 1.0 ) return
      do i=1,n
	x(i) = y(i) * c
      end do
	return
	end
C*BNU ... Planck function
C+
	FUNCTION BNU( WAVE, TEMP )
*
* input:
*	WAVE	R4 wavelength in Angstroms
*	TEMP	R4 Kelvin temperature
* output:
*	BNU	R4 blackbody intensity (erg/cm2/s/Hz/ster)
C--
* 1989 Jun Keith Horne @ STScI - clean up old subroutine
* 2002 Aug KDH @ St.And - better long-wavelength limit
	DATA C1/1.43883E8/
	DATA C2/1.95722E5/
	BNU = 0.
	X = WAVE * TEMP
	IF( X.LE.0. ) RETURN
	X = C1 / X
      IF( X.LT.1.E-4 ) THEN
	FACTOR = 2. / ( X * ( X + 2. ) )
      ELSE IF( X.LT.85. ) THEN
	FACTOR = 1. / ( EXP( X ) - 1. )
      ELSE
	bnuln = 3. * alog( ( c1 / wave ) / c2 ) - x
	bnu = exp( bnuln )
	RETURN
      END IF
	X = X * TEMP / C2
	BNU = FACTOR * X**3
	RETURN
	END 
*----------------------------------------
	function dl_zoo( z, om, ol )
* luminosity diameter distance in units of c/H_0
* 2004 Mar Keith Horne @ St-And
	dl_zoo = dp_zoo( z, om, ol ) * ( 1. + z )
	return
	end
*----------------------------------------
	function da_zoo( z, om, ol )
* angular diameter distance in units of c/H_0
* 2004 Mar Keith Horne @ St-And
	da_zoo = dp_zoo( z, om, ol ) / ( 1. + z )
	return
	end
*----------------------------------------
	function dp_zoo( z, om, ol )
* proper distance (with curvature) in units of c/H_0
* 2004 Mar Keith Horne @ St-And
* 2005 May KDH @ St-And -- 
	d = d_zoo( z, om, ol )
	o = om + ol
      if( o .gt. 1. ) then
	r0 = 1. / sqrt( o - 1. )
	dp = r0 * sin( d / r0 )
      else if( o .lt. 1. ) then
	r0 = 1. / sqrt( 1. - o )
	dp = r0 * sinh( d / r0 )
      else
	dp = d
      end if
	dp_zoo = dp
	return
	end
*----------------------------------------
	function d_zoo( z, om, ol )
* co-moving coordinate distance in units of c/H_0
* 2005 Dec KDH @ St-And
	d = cosmic( 0., z, 0., -1., 0., om, ol )
	d_zoo = d
	return
	end
*----------------------------------------
	function t_oo( om, ol )
* age in units of 1/H_0
* 2004 Feb Keith Horne @ St-And
	zinf = 1e3
	t_oo = t_zzooo( 0., zinf, 0., om, ol )
	return
	end
*----------------------------------------
	function t_zoo( z, om, ol )
* look-back time in units of 1/H_0
* 2004 Feb Keith Horne @ St-And
	t_zoo = t_zzooo( 0., z, 0., om, ol )
	return
	end
*----------------------------------------
	function t_zzooo( z1, z2, or, om, ol )
* elapsed time in units of 1/H_0
* 2004 Jul Keith Horne @ St-And
	t_zzooo = cosmic( z1, z2, -1., -1., or, om, ol )
	return
	end
*----------------------------------------
	function h_zooo( z, or, om, ol )
* Hubble parameter in units of H_0
* 2004 Jul Keith Horne @ St-And
	x = 1. + z
	x2 = x * x
	ok = 1. - (or + om + ol)
	hh = ol + x2 * ( ok + x * (om + x * or ) )
	h_zooo = sqrt( hh )
	return
	end
*----------------------------------------
	function cosmic( zlo, zhi, powx, powh, or, om, ol )
* integrate (1+z) ** powx * H(z) ** powh * dz
* 2004 Feb Keith Horne @ St-And
* 2004 Jul KDH @ St-And - implicit real*8
	implicit real*8 (a-h,o-z)
	real*4 cosmic, zlo, zhi, powx, powh, or, om, ol
	logical logint/.false./

	xlo = 1. + zlo
	xhi = 1. + zhi

      if( logint ) then
	xloglo = dlog( xlo )
	xloghi = dlog( xhi )
	powxp = powx + 1.
      end if

	powh2 = powh / 2.
	ok = 1.d0 - ( or + om + ol )

	n = abs( zhi - zlo ) * 1000
	n = max( 10, min( 10000, n ) )

	sum = 0.d0
      do i = 1, n
	part = ( i - 0.5 ) / n
      if( logint ) then
	xlog = xloglo * ( 1. - part ) + xloghi * part
	x = exp( xlog )
      else
	x = xlo * ( 1. - part ) + xhi * part
      end if

	h2 = ol + x*x * ( ok + x * ( om + x * or)  )

	add = 1.d0
      if( logint ) then
	if( powxp .ne. 0. ) add = x ** powxp
	if( powh2 .ne. 0. ) add = add * h2 ** powh2
      else
	if( powx .ne. 0. ) add = x ** powx
	if( powh2 .ne. 0. ) add = add * h2 ** powh2
      end if
	sum = sum + add
      end do

      if( logint ) then
	dlogx = ( xloghi - xloglo ) / n
	cosmic = dlogx * sum
      else
	dx = ( xhi - xlo ) / n
	cosmic = dx * sum
      end if
	return
	end
C*LOAD1D .. load up to 140 columns of data from multi-column ascii file
C+
      SUBROUTINE LOAD1D( NAME, NCOL, NROW,
     #	C1, C2, C3, C4, C5, C6, C7, C8, C9,
     #	C10, C11, C12, C13, C14, C15, C16, C17, C18, C19,
     #	C20, C21, C22, C23, C24, C25, C26, C27, C28, C29,
     #	C30, C31, C32, C33, C34, C35, C36, C37, C38, C39,
     #	C40, C41, C42, C43, C44, C45, C46, C47, C48, C49,
     #	C50, C51, C52, C53, C54, C55, C56, C57, C58, C59,
     #	C60, C61, C62, C63, C64, C65, C66, C67, C68, C69,
     #	C70, C71, C72, C73, C74, C75, C76, C77, C78, C79,
     #	C80, C81, C82, C83, C84, C85, C86, C87, C88, C89,
     #	C90, C91, C92, C93, C94, C95, C96, C97, C98, C99,
     #	C100, C101, C102, C103, C104, C105, C106, C107, C108, C109,
     #	C110, C111, C112, C113, C114, C115, C116, C117, C118, C119,
     #	C120, C121, C122, C123, C124, C125, C126, C127, C128, C129,
     #	C130, C131, C132, C133, C134, C135, C136, C137, C138, C139,
     #	C140 )
*
* Loads up to 140 columns of data from a free-format ASCII data file
* Use as follows:
*       parameter ( maxdat = 10000 )
*       real*4 c1(maxdat), c2(maxdat), c3(maxdat)
*       character*256 file
*       file = 'file.dat'
*       ndat = maxdat
*       call load1d( file, 3, ndat, c1, c2, c3 )
*
*  Input:
*	NAME	= C* Name of disk file, ' ' for interactive selection.
*	NCOL	= I4 Number of columns to be read (if NCOL > 0)
*		     (if NCOL < 0, then only column -NCOL is read)
*	NROW	= I4 Maximum allowed number of rows
*  Output:
*	NAME	= C* Name of disk file (only if selected interactively)
*	NROW	= I4 Number of data points (0 if none loaded)
*	C1(NROW)= R4 first column
*	C2(NROW)= R4 second column
*	  ...
*	C140(NROW)= R4 140th column
C--
* May 1987 Keith Horne @ STScI
* Sep 1987 KDH - add NCOL<0 option to read specified column.
* May 1988 KDH - add final read to pick up eof after reading all rows.
* Jul 1988 KDH @ STScI - change name from LOAD_ASCII to LOAD_DAT
* Sep 1991 KDH @ STScI - skip words when COL < 0.
* Jan 1994 KDH @ Utrecht - split LOAD_DAT to LOAD1D,LOAD1D10,...
* Mar 1994 KDH @ Utrecht - entry points LOAD1DX, X=1...20
* Mar 1998 KDH @ St.Andrews - MAXCOL=20->140
* 1999 Oct KDH @ St.And - change OPEN for unix compatibility
* 2001 Jul KDH @ St.And -- no // of character*(*)
* 2001 Aug KDH @ St.And -- use OPENFILE
* 2010 Mar KDH @ St.And -- continue to avoid err=n branch to same line
* 2018 Jul KDH @ St.And -- entry points LOAD1DX, X=21...35
	REAL*4 C1(*), C2(*), C3(*), C4(*), C5(*)
	REAL*4 C6(*), C7(*), C8(*), C9(*), C10(*)
	REAL*4 C11(*), C12(*), C13(*), C14(*), C15(*)
	REAL*4 C16(*), C17(*), C18(*), C19(*), C20(*)
	REAL*4 C21(*), C22(*), C23(*), C24(*), C25(*)
	REAL*4 C26(*), C27(*), C28(*), C29(*), C30(*)
	REAL*4 C31(*), C32(*), C33(*), C34(*), C35(*)
	REAL*4 C36(*), C37(*), C38(*), C39(*), C40(*)
	REAL*4 C41(*), C42(*), C43(*), C44(*), C45(*)
	REAL*4 C46(*), C47(*), C48(*), C49(*)
	REAL*4 C50(*), C51(*), C52(*), C53(*), C54(*)
	REAL*4 C55(*), C56(*), C57(*), C58(*), C59(*)
	REAL*4 C60(*), C61(*), C62(*), C63(*), C64(*)
	REAL*4 C65(*), C66(*), C67(*), C68(*), C69(*)
	REAL*4 C70(*), C71(*), C72(*), C73(*), C74(*)
	REAL*4 C75(*), C76(*), C77(*), C78(*), C79(*)
	REAL*4 C80(*), C81(*), C82(*), C83(*), C84(*)
	REAL*4 C85(*), C86(*), C87(*), C88(*), C89(*)
	REAL*4 C90(*), C91(*), C92(*), C93(*), C94(*)
	REAL*4 C95(*), C96(*), C97(*), C98(*), C99(*)
	REAL*4 C100(*), C101(*), C102(*), C103(*), C104(*)
	REAL*4 C105(*), C106(*), C107(*), C108(*), C109(*)
	REAL*4 C110(*), C111(*), C112(*), C113(*), C114(*)
	REAL*4 C115(*), C116(*), C117(*), C118(*), C119(*)
	REAL*4 C120(*), C121(*), C122(*), C123(*), C124(*)
	REAL*4 C125(*), C126(*), C127(*), C128(*), C129(*)
	REAL*4 C130(*), C131(*), C132(*), C133(*), C134(*)
	REAL*4 C135(*), C136(*), C137(*), C138(*), C139(*)
	REAL*4 C140(*)
	CHARACTER*(*) NAME
	CHARACTER*256 FILE	! file name
	PARAMETER (INUNIT=49)	! input unit number
	PARAMETER (MAXCOL=140)	! maximum number of columns that can be read
	CHARACTER*(15*MAXCOL) RECORD ! buffer into which data are read
	REAL*4 ROW( MAXCOL )

* entry points since UNIX requires number of arguments to match exactly


	ENTRY LOAD1D35( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29, C30, C31, C32, C33, C34, C35 )
        GOTO 100
	ENTRY LOAD1D34( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29, C30, C31, C32, C33, C34 )
        GOTO 100
	ENTRY LOAD1D33( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29, C30, C31, C32, C33 )
        GOTO 100
	ENTRY LOAD1D32( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29, C30, C31, C32 )
        GOTO 100
	ENTRY LOAD1D31( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29, C30, C31 )
        GOTO 100

	ENTRY LOAD1D30( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29, C30 )
        GOTO 100
	ENTRY LOAD1D29( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28, C29 )
        GOTO 100
	ENTRY LOAD1D28( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27, C28 )
        GOTO 100
	ENTRY LOAD1D27( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26, C27 )
        GOTO 100
	ENTRY LOAD1D26( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25,
     #	C26 )
        GOTO 100

	ENTRY LOAD1D25( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24, C25 )
        GOTO 100
	ENTRY LOAD1D24( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23, C24 )
        GOTO 100
	ENTRY LOAD1D23( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22, C23 )
        GOTO 100
	ENTRY LOAD1D22( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21, C22 )
        GOTO 100
	ENTRY LOAD1D21( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20, C21 )
        GOTO 100

	ENTRY LOAD1D20( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19, C20 )
        GOTO 100
	ENTRY LOAD1D19( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18, C19 )
        GOTO 100
	ENTRY LOAD1D18( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17, C18 )
        GOTO 100
	ENTRY LOAD1D17( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16, C17 )
        GOTO 100
	ENTRY LOAD1D16( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15,
     #	C16 )
        GOTO 100
        ENTRY LOAD1D15( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14, C15 )
        GOTO 100
        ENTRY LOAD1D14( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13, C14 )
        GOTO 100
        ENTRY LOAD1D13( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12, C13 )
        GOTO 100
        ENTRY LOAD1D12( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11, C12 )
        GOTO 100
        ENTRY LOAD1D11( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10, C11 )
        GOTO 100
        ENTRY LOAD1D10( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9, C10 )
        GOTO 100
        ENTRY LOAD1D9( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8, C9 )
        GOTO 100
        ENTRY LOAD1D8( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7, C8 )
        GOTO 100
        ENTRY LOAD1D7( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6, C7 )
        GOTO 100
        ENTRY LOAD1D6( NAME, NCOL, NROW, C1, C2, C3, C4, C5,
     #	C6 )
        GOTO 100
        ENTRY LOAD1D5( NAME, NCOL, NROW, C1, C2, C3, C4, C5 )
        GOTO 100
        ENTRY LOAD1D4( NAME, NCOL, NROW, C1, C2, C3, C4 )
        GOTO 100
        ENTRY LOAD1D3( NAME, NCOL, NROW, C1, C2, C3 )
        GOTO 100
        ENTRY LOAD1D2( NAME, NCOL, NROW, C1, C2 )
        GOTO 100
        ENTRY LOAD1D1( NAME, NCOL, NROW, C1 )
        GOTO 100

  100 IF( NROW.LE.0. OR. NCOL.EQ.0 ) THEN
	WRITE(*,*) '** INVALID DIMENSIONS.', NROW, NCOL
	GOTO 999
      ELSE IF( NCOL.GT.MAXCOL ) THEN
	WRITE(*,*) '** WARNING: Too many columns for LOAD1D.'
	WRITE(*,*) '** Columns', MAXCOL, ' thru', NCOL, ' skipped.'
      END IF

	MAXROW = NROW

* get name of input file

	FILE = NAME
c        WRITE(*,*) 'FILE"', NAME( :LEN1(NAME) ), '"'
      IF( FILE.EQ.' ' ) THEN
	WRITE(*,'(A,I3,A)')
     #	' Enter name of', IABS(NCOL),'-column free-format ASCII file :'
	READ(*,'(A)') FILE
	CALL WORD1( FILE, L1, L2 )
	IF( L1 .LE. 0 ) GOTO 99
	NAME = FILE( :LEN1( FILE ) )
      END IF

* get full path to file
c	CALL FULLPATH( NAME, FILE )
c        WRITE(*,*) 'PATH"', FILE( :LEN1(FILE) ), '"'

* open input disk file

C	OPEN( UNIT=INUNIT, NAME=FILE, FORM='FORMATTED', TYPE='OLD',
C    &	ERR=99, READONLY )
* this works on unix machines
C	OPEN( UNIT=INUNIT, NAME=FILE, FORM='FORMATTED', STATUS='OLD',
C     &	ERR=99 )
	CALL OPENFILE( INUNIT, FILE, 'FORMATTED', 'OLD', IER )
	IF( IER .NE. 0 ) GOTO 99

* read desired column by skipping unwanted columns that precede it.

      IF( NCOL.LE.-2 ) THEN
	NSKIP = IABS(NCOL)-1
      DO I=1,MAXROW
   51	continue
	READ( INUNIT, '(A)', ERR=51, END=100 ) RECORD
      DO J=1,NSKIP
	CALL WORD1( RECORD, L1, L2 )
	IF( L1.LT.1 ) GOTO 51
	RECORD = RECORD(L2+1:)
      END DO
	READ( RECORD, *, ERR=51, END=51 ) C1(I)
C	READ( RECORD, *, ERR=51, END=51 ) 
C     #	(DUMMY,J=1,NSKIP), C1(I)
      END DO

* read desired consecutive columns of data from the file

      ELSE IF( IABS(NCOL).EQ.1 ) THEN
      DO I=1,MAXROW
    1   continue
	READ( INUNIT, '(A)', ERR=1, END=200 ) RECORD
	READ( RECORD, *, ERR=1, END=1 ) 
     #	C1(I)
      END DO

      ELSE IF( NCOL.EQ.2 ) THEN
      DO I=1,MAXROW
    2	continue
	READ( INUNIT, '(A)', ERR=2, END=200 ) RECORD
	READ( RECORD, *, ERR=2, END=2 ) 
     #	C1(I), C2(I)
      END DO

      ELSE IF( NCOL.EQ.3 ) THEN
      DO I=1,MAXROW
    3	continue
	READ( INUNIT, '(A)', ERR=3, END=200 ) RECORD
	READ( RECORD, *, ERR=3, END=3 ) 
     #	C1(I), C2(I), C3(I)
      END DO

      ELSE IF( NCOL.EQ.4 ) THEN
      DO I=1,MAXROW
    4	continue
	READ( INUNIT, '(A)', ERR=4, END=200 ) RECORD
	READ( RECORD, *, ERR=4, END=4 ) 
     #	C1(I), C2(I), C3(I), C4(I)
      END DO

      ELSE IF( NCOL.EQ.5 ) THEN
      DO I=1,MAXROW
    5	continue
	READ( INUNIT, '(A)', ERR=5, END=200 ) RECORD
	READ( RECORD, *, ERR=5, END=5 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I)
      END DO

      ELSE IF( NCOL.EQ.6 ) THEN
      DO I=1,MAXROW
    6   continue
	READ( INUNIT, '(A)', ERR=6, END=200 ) RECORD
        READ( RECORD, *, ERR=6, END=6 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I)
      END DO


      ELSE IF( NCOL.EQ.7 ) THEN
      DO I=1,MAXROW
    7   continue
	READ( INUNIT, '(A)', ERR=7, END=200 ) RECORD
	READ( RECORD, *, ERR=7, END=7 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I)
      END DO


      ELSE IF( NCOL.EQ.8 ) THEN
      DO I=1,MAXROW
    8   continue
	READ( INUNIT, '(A)', ERR=8, END=200 ) RECORD
	READ( RECORD, *, ERR=8, END=8 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I)
      END DO


      ELSE IF( NCOL.EQ.9 ) THEN
      DO I=1,MAXROW
    9   continue
	READ( INUNIT, '(A)', ERR=9, END=200 ) RECORD
	READ( RECORD, *, ERR=9, END=9 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I)
      END DO


      ELSE IF( NCOL.EQ.10 ) THEN
      DO I=1,MAXROW
   10   continue
	READ( INUNIT, '(A)', ERR=10, END=200 ) RECORD
	READ( RECORD, *, ERR=10, END=10 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I)
      END DO

      ELSE IF( NCOL.EQ.11 ) THEN
      DO I=1,MAXROW
   11   continue
	READ( INUNIT, '(A)', ERR=11, END=200 ) RECORD
	READ( RECORD, *, ERR=11, END=11 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I)
      END DO


      ELSE IF( NCOL.EQ.12 ) THEN
      DO I=1,MAXROW
   12   continue
	READ( INUNIT, '(A)', ERR=12, END=200 ) RECORD
	READ( RECORD, *, ERR=12, END=12 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I)
      END DO


      ELSE IF( NCOL.EQ.13 ) THEN
      DO I=1,MAXROW
   13   continue
	READ( INUNIT, '(A)', ERR=13, END=200 ) RECORD
	READ( RECORD, *, ERR=13, END=13 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I)
      END DO


      ELSE IF( NCOL.EQ.14 ) THEN
      DO I=1,MAXROW
   14   continue
	READ( INUNIT, '(A)', ERR=14, END=200 ) RECORD
	READ( RECORD, *, ERR=14, END=14 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I)
      END DO


      ELSE IF( NCOL.EQ.15 ) THEN
      DO I=1,MAXROW
   15   continue
	READ( INUNIT, '(A)', ERR=15, END=200 ) RECORD
	READ( RECORD, *, ERR=15, END=15 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I), C15(I)
      END DO


      ELSE IF( NCOL.EQ.16 ) THEN
      DO I=1,MAXROW
   16   continue
	READ( INUNIT, '(A)', ERR=16, END=200 ) RECORD
	READ( RECORD, *, ERR=16, END=16 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I), C15(I),
     #	C16(I)
      END DO

      ELSE IF( NCOL.EQ.17 ) THEN
      DO I=1,MAXROW
   17   continue
	READ( INUNIT, '(A)', ERR=17, END=200 ) RECORD
	READ( RECORD, *, ERR=17, END=17 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I), C15(I),
     #	C16(I), C17(I)
      END DO

      ELSE IF( NCOL.EQ.18 ) THEN
      DO I=1,MAXROW
   18   continue
	READ( INUNIT, '(A)', ERR=18, END=200 ) RECORD
	READ( RECORD, *, ERR=18, END=18 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I), C15(I),
     #	C16(I), C17(I), C18(I)
      END DO

      ELSE IF( NCOL.EQ.19 ) THEN
      DO I=1,MAXROW
   19   continue
	READ( INUNIT, '(A)', ERR=19, END=200 ) RECORD
	READ( RECORD, *, ERR=19, END=19 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I), C15(I),
     #	C16(I), C17(I), C18(I), C19(I)
      END DO

      ELSE IF( NCOL.EQ.20 ) THEN
      DO I=1,MAXROW
   20   continue
	READ( INUNIT, '(A)', ERR=20, END=200 ) RECORD
	READ( RECORD, *, ERR=20, END=20 ) 
     #	C1(I), C2(I), C3(I), C4(I), C5(I),
     #	C6(I), C7(I), C8(I), C9(I), C10(I),
     #	C11(I), C12(I), C13(I), C14(I), C15(I),
     #	C16(I), C17(I), C18(I), C19(I), C20(I)
      END DO

      ELSE IF( NCOL.LE.MAXCOL ) THEN
      DO I=1,MAXROW
   50   continue
	READ( INUNIT, '(A)', ERR=50, END=200 ) RECORD
	READ( RECORD, *, ERR=50, END=50 ) 
     #		(ROW(J),J=1,NCOL)
	IF( NCOL.GE.1 ) C1(I) = ROW(1)
	IF( NCOL.GE.2 ) C2(I) = ROW(2)
	IF( NCOL.GE.3 ) C3(I) = ROW(3)
	IF( NCOL.GE.4 ) C4(I) = ROW(4)
	IF( NCOL.GE.5 ) C5(I) = ROW(5)
	IF( NCOL.GE.6 ) C6(I) = ROW(6)
	IF( NCOL.GE.7 ) C7(I) = ROW(7)
	IF( NCOL.GE.8 ) C8(I) = ROW(8)
	IF( NCOL.GE.9 ) C9(I) = ROW(9)
	IF( NCOL.GE.10 ) C10(I) = ROW(10)
	IF( NCOL.GE.11 ) C11(I) = ROW(11)
	IF( NCOL.GE.12 ) C12(I) = ROW(12)
	IF( NCOL.GE.13 ) C13(I) = ROW(13)
	IF( NCOL.GE.14 ) C14(I) = ROW(14)
	IF( NCOL.GE.15 ) C15(I) = ROW(15)
	IF( NCOL.GE.16 ) C16(I) = ROW(16)
	IF( NCOL.GE.17 ) C17(I) = ROW(17)
	IF( NCOL.GE.18 ) C18(I) = ROW(18)
	IF( NCOL.GE.19 ) C19(I) = ROW(19)
	IF( NCOL.GE.20 ) C20(I) = ROW(20)
	IF( NCOL.GE.21 ) C21(I) = ROW(21)
	IF( NCOL.GE.22 ) C22(I) = ROW(22)
	IF( NCOL.GE.23 ) C23(I) = ROW(23)
	IF( NCOL.GE.24 ) C24(I) = ROW(24)
	IF( NCOL.GE.25 ) C25(I) = ROW(25)
	IF( NCOL.GE.26 ) C26(I) = ROW(26)
	IF( NCOL.GE.27 ) C27(I) = ROW(27)
	IF( NCOL.GE.28 ) C28(I) = ROW(28)
	IF( NCOL.GE.29 ) C29(I) = ROW(29)
	IF( NCOL.GE.30 ) C30(I) = ROW(30)
	IF( NCOL.GE.31 ) C31(I) = ROW(31)
	IF( NCOL.GE.32 ) C32(I) = ROW(32)
	IF( NCOL.GE.33 ) C33(I) = ROW(33)
	IF( NCOL.GE.34 ) C34(I) = ROW(34)
	IF( NCOL.GE.35 ) C35(I) = ROW(35)
	IF( NCOL.GE.36 ) C36(I) = ROW(36)
	IF( NCOL.GE.37 ) C37(I) = ROW(37)
	IF( NCOL.GE.38 ) C38(I) = ROW(38)
	IF( NCOL.GE.39 ) C39(I) = ROW(39)
	IF( NCOL.GE.40 ) C40(I) = ROW(40)
	IF( NCOL.GE.41 ) C41(I) = ROW(41)
	IF( NCOL.GE.42 ) C42(I) = ROW(42)
	IF( NCOL.GE.43 ) C43(I) = ROW(43)
	IF( NCOL.GE.44 ) C44(I) = ROW(44)
	IF( NCOL.GE.45 ) C45(I) = ROW(45)
	IF( NCOL.GE.46 ) C46(I) = ROW(46)
	IF( NCOL.GE.47 ) C47(I) = ROW(47)
	IF( NCOL.GE.48 ) C48(I) = ROW(48)
	IF( NCOL.GE.49 ) C49(I) = ROW(49)
	IF( NCOL.GE.50 ) C50(I) = ROW(50)
	IF( NCOL.GE.51 ) C51(I) = ROW(51)
	IF( NCOL.GE.52 ) C52(I) = ROW(52)
	IF( NCOL.GE.53 ) C53(I) = ROW(53)
	IF( NCOL.GE.54 ) C54(I) = ROW(54)
	IF( NCOL.GE.55 ) C55(I) = ROW(55)
	IF( NCOL.GE.56 ) C56(I) = ROW(56)
	IF( NCOL.GE.57 ) C57(I) = ROW(57)
	IF( NCOL.GE.58 ) C58(I) = ROW(58)
	IF( NCOL.GE.59 ) C59(I) = ROW(59)
	IF( NCOL.GE.60 ) C60(I) = ROW(60)
	IF( NCOL.GE.61 ) C61(I) = ROW(61)
	IF( NCOL.GE.62 ) C62(I) = ROW(62)
	IF( NCOL.GE.63 ) C63(I) = ROW(63)
	IF( NCOL.GE.64 ) C64(I) = ROW(64)
	IF( NCOL.GE.65 ) C65(I) = ROW(65)
	IF( NCOL.GE.66 ) C66(I) = ROW(66)
	IF( NCOL.GE.67 ) C67(I) = ROW(67)
	IF( NCOL.GE.68 ) C68(I) = ROW(68)
	IF( NCOL.GE.69 ) C69(I) = ROW(69)
	IF( NCOL.GE.70 ) C70(I) = ROW(70)
	IF( NCOL.GE.71 ) C71(I) = ROW(71)
	IF( NCOL.GE.72 ) C72(I) = ROW(72)
	IF( NCOL.GE.73 ) C73(I) = ROW(73)
	IF( NCOL.GE.74 ) C74(I) = ROW(74)
	IF( NCOL.GE.75 ) C75(I) = ROW(75)
	IF( NCOL.GE.76 ) C76(I) = ROW(76)
	IF( NCOL.GE.77 ) C77(I) = ROW(77)
	IF( NCOL.GE.78 ) C78(I) = ROW(78)
	IF( NCOL.GE.79 ) C79(I) = ROW(79)
	IF( NCOL.GE.80 ) C80(I) = ROW(80)
	IF( NCOL.GE.81 ) C81(I) = ROW(81)
	IF( NCOL.GE.82 ) C82(I) = ROW(82)
	IF( NCOL.GE.83 ) C83(I) = ROW(83)
	IF( NCOL.GE.84 ) C84(I) = ROW(84)
	IF( NCOL.GE.85 ) C85(I) = ROW(85)
	IF( NCOL.GE.86 ) C86(I) = ROW(86)
	IF( NCOL.GE.87 ) C87(I) = ROW(87)
	IF( NCOL.GE.88 ) C88(I) = ROW(88)
	IF( NCOL.GE.89 ) C89(I) = ROW(89)
	IF( NCOL.GE.90 ) C90(I) = ROW(90)
	IF( NCOL.GE.91 ) C91(I) = ROW(91)
	IF( NCOL.GE.92 ) C92(I) = ROW(92)
	IF( NCOL.GE.93 ) C93(I) = ROW(93)
	IF( NCOL.GE.94 ) C94(I) = ROW(94)
	IF( NCOL.GE.95 ) C95(I) = ROW(95)
	IF( NCOL.GE.96 ) C96(I) = ROW(96)
	IF( NCOL.GE.97 ) C97(I) = ROW(97)
	IF( NCOL.GE.98 ) C98(I) = ROW(98)
	IF( NCOL.GE.99 ) C99(I) = ROW(99)
	IF( NCOL.GE.100 ) C100(I) = ROW(100)
	IF( NCOL.GE.101 ) C101(I) = ROW(101)
	IF( NCOL.GE.102 ) C102(I) = ROW(102)
	IF( NCOL.GE.103 ) C103(I) = ROW(103)
	IF( NCOL.GE.104 ) C104(I) = ROW(104)
	IF( NCOL.GE.105 ) C105(I) = ROW(105)
	IF( NCOL.GE.106 ) C106(I) = ROW(106)
	IF( NCOL.GE.107 ) C107(I) = ROW(107)
	IF( NCOL.GE.108 ) C108(I) = ROW(108)
	IF( NCOL.GE.109 ) C109(I) = ROW(109)
	IF( NCOL.GE.110 ) C110(I) = ROW(110)
	IF( NCOL.GE.111 ) C111(I) = ROW(111)
	IF( NCOL.GE.112 ) C112(I) = ROW(112)
	IF( NCOL.GE.113 ) C113(I) = ROW(113)
	IF( NCOL.GE.114 ) C114(I) = ROW(114)
	IF( NCOL.GE.115 ) C115(I) = ROW(115)
	IF( NCOL.GE.116 ) C116(I) = ROW(116)
	IF( NCOL.GE.117 ) C117(I) = ROW(117)
	IF( NCOL.GE.118 ) C118(I) = ROW(118)
	IF( NCOL.GE.119 ) C119(I) = ROW(119)
	IF( NCOL.GE.120 ) C120(I) = ROW(120)
	IF( NCOL.GE.121 ) C121(I) = ROW(121)
	IF( NCOL.GE.122 ) C122(I) = ROW(122)
	IF( NCOL.GE.123 ) C123(I) = ROW(123)
	IF( NCOL.GE.124 ) C124(I) = ROW(124)
	IF( NCOL.GE.125 ) C125(I) = ROW(125)
	IF( NCOL.GE.126 ) C126(I) = ROW(126)
	IF( NCOL.GE.127 ) C127(I) = ROW(127)
	IF( NCOL.GE.128 ) C128(I) = ROW(128)
	IF( NCOL.GE.129 ) C129(I) = ROW(129)
	IF( NCOL.GE.130 ) C130(I) = ROW(130)
	IF( NCOL.GE.131 ) C131(I) = ROW(131)
	IF( NCOL.GE.132 ) C132(I) = ROW(132)
	IF( NCOL.GE.133 ) C133(I) = ROW(133)
	IF( NCOL.GE.134 ) C134(I) = ROW(134)
	IF( NCOL.GE.135 ) C135(I) = ROW(135)
	IF( NCOL.GE.136 ) C136(I) = ROW(136)
	IF( NCOL.GE.137 ) C137(I) = ROW(137)
	IF( NCOL.GE.138 ) C138(I) = ROW(138)
	IF( NCOL.GE.139 ) C139(I) = ROW(139)
	IF( NCOL.GE.140 ) C140(I) = ROW(140)
      END DO

      END IF
	I = MAXROW + 1

* final read to pick up the eof

	READ( INUNIT, '(A)', ERR=200, END=200 ) RECORD

* give warning if there was data past the final row

	WRITE(*,*)  '** WARNING: Data beyond row', MAXROW, ' was not read.'
  200	NROW = I - 1

* close input file

	CLOSE( UNIT = INUNIT )

* normal return

 1000	RETURN

* Error return

   99	WRITE(*,*) CHAR(7),'** Error opening file. '
	WRITE(*,*) '** ', FILE(:LEN1(FILE))
	GOTO 999

  999   WRITE(*,*) '** LOAD1D aborted.'
	NROW = 0
	RETURN
	END
C* DUMP1D ... Dump NCOL arrays to a free-format ascii disk file.
C+
      SUBROUTINE DUMP1D( FILE, NCOL, NROW
     &	, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &	, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20 )
* Input:
*	FILE	C* file name or ' ' for interactive selection
*	NCOL	I4 Number of data arrays
*	NROW	I4 Number of data values in each array
*	X1-XNCOL	R4(NROW) The data arrays
C--
* Oct 1985 by KDH at STScI
* Apr 1988 KDH @ STScI - generalized from DUMPMONGO
* May 1991 KDH @ STScI - install in XCAL
* Jan 1994 KDH @ Utrecht - maxcol=5, DUMP_DAT -> DUMP1D, improve error msgs
* Mar 1994 KDH @ Utrecht - entry points for 1-5 images
* May 1995 KDH @ St.And - fix bug causing no output when file specified
* 1999 Nov KDH @ St.And - maxcol=5->6
* 2000 Nov KDH @ Austin - open with type='unknown'      
* 2001 Aug KDH @ St.And - use OPENFILE
* 2002 Jul KDH @ St.And - add array 7
* 2007 Oct KDH @ Christchurch - arrays 8-20
* 2007 Oct KDH @ Christchurch - PARAMETER( IOUT = 93 )
* 2007 Oct KDH @ Christchurch -- WRITE(RECORD,*), then WRITE(IOUT,'(A)')
*	       since WRITE(IOUT,*) sometimes writes more than 1 line
* 2012 Oct KDH @ StA - extend 10 to 20 columns
	REAL*4 X1(*), X2(*), X3(*), X4(*), X5(*)
	REAL*4 X6(*), X7(*), X8(*), X9(*), X10(*)
	REAL*4 X11(*), X12(*), X13(*), X14(*), X15(*)
	REAL*4 X16(*), X17(*), X18(*), X19(*), X20(*)
	CHARACTER*(*) FILE
	CHARACTER*128 OFILE
	CHARACTER*512 RECORD
	CHARACTER*10 REPLY
	PARAMETER ( MAXCOL = 20 )
	PARAMETER ( IOUT = 93 )
	GOTO 100
* multiple entries because unix/linux can't handle more arrays in
* the subroutine argument list than in the calling program
	ENTRY DUMP1D20( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20 )
	GOTO 100
	ENTRY DUMP1D19( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14, X15, X16, X17, X18, X19 )
	GOTO 100
	ENTRY DUMP1D18( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14, X15, X16, X17, X18 )
	GOTO 100
	ENTRY DUMP1D17( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14, X15, X16, X17 )
	GOTO 100
	ENTRY DUMP1D16( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14, X15, X16 )
	GOTO 100
	ENTRY DUMP1D15( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14, X15 )
	GOTO 100
	ENTRY DUMP1D14( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13, X14 )
	GOTO 100
	ENTRY DUMP1D13( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12, X13 )
	GOTO 100
	ENTRY DUMP1D12( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11, X12 )
	GOTO 100
	ENTRY DUMP1D11( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10
     &		, X11 )
	GOTO 100
	ENTRY DUMP1D10( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10 )
	GOTO 100
	ENTRY DUMP1D9( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8, X9 )
	GOTO 100
	ENTRY DUMP1D8( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7, X8 )
	GOTO 100
	ENTRY DUMP1D7( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6, X7 )
	GOTO 100
	ENTRY DUMP1D6( FILE, NCOL, NROW
     &		, X1, X2, X3, X4, X5, X6 )
	GOTO 100
	ENTRY DUMP1D5( FILE, NCOL, NROW, X1, X2, X3, X4, X5 )
	GOTO 100
	ENTRY DUMP1D4( FILE, NCOL, NROW, X1, X2, X3, X4 )
	GOTO 100
	ENTRY DUMP1D3( FILE, NCOL, NROW, X1, X2, X3 )
	GOTO 100
	ENTRY DUMP1D2( FILE, NCOL, NROW, X1, X2 )
	GOTO 100
	ENTRY DUMP1D1( FILE, NCOL, NROW, X1 )

* check dimensions
  100 IF( NCOL .LT. 1 ) RETURN
      IF( NROW .LT. 1 ) RETURN
      NC = MIN( NCOL, MAXCOL )
      IF( NC .NE. NCOL ) THEN
        WRITE(*,*) '**DUMP1D: skipping columns', NC + 1,'-', NCOL
      END IF
      OFILE = FILE
* get filename
   10 IF( OFILE.EQ.' ' ) THEN
        WRITE(*,'(A,$)') ' Output filename : '
	READ(*,'(A)') OFILE
	REPLY = OFILE
	CALL UPCASE( REPLY, REPLY )
	IF( REPLY .EQ. ' ') RETURN
	IF( REPLY .EQ. 'NONE' ) RETURN
      END IF
* open file
	CALL OPENFILE( IOUT, OFILE, 'FORMATTED', 'UNKNOWN', IER )
	IF( IER .NE. 0 ) GOTO 99
      WRITE(*,*) 'File open ', OFILE(:LEN1(OFILE))
* write data
      NR = 0
      DO I=1,NROW
      IF( NC .EQ. 1 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I)
      ELSE IF( NC .EQ. 2 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I)
      ELSE IF( NC .EQ. 3 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I)
      ELSE IF( NC .EQ. 4 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I)
      ELSE IF( NC .EQ. 5 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
      ELSE IF( NC .EQ. 6 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I)
      ELSE IF( NC .EQ. 7 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I)
      ELSE IF( NC .EQ. 8 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I), X8(I)
      ELSE IF( NC .EQ. 9 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I), X8(I), X9(I)
      ELSE IF( NC .EQ. 10 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I), X8(I), X9(I), X10(I)
      ELSE IF( NC .EQ. 11 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I), X8(I), X9(I), X10(I), x11(i)
      ELSE IF( NC .EQ. 12 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I), X8(I), X9(I), X10(I), x11(i), x12(i)
      ELSE IF( NC .EQ. 13 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X6(I), X7(I), X8(I), X9(I), X10(I), x11(i), x12(i), x13(i)
      ELSE IF( NC .EQ. 14 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I)
      ELSE IF( NC .EQ. 15 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I), x15(i)
      ELSE IF( NC .EQ. 16 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I), x15(i), x16(i)
      ELSE IF( NC .EQ. 17 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I), x15(i), x16(i), x17(i)
      ELSE IF( NC .EQ. 18 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I), x15(i), x16(i), x17(i), x18(i)
      ELSE IF( NC .EQ. 19 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I), x15(i), x16(i), x17(i), x18(i), x19(i)
      ELSE IF( NC .EQ. 20 ) THEN
	WRITE( RECORD, *, ERR=98 ) X1(I), X2(I), X3(I), X4(I), X5(I)
     &	, X14(I), x15(i), x16(i), x17(i), x18(i), x19(i), x20(i)
      ELSE
	WRITE(*,*) '** DUMP1D TOO MANY COLUMNS', NC
	STOP ':(((((.'
      END IF
	WRITE( IOUT, '(A)', ERR=98 ) RECORD( : LEN1( RECORD ) )
        NR = I
      END DO

* close file
   50 WRITE(*,*) 'Dump columns:', NC, ' rows:', NR, ' of', NROW
      CLOSE( UNIT = IOUT )
* normal return
 1000 RETURN

* file write error
 98   WRITE(*,*) '**DUMP1D: error writing to file ', OFILE(:LEN1(OFILE))
      GOTO 50
       
* file open error
 99   WRITE(*,*) '**DUMP1D: error opening file ', OFILE(:LEN1(OFILE))
      OFILE = ' '
      GOTO 10
      END
	subroutine openfile( iunit, name, form, status, ier )
* Input:
*	iunit	i4 unit number
*	name	c* file name
*	form	c* 'formatted' or 'unformatted'
*	status	c* 'old', 'new', ...
* Output:
*	ier	i4 0 if successful
* 2001 Jul Keith Horne @ St.Andrews -- use FULLPATH
* 2001 Aug KDH @ St-And -- silence
* 2001 Aug KDH @ St-And -- default extension .dat

	character*(*)	name, form, status
	character*256	path
	logical silence
	common/silence/ silence

* get full path
	call fullpath( name, path )
      if( .not. silence ) then
      if( path .ne. name ) then
	write(*,*) 'Path( ', name(:len1(name)), ' ) = ', 
     &		path(:len1(path)) 
      end if
      end if

* open file
	open( unit=iunit, file=path, form=form, status=status, err=2 )
	goto 10

* try again with default extension
    2	path = path( :len1( path ) ) // '.dat'
	open( unit=iunit, file=path, form=form, status=status, err=3 )
	goto 10

* normal return

   10	ier = 0
      if( .not. silence ) then
	write(*,*) 'Open(', iunit, ', ', path(:len1(path)),
     &	', ', form(:len1(form)), ', ', status(:len1(status)), ' )' 
      end if
	return

* error return

    3	ier = 1
	write(*,*) '** OPENFILE FAILED :( '
	write(*,*) '** Path( ', name(:len1(name)), ' ) = ', 
     &		path(:len1(path)) 
	return
	end
	subroutine fullpath( s1, s2 )
* evaluate unix path by translating environment variables
* 2001 Jul Keith Horne @ St.Andrews
	character*(*) s1, s2
	parameter ( maxchar = 256 )
	character*( maxchar ) t1, t2
* null string
	call word1( s1, l1, l2 )
      if( l1.lt.1 ) then
	s2 = s1
	return
      end if
* copy to a local string
	t1 = s1( l1:l2 )
      if( l2-l1+1 .gt. maxchar ) then
	write(*,*) '** fullpath truncation:', l1, l2, l2-l1+1, maxchar
	write(*,*) '"', s1( :len1(s1) ), '"'
	write(*,*) '"', t1( :len1(t1) ), '"'
      end if
* isolate the environment variable to be translated
	ldollar = index( t1, '$' )
      do while( ldollar .gt. 0 )
	l1 = ldollar + 1
	lslash = ldollar + index( t1( l1:len1(t1) ) // '/', '/' )
	l2 = lslash - 1
* translate the environment variable
	call getenv( t1( l1:l2 ), t2 )
      if( len1( t2 ) .ge. maxchar ) then
	write(*,*) '** fullpath possible truncation:'
	write(*,*) '"', t1( :len1(t1) ), '"'
	write(*,*) '"', t2( :len1(t2) ), '"'
      end if

	if( ldollar .gt. 1 ) t2 = t1(:ldollar-1) // t2(:len1(t2) )

	t1 = t2(:len1(t2) ) // t1(l2+1:)

	ldollar = index( t1, '$' )
      end do
* copy to output string
	s2 = t1
	l1 = len1( t1 )
	l2 = len1( s2 )
      if( l1 .gt. l2 ) then
	write(*,*) '** fullpath truncation:', l1, l2
	write(*,*) '"', t1(:l1), '"'
	write(*,*) '"', s2(:l2), '"'
      end if
	return
	end
************************      
      subroutine cpu( tcpu )
*
* cpu time in seconds since previous call
*
* 1994 Mar kdh @ Utrecht - use dtime() to get the answer
* 2010 Mar kdh @ LCO - try cpu_time from gfortran
	real*4 told/-1./
c	real*4 tarray(2)
c	external dtime

	call cpu_time( tnow )
	tcpu = tnow - told
	if( told .lt. 0. ) tcpu = 0.
	told = tnow

c	tcpu = dtime( tarray )
c	tcpu = etime( tarray )
c	tuser = tarray(1)
c	tsystem = tarray(2)
	return
	end
************************
	subroutine spawn( task )
*
* spawn a task to the operating system
*
* 1994 Jan KDH @ Utrecht - dummy version
* 1994 Mar KDH @ Utrecht - call system( task )
* 1994 Mar KDH @ Utrecht - call system( spawn task )
* 2000 Nov KDH @ Austin - call system( task )
* 2001 Jul KDH @ St-And - linux g77 rejects char*(*) concatenations
* 2015 Aug KDH @ St-And - comment out external system
	character*(*) task

* external seems not to be needed and causes problems on some compilers
c	external system

* at Utrecht, this works:
c	write(*,*) ' call system( spawn ' // TASK(:LEN1(TASK)) // ' )'
c	write(*,*) ' '
c	CALL SYSTEM( 'spawn ' // TASK(:LEN1(TASK)) )
* Note: $HOME/bin/spawn is an executable script as follows :
*  #!/bin/csh
*  source $HOME/.login ; $*
* at Austin, this works:
c	write(*,*) ' call system( ' // task(:len1(task)) // ' )'
c	call system( task(:len1(task)) )

* on linux, this works
* Note: linux g77 rejects character*(*) concatenation
	l1 = len1( task )
	write(*,*) 'call system( ', task(:l1), ' ) '
	call system( task( :l1 ) )
	return
	end
C* PGSTART ... open a plot
C+
      SUBROUTINE PGSTART( NX, NY, NEW )
*
* PGSTART opens and PGSTOP closes a plot.
* Like PGBEGIN and PGEND, but keeping track of previous plot device
* and local plot queues so plots can be spooled at run time.
*
* Input:
*	NX,NY	I4 number of plots per page in X and Y direction
*		selected interactively if input values are 0
* Output:
*	NEW	I4 1 if opening a new plot, otherwise
*	NX,NY	I4 if changed by interaction
*
C--
* 1980 KDH at CIT - Original version
* 1985 Nov KDH @ STScI - added QMS devices
* 1986 Apr KDH @ STScI - added laser printer
* 1987 Jul KDH @ STScI - added Steve Baron's DeAnza driver
* 1987 Nov AWC @ STScI - option not to plot named Q,Z,P,V files.
* 1987 Dec KDH @ STScI - add retrographics option
* 1987 Dec KDH @ STScI - add metafile options
* 1988 Feb KDH @ STScI - call PGIDEN in PGSTOP
* 1988 Nov KDH @ STSCI - PRINT/NOFLAG
* 1989 Apr KDH @ STSCI - CHANGE QUEUE NAMES
* 1989 Jun KDH @ STSCI - put flag page back in
* 1989 Oct KDH @ STSCI - change queue 3C_T800 to 3C_T2400
* 1989 Nov KDH @ STSCI - add VAX WS option, also ?
* 1989 Dec KDH @ STSCI - update STScI printer queues
* 1990 Jan KDH @ STScI - interactive setting of NX, NY if 0 on input
* 1990 Feb KDH @ STScI - pgbegin option, add 3CLPS_POST queue
* 1991 Apr KDH @ STScI - Ant queue
* 1992 Mar Carole Hawell @ STScI - LP5 queue added
* 1992 Mar KDH @ STScI - distinguish VT125 and VT340 drivers
* 1992 Jul KDH @ Amsterdam - add amsterdam postscript queues
* 1992 Jul KDH @ cavad - add cavad canon queues
* 1992 Jul KDH @ cavad - add pericom_mg
* 1992 Oct KDH @ STScI - update STScI postscript queues
* 1993 Feb CAH @ STScI - update P&A postscript queues
* 1993 Apr KDH @ STScI - change SOL331 and ANT423S to 331PS1 and S423PS
* 1994 Jan KDH @ Utrecht - write(*,'(A)') '$...' -> write(*,'(A,$)') ' ...'
* 1994 Jan KDH @ Utrecht - add lpr for use on unix versions
* 1994 Jan KDH @ Utrecht - add XW devicetype
* 1994 Mar KDH @ Utrecht - eliminate de-anzas
* 1994 Mar KDH @ Utrecht - change LIB$SPAWN to SPAWN
* 1994 Mar KDH @ Utrecht - don't upcase filenames
* 1994 May KDH @ STScI - import Utrecht version, set OPSYS->VMS, SITE->STSCI
* 1995 Jan KDH @ St.Andrews - implement St-Andrews queues
* 1999 Sep KDH @ St.And - remove PGIDEN from PGSTOP
* 2000 Nov KDH @ Austin - implement Austin queues
* 2001 Jul KDH @ St.And - new St-Andrews queues
* 2001 Jul KDH @ StA - linux port: must use g77 -fno-automatic
* 2001 Jul KDH @ StA - simplify, strip vms & multi-sites
* 2001 Aug KDH @ StA - FULLPATH plot files, ASPECT
* 2004 Jan KDH @ StA - add err= to read (linux)
* 2005 Feb KDH @ StA - strip @ from plotfile
* 2006 May KDH @ StA - add farside printer
* 2010 Jun KDH @ StA - update StA printers
* 2014 Mar KDH @ StA - kdh -> kdh1
* 2014 Sep KDH @ StA - kdh1 -> kdh
* 2022 May KDH @ StA - update StA printers (kdh1,kdh1_uniprint)

        character*20 site
	character*1 aspect, queue, mode, oldmode
	character*256 reply
	character*64 terminal, devtype
	character*128 prtcom, plotdir, plotfile 
	character*128 psfile, cpsfile, giffile, pltfile, jfile
	character*256 task, path
	logical soft, spool, auto
	integer pgbegin
	common /pgstartmode/ mode
	data icall/0/
	icall = icall + 1

* first call initializations
      if( icall .eq. 1 ) then
	aspect = 'L'
	mode = 'H'
	queue = '?'
	iopen = 0
* page format
	nx0 = 1
	ny0 = 1

* site
	site = 'STANDREWS'
	write(*,*) 'PGSTART at ', SITE( :LEN1(SITE) )

* default directory
	plotdir='$SCR/'
* if system doesn't translate the environment variable $SCR
* correctly try the current directory as the default
c	plotdir = './'
	l1 = len1( plotdir )
        write(*,*) 'Default plot dir: ', plotdir(:l1)
	call fullpath( plotdir, path )
      if( plotdir(:l1) .ne. path(:len1(path)) ) then
	write(*,*) 'Default path: ', path(:len1(path))
      end if

* default plot files
	psfile = plotdir(:l1) // 'pgplot.ps'
	cpsfile = plotdir(:l1) // 'pgplot.cps'
	giffile = plotdir(:l1) // 'pgplot.gif'
	pltfile = plotdir(:l1) // 'pgplot.plt'

* end of first call initializations
      end if

* every call initializations
	oldmode = mode
	new = 0
	goto 1

* menu of device types
   10	write(*,*) 'PGSTART device types:'
	write(*,*) ' X=Xserve W=xWin T=Tektronix G=Gif H=Help'
	write(*,*) ' P=PostScript C=ColourPS N=null B=pgBegin'

* select device type
    1	write(*,'(A,$)')
     &	' PGSTART device type (H=help) [' // mode // '] '
	read(*,'(A)', err=1 ) reply
	call upcase( reply, reply )
	if( reply .eq. 'H' ) goto 10
	if( reply .eq. ' ' ) reply = mode

* soft devices (terminals) 
      if(    reply .eq. 'T'
     &	.or. reply .eq. 'X'
     &	.or. reply .eq. 'W'
     &	.or. reply .eq. 'B'
     &	.or. reply .eq. 'N'
     &	) then
	mode = reply
	soft = .true.
	spool = .false.
	auto = .false.

* hardcopy devices require a disk file
      else if( reply .eq. 'P'
     &	.or.   reply .eq. 'C'
     &	.or.   reply .eq. 'G' ) then
	mode = reply
	soft = .false.
	spool = mode .eq. 'P' .or. mode .eq. 'C'
	auto = .false.

* trap invalid replies
      else if( reply .ne. ' ' ) then
	goto 10
      else if( mode .eq. ' ' .or. mode .eq. 'H' ) then
	goto 10
      end if

* close previous plot if one was open
      if( mode .ne. oldmode .and. iopen .ne. 0 ) then
	call pgend
	iopen = 0
      end if

* terminals
      if( soft ) then
	if( mode .eq. 'N' ) terminal = '/NU'
	if( mode .eq. 'T' ) terminal = '/TE'
	if( mode .eq. 'W' ) terminal = '/XW'
	if( mode .eq. 'X' ) terminal = '/XSERV'
	if( mode .eq. 'B' ) terminal = '?'

* plot file 
      else

	jfile = pltfile
	if( mode .eq. 'P' ) jfile = psfile
	if( mode .eq. 'C' ) jfile = cpsfile
	if( mode .eq. 'G' ) jfile = giffile
      if( mode .ne. oldmode ) then
	plotfile = pltfile
	if( mode .eq. 'P' ) plotfile = psfile
	if( mode .eq. 'C' ) plotfile = cpsfile
	if( mode .eq. 'G' ) plotfile = giffile
      end if

	if( plotfile .ne. ' ' )
     &	write(*,*) '<CR> uses plotfile ', plotfile( :len1(plotfile) )
	write(*,*) 'JUNK autoprints ', jfile( :len1(jfile) )
	if( mode .eq. 'P' ) write(*,'(A,$)') ' PostScript file : '
	if( mode .eq. 'C' ) write(*,'(A,$)') ' ColourPostScript file : '
	if( mode .eq. 'G' ) write(*,'(A,$)') ' GIF file : '
 	read(*,'(A)') reply
      if( reply .eq. ' ' ) then
	reply = plotfile
      else
	if( reply(1:1) .eq. '@' ) reply = reply(2:)
	plotfile = reply
	call upcase( reply, reply )
	auto = .false.
      if( reply .eq. 'JUNK') then
	plotfile = jfile
	auto = .true.
      end if
	if( plotfile .eq. ' ' ) goto 10
      end if

      if( spool ) then

* print queues at each site
      if( site .eq. 'STANDREWS' ) then
   16	write(*,'(A,$)')
     &	' Which Queue (?,U,K,D,T,C) ? [' // queue // '] '
	read(*,'(A)') reply
	if( reply .eq. ' ' ) reply = queue
	call upcase( reply, reply )
	prtcom = ' '
      if( reply .eq. 'K') then
	prtcom = 'lpr -Pkdh1'
c	prtcom = 'lpr -Pkdh'
      else if( reply .eq. 'D') then
	prtcom = 'lpr -Pduplex'
      else if( reply .eq. 'T') then
	prtcom = 'lpr -PT644N'
      else if( reply .eq. 'C') then
	prtcom = 'lpr -PHPCP2025'
      else if( reply .eq. 'U') then
c	prtcom = 'lpr -PUniPrint_Ricoh'
	prtcom = 'lpr -Pkdh1_uniprint'
      else
        write(*,*) 'St.Andrews Printer queues :'
        write(*,*) ' U  kdh1_uniprint'
        write(*,*) ' K  kdh1    (Keith Horne''s office)'
        write(*,*) ' D  duplex  double-sided (Starlink room)'
        write(*,*) ' T  T644N   double-sided (south corridor)'
        write(*,*) ' C HPCP2025 colour (Starlink room)'
	goto 16
      end if
	queue = reply

* last site
      else
	write(*,*) '** PGSTART: no queues at ', site( :len1(site) )
	spool = .false.
	auto = .false.
      end if

      end if

* portrait vs landscape
	write(*,'(A,$)')
     &	' P(ortrait) or L(andscape) ? [' // aspect // '] '
	read(*,'(A)') reply
	call upcase( reply, reply )
	if( reply .eq. ' ' ) reply = aspect
	if( reply .eq. 'P' .or. reply .eq. 'L' ) aspect = reply

* device type
      if( aspect .eq. 'L' ) then
	if( mode .eq. 'P' ) devtype  =  'PS'
	if( mode .eq. 'C' ) devtype  =  'CPS'
	if( mode .eq. 'G' ) devtype  =  'GIF'
      else if( aspect .eq. 'P' ) then
	if( mode .eq. 'P' ) devtype  =  'VPS'
	if( mode .eq. 'C' ) devtype  =  'VCPS'
	if( mode .eq. 'G' ) devtype  =  'VGIF'
      end if

* end of hard devices
      end if

* need to open a new plot ?
      if( iopen .le. 0 ) then

* page format
      if( nx .le. 0 .or. ny .le. 0 ) then
   40	write(*,'(A,I3,A,I3,A,$)')
     #	' Plot page format (NX,NY) [', nx0, ' ,', ny0, ' ] '
	read(*,'(A)') reply
      if( reply .eq. ' ' ) then
	nx = nx0
	ny = ny0
      else
	read( reply,*,err=40, end=40 ) nx, ny
	IF( nx .le. 0 .or. ny .le. 0 ) goto 40
      end if
      end if
	nx0 = nx
	ny0 = ny

* open the plot

* terminals
      if( soft ) then
	write(*,*) 'PGBEGIN( 0, ',
     &	'''', terminal( :len1( terminal ) ), ''',',
     &	nx, ',', ny, ' )'
	ierr = pgbegin( 0, terminal, nx, ny )

* hardcopy devices
      else
* get full path
	call fullpath( plotfile, path )
* call pgbegin
	task = '"' // path( :len1( path ) ) // '"'
     &		// '/' // devtype( :len1( devtype ) )
	write(*,*) 'PGBEGIN( 0, ', task( :len1(task) ),
     &	', ', nx, ', ', ny, ' )'
	ierr = pgbegin( 0, task( :len1(task) ), nx, ny )
      end if
* failure
      if( ierr .ne. 1 ) then
	write(*,*) '** FAILED TO OPEN PLOT.'
	goto 1
      end if
* success
	iopen = 1
	new = 1
      end if

c	write(*,*) 'PGSTART spool', spool, ' auto', auto
	return

C* PGSTOP ... close plot, send output to printer queue
C+
	entry pgstop
C--
* write time and username at bottom of plot
c	call pgiden

* close the plot
	call pgend
	iopen = 0

* how to print the plot file
	write(*,*) 'PGSTOP spool', spool, ' auto', auto
      if( spool ) then
	call word1( prtcom, l1, l2 )
	l2 = len1( prtcom )
	write(*,*) 'COMMAND: ', prtcom( l1:l2 )
	call word1( plotfile, l3, l4 )
	write(*,*) 'PLOTFILE: ', plotfile( l3:l4 )
	call fullpath( plotfile, path )
	call word1( path, l5, l6 )
      if( plotfile( l3:l4 ) .ne. path( l5:l6 ) ) then
	write(*,*) 'PATH: ', path( l5:l6 )
      end if
	task = prtcom( l1:l2 ) // ' ' // path( l5:l6 )
	write(*,*) 'TASK: ', task( :len1( task ) )

* decide to spool, or not to spool
      if( .not. auto ) then
	write(*,'(A,$)') ' Queue file for plotting? [Y] '
	read(*,'(A)') reply
	call upcase( reply, reply )
	if( reply .eq. ' ' ) reply = 'Y'
      else
	reply = 'Y'
      end if

* print the plot file
      if( reply .eq. 'Y' ) then
	call spawn( task( : len1( task ) ) )
      end if
      end if

	return
	end
*=====================================
	logical function pgagain()
* plot again ?
* 1997 Aug Keith Horne @ St.Andrews
* 1999 Oct KDH @ St.And - make logical
	character*10 reply
	write(*,'(A,$)') ' Plot again ? [N] '
	read(*,'(A)') reply
	call upcase( reply, reply )
	pgagain = .false.
	if( reply.eq.'Y') pgagain = .true.
	return
	end
*----------------------------------------
	subroutine pgbgw
* set background/foreground colours to white/black
* 2010 Oct Keith Horne @ St-And
	call pgscr( 0, 1., 1., 1. )
	call pgscr( 1, 0., 0., 0. )
	return
	end
