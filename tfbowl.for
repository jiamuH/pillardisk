! --------------------------------------------------------------
	subroutine tfb_inclined( nr, r_r, h_r, T_r, 
     &		hlamp, cosi, wrest, ntau, tau, psi )
* 2024 Jan KDH @ StA - retain for backwards compatibility.
	real*4 r_r( nr ), h_r( nr ), t_r( nr ), tau( ntau ), psi( ntau )
	call tfbowl( nr, r_r, h_r, T_r,
     &		hlamp, cosi, 1., wrest, ntau, tau, psi )
	return
	end
! --------------------------------------------------------------
	subroutine tfbowl( nr, r_r, h_r, T_r, 
     &		hlamp, cosi, fcol, wrest, ntau, tau, psi )
! normalised delay distribution for inclined blackbody disk
! assumes no self-occultations
! In:	nr	i4 number of annuli
!	r_r(nr)	r4 radius at inner edge (light days)
!	h_r(nr)	r4 height at inner edge (light days)
!	T_r(nr)	r4 Teff at inner edge (Kelvin)
!	hlamp	r4 lamppost height (light days)
!	cosi	r4 cos(inclination) (1=face-on)
!	fcol	r4 colour temperature boost fcol = Tc/Teff
!	wrest	r4 rest-frame wavelength (Angstroms)
!	ntau		i4 number of delays
!	tau(ntau)	r4 time delays (days) (even spacing)
! Out:	psi(ntau)	r4 delay distribution (Int psi dtau = 1)
!
! 2021 Aug Keith Horne @ Crail - adapt from David Starkey's tfb subroutine
! 2021 Oct KDH @ Crail - discshad and i_r
! 2024 Jan KDH @ StA - implement fcol  and F90 => F77
! 2024 Mar KDH @ StA - trap nan psi(i), invalid sum, ...
! 2024 Mar KDH @ StA - cold T0/15 => T0/100 fixes Psi(tau)=0 at short lam.
! 2024 Apr KDH @ StA - use Xref to prevent Psi(tau)=0 at short lam.
! 2024 Apr KDH @ StA - report dark annuli
! 2024 Apr KDH @ StA - Loop 2 avoids panic crash when psi(tau)=0
! 2025 Jan KDH @ StA - enumerate PANIC messages

	real*4 r_r( nr ), h_r( nr ), T_r( nr )
	real*4 tau( ntau ), psi( ntau )

! irradiation (1) or not (0)
	parameter( maxr = 10000 )
	integer*4 i_r( maxr )

	double precision sum, top, bot
	logical verbose, ok, yes, yesr, yesaz, panic
	logical firstcall / .true. /


! toggle diagnostics
	panic = .false.
	verbose = .true.
	verbose = .false.


      if( verbose ) then
	write(*,*) '--- tfbowl --------------', nr, ntau, nint(wrest)
      end if

! sanity checks
	ok = nr .gt. 1 .and. ntau .gt. 1 .and. wrest .gt. 0.
	ok = ok .and. tau( 1 ) .lt. tau( ntau )
	ok = ok .and. r_r( 1 ) .lt. r_r( nr )
      if( .not. ok .or. verbose ) then
	if( .not. ok ) write(*,*) '** PANIC 1 in tfbowl. INPUT ERRORS :(.'
	write(*,*) ' radii', nr, r_r( 1 ), r_r( nr )
	write(*,*) 'delays', ntau, tau( 1 ), tau( ntau )
	write(*,*) 'rest wavelength', wrest
	if( .not. ok ) stop
      end if

      if( nr .gt. maxr ) then
	write(*,*) '** PANIC 2 in tfbowl. Nr', nr, ' max', maxr
	stop
      end if

! need tmn and tmx for xmn, xmx below
c 	call mnmx(   nr, t_r, tmn, tmx )
	call imnmx( nr, t_r, itmn, itmx )
	tmn = t_r( itmn )
	tmx = t_r( itmx )

      if( verbose ) then
	call mnmx(   nr, r_r, rmn, rmx )
	call mnmx( ntau, tau, dmn, dmx )
	write(*,*) 'Nr', nr, ' r ', rmn, rmx, ' T', tmn, tmx
	write(*,*) ' max T(', itmx, ' )', nint( t_r( itmx ) ), nint( tmx )
	write(*,*) ' min T(', itmn, ' )', nint( t_r( itmn ) ), nint( tmn )
	write(*,*) 'tau', dmn, dmx, ' lam', wrest
      end if

! read clock
	call system_clock( istartclock )

! random number seed
      if( firstcall ) then
	iseed = 926056
	firstcall = .false.
      end if
! constants
	twopi = 8. * atan2( 1., 1. )
! physics constants (cgs units)
	c  = cgs( 'C' )	! light speed
	bk = cgs( 'K' )	! Boltzmann k
	hp = cgs( 'H' )	! Planck h
	sb = cgs( 'SIGMA' )	! Stefan-Boltzman 
! astro constants (cgs units)
	angst = cgs( 'ANGSTROM' )
! Wien temperature (x = hnu/kT  = T0/T)
	T0 = hp * c / ( bk * wrest * angst )

! range of dimensionless frequency x and dB/dx
	t = t0 / fcol
	xmn = t / tmx
	xmx = t / tmn
	dBmx = xmn * xmn / ( cosh( xmn ) - 1. )
	dBmn = xmx * xmx / ( cosh( xmx ) - 1. )
	xref = 1.

! detect underflow problems
	panic = dBmx .le. 0.

      if( panic ) then
	call mnmx( nr, r_r, rmn, rmx )
	write(*,*) '** PANIC 3 in tfbowl. wrest', wrest, ' Nr', nr
	write(*,*) '(',  rmn, ' <   R    <',  rmx, ' )'
	write(*,*) '(',  tmn, ' <   T    <',  tmx, ' )'
	write(*,*) '(',  xmn, ' < hnu/kT <',  xmx, ' )'	
	write(*,*) '(', dBmn, ' < dB/dT  <', dBmx, ' )'	

* try normalising to xref
* KDH : OK BECAUSE DELAY MAP IS NORMALISED ?
c	xref = sqrt( xmx * xmn )
	xref = xmx
	top = cosh( xref ) - 1.
	bot = cosh( xmn ) - 1.
	wien = top / bot
	if( wien .eq. 0. .or. wien .ne. wien ) wien = exp( xref - xmn )
	x1 = xmn
	wien1 = wien
	dBmx = x1 * x1 * wien1

	bot = cosh( xmx ) - 1.
	wien = top / bot
	if( wien .eq. 0. .or. wien .ne. wien ) wien = exp( xref - xmx )
	x2 = xmn
	wien2 = wien
	dBmx = x2 * x2 * wien2

	write(*,*) 'Normalise X to', xref, ' wien', wien1, wien2
	write(*,*) '(',  x1,  ' < hnu/kT <',  x2, ' )'	
	write(*,*) '(', dBmn, ' < dB/dT  <', dBmx, ' )'	
	if( dBmx .le. 0. ) stop
	write(*,*) 'Should now be OK.'
* KDH : SWITCH ON DEBUG DIAGNOSTICS
c	verbose = .true.
      end if

! neglect blackbody when T < cold
c	cold = T0 / 15.
	cold = T0 / 100.
	cold = T0 / 200.
	cold = min( cold, tmn / 10. )
* KDH : INCREASE COLD FOR SPEED, BUT MAY INDUCE PANIC CRASH
! Use Bnu( fcol * Teff ) / fcol**4
	f4 = fcol ** 4

      if( verbose .or. panic ) then
	write(*,*) 'rest wavelength', wrest, ' Tcold', cold
      end if

! mean delay step
	dtau = ( tau( ntau ) - tau( 1 ) ) / max( 1,  ntau - 1 )

! check uniform spacing
	dt1 = dtau
	dt2 = dt1
      do i = 2, ntau
	dt = tau( i ) - tau( i - 1 )
	dt1 = min( dt1, dt )
	dt2 = max( dt2, dt )
      end do

      if( verbose ) then
	write(*,*) 'Delays', ntau, tau(1), tau(ntau), dtau
	write(*,*) 'Delay spacing', dt1, dt2
      end if

! report and abort if not uniform
      if( dt2 - dt1 .gt. dtau / 100. ) then
	write(*,*) '** PANIC 4 in tfbowl. NON-UNIFORM DELAY GRID.'
	write(*,*) '** delays', ntau, tau(1), tau(ntau), dtau
	write(*,*) '** delay range', dt1, dt2
	stop
      end if 

! rms of delay blur (prevents  undersampling)
	sigtau = dtau
	chimax = 5.0

! cos and sin of inclination
	ci = min( 1., max( 0., cosi ) )
	si = sqrt( 1. - ci * ci )

! KDH : MAY NEED A DOUBLE PRECISION SUM
      do i = 1, ntau
	psi( i ) = 0.0
      end do

! disc shadows 
c	allocate ( i_r( nr ) )
	call diskshad( nr, r_r, h_r, hlamp, i_r )
	n = 0
      do i = 1, nr
	if( i_r( i ) .eq. 0 ) n = n + 1
      end do
      if( n .gt. 0 ) then
	p = n * 100. / nr
      if( verbose ) then
	write(*,*) 'Shadowed annuli', n, ' of', nr, ' =', p, ' %'
      end if
      end if

* SECOND LOOP OCCURS IF TRANSFER FUNCTION VANISHES
      DO LOOP = 1, 2

* all annuli
	ir1 = 1
	ir2 = nr - 1
* KDH : requires nr > 1

* hot annulus
      if( loop .eq. 2 ) then
	verbose = .true.
	ir1 = max( 1, itmx - 1 )
	ir2 = ir1 + 1
	write(*,*) 'LOOP 2 : Nr', nr, ' Using', ir1, ir2
	write(*,*) 'Min T(', itmn, ' )=', nint( t_r( itmn ) )
	write(*,*) 'Max T(', itmx, ' )=', nint( t_r( itmx ) )
	write(*,*) ' T(', ir1, ' )=', nint( t_r( ir1 ) )
	write(*,*) ' T(', ir2, ' )=', nint( t_r( ir2 ) )
      do i = ir1, ir2
	write(*,*) ' T(', i, ' )=', nint( t_r( i ) )
      end do
      end if

* dark annulus counter
	ndark = 0
      do ir = ir1, ir2
* compute dpsi0 for this annulus
! skip shadowed annuli
	if( i_r( ir ) .le. 0 ) cycle
	i = ir + 1
	yesr = ir .le. 10 .or. i .eq. nr .or. mod(ir, max(1,nr-1)) .eq. 0
	yesr = .true.
	if( yesr .and. verbose ) write(*,*) '----------------------'

* T_R(I) IS TEMP AT OUTER EDGE OF ANNULUS I
! mean temperature
	T = ( T_r( ir ) + T_r( i ) ) / 2.
* KDH : COULD INTERPOLATE TO RANDOM RADII ON THE ANNULUS
! skip if cold or invalid
      if( T .le. cold  .or. T .ne. T ) then
      if( verbose .or. panic .or. .true. ) then
	write(*,*) ir, ' R', r_r(ir), ' T', nint(T), '< cold', nint(cold)	
      end if
	cycle
      end if
! inner radius and height 
	rin = r_r( ir )
	hin = h_r( ir )
! outward change
	dr = r_r( i ) - rin
	dh = h_r( i ) - hin
! midpoint
	rmid = rin + dr / 2.
! inward tilt
	dl = sqrt( dr * dr + dh * dh )
	sintilt = dh / dl
	costilt = dr / dl
! azimuth step (approx square surface elements)
	daz = dr / rmid
! number of azimuths
	naz = nint( twopi / daz )
	naz = max( 100, naz )
	daz = twopi / naz

! pre-evaluate azimuth-independent quantities
	dsa0 = rmid * dl * daz
! boost Teff by factor fcol to obtain colour temperature
	Tcol = fcol * t
	x = T0 / Tcol
* KDH : dBdT tiny for large x, dpsi0 then underflows :-(
	top = cosh( xref ) - 1.
	bot = cosh( x ) - 1.
	wien = top / bot
      if( wien .le. 0. .or. wien .ne. wien ) then
	old = wien
	wien = exp( xref - x )
c	write(*,*) ' Tcol', nint(Tcol), ' X', x, ' Wien', old, ' =>', wien
      end if
c	x = x / xref
	dBdT = x * x * wien
* scale to max ( hope to avoid a null response )
	dBdT = dBdT / dBmx
	dpsi0 = dBdT / T / T / T

* possible alternative
      if( dpsi0 .eq. 0. ) then
	botln = alog( cosh( x ) - 1. )
	x_t = x / t
	top = x_t * x_t / t
      end if

      if( verbose .and. yesr .or. ( x .ne. x .or. dBdT .ne. dBdT ) ) then
	write(*,*) 'R', rin, ' dr', dr, ' Naz', naz, ' dsa0', dsa0
	write(*,*) 'H', hin, ' dh', dh, ' dl', dl, ' sintilt', sintilt
	write(*,*) 'T', nint( t ) , ' X', x, ' dBdT', dbdt, ' dpsi0', dpsi0
      end if
* end calculation of dpsi0

* force dpsi = 1 on second loop
      if( loop .gt. 1 ) then
	write(*,*) 'LOOP', loop, ' dpsi', dpsi, ' =>', 1.
	dpsi0 = 1.
      end if

* azimuths
      do iaz = 1, naz
	yesaz = iaz .eq. 1 .or. iaz .eq. naz/2 .or. iaz .eq. naz/4
	azlo = ( iaz - 1 ) * daz

! random samples on each surface element
! helps to prevent beating with the delay grid
! number of random samples (increase if needed)
	nran = 10
      do iran = 1, nran
! random radius and corresponding height
	u = ran3( iseed )
	r = ( rin + u * dr )
	h = ( hin + u * dh )
! random azimuth
	u = ran3( iseed )
	az = azlo + u * daz
! KDH : assumes dr/r is small
	caz = cos( az )
!	saz = sin( az )  ! not used
! cartesian coords relative to lamp
	dx = r * caz
!	dy = r * saz     ! not used
	dz = h - hlamp
! distance from lamp
	ddlamp = dz * dz + r * r
	dlamp = sqrt( ddlamp )
! p = unit vector toward lamp
	pz = -dz / dlamp
	pr = - r / dlamp 
	costheta = -pr * sintilt + pz * costilt

! skip if tilted away from the lamp
      if( costheta .le. 0. ) then
      if( verbose ) then
	deg = 360. * az / twopi
	write(*,*) '** WARNING: In shadow at r=', r, ' deg', deg
	write(*,*) 'dr', dr, ' dh', dh, 'costilt', costilt, ' sintilt', sintilt
	write(*,*) 'dz', dz, ' dlamp', dlamp, ' pz', pz, ' pr', pr
	write(*,*) 'costheta', costheta
C	stop
      end if
	cycle
      end if

! foreshorten by e dot n (n=normal vector, e=earth vector)
!	nx = - sintilt * caz   ex = si
!	ny = - sintilt * saz   ey = 0
!	nz =   costilt         ez = ci
	edotn = - sintilt * caz * si + costilt * ci

! skip if facing away from Earth
	if( edotn .le. 0. ) cycle

! KDH : test for self-occultation could go here
! code currently assumes no self-occultations
! could check for occultation by disk rim

! area presented to Earth (solid angle is dsa/(4 pi D^2)
!	dsa0 = rmid * dl * daz   ! dsa0 is pre-calculated
	dsa = dsa0 * edotn

! time delay in days
!	tdelay = dlamp - r * caz * si - dz * ci
	tdelay = dlamp - dx * si - dz * ci

!  spectrum : Fnu = Int Inu dOmega
!	Inu = Bnu( f T ) / f^4
!  response : psi(tau) = Int dI/dB dB/dT dT/dL delta(tau) dOmega
!
!  x = (hc/kfTw)  dx/dT = -x/T = -(kfw/hc) x^2
!  Bnu = (2kT/w^2)x/(exp(x)-1)
!      = (2hc/w^3f)/(exp(x)-1)
!  dB/dx = -(2hc/w^3f) exp(x)/(exp(x)-1)^2
!        = -(2hc/w^3f)/(exp(x)-2+exp(-x))
!        = -(hc/w^3f)/(cosh(x)-1)
!  dB/dT = (dB/dx)(dx/dT)
!        = (kw/hc)(hc/w^3f) x^2/(cosh(x)-1)
!        = (k/w^2f) x^2/(cosh(x)-1)
!    T^4 = T_v^4 + L (1-a) cos(theta) / ( sig 4 pi dlamp^2 )
!   d(T^4)/dL = (1-a) cos(theta) / ( sig 4 pi dlamp^2 )
!   d(T^4)/dT = 4 T^3
!   dT/dL = (d(T^4)/dL)(dT/dT^4)
!         = (1-a)/(sig pi) cos(theta) / (T^3 dlamp^2)
!
! For normalised Psi, omit all constant factors:
!    dIdB => 1  (assume fcol independent of T)
!    dBdT => x^2 / ( cosh( x ) - 1. )
!    dTdL => cos(theta) / ( T^3 dlamp^2 )
!  dOmega => dsa
!    dpsi => dBdT * dTdL * dsa / f^4

! dpsi0 = dBdT / T^3 is pre-calculated
      if( dpsi0 .ne. 0. ) then
	dpsi = dpsi0 * costheta / ddlamp * dsa
      else
	dpsi = top * costheta / ddlamp * dsa / exp( botln )
      end if

! average over nran random samples
	dpsi = dpsi / nran

! trap NAN 
	yes = mod( iaz, max(1,naz/12) ) .eq. 0 .or. iaz .eq. 1 .or. iaz .eq. naz 
	yes = verbose .and. yes
      if( dpsi .ne. dpsi .or. yes ) then
	deg = 360. * az / twopi
	write(*,*) '** PANIC 5 in tfbowl. wrest', nint( wrest ), ' dPsi', dpsi
	write(*,*) 'ir', ir, ' deg', nint(deg), ' iran', iran
	write(*,*) 'dB/dT', dBdT, ' dPsi0', dpsi0, ' X', t0 / Tcol
	write(*,*) 'T', T, 'T^3', T*T*T, ' dlamp^2', ddlamp
	write(*,*) 'dOmega', dsa, ' cos(theta)', costheta
	write(*,*) 
      if( dpsi .ne. dpsi ) then
	dpsi = 0.
	WRITE(*,*) '** TO AVERT CRASH, SET DPSI =', DPSI
      end if
	if( dpsi .ne. dpsi ) stop
      end if

* largest contribution at this annulus
      if( iran .eq. 1 .and. iaz .eq. 1 ) then
	dbdtmax = dbdt
	dpsimax = dpsi
      else
	dbdtmax = max( dbdtmax, dbdt )
	dpsimax = max( dpsimax, dpsi )
      end if

* KDH : SKIP IF DPSI = 0 ?
      IF( DPSI .EQ. 0. ) THEN
* 2024-04-02 KDH : DEBUG REPORT
      if( ( iran .eq. 1 .or. iran .eq. nran ) .and. verbose ) then
	write(*,*)
	write(*,*) ir, iaz, iran, ' T', nint(t), dpsi
	write(*,*) 'dB/dT', dBdT, ' dPsi0', dpsi0, ' X', t0 / Tcol
	write(*,*) 'T', T, 'T^3', T*T*T, ' dlamp^2', ddlamp
	write(*,*) 'dOmega', dsa, ' cos(theta)', costheta
	write(*,*)
      end if
      ELSE

! deposit dpsi in the delay grid
! gaussian approximates delta( tau )
! pixel range - may extrapolate (1,ntau)

* KDH : tau0 extrapolates tau(i) to i=0
	tau0 = tau( 1 ) - dtau
	dt = tdelay - tau0
	add = sigtau * chimax
* KDH : FLOOR ROUNDS DOWN AND CEILING ROUNDS UP ?
	i1 =   floor( ( dt - add ) / dtau ) 
	i2 = ceiling( ( dt + add ) / dtau )

* KDH : (i1,i2) can extend outside (1,ntau)
C      IF( I1 .LT. 1 .OR. I2 .GT. NTAU ) THEN
C	WRITE(*,*) ' I1', I1, ' I2', I2, ' NTAU', NTAU
C	WRITE(*,*) ' DT', DT, ' ADD', ADD, ' DTAU', DTAU
C      END IF

	sum = 0.d0
      do k = 1, 2
! k=1 : sum gaussian samples for normalisation
! k=2 : stow in the psi(tau) grid
      do i = i1, i2
! delay samples
      if( i .ge. 1 .and. i .le. ntau ) then
	tgrid = tau( i ) 
! extrapolate the tau grid
      else
	tgrid = tau0 + i * dtau
      end if
! gaussian weights
	chi = ( tdelay - tgrid ) / sigtau
	g = exp( -0.5 * chi * chi )
! normalisation sum
      if( k .eq. 1 ) then
	sum = sum + g
! increment psi
      else if( i .ge. 1 .and. i .le. ntau ) then
	add = dpsi * g
	psi( i ) = psi( i ) + add

* trap NAN
      if( psi( i ) .ne. psi( i ) ) then
	write(*,*) '** PANIC 6 in tfbowl. psi(i)', psi(i)
	write(*,*) 'tdelay', tdelay, ' tgrid', tgrid, ' sigtau', sigtau
	write(*,*) ' chi', chi, ' g', g, ' dpsi', dpsi
	psi( i ) = 0.
	write(*,*) '** TO AVERT CRASH, SET PSI(', I, ' )= ', psi( i )
c	stop
      end if

! KDH : COULD SUM THE EXTRAPOLATED PART HERE
      else

      end if
! next delay i
      end do

* trap NAN
c      IF( SUM .EQ. 0.D0 .OR. SUM .NE. SUM ) THEN
      IF( SUM .NE. SUM ) THEN
	WRITE(*,*) '** PANIC 7 in tfbowl. sum', sum
	WRITE(*,*) 'i1', i1, ' i2', i2, ' dpsi', dpsi
	SUM = 0.
	write(*,*) '** TO AVERT CRASH, SET SUM = ', sum
	STOP
      END IF

* exit if no contribution
	if( sum .eq. 0.d0 ) exit
! divide by normalisation sum
	dpsi = dpsi / sum
! next k = 1, 2
      end do

* report
      if( verbose .and. yesr .and. yesaz ) then
	deg = 360. * az / twopi
	write(*,*) 'R', r, ' deg', deg, ' tau', tdelay, ' dpsi', dpsi
	write(*,*) 'tau(', i1, i2, ') add', add, ' sum', sum
      end if

* END IF DPSI > 0
      END IF

! next random sample iran
      end do
! next azimuth iaz
      end do 

* count and report dark annuli
      if( dbdtmax .le. 0. .or. dpsimax .le. 0. ) then
      if( ndark .eq. 0 ) then
	t1 = t
	t2 = t
      end if
	ndark = ndark + 1
	t1 = min( t1, t )
	t2 = max( t2, t )
      if( verbose ) then
	write(*,*) '** annulus', ir, ' R', r_r(ir),
     &		' T', nint(t_r(ir)), ' dark at', nint( wrest )
      end if
      end if

! next radius ir
      end do

* report number of dark annuli
      if( verbose ) then
      if( ndark .gt. 0 ) then
	write(*,*) 'W', nint( wrest ), ' T', nint( t1 ), nint( t2 ),
     &		' dark', ndark, ' of', nr
      end if
      end if

! report
      if( verbose ) then
	call system_clock( iendclock )
	write(*,*)'Run time...', iendclock - istartclock
! mean delay
	top = 0.d0
	bot = 0.d0
      do i = 1, ntau 
	top = top + psi( i ) * tau( i )
	bot = bot + psi( i )
      end do
	write(*,*) 'Rest', wrest, ' mean delay', top / bot
	write(*,*) 'top', top, ' bot', bot
! KDH : COULD REPORT LOST PART OUTSIDE (1,NTAU)
!	stop
      endif

! sum response over the delay grid
	sum = 0.d0
      do i = 1, ntau
* trap NAN
	if( psi( i ) .ne. psi( i ) ) psi( i ) = 0.
	sum = sum + psi( i )
      end do
	sum = sum * dtau
! KDH : COULD TRY TO ADD RESPONSE OUTSIDE (1,NTAU)


	if( loop .gt. 1 ) write(*,*) 'LOOP', loop, ' SUM', sum

* normalise as probability distribution
      if( sum .gt. 0.d0 ) then
	top = 0.d0
	bot = 0.d0
      do i = 1, ntau
	psi( i ) = psi( i ) / sum
* trap NAN
	if( psi( i ) .ne. psi( i ) ) psi( i ) = 0.
	top = top + psi( i ) * tau( i )
	bot = bot + psi( i )
      end do
	avg = top / bot
* successful return
	if( avg .eq. avg ) return
      end if
	
* TRAP SUM .le. 0 or mean lag NAN
* THIS CAN OCCUR FOR SHORT WAVELENGTHS AND LOW TEMPERATURES
	WRITE(*,*) '** LOOP', loop
	WRITE(*,*) '** PANIC 8 in tfbowl: sum', sum, ' not > 0'
	WRITE(*,*) '** top', top, ' bot', bot, ' avg', avg
      DO I = 1, NTAU, MAX( 1, NTAU / 20 )
	WRITE(*,*) I, TAU( I ), PSI( I )
      END DO
	WRITE(*,*) 'NTAU', NTAU, ' DTAU', DTAU, ' DPSI', DPSI
	write(*,*) 'Nr', nr, ' wrest', wrest, ' fcol', fcol
	write(*,*) 'Hx', hlamp, ' cosi', cosi
	i = max( 1, nr / 2 )
	write(*,*) ' R:', r_r( 1 ), r_r( i ),  r_r( nr )
	write(*,*) ' T:', t_r( 1 ), t_r( i ), t_r( nr )
	write(*,*) ' H:', h_r( 1 ), h_r( i ), h_r( nr )

* next loop  ( TRY AGAIN )
      END DO

	write(*,*) '** LOOP 2 failed to fix PSI(TAU)=0.'
	r = r_r( nr )
	call lookup( ntau, tau, r, i1, i2, p )
	write(*,*) 'Nr', nr, ' R', r_r(1), r_r( nr ), ' dtau', dtau
	write(*,*) ' tau(', 1, ' )=', tau( 1 ), ' tau(', ntau, ' )=', tau( ntau )
	write(*,*) ' tau(', i1, ' )=', tau( i1 ), ' tau(', i2, ' )=', tau( i2 )
	write(*,*) ' psi(', i1, ' )=', psi( i1 ), ' psi(', i2, ' )=', psi( i2 )
	psi( i2 ) = 1.
	if( dtau .gt. 0. ) psi( i2 ) = 1. / dtau
	write(*,*) 'TO AVERT CRASH, PSI(', i2, ' ) <=', psi( i2 )
	write(*,*) ' psi(', i1, ' )=', psi( i1 ), ' psi(', i2, ' )=', psi( i2 )

* second loop failed
c	STOP 'psi(tau) still vanishes :(.'
	return
	end subroutine

