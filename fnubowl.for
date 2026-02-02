
! -----------------------------------------------------------
	subroutine fnudisk_inclined( nr, r, h, t, dmpc, cosi, nw, w, fmjy )
! 2024-Jan KDH @ StA : retain for backwards compatibility
	real*4 w( nw ), fmjy( nw ), r( nr ), h( nr ), t( nr )	
	call fnubowl( nr, r, h, t, dmpc, cosi, 1., nw, w, fmjy )
	return
	end
! -----------------------------------------------------------
	subroutine fnubowl( nr, r, h, t, dmpc, cosi, fcol, nw, w, fmjy )
! compute spectrum of inclined axi-symmetric blackbody disk
! assumes all disk elements visible from earth - no self-occultation
! In: 	nr	i4 number of radii
!	r(nr)	r4 radii (light days)
!	h(nr)	r4 height above disk plane (light days)
!	t(nr)	r4 effective temperature (Kelvin)
!	dmpc	r4 luminosity distance (Mpc)
!	cosi	r4 cosine of disk inclination (1=face on)
!	fcol	r4 colour temperature Tcol = fcol Teff
!	nw	i4 number of wavelengths
!	w(nw)	r4 wavelengths (Angstroms)
! Out:	fmjy(nw) r4 disk spectrum (mJy)
!
! 2021-Aug KDH @ Crail - adapt from D.Starkey subroutine fnudisk.
! 2021 Nov KDH @ Crail - verbose
! 2024-Jan KDH @ StA - implement fcol and f90->f77
c	integer, intent(in) :: nw, nr
c	real,    intent(in) :: w( nw ), r( nr ), h( nr ), t( nr ), dmpc, cosi
c	real,   intent(out) :: fmjy( nw )
c	logical ::  ok , verbose
	real*4 w( nw ), fmjy( nw )
	real*4 r( nr ), h( nr ), t( nr )
	logical ok, verbose
	double precision sum

	verbose = .false.

! sanity check
	ok = nr .ge. 2 .and. nw .ge. 1 .and. dmpc .gt. 0.
      if( .not. ok ) then
	write(*,*) '** FNUBOWL INPUT ERRORS :(.'
	write(*,*) '** Nr', nr, ' Nw', nw, ' D(Mpc)', dmpc
	stop
      end if

! report
      if( verbose ) then
	write(*,*) 'R:', nr, r(1), r(nr)
	write(*,*) 'H:', nr, h(1), h(nr)
	write(*,*) 'T:', nr, t(1), t(nr)
	write(*,*) 'W:', nw, w(1), w(nw)
	write(*,*) 'Dmpc', dmpc, ' cosi', cosi, ' fcol', fcol
      end if

! constants
	twopi = 8. * atan2( 1., 1. )
	pc = cgs( 'PC' )
	clight = cgs( 'C' )
	day = 24. * 3600.
! light days per Mpc
	pcM2ld = 1.e6 * pc / clight / day
! distance in light days
	d = dmpc * pcM2ld
! mJy per (erg/cm2/s/Hz)
	fnu2mjy = 1.e26
! units factor for bnu/cgs, area/ld^2, and fnu/mJy
	units = fnu2mjy / d / d
! divide Bnu( fcol Teff ) by factor f4
	f4 = fcol ** 4

      if( verbose ) then
	write(*,*) 'Mpc/ld', pcm2ld
	write(*,*) 'D(Mpc)', dmpc, ' D(cm)', d
	write(*,*) 'fnu/mJy', fnu2mjy, ' units', units
	write(*,*) 'fcol', fcol, ' fcol^4', f4
      end if

! unit vector pointing to Earth (in the x-z plane)
	ci = max( 0., min( 1., cosi ) )
	si = sqrt( 1. - ci * ci )

	if( verbose ) write(*,*) 'sini', si, ' cosi', ci

	ex = si
!	ey = 0.
	ez = ci
! number of azimuths (1 degree steps)
	naz = 365
! full ring if face on
	if( abs( ci ) .eq. 1. ) naz = 1
! azimuth step in radians
	daz = twopi / naz

	if( verbose ) write(*,*) 'Azimuths', naz, daz, ' 2pi', naz*daz
! wavelengths
      do iw = 1, nw
	wnow = w( iw )
! zero the sum
	sum = 0.d0
! radius grid
      do ir = 1, nr - 1
! radius and height above disk plane (light days)
	rnow = r( ir )
	hnow = h( ir )
! area of tilted surface element (light day)^2
	ip = ir + 1
	dr = r( ip ) - rnow
	dh = h( ip ) - hnow
	ds = sqrt( dr * dr + dh * dh )
	da = ds * rnow * daz
! inward tilt
	sintilt = dh / ds
	costilt = dr / ds
! kelvin temperature and planck function (Bnu in cgs units)
	tnow = t( ir )
c	bb = bnu( wnow, tnow ) 
	tc = fcol * tnow
	bb = bnu( wnow, tc ) / f4
! azimuths (0 = x axis)
      do iaz = 1, naz
	az = iaz * daz
	cosaz = cos( az )
	sinaz = sin( az )
! unit vector outward normal to tilted surface (py not used)
	px = -sintilt * cosaz
!	py = -sintilt * sinaz
	pz =  costilt
! foreshortening factor
!	dot = ex * px + ey * py + ez * pz
! omit ey * py since ey=0
	dot = ex * px + ez * pz
! sum blackbody times projected area if facing toward Earth
! assume no self-occultations
      if( dot .gt. 0. ) then
	sum = sum + bb * da * dot
      end if
! next azimuth
      end do
! next radius ir
      end do
! convert fnu to mjy and stow
	fmjy( iw ) = sum * units
! next wavelength iw
      end do
	return
	end subroutine
	
