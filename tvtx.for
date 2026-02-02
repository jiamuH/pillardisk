*==================================================
	subroutine tvtx( nr, r_r, h_r, hlamp, risco, r1, alp, qin,
     &		tv1, tx1, tv_r, tx_r, tq_r, t_r )
* Teff(r) for lamp-irradiated viscous curved accretion disc
* in:	nr	i4 number of annuli
*	r_r(nr)	r4 radii (light days)
*	h_r(nr)	r4 height above disc plane (light days)
*	hlamp	r4 height of lamppost (light days)	
*	risco	r4 inner disc radius (light days)
*	r1	r4 reference radius (light days)
*	alp	r4 power-law index Tv = Tv1 (R1/R)^alp
*	qin	r4 inner torque
*	tv1	r4 (Kelvin) Tv1^4 = ( 3 G M dM/dt ) / ( 8 pi sigma R1^3 )
*	tx1	r4 (Kelvin) Tx1^4 =  L (1-a) / (4 pi sigma R1^2 )
* out:	tv_r(nr) r4 (Kelvin) viscous temperature
*	tx_r(nr) r4 (Kelvin) irradiation temperature
*	tq_r(nr) r4 (Kelvin) torque temperature
*	t_r(nr)	r4 (Kelvin) effective temperature profile 
* 2021 Aug Keith Horne @ St Andrews
* 2022 Aug KDH @ St.A - prevent tv = NAN
* 2023 Jun KDH @ Crail - fix tv bug
* 2023 Nov KDH @ StA - torque parameter qin
* 2023-Dec-09 KDH @ StA - torque heating tq_r(nr)
	real*4 r_r( nr )
	real*4 h_r( nr )
	real*4 tv_r( nr )
	real*4 tx_r( nr )
	real*4 tq_r( nr )
	real*4 t_r( nr )
*
c	write(*,*) '** tvtx(', nr, risco, r1, alp, nint( tv1 ), nint( tx1 ), ' )'

* reference temperature (avoids T^4 overflow) 
	t0 = 1.e4
* floor temperature
	tcmb = 2.7
* temperature profile
      do i = 1, nr
	r = r_r( i )
* TRAP DIVIDE BY 0
      if( r .le. 0. ) then
	write(*,*) 'R(', i, ') = 0 in tvtx.'
	stop
      end if
* power-law temperature profile (thin,steady,keplerian disc)
	tv = tv1 * ( r1 / r ) ** alp
* torque at inner edge heats Q sqrt(Risco/R)
	root = sqrt( risco / r )
	torque = max( 0., qin * root )
* sub-Keplerian shear at inner edge cools (less viscous heating)
c	cool1 = max( 0., 1. + ( qin - 1. ) * root )
	cool1 = max( 0., 1. - root + torque )
	coolq = max( 0., qin * root )
* curved surface has larger surface area
	h = h_r( i )
	ip = min( i + 1, nr )
	im = max( ip- 1, 1 )
	n = ip - im
* KDH : SHOULD R(I) BE INNER RADIUS OF ANNULUS ?
	dr = ( r_r( ip ) - r_r( im ) ) / n
	dh = ( h_r( ip ) - h_r( im ) ) / n
	ds = sqrt( dh * dh + dr * dr )
	cool2 = dr / ds
* viscous heating
* KDH: 2023 revise to include cool1 * cool2
* Qv = (Tv1/T0)^4(R1/R)^(4alpha) (1+(Q-1)sqrt(Risco/R)) ( dr /ds )
	tv4 = ( tv / t0 ) ** 4
	qv = tv4 * cool1 * cool2
* viscous temp
	tv = t0 * qv ** 0.25
	tv = max( tv, tcmb )
* torque temp
* KDH: 2023-Dec-09
	tq = t0 * ( tv4 * coolq * cool2 ) ** 0.25
	tq = max( tq, tcmb )
c	write(*,*) '** tvtx : t0', t0, ' tv4', ' coolq', coolq, ' tq', tq
* TRAP NAN
      if( tv .ne. tv ) then
	write(*,*) '** WARNING : Tv', tv, ' qv', qv, ' cool1', cool1, ' cool2', cool2
	write(*,*) 'r1', r1, ' tv1', tv1, ' alp', alp
	write(*,*) 'r ', r, ' risco', risco
	stop
	tv = 100.
      end if
* lamp heating
	tx = 0.
	qx = 0.
      if( tx1 .gt. 0. ) then
* distance to lamp
	dz = h - hlamp
	dlamp = sqrt( dz * dz + r * r )
* n = normal to surface (nr,nz) = (-dh,dr) / ds
* p = unit toward lamp  (pr,pz) = (-r,-dz) / dlamp
* cos( incidence angle ) = dot(p,n) = nz*pz + nr*pr
	pdotn = ( r * dh - dz * dr ) / dlamp / ds
* Qx = pdotn (Tx/T0)^4 (R1/dlamp)^2
      if( pdotn .gt. 0. ) then
	t = tx1 / t0
	t = t * t * r1 / dlamp
	qx = t * t * pdotn
	tx = t0 * qx ** 0.25
      end if
      end if

* teff
      if( tx .le. 0. ) then
	t = tv
      else if( tv .le. 0. ) then
	t = tx
      else
	t  = t0 * ( qv + qx ) ** 0.25	
      end if
* stow
	 t_r( i ) = t
	tv_r( i ) = tv
	tx_r( i ) = tx
	tq_r( i ) = tq
      end do
	return
	end subroutine
