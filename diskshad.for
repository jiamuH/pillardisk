! =======================================================
	subroutine diskshad( n, r_r, h_r, hx, i_r )
! find shadows on lamp-irradiated axi-symmetric disk surface
! in:	n	i4 number of radii
!	r_r(n)	r4 radii (light days)
!	h_r(n)	r4 height above midplane (light days)
!	hx	r4 lamp height above origin (light days)
! out:	i_r(n)	i4 1 if irradiated, 0 if shadowed
!
! 2918 Aug (D.Starkey version of KDH algorithm)
! 2021 Sep Keith Horne @ St.Andrews (revised)
! 2021 Nov KDH @ St.A - test for new shadow crest at each radius
! 2021 Nov KDH @ St.A - parabola interpolation of shadow crest

	real*4 r_r( n )
	real*4 h_r( n )
	integer*4 i_r( n )
	logical verbose / .false. /
	logical interpolate / .true. /

      if( n .le. 0 ) then
	write(*,*) '** FATAL ERROR DISKSHAD N =', n
	stop
      end if 

* crest counter
	nc = 0
* shadow threshold (sin( lamp elevation ) = 1 at r = 0)
	smn = 1.
* loop over annuli
      do i = 1, n
* radius and height
	r = r_r( i )
	h = h_r( i )
* e = (er,ez) surface normal 
	im = max( i - 1, 1 )
	ip = min( i + 1, n )
	dr =  r_r( ip ) - r_r( im )
	dh =  h_r( ip ) - h_r( im )
	div = sqrt( dr * dr + dh * dh )
	er = -dh / div
	ez =  dr / div
* u = (ur,uz) direction to lamp
	dr =    - r
	dh = hx - h
	div = sqrt( dr * dr + dh * dh )
	ur = dr / div
	uz = dh / div
* sine of lamp elevation
	s = uz	
* self-shadowed (faces away from the lamp)
	dot = ur * er +  uz * ez
      if( dot .le. 0. ) then
	i_r( i ) = 0
* just past a crest - parabolic interpolation last 3 points
      if( i_r( im ) .eq. 1 ) then
	nc = nc + 1
* report
      if( verbose ) then
	write(*,*) '--- Crest', nc, ' annulus', i, ' ---------'
	imm = max( 1, im - 1 )
	rmm = r_r( imm )
	rm  = r_r( im )
	hmm = h_r( imm )
	hm  = h_r( im )
	write(*,*) '  i ', imm, im, i
	write(*,*) '  R ', rmm, rm, r
	write(*,*) '  H ', hmm, hm, h
	write(*,*) ' irr', i_r( imm ), i_r( im ), i_r( i )
      end if

* interpolate sin( lamp elevation ) with a parabola
      if( interpolate .and. i .ge. 3 ) then
*	s(x) = a + b x + c x^2
*	s   = a + b + c
*	sm  = a
*	smm = a - b + c
	a  = sm
	b  = ( s - smm ) / 2.
	c2 = ( smm + s ) / 2. - sm
* crest at a minimum of s = sin( lamp elevation )
      if( c2 .gt. 0. ) then
	x = - b / c2
	sc = a + ( b + c * x ) * x
* report
      if( verbose ) then
	write(*,*) '  S ', smm, sm, s
	write(*,*) ' abc', a, b, c2/2.
	write(*,*) '  X ', x, ' Sc', sc, ' Smn', smn
      end if
	s = min( s, sc )
      end if
* new crest may lower sin( lamp elevation )
	smn = min( smn, s )
      end if
      end if
* not self-shadowed, but shadowed by previous crest
      else if( s .ge. smn ) then
	i_r( i ) = 0
* exposed to irradiation
      else
	i_r( i ) = 1
      end if
* new crest may lower sin( lamp elevation)
	smn = min( smn, s )
* stow 2 values for parabolic interpolation
      if( interpolate ) then
	smm = sm
	sm = s
      end if
* next radius i
      end do

	if( verbose ) write(*,*) 'Ripple crests', nc

	return
	end subroutine


