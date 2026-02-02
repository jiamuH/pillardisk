* BOWL - fit continuum lags and faint/bright disc SEDs
*        for lamp-post irradiating a bowl-shaped accretion disc
*
* 2021 Aug Keith Horne @ St Andrews - test inclined disc routines
* 2021 Sep KDH @ Crail - rippled surface
* 2021 Oct KDH @ StA - chi2 evaluation
* 2021 Oct KDH @ StA - modules for function bofcalc
* 2021 Oct KDH @ StA - fix bug in tfb (no response from shadows)
* 2021 Oct KDH @ StA - plot models with and without irradiation
* 2021 Oct KDH @ StA - ripple amplitude control - wamp, wpow
* 2021 Nov KDH @ StA - ripple wavelength control wnum, wnpow
* 2021 Nov KDH @ StA - dexflx, chi2brt, chi2fnt, sumlnvar
* 2021 Nov KDH @ StA - maxpar, p_p, dp_p, bof_p
* 2021 Nov KDH @ StA - lagcalc, flxcalc, needflx, needlag
* 2021 Nov KDH @ StA - plot menu
* 2021 Nov KDH @ StA - MCMC hardwired
* 2021 Nov KDH @ Crail - MCMC iteration and corner plots
* 2021 Nov KDH @ Crail - cpu timer, crest
* 2021 Nov KDH @ Crail - mcmcfit -> subroutine
* 2021 Nov KDH @ StA - auto-adjust step sizes
* 2022 Apr KDH @ StA - Tmax
* 2022 Jun KDH @ StA - select files for lag and flux data
* 2022 Jun KDH @ Crail - E(B-V) for MW dust
* 2022 Jun KDH @ Crail - known agn list (5548, 7469)
* 2022 Jun KDH @ StA - rms bandwidths -> fwhm bandwidths
* 2022 Jul KDH @ Crail - tweak plots
* 2022 Jul KDH @ Crail - parameter ranges bot_p top_p
* 2022 Jul KDH @ Crail - logical fit_parameters
* 2022 Jul KDH @ Crail - fitlist, codes, index i_p
* 2022 Jul KDH @ StA - Mdot replaces L/Led as primary
* 2022 Jul KDH @ StA - Rin replaces eta as primary
* 2022 Jul KDH @ StA - units, units_p, symbol, symbol_p
* 2022 Jul KDH @ StA - add M/Msun and D/mpc to primary
* 2022 Jul KDH @ StA - H/R now in percent
* 2022 Aug KDH @ StA - pgmdot
* 2022 Aug KDH @ StA - Lx on T(r) plot
* 2023 May KDH @ StA - add MCG 8-11-11
* 2023-Jun-14 KDH @ StA - needgeom
* 2023-Jun-25 KDH @ StA - fix bug in getpar (cool1*cool2)
* 2023-Jun-27 KDH @ StA - output plot option
* 2023-Jul-13 KDH @ StA - entry of M or logM
* 2023-Aug KDH @ StA - simpleplot option
* 2023-Aug-30 KDH @ StA - implement syslag
* 2023-Aug-30 KDH @ StA - flxdex => dexflx, sysfac => sysflx
* 2023-Sep-06 KDH @ StA - Change units from Rs to Rg
* 2023-Sep-07 KDH @ StA - dexflx => dexsed, sysflx => sigsed
* 2023-Nov-07 KDH @ StA - qin - inner disc torque parameter
* 2023-Nov-15 KDH @ StA - Dkappa = (Hx/rg) (1-A)(Lx/Mdot c^2)
* 2023-Dec KDH @ StA - DlogT = log ( T / Tv )
* 2023-Dec KDH @ StA - fiducial temperatures on T(R) plot
* 2023-Dec KDH @ StA - denser grid near Rin and Rout
* 2023-Dec KDH @ StA - mcmcload
* 2023-Dec KDH @ StA - flxcalc,needflx -> sedcalc,needsed
* 2024-Jan KDH @ StA - 3c273 from data Jiamu Huang
* 2024-Jan KDH @ StA - implement fcol and tc_r
* 2024-Jan KDH @ StA - fnudisk_inclined => fnubowl
* 2024-Jan KDH @ StA -     tfb_inclined =>  tfbowl
* 2024-Feb KDH @ StA - NGC 6814
* 2024-Feb KDH @ StA - NGC 4051
* 2024-Mar-30 KDH @ StA - histogram fiducial temperatures
* 2024-Mar-31 KDH @ StA - histogram Rout and Hout
* 2024-Apr-01 KDH @ StA - rim and no-rim SED and lags
* 2024-Apr-18 KDH @ StA - divide SED wavelengths by (1+z)
* 2024-Apr-24 KDH @ StA - Dmpc=16.6 for NGC 4051
* 2024-Apr-30 KDH @ StA - bhrg was actually bhrs (bug reported by Jiamu Huang)
* 2024-May-09 KDH @ StA - toggle step size adjustments
* 2025-Jan-01 KDH @ Crail - avoid Nr=1 crash in tfbowl
* 2025-Jan-28 KDH @ StA - tidy
* 2025-Mar-06 KDH @ StA - SMC dust on input
* 2025-Mar-06 KDH @ StA - w_b => wobs_b
* 2025-Mar-12 KDH @ StA - keep dlogt in range (0,1)
* 2025-Apr-04 KDH @ StA - line thickness control
* 2025-Apr-28 KDH @ StA - I Zw 1

*------------------------------------------------
	module all
* radii
	parameter( maxr = 10000 )
* wavelengths
	parameter( maxw = 10000 )
* delays
	parameter( maxtau = 10000 )
* bands
	parameter( maxb = 1000 )
* parameters
	parameter( maxpar = 100 )
* iterations
	parameter( maxiter = 100000 )
	parameter( maxplot = maxiter )
* max samples
	parameter( maxi = 1000 )
* disk structure
	real*4 r_r( maxr )
	real*4 h_r( maxr )
	real*4 t_r( maxr )
	real*4 tc_r( maxr )
	real*4 tx_r( maxr )
	real*4 tv_r( maxr )
	real*4 tq_r( maxr )

	real*4 h0_r( maxr )
	real*4 t0_r( maxr )
	real*4 tau_r( maxr )
	real*4 taulo_r( maxr )
	real*4 tauhi_r( maxr )
	integer*4 i_r( maxr )


* flux data (mJy)
	real*4 wdat_b( maxb )
	real*4 fwhm_b( maxb )
	real*4 datflx_b( maxb )
	real*4 sigflx_b( maxb )
	real*4 datflx0_b( maxb )
	real*4 sigflx0_b( maxb )
* lag data (days)
	real*4 wlag0
	real*4 wlag_b( maxb )
	real*4 wlagfwhm_b( maxb )
	real*4 datlag_b( maxb )
	real*4 siglaglo_b( maxb )
	real*4 siglaghi_b( maxb )

* spectra (bright,faint,norim,rim) = (f,f0,f1,f2)
* KDH: w_w = rest-frame wavelengths ?
	real*4 w_w( maxw )
	real*4 f_w( maxw )
	real*4 f0_w( maxw )
	real*4 f1_w( maxw )
	real*4 f2_w( maxw )
* bands  wobs_b = observed-frame
	real*4 wobs_b( maxb )
	integer*4 kolor_b( maxb )

* lags (bright,faint,norim,rim) = (tau,tau0,tau1,tau2)
	real*4 tau_b( maxb )
	real*4 tau0_b( maxb )
	real*4 tau1_b( maxb )
	real*4 tau2_b( maxb )

* delay maps (bright,faint,norim,rim) = (psi,psi0,psi1,psi2)
	real*4  tau_d( maxtau )
	real*4  psi_d( maxtau )
	real*4 psi0_d( maxtau )
	real*4 psi1_d( maxtau )
	real*4 psi2_d( maxtau )

	real*4 psi_db( maxtau, maxb )
	real*4 psi0_db( maxtau, maxb )
	real*4 psi1_db( maxtau, maxb )
	real*4 psi2_db( maxtau, maxb )

* MCMC subset for plots
	integer*4 nsp/ -1 /
	integer*4 nshow/ 50 /
	integer*4 nsamp/ 0 /
	integer*4 nthin/ 1 /
	logical	plot_i( maxi )
	integer*4 iter_i( maxi )

* geometry
	integer*4 nr_i( maxi )
	real*4  r_ri( maxr, maxi )
	real*4  h_ri( maxr, maxi )
	real*4 h0_ri( maxr, maxi )
	real*4 hlamp_i( maxi )
* temperatures
	real*4  t_ri( maxr, maxi )
	real*4 t0_ri( maxr, maxi )
	real*4 tc_ri( maxr, maxi )
	real*4 tx_ri( maxr, maxi )
	real*4 tv_ri( maxr, maxi )
	real*4 tq_ri( maxr, maxi )
* seds
	integer*4 nw_i( maxi )
	real*4  w_wi( maxw, maxi )
	real*4  f_wi( maxw, maxi )
	real*4 f0_wi( maxw, maxi )
	real*4 f1_wi( maxw, maxi )
	real*4 f2_wi( maxw, maxi )
* lags
	integer*4 nb_i( maxi )
	real*4 wobs_bi( maxb, maxi )
	real*4  tau_bi( maxb, maxi )
	real*4 tau0_bi( maxb, maxi )
	real*4 tau1_bi( maxb, maxi )
	real*4 tau2_bi( maxb, maxi )

* real*8 variables
	real*8 sum, top, bot
	real*8 chi2lag, chi2brt, chi2fnt
	real*8 sumlnvarlag, sumlnvarsed
	real*8 bof, chi2, sumlnvar
	real*8 bold, bofold, bofmed

* parameter values, steps, names, codes
	real*4 p_p( maxpar )
	real*4 dp_p( maxpar )
	logical log_p( maxpar )
	logical fit_p( maxpar )
	character*75 units, units_p( maxpar )
	character*75 symbol, symbol_p( maxpar )
	character*75 name, name_p( maxpar )
	character*1 code, code_p( maxpar )
	character*(maxpar) fitlist
* prior
	real*4 bot_p( maxpar )
	real*4 top_p( maxpar )
	real*8 sum_p( maxpar )
* index of fit parameters
	real*4 ifit_p( maxpar )
	real*4 i_p( maxpar )
	real*4 j_p( maxpar )
* simplex (for use by Amoeba)
	parameter( maxvert = 1 + maxpar )
	real*4 p_pv( maxpar, maxvert )
	real*4 bof_pv( maxpar, maxvert )
* mcmc
	integer*4 ibest, iburn1
	real*8 bof_i( maxiter )
	real*8 chi2_i( maxiter )
	real*4 cpu_i( maxiter )

	real*4 avg_p( maxpar )
	real*4 rms_p( maxpar )
	real*4 pmed_p( maxpar )
	real*4 pmad_p( maxpar )

	real*4  p_pi( maxpar, maxiter )
	real*4 dp_pi( maxpar, maxiter )
	integer*4 keep_p( maxpar )
	integer*4 keep_pi( maxpar, maxiter )
	real*4 acc_pi( maxpar, maxiter )

	real*4 old_p( maxpar )
	real*4 save_p( maxpar )
	real*4 best_p( maxpar )

* KDH : ARRAYS FOR SECONDARY PARAMETERS
	real*4 p_si( maxpar, maxiter )
	real*4 avg_s( maxpar )
	real*4 rms_s( maxpar )
	real*4 pmed_s( maxpar )
	real*4 pmad_s( maxpar )

* number of parameters
	integer*4 npar/ 0 /
* number of fit parameters
	integer*4 nfit/ 0 /
* number of secondary parameters
	integer*4 nsec/ 0 /
* number of iterations
	integer*4 nits/ 0 /
* new iterations
	integer*4 newits/ 5 /
* iterations between bof plots and corner plots
	integer*4 itfit/ 0 /
	integer*4 itbof/ 0 /
	integer*4 itcorner/ 0 /
	integer*4 itsave/ 0 /
* adjust step sizes
	logical adjust /.true./
* interval between adjustments
	integer*4 nfix / 30 /

* plot arrays
	real*4 plot0( maxplot )
	real*4 plot1( maxplot )
	real*4 plot2( maxplot )
	real*4 plot3( maxplot )
	real*4 plot4( maxplot )
	real*4 plot5( maxplot )
	real*4 plot6( maxplot )
	real*4 plot7( maxplot )
	real*4 plot8( maxplot )
	real*4 plot9( maxplot )

	real*4 plot10( maxplot )
	real*4 plot11( maxplot )
	real*4 plot12( maxplot )
	real*4 plot13( maxplot )
	real*4 plot14( maxplot )
	real*4 plot15( maxplot )
	real*4 plot16( maxplot )
	real*4 plot17( maxplot )
	real*4 plot18( maxplot )
	real*4 plot19( maxplot )

	real*4 plot20( maxplot )
	real*4 plot21( maxplot )
	real*4 plot22( maxplot )
	real*4 plot23( maxplot )
	real*4 plot24( maxplot )
	real*4 plot25( maxplot )
	real*4 plot26( maxplot )
	real*4 plot27( maxplot )
	real*4 plot28( maxplot )
	real*4 plot29( maxplot )

	real*4 plot30( maxplot )
	real*4 plot31( maxplot )
	real*4 plot32( maxplot )
	real*4 plot33( maxplot )
	real*4 plot34( maxplot )
	real*4 plot35( maxplot )
	real*4 plot36( maxplot )
	real*4 plot37( maxplot )
	real*4 plot38( maxplot )
	real*4 plot39( maxplot )

	integer*4 iplot( maxplot )
* 2d arrays
	parameter( maxbin = 1024 )
	parameter( maxpic = maxbin * maxbin )
	real*4 pic( maxpic )
	integer*4 i2d( maxpic )
* contours
	parameter( maxc = 100 )
	real*4 clev( maxc )

	logical simpleplot/ .false. /
	logical bestfitplot/ .false. /
	logical medianplot/ .false. /

* plot window
	integer*4 logx, logy
	real*4 xmn, xmx, ymn, ymx
* plot labels
	character*100 reply, label, xlabel, ylabel, title
	character*75 plotoption / '6' /
* files
	character*75 file, file_lag, file_sed

* dimensions
* radii
c	integer*4 nr
* bands (for lag spectra)
	integer*4 nb
* wavelengths (for flux spectra)
	integer*4 nw
* delays
	integer*4 ntau
* lags (for lag data)
	integer*4 nlag
* flux data
	integer*4 ndat
* total data
	integer*4 ntot

* reference wavelength
	real*4 wobs, wrest
* lag with lamp on and off
	real*4 tauavg, tauavg0

	save

	end module all
*----------------------------------------------
	module pgplot
* pgplot colour indices
	integer kblack / 1 /
	integer	kred / 2 /
	integer kgreen / 3 /
	integer kblue / 4 /
	integer kcyan / 5 /
	integer kmagenta / 6 /
	integer kyellow / 7 /
	integer korange / 8 /
	integer kgreenyellow / 9 /
	integer kgreencyan / 10 /
	integer kbluecyan / 11 /
	integer kbluemagenta / 12 /
	integer kredmagenta / 13 /
	integer kdarkgrey / 14 /
	integer kgrey / 15 /
* pgplot line styles
	integer ldash / 2 /
	integer ldashdot / 3 /
	integer ldot / 4 /
* pgplot symbols
	integer idot / 1 /
	integer iopencircle / 4 /
	integer ifilledcircle / 17 /
	integer ismallcircle / 20 /
	integer istar / 18 /
	real*4 tr( 6 )
	save
	end module pgplot
*----------------------------------------------
	module constants
	real*4 pi, dtr
	real*4 clight, gnewt, stefan, sigmat, emp
	real*4 pc, year, day, vkms
	real*4 sunm, sunl, sunrs, sunrg, sunledd
	save
	end module constants
*----------------------------------------------
	module cosmology
* cosmology parameters
	real*4	h0 / 70. /
	real*4	omat / 0.3 /
	real*4	olam / 0.7 /

c	real*4	h0 / 70. /
c	real*4	omat / 0.28 /
c	real*4	olam / 0.72 /
* Planck collaboration 2015
c	real*4	h0 / 67.8. /
c	real*4	omat / 0.31 /
c	real*4	olam / 0.69 /
	end module cosmology
*----------------------------------------------
* agn parameters (with values for primaries)
	module agn
* 2025-01-27 Keith Horne @ St Andrews -- logical primary_xlamp
* 2025-04-14 KDH @ StA -- logical primary_dh

* known agn
	parameter (maxa = 20 )
	integer nagn, iagn
	real*4 redshift_a( maxa )
	real*4 ebmv_a( maxa )
	real*4 embh_a( maxa )
	character*75 name_a( maxa )
	character*75 file_lag_a( maxa )
	character*75 file_sed_a( maxa )

* agn-specific parameters
* redshift (default values for NGC 5548)
	real*4 redshift / 0.0172 /
	real*4 dmpc, dmpcz
* M_BH (solar)
	real*4 bhM_Msun / 7.e7 /
	real*4 bhRg
* L/Ledd
	real*4 bhL_Ledd/ 0.0022 /
* L = eta Mdot c^2
	real*4 eta/ 0.0833 /

* dust extinction (applied on input)
	real*4 ebmv_mw/ 0. /
	real*4 ebmv_smc/ 0. /
* more SMC-like dust
	real*4 ebmv_agn/ 0. /

	real*4 bhL_Lsun
	real*4 bhLedd_Lsun
	real*4 bhMsun_yr/ 0.0014 /

* disc geometry parameters
* number of radii
	integer*4 nr / 500 /
* R1 reference radius (light days)  (use rout if r1 = 0 )
	real*4 r1 / 0. /
	real*4 r1_rg
	real*4 r0
* Rout outer radius (light days)
	real*4 rout / 6. /
	real*4 rout_rg
* Risco / Rg (6 for zero spin, 3 for max prograde spin.
	real*4 risco_rg / 6. /
	real*4 risco

* power-law disk parameters
* h(r) = h1 (r/r1)**bet (light days)
	real*4 h1
	real*4 h1_rg
* H/R
	real*4 h1_r1_pct / 1. /
	real*4 bet / 1.4 /
* primary dh = h1-hx to de-correlate Hx and H/R
	logical primary_dh/ .false. /
	real*4 dh_r1_pct/ 1. /

* wave amplitude and power-law index
	real*4 wamp/ 0. /
	real*4 wpow/ 0. /
* wave number and power-law index
	real*4 wnum/ 4. /
	real*4 wnpow/ 0. /
* wave crest shape
	real*4 crest/ 1. /

* power-law temperature profile
* T(r) = T1 (r/r1)**(-alp) (K)
	real*4 t1 / 4.e4 /
	real*4 alp / 0.75 /

	real*4 tv1, tx1, tvmx, tvmn

* inner disc torque parameter and temperature
	real*4 qin / 0. /
	real*4 tqin, tqmx

* Hlamp / Rg
	real*4 hlamp_rg / 6. /
	real*4 hlamp
* albedo
*	real*4 albedo / 0.5 /
	real*4 albedo / 0.0 /
* colour temperature boost factor
	real*4 fcol / 1.0 /

* primary xlamp vs dlogT
	logical primary_xlamp/ .false. /
* lamp boost xlamp = Lx/( Mdot c^2)
	real*4 xlamp / 5. /
* dlogt = log( T / Tv )
	real*4 dlogt/ 0.1 /
* dkappa = (Hx/Rg) (Lx/Mdot c^2) (1-A)
	real*4 dkappa

* inclination (deg)
	real*4 dinc / 45. /
	real*4 dincmn/ 0. /
	real*4 dincmx/ 60. /
	real*4 cosi, sini

* flux model uncertainties (dex)
	real*4 dexsed / 0.005 /
	real*4 sigsed
* lag model uncertainties (day)
	real*4 syslag/ 1e-5 /


	logical needsed / .true. /
	logical needlag / .true. /
	logical needgeom / .true. /
	save
	end module agn

* ===========================================
* modules:
	use all
	use constants
	use cosmology
	use agn
	use pgplot
	logical pgagain

* toggle verbose diagnostics
	logical verbose/ .false. /
	logical yes

	r0 = r1
	if( r0 .le. 0. ) r0 = rout
      if( primary_dh ) then
	h1_r1_pct = hlamp / r0 * 100. + dh_r1_pct
      else
	dh_r1_pct = h1_r1_pct - hlamp / r0 * 100.
      end if


* preamble
	write(*,*)
	write(*,*) '   ************************************************'
	write(*,*) '   *  BOWL - fits lags and faint/bright disc SEDs *'
	write(*,*) '   * with a lamp post irradiated accretion disc   *'
	write(*,*) '   * with a flared and/or rippled H(r) thickness. *'
	write(*,*) '   *                                              *'
	write(*,*) '   *  Author: Keith Horne, SUPA St Andrews        *'
	write(*,*) '   *  e-mail: kdh1@st-andrews.ac.uk               *'
	write(*,*) '   *  2021-Aug - First version                    *'
	write(*,*) '   *  2021-Nov - MCMC fits to NGC 5548            *'
	write(*,*) '   *  2023-Aug - Simple plots for F9 Swift paper  *'
	write(*,*) '   *  2023-Nov - Plots for last,best,median parms *'
	write(*,*) '   *  2023-Dec - Save and Load MCMC chain         *'
	write(*,*) '   ************************************************'
	write(*,*)
	write(*,*) 'redshift', redshift
	write(*,*) 'Mbh', bhM_Msun, ' Mdot', bhMsun_yr
	write(*,*) 'Rin/Rg', risco_rg, ' Rout', rout, ' ltd'
	write(*,*) '     H/R', h1_r1_pct, ' %', ' R1', r1
	write(*,*) '(H-Hx)/R', dh_r1_pct, ' % ', primary_dh
	write(*,*) 'beta', bet, ' dH/dR', bet * h1_r1_pct, ' %'
	write(*,*) 'Hx/Rg', hlamp_rg, ' Xlamp', xlamp, ' DlogT', dlogt
	write(*,*) 'albedo', albedo, ' fcol', fcol, ' primary_xlamp', primary_xlamp
* initialise constants
	call setconstants

	iseed = 729582

* known agn
	n = 0
	n = n + 1
	pg0844 = n
	name_a( n ) = 'PG 0844+349'
	redshift_a( n ) = 0.064
	embh_a( n ) = 4e7
	ebmv_a( n ) = 0.032
	file_lag_a( n ) = 'pg0844_lag.dat'
	file_sed_a( n ) = 'pg0844_mjy.dat'

	n = n + 1
	n5548 = n
	name_a( n ) = 'NGC 5548'
	redshift_a( n ) = 0.0172
	embh_a( n ) = 7e7
	ebmv_a( n ) = 0.0
	file_lag_a( n ) = 'ngc5548_lag.dat'
	file_sed_a( n ) = 'ngc5548_mjy.dat'

	n = n + 1
	n7469 = n
	name_a( n ) = 'NGC 7469'
	redshift_a( n ) = 0.01656
	embh_a( n ) = 9e6
	ebmv_a( n ) = 0.069
	file_lag_a( n ) = 'ngc7469_lag.dat'
	file_sed_a( n ) = 'ngc7469_mjy.dat'

	n = n + 1
	j1119 = n
	name_a( n ) = 'PG 1119+120'
	redshift_a( n ) = 0.050201
	embh_a( n ) = 10. ** 7.73
	embh_a( n ) = 1e7
	ebmv_a( n ) = 0.033
	file_lag_a( n ) = 'pg1119_lag.dat'
	file_sed_a( n ) = 'pg1119_mjy.dat'


c	n = n + 1
c	jf9 = n
c	name_a( n ) = 'Fairall 9'
c	redshift_a( n ) = 0.046145
c	embh_a( n ) = 10. ** 8.3
c	ebmv_a( n ) = 0.026
c	file_lag_a( n ) = 'f9_lag.dat'
c	file_sed_a( n ) = 'f9_mjy.dat'


	n = n + 1
	n4395 = n
	name_a( n ) = 'NGC 4395'
	redshift_a( n ) = 0.001064
	embh_a( n ) = 10. ** 5.44
	ebmv_a( n ) = 0.017
	file_lag_a( n ) = 'ngc4395_lag.dat'
	file_sed_a( n ) = 'ngc4395_mjy.dat'

	n = n + 1
	mrk876 = n
	name_a( n ) = 'Mrk 876'
	redshift_a( n ) = 0.12109
	embh_a( n ) = 10. ** 8.34
	ebmv_a( n ) = 0.027
	file_lag_a( n ) = 'mrk876_lag.dat'
	file_sed_a( n ) = 'mrk876_mjy.dat'

	n = n + 1
	mrk817 = n
	name_a( n ) = 'Mrk 817'
	redshift_a( n ) = 0.031328
	embh_a( n ) = 10. ** 7.6
	ebmv_a( n ) = 0.007
	ebmv_a( n ) = 0.0
	file_lag_a( n ) = 'Mrk817_lag.dat'
	file_sed_a( n ) = 'Mrk817_mjy.dat'

	n = n + 1
	if9 = n
	name_a( n ) = 'Fairall 9'
	redshift_a( n ) = 0.046145
	embh_a( n ) = 10. ** 8.3
	ebmv_a( n ) = 0.026
	ebmv_a( n ) = 0.
	file_lag_a( n ) = 'f9_lag.dat'
	file_sed_a( n ) = 'f9_mjy.dat'

	n = n + 1
	mcg = n
	name_a( n ) = 'MCG 8-11-11'
	redshift_a( n ) = 0.020457
	embh_a( n ) = 10. ** 7.19
	ebmv_a( n ) = 0.0
	file_lag_a( n ) = 'mcg_lag.dat'
	file_sed_a( n ) = 'mcg_mjy.dat'

	n = n + 1
	i273 = n
	name_a( n ) = '3C273'
	redshift_a( n ) = 0.158339
	embh_a( n ) = 10. ** 8.84
	ebmv_a( n ) = 0.0
	file_lag_a( n ) = '3c273_lag.dat'
	file_sed_a( n ) = '3c273_mjy.dat'

	n = n + 1
	n6814 = n
	name_a( n ) = 'NGC 6814'
	redshift_a( n ) = 0.00522
	embh_a( n ) = 10. ** 6.42
	ebmv_a( n ) = 0.185
	ebmv_a( n ) = 0.0
	file_lag_a( n ) = 'ngc6814_lag.dat'
	file_sed_a( n ) = 'ngc6814_mjy.dat'

	n = n + 1
	n6814 = n
	name_a( n ) = 'NGC 4051'
	redshift_a( n ) = 0.002366
	embh_a( n ) = 10. ** 5.89
	ebmv_a( n ) = 0.013
	ebmv_a( n ) = 0.0
	file_lag_a( n ) = 'ngc4051_lag.dat'
	file_sed_a( n ) = 'ngc4051_mjy.dat'


	n = n + 1
	n3783 = n
	name_a( n ) = 'NGC 3783'
	redshift_a( n ) = 0.00973
	embh_a( n ) = 10. ** 7.37
	ebmv_a( n ) = 0.121
	ebmv_a( n ) = 0.127
	file_lag_a( n ) = 'ngc3783_yr1_lag.dat'
	file_sed_a( n ) = 'ngc3783_yr1_mjy.dat'

	n = n + 1
	n4593 = n
	name_a( n ) = 'NGC 4593'
	redshift_a( n ) = 0.008312
	embh_a( n ) = 10. ** 6.9
	ebmv_a( n ) = 0.021
	file_lag_a( n ) = 'ngc4593_lag.dat'
	file_sed_a( n ) = 'ngc4593_mjy.dat'

	n = n + 1
	n1302 = n
	name_a( n ) = 'PG 1302-102'
	redshift_a( n ) = 0.2784
	embh_a( n ) = 5e8
	ebmv_a( n ) = 0.037
	file_lag_a( n ) = 'pg1302_lag.dat'
	file_sed_a( n ) = 'pg1302_mjy.dat'

	n = n + 1
	izw1 = n
	name_a( n ) = 'I Zw 1'
	redshift_a( n ) = 0.061169
	embh_a( n ) = 1e7
	ebmv_a( n ) = 0.057
	ebmv_a( n ) = 0.00
	file_lag_a( n ) = 'izw1_lag.dat'
	file_sed_a( n ) = 'izw1_mjy.dat'

	nagn = n

      if( nagn .gt. maxa ) then
	write(*,*) '** AGN list overflow', nagn, ' max', maxa
	stop
      end if

* default agn (last one holds)
	i = j1119
	i = n7469
	i = n4395
	i = mcg
	i = if9
	i = n6814
	i = n5548

	redshift = redshift_a( i )
	bhM_Msun = embh_a( i )
	ebmv_mw = ebmv_a( i )
	file_lag = file_lag_a( i )
	file_sed = file_sed_a( i )
	iagn = i

      if( iagn .eq. if9 ) then
	bhM_Msun = 2.6e8
	bhMsun_yr = 0.03
	bet = 100
	rout = 11
	qin = 1.
	risco_rg = 1.5
      end if

* 5 bands for the model
	nb = 5
	wobs_b( 1 ) = 1000.
	wobs_b( 2 ) = 1500.
	wobs_b( 3 ) = 3000.
	wobs_b( 4 ) = 5000.
	wobs_b( 5 ) = 9000.
	kolor_b( 1 ) = kblack
	kolor_b( 2 ) = kmagenta
	kolor_b( 3 ) = kblue
	kolor_b( 4 ) = kgreen
	kolor_b( 5 ) = kred
* 10 bands for the model
	nb = 10
	wobs_b( 1 ) = 1000.
	wobs_b( 2 ) = 1200.
	wobs_b( 3 ) = 1500.
	wobs_b( 4 ) = 2000.
	wobs_b( 5 ) = 3000.
	wobs_b( 6 ) = 4000.
	wobs_b( 7 ) = 5000.
	wobs_b( 8 ) = 6000.
	wobs_b( 9 ) = 8000.
	wobs_b( 10 ) = 10000.
	kolor_b( 1 ) = kblack
	kolor_b( 2 ) = kblack
	kolor_b( 3 ) = kmagenta
	kolor_b( 4 ) = kmagenta
	kolor_b( 5 ) = kblue
	kolor_b( 6 ) = kblue
	kolor_b( 7 ) = kgreen
	kolor_b( 8 ) = kgreen
	kolor_b( 9 ) = kred
	kolor_b( 10 ) = kred


* number of radii
	nr = 500
* smaller for speed
	nr = 50

* parameters to fit
	fitlist = 'AXIOHB*FL'
	nfit = len1(fitlist)

* load faint and bright disc SED data
    1	write(*,*)
	write(*,*) 'Input faint/bright SED file = 6-column ascii file:'
	write(*,*) ' Col 1 = pivot wavelength (A) observed frame'
	write(*,*) ' Col 2 =  fwhm  bandwidth (A) observed frame'
	write(*,*) ' Col 3 = faint  disc flux (mJy)'
	write(*,*) ' Col 4 = 1 sigma error bar on Col 3'
	write(*,*) ' Col 5 = bright disc flux (mJy)'
	write(*,*) ' Col 6 = 1 sigma error bar on Col 5'
	write(*,*)


	write(*,*)
	write(*,*) 'AGN known to Bowl:'
* internal options
      do i = 1, nagn
	name = name_a( i )
	call nospace( name )
	file = file_sed_a( i )
	z = redshift_a( i )
	e = ebmv_a( i )
	b = alog10( embh_a( i ) )
	title = ' '
	call append_i( title, i )
	call append( title, name(:len1(name)) )
	call append( title, 'file:' )
	call append( title, file(:len1(file)) )
c	write(*,*) i, ' ', name(:len1(name)), ' file: ', file(:len1(file))
	write(*,*) title(:len1(title) )
	write(*,*) '      z', z, ' E(B-V)', e, ' logM', b
      end do

	file = file_sed
	write(*,*) 'Pick an AGN (by number, name, or file name from above list)'
	write(*,*) ' or enter the name of the disc SED file.'
* external options
	write(*,*)
	write(*,*) 'Candidate SED input files:'
	call spawn( 'ls -l *mjy.dat' )
	write(*,*)
	call inq( 'Disc SED file', file, file )

* known agn?
c	iagn = 0
	match = 0
      do i = 1, nagn

	title = file
	call upcase( title, title )
	call nospace( title )
	if( index( title, '.DAT' ) .gt. 0 ) exit

* match the AGN name
      if( len1( title ) .ge. 3 ) then
	name = name_a( i )
	call upcase( name, name )
	call nospace( name )
	m = index( name(:len1(name)) , title(:len1(title)) )
	write(*,*) 'index( ', name(:len1(name)), ' , ', title(:len1(title)), ' ) = ', m
c	yes = index( name, title ) .gt. 0
	yes = m .gt. 0
	write(*,*) 'Does ', title(:len1(title)), ' match ', name(:len1(name)), ' ? ', yes
* match the sequence number
      else
	write( title, * ) i
	call word1( title, l1, l2 )
	call word1(  file, i1, i2 )
	yes = l2 - l1 .eq. i2 - i1 .and. title(l1:l2) .eq. file(i1:i2)
	write(*,*) 'Does ', file(i1:i2), ' match ', title(l1:l2), ' ? ', yes
      end if

* match achieved
      if( yes ) then
	match = i
	iagn = i
	name = name_a( i )
	write(*,*) 'MATCH with ', name(:len1(name))
	redshift = redshift_a( i )
	bhM_Msun = embh_a( i )
	ebmv_mw = ebmv_a( i )
	file = file_sed_a( i )
	file_sed = file_sed_a( i )
	file_lag = file_lag_a( i )
	write(*,*) 'AGN parameters set for ', name(:len1(name))
	write(*,*) '     z ', redshift
	write(*,*) ' E(B-V)', ebmv_mw, ' MW'
	write(*,*) ' M/Msun', bhM_Msun, ' log(M/Msun)', alog10( bhM_Msun )
* exit the loop
	exit
      end if
* next AGN i
      end do

* load SED data
	n = maxb
	call load1d6( file, 6, n, wdat_b, fwhm_b,
     &		datflx0_b, sigflx0_b, datflx_b, sigflx_b )
	write(*,*) 'Disc SED data', n
      if( n .le. 0 ) then
	write(*,*) 'No SED data.'
	goto 20
      end if

	ndat = n
      do i = 1, ndat
	write(*,*) nint( wdat_b(i) ), nint( fwhm_b(i) ),
     &		datflx0_b(i), datflx_b(i)
      end do
	file_sed = file

* dust extinction correction
	write(*,*) 'Enter E(B-V) for MW and host (SMC-like) dust corrections.'
	val1 = ebmv_mw
	val2 = ebmv_smc
	call inq2r( 'E(B-V) for MW and SMC dust', val1, val2 )
	ebmv_mw = val1
	ebmv_smc = val2

	label = '      lam(A)     fwhm(A)      factor'
	label = label(:len1(label)) // '    faint(mJy)    bright(mJy)'
* apply MW dust extinction (observed frame)
	ebmv = ebmv_mw
      if( ebmv .eq. 0. ) then
	write(*,*) 'No correction here for observed-frame MW dust extinction.'
	ebmv_mw = 0.
      else
	ebmv_mw = ebmv
	write(*,*) 'Apply observed-frame MW dust extinction for E(B-V)', ebmv
	write(*,*) label(:len1(label))
      do i = 1, ndat
	w = wdat_b( i )
	fac = 10. ** ( 0.4 * extmag_mw( w, ebmv ) )
	datflx0_b( i ) = datflx0_b( i ) * fac
	 datflx_b( i ) = datflx_b( i ) * fac
	write(*,*) nint( w ), nint( fwhm_b(i) ), fac,
     &		datflx0_b(i), datflx_b(i)
      end do
	write(*,*) label(:len1(label))
      end if

* apply SMC-like dust extinction (rest frame)
	ebmv = ebmv_smc
      if( ebmv .eq. 0. ) then
	write(*,*) 'No correction here for rest-frame SMC-like dust extinction.'
	ebmv_smc = 0.
      else
	ebmv_smc = ebmv
	write(*,*) 'Apply rest-frame SMC-like dust extinction for E(B-V)', ebmv
	opz = 1. + redshift
	write(*,*) 'redshift', redshift, ' 1+z', opz
	write(*,*) label(:len1(label))
      do i = 1, ndat
	w = wdat_b( i ) / opz
	fac = 10. ** ( 0.4 * extmag_smc( w, ebmv ) )
	datflx0_b( i ) = datflx0_b( i ) * fac
	 datflx_b( i ) = datflx_b( i ) * fac
	write(*,*) nint( w ), nint( fwhm_b(i) / opz ), fac,
     &		datflx0_b(i), datflx_b(i)
      end do
	write(*,*) label(:len1(label))
      end if

* adjust eta
      if( iagn .eq. n7469 ) then
	eta = 0.1
      else if( i .eq. n4051 ) then
	dmpc = 16.6
	write(*,*) 'NGC 4051 distance set to', dmpc, ' Mpc'
      end if

* what ends here?
c      end if

	if( reply(1:2) .eq. 'DS' ) goto 4

* load lag data
    2	write(*,*)
	write(*,*) 'Disk LAG file = 5-column ascii file:'
	write(*,*) ' Col 1 = pivot wavelength (A)'
	write(*,*) ' Col 2 =  fwhm  bandwidth (A)'
	write(*,*) ' Col 3 = rest-frame lag (days)'
	write(*,*) ' Col 4 = -1 sigma error bar on Col 3'
	write(*,*) ' Col 5 = +1 sigma error bar on Col 3'
	write(*,*)
	i1 = 1
	i2 = nagn
      if( match .ge. 1 .and. match .le. nagn ) then
	i1 = match
	i2 = i1
      end if
      do i = i1, i2
	name = name_a( i )
	file = file_lag_a( i )
	write(*,*) i, ' ', file(:len1(file))
      end do

	write(*,*)
	write(*,*) 'Candidate LAG input files:'
	call spawn( 'ls -l *lag.dat' )
	write(*,*)
	file = file_lag
	write(*,*) 'Enter file name or AGN number from above list.'
	call inq( 'Disc LAG file', file, file )
	n = maxb
	call load1d5( file, 5, n, wlag_b, wlagfwhm_b,
     &		datlag_b, siglaglo_b, siglaghi_b )
      if( n .le. 0 ) then
	write(*,*) 'No lag data.'
	goto 20
      else
	write(*,*) 'Disc LAG data', n
	nlag = n
	i1 = 1
	i0 = 0
      do i = 1, nlag
	write(*,*) nint( wlag_b(i) ), nint( wlagfwhm_b(i) ),
     &		datlag_b(i), siglaglo_b(i), siglaghi_b(i)
* keep lower error bar positive
	s = siglaglo_b( i )
	if( s .lt. 0. ) siglaglo_b(i) = abs( s )
* smallest lag
	if( abs( datlag_b( i ) ) .lt. abs( datlag_b( i1 ) ) ) i1 = i
* negative lam
	if( wlag_b(i) .lt. 0. ) i0 = i
      end do
* smallest lag
	write(*,*) 'Smallest lag', datlag_b(i1), ' at', nint( wlag_b(i1) )
	wlag0 = wlag_b( i1 )
* negative wavelength
      if( i0 .gt. 0 ) then
	write(*,*) 'Negative lam', datlag_b(i0), ' at', nint( wlag_b(i0) )
	wlag0 = wlag_b( i0 )
      end if

	file_lag = file

* reference wavelength for CCF lag = 0
* the wavelengths are observed frame. The lags have been corrected.
* KDH : NEED TO SET THIS IN THE INPUT FILE IN SOME WAY  (NEGATIVE LAMBDA?)
c	wlag0 = 1367

* KDH : WHY NOT DO 3 HERE AS WELL?
C	if( reply(1:2) .eq. 'DL' .or. reply(1:2) .eq. 'DS' ) goto 4
c	goto 4
      end if


* initialise MCMC framework
    3	nfit = 0
	npar = 0
	call initpar
	write(*,*) 'Npar', npar, ' Nfit', nfit

* branch to here if parameters changed
    4	needsed = .true.
	needlag = .true.

* KDH : SET PARAMETERS IN P_P TO AVOID REVERSION IN CALL TO BOFCALC
	call setpar

* KDH : WHY SET REPLY TO 0 ? -- AVOIDS INFINITE LOOP
	reply = '0'

* branch to here if parameters unchanged
    5	call upcase( reply, reply )
* compute secondary parameters
	call secondary
* set parameters in p_p
c	call setpar

* update model
* KDH : THESE NOT INVOKED SINCE REPLY = '0' ON BRANCH TO 4
      if( reply(1:1) .eq. 'U'
     &	.or. reply(1:1) .eq. 'P'
     &	.or. reply(1:1) .eq. 'N'
     &	.or. needlag
     &	.or. needsed
     &	.or. needgeom
     &	 ) then

	bof = bofcalc( p_p )
* plot
	if( reply(1:1) .eq. 'P' ) goto 90
      end if

* command menu  ----------------------------------
* report parameters
   10	write(*,*)
	write(*,*) ' BOWL parameters:'
	name = name_a( iagn )
	write(*,*) 'AGN ', name(:len1(name))
	write(*,*) '-- Cosmology: --------------------'
	write(*,*) '       H0', h0, ' km/s/Mpc'
	write(*,*) '  Omega_M', omat, ' Omega_L', olam
	write(*,*) '-- AGN: --------------------------'
	write(*,*) 'Z       z', redshift, ' D_L(z)', dmpcz, ' Mpc'
	write(*,*) 'D     D_L', dmpc, ' Mpc'
	write(*,*) 'M     Mbh', bhM_Msun, ' Msun', ' log(M/Msun)', alog10( bhM_Msun )
	write(*,*) '       Rg', bhRg, ' ltd'
	write(*,*) '     Ledd', bhLedd_Lsun, ' Lsun'

	write(*,*) '-- Disk Geometry: ----------------'
	write(*,*) 'I       i', dinc, ' deg', ' cosi', cosi
	write(*,*) 'IM mnmx(i)=(', dincmn, dincmx, ') deg'
	write(*,*) 'N    Nrad', nr, ' max', maxr
	write(*,*) 'S   Risco', risco_rg, ' Rg =', risco, ' ltd'
	write(*,*) 'O    Rout', rout_rg, ' Rg =', rout, ' ltd'
      if( r1 .ne. 0. ) then
	write(*,*) '1      R1', r1_rg, ' Rg =', r1, ' ltd'
      else
	write(*,*) '1      R1 => Rout'
      end if
	h1_r1 = h1_r1_pct / 100.
	write(*,*) 'H      H/R', h1_r1_pct, '% at R1'
	write(*,*) 'H1   H(R1)', h1/bhrg, ' Rg =', h1, ' ltd'
	hout = h1
	if( r1 .ne. 0. ) hout = hout * ( rout / r1 ) ** bet
	write(*,*) '     H/R', hout / rout * 100., ' % at Rout'
	write(*,*) '    (H-Hx)/R', dh_r1_pct, ' %', primary_dh
	write(*,*) 'B    beta', bet, ' = dlnH/dlnR  dH/dR', bet * h1_r1_pct, ' %'

	write(*,*) '-- Ripples: ----------------------'
	write(*,*) 'W    Wamp', wamp, ' Wpow', wpow
	write(*,*) 'K    K(R1)', wnum, ' dlnK/dlnR', wnpow
	write(*,*) 'C    Crest', crest, ' crest peakiness'

	write(*,*) '-- Accretion: --------------------'
	write(*,*) 'MD   Mdot', bhMsun_yr, ' Msun/yr'
	write(*,*) '   M Mdot', bhMsun_yr*bhM_Msun, ' Msun^2/yr'
	write(*,*) 'L   Ldisc', bhL_Lsun, ' Lsun'
	write(*,*) 'ED  Ldisc', bhL_Ledd, ' Ledd'
	write(*,*) 'E     eta', eta, ' = L/(Mdot c^2)'
     &		, ' = Rg/2Risco=', 0.5/risco_rg
	write(*,*) 'A   alpha', alp, ' = dlnT/dlnR'
	write(*,*) '   Tv(R1)', nint(tv1), ' K', '  Tv(max)', nint(tvmx), ' K'
	write(*,*) '   Tv(min)', nint(tvmn), ' K'
	write(*,*) '^    fcol', fcol, ' colour temp boost Tc/Teff'
	write(*,*) '@     Qin', qin, ' inner disc torque parameter'
	write(*,*) '   Tq(Rin)', nint(tqin), ' K', '  Tq(max)', nint(tqmx), ' K'
	write(*,*) '-- Lamp: -------------------------'
	write(*,*) '*   Hlamp', hlamp_rg, ' Rg =', hlamp, ' ltd'
	write(*,*) 'X    lamp', xlamp, ' (Lx/Mdot c^2) (2 isotropic lamps)'
	write(*,*) 'Y    Flat', ylamp, ' (2 flat lamps)'
	write(*,*) 'DT  DlogT', dlogt, ' (dex)'
	write(*,*) 'AO albedo', albedo
	write(*,*) 'DK	DK =', dkappa, ' (Hx/Rg)(1-A)(Lx/Mdot c^2)'
	write(*,*) '   Tx(R1)', tx1, ' K'

	write(*,*) '-- Lag and SED data -------------'
	write(*,*) 'DL   disc Lag data ', file_lag(:len1(file_lag))
	write(*,*) 'DS   disc SED data ', file_sed(:len1(file_sed))
	write(*,*) 'DUst E(B-V) MW', ebmv_mw, ' SMC', ebmv_smc, ' AGN', ebmv_agn

	write(*,*) '-- Fit Statistics: ---------------'
c	call setpar
	write(*,*) 'EL  lag sys err', syslag, ' d'
	write(*,*) 'EF flux sys err', dexsed, ' dex =', 100.*sigsed, ' %'

* lag data
	write(*,*) nlag, ' Lag data'
	n = max( 1, nlag )
	x = chi2lag / n
	write(*,*) '   Lag Chi2', real( chi2lag ), ' /', n, ' =', x
	x = sumlnvarlag / n
	write(*,*) '   SumLnVar', real( sumlnvarlag ), ' /', n, ' =', x
* flux data
	write(*,*) ndat, ' SED data'
	n = max( 1, ndat )
* bright
	x = chi2brt / n
	write(*,*) 'Bright Chi2', real( chi2brt ), ' /', n, ' =', x
* faint
	x = chi2fnt / n
	write(*,*) ' Faint Chi2', real( chi2fnt ), ' /', n, ' =', x
	n = max( 1, 2 * ndat )
	x = sumlnvarsed / n
	write(*,*) '   SumLnVar', real( sumlnvarsed ), ' /', n, ' =', x

* total data and fit parameters
	write(*,*) ntot, ' Total data',  nfit, ' fit parameters'
* degrees of freedom
	ndof = ntot - nfit
	write(*,*) ndof, ' Degrees of Freedom'
	n = max( 1, ndof )
	x = chi2 / n
	write(*,*) ' Total Chi2', real( chi2 ), ' /', n, ' =', x
	x = sumlnvar / n
	write(*,*) '   SumLnVar', real( sumlnvar ), ' /', n, ' =', x
	x = bof / n
	write(*,*) 'BoF = Chi2 + SumLnVar =', bof, ' /', n, ' =', x

* warning  (JUST DO IT)
      if( needlag .or. needsed ) then
	write(*,*) '*** MODEL NEEDS UPDATED ***'
c	reply = 'U'
c	goto 25
      else
	write(*,*) '*** MODEL IS UP TO DATE ***'
      end if
* actions
	write(*,*)
	write(*,*) '-- Actions: ----------------------'
	write(*,*) 'U ... Update model'
	write(*,*) 'P ... Plot ( Panels: F L T H D W 2 4 6 )'
      if( nits .gt. 0 ) then
	write(*,*) ' PB .. Best parameters  ', bestfitplot
	write(*,*) ' PM .. Median parameters', medianplot
	write(*,*) ' PS .. Simple plot      ', simpleplot
	write(*,*) ' P- .. reset options (BMS) => ', .false.
	write(*,*) ' PO .. Plot with Output to files'
      end if
	write(*,*) 'F ... Fit iterations', nits, ' new', newits
      if( nits .gt. 0 ) then
	write(*,*) '      CPU:', cpumcmc, ' CPU/it:', cpumcmc / max(1,nits)
	write(*,*) '      MCMC(', 1,     ' ) BoF', bof_i( 1 )
	write(*,*) '      MCMC(', nits,  ' ) BoF', bof_i( nits )
	write(*,*) '      MCMC(', ibest, ' ) BoF', bof_i( ibest ), ' best'
      end if

	write(*,*) 'Q ... Quit'

* command node  ----------------------------------

   20	if( reply .ne. '?' ) reply = ' '
	call inq( 'BOWL: What next to do or change >', reply, reply )
	call upcase( reply, reply )

	write(*,*)

      if( reply(1:1) .eq. '?' ) then
	reply = ' '
	goto 5
      end if

* Note: 3 and 2-letter commands must precede corresponding 1-letter commands
      if( reply(1:3) .eq. 'AGN' ) then
	name = ' '
	if( iagn .gt. 0 ) name = name_a( iagn )
	call inq( 'AGN name :', name, reply )
	name = reply
* KDH: SHOULD WE ADD A NEW AGN TO THE LIST HERE?
	name_a( iagn ) = name
	goto 10
      end if

	if( reply(1:2) .eq. 'DS' ) goto 1
	if( reply(1:2) .eq. 'DL' ) goto 2

	if( reply(1:2) .eq. 'MD' ) goto 43
	if( reply(1:3) .eq. 'ERR' ) goto 39
	if( reply(1:3) .eq. 'EF' ) goto 39
	if( reply(1:3) .eq. 'EL' ) goto 391
	if( reply(1:3) .eq. 'ET' ) goto 391
	if( reply(1:2) .eq. 'ED' ) goto 40
	if( reply(1:1) .eq. 'H1' ) goto 30
	if( reply(1:1) .eq. 'H' ) goto 31
	if( reply(1:1) .eq. 'B' ) goto 32
	if( reply(1:1) .eq. '^' ) goto 33
	if( reply(1:2) .eq. 'IM' ) goto 35
	if( reply(1:1) .eq. 'I' ) goto 34
	if( reply(1:1) .eq. 'Z' ) goto 36
	if( reply(1:2) .eq. 'DT' ) goto 46
	if( reply(1:2) .eq. 'DU' ) goto 461
	if( reply(1:1) .eq. 'D' ) goto 37
	if( reply(1:1) .eq. 'M' ) goto 38
	if( reply(1:1) .eq. 'E' ) goto 42
	if( reply(1:1) .eq. 'L' ) goto 41
	if( reply(1:2) .eq. 'AO' ) goto 44
	if( reply(1:2) .eq. 'A' ) goto 45
	if( reply(1:1) .eq. 'X' ) goto 46
	if( reply(1:1) .eq. 'Y' ) goto 56
	if( reply(1:1) .eq. '@' ) goto 47
	if( reply(1:1) .eq. 'U' ) goto 25
	if( reply(1:1) .eq. '*' ) goto 26
	if( reply(1:1) .eq. 'S' ) goto 28
	if( reply(1:1) .eq. 'O' ) goto 48
	if( reply(1:1) .eq. '1' ) goto 49
	if( reply(1:1) .eq. 'N' ) goto 50
	if( reply(1:1) .eq. 'W' ) goto 52
	if( reply(1:1) .eq. 'K' ) goto 53
	if( reply(1:1) .eq. 'C' ) goto 54
	if( reply(1:1) .eq. 'F' ) goto 60



* ACTIONS ------------------------

      if( reply(1:1) .eq. 'Q' ) then
	reply = 'N'
	call inq( 'Do you REALLY want to quit BOWL?', reply, reply )
	call upcase( reply, reply )
	if( reply(1:1) .ne. 'Y' ) goto 20
	write(*,*)
	write(*,*) 'Well. BOWL me over.'
	stop 'Fini.'
	goto 20
      end if

* plot options
      if( reply(1:1) .eq. 'P' ) then
	call nospace( reply )
	goto 90
      end if

* polite
	write(*,*) '** Command not recognised. Please try again.'
	reply = '?'
	goto 20

* update model -----------------------------------------
   25	write(*,*) 'Updating model'
	call initpar
	call setpar
	bof = bofcalc( p_p )
	goto 4
c	goto 5

* change parameter values -------------------------------

* Hlamp
   26	write(*,*) 'Lamp height Hlamp =', hlamp_rg, ' Rg =', hlamp, ' ltd'
	write(*,*) '(Enter negative value for Rg units)'
	old = hlamp
	call inqr( 'Hlamp (light days)', old, val )
      if( val .ge. 0. ) then
	hlamp = val
	hlamp_rg = hlamp / bhRg
      else
	hlamp_rg = -val
	hlamp = hlamp_rg * bhRg
      end if
	write(*,*) old, ' =>', hlamp
	if( old .ne. hlamp ) goto 4
c	goto 5
	goto 4

* isco radius
   28	write(*,*) 'ISCO radius =', risco_rg, ' Rg =', risco, ' ltd'
	write(*,*) '(Enter negative value for Rg units)'
	old = risco
	call inqr( 'Risco (light days)', old, val )
      if( val .ge. 0. ) then
	risco = val
	risco_rg = risco / bhRg
      else
	risco_rg = -val
	risco = bhRg * risco_rg
      end if
	write(*,*) old, ' =>', risco
	if( old .ne. risco ) goto 4
c	goto 5
	goto 4

* H1
   30	write(*,*) 'Disc height H(R1) =', h1_rg, ' Rg =', h1, ' ltd'
      if( r1 .ne. 0. ) then
	write(*,*) 'at ref radius R1  =', r1_rg, ' Rg =', r1, ' ltd'
	r0 = r1
      else
	write(*,*) 'at ref radius R1 = Rout =', rout_rg, ' Rg =', rout, ' ltd'
	r0 = rout
      end if
	write(*,*) 'Lamp height', hlamp/bhRg, ' Rg =', hlamp, ' ltd'
	write(*,*) '(Enter negative value for Rg units)'
	old = h1
	call inqr( 'H1 (light days)', old, val )
      if( val .ge. 0. ) then
	h1 = val
	h1_rg = h1 / bhRg
      else
	h1_rg = - val
	h1 = bhRg * h1_rg
      end if
	h1_r1 = h1 / r0
	write(*,*) old, ' =>', h1
      if( primary_dh ) then
	dh_r1_pct = ( h1 - hlamp ) / r0 * 100.
	write(*,*) '(H-Hx)/R =', dh_r1_pct, ' %'
      end if
	if( old .ne. h1 ) goto 4
c	goto 5
	goto 4

* H/R at ref radius r0
   31	write(*,*) 'Disc H/R =', h1_r1_pct, ' %'
	write(*,*) '(H-Hx)/R =', dh_r1_pct, ' %'
      if( r1 .ne. 0. ) then
	write(*,*) 'at ref radius R1  =', r1_rg, ' Rg =', r1, ' ltd'
	r0 = r1
      else
	write(*,*) 'at ref radius R1 = Rout =', rout_rg, ' Rg =', rout, ' ltd'
	r0 = rout
      end if
	old = h1_r1_pct
	call inqr( 'H/R (%)', old, val )
      if( val .ge. 0. ) then
	h1_r1_pct = val
	h1_r1 = h1_r1_pct / 100.
	h1 = h1_r1 * r0
	h1_rg = h1 / bhRg
      end if
	write(*,*) ' H/R =', old, ' =>', h1_r1_pct, ' %'
      if( primary_dh ) then
	pct = dh_r1_pct
	dh_r1_pct = h1_r1_pct - hlamp / r0 * 100.
	write(*,*) '(H-Hx)/R =', pct, ' =>', dh_r1_pct, ' %'
      end if

	if( old .ne. h1_r1_pct ) goto 4
c	goto 5
	goto 4

* beta
   32	old = bet
	call inqr( 'Beta = dlnH/dlnR :', old, val )
	if( val .ge. 0. ) bet = val
	write(*,*) old, ' =>', bet
	if( old .ne. bet ) goto 4
c	goto 5
	goto 4

* fcol
   33	old = fcol
	call inqr( 'fcol = Tc/Teff', old, val )
	if( val .ge. 0. ) fcol = val
	write(*,*) old, ' =>', fcol
	if( old .eq. fcol ) goto 5
	goto 4

* inclination
   34	old = dinc
	write(*,*) 'Current inclination', dinc, ' (deg)'
	write(*,*) 'Min and Max inclination (deg)', dincmn, dincmx
	call inqr( 'Inclination (deg)', old, val )
	if( val .ge. dincmn .and. val .le. dincmx ) dinc = val
	write(*,*) old, ' =>', dinc
	if( old .ne. dinc ) goto 4
c	goto 5
	goto 4

* max inclination
   35	old1 = dincmn
	old2 = dincmx
	val1 = dincmn
	val2 = dincmx
	write(*,*) 'Current inclination', dinc, ' (deg)'
	call inq2r( 'Min and Max  inclination (deg)', val1, val2 )
      if( val1 .ge. 0. .and. val1 .le. 90.
     &	.and. val2 .ge. 0. .and. val2 .le. 90. ) then
	dincmn = min( val1, val2 )
	dincmx = max( val1, val2 )
      end if
	write(*,*) old1, ' =>', dincmn
	write(*,*) old2, ' =>', dincmx
	if( old1 .ne. dincmn .or. old2 .ne. dincmx ) goto 4
c	goto 5
	goto 4

* redshift
   36	old = redshift
	call inqr( 'Redshift', old, val )
	if( val .ge. 0. ) redshift = val
	write(*,*) old, ' =>', redshift
	if( old .eq. redshift ) goto 5
* force re-evaluation of distance when secondary is called
	dmpc = -1.
	goto 4

* distance
   37	old = dmpc
	write(*,*) 'Enter 0 for D_L(z)=', dmpcz
	call inqr( 'D_L / Mpc', old, val )
	if( val .ge. 0. ) dmpc = val
	if( dmpc .eq. 0. ) dmpc = dmpcz
	write(*,*) old, ' =>', dmpc
	if( old .eq. dmpc ) goto 5
	goto 4

* mass
   38	old = bhM_Msun
	write(*,*) 'Enter <0 for log( Mbh / Msun ), currently', alog10( old )
	call inqr( 'Mbh / Msun', old, val )
      if( val .gt. 0. ) then
	bhM_Msun = val
      else if( val .lt. 0. ) then
	bhM_Msun = 10. ** abs( val )
      end if
	write(*,*) old, ' =>', bhM_Msun
	if( old .eq. bhM_Msun ) goto 5
	goto 4

* systematic flux error
   39	old = dexsed
	write(*,*) '10 ** (', dexsed, ' ) =', 100. * sigsed, ' %'
	write(*,*) 'Enter negative value for percent.'
	call inqr( 'Systematic Flux Error (dex)', old, val )
	if( val .eq. old ) goto 5
      if( val .ge. 0. ) then
	dexsed = val
	write(*,*) old, ' =>', dexsed, ' (dex)'
	old = sigsed * 100.
	sigsed = 10.d0 ** dexsed - 1.d0
	write(*,*) old, ' =>', 100. * sigsed, ' %'
      else
	old = sigsed
	sigsed = abs( val ) / 100.
	write(*,*) 100. * old, ' =>', 100. * sigsed, ' %'
	old = dexsed
	dexsed = alog10( 1. + sigsed )
	write(*,*) old, ' =>', dexsed, ' (dex)'
      end if
	goto 4

* systematic lag
  391	old = syslag
	call inqr( 'Systematic Lag Error (d)', old, val )
	if( val .ge. 0. ) syslag = val
	write(*,*) old, ' =>', syslag
	if( old .eq. syslag ) goto 5
	goto 4


* L/Ledd
   40	old = bhL_Ledd
	call inqr( 'Ldisc / Ledd', old, val )
      if( val .ge. 0. ) then
	write(*,*) 'Ldisc/Ledd', bhL_Ledd, ' =>', val
	bhL_Ledd = val
* revise primary parameter bhMsun_yr
	val = bhL_Ledd * bhLedd_Lsun
	write(*,*) 'Ldisc/Lsun', bhL_Lsun, ' =>', val
	bhL_Lsun = val
	val = bhL_Lsun / eta / clight / clight * ( sunL / sunM ) * year
	write(*,*) '      Mdot', bhMsun_yr, ' =>', val
	bhMsun_yr = val
      end if
	write(*,*) old, ' =>', bhL_Ledd
	if( old .eq. bhL_Ledd ) goto 5
	goto 4

* L/Lsun
   41	old = bhL_Lsun
	call inqr( 'Ldisc / Lsun', old, val )
      if( val .ge. 0. ) then
	bhL_Lsun = val
* revise primary parameter L/Ledd
c	bhL_Ledd = bhL_Lsun / bhLedd_Lsun
* revise primary parameter bhMsun_yr
	bhMsun_yr = bhL_Lsun / eta / clight / clight * ( sunL / sunM ) * year
      end if
	write(*,*) old, ' =>', bhL_Lsun
	if( old .eq. bhL_Lsun ) goto 5
	goto 4

* eta = disc radiative efficiency, eta = (1/4)(Rs/Rin) < 1/2
   42	old = eta
	val = 0.5 / risco_rg
	write(*,*) 'eta = Ldisc/(Mdot c^2) = Rg/Rin/2 =', val, ' < 1/2'
	call inqr( 'eta = ', old, val )
      if( val .le. 0.5 ) then
	eta = val
* revise primary parameter Risco_rg
	risco_rg = 0.5 / eta
	risco = risco_rg * bhRg
      end if
	write(*,*) old, ' =>', eta
	if( old .eq. eta ) goto 5
	goto 4

* Mdot = accretion rate
   43	old = bhMsun_yr
	call inqr( 'Mdot (Msun/yr)  = ', old, val )
      if( val .gt. 0. ) then
	bhMsun_yr = val
      end if
	write(*,*) old, ' =>', bhMsun_yr
	if( old .eq. bhMsun_yr ) goto 5
	goto 4


* albedo
   44	old = albedo
	call inqr( 'albedo', old, val )
	if( val .ge. 0. .and. val .le. 1. ) albedo = val
	write(*,*) old, ' =>', albedo
	if( old .eq. albedo ) goto 5
	goto 4

* alpha
   45	old = alp
	call inqr( 'alpha = -dlnT/dlnR ', old, val )
	if( val .ge. 0. ) alp = val
	write(*,*) old, ' =>', alp
	if( old .eq. alp ) goto 5
	goto 4

* Lamp boost
   46 if( primary_xlamp ) then
	write(*,*) 'Xlamp', xlamp, ' primary'
	write(*,*) 'DlogT', dlogt, ' secondary'
      else
	write(*,*) 'Xlamp', xlamp, ' secondary'
	write(*,*) 'DlogT', dlogt, ' primary'
	goto 66
      end if
	old = xlamp
	call inqr( 'X = L_lamp/(eta Mdot c^2)', old, val )
	if( val .ge. 0. ) xlamp = val
	write(*,*) old, ' =>', xlamp
	if( old .eq. xlamp ) goto 5
	goto 4

* Dust (rest-frame SMC-like) E(B-V)
  461	write(*,*) 'E(B-V) for rest-frame SMC-like AGN dust.'
	ebmv = ebmv_smc
	if( ebmv .gt. 0. ) write(*,*) '( Applied already', ebmv, ' )'
	old = ebmv_agn
	call inqr( 'E(B-V)', old, val )
	ebmv_agn = val
	write(*,*) old, ' =>', ebmv_agn
	if( old .eq. ebmv_agn ) goto 5
	goto 4

* Flat Lamp boost
   56	old = ylamp
	call inqr( 'Flat lamp F = L_lamp/(eta Mdot c^2)', old, val )
	if( val .ge. 0. ) ylamp = val
	write(*,*) old, ' =>', ylamp
	if( old .eq. ylamp ) goto 5
	goto 4

* DlogT
   66	old = dlogt
	write(*,*) '** Keep DlogT in range (0,1) **'
	call inqr( 'DlogT (dex)', old, val )
	dlogt = max( 0., min( 1., val ) )
	if( dlogt .ne. val ) write(*,*) '** Requested DlogT =', val, ' =>', dlogt
	write(*,*) old, ' =>', dlogt
	if( old .eq. dlogt ) goto 5
	goto 4

* toggle inner disc torque
   47	old = qin
	call inqr( 'Inner disc torque ? (>0)', old, val )
	if( val .ge. 0. ) qin = val
	write(*,*) old, ' =>', qin
	if( old .eq. qin ) goto 5
	goto 4

* outer disc radius Rout
   48	old = rout
	write(*,*) ' Rin =', risco_rg, ' Rg'
	write(*,*) 'Rout =', rout_rg, ' Rg =', rout, ' ltd'
	write(*,*) '(Enter Rout<0 for Rg units)'
	call inqr( 'Rout (light days)', old, val )
      if( val .ge. 0. ) then
	rout = val
	rout_rg = rout / bhRg
      else
	rout_rg = -val
	rout = bhRg * rout_rg
      end if

* trap Rout < Risco
      if( rout_rg .le. risco_rg ) then
	write(*,*) '** INVALID OUTER RADIUS'
	write(*,*) '** Rout', rout_rg, ' < Risco', risco_rg
	rout = old
	rout_rg = rout / bhRg
      end if
	write(*,*) old, ' =>', rout
	if( old .eq. rout ) goto 5
	goto 4

* reference radius R1
   49	write(*,*) 'Reference radius R1  =', r1_rg, ' Rg =', r1, ' ltd'
        write(*,*) ' < 0 for Rg units '
	write(*,*) ' = 0 for R1 = Rout =', rout_rg, ' Rg =', rout, ' ltd'
	old = r1
        call inqr( 'R1 (light days)', old, val )
      if( val .ge. 0. ) then
        r1 = val
        r1_rg = r1 / bhRg
      else
        r1_rg = - val
        r1 = bhRg * r1_rg
      end if
	write(*,*) old, ' =>', r1
	if( old .eq. r1 ) goto 5
        goto 4

* number of radii
   50	write(*,*) 'Max radii', maxr
	n = nr
	call inqi( 'Number of radii Nr', n, i )
	if( i .gt. 1 ) nr = min( maxr, i )
	write(*,*) n, ' =>', nr
* no change
	if( n .eq. i ) goto 5
* change
	write(*,*) 'Radius grid changes. Nr=', n, ' =>', nr
	needgeom = .true.
	goto 4

* ripple amplitude
   52	val1 = wamp
	val2 = wpow
	old1 = val1
	old2 = val2
	call inq2r( 'Wave Amplitude (-1,1) and power-law index', val1, val2 )
      if( abs( val1 ) .gt. 1. ) then
	write(*,*) '** Invalid wave amplitude', val1, ' nothing changed.'
      else
	wamp = min( max( val1, -1. ), 1. )
	wpow = val2
      end if
	write(*,*) old1, ' =>', wamp
	write(*,*) old2, ' =>', wpow
	if( old1 .eq. wamp .and. old2 .eq. wpow ) goto 5
	goto 4

* ripple wave number
   53	val1 = wnum
	val2 = wnpow
	old1 = val1
	old2 = val2
	write(*,*) 'Ripple (cycles / light day) K = K1(R/R1)^powk'
	call inq2r( 'K1 and powk', val1, val2 )
      if( val1 .eq. 0. ) then
	write(*,*) '** Invalid wavenumber', val1, ' nothing changed.'
      else
	wnum = val1
	wnpow = val2
      end if
	write(*,*) old1, ' =>', wnum
	write(*,*) old2, ' =>', wnpow
	if( old1 .eq. wnum .and. old2 .eq. wnpow ) goto 5
	goto 4

* ripple crest peakiness
   54	write(*,*) 'Ripple crest peakiness : 2(cos/2)^c'
	old = crest
	call inqr( 'Crest peakiness (>0) ', old, val )
	if( val .gt. 0. ) crest = val
	write(*,*) old, ' =>', crest
	if( old .eq. crest ) goto 5
	goto 4

* MCMC fit  ----------------------------------------------
   60	write(*,*)
	write(*,*) '*** MCMC parameter sampling ***'
	call mcmcfit
	goto 4


* BOWL fit plot ----------------------------------------------
   90	write(*,*)
	call bowlplot

	goto 10
	end


* =================================================
	subroutine bowlplot
* multi-panel plot of bowl fit
* In:	reply	c* Plot control option list
*
* 2022 Aug Keith Horne @ St Andrews
* 2022 Aug KDH @ Prague - lightsec
* 2023-06-14 KDH @ StA - needgeom
* 2023-06-27 KDH @ StA - writeout
* 2023-07-13 KDH @ StA - simpleplot - omit fiducial T,R,w
* 2023-07-13 KDH @ StA - omit chi^2 and BoF if simpleplot
* 2023-07-13 KDH @ StA - revise SED window if simpleplot
* 2023-08-30 KDH @ StA - report reflag on lag spectrum plot
* 2023-12-01 KDH @ StA - best fit option
* 2023-12-02 KDH @ StA - toggle off options
* 2023-12-09 KDH @ StA - median parameters
* 2023-12-09 KDH @ StA - dlogt
* 2024-03-07 KDH @ StA - MCMC sub-samples
* 2024-03-07 KDH @ StA - divide SED wavelengths wobs_b by (1+z)
* 2024-05-07 KDH @ StA - adjust histogram scaling to avoid pile-up
* 2025-04-04 KDH @ StA - control lwide

	use all
	use agn
	use constants
	use cosmology
	use pgplot
	logical pgagain

	logical lightsec
	logical writeout
	logical yes

* colour indices for various plot features
	ktvmin = kred
	ktvmax = kblack

	ktin = kblue
	ktmax = kmagenta
	ktmin = korange
	ktout = kgreen

	krin = kblue
	krout = kblue

	khout = kblue
	klamp = kmagenta

	kfnt = kred
	kbrt = kblue
	krim = kgreen

	if( lwide .le. 0 ) lwide = 2
	if( csize .le. 0. ) csize = 2


 	write(*,*)
	write(*,*) '-- bowlplot ----------------'


* time units

	lightsec = rout .lt. 0.1
	write(*,*) 'Rout/c', rout, ' days', rout/day, ' seconds'

* compute BoF
	bof = bofcalc( p_p )

* KDH : IS REPLY KNOWN?

   10	write(*,*) ' reply : ', reply(:len1(reply))

      if( len1( reply ) .eq. 1 .and. reply(1:1) .eq. 'P' ) then
	reply = plotoption
      else
	reply = reply(2:)
      end if
	write(*,*) ' reply : ', reply(:len1(reply))

* plot options
	call nospace( reply )
	call upcase( reply, reply )

* toggle off plot options
	i = index( reply, '-' )
	write(*,*) 'index( ', reply(:len1(reply)), ', - ) =', i
      if( i .gt. 0 ) then
	write(*,*) '  Simpleplot',  simpleplot, ' => ', .false.
	write(*,*) ' Bestfitplot', bestfitplot, ' => ', .false.
	write(*,*) '  Medianplot',  medianplot, ' => ', .false.
	simpleplot = .false.
	bestfitplot = .false.
	medianplot = .false.
	reply(i:i) = ' '
      end if

* simple plot
	i =  index( reply, 'S' )
	write(*,*) 'index( ', reply(:len1(reply)), ', S ) =', i
      if( i .gt. 0 ) then
	simpleplot = .true.
	reply(i:i) = ''
	write(*,*) 'Simpleplot', simpleplot, ' reply : ', reply(:len1(reply))
      end if

* median plot
	i =  index( reply, 'M' )
	write(*,*) 'index( ', reply(:len1(reply)), ', M ) =', i
      if( i .gt. 0 ) then
	medianplot = .true.
	bestfitplot = .false.
	reply(i:i) = ''
	write(*,*) 'Medianplot', medianplot, ' reply : ', reply(:len1(reply))
      end if

* bestfit plot
	i =  index( reply, 'B' )
	write(*,*) 'index( ', reply(:len1(reply)), ', B ) =', i
      if( i .gt. 0 ) then
	bestfitplot = .true.
	reply(i:i) = ''
	write(*,*) 'Bestfitplot', bestfitplot, ' reply : ', reply(:len1(reply))
      end if

* output files
	i = index( reply, 'O' )
	write(*,*) 'index( ', reply(:len1(reply)), ', O ) =', i
      if( i .gt. 0 ) then
	writeout = .true.
	reply(i:i) = ''
	write(*,*) 'writeout', writeout, ' reply : ', reply(:len1(reply))
      end if


* random samples to plot
	if( nsp .lt. 0 ) nsp = nshow
	nsp = max( 0, min( nsamp, nsp ) )

* plot menu
	write(*,*) 'Plot Panel Options:'
	write(*,*) ' F ... Fnu(lam)'
	write(*,*) ' L ... Lag(lam)'
	write(*,*) ' T ... T(r)'
	write(*,*) ' H ... H(r)'
	write(*,*) ' D ... Delay maps at', nb, ' wavelengths'
	write(*,*) ' W ... Delay map at', nint(wrest), ' A'
	write(*,*) ' 2 ... 2-panel FT'
	write(*,*) ' 4 ... 4-panel FLTH'
	write(*,*) ' 6 ... 6-panel FLTHWD'
	write(*,*)
	write(*,*) 'Plot Options:'
      if( nfit .gt. 0 .and. nits .gt. 0 ) then
	write(*,*) '      MCMC(', 1,     ' ) BoF', bof_i( 1 )
	write(*,*) '      MCMC(', nits,  ' ) BoF', bof_i( nits )
      if( ibest .gt. 0 ) then
	write(*,*) '      MCMC(', ibest, ' ) BoF', bof_i( ibest ), ' best'
      end if
	write(*,*) ' B ... Best parameters', bestfitplot
	write(*,*) ' M ... Median parameters', medianplot
	write(*,*) ' R ... Random samples', nsp, ' of', nsamp, ' thin', nthin
      end if
	write(*,*) ' S ... Simple plot', simpleplot
	write(*,*) ' - ... reset (BMS) => ', .false.
	write(*,*) ' O ... Output to files', writeout
	write(*,*) ' C ... Custom LineWidth', lwide, ' CharHght', csize
	write(*,*)
	call inq( 'What to plot? ', reply, reply )

	call nospace( reply )
	call upcase( reply, reply )

* toggle options and strip from the reply
	i = index( reply, '-' )
      if( i .gt. 0 ) then
	simpleplot = .false.
	bestfitplot = .false.
	medianplot = .false.
	reply = reply(:i-1) // reply(i+1:)
      end if

* random samples
	i = index( reply, 'R' )
      if( i .gt. 0 ) then
	reply = reply(:i-1) // reply(i+1:)
	call inqi( 'Random samples to plot', nsp, j )
	if( j .ge. 0 ) nsp = min( nsamp, j )
      end if


* custom line width and character height
	i = index( reply, 'C' )
      if( i .gt. 0 ) then
	reply = reply(:i-1) // reply(i+1:)
	w = lwide
	c = csize
	call inq2r( 'Custom LineWidth and CharHgt', w, c )
	i = nint( w )
	if( i .gt. 0 ) lwide = i
	if( c .gt. 0. ) csize = c
      end if


* simpler plot
	i = index( reply, 'S' )
      if( i .gt. 0 ) then
	simpleplot = .true.
	reply = reply(:i-1) // reply(i+1:)
      end if
* write-out
	i = index( reply, 'O' )
      if( i .gt. 0 ) then
	writeout = .true.
	reply = reply(:i-1) // reply(i+1:)
      end if

* best parameters
      if( nfit .gt. 0 .and. nits .gt. 0 ) then
	i = index( reply, 'B' )
      if( i .gt. 0 ) then
	bestfitplot = .true.
	medianplot = .false.
	reply = reply(:i-1) // reply(i+1:)
      end if
* median parameters
	i = index( reply, 'M' )
      if( i .gt. 0 ) then
	medianplot = .true.
	bestfitplot = .false.
	reply = reply(:i-1) // reply(i+1:)
      end if

      end if

* plot panel options
* 6-panel plot
	i = index( reply, '6' )
	if( i .gt. 0 ) reply = reply(:i-1) // 'FLTHDW'
* 4-panel plot
	i = index( reply, '4' )
	if( i .gt. 0 ) reply = reply(:i-1) // 'FLTH'
* 2-panel plot
	i = index( reply, '2' )
	if( i .gt. 0 ) reply = reply(:i-1) // 'FT'

	if( reply(1:1) .ne. 'P' ) reply = 'P' // reply(:len1(reply))
	ipanh = 0
	ipant = 0
	ipanf = 0
	ipantau = 0
	ipanpsi = 0
	ipanpsib = 0
* assign panels
	call word1( reply, i1, i2 )
	n = 0
      if( i2 .ge. i1 ) then
      do i = i1, i2
	n = n + 1
* H(R)
      if(      reply(i:i) .eq. 'H' ) then
	write(*,*) 'Panel', n, ' H(R)'
	ipanh = n
* F(lam)
      else if( reply(i:i) .eq. 'F' ) then
	write(*,*) 'Panel', n, ' Flux(lam)'
	ipanf = n
* T(R)
      else if( reply(i:i) .eq. 'T' ) then
	write(*,*) 'Panel', n, ' T(R)'
	ipant = n
* tau(lam)
      else if( reply(i:i) .eq. 'L' ) then
	write(*,*) 'Panel', n, ' Lag(lam)'
	ipantau = n
* delay map (1-band)
      else if( reply(i:i) .eq. 'D' ) then
	write(*,*) 'Panel', n, ' Delay Maps at', nb, ' wavelengths'
	ipanpsib = n
* delay maps (multi-band)
      else if( reply(i:i) .eq. 'W' ) then
	write(*,*) 'Panel', n, ' Delay Map at', nint(wrest), ' A'
	ipanpsi = n
      else
	write(*,*) 'Panel ', reply(i:i), ' not recognised.'
	n = n - 1
      end if
      end do
* no valid panels
      if( n .eq. 0 ) then
	write(*,*) '** NOTHING TO PLOT - REVERT TO DEFAULT'
	reply = 'P'
	if( simpleplot ) reply = 'PS'
	goto 10
      end if
* panel layout
	nx = 1
	ny = 1
	if( n .eq. 2 ) ny = 2
	if( n .eq. 3 ) ny = 3
	if( n .ge. 4 ) nx = 2
	if( n .eq. 4 ) ny = 2
	if( n .ge. 5 ) ny = 3
	if( n .ge. 7 ) ny = 4
	if( n .ge. 9 ) nx = 3
* KDH : could use algorithm
	npan = n
	write(*,*) 'Panels', npan, ' Nx', nx, ' Ny', ny
      end if


* pluck best fit parameters
	nstow = 0
      if( ( bestfitplot .or. medianplot ) .and. ibest .gt. 0 ) then
	nstow = 1
* save current fit parameters
	bofsave = bof
	write(*,*) 'Stow current parameters. BoF=', bof
	call setpar
      do i = 1, npar
	save_p( i ) = p_p( i )
      end do
* best parameters
      if( bestfitplot ) then
	write(*,*) '** Pluck best fit parameters from MCMC step', ibest
	call pluck( ibest )
* median parameters
      else if( medianplot ) then
	write(*,*) '** Pluck median parameters from MCMC chain'
	call pluck( 0 )
      end if
      end if

* open
  100	call pgstart( nx, ny, nouse )
	call pgbgw
	if( csize .le. 0. ) csize = 2.
	if( lwide .le. 0 ) lwide = 2
	lwider = lwide * 3
	call pgsch( csize )
	call pgslw( lwide )

* select samples to plot
      do i = 1, nsamp
	plot_i( i ) = .false.
      end do
* toggle on nsp samples
      if( nsp .ge. 1 .and. nsp .le. nsamp ) then
	write(*,*) 'Plot', nsp, ' of', nsamp, ' samples.'
      do i = 1, nsp
	p = ( i - 1. ) / max( 1, nsp - 1 )
	p = 0.5001 * ( 1. - p ) + p * ( nsamp + 0.4999 )
	ip = nint( p )
	plot_i( ip ) = .true.
	write(*,*) i, ip, iter_i( ip )
      end do
      end if

* panels
      do ipan = 1, npan
* defaults
	logx = 1
	logy = 1
	np = nr

* compute model faint and bright flux spectra
	call sedcalc
* compute model lag spectrum
	call lagcalc

* ref radius
	r0 = r1
	if( r0 .le. 0. ) r0 = rout
* delay range of annuli
      do i = 1, nr
	r = r_r( i )
	h = h_r( i )
	x = r / r0
	ip = min( i + 1, nr )
	im = max( ip- 1, 1 )
	n = ip - im
* KDH : SHOULD R(I) BE INNER RADIUS OF ANNULUS?
	dr = ( r_r( ip ) - r_r( im ) ) / n
	dh = ( h_r( ip ) - h_r( im ) ) / n
	ds = sqrt( dh * dh + dr * dr )
* height above lamp
	dz = h - hlamp
* distance to lamp
	ddlamp = dz * dz + r * r
	dlamp = sqrt( ddlamp )
* delay range
	tau = dlamp - dz * cosi
	add = r * sini
	taulo = tau - add
	tauhi = tau + add
* stow
	tau_r( i ) = tau
	taulo_r( i ) = taulo
	tauhi_r( i ) = tauhi
      end do

* AGN name
	name = name_a( iagn )
	title = name( : len1( name ) )

* geometry
      if( ipan .eq. ipanh ) then
	call append( title, 'Disk Geometry' )
	xlabel = 'radius R/c (d)'
	ylabel = 'height H/c (d)'
	u  = 1.
      if( lightsec ) then
	u = day
	xlabel = 'radius R/c (s)'
	ylabel = 'height H/c (s)'
      end if
	logy = 0
	logx = 0
      do i = 1, np
	plot1( i ) = r_r( i ) * u
	plot2( i ) = h_r( i ) * u
	plot3( i ) = h0_r( i ) * u
      end do

* temperature profile
      else if( ipan .eq. ipant ) then

	call append( title, 'Disk Temperature Profile' )
c	title = 'Disk Temperature Profile'
	ylabel = 'T\deff\u (10\u4\dK)'
	powtemp = 4
	add = - powtemp
	xlabel = 'radius R/c (d)'
	u = 1.
      if( lightsec ) then
	u = day
	xlabel = 'radius R/c (s)'
      end if

      do i = 1, nr
	plot1( i ) = alog10( r_r( i ) * u )
* total temp
	plot2( i ) = alog10( t_r( i ) ) + add
* viscous temp
	plot3( i ) = alog10( tv_r( i ) ) + add
* power-law temp
	plot4( i ) = alog10( t0_r( i ) ) + add
* irradiation temp
	plot5( i ) = alog10( tx_r( i ) ) + add
* torque temp
	plot6( i ) = alog10( tq_r( i ) ) + add
* colour temp
	plot7( i ) = alog10( tc_r( i ) ) + add
      end do

* REPORT FOR TESTING
	r  = 10. ** plot1( 1 )
	t  = 10. ** ( plot2( 1 ) - add )
	tv = 10. ** ( plot3( 1 ) - add )
	t0 = 10. ** ( plot4( 1 ) - add )
	tx = 10. ** ( plot5( 1 ) - add )
	tq = 10. ** ( plot6( 1 ) - add )
	tc = 10. ** ( plot7( 1 ) - add )
	write(*,*) 'Rin', r, ' T', nint( t )
	write(*,*) ' Tv', nint( tv ), ' T0', nint( t0 )
	write(*,*) ' Tq', nint( tq ), ' Tx', nint( tx )
	write(*,*) ' Tq(1)', nint( tq_r(1) ), ' Tq(Nr)', nint( tq_r(nr) )
	write(*,*) ' plot6(1)', plot6(1), ' plot6(Nr)', plot6( nr )
	write(*,*) ' Tc(1)', nint( tc_r(1) ), ' Tc(Nr)', nint( tc_r(nr) )
	write(*,*) ' plot7(1)', plot7(1), ' plot7(Nr)', plot7( nr )

* fnu spectrum
      else if( ipan .eq. ipanf ) then

	call append( title, 'Disk Spectrum' )
	xlabel = 'rest wavelength \gl (\A)'
	ylabel = 'F\d\gn\u (mJy)'

* model SEDs
	np = nw
	call mnmx( nw, f_w, fmn, fmx )
	floor = max( fmn / 10., fmx * 1e-10 )
	if( floor .le. 0. ) floor = 0.1
      do i = 1, nw
	w = w_w( i )
	f = f_w( i )
	f0 = f0_w( i )
	f1 = f1_w( i )
	f2 = f2_w( i )
	plot1( i ) = alog10( w )
	plot2( i ) = alog10( max( floor, f ) )
	plot3( i ) = alog10( max( floor, f0 ) )
	plot4( i ) = alog10( max( floor, f1 ) )
	plot5( i ) = alog10( max( floor, f2 ) )
      end do

* 2024-04-01 KDH : DEBUG REPORT
      IF( .FALSE. ) THEN
	call mnmx( nw, f_w, fmn, fmx )
	call mnmx( nw, f0_w, f0mn, f0mx )
	call mnmx( nw, f1_w, f1mn, f1mx )
	call mnmx( nw, f2_w, f2mn, f2mx )
	write(*,*) ' Fnu/mJy',  fmn,  fmx, ' bright disk+rim'
	write(*,*) ' Fnu/mJy', f1mn, f1mx, ' bright disk'
	write(*,*) ' Fnu/mJy', f2mn, f2mx, ' bright rim'
	write(*,*) ' Fnu/mJy', f0mn, f0mx, '  faint disk'
	call mnmx( nw, plot2, fmn, fmx )
	call mnmx( nw, plot3, f0mn, f0mx )
	call mnmx( nw, plot4, f1mn, f1mx )
	call mnmx( nw, plot5, f2mn, f2mx )
	write(*,*) 'Log(Fnu/mJy)',  fmn,  fmx, ' bright disk+rim'
	write(*,*) 'Log(Fnu/mJy)', f1mn, f1mx, ' bright disk'
	write(*,*) 'Log(Fnu/mJy)', f2mn, f2mx, ' bright rim'
	write(*,*) 'Log(Fnu/mJy)', f0mn, f0mx, '  faint disk'
      END IF

* lag spectrum
      else if( ipan .eq. ipantau ) then

	call append( title, 'Lag Spectrum' )
	xlabel = 'rest wavelength \gl (\A)'
	ylabel = 'delay \gt (d)'

	u = 1.
      if( lightsec ) then
	u = day
	ylabel = 'delay \gt (s)'
      end if

	np = nb
	call mnmx( nb, tau_b, tmn, tmx )
	logy = 0
	logx = 1
	if( tmx - tmn .gt. tmx / 15. ) logx = 0
	floor = max( tmn, tmx * 1e-10 )
	if( floor .le. 0. ) floor = 0.1
	floor = floor * u
      do i = 1, nb
	 x = wobs_b( i )
	 y =  tau_b( i ) * u
	y0 = tau0_b( i ) * u
	y1 = tau1_b( i ) * u
	y2 = tau2_b( i ) * u
	write(*,*) nint( x ), ' tau', y, y0, y1, y2
	if( logx .ne. 0 )  x = alog10( x )
	if( logy .ne. 0 )  y = alog10( max( floor, y  ) )
	if( logy .ne. 0 ) y0 = alog10( max( floor, y0 ) )
	if( logy .ne. 0 ) y1 = alog10( max( floor, y1 ) )
	if( logy .ne. 0 ) y2 = alog10( max( floor, y2 ) )
	plot1( i ) = x
	plot2( i ) = y
	plot3( i ) = y0
	plot4( i ) = y1
	plot5( i ) = y2
      end do

* delay map
      else if( ipan .eq. ipanpsi ) then
	lam = nint( wrest )
	title = ' '
	call append_i( title, lam )
	call append( title, '\A Delay Map' )
	xlabel = 'delay \gt (d)'
	ylabel = 'response \gQ(\gt) (d\u-1\d)'
	u = 1.
      if( lightsec ) then
	u = day
	xlabel = 'delay \gt (s)'
	ylabel = 'response \gQ(\gt) (s\u-1\d)'
      end if
	logx = 0
	logy = 0
	np = ntau
      do i = 1, ntau
	plot1( i ) = tau_d( i ) * u
	plot2( i ) = psi_d( i ) / u
	plot3( i ) = psi0_d( i ) / u
	plot4( i ) = psi1_d( i ) / u
      end do

* multi-band delay maps
      else if( ipan .eq. ipanpsib ) then
	title = 'Delay Maps'
	xlabel = 'delay \gt (d)'
	ylabel = 'response \gQ(\gt|\gl) (d\u-1\d)'
	u = 1.
      if( lightsec ) then
	u = day
	xlabel = 'delay \gt (s)'
	ylabel = 'response \gQ(\gt|\gl) (s\u-1\d)'
      end if

	logx = 1
	logy = 1
	np = nb
      do i = 1, nb
	call mnmx( ntau, psi_db( 1, i ), pmn, pmx )
	plot1( i ) = tau_b( i ) * u
	plot2( i ) = pmx / u
c	write(*,*) 'W', w_w(i), ' tau', tau_b(i), ' psi', pmx
      end do

* last panel option
      end if

* report
	write(*,*)
	write(*,*) 'Panel', ipan, ' Np', np
	write(*,*) ' X : ', xlabel(: len1( xlabel ) )
	write(*,*) ' Y : ', ylabel(: len1( ylabel ) )
	write(*,*) ' T : ', title (: len1( title  ) )

* viewport
	call mnmx( np, plot1, xmn, xmx )
	call mnmx( np, plot2, ymn, ymx )
	call mnmx( np, plot3, p, y )
	ymx = max( ymx, y )

* revise viewport if needed =======================
* multi-band delay maps
      if( ipan .eq. ipanpsib ) then
	xmx = alog10( xmx ) + 1.0
	xmn = alog10( xmn ) - 0.5
	ymx = alog10( ymx ) + 0.3
c	ymn = alog10( ymn ) - 2.5
	ymn = alog10( ymn ) - 3.0
      end if

* flux data
      if( ipan .eq. ipanf ) then

      if( ndat .gt. 0 ) then
	call mnmx( ndat,  datflx_b,  fmn,  fmx )
	call mnmx( ndat, datflx0_b, f0mn, f0mx )
      if( logy .ne. 0 ) then
	add = 0.1
	fmn = alog10( fmn ) - add
	fmx = alog10( fmx ) + add
	f0mn = alog10( f0mn ) - add
	f0mx = alog10( f0mx ) + add

* revise plot window for simple plots
      if( simpleplot ) then
	call mnmx( ndat, wdat_b, xmn, xmx )
	opz = 1. + redshift
	xmn = xmn / opz
	xmx = xmx / opz
	x1 = alog10( xmn )
	x2 = alog10( xmx )
	call lookup( np, plot1, x1, i1, i, p )
	call lookup( np, plot1, x2, i, i2, p )
	y1 = plot2( i1 )
	y2 = plot2( i2 )
	y3 = plot3( i1 )
	y4 = plot3( i2 )
	xmn = min( xmn, 1e3 )
	xmx = max( xmx, 1e4 )
	add = 0.1
	xmn = alog10( xmn ) - add
	xmx = alog10( xmx ) + add
	ymx = max( fmx, f0mx, y1, y2, y3, y4 )
	ymn = min( fmn, f0mn, y1, y2, y3, y4 )
      end if
	ymx = max( ymx, fmx, f0mx )
	ymn = min( ymn, fmn, f0mn )
      end if

* no rim SED
	call mnmx( np, plot4, f4mn, f4mx )
* rim SED
	call mnmx( np, plot5, f5mn, f5mx )
      if( f5mx .gt. f5mn ) then
	if( logy .ne. 0 ) y = f5mx - 0.1
	if( logy .eq. 0 ) y = f5mx * 0.8
* include peak of rim SED
	ymn = min( ymn, y )
      end if
      end if
      end if

* lag data
      if( ipan .eq. ipantau ) then

      if( nlag .gt. 0 ) then
* reference lag = model lag at reference wavelength
	opz = 1. + redshift
	w = wlag0 / opz
	call lookup( nb, wobs_b, w, i1, i2, p )
	reflag = terp1( nb, tau_b, i1, i2, p )
	write(*,*) 'tau(', nint( w ), ' ) =', reflag

	u = 1
	if( lightsec ) u = day
      do i = 1, nlag
	x = wlag_b( i )
	dx = wlagfwhm_b( i ) / 2.
	xp = x + dx
	xm = x - dx
* add reference lag to the data
	dat   =   datlag_b( i ) + reflag
	sighi = siglaghi_b( i )
	siglo = siglaglo_b( i )
	y = dat * u
	yp = y + sighi * u
	ym = y - siglo * u
      if( logx .ne. 0 ) then
	xp = alog10( xp )
	xm = alog10( xm )
      end if
      if( logy .ne. 0 ) then
	y = 10. ** ( ymn - 0.5 )
	yp = alog10( max( y, yp ) )
	ym = alog10( max( y, ym ) )
      end if

* revise window if simpleplot
      if( i .eq. 1 .and. simpleplot ) then
	xmn = min( 0., xm )
	xmx = xp
	ymn = min( 0., ym )
	ymx = yp
      end if

	xmx = max( xmx, xp )
	ymx = max( ymx, yp )
	xmn = min( xmn, xm )
	ymn = min( ymn, ym )
* next lag
      end do

* include rim lags in plot5
	call mnmx( nb, plot4,  tmn,  tmx )
	call mnmx( nb, plot5, t1mn, t1mx )
	write(*,*) ' ymx', ymx, tmx, t1mx
	y = max( ymx, tmx, t1mx )
c      if( y .lt. 3. * ymx ) then
	write(*,*) ' ymx', ymx, ' =>', y
	ymx = y
c      end if

      end if
      end if

* focus on relevant radius or delay  range
      if( logx .eq. 0 ) then
      if( ipan .eq. ipanh .or. ipan .eq. ipanpsi ) then
	tau = tauavg
	call mnmx( nb, tau_b, taumn, taumx )
      if( ipan .eq. ipanpsi ) then
	xmx = max( rout * ( 1. + sini ) , 2. * tau ) * u
	x = tau * 10. * u
	xmx = min( x, xmx )
      else
	tau = max( tau, taumx )
	call lookup( nr, tau_r, tau, i1, i2, p )
	r = terp1( nr, r_r, i1, i2, p )
	xmx = 1.5 * r * u
* change to rout
	xmx = rout * u
      end if
	if( xmx .le. 0. ) xmx = rout * u

	ymx = 0.
      do i = 1, nr
	x = plot1( i )
	if( x .le. xmx ) ymx = max( ymx, plot2( i ) )
      end do
      end if
      end if

* expand x range
	write(*,*) 'X range', xmn, xmx
	if( ipan .eq. ipanh ) xmn = - xmx
	xpand = 0.05
	if( ipan .eq. ipanh ) xpand = 0.3
	if( ipan .eq. ipant ) xpand = 0.1
	add = ( xmx - xmn ) * xpand
	xmn = xmn - add

	if( ipan .eq. ipant .and. logx .ne. 0 )
     &		add = add + 0.3
	if( ipan .eq. ipanh .and. logx .eq. 0 )
     &		add = add + ( xmx - xmn ) / 4.
	xmx = xmx + add
	write(*,*) '   ==> ', xmn, xmx

* crop log span
      if( logy .eq. 1 ) then
	dex = 3.5
	ymn = max( ymn, ymx - dex )
* include peak of Tx
      if( ipan .eq. ipant ) then
	call mnmx( np, plot5, p, y )
	ymn = min( ymn, y - 1. )
* include Tmin=1000K < 1500K dust sublimation temp
	y = 3. - powtemp
	ymn = min( ymn, y )
* go lower (500K)
	y = alog10( 500. ) - powtemp
	ymn = min( ymn, y )
* go below the power-law T(r) curve
	y = plot4( nr ) - 0.4
	ymn = min( ymn, y )
* max Tc
      if( fcol .ne. 1 ) then
	call mnmx( np, plot7, p, y )
	ymx = max( ymx, y )
      end if
	ymx = ymx + 0.2
      end if

      end if

* include lamp post
* KDH : sometimes best hlamp is off top of plot
      if( ipan .eq. ipanh .and. logy .eq. 0 ) then
	y = hlamp * 1.5 * u
	ymx = max( ymx, y )
* expand y range
	ymx = 1.5 * ymx
	ymn = 1.5 * ymn

      end if
* positive and negative r
      if( ipan .eq. ipanh .and. logx .eq. 0 ) then
	ymn = min( ymn, 0. )
c	xmn = - xmx
	xpand = 1.2
      end if

	add = ( ymx - ymn ) * xpand
	if( add .le. 0. ) add = max( abs( ymn ), abs( ymx ) )
	if( add .le. 0. ) add = 0.1
	ymn = ymn - add
	ymx = ymx + add

* space at the top
      if( ipan .eq. ipantau .or. ipan .eq. ipanf ) then
	p = 1.15
	ymx = ymn * ( 1.-p ) + p * ymx
      end if

* window
	logxy = logx * 10 + logy * 20
	write(*,*) 'pgenv(', xmn, xmx, ymn, ymx, logxy, ' )'
	call pgenv( xmn, xmx, ymn, ymx, 0, logxy )

* grid (dotted)
	call pgsls( ldot )
      if( .false. ) then
	call pgbox( 'G', 0., 0, 'G', 0., 0 )
      else
* y axis
      if( logx .eq. 0 ) then
	call pgmove( 0., ymx )
	call pgdraw( 0., ymn )
	call pgbox( 'G', 0., 0, ' ', 0., 0 )
      else
	call pgbox( 'G', 1., 0, ' ', 0., 0 )
      end if
* x axis
      if( logy .eq. 0 ) then
	call pgmove( xmn, 0. )
	call pgdraw( xmx, 0. )
	call pgbox( ' ', 0., 0, 'G', 0., 0 )
      else
	call pgbox( ' ', 0., 0, 'G', 1., 0 )
      end if
      end if
	call pgsls( 1 )

* labels
	call pglabel( xlabel, ylabel, title )

c      if( .not. simpleplot ) then
      if( nits .gt. 0 ) then
      if( bestfitplot ) then
	label = 'best('
      else if( medianplot ) then
	label = 'median('
      else
	label = 'MCMC('
	label = 'mcmc('
      end if
	call append_i( label, nits )
	call append( label, ')' )
	call nospace( label )
	call pgmtxt( 'R', 1.5, 0., 0., label )
      end if
c      end if


* stowed samples
	write(*,*) 'MCMC samples', nsamp, ' thin', nthin
      if( nsamp .ge. 1 ) then
      do is = 1, nsamp
	n = 0
* colours
	k12 = kcyan
	k13 = kyellow
	k14 = kyellow
	k15 = kyellow

* line styles
	l12 = 1
	l13 = 1
	l14 = 1
	l15 = 1
* SEDs ----------------------------
      if( ipan .eq. ipanf ) then
	n = nw_i( is )
      do i = 1, n
	plot11( i ) = alog10( w_wi( i, is ) )
	plot12( i ) = alog10( max( floor,  f_wi( i, is ) ) )
	plot13( i ) = alog10( max( floor, f0_wi( i, is ) ) )
	plot14( i ) = alog10( max( floor, f1_wi( i, is ) ) )
	plot15( i ) = alog10( max( floor, f2_wi( i, is ) ) )
      end do
* lags ----------------------------
      else if( ipan .eq. ipantau ) then
	n = nb_i( is )
      do i = 1, n
	x = wobs_bi( i, is )
	y =  tau_bi( i, is ) * u
	y0 = tau0_bi( i, is ) * u
	y1 = tau1_bi( i, is ) * u
	y2 = tau2_bi( i, is ) * u
	if( logx .ne. 0 ) x = alog10( x )
      if( logy .ne. 0 ) then
	y  = alog10( max( floor, y  ) )
	y0 = alog10( max( floor, y0 ) )
	y1 = alog10( max( floor, y1 ) )
	y2 = alog10( max( floor, y2 ) )
      end if
	plot11( i ) = x
	plot12( i ) = y
	plot13( i ) = y0
	plot14( i ) = y1
	plot15( i ) = y2
      end do

* temperatures ---------------------
      else if( ipan .eq. ipant ) then
	n = nr_i( is )
	powtemp = 4
	add = - powtemp
      do i = 1, n
	plot11( i ) = alog10( r_ri( i, is ) * u )
* total temp
	plot12( i ) = alog10( t_ri( i, is ) ) + add
* viscous temp
	plot13( i ) = alog10( tv_ri( i, is ) ) + add
* power-law temp
	plot14( i ) = alog10( t0_ri( i, is ) ) + add
* irradiation temp
	plot15( i ) = alog10( tx_ri( i, is ) ) + add
* torque temp
	plot16( i ) = alog10( tq_ri( i, is ) ) + add
* colour temp
	plot17( i ) = alog10( tc_ri( i, is ) ) + add
      end do

* fiducial temperatures and radii
* T(in)
	call imnmx( n, plot12, imn, imx )
	plot21( is ) = plot12( 1 )
	k21 = ktin
* T(max)
	plot22( is ) = plot12( imx )
	k22 = ktmax
* T(min)
	plot23( is ) = plot12( imn )
	k23 = ktmin
* T(out)
	plot24( is ) = plot12( n )
	k24 = ktout
* Tv(mx)
	call imnmx( n, plot13, imn, imx )
	plot25( is ) = plot13( imx )
	k25 = ktvmax
* Tv(out)
	plot26( is ) = plot14( n )
	k26 = ktvmin
	ky1 = 21
	ky2 = 26
* R(in)
	k = kblue
	plot31( is ) = plot11( 1 )
	k31 = krin
* R(out)
	plot32( is ) = plot11( n )
	k32 = krout
	kx1 = 31
	kx2 = 32

* geometry --------------------
      else if( ipan .eq. ipanh ) then
	n = nr_i( is )
      do i = 1, n
	plot11( i ) =  r_ri( i, is ) * u
	plot12( i ) =  h_ri( i, is ) * u
	plot13( i ) = h0_ri( i, is ) * u
      end do
	hlamp = hlamp_i( is ) * u

* R(out)
	x = plot11( n )
	plot31( is ) = +x
	plot32( is ) = -x
	k31 = krout
	k32 = krout
	kx1 = 31
	kx2 = 32
* H(out)
	y = plot12( n )
	plot21( is ) = +y
	plot22( is ) = -y
	k21 = khout
	k22 = khout
	ky1 = 21
	ky2 = 22
* Hx
	plot23( is ) = +hlamp
	plot24( is ) = -hlamp
	k23 = klamp
	k24 = klamp

	ky2 = 24
      end if

* plot selected MCMC sample curves
      if( n .gt. 0 .and. plot_i( is ) ) then
	call pgslw( 1 )

* no-rim and rim curves  SED and lag
      if( ipan .eq. ipanf .or. ipan .eq. ipantau ) then
      if( k14 .ne. 0 ) then
	call pgsci( k14 )
	call pgsls( l14 )
	call pgline( n, plot11, plot14 )
      end if
      if( k15 .ne. 0 ) then
	call pgsci( k15 )
	call pgsls( l15 )
	call pgline( n, plot11, plot15 )
      end if
* debug report
      IF( .FALSE. ) THEN
	call mnmx( n, plot14, pmn, pmx )
	write(*,*) 'plot14', pmn, pmx
	call mnmx( n, plot15, pmn, pmx )
	write(*,*) 'plot15', pmn, pmx
      END IF
      end if

* faint and bright curves
	kk = 1
      if( ipan .eq. ipanh ) then
	kk = 4
* lamp
	x = 0.
	y = hlamp * u
	call pgsci( klamp )
	call pgpoint( 1, x, +y, istar )
	call pgpoint( 1, x, -y, istar )
      end if
	x = ( xmn + xmx )
      do k = 1, kk
* reflect geometry across x and y axes
	if( k .eq. 2 ) call pgswin(   xmn,   xmx, -ymn, -ymx )
	if( k .eq. 3 ) call pgswin( xmx-x, xmn-x, -ymn, -ymx )
	if( k .eq. 4 ) call pgswin( xmx-x, xmn-x,  ymn,  ymx )
	call pgsci( k12 )
	call pgsls( l12 )
	call pgline( n, plot11, plot12 )
	call pgsci( k13 )
	call pgsls( l13 )
	call pgline( n, plot11, plot13 )
	call pgswin( xmn, xmx, ymn, ymx )
      end do

* debug report
      IF( .FALSE. ) THEN
	call mnmx( n, plot11, pmn, pmx )
	write(*,*) 'plot11', pmn, pmx
	call mnmx( n, plot12, pmn, pmx )
	write(*,*) 'plot12', pmn, pmx
	call mnmx( n, plot13, pmn, pmx )
	write(*,*) 'plot13', pmn, pmx
      END IF
	call pgslw( lwide )
      end if

* next sample is
      end do
      end if

* KDH : MAKE A SUBROUTINE
* histograms on x and y axes
      if( ipan .eq. ipant .or. ipan .eq. ipanh ) then
	ns = nsamp
      if( ns .gt. 0 ) then
	iy = 1
	ix = 2
      do ixy = 1, 2
      if( ixy .eq. iy ) then
	k1 = ky1
	k2 = ky2
      else
	k1 = kx1
	k2 = kx2
      end if
      do k = k1, k2
	if( k .eq. 21 ) kolor = k21
	if( k .eq. 22 ) kolor = k22
	if( k .eq. 23 ) kolor = k23
	if( k .eq. 24 ) kolor = k24
	if( k .eq. 25 ) kolor = k25
	if( k .eq. 26 ) kolor = k26
	if( k .eq. 31 ) kolor = k31
	if( k .eq. 32 ) kolor = k32
      do i = 1, ns
	if( k .eq. 22 ) plot21( i ) = plot22( i )
	if( k .eq. 23 ) plot21( i ) = plot23( i )
	if( k .eq. 24 ) plot21( i ) = plot24( i )
	if( k .eq. 25 ) plot21( i ) = plot25( i )
	if( k .eq. 26 ) plot21( i ) = plot26( i )
	if( k .eq. 31 ) plot21( i ) = plot31( i )
	if( k .eq. 32 ) plot21( i ) = plot32( i )
      end do
	call avgrms( ns, plot21, avg, rms )
	write(*,*) 'k', k, ' ns', ns, ' avg', avg, ' rms', rms
* histogram bins
	n = 300
	n = 200
      if( ixy .eq. iy ) then
	p1 = ymn
	p2 = ymx
      else
	p1 = xmn
	p2 = xmx
      end if
	if( ipan .eq. ipanh ) p1 = min( p1, -p2 )
      do i = 1, n
	p = ( i - 1. ) / max( 1, n - 1 )
	p = p1 * ( 1. - p ) + p * p2
	plot20( i ) = p
	plot19( i ) = 0.
      end do
      do i = 1, ns
	val = plot21( i )
	call lookup( n, plot20, val, j1, j2, p )
	a = j1 * ( 1. - p ) + p * j2
	j = max( 1, min( n, nint( a ) ) )
	plot19( j ) = plot19( j ) + 1
      end do
* min and max values (use pmx to scale histogram on plot )
	call mnmx( n, plot19, pmn, pmx )
	write(*,*) 'Histogram counts', pmn, pmx
* min and max values (omit first and last to avoid pile-up )
	j = min( n, 2 )
	i = max( 1, n-j )
	call mnmx( i, plot19( j ), pmn, pmx )
	write(*,*) 'Histogram counts', pmn, pmx
	p = pmx * 20.
	p = pmx * 15.
	call pgsci( kolor )
      if( ixy .eq. iy ) then
	call pgswin( p, 0., ymn, ymx )
	call pgbin( n, plot19, plot20, .true. )
      else
	call pgswin( xmn, xmx, 0., p )
	call pgbin( n, plot20, plot19, .true. )
      end if
	call pgswin( xmn, xmx, ymn, ymx )
	call pgsci( kblack )
* next k
      end do
* next ixy
      end do
      end if
      end if

* extra T(r) curves:
      if( ipan .eq. ipant ) then

* power-law T(r)
	call pgsci( kfnt )
	call pgsls( ldot )
	call pgline( np, plot1, plot4 )
	call pgsls( 1 )

* filled circle at T(Rout) and T(max)
	call pgsci( ktvmin )
	call pgpoint( 1, plot1( np ), plot4( np ), ifilledcircle )

	call imnmx( np, plot3, imn, imx )
	call pgsci( ktvmax )
	call pgpoint( 1, plot1( imx ), plot3( imx ), ifilledcircle )

* Tx(r)
	call pgsci( kblue )
	call pgsls( ldash )
	call pgline( np, plot1, plot5 )
* Tq(r)
      IF( .FALSE. ) THEN
      if( qin .gt. 0. ) then
	call pgsci( kgreen )
	call pgsls( ldot )
	call pgline( np, plot1, plot6 )
      end if
      END IF
* Tc(r)
      if( fcol .ne. 1. ) then
	call pgsci( kgreen )
	call pgsls( ldash )
	call pgline( np, plot1, plot7 )
      end if

	call pgsci( kblack )
	call pgsls( 1 )
      end if

* faint-state fiducial curve
      if( ipan .ne. ipanpsib ) then
	call pgsci( kred )
	call pgsls( 1 )
	call pgslw( lwider )
	call pgline( np, plot1, plot3 )

* no-rim lag spectrum
      if( ipan .eq. ipantau ) then
	call pgsci( krim )
	call pgslw( lwider )
	call pgline( np, plot1, plot4 )
* rim lag
	call pgsci( krim )
	call pgslw( lwider )
	call pgline( np, plot1, plot5 )
	call pgsls( 1 )
	call pgslw( lwide )
      end if

* reflect across symmetry axis
      if( ipan .eq. ipanh ) then
* KDH : FASTER (BUT TRICKY) WITH PGSWIN
	x = ( xmx + xmn )
	call pgswin( xmx-x, xmn-x, ymn, ymx )
	call pgline( np, plot1, plot3 )
	call pgswin( xmn, xmx, ymn, ymx )
      end if

	call pgslw( lwide )

* filled dot at T(max)
      if( ipan .eq. ipant ) then
	call imnmx( np, plot3, imn, imx )
	call pgsci( ktvmax )
	call pgpoint( 1, plot1( imx ), plot3( imx ), ifilledcircle )


* histogram and dots (TO TEST ADEQUACY OF RADIAL SAMPLING)
      if( .not. simpleplot ) then
	call pgsci( kfnt )
	call pgbin( np, plot1, plot3, .true. )
* call pgpoint( np, plot1, plot3, 20 )
	do i = 1, np
     call pgpoint( 1, plot1(i), plot3(i), 20 )
	end do

      end if

      end if

* rim and no-rim sed
      if( ipan .eq. ipanf ) then
	call pgslw( lwider )
	call pgsci( krim )
	call pgline( np, plot1, plot4 )
	call pgsci( krim )
	call pgline( np, plot1, plot5 )
	call pgsls( 1 )
	call pgslw( lwide )

* shift log flux spectrum up and down by dexsed
      if( logy .ne. 0 .and. dexsed .ne. 0 ) then
	add = dexsed
	call pgslw( lwide )
	call pgsci( kfnt )
	call pgsls( ldash )
	call pgswin( xmn, xmx, ymn + add, ymx + add )
	call pgline( np, plot1, plot3 )
	call pgswin( xmn, xmx, ymn - add, ymx - add )
	call pgline( np, plot1, plot3 )
	call pgswin( xmn, xmx, ymn, ymx )
      end if
	call pgsls( 1 )
	call pgsci( kfnt )
	call pgslw( lwider )
	call pgline( np, plot1, plot3 )
	call pgslw( lwide )

      end if


* mark mean delay
      if( ipan .eq. ipanpsi ) then
	x = tauavg0 * u
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
      end if
	call pgsls( 1 )
	call pgsci( kblack )

* end faint-state fiducial curve
      end if

* multi-band delay maps
      if( ipan .eq. ipanpsib ) then

* aim for 5 cases to avoid clutter
	ibstep = 1 + (nb-3)/4
	j = 0

* nb = 10, j = 0 1 2 3 4
*	1  1+9/4 1+18/4 1+27/4 1+36/4
*       1  3.25  5.5    7.75   10
*       1  3     5      8      10
*     1000 1500  3000   6000   10000

      do ib = 1, nb

	iwant = nint( 1 + ( nb - 1 ) * j / 4. - 0.0001 )
      if( ib .eq. iwant ) then

c      if( ib .eq. 1 .or. ib .eq. nb
c     &	.or. mod( ib-1, ibstep ) .eq. 0 ) then

* k=0 for fiducial model
	kkfnt = 0
	kkbrt = 1
	kkdsk = 2
	kkrim = 3
      do kk = 0, 2
	n = 0
      do i = 1, ntau
	tau = tau_d( i )
	if( tau .le. 0. ) cycle
	if( kk .eq. kkbrt ) psi =  psi_db( i, ib )
	if( kk .eq. kkfnt ) psi = psi0_db( i, ib )
	if( kk .eq. kkdsk ) psi = psi1_db( i, ib )
	if( kk .eq. kkrim ) psi = psi2_db( i, ib )
	if( psi .le. 0. ) cycle
	n = n + 1
	plot1( n ) = alog10( tau * u )
	plot2( n ) = alog10( psi / u )
      end do
      if( kk .eq. kkbrt ) then
	tau = tau_b( ib )
	call pgsls( 1 )
      else if( kk .eq. kkfnt ) then
	tau = tau0_b( ib )
	call pgsls( ldashdot )
      else if( kk .eq. kkdsk ) then
	tau = tau1_b( ib )
	call pgsls( ldashdot )
      else if( kk .eq. kkrim ) then
	tau = tau2_b( ib )
	call pgsls( ldashdot )
      end if
	call pgsci( kolor_b( ib ) )
	call pgline( n, plot1, plot2 )
	x = alog10( tau * u )
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
	call pgsls( 1 )
* next kk = 0, 2
      end do
* wavelength
	w = wobs_b( ib )
	write( label, * ) nint( w )
	call nospace( label )
	x = 0.95
	j = j + 1
	y = -0.5 -1.5 * j
	call pgmtxt( 'T', y, x, 1., label )
	write(*,*) ib, j, nint( w )
* next band ib
      end if
      end do
	call pgsci( kblack )
	call pgmtxt( 'T', 0.5, x, 1., '\gl/\A' )

	np = 0
      end if

* main curve
c	call pgsci( kblue )
	call pgsci( kbrt )
	call pgslw( lwider )
	call pgline( np, plot1, plot2 )
	call pgslw( lwide )

* filled circle at T(max) T(min) and T(out)
      if( ipan .eq. ipant ) then
	call imnmx( np, plot2, imn, imx )
	call pgsci( ktmax )
	call pgpoint( 1, plot1( imx ), plot2( imx ), ifilledcircle )
	call pgsci( ktmin )
	call pgpoint( 1, plot1( imn ), plot2( imn ), ifilledcircle )
	call pgsci( ktout )
	call pgpoint( 1, plot1(  np ), plot2(  np ), ifilledcircle )
	call pgsci( kbrt )

* small dot for testing
      if( .not. simpleplot ) then
	call pgsci( kbrt )
	call pgbin( np, plot1, plot2, .true. )
* call pgpoint( np, plot1, plot2, 20 )
	do i = 1, np
     call pgpoint( 1, plot1(i), plot2(i), 20 )
	end do

      end if
* end if ipant
      end if

* shift up and down by dexsed
      if( logy .ne. 0 .and. ipan .eq. ipanf .and. dexsed .ne. 0. ) then
	add = dexsed
c	call pgslw( lwider )
	call pgsci( kbrt )
	call pgsls( ldash )
	call pgswin( xmn, xmx, ymn + add, ymx + add )
	call pgline( np, plot1, plot2 )
	call pgswin( xmn, xmx, ymn - add, ymx - add )
	call pgline( np, plot1, plot2 )
	call pgswin( xmn, xmx, ymn, ymx )
	call pgsls( 1 )
	call pgsci( kbrt )
	call pgline( np, plot1, plot2 )
	call pgslw( lwide )
      end if

* shift lag spectrum up and down by syslag
      if( ipan .eq. ipantau .and. logy .eq. 0 .and. syslag .gt. 0. ) then
	add = syslag
c	call pgslw( lwider )
	call pgsci( kbrt )
	call pgsls( ldash )
	call pgswin( xmn, xmx, ymn + add, ymx + add )
	call pgline( np, plot1, plot2 )
	call pgswin( xmn, xmx, ymn - add, ymx - add )
	call pgline( np, plot1, plot2 )
	call pgswin( xmn, xmx, ymn, ymx )
	call pgsci( kbrt )
	call pgsls( 1 )
	call pgline( np, plot1, plot2 )
	call pgslw( lwide )
      end if

	kolor = kblue

* delays for each band

      if( ipan .eq. ipantau ) then
      if( .not. simpleplot ) then
	call pgqci( kolor )
	call pgsch( csize * 1.5 )
      do ib = 1, nb
	x = plot1( ib )
	y = plot2( ib )
	call pgsci( kolor_b( ib ) )
	call pgpoint( 1, x, y, ifilledcircle )
	y = plot3( ib )
	call pgpoint( 1, x, y, iopencircle )
      end do
	call pgsch( csize )
	call pgsci( kolor )
      end if

* lag data
      if( nlag .gt. 0 ) then
* reference lag
	opz = 1. + redshift
	w = wlag0 / opz
	call lookup( nb, wobs_b, w, i1, i2, p )
	reflag = terp1( nb, tau_b, i1, i2, p )
	write(*,*) 'tau(', nint( w ), ' ) =', reflag

* mark reference wavelength
	x = w
	if( logx .ne. 0 ) x = alog10( x )
c	call pgsci( kbrt )
	call pgsci( kblack )
	call pgsls( ldash )
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
* mark reference lag
	y = reflag
	call pgmove( 0., y )
	call pgdraw( x, y )
	call pgsls( 1 )

* reference wavelength
	label = '\gl\d0\u='
      if( w .gt. 10 ) then
	call append_i( label, nint( w ) )
      else
	call append_n( label, w, 3 )
      end if
	call append( label, '\A' )
	call nospace( label )
	x = 0.05
	call pgmtxt( 'T', -1.5, x, 0.0, label )
* reference lag
	label = '\gt\d0\u='
	call append_n( label, reflag, 3 )
	call append( label, 'd' )
	call nospace( label )
	call pgmtxt( 'T', -3.0, x, 0.0, label )

* lag systematic error (only if large enough)
c      if( .not. simpleplot .and. syslag .gt. 1.d-4 ) then
* threshold is smallest lag error bar / 10
	call mnmx( nlag, siglaglo_b, e1, e2 )
	call mnmx( nlag, siglaghi_b, e3, e4 )
	e = min( e1, e2, e3, e4 ) / 10.
      if( syslag .lt. e ) then
	write(*,*) 'syslag', syslag, ' <', e
      else
	write(*,*) 'syslag', syslag, ' >', e
	label = '\gs\d\gt\u='
	call append_n( label, syslag, 3 )
      if( index( label, 'E' ) .le. 0 ) then
	call tidy( label )
      end if
	call append( label, 'd' )
	call nospace( label )
* top right
c	call pgmtxt( 'T', -1.5, 0.95, 1., label )
	call pgmtxt( 'T', -4.5, x, 0., label )
      end if

	call pgsci( kblack )


* lags
	opz = 1. + redshift
      do i = 1, nlag
	x =   wlag_b( i ) / opz
	xp = wlagfwhm_b( i ) / opz / 2.
	xm = xp
	dat   =   datlag_b( i ) + reflag
	sighi = siglaghi_b( i )
	siglo = siglaglo_b( i )
	y = dat * u
	yp = sighi * u
	ym = siglo * u
	write(*,*) nint(x), nint( xp ), y, yp, ym
* log x
      if( logx .ne. 0 ) then
	xp = x + xp
	xm = x - xm
	x  = alog10( x  )
	xp = alog10( xp ) - x
	xm = x - alog10( xm )
      end if
* log y
      if( logy .ne. 0 ) then
	floor = 10. ** ( ymn - 0.5 )
	yp = y + yp
	ym = y - ym
	y  = alog10( max( floor, y  ) )
	yp = alog10( max( floor, yp ) ) - y
	ym = y - alog10( max( floor, ym ) )
      end if
* error bars
	call pgsci( kblack )
	call pgerrb( 1, 1, x, y, xp, 1. )
	call pgerrb( 2, 1, x, y, yp, 1. )
	call pgerrb( 3, 1, x, y, xm, 1. )
	call pgerrb( 4, 1, x, y, ym, 1. )
	call pgsci( kblack )
* black star
	call pgpoint( 1, x, y, istar )

c	write(*,*) '  x,+,-', x, xp, xm
c	write(*,*)  ' y,+,-', y, yp, ym

* increment lag chi2
	call lookup( np, plot1, x, i1, i2, p )
	fit = terp1( np, plot2,    i1, i2, p ) / u
	sig = sighi
	if( dat .lt. fit ) sig = siglo
	chi = ( dat - fit ) / sig
	if( i .eq. 1 ) chi2lag = 0.d0
	chi2lag = chi2lag + chi * chi

* next lag
      end do

* chi2 label
      if( .not. simpleplot .or. .true. ) then
	label = '\gx\u2\d='
	c = chi2lag
	x = c / max( 1, nlag )
	call append_f( label, c, '(f15.1)' )
	call append( label, '/' )
	call append_i( label, nlag )
	call append( label, '=' )
	call append_f( label, x, '(f15.2)' )
	call nospace( label )
	call pgmtxt( 't', 0.5, 1., 1., label )
      end if

* end nlag > 0
      end if
* end lag panel
      end if

* flux data (k=0,1 for faint,bright)
      if( ipan .eq. ipanf ) then
      if( ndat .gt. 0 ) then
* 2 cases (bright and faint)
      do k = 0, 1
* data
	opz = 1. + redshift
      do i = 1, ndat
* wavelength
        x = wdat_b( i )
	dx = fwhm_b( i ) / 2.
	x = x / opz
	dx = dx / opz
* flux
      if( k .eq. 1 ) then
	dat = datflx_b( i )
	sig = sigflx_b( i )
      else
	dat = datflx0_b( i )
	sig = sigflx0_b( i )
      end if
	y = dat
	yp = sig
	ym = yp
c        write(*,*) k, nint(x), y, yp

* log x
      if( logx .ne. 0 ) then
	xp = x + dx
	xm = x - dx
        x = alog10( x )
	xp = alog10( xp ) - x
	xm = x - alog10( xm )
      end if
* log y
      if( logy .ne. 0 ) then
        floor = 10. ** ( ymn - 0.5 )
        yp = y + yp
        ym = y - ym
        y  = alog10( max( floor, y ) )
        yp = alog10( max( floor, yp ) ) - y
        ym = y - alog10( max( floor, ym ) )
      end if
        call pgsci( kblack )
        call pgerrb( 1, 1, x, y, xp, 1. )
        call pgerrb( 2, 1, x, y, yp, 1. )
        call pgerrb( 3, 1, x, y, xm, 1. )
        call pgerrb( 4, 1, x, y, ym, 1. )
        call pgsci( kblack )
        call pgpoint( 1, x, y, istar )
c        write(*,*)  '  x,y,+,-', x, y, yp, ym

* next datum i
      end do

* next case k
      end do

* systematic error label
      if( simpleplot ) then
	label = '\gs\d0\u='
	p = 100. * ( 10. ** dexsed - 1. )
      if( p .ge. 30. ) then
	call append_n( label, dexsed, 2 )
	call nospace( label )
	call append( label, 'dex' )
      else
      if( p .ge. 10. ) then
	call append_i( label, nint( p ) )
      else
	call append_n( label, p, 2 )
      end if
	call append( label, '%' )
	call nospace( label )
      end if
	call pgmtxt( 't', -1.5, 0.95, 1., label )

      else
	dy = -1.5
	y = -1.0 - dy
      if( dexsed .ne. 0. ) then
	label = '\gSln\gs\u2\d='
	call append_f( label, real(sumlnvar), '(f15.1)' )
	call nospace( label )
	call word1( label, l1, l2 )
	call append( label, '\gs\d0\u=' )
	call append_f( label, dexsed, '(f15.3)' )
	call nospace( label( l2+2:) )
	call append( label, 'dex' )
	y = y + dy
	call pgmtxt( 'b', y, 0.5, 0.5, label )
      end if


      end if

* chi2 label
      if( .not. simpleplot ) then
	kkfnt = 1
	kkbrt = 2
	kkbof = 3
      do kk = 1, 3
	if( kk .eq. kkbrt ) c = chi2brt
	if( kk .eq. kkfnt ) c = chi2fnt
	if( kk .eq. kkbof ) c = chi2brt + chi2fnt + sumlnvar
	n = ndat
* double n since n bright + n faint data
	if( kk .eq. kkbof ) n = n + n
	x = c / max( 1, n )
	label = '\gx\u2\d='
	if( kk .eq. kkbof .and. sumlnvar .ne. 0.d0 ) label = 'BoF='
	call append_f( label, c, '(f15.1)' )
	call append( label, '/' )
	call append_i( label, n )
	call append( label, '=' )
	call append_f( label, x, '(f15.2)' )
	call nospace( label )
	y = y + dy
      if( kk .eq. kkbrt ) then
	call append( label, 'Bright' )
	call pgmtxt( 'b', y, 0.5, 0.5, label )
      else if( kk .eq. kkfnt ) then
	call append( label, 'Faint' )
	call pgmtxt( 'b', y, 0.5, 0.5, label )
      else if( kk .eq. kkbof ) then
	call pgmtxt('t', 0.5, 1.0, 1.0, label )
      end if
* next kk = 1, 3
      end do
      end if

* end flux data ndat > 0
      end if
* end flux panel ipan = ipanf
      end if

* mark radii for each band
      if( ipan .eq. ipant ) then
      if( .not. simpleplot ) then
	call pgqci( kolor )
	call pgsch( csize * 1.5 )
      do ib = 1, nb
* first exposed annulus outside this lag
	tau = tau_b( ib )
	call lookup( nr, tau_r, tau, i1, i2, p )
      if( i_r( i1 ) .le. 0 ) then
      do i = i1, nr
	k = i
	if( i_r( i ) .gt. 0 ) exit
      end do
	i1 = k
	i2 = k
      end if
	x = terp1( nr, plot1, i1, i2, p )
	y = terp1( nr, plot2, i1, i2, p )
      if( .not. simpleplot ) then
	call pgsci( kolor_b( ib ) )
	call pgpoint( 1, x, y, ifilledcircle )
	y = terp1( nr, plot3, i1, i2, p )
	call pgpoint( 1, x, y, iopencircle )
      end if
      end do
	call pgsch( csize )
	call pgsci( kolor )
      end if
      end if

* mark flux for each band
      if( ipan .eq. ipanf ) then
      if( .not. simpleplot ) then
	call pgsch( csize * 1.5 )
      do ib = 1, nb
	w = wobs_b( ib )
	if( logx .ne. 0 ) w = alog10( w )
	call lookup( np, plot1, w, i1, i2, p )
	x = terp1( np, plot1, i1, i2, p )
	y = terp1( np, plot2, i1, i2, p )
	call pgsci( kolor_b( ib ) )
	call pgpoint( 1, x, y, ifilledcircle )
	y = terp1( np, plot3, i1, i2, p )
	call pgpoint( 1, x, y, iopencircle )
      end do
	call pgsch( csize )
      end if
      end if
	call pgsci( kblack )

* reflect disk profile thru symmetry axes
      if( ipan .eq. ipanh .and. logxy .eq. 0 ) then
	call pgslw( lwider )
	x = ( xmx + xmn )
      do k = 1, 3
	if( k .eq. 3 ) call pgswin(   xmn,   xmx, -ymn, -ymx )
	if( k .eq. 2 ) call pgswin( xmx-x, xmn-x, -ymn, -ymx )
	if( k .eq. 1 ) call pgswin( xmx-x, xmn-x,  ymn,  ymx )
	call pgline( np, plot1, plot2 )
      end do
	call pgswin( xmn, xmx, ymn, ymx )
	call pgslw( lwide )

* highlight irradiated crests
	call pgsci( kblack )
	call pgslw( lwider )
* reflect plot window across x=0
	x = ( xmx + xmn )
	call pgswin( xmx-x, xmn-x, ymn, ymx )
      do k = 0, 1
	i = 1
      do j = 1, np
* edge of a shadow
      if( i_r( i ) .ne. i_r( j ) ) then
* inner edge of exposed region
      if( i_r( i ) .eq. 0 ) then
* plot inner (hotter) half segment in blue
	jm = max( 1, j - 1 )
	x = ( plot1( jm ) + plot1( j ) ) / 2.
	y = ( plot2( jm ) + plot2( j ) ) / 2.
	call pgsci( kblue )
	call pgmove( x, y )
	call pgdraw( plot1( j ), plot2( j ) )
	call pgsci( kblack )
* outer edge
      else if( i_r( i ) .eq. 1 ) then
* plot outer (cooler) half segment in red
	jm = max( 1, j - 1 )
	x = ( plot1( jm ) + plot1( j ) ) / 2.
	y = ( plot2( jm ) + plot2( j ) ) / 2.
	call pgsci( kred )
	call pgmove( x, y )
	call pgdraw( plot1( jm ), plot2( jm ) )
	call pgsci( kblack )

	n = j - i
	call pgline( n , plot1( i ), plot2( i ) )
cc	if( n .eq. 1 ) call pgpoint( n, plot1( i ), plot2( i ), ifilledcircle )

      end if
	i = j
      end if

* next j
      end do

* final segment (if needed)
	n = np - i + 1
      if( i_r( i ) .eq. 1 .and. n .gt. 0 ) then
	call pgline( n , plot1( i ), plot2( i ) )
c	if( n .eq. 1 ) call pgpoint( n, plot1( i ), plot2( i ), ifilledcircle )
      end if

* mark all exposed pixels
      do i = 1, np
      if( i_r( i ) .gt. 0 ) then
	call pgline( 1, plot1( i ), plot2( i ) )
      end if
      end do

* restore plot window for next pass
	call pgswin( xmn, xmx, ymn, ymx )
* next k
      end do

	call pgslw( lwide )

* photon trajectories
	p = max( xmx - xmn, ymx - ymn )
	xe = p * sini
	ye = p * cosi
      do ib = 1, nb
	kolor = kolor_b( ib )
	call pgsci( kolor )
	call pgsls( ldot )
	call pgsls( 1 )

* first exposed annulus outside this delay
	tau = tau_b( ib )
	call lookup( nr, tau_r, tau, i1, i2, p )
      if( i_r( i1 ) .le. 0 ) then
      do i = i1, nr
	k = i
	if( i_r( i ) .gt. 0 ) exit
      end do
	i1 = k
	i2 = k
      end if

	klo = 1
	khi = 2
      do k = 1, 2
	x = terp1( nr, r_r, i1, i2, p ) * u
	y = terp1( nr, h_r, i1, i2, p ) * u
	if( k .eq. khi ) x = - x
c	write(*,*) 'tau', tau, ' R', x, ' H', y
	call pgmove( 0., hlamp * u )
	call pgdraw( x, y )
	call pgsch( csize * 1.5 )
	call pgpoint( 1, x, y, ifilledcircle )
	call pgsch( csize )
	x = x + xe
	y = y + ye
	call pgdraw( x, y )
* next k
      end do
* next band ib
      end do
	call pgsls( 1 )
	call pgsci( 1 )

* inclined line of sight from lamp
	call pgsci( 2 )
	call pgsls( ldot )
	x = 0.
	y = hlamp * u
	call pgmove( x, y )
	x = x + ( ymx - y ) * sini / cosi
	y = ymx
	call pgdraw( x, y )
	call pgsls( 1 )
	call pgsci( 1 )
      if( .false. ) then
	x = ( x - xmn ) / ( xmx - xmn )
	label = 'i='
	call append_i( label, nint( dinc ) )
	call append( label, '\uo\d' )
	call nospace( label )
	call pgmtxt( 'T', 0.5, x, 0.5, label )
      end if

* lamp post (star)
	call pgsci( kblack )
	call pgsch( 2 * csize )
	y = hlamp * u
	call pgpoint( 1, 0., y, istar )
	call pgpoint( 1, 0., -y, istar )
	call pgsch( csize )

      end if

* mean lag
      if( ipan .eq. ipanpsi ) then
	x = tauavg * u
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
	call pgsls( 1 )
	call pgsci( 1 )

	label = '<\gt>='
	call append_n( label, tauavg0 * u, 3 )
	call append( label, '=>' )
	call append_n( label, tauavg * u, 3 )
	call nospace( label )
	x = ( tauavg0 + tauavg ) / 2. * u
	x = ( x - xmn ) / ( xmx - xmn )
	call pgmtxt( 'T', 0.5, x, x, label )
      end if

	call pgsci( 1 )



* parameter legend ======================
* KDH : MAKE A SUBROUTINE ?
* KDH : DO THESE LAST

      if( ipan .eq. ipanh ) then

* lamp height
	label = 'H\dx\u='
	call append_n( label, hlamp_rg, 2 )
	call tidy( label )
	call append( label, 'R\dg\u' )
	call nospace( label )
	x = ( 0. - xmn ) / ( xmx - xmn )
	call pgmtxt( 'b', -1.0, x, 0.5, label )


* outer radius
	label = 'R\dout\u='
	call append_n( label, rout, 3 )
	call tidy( label )
      if( lightsec ) then
	call append( label, 's' )
      else
	call append( label, 'd' )
      end if
	call nospace( label )
	call pgmtxt( 't', -1.5, 0.95, 1.0, label )

* H/R
	label = 'H/R='
	hout = h1
	if( r1 .ne. 0. ) hout = hout * ( rout / r1 ) ** bet
	pct = 100. * hout / rout
* KDH : USE h1_r1_pct ?
      if( pct .lt. 1.e-3 ) then
	call append_i( label, 0 )
	call nospace( label )
      else
	call append_n( label, pct, 3 )
	call tidy( label )
	call append( label, '%' )
	call nospace( label )
	call pgmtxt( 't', -3.0, 0.95, 1.0, label )

* beta = dlnH/dlnR
	label = '\gb='
	i = 1
c	i = len1( label ) + 2
c	call append( label, '\gb=' )
      if( bet .gt. 1e4 ) then
	label = 'log\gb='
	call append_f( label, alog10( bet ), '(f10.2)' )
      else if( bet .ge. 100. ) then
	call append_i( label, nint( bet ) )
      else if( bet .gt. 1.e-3 ) then
	call append_n( label, bet, 3 )
	call noeqbl( label )
	call tidy( label )
      else
	call append_i( label, 0 )
      end if
	call nospace( label(i:) )
      end if
	call pgmtxt( 't', -4.5, 0.95, 1.0, label )
* top left
c	call pgmtxt( 't', 0.5, 0., 0., label )


* inclination
	label = 'i='
	call append_i( label, nint( dinc ) )
	call append( label, '\uo\d' )
	call nospace( label )
c	call pgmtxt( 't', 0.5, 1., 1., label )
	call pgmtxt( 't', 0.5, 0.5, 0.5, label )

      else if( ipan .eq. ipant ) then
* Mbh
	x = bhM_Msun
      IF( .FALSE. ) THEN
	label = 'logM='
	x = alog10( bhM_Msun )
	call append_f( label, x, '(f5.2)' )
      ELSE
	label = 'M='
	i = int( alog10( x ) )
	if( i .le. 0 ) i = i - 1
      if( iabs( i ) .le. 3 ) then
	call append_n( label, x, 3 )
      else
	x = x / 10. ** i
	call append_n( label, x, 3 )
	call tidy( label )
	call append( label, 'E' )
	call append_i( label, i )
      end if
      END IF
	call nospace( label )
	call append( label, 'M\d\(2281)\u' )
	call pgmtxt( 't', 0.5, 0., 0., label )

* DlogT
      IF( .TRUE. ) THEN
	label = '\gDT='
	call append_n( label, dlogt, 3 )
	call noeqbl( label )
	call append( label, 'dex' )
	call nospace( label )
	call pgsci( kblue )
C	call pgmtxt( 't', -1.5, 0.95, 1., label )
	x = 0.7
c	call pgmtxt( 't', -1.5, 0.7, 0.7, label )
	call pgmtxt( 't', -1.5, 0.55, 0.0, label )
	call pgsci( kblack )

* colour correction
      if( fcol .gt. 0. .and. fcol .ne. 1. ) then
	label = '\gDT\dcol\u='
	dlogtc = alog10( fcol )
	call append_n( label, dlogtc, 3 )
	call noeqbl( label )
	call append( label, 'dex' )
	call nospace( label )
	call pgsci( kgreen )
      if( simpleplot ) then
c	call pgmtxt( 't', -3.0, 0.95, 1., label )
	call pgmtxt( 't', -1.5, 0.5, 1., label )
      else
	call pgmtxt( 't', -1.5, 0.5, 1., label )
      end if
	call pgsci( kblack )
      end if

* Delta(LxHx)
      ELSE IF( .FALSE. ) THEN
	label = '\gD(L\dx\uH\dx\u'
	if( albedo .ne. 0 ) call append( label, '(1-A)' )
	call append( label, ')=' )
	call nospace( label )
	call append_n( label, dkappa, 3 )
	call noeqbl( label )
	call tidy( label )
	l = len1( label ) + 1
	call pgmdot( label(l:) )
	l = len1( label )
	call append( label, 'c\u2\dR\dg\u' )
	call nospace( label(l:) )
	call pgmtxt( 't', -1.5, 0.95, 1., label )
      END IF

* Lx/Mdotc^2
	label = '\gDL\dx\u='
	call append_n( label, xlamp, 3 )
	call nospace( label )
	call tidy( label )

	call pgmdot( label(len1(label)+1:) )
	l = len1( label )
	call append( label, 'c\u2\d' )
	call nospace( label(l:) )
c	call pgmtxt( 't', -3.0, 0.95, 1., label )
	call pgmtxt( 't', -3.0, 0.55, 0., label )


* Hx/Rg
      if( .not. simpleplot ) then
	label = 'H\dx\u='
	call append_n( label, hlamp_rg , 3 )
	call noeqbl( label )
	call tidy( label )
	call append( label, 'R\dg\u' )
	call nospace( label )
	call pgmtxt( 't', -4.5, 0.95, 1., label )
      end if

* L/Ledd
	label = 'L/L\dEdd\u='
	x = bhL_Ledd
	yes = x .gt. 0.5
      if( yes ) then
	call append_n( label, x, 3 )
      else
	call append_n( label, 100 * x, 3 )
      end if
	call nospace( label )
	call tidy( label )
	if( .not. yes ) call append( label, '%' )
* eta
      if( eta .ne. 1. ) then
	call append( label, '(\gy/' )
	call append_n( label, eta, 3 )
	call nospace( label )
	call tidy( label )
	call append( label, ')' )
      end if
	call nospace( label )
	call pgmtxt( 't', 0.5, 1., 1., label )

* Rin/Rg
      IF( .FALSE. ) THEN
	label = 'R\din\u='
	call append_n( label, risco_rg, 3 )
	call nospace( label )
	call tidy( label )
	call append( label, 'R\dg\u' )
	call nospace( label )
* eta
	call word1( label, l1, l2 )
	call append( label, ' (\gy=' )
	call append_n( label, 0.5/risco_rg, 3 )
	call tidy( label )
	call append( label, ')' )
	call nospace( label(l2+2:) )

* R=(Rin,Rout)
      ELSE
	label = 'R=('
	call append_n( label, risco_rg, 3 )
	call tidy( label )
	call append( label, ',' )
	i = len1( label )
	call append_n( label, rout_rg, 3 )
	call tidy( label )
	call append( label, ')' )
	call nospace( label )
	call append( label, 'R\dg\u' )
	call nospace( label )
      END IF
	x = 0.15
	call pgmtxt( 'b', -2.5, x, 0., label )

* Mdot
	call pgmdot( label )
	call append( label, '<' )
	x = bhMsun_yr
	i = int( alog10( x ) )
	if( i .le. 0 ) i = i - 1
      if( iabs( i ) .le. 3 ) then
	call append_n( label, x, 3 )
	call tidy( label )
      else
	x = x / 10. ** i
	call append_n( label, x, 3 )
	call tidy( label )
	call append( label, 'E' )
	call append_i( label, i )
      end if
	call append( label, 'M\d\(2281)\u/yr' )
	call nospace( label )

* torque parameter
      if( .not. simpleplot .and. qin .gt. 1.e-5 ) then
	call append( label, ' Q=' )
	i = len1( label )
	p = 1.
	if( qin .gt. 0. .and. qin .lt. 1. ) p = 100.
	call append_n( label, qin * p, 3 )
	call nospace( label(i:) )
	call tidy( label(i:) )
	if( p .gt. 1. ) call append( label, '%' )
	call nospace( label(i:) )
      end if

	x = 0.15
	call pgmtxt( 'b', -1.0, x, 0.0, label )

* temperature labels
	call imnmx( nr, t_r, imn, imx )

	i1 = 1
	if( imn .lt. nr .and. simpleplot ) i1 = 0
      do i = i1, 4
* T(min)
	if( i .eq. 0 ) t = t_r( imn )
* T(max)
	if( i .eq. 1 ) t = t_r( imx )
* T(out)
	if( i .eq. 2 ) t = t_r( nr )
* Tv(min)
	if( i .eq. 3 ) t = tvmn
* Tq(max)
	if( i .eq. 4 ) t = tqmx
* temperature label
      if( t .gt. 0. ) then
	label = ' '
	k = int( alog10( t ) )
      if( k .lt. 5 ) then
	call append_i( label, nint( t ) )
      else
	call append_n( label, t / 10. ** k, 3 )
	call append( label, 'E' )
	call append_i( label, k )
      end if
	call append( label, 'K' )
	call nospace( label )
	x = 0.15
	xx = 0.90
	y = -6.0
c	if( simpleplot ) y = -5.0
* T(max)
	call pgsci( kblue )
      if( i .eq. 1 ) then
	call pgmtxt( 'b', -5.5, x , 0.0, label )
      end if
* T(min)
      if( i .eq. 0 ) then
	call pgsci( korange )
	call pgmtxt( 't', y, xx, 1.0, label )
      end if
* T(out)
      if( i .eq. 2 ) then
	call pgsci( kgreen )
	if( i1 .eq. 0 ) y = y + 1.5
	call pgmtxt( 't', y, xx, 1.0, label )
      end if
* Tv(min)
	call pgsci( kred )
	if( i .eq. 3 ) call pgmtxt( 'b', -1.0, xx, 1.0, label )
* Tq(max)
	if( i .eq. 4 ) call pgmtxt( 'b', -4.0, x , 0.0, label )

	call pgsci( kblack )
* next temp
      end if
      end do

      else if( ipan .eq. ipanf ) then
	label = 'z='
	call append_n( label, redshift, 3 )
	call nospace( label )
	call pgmtxt( 't', 0.5, 0., 0., label )

	label = 'D\dL\u='
	call append_n( label, dmpc, 3 )
	call nospace( label )
	call append( label, 'Mpc' )
	call nospace( label )
	x = 0.05
	call pgmtxt( 't', -1.5, x, 0., label )

      if( ebmv_agn .ne. 0. ) then
	label = 'E\uSMC\b\b\b\d\dB-V\u='
	label = 'E\dB-V\u='
	call append_n( label, ebmv_agn, 3 )
	call nospace( label )
	x = 0.05
	call pgmtxt( 't', -3.0, x, 0., label )
      end if

* cosmic variance uncertainty
      if( .not. simpleplot ) then
      if( dmpc .eq. dmpcz ) then
	e = 300. / h0
	e = 500. / h0
	e = 2. * alog10 ( 1. + e / dmpc )
	x = xmn + 0.9 * ( xmx - xmn )
	y = ymn + 0.9 * ( ymx - ymn )
	call pgsci( 2 )
	write(*,*) 'Cosmic uncertainty', e, ' (dex)'
	write(*,*) 'PGERRB(', x, y, e, ' )'
	call pgerrb( 6, 1, x, y, e, 1. )
	call pgsci( 1 )
      end if
      end if
      end if

* re-draw the box
	logxy = logx * 10 + logy * 20
	if( logx .eq. 0 ) call pgbox( 'BCST',  0., 0, ' ', 0., 0 )
	if( logx .ne. 0 ) call pgbox( 'BCSTL', 0., 0, ' ', 0., 0 )
	if( logy .eq. 0 ) call pgbox( ' ', 0., 0, 'BCST',  0., 0 )
	if( logy .ne. 0 ) call pgbox( ' ', 0., 0, 'BCSTL', 0., 0 )

* writeout
      if( writeout ) then
      if( ipan .eq. ipanf ) then
	file = 'bowl_sed.dat'
	call inq( 'SED output file:', file, file )
	call dump1d5( file, 5, np, plot1, plot2, plot3, plot4, plot5 )
  		end if
  		if( ipan .eq. ipantau ) then
	file = 'bowl_lag.dat'
	call inq( 'LAG output file:', file, file )
	call dump1d3( file, 3, np, plot1, plot2, plot3 )
      end if
			if( ipan .eq. ipanh ) then
	file = 'bowl_aspect.dat'
	call inq( 'ASPECT output file:', file, file )
	call dump1d3( file, 3, np, plot1, plot2, plot3 )
      end if
			if( ipan .eq. ipant ) then
	file = 'bowl_temp.dat'
	call inq( 'ASPECT output file:', file, file )
	call dump1d7( file, 7, np, plot1, plot2, plot3, plot4, plot5, plot6, plot7 )
			end if
      end if

* next panel ipan
      end do
* close
	call pgstop
* again
	if( pgagain() ) goto 100

* restore previous parameters
      if( nstow .gt. 0 .and. nfit .gt. 0 .and. nits .gt. 0 ) then
	write(*,*) 'Restore saved parameters.  BoF', bofsave
      do i = 1, npar
	old = p_p( i )
	p = save_p( i )
	name = name_p( i )
	write(*,*) i, ' ', name(:len1(name)), old, ' =>', p
	p_p( i ) = p
      end do
* extract primary parameters from p_p
	call getpar
* new BoF
	needsed = .true.
	needlag = .true.
	bof = bofcalc( p_p )
	write(*,*) 'Re-evaluated BoF', bof, ' was', bofsave

      end if

	return
	end

*================================================
	function bofcalc( params )
* compute the badness of fit for given parameters
* 2021 Nov Keith Horne @ Crail
* 2021 Nov KDH @ Crail - verbose
* 2022 Jul KDH @ StA - omit argument params (not used)
* 2023 Jun KDH @ StA - fix bug in getpar change to needsed and needlag
* 2025-Mar-07 KDH @ StA - rest-frame SMC-like dust extinction
	use all
	use constants
	use agn

	real*4 params( * )
* KDH : PARAMS SEEMS NOT TO BE USED

	logical verbose/ .false. /

	write(*,*) '--- geometry ------------ Nr', nr
	write(*,*) '--- bofcalc ----------- Npar', npar, ' Nfit', nfit

* set up the parameter list
      if( npar .le. 0 ) then
	call initpar
	call setpar
      end if

* extract nfit parameters
      if( nfit .gt. 0 ) then
* extract from p_p
	call getpar
      end if
* update seds
	call sedcalc
* update lags
	call lagcalc

* zero sums
	chi2lag = 0.d0
	chi2brt = 0.d0
	chi2fnt = 0.d0
	chi2 = 0.d0
	sumlnvar = 0.d0
	sumlnvarlag = 0.d0
	sumlnvarsed = 0.d0
	bof = 0.d0

* delays for each band

* lag data
      if( nlag .gt. 0 ) then

	logy = 0
	logx = 1
	call mnmx( nb, tau_b, tmn, tmx )
	floor = max( tmn, tmx * 1e-10 )
	if( floor .le. 0. ) floor = 0.1
* report
      if( verbose ) then
	write(*,*) 'Lag data', nlag
	write(*,*) 'Lag model wavelengths', nb
	write(*,*) 'LogX', logx, ' LogY', logy, ' floor', floor
      end if

      do i = 1, nb
	x = wobs_b( i )
	y = tau_b( i )
	if( verbose ) write(*,*) i, ' lam', nint(x), ' lag', y
	if( logx .ne. 0 ) x = alog10( x )
	if( logy .ne. 0 ) y = alog10( max( floor, y ) )
	plot1( i ) = x
	plot2( i ) = y
      end do
	np = nb

* reference lag
	opz = 1. + redshift
	w = wlag0 / opz

	call lookup( nb, wobs_b, w, i1, i2, p )
	reflag = terp1( nb, tau_b, i1, i2, p )
* report
      if( verbose ) then
	write(*,*) 'redshift', redshift
	write(*,*) 'reference wavelength', w
	write(*,*) 'tau(', nint( w ), ' ) =', reflag
	write(*,*) 'Chi2lag', chi2lag
      end if
* lags
      do i = 1, nlag
	w = wlag_b( i ) / opz
	dw = wlagfwhm_b( i ) / opz / 2.
	dat   =   datlag_b( i ) + reflag
	sighi = siglaghi_b( i )
	siglo = siglaglo_b( i )
	x =   w
	xp = dw
	xm = xp
	y = dat
	yp = sighi
	ym = siglo

* log x
      if( logx .ne. 0 ) then
	xp = x + xp
	xm = x - xm
	x  = alog10( x  )
	xp = alog10( xp ) - x
	xm = x - alog10( xm )
      end if

* log y (NOT FULLY IMPLEMENTED.  ABORT)
      if( logy .ne. 0 ) then
	write(*,*) '** ERROR : LOG(TAU) DATA NOT FULLY IMPLEMENTED.'
	stop
	floor = 10. ** ( ymn - 0.5 )
	yp = y + yp
	ym = y - ym
	y  = alog10( max( floor, y  ) )
	yp = alog10( max( floor, yp ) ) - y
	ym = y - alog10( max( floor, ym ) )
      end if

* fit interpolates model lags
	call lookup( np, plot1, x, i1, i2, p )
	fit = terp1( np, plot2,    i1, i2, p )

* extra variance
	ss = syslag * syslag
* asymmetric error bar
      if( dat .ge. fit ) then
	old = yp
	sig = sqrt( yp * yp + ss )
      else
	old = ym
	sig = sqrt( ym * ym + ss )
      end if

* increment lag chi2
	chi = ( dat - fit ) / sig
	add = chi * chi
	chi2lag = chi2lag + add

* error bar expansion penalty
      if( sig .ne. old ) then
	sumlnvarlag = sumlnvarlag + 2. * alog( sig / old )
      end if

* report
      if( verbose ) then
	write(*,*) '-------------------------------'
	write(*,*) i, 'lam', nint( w ), ' +/-', nint( dw )
	write(*,*) i, 'lag dat', dat, ' fit', fit
	write(*,*) i, 'lag sig', sig, ' chi', chi
	write(*,*) i, 'lag Chi^2', chi2lag, ' sum(ln(var))', sumlnvarlag
      end if
* next lag
      end do
* end nlag > 0
      end if

* SED data
      if( ndat .gt. 0 ) then

	logx = 1
	logy = 1

	call mnmx( nw, f_w, fmn, fmx )
	floor = max( fmn / 10. , fmx * 1e-10 )
	if( floor .le. 0. ) floor = 1.e-5
* report
      if( verbose ) then
	write(*,*) 'Wavelengths:', nw, nint(w_w(1)), nint(w_w(nw))
	write(*,*) 'LogX', logx, ' LogY', logy, ' floor', floor
	write(*,*) 'Model Fluxes', nw
      end if
      do i = 1, nw
	w  =  w_w( i )
	f  =  f_w( i )
	f0 = f0_w( i )
	f1 = f1_w( i )
	f2 = f2_w( i )
	if( logx .ne. 0 )  w = alog10( w )
	if( logy .ne. 0 )  f = alog10( max( floor, f  ) )
	if( logy .ne. 0 ) f0 = alog10( max( floor, f0 ) )
	if( logy .ne. 0 ) f2 = alog10( max( floor, f2 ) )
	plot1( i ) = w
	plot2( i ) = f
	plot3( i ) = f0
	plot4( i ) = f2
      if( verbose .and. mod( i, max(1,nw/10) ) .eq. 0 ) then
	write(*,*) i, ' lam', nint(w), ' f/mjy', f, f0, f1, f2
      end if
      end do
	np = nw

	sigsed = 10.d0 ** dexsed - 1.d0
      if( verbose ) then
	write(*,*) 'Flux data', ndat
	write(*,*) 'systematics', 100.*sigsed, ' %', dexsed, ' (dex)'
      end if

* floor
	call mnmx( ndat, datflx0_b, ymn, ymx )
	floor = ymn * 1.e-6
	if( floor .lt. 0. ) floor = 1.e-10
* bright/faint loop
	ibrt = 1
	ifnt = 2
      do ibf = 1, 2
* report
      if( verbose ) then
	if( ibf .eq. ibrt ) write(*,*) 'Bright Chi2', chi2brt
	if( ibf .eq. ifnt ) write(*,*) 'Faint Chi2', chi2fnt
      end if
      do i = 1, ndat
* wavelength and bandwidth
        x = wdat_b( i )
	dx = fwhm_b( i ) / 2
	opz = 1. + redshift
	x = x / opz
	dx = dx / opz
* report
      if( verbose ) then
	write(*,*) '-------------------------------------'
	write(*,*) i, ' lam', nint(x), ' +/-', nint( dx )
      end if
* bright flux data
      if( ibf .eq. ibrt ) then
	dat = datflx_b( i )
	sig = sigflx_b( i )
* faint flux data
      else if( ibf .eq. ifnt ) then
	dat = datflx0_b( i )
	sig = sigflx0_b( i )
      end if

* augment flux error bar
	sigold = sig
      if( dexsed .ne. 0. ) then
	snr = dat / sig
	s = dat * sigsed
	sig = sqrt( sig * sig + s * s )
* report
      if( verbose ) then
	write(*,*) i, ' S/N', snr, ' ->', dat / sig
      end if
      end if
* report
      if( verbose ) then
	write(*,*) i, ' flx', dat, ' +/-', sig, ' S/N', dat / sig
      end if
	y = dat
* symmetric error bars
	yp = sig
	ym = yp
* log x
      if( logx .ne. 0 ) then
	xp = x + dx
	xm = x - dx
        x = alog10( x )
	xp = alog10( xp ) - x
	xm = x - alog10( xm )
      end if
* log y
      if( logy .ne. 0 ) then
c        floor = 10. ** ( ymn - 0.5 )
        yp = y + yp
        ym = y - ym
        y  = alog10( max( floor, y ) )
        yp = alog10( max( floor, yp ) ) - y
        ym = y - alog10( max( floor, ym ) )
      end if

* fit interpolates model fluxes
	call lookup( np, plot1, x, i1, i2, p )
      if( ibf .eq. ibrt ) then
	fit = terp1( np, plot2, i1, i2, p )
      else if( ibf .eq. ifnt ) then
	fit = terp1( np, plot3, i1, i2, p )
      end if

* asymmetric error bar
	dat = y
	sig = yp
	if( dat .lt. fit ) sig = ym

* chi
	chi = ( dat - fit ) / sig
* increment
	add = chi * chi
	if( ibf .eq. ibrt ) chi2brt = chi2brt + add
	if( ibf .eq. ifnt ) chi2fnt = chi2fnt + add
* error bar expansion penalty
      if( sig .ne. sigold ) then
	sumlnvarsed = sumlnvarsed + 2. * alog( sig / sigold )
      end if
* report
      if( verbose ) then
	write(*,*) i, ' dat', dat, ' fit', fit
	write(*,*) i, ' sig', sig, ' chi', chi
	if( ibf .eq. ibrt ) write(*,*) i, ' Bright Chi2', chi2brt
	if( ibf .eq. ifnt ) write(*,*) i, ' Faint Chi2', chi2fnt
      end if
* next datum i
      end do
* next bright/faint ibf
      end do

* end flux data ndat > 0
      end if

	ntot = nlag + ndat + ndat
	chi2 = chi2lag + chi2brt + chi2fnt
	sumlnvar = sumlnvarlag + sumlnvarsed
	bofcalc = chi2 + sumlnvar

* report
	write(*,*) 'Nlag', nlag, '    Lag Chi2', real( chi2lag )
     &		, ' sum(ln(var))', real( sumlnvarlag )
	write(*,*) 'Nbrt', ndat, ' Bright Chi2', real( chi2brt )
	write(*,*) 'Nfnt', ndat, '  Faint Chi2', real( chi2fnt )
     &		, ' sum(ln(var))', real( sumlnvarsed )
	write(*,*) 'Ntot', ntot, '  Total Chi2', real( chi2 )
     &		, ' sum(ln(var))', real( sumlnvar )
	write(*,*) ' BoF =', real(chi2), ' +', real(sumlnvar),' =', bofcalc
	return
	end
* ------------------------------
	subroutine setconstants
	use constants
* constants
	pi = 4. * atan2( 1., 1. )
	dtr = pi / 180.

* physics constants
	clight = cgs( 'C' )
	gnewt = cgs( 'G' )
	stefan = cgs( 'SIGMA' )
	sigmat = cgs( 'SIGMAT' )
	emp = cgs( 'MP' )
* astro constants
	pc = cgs( 'PC' )
	year = cgs( 'YEAR' )
	vkms = 1.e5
	day = 24. * 3600.

	sunM = cgs( 'MSUN' )
	sunL = cgs( 'LSUN' )
	sunRs = 2. * gnewt * ( sunM / clight / clight )
	sunRg = gnewt * ( sunM / clight / clight )
	sunLedd = 4. * pi * gnewt * (sunM * emp) * clight / sigmat

	write(*,*) 'Sun: M', sunM, ' L', sunL
	write(*,*) ' Rg/km', sunRg/vkms, ' Ledd', sunLedd, ' =', sunLedd/sunL, ' Lsun'
	return
	end

* -----------------------------------------------------
	subroutine geometry( nr, r_r, h_r, risco, rout, r1, h1, bet,
     &		wamp, wpow, wnum, wnpow, crest )
* compute the disc surface geometry
* In:	nr	i4 number of radii
*	risco	r4 inner radius (light days)
*	rout	r4 outer radius (light days)
*	r1	r4 reference radius (light days) (0 for r1=rout)
*	h1	r4 height above midplane at r1 (light days)
*	bet	r4 power-law index H(r) = H1 (r/r1)^bet
*	wamp	r4 ripple amplitude (fraction of h)
*	wpow	r4 power-law index A(r) = A (r/r1)^wpow
*	wnum	r4 ripple wave number k1 at r=r1
*	wnpow	r4 power-law index k(r) = k1 (r/r1)^wpow
*	crest	r4 ripple crest shape ( >1 for peakier)
* Out:	r_r(nr)	r4 radius grid - uniform in log(r)
*	h_r(nr)	r4 height grid h(r) = H(r) ( 1 + A(r) cos( phi(r) ) )
* 2021 Oct Keith Horne @ St.Andrews
* 2021 Nov KDH @ StA - mix linear and log spacing
* 2021 Nov KDH @ StA - ripple wave number wnum, wpow
* 2021 Nov KDH @ Crail - crest shape
* 2021 Nov KDH @ Crail - cosine to better sample near inner/outer edges
* 2023 Dec KDH @ StA - finer sampling near Rin
	real*4 r_r( nr )
	real*4 h_r( nr )
	pi = 4. * atan2( 1., 1. )

	write(*,*) '--- geometry ------------ Nr', nr


* KDH : SHOULD CHECK RIPPLE RESOLUTION
* radius grid -uniform in log(r)
* grid spacing is min of U(logr) and U(r)
	rlo = risco
	rhi = rout
	rlog1 = alog10( rlo )
	rlog2 = alog10( rhi )
	frac = 2.
	frac = 1.5

* finer sampling near rin
	dex = alog10( rhi / rlo )
	dpin = alog10( 2. ) / dex
* finer sampling near rout
	dpout = dpin
c	dpout = dpout * max( 1., bet - 1. )

      do i = 1, nr
	p = ( i - 1. ) / max( 1, nr - 1 )

* cosine to give more points near inner/outer edge
* KDH : could tune this more carefully, e.g. depending on beta.
c	p = ( 1. - cos( pi * p ) ) / 2.

* finer sampling near Rout
	pc = 1. - p
	pout = 1. - exp( - pc / dpout )
	pc = pc * pout
	p = 1. - pc

* finer sampling near Rin
	pin = 1. - exp( - p / dpin )
	p = p * pin


* log spacing
	r = 10. ** ( rlog1 * ( 1. - p ) + p * rlog2 )
      IF( .FALSE. ) THEN
* linear spacing
	rlin = rlo * ( 1. - p ) + p * rhi
* mix linear with stretched log sampling
	plog = p * frac
	rlog = rlog1 * ( 1. - plog ) + plog * rlog2
	rlog = 10. ** rlog
	rmix = min( rlog, rlin )
      END IF
	r_r( i ) = r
      end do

* reference radius r0
	r0 = r1
	if( r0 .le. 0. ) r0 = rout

* disk geometry
      do i = 1, nr
	r = r_r( i )
	x = r / r0
	h_r( i ) = h1 * x ** bet
      end do

* ripples
      if( wamp .ne. 0. ) then
* phase at r=0
	radian0 = 0.
      do i = 1, nr
	r = r_r( i )
	h = h_r( i )
	x = r / r0
	amp = wamp * x ** wpow
	amp = max( -1., min( 1., amp ) )
	pow = wnpow + 1.
* phase 0 at r=r0
c	radian = phi1 + wnum * r0 * ( x ** pow - 1. ) / pow
* phase 0 at r=0
	old = radian
	radian = radian0 + wnum * r0 * x ** pow / pow
	ripple = sin( radian )
* crest shape
	if( crest .ne. 1. ) ripple = 2. * ( ( ripple + 1. ) / 2. ) ** crest
	h = h * ( 1. + amp * ripple )
	h_r( i ) = h
	if( i .eq. 1 ) big = 0.
	if( i .gt. 1 ) big = max( big, radian - old )
      end do

      if( big .gt. 0.1 ) then
	write(*,*) 'Ripple resolution K dr > ', 2. * acos( 0. ) / big
      end if
      end if
	return
	end

* -----------------------------------------------
	subroutine lagcalc
* compute delay maps and lag spectra
* In:	nb	i4 number of bands
*	wobs_b(nb)	r4 pivot wavelengths (Angstroms)
*	nr	i4 number of radii
*	r_r(nr)	r4 radii of annuli (light days)
*	h_r(nr)	r4 heights above midplane (light days)
*	tv_r(nr) r4 Kelvin temperature - viscous heating
*	t_r(nr)	r4 Kelvin temperature - lamp+viscous heating
*	hlamp	r4 lamp height (light days)
*	cosi	r4 cos(inclination)
*	redshift r4 redshift
*	needlag	l4 .true. if lags needed
* Out:	ntau	i4 number of delays
*	tau_d(ntau)	r4 delay grid (days)
*	psi_db(ntau,nb)	r4 delay maps lamp on (prob/day)
*	psi0_db(ntau,nb) r4 delay maps lamp off (prob/day)
*	psi1_db(ntau,nb) r4 delay maps lamp no rim (prob/day)
*	tau_b(nb)	r4 mean delay lamp on (days)
*	tau0_b(nb)	r4 mean delay lamp off (days)
*	tau1_b(nb)	r4 mean delay lamp on no rim (days)
*	tau2_b(nb)	r4 mean delay lamp on rim (days)
*
*	wobs	r4 fiducial observed wavelength (Angstroms)
*	wrest	r4 fiducial rest wavelength (Angstroms)
*	psi_d(ntau)	r4 delay map lamp on (prob/day)
*	psi0_d(ntau)	r4 delay map lamp off (prob/day)
*	psi1_d(ntau)	r4 delay map lamp on no rim (prob/day)
*	psi2_d(ntau)	r4 delay map lamp on rim (prob/day)
*	tauavg	r4 mean delay lamp on (days)
*	tauavg0	r4 mean delay lamp off (days)
*	needlag	l4 .false.
* 2021 Nov Keith Horne @ Crail
* 2021 Nov KDH @ Crail - needlag
* 2021 Nov KDH @ Crail - verbose
* 2024 Apr KDH @ StA - psi1_db, tau1_db
	use all
	use agn

	logical verbose/ .false. /
	logical yes

	write(*,*) '--- lagcalc ------------- Nr', nr, ' needlag', needlag
	if( .not. needlag ) return
	needlag = .false.

* delay step
	drmn = r_r( nr )
      do i = 2, nr
	dr = abs( r_r( i ) - r_r( i - 1 ) )
	drmn = min( drmn, dr )
      end do
	dtau = 1. / 24.

* delay grid
	taumn = 0.
	taumx = 1000.
	dtau = min( dtau, dr )
	taumx = min( taumx, 3. * r_r( nr ) )

	ntau = 1 + ( taumx - taumn ) / dtau
	ntau = min( ntau, maxtau )
	dtau = ( taumx - taumn ) / max( 1, ntau - 1 )
	tau0 = taumn - dtau
      do i = 1, ntau
	tau_d( i ) = tau0 + i * dtau
	psi_d( i ) = 0.
      end do
* report
      if( verbose ) then
	write(*,*) 'Delay grid:', ntau, taumn, taumx, dtau
      end if

* rim present
	call imnmx( nr, t_r, imn, imx )
	yes = imn .lt. nr

* bands
	opz = 1. + redshift
      do ib = 1, nb
	wobs = wobs_b( ib )
	wrest = wobs / opz
* toggle irradiation
      do irr = 1, 0, -1
	if( verbose ) write(*,*) 'Irradiation', irr
* irradiation off
      if( irr .eq. 0 ) then
	call tfbowl( nr, r_r, h_r, tv_r,
     &		hlamp, cosi,   1., wrest, ntau, tau_d, psi_d )
* irradiation on
      else
	call tfbowl( nr, r_r, h_r, t_r,
     &		hlamp, cosi, fcol, wrest, ntau, tau_d, psi_d )
      end if

* mean delay
	top = 0.d0
	bot = 0.d0
      do i = 1, ntau
	psi = psi_d( i )
	top = top + psi * tau_d( i )
	bot = bot + psi
* stow
      if( irr .eq. 1 ) then
	 psi_db( i, ib ) = psi
* default to no rim
	psi1_db( i, ib ) = psi
	psi2_db( i, ib ) = 0.
      else
	psi0_db( i, ib ) = psi
      end if
      end do
	tau = top / bot

* KDH : TRAP NAN
      if( tau .ne. tau ) then
	write(*,*) '** INVALID LAG', top, ' /', bot, ' =', tau
      do i = 1, ntau
	write(*,*) i, tau_d( i ), psi_d( i )
      end do
	write(*,*) 'Irradiation', irr, ' Fcol', fcol
	stop
      end if

* stow
      if( irr .eq. 0 ) then
	tau0_b( ib ) = tau
      else
	 tau_b( ib ) = tau
	tau1_b( ib ) = tau
	tau2_b( ib ) = tau

* disc response omitting the rim
      if( yes ) then
	i = imn
	n = nr - i + 1

* KDH : 2025-Jan-01 INCLUDE IMN IN THE RIM TO FIX NR=1 CRASH IN TFBOWL
      if( n .le. 1 ) then
	write(*,*) '** Rim annuli', i, ' to', nr, ' total', n, ' :('
	i = max( 1, imn - 1 )
	n = nr - i + 1
	write(*,*) '** Rim annuli', i, ' to', nr, ' total', n, ' OK?'
      end if

      if( i .le. 1 ) then
	write(*,*) '** Rim annuli', i, ' to', nr, ' total', n, ' :('
	i = min( nr, imn + 1 )
	n = nr - i + 1
	write(*,*) '** Rim annuli', i, ' to', nr, ' total', n, ' OK?'
      end if

* disc response omitting the rim
	call tfbowl( i, r_r, h_r, t_r,
     &		hlamp, cosi, fcol, wrest, ntau, tau_d, psi1_d )

      if( verbose ) then
	write(*,*) 'tfbowl(', n, i, nint( t_r( i ) ), nint( t_r( nr ) ), ' )'
      end if
* rim response
	call tfbowl( n, r_r( i ), h_r( i ), t_r( i ),
     &		hlamp, cosi, fcol, wrest, ntau, tau_d, psi2_d )

* mean lags
      do k = 1, 2
	top = 0.d0
	bot = 0.d0
      do i = 1, ntau
	tau = tau_d( i )
	if( k .eq. 1 ) psi = psi1_d( i )
	if( k .eq. 2 ) psi = psi2_d( i )
	top = top + psi * tau
	bot = bot + psi
      end do
	tau = top / bot
* panic
      if( bot .le. 0.d0 ) then
	write(*,*) '** PANIC in lagcalc'
      do i = 1, ntau
	write(*,*) i, tau_d( i ), psi_d( i ), psi1_d( i ), psi2_d( i )
      end do
	write(*,*) k, ' top', top, ' bot', bot, tau
	tau = 0.
      end if
	if( k .eq. 1 ) tau1_b( ib ) = tau
	if( k .eq. 2 ) tau2_b( ib ) = tau
      end do
      end if

      end if
* next irr
      end do

* report
      if( verbose ) then
	tau = tau_b( ib )
* KDH : WARNING IS TAU0 THE LAG AT THE REFERENCE WAVELENGTH ?
	tau0 = tau0_b( ib )
	tau1 = tau1_b( ib )
	tau2 = tau2_b( ib )
	call mnmx( ntau,  psi_db(1,ib), p, psi )
	call mnmx( ntau, psi0_db(1,ib), p, psi0 )
	call mnmx( ntau, psi1_db(1,ib), p, psi1 )
	call mnmx( ntau, psi2_db(1,ib), p, psi2 )
	write(*,*) 'Wrest', nint( wrest )
	write(*,*) ' Flat,viscous,face-on Psi0', psi0, ' tau0', tau0
	write(*,*) ' Curved,flashed,tilted Psi', psi,  ' tau ', tau
	write(*,*) '     Ditto but no rim Psi1', psi1, ' tau1', tau1
	write(*,*) '     Ditto but    rim Psi2', psi2, ' tau2', tau2
      end if

* stow
      do i = 1, ntau
	psi1_db( i, ib ) = psi1_d( i )
	psi2_db( i, ib ) = psi2_d( i )
      end do

* next band ib
      end do

c	wobs = 4132
c	wrest = wobs / ( 1. + redshift )
	wrest = 5000.
	wrest = 10000.
	wobs = wrest * ( 1. + redshift )
	if( verbose ) write(*,*) 'Wobs', wobs, ' wrest', wrest

* fiducial disk (flat, face-on, viscous) ------------------------
c	call tfb_inclined( nr, r_r, h0_r, t0_r,
c     &		hlamp, 1.0, wrest, ntau, tau_d, psi0_d )
* fiducial disk (no irradiation)
c	call tfb_inclined( nr, r_r, h0_r, t0_r,
c     &		hlamp, cosi, wrest, ntau, tau_d, psi0_d )
* curved rippled disk (no irradiation)
c	call tfb_inclined( nr, r_r, h_r, tv_r,
c     &		hlamp, cosi, wrest, ntau, tau_d, psi0_d )
* no colour correction fcol = Tc/Teff
	call tfbowl( nr, r_r, h_r, tv_r,
     &		hlamp, cosi, 1., wrest, ntau, tau_d, psi0_d )

	call mnmx( ntau, psi0_d, p0mn, p0mx )

* mean delay
	top = 0.d0
	bot = 0.d0
      do i = 1, ntau
	psi = psi0_d( i )
	tau = tau_d( i )
	top = top + psi * tau
	bot = bot + psi
      end do
	tauavg0 = top / bot

* curved, tilted, irradiated disk ----------------------------
c	call tfb_inclined( nr, r_r, h_r, t_r,
c     &		hlamp, cosi, wrest, ntau, tau_d, psi_d )
* include colour correction fcol = Tc/Teff
	call tfbowl( nr, r_r, h_r, t_r,
     &		hlamp, cosi, fcol, wrest, ntau, tau_d, psi_d )
* mean delay
	top = 0.d0
	bot = 0.d0
      do i = 1, ntau
	tau = tau_d( i )
	psi = psi_d( i )
	top = top + psi * tau
	bot = bot + psi
      end do
	tauavg = top / bot

* report
      if( verbose ) then
	call mnmx( ntau, tau_d, tmn, tmx )
	call mnmx( ntau, psi0_d, p0mn, p0mx )
	call mnmx( ntau, psi_d, pmn, pmx )
	write(*,*) 'Delays', nd, tmn, tmx, ( tmx - tmn ) / max( 1, nd-1 )
	write(*,*)
	write(*,*) 'Psi0(tau)', pmn, pmx, ' <tau>', tauavg0
	write(*,*) ' Psi(tau)', pmn, pmx, ' <tau>', tauavg
	write(*,*)
	write(*,*) 'D(Mpc)', dmpc, ' i(deg)', dinc, ' cosi', cosi
	write(*,*) 'f_nu(mJy)', f_w( 1 ), f_w( nw )
      end if

	return
	end

* ----------------------------------------------
	subroutine sedcalc
* compute the disc spectra with and without lamp heating
* In:	nr	i4 number of annuli
*	r_r(nr)	r4 radii of annuli (light days)
*	h_r(nr)	r4 heights above midplane (light days)
*	hlamp	r4 lamp height (light days)
*	risco	r4 ISCO radius (light days)
*	rout	r4 outer radius (light days)
*	r1	r4 reference radius (light days) (rout if r1=0)
*	alp	r4 viscous heating index (3/4)
*	tv1	r4 Kelvin temperature at r1 - viscosity
*	tx1	r4 Kelvin temperature at r1 - face-on irradiation
*	dmpc	r4 distance (Mpc)
*	cosi	r4 cos(inclination)
*	ebmv_agn r4 E(B-V) rest-frame SMC-like dust
*	needsed l4 .true. if fluxes needed
* Out:	tv_r(nr) r4 Kelvin temperature - viscosity
*	tx_r(nr) r4 Kelvin temperature - irradiation
*	t_r(nr)	r4 Kelvin temperature - both
*	tc_r(nr) r4 Kelvin temperature - scaled by fcol
*	i_r(nr)	i4 0 (in shadow) or 1 (exposed to lamp)
*	nw	i4 number of wavelengths
*	w_w(nw)	r4 wavelengths (Angstrom)
*	f_w(nw)	r4  sed with lamp on  (mJy)
*	f0_w(nw) r4 sed with lamp off (mJy)
*	f1_w(nw) r4 sed with lamp on, no rim (mJy)
*	f2_w(nw) r4 sed with lamp on, rim (mJy)
*	needsed	l4 .false.
* 2021 Nov Keith Horne @ Crail
* 2021 Nov KDH @ Crail - needflx
* 2021 Nov KDH @ Crail - verbose
* 2023-Dec KDH @ StA - flxcalc,needflx -> sedcalc,needsed
* 2024-Apr KDH @ StA - f1_w, f2_w
* 2025-Mar-07 KDH @ StA - ebmv_agn

	use all
	use agn

	logical verbose/ .false. /

	write(*,*) '--- sedcalc ------------- Nr', nr, ' needsed', needsed
	if( .not. needsed ) return
	needsed = .false.

      if( nr .le. 0 ) then
	write(*,*) '** FATAL ERROR.  NR=', nr
	stop
      end if

* reference radius
	r0 = r1
	if( r0 .le. 0. ) r0 = rout
	if( r0 .le. 0. ) r0 = r_r( nr )

* Teff profile for curved disc with viscous+lamp heating
	call tvtx( nr, r_r, h_r, hlamp, risco, r0, alp, qin,
     &		tv1, tx1, tv_r, tx_r, tq_r, t_r )

* power-law Teff profile
      do i = 1, nr
	r = r_r( i )
	t0_r( i ) = tv1 * ( r0 / r ) ** alp
      end do

* shadows
	call diskshad( nr, r_r, h_r, hlamp, i_r )
      do i = 1, nr
      if( i_r( i ) .eq. 0 ) then
	t_r( i ) = tv_r( i )
* may cause problems on log-log plots
	tx_r( i ) = 0.
      end if
      end do

* colour temperature
      do i = 1, nr
	tc_r( i ) = t_r( i ) * fcol
      end do


* report
      if( verbose ) then
	write(*,*) 'R1', r0, ' H1', h1, ' beta', bet
	write(*,*) 'Tv1', tv1, ' alpha', alp
	write(*,*) 'Tx1', tx1, ' albedo', albedo
	write(*,*) ' Xlamp', X, ' Ylamp', Y, ' fcol', fcol
	rlog1 = alog10( risco )
	rlog2 = alog10( rout )
	write(*,*) 'Nr', nr, ' logR', rlog1, rlog2
	write(*,*) 'Radii ', r_r( 1 ), r_r( nr )
	write(*,*) 'Height', h_r( 1 ), h_r( nr )
	write(*,*) '  Tv  ', tv_r( 1 ), tv_r( nr )
	write(*,*) '  Tx  ', tx_r( 1 ), tx_r( nr )
	write(*,*) '  Tq  ', tq_r( 1 ), tq_r( nr )
	write(*,*) '  T0  ', t0_r( 1 ), t0_r( nr )
	write(*,*) ' Teff ',  t_r( 1 ),  t_r( nr )
	write(*,*) '  Tc  ', tc_r( 1 ), tc_r( nr )
      end if


* wavelength grid (Angstroms)
	wlog1 = 2.
	wlog2 = 5.5
	nw = 1000
	nw = min( nw, maxw )
      do i = 1, nw
	p = ( i - 1. ) / max( 1, nw - 1 )
	wlog = wlog1 * ( 1. - p ) + p * wlog2
	w = 10. ** wlog
	w_w( i ) = w
      end do

* report
      if( verbose ) then
	write(*,*) 'Nw', nw, ' log(w)', wlog1, wlog2
	write(*,*) 'Wavelengths', w_w( 1 ), w_w( nw )
      end if

* disk spectrum
* flat face-on viscous disk with no irradiation
c	call fnudisk_inclined( nr, r_r, h0_r, t0_r, dmpc,  1.0, nw, w_w, f0_w )
* curved rippled inclined disk, no irradiation
c	call fnudisk_inclined( nr, r_r,  h_r, tv_r, dmpc, cosi, nw, w_w, f0_w )

* include colour correction fcol = Tc/Teff
	call fnubowl( nr, r_r,  h_r, tv_r, dmpc, cosi, fcol, nw, w_w, f0_w )


* KDH : TESTING
      IF( .FALSE. ) THEN
c      IF( .TRUE. ) THEN
      do i = 1, nr
      if( i .le. 5 .or. i .eq. nr .or. mod( i, max(1,nr/5) ) .eq. 0 ) then
	write(*,*) 'tv_r(', r_r(i), ' ) =', tv_r( i )
      end if
      end do
      do i = 1, nw
      if( i .eq. 1 .or. i .eq. nw .or. mod( i, max(1,nw/5) ) .eq. 0 ) then
	write(*,*) 'f0_w(', nint( w_w(i) ), ' ) =', f0_w( i )
      end if
      end do
      END IF

* curved rippled inclined disk yes irradiation
c	call fnudisk_inclined( nr, r_r, h_r, t_r, dmpc, cosi, nw, w_w, f_w )

* include colour correction fcol = Tc/Teff
	call fnubowl( nr, r_r, h_r, t_r, dmpc, cosi, fcol, nw, w_w, f_w )

* repeat f1 = norim, f2 = rim
* KDH : FASTER TO SUM OVER THE RIM?
* KDH : THIS CAN FAIL IF IMN = 1 due to zero torque
c	call imnmx( nr, t_r, imn, imx )
	i = max( 1, nr / 2 )
	n = nr - i + 1
	call imnmx( n, t_r( i ), imn, imx )
	imn = imn + i - 1
	imx = imx + i - 1

      if( imn .lt. nr ) then
	call fnubowl( imn, r_r, h_r, t_r, dmpc, cosi, fcol, nw, w_w, f1_w )
      do i = 1, nw
	f2_w( i ) = f_w( i ) - f1_w( i )
      end do
      else
      do i = 1, nw
	f1_w( i ) = f_w( i )
	f2_w( i ) = 0.
      end do
      end if

* apply rest-frame SMC-like dust extinction
	ebmv = ebmv_agn
      if( ebmv .ne. 0. ) then
	opz = 1. + redshift
      do i = 1, nw
	w = w_w( i ) / opz
	fac = 10. ** ( -0.4 * extmag_smc( w, ebmv ) )
	 f_w( i ) =  f_w( i ) * fac
	f0_w( i ) = f0_w( i ) * fac
	f1_w( i ) = f1_w( i ) * fac
	f2_w( i ) = f2_w( i ) * fac
      end do
      end if


* report
      if( verbose ) then
	call mnmx( nw,  f_w,  fmn,  fmx )
	call mnmx( nw, f0_w, f0mn, f0mx )
	call mnmx( nw, f1_w, f1mn, f1mx )
	call mnmx( nw, f2_w, f2mn, f2mx )
	write(*,*)
c	write(*,*) ' Fnu/mJy', fmn, fmx, ' face-on,flat,viscous disk'
	write(*,*) ' Fnu/mJy',  fmn,  fmx, ' Bright disk(+rim)'
	write(*,*) ' Fnu/mJy', f1mn, f1mx, ' Bright no rim'
	write(*,*) ' Fnu/mJy', f2mn, f2mx, ' Bright rim'
	write(*,*) ' Fnu/mJy', f0mn, f0mx, ' Faint disk'
	write(*,*)
      end if
	return
	end

* ------------------------------
	subroutine initpar
* Install BOWL parameters in the MCMC framework
* Out:	nfit    	i4 number of MCMC parameters
*	code_p(nfit)	c1 parameter 1-character codes
*	name_p(nfit)	c* parameter names
*	units_p(nfit)	c* parameter units
*	symbol_p(nfit)	c* parameter pgplot symbol
*	 log_p(nfit)	l  .true. for log10(parameter)
*	 bot_p(nfit)	r4 parameter lower limit
*	 top_p(nfit)	r4 parameter upper limit
*	   p_p(nfit)	r4 parameter values
*	  dp_p(nfit)	r4 parameter steps (set if 0 on input)
* 2022 Jul KDH @ Crail
* 2022 Jul KDH @ StA - add units and symbols
* 2022 Jul KDH @ StA - add M/Msun and D/Mpc
* 2022 Jul KDH @ StA - h1_r1_pct
* 2023 Jul KDH @ Warsaw - re-order parameters
* 2024 Jan KDH @ StA - max Rout = 100 => 1000 d (for 3C273)
* 2024 Jan KDH @ StA - call secondary (since dmpc = 0 on first call)
* 2024 Jan KDH @ StA - fcol
* 2025 Jan KDH @ StA - primary_xlamp for xlamp vs dlogt
* 2025-Mar-07 KDH @ StA - add ebmv_agn
* 2025-Mar-12 KDH @ StA - keep log( dlogt ) < 0
* 2025-Apr-14 KDH @ StA - primary_dh, dh_r1_pct

	use all
	use constants
	use agn

* compute secondary parameters
	call secondary

* parameter count
	n = 0

* log( D/Mpc )
	n = n + 1
	code_p( n ) = 'D'
	name_p( n ) = 'D'
	units_p( n ) = 'Mpc'
	symbol_p( n ) = 'D\dL\u'
	log_p( n ) = .true.
	bot_p( n ) = -5
	top_p( n ) = +10
	p_p( n ) = alog10( Dmpc )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05

	write(*,*) '** initpar Dmpc', dmpc

* E(B-V)
	n = n + 1
	code_p( n ) = 'E'
	name_p( n ) = 'E(B-V)'
	units_p( n ) = ' '
	symbol_p( n ) = 'E(B-V)'
	p = ebmv_agn
	log_p( n ) = .false.
	bot_p( n ) = 0.
	top_p( n ) = 2.
	p_p( n ) = p
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.01
	write(*,*) '** initpar E(B-V)', ebmv_agn

* log( M/Msun )
	n = n + 1
	code_p( n ) = 'M'
	name_p( n ) = 'Mbh'
	units_p( n ) = 'M\d\(2281)\u'
	symbol_p( n ) = 'M'
	log_p( n ) = .true.
	bot_p( n ) = -2
	top_p( n ) = +13
	p_p( n ) = alog10( bhM_Msun )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05
	write(*,*) '** initpar Mbh', bhM_Msun

* log( mdot )
	n = n + 1
	code_p( n ) = 'A'
	name_p( n ) = 'Mdot'
	units_p( n ) = 'M\d\(2281)\u/yr'
	symbol_p( n ) = ' \u\(828)\b\b\d\bM'
	call pgmdot( symbol_p( n ) )
	log_p( n ) = .true.
	bot_p( n ) = -10
	top_p( n ) = +5
	p_p( n ) = alog10( bhMsun_yr )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05
	write(*,*) '** initpar Mdot', bhMsun_yr

* cos( i )
	n = n + 1
	code_p( n ) = 'C'
	name_p( n ) = 'cos(i)'
	units_p( n ) = ' '
	symbol_p( n ) = 'cos(i)'

	bot_p( n ) = cos( dincmx * dtr )
	top_p( n ) = cos( dincmn * dtr )
	p_p( n ) = cosi
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05
	write(*,*) '** initpar cosi', cosi, ' Dinc', dinc

* log( xlamp )
      if( primary_xlamp ) then
	n = n + 1
	code_p( n ) = 'X'
	name_p( n ) = 'Xlamp'
	units_p( n ) = ' \u\(828)\b\b\d\bMc\u2\d'
	symbol_p( n ) = '\gDL\dx\u'
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = +3
	 p_p( n ) = alog10( xlamp )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.2
	write(*,*) '** initpar Xlamp', xlamp

* log( dlogT )
      else
	n = n + 1
	code_p( n ) = 'X'
	name_p( n ) = 'DlogT'
	units_p( n ) = 'dex'
	symbol_p( n ) = '\gDlogT'
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = 0.
	 p_p( n ) = alog10( dlogt )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.2
	write(*,*) '** initpar DlogT', dlogt
      end if

* log( hlamp_rg )
	n = n + 1
	code_p( n ) = '*'
	name_p( n ) = 'Hx/Rg'
	units_p( n ) = 'R\dg\u'
	symbol_p( n ) = 'H\dx\u'
	log_p( n ) = .true.
	bot_p( n ) = 0.
	top_p( n ) = +10
	p_p( n ) = alog10( hlamp_rg )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05
	write(*,*) '** initpar Hx/Rg', hlamp_rg

* log( rin )
	n = n + 1
	code_p( n ) = 'I'
	name_p( n ) = 'Rin'
c	units_p( n ) = 'GM/c\u2\d'
	units_p( n ) = 'R\dg\u'
	symbol_p( n ) = 'R\din\u'
	log_p( n ) = .true.
	bot_p( n ) = alog10( 1. )
	top_p( n ) = alog10( 100. )
	p_p( n ) = alog10( risco_rg )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05
	write(*,*) '** initpar Rin/Rg', risco_rg

* log( rout )
	n = n + 1
	code_p( n ) = 'O'
	name_p( n ) = 'Rout'
	log_p( n ) = .true.
	bot_p( n ) = -3
* KDH : may need +3 for quasars like 3C273
	top_p( n ) = +3
c	top_p( n ) = +2
	p_p( n ) = alog10( rout )
	units_p( n ) = 'lt.d'
	symbol_p( n ) = 'R\dout\u'
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.1
	write(*,*) '** initpar Rout/c', rout, ' d'

* H/R ( h1_r1_pct )
	n = n + 1
	code_p( n ) = 'H'
	name_p( n ) = 'H/R'
	units_p( n ) = '%'
	symbol_p( n ) = 'H/R'
	p = h1_r1_pct
* (H-Hx)/R ( dh_r1_pct )
      if( primary_dh ) then
	name_p( n ) = '(H-Hx)/R'
	symbol_p( n ) = '\gDH/R'
	p = dh_r1_pct
      end if
* log
      if( p .gt. 0. ) then
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = 1.
	p_p( n ) = alog10( p )
* linear
      else
	log_p( n ) = .false.
	bot_p( n ) = 0.
	top_p( n ) = 10.
	p_p( n ) = p
      end if
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.05
	write(*,*) '** initpar H/R', h1_r1_pct, ' %'
	write(*,*) '      (H-Hx)/R', dh_r1_pct, ' %', primary_dh

* beta = dlnT/dlnR
	n = n + 1
	code_p( n ) = 'B'
	name_p( n ) = 'RdH/HdR'
	name_p( n ) = 'dlnH/dlnR'
	units_p( n ) = ' '
	symbol_p( n ) = '\gb'
	p = bet
* log
      if( p .gt. 0. ) then
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = +3
	p_p( n ) = alog10( p )
* linear
      else
	log_p( n ) = .false.
	bot_p( n ) = -100
	top_p( n ) = +100
	p_p( n ) = p
      end if
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.1
	write(*,*) '** initpar beta', bet


* torque parameter
	n = n + 1
	code_p( n ) = '@'
	name_p( n ) = 'Q'
	units_p( n ) = ' '
	symbol_p( n ) = 'Q'
	p = qin
* log
      if( .TRUE. .and. p .ge. 0. ) then
* KDH : P CAN BE OUT OF RANGE HERE ?
	p = max( p, 1.e-5 )
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = +3
	p_p( n ) = alog10( p )
* linear
      else
	log_p( n ) = .false.
	bot_p( n ) = 0.
	top_p( n ) = 100.
	p_p( n ) = p
      end if
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.1
	write(*,*) '** initpar Qin', qin

* fcol = Tc/Teff
	n = n + 1
	code_p( n ) = '^'
	name_p( n ) = 'fcol'
	units_p( n ) = ' '
	symbol_p( n ) = 'f\dcol\u'
	p = fcol
      if( p .gt. 0. ) then
	log_p( n ) = .true.
	bot_p( n ) = -3.
	top_p( n ) = +3.
	p_p( n ) = alog10( p )

* KDH could let fcol < 0 implement Done+Kobota
      else
	log_p( n ) = .false.
	bot_p( n ) = 1.
	top_p( n ) = 1.
	p_p( n ) = p
      end if

* log( dexsed )
	n = n + 1
	code_p( n ) = 'F'
	name_p( n ) = 'dexsed'
	units_p( n ) = 'dex'
	symbol_p( n ) = '\gs\d0\u'
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = 3
	p_p( n ) = alog10( dexsed )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.2
	write(*,*) '** initpar sig_sed', dexsed, ' dex'

* log( syslag )
	n = n + 1
	code_p( n ) = 'L'
	name_p( n ) = 'syslag'
	units_p( n ) = 'd'
	symbol_p( n ) = '\gs\d\gt\u'
	log_p( n ) = .true.
	bot_p( n ) = -3
	top_p( n ) = 3
	p_p( n ) = alog10( syslag )
	if( dp_p( n ) .le. 0. ) dp_p( n ) = 0.2
	write(*,*) '** initpar sig_lag', syslag, ' d'


* total parameters
	npar = n

* report
	write(*,*)
	write(*,*) '--- initpar ------', npar
      do i = 1, npar
	code = code_p( i )
	name = name_p( i )
	p = p_p( i )
	dp = dp_p( i )
c	if( log_p( i ) ) name = 'log[' // name(:len1(name)) // ']'
	if( log_p( i ) ) name = 'log(' // name(:len1(name)) // ')'
	write(*,*) i, code, p, dp, ' ', name(:len1(name))
      end do

	return
	end


* ------------------------------
	subroutine pluck( iter )
* extract fit parameters from MCMC chain
* In:	iter	i4 MCMC sample (0 for median)
* In from module all:
*	nits	i4 number of iterations
*	npar	i4 number of parameters
*	nfit	i4 number of fit parameters
*	p_pi(nfit,nits)	r4 fit parameters from all iterations
*	bof_i(nits)	r8 Badness-of-fit values
* Out to module all:
*	p_p(nfit)	r4 selected parameters
* 2023 Nov Keith Horne @ St-Andrews
	use constants
	use all
	use agn

	imed = 0
      if( iter .gt. 0 ) then
	write(*,*) '** pluck MCMC sample', iter, ' of ', nits
      else
	write(*,*) '** pluck MEDIAN of ', nits, ' MCMC samples.'
      end if

	write(*,*) 'Fit parameters', nfit, ' of', npar
      if( nfit .le. 0 ) then
	write(*,*) 'No fit parameters.', nfit
	return
      end if

      if( iter .lt. imed ) then
	write(*,*) '*** ERROR: MCMC sample', iter, ' <', imed
	return
      else if( iter .gt. nits ) then
	write(*,*) '*** ERROR: MCMC sample', iter, ' >', nits
	return
      end if
* report
	write(*,*) 'BoF(', 1, ' ) =', bof_i( 1 )
	write(*,*) 'BoF(', ibest, ' ) =', bof_i( ibest ), ' best'
      if( iter .gt. 0 ) then
	bold = bof_i( iter )
	write(*,*) 'BoF(', iter, ' ) =', bold
      else
	bold = bofmed
	write(*,*) 'BoF( median ) =', bofmed
      end if
	write(*,*) 'BoF(', nits, ' ) =', bof_i( nits )

* extract parameters
      do i = 1, npar
	name = name_p( i )
	old = p_p( i )
* to p_p from p_pi
      if( iter .gt. 0 ) then
	p = p_pi( i, iter )
* to p_p from pmed_p
      else if( iter .eq. imed ) then
	p = pmed_p( i )
      end if
	write(*,*) i, name(:len1(name)), old, ' =>', p
	p_p( i ) = p
      end do
* extract primary parameters from p_p
	call getpar
* update secondaries ( done in bofcalc )
	call secondary
* new BoF
	needsed = .true.
	needlag = .true.
	bof = bofcalc( p_p )
	write(*,*) 'Re-evaluated BoF', bof, ' vs', bold
	return
	end
* ------------------------------
	subroutine setpar
* stow MCMC-safe parameters and steps for fitlist parameters
* In:	npar    	i4 number of parameters
*	code_p(npar)	c1 parameter codes
*	log_p(npar)	l4 .true. for log-uniform prior
*	bot_p(npar)	r4 lower limit
*	top_p(npar)	r4 upper limit
*	fitlist	c*	c* codes of parameters to fit
* Out:	nfit    	i4 number of fit parameters
*	i_p(nfit)	i4 parameter index
*	p_p(i_p(nfit))	r4 parameter values
*	dp_p(i_p(nfit))	r4 parameter steps
*
* 2021 Nov Keith Horne @ Crail - hardwire parameters
*      xlamp, dexsed, h1_hs, bet, rout
* 2022 Jul KDH @ Crail - bot_p, top_p
* 2022 Jul KDH @ Crail - common/mcmc/ -> module agn
* 2022 Jul KDH @ Crail - recode to use fitlist codes and index i_p
* 2022 Jul KDH @ StA - replace h1_rg by h1_r1
* 2022 Jul KDH @ StA - add M/Msun and D/Mpc
* 2022 Jul KDH @ StA - replace h1_r1 by h1_r1_pct
* 2023 Nov KDH @ StA - add qin
* 2024 Jan KDH @ StA - add fcol
* 2025 Jan KDH @ StA - primary_xlamp for xlamp vs dlogt
* 2025-Mar-06 KDH @ StA - ebmv_agn
* 2025-Apr-14 KDH @ StA - primary_dh, dh_r1_pct

	use all
	use constants
	use agn

	logical yes

* sanity checks
* no parameters
      if( npar .le. 0 ) then
	write(*,*) '** ERROR IN SETPAR : NO PARAMETERS ', npar
	stop
      end if
* no parameters to fit
	nfit = len1( fitlist )
      if( nfit .le. 0 ) then
	write(*,*) '** ERROR IN SETPAR : EMPTY FITLIST', nfit
	return
      end if

* use fitlist to construct index to MCMC fit parameters
	n = 0
      do i = 1, npar
	code = code_p( i )
	yes = index( fitlist, code ) .gt. 0
	fit_p( i ) = yes
      if( yes ) then
	n = n + 1
	i_p( n ) = i
      end if
      end do
	nfit = n

* stow fit parameter values
* KDH : MUST MATCH GETPAR
      do ifit = 1, nfit
	i = i_p( ifit )
	code = code_p( i )
* A = accretion rate Mdot
      if( code .eq. 'A' ) then
	p = bhMsun_yr
* X = lamp-post boost factor X
      else if( code .eq. 'X' ) then
	p = dlogt
	if( primary_xlamp ) p = xlamp
* F = systematic flux uncertainty dexsed
      else if( code .eq. 'F' ) then
	p = dexsed
* L = systematic lag uncertainty
      else if( code .eq. 'L' ) then
	p = syslag
* I = Inner radius Rin
      else if( code .eq. 'I' ) then
	p = risco_rg
* O = Outer radius Rout
      else if( code .eq. 'O' ) then
	p = rout
* H = H(R1)
      else if( code .eq. 'H' ) then
c	p = H1_rg
	p = h1_r1_pct
      if( primary_dh ) then
	p = dh_r1_pct
      end if
* B = beta = dlnH/dlnR
      else if( code .eq. 'B' ) then
	p = bet
* C = cos(i)
      else if( code .eq. 'C' ) then
	p = cosi
* KDH : WHY DO WE NEED TO ADJUST THESE HERE?
	sini = 1. - p * p
	dinc = acos( p ) / dtr
* * = Hlamp
      else if( code .eq. '*' ) then
	p = hlamp_rg
* D = D/Mpc
      else if( code .eq. 'D' ) then
	p = dmpc
* E(B-V) for SMC dust
      else if( code .eq. 'E' ) then
	p = ebmv_agn
* M = M/Msun
      else if( code .eq. 'M' ) then
	p = bhM_Msun
* qin torque parameter
      else if( code .eq. '@' ) then
	p = qin
* fcol = Tc/Teff
      else if( code .eq. '^' ) then
	p = fcol
      else
	write(*,*) '** ERROR : INVALID PARAMETER CODE ', CODE
	stop
      end if

* default parameter step
	dex = 0.05
* log-uniform prior
      if( log_p( i ) ) then
	p = alog10( p )
	dp = dex
      else
	dp = p * ( 10. ** dex - 1. )
      end if
* clip
	bot = bot_p( i )
	top = top_p( i )
	p = max( bot, min( top, p ) )
* stow
	p_p( i ) = p
	if( dp_p( i ) .le. 0. ) dp_p( i ) = dp

* next fit parameter ifit
      end do

* report
	write(*,*) '--- setpar ------', nfit, ' ', fitlist(:len1(fitlist))
      do ifit = 1, nfit
	i = i_p( ifit )
	name = name_p(i)
	code = code_p(i)
c	if( log_p( i ) ) name = 'log[' // name(:len1(name)) // ']'
	if( log_p( i ) ) name = 'log(' // name(:len1(name)) // ')'
	p = p_p( i )
	dp = dp_p( i )
	write(*,*) i, code, p, dp, ' ', name(:len1(name))
      end do

	return
	end
* -------------------------------
	subroutine getpar
* extract primary parameters from p_p
* In:	nfit		i4 number of fit parameters
*	i_p(nfit)	i4 index to Nfit of Npar parameters
*	p_p(nfit)	i4 parameter values
*	log_p(nfit)	.true. for log-uniform prior
*	name_p(nfit)	c* parameter name
*	code_p(nfit)	c* parameter code
* Out:
*
* 2021 Nov Keith Horne @ Crail - hardwire parameters
*      xlamp, dexsed, h1_hs, bet
* 2021 Nov KDH @ Crail - test for changes
* 2021 Nov KDH @ Crail - log_p
* 2022 Jul KDH @ Crail - common/mcmc/ -> module agn
* 2022 Jul KDH @ Crail - recode to use i_p(nfit)
* 2022 Jul KDH @ STA - add M/Msun and D/Mpc
* 2022 Jul KDH @ STA - h1_r1_pct replaces h1_r1 as primary
* 2023 Jun KDH @ STA - fix bug affecting needsed, needlag
* 2023 Nov KDH @ STA - add qin
* 2023 Dec KDH @ STA - extract npar primary parameters, not just nfit
* 2024 Jan KDH @ STA - add fcol
* 2025 Jan KDH @ StA - primary_xlamp for xlamp vs dlogt
* 2025-Mar-06 KDH @ StA - ebmv_agn
* 2025-Apr-14 KDH @ StA - primary_dh, dh_r1_pct
	use all
	use constants
	use agn

	logical verbose, change, changes
	verbose = .true.

* nothing to do
	if( nfit .le. 0 ) return

* MCMC parameter changes
	changes = .false.

* KDH : MOVING ALL PARAMETERS CAUSES COSI=0 AND DMPC=0
*        MCMCLOAD NEEDS TO MOVE ALL PARAMETERS
c      do i = 1, npar

      do ifit = 1, nfit
	i = i_p( ifit )
	code = code_p( i )
	name = name_p( i )
	p = p_p( i )
	if( log_p( i ) ) p  = 10. ** p

* extract primary parameters for MCMC fitting
	change = .false.
      if( code .eq. 'A' ) then
	old = bhMsun_yr
	bhMsun_yr = p
      else if( code .eq. 'X' ) then
      if( primary_xlamp ) then
	old = xlamp
	xlamp = p
      else
	old = dlogt
	dlogt = p
      end if
      else if( code .eq. 'F' ) then
	old = dexsed
	dexsed = p
      else if( code .eq. 'L' ) then
	old = syslag
	syslag = p
      else if( code .eq. 'I' ) then
	old = risco_rg
	risco_rg = p
      else if( code .eq. 'O' ) then
	old = rout
	rout = p
      else if( code .eq. '@' ) then
	old = qin
	qin = p
      else if( code .eq. 'E' ) then
	old = ebmv_agn
	ebmv_agn = p
      else if( code .eq. '^' ) then
	old = fcol
	fcol = p
      else if( code .eq. 'H' ) then
      if( primary_dh ) then
	old = dh_r1_pct
	dh_r1_pct = p
      else
	old = h1_r1_pct
	h1_r1_pct = p
      end if
      else if( code .eq. 'B' ) then
	old = bet
	bet = p
      else if( code .eq. 'C' ) then
	old = cosi
	cosi = p
* dinc is the primary parameter
	sini = sqrt( 1 - p * p )
	dinc = acos( p ) / dtr
      else if( code .eq. '*' ) then
	old = hlamp_rg
	hlamp_rg = p
      else if( code .eq. 'M' ) then
	old = bhM_Msun
	bhM_Msun = p
      else if( code .eq. 'D' ) then
	old = dmpc
	dmpc = p
      else
	write(*,*) '** ERROR IN GETPAR : INVALID CODE ', code
      end if
* report changes
	change = old .ne. p

      if( change ) then
	write(*,*) '     Model ', fitlist(:len1(fitlist))
     &	, ' Code ',  code, ' ', name(:len1(name)), old, '->', p
c	write(*,*) name(:len1(name)), old, ' ->', p
      end if
	changes = changes .or. change

* next fit parameter i = ifit_p( ifit )
      end do

* toggle need for updates
c	needlag = changes
c	needsed = changes
* revise input values
	needlag = needlag .or. changes
	needsed = needsed .or. changes

* update secondary parameters
	if( changes ) call secondary

	return
	end

* -------------------------------------
	subroutine showprimary
* report primary parameters
* 2023 Dec Keith Horne @ St Andrews
* 2025-Mar-06 KDH @ StA - ebmv_agn
* 2025-Apr-14 KDH @ StA - primary_dh, dh_r1_pct
c	use all
	use cosmology
	use agn
	write(*,*)
	write(*,*) '-- primary BOWL parameters -------------'
	write(*,*) ' cosmology: H0', h0, ' Omat', omat, ' Olam', olam
	write(*,*) '            z', redshift, ' D/Mpc', dmpc
	write(*,*) 'black hole:', bhM_Msun, ' Msun',  bhMsun_yr, ' Msun/yr'
	write(*,*) '  geometry: Nr', nr, ' Rin/Rg', risco_rg, ' Rout/c', rout, ' d'
	write(*,*) '           H/R', h1_r1_pct, '% R1', R1, ' beta', bet
	write(*,*) '     (H1-Hx)/R', dh_r1_pct, '% ', primary_dh
      if( wamp .ne. 0. ) then
	write(*,*) '   ripples: amp', wamp, ' pow', wpow
	write(*,*) '        wavenum', wnum, ' pow', wnpow, ' crest', crest
      end if
	write(*,*) '      tilt:', dinc, ' deg  (', dincmn, dincmx, ' )'
	write(*,*) '      lamp:', xlamp, ' dlogT', dlogt, ' Hx/Rg', hlamp_rg
	write(*,*) '    albedo:', albedo, ' fcol:', fcol, ' primary_xlamp', primary_xlamp
	write(*,*) '       fit: sig(lag)', syslag, 'd  sig(SED)', dexsed, ' dex'
	write(*,*) '      dust: E(B-V)', ebmv_agn
	write(*,*)
	return
	end
* -------------------------------------
	subroutine secondary
* compute secondary parameters from primaries
* In : primary parameters
*     cosmology: h0, omat, olam, redshift
*    black hole: dmpc, bhM_Msun, bhMsun_yr
* disk geometry: nr, risco_rg, rout, r1, h1_r1_pct (or dh_r1_pct), bet
*        ripple: wamp, wpow, wnum, wnpow, crest
*   inclination: dinc  mnmx=(dincmn, dincmx )
*          lamp: xlamp (or dlogt), hlamp_rg, albedo, fcol
*       fitting: dexsed, syslag
*      geometry: needgeom
*          dust: ebmv_agn
* ------------------------------------------
* Out: secondary parameters (derived from primaries)
*	dmpcz, bhRg
*	geometry ( r_r, h_r, h0_r )
*	r0 = r1 or r0 = rout (if r1=0)
*	risco, eta, rout_rg, hlamp, h1_r1, h1
*	r_rg, bhLedd_Lsun, bhL_Lsun
*	bhMsun_yr, tv1, tx1
*	cosi, sini, sigsed
*	dkappa, dlogt (or xlamp), tqin, tvmax, tvmin
* 2021 Nov Keith Horne @ Crail
* 2022 Jul KDH @ Crail - Rin replaces eta as primary
* 2022 Jul KDH @ StA - Mdot replaces L/Led as primary
* 2022-07-01 KDH @ StA - demote eta = M/2Rin to secondary
* 2022-07-05 KDH @ StA - demote L/Led = eta Mdot c^2/Led, promote Mdot
* 2022-07-05 KDH @ StA - demote h1_rg, promote h1_r1
* 2022-07-11 KDH @ StA - promote Dmpc to primary, dmpcz to secondary
* 2022-07-29 KDH @ StA - h1_r1_pct replaces h1_r1 as primary
* 2023-06-14 KDH @ StA - needgeom
* 2023-08-30 KDH @ StA - syslag
* 2023-11-15 KDH @ StA - dkappa
* 2023-12-09 KDH @ StA - dlogt, tqin
* 2023-12-10 KDH @ StA - tvmx, tvmn, tqmx, tqmn
* 2025-01-27 KDH @ StA - primary_xlamp for xlamp vs dlogt
* 2025-04-14 KDH @ StA - primary_dh for dh_r1_pct vs h1_r1_pct
*
* KDH : L_disc = eta Mdot c^2 (not yet done)

	use all
	use constants
	use cosmology
	use agn
	logical verbose/ .false. /


	if( verbose ) call showprimary
	write(*,*) '--- secondary BOWL parameters ----------'

* reference radius (light days)
	r0 = r1
	if( r0 .le. 0. ) r0 = rout

* Report for debugging
	write(*,*) 'R1', r1, ' Rout', rout, ' R0', r0

* D/Mpc (from redshift, h0, omat, olam)
      if( dmpc .le. 0. .or. dmpcz .le. 0. ) then
	dmpcz = clight / vkms / h0 * dl_zoo( redshift, omat, olam )
      end if
	if( dmpc .le. 0. ) dmpc = dmpcz

* Rg (light days)
c	bhRg = bhM_Msun * ( sunRs / clight / day )
* 2024-04-30 FIX BUG REPORTED BY JIAMU HUANG
	bhRg = bhM_Msun * ( sunRg / clight / day )

* ref radius R1 / Rg  KDH : SHOULD THIS BE R0?
	r1_rg = r0 / bhRg

* scale geometry from Rg to light days
* Risco (light days)
	risco = bhRg * risco_rg
* disc radiative efficiency
	eta = 0.5 / risco_rg
* Rout / Rg
	rout_rg = rout / bhRg
* Hlamp (light days)
	hlamp = bhRg * hlamp_rg
	hlamp_r1 = hlamp / r0
* H/R
      if( primary_dh ) then
	h1_r1_pct = hlamp_r1 * 100. + dh_r1_pct
      end if
	h1_r1 = h1_r1_pct / 100.
* disk thickness h1 (light days)
	h1 = h1_r1 * r0
	h1_rg = h1 / bhRg
      if( .not. primary_dh ) then
	dh_r1_pct = ( h1_r1 - hlamp_r1 ) * 100.
      end if

* Ledd/Lsun
	bhLedd_Lsun = bhM_Msun * ( sunLedd / sunL )

* 2022-07-05 KDH : bhMsun_yr => bhL_Lsun => bhL_Ledd
	bhL_Lsun = bhMsun_yr * eta * clight * clight * ( sunM / sunL ) / year
	bhL_Ledd = bhL_Lsun / bhLedd_Lsun
* Lbh/Lsun
c	bhL_Lsun = bhL_Ledd * bhLedd_Lsun
* Mdot/(Msun/yr)
c	bhMsun_yr = ( bhL_Lsun / eta / clight / clight ) * ( sunL / sunM ) * year


* Tv^4 = ( 3 G M Mdot ) / ( 8 pi sigma R^p )
	r = r0 * clight * day
	x = alog( sunM ) - ( 4. * alp ) * alog( r )
	t4unit = 3./8./pi * (gnewt/stefan) * exp(x) * (sunM/year)
c	t4unit = 3./8./pi * (gnewt/stefan) * (sunM/r/r/r) * (sunM/year)
	tv4 = t4unit * bhm_Msun * bhMsun_yr
* Tv1 = Tv(r0)
	tv1 = tv4 ** 0.25

      if( verbose ) then
	write(*,*) ' alp = -dlnT/dlnr', alp
	write(*,*) 't4unit', t4unit, ' Tv4', tv4, ' Tv1', tv1
	write(*,*) 'r1/rg', r1_rg, ' risco_rg', risco_rg
      end if
* Tvin = Tv(Rin)
	ref = r1_rg
	if( ref .le. 0. ) ref = rout_rg
	x = ref / risco_rg
	tvin = tv1 * x ** alp
* Tvmx = Tv(peak)   for alpha = 3/4 and qin = 0
	tvmx = tvin * ( 6. ** 1.5 / 7. ** 1.75 )
* Tvmn = Tv(Rout)
	x = ref / rout_rg
	tvmn = tv1 * x ** alp
      if( verbose ) then
	write(*,*) 'Tv(Risco)', tvin, ' Tv(max)', tvmx, ' Tv(min)', nint(tvmn)
      end if
* torque corrections
*      T^4(x) = T1^4 x^{-4a-1/2) ( sqrt(x) + (Q-1) )
*  0 = (T^4)' = T1^4 x^{-4a-3/2) ( (1-Q)(4a+1/2) - 4a sqrt(x) )
*     sqrt(x) = (1-Q) ( 1 + (1/8a) )
* Tqin = Tq(Rin)
	tqin = tvin * qin ** 0.25
* Tqmx = Tq(max)
	omq = 1. - qin
	rootx = omq * ( 1. + 1. / ( 8. * alp ) )
	rootx = max( 1., rootx )
	x = rootx * rootx
	tqmx = tvin / x ** alp * ( 1. - omq / rootx ) ** 0.25
	rqmx = x * risco
      if( verbose ) then
	write(*,*) 'qin', qin, ' Tq(Risco)', tqin, ' Tq(max)', nint(tqmx)
      end if

	t4unit = ( 1. - albedo ) * (sunL/r/r) / ( 4. * pi * stefan )

* xlamp is primary
      if( primary_xlamp ) then
* Tx^4 = (1-A) Lx / ( 4 pi sigma R^2 )
	tx4 = t4unit * xlamp * bhL_Lsun
* DlogT
	dlogt = alog10( ( tv4 + tx4 * hlamp / r0 ) ** 0.25 / tv1 )
	write(*,*) primary_xlamp, ' X', xlamp, ' => Tx4', tx4, ' => DlogT', dlogt
* Dlogt is primary
      else
	tx4 = tv4 * ( 10. ** ( 4. * dlogt ) - 1. ) * r0 / hlamp
	xlamp = tx4 / t4unit / bhL_Lsun
	write(*,*) primary_xlamp, ' DlogT', dlogt, ' => Tx4', tx4, ' => X', xlamp
      end if

* Tx(R1)  lamp flux at distance R1 (face-on irradiation)
	tx1 = tx4 ** 0.25
* delta Kappa = ( Lx Hx (1-A) / G M Mdot ) = (Hx/Rg)( Lx (1-A) / Mdot c^2 )
* Lx is bolometric luminosity of 2 lamps
	dkappa = hlamp_rg * xlamp * ( 1. - albedo )

      if( verbose ) then
	write(*,*) 't4unit', t4unit, ' Tx4', tx4, ' Tx1', tx1
	write(*,*) ' Hx/R', hlamp/r0, ' log( T / Tv )', dlogt
      end if


* KDH: CALL GEOMETRY SEPARATELY FROM SECONDARY ?
* fix invalid nr
      if( nr .le. 1 ) then
	nr = min( maxr, 500 )
	needgeom = .true.
      end if

* radius and height grids
c      if( needgeom ) then
	call geometry( nr, r_r, h_r, risco, rout, r1, h1, bet,
     &		wamp, wpow, wnum, wnpow, crest )
	needgeom = .false.
c      end if

* flat disc
      do i = 1, nr
	h0_r( i ) = 0.d0
      end do

* inclination
	radian = dtr * dinc
	cosi = cos( radian )
	sini = sin( radian )

* systematic flux error
	sigsed = 10.d0 ** dexsed - 1.d0

	return
	end

*--------------------------------
	subroutine bofplot
* MCMC iteration summary plot
* 2021 Keith Horne @ Crail
* 2021 Jul KDH @ Crail - clip to (plo,phi)
* 2021 Jul KDH @ Crail - pgenv -> pgvstd,pgqvp,pgsvp,pgswin,pgbox
* 2023 Jun KDH @ StA - Peff on plot
* 2023 Jun KDH @ StA - logy decade labels (small)
* 2024-05-09 KDH @ StA - revise Peff report to med(DBoF)/Np=x/y=z.
	use all
	use agn
	use pgplot
	logical pgagain

	write(*,*)
	write(*,*) '---- bofplot ------------'
	write(*,*) 'Fit parameters', nfit, ' of', npar
	write(*,*) 'MCMC iterations', nits
	if( npar .le. 0 ) return
	if( nfit .le. 0 ) return
	if( nits .le. 0 ) return


* nothing to do
	if( nits .le. 0 ) return
	if( nfit .le. 0 ) return
	if( npar .le. 0 ) return

  100	nx = 1
	npan = nfit + 1
	ipanb = 1
	ny = npan

* open
	call pgstart( 1, 1, nouse )
	csize = 1.
	call pgsch( csize )
	call pgbgw

* standard viewport
	call pgvstd
	call pgqvp( 0, x1, x2, y1, y2 )
	write(*,*) 'pgqvp(', x1, x2, y1, y2, ' )'


* panels
      do ipan = 1, npan
	iy = ipan
* sub-panel viewport
	dy = ( y2 - y1 ) / max( 1, ny )
	yb = y2 - dy * iy
	yt = yb + dy
	write(*,*) 'pgsvp(', x1, x2, yb, yt, ' )'
	call pgsvp( x1, x2, yb, yt )
* verify
c	call pgqvp( 0, x1, x2, y1, y2 )
c	write(*,*) 'pgqvp(', x1, x2, y1, y2, ' )'

* defaults
	logx = 0
	logy = 0

	xlabel = 'MCMC samples'
	call append_i( xlabel, nits )

	title = 'BOWL : MCMC Iterations'

	np = nits
      do i = 1, np
	plot1( i ) = i
      end do

* find best BoF and burn-in stage 1
	ibest = 1
	iburn1 = 0
      do i = 1, np
	b = bof_i( i )
* best BoF
	if( bof_i( i ) .lt. bof_i( ibest ) ) ibest = i
	bmn = bof_i( ibest )
* median of previous BoFs
	plot2( i ) = b
	bmed = xmedian( i, plot2 )
* detect end of burn-in
	if( iburn1 .eq. 0 .and. b .gt. bmed ) iburn1 = i
	plot3( i ) = bmn
	plot4( i ) = bmed
      end do
	write(*,*) 'End burn-in phase 1 at', iburn1
	write(*,*) 'Best BoF at', ibest, bof_i( ibest )

* KDH : IS THIS THE BEST THRESHOLD ?
* identify the core
	bcore = ( bmed + bmn ) / 2.
	nc = 0
      do i = 1, np
	b = plot2( i )
	p = 0.
      if( b .le. bcore ) then
	p = 1.
	nc = nc + 1
      end if
	plot5( i ) = p
      end do
	write(*,*) 'Bmed', bmed, ' Bmin', bmn, ' Bcore', bcore
	write(*,*) 'Core samples', nc, ' of', np, ' =', nc * 100. / np, ' %'

* BoF panel
      if( ipan .eq. ipanb ) then
	ylabel = 'BoF'
* log BoF
	i1 = max( 1, iburn1 )
	n1 = np - i1 + 1
	call mnmx( n1, plot2( i1 ), bmn, bmx )
	logy = 0
      if( bmx .gt. bmn ) then
      if( bmn .gt. 0. .and. bmx / bmn .gt. 10. ) then
	logy = 1
      do i = 1, np
	plot2( i ) = alog10( plot2( i ) )
      end do
      end if
      end if

* parameter panels
      else
	ifit = ipan
	if( ipanb .eq. 1 ) ifit = ipan - 1
	ip = i_p( ifit )
	name = name_p( ip )
	symbol = symbol_p( ip )
	units = units_p( ip )
	dp = dp_p( ip )
	pbest = p_pi( ip, ibest )
	pavg = avg_p( ip )
	prms = rms_p( ip )
	pmed = pmed_p( ip )
	pmad = pmad_p( ip )
	keep = keep_p( ip )
	logy = 0
	if( log_p( ip ) ) logy = 1
	ylabel = symbol
	ybot = bot_p( ip )
	ytop = top_p( ip )
      do i = 1, np
	plot2( i ) =  p_pi( ip, i )
	plot3( i ) = dp_pi( ip, i )
      end do
* last panel option
      end if

	write(*,*) 'Panel', ipan, ' ', ylabel(:len1(ylabel))

* window
	call mnmx( np, plot1, xmn, xmx )
* omit first 20%
	i1 = max( 1, nint( np * 0.2 ) )
	n1 = np - i1 + 1
	call mnmx( n1, plot2( i1 ), ymn, ymx )

	xmn = xmn - 0.5
	xmx = xmx + 0.5
	add = ( ymx - ymn ) * 0.05
	if( add .eq. 0. ) add = 0.1
	ymn = ymn - add
	ymx = ymx + add

      if( logy .ne. 0 ) then
	y = ( ymx + ymn ) / 2.
	add = 0.5
	ymx = max( ymx, y + add )
	ymn = min( ymn, y - add )
      end if


* clip to parameter domain
      if( ipan .ne. ipanb ) then
	ymn = max( ymn, ybot )
	ymx = min( ymx, ytop )
      end if

	yspan = abs( ymx - ymn )
	xspan = abs( xmx - xmn )
	logxy = 0
	if( logx .ne. 0 .and. xspan .le. 10. ) logxy = logxy + 10 * logx
	if( logy .ne. 0 .and. yspan .le. 10. ) logxy = logxy + 20 * logy
	write(*,*) 'logx', logx, ' logy', logy, ' logxy', logxy

* replace pgenv to eliminate white space between panels
c	write(*,*) 'pgenv(', xmn, xmx, ymn, ymx, logxy, ' )'
c	call pgenv( xmn, xmx, ymn, ymx, 0, logxy )

* window
	write(*,*) 'pgswin(', xmn, xmx, ymn, ymx, ' )'
	call pgswin( xmn, xmx, ymn, ymx )


* burn-in
	write(*,*) 'Burn-in at sample', iburn1
      if( iburn1 .ge. 1 .and. iburn1 .le. np ) then
	call pgsci( korange )
	call pgsls( ldash )
	x = plot1( iburn1 )
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
	call pgsls( 1 )
	call pgsci( kblack )
      end if



* y label
	label = 'BCST'
	if( nfit .le. 5 ) label = label(:len1(label)) // 'N'
	if( logxy .ge. 20 .and. logy .gt. 0 )
     &		label = label(:len1(label)) // 'L'
	call pgbox( ' ', 0., 0, label, 0., 0 )

* logy decade labels (small)
      if( logy .gt. 0 ) then
	s = 1.
	call pgsch( 0.7 * csize )
	call pgbox( ' ', 0., 0, 'NLP', s, 0 )
	call pgsch( csize )
      end if

* x label
	label = 'BCST'
	if( iy .eq. ny ) label = 'BCNST'
	call pgbox( label, 0., 0, ' ', 0., 0 )

* title
      if( ipan .eq. ipanb ) then
	title = 'BOWL fit to'
	call append( title, name_a( iagn ) )
      else
* parameter avg and rms
	title = ylabel
	call append( title, '=' )
      if( logy .ne. 0 .and. logxy .le. 20 ) then
	call append( title, '10\u' )
	call append_datsig( title, pavg, prms, 2 )
	if( logxy .le. 20 ) call append( title, '\d' )
      else
	call append_datsig( title, pavg, prms, 2 )
      end if
	write(*,*) title(:len1(title))

* med(mad) label
	label = '10\umed(mad)\d='
	label = '='
      if( logy .eq. 0 ) then
	call append_datsig( label, pmed, pmad, 2 )
      else
	call append( label, '10\u' )
	call append_datsig( label, pmed, pmad, 2 )
	call append( label, '\d' )
      end if

* logy treatment
      if( logy .ne. 0 ) then
	p = 10. ** pavg
	s = p * sqrt( ( 10. ** prms - 1. ) * ( 1. - 10. ** ( -prms ) ) )
	call append( title, '=' )
	call append_datsig( title, p, s, 2 )
	call nospace( title )
	p = 10. ** pmed
	s = p * sqrt( ( 10. ** prms - 1. ) * ( 1. - 10. ** ( -prms ) ) )
	call append( label, '=' )
	call append_datsig( label, p, s, 2 )
      end if
	call nospace( label )
	x = 0.9
	x = 0.05
	call pgsci( kblue )
	call pgmtxt( 'b', -0.5, x, 0., label )
	write(*,*) label(:len1(label)) // ' ' // label(:len1(label))
	call pgsci( kblack )

* best label (red)
	label = ' '
	p = pbest
	s = pmad
      if( logy .eq. 0 ) then
	call append_datsig( label, p, s, 2 )
      else
	p = 10. ** p
	s = p * sqrt( ( 10. ** s - 1. ) * ( 1. - 10. ** ( -s ) ) )
	call append_datsig( label, p, s, 2 )
      end if
	call pgsci( kred )
	call pgmtxt( 'T', -1.5, 0.95, 1., label )
	call pgsci( kblack )

      end if


* x label
	if( ipan .eq. npan ) call pglabel( xlabel, ' ', ' ' )
* y label
	call pglabel( ' ', ylabel, ' ' )
* title
      if( ipan .eq. 1    ) then
	call pglabel( ' ', ' ', title )
      else
	call pgmtxt( 'T', -1.5, x, x, title )
      end if


* bof panel
      if( ipan .eq. ipanb ) then

* effective number of parameters
	ybest = plot2( ibest )
	ymed = plot4( np )
	peff = ymed - ybest
	label = 'med(\gDBoF)/N\dp\u='
	call append_f( label, peff, '(f15.2)' )
	call append( label, '/' )
	n = npan - 1
	call append_i( label, n )
	call append( label, '=' )
	call append_n( label, peff / n, 3 )
	call nospace( label )
	call pgmtxt( 'T', 0.5, 1., 1., label )

* min BoF
	call pgsci( kred )
	call pgsls( 4 )
	call mnmx( np, plot2, bmn, bmx )
	y = bmn
	write(*,*) 'BoF min', y
	call pgsls( 2 )
	call pgmove( xmn, y )
	call pgdraw( xmx, y )
	call pgsls( 1 )


	call pgsci( kblack )

* parameter panels
      else
* mean and rms
	call pgsls( 4 )
      do k = -1, 1
	y = pavg + k * prms
	call pgmove( xmn, y )
	call pgdraw( xmx, y )
      end do
* med and mad
	call pgsls( 2 )
	call pgsci( kblue )
      do k = -1, 1
	y = pmed + k * pmad
	call pgmove( xmn, y )
	call pgdraw( xmx, y )
      end do
	call pgsls( 1 )

* keep percentage
	label = 'keep'
	label = ' '
	n = nint( keep * 100. / nits )
	call append_i( label, n )
	call append( label, '%' )
c	call pgmtxt( 'T', 0.5, 1., 1., label )
	x = 0.95
	call pgsci( korange )
	call nospace( label )
	call pgmtxt( 'B', -0.5, 1., 0., ' ' // label )

* end xtras
      end if

* best
	call pgsci( kred )
	x = plot1( ibest )
	y = plot2( ibest )
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
	call pgsch( 2. * csize )
	call pgpoint( 1, x, y, istar )
	call pgsch( csize )

* minimum BoF
      if( ipan .eq. ipanb ) then
	call pgsci( kred )
	call pgsls( 2 )
	call pgbin( np, plot1, plot3, .true. )
	call pgsls( 1 )

* median of previous BoF
	call pgsci( korange )
	call pgsls( ldashdot )
	call pgbin( np, plot1, plot4, .true. )

* end of burnin phase 1
      if( iburn1 .gt. 0 ) then
	x = plot1( iburn1 )
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
      end if
	call pgsls( 1 )
      end if

* bin histogram
	call pgsci( kblue )
	call pgbin( np, plot1, plot2, .true. )

* core segments in red
	call pgsci( kred )
	i1 = np
	i2 = 0
      do i = 1, np
      if( plot5( i ) .gt. 0.5 ) then
	i1 = min( i1, i )
	i2 = max( i2, i )
      else if( i2 .ge. i1 ) then
	n = i2 - i1 + 1
	call pgbin( n, plot1( i1 ), plot2( i1 ), .true. )
c	write(*,*) 'Core segment', i1, i2, n
	i1 = np
	i2 = 0
      end if
      end do
	n = i2 - i1 + 1
      if( n .gt. 0 ) then
	call pgbin( n, plot1( i1 ), plot2( i1 ), .true. )
c	write(*,*) 'Core segment', i1, i2, n
      end if

* parameter steps
      if( ipan .ne. ipanb ) then
	call pgswin( xmn, xmx, 0., ymx - ymn )
	call pgsci( korange )
	call pgbin( np, plot1, plot3, .true. )
	call pgswin( xmn, xmx, ymn, ymx )
      end if
	call pgsci( kblack )

* next panel ipan
      end do


	call pgstop
	if( pgagain() ) goto 100

	return
	end
*--------------------------------
	subroutine mcmcreport
* summary of MCMC iterations
* 2021 Keith Horne @ Crail
	use all

	write(*,*)
	write(*,*) 'MCMC iterations', nits

      if( nits .gt. 0 ) then
* skip burn-in
	pburn = 0.2
	j1 = nits * pburn
	write(*,*) 'Burn-in', j1, ' iterations', j1 * 100. / nits, ' %'

* iteration summary
	n = 0
	write(*,*) '     iteration', '   BoF',
     &	'          BoFmin', '         BoF-BoFmin'
      do i = 1, nits
	b = bof_i( i )
	if( i .eq. 1 ) bmn = b
	bmn = min( bmn, b )
      if( i .ge. nits-5 .or. mod( i, max( 1, nits/5 ) ) .eq. 0 ) then
	write(*,*) i, b, bmn, b - bmn
      end if
      if( i .gt. j1 ) then
	n = n + 1
	plot1( n ) = b
      end if
      end do
	call medmad( n, plot1, bmed, bmad )
	write(*,*) 'BoF med', bmed, ' min', bmn
	write(*,*) '    med - min', bmed - bmn, ' Npar', npar

* stats
	write(*,*)
	write(*,*) 'Parameter Statistics:'
c      do ifit = 1, nfit
c	i = i_p( ifit )
      do i = 1, npar
	name = name_p( i )
	code = code_p( i )
	n = 0
      do j = j1+1, nits
	n = n + 1
	plot1( n ) = p_pi( i, j )
      end do
	call medmad( n, plot1, pmed, pmad )
	call avgrms( n, plot1, avg, rms )
	avg_p( i ) = avg
	rms_p( i ) = rms
	pmed_p( i ) = pmed
	pmad_p( i ) = pmad
	p = avg
	if( log_p( i ) ) p = 10. ** p
	pkeep = keep_p( i ) * 100. / nits
	write(*,*) i, ' ', code, ' -------- ', name(:len1(name)), ' --------', p
	write(*,*) '    avg', avg, ' med', pmed, ' keep', pkeep, ' %'
	write(*,*) '    rms', rms, ' mad', pmad, ' step', dp_p( i )
* next parameter i
      end do
      end if
	write(*,*)

	return
	end

* -------------------------------------
	subroutine cornerplot
* corner plot of MCMC parameters
* 2021 Nov Keith Horne @ Crail
* 2021 Jul KDH @ Crail - clip to (plo,phi)
* 2021 Jul KDH @ Crail - i_p( nfit )
* 2022 Aug KDH @ Crail - pggray
* 2022 Aug KDH @ Crail - contours
* 2022 Aug KDH @ StA - report fixed parameters
* 2023 Jun KDH @ StA - small x and y tick labels
* 2023 Jun KDH @ StA - codex, codey
* 2023 Jun KDH @ StA - expected correlations
* 2024 Jan KDH @ StA - fcol correlations with Mdot, Hx, Lx
* 2025 Mar 17 KDH @ StA - add Nr to legend
	use all
	use agn
	use constants
	use pgplot
	logical pgagain
	logical yes

	character*50 codex, codey

	write(*,*)
	write(*,*) '---- cornerplot ------------'
	write(*,*) 'Fit parameters', nfit, ' of', npar
	write(*,*) 'MCMC iterations', nits
	if( npar .le. 0 ) return
	if( nfit .le. 0 ) return
	if( nits .le. 0 ) return

  100	nx = 1
	ny = 1

* open
	call pgstart( nx, ny, new )
	csize = 1.
	call pgsch( csize )
	call pgbgw

* standard viewport
	call pgvstd
	call pgqvp( 0, x1, x2, y1, y2 )
	write(*,*) 'pgqvp(', x1, x2, y1, y2, ' )'

	burnin = 0.2
	i1 = max( 1, nint( nits * burnin ) )
	n = nits - i1 + 1

	ibest = 1
	bmn = bof_i( 1 )
      do i = 1, nits
	b = bof_i( i )
	if( b .lt. bmn ) ibest = i
	bmn = min( b, bmn )
      end do
	write(*,*) 'Best BoF at', ibest, bof_i( ibest )

* y panel
c      do iy = 1, npar
      do iy = 1, nfit
	iyp = i_p( iy )

* x panel
      do ix = 1, iy
	ixp = i_p( ix )

* viewport
	dx = ( x2 - x1 ) / nfit
	xr = x1 + ix * dx
	xl = xr - dx
	dy = ( y2 - y1 ) / nfit
	yb = y2 - iy * dy
	yt = yb + dy
	write(*,*) 'pgsvp(', xl, xr, yb, yt, ' )'
	call pgsvp( xl, xr, yb, yt )

* x values
	xmed = pmed_p( ixp )
	xmad = pmad_p( ixp )
	xavg = avg_p( ixp )
	xrms = rms_p( ixp )
	xsig =  dp_p( ixp )
	xbot = bot_p( ixp )
	xtop = top_p( ixp )
      do i = 1, nits
	plot1( i ) = p_pi( ixp, i )
      end do
	logx = 0
	if( log_p( ixp ) ) logx = 1
	name = name_p( ixp )
	code = code_p( ixp )
	symbol = symbol_p( ixp )
	units = units_p( ixp )
	units = units_p( ixp )
	write(*,*) 'X avg rms', xavg, xrms, ' ', name(:len1(name))
	write(*,*) '  med mad', xmed, xmad, ' sig', xsig

	codex = code

* y values
	ymed = pmed_p( iyp )
	ymad = pmad_p( iyp )
	yavg = avg_p( iyp )
	yrms = rms_p( iyp )
	ysig =  dp_p( iyp )
	ybot = bot_p( iyp )
	ytop = top_p( iyp )
      do i = 1, nits
	plot2( i ) = p_pi( iyp, i )
      end do
	logy = 0
	if( log_p( iyp ) ) logy = 1
	name = name_p( iyp )
	code = code_p( iyp )
	write(*,*) 'Y avg rms', yavg, yrms, ' ', name(:len1(name))
	write(*,*) '  med mad', ymed, ymad, ' sig', ysig

	codey = code

* x window
	wide = 4.
	add = xmad * wide
	add = xrms * wide
	if( add .le. 0. ) add = 0.1
	x = xmed
	x = xavg
	xmn = x - add
	xmx = x + add

* clip x to parameter domain
	xmn = max( xmn, xbot )
	xmx = min( xmx, xtop )

* y window
      if( ix .ne. iy ) then
	add = ymad * wide
	add = yrms * wide
	if( add .le. 0. ) add = 0.1
	y = ymed
	y = yavg
	ymn = y - add
	ymx = y + add

* clip y to parameter domain
	ymn = max( ymn, ybot )
	ymx = min( ymx, ytop )


* binned histogram
      else
* bins
	nbin = 50
	nbin = 31
      do i = 1, nbin
* bin centre from 0.5 to nbin-0.5
	b = i - 0.5
* distribute bins from xmn to xmx
	p = b / max( 1, nbin )
	plot1( i ) = xmn * ( 1. - p ) + xmx * p
* zero bin counts
	plot3( i ) = 0.
      end do
* sort
	call quicksort( nits, plot2, iplot )
	add = 1.
      do ir = 1, nits
	i = iplot( ir )
	call lookup( nbin, plot1, plot2(i), j1, j2, p )
	b = j1 * ( 1. - p ) + p * j2
	j = nint( b )
	plot3( j ) = plot3( j ) + add
      end do
	call mnmx( nbin, plot3, ymn, ymx )
	write(*,*) 'Bins', nbin, ' counts', ymn, ymx
* normalise
	div = ymx * 1.1
      do i = 1, nbin
	plot3( i ) = plot3( i ) / div
      end do
	ymn = 0.
	ymx = 1.
      end if

	write(*,*) 'pgswin(', xmn, xmx, ymn, ymx, ' )'
	call pgswin( xmn, xmx, ymn, ymx )


* 2d greyscale image
* KDH : do this here so that pgbox is called later
	nitmn = 100
      if( ix .ne. iy .and. nits .lt. nitmn ) then
	write(*,*) 'Have', nits, ' need', nitmn, ' iterations for greyscale.'
      else if( ix .ne. iy ) then

	n = sqrt( max( 1., nits / 10. ) )
	n = max( 9, min( 31, n ) )
	n2d = n * n
      do i = 1, n2d
	pic( i ) = 0.
      end do
      do i = 1, nits
	x = p_pi( ixp, i )
	y = p_pi( iyp, i )
	px = ( x - xmn ) / ( xmx - xmn )
	py = ( y - ymn ) / ( ymx - ymn )
	jx = nint ( px * n + 0.5 )
	jy = nint ( py * n + 0.5 )
      if(     jx .ge. 1 .and. jx .le. n
     &	.and. jy .ge. 1 .and. jy .le. n ) then
	j = jx + n * ( jy - 1 )
	pic( j ) = pic( j ) + 1.
      end if
      if( i .lt. 10 ) then
	write(*,*) i, x, y, jx, jy
      end if
      end do
	call mnmx( n2d, pic, pmn, pmx )
	write(*,*) 'n', n, ' n2d', n2d, ' pmn,mx', pmn, pmx
	dx = ( xmx - xmn ) / n
	dy = ( ymx - ymn ) / n
	tr( 6 ) = dy
	tr( 5 ) = 0.
	tr( 4 ) = ymn - dy / 2.
	tr( 3 ) = 0.
	tr( 2 ) = dx
	tr( 1 ) = xmn - dx / 2.
	write(*,*) ' tr', tr
	bg = 0.
	fg = pmx * 1.5
	write(*,*) 'pggray(', n, n, fg, bg, ' )'
	call pggray( pic, n,n, 1,n, 1,n, fg, bg, tr )

* enclosed probability levels
	clev( 1 ) = 0.68
	clev( 2 ) = 0.95
	nc = 2
* require sufficient counts per pixel
* KDH : threshold could depend on nits
	p = nc * 2.
	write(*,*) 'Have', nint(pmx), ' need', nint(p), ' samples/pixel for contours.'
      if( pmx .ge. p ) then
	write(*,*) 'quicksort(', n2d, ' )'
* sort by increasing MCMC samples / pixel
	call quicksort( n2d, pic, i2d )
* total samples should be nits
	sum = 0.d0
      do i = 1, n2d
	sum = sum + pic( i )
      end do
      if( nint( sum ) .ne. nits ) then
	write(*,*) '*** WARNING: Sum', sum, ' Nits', nits, ' NOT THE SAME.'
      end if
* contour levels
      do ic = 1, nc
	if( ic .eq. 1 ) call pgsci( kblack )
	if( ic .eq. 2 ) call pgsci( kblue )
	if( ic .eq. 3 ) call pgsci( kgreen )
* enclose fraction p of nits
	p = clev( ic )
	e = p * nits
	write(*,*) nint(p*100), ' % contour', nint(e), ' of', nits, ' samples'
* accumulate sorted pixels (high density first)
	sum = 0.d0
	old = pmx
      do ir = n2d, 1, -1
      if( sum .le. e ) then
	add = pic( i2d( ir ) )
	sum = sum + add
      if( sum .gt. e ) then
* report
	write(*,*) nint(p*100), ' % contour', nint(sum), ' of', nits, ' samples'
	write(*,*) ' between density', add, ' and', old, ' samples / pixel'
* 2 contours bracketing the confidence interval
c	call pgcont( pic, n, n, 1, n, 1, n, old, 1, tr )
c	call pgcont( pic, n, n, 1, n, 1, n, add, 1, tr )
* 1 contour at mean
	c = ( old + add ) / 2.
	call pgcont( pic, n, n, 1, n, 1, n, c, 1, tr )
	goto 50
      end if
	old = add
* next pixel
      end if
      end do
   50	write(*,*) 'Sum', sum, ' target', e
	call pgsci( kblack )
* next contour ic
      end do
* end contouring
      end if
* end greyscale
      end if


	xspan = abs( xmx - xmn )
	yspan = abs( ymx - ymn )
c	big = 3.
* x axis
	xlabel = 'bcst'
	nup = 5
	if( iy .eq. nfit .and. nfit .le. nup ) xlabel = 'n' // xlabel
* y axis
	ylabel = 'bcst'
	if( ix .eq. 1 .and. nfit .le. nup )    ylabel = 'n' // ylabel
	call pgbox( ' ', 0., 0, ylabel, 0., 0 )
	write(*,*) 'pgbox( ', xlabel(:len1(xlabel)),
     &		' ', ylabel(:len1(ylabel)), ' )'
	call pgbox( xlabel, 0., 0, ylabel, 0., 0 )


* x tick labels (small)
      if( nfit .gt. nup ) then
	call pgsch( 0.7 * csize )
      if( iy .eq. nfit ) then
	span = xspan
	s = 1
	if( span .lt. 1.5 ) s = 0.5
	if( span .lt. 0.6 ) s = 0.2
	if( span .lt. 0.3 ) s = 0.1
	if( span .lt. 0.15 ) s = 0.05
	call pgbox( 'TNP', s, 0, ' ', 0., 0 )
      end if
* y tick labels (small)
      if( ix .eq. 1 ) then
	span = yspan
	s = 1
	if( span .lt. 1.5 ) s = 0.5
	if( span .lt. 0.6 ) s = 0.2
	if( span .lt. 0.3 ) s = 0.1
	if( span .lt. 0.15 ) s = 0.05
	call pgbox( ' ', 0., 0, 'TNP', s, 0 )
      end if
	call pgsch( csize )
      end if



* x label
	label = name_p( ixp )
	label = symbol_p( ixp )
	yes = logx .ne. 0
	if( yes ) label = 'log(' // label(:len1(label)) // ')'
	xlabel = ' '
      if( iy .eq. nfit ) then
	xlabel = label
      end if
* y label
	ylabel = ' '
      if( ix .eq. 1 .and. iy .gt. ix ) then
	ylabel = name_p( iyp )
	ylabel = symbol_p( iyp )
	yes = logy .ne. 0
	if( yes ) ylabel = 'log(' // ylabel(:len1(ylabel)) // ')'
      end if
	call pglabel( xlabel, ylabel, ' ' )
* title
	title = ' '
      if( ix .eq. iy ) then
	title = symbol_p( ixp )
	if( yes ) title = 'log'
	call append( title, '=' )
	l = len1( title ) - 1
	call append_datsig( title, xavg, xrms, 2 )
	call nospace( title(l:) )
	if( .not. yes ) title = title (:len1(title)) // units

	y = 0.5
	if( logx .ne. 0 ) y = 2.
	xt = 0.5
	if( nfit .gt. 3 ) xt = 0.0
	call pgmtxt( 't', y, xt, xt, ' ' // title )

* anti-logged results
      if( logx .ne. 0 ) then
	x = 10. ** xavg
	s1 = 10. ** xrms  - 1.
	s2 = 1. - 10. ** ( - xrms )
	s = x * sqrt( s1 * s2 )
	write(*,*) x, ' +/-', s
	write(*,*) xavg, ' +/-', xrms, ' =>', x, '+/-', s
	title = ' '
	call append( title, symbol_p( ixp ) )
      if( name_p( ixp ) .eq. 'Mdot' ) then
	call append( title, '<' )
      else
	call append( title, '=' )
      end if
	l = len1( title ) - 1
	call append_datsig( title, x, s, 2 )
	call nospace( title(l:) )
	title = title(:len1(title)) // units(:len1(units))
c	call append( title, units )
	call pgmtxt( 't', 0.5, xt, xt, ' ' // title )
      end if
* cosi
      if( index( name_p( ixp ), 'cos' ) .gt. 0 ) then
	cosi = p_p( ixp )
	sini = sqrt( 1. - cosi * cosi )
	dinc = atan2( sini, cosi ) * 180. / pi
	write(*,*) 'Cos(i)', cosi, ' Sin(i)', sini, ' i', dinc
* confidence limits
	c = xavg
	s = sqrt( 1. - c * c )
	d = atan2( s, c ) * 180. / pi
	e = xrms / abs( s ) * 180. / pi
	title = 'i='
c	call append_i( title, nint( dinc ) )
c	call append_datsig( title, d, e, 2 )
	call append_i( title, nint( d ) )
	call append( title, '(' )
	call append_i( title, nint( e ) )
	call append( title, ')' )
	call append( title, '\uo\d' )
	call nospace( title )
	call pgmtxt( 't', 0.5, xt, xt, ' ' // title )
      end if
      end if


* horizontal
	call pgsls( 4 )
      if( ix .ne. iy ) then
      do k = -1, 1
	y = ymed + k * ymad
	y = yavg + k * yrms
	call pgmove( xmn, y )
	call pgdraw( xmx, y )
      end do
      end if
* vertical
      do k = -1, 1
	x = xmed + k * xmad
	x = xavg + k * xrms
	call pgmove( x, ymn )
	call pgdraw( x, ymx )
      end do
	call pgsls( 1 )

* MCMC posterior
      if( ix .eq. iy ) then

	call pgsci( kblue )
	call pgbin( nbin, plot1, plot3, .true. )
* staircase cumulative
	call pgsci( kred )
	x = xmn
	y = 0.
	add = ( ymx - ymn ) / max(1,nits)
	call pgmove( x, y )
      do ir = 1, nits
	i = iplot( ir )
* step right
	x = plot2( i )
	call pgdraw( x, y )
* mark best with a star
      if( i .eq. ibest ) then
	call pgsch( 2. * csize )
	call pgpoint( 1, x, y + add/2., istar )
	call pgsch( csize )
	call pgmove( x, y )
      end if
* step up
	y = y + add
	call pgdraw( x, y )
      end do
* last point
	x = xmx
	call pgdraw( x, y )
	call pgsci( kblack )

      else

* MCMC trajectory
	nmax = 500
	nmax = 200
      if( nits .le. nmax ) then
* trajectory in fifths
	i1 = max( 1, nint( nits * 0.2 ) )
	i2 = max( 1, nint( nits * 0.4 ) )
	i3 = max( 1, nint( nits * 0.6 ) )
	i4 = max( 1, nint( nits * 0.8 ) )
	call pgsci( korange )
	call pgline( i1, plot1, plot2 )
	call pgsci( kgreen )
	n = i2 - i1 + 1
	call pgline( n, plot1( i1 ), plot2( i1 ) )
	call pgsci( kblue )
	n = i3 - i2 + 1
	call pgline( n, plot1( i2 ), plot2( i2 ) )
	call pgsci( kred )
	n = i4 - i3 + 1
	call pgline( n, plot1( i3 ), plot2( i3 ) )
	call pgsci( kblack )
	n = nits - i4 + 1
	call pgline( n, plot1( i4 ), plot2( i4 ) )
      end if

* last point
	x = plot1( nits )
	y = plot2( nits )
	call pgpoint( 1, x, y, 23 )

* best
	x = plot1( ibest )
	y = plot2( ibest )
	call pgsci( kred )
	call pgsch( csize * 2. )
	call pgpoint( 1, x, y, istar )
	call pgsch( csize )
	call pgsci( kblack )

* expected power-law slopes
      if( logx .ne. 0 .and. logy .ne. 0 ) then
	call pgsci( kred )
	call pgsls( 2 )
	slope = 0.

* fcol vs Mdot
	if( codex .eq. 'A' .and. codey .eq. '^' ) slope = 0.25
	if( codex .eq. '^' .and. codey .eq. 'A' ) slope = 4.
      IF( .FALSE. ) THEN
* fcol vs Hx
	if( codex .eq. 'H' .and. codey .eq. '^' ) slope = 0.25
	if( codex .eq. '^' .and. codey .eq. 'H' ) slope = 4.
* fcol vs Lx
      if( primary_xlamp ) then
	if( codex .eq. 'X' .and. codey .eq. '^' ) slope = 0.25
	if( codex .eq. '^' .and. codey .eq. 'X' ) slope = 4.
      end if
      END IF

* Lx vs Hx
      if( primary_xlamp ) then
	if( codex .eq. 'X' .and. codey .eq. '*' ) slope = -1.
	if( codex .eq. '*' .and. codey .eq. 'X' ) slope = -1.
* Lx vs H/R
	if( codex .eq. 'X' .and. codey .eq. 'H' ) slope = -1.
	if( codex .eq. 'H' .and. codey .eq. 'X' ) slope = -1.
      end if
* H/R vs R
	if( codex .eq. 'H' .and. codey .eq. 'O' ) slope = -1.
	if( codex .eq. 'O' .and. codey .eq. 'H' ) slope = -1.
* M vs Mdot
	if( codex .eq. 'M' .and. codey .eq. 'A' ) slope = -1.
	if( codex .eq. 'A' .and. codey .eq. 'M' ) slope = -1.
* H/R vs Hx
	if( codex .eq. 'H' .and. codey .eq. '*' ) slope = 1.
	if( codex .eq. '*' .and. codey .eq. 'H' ) slope = 1.

* D vs M
	if( codex .eq. 'D' .and. codey .eq. 'M' ) slope = 3.
	if( codex .eq. 'M' .and. codey .eq. 'D' ) slope = 1./3.
* D vs Mdot
	if( codex .eq. 'D' .and. codey .eq. 'A' ) slope = 3.
	if( codex .eq. 'A' .and. codey .eq. 'D' ) slope = 1./3.
* D vs cos(i)
c	if( codex .eq. 'D' .and. codey .eq. 'C' ) slope = 2.
c	if( codex .eq. 'C' .and. codey .eq. 'D' ) slope = 0.5

* draw the power-law slope
      if( abs( slope ) .gt. 0.1 ) then
	x = plot1( ibest )
	y = plot2( ibest ) + slope * ( xmn - x )
	call pgmove( xmn, y )
	y = plot2( ibest ) + slope * ( xmx - x )
	call pgdraw( xmx, y )
* label the slope
	label = '+'
	if( slope .lt. 0. ) label = '-'
	i = nint( slope )
	ii = nint( 1. / slope )
      if( abs( slope - i ) .lt. 0.01 ) then
	call append_i( label, iabs( i ) )
      else if( abs( 1/slope - ii ) .lt. 0.01 ) then
	call append( label, '1/' )
	call append_i( label, iabs( ii ) )
      else
	call append_n( label, slope, 3 )
      end if
	call nospace( label )
      if( len1( label ) .gt. 1 ) then
	call pgsch( 0.7 * csize )
      if( slope .gt. 0. ) then
	call pgmtxt( 'T', -1.5, 0.25, 0.5, label )
      else
	call pgmtxt( 'T', -1.5, 0.75, 0.5, label )
      end if
	call pgsch( csize )
      end if

      end if


	call pgsls( 1 )
	call pgsci( kblack )
      end if

* step sigmas

	p = 0.15
	x = xmn * ( 1. - p ) + p * xmx
	y = ymn * ( 1. - p ) + p * ymx
	write(*,*) 'Steps:', x, xsig, y, ysig
	call pgsci( korange )
	call pgmove( x + xsig, y )
	call pgdraw( x, y )
	call pgdraw( x, y + ysig )
	call pgsci( kblack )

      end if

* next ix
      end do
* next iy
      end do

* report fixed parameters
	call pgsvp( x1, x2, y1, y2 )
	x = 1.
	y = 2.
	dy = - 1.5
c	call pgmtxt( 'T', y, x, x, 'Fixed parameters:' )
	title = 'BOWL fit to'
	call append( title, name_a( iagn ) )
	call pgmtxt( 'T', y, x, x, title )

* iterations
	title = 'MCMC'
	call append_i( title, nits )
	call append( title, 'N\dr\u' )
	call append_i( title, nr )
	y = y + dy
	call pgmtxt( 'T', y, x, x, title )

* fixed parameters
      do i = 1, npar
      if( .not. fit_p( i ) ) then
	y = y + dy
	title = symbol_p( i )
	avg = avg_p( i )
	rms = rms_p( i )

* log parameter
      if( log_p( i ) ) then
      if( abs( avg ) .gt. 5. ) then
	title = 'log(' // title(:len1(title)) // ')'
      else
	avg = 10. ** avg
	s1 = 10. ** rms  - 1.
	s2 = 1. - 10. ** ( - rms )
	rms = avg * sqrt( s1 * s2 )
      end if
      end if

      if( name_p( i ) .eq. 'Mdot' ) then
	call append( title, '<' )
      else
	call append( title, '=' )
      end if

* sig figs
	nsf = 2
	if( rms .le. 0. ) nsf = 3
* avg and rms
	l = len1( title ) - 1
	call append_datsig( title, avg, rms, nsf )

* trim fixed number
      if( rms .le. 0. ) then
* locate . at ld
	ld = index( title(l:), '.' )
	ld = l + ld - 1
* . present
      if( ld .gt. l ) then
* no E present
	le = index( title(l:), 'E' )
      if( le .le. 0 ) then
* snip tailing 0 and .
	call tidy( title )
* E present at le > ld
      else
	le = l + le - 1
* snip tailing 0 between . and E
      do k = le - 1, ld + 1, -1
	if( title(k:k) .ne. '0' ) exit
	title = title(:k-1) // title(k+1:)
	le = le - 1
      end do
* snip . before E
	if( le .eq. ld + 1) title = title(:ld-1) // title(le:)
      end if
      end if
      end if

	call nospace( title(l:) )
* units inside : log(val/units)
      if( index( title, 'log(' ) .gt. 0 ) then
	l = index( title, ')' )
	lu = len1( units_p( i ) )
	if( lu .gt. 0 ) title = title(:l-1) // '/' // units_p( i )(:lu) // title(l:)
	call nospace( title )
      else
	call append( title, units_p( i ) )
      end if
* report both cosi and i
      if( index( name_p( i ), 'cos' ) .gt. 0 ) then
	cosi = p_p( i )
	sini = sqrt( 1. - cosi * cosi )
	dinc = atan2( sini, cosi ) * 180. / pi
	write(*,*) 'Cos(i)', cosi, ' Sin(i)', sini, ' i', dinc
	l = len1( title )
	call append( title, 'i=' )
	call append_f( title, dinc, '(f15.1)' )
	call noeqbl( title )
	call tidy( title )
	call append( title, '\uo\d' )
	call nospace( title(l+2:) )
      end if
	call pgmtxt( 'T', y, x, x, title )
* next fixed parameter
      end if
      end do

* KDH : NEED TO STOW SECONDARY QUANTITIES
* secondary quantities:
*  Lx Hx
*  Tv(Rout)
*  Tmx
* Xlamp or DlogT



* close
	call pgstop
* again
	if( pgagain() ) goto 100

	return
	end

*------------------------------------------
	subroutine mcmcsave
* save MCMC iterations
* 2021 Nov Keith Horne @ Crail
* 2023 Dec KDH @ StA - write out all npar parameters, not just nfit
	use all

	character*75 savefile/ 'bowl_mcmc.dat' /

	write(*,*)
	write(*,*) '--- mcmcsave ------------'

	write(*,*) 'MCMC iterations:', nits
	if( nits .le. 0 ) return
	write(*,*) 'MCMC parameters:', nfit, ' of', npar
	if( npar .le. 0 ) return
	if( nfit .le. 0 ) return


	savefile = 'bowl_mcmc_'
	call append_i( savefile, nits )
	call append( savefile, '.dat' )
	call nospace( savefile )

* select save file
	write(*,*)
	write(*,*) 'Extant MCMC files:'
	call spawn( 'ls -l *mcmc*' )
	write(*,*)
	write(*,*) 'Avoid accidental over-write.'
	write(*,*)
	write(*,*) 'Enter none to abort.'
	file = savefile
	call inq( 'Save file', savefile, file )

* abort
	if( file .eq. 'none' .or. file .eq. 'NONE' ) return

	write(*,*)
	call spawn( 'ls -l ' // file(:len1(file)) )
	call spawn( 'head ' // file(:len1(file)) )

	write(*,*)

* open file
	iu = 50
	call openfile( iu, file, 'formatted', 'unknown', ierr )
      if( ierr .ne. 0 ) then
	write(*,*) '** OPENFILE FAILED.'
	return
      end if


* parameter details
	write( iu, * ) npar, ' MCMC parameters'
      do i = 1, npar
c      do ifit = 1, nfit
c	i = i_p( ifit )
	avg = avg_p( i )
	rms = rms_p( i )
	name = name_p( i )
	code = code_p( i )
	if( log_p( i ) ) name = 'log(' // name(:len1(name)) // ')'
	write( iu, * ) i, code, avg, rms, ' ', name(:len1(name))
      end do
* iterations
	write( iu, *) nits, ' MCMC iterations'
      do i = 1, nits
	write( iu, * ) i, bof_i( i )
	write( iu, * ) ( p_pi( j, i ), j = 1, npar )
      end do
* close
	close( iu )
	write(*,*) 'File closed : ', file(:len1(file))
* head/tail
	write(*,*)
	call spawn( 'ls -l ' // file(:len1(file)) )
	call spawn( 'head ' // file(:len1(file)) )
	call spawn( 'tail ' // file(:len1(file)) )

	return
	end
*------------------------------------------
	subroutine mcmcload
* load MCMC iterations
* 2021 Nov Keith Horne @ Crail - not yet implemented
* 2023 Dec KDH @ StA - first version
* 2023 Dec KDH @ StA - unlog and relog parameters as rquired
* 2023 Dec KDH @ StA - keep parameters in range
* 2023 Dec KDH @ StA - keep_pi counts parameter changes
* 2023 Dec KDH @ StA - acc_pi acceptance probability
* 2023 Dec KDH @ StA - dp_pi rms of previous parameter steps
	use all
	character*75 loadfile/ 'bowl_mcmc.dat' /
	character*256 record
	logical yes, inlog_p( maxpar )
	logical verbose / .true. /

	write(*,*)
	write(*,*) '--- mcmcload ------------'

* select load file
	write(*,*)
	write(*,*) 'Extant MCMC files:'
	call spawn( 'ls -l *mcmc*' )
    1	write(*,*)
	write(*,*) 'Enter none to abort'
	write(*,*) 'Enter ? for list of .dat files'
	file = loadfile
	call inq( 'Load file', loadfile, file )

      if( file .eq. '?' ) then
	call spawn( 'ls -l *.dat' )
	goto 1
      end if

* abort
	if( file .eq. 'none' .or. file .eq. 'NONE' ) return

	write(*,*)
* head and tail
	call spawn( 'ls -l ' // file(:len1(file)) )
	call spawn( 'head ' // file(:len1(file)) )
	call spawn( 'tail ' // file(:len1(file)) )

* open file
	iu = 50
	call openfile( iu, file, 'formatted', 'unknown', ierr )
	write(*,*) 'File open', ierr
      if( ierr .ne. 0 ) then
	write(*,*) '** OPENFILE FAILED.'
	return
      end if


* number of parameters
	read( iu, '(A)', err=998, end=999 ) record
	write(*,*) record( : len1( record ) )
	read( record, *, err=998, end=997 ) ipar
	write(*,*) ' MCMC parameters', ipar, ' should be', npar

      if( ipar .ne. npar ) then
	write(*,*) '** ERROR: THESE SHOULD MATCH. :(((.'
	goto 997
      end if

* parameters
      do i = 1, npar
* read record from file
	read( iu, '(A)', err=998, end=999 ) record
	write(*,*)
	write(*,*) record( : len1(record) )
* decode record  ( in, code, avg, rms, name )
      do k = 1, 5
	call word1( record, l1, l2 )
	if( l1 .le. 0 .or. l2 .lt. l1 ) goto 997
	if( k .eq. 1 ) read( record, *) in
	if( k .eq. 2 ) code = record(l1:l2)
	if( k .eq. 3 ) read( record, * ) avg
	if( k .eq. 4 ) read( record, * ) rms
	if( k .eq. 5 ) name = record(l1:l2)
	if( k .lt. 5 ) record = record( l2+1: )
      end do

* sequence number
	write(*,*) 'Parameter', in, ' should be',  i, ' of', npar
	if( in .ne. i ) goto 997

* code
	l1 = len1( code )
	l2 = len1( code_p( i ) )
	write(*,*) ' Code ', code(:l1), ' should be ', code_p(i)(:l2)
	if( code(:l1) .ne. code_p(i)(:l2) ) goto 997

* parameter name
	call word1( name, l1, l2 )
	call word1( name_p( i ), l3, l4 )
* log parameter
	logp = index( name, 'log' )
	yes = logp .gt. 0
	inlog_p( i ) = yes
      if( yes ) then
	l1 = logp + 4
	l2 = l2 - 1
      end if
	write(*,*) ' Log parameter ', yes, ' ', name(:len1(name))
	write(*,*) ' Name (', name(l1:l2), ') should be (', name_p(i)(l3:l4), ')'
	if( name(l1:l2) .ne. name_p(i)(l3:l4) ) goto 997

* stats
	write(*,*) '  avg:', avg, ' rms', rms

* next parameter
      end do
	write(*,*)
	write(*,*) 'Parameters in the file match those in the code. :)'
	write(*,*)
	call inq( 'C(ontinue) or A(bort) ? ', 'C', reply )
	call upcase( reply, reply )
	if( reply .ne. 'C' ) return

* iterations
	write(*,*)
	read( iu, '(A)', err=998, end=999 ) record
	write(*,*) record( : len1( record ) )

* number of new MCMC samples expected
	read( record, *, err=998, end=997 ) new
	write(*,*) ' MCMC samples in code', nits, ' of max', maxiter
	maxnew = maxiter - nits
	lost = max( 0, nits - maxnew )
	write(*,*) ' MCMC samples in file', new, ' of max', maxnew
	write(*,*) ' MCMC samples lost', lost

* new MCMC samples from the file
	new = min( new, maxnew )
	i1 = nits + 1
	i2 = nits + new
      do i = i1, i2
* read BoF
	read( iu, *, err=998, end=999 ) in, bof
* stow BoF
	bof_i( i ) = bof
	write(*,*)
	write(*,*) 'Iteration', in, ' =>', i, ' BoF', bof

* read parameters
	read( iu, * ) ( p_p( j ), j = 1, npar )
	write(*,*) 'Parameters:', ( p_p( j ), j = 1, npar )
      do j = 1, npar
	name = name_p( j )
	bot = bot_p( j )
	top = top_p( j )
* input parameter
	p = p_p( j )
	yes = inlog_p( j )
	if( yes) p = 10. ** p
* report
	write(*,*) p, ' log', yes, ' =>', log_p( j ), ' ', name(:len1(name))
* log parameter
	if( log_p( j ) ) p = alog10( max( p, real( 10. ** bot ) ) )
* parameter range (bot,top)
      if( p .le. bot .or. p .ge. top ) then
	write(*,*) '** Parameter', p, ' confined to range (', bot, top, ') **'
	p = max( bot, min( top, p ) )
      end if

* recover cumulative keep_pi(j,i) by noting parameter changes
	keep = 0
	im = i - 1
      if( im .gt. 0 ) then
	keep = keep_pi( j, im )
	change = p - p_pi( j, im )
      if( change .ne. 0. ) then
	keep = keep + 1
      end if
      end if
	keep_p( j ) = keep

* acceptance probability (%) over past nfix iterations
* KDH : assume no change of algorithm since previous mcmc run
* KDH : should make this a subroutine
      if( nfix .le. 0 ) then
	nfix = 50
	nfix = 30
      end if
	k1 = max( 1, i - nfix + 1 )
	keep = keep_p( j ) - keep_pi( j, k1 )
	n = max( 1, i - k1 )
	acc  = keep * 100. / n
	write(*,*) 'MCMC(', k1, i, ' ) keep', keep, ' of', n, '=', acc, ' %'

* approximately recover previously-used parameter steps.
* use rms of recently kept steps
* KDH : could use nMAD
	kback = 100
	k1 = max( 1, i - kback + 1 )
	k2 = i
	n = 0
	sum = 0.d0
      do k = k1, k2
	km = max( 1, k - 1 )
	dp = p_pi( j, k ) - p_pi( j, km )
      if( dp .ne. 0 ) then
	n = n + 1
	sum = sum + dp * dp
      end if
      end do
	dp = dp_p( j )
	if( n .gt. 3 ) dp = dsqrt( sum / n )
* stow
	p_pi( j, i ) = p
	dp_pi( j, i ) = dp
	acc_pi( j, i ) = acc
	keep_pi( j, i ) = keep_p( j )
* next parameter j
      end do

* next MCMC sample i
	nits = nits + 1
      end do
	write(*,*) 'New MCMC samples', new, ' Total', nits

      do i = 1, npar
	p_p( i ) = p_pi( i, nits )
	dp_p( i ) = dp_pi( i, nits )
      end do

* restore all parameters (both fixed and free)
* KDH : BELOW KLUDGE MAKES GETPAR THINK ALL PARAMETERS ARE FREE
* stow fit parameters
      do i = 1, nfit
	j_p( i ) = i_p( i )
      end do
	n = nfit
* claim all fit parameters are free
      do i = 1, npar
	i_p( i ) = i
      end do
	nfit = npar
* extract parameters from p_p
	call getpar
* restore previous fit parameters
	nfit = n
      do i = 1, nfit
	i_p( i ) = j_p( i )
      end do

* close
  100	close( iu )
	write(*,*) 'File closed : ', file(:len1(file))
	write(*,*)
* normal return
	return

  997	write(*,*) '** PARAMETER MISMATCH.'
	goto 100
  998	write(*,*) '** READ ERROR.'
	goto 100
  999	write(*,*) '** PREMATURE END OF FILE.'
	goto 100

	return
	end
* --------------------------------------
	subroutine mcmcfit
* Markov Chain Monte Carlo fit
* In:	from module all:
*	maxpar	i4 max number of parameters
*	maxiter	i4 max number of iterations
*	npar	i4 number of model parameters
*	nfit	i4 number of fit parameters
*	name_p(npar)	c* parameter names
*	code_p(npar)	c1 parameter code letter
*	units_p(npar)	c1 parameter units
*	symbol_p(npar)	c1 parameter pgplot symbol
*	p_p(npar) r4 parameters
*	dp_p(npar)	r4 parameter steps (gaussian sigma)
*	log_p(npar)	l4 .true. for log-uniform prior
*	bot_p(npar)	r4 lower limit
*	top_p(npar)	r4 upper limit
*	chi2	r8 chi-squared
*	sumlnv	r8 sum(ln(var))
*	bof	r8 badness of fit
*	nits	i4 iterations done
*	newits	i4 new iterations to be done
*	bof_i(nits)	r8 badness-of-fit per iteration
*	chi2_i(nits)	r8 chi-squared per iteration
*	cpu_i(nits)	r4 cpu used
*	p_pi(maxpar,nits)	r8 parameters
* Out:	more MCMC iterations updating most of above
*	ibest		i4 best mcmc sample
* KDH : THINNED SAMPLES OF FUNCTIONS FOR LATER PLOTTING
*	nsamp		i4 number of stowed sample functions
*	nthin		i4 thining interval (1,2,4,...)
*	iter_i(nsamp)	i4 iteration number
*	nr_i(nsamp)	i4 number of radii
*	nw_i(nsamp)	i4 number of wavelengths
*	nb_i(nsamp)	i4 number of bands
*	r_ri, h_ri, h0_ri	r4 geometry
*	t_ri, t0_ri, tc_ri	r4 temps
*	tx_ri, tv_ri, tq_ri	r4 temps
*	w_ri, f_ri, f0_ri	r4 seds
*	wb_ri, tau_ri, tau0_ri	r4 lags
* Use:
*	all	module supplies all variables needed
*	agn	module supplies agn variables
*	setpar	subroutine sets parameters into p_p
*	getpar	subroutine extracts parameters from p_p
*	bofcalc		r8 function computes badness-of-fit
*	mcmcreport	subroutine MCMC summary
*	cornerplot	subroutine corner plot
*	bofplot		subroutine iteration summary plot
*	mcmcsave	subroutine save MCMC iterations
*	mcmcload	subroutine load MCMC iterations
* 2021 Nov Keith Horne @ Crail
* 2022 Jul KDH @ Crail - nfit parameters i_p( ifit )
* 2022 Jul KDH @ StA - parameter symbols and units
* 2022 Aug KDH @ StA - plot or save after n iterations
* 2022 Aug KDH @ StA - bowlplot
* 2023 Dec KDH @ StA - mcmcload
* 2023 Dec KDH @ StA - verbose
* 2024 Mar KDH @ StA - stow thinned samples of functions for plotting
* 2024 Apr KDH @ StA - change nsp only if all samples requested
* 2024-May-09 KDH @ StA - adjust toggles step size adjustments
	use all
	use agn

	logical verbose/ .false. /
	logical yes

	write(*,*)
	write(*,*) '--- mcmcfit ---------------', nits
	write(*,*)


* initialise MCMC
  200 if( nits .le. 0 ) then
	nits = 0
	nsamp = 0
	nthin = 1
	call cpu( cpumcmc )
      do i = 1, npar
	keep_p( i ) = 0
	old_p( i ) = 0
      end do
	if( newits .le. 0 ) newits = 5
	call initpar
	call setpar
	write(*,*) 'MCMC initialised.'
      end if


* best fit MCMC sample
	ibest = 0
      if( nits .gt. 0 ) then
	ibest = 1
      do i = 1, nits
	if( bof_i( i ) .le. bof_i( ibest ) ) ibest = i
      end do
	write(*,*) '     MCMC(', 1,     ' ) BoF', bof_i( 1 )
	write(*,*) '     MCMC(', nits,  ' ) BoF', bof_i( nits )
	write(*,*) 'Best MCMC(', ibest, ' ) BoF', bof_i( ibest )
      end if

* initial BoF = Badness of Fit
	bof = bofcalc( p_p )

* report
	write(*,*) 'MCMC parameters', npar
	write(*,*)
     &	'      index',
     &	'    parameter',
     &	'          step'
      do i = 1, npar
	name = name_p( i )
	code = code_p( i )
	if( log_p(i) ) name = 'log(' // name(:len1(name)) // ')'
	write(*,*) i, code, p_p(i), dp_p(i), ' ', name(:len1(name))
      end do

  203	call mcmcreport

* MCMC menu -------------------------------------
  205	write(*,*) '------------------------------------------------'
	label = ' '
      do i = 1, npar
	label(i:i) = code_p(i)
      end do
	write(*,*) 'MCMC parameters:', npar, ' ', label(:len1(label))
	n = len1( fitlist )
      do i = 1, npar
	code = code_p(i)
	label(i:i) = ' '
	if( index( fitlist, code ) .gt. 0 ) label(i:i) = code
      end do
	write(*,*) ' BOWL fit parameters:', n, ' ', label(:len1(label))
	write(*,*) ' BOWL fit parameters:', n, ' ', fitlist(:len1(fitlist))

	if( n .le. 0 ) write(*,*) '** WARNING : NO FIT PARAMETERS.'

* parameter menu
      do i = 1, npar
	code = code_p( i )
	name = name_p( i )
	p = p_p( i )
	dp = 0.
	if( fit_p( i ) ) dp = dp_p( i )
	if( log_p( i ) ) name = 'log(' // name(:len1(name)) // ')'
	label = '---'
      do j = 1, n
	if( code .eq. fitlist(j:j) ) label = 'FIT'
      end do
	write(*,*) i, fit_p(i), ' ', label(:len1(label)), ' ', code,
     &	p, dp,	' ', name(:len1(name))
      end do


	write(*,*)
	write(*,*) ' + ... Add fit parameters (e.g. +DMA or +123)'
	write(*,*) ' - ... Cut fit parameters (e.g. -C* or -1011)'
	write(*,*) ' N ... New fit parameter list'
  208	write(*,*)
      if( nits .gt. 0 ) then
      if( iburn1 .gt. 0 ) then
	write(*,*) ' R ... Reset MCMC.  (Burnin stage 1 at', iburn1, '  )'
      else
	write(*,*) ' R ... Reset MCMC (Burnin in progress).'
      end if
      end if
	write(*,*) ' A ... Adjust step sizes every', nfix, ' iterations :', adjust
	write(*,*) ' I ... Iterate  done', nits, ' new', newits
      if( nits .gt. 0 ) then
	write(*,*) '      CPU:', cpumcmc, ' CPU/it:', cpumcmc / max(1,nits)
	write(*,*) '     MCMC(', nits,  ' ) BoF', bof_i( nits )
	write(*,*) '     MCMC(', ibest, ' ) BoF', bof_i( ibest ), ' best'
      end if
	write(*,*) '  (IP,IB, IC, IS iterates with plots or saves)'
	write(*,*) ' P ... Plot fit     every', itfit,    ' iterations.'
	write(*,*) ' B ... BoF plot     every', itbof,    ' iterations.'
	write(*,*) ' C ... Corner plot  every', itcorner, ' iterations.'
	write(*,*) ' S ... Save results every', itsave,   ' iterations.'
	write(*,*) '  (Pn, Bn, Cn, Sn to plot or save after n iterations.)'
	write(*,*) ' L ... Load MCMC iterations'
	write(*,*) ' Q ... Quit MCMC fitting'

* MCMC command
  210	if( reply .ne. '?' ) reply = ' '
	call inq( 'MCMC>', reply, reply )
	call upcase( reply, reply )

      if( reply(1:1) .eq. '?' ) then
	reply = ' '
	goto 205
      end if
	if( reply(1:1) .eq. 'Q' ) return

      if( nits .gt. 0 ) then
	if( reply(1:1) .eq. 'R' ) goto 215
      end if
	if( reply(1:1) .eq. 'P' ) goto 218
	if( reply(1:1) .eq. 'B' ) goto 220
	if( reply(1:1) .eq. 'C' ) goto 225
	if( reply(1:1) .eq. 'S' ) goto 230
	if( reply(1:1) .eq. 'L' ) goto 235
	if( reply(1:1) .eq. 'I' ) goto 240
	if( reply(1:1) .eq. 'N' ) goto 250
	if( reply(1:1) .eq. '+' ) goto 251
	if( reply(1:1) .eq. '-' ) goto 251
	if( reply(1:1) .eq. 'A' ) goto 214

* is the reply a valid number of iterations?
	call word1( reply, i1, i2 )
      if( i1 .gt. 0 .and. i2 .ge. i1 ) then
	read( reply(i1:i2), *, err=213, end=213 ) i
      if( i .ge. 1 ) then
	newits = i
	goto 240
      end if
      end if

* polite
  213	write(*,*) '** Command not recognised.  Please try again.'
	reply = '?'
	goto 210

* adjust step size interval
  214	write(*,*) 'Adjust step sizes every NFIX iterations.'
	write(*,*) 'Enter -Nfix to suspend step size adjustments.'
	call inqi( 'Adjust interval NFIX', nfix, i )
	adjust = i .gt. 0
	if( i .ne. 0 ) nfix = iabs( i )
	goto 208

* reset
  215	write(*,*) 'MCMC Iterations', nits
	call inq( 'Do you REALLY want to reset?', 'N', reply )
	call upcase( reply, reply )
	if( reply(1:1) .ne. 'Y' ) goto 200
	write(*,*) 'MCMC Reset.'
* KDH : NEED OPTION TO CUT ONLY FIRST PART OF THE CHAIN
	nits = 0
	nsamp = 0
	nthin = 1
	goto 200


* bowl fit plot
  218 if( len1( reply ) .eq. 1 ) then
	call bowlplot
* iterations between plots
      else
	read( reply(2:), *, end=200, err=200 ) i
	itfit = i
	write(*,*) 'Plot fit every', i, ' iterations.'
      end if
	goto 208


* bof vs iteration plot
  220 if( len1( reply ) .eq. 1 ) then
	call bofplot
* iterations between plots
      else
	read( reply(2:), *, end=200, err=200 ) i
	itbof = i
	write(*,*) 'Bof plot every', i, ' iterations.'
      end if
	goto 208

* corner plot
  225 if( len1( reply ) .eq. 1 ) then
	call cornerplot
* iterations between plots
      else
	read( reply(2:), *, end=200, err=200 ) i
	itcorner = i
	write(*,*) 'Corner plot every', i, ' iterations.'
      end if
	goto 208

* save MCMC chain
  230 if( len1( reply ) .eq. 1 ) then
	call mcmcsave
* iterations between saves
      else
	read( reply(2:), *, end=200, err=200 ) i
	itsave = i
	write(*,*) 'Save results every', i, ' iterations.'
      end if
C	goto 200
	goto 208

* load iterations
  235	n = nits
	call mcmcload
	if( n .ne. nits ) call mcmcreport
c	goto 210
	goto 200

* new iterations
  240	iter = nits
	if( newits .eq. 0 ) newits = 5
	maxnewits = maxiter - nits

	write(*,*) 'BOWL fit parameters', nfit, ' of', npar
      if( nfit .le. 0 ) then
	write(*,*) '** NO FIT PARAMETERS SELECTED.'
	goto 210
      end if

	write(*,*) 'Iterations done', nits, ' of max', maxiter
	if( index( reply, 'P' ) .gt. 0 .and. itfit .gt. 0 )
     &		write(*,*) 'Plot fit every', itfit, ' iterations.'
	if( index( reply, 'B' ) .gt. 0 .and. itbof .gt. 0 )
     &		write(*,*) 'BoF plot every', itbof, ' iterations.'
	if( index( reply, 'C' ) .gt. 0 .and. itcorner .gt. 0 )
     &		write(*,*) 'Corner plot every', itcorner, ' iterations.'
	if( index( reply, 'S' ) .gt. 0 .and. itsave .gt. 0 )
     &		write(*,*) 'Save every', itsave, ' iterations.'

	write(*,*) 'Max new iterations:', maxnewits
	call inqi( 'New iterations', newits, i )
c	if( newits .le. 0 ) goto 203
	if( i .le. 0 ) goto 203
	newits = min( i, maxnewits )

	last = nits + newits
	write(*,*) 'New iterations', newits, ' last', last

	call cpu( t )
      do it = 1, newits
	iter = iter + 1
	nits = iter
	ntogo = newits - it
	write(*,*)
	write(*,*) '==== MCMC ', nits, ' of', last, ' TO GO', ntogo, ' ===='

* gibbs sampling
      do ifit = 1, nfit
	i = i_p( ifit )
	name = name_p( i )
	code = code_p( i )
	write(*,*)
	itogo = nfit - ifit + 1
	write(*,*) '==== MCMC iter', nits, ' of', last, ' TO GO', ntogo, ' ===='
	write(*,*) '==== Gibbs par', ifit, ' of', nfit, ' TO GO', itogo, ' ===='
* stow
	bold = bof
	pold = p_p( i )
* step parameter i
	dp = dp_p( i )
	pnew = rang( pold, dp, iseed )

* report trial sample
      if( log_p( i ) ) then
	write(*,*) '     Model ', fitlist(:len1(fitlist))
     &	, ' Code ',  code, ' log(', name(:len1(name)), ')', pold, '->', pnew
      else
	write(*,*) '     Model ', fitlist(:len1(fitlist))
     &	, ' Code ',  code, ' ', name(:len1(name)), pold, '->', pnew
      end if

* reject if out of range
	top = top_p( i )
	bot = bot_p( i )
      if( pnew .lt. bot .or. pnew .gt. top ) then
	write(*,*) '** OUT OF RANGE (', bot, top, ' )'
	p_p( i ) = pold
	keep = 0
* ok if inside range or on boundary
      else
	p_p( i ) = pnew
	keep = 1

* KDH : MOST PARAMETER STEPS DON'T NEED FULL SECONDARY UPDDATE
*       COULD STREAMLINE FOR SPEED
* update secondary parameters 	NOW DONE INSIDE BOFCALC
C	call secondary
* DEBUG REPORT
      if( verbose ) then
	write(*,*) 'Xlamp', xlamp, ' Tv1', tv1, ' Tx1', tx1
	write(*,*) 'H1/Rg', h1_rg, ' bet = dlnH/dlnR', bet
	write(*,*) 'Rout', rout
	write(*,*) 'sigsed', dexsed, ' dex =', 100.*sigsed, ' %'
	write(*,*) 'syslag', syslag, ' d'
      end if
* new BoF
	needsed = .true.
	needlag = .true.
	bof = bofcalc( p_p )
* relative probability
	dbof = bof - bold
	p = exp( - dbof / 2. )
	write(*,*) 'dBoF =', real( bof ), ' -', real( bold )
     &		, ' =', real( dbof ), ' Prel', p
* good move
c	keep = 1
      if( bof .lt. bold ) then
	write(*,*) '*** GOOD MOVE', p, ' KEEP !'
* keep bad move
      else
	u = ran3( iseed )
      if( u .lt. p ) then
	write(*,*) '*** BAD MOVE', p, ' KEEP ANYWAY'
* retract bad move
      else
	write(*,*) '*** BAD MOVE', P, ' RETRACT'
	keep = 0
	bof = bold
	p_p( i ) = pold
      end if
      end if
* count keeps for this parameter
	keep_p( i ) = keep_p( i ) + keep
      end if

* acceptance probability (%) over past nfix iterations
      if( nfix .le. 0 ) then
	nfix = 50
	nfix = 30
      end if
	k = max( 1, iter - nfix + 1 )
	keep = keep_p( i ) - keep_pi( i, k )
	k = max( 1, iter - k )
	acc  = keep * 100. / k
	write(*,*) 'Kept over last', k, ' iterations:', nint(acc), ' %'
* stow
	  dp_pi( i, iter ) = dp
	keep_pi( i, iter ) = keep_p( i )
	 acc_pi( i, iter ) = acc
* next fit parameter i = i_p( ifit )
      end do

* adjust step size every nfix iterations
      if( adjust .and. mod( iter, nfix ) .eq. 0 ) then
	boost = 1.5
	write(*,*)
	write(*,*) '** Revise MCMC parameter steps', iter, ' **'
* parameters
      do ifit = 1, nfit
	i = i_p( ifit )
	code = code_p( i )
	name = name_p( i )
	if( log_p( i ) ) name = 'log(' // name(:len1(name)) // ')'
	write(*,*) ' MCMC parameter', i, ' ', code, ' ', name(:len1(name))
* percent kept
	acc = acc_pi( i, iter )
* boost if > 40% accepted
	old = dp_p( i )
      if( acc .gt. 80. ) then
	dp_p( i ) = old * boost * boost * boost
	write(*,*) 'Keep', nint(acc), ' > 80%  step', old, ' =>', dp_p(i)
      else if( acc .gt. 60. ) then
	dp_p( i ) = old * boost * boost
	write(*,*) 'Keep', nint(acc), ' > 60%  step', old, ' =>', dp_p(i)
      else if( acc .gt. 40. ) then
	dp_p( i ) = old * boost
	write(*,*) 'Keep', nint(acc), ' > 40%  step', old, ' =>', dp_p(i)
* stifle if < 20% accepted
      else if( acc .lt. 20. ) then
	dp_p( i ) = old / boost
	write(*,*) 'Keep', nint(acc), ' < 20%  step', old, ' =>', dp_p(i)
* stifle if < 10% accepted
      else if( acc .lt. 10. ) then
	dp_p( i ) = old / boost / boost
	write(*,*) 'Keep', nint(acc), ' < 10%  step', old, ' =>', dp_p(i)
* otherise ok
      else
	write(*,*) 'Keep', nint(acc), ' %  step', old, ' OK'
      end if
* cap log steps
	big = 0.3
	big = 0.5
      if( log_p( i ) .and. dp_p( i ) .gt. big ) then
	write(*,*) '** cap the step at', big, ' dex.'
	dp_p( i ) = big
      end if
* next parameter i = i_p( ifit )
      end do
	write(*,*) 'End of step size adjustments.'
	write(*,*)
      end if

* stow all parameter values
      do i = 1, npar
	p_pi( i, iter ) = p_p( i )
      end do
	bof_i( iter ) = bof

* stow functions (for later plotting)
* KDH : MAKE THIS A SUBROUTINE
* KDH : THIN ON STOW TO AVOID WRAPAROUND
      if( nsamp .ge. maxi ) then
	i = 0
      do ii = 2, nsamp, 2
	i = i + 1
	iter_i( i ) = iter_i( ii )
	  nr_i( i ) =   nr_i( ii )
      do ir = 1, nr_i( i )
	 r_ri( ir, i ) =  r_ri( ir, ii )
	 h_ri( ir, i ) =  h_ri( ir, ii )
	h0_ri( ir, i ) = h0_ri( ir, ii )
	 t_ri( ir, i ) =  t_ri( ir, ii )
	t0_ri( ir, i ) = t0_ri( ir, ii )
	tc_ri( ir, i ) = tc_ri( ir, ii )
	tx_ri( ir, i ) = tx_ri( ir, ii )
	tv_ri( ir, i ) = tv_ri( ir, ii )
	tq_ri( ir, i ) = tq_ri( ir, ii )
      end do
	hlamp_i( i ) = hlamp_i( ii )
	nb_i( i ) = nb_i( ii )
      do ib = 1, nb_i( i )
	   wobs_bi( ib, i ) =    wobs_bi( ib, ii )
	 tau_bi( ib, i ) =  tau_bi( ib, ii )
	tau0_bi( ib, i ) = tau0_bi( ib, ii )
	tau1_bi( ib, i ) = tau1_bi( ib, ii )
	tau2_bi( ib, i ) = tau2_bi( ib, ii )
      end do
	nw_i( i ) = nw_i( ii )
      do iw = 1, nw_i( i )
	 w_wi( iw, i ) =  w_wi( iw, ii )
	 f_wi( iw, i ) =  f_wi( iw, ii )
	f0_wi( iw, i ) = f0_wi( iw, ii )
	f1_wi( iw, i ) = f1_wi( iw, ii )
	f2_wi( iw, i ) = f2_wi( iw, ii )
      end do
* next ii
      end do
	write(*,*) 'Nsamp', nsamp, ' =>', i
	write(*,*) 'Nthin', nthin, ' =>', nthin * 2
	nsamp = i
	nthin = nthin * 2
      end if

* where to stow
	yes = nsp .eq. nsamp
	i = ( iter + nthin - 1 ) / nthin
	nsamp = min( maxi, i )
	if( yes ) nsp = nsamp
* iteration number
	iter_i( i ) = iter
* report
	write(*,*) 'Stow iteration', iter, ' at sample', i, ' of', nsamp
	write(*,*) 'Nr', nr, ' Nw', nw, ' Nb', nb
* geometry
	nr_i( i ) = nr
      do ir = 1, nr
	 r_ri( ir, i ) =  r_r( ir )
	 h_ri( ir, i ) =  h_r( ir )
	h0_ri( ir, i ) = h0_r( ir )
	hlamp_i( i ) = hlamp
* temps T(r)
	 t_ri( ir, i ) =  t_r( ir )
	t0_ri( ir, i ) = t0_r( ir )
	tc_ri( ir, i ) = tc_r( ir )
	tx_ri( ir, i ) = tx_r( ir )
	tv_ri( ir, i ) = tv_r( ir )
	tq_ri( ir, i ) = tq_r( ir )
      end do
* lags
	nb_i( i ) = nb
      do ib = 1, nb
	   wobs_bi( ib, i ) =    wobs_b( ib )
	 tau_bi( ib, i ) =  tau_b( ib )
	tau0_bi( ib, i ) = tau0_b( ib )
	tau1_bi( ib, i ) = tau1_b( ib )
	tau2_bi( ib, i ) = tau2_b( ib )
      end do
* SEDs Fnu(lam)
	nw_i( i ) = nw
      do iw = 1, nw
	 w_wi( iw, i ) =  w_w( iw )
	 f_wi( iw, i ) =  f_w( iw )
	f0_wi( iw, i ) = f0_w( iw )
	f1_wi( iw, i ) = f1_w( iw )
	f2_wi( iw, i ) = f2_w( iw )
      end do

* report
	call mcmcreport

* cpu
	call cpu( t )
	cpumcmc = cpumcmc + t
	write(*,*) 'CPU:', cpumcmc, ' CPU/iter:', t, cpumcmc / max(1,iter)

* make plots or save
	ip = index( reply, 'P' )
	ib = index( reply, 'B' )
	ic = index( reply, 'C' )
	is = index( reply, 'S' )
* retain the requested order
      do i = 1, max( ib, ic, is, ip )
* fit plot
      if( i .eq. ip .and. itfit .gt. 0 ) then
	if( mod( it, itfit ) .eq. 0 ) call bowlplot
      end if

* bof plot
      if( i .eq. ib .and. itbof .gt. 0 ) then
	if( mod( it, itbof ) .eq. 0 ) call bofplot
      end if

* corner plot
      if( i .eq. ic .and. itcorner .gt. 0 ) then
	if( mod( it, itcorner ) .eq. 0 ) call cornerplot
      end if
* save
      if( i .eq. is .and. itsave .gt. 0 ) then
	if( mod( it, itsave ) .eq. 0 ) call mcmcsave
      end if
* next plot or save i
      end do


* next iteration it
      end do

	goto 200


* select parameters to fit
  250	write(*,*)
	write(*,*) '+XYZ appends XYZ and -XYZ removes XYZ'

	reply = fitlist
	call inq( 'MCMC fit parameter(s)', reply, reply )
	call upcase( reply, reply )

  251	write(*,*) '  Reply ',   reply( : len1(reply) )
	write(*,*)' Fitlist ', fitlist( : len1(fitlist) )

* stow for later reference
	label = fitlist

* Let +/ - append/strip codes from the list

	ip = index( reply, '+' )
	im = index( reply, '-' )
      if( ip .gt. 0 .and. im .gt. 0 ) then
	write(*,*) '+ index', ip, ' - index', im
	write(*,*) '** CANNOT HAVE BOTH + and -'
      end if

* replace numbers by codes  (crude but should work up to 11)
	i2 = len1( reply )
	i1 = i2
	if( ip .gt. 0 ) i1 = min( i1, ip )
	if( im .gt. 0 ) i1 = min( i1, im )
      if( i1 .gt. 0 .and. i1 .lt. i2 ) then
      do i = i1, i2
	if( reply(i:i+1) .eq. '10' ) reply(i:i+1) = code_p( 10 )
	if( reply(i:i+1) .eq. '11' ) reply(i:i+1) = code_p( 11 )
* KDH : 12 is ambiguous with 1 and 2, and similarly for 13, 14, ...
	if( reply(i:i) .eq. '1' ) reply(i:i) = code_p( 1 )
	if( reply(i:i) .eq. '2' ) reply(i:i) = code_p( 2 )
	if( reply(i:i) .eq. '3' ) reply(i:i) = code_p( 3 )
	if( reply(i:i) .eq. '4' ) reply(i:i) = code_p( 4 )
	if( reply(i:i) .eq. '5' ) reply(i:i) = code_p( 5 )
	if( reply(i:i) .eq. '6' ) reply(i:i) = code_p( 6 )
	if( reply(i:i) .eq. '7' ) reply(i:i) = code_p( 7 )
	if( reply(i:i) .eq. '8' ) reply(i:i) = code_p( 8 )
	if( reply(i:i) .eq. '9' ) reply(i:i) = code_p( 9 )
	write(*,*) i, ' ', reply(:len1(reply))
      end do
      end if

* add to fitlist
      if( reply(1:1) .eq. '+' ) then
	fitlist = fitlist(:len1(fitlist)) // reply(2:)
* reply is new fitlist
      else if( reply(1:1) .ne. '-' ) then
	fitlist = reply
* remove from fitlist
      else
      do while ( len1(reply) .gt. 0 )
	write(*,*) 'Remove ', reply(:len1(reply))
	code = reply(1:1)
	i = index( fitlist, code )
      if( i .eq. 1 ) then
	fitlist = fitlist(2:)
      else if( i .gt. 1 ) then
	fitlist = fitlist(:i-1) // fitlist(i+1:)
      end if
	reply = reply(2:)
      end do
      end if

* verify fitlist
	write(*,*) 'Validate ', fitlist(:len1(fitlist))
	reply = ' '
      do i = 1, len1( fitlist )
	code = fitlist( i : i )
      do j = 1, npar
      if( code .eq. code_p( j ) ) then
	call append( reply, code )
	call nospace( reply )
	write(*,*) i, j, code, ' ', reply(:len1(reply) )
      end if
      end do
      end do
	fitlist = reply(:len1(reply))

* reset if changed
      if( fitlist .ne. label ) then
c	call initpar
	call setpar
      end if

	goto 203


* end subroutine mcmcfit
	end
*===========================================
	subroutine pgmdot( string )
* pgplot character string for Mdot (a dot over an M)
* 2022 Aug Keith Horne @ St.Andrews
* 2023 Sep KDH @ Crail - try to avoid starting with a space
	character*(*) string
	string =' \u\(828)\b\b\d\bM'
* possible space characters 697, 698, 699
c	string ='M\u\(828)\d'
c	string ='M\u\b\b\(828)\b\b\dM'
	return
	end

*===========================================
	subroutine tidy( string )
* remove trailing 0 and decimal point
* 2025 Jan Keith Horne @ StAndrews
	character*(*) string
	i = index( string, '.' ) - 1
      if( i .gt. 0 ) then
	l = len1( string )
      do k = i, l
	call notail( string(i:), '0' )
      end do
	call notail( string(i:), '.' )
      end if

	return
	end
