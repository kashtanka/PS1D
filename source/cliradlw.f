
c**********************   June 2003  *****************************
 
      subroutine irrad (m,np,pl,ta,wa,oa,tb,co2,
     *                  high,trace,n2o,ch4,cfc11,cfc12,cfc22,
     *                  vege,ns,fs,tg,eg,tv,ev,rv,
     *                  overcast,cldwater,cwc,taucl,fcld,ict,icb,
     *                  aerosol,na,taual,ssaal,asyal,
     *                  flx,flc,dfdts,sfcem)
 
c*********************************************************************
 
c   THE EQUATION NUMBERS noted in this code follows the latest  
c    version (May 2003) of the NASA Tech. Memo. (2001), which can 
c    be accessed at ftp://climate.gsfc.nasa.gov/pub/chou/clirad_lw/
 
c*********************************************************************
 
c  CHANGE IN MAY 2003
 
c    The effective size of ice particles is replaced changed to the  
c    effective radius, and the definition of the effective radius 
c    follows that given in Chou, Lee and Yang (JGR, 2002).
c    The reff is no longer an input parameter.
 
c  CHANGE IN DECEMBER 2002
 
c    Do-loop 1500 is created inside the do-loop 1000 to compute 
c    the upward and downward emission of a layer
 
c   CHANGE IN JULY 2002
 
c    The effective Planck functions of a layer are separately
c    computed for the upward and downward emission (bu and bd).
c    For a optically thick cloud layer, the upward emission will be
c    at the cloud top temperature, and the downward emission will
c    at the cloud base temperature.
 
c   RECENT CHANGES:
c
c    Subroutines for planck functions
c    Subroutines for cloud overlapping
c    Eliminate "rflx" and "rflc". Fold the flux calculations in
c      Band 10 to that of the other bands.
c    Return the calculations when ibn=10 and trace=.false.
c    The number of aerosol types is allowed to be more than one.
c    Include sub-grid surface variability and vegetation canopy.
c    Include the CKD continuum absorption coefficient as an option.
c    
c********************************************************************
 
c Ice and liquid cloud particles are allowed to co-exist in each of the
c  np layers. 
c
c The maximum-random assumption is applied for cloud overlapping. 
c  Clouds are grouped into high, middle, and low clouds separated 
c  by the level indices ict and icb.  Within each of the three groups,
c  clouds are assumed maximally overlapped.  Clouds among the three 
c  groups are assumed randomly overlapped. The indices ict and icb 
c  correspond approximately to the 400 mb and 700 mb levels.
c
c Various types of aerosols are allowed to be in any of the np layers. 
c  Aerosol optical properties can be specified as functions of height  
c  and spectral band.
c
c The surface can be divided into a number of sub-regions either with or 
c  without vegetation cover. Reflectivity and emissivity can be 
c  specified for each sub-region.
c
c There are options for computing fluxes:
c
c   If high = .true., transmission functions in the co2, o3, and the
c   three water vapor bands with strong absorption are computed using
c   table look-up.  cooling rates are computed accurately from the
c   surface up to 0.01 mb.
c   If high = .false., transmission functions are computed using the
c   k-distribution method with linear pressure scaling for all spectral
c   bands except Band 5.  cooling rates are not accurately calculated 
c   for pressures less than 10 mb. the computation is faster with
c   high=.false. than with high=.true.
c
c   If trace = .true., absorption due to n2o, ch4, cfcs, and the 
c   two minor co2 bands in the window region is included.
c   Otherwise, absorption in those minor bands is neglected.
c
c   If vege=.true., a vegetation layer is added, and the emission and 
c   reflectivity are computed for the ground+vegetation surface.
c   Otherwise, only ground and ocean surfaces are considered.
c
c   If overcast=.true., the layer cloud cover is either 0 or 1.
c   If overcast=.false., the cloud cover can be anywhere between 0 and 1.
c   Computation is faster for the .true. option than the .false. option.
c
c   If cldwater=.true., taucl is computed from cwc and reff as a
c   function of height and spectral band. 
c   If cldwater=.false., taucl must be given as input to the radiation
c   routine. For this case, taucl is independent of spectral band.
c
c   If aerosol = .true., aerosols are included in calculating transmission
c   functions. Otherwise, aerosols are not included.
c   
c
c The IR spectrum is divided into nine bands:
c   
c   band     wavenumber (/cm)   absorber
c
c    1           0 - 340           h2o
c    2         340 - 540           h2o
c    3         540 - 800       h2o,cont,co2
c    4         800 - 980       h2o,cont
c                              co2,f11,f12,f22
c    5         980 - 1100      h2o,cont,o3
c                              co2,f11
c    6        1100 - 1215      h2o,cont
c                              n2o,ch4,f12,f22
c    7        1215 - 1380      h2o,cont
c                              n2o,ch4
c    8        1380 - 1900          h2o
c    9        1900 - 3000          h2o
c
c In addition, a narrow band in the 17 micrometer region (Band 10) is added
c    to compute flux reduction due to n2o
c
c    10        540 - 620       h2o,cont,co2,n2o
c
c Band 3 (540-800/cm) is further divided into 3 sub-bands :
c
c   subband   wavenumber (/cm)
c
c    3a        540 - 620
c    3b        620 - 720
c    3c        720 - 800
c
c---- Input parameters                               units    size
c
c   number of soundings (m)                            --      1
c   number of atmospheric layers (np)                  --      1
c   level pressure (pl)                               mb      m*(np+1)
c   layer temperature (ta)                            k       m*np
c   layer specific humidity (wa)                      g/g     m*np
c   layer ozone mixing ratio by mass (oa)             g/g     m*np
c   surface air temperature (tb)                      k        m
c   co2 mixing ratio by volume (co2)                  pppv     1
c   option (high) (see explanation above)              --      1
c   option (trace) (see explanation above)             --      1
c   n2o mixing ratio by volume (n2o)                  pppv     1
c   ch4 mixing ratio by volume (ch4)                  pppv     1
c   cfc11 mixing ratio by volume (cfc11)              pppv     1
c   cfc12 mixing ratio by volume (cfc12)              pppv     1
c   cfc22 mixing ratio by volume (cfc22)              pppv     1
c   option for including vegetation cover (vege)       --      1
c   number of sub-grid surface types (ns=2)            --      m
c   fractional cover of sub-grid regions (fs)       fraction  m*ns
c   land or ocean surface temperature (tg)            k       m*ns
c   land or ocean surface emissivity (eg)           fraction  m*ns*9
c   vegetation temperature (tv)                       k       m*ns
c   vegetation emissivity (ev)                      fraction  m*ns*9
c   vegetation reflectivity (rv)                    fraction  m*ns*9
c   option for cloud fractional cover                  --      1
c      (overcast)   (see explanation above)
c   option for cloud optical thickness                 --      1
c      (cldwater)   (see explanation above)
c   cloud water mixing ratio (cwc)                   gm/gm   m*np*3
c       index 1 for ice particles
c       index 2 for liquid drops
c       index 3 for rain drops
c   cloud optical thickness (taucl)                    --    m*np*3
c       index 1 for ice particles
c       index 2 for liquid drops
c       index 3 for rain drops
c   cloud amount (fcld)                             fraction  m*np
c   level index separating high and middle             --      1
c       clouds (ict)
c   level index separating middle and low              --      1
c       clouds (icb)
c   option for including aerosols (aerosol)            --      1
c   number of aerosol types (na)                       --      1
c   aerosol optical thickness (taual)                  --   m*np*10*na
c   aerosol single-scattering albedo (ssaal)           --   m*np*10*na
c   aerosol asymmetry factor (asyal)                   --   m*np*10*na
c
c---- output parameters
c
c   net downward flux, all-sky   (flx)             w/m**2  m*(np+1)
c   net downward flux, clear-sky (flc)             w/m**2  m*(np+1)
c   sensitivity of net downward flux  
c       to surface temperature (dfdts)            w/m**2/k m*(np+1)
c   emission by the surface (sfcem)                w/m**2     m
c
c Data used in table look-up for transmittance calculations:
c
c   c1 , c2, c3: for co2 (band 3)
c   o1 , o2, o3: for  o3 (band 5)
c   h11,h12,h13: for h2o (band 1)
c   h21,h22,h23: for h2o (band 2)
c   h81,h82,h83: for h2o (band 8)
c 
c Notes: 
c
c   (1) Scattering is parameterized for clouds and aerosols.
c   (2) Diffuse cloud and aerosol transmissions are computed
c       from exp(-1.66*tau).
c   (3) If there are no clouds, flx=flc.
c   (4) plevel(1) is the pressure at the top of the model atmosphere,
c        and plevel(np+1) is the surface pressure.
c   (5) Downward flux is positive and upward flux is negative.
c   (6) sfcem and dfdts are negative because upward flux is defined as negative.
c   (7) For questions and coding errors, please contact Ming-Dah Chou,
c       Code 913, NASA/Goddard Space Flight Center, Greenbelt, MD 20771.
c       Phone: 301-614-6192, Fax: 301-614-6307,
c       e-mail: chou@climate.gsfc.nasa.gov
c
c***************************************************************************
 
      implicit none
 
c---- input parameters ------
 
      integer m,np,na,ns,ict,icb
      real pl(m,np+1),ta(m,np),wa(m,np),oa(m,np),tb(m)
      real co2,n2o,ch4,cfc11,cfc12,cfc22
      real fs(m,ns),tg(m,ns),eg(m,ns,10)
      real tv(m,ns),ev(m,ns,10),rv(m,ns,10)
      real cwc(m,np,3),taucl(m,np,3),fcld(m,np)
      real taual(m,np,10,na),ssaal(m,np,10,na),asyal(m,np,10,na)
      logical high,vege,trace,overcast,cldwater,aerosol
 
c---- output parameters ------
 
      real flx(m,np+1),flc(m,np+1),dfdts(m,np+1),sfcem(m)
 
c---- static data -----
 
      real xkw(9),xke(9),aw(9),bw(9),pm(9),fkw(6,9),gkw(6,3)
      real aib(3,10),awb(4,10),aiw(4,10),aww(4,10),aig(4,10),awg(4,10)
      integer mw(9)
 
c-----parameters defining the size of the pre-computed tables for
c     transmittance using table look-up.
 
c     "nx" is the number of intervals in pressure
c     "no" is the number of intervals in o3 amount
c     "nc" is the number of intervals in co2 amount
c     "nh" is the number of intervals in h2o amount
 
      integer nx,no,nc,nh
      parameter (nx=26,no=21,nc=30,nh=31)
 
      real c1 (nx,nc),c2 (nx,nc),c3 (nx,nc)
      real o1 (nx,no),o2 (nx,no),o3 (nx,no)
      real h11(nx,nh),h12(nx,nh),h13(nx,nh)
      real h21(nx,nh),h22(nx,nh),h23(nx,nh)
      real h81(nx,nh),h82(nx,nh),h83(nx,nh)
 
c---- temporary arrays -----
 
      real pa(m,np),dt(m,np),tx(m),xlayer(m),reff(m,np,3)
      real x1(m),x2(m),x3(m)
      real dh2o(m,np),dcont(m,np),dco2(m,np),do3(m,np)
      real dn2o(m,np),dch4(m,np)
      real df11(m,np),df12(m,np),df22(m,np)
      real th2o(m,6),tcon(m,3),tco2(m,6,2)
      real tn2o(m,4),tch4(m,4),tcom(m,6)
      real tf11(m),tf12(m),tf22(m)
      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
      real n2oexp(m,np,4),ch4exp(m,np,4),comexp(m,np,6)
      real f11exp(m,np),f12exp(m,np),f22exp(m,np)
      real blayer(m,0:np+1),blevel(m,np+1)
      real bd(m,0:np+1),bu(m,0:np+1)
      real bs(m),dbs(m),rflxs(m)
      real dp(m,np),cwp(m,np,3)
      real trant(m,np+1),tranal(m),transfc(m,np+1),trantcr(m,np+1)
      real flxu(m,np+1),flxd(m,np+1),flcu(m,np+1),flcd(m,np+1)
      real taua(m,np),ssaa(m,np),asya(m,np),taerlyr(m,np)
 
      integer it(m),im(m),ib(m),itx(m,np),imx(m,np),ibx(m,np)
      real cldhi(m),cldmd(m),cldlw(m),tcldlyr(m,np),fclr(m,np+1)
 
      integer i,j,k,ip,iw,ibn,ik,iq,isb,k1,k2,ne
      real x,xx,yy,p1,dwe,dpe,a1,b1,fk1,a2,b2,fk2
      real w1,w2,w3,g1,g2,g3,ww,gg,ff,tauc
 
      logical oznbnd,co2bnd,h2otbl,conbnd,n2obnd
      logical ch4bnd,combnd,f11bnd,f12bnd,f22bnd,b10bnd
 
c-----coefficients for computing effective particle size following
c     McFarquhar (QJRMS, 2000)
 
      real ai(5),bi(5),ci(5),di(5)
 
      data ai/    2.076,    2.054,    2.035,    2.019,    2.003/
      data bi/    0.148,    0.130,    0.119,    0.111,    0.102/
      data ci/  -0.0453,  -0.0491,  -0.0507,  -0.0517,  -0.0532/
      data di/ -0.00686, -0.00711, -0.00716, -0.00717, -0.00725/
 
c-----xkw is the absorption coefficient for the first k-distribution
c     interval due to water vapor line absorption (Table 4)
c     Units are cm**2/g    
 
      data xkw / 29.55  , 4.167e-1, 1.328e-2, 5.250e-4,
     *           5.25e-4, 9.369e-3, 4.719e-2, 1.320e-0, 5.250e-4/
 
c-----xke is the absorption coefficient for the first k-distribution
c     function due to water vapor continuum absorption (Table 9).
c     Units are cm**2/g
 
c-----Roberts et al's continuum k data
 
c     data xke /  0.00,   339.00,  27.40,   15.8,
c    *            9.40,   7.75,    7.70,    0.0,   0.0/
 
c-----CKD (Version 2.3) continuum k data
 
      data xke /  0.0,    271.,    25.00,   16.8,
     *            8.31,   6.52,    12.7,    0.0,  0.0/
 
c-----mw is the ratio between neighboring absorption coefficients
c     for water vapor line absorption (Table 4).
 
      data mw /6,6,8,6,6,8,9,6,16/
 
c-----aw and bw (Table 3) are the coefficients for temperature scaling
c     for water vapor in Eq. (4.2).
 
      data aw/ 0.0021, 0.0140, 0.0167, 0.0302,
     *         0.0307, 0.0195, 0.0152, 0.0008, 0.0096/
      data bw/ -1.01e-5, 5.57e-5, 8.54e-5, 2.96e-4,
     *          2.86e-4, 1.108e-4, 7.608e-5, -3.52e-6, 1.64e-5/
 
c-----pm is the pressure-scaling parameter for water vapor absorption
c     Eq. (4.1) and Table 3.
 
      data pm/ 1.0, 1.0, 1.0, 1.0, 1.0, 0.77, 0.5, 1.0, 1.0/
 
c-----fkw is the planck-weighted k-distribution function due to h2o
c     line absorption (Table 4).
c     The k-distribution function for Band 3, fkw(*,3), 
c     is not used (see the parameter gkw below).
 
      data fkw / 0.2747,0.2717,0.2752,0.1177,0.0352,0.0255,
     2           0.1521,0.3974,0.1778,0.1826,0.0374,0.0527,
     3           6*1.00,
     4           0.4654,0.2991,0.1343,0.0646,0.0226,0.0140,
     5           0.5543,0.2723,0.1131,0.0443,0.0160,0.0000,
     6           0.5955,0.2693,0.0953,0.0335,0.0064,0.0000,
     7           0.1958,0.3469,0.3147,0.1013,0.0365,0.0048,
     8           0.0740,0.1636,0.4174,0.1783,0.1101,0.0566,
     9           0.1437,0.2197,0.3185,0.2351,0.0647,0.0183/
 
c-----gkw is the planck-weighted k-distribution function due to h2o
c     line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
c     of band 3 (Table 10).  Note that the order of the sub-bands
c     is reversed.
 
      data gkw/  0.1782,0.0593,0.0215,0.0068,0.0022,0.0000,
     2           0.0923,0.1675,0.0923,0.0187,0.0178,0.0000,
     3           0.0000,0.1083,0.1581,0.0455,0.0274,0.0041/
 
 
c-----Coefficients for computing the extinction coefficient
c     for cloud ice particles (Table 11a, Eq. 6.4a).
c
      data aib /  -0.44171,    0.61222,   0.06465,
     2            -0.13727,    0.54102,   0.28962,
     3            -0.01878,    1.19270,   0.79080,
     4            -0.01896,    0.78955,   0.69493,
     5            -0.04788,    0.69729,   0.54492,
     6            -0.02265,    1.13370,   0.76161,
     7            -0.01038,    1.46940,   0.89045,
     8            -0.00450,    1.66240,   0.95989,
     9            -0.00044,    2.01500,   1.03750,
     *            -0.02956,    1.06430,   0.71283/
c
c-----coefficients for computing the extinction coefficient
c     for cloud liquid drops. (Table 11b, Eq. 6.4b)
c
      data awb /   0.08641,    0.01769,    -1.5572e-3,   3.4896e-5,
     2             0.22027,    0.00997,    -1.8719e-3,   5.3112e-5,
     3             0.38074,   -0.03027,     1.0154e-3,  -1.1849e-5,
     4             0.15587,    0.00371,    -7.7705e-4,   2.0547e-5,
     5             0.05518,    0.04544,    -4.2067e-3,   1.0184e-4,
     6             0.12724,    0.04751,    -5.2037e-3,   1.3711e-4,
     7             0.30390,    0.01656,    -3.5271e-3,   1.0828e-4,
     8             0.63617,   -0.06287,     2.2350e-3,  -2.3177e-5,
     9             1.15470,   -0.19282,     1.2084e-2,  -2.5612e-4,
     *             0.34021,   -0.02805,     1.0654e-3,  -1.5443e-5/
c
c-----coefficients for computing the single-scattering albedo
c     for cloud ice particles. (Table 12a, Eq. 6.5)
c
      data aiw/    0.17201,    1.8814e-2,  -3.5117e-4,   2.1127e-6,
     2             0.81470,   -4.1989e-3,   2.3152e-7,   2.0992e-7,
     3             0.54859,   -7.4266e-4,   1.2865e-5,  -5.7092e-8,
     4             0.39218,    6.4180e-3,  -1.1567e-4,   6.9710e-7,
     5             0.71773,   -5.1754e-3,   4.6658e-5,  -1.2085e-7,
     6             0.77345,   -8.4966e-3,   1.1451e-4,  -5.5170e-7,
     7             0.74975,   -8.7083e-3,   1.3367e-4,  -7.1603e-7,
     8             0.69011,   -6.9766e-3,   1.1674e-4,  -6.6472e-7,
     9             0.83963,   -1.0347e-2,   1.4651e-4,  -7.5965e-7,
     *             0.64860,   -4.4142e-3,   6.5458e-5,  -3.2655e-7/
 
c-----coefficients for computing the single-scattering albedo
c     for cloud liquid drops. (Table 12b, Eq. 6.5)
c
      data aww/   -7.8566e-2,  8.0875e-2,  -4.3403e-3,   8.1341e-5,
     2            -1.3384e-2,  9.3134e-2,  -6.0491e-3,   1.3059e-4,
     3             3.7096e-2,  7.3211e-2,  -4.4211e-3,   9.2448e-5,
     4            -3.7600e-3,  9.3344e-2,  -5.6561e-3,   1.1387e-4,
     5             0.40212,    7.8083e-2,  -5.9583e-3,   1.2883e-4,
     6             0.57928,    5.9094e-2,  -5.4425e-3,   1.2725e-4,
     7             0.68974,    4.2334e-2,  -4.9469e-3,   1.2863e-4,
     8             0.80122,    9.4578e-3,  -2.8508e-3,   9.0078e-5,
     9             1.02340,   -2.6204e-2,   4.2552e-4,   3.2160e-6,
     *             0.05092,    7.5409e-2,  -4.7305e-3,   1.0121e-4/ 
c
c-----coefficients for computing the asymmetry factor for cloud ice 
c     particles. (Table 13a, Eq. 6.6)
c
      data aig /   0.57867,    1.5592e-2,  -2.6372e-4,   1.5125e-6,
     2             0.72259,    4.7922e-3,  -4.7164e-5,   2.0400e-7,
     3             0.76109,    6.9922e-3,  -1.0935e-4,   5.9885e-7,
     4             0.86934,    4.2268e-3,  -7.4085e-5,   4.3547e-7,
     5             0.89103,    2.8482e-3,  -3.9174e-5,   2.0098e-7,
     6             0.86325,    3.2935e-3,  -3.9872e-5,   1.8015e-7,
     7             0.85064,    3.8505e-3,  -4.9259e-5,   2.3096e-7,
     8             0.86945,    3.7869e-3,  -5.6525e-5,   3.0016e-7,
     9             0.80122,    4.9086e-3,  -5.8831e-5,   2.6367e-7,
     *             0.73290,    7.3898e-3,  -1.0515e-4,   5.4034e-7/
c
c-----coefficients for computing the asymmetry factor for cloud liquid 
c     drops. (Table 13b, Eq. 6.6)
c
      data awg /  -0.51930,    0.20290,    -1.1747e-2,   2.3868e-4,
     2            -0.22151,    0.19708,    -1.2462e-2,   2.6646e-4,
     3             0.14157,    0.14705,    -9.5802e-3,   2.0819e-4,
     4             0.41590,    0.10482,    -6.9118e-3,   1.5115e-4,
     5             0.55338,    7.7016e-2,  -5.2218e-3,   1.1587e-4,
     6             0.61384,    6.4402e-2,  -4.6241e-3,   1.0746e-4,
     7             0.67891,    4.8698e-2,  -3.7021e-3,   9.1966e-5,
     8             0.78169,    2.0803e-2,  -1.4749e-3,   3.9362e-5,
     9             0.93218,   -3.3425e-2,   2.9632e-3,  -6.9362e-5,
     *             0.01649,    0.16561,    -1.0723e-2,   2.3220e-4/ 
c
c-----include tables used in the table look-up for co2 (band 3), 
c     o3 (band 5), and h2o (bands 1, 2, and 8) transmission functions.
c     "co2.tran4" is the co2 transmission table applicable to a large
c     range of co2 amount (up to 100 times of the present-time value).
 
      include "optics/h2o.tran3"
      include "optics/co2.tran4"
      include "optics/o3.tran3"
 
c-----compute layer pressure (pa) and layer temperature minus 250K (dt)
 
      do k=1,np
       do i=1,m
         pa(i,k)=0.5*(pl(i,k)+pl(i,k+1))
         dt(i,k)=ta(i,k)-250.0
       enddo
      enddo
 
c-----compute layer absorber amount
 
c     dh2o : water vapor amount (g/cm**2)
c     dcont: scaled water vapor amount for continuum absorption
c            (g/cm**2)
c     dco2 : co2 amount (cm-atm)stp
c     do3  : o3 amount (cm-atm)stp
c     dn2o : n2o amount (cm-atm)stp
c     dch4 : ch4 amount (cm-atm)stp
c     df11 : cfc11 amount (cm-atm)stp
c     df12 : cfc12 amount (cm-atm)stp
c     df22 : cfc22 amount (cm-atm)stp
c     the factor 1.02 is equal to 1000/980
c     factors 789 and 476 are for unit conversion
c     the factor 0.001618 is equal to 1.02/(.622*1013.25) 
c     the factor 6.081 is equal to 1800/296
 
      do k=1,np
       do i=1,m
 
         dp   (i,k) = pl(i,k+1)-pl(i,k)
 
         dh2o (i,k) = 1.02*wa(i,k)*dp(i,k)
         dh2o (i,k) = max(dh2o (i,k),1.e-10)
         do3  (i,k) = 476.*oa(i,k)*dp(i,k)
         do3 (i,k) = max(do3 (i,k),1.e-6)
         dco2 (i,k) = 789.*co2*dp(i,k)
         dco2 (i,k) = max(dco2 (i,k),1.e-4)
 
         dch4 (i,k) = 789.*ch4*dp(i,k)
         dn2o (i,k) = 789.*n2o*dp(i,k)
         df11 (i,k) = 789.*cfc11*dp(i,k)
         df12 (i,k) = 789.*cfc12*dp(i,k)
         df22 (i,k) = 789.*cfc22*dp(i,k)
 
c-----compute scaled water vapor amount for h2o continuum absorption
c     following eq. (4.21).
 
         xx=pa(i,k)*0.001618*wa(i,k)*wa(i,k)*dp(i,k)
         dcont(i,k) = xx*exp(1800./ta(i,k)-6.081)
 
       enddo
      enddo
 
c-----Set default values for reff.
c     Index is 1 for ice, 2 for waterdrops and 3 for raindrops.
 
        do k=1,np
         do i=1,m
          reff(i,k,1)=40.0
          reff(i,k,2)=10.0
        enddo
       enddo
 
c-----compute layer cloud water amount (gm/m**2)
 
       if (cldwater) then
        do k=1,np
         do i=1,m
             xx=1.02*10000.*(pl(i,k+1)-pl(i,k))
             cwp(i,k,1)=xx*cwc(i,k,1)
             cwp(i,k,2)=xx*cwc(i,k,2)
             cwp(i,k,3)=xx*cwc(i,k,3)
         enddo
        enddo
 
      do k=1,np
       do i=1,m
 
        if (cwp(i,k,1) .gt. 0.000001) then
 
c-----Compute effective radius of ice cloud particles following Equation (6.9)
 
          j=(ta(i,k)-193.)*0.1
          if (j.lt.1) j=1
          if (j.gt.5) j=5
 
c-----Conversion of the unit of cwc in g/g to the unit of x in g/m^3.
c     The constant 348.43 is equal to (100/0.287), 
c     where the constant 0.287 is related to the gas constant of dry air.
 
          x=cwc(i,k,1)*348.43*pa(i,k)/ta(i,k)
          x=log10(x)
          reff(i,k,1)=0.65*10.**(ai(j)+bi(j)*x+ci(j)*x*x+di(j)*x*x*x)
          reff(i,k,1)=max(reff(i,k,1),10.0)
          reff(i,k,1)=min(reff(i,k,1),70.0)
 
        endif
 
        if(cwp(i,k,2) .gt. 0.000001) then
 
c-----Effective radius of water cloud particles following Equation (6.13).
 
          x=cwc(i,k,2)*348.43*pa(i,k)/ta(i,k)
          reff(i,k,2)=14.3*x**0.1667
          reff(i,k,2)=max(reff(i,k,2),4.0)
          reff(i,k,2)=min(reff(i,k,2),20.0)
 
        endif
 
       enddo
      enddo
 
      endif
 
c-----the surface (np+1) is treated as a layer filled with black clouds.
c     transfc is the transmittance between the surface and a pressure level.
c     trantcr is the clear-sky transmittance between the surface and a
c     pressure level.
 
      do i=1,m
        sfcem(i)       =0.0
        transfc(i,np+1)=1.0
        trantcr(i,np+1)=1.0
      enddo
 
c-----initialize fluxes
 
      do k=1,np+1
       do i=1,m
         flx(i,k)  = 0.0
         flc(i,k)  = 0.0
         dfdts(i,k)= 0.0
       enddo
      enddo
 
c-----integration over spectral bands
 
      do 1000 ibn=1,10
 
       if (ibn.eq.10 .and. .not.trace) return
 
c-----if h2otbl, compute h2o (line) transmittance using table look-up.
c     if conbnd, compute h2o (continuum) transmittance in bands 2-7.
c     if co2bnd, compute co2 transmittance in band 3.
c     if oznbnd, compute  o3 transmittance in band 5.
c     if n2obnd, compute n2o transmittance in bands 6 and 7.
c     if ch4bnd, compute ch4 transmittance in bands 6 and 7.
c     if combnd, compute co2-minor transmittance in bands 4 and 5.
c     if f11bnd, compute cfc11 transmittance in bands 4 and 5.
c     if f12bnd, compute cfc12 transmittance in bands 4 and 6.
c     if f22bnd, compute cfc22 transmittance in bands 4 and 6.
c     if b10bnd, compute flux reduction due to n2o in band 10.
 
       h2otbl=high.and.(ibn.eq.1.or.ibn.eq.2.or.ibn.eq.8)
       conbnd=ibn.ge.2.and.ibn.le.7
       co2bnd=ibn.eq.3
       oznbnd=ibn.eq.5
       n2obnd=ibn.eq.6.or.ibn.eq.7
       ch4bnd=ibn.eq.6.or.ibn.eq.7
       combnd=ibn.eq.4.or.ibn.eq.5
       f11bnd=ibn.eq.4.or.ibn.eq.5
       f12bnd=ibn.eq.4.or.ibn.eq.6
       f22bnd=ibn.eq.4.or.ibn.eq.6
       b10bnd=ibn.eq.10
 
c-----blayer is the spectrally integrated planck flux of the mean layer
c     temperature derived from eq. (3.11)
c     The fitting for the planck flux is valid for the range 160-345 K.
 
       do k=1,np
 
        do i=1,m
          tx(i)=ta(i,k)
        enddo
            call planck(ibn,m,tx,xlayer)
 
        do i=1,m
          blayer(i,k)=xlayer(i)
        enddo
 
       enddo
 
c-----Index "0" is the layer above the top of the atmosphere.
 
        do i=1,m
          blayer(i,0)=0.0
        enddo
 
c-----Surface emission and reflectivity. See Section 9.
c     bs and dbs include the effect of surface emissivity.
 
          call sfcflux (ibn,m,ns,fs,tg,eg,tv,ev,rv,vege,bs,dbs,rflxs) 
 
         do i=1,m
          blayer(i,np+1)=bs(i)
         enddo
 
c------interpolate Planck function at model levels (linear in p)
 
       do k=2,np
        do i=1,m
         blevel(i,k)=(blayer(i,k-1)*dp(i,k)+blayer(i,k)*dp(i,k-1))/
     *               (dp(i,k-1)+dp(i,k))
        enddo
       enddo
 
c-----Extrapolate blevel(i,1) from blayer(i,2) and blayer(i,1)
 
       do i=1,m
         blevel(i,1)=blayer(i,1)+(blayer(i,1)-blayer(i,2))*dp(i,1)/
     *               (dp(i,1)+dp(i,2))
       enddo
 
c-----If the surface air temperature tb is known, compute blevel(i,np+1)
 
         call planck(ibn,m,tb,xlayer)
       do i=1,m
         blevel(i,np+1)=xlayer(i)
       enddo
 
c-----Otherwise, extrapolate blevel(np+1) from blayer(np-1) and blayer(np)
 
c      do i=1,m
c        blevel(i,np+1)=blayer(i,np)+(blayer(i,np)-blayer(i,np-1))
c    *                 *dp(i,np)/(dp(i,np)+dp(i,np-1))
c      enddo
 
c-----Compute cloud optical thickness following Eqs. (6.4a,b) and (6.7)
c     Rain optical thickness is set to 0.00307 /(gm/m**2).
c     It is for a specific drop size distribution provided by Q. Fu.
 
      if (cldwater) then
       do k=1,np
        do i=1,m
          taucl(i,k,1)=cwp(i,k,1)*(aib(1,ibn)+aib(2,ibn)/
     *      reff(i,k,1)**aib(3,ibn))
          taucl(i,k,2)=cwp(i,k,2)*(awb(1,ibn)+(awb(2,ibn)+
     *      (awb(3,ibn)+awb(4,ibn)*reff(i,k,2))*reff(i,k,2))
     *      *reff(i,k,2))
          taucl(i,k,3)=0.00307*cwp(i,k,3)
        enddo
       enddo
      endif
 
c-----Compute cloud single-scattering albedo and asymmetry factor for
c     a mixture of ice particles and liquid drops following 
c     Eqs. (6.5), (6.6), (6.15) and (6.16).
c     Single-scattering albedo and asymmetry factor of rain are set
c     to 0.54 and 0.95, respectively, based on the information provided
c     by Prof. Qiang Fu.
 
       do k=1,np
        do i=1,m
 
           tcldlyr(i,k) = 1.0
           tauc=taucl(i,k,1)+taucl(i,k,2)+taucl(i,k,3)
 
          if (tauc.gt.0.02 .and. fcld(i,k).gt.0.01) then
 
           w1=taucl(i,k,1)*(aiw(1,ibn)+(aiw(2,ibn)+(aiw(3,ibn)
     *       +aiw(4,ibn)*reff(i,k,1))*reff(i,k,1))*reff(i,k,1))
           w2=taucl(i,k,2)*(aww(1,ibn)+(aww(2,ibn)+(aww(3,ibn)
     *       +aww(4,ibn)*reff(i,k,2))*reff(i,k,2))*reff(i,k,2))
           w3=taucl(i,k,3)*0.54
           ww=(w1+w2+w3)/tauc
 
           g1=w1*(aig(1,ibn)+(aig(2,ibn)+(aig(3,ibn)
     *      +aig(4,ibn)*reff(i,k,1))*reff(i,k,1))*reff(i,k,1))
           g2=w2*(awg(1,ibn)+(awg(2,ibn)+(awg(3,ibn)
     *      +awg(4,ibn)*reff(i,k,2))*reff(i,k,2))*reff(i,k,2))
           g3=w3*0.95
 
           gg=(g1+g2+g3)/(w1+w2+w3)
 
c-----Parameterization of LW scattering following Eqs. (6.11) and (6.12). 
 
           ff=0.5+(0.3739+(0.0076+0.1185*gg)*gg)*gg
           tauc=(1.-ww*ff)*tauc
 
c-----compute cloud diffuse transmittance. It is approximated by using 
c     a diffusivity factor of 1.66.
 
           tcldlyr(i,k)=exp(-1.66*tauc)
 
          endif
 
        enddo
       enddo
 
c-----Compute optical thickness, single-scattering albedo and asymmetry
c     factor for a mixture of "na" aerosol types. Eqs. (7.1)-(7.3)
 
       if (aerosol) then
 
        do k=1,np
 
         do i=1,m
          taua(i,k)=0.0
          ssaa(i,k)=0.0
          asya(i,k)=0.0
         enddo
 
         do j=1,na
          do i=1,m
           taua(i,k)=taua(i,k)+taual(i,k,ibn,j)
           w1=ssaal(i,k,ibn,j)*taual(i,k,ibn,j)
           ssaa(i,k)=ssaa(i,k)+w1
           asya(i,k)=asya(i,k)+asyal(i,k,ibn,j)*w1
          enddo
         enddo
 
c-----taerlyr is the aerosol diffuse transmittance
 
         do i=1,m
           taerlyr(i,k)=1.0
 
          if (taua(i,k) .gt. 0.001) then 
          if (ssaa(i,k) .gt. 0.001) then
            asya(i,k)=asya(i,k)/ssaa(i,k)
            ssaa(i,k)=ssaa(i,k)/taua(i,k)
 
c-----Parameterization of aerosol scattering following Eqs. (6.11) and (6.12). 
 
           ff=0.5+(0.3739+(0.0076+0.1185*asya(i,k))*asya(i,k))*asya(i,k)
           taua(i,k)=taua(i,k)*(1.-ssaa(i,k)*ff)
 
          endif
           taerlyr(i,k)=exp(-1.66*taua(i,k))
          endif
 
         enddo
        enddo
 
       endif
 
c-----Compute the exponential terms (Eq. 8.21) at each layer due to
c     water vapor line absorption when k-distribution is used
 
      if (.not.h2otbl .and. .not.b10bnd) then
        call h2oexps(ibn,m,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
      endif
 
c-----compute the exponential terms (Eq. 4.24) at each layer due to
c     water vapor continuum absorption.
c     ne is the number of terms used in each band to compute water 
c     vapor continuum transmittance (Table 9).
 
        ne=0
      if (conbnd) then
 
        ne=1
        if (ibn.eq.3) ne=3
 
        call conexps(ibn,m,np,dcont,xke,conexp)
 
      endif
 
c-----compute the exponential terms (Eq. 8.21) at each layer due to
c     co2 absorption
 
      if (.not.high .and. co2bnd) then
        call co2exps(m,np,dco2,pa,dt,co2exp)
      endif
 
c***** for trace gases *****
 
      if (trace) then
 
c-----compute the exponential terms at each layer due to n2o absorption
 
       if (n2obnd) then
        call n2oexps(ibn,m,np,dn2o,pa,dt,n2oexp)
       endif
 
c-----compute the exponential terms at each layer due to ch4 absorption
 
       if (ch4bnd) then
        call ch4exps(ibn,m,np,dch4,pa,dt,ch4exp)
       endif
 
c-----Compute the exponential terms due to co2 minor absorption
 
       if (combnd) then
        call comexps(ibn,m,np,dco2,dt,comexp)
       endif
 
c-----Compute the exponential terms due to cfc11 absorption.
c     The values of the parameters are given in Table 7.
 
       if (f11bnd) then
            a1  = 1.26610e-3
            b1  = 3.55940e-6
            fk1 = 1.89736e+1
            a2  = 8.19370e-4
            b2  = 4.67810e-6
            fk2 = 1.01487e+1
        call cfcexps(ibn,m,np,a1,b1,fk1,a2,b2,fk2,df11,dt,f11exp)
       endif
 
c-----Compute the exponential terms due to cfc12 absorption.
 
       if (f12bnd) then
            a1  = 8.77370e-4
            b1  =-5.88440e-6
            fk1 = 1.58104e+1
            a2  = 8.62000e-4
            b2  =-4.22500e-6
            fk2 = 3.70107e+1
        call cfcexps(ibn,m,np,a1,b1,fk1,a2,b2,fk2,df12,dt,f12exp)
       endif
 
c-----Compute the exponential terms due to cfc22 absorption.
 
       if (f22bnd) then
            a1  = 9.65130e-4
            b1  = 1.31280e-5
            fk1 = 6.18536e+0
            a2  =-3.00010e-5 
            b2  = 5.25010e-7
            fk2 = 3.27912e+1
        call cfcexps(ibn,m,np,a1,b1,fk1,a2,b2,fk2,df22,dt,f22exp)
       endif
 
c-----Compute the exponential terms at each layer in band 10 due to
c     h2o line and continuum, co2, and n2o absorption
 
       if (b10bnd) then
        call b10exps(m,np,dh2o,dcont,dco2,dn2o,pa,dt
     *              ,h2oexp,conexp,co2exp,n2oexp)
       endif
 
      endif
 
c-----blayer(i,np+1) includes the effect of surface emissivity.
 
      do i=1,m
        bd(i,0)=0.0
        bu(i,np+1)=blayer(i,np+1)
      enddo
 
c-----do-loop 1500 is for computing upward (bu) and downward (bd)
c     emission of a layer following Eqs. (8.17), (8.18), (8.19).
c     Here, trant(i,k2) is the transmittance of the layer k2-1.
 
      do 1500 k2=2,np+1
 
c-----for h2o line transmission
 
      if (.not. h2otbl) then
        do ik=1,6
         do i=1,m
           th2o(i,ik)=1.0
         enddo
        enddo
      endif
 
c-----for h2o continuum transmission
 
         do iq=1,3
          do i=1,m
            tcon(i,iq)=1.0
          enddo
         enddo
 
c-----for co2 transmission using k-distribution method.
c     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
c     are combined in computing the co2 transmittance.
 
       if (.not.high .and. co2bnd) then
         do isb=1,2
          do ik=1,6
           do i=1,m
             tco2(i,ik,isb)=1.0
           enddo
          enddo
         enddo
       endif
 
 
       do i=1,m
        x1(i)=0.0
        x2(i)=0.0
        x3(i)=0.0
        trant(i,k2)=1.0
       enddo
 
      if (h2otbl) then
 
c-----Compute water vapor transmittance using table look-up.
c     The following values are taken from Table 8.
 
          w1=-8.0
          p1=-2.0
          dwe=0.3
          dpe=0.2
 
          if (ibn.eq.1) then
           call tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3,
     *                 w1,p1,dwe,dpe,h11,h12,h13,trant)
 
          endif
          if (ibn.eq.2) then
           call tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3,
     *                 w1,p1,dwe,dpe,h21,h22,h23,trant)
 
          endif
          if (ibn.eq.8) then
           call tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3,
     *                 w1,p1,dwe,dpe,h81,h82,h83,trant)
          endif
 
c-----for water vapor continuum absorption
 
           if (conbnd) then
            do i=1,m
             tcon(i,1)=tcon(i,1)*conexp(i,k2-1,1)
             trant(i,k2)=trant(i,k2)*tcon(i,1)
            enddo
           endif
 
      else
 
c-----compute water vapor transmittance using k-distribution
 
       if (.not.b10bnd) then
        call h2okdis(ibn,m,np,k2-1,fkw,gkw,ne,h2oexp,conexp,
     *               th2o,tcon,trant)
 
       endif
 
      endif
 
      if (co2bnd) then
 
        if (high) then
 
c-----Compute co2 transmittance using table look-up method.
c     The following values are taken from Table 8.
 
          w1=-4.0
          p1=-2.0
          dwe=0.3
          dpe=0.2
          call tablup(k2,m,np,nx,nc,dco2,pa,dt,x1,x2,x3,
     *                w1,p1,dwe,dpe,c1,c2,c3,trant)
        else
 
c-----compute co2 transmittance using k-distribution method
          call co2kdis(m,np,k2-1,co2exp,tco2,trant)
 
        endif
 
      endif
 
c-----Always use table look-up to compute o3 transmittance.
c     The following values are taken from Table 8.
 
      if (oznbnd) then
          w1=-6.0
          p1=-2.0
          dwe=0.3
          dpe=0.2
          call tablup(k2,m,np,nx,no,do3,pa,dt,x1,x2,x3,
     *                w1,p1,dwe,dpe,o1,o2,o3,trant)
      endif
 
c-----include aerosol effect
 
      if (aerosol) then
       do i=1,m
         trant(i,k2)=trant(i,k2)*taerlyr(i,k2-1)
       enddo
      endif
 
c-----Compute upward and downward emission of the layer k2-1
 
        do i=1,m
 
          xx=(blayer(i,k2-1)-blevel(i,k2-1))*
     *       (blayer(i,k2-1)-blevel(i,k2))
 
         if (xx.gt.0.0) then
 
c-----If xx>0, there is a local temperature minimum or maximum.
c     Computations of bd and bu follow Eq. (8.20).
 
          bd(i,k2-1)=.5*blayer(i,k2-1)+.25*(blevel(i,k2-1)+blevel(i,k2))
          bu(i,k2-1)=bd(i,k2-1)
 
         else
 
c-----Computations of bd and bu following Eqs.(8.17) and (8.18).
c     The effect of clouds on the transmission of a layer is taken
c     into account, following Eq. (8.19).
 
          xx=(fcld(i,k2-1)*tcldlyr(i,k2-1)+(1.-fcld(i,k2-1)))
     *       *trant(i,k2)
 
          yy=min(0.9999,xx)
          yy=max(0.00001,yy)
          xx=(blevel(i,k2-1)-blevel(i,k2))/alog(yy)
          bd(i,k2-1)=(blevel(i,k2)-blevel(i,k2-1)*yy)/(1.0-yy)-xx
          bu(i,k2-1)=(blevel(i,k2-1)+blevel(i,k2))-bd(i,k2-1)
 
         endif
 
        enddo
 
 1500 continue
 
c-----initialize fluxes
 
      do k=1,np+1
       do i=1,m
         flxu(i,k) = 0.0
         flxd(i,k) = 0.0
         flcu(i,k) = 0.0
         flcd(i,k) = 0.0
       enddo
      enddo
 
c-----
 
      do 2000 k1=1,np
 
c-----initialization
c
c     it, im, and ib are the numbers of cloudy layers in the high,
c     middle, and low cloud groups between levels k1 and k2.
c     cldlw, cldmd, and cldhi are the equivalent black-cloud fractions
c     of low, middle, and high troposphere.
c     tranal is the aerosol transmission function
 
        do i=1,m
          it(i) = 0
          im(i) = 0
          ib(i) = 0
          cldlw(i) = 0.0
          cldmd(i) = 0.0
          cldhi(i) = 0.0
          tranal(i)= 1.0
        enddo
 
c-----for h2o line transmission
 
      if (.not. h2otbl) then
        do ik=1,6
         do i=1,m
           th2o(i,ik)=1.0
         enddo
        enddo
      endif
 
c-----for h2o continuum transmission
 
         do iq=1,3
          do i=1,m
            tcon(i,iq)=1.0
          enddo
         enddo
 
c-----for co2 transmission using k-distribution method.
c     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
c     are combined in computing the co2 transmittance.
 
       if (.not.high .and. co2bnd) then
         do isb=1,2
          do ik=1,6
           do i=1,m
             tco2(i,ik,isb)=1.0
           enddo
          enddo
         enddo
       endif
 
c***** for trace gases *****
 
      if (trace) then
 
c-----for n2o transmission using k-distribution method.
 
       if (n2obnd) then
          do ik=1,4
           do i=1,m
             tn2o(i,ik)=1.0
           enddo
          enddo
       endif
 
c-----for ch4 transmission using k-distribution method.
 
       if (ch4bnd) then
          do ik=1,4
           do i=1,m
             tch4(i,ik)=1.0
           enddo
          enddo
       endif
 
c-----for co2-minor transmission using k-distribution method.
 
       if (combnd) then
          do ik=1,6
           do i=1,m
             tcom(i,ik)=1.0
           enddo
          enddo
       endif
 
c-----for cfc-11 transmission using k-distribution method.
 
       if (f11bnd) then
           do i=1,m
             tf11(i)=1.0
           enddo
       endif
 
c-----for cfc-12 transmission using k-distribution method.
 
       if (f12bnd) then
           do i=1,m
             tf12(i)=1.0
           enddo
       endif
 
c-----for cfc-22 transmission when using k-distribution method.
 
       if (f22bnd) then
           do i=1,m
             tf22(i)=1.0
           enddo
       endif
 
c-----for the transmission in band 10 using k-distribution method.
 
       if (b10bnd) then
          do ik=1,5
           do i=1,m
              th2o(i,ik)=1.0
           enddo
          enddo
 
          do ik=1,6
           do i=1,m
              tco2(i,ik,1)=1.0
           enddo
          enddo
 
          do i=1,m
             tcon(i,1)=1.0
          enddo
 
          do ik=1,2
            do i=1,m
              tn2o(i,ik)=1.0
           enddo
          enddo
       endif
 
      endif
 
c***** end trace gases *****
 
       do i=1,m
        x1(i)=0.0
        x2(i)=0.0
        x3(i)=0.0
       enddo
 
c-----trant is the total transmittance between levels k1 and k2.
c     fclr is the clear line-of-sight  between levels k1 and k2.
 
       do k=1,np+1
        do i=1,m
          trant(i,k)=1.0
          fclr(i,k) =1.0
        enddo
       enddo
 
c-----do-loop 3000 are for computing (a) transmittance, trant(i,k2),
c     and (b) clear line-of-sight, fclr(i,k2), between levels k1 and k2.
 
      do 3000 k2=k1+1,np+1
 
      if (h2otbl) then
 
c-----Compute water vapor transmittance using table look-up.
c     The following values are taken from Table 8.
 
          w1=-8.0
          p1=-2.0
          dwe=0.3
          dpe=0.2
 
          if (ibn.eq.1) then
           call tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3,
     *                 w1,p1,dwe,dpe,h11,h12,h13,trant)
 
          endif
          if (ibn.eq.2) then
           call tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3,
     *                 w1,p1,dwe,dpe,h21,h22,h23,trant)
 
          endif
          if (ibn.eq.8) then
           call tablup(k2,m,np,nx,nh,dh2o,pa,dt,x1,x2,x3,
     *                 w1,p1,dwe,dpe,h81,h82,h83,trant)
          endif
 
           if (conbnd) then
            do i=1,m
             tcon(i,1)=tcon(i,1)*conexp(i,k2-1,1)
             trant(i,k2)=trant(i,k2)*tcon(i,1)
            enddo
           endif
 
      else
 
c-----compute water vapor transmittance using k-distribution
 
       if (.not.b10bnd) then
        call h2okdis(ibn,m,np,k2-1,fkw,gkw,ne,h2oexp,conexp,
     *               th2o,tcon,trant)
       endif
 
      endif
 
      if (co2bnd) then
 
        if (high) then
 
c-----Compute co2 transmittance using table look-up method.
c     The following values are taken from Table 8.
 
          w1=-4.0
          p1=-2.0
          dwe=0.3
          dpe=0.2
          call tablup(k2,m,np,nx,nc,dco2,pa,dt,x1,x2,x3,
     *                w1,p1,dwe,dpe,c1,c2,c3,trant)
       else
 
c-----compute co2 transmittance using k-distribution method
          call co2kdis(m,np,k2-1,co2exp,tco2,trant)
 
        endif
 
      endif
 
c-----Always use table look-up to compute o3 transmittance.
c     The following values are taken from Table 8.
 
      if (oznbnd) then
          w1=-6.0
          p1=-2.0
          dwe=0.3
          dpe=0.2
          call tablup(k2,m,np,nx,no,do3,pa,dt,x1,x2,x3,
     *                w1,p1,dwe,dpe,o1,o2,o3,trant)
      endif
 
c***** for trace gases *****
 
      if (trace) then
 
c-----compute n2o transmittance using k-distribution method
 
       if (n2obnd) then
          call n2okdis(ibn,m,np,k2-1,n2oexp,tn2o,trant)
       endif
 
c-----compute ch4 transmittance using k-distribution method
 
       if (ch4bnd) then
          call ch4kdis(ibn,m,np,k2-1,ch4exp,tch4,trant)
       endif
 
c-----compute co2-minor transmittance using k-distribution method
 
       if (combnd) then
          call comkdis(ibn,m,np,k2-1,comexp,tcom,trant)
       endif
 
c-----compute cfc11 transmittance using k-distribution method
 
       if (f11bnd) then
          call cfckdis(m,np,k2-1,f11exp,tf11,trant)
       endif
 
c-----compute cfc12 transmittance using k-distribution method
 
       if (f12bnd) then
          call cfckdis(m,np,k2-1,f12exp,tf12,trant)
       endif
 
c-----compute cfc22 transmittance using k-distribution method
 
       if (f22bnd) then
          call cfckdis(m,np,k2-1,f22exp,tf22,trant)
       endif
 
c-----Compute transmittance in band 10 using k-distribution method.
c     For band 10, trant is the change in transmittance due to n2o 
c     absorption.
 
       if (b10bnd) then
          call b10kdis(m,np,k2-1,h2oexp,conexp,co2exp,n2oexp
     *                ,th2o,tcon,tco2,tn2o,trant)
 
       endif
 
      endif
 
c*****   end trace gases  *****
 
c-----include aerosol effect
 
      if (aerosol) then
       do i=1,m
         tranal(i)=tranal(i)*taerlyr(i,k2-1)
         trant(i,k2)=trant(i,k2) *tranal(i)
       enddo
      endif
 
c***** cloud overlapping *****
 
      if (.not. overcast) then
        call cldovlp (m,np,k2,ict,icb,it,im,ib,itx,imx,ibx,
     *           cldhi,cldmd,cldlw,fcld,tcldlyr,fclr)
 
      else
 
       do i=1,m
        fclr(i,k2)=fclr(i,k2)*tcldlyr(i,k2-1)
       enddo
 
      endif
 
 3000 continue
 
c-----do-loop 4000 is for computing upward and downward fluxes
c     for each spectral band
c     flcu, flcd: clear-sky upward and downward fluxes
c     flxu, flxd: all-sky   upward and downward fluxes
 
       do 4000 k2=k1+1,np+1
 
        if (k2.eq.k1+1 .and. ibn .ne. 10) then
 
c-----The first terms on the rhs of Eqs. (8.15) and (8.16)
 
         do i=1,m
          flcu(i,k1)=flcu(i,k1)-bu(i,k1)
          flcd(i,k2)=flcd(i,k2)+bd(i,k1)
          flxu(i,k1)=flxu(i,k1)-bu(i,k1)
          flxd(i,k2)=flxd(i,k2)+bd(i,k1)
         enddo
 
        endif
 
c-----The summation terms on the rhs of Eqs. (8.15) and (8.16).
c     Also see Eqs. (5.4) and (5.5) for Band 10.
 
         do i=1,m
          xx=trant(i,k2)*(bu(i,k2-1)-bu(i,k2))
          flcu(i,k1) =flcu(i,k1)+xx
          flxu(i,k1) =flxu(i,k1)+xx*fclr(i,k2)
          xx=trant(i,k2)*(bd(i,k1-1)-bd(i,k1))
          flcd(i,k2) =flcd(i,k2)+xx
          flxd(i,k2) =flxd(i,k2)+xx*fclr(i,k2)
         enddo
 
 4000 continue
 
c-----Here, fclr and trant are, respectively, the clear line-of-sight 
c     and the transmittance between k1 and the surface.
 
       do i=1,m
         trantcr(i,k1) =trant(i,np+1)
         transfc(i,k1) =trant(i,np+1)*fclr(i,np+1)
       enddo
 
c-----compute the partial derivative of fluxes with respect to
c     surface temperature (Eq. 3.12). 
c     Note: upward flux is negative, and so is dfdts.
 
       do i=1,m
         dfdts(i,k1) =dfdts(i,k1)-dbs(i)*transfc(i,k1)
       enddo
 
 2000 continue
 
      if (.not. b10bnd) then
 
c-----For surface emission.
c     Note: blayer(i,np+1) and dbs include the surface emissivity effect.
c     Both dfdts and sfcem are negative quantities.
 
        do i=1,m
          flcu(i,np+1)=-blayer(i,np+1)
          flxu(i,np+1)=-blayer(i,np+1)
          sfcem(i)=sfcem(i)-blayer(i,np+1)
          dfdts(i,np+1)=dfdts(i,np+1)-dbs(i)
        enddo
 
 
c-----Add the flux reflected by the surface. (Second term on the
c     rhs of Eq. 8.16)
 
        do k=1,np+1
         do i=1,m
           flcu(i,k)=flcu(i,k)-
     *          flcd(i,np+1)*trantcr(i,k)*rflxs(i)
           flxu(i,k)=flxu(i,k)-                         
     *          flxd(i,np+1)*transfc(i,k)*rflxs(i)   
         enddo
        enddo
 
      endif
 
c-----Summation of fluxes over spectral bands
 
       do k=1,np+1
        do i=1,m
          flc(i,k)=flc(i,k)+flcd(i,k)+flcu(i,k)
          flx(i,k)=flx(i,k)+flxd(i,k)+flxu(i,k)
        enddo
       enddo
 
 1000 continue
 
      return
      end
 
c***********************************************************************
      subroutine planck(ibn,m,t,xlayer)
c***********************************************************************
c
c-----Compute spectrally integrated Planck flux
c
      implicit none
 
      integer ibn                   ! spectral band index
      integer m                     ! no of points
      real t(m)                     ! temperature (K)
      real xlayer(m)                ! planck flux (w/m2)
      real cb(6,10)
      integer i
 
c-----the following coefficients are given in Table 2 for computing  
c     spectrally integrated planck fluxes using Eq. (3.11)
 
       data cb/
     1      5.3443e+0,  -2.0617e-1,   2.5333e-3,
     1     -6.8633e-6,   1.0115e-8,  -6.2672e-12,
     2      2.7148e+1,  -5.4038e-1,   2.9501e-3,
     2      2.7228e-7,  -9.3384e-9,   9.9677e-12,
     3     -3.4860e+1,   1.1132e+0,  -1.3006e-2,
     3      6.4955e-5,  -1.1815e-7,   8.0424e-11,
     4     -6.0513e+1,   1.4087e+0,  -1.2077e-2,
     4      4.4050e-5,  -5.6735e-8,   2.5660e-11,
     5     -2.6689e+1,   5.2828e-1,  -3.4453e-3,
     5      6.0715e-6,   1.2523e-8,  -2.1550e-11,
     6     -6.7274e+0,   4.2256e-2,   1.0441e-3,
     6     -1.2917e-5,   4.7396e-8,  -4.4855e-11,
     7      1.8786e+1,  -5.8359e-1,   6.9674e-3,
     7     -3.9391e-5,   1.0120e-7,  -8.2301e-11,
     8      1.0344e+2,  -2.5134e+0,   2.3748e-2,
     8     -1.0692e-4,   2.1841e-7,  -1.3704e-10,
     9     -1.0482e+1,   3.8213e-1,  -5.2267e-3,
     9      3.4412e-5,  -1.1075e-7,   1.4092e-10,
     *      1.6769e+0,   6.5397e-2,  -1.8125e-3,
     *      1.2912e-5,  -2.6715e-8,   1.9792e-11/
c
      do i=1,m
         xlayer(i)=t(i)*(t(i)*(t(i)*(t(i)*(t(i)*cb(6,ibn)+cb(5,ibn))
     *            +cb(4,ibn))+cb(3,ibn))+cb(2,ibn))+cb(1,ibn)
      enddo
 
      return
      end
 
c***********************************************************************
      subroutine plancd(ibn,m,t,dbdt) 
c***********************************************************************
c
c-----Compute the derivative of Planck flux wrt temperature
c
      implicit none
 
      integer ibn               ! spectral band index
      integer m                 ! no of points
      real t(m)                 ! temperature (K)
      real dbdt(m)              ! derivative of Planck flux wrt temperature
      real dcb(5,10)
      integer i
 
c-----Coefficients for computing the derivative of Planck function
c     with respect to temperature (Eq. 3.12).
c     dcb(1)=1*cb(2), dcb(2)=2*cb(3), dcb(3)=3*cb(4) ...  etc
 
       data dcb/
     1  -2.0617E-01, 5.0666E-03,-2.0590E-05, 4.0460E-08,-3.1336E-11,
     2  -5.4038E-01, 5.9002E-03, 8.1684E-07,-3.7354E-08, 4.9839E-11,
     3   1.1132E+00,-2.6012E-02, 1.9486E-04,-4.7260E-07, 4.0212E-10,
     4   1.4087E+00,-2.4154E-02, 1.3215E-04,-2.2694E-07, 1.2830E-10,
     5   5.2828E-01,-6.8906E-03, 1.8215E-05, 5.0092E-08,-1.0775E-10,
     6   4.2256E-02, 2.0882E-03,-3.8751E-05, 1.8958E-07,-2.2428E-10,
     7  -5.8359E-01, 1.3935E-02,-1.1817E-04, 4.0480E-07,-4.1150E-10,
     8  -2.5134E+00, 4.7496E-02,-3.2076E-04, 8.7364E-07,-6.8520E-10,
     9   3.8213E-01,-1.0453E-02, 1.0324E-04,-4.4300E-07, 7.0460E-10,
     *   6.5397E-02,-3.6250E-03, 3.8736E-05,-1.0686E-07, 9.8960E-11/
c
      do i=1,m
         dbdt(i)=t(i)*(t(i)*(t(i)*(t(i)*dcb(5,ibn)+dcb(4,ibn))
     *          +dcb(3,ibn))+dcb(2,ibn))+dcb(1,ibn)
      enddo
 
      return
      end
 
c**********************************************************************
      subroutine h2oexps(ib,m,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
c**********************************************************************
c   Compute exponentials for water vapor line absorption
c   in individual layers using Eqs. (8.21) and (8.22).
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer water vapor amount for line absorption (dh2o) 
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c  absorption coefficients for the first k-distribution
c     function due to h2o line absorption (xkw)
c  coefficients for the temperature and pressure scaling (aw,bw,pm)
c  ratios between neighboring absorption coefficients for
c     h2o line absorption (mw)
c
c---- output parameters
c  6 exponentials for each layer  (h2oexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k,ik
 
c---- input parameters ------
 
      real dh2o(m,np),pa(m,np),dt(m,np)
 
c---- output parameters -----
 
      real h2oexp(m,np,6)
 
c---- static data -----
 
      integer mw(9)
      real xkw(9),aw(9),bw(9),pm(9)
 
c---- temporary arrays -----
 
      real xh
 
c**********************************************************************
c    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
c    and bw,  therefore, h2oexp for these sub-bands are identical.
c**********************************************************************
 
        do k=1,np
         do i=1,m
 
c-----xh is the scaled water vapor amount for line absorption
c     computed from Eq. (4.4).
 
           xh = dh2o(i,k)*(pa(i,k)/500.)**pm(ib)
     1        * ( 1.+(aw(ib)+bw(ib)* dt(i,k))*dt(i,k) )
 
c-----h2oexp is the water vapor transmittance of the layer k
c     due to line absorption
 
           h2oexp(i,k,1) = exp(-xh*xkw(ib))
 
         enddo
        enddo
 
c-----compute transmittances from Eq. (8.22)
 
        do ik=2,6
 
         if (mw(ib).eq.6) then
 
          do k=1,np
           do i=1,m
             xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             h2oexp(i,k,ik) = xh*xh*xh
           enddo
          enddo
 
        elseif (mw(ib).eq.8) then
 
          do k=1,np
           do i=1,m
             xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             xh = xh*xh
             h2oexp(i,k,ik) = xh*xh
           enddo
          enddo
 
        elseif (mw(ib).eq.9) then
 
          do k=1,np
           do i=1,m
             xh=h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             h2oexp(i,k,ik) = xh*xh*xh
           enddo
          enddo
 
        else
 
          do k=1,np
           do i=1,m
             xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             xh = xh*xh
             xh = xh*xh
             h2oexp(i,k,ik) = xh*xh
           enddo
          enddo
 
        endif
       enddo
 
      return
      end
 
c**********************************************************************
      subroutine conexps(ib,m,np,dcont,xke,conexp)
c**********************************************************************
c   compute exponentials for continuum absorption in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer scaled water vapor amount for continuum absorption (dcont) 
c  absorption coefficients for the first k-distribution function
c     due to water vapor continuum absorption (xke)
c
c---- output parameters
c  1 or 3 exponentials for each layer (conexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters ------
 
      real dcont(m,np)
 
c---- updated parameters -----
 
      real conexp(m,np,3)
 
c---- static data -----
 
      real xke(9)
 
c****************************************************************
 
        do k=1,np
         do i=1,m
           conexp(i,k,1) = exp(-dcont(i,k)*xke(ib))
         enddo
        enddo
 
       if (ib .eq. 3) then
 
c-----The absorption coefficients for sub-bands 3b and 3a are, respectively,
c     two and four times the absorption coefficient for sub-band 3c (Table 9).
c     Note that conexp(i,k,3) is for sub-band 3a. 
 
         do k=1,np
          do i=1,m
            conexp(i,k,2) = conexp(i,k,1) *conexp(i,k,1)
            conexp(i,k,3) = conexp(i,k,2) *conexp(i,k,2)
          enddo
         enddo
 
       endif
 
      return
      end
 
c**********************************************************************
      subroutine co2exps(m,np,dco2,pa,dt,co2exp)
c**********************************************************************
c   Compute co2 exponentials for individual layers.
c
c---- input parameters
c  number of grid intervals (m)
c  number of layers (np)
c  layer co2 amount (dco2)
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  6 exponentials for each layer (co2exp)
c**********************************************************************
      implicit none
      integer m,np,i,k
 
c---- input parameters -----
 
      real dco2(m,np),pa(m,np),dt(m,np)
 
c---- output parameters -----
 
      real co2exp(m,np,6,2)
 
c---- temporary arrays -----
 
      real xc
 
c**********************************************************************
 
        do k=1,np
         do i=1,m
 
c-----The scaling parameters are given in Table 3, and values of
c     the absorption coefficient are given in Table 10.
 
c     Scaled co2 amount for band-wings (sub-bands 3a and 3c)
 
           xc = dco2(i,k)*(pa(i,k)/300.0)**0.5
     1             *(1.+(0.0182+1.07e-4*dt(i,k))*dt(i,k))
 
c-----six exponentials by powers of 8 (See Eqs. 8.21, 8.22 and Table 10).
 
           co2exp(i,k,1,1)=exp(-xc*2.656e-5)
 
           xc=co2exp(i,k,1,1)*co2exp(i,k,1,1)
           xc=xc*xc
           co2exp(i,k,2,1)=xc*xc
 
           xc=co2exp(i,k,2,1)*co2exp(i,k,2,1)
           xc=xc*xc
           co2exp(i,k,3,1)=xc*xc
 
           xc=co2exp(i,k,3,1)*co2exp(i,k,3,1)
           xc=xc*xc
           co2exp(i,k,4,1)=xc*xc
 
           xc=co2exp(i,k,4,1)*co2exp(i,k,4,1)
           xc=xc*xc
           co2exp(i,k,5,1)=xc*xc
 
           xc=co2exp(i,k,5,1)*co2exp(i,k,5,1)
           xc=xc*xc
           co2exp(i,k,6,1)=xc*xc
 
c-----For band-center region (sub-band 3b)
 
           xc = dco2(i,k)*(pa(i,k)/30.0)**0.85
     1             *(1.+(0.0042+2.00e-5*dt(i,k))*dt(i,k))
 
           co2exp(i,k,1,2)=exp(-xc*2.656e-3)
 
           xc=co2exp(i,k,1,2)*co2exp(i,k,1,2)
           xc=xc*xc
           co2exp(i,k,2,2)=xc*xc
 
           xc=co2exp(i,k,2,2)*co2exp(i,k,2,2)
           xc=xc*xc
           co2exp(i,k,3,2)=xc*xc
 
           xc=co2exp(i,k,3,2)*co2exp(i,k,3,2)
           xc=xc*xc
           co2exp(i,k,4,2)=xc*xc
 
           xc=co2exp(i,k,4,2)*co2exp(i,k,4,2)
           xc=xc*xc
           co2exp(i,k,5,2)=xc*xc
 
           xc=co2exp(i,k,5,2)*co2exp(i,k,5,2)
           xc=xc*xc
           co2exp(i,k,6,2)=xc*xc
 
         enddo
        enddo
 
      return
      end
 
c**********************************************************************
      subroutine n2oexps(ib,m,np,dn2o,pa,dt,n2oexp)
c**********************************************************************
c   Compute n2o exponentials for individual layers 
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer n2o amount (dn2o)
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  2 or 4 exponentials for each layer (n2oexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters -----
 
      real dn2o(m,np),pa(m,np),dt(m,np)
 
c---- output parameters -----
 
      real n2oexp(m,np,4)
 
c---- temporary arrays -----
 
      real xc,xc1,xc2
 
c-----Scaling and absorption data are given in Table 5.
c     Transmittances are computed using Eqs. (8.21) and (8.22).
 
       do k=1,np
        do i=1,m
 
c-----four exponential by powers of 21 for band 6.
 
          if (ib.eq.6) then
 
           xc=dn2o(i,k)*(1.+(1.9297e-3+4.3750e-6*dt(i,k))*dt(i,k))
           n2oexp(i,k,1)=exp(-xc*6.31582e-2)
 
           xc=n2oexp(i,k,1)*n2oexp(i,k,1)*n2oexp(i,k,1)
           xc1=xc*xc
           xc2=xc1*xc1
           n2oexp(i,k,2)=xc*xc1*xc2
 
c-----four exponential by powers of 8 for band 7
 
          else
 
           xc=dn2o(i,k)*(pa(i,k)/500.0)**0.48
     *        *(1.+(1.3804e-3+7.4838e-6*dt(i,k))*dt(i,k))
           n2oexp(i,k,1)=exp(-xc*5.35779e-2)
 
           xc=n2oexp(i,k,1)*n2oexp(i,k,1)
           xc=xc*xc
           n2oexp(i,k,2)=xc*xc
           xc=n2oexp(i,k,2)*n2oexp(i,k,2)
           xc=xc*xc
           n2oexp(i,k,3)=xc*xc
           xc=n2oexp(i,k,3)*n2oexp(i,k,3)
           xc=xc*xc
           n2oexp(i,k,4)=xc*xc
 
          endif
 
        enddo
       enddo
 
      return
      end
 
c**********************************************************************
      subroutine ch4exps(ib,m,np,dch4,pa,dt,ch4exp)
c**********************************************************************
c   Compute ch4 exponentials for individual layers
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer ch4 amount (dch4)
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  1 or 4 exponentials for each layer (ch4exp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters -----
 
      real dch4(m,np),pa(m,np),dt(m,np)
 
c---- output parameters -----
 
      real ch4exp(m,np,4)
 
c---- temporary arrays -----
 
      real xc
 
c*****  Scaling and absorption data are given in Table 5  *****
 
       do k=1,np
        do i=1,m
 
c-----four exponentials for band 6
 
          if (ib.eq.6) then
 
           xc=dch4(i,k)*(1.+(1.7007e-2+1.5826e-4*dt(i,k))*dt(i,k))
           ch4exp(i,k,1)=exp(-xc*5.80708e-3)
 
c-----four exponentials by powers of 12 for band 7
 
          else
 
           xc=dch4(i,k)*(pa(i,k)/500.0)**0.65
     *       *(1.+(5.9590e-4-2.2931e-6*dt(i,k))*dt(i,k))
           ch4exp(i,k,1)=exp(-xc*6.29247e-2)
 
           xc=ch4exp(i,k,1)*ch4exp(i,k,1)*ch4exp(i,k,1)
           xc=xc*xc
           ch4exp(i,k,2)=xc*xc
 
           xc=ch4exp(i,k,2)*ch4exp(i,k,2)*ch4exp(i,k,2)
           xc=xc*xc
           ch4exp(i,k,3)=xc*xc
 
           xc=ch4exp(i,k,3)*ch4exp(i,k,3)*ch4exp(i,k,3)
           xc=xc*xc
           ch4exp(i,k,4)=xc*xc
 
          endif
 
        enddo
       enddo
 
      return
      end
 
c**********************************************************************
      subroutine comexps(ib,m,np,dcom,dt,comexp)
c**********************************************************************
c   Compute co2-minor exponentials for individual layers using 
c   Eqs. (8.21) and (8.22).
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  layer co2 amount (dcom)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  6 exponentials for each layer (comexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k,ik
 
c---- input parameters -----
 
      real dcom(m,np),dt(m,np)
 
c---- output parameters -----
 
      real comexp(m,np,6)
 
c---- temporary arrays -----
 
      real xc
 
c*****  Scaling and absorpton data are given in Table 6  *****
 
       do k=1,np
        do i=1,m
 
          if (ib.eq.4) then
           xc=dcom(i,k)*(1.+(3.5775e-2+4.0447e-4*dt(i,k))*dt(i,k))
          endif
 
          if (ib.eq.5) then
           xc=dcom(i,k)*(1.+(3.4268e-2+3.7401e-4*dt(i,k))*dt(i,k))
          endif
 
           comexp(i,k,1)=exp(-xc*1.922e-7)
 
          do ik=2,6
           xc=comexp(i,k,ik-1)*comexp(i,k,ik-1)
           xc=xc*xc
           comexp(i,k,ik)=xc*comexp(i,k,ik-1)
          enddo
 
        enddo
       enddo
 
      return
      end
 
c**********************************************************************
      subroutine cfcexps(ib,m,np,a1,b1,fk1,a2,b2,fk2,dcfc,dt,cfcexp)
c**********************************************************************
c   compute cfc(-11, -12, -22) exponentials for individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of layers (np)
c  parameters for computing the scaled cfc amounts
c             for temperature scaling (a1,b1,a2,b2)
c  the absorption coefficients for the
c     first k-distribution function due to cfcs (fk1,fk2)
c  layer cfc amounts (dcfc)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  1 exponential for each layer (cfcexp)
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters -----
 
      real dcfc(m,np),dt(m,np)
 
c---- output parameters -----
 
      real cfcexp(m,np)
 
c---- static data -----
 
      real a1,b1,fk1,a2,b2,fk2
 
c---- temporary arrays -----
 
      real xf
 
c**********************************************************************
 
       do k=1,np
        do i=1,m
 
c-----compute the scaled cfc amount (xf) and exponential (cfcexp)
 
          if (ib.eq.4) then
           xf=dcfc(i,k)*(1.+(a1+b1*dt(i,k))*dt(i,k))
           cfcexp(i,k)=exp(-xf*fk1)
          else
           xf=dcfc(i,k)*(1.+(a2+b2*dt(i,k))*dt(i,k))
           cfcexp(i,k)=exp(-xf*fk2)
          endif
 
        enddo
       enddo
 
      return
      end
 
c**********************************************************************
      subroutine b10exps(m,np,dh2o,dcont,dco2,dn2o,pa,dt
     *          ,h2oexp,conexp,co2exp,n2oexp)
c**********************************************************************
c   Compute band3a exponentials for individual layers
c
c---- input parameters
c  number of grid intervals (m)
c  number of layers (np)
c  layer h2o amount for line absorption (dh2o)
c  layer h2o amount for continuum absorption (dcont)
c  layer co2 amount (dco2)
c  layer n2o amount (dn2o)
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c
c  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
c**********************************************************************
      implicit none
      integer m,np,i,k
 
c---- input parameters -----
 
      real dh2o(m,np),dcont(m,np),dn2o(m,np)
      real dco2(m,np),pa(m,np),dt(m,np)
 
c---- output parameters -----
 
      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
     *    ,n2oexp(m,np,4)
 
c---- temporary arrays -----
 
      real xx,xx1,xx2,xx3
 
c**********************************************************************
 
        do k=1,np
         do i=1,m
 
c-----Compute scaled h2o-line amount for Band 10 (Eq. 4.4 and Table 3).
 
           xx=dh2o(i,k)*(pa(i,k)/500.0)
     1           *(1.+(0.0149+6.20e-5*dt(i,k))*dt(i,k))
 
c-----six exponentials by powers of 8
 
           h2oexp(i,k,1)=exp(-xx*0.10624)
 
           xx=h2oexp(i,k,1)*h2oexp(i,k,1)
           xx=xx*xx
           h2oexp(i,k,2)=xx*xx
 
           xx=h2oexp(i,k,2)*h2oexp(i,k,2)
           xx=xx*xx
           h2oexp(i,k,3)=xx*xx
 
           xx=h2oexp(i,k,3)*h2oexp(i,k,3)
           xx=xx*xx
           h2oexp(i,k,4)=xx*xx
 
           xx=h2oexp(i,k,4)*h2oexp(i,k,4)
           xx=xx*xx
           h2oexp(i,k,5)=xx*xx
 
c-----one exponential of h2o continuum for sub-band 3a (Table 9).
 
           conexp(i,k,1)=exp(-dcont(i,k)*109.0)
 
c-----Scaled co2 amount for the Band 10 (Eq. 4.4, Tables 3 and 6).
 
           xx=dco2(i,k)*(pa(i,k)/300.0)**0.5
     1           *(1.+(0.0179+1.02e-4*dt(i,k))*dt(i,k))
 
c-----six exponentials by powers of 8
 
           co2exp(i,k,1,1)=exp(-xx*2.656e-5)
 
           xx=co2exp(i,k,1,1)*co2exp(i,k,1,1)
           xx=xx*xx
           co2exp(i,k,2,1)=xx*xx
 
           xx=co2exp(i,k,2,1)*co2exp(i,k,2,1)
           xx=xx*xx
           co2exp(i,k,3,1)=xx*xx
 
           xx=co2exp(i,k,3,1)*co2exp(i,k,3,1)
           xx=xx*xx
           co2exp(i,k,4,1)=xx*xx
 
           xx=co2exp(i,k,4,1)*co2exp(i,k,4,1)
           xx=xx*xx
           co2exp(i,k,5,1)=xx*xx
 
           xx=co2exp(i,k,5,1)*co2exp(i,k,5,1)
           xx=xx*xx
           co2exp(i,k,6,1)=xx*xx
 
c-----Compute the scaled n2o amount for Band 10 (Table 5).
 
           xx=dn2o(i,k)*(1.+(1.4476e-3+3.6656e-6*dt(i,k))*dt(i,k))
 
c-----Two exponentials by powers of 58
 
           n2oexp(i,k,1)=exp(-xx*0.25238)
 
           xx=n2oexp(i,k,1)*n2oexp(i,k,1)
           xx1=xx*xx
           xx1=xx1*xx1
           xx2=xx1*xx1
           xx3=xx2*xx2
           n2oexp(i,k,2)=xx*xx1*xx2*xx3
 
         enddo
        enddo
 
      return
      end
 
c**********************************************************************
      subroutine tablup(k2,m,np,nx,nh,dw,p,dt,s1,s2,s3,w1,p1,
     *                  dwe,dpe,coef1,coef2,coef3,tran)
c**********************************************************************
c   Compute water vapor, co2 and o3 transmittances between level
c   k1 and and level k2 for m soundings, using table look-up.
c
c   Calculations follow Eq. (4.16).
c
c---- input ---------------------
c
c  index for level (k2)
c  number of grid intervals (m)
c  number of atmospheric layers (np)
c  number of pressure intervals in the table (nx)
c  number of absorber amount intervals in the table (nh)
c  layer absorber amount (dw)
c  layer pressure in mb (p)
c  deviation of layer temperature from 250K (dt)
c  first value of absorber amount (log10) in the table (w1) 
c  first value of pressure (log10) in the table (p1) 
c  size of the interval of absorber amount (log10) in the table (dwe)
c  size of the interval of pressure (log10) in the table (dpe)
c  pre-computed coefficients (coef1, coef2, and coef3)
c
c---- updated ---------------------
c
c  column integrated absorber amount (s1)
c  absorber-weighted column pressure (s2)
c  absorber-weighted column temperature (s3)
c  transmittance (tran)
c
c  Note: Units of s1 are g/cm**2 for water vapor and
c       (cm-atm)stp for co2 and o3.
c   
c**********************************************************************
      implicit none
      integer k2,m,np,nx,nh,i
 
c---- input parameters -----
 
      real w1,p1,dwe,dpe
      real dw(m,np),p(m,np),dt(m,np)
      real coef1(nx,nh),coef2(nx,nh),coef3(nx,nh)
 
c---- update parameter -----
 
      real s1(m),s2(m),s3(m),tran(m,np+1)
 
c---- temporary variables -----
 
      real we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2,x1,x2,x3
      integer iw,ip
 
c-----Compute effective pressure (x2) and temperature (x3) following 
c     Eqs. (8.28) and (8.29)
 
      do i=1,m
 
        s1(i)=s1(i)+dw(i,k2-1)
        s2(i)=s2(i)+p(i,k2-1)*dw(i,k2-1)
        s3(i)=s3(i)+dt(i,k2-1)*dw(i,k2-1)
 
        x1=s1(i)
        x2=s2(i)/s1(i)
        x3=s3(i)/s1(i)
 
c-----normalize we and pe
 
        we=(log10(x1)-w1)/dwe
        pe=(log10(x2)-p1)/dpe
 
c-----restrict the magnitudes of the normalized we and pe.
 
        we=min(we,float(nh-1))
        pe=min(pe,float(nx-1))
 
c-----assign iw and ip and compute the distance of we and pe 
c     from iw and ip.
 
        iw=int(we+1.0)
        iw=min(iw,nh-1)
        iw=max(iw, 2)
        fw=we-float(iw-1)
 
        ip=int(pe+1.0)
        ip=min(ip,nx-1)
        ip=max(ip, 1)
        fp=pe-float(ip-1)
 
c-----linear interpolation in pressure
 
        pa = coef1(ip,iw-1)*(1.-fp)+coef1(ip+1,iw-1)*fp
        pb = coef1(ip,  iw)*(1.-fp)+coef1(ip+1,  iw)*fp
        pc = coef1(ip,iw+1)*(1.-fp)+coef1(ip+1,iw+1)*fp
 
c-----quadratic interpolation in absorber amount for coef1
 
        ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)
 
c-----linear interpolation in absorber amount for coef2 and coef3
 
        ba = coef2(ip,  iw)*(1.-fp)+coef2(ip+1,  iw)*fp
        bb = coef2(ip,iw+1)*(1.-fp)+coef2(ip+1,iw+1)*fp
        t1 = ba*(1.-fw) + bb*fw
 
        ca = coef3(ip,  iw)*(1.-fp)+coef3(ip+1,  iw)*fp
        cb = coef3(ip,iw+1)*(1.-fp)+coef3(ip+1,iw+1)*fp
        t2 = ca*(1.-fw) + cb*fw
 
c-----update the total transmittance between levels k1 and k2
 
        tran(i,k2)= (ax + (t1+t2*x3) * x3)*tran(i,k2)
        tran(i,k2)=min(tran(i,k2),0.9999999)
        tran(i,k2)=max(tran(i,k2),0.0000001)
 
      enddo
 
      return
      end
 
c**********************************************************************
      subroutine h2okdis(ib,m,np,k,fkw,gkw,ne,h2oexp,conexp,
     *                   th2o,tcon,tran)
c**********************************************************************
c   compute water vapor transmittance between levels k1 and k2 for
c   m soundings, using the k-distribution method.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals (m)
c  number of levels (np)
c  current level (k)
c  planck-weighted k-distribution function due to
c    h2o line absorption (fkw)
c  planck-weighted k-distribution function due to
c    h2o continuum absorption (gkw)
c  number of terms used in each band to compute water vapor
c     continuum transmittance (ne)
c  exponentials for line absorption (h2oexp) 
c  exponentials for continuum absorption (conexp) 
c
c---- updated parameters
c  transmittance between levels k1 and k2 due to
c    water vapor line absorption (th2o)
c  transmittance between levels k1 and k2 due to
c    water vapor continuum absorption (tcon)
c  total transmittance (tran)
c
c**********************************************************************
      implicit none
 
c---- input parameters ------
 
      integer ib,m,np,k,ne
      real conexp(m,np,3),h2oexp(m,np,6)
      real  fkw(6,9),gkw(6,3)
 
c---- updated parameters -----
 
      real th2o(m,6),tcon(m,3),tran(m,np+1)
 
c---- temporary arrays -----
 
      real trnth2o
      integer i
 
c-----tco2 are the six exp factors between levels k1 and k2 
c     tran is the updated total transmittance between levels k1 and k2
 
c-----th2o is the 6 exp factors between levels k1 and k2 due to
c     h2o line absorption. 
 
c-----tcon is the 3 exp factors between levels k1 and k2 due to
c     h2o continuum absorption.
 
c-----trnth2o is the total transmittance between levels k1 and k2 due
c     to both line and continuum absorption.
 
 
 
c-----Compute th2o following Eq. (8.23).
 
         do i=1,m
           th2o(i,1) = th2o(i,1)*h2oexp(i,k,1)
           th2o(i,2) = th2o(i,2)*h2oexp(i,k,2)
           th2o(i,3) = th2o(i,3)*h2oexp(i,k,3)
           th2o(i,4) = th2o(i,4)*h2oexp(i,k,4)
           th2o(i,5) = th2o(i,5)*h2oexp(i,k,5)
           th2o(i,6) = th2o(i,6)*h2oexp(i,k,6)
         enddo
 
 
      if (ne.eq.0) then
 
c-----Compute trnh2o following Eq. (8.25). fkw is given in Table 4.
 
         do i=1,m
 
           trnth2o      =(fkw(1,ib)*th2o(i,1)
     *                  + fkw(2,ib)*th2o(i,2)
     *                  + fkw(3,ib)*th2o(i,3)
     *                  + fkw(4,ib)*th2o(i,4)
     *                  + fkw(5,ib)*th2o(i,5)
     *                  + fkw(6,ib)*th2o(i,6))
 
          tran(i,k+1)=tran(i,k+1)*trnth2o
 
         enddo
 
      elseif (ne.eq.1) then
 
c-----Compute trnh2o following Eqs. (8.25) and (4.27).
 
         do i=1,m
 
           tcon(i,1)= tcon(i,1)*conexp(i,k,1)
 
           trnth2o      =(fkw(1,ib)*th2o(i,1)
     *                  + fkw(2,ib)*th2o(i,2)
     *                  + fkw(3,ib)*th2o(i,3)
     *                  + fkw(4,ib)*th2o(i,4)
     *                  + fkw(5,ib)*th2o(i,5)
     *                  + fkw(6,ib)*th2o(i,6))*tcon(i,1)
 
          tran(i,k+1)=tran(i,k+1)*trnth2o
 
         enddo
 
      else
 
c-----For band 3. This band is divided into 3 subbands.
 
         do i=1,m
 
           tcon(i,1)= tcon(i,1)*conexp(i,k,1)
           tcon(i,2)= tcon(i,2)*conexp(i,k,2)
           tcon(i,3)= tcon(i,3)*conexp(i,k,3)
 
c-----Compute trnh2o following Eqs. (4.29) and (8.25).
 
           trnth2o      = (  gkw(1,1)*th2o(i,1)
     *                     + gkw(2,1)*th2o(i,2)
     *                     + gkw(3,1)*th2o(i,3)
     *                     + gkw(4,1)*th2o(i,4)
     *                     + gkw(5,1)*th2o(i,5)
     *                     + gkw(6,1)*th2o(i,6) ) * tcon(i,1)
     *                  + (  gkw(1,2)*th2o(i,1)
     *                     + gkw(2,2)*th2o(i,2)
     *                     + gkw(3,2)*th2o(i,3)
     *                     + gkw(4,2)*th2o(i,4)
     *                     + gkw(5,2)*th2o(i,5)
     *                     + gkw(6,2)*th2o(i,6) ) * tcon(i,2)
     *                  + (  gkw(1,3)*th2o(i,1)
     *                     + gkw(2,3)*th2o(i,2)
     *                     + gkw(3,3)*th2o(i,3)
     *                     + gkw(4,3)*th2o(i,4)
     *                     + gkw(5,3)*th2o(i,5)
     *                     + gkw(6,3)*th2o(i,6) ) * tcon(i,3)
 
          tran(i,k+1)=tran(i,k+1)*trnth2o
 
         enddo
 
      endif
 
      return
      end
 
c**********************************************************************
      subroutine co2kdis(m,np,k,co2exp,tco2,tran)
c**********************************************************************
c   compute co2 transmittances between levels k1 and k2 for
c    m soundings, using the k-distribution method with linear
c    pressure scaling.
c
c---- input parameters
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for co2 absorption (co2exp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2 absorption
c     for the various values of the absorption coefficient (tco2)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,np,i,k
 
c---- input parameters -----
 
      real co2exp(m,np,6,2)
 
c---- updated parameters -----
 
      real tco2(m,6,2),tran(m,np+1)
 
c---- temporary arrays -----
 
      real xc
 
c-----tco2 is the 6 exp factors between levels k1 and k2 computed
c     from Eqs. (8.23) and (8.25). Also see Eq. (4.30).
c     The k-distribution functions are given in Table 10.
 
         do i=1,m
 
c-----band-wings
 
           tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
           xc=   0.1395 *tco2(i,1,1)
 
           tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
           xc=xc+0.1407 *tco2(i,2,1)
 
           tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
           xc=xc+0.1549 *tco2(i,3,1)
 
           tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
           xc=xc+0.1357 *tco2(i,4,1)
 
           tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
           xc=xc+0.0182 *tco2(i,5,1)
 
           tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
           xc=xc+0.0220 *tco2(i,6,1)
 
c-----band-center region
 
           tco2(i,1,2)=tco2(i,1,2)*co2exp(i,k,1,2)
           xc=xc+0.0766 *tco2(i,1,2)
 
           tco2(i,2,2)=tco2(i,2,2)*co2exp(i,k,2,2)
           xc=xc+0.1372 *tco2(i,2,2)
 
           tco2(i,3,2)=tco2(i,3,2)*co2exp(i,k,3,2)
           xc=xc+0.1189 *tco2(i,3,2)
 
           tco2(i,4,2)=tco2(i,4,2)*co2exp(i,k,4,2)
           xc=xc+0.0335 *tco2(i,4,2)
 
           tco2(i,5,2)=tco2(i,5,2)*co2exp(i,k,5,2)
           xc=xc+0.0169 *tco2(i,5,2)
 
           tco2(i,6,2)=tco2(i,6,2)*co2exp(i,k,6,2)
           xc=xc+0.0059 *tco2(i,6,2)
 
           tran(i,k+1)=tran(i,k+1)*xc
 
         enddo
 
      return
      end
 
c**********************************************************************
      subroutine n2okdis(ib,m,np,k,n2oexp,tn2o,tran)
c**********************************************************************
c   compute n2o transmittances between levels k1 and k2 for
c    m soundings, using the k-distribution method with linear
c    pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for n2o absorption (n2oexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to n2o absorption
c     for the various values of the absorption coefficient (tn2o)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters -----
 
      real n2oexp(m,np,4)
 
c---- updated parameters -----
 
      real tn2o(m,4),tran(m,np+1)
 
c---- temporary arrays -----
 
      real xc
 
c-----tn2o is computed from Eq. (8.23). 
c     xc is the total n2o transmittance computed from (8.25)
c     The k-distribution functions are given in Table 5.
 
        do i=1,m
 
c-----band 6
 
          if (ib.eq.6) then
 
           tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
           xc=   0.940414*tn2o(i,1)
 
           tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
           xc=xc+0.059586*tn2o(i,2)
 
c-----band 7
 
          else
 
           tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
           xc=   0.561961*tn2o(i,1)
 
           tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
           xc=xc+0.138707*tn2o(i,2)
 
           tn2o(i,3)=tn2o(i,3)*n2oexp(i,k,3)
           xc=xc+0.240670*tn2o(i,3)
 
           tn2o(i,4)=tn2o(i,4)*n2oexp(i,k,4)
           xc=xc+0.058662*tn2o(i,4)
 
          endif
 
           tran(i,k+1)=tran(i,k+1)*xc
 
        enddo
 
      return
      end
 
c**********************************************************************
      subroutine ch4kdis(ib,m,np,k,ch4exp,tch4,tran)
c**********************************************************************
c   compute ch4 transmittances between levels k1 and k2 for
c    m soundings, using the k-distribution method with
c    linear pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for ch4 absorption (ch4exp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to ch4 absorption
c     for the various values of the absorption coefficient (tch4)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters -----
 
      real ch4exp(m,np,4)
 
c---- updated parameters -----
 
      real tch4(m,4),tran(m,np+1)
 
c---- temporary arrays -----
 
      real xc
 
c-----tch4 is computed from Eq. (8.23). 
c     xc is the total ch4 transmittance computed from (8.25)
c     The k-distribution functions are given in Table 5.
 
        do i=1,m
 
c-----band 6
 
          if (ib.eq.6) then
 
           tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
           xc= tch4(i,1)
 
c-----band 7
 
          else
 
           tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
           xc=   0.610650*tch4(i,1)
 
           tch4(i,2)=tch4(i,2)*ch4exp(i,k,2)
           xc=xc+0.280212*tch4(i,2)
 
           tch4(i,3)=tch4(i,3)*ch4exp(i,k,3)
           xc=xc+0.107349*tch4(i,3)
 
           tch4(i,4)=tch4(i,4)*ch4exp(i,k,4)
           xc=xc+0.001789*tch4(i,4)
 
          endif
 
           tran(i,k+1)=tran(i,k+1)*xc
 
        enddo
 
      return
      end
 
c**********************************************************************
      subroutine comkdis(ib,m,np,k,comexp,tcom,tran)
c**********************************************************************
c  compute co2-minor transmittances between levels k1 and k2
c   for m soundings, using the k-distribution method
c   with linear pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for co2-minor absorption (comexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2-minor absorption
c     for the various values of the absorption coefficient (tcom)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,np,i,k
 
c---- input parameters -----
 
      real comexp(m,np,6)
 
c---- updated parameters -----
 
      real tcom(m,6),tran(m,np+1)
 
c---- temporary arrays -----
 
      real xc
 
c-----tcom is computed from Eq. (8.23). 
c     xc is the total co2 transmittance computed from (8.25)
c     The k-distribution functions are given in Table 6.
 
         do i=1,m
 
c-----band 4
 
           if (ib.eq.4) then
 
            tcom(i,1)=tcom(i,1)*comexp(i,k,1)
            xc=   0.12159*tcom(i,1)
            tcom(i,2)=tcom(i,2)*comexp(i,k,2)
            xc=xc+0.24359*tcom(i,2)
            tcom(i,3)=tcom(i,3)*comexp(i,k,3)
            xc=xc+0.24981*tcom(i,3)
            tcom(i,4)=tcom(i,4)*comexp(i,k,4)
            xc=xc+0.26427*tcom(i,4)
            tcom(i,5)=tcom(i,5)*comexp(i,k,5)
            xc=xc+0.07807*tcom(i,5)
            tcom(i,6)=tcom(i,6)*comexp(i,k,6)
            xc=xc+0.04267*tcom(i,6)
 
c-----band 5
 
           else
 
            tcom(i,1)=tcom(i,1)*comexp(i,k,1)
            xc=   0.06869*tcom(i,1)
            tcom(i,2)=tcom(i,2)*comexp(i,k,2)
            xc=xc+0.14795*tcom(i,2)
            tcom(i,3)=tcom(i,3)*comexp(i,k,3)
            xc=xc+   0.19512*tcom(i,3)
            tcom(i,4)=tcom(i,4)*comexp(i,k,4)
            xc=xc+   0.33446*tcom(i,4)
            tcom(i,5)=tcom(i,5)*comexp(i,k,5)
            xc=xc+   0.17199*tcom(i,5)
            tcom(i,6)=tcom(i,6)*comexp(i,k,6)
            xc=xc+   0.08179*tcom(i,6)
           endif
 
            tran(i,k+1)=tran(i,k+1)*xc
 
         enddo
 
      return
      end
 
c**********************************************************************
      subroutine cfckdis(m,np,k,cfcexp,tcfc,tran)
c**********************************************************************
c  compute cfc-(11,12,22) transmittances between levels k1 and k2
c   for m soundings, using the k-distribution method with
c   linear pressure scaling.
c
c---- input parameters
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for cfc absorption (cfcexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to cfc absorption
c     for the various values of the absorption coefficient (tcfc)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,np,i,k
 
c---- input parameters -----
 
      real cfcexp(m,np)
 
c---- updated parameters -----
 
      real tcfc(m),tran(m,np+1)
 
c-----tcfc is the exp factors between levels k1 and k2. 
 
         do i=1,m
 
            tcfc(i)=tcfc(i)*cfcexp(i,k)
            tran(i,k+1)=tran(i,k+1)*tcfc(i)
 
         enddo
 
      return
      end
 
c**********************************************************************
      subroutine b10kdis(m,np,k,h2oexp,conexp,co2exp,n2oexp
     *          ,th2o,tcon,tco2,tn2o,tran)
c**********************************************************************
c
c   compute h2o (line and continuum),co2,n2o transmittances between
c   levels k1 and k2 for m soundings, using the k-distribution
c   method with linear pressure scaling.
c
c---- input parameters
c   number of grid intervals (m)
c   number of levels (np)
c   current level (k)
c   exponentials for h2o line absorption (h2oexp)
c   exponentials for h2o continuum absorption (conexp)
c   exponentials for co2 absorption (co2exp)
c   exponentials for n2o absorption (n2oexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to h2o line absorption
c     for the various values of the absorption coefficient (th2o)
c   transmittance between levels k1 and k2 due to h2o continuum
c     absorption for the various values of the absorption
c     coefficient (tcon)
c   transmittance between levels k1 and k2 due to co2 absorption
c     for the various values of the absorption coefficient (tco2)
c   transmittance between levels k1 and k2 due to n2o absorption
c     for the various values of the absorption coefficient (tn2o)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,np,i,k
 
c---- input parameters -----
 
      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
     *    ,n2oexp(m,np,4)
 
c---- updated parameters -----
 
      real th2o(m,6),tcon(m,3),tco2(m,6,2),tn2o(m,4)
     *    ,tran(m,np+1)
 
c---- temporary arrays -----
 
      real xx
 
c-----For h2o line. The k-distribution functions are given in Table 4.
 
        do i=1,m
 
           th2o(i,1)=th2o(i,1)*h2oexp(i,k,1)
           xx=   0.3153*th2o(i,1)
 
           th2o(i,2)=th2o(i,2)*h2oexp(i,k,2)
           xx=xx+0.4604*th2o(i,2)
 
           th2o(i,3)=th2o(i,3)*h2oexp(i,k,3)
           xx=xx+0.1326*th2o(i,3)
 
           th2o(i,4)=th2o(i,4)*h2oexp(i,k,4)
           xx=xx+0.0798*th2o(i,4)
 
           th2o(i,5)=th2o(i,5)*h2oexp(i,k,5)
           xx=xx+0.0119*th2o(i,5)
 
           tran(i,k+1)=xx
 
        enddo
 
c-----For h2o continuum. Note that conexp(i,k,3) is for subband 3a.
 
        do i=1,m
 
           tcon(i,1)=tcon(i,1)*conexp(i,k,1)
           tran(i,k+1)=tran(i,k+1)*tcon(i,1)
 
        enddo
 
c-----For co2 (Table 6)
 
        do i=1,m
 
           tco2(i,1,1)=tco2(i,1,1)*co2exp(i,k,1,1)
           xx=    0.2673*tco2(i,1,1)
 
           tco2(i,2,1)=tco2(i,2,1)*co2exp(i,k,2,1)
           xx=xx+ 0.2201*tco2(i,2,1)
 
           tco2(i,3,1)=tco2(i,3,1)*co2exp(i,k,3,1)
           xx=xx+ 0.2106*tco2(i,3,1)
 
           tco2(i,4,1)=tco2(i,4,1)*co2exp(i,k,4,1)
           xx=xx+ 0.2409*tco2(i,4,1)
 
           tco2(i,5,1)=tco2(i,5,1)*co2exp(i,k,5,1)
           xx=xx+ 0.0196*tco2(i,5,1)
 
           tco2(i,6,1)=tco2(i,6,1)*co2exp(i,k,6,1)
           xx=xx+ 0.0415*tco2(i,6,1)
 
           tran(i,k+1)=tran(i,k+1)*xx
 
        enddo
 
c-----For n2o (Table 5)
 
        do i=1,m
 
           tn2o(i,1)=tn2o(i,1)*n2oexp(i,k,1)
           xx=   0.970831*tn2o(i,1)
 
           tn2o(i,2)=tn2o(i,2)*n2oexp(i,k,2)
           xx=xx+0.029169*tn2o(i,2)
           tran(i,k+1)=tran(i,k+1)*(xx-1.0)
 
        enddo
 
      return
      end
 
c***********************************************************************
      subroutine cldovlp (m,np,k2,ict,icb,it,im,ib,itx,imx,ibx,
     *               cldhi,cldmd,cldlw,fcld,tcldlyr,fclr)
 
c***********************************************************************
c     compute the fractional clear line-of-sight between levels k1
c     and k2 following Eqs.(6.18)-(6.21).
c
c input parameters
c
c  m:       number of soundings
c  np:      number of layers
c  k2:      index for the level
c  ict:     the level separating high and middle clouds
c  icb:     the level separating middle and low clouds
c  it:      number of cloudy layers in the high-cloud group
c  im:      number of cloudy layers in the middle-cloud group
c  ib:      number of cloudy layers in the low-cloud group
c  fcld:    fractional cloud cover of a layer
c  tcldlyr: transmittance of a cloud layer
c  
c output parameter
c
c  fclr:    clear line-of-sight between levels k1 and k2
c***********************************************************************
 
      implicit none
      integer m,np,k2,ict,icb
      integer i,j,k,ii,it(m),im(m),ib(m),itx(m,np),imx(m,np),ibx(m,np)
      real cldhi(m),cldmd(m),cldlw(m)
      real fcld(m,np),tcldlyr(m,np),fclr(m,np+1)
 
c***********************************************************************
       do i=1,m
 
c-----For high clouds
c     "it" is the number of high-cloud layers
 
        if (k2.le.ict) then
         if(fcld(i,k2-1).gt.0.001) then
 
          it(i)=it(i)+1
          ii=it(i)
          itx(i,ii)=k2-1
 
         if (ii .eq. 1) go to 11
 
c-----Rearrange the order of cloud layers with increasing cloud amount
 
         do k=1,ii-1
           j=itx(i,k)
          if(fcld(i,j).gt.fcld(i,k2-1)) then
           do j=ii-1,k,-1
            itx(i,j+1)=itx(i,j)
           enddo
            itx(i,k)=k2-1
            go to 11
          endif
         enddo
 
   11   continue
 
c-----compute equivalent black-body high cloud amount
 
           cldhi(i)=0.0
          do k=1,ii
           j=itx(i,k)
           cldhi(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldhi(i))
          enddo
 
        endif
       endif
 
c-----For middle clouds
c     "im" is the number of middle-cloud layers
 
       if (k2.gt.ict .and. k2.le.icb) then
        if(fcld(i,k2-1).gt.0.001) then
 
         im(i)=im(i)+1
         ii=im(i)
         imx(i,ii)=k2-1
 
        if (ii .eq. 1) go to 21
 
c-----Rearrange the order of cloud layers with increasing cloud amount
 
         do k=1,ii-1
            j=imx(i,k)
           if(fcld(i,j).gt.fcld(i,k2-1)) then
            do j=ii-1,k,-1
             imx(i,j+1)=imx(i,j)
            enddo
             imx(i,k)=k2-1
             go to 21
           endif
          enddo
 
   21   continue
 
c-----compute equivalent black-body middle cloud amount
 
           cldmd(i)=0.0
          do k=1,ii
           j=imx(i,k)
           cldmd(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldmd(i))
          enddo
 
        endif
       endif
 
c-----For low clouds
c     "ib" is the number of low-cloud layers
 
       if (k2.gt.icb) then
        if(fcld(i,k2-1).gt.0.001) then
 
         ib(i)=ib(i)+1
         ii=ib(i)
         ibx(i,ii)=k2-1
 
        if (ii .eq. 1) go to 31
 
c-----Rearrange the order of cloud layers with increasing cloud amount
 
         do k=1,ii-1
          j=ibx(i,k)
           if(fcld(i,j).gt.fcld(i,k2-1)) then
            do j=ii-1,k,-1
             ibx(i,j+1)=ibx(i,j)
            enddo
             ibx(i,k)=k2-1
             go to 31
           endif
          enddo
 
   31    continue
 
c-----compute equivalent black-body low cloud amount
 
           cldlw(i)=0.0
          do k=1,ii
           j=ibx(i,k)
           cldlw(i)=fcld(i,j)-tcldlyr(i,j)*(fcld(i,j)-cldlw(i))
          enddo
 
        endif
       endif
 
c-----fclr is the equivalent clear fraction between levels k1 and k2
c     assuming the three cloud groups are randomly overlapped.
c     It follows Eqs. (6.20) and (6.21).
 
        fclr(i,k2)=(1.0-cldhi(i))*(1.0-cldmd(i))*(1.0-cldlw(i))   
 
      enddo
 
      return
      end
 
c***********************************************************************
      subroutine sfcflux (ibn,m,ns,fs,tg,eg,tv,ev,rv,vege,
     *                    bs,dbs,rflxs)
c***********************************************************************
c Compute emission and reflection by an homogeneous/inhomogeneous 
c  surface with vegetation cover.
c
c-----Input parameters
c  index for the spectral band (ibn)
c  number of grid box (m)
c  number of sub-grid box (ns)
c  fractional cover of sub-grid box (fs)
c  sub-grid ground temperature (tg)
c  sub-grid ground emissivity (eg)
c  sub-grid vegetation temperature (tv)
c  sub-grid vegetation emissivity (ev)
c  sub-grid vegetation reflectivity (rv)
c  if there is vegetation cover, vege=.true.
c
c-----Output parameters
c  Emission by the surface (ground+vegetation) (bs)
c  Derivative of bs rwt temperature (dbs)
c  Reflection by the surface (rflxs)
 
c**********************************************************************
       implicit none
 
c---- input parameters -----
       integer ibn,m,ns
       real fs(m,ns),tg(m,ns),eg(m,ns,10)
       real tv(m,ns),ev(m,ns,10),rv(m,ns,10)
       logical vege
 
c---- output parameters -----
       real bs(m),dbs(m),rflxs(m)
 
c---- temporary arrays -----
 
       integer i,j
       real bg(m),dbg(m),bv(m),dbv(m),tx(m),ty(m),xx,yy,zz
 
c*********************************************************
       write(0,*) 'tg =', tg
        if (ns.eq.1) then
 
         if (.not.vege) then
 
c-----for homogeneous surface without vegetation
c     following Eqs. (9.4), (9.5), and (3.13)
 
          do i=1,m
           tx(i)=tg(i,1)
          enddo
 
           call planck(ibn,m,tx,bg)
           call plancd(ibn,m,tx,dbg)
 
          do i=1,m
           bs(i) =eg(i,1,ibn)*bg(i)
           dbs(i)=eg(i,1,ibn)*dbg(i)
           rflxs(i)=1.0-eg(i,1,ibn)
          enddo
 
         else
 
c-----With vegetation, following Eqs. (9.1), (9.3), and (9.13)
 
          do i=1,m
           tx(i)=tg(i,1)
           ty(i)=tv(i,1)
          enddo
 
           call planck(ibn,m,tx,bg)
           call planck(ibn,m,ty,bv)
           call plancd(ibn,m,tx,dbg)
           call plancd(ibn,m,ty,dbv)
 
          do i=1,m
           xx=ev(i,1,ibn)*bv(i)
           yy=1.0-ev(i,1,ibn)-rv(i,1,ibn)
           zz=1.0-eg(i,1,ibn)
           bs(i)=yy*(eg(i,1,ibn)*bg(i)+zz*xx)+xx
 
           xx=ev(i,1,ibn)*dbv(i)
           dbs(i)=yy*(eg(i,1,ibn)*dbg(i)+zz*xx)+xx
 
           rflxs(i)=rv(i,1,ibn)+zz*yy*yy/(1.0-rv(i,1,ibn)*zz)
          enddo
 
         endif
 
        else
 
c-----for nonhomogeneous surface
 
          do i=1,m
           bs(i)=0.0
           dbs(i)=0.0
           rflxs(i)=0.0
          enddo
 
         if(.not.vege) then
 
c-----No vegetation, following Eqs. (9.9), (9.10), and (9.13)
 
          do j=1,ns
 
           do i=1,m
            tx(i)=tg(i,j)
           enddo
 
            call planck(ibn,m,tx,bg)
            call plancd(ibn,m,tx,dbg)
 
           do i=1,m
            bs(i)=bs(i)+fs(i,j)*eg(i,j,ibn)*bg(i)
            dbs(i)=dbs(i)+fs(i,j)*eg(i,j,ibn)*dbg(i)
            rflxs(i)=rflxs(i)+fs(i,j)*(1.0-eg(i,j,ibn))
           enddo
 
          enddo
 
         else
 
c-----With vegetation, following Eqs. (9.6), (9.7), and (9.13)
 
          do j=1,ns
           do i=1,m
            tx(i)=tg(i,j)
            ty(i)=tv(i,j)
           enddo
 
            call planck(ibn,m,tx,bg)
            call planck(ibn,m,ty,bv)
            call plancd(ibn,m,tx,dbg)
            call plancd(ibn,m,ty,dbv)
 
           do i=1,m
            xx=ev(i,j,ibn)*bv(i)
            yy=1.0-ev(i,j,ibn)-rv(i,j,ibn)
            zz=1.0-eg(i,j,ibn)
            bs(i)=bs(i)+fs(i,j)*(yy*(eg(i,j,ibn)*bg(i)+zz*xx)+xx)
 
            xx=ev(i,j,ibn)*dbv(i)
            dbs(i)=dbs(i)+fs(i,j)*(yy*(eg(i,j,ibn)*dbg(i)+zz*xx)+xx)
 
            rflxs(i)=rflxs(i)+fs(i,j)*(rv(i,j,ibn)+zz*yy*yy
     *             /(1.0-rv(i,j,ibn)*zz))
           enddo
          enddo
 
         endif
 
        endif
 
      return
      end
 
