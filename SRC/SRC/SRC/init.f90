! Initialize some arrays from vars module:

subroutine init()

  use vars
  use tracers

  implicit none

  nstep = 0
  fzero = 0.
 
  ttend = 0.
  qtend = 0.
  wsub = 0.
  unudge = 0.
  vnudge = 0.
  tnudge = 0.
  qnudge = 0.
  qlsvadv = 0.
  tlsvadv = 0.
  ulsvadv = 0.
  vlsvadv = 0.
  qstor = 0.
  tstor = 0.
  ustor = 0.
  vstor = 0.
  raf = 0.
  tg0(:)=0.
  qg0(:)=0.
  ug0(:)=0.
  vg0(:)=0.

  radlwup = 0.
  radlwdn = 0.
  radswup = 0.
  radswdn = 0.
  radqrlw = 0.
  radqrsw = 0.
  radqrclw = 0.
  radqrcsw = 0.
 
  tlat = 0. 
  tdiss = 0. 
  tpdiss = 0. 
  tbuoy = 0. 
  tlatqi = 0.
  tlatqc = 0.
  tadv = 0.
  tdiff = 0.
  qifall = 0.
  qcfall = 0.
  qpfall = 0.

  twle = 0.
  twsb = 0.
  uwle = 0.
  uwsb = 0.
  vwle = 0.
  vwsb = 0.

  trwle = 0.
  trwsb = 0.
  tradv = 0.
  trdiff = 0.
  trphys = 0.

  t2leprec = 0.
  q2leprec = 0.
  twleprec = 0.
  qwleprec = 0.

  tkelediff = 0.
  t2lediff = 0.
  q2lediff = 0.
  twlediff = 0.
  qwlediff = 0.
  momlediff = 0.

  tkelebuoy = 0.
  twlebuoy = 0.
  qwlebuoy = 0.
  momlebuoy = 0.

  tkelediss = 0.
  t2lediss = 0.
  q2lediss = 0.


  gamt0 = 0.
  gamq0 = 0.

  precsfc = 0.
  precsfcsnow = 0.
  precflux = 0.
  prectot = 0.
  wvp = 0.
  wvpobs = 0.
  cwp = 0.
  iwp = 0.
  rwp = 0.
  gwp = 0.
  swp = 0.

  fluxbu = 0.
  fluxbv = 0.
  fluxbt = 0.
  fluxbq = 0.
  fluxtu = 0.
  fluxtv = 0.
  fluxtt = 0.
  fluxtq = 0.

  precinst = 0.
  shf_ocn(:,:) = 0.
  lhf_ocn(:,:) = 0.
  taux_ocn(:,:) = 0.
  tauy_ocn(:,:) = 0.
  prec_ocn(:,:) = 0.
  roff_ocn(:,:) = 0.
  lw_ocn(:,:) = 0.
  sw_ocn(:,:) = 0.
  u_ocn(:,:) = 0.
  v_ocn(:,:) = 0.

  shf_ocean(:,:) = 0.
  lhf_ocean(:,:) = 0.
  shf_land(:,:) = 0.
  lhf_land(:,:) = 0.
  shf_all(:,:) = 0.
  lhf_all(:,:) = 0.
  evp_all(:,:) = 0.
  shf_top(:,:) = 0.
  lhf_top(:,:) = 0.

  s_acldisccp = 0.
  s_acldlisccp = 0.
  s_acldmisccp = 0.
  s_acldhisccp = 0.
  s_ptopisccp = 0.
  s_acldmodis = 0.
  s_acldlmodis = 0.
  s_acldmmodis = 0.
  s_acldhmodis = 0.
  s_ptopmodis = 0.
  s_acldmisr = 0.
  s_ztopmisr = 0.
  s_relmodis = 0.
  s_reimodis = 0.
  s_lwpmodis = 0.
  s_iwpmodis = 0.
  s_tbisccp = 0.
  s_tbclrisccp = 0.
  s_acldliqmodis = 0.
  s_acldicemodis = 0.
  s_cldtauisccp = 0.
  s_cldtaumodis = 0.
  s_cldtaulmodis = 0.
  s_cldtauimodis = 0.
  s_cldalbisccp = 0.

  dtstat = 0.

  dudtd(:,:,:) = 0.
  dvdtd(:,:,:) = 0.
  dwdtd(:,:,:) = 0.

end subroutine init



