;=============================================================
; Composite radar reflectivity from gSAM 3D file + CONUS map
;=============================================================
begin
  ;--------------------------
  ; User input
  ;--------------------------
  fpath = "../../OUT_3D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-00-00-00_0000008640.3D_atm.nc"
  it    = 0          ; time index
  dmin  = -30.0      ; dBZ floor (display)
  dmax  = 75.0       ; dBZ cap (display)
  tinyLWC = 1.0e-8   ; kg m^-3 guard
  Zfloor  = 1.0e-6   ; mm^6 m^-3 to avoid log(0)

  ; Microphysical/MP distribution parameters (match your Fortran "micro_params")
  ; Adjust if your scheme uses different values.
  N0r  = 8.0e6       ; m^-4  (rain intercept)
  N0s  = 3.0e6       ; m^-4  (snow intercept)
  N0g  = 4.0e6       ; m^-4  (graupel/hail intercept)
  rhor = 1000.0      ; kg m^-3 (rain density)
  rhos = 100.0       ; kg m^-3 (snow density, dry)
  rhog = 400.0       ; kg m^-3 (graupel density)

  ; Electromagnetic (S-band) and melting-layer blend (as in your Fortran)
  K2w   = 0.93       ; |K|^2, liquid water
  K2i   = 0.197      ; |K|^2, dry ice/snow/graupel
  Tm0   = 273.15 - 2.0   ; start blend at -2 C
  Tm1   = 273.15         ; fully wet at 0 C
  pi    = 3.141592653589793
  mm6m3_from_m = 1.0e18  ; convert m^6 to mm^6

  ;--------------------------
  ; Read data
  ;--------------------------
  f    = addfile(fpath,"r")
  lon  = f->lon               ; (lon)
  lat  = f->lat               ; (lat)
  rho  = f->rho               ; (z) kg m^-3
  TABS = f->TABS(it,:,:,:)    ; (z,lat,lon) or (time,z,lat,lon) -> here sliced to (z,lat,lon)
  QR   = f->QR(it,:,:,:)      ; g/kg
  QS   = f->QS(it,:,:,:)      ; g/kg
  QG   = f->QG(it,:,:,:)      ; g/kg

  ; Ensure coordinate order is (z,lat,lon)
  ; (gSAM file already stores as (time,z,lat,lon), so above slicing is fine.)

  ; Convert mixing ratios to kg/kg
  QR = QR*0.001
  QS = QS*0.001
  QG = QG*0.001

  ; If longitudes are 0..360, shift to -180..180 for CONUS plotting
  lon_plot = lon
  lon_plot = where(lon_plot.gt.180.,lon_plot-360.,lon_plot)

  print("lon_plot = "+min(lon_plot)+" "+max(lon_plot))
  print("lat = "+min(lat)+" "+max(lat))

  nz = dimsizes(TABS(:,0,0))
  ny = dimsizes(TABS(0,:,0))
  nx = dimsizes(TABS(0,0,:))
  print("nx ny nz: "+nx+" "+ny+" "+nz)

  ; Allocate outputs
  Zcomp = new((/ny,nx/), "float")
  Zcomp = dmin

  ;--------------------------
  ; Level-by-level computation (Fortran-aligned)
  ; lam = [(pi * rho_p * N0) / LWC]^(1/4)
  ; Z_species = 720 * N0 / lam^7 * 1e18    (mm^6 m^-3)
  ; Wet-ice dielectric blend across [Tm0, Tm1]
  ; Ze_lin = ( Zr*K2w + (Zs+Zg)*K2ice_eff ) / K2w
  ; dBZ = 10*log10(max(Ze_lin, Zfloor)), then clamp to [dmin,dmax]
  ;--------------------------
  do k = 0, nz-1
    Tk  = TABS(k,:,:)
    rk  = rho(k)             ; scalar

    lwc_r = rk * QR(k,:,:)   ; kg m^-3
    lwc_s = rk * QS(k,:,:)
    lwc_g = rk * QG(k,:,:)

    ; Avoid division by 0: set lam only where LWC > tiny
    valid_r = where(lwc_r.gt.tinyLWC, 1, 0)
    valid_s = where(lwc_s.gt.tinyLWC, 1, 0)
    valid_g = where(lwc_g.gt.tinyLWC, 1, 0)

    lam_r = new((/ny,nx/),"double")
    lam_s = new((/ny,nx/),"double")
    lam_g = new((/ny,nx/),"double")

    lam_r = ((pi*rhor*N0r) / (lwc_r+tinyLWC))^(0.25)
    lam_s = ((pi*rhos*N0s) / (lwc_s+tinyLWC))^(0.25)
    lam_g = ((pi*rhog*N0g) / (lwc_g+tinyLWC))^(0.25)

    ; Species reflectivities (linear)
    Zr = 720.0*N0r / (lam_r^7) * mm6m3_from_m
    Zs = 720.0*N0s / (lam_s^7) * mm6m3_from_m
    Zg = 720.0*N0g / (lam_g^7) * mm6m3_from_m

    ; Zero them where invalid
    Zr = Zr*valid_r
    Zs = Zs*valid_s
    Zg = Zg*valid_g

    ; Dielectric factor for (snow+graupel): linear blend across melting band
    K2ice_eff = new((/ny,nx/),"double")
    K2ice_eff = K2i
    K2ice_eff = where(Tk.ge.Tm1, K2w, K2ice_eff)
    blendMask = where((Tk.gt.Tm0).and.(Tk.lt.Tm1), 1, 0)
    frac = (Tk - Tm0)/(Tm1 - Tm0)
    K2ice_eff = where(blendMask.eq.1, K2i + (K2w - K2i)*frac, K2ice_eff)

    Ze_lin = ( Zr*K2w + (Zs+Zg)*K2ice_eff )/K2w
    Ze_lin = where(Ze_lin.gt.Zfloor,Ze_lin,Zfloor)

    dBZ = tofloat(10.0*log10(Ze_lin))
    dBZ = where(dBZ.gt.dmin,dBZ,dmin)
    dBZ = where(dBZ.lt.dmax,dBZ,dmax)

    Zcomp       = where(Zcomp.gt.dBZ,Zcomp,dBZ)
  end do

  ; Attach simple coords to the 2-D field
  Zcomp!0   = "lat"
  Zcomp!1   = "lon"
  Zcomp&lat = lat
  Zcomp&lon = lon_plot
  Zcomp@long_name = "Composite Reflectivity"
  Zcomp@units     = "dBZ"

  ;--------------------------
  ; Plot
  ;--------------------------
;  wks = gsn_open_wks("png","comp_reflectivity_CONUS")  ; output: comp_reflectivity_CONUS.png
  wks = gsn_open_wks("ps","ps")  ; output: comp_reflectivity_CONUS.png
  res                       = True
  res@cnFillOn              = True
  res@cnFillMode            = "RasterFill"
  res@cnLinesOn             = False
;  res@cnLevelSelectionMode  = "ManualLevels"
;  res@cnMinLevelValF        = -10.0
;  res@cnMaxLevelValF        = 60.0
;  res@cnLevelSpacingF       = 2.0
  res@lbTitleOn             = True
  res@lbTitleString         = "Composite Reflectivity (dBZ)"
  res@pmTickMarkDisplayMode = "Always"

  ; Colormap (any perceptual map works; change if you prefer a radar palette)
;  res@cnFillPalette         = "MPL_viridis"

  ; Map settings
  res@mpProjection                = "CylindricalEquidistant"
  res@mpDataBaseVersion          = "MediumRes"
  res@mpOutlineBoundarySets      = "GeophysicalAndUSStates"
  res@mpUSStateLineThicknessF    = 2.0
  res@mpGridAndLimbOn            = False
  res@gsnAddCyclic               = False
  res@mpMinLatF                  = 20.0      ; crop to CONUS box (adjust if needed)
  res@mpMaxLatF                  = 55.0
  res@mpMinLonF                  = -130.0
  res@mpMaxLonF                  = -60.0

  ; Tell NCL lat/lon are 1D coords of a rectilinear grid
  res@sfXArray = lon_plot
  res@sfYArray = lat

  colors = (/ (/255,255,255/), (/77,230,233/), \
              (/43,150,235/), (/0,35,232/), (/91,247,66/), \
              (/68,186,49/), (/47,128,34/), \
              (/255,252,72/),(/224,182,54/),(/245,136,45/), \
              (/239,47,29/),(/187,57,26/),(/134,27,16/), \
              (/237,60,245/),(/133,59,190/)/) / 255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,70.,75./)
  res@cnFillPalette        = colors


  print("Zcomp: "+min(Zcomp)+"  "+max(Zcomp))
  plot = gsn_csm_contour_map(wks, Zcomp, res)

end

