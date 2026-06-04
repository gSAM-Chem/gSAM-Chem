load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; SAM datset
;filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nonudge_4608_0000057120"
filenameroot = "/glade/campaign/univ/unsb0021/SAM7.3_SOCRATES/OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nonudge_SAM1MOM_gmg6_4608_0000057600"
;-----------------------------------------------------------------------------------------
; ERA5 datset

;date = 2018017023
date = 2018022100

;-------------------------------------------------------------------
; plot resolution

lowres_nx = 128
lowres_ny = 65

;-----------------------------------------------------------------------------------------
; find corresponding ERA5 file names:

year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
day = floattointeger((date-year*1000000-month*10000)/100.)
hour = date-year*1000000-month*10000-day*100
date_obs = (year*1000000+month*10000+day*100)+"_"+(year*1000000+month*10000+day*100+23)
month_obs = ""+(year*100+month)
print(""+date_obs)
dirpath = "/glade/collections/rda/data/ds633.0/e5.oper.an.pl/"+month_obs+"/"
f_u = addfile(dirpath+"e5.oper.an.pl.128_131_u.ll025uv."+date_obs+".nc","r")
f_v = addfile(dirpath+"e5.oper.an.pl.128_132_v.ll025uv."+date_obs+".nc","r")
f_w = addfile(dirpath+"e5.oper.an.pl.128_135_w.ll025sc."+date_obs+".nc","r")
f_z = addfile(dirpath+"e5.oper.an.pl.128_129_z.ll025sc."+date_obs+".nc","r")
f_t = addfile(dirpath+"e5.oper.an.pl.128_130_t.ll025sc."+date_obs+".nc","r")
f_q = addfile(dirpath+"e5.oper.an.pl.128_133_q.ll025sc."+date_obs+".nc","r")
f_clwc = addfile(dirpath+"e5.oper.an.pl.128_246_clwc.ll025sc."+date_obs+".nc","r")
f_ciwc = addfile(dirpath+"e5.oper.an.pl.128_247_ciwc.ll025sc."+date_obs+".nc","r")
times = calendar_decode2(f_u->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(""+dates)
ntime = dimsizes(f_u->time)
itime = -1
do i=0,ntime-1
 if(dates(i).eq.date) then
  if(itime.eq.-1) then
     itime = i
  end if
 end if
end do
if(itime.eq.-1) then
  print("could not find date "+date)
  exit
end if
time = day_of_year(times(itime,0),times(itime,1),times(itime,2))+times(itime,3)/24.
print("  time="+time)
mmm = itime
print("itime="+itime)
;-----------------------------------------------------------------------------------------
; prepare grid for comparison (ERA5 native grid):

nlati = dimsizes(f_u->latitude)
nloni = dimsizes(f_u->longitude)
nlevi = dimsizes(f_u->level)
grid = nloni+"x"+nlati
print("data grid:"+grid+"x"+nlevi)

; determine pressure levels at z:

presi = doubletofloat(f_u->level(:))
zi = f_z->Z(mmm,nlevi-1:0:1,nlati-1:0:1,:)
lati = doubletofloat(zi&latitude)
lati@units = zi&latitude@units
loni = doubletofloat(zi&longitude)
loni@units = zi&longitude@units
wgt = cos(lati*acos(-1.)/180.)
zim = dim_avg_Wrap(dim_avg_wgt_n_Wrap(zi,wgt,0,1))/9.80665
print("zim: "+zim+"  presi: "+presi(nlevi-1:0:1))

;-----------------------------------------------------------------------------------------
; prepare initial settings for plots:

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnLinesOn       = False
  res@cnFillPalette        = "BlueYellowRed"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse        = True

;controlling label bar
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = True         ; labelbar box lines.
  res@lbTitlePosition     = "Bottom"
  res@lbTitleFontHeightF = 0.02
  res@lbLabelFontHeightF = 0.015
  res@pmLabelBarHeightF = 0.1
  res@lbBoxMinorExtentF = 0.2
  res@lbBoxSeparatorLinesOn=False
  res@gsnStringFontHeightF = 0.02

  res@gsnDraw         = False               ; Do not draw plot
  res@gsnFrame        = False               ; Do not advance frame

plot = new(3,graphic)

;-----------------------------------------------------------------------------------------

; ------ U
f=addfile(filenameroot+"_U.nc","r")
pres = f->p
nlat = dimsizes(f->lat)
fld = dim_avg_Wrap(f->U(0,:,:,:))
fld&z = pres
fld&lat@units = "degrees_north"
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
f=addfile(filenameroot+"_TERRA.nc","r")
terr = dim_avg_Wrap(f->TERRA(0,:,:,:))
terr&z = pres
terr&lat@units = "degrees_north"
terr_lowres = linint1_Wrap(yi,terr,True,yo,0)
print("terr_lowre:"+min(terr_lowres)+" "+max(terr_lowres))
fld_lowres = linint1_Wrap(yi,fld,True,yo,0)
print("u: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs=dim_avg(f_u->U(mmm,nlevi-1:0:1,nlati-1:0:1,:))
fld_obs_lowres = linint1(lati,fld_obs,True,yo,0)
fld_obs_lowres!0 = "lev"
fld_obs_lowres&lev = presi(nlevi-1:0:1)
fld_obs_lowres!1 = "lat"
fld_obs_lowres&lat = yo
print("u_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld0 = int2p_n_Wrap(presi(nlevi-1:0:1),fld_obs_lowres,pres,-1,0)
print("u_obs0:"+min(fld0)+"  "+max(fld0))
fld_lowres@units=" SAM"
res@cnLevelSpacingF      = 5.                  ; contour spacing
res@cnMinLevelValF       = -50.                ; min level
res@cnMaxLevelValF       = 50.                 ; max level
plot(0)  = gsn_csm_contour(wks, fld_lowres, res )   ; place holder
fld0@units="ERA5"
plot(1)  = gsn_csm_contour(wks, fld0, res )   ; place holder
fld0 = fld_lowres-fld0
fld0@units="SAM-ERA5"
;symMinMaxPlt (fld0,20,False,res)
res@cnLevelSpacingF      = 1.                  ; contour spacing
res@cnMinLevelValF       = -5.                ; min level
res@cnMaxLevelValF       = 5.                 ; max level
plot(2)  = gsn_csm_contour(wks, fld0, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

; ------ V
f=addfile(filenameroot+"_V.nc","r")
pres = f->p
nlat = dimsizes(f->lat)
fld = dim_avg_Wrap(f->V(0,:,:,:))
fld&z = pres
fld&lat@units = "degrees_north"
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
fld_lowres = linint1_Wrap(yi,fld,True,yo,0)
print("v: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs=dim_avg(f_v->V(mmm,nlevi-1:0:1,nlati-1:0:1,:))
fld_obs_lowres = linint1(lati,fld_obs,True,yo,0)
fld_obs_lowres!0 = "lev"
fld_obs_lowres&lev = presi(nlevi-1:0:1)
fld_obs_lowres!1 = "lat"
fld_obs_lowres&lat = yo
print("v_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld0 = int2p_n_Wrap(presi(nlevi-1:0:1),fld_obs_lowres,pres,-1,0)
print("v_obs0:"+min(fld0)+"  "+max(fld0))
fld_lowres@units=" SAM"
res@cnLevelSpacingF      = 0.5                  ; contour spacing
res@cnMinLevelValF       = -5.5                ; min level
res@cnMaxLevelValF       = 5.5                 ; max level
plot(0)  = gsn_csm_contour(wks, fld_lowres, res )   ; place holder
fld0@units="ERA5"
plot(1)  = gsn_csm_contour(wks, fld0, res )   ; place holder
fld0 = fld_lowres-fld0
fld0@units="SAM-ERA5"
;symMinMaxPlt (fld0,20,False,res)
res@cnLevelSpacingF      = 0.2                  ; contour spacing
res@cnMinLevelValF       = -1.                ; min level
res@cnMaxLevelValF       = 1.                 ; max level
plot(2)  = gsn_csm_contour(wks, fld0, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Meridional Wind"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

; ------ T
f=addfile(filenameroot+"_TABS.nc","r")
pres = f->p
nlat = dimsizes(f->lat)
fld = dim_avg_Wrap(f->TABS(0,:,:,:))
fld&z = pres
fld&lat@units = "degrees_north"
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
fld_lowres = linint1_Wrap(yi,fld,True,yo,0)
print("t: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs=dim_avg(f_t->T(mmm,nlevi-1:0:1,nlati-1:0:1,:))
fld_obs_lowres = linint1(lati,fld_obs,True,yo,0)
fld_obs_lowres!0 = "lev"
fld_obs_lowres&lev = presi(nlevi-1:0:1)
fld_obs_lowres!1 = "lat"
fld_obs_lowres&lat = yo
print("t_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld0 = int2p_n_Wrap(presi(nlevi-1:0:1),fld_obs_lowres,pres,-1,0)
print("t_obs0:"+min(fld0)+"  "+max(fld0))
fld_lowres@units=" SAM"
res@cnLevelSpacingF      = 5.                  ; contour spacing
res@cnMinLevelValF       = 190.                ; min level
res@cnMaxLevelValF       = 300.                 ; max level
plot(0)  = gsn_csm_contour(wks, fld_lowres, res )   ; place holder
fld0@units="ERA5"
plot(1)  = gsn_csm_contour(wks, fld0, res )   ; place holder
fld0 = fld_lowres-fld0
fld0@units="SAM-ERA5"
;symMinMaxPlt (fld0,20,False,res)
res@cnLevelSpacingF      = 0.4                  ; contour spacing
res@cnMinLevelValF       = -2.                ; min level
res@cnMaxLevelValF       = 2.                 ; max level
plot(2)  = gsn_csm_contour(wks, fld0, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

; ------ Q
f=addfile(filenameroot+"_QV.nc","r")
pres = f->p
nlat = dimsizes(f->lat)
fld = dim_avg_Wrap(f->QV(0,:,:,:))
fld&z = pres
fld&lat@units = "degrees_north"
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
fld_lowres = linint1_Wrap(yi,fld,True,yo,0)
print("q: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs=dim_avg(f_q->Q(mmm,nlevi-1:0:1,nlati-1:0:1,:))
fld_obs = fld_obs*1.e3
fld_obs_lowres = linint1(lati,fld_obs,True,yo,0)
fld_obs_lowres!0 = "lev"
fld_obs_lowres&lev = presi(nlevi-1:0:1)
fld_obs_lowres!1 = "lat"
fld_obs_lowres&lat = yo
print("q_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld0 = int2p_n_Wrap(presi(nlevi-1:0:1),fld_obs_lowres,pres,-1,0)
print("q_obs0:"+min(fld0)+"  "+max(fld0))
fld_lowres@units=" SAM"
res@cnLevelSpacingF      = 1.                  ; contour spacing
res@cnMinLevelValF       = 0.                ; min level
res@cnMaxLevelValF       = 18.                 ; max level
plot(0)  = gsn_csm_contour(wks, fld_lowres, res )   ; place holder
fld0@units="ERA5"
plot(1)  = gsn_csm_contour(wks, fld0, res )   ; place holder
fld0 = fld_lowres-fld0
fld0@units="SAM-ERA5"
;symMinMaxPlt (fld0,20,False,res)
res@cnLevelSpacingF      = 0.4                  ; contour spacing
res@cnMinLevelValF       = -2.                ; min level
res@cnMaxLevelValF       = 2.                 ; max level
plot(2)  = gsn_csm_contour(wks, fld0, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Water Vapor"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down


; ------ QC/QI
f=addfile(filenameroot+"_QC.nc","r")
f1=addfile(filenameroot+"_QI.nc","r")
pres = f->p
nlat = dimsizes(f->lat)
fld = dim_avg_Wrap(f->QC(0,:,:,:))
fld = fld+dim_avg_Wrap(f1->QI(0,:,:,:))
fld&z = pres
fld&lat@units = "degrees_north"
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
fld_lowres = linint1_Wrap(yi,fld,True,yo,0)
print("q: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs=dim_avg(f_clwc->CLWC(mmm,nlevi-1:0:1,nlati-1:0:1,:))
fld_obs=fld_obs+dim_avg(f_ciwc->CIWC(mmm,nlevi-1:0:1,nlati-1:0:1,:))
fld_obs = fld_obs*1.e3
fld_obs_lowres = linint1(lati,fld_obs,True,yo,0)
fld_obs_lowres!0 = "lev"
fld_obs_lowres&lev = presi(nlevi-1:0:1)
fld_obs_lowres!1 = "lat"
fld_obs_lowres&lat = yo
print("q_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld0 = int2p_n_Wrap(presi(nlevi-1:0:1),fld_obs_lowres,pres,-1,0)
print("q_obs0:"+min(fld0)+"  "+max(fld0))
fld_lowres@units=" SAM"
res@cnLevelSpacingF      = 0.01                  ; contour spacing
res@cnMinLevelValF       = 0.                ; min level
res@cnMaxLevelValF       = 0.1                 ; max level
plot(0)  = gsn_csm_contour(wks, fld_lowres, res )   ; place holder
fld0@units="ERA5"
plot(1)  = gsn_csm_contour(wks, fld0, res )   ; place holder
fld0 = fld_lowres-fld0
fld0@units="SAM-ERA5"
symMinMaxPlt (fld0,20,False,res)
plot(2)  = gsn_csm_contour(wks, fld0, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Cloud Water/Ice"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down



