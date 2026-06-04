load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; SAM datset
filenameroot = "../..//OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nudge24h_SAM1MOM_gmg6_newcanopyinit_4608_0000000001"
;-----------------------------------------------------------------------------------------
; ERA5 datset

;date = 2018021700
date = 2018021600

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
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
plot = new(3,graphic) 
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frame
res@cnLevelSelectionMode = "ManualLevels" 
res1 = res
;-----------------------------------------------------------------------------------------

g=addfile(filenameroot+"_TERRA.nc","r")

; ------ U950
f=addfile(filenameroot+"_U.nc","r")
pres = f->p
print(pres)
p0=950.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->U(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
xi = fld&lon
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlon=lowres_nx
dlon = 360. / mlon
xo  = fspan (0, (mlon - 1) * dlon, mlon)
xo@units = "degrees east"
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
print("terr_lowre:"+min(terr_lowres)+" "+max(terr_lowres))
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("u: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_u->U(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_u->U(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("u_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res) 
fld_lowres@units=" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
;symMinMaxPlt (dfld,20,False,res) 
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind 950 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down


; ------ V850
f=addfile(filenameroot+"_V.nc","r")
p0=850.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->V(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("v: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_v->V(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_v->V(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("v_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Meridional Wind 850 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down


; ------ U500
f=addfile(filenameroot+"_U.nc","r")
p0=500.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->U(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("u: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_u->U(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_u->U(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("u_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind 500 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down


; ------ V500
f=addfile(filenameroot+"_V.nc","r")
p0=500.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->V(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("v: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_v->V(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_v->V(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("v_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units="ERA5"
;plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Meridional Wind 500 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

; ------ U200
f=addfile(filenameroot+"_U.nc","r")
p0=200.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->U(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("u: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_u->U(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_u->U(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("v_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind 200 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

; ------ V200
f=addfile(filenameroot+"_V.nc","r")
p0=200.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->V(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("v: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_v->V(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_v->V(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("v_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Meridional Wind 200 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

; ------ T1000
f=addfile(filenameroot+"_TABS.nc","r")
p0=1000.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("t: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
print(presi(ilevi)+"  "+ilevi)
fld_obs=f_t->T(mmm,ilevi,nlati-1:0:1,:)
print("t_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  200.                 ; set min contour level
res1@cnMaxLevelValF       =  320.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -5.                 ; set min contour level
res@cnMaxLevelValF       =  5.                  ; set max contour level
res@cnLevelSpacingF      =   0.5                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature 1000 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ T850
f=addfile(filenameroot+"_TABS.nc","r")
p0=850.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("t: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_t->T(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_t->T(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("t_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  200.                 ; set min contour level
res1@cnMaxLevelValF       =  320.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -5.                 ; set min contour level
res@cnMaxLevelValF       =  5.                  ; set max contour level
res@cnLevelSpacingF      =   0.5                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature 850 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ T500
f=addfile(filenameroot+"_TABS.nc","r")
p0=500.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("t: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_t->T(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_t->T(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("t_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  210.                 ; set min contour level
res1@cnMaxLevelValF       =  280.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -5.                 ; set min contour level
res@cnMaxLevelValF       =  5.                  ; set max contour level
res@cnLevelSpacingF      =   0.5                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature 500 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ T200
f=addfile(filenameroot+"_TABS.nc","r")
p0=200.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("t: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_t->T(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_t->T(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
print("t_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  180.                 ; set min contour level
res1@cnMaxLevelValF       =  240.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -5.                 ; set min contour level
res@cnMaxLevelValF       =  5.                  ; set max contour level
res@cnLevelSpacingF      =   0.5                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature 200 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ Q1000
f=addfile(filenameroot+"_QV.nc","r")
p0=1000.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->QV(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("q: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
print(presi(ilevi)+"  "+ilevi)
fld_obs=f_q->Q(mmm,ilevi,nlati-1:0:1,:)
fld_obs = fld_obs*1.e3
print("q_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  20.                  ; set max contour level
res1@cnLevelSpacingF      =   1.                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -3.                 ; set min contour level
res@cnMaxLevelValF       =  3.                  ; set max contour level
res@cnLevelSpacingF      =   0.3                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Vapor 1000 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ Q850
f=addfile(filenameroot+"_QV.nc","r")
p0=850.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->QV(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("q: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_q->Q(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_q->Q(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
fld_obs = fld_obs*1.e3
print("q_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  17.                  ; set max contour level
res1@cnLevelSpacingF      =   1.                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -5.                 ; set min contour level
res@cnMaxLevelValF       =  5.                  ; set max contour level
res@cnLevelSpacingF      =   0.5                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Vapor 850 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ Q500
f=addfile(filenameroot+"_QV.nc","r")
p0=500.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->QV(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("q: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_q->Q(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_q->Q(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
fld_obs = fld_obs*1.e3
print("q_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  7.                  ; set max contour level
res1@cnLevelSpacingF      =  0.5                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -3.                 ; set min contour level
res@cnMaxLevelValF       =  3.                  ; set max contour level
res@cnLevelSpacingF      =   0.3                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Vapor 500 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ Q200
f=addfile(filenameroot+"_QV.nc","r")
p0=200.
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->QV(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("q: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_q->Q(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_q->Q(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
fld_obs = fld_obs*1.e3
print("q_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  0.1                  ; set max contour level
res1@cnLevelSpacingF      =   0.01                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
;symMinMaxPlt (dfld,20,False,res)
res@cnMinLevelValF       =  -0.05                 ; set min contour level
res@cnMaxLevelValF       =  0.05                  ; set max contour level
res@cnLevelSpacingF      =   0.005                  ; set contour spacing
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Vapor 200 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ QCL900
p0=900.
f=addfile(filenameroot+"_QC.nc","r")
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->QC(0,ilev,:,:)
f=addfile(filenameroot+"_QI.nc","r")
fld = fld+f->QI(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:) 
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"                
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("qc: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)   
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_clwc->CLWC(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_clwc->CLWC(mmm,ilevi1,nlati-1:0:1,:)
fld_obs=fld_obs+f_ciwc->CIWC(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=fld_obs1+f_ciwc->CIWC(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
fld_obs = fld_obs*1.e3
print("qc_obs:"+min(fld_obs)+"  "+max(fld_obs))     
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  0.1                  ; set max contour level
res1@cnLevelSpacingF      =  0.01                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
symMinMaxPlt (dfld,20,False,res)
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Cloud Water/Ice 900 mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ QCL200
p0=200.
f=addfile(filenameroot+"_QC.nc","r")
ilev = closest_val(p0,pres)
print(pres(ilev)+"  "+ilev)
nlat = dimsizes(f->lat)
fld = f->QC(0,ilev,:,:)
f=addfile(filenameroot+"_QI.nc","r")
fld = fld+f->QI(0,ilev,:,:)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees_east"
terr&lat@units = "degrees_north"
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
print("qc: "+min(fld_lowres)+"  "+max(fld_lowres))
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
ilevi = closest_val(pres(ilev),presi)
if(presi(ilevi).ge.p0) then
 ilevi1 = ilevi+1
else
 ilevi1 = ilevi-1
end if
print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
fld_obs=f_clwc->CLWC(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=f_clwc->CLWC(mmm,ilevi1,nlati-1:0:1,:)
fld_obs=fld_obs+f_ciwc->CIWC(mmm,ilevi,nlati-1:0:1,:)
fld_obs1=fld_obs1+f_ciwc->CIWC(mmm,ilevi1,nlati-1:0:1,:)
;fld_obs = fld_obs1+(fld_obs-fld_obs1)*(pres(ilev)-presi(ilevi1))/(presi(ilevi)-presi(ilevi1))
fld_obs = fld_obs1+(fld_obs-fld_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
fld_obs = fld_obs*1.e3
print("qc_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  0.1                  ; set max contour level
res1@cnLevelSpacingF      =  0.01                  ; set contour spacing
fld_lowres@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
symMinMaxPlt (dfld,20,False,res)
dfld@units="SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Cloud Water/Ice 200 mb"
gsn_panel(wks,plot,(/3,1/),resp)


