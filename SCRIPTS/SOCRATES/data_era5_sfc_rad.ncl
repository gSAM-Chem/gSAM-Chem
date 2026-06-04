load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

;-----------------------------------------------------------------------------------------
; ERA5 datset

;date = 2018021700
date = 2018021612

date_obs = "2018020100_2018022823"
rad = "../../GLOBAL_DATA/NC/era5_rad_20180216.nc"
nrad = 2
;-------------------------------------------------------------------
; plot resolution

resol="576x288_dyvar"

lowres_nx = 576
lowres_ny = 288
;lowres_nx = 128
;lowres_ny = 65

system("rm ../../GLOBAL_DATA/NC/era5_sfc_"+resol+"_"+date+".nc")
g=addfile("../../GLOBAL_DATA/NC/era5_sfc_"+resol+"_"+date+".nc","c")
frad = addfile(rad,"r")

f1 = addfile("../../GLOBAL_DATA/NC/landmask_"+resol+".nc","r")

;-----------------------------------------------------------------------------------------
; find corresponding ERA5 file names:

year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
day = floattointeger((date-year*1000000-month*10000)/100.)
hour = date-year*1000000-month*10000-day*100
;date_obs = (year*1000000+month*10000+day*100)+"_"+(year*1000000+month*10000+day*100+23)
month_obs = ""+(year*100+month)
print(""+date_obs)
dirpath = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/"+month_obs+"/"
f_skt = addfile(dirpath+"e5.oper.an.sfc.128_235_skt.ll025sc."+date_obs+".nc","r")
f_2t = addfile(dirpath+"e5.oper.an.sfc.128_167_2t.ll025sc."+date_obs+".nc","r")
f_10u = addfile(dirpath+"e5.oper.an.sfc.128_165_10u.ll025sc."+date_obs+".nc","r")
f_10v = addfile(dirpath+"e5.oper.an.sfc.128_166_10v.ll025sc."+date_obs+".nc","r")
f_sd = addfile(dirpath+"e5.oper.an.sfc.128_141_sd.ll025sc."+date_obs+".nc","r")
f_stl1 = addfile(dirpath+"e5.oper.an.sfc.128_139_stl1.ll025sc."+date_obs+".nc","r")
f_swvl1 = addfile(dirpath+"e5.oper.an.sfc.128_039_swvl1.ll025sc."+date_obs+".nc","r")
f_laihv = addfile(dirpath+"e5.oper.an.sfc.128_067_laihv.ll025sc."+date_obs+".nc","r")
f_lailv = addfile(dirpath+"e5.oper.an.sfc.128_066_lailv.ll025sc."+date_obs+".nc","r")
f_shf = addfile(dirpath+"e5.oper.an.sfc.128_231_ishf.ll025sc."+date_obs+".nc","r")
f_lhf = addfile(dirpath+"e5.oper.an.sfc.128_232_ie.ll025sc."+date_obs+".nc","r")
f_taux = addfile(dirpath+"e5.oper.an.sfc.128_229_iews.ll025sc."+date_obs+".nc","r")
f_tauy = addfile(dirpath+"e5.oper.an.sfc.128_230_inss.ll025sc."+date_obs+".nc","r")
times = calendar_decode2(f_skt->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(""+dates)
ntime = dimsizes(f_skt->time)
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

nlati = dimsizes(f_skt->latitude)
nloni = dimsizes(f_skt->longitude)
grid = nloni+"x"+nlati
print("data grid:"+grid)

; determine pressure levels at z:

lati = doubletofloat(f_skt->latitude(nlati-1:0:1))
lati@units = f_skt->latitude@units
loni = doubletofloat(f_skt->longitude)
loni@units = f_skt->longitude@units
wgt = cos(lati*acos(-1.)/180.)

;-----------------------------------------------------------------------------------------
; prepare initial settings for plots:

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
res@cnLinesOn            = False
res@cnLevelSelectionMode = "ManualLevels" 
res@cnRasterModeOn = True
res1 = res
;-----------------------------------------------------------------------------------------

; ------ SKT
mlon=lowres_nx
dlon = 360. / mlon
xo  = fspan (0, (mlon - 1) * dlon, mlon)
xo@units = "degrees east"
mlat=lowres_ny
;yo  = fspan (-89., 89., mlat)
yo  = f1->lat
yo@units = "degrees north"

fld_obs=f_skt->SKT(mmm,nlati-1:0:1,:)
fld_obs = fld_obs-273.15
print("skt_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->SKT = fld_obs_lowres

fld_obs=f_2t->VAR_2T(mmm,nlati-1:0:1,:)
fld_obs = fld_obs-273.15
print("t2m_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->T2m = fld_obs_lowres

fld_obs=f_10u->VAR_10U(mmm,nlati-1:0:1,:)
print("u10_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->U10m = fld_obs_lowres


fld_obs=f_10v->VAR_10V(mmm,nlati-1:0:1,:)
print("v10_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->V10m = fld_obs_lowres

fld_obs=f_sd->SD(mmm,nlati-1:0:1,:)
print("sd_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1.                  ; set max contour level
res@cnLevelSpacingF      =  0.5                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->SD = fld_obs_lowres

fld_obs=f_stl1->STL1(mmm,nlati-1:0:1,:)
fld_obs = fld_obs - 273.15
print("stl1_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =  5.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->STL1 = fld_obs_lowres

fld_obs=f_swvl1->SWVL1(mmm,nlati-1:0:1,:)
print("swl1_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1.                  ; set max contour level
res@cnLevelSpacingF      =  0.1                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->SWVL1 = fld_obs_lowres

fld_obs=f_laihv->LAIHV(mmm,nlati-1:0:1,:)
print("laihv_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  7.                  ; set max contour level
res@cnLevelSpacingF      =  0.5                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->LAIHV = fld_obs_lowres


fld_obs=f_shf->ISHF(mmm,nlati-1:0:1,:)
fld_obs = -fld_obs
print("shf_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -400.                 ; set min contour level
res@cnMaxLevelValF       =  400.                  ; set max contour level
res@cnLevelSpacingF      =  10.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
g->SHF = fld_obs_lowres

fld_obs=f_lhf->IE(mmm,nlati-1:0:1,:)
fld_obs = -2.5e6*fld_obs
print("lhf_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -400.                 ; set min contour level
res@cnMaxLevelValF       =  400.                  ; set max contour level
res@cnLevelSpacingF      =  10.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "ilhf"
fld_obs_lowres@units = "W/m2"
g->LHF = fld_obs_lowres

fld_obs=f_taux->IEWS(mmm,nlati-1:0:1,:)
fld_obs = -fld_obs
print("taux_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =  1.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "taux"
fld_obs_lowres@units = "N/m2"
g->TAUX = fld_obs_lowres

fld_obs=f_tauy->INSS(mmm,nlati-1:0:1,:)
fld_obs = -fld_obs
print("tauy_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =  1.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "tauy"
fld_obs_lowres@units = "N/m2"
g->TAUY = fld_obs_lowres




delete(fld_obs)
delete(fld_obs_lowres)

fld_obs=short2flt(frad->tisr(nrad,nlati-1:0:1,:))
fld_obs = fld_obs/3600.
print("solin_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =  50.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "solin"
fld_obs_lowres@units = "W/m2"
g->SOLIN = fld_obs_lowres


fld_obs=short2flt(frad->tsr(nrad,nlati-1:0:1,:))
fld_obs = fld_obs/3600.
print("swnt_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =  50.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "swnt"
fld_obs_lowres@units = "W/m2"
g->SWNT = fld_obs_lowres

fld_obs=short2flt(frad->tsrc(nrad,nlati-1:0:1,:))
fld_obs = fld_obs/3600.
print("swnt_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =  50.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "swntc"
fld_obs_lowres@units = "W/m2"
g->SWNTC = fld_obs_lowres


fld_obs=short2flt(frad->ssr(nrad,nlati-1:0:1,:))
fld_obs = fld_obs/3600.
print("swns_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =  50.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "swns"
fld_obs_lowres@units = "W/m2"
g->SWNS = fld_obs_lowres

fld_obs=short2flt(frad->ttr(nrad,nlati-1:0:1,:))
fld_obs = -fld_obs/3600.
print("lwnt_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  500.                  ; set max contour level
res@cnLevelSpacingF      =  25.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "lwnt"
fld_obs_lowres@units = "W/m2"
g->LWNT = fld_obs_lowres

fld_obs=short2flt(frad->ttrc(nrad,nlati-1:0:1,:))
fld_obs = -fld_obs/3600.
print("lwntc_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  500.                  ; set max contour level
res@cnLevelSpacingF      =  25.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "lwntc"
fld_obs_lowres@units = "W/m2"
g->LWNTC = fld_obs_lowres

fld_obs=short2flt(frad->str(nrad,nlati-1:0:1,:))
fld_obs = -fld_obs/3600.
print("lwns_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  500.                  ; set max contour level
res@cnLevelSpacingF      =  25.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "swns"
fld_obs_lowres@units = "W/m2"
g->LWNS = fld_obs_lowres

fld_obs=short2flt(frad->alnid(nrad,nlati-1:0:1,:))
print("albnirdiff_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1.                  ; set max contour level
res@cnLevelSpacingF      =  0.1                 ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "alnid"
fld_obs_lowres@units = ""
g->ALBNIRDIF = fld_obs_lowres

fld_obs=short2flt(frad->alnip(nrad,nlati-1:0:1,:))
print("albnirdir:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1.                  ; set max contour level
res@cnLevelSpacingF      =  0.1                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "alnid"
fld_obs_lowres@units = ""
g->ALBNIRDIR = fld_obs_lowres

fld_obs=short2flt(frad->flsr(nrad,nlati-1:0:1,:))
print("flsr:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
fld_obs_lowres@_FillValue = -999.
res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  1.                  ; set max contour level
res@cnLevelSpacingF      =  1.                  ; set contour spacing
plot = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
fld_obs_lowres@short_name = "logz0t"
fld_obs_lowres@units = ""
g->LOGZ0T = fld_obs_lowres




