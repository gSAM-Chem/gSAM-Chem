load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; SAM datset
filenameroot = "/glade/p/univ/unsb0007/SAM7.3_SOCRATES/OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nonudge_SAM1MOM_gmg6_4608_0000023040"
;-----------------------------------------------------------------------------------------
; ERA5 datset

;date = 2018021700
date = 2018021800
;date = 2018022100

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

;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("ps","ps")             ; send graphics to PNG file
  gsn_define_colormap(wks, "wgne15")

  res            = True                         ; plot mods desired
  res@gsnPolar   = "NH"                         ; specify the hemisphere
  res@mpMinLatF  = 60                           ; minimum lat to plot
  res@mpFillOn   = False

  res@cnFillOn          = True                  ; color fill
  res@cnLevelSpacingF   = 3.                    ; interval spacing
  res@cnLinesOn         = False                 ; no contour lines

;  res@gsnSpreadColorEnd     = -3
  res@gsnScalarContour      = True              ; vectors over scalars
  res@vcRefMagnitudeF       = 20.0              ; vec magnitude
  res@vcRefLengthF          = 0.04              ; size of ref vector
  res@vcLineArrowColor      = "white"           ; white vectors
  res@vcMinDistanceF        = 0.02              ; thins arrows near pole  
  res@vcGlyphStyle          = "CurlyVector"     ; turns on curly vectors

  res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
  res@vcRefAnnoArrowUseVecColor = False
  res@cnFillPalette        = "BkBlAqGrYeOrReViWh200"
;  res@gsnLeftString       = "Wind Velocity over Surface Temperature"
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frame
res@cnFillMode = "RasterFill"


 plot = new(2,graphic)
; because this is ice data, which has a gap in the tropics, we need to
; explicitly pass the range of the data to plot.  Since we are coloring 
; the vectors, this range should also match the MinLatF above, since the
; range for the colors is chosen over the full data passed, and not the
; map limits.
 f=addfile(filenameroot+"_U.nc","r")
 pres = f->p
 p0=850.
 ilev = closest_val(p0,pres)
 print(pres(ilev)+"  "+ilev)
 u = f->U(0,ilev,:,:)
 f=addfile(filenameroot+"_V.nc","r")
 v = f->V(0,ilev,:,:)
 f=addfile(filenameroot+"_TABS.nc","r")
 t = f->TABS(0,ilev,:,:)
 u&lat@units = "degrees north"
 v&lat@units = "degrees north"
 t&lat@units = "degrees north"
 u&lon@units = "degrees east"
 v&lon@units = "degrees east"
 t&lon@units = "degrees east"
 res@cnLevelSelectionMode = "ManualLevels"
 res@cnMinLevelValF       =  220.                 ; set min contour level
 res@cnMaxLevelValF       =  280.                  ; set max contour level
 res@cnLevelSpacingF      =   1.     
 plot(0) = gsn_csm_vector_scalar_map_polar(wks,u({62.:90.},:),v({62.:90.},:), \
         t({62.:90.},:),res)

 ilevi = closest_val(pres(ilev),presi)
 if(presi(ilevi).ge.p0) then
  ilevi1 = ilevi+1
 else
  ilevi1 = ilevi-1
 end if
 print(presi(ilevi)+"  "+ilevi+"  "+presi(ilevi1)+"  "+ilevi1)
 u_obs=f_u->U(mmm,ilevi,nlati-1:0:1,:)
 u_obs1=f_u->U(mmm,ilevi1,nlati-1:0:1,:)
 u_obs = u_obs1+(u_obs-u_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
 print("u_obs:"+min(u_obs)+"  "+max(u_obs))
 v_obs=f_v->V(mmm,ilevi,nlati-1:0:1,:)
 v_obs1=f_v->V(mmm,ilevi1,nlati-1:0:1,:)
 v_obs = v_obs1+(v_obs-v_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
 print("v_obs:"+min(v_obs)+"  "+max(v_obs))
 t_obs=f_t->T(mmm,ilevi,nlati-1:0:1,:)
 t_obs1=f_t->T(mmm,ilevi1,nlati-1:0:1,:)
 t_obs = t_obs1+(t_obs-t_obs1)*log(pres(ilev)/presi(ilevi1))/log(presi(ilevi)/presi(ilevi1))
 print("t_obs:"+min(t_obs)+"  "+max(t_obs))
 plot(1) = gsn_csm_vector_scalar_map_polar(wks,u_obs({62.:90.},:),v_obs({62.:90.},:), \
         t_obs({62.:90.},:),res)
 
 resp             = True                    ; panel mods desired
 resp@gsnMaximize = True                    ; fill the page
 resp@txString   = ""
 gsn_panel(wks,plot,(/1,2/),resp)           ; draw: 2-across, 2-down

