load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

;-----------------------------------------------------------------------------------------
; ERA5 datset

date = "20180217"
filename = "../../GLOBAL_DATA/NC/era5_rad_20180217_daily.nc"
;-------------------------------------------------------------------
; plot resolution


lowres_nx = 576
lowres_ny = 288

resol = lowres_nx+"x"+lowres_ny+"_dyvar"
system("rm ../../GLOBAL_DATA/NC/era5_sfc_"+resol+"_"+date+".nc")
g=addfile("../../GLOBAL_DATA/NC/era5_sfc_"+resol+"_"+date+".nc","c")

f1 = addfile("../../GLOBAL_DATA/NC/landmask_"+resol+".nc","r")

f = addfile(filename,"r")

;-----------------------------------------------------------------------------------------
; prepare grid for comparison (ERA5 native grid):

nlati = dimsizes(f->latitude)
nloni = dimsizes(f->longitude)
grid = nloni+"x"+nlati
print("data grid:"+grid)

lati = doubletofloat(f->latitude(nlati-1:0:1))
lati@units = f->latitude@units
loni = doubletofloat(f->longitude)
loni@units = f->longitude@units
wgt = cos(lati*acos(-1.)/180.)



;-----------------------------------------------------------------------------------------
; prepare initial settings for plots:

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
res@cnLinesOn            = False
res@cnLevelSelectionMode = "AutomaticLevels" 
res@cnRasterModeOn = True
;-----------------------------------------------------------------------------------------
names = getfilevarnames(f)

xo  = f1->lon
yo  = f1->lat
xo@units = "degrees east"
yo@units = "degrees north"

nt = dimsizes(f->time)

do i=0,dimsizes(names)-1

  if(dimsizes(dimsizes(f->$names(i)$)).eq.3) then
    print(""+names(i))
    fld_obs=dble2flt(f->$names(i)$(:,nlati-1:0:1,:))
    print("obs:"+min(fld_obs)+"  "+max(fld_obs))
    fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
    plot = gsn_csm_contour_map(wks,fld_obs_lowres(nt-1,:,:),res)        ; create plot
    g->$names(i)$ = fld_obs_lowres
   end if

end do




