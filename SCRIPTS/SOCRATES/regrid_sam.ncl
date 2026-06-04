load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

;-----------------------------------------------------------------------------------------
; ERA5 datset

filename_root = "../../OUT_2D/SOCRATES_9216x4608x74_THOM_4608_20180219020000_000"
N1 = 26640
N2 = 26640
NN = 360

names = (/"Preca","LWNTA","SWNTA","LHF","SHF","SOLINA","T2m","U10m","V10m","SST"/)

doplot = False

do N = N1,N2,NN

;if(N < 10) then
; M = "000000"
;else if(N < 100) then
; M = "00000"
;else if(N < 1000) then
; M = "0000"
;else if(N < 10000) then
; M = "000"
;else if(N < 100000) then
; M = "00"
;else if(N < 1000000) then
; M = "0"
;end if
;end if
;end if
;end if
;end if
;end if

M="00"
filename = filename_root+M+N
print(""+filename)


;-------------------------------------------------------------------
; resolution


lowres_nx = 576
lowres_ny = 288

resol = lowres_nx+"x"+lowres_ny+"_dyvar"
fileout = filename+"_"+resol+".nc"
system("rm "+fileout)
g=addfile(fileout,"c")

f1 = addfile("../../GLOBAL_DATA/NC/landmask_"+resol+".nc","r")

f = addfile(filename+".2D_1.nc","r")

;-----------------------------------------------------------------------------------------
; prepare grid for comparison (ERA5 native grid):

nlati = dimsizes(f->lat)
nloni = dimsizes(f->lon)
grid = nloni+"x"+nlati
print("data grid:"+grid)

lati = f->lat
lati@units = f->lat@units
loni = f->lon
loni@units = f->lon@units
wgt = cos(lati*acos(-1.)/180.)



;-----------------------------------------------------------------------------------------
; prepare initial settings for plots:

if(doplot) then
wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
res@cnLinesOn            = False
res@cnLevelSelectionMode = "AutomaticLevels" 
res@cnRasterModeOn = True
end if
;-----------------------------------------------------------------------------------------
;names = getfilevarnames(f)

xo  = f1->lon
yo  = f1->lat
xo@units = "degrees east"
yo@units = "degrees north"

do i=0,dimsizes(names)-1

  if(dimsizes(dimsizes(f->$names(i)$)).eq.3) then
    print(""+names(i))
    fld_obs=f->$names(i)$(:,:,:)
    print("obs:"+min(fld_obs)+"  "+max(fld_obs))
    fld_obs_lowres = area_hi2lores_Wrap(loni,lati,fld_obs,True,wgt,xo,yo,False)
    if(doplot) then
     plot = gsn_csm_contour_map(wks,fld_obs_lowres(0,:,:),res)        ; create plot
    end if
    g->$names(i)$ = fld_obs_lowres
  end if

end do

system("ncks -O --mk_rec_dmn time "+fileout+" tmp.nc") 
system("mv tmp.nc "+fileout)

end do ; NN

