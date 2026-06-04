f1 = addfile("/glade/campaign/univ/unsb0027/SAM7.7_DYAMOND2/OUT_EXTR/SAM_20200121-20200125_atm_2D.daily_360x180_davg.nc","r")
f2 = addfile("../../OBS/CERES.SYN1deg.20200121-20200125.davg.nc","r")
f3 = addfile("../../OBS/e5.jan21-25.2020_mean.nc","r")

deg2rad = acos(-1.)/180.

wks   = gsn_open_wks ("ps","ps")
res                  = True                     ; plot mods desired

fld1 = f1->LWNT(:,:)
print("LWNT SAM: "+min(fld1)+" "+max(fld1))
; compute averaging weights
lat1 = f1->lat
g1 = cos(lat1*deg2rad)
g1 = g1/sum(g1)
fld1_zonal = dim_avg(fld1)
print(sum(fld1_zonal*g1))

fld2 = f2->toa_lw_all_daily(0,:,:)
print("olr : "+min(fld2)+" "+max(fld2))
; compute averaging weights
lat2 = f2->lat(:)
g2 = cos(lat2*deg2rad)
g2 = g2/sum(g2)
fld2_zonal = dim_avg(fld2)
print(sum(fld2_zonal*g2))

fld3 = f3->olr(0,:,:)
print("olr : "+min(fld3)+" "+max(fld3))
; compute averaging weights
lat3 = f3->latitude(:)
print(""+dimsizes(lat3))
g3 = cos(lat3*deg2rad)
g3 = g3/sum(g3)
fld3_zonal = dim_avg(fld3)
print(sum(fld3_zonal*g3))
fld4=fld1_zonal
nlat = dimsizes(lat1)
do m=0,nlat-1
  fld4(m) = avg(fld3_zonal(m*4:m*4+3))
end do

fld = new((/3,nlat/),float)

res@tiMainString     = "OLR"          ; add title
colors = (/"blue","black","red"/)
dashes = (/0,0,0,0,0,0,0,0,0/)
res@xyLineThicknessF = 3.
res@xyLineColors = colors
res@xyDashPatterns = dashes
res@vpWidthF            = 0.7           ; change aspect ratio of plot
res@vpHeightF           = 0.35
fld(0,:) = fld1_zonal
fld(1,:) = fld2_zonal
fld(2,:) = -fld4(nlat-1:0:1)

plot  = gsn_csm_xy (wks,lat1,fld,res) ; create plot

fld1 = f1->SWNT(:,:)
print("SWNT SAM: "+min(fld1)+" "+max(fld1))
; compute averaging weights
lat1 = f1->lat
g1 = cos(lat1*deg2rad)
g1 = g1/sum(g1)
fld1_zonal = dim_avg(fld1)
print(sum(fld1_zonal*g1))

fld2 = f2->toa_sw_all_daily(0,:,:)
print("asr : "+min(fld2)+" "+max(fld2))
; compute averaging weights
lat2 = f2->lat(:)
g2 = cos(lat2*deg2rad)
g2 = g2/sum(g2)
fld2_zonal = dim_avg(fld2)
print(sum(fld2_zonal*g2))

fld3 = f3->ntsrf(0,:,:)
print("asr : "+min(fld3)+" "+max(fld3))
; compute averaging weights
lat3 = f3->latitude(:)
print(""+dimsizes(lat3))
g3 = cos(lat3*deg2rad)
g3 = g3/sum(g3)
fld3_zonal = dim_avg(fld3)
print(sum(fld3_zonal*g3))
do m=0,nlat-1
  fld4(m) = avg(fld3_zonal(m*4:m*4+3))
end do

res@tiMainString     = "ASR"          ; add title
res@xyLineThicknessF = 3.
res@xyLineColors = colors
res@xyDashPatterns = dashes
res@vpWidthF            = 0.7           ; change aspect ratio of plot
res@vpHeightF           = 0.35
fld(0,:) = fld1_zonal
fld(1,:) = fld2_zonal
fld(2,:) = fld4(nlat-1:0:1)
fld = where(fld.lt.0.,0.,fld)

plot  = gsn_csm_xy (wks,lat1,fld,res) ; create plot
