f1 = addfile("/glade/campaign/univ/unsb0027/SAM7.7_DYAMOND2/OUT_EXTR/SAM_20200121-20200125_atm_2D.daily_1440x720_davg.nc","r")
f2 = addfile("../../OBS/e5.jan21-25.2020_mean.nc","r")
f3 = addfile("../../OBS/3B-HHR.MS.MRG.3IMERG.20200121-20200125.davg_1440x720.nc","r")

deg2rad = acos(-1.)/180.

prec1 = f1->Preca(:,:)
print("prec SAM: "+min(prec1)+" "+max(prec1))
; compute averaging weights
lat1 = f1->lat
g1 = cos(lat1*deg2rad)
g1 = g1/sum(g1)
prec1_zonal = dim_avg(prec1)
print(sum(prec1_zonal*g1))


prec2 = f2->pr(0,1:,:)
print("prec ERA5: "+min(prec2)+" "+max(prec2))
; compute averaging weights
lat2 = f2->latitude(1:)
g2 = cos(lat2*deg2rad)
g2 = g2/sum(g2)
prec2_zonal = dim_avg(prec2)
print(sum(prec2_zonal*g2))

prec3 = f3->IRprecipitation(:,:)
prec3 = prec3*24.
print("prec GPM: "+min(prec3)+" "+max(prec3))
lat3 = f3->lat
prec3_zonal = dim_avg(prec3)

nlat = dimsizes(lat1)
fld = new((/3,nlat/),float)
fld(0,:) = prec1_zonal
fld(1,:) = prec2_zonal(nlat-1:0:1)
fld(2,:) = prec3_zonal
print(prec2_zonal)
print(lat1+"    "+fld(0,:)+"    "+fld(1,:)+"    "+fld(2,:))

wks   = gsn_open_wks ("ps","ps")
res                  = True                     ; plot mods desired
res@tiMainString     = "Precipitation"          ; add title
colors = (/"blue","red","black"/)
dashes = (/0,0,0,0,0,0,0,0,0/)
res@xyLineThicknessF = 3.
res@xyLineColors = colors
res@xyDashPatterns = dashes
res@vpWidthF            = 0.7           ; change aspect ratio of plot
res@vpHeightF           = 0.35

plot  = gsn_csm_xy (wks,lat1,fld,res) ; create plot

fld1 = f1->T2m(:,:)
fld1_zonal = dim_avg(fld1)
print(sum(fld1_zonal*g1))
fld2 = f2->t2m(0,1:,:)
fld2_zonal = dim_avg(fld2)
print(sum(fld2_zonal*g2))

res@tiMainString     = "2-m Temperature"          ; add title
colors = (/"blue","red","black"/)
dashes = (/0,0,0,0,0,0,0,0,0/)
res@xyLineThicknessF = 3.
res@xyLineColors = colors
res@xyDashPatterns = dashes
res@vpWidthF            = 0.75           ; change aspect ratio of plot
res@vpHeightF           = 0.35
fld(0,:) = fld1_zonal
fld(1,:) = fld2_zonal(nlat-1:0:1)
plot  = gsn_csm_xy (wks,lat1,fld(0:1,:),res) ; create plot

fld1 = f1->LHF(:,:)
fld1_zonal = dim_avg(fld1)
print(sum(fld1_zonal*g1))
fld2 = f2->lhf(0,1:,:)
fld2_zonal = dim_avg(fld2)
print(sum(fld2_zonal*g2))

res@tiMainString     = "LHF"          ; add title
colors = (/"blue","red","black"/)
dashes = (/0,0,0,0,0,0,0,0,0/)
res@xyLineThicknessF = 3.
res@xyLineColors = colors
res@xyDashPatterns = dashes
res@vpWidthF            = 0.75           ; change aspect ratio of plot
res@vpHeightF           = 0.35
fld(0,:) = fld1_zonal
fld(1,:) = -fld2_zonal(nlat-1:0:1)
plot  = gsn_csm_xy (wks,lat1,fld(0:1,:),res) ; create plot

