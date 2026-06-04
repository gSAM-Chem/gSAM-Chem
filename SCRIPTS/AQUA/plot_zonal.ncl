load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;=================================================;
; open file and read in data
;=================================================;

filename = (/ \
"../../OUT_2D/AQUA_2048x512x38_za30_seaice1m_1024_AVG.2D_atm.nc", \
"../../OUT_2D/AQUA_2048x512x38_za30_seaice1m_3xco2_1024_AVG.2D_atm.nc", \ 
"../../OUT_2D/AQUA_2048x512x38_za50_seaice1m_radhomoz_1024_AVG.2D_atm.nc", \
"../../OUT_2D/AQUA_2048x512x38_za50_seaice1m_radhomoz_3xco2_1024_AVG.2D_atm.nc" \
/)

lat_max = 90.

  wks  = gsn_open_wks("ps","ps")    ; open a ncgm file
  gsn_define_colormap(wks, "ViBlGrWhYeOrRe")   ; generate new color map

  res                     = True          ; plot mods desired

 res@pmLegendDisplayMode    = "Always"            ; turn on legend

 res@pmLegendSide           = "Top"               ; Change location of
 res@pmLegendParallelPosF   = .5                 ; move units right
 res@pmLegendOrthogonalPosF = -1.2                ; move units down

 res@pmLegendWidthF         = 0.15                ; Change width and
 res@pmLegendHeightF        = 0.15                ; height of legend.
 res@lgLabelFontHeightF     = .015                 ; change font height
 res@lgTitleOn              = True                ; turn on legend title
 res@lgTitleString          = ""           ; create legend title
 res@lgTitleFontHeightF     = .025                ; font of legend title
 res@lgPerimOn              = False
 res@xyExplicitLegendLabels = (/"1xCO2"," 3xCO2","1xCO2+RADHONOZ","3xCO2+RADHONOZ"/) 

 colors = (/"blue","red","blue","red","cyan","burlywood","green","orange"/)
dashes = (/0,0,1,1,0,0,0,0,0/)
 res@xyLineThicknessF = 3.
 res@xyLineColors = colors
 res@xyDashPatterns = dashes
; res@trXMinF=-45.
; res@trXMaxF=45.
 res@vpWidthF            = 0.7           ; change aspect ratio of plot
 res@vpHeightF           = 0.35


res@pmLegendDisplayMode    = "Never"            ; turn on legend
nfiles = dimsizes(filename)
fm = new((/nfiles,nlat1/),float)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat) 
 if(i.eq.0) then
  fld = new((/nfiles,nlat/),float)
  lat = f->lat
  dlat = lat
  nlat = dimsizes(lat)
  dlat(0)=0.5*(lat(0)+lat(1))+90.
  do j=1,nlat-2
   dlat(j) = 0.5*(lat(j+1)-lat(j-1))
  end do
  dlat(nlat-1)=dlat(nlat-2)
  dlat(0)=dlat(1)
  print(lat+"  "+dlat)
  wgt = dlat*cos(lat*acos(-1.)/180.)
  wgt = wgt/sum(wgt)
  wgt!0 = "lat"
  wgt&lat = lat
 end if
 fld(i,:) = dim_avg_Wrap(f->Preca(0,:,:))
 print("PREC:"+sum(fld(i,:)*wgt(:)))
 do j=0,nlat1-1
    fm(i,j) = dim_avg(fld(i,{lat1(j)-0.5:lat1(j)+0.5}))
 end do
 fm!1="lat"
 fm&lat = lat1
end do
res@tiMainString          = "precipitation" 
plot = gsn_csm_xy (wks,lat1,fm,res)

res@pmLegendDisplayMode    = "Always"            ; turn on legend
do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->Preca(0,:,:))
 print("Preca:"+sum(fld(i,:)*wgt(:)))
end do
res@tiMainString          = "precipitation"
plot = gsn_csm_xy (wks,lat,fld,res)
res@pmLegendDisplayMode    = "Never"            ; turn on legend


res@pmLegendDisplayMode    = "Always"            ; turn on legend
do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->PW(0,:,:))
 print("PW:"+sum(fld(i,:)*wgt(:)))
end do
res@tiMainString          = "precipitable water" 
plot = gsn_csm_xy (wks,lat,fld,res)
res@pmLegendDisplayMode    = "Never"            ; turn on legend

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->CLD(0,:,:))
 print("CLD:"+sum(fld(i,:)*wgt(:)))
end do
res@tiMainString          = "Cloud Fraction"
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
; if(i.eq.0.) then
;  fld(i,:) = 27.
; else
  f=addfile(filename(i),"r")
  nlat = dimsizes(f->lat)
  fld(i,:) = dim_avg(f->SKT(0,:,:))
; end if
end do
res@tiMainString          = "SST"
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->LHF(0,:,:))
 fld(i,:) = fld(i,:) + dim_avg(f->SHF(0,:,:))
 print("LHF+SHF:"+sum(fld(i,:)*wgt(:)))
end do
res@tiMainString          = "lhf+shf" 
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->U200(0,:,:))
end do
res@tiMainString          = "u200" 
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->U850(0,:,:))
end do
res@tiMainString          = "u850" 
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->SOLIN(0,:,:))
 print("SOLIN:"+sum(fld(i,:)*wgt(:)))
 print("SOLIN (tropics):"+res@xyExplicitLegendLabels(i)+" "+sum(fld(i,{-lat_max:lat_max})*wgt({-lat_max:lat_max}))/sum(wgt({-lat_max:lat_max})))
end do
res@tiMainString          = "SOLIN"
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->LWNT(0,:,:))
 print("OLR:"+sum(fld(i,:)*wgt(:)))
 print("OLR (tropics):"+res@xyExplicitLegendLabels(i)+" "+sum(fld(i,{-lat_max:lat_max})*wgt({-lat_max:lat_max}))/sum(wgt({-lat_max:lat_max})))
end do
res@tiMainString          = "OLR" 
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->SWNT(0,:,:))
 print("ASR:"+sum(fld(i,:)*wgt(:)))
 print("ASR (tropics):"+res@xyExplicitLegendLabels(i)+" "+sum(fld(i,{-lat_max:lat_max})*wgt({-lat_max:lat_max}))/sum(wgt({-lat_max:lat_max})))
end do
res@tiMainString          = "ASR" 
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->SWNT(0,:,:))
 fld(i,:) = fld(i,:) - dim_avg(f->LWNT(0,:,:))
 print("TOP:"+sum(fld(i,:)*wgt(:)))
 print("TOP (tropics):"+res@xyExplicitLegendLabels(i)+" "+sum(fld(i,{-lat_max:lat_max})*wgt({-lat_max:lat_max}))/sum(wgt({-lat_max:lat_max})))
end do
res@tiMainString          = "TOP" 
plot = gsn_csm_xy (wks,lat,fld,res)

do i=0,nfiles-1
 f=addfile(filename(i),"r")
 nlat = dimsizes(f->lat)
 fld(i,:) = dim_avg(f->SWNS(0,:,:))
 fld(i,:) = fld(i,:) - dim_avg(f->LWNS(0,:,:))
 fld(i,:) = fld(i,:) - dim_avg(f->LHF(0,:,:))
 fld(i,:) = fld(i,:) - dim_avg(f->SHF(0,:,:))
 print("SFC:"+sum(fld(i,:)*wgt(:)))
end do
res@tiMainString          = "SFC"
plot = gsn_csm_xy (wks,lat,fld,res)





