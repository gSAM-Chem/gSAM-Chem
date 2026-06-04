; modle output:
filenames = (/"../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR12_2019", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR11_2019", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_2019"/)

; ERA5 dataset
fobs = addfile("/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/201904/e5.oper.an.pl.128_129_z.ll025sc.2019041400_2019041423.nc","r")

nfilenames = dimsizes(filenames)

wks = gsn_open_wks("ps","ps")

res = True
res@gsnDraw             = False
res@gsnFrame            = False
res@cnFillOn            = False
res@cnLinesOn           = True
res@cnLineThicknessF    = 4.0
res@cnInfoLabelOn       = False
res@lbLabelBarOn        = False
res@cnLevelSelectionMode = "ManualLevels"


color = (/"black","red","blue","green","magenta"/)

;---Base map (empty), so we can overlay multiple contours of different grids
res@gsnMaximize          = True
res@mpOutlineOn          = True
res@mpGeophysicalLineThicknessF = 1.5


; plot observations first:


do j=0,nfilenames-1

 files = systemfunc("ls -1 "+filenames(j)+"*.2D_atm.nc")
 nfiles = dimsizes(files)
 f = addfiles(files,"r")
 
 do i = 0,nfiles-1 
 fld = f[i]->PHI500(0,:,:)
 print("fld: "+min(fld)+" "+max(fld))
 
 if(i.eq.0.and.j.eq.0.) then
  latmin = min(f[i]->lat)
  latmax = max(f[i]->lat)
  lonmin = min(f[i]->lon)
  lonmax = max(f[i]->lon)
  res@mpMinLatF   = latmin
  res@mpMaxLatF   = latmax
  res@mpMinLonF   = lonmin
  res@mpMaxLonF   = lonmax
 end if

; observation:
 if(i.eq.0.and.j.eq.0) then
   fldobs1 = fobs->Z(0,{500},{latmin:latmax},{lonmin:lonmax}) ; 500-mb geopotential
   fldobs = fldobs1(::-1,:)
   delete(fldobs1)
   fldobs = fldobs/9.797  ; from geopotential to height
   print("fldobs: "+min(fldobs)+" "+max(fldobs))
   ; find domain average:
   latobs = fldobs&latitude
   lonobs = fldobs&longitude
   fldobs1 = dim_avg_Wrap(fldobs)
   wobs = cos(latobs * 0.0174533)
   w2Dobs = conform(fldobs1, wobs, 0)
   fldobs_mean = dim_avg_wgt_n(fldobs1,w2Dobs,0,0)
   fldobs = fldobs - fldobs_mean
   print("fldobs perturb: "+min(fldobs)+" "+max(fldobs))
 end if


 ww = cos(fld&lat*0.0174533)
 fld_lo = area_hi2lores_Wrap(fld&lon, fld&lat, fld, False, ww, lonobs, latobs, False)
 ; find domain average:
 lat = fld_lo&latitude
 fld1 = dim_avg_Wrap(fld_lo)
 if(i.eq.0) then
  w = cos(lat * 0.0174533)
  w2D = conform(fld_lo, w, 0)
 end if
 fld_mean = dim_avg_wgt(fld1,w,0)
 fld_lo = fld_lo - fld_mean

; corr = corr + w2D*fld_lo*fldobs
; var = 


 print("fld perturb: "+min(fld_lo)+" "+max(fld_lo))

 

 end do
end do

