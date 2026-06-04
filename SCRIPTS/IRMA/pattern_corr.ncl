

files = systemfunc("ls -1 ../../OUT_2D/IRMA_2304x1152x74_15s_01_mom2_impldamp_nofudge_2304_20170*0.2D_atm.nc")
nfiles = dimsizes(files)
f=addfiles(files,"r")

;g = addfile("/glade/work/marat/ERA5/out.nc","r")
g = addfile("/glade/work/marat/ERA5/ERA5_2017-05-11Sep_sfc_hourly_2304x1152.nc","r")


lat = f[0]->lat
pi = acos(-1.)
gw = lat
gw = cos(pi/180.*lat)
gw = gw/sum(gw)
land = f[0]->LANDMASK(0,:,:)

do i=0,nfiles-1

  fld1 = g->tcwv(i,:,:)
  fld2 = f[i]->PW(0,:,:)
  corr = pattern_cor(fld1,fld2,gw,0)
  rmse = wgt_arearmse(fld1,fld2, gw, 1.0, 0)
  mean1 = wgt_areaave(fld1, gw, 1.0, 1)
  mean2 = wgt_areaave(fld2, gw, 1.0, 1)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do

do i=0,nfiles-1

  fld1 = g->t2m(i,:,:)
  fld2 = f[i]->T2m(0,:,:)
  fld1 = where(land.gt.0.5,fld1,-999.) 
  fld2 = where(land.gt.0.5,fld2,-999.) 
  fld1@_FillValue = -999.
  fld2@_FillValue = -999.
  corr = pattern_cor(fld1,fld2,gw,0)
  rmse = wgt_arearmse(fld1,fld2, gw, 1.0, 0)
  mean1 = wgt_areaave(fld1, gw, 1.0, 0)
  mean2 = wgt_areaave(fld2, gw, 1.0, 0)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do

do i=0,nfiles-1

  fld1 = g->u10(i,:,:)
  fld2 = f[i]->U10m(0,:,:)
  corr = pattern_cor(fld1,fld2,gw,0)
  rmse = wgt_arearmse(fld1,fld2, gw, 1.0, 0)
  mean1 = wgt_areaave(fld1, gw, 1.0, 1)
  mean2 = wgt_areaave(fld2, gw, 1.0, 1)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do


