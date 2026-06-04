

files = systemfunc("ls -1 ../../OUT_2D/IRMA_2304x1152x74_15s_mom2_dampww6u4_00_2304_20170*0.2D_atm_1440x721.nc")
nfiles = dimsizes(files)
f=addfiles(files,"r")

g = addfile("/glade/work/marat/ERA5/ERA5_2017-05-11Sep_sfc_hourly_unlim.nc","r")
nlat = dimsizes(g->latitude)

latmax = 90.

lat = f[0]->lat({-latmax:latmax})
pi = acos(-1.)
gw = lat
gw = cos(pi/180.*lat)
gw = gw/sum(gw)
gw!0 = "lat"
gw&lat = lat
land = f[0]->LANDMASK(0,:,:)

do i=0,nfiles-1

  fld1 = short2flt(g->t2m(i,nlat-1:0:1,:))
  fld2 = f[i]->T2m(0,:,:)
  fld1 = where(land.lt.0.1,fld1,-999.) 
  fld2 = where(land.lt.0.1,fld2,-999.) 
  fld1@_FillValue = -999.
  fld2@_FillValue = -999.
  corr = pattern_cor(fld1({-latmax:latmax},:),fld2({-latmax:latmax},:),gw,0)
  rmse = wgt_arearmse(fld1({-latmax:latmax},:),fld2({-latmax:latmax},:), gw, 1.0, 0)
  mean1 = wgt_areaave(fld1({-latmax:latmax},:), gw, 1.0, 0)
  mean2 = wgt_areaave(fld2({-latmax:latmax},:), gw, 1.0, 0)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do

do i=0,nfiles-1

  fld1 = short2flt(g->tcwv(i,nlat-1:0:1,:))
  fld2 = f[i]->PW(0,:,:)
  corr = pattern_cor(fld1({-latmax:latmax},:),fld2({-latmax:latmax},:),gw,0)
  rmse = wgt_arearmse(fld1({-latmax:latmax},:),fld2({-latmax:latmax},:), gw, 1.0, 0)
  mean1 = wgt_areaave(fld1({-latmax:latmax},:), gw, 1.0, 0)
  mean2 = wgt_areaave(fld2({-latmax:latmax},:), gw, 1.0, 0)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do

do i=0,nfiles-1

  fld1 = short2flt(g->u10(i,nlat-1:0:1,:))
  fld2 = f[i]->U10m(0,:,:)
  corr = pattern_cor(fld1,fld2,gw,0)
  rmse = wgt_arearmse(fld1,fld2, gw, 1.0, 0)
  mean1 = wgt_areaave(fld1, gw, 1.0, 0)
  mean2 = wgt_areaave(fld2, gw, 1.0, 0)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do

do i=0,nfiles-1

  fld1 = short2flt(g->v10(i,nlat-1:0:1,:))
  fld2 = f[i]->V10m(0,:,:)
  corr = pattern_cor(fld1,fld2,gw,0)
  rmse = wgt_arearmse(fld1,fld2, gw, 1.0, 0)
  mean1 = wgt_areaave(fld1, gw, 1.0, 0)
  mean2 = wgt_areaave(fld2, gw, 1.0, 0)
  print(i+" "+corr+"  "+rmse+"  "+mean1+"  "+mean2)

end do


