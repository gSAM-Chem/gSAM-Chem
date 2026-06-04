load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
; IAN
;date = 2022092300
;dates_obs = "2022092300_2022092323"
;month = "202209"
; IRMA
;date = 2017090200
;dates_obs = "2017090200_2017090223"
;month = "201709"
;date = 2020052000
;dates_obs = "2020052000_2020052023"
;month = "202005"
;date = 2020050100
;dates_obs = "2020050100_2020050123"
;month = "202005"
; CA 
date = 2021102112
dates_obs = "2021102100_2021102123"
month = "202110"


nodyn = True ; don't add perturbation to wind


outdir = "BIN_D"

netcdf_out = True
outdir_nc = "NC_DY"

; difference between future and present climates from CESM:
cc_file = "/glade/work/marat/CESM/RCP8.5/DIFF_2090s-2010s_10.nc_era5.nc"
;cc_file = "/glade/work/marat/CESM/RCP8.5/DIFF_2090s-2000s.rcp8_5.hyb.05_era5.nc"
f_cc = addfile(cc_file,"r")

dataset = "era5_CESM2090s"
dirpath = "/glade/campaign/collections/rda/data/ds633.0/e5.oper.an.pl/"+month+"/"
f_u = addfile(dirpath+"e5.oper.an.pl.128_131_u.ll025uv."+dates_obs+".nc","r")
f_v = addfile(dirpath+"e5.oper.an.pl.128_132_v.ll025uv."+dates_obs+".nc","r")
f_w = addfile(dirpath+"e5.oper.an.pl.128_135_w.ll025sc."+dates_obs+".nc","r")
f_z = addfile(dirpath+"e5.oper.an.pl.128_129_z.ll025sc."+dates_obs+".nc","r")
f_t = addfile(dirpath+"e5.oper.an.pl.128_130_t.ll025sc."+dates_obs+".nc","r")
f_q = addfile(dirpath+"e5.oper.an.pl.128_133_q.ll025sc."+dates_obs+".nc","r")
f_clwc = addfile(dirpath+"e5.oper.an.pl.128_246_clwc.ll025sc."+dates_obs+".nc","r")
f_ciwc = addfile(dirpath+"e5.oper.an.pl.128_247_ciwc.ll025sc."+dates_obs+".nc","r")

times = calendar_decode2(f_u->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(dates)
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
time = times(itime,0)*365+day_of_year(times(itime,0),times(itime,1),times(itime,2))+times(itime,3)/24.
print("  time="+time)

nlat = dimsizes(f_u->latitude)
nlon = dimsizes(f_u->longitude)
nzi = dimsizes(f_u->level)
grid = nlon+"x"+nlat
print("data grid:"+grid+"x"+nzi)

fname = "./"+outdir+"/init_"+dataset+"_"+date+"_"+grid+"x"+nzi+".bin"
system("rm -f "+fname)
inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)

mmm = itime

; determine pressure levels at z:

nlati = dimsizes(f_u->latitude)
presi = doubletofloat(f_u->level(nzi-1:0:1))
z3 = f_z->Z(mmm,nzi-1:0:1,nlati-1:0:1,:)
diff = f_cc->Z3(nzi-1:0:1,:,:)
z3 = z3 + diff
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
wgt = cos(z3&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z3,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)

print("reading U...")
fld3=f_u->U(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(.not.nodyn) then
 diff = f_cc->U(nzi-1:0:1,:,:)
 fld3 = fld3 + diff
end if
print("u: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" U")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_U.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->u = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading V...")
fld3=f_v->V(mmm,nzi-1:0:1,nlati-1:0:1,:)
diff = f_cc->V(nzi-1:0:1,:,:)
if(.not.nodyn) then
 diff = f_cc->V(nzi-1:0:1,:,:)
 fld3 = fld3 + diff
end if
print("v: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" V")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_V.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->v = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading W...")
fld3=f_w->W(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("w: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" OMEGA")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_W.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->omega = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading T...")
fld3=f_t->T(mmm,nzi-1:0:1,nlati-1:0:1,:)
diff = f_cc->T(nzi-1:0:1,:,:)
fld3 = fld3 + diff
print("t: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" T")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_T.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->t = fld3
 gg->z = z3m
end if
delete(fld3)


print("reading Q...")
fld3=f_q->Q(mmm,nzi-1:0:1,nlati-1:0:1,:)
diff = f_cc->Q(nzi-1:0:1,:,:)
fld3 = fld3 + diff
print("q: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QV")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_Q.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qv = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QC...")
fld3=f_clwc->CLWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("qc: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QC")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_QC.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qc = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QCI..")
fld3=f_ciwc->CIWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("qi: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QI")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"x"+nzi+"_QI.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qi = fld3
 gg->z = z3m
end if

fld3 = 0.
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QR")

fld3 = 0.
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QS")
delete(fld3)

