; create initial datasets from ERA5 ensembles downloaded from Copernicus site.
; Note that the pressure levels in the data is reverse than analysis, that it from 1000
; down to small. Always check the order of pressure levels, as it has changed for
; some reason.
; Also note that the ERA5 ensemble data are on a twice as coarse a grid as the reanalysis data.
; that is 0.5 deg instead of 0.25 deg.

; This version adds the climate perturbation. It is on ERA5 grid, so needs to be degrated
; in resolution.

; Coded by Marat Khairoutdinov, July 2025

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

dataset = "era5_DELTA24-03"
;========================================================================================
; file with climate perturbation (delta):

delta_file = "/glade/work/marat/CLIMO_PERT/ISABEL/TQ_2015-2024--1994-2003_Sep_DELTA.nc"
;delta_file = "/glade/work/marat/CLIMO_PERT/TQ_ERA5_June_2020-2024--1968-1972_DELTA.nc"

;========================================================================================
; ISABEL
date = 2003091600
filedata = "/glade/work/marat/ERA5_ENS/ISABEL/era5_isabel_ens_091600.nc"
identifier = "ISABEL"
; AGNES
;date = 1972061812
;filedata = "/glade/work/marat/ERA5_ENS/AGNES/era5_agnes_ens_061812.nc"
;identifier = "AGNES"
; IRMA
;date = 2017083000
;filedata = "/glade/derecho/scratch/marat/ERA5_ENS/IRMA/era5_irma_ens_0830.nc"
;identifier = "IRMA"
; HELENE
;date = 2024092400
;filedata = "/glade/derecho/scratch/marat/ERA5_ENS/HELENE/0924/era5_helene_ens_0924.nc"
;identifier = "HELENE"
; MILTON
;date = 2024100500
;filedata = "/glade/derecho/scratch/marat/ERA5_ENS/MILTON/1005/data_stream-enda_stepType-instant.nc"
;identifier = "MILTON"

add360 = True ; change Copernicus based data from degrees W (negative nubers) to all-positive degrees E
              ; as the case in NCAR based ERA5 datasets

outdir = "BIN_D"

netcdf_out = True
outdir_nc = "NC_D"

;========================================================================================

f = addfile(filedata,"r")

nlat = dimsizes(f->latitude)
nlon = dimsizes(f->longitude)
nzi = dimsizes(f->pressure_level)

do nens = 0,9

fname = "./"+outdir+"/init_"+dataset+"_"+date+"_"+identifier+"_0"+nens+".bin"
system("rm -f "+fname)
inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)

mmm = 0

; determine pressure levels at z:

nlati = dimsizes(f->latitude)
presi = doubletofloat(f->pressure_level(:))
z3 = f->z(nens,mmm,:,nlati-1:0:1,:)
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
if(add360) then
 lon = lon + 360.
end if
print("lat:"+min(lat)+"  "+max(lat))
print("lon:"+min(lon)+"  "+max(lon))
;print("lon:"+lon)

;========================================================================================
; extract the DELTA perturbation data and regrid to coarser grid:

if(nens.eq.0) then

 f_delta = addfile(delta_file,"r")
 ; prepare delta-fields:
 nlati1 = dimsizes(f_delta->latitude)

 fld = f_delta->z(:,nlati1-1:0:1,:)
 if(add360) then
  fld&longitude=fld&longitude + todouble(360.)
 end if
 z_delta1 = fld
 lat1 = doubletofloat(z_delta1&latitude)
 lon1 = doubletofloat(z_delta1&longitude)
 print("lat DELTA:"+min(lat1)+"  "+max(lat1))
 print("lon DELTA:"+min(lon1)+"  "+max(lon1))
; print("lon1:"+lon1)
 print("z_delta orig: "+min(z_delta1)+"  "+max(z_delta1))
 z_delta = area_hi2lores_Wrap(lon1,lat1,z_delta1,False,1,lon,lat,False)
 print("z_delta regrided: "+min(z_delta)+"  "+max(z_delta))

 fld = f_delta->t(:,nlati1-1:0:1,:)
 if(add360) then
  fld&longitude=fld&longitude + todouble(360.)
 end if
 t_delta1 = fld
 print("t_delta orig: "+min(t_delta1)+"  "+max(t_delta1))
 t_delta = area_hi2lores_Wrap(lon1,lat1,t_delta1,False,1,lon,lat,False)
 print("t_delta regrided: "+min(t_delta)+"  "+max(t_delta))
 do k=0,dimsizes(f_delta->pressure_level)-1
   print(k+"  "+f_delta->pressure_level(k)+" t_delta = "+min(t_delta(k,:,:))+"  "+max(t_delta(k,:,:)))
 end do

 fld = f_delta->q(:,nlati1-1:0:1,:)
 if(add360) then
  fld&longitude=fld&longitude + todouble(360.)
 end if
 q_delta1 = fld
 print("q_delta orig: "+min(q_delta1)+"  "+max(q_delta1))
 q_delta = area_hi2lores_Wrap(lon1,lat1,q_delta1,False,1,lon,lat,False)
 print("q_delta regrided: "+min(q_delta)+"  "+max(q_delta))
 do k=0,dimsizes(f_delta->pressure_level)-1
   print(k+"  "+f_delta->pressure_level(k)+" q_delta = "+min(q_delta(k,:,:))+"  "+max(q_delta(k,:,:)))
 end do

; done

end if ; ens.eq.0

z3 = z3 + z_delta
wgt = cos(z3&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z3,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)

;========================================================================================
print("reading U...")
fld3=f->u(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_U.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->u = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading V...")
fld3=f->v(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+identifier+"_0"+nens+"_V.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->v = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading W...")
fld3=f->w(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_W.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->omega = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading T...")
fld3=f->t(nens,mmm,:,nlati-1:0:1,:)
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
  fbinrecwrite(fname,-1,fld3(k,:,:)+t_delta(k,:,:))
end do
print("done"+" T")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_T.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 fff = fld3+t_delta
 copy_VarAtts(fld3, fff)
 copy_VarCoords(fld3, fff)
 gg->t = fff
 gg->z = z3m
end if
delete(fld3)


print("reading Q...")
fld3=f->q(nens,mmm,:,nlati-1:0:1,:)
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
  fbinrecwrite(fname,-1,fld3(k,:,:)+q_delta(k,:,:))
end do
print("done"+" QV")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_Q.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 fff = fld3+q_delta
 copy_VarAtts(fld3, fff)
 copy_VarCoords(fld3, fff)
 gg->qv = fff
 gg->z = z3m
end if
delete(fld3)

print("reading QC...")
fld3=f->clwc(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_QC.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qc = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QCI..")
fld3=f->ciwc(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_QI.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qi = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QR...")
fld3=f->crwc(nens,mmm,:,nlati-1:0:1,:)
print("qr: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QR")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_QR.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qr = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QS...")
fld3=f->cswc(nens,mmm,:,nlati-1:0:1,:)
print("qs: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QS")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"+DELTA_"+date+"_"+"_0"+nens+"_QS.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qs = fld3
 gg->z = z3m
end if
delete(fld3)


end do ; nens
