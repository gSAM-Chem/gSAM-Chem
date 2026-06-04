;              W A R N I N G !;
; This script use the convention that day is year*365+calday of that year.
;
;==========================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

;grid = "768x768_EUSA"
grid = "384x384_EUSA"

; EUSA
date=2024021618

indir = "./NC_D"
outdir = "BIN_D"

netcdf_out = True
outdir_nc = "NC_D"

fobs = addfile("/glade/work/marat/GFS/gfs.0p25."+date+".f000.grib2","r")

;=========================================================

fland=addfile(indir+"/landmask_"+grid+".nc","r")

dataset="GFS"

lat = fland->lat
lon = fland->lon
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("nlat="+nlat)
print("nlon="+nlon)
land = fland->LANDMASK
print("landmask>>>"+min(land)+"  "+max(land))

nnn=dimsizes(fobs->lat_0)-1


fld = new((/1,nlat,nlon/),float)

fname = "./"+outdir+"/snow_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
fbinrecwrite(fname,-1,1)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,0.)

  fld1= fobs->SNOD_P0_L1_GLL0(nnn:0:1,:)
  poisson_grid_fill(fld1, True, 1, 2000, 1.e-2, 0.6, 0)  
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  fld1 = fld1 * 0.1  ; convert to liquid water equivalent depth as GFS has actual depth
  print(date+" original snow depth>>>"+min(fld1)+"  "+max(fld1))
  fld1 = where(fld1.lt.0.,0.,fld1)
  if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
  end if
  fld(0,:,:) = linint2_Wrap(fld1&lon_0,fld1&lat_0,fld1(:,:),True,land&lon,land&lat,0)
  fld(0,:,:) = where(land.ne.1,0.,fld(0,:,:))
  fld!0="time"
  fld&time = 0.
  fld&time@units = "day"
  fld!1="lat"
  fld!2="lon"
  fld&lon = lon
  fld&lat = lat
  delete(fld@_FillValue)
  print("snowd:"+min(fld)+"  "+max(fld))
  fbinrecwrite(fname,-1,fld(0,:,:))
  if(netcdf_out) then
   system("rm -f "+outdir_nc+"/snow_"+date+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdir_nc+"/snow_"+date+"_"+grid+"_"+dataset+".nc","c")
   out->snow = fld
  end if

fname = "./"+outdir+"/snowt_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
fbinrecwrite(fname,-1,1)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,0.)
  fld1= fobs->TMP_P0_L1_GLL0(nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(date+"  original snow temperature>>>"+min(fld1)+"  "+max(fld1))
  if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
  end if
  fld(0,:,:) = linint2_Wrap(fld1&lon_0,fld1&lat_0,fld1(:,:),True,land&lon,land&lat,0)
  fld!0="time"
  fld&time = 0.
  fld&time@units = "day"
  fld!1="lat"
  fld!2="lon"
  fld&lon = lon
  fld&lat = lat
  delete(fld@_FillValue)
  print("sbowt:"+min(fld)+"  "+max(fld))
  fbinrecwrite(fname,-1,fld(0,:,:))
  if(netcdf_out) then
   system("rm -f "+outdir_nc+"/snowt_"+date+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdir_nc+"/snowt_"+date+"_"+grid+"_"+dataset+".nc","c")
   out->snowt = fld
  end if



