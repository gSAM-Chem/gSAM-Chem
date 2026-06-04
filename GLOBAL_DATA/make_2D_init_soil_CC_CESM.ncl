load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
grid = "2304x1152_dyvar"

netcdf_out = True

indir = "NC_D"
outdir = "BIN_D"

dataset= "era5_CESM2090s"

f_grid = "BIN_D/soil_init_2020050100_2304x1152_dyvar_era5.bin"
f_grid_new = "BIN_D/soil_init_2020050100_2304x1152_dyvar_era5_CESM2090s.bin"
system("rm -f f_grid_new")

fland=addfile("./"+indir+"/landmask_"+grid+".nc","r")

f1 = addfile("/glade/work/marat/CESM/RCP8.5/AVG_2000-2009.rcp8.5.TS_30dayavg.nc","r")
f2 = addfile("/glade/work/marat/CESM/RCP8.5/AVG_2090-2099.rcp8.5.TS_30dayavg.nc","r")

nday = 122

fland=addfile("./NC_D/landmask_"+grid+".nc","r")
dts_grid = fland->LANDMASK
dts_grid = -999.

lm = addfile("/glade/work/marat/CESM/landfrac.nc","r")
landfrac = lm->LANDFRAC
landmask = where(landfrac.eq.0.,0.,1.)

ts1 = f1->TS(nday*4-1,:,:)
ts1 = where(landmask.gt.0.,ts1,-999.)
ts1@_FillValue = -999.
poisson_grid_fill(ts1, True, 0, 1500, 1.e-2, 0.6, 0)
delete(ts1@_FillValue)
ts2 = f2->TS(nday*4-1,:,:)
ts2 = where(landmask.gt.0.,ts2,-999.)
ts2@_FillValue = -999.
poisson_grid_fill(ts2, True, 0, 1500, 1.e-2, 0.6, 0)
delete(ts2@_FillValue)
dts = ts2
dts = dts - ts1
dts_grid = linint2_Wrap(dts&lon,dts&lat,dts,True,dts_grid&lon,dts_grid&lat,0)

if(netcdf_out) then
  fld = dts_grid
  system("rm -f ./NC_D/dts_"+grid+"_"+dataset+".nc")
  out = addfile("./NC_D/dts_"+grid+"_"+dataset+".nc","c")
  fld@long_name = "DTS"
  fld@units = "K"
  out->dts_grid = fld
end if

nsoil = fbinrecread(f_grid,0,1,"integer")
nlon = fbinrecread(f_grid,1,1,"integer")
nlat = fbinrecread(f_grid,2,1,"integer")
zsoil = new((/nsoil/),float)
zsoil = fbinrecread(f_grid,3,(/nsoil/),"float")

fbinrecwrite(f_grid_new,0,nsoil)
fbinrecwrite(f_grid_new,1,nlon)
fbinrecwrite(f_grid_new,2,nlat)
fbinrecwrite(f_grid_new,3,zsoil)
print("nlon x nlat "+nlon+" "+nlat)
print("nsoil "+nsoil)
print("zsoil "+zsoil)
rec = 3
do i=1,nsoil
  rec = rec+1
  fld = fbinrecread(f_grid,rec,(/nlat,nlon/),"float")
  fld = fld + dts_grid
  print("tsoil:"+" "+i+" "+min(fld)+" "+max(fld))
  fbinrecwrite(f_grid_new,rec,fld)
end do
do i=1,nsoil
  rec = rec+1
  fld = fbinrecread(f_grid,rec,(/nlat,nlon/),"float")
  print("swsoil:"+" "+i+" "+min(fld)+" "+max(fld))
  fbinrecwrite(f_grid_new,rec,fld)
end do



