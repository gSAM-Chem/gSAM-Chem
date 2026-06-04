load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

grid = "2304x1152_dyvar"

netcdf_out = True
dataset = "CESM2090s"

f_sst = "BIN_D/sst_2020050100-2020123100_2304x1152_dyvar_era5_7dayrm.bin"
f_ci = "BIN_D/icemask_2020050100-2020123100_2304x1152_dyvar_era5_7dayrm.bin"
f_sst_new = "BIN_D/sst_2020050100-2020123100_2304x1152_dyvar_era5_7dayrm_CESM2090s.bin"
f_ci_new = "BIN_D/icemask_2020050100-2020123100_2304x1152_dyvar_era5_7dayrm_CESM2090s.bin"

f1 = addfile("/glade/work/marat/CESM/RCP8.5/AVG_2000-2009.rcp8.5.TS_30dayavg.nc","r")
f2 = addfile("/glade/work/marat/CESM/RCP8.5/AVG_2090-2099.rcp8.5.TS_30dayavg.nc","r")

fland=addfile("./NC_D/landmask_"+grid+".nc","r")
dsst_grid = fland->LANDMASK
dsst_grid = -999.
ci_grid = dsst_grid

lm = addfile("/glade/work/marat/CESM/landfrac.nc","r")
landfrac = lm->LANDFRAC
landmask = where(landfrac.eq.0.,0.,1.)

fname = f_sst
ntime = fbinrecread(fname,0,1,"integer")
print("ntime="+ntime)
nlon = fbinrecread(fname,1,1,"integer")
nlat = fbinrecread(fname,2,1,"integer")
print("nlon x nlat:"+nlon+" "+nlat)
times = fbinrecread(fname,3,(/ntime/),"float")
print("times:"+times)
rec = 3
fname1 = f_ci
times(0) = times(0)-1.

fname_sst = f_sst_new
system("rm -f "+fname_sst)
fbinrecwrite(fname_sst,0,ntime)
fbinrecwrite(fname_sst,-1,nlon)
fbinrecwrite(fname_sst,-1,nlat)
fbinrecwrite(fname_sst,-1,times(0:ntime-1))

fname_ci = f_ci_new
system("rm -f "+fname_ci)
fbinrecwrite(fname_ci,0,ntime)
fbinrecwrite(fname_ci,-1,nlon)
fbinrecwrite(fname_ci,-1,nlat)
fbinrecwrite(fname_ci,-1,times(0:ntime-1))

do m=0,ntime-1
 rec = rec+1
 sst = fbinrecread(fname,rec,(/nlat,nlon/),"float")
 ci = fbinrecread(fname1,rec,(/nlat,nlon/),"float")
 n = tointeger(times(m))*4-1
 sst1 = f1->TS(n,:,:)
 ci1 = where(landmask.eq.0..and.sst1.lt.271.,1.,0.)
 sst1 = where(landmask.gt.0.or.ci1.gt.0.,271,sst1) 
 sst2 = f2->TS(n,:,:)
 ci2 = where(landmask.eq.0..and.sst2.lt.271.,1.,0.)
 sst2 = where(landmask.gt.0..or.ci2.gt.0.,271,sst2) 
 dsst = sst2 - sst1
 copy_VarCoords(sst1,dsst)
 copy_VarCoords(sst1,ci2)

 dsst = where(landmask.eq.0.,dsst,-999.)
 dsst@_FillValue = -999.
 poisson_grid_fill(dsst, True, 0, 1500, 1.e-2, 0.6, 0) 
 delete(dsst@_FillValue)
 dsst_grid = linint2_Wrap(dsst&lon,dsst&lat,dsst,True,dsst_grid&lon,dsst_grid&lat,0)
 ci_grid = linint2_Wrap(ci2&lon,ci2&lat,ci2,True,ci_grid&lon,ci_grid&lat,0)
 ci_grid = where(ci_grid.lt.0.9,0.,1.)
; print(times(m)+" sst(ERA5): "+min(sst)+" "+max(sst)+" sst1 "+min(sst1)+" "+max(sst1)+" sst2 "+min(sst2)+" "+max(sst2)+" dsst: "+min(dsst)+" "+max(dsst)+" dsst_g "+min(dsst_grid)+" "+max(dsst_grid)) 
 print(times(m)+" dsst_g "+min(dsst_grid)+" "+max(dsst_grid)+" ci_g "+min(ci_grid)+" "+max(ci_grid)) 
 
 fld = sst + dsst_grid
 fbinrecwrite(fname_sst,-1,fld)
 if(netcdf_out) then
   system("rm -f ./NC_D/sst_"+times(m)+"_"+grid+"_"+dataset+".nc")
   out = addfile("./NC_D/sst_"+times(m)+"_"+grid+"_"+dataset+".nc","c")
   fld@long_name = "SST"
   fld@units = "K"
   out->sst = fld
 end if
 fld = ci_grid
 fbinrecwrite(fname_ci,-1,fld)
 if(netcdf_out) then
   system("rm -f ./NC_D/icemask_"+times(m)+"_"+grid+"_"+dataset+".nc")
   out = addfile("./NC_D/icemask_"+times(m)+"_"+grid+"_"+dataset+".nc","c")
   fld@long_name = "Ice Mask"
   out->ci = fld
 end if

end do

