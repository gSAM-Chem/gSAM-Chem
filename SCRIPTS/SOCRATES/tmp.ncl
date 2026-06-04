
filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_test_4608_0000000001"
filename_obs = "../../GLOBAL_DATA/BIN/init_2018021600_9216x4608_dyvarx74_era5.bin"


g=addfile(filenameroot+"_TERRA.nc","r")

f=addfile(filenameroot+"_TABS.nc","r")
pres = f->p
nx = dimsizes(f->lon)
ny = dimsizes(f->lat)
nz = dimsizes(pres)
i = nx/2
j = ny/2
tabs = f->TABS(0,:,j,i)
f=addfile(filenameroot+"_QV.nc","r")
qv = f->QV(0,:,j,i)
f=addfile(filenameroot+"_QC.nc","r")
qc = f->QC(0,:,j,i)
f=addfile(filenameroot+"_QI.nc","r")
qi = f->QI(0,:,j,i)
tabs_obs = tabs
do ilev=1,nz
 fld = fbinrecread(filename_obs,3+3*(nz+1)+ilev,(/ny,nx/),"float")
 tabs_obs(ilev-1) = fld(j,i) 
end do

print(pres+" "+tabs+" "+tabs_obs)
