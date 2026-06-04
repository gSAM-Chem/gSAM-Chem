load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

grid="3072x3072_terr4km_CA"

indir = "./NC_D"
outdir = "./BIN_D"

f=addfile(indir+"/terrain_"+grid+".nc","r")
fld = f->terr_smooth
lat = f->lat
lon = f->lon
delete(fld&lon)
delete(fld&lat)
fld!0="lat"
fld!1="lon"
fld&lat = lat
fld&lon = lon
nlon=dimsizes(lon)
nlat=dimsizes(lat)
print(nlon+"  "+nlat)
print("terr:"+min(fld)+"  "+max(fld))
system("rm -f "+outdir+"/terrain_"+grid+".bin")
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,1)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,nlon)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,nlat)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,0.)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,fld)

delete(fld)

