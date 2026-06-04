filename = "./BIN_D/out_terr.bin"
filename_out = "./NC_D/out_terr.nc"
nlon = fbinrecread(filename,0,1,"integer")
nlat = fbinrecread(filename,1,1,"integer")
lon = fbinrecread(filename,2,nlon,"float")
lat = fbinrecread(filename,3,nlat,"float")
terr = fbinrecread(filename,4,(/nlat,nlon/),"float")

system("rm -f "+filename_out)
f=addfile(filename_out,"c")
terr!0="lat"
terr!1="lon"
terr&lat = lat
terr&lon = lon
terr&lat@units = "degrees north"
terr&lon@units = "degrees east"
f->terr = terr
terr = fbinrecread(filename,5,(/nlat,nlon/),"float")
f->terr_smooth = terr
terr = fbinrecread(filename,6,(/nlat,nlon/),"float")
f->lake_mask = terr


