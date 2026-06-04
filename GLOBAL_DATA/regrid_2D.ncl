in_lr = addfile("NC_DY/terrain_1152x576_dyvar.nc","r")
in_hr = addfile("BIN_D/terrain_9216x4608_dyvar.nc","r")

file_out_nc = "NC_DY/terrain_9216x4608_dyvar_smooth8.nc"
out_bin = "BIN_D/terrain_9216x4608_dyvar_smooth8.bin"
system("rm -f "+file_out_nc+"  "+out_bin)

fld_lr = in_lr->terr
fld_hr = in_hr->terr

print("fld_lr:"+min(fld_lr)+" "+max(fld_lr))

fld_hr = linint2_Wrap(fld_lr&lon,fld_lr&lat,fld_lr,True,fld_hr&lon,fld_hr&lat,0)

print("fld_hr:"+min(fld_hr)+" "+max(fld_hr))

out_nc = addfile(file_out_nc,"c")
out_nc->terr = fld_hr

fbinrecwrite(out_bin,-1,1)
fbinrecwrite(out_bin,-1,dimsizes(fld_hr&lon))
fbinrecwrite(out_bin,-1,dimsizes(fld_hr&lat))
fbinrecwrite(out_bin,-1,0.)
fbinrecwrite(out_bin,-1,fld_hr)

