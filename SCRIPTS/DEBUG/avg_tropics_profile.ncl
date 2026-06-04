print("filename="+filename)
f=addfile(filename+".nc","r")
tabs = dim_avg_Wrap(dim_avg_Wrap(f->TABS(:,:,{-20.:20.},:)))
print("tabs:"+min(tabs)+" "+max(tabs))
system("rm -f "+filename+"_TROPAVG.nc")
g=addfile(filename+"_TROPAVG.nc","c")
g->TABS = tabs
