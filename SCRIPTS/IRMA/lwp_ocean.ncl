f=addfile("../../OUT_2D/IRMA_9216x4608x72_10s_u25w30_KKau_qi5_2017-09-06-17-00-00_0000014760.2D_atm.nc","r")
system("rm -f ../../OUT_2D/out_lwp.nc")
g = addfile("../../OUT_2D/out_lwp.nc","c")

landmask = round(f->LANDMASK,3)
print(""+min(landmask)+" "+max(landmask))
cwp = f->CWP
print(""+min(cwp)+" "+max(cwp))
cwp = where(landmask.eq.0,cwp,-999.)
print(""+min(cwp)+" "+max(cwp))
cwp@_FillValue = -999.

g->CWP = cwp
g->LANDMASK = landmask
