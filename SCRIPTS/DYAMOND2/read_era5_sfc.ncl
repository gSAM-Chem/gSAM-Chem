filename1 = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/202001/e5.oper.an.sfc.128_235_skt.ll025sc.2020010100_2020013123.nc"
filename2 = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/202001/e5.oper.an.sfc.128_167_2t.ll025sc.2020010100_2020013123.nc"
filename3 = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/202001/e5.oper.an.sfc.128_139_stl1.ll025sc.2020010100_2020013123.nc"
filename4 = "/glade/collections/rda/data/ds633.0/e5.oper.fc.sfc.accumu/202001/e5.oper.fc.sfc.accumu.128_178_tsr.ll025sc.2020011606_2020020106.nc"
f=addfile(filename1,"r")

ntime = 696
date=f->utc_date(ntime)

fout = "../../OUT_2D/era5_"+date+".nc"
system("rm -f "+fout)
g=addfile(fout,"c")

fld = f->SKT(ntime,:,:)
fld = fld - 273.16
g->SKT = fld

f=addfile(filename2,"r")
fld = f->VAR_2T(ntime,:,:)
g->T2m = fld

f=addfile(filename3,"r")
fld = f->STL1(ntime,:,:)
g->TGR = fld


