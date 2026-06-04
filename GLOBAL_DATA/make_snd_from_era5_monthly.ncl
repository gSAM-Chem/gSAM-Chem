; make snd file from monthly era5 data at NCAR


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
year = 2019
month = 12

; chose coordinate of region for acaraging:
lat1 = -9.
lat2 = -4.
lon1 = 298.
lon2 = 304.

dirpath = "/glade/campaign/collections/rda/data/ds633.1/e5.moda.an.pl/"+year+"/e5.moda.an.pl.128_"
endpath = "ll025sc."+year+"010100_"+year+"120100.nc"
endpathuv = "ll025uv."+year+"010100_"+year+"120100.nc"

f_u = addfile(dirpath+"131_u."+endpathuv,"r")
f_v = addfile(dirpath+"132_v."+endpathuv,"r")
f_z = addfile(dirpath+"129_z."+endpath,"r")
f_t = addfile(dirpath+"130_t."+endpath,"r")
f_q = addfile(dirpath+"133_q."+endpath,"r")

filepath = "/glade/campaign/collections/rda/data/ds633.1/e5.moda.an.sfc/"+year+"/e5.moda.an.sfc.128_134_sp.ll025sc."+year+"010100_"+year+"120100.nc"
f_sp = addfile(filepath,"r")

pres = f_u->level

z = dim_avg_Wrap(dim_avg_Wrap(f_z->Z(month-1,:,{lat1:lat2},{lon1:lon2})))
u = dim_avg_Wrap(dim_avg_Wrap(f_u->U(month-1,:,{lat1:lat2},{lon1:lon2})))
v = dim_avg_Wrap(dim_avg_Wrap(f_v->V(month-1,:,{lat1:lat2},{lon1:lon2})))
t = dim_avg_Wrap(dim_avg_Wrap(f_t->T(month-1,:,{lat1:lat2},{lon1:lon2})))
q = dim_avg_Wrap(dim_avg_Wrap(f_q->Q(month-1,:,{lat1:lat2},{lon1:lon2})))
pres0 = dim_avg_Wrap(dim_avg_Wrap(f_sp->SP(month-1,{lat1:lat2},{lon1:lon2})))
pres0 = pres0/100.
q = q*1.e3
t = -t
z = z/9.81

time = 0.
nz = dimsizes(pres)

system("rm -f snd")
header = (/"  z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]"/)
write_table("snd", "a", [/header/], "%s")
write_table("snd", "a", [/time,nz,pres0/], "%f  %i  %f")
write_table("snd", "a", [/z(::-1),pres(::-1),t(::-1),q(::-1),u(::-1),v(::-1)/],"%f  %f  %f  %f  %f  %f")
time = 1000.
write_table("snd", "a", [/time,nz,pres0/], "%f  %i  %f")
write_table("snd", "a", [/z(::-1),pres(::-1),t(::-1),q(::-1),u(::-1),v(::-1)/],"%f  %f  %f  %f  %f  %f")

