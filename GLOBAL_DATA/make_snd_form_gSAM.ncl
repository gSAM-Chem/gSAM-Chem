; make standard SAM sounding from 3D gSAM file at specified point.
f=addfile("../OUT_3D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_TR_2025-09-03-18-00-00_0000010800.3D_atm.nc","r")
g=addfile("../OUT_2D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_TR_2025-09-03-18-00-00_0000010800.2D_atm.nc","r")
lat = 40.48
lon = 286.199

z = f->z
p = f->p
t = f->TABS(0,:,{lat},{lon})
q = f->QV(0,:,{lat},{lon})
u = f->U(0,:,{lat},{lon})
v = f->V(0,:,{lat},{lon})

sst = g->SKT(0,{lat},{lon})
print("sst= "+sst)

time = f->day(0)
pres0 = f->pi(0)
nzi = dimsizes(z)
write_table("snd", "w", [/"  z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]"/], "%s")
write_table("snd", "a", [/time,nzi,pres0/], "%f  %i  %f")
write_table("snd", "a", [/z,p,(t*(1000./p)^(287./1004.)),q,u,v/],"%f  %f  %f  %f  %f  %f")
time =time + 100.
write_table("snd", "a", [/time,nzi,pres0/], "%f  %i  %f")
write_table("snd", "a", [/z,p,(t*(1000./p)^(287./1004.)),q,u,v/],"%f  %f  %f  %f  %f  %f")
print(z+" "+p+" "+(t*(1000./p)^(287./1004.))+" "+q+" "+u+" "+v)



