f=addfile("../.././OUT_3D/EUROPE_DYAMOND2_2592x1296_1000s_2592_20200121060000_0000010800.atm.3D.nc","r")
i=463
j=1268

u = f->U(0,:,j,i)
v = f->V(0,:,j,i)
t = f->TABS(0,:,j,i)
q = f->QV(0,:,j,i)
z = f->z
p = f->p
t = t*(1000./p)^(287./1004.)

p0 = p(0)-(p(1)-p(0))/(z(1)-z(0))*z(0)
print(0.+"  "+dimsizes(z)+"  "+p0)
print(z+"    "+p+"    "+t+"    "+q+"    "+u+"    "+v)
