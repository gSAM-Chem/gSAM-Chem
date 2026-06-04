; make initial conditions from gSAM run in a chamber.
; make sure that velocity components at the walls are zro, and 
; also that the dataset domain is larger than actual domain
; as needed for interpolation by gSAM.
; all arrays should be on the same grid, which is a scalar grid
; which means that the velocity components should be interpolated to scalar positions.

filename="../OUT_3D/CHAMBER_288x288x144_7mm_0.05s_dry__0000072000.3D_atm.nc"

dataset = "SAM"

dirout = "BIN_D"

f=addfile(filename,"r")
dims = dimsizes(f->U(0,:,:,:))
nz = dims(0)
ny = dims(1)
nx = dims(2)
print("SAM fields dimensions:"+nx+"x"+ny+"x"+nz)
x = new((/nx+1/),float)
y = new((/ny+1/),float)
z = new((/nz+1/),float)
p = z
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))

fname = dirout+"/init_"+dataset+"_"+nx+"x"+ny+"x"+nz+".bin"
print(""+fname)
system("rm -f "+fname)
inds = (/nx+1,ny+1,nz+1/)
print("inds:"+inds)
fbinrecwrite(fname,-1,inds)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)

print("reading U...")
fld = new((/nz+1,ny+1,nx+1/),float)
fld(0:nz-1,0:ny-1,0:nx-1) = f->U(0,:,:,:)
fld(:,:,nx) = 0.
fld(:,ny,:) = fld(:,ny-1,:)
fld(nz,:,:) = fld(nz-1,:,:)
print(min(fld)+" "+max(fld))
x(0:nx-1) = f->x
dx = x(2)-x(1)
x(nx) = x(nx-1)+dx
y(0:ny-1) = f->y
dy = y(2)-y(1)
y(0:ny-1) = y(0:ny-1)+0.5*dy
y(ny) = y(ny-1)+dy
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
;print("x:"+x)
;print("y:"+y)
;print("z:"+z)
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
delete(fld)
print("done"+" U")

print("reading V...")
fld = new((/nz+1,ny+1,nx+1/),float)
fld(0:nz-1,0:ny-1,0:nx-1) = f->V(0,:,:,:)
fld(:,ny,:) = 0.
fld(:,:,nx) = fld(:,:,nx-1)
fld(nz,:,:) = fld(nz-1,:,:)
print(min(fld)+" "+max(fld))
y(0:nx-1) = f->y
dy = y(2)-y(1)
y(nx) = y(nx-1)+dy
x(0:nx-1) = f->x
dx = x(2)-x(1)
x(0:nx-1) = x(0:nx-1)+0.5*dx
x(nx) = x(nx-1)+dx
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
delete(fld)
print("done"+" V")

print("reading W...")
fld = new((/nz+1,ny+1,nx+1/),float)
fld(0:nz-1,0:ny-1,0:nx-1) = f->W(0,:,:,:)
fld(:,ny,:) = fld(:,ny-1,:)
fld(:,:,nx) = fld(:,:,nx-1)
fld(nz,:,:) = 0.
print(min(fld)+" "+max(fld))
x(0:nx-1) = f->x
dx = x(2)-x(1)
x(0:nx-1) = x(0:nx-1)+0.5*dx
x(nx) = x(nx-1)+dx
y(0:ny-1) = f->y
dy = y(2)-y(1)
y(0:ny-1) = y(0:ny-1)+0.5*dy
y(ny) = y(ny-1)+dy
z(0:nz-1) = f->zi
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->pi
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
delete(fld)
print("done"+" W")

print("reading T...")
fld = new((/nz+1,ny+1,nx+1/),float)
fld(0:nz-1,0:ny-1,0:nx-1) = f->TABS(0,:,:,:)
fld(:,ny,:) = fld(:,ny-1,:)
fld(:,:,nx) = fld(:,:,nx-1)
fld(nz,:,:) = fld(nz-1,:,:)
print(min(fld)+" "+max(fld))
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
print("done"+" T")

print("reading Q...")
fld(0:nz-1,0:ny-1,0:nx-1) = f->QV(0,:,:,:)*1.e-3
fld(:,ny,:) = fld(:,ny-1,:)
fld(:,:,nx) = fld(:,:,nx-1)
fld(nz,:,:) = fld(nz-1,:,:)
print(min(fld)+" "+max(fld))
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
print("done"+" Q")

print("reading QC...")
fld(0:nz-1,0:ny-1,0:nx-1) = f->QC(0,:,:,:)*1.e-3
fld(:,ny,:) = fld(:,ny-1,:)
fld(:,:,nx) = fld(:,:,nx-1)
fld(nz,:,:) = fld(nz-1,:,:)
print(min(fld)+" "+max(fld))
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
print("done"+" QC")

print("reading QI...")
fld(0:nz-1,0:ny-1,0:nx-1) = f->QI(0,:,:,:)*1.e-3
fld(:,ny,:) = fld(:,ny-1,:)
fld(:,:,nx) = fld(:,:,nx-1)
fld(nz,:,:) = fld(nz-1,:,:)
print(min(fld)+" "+max(fld))
z(0:nz-1) = f->z
z(nz) = z(nz-1)+(z(2)-z(1))
p(0:nz-1) = f->p
p(nz) = p(nz-1)+(p(nz-1)-p(nz-2))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
print("done"+" QI")

print("reading QR...")
fld = 0.
print(min(fld)+" "+max(fld))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
print("done"+" QR")

print("reading QS...")
fld = 0.
print(min(fld)+" "+max(fld))
fbinrecwrite(fname,-1,x)
fbinrecwrite(fname,-1,y)
fbinrecwrite(fname,-1,z)
fbinrecwrite(fname,-1,p)
do k=0,nz
  fbinrecwrite(fname,-1,fld(k,:,:))
end do
print("done"+" QS")


