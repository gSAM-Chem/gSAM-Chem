
filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_4608_0000011520"

lowres_nx = 128
lowres_ny = 65

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file

res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill

; ------ P
f=addfile(filenameroot+"_PP.nc","r")
pres = f->p
p0=200.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->PP(0,ilev,:,:)
print(min(fld)+"  "+max(fld))
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlon=lowres_nx
dlon = 360. / mlon
xo  = fspan (0, (mlon - 1) * dlon, mlon) 
xo@units = "degrees east"
mlat=lowres_ny
yo  = fspan (-89., 89., mlat) 
yo@units = "degrees north"
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res) 
plot = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot

; ------ P
f=addfile(filenameroot+"_PP.nc","r")
pres = f->p
p0=1000.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->PP(0,ilev,:,:)
print(min(fld)+"  "+max(fld))
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlon=lowres_nx
dlon = 360. / mlon
xo  = fspan (0, (mlon - 1) * dlon, mlon)
xo@units = "degrees east"
mlat=lowres_ny
yo  = fspan (-89., 89., mlat)
yo@units = "degrees north"
fld_lowres = area_hi2lores_Wrap(xi,yi,fld,True,wy,xo,yo,False)
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
plot = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot

