
;filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_new_nudge6h1_4608_0000011520"
;filename_obs = "../../../SAM7.3/GLOBAL_DATA/BIN/init_2018021700_9216x4608_dyvarx74_era5.bin"
filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_4608_0000031680"
filename_obs = "../../../SAM7.3/GLOBAL_DATA/BIN/init_2018021818_9216x4608_dyvarx74_era5.bin"
;filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_test_4608_0000000001"
;filename_obs = "../../GLOBAL_DATA/BIN/init_2018021600_9216x4608_dyvarx74_era5.bin"

lowres_nx = 128
lowres_ny = 65

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file

res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
plot = new(3,graphic) 
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frame
res1 = res

g=addfile(filenameroot+"_TERRA.nc","r")

; ------ U850
f=addfile(filenameroot+"_U.nc","r")
pres = f->p
p0=850.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->U(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
terr = g->TERRA(0,ilev,:,:)
terr&lon@units = "degrees east"
terr&lat@units = "degrees north"
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
print("terr_lowre:"+min(terr_lowres)+" "+max(terr_lowres))
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+ilev,(/ny,nx/),"float")
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res) 
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
symMinMaxPlt (dfld,20,False,res) 
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind 850 mb"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down


; ------ V850
f=addfile(filenameroot+"_V.nc","r")
pres = f->p
p0=850.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->V(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+nz+1+ilev,(/ny,nx/),"float")
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+" ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
symMinMaxPlt (dfld,20,False,res) 
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Meridional Wind 850mb"
gsn_panel(wks,plot,(/3,1/),resp)  


; ------ U200
f=addfile(filenameroot+"_U.nc","r")
pres = f->p
p0=200.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->U(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+ilev,(/ny,nx/),"float")
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind 200mb"
gsn_panel(wks,plot,(/3,1/),resp)


; ------ V200
f=addfile(filenameroot+"_V.nc","r")
pres = f->p
p0=200.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->V(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+nz+1+ilev,(/ny,nx/),"float")
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res@cnFillPalette        = "BlWhRe"              ; choose colormap
res@cnLinesOn            = False
symMinMaxPlt (fld,20,False,res)
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res)        ; create plot
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Meridional Wind 200mb"
gsn_panel(wks,plot,(/3,1/),resp)


; ------ T850
f=addfile(filenameroot+"_TABS.nc","r")
pres = f->p
p0=850.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+3*(nz+1)+ilev,(/ny,nx/),"float")
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  200.                 ; set min contour level
res1@cnMaxLevelValF       =  320.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature 850mb"
gsn_panel(wks,plot,(/3,1/),resp)


; ------ T1000
f=addfile(filenameroot+"_TABS.nc","r")
pres = f->p
p0=1000.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+3*(nz+1)+ilev,(/ny,nx/),"float")
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  200.                 ; set min contour level
res1@cnMaxLevelValF       =  320.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)        ; create plot
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Temperature 1000mb"
gsn_panel(wks,plot,(/3,1/),resp)


; ------ T200
f=addfile(filenameroot+"_TABS.nc","r")
pres = f->p
p0=200.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->TABS(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+3*(nz+1)+ilev,(/ny,nx/),"float")
print("tabs_obs:"+fld_obs(ny/2,nx/2)+" tabs:  "+fld(ny/2,nx/2))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
print("tabs_obs_lowres:"+fld_obs_lowres(lowres_ny/2,lowres_nx/2)+" tabs_lowres:  "+fld_lowres(lowres_ny/2,lowres_nx/2))
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  180.                 ; set min contour level
res1@cnMaxLevelValF       =  240.                  ; set max contour level
res1@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1) 
res@cnFillPalette        = "BlWhRe"              ; choose colormap
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "Temperature 200 mb "
gsn_panel(wks,plot,(/3,1/),resp)


; ------ Q850
f=addfile(filenameroot+"_QV.nc","r")
pres = f->p
p0=850.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->QV(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+4*(nz+1)+ilev,(/ny,nx/),"float")
fld_obs = fld_obs*1.e3
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  17.                  ; set max contour level
res1@cnLevelSpacingF      =   1.                  ; set contour spacing
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1) 
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Vapor 850mb"
gsn_panel(wks,plot,(/3,1/),resp)


; ------ Q1000
f=addfile(filenameroot+"_QV.nc","r")
pres = f->p
p0=1000.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->QV(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+4*(nz+1)+ilev,(/ny,nx/),"float")
fld_obs = fld_obs*1.e3
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  17.                  ; set max contour level
res1@cnLevelSpacingF      =   1.                  ; set contour spacing
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1) 
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Vapor 1000mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ QC900
f=addfile(filenameroot+"_QC.nc","r")
pres = f->p
p0=900.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->QC(0,ilev,:,:)
f=addfile(filenameroot+"_QI.nc","r")
fld = fld+f->QI(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
terr_lowres = area_hi2lores_Wrap(xi,yi,terr,True,wy,xo,yo,False)
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+5*(nz+1)+ilev,(/ny,nx/),"float")
fld_obs = fld_obs+fbinrecread(filename_obs,4+6*(nz+1)+ilev,(/ny,nx/),"float")
fld_obs = fld_obs*1.e3
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
fld_lowres = where(terr_lowres.lt.0.1,-999.,fld_lowres)
fld_lowres@_FillValue = -999.
fld_obs_lowres = where(terr_lowres.lt.0.1,-999.,fld_obs_lowres)
fld_obs_lowres@_FillValue = -999.
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  0.1                  ; set max contour level
res1@cnLevelSpacingF      =  0.01                  ; set contour spacing
fld_lowres@long_name="CloudWater/Ice"
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
fld_obs_lowres@long_name="CloudWater/Ice"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Cloud Water/Ice 900mb"
gsn_panel(wks,plot,(/3,1/),resp)

; ------ QC200
f=addfile(filenameroot+"_QC.nc","r")
pres = f->p
p0=200.
ilev = closest_val(p0,pres)
print(p0+"  "+ilev)
fld = f->QC(0,ilev,:,:)
f=addfile(filenameroot+"_QI.nc","r")
fld = fld+f->QI(0,ilev,:,:)
fld&lon@units = "degrees east"
fld&lat@units = "degrees north"
xi = fld&lon
yi = fld&lat
nx = dimsizes(xi)
ny = dimsizes(yi)
nz = dimsizes(pres)
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
fld_obs = fld
fld_obs = fbinrecread(filename_obs,4+5*(nz+1)+ilev,(/ny,nx/),"float")
fld_obs = fld_obs+fbinrecread(filename_obs,4+6*(nz+1)+ilev,(/ny,nx/),"float")
fld_obs = fld_obs*1.e3
print("fld_obs:"+min(fld_obs)+"  "+max(fld_obs))
fld_obs_lowres = area_hi2lores_Wrap(xi,yi,fld_obs,True,wy,xo,yo,False)
dfld = fld_lowres
dfld = fld_lowres - fld_obs_lowres
res1@cnFillPalette        = "BlAqGrYeOrReVi200"              ; choose colormap
res1@cnLinesOn            = False
res1@cnLevelSelectionMode = "ManualLevels"        ; set manual contour levels
res1@cnMinLevelValF       =  0.                 ; set min contour level
res1@cnMaxLevelValF       =  0.1                  ; set max contour level
res1@cnLevelSpacingF      =  0.01                  ; set contour spacing
fld_lowres@long_name="CloudWater/Ice"
fld_lowres@units=fld_lowres@units+" SAM"
plot(0) = gsn_csm_contour_map(wks,fld_lowres,res1)        ; create plot
fld_obs_lowres@units=fld_obs_lowres@units+"ERA5"
fld_obs_lowres@long_name="CloudWater/Ice"
plot(1) = gsn_csm_contour_map(wks,fld_obs_lowres,res1)
symMinMaxPlt (dfld,20,False,res)
dfld@units=dfld@units+" SAM-ERA5"
plot(2) = gsn_csm_contour_map(wks,dfld,res)        ; create plot
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Cloud Water/Ice 200mb"
gsn_panel(wks,plot,(/3,1/),resp)

