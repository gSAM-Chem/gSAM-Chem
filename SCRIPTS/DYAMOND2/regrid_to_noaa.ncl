load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; SAM datset
file_in = "/glade/campaign/univ/unsb0027/SAM7.7_DYAMOND2/OUT_EXTR/SAM_202002_monthly.nc"
file_out = "/glade/campaign/univ/unsb0027/SAM7.7_DYAMOND2/OUT_EXTR/SAM_202002_monthly_144x72.nc"
;-----------------------------------------------------------------------------------------
; ERA5 datset

;-------------------------------------------------------------------
; target resolution

lowres_nx = 144
lowres_ny = 72
latmin=-88.
latmax= 88.

;-----------------------------------------------------------------------------------------
; prepare initial settings for plots:

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
plot = new(3,graphic) 
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frame
res@cnLevelSelectionMode = "ManualLevels" 
res1 = res
;-----------------------------------------------------------------------------------------

f=addfile(file_in,"r")
system ("rm -f "+file_out)
g=addfile(file_out,"c")

vars = (/"LWNT","SWNT","Preca","LHF","SHF","LWNS","SWNS","T2m"/)

do ivar=0,dimsizes(vars)-1

fld = f->$vars(ivar)$(0,:,:)
print(""+fld@long_name)
fld&lon@units = "degrees_east"
fld&lat@units = "degrees_north"
xi = fld&lon
yi = fld&lat
pi = acos(-1.)
wy = cos(yi*pi/180.)
mlon=lowres_nx
dlon = 360. / mlon
xo  = fspan (0, (mlon - 1) * dlon, mlon)
xo@units = "degrees east"
mlat=lowres_ny
yo  = fspan (latmin,latmax, mlat)
yo!0 = "lat"
yo&lat = yo
yo@units = "degrees north"
print("high res: "+min(fld)+"  "+max(fld))
fld1 = fld({-90.:-65.},:)
fld_lowres_poles1 = linint2_Wrap(xi,fld1&lat,fld1,True,xo,yo({-90.:-65}),0)
print("low res: "+min(fld_lowres_poles1)+"  "+max(fld_lowres_poles1))
fld2 = fld({65.:90.},:)
fld_lowres_poles2 = linint2_Wrap(xi,fld2&lat,fld2,True,xo,yo({65.:90}),0)
print("low res: "+min(fld_lowres_poles2)+"  "+max(fld_lowres_poles2))
fld3 = fld({-65.:65.},:)
wy3 = cos(fld3&lat*pi/180.)
fld_lowres_in = area_hi2lores_Wrap(xi,fld3&lat,fld3,True,wy3,xo,yo({-65.:65}),False)
print("low res: "+min(fld_lowres_in)+"  "+max(fld_lowres_in))
;print(fld_lowres_poles1&LAT)
;print(fld_lowres_poles2&LAT)
;rint(fld_lowres_in&LAT)
ny1=dimsizes(fld_lowres_poles1&LAT)
ny2=dimsizes(fld_lowres_poles2&LAT)
ny3=dimsizes(fld_lowres_in&LAT)
print(ny1+" "+ny2+" "+ny3)
print(""+(ny1+ny2+ny3))
field = new((/lowres_ny,lowres_nx/),float)
lat = new((/lowres_ny/),float)
lon = new((/lowres_nx/),float)
field(0:ny1-1,:) = fld_lowres_poles1
field(ny1:ny1+ny3-1,:) = fld_lowres_in
field(ny1+ny3:ny1+ny3+ny2-1,:) = fld_lowres_poles2
lat(0:ny1-1) = fld_lowres_poles1&LAT
lat(ny1:ny1+ny3-1) = fld_lowres_in&LAT
lat(ny1+ny3:ny1+ny3+ny2-1) = fld_lowres_poles2&LAT
;print(field(:,0))
lon = fld_lowres_poles2&LON
field!0="lat"
field!1="lon"
field&lat = lat
field&lon = lon
field&lat@units = "degrees north"
field&lon@units = "degrees east"
g->$vars(ivar)$ = field

end do

