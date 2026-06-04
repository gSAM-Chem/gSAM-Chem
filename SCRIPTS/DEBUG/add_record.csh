foreach f (../../OUT_3D/*TROP*.nc)
 ncks --mk_rec_dmn time $f ../../OUT_3D/out.nc
 mv ../../OUT_3D/out.nc $f
end
