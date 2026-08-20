cd OUT_3D

for f in goAMAZON_deep_convection_gSAM_wet_season_3D_128_128_256_2014-03-07*.3D; do
    echo "Converting $f"
    ../UTIL/3D2nc "$f"
done
