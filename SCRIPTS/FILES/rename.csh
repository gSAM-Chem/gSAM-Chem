#! /bin/csh -f
set echo
foreach f(../../OUT_*/USA_6144x3072x72_2.5km_10s_APR10*)
set g = `echo $f | sed s/"USA_6144x3072x72_2.5km_10s_APR10"/"USA_6144x3072x72_1.25km_10s_APR10"/`
mv $f $g
end
