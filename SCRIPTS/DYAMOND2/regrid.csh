#! /bin/csh -f

set echo

set filename = ../../OUT_2D/DYAMOND2_9216x4608x74_10s_4608
set filenameM = ../../OUT_MISC/DYAMOND2_9216x4608x74_10s_4608
set filenameL = ../../OUT_2DL/DYAMOND2_9216x4608x74_10s_4608
set filenameI = ../../OUT_INV/DYAMOND2_9216x4608x74_10s_4608_INV.2D_lnd.nc
set outdir = ../../OUT_COARSE
set N1 = 360
set N2 = 354240
set NN = 360
set MM = 0

while ($N1 <= $N2)

 set M = ""
 if($N1 < 10) then
  set M = "000000"
 else if($N1 < 100) then
  set M = "00000"
 else if($N1 < 1000) then
  set M = "0000"
 else if($N1 < 10000) then
  set M = "000"
 else if($N1 < 100000) then
  set M = "00"
 else if($N1 < 1000000) then
  set M = "0"
 endif

 set K = ""
 if($MM < 10) then
  set K = "000"
 else if($MM < 100) then
  set K = "00"
 else if($MM < 1000) then
  set K = "0"
 endif
 rm -f input.nc
 set f = `echo $filename*_000$M$N1.2D_atm.nc`
 ls $f
 ln -s $f input.nc
 ls -l input.nc
 rm -f inputM.nc
 set f = `echo $filenameM*_000$M$N1.2D_atm.nc`
 ls $f
 ln -s $f inputM.nc
 ls -l inputM.nc
 rm inputL.nc
 ln -s ${filenameL}*000$M$N1.2D_lnd.nc inputL.nc
 ls -l inputL.nc
 rm inputI.nc
 ln -s $filenameI inputI.nc
 ls -l inputI.nc
 ncl regrid.ncl
 set f = `echo ${filename}_000$M$N1.2D_atm.1440x720.nc` 
 echo $f
 mv ${outdir}/data.nc $f

 @ N1 = $N1 + $NN
 @ MM = $MM + 1

end
