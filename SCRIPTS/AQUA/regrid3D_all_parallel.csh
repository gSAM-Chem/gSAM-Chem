#! /bin/csh -f
set filename = ../../OUT_3D/AQUA_2048x512x36_300K_SAM1MOM
set N1 = 259200
set N2 = 3110400
set NN = 4320

set NC  = 1  # number channels
#@ N2 = $N2 - $NN * ($NC - 1)
echo $N2

source ~/.cshrc

while ($N1 < $N2)

  set CH = 0

  while ($CH < $NC)

   @ CH = $CH + 1
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

set f = `echo $filename*_000$M$N1.3D_atm.nc`
set g = `echo $f | sed s/"3D_atm.nc"/".3D_1deg.nc"/`

   if($CH == $NC) then
    ncl filename_in=\"{$f}\"  filename_out=\"{$g}\" regrid3D.ncl  #>& out3D$CH
    wait
   else
    ncl filename_in=\"{$f}\"  filename_out=\"{$g}\" regrid3D.ncl  >& out3D$CH &
   endif
  @ N1 = $N1 + $NN


  end

#  set COUNT=1
#  while ($COUNT != 0)
   sleep 5
#   set COUNT = `ps -ef | grep $filename | wc -l`
#  end


end



