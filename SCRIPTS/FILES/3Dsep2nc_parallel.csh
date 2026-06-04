#! /bin/csh -f

# a script to convert many 3D files to netcdf simultaneously - MK

#Set the common file root-name:

set filename = ../../OUT_3D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_20
set N1 = 6000
set N2 = 16000
set NN = 400
#set filename = ../../OUT_3D/HUMBERTO_2048x2048x72_1km_10s_ICEF8_20
#set N1 = 180
#set N2 = 17280
#set NN = 180

# set the utility to execute (only for *.3D_* files, not *3D files)

#set UTILITY = ../../UTIL/3Dsep2nc_sepfields
set UTILITY = ../../UTIL/3Dsep2nc

# set number of parallel channels. Use batch script (like run.script.2D)
# for interactive jobs, set NC to 1
set NC  = 24  # number channels

# Note that the print-out for each of the channels is writen to out2D_* files

#=====================================================================
# don't need to edit below this line

#@ N2 = $N2 - $NN * ($NC - 1)
echo $N2

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


set g =  `echo $filename*_000$M$N1.3D_0`
set f = `echo $g | sed s/"3D_0"/"3D"/`



   if($CH == $NC) then
    $UTILITY  $f #>& out2D$CH
    wait
   else
    $UTILITY  $f >& out2D$CH &
   endif
  @ N1 = $N1 + $NN


  end

#  set COUNT=1
#  while ($COUNT != 0)
   sleep 5
#   set COUNT = `ps -ef | grep $filename | wc -l`
#  end


end



