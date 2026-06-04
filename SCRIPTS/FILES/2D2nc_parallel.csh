#! /bin/csh -f

#Set the common file root-name:

set filename = ../../OUT_2D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_TR_TKEFIX
set N1 = 200
set N2 = 24400
set NN = 200
#set filename = ../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_NEWNUDGE0.02_SPNUDGE20_20
#set N1 = 360
#set N2 = 25920
#set NN = 360

# set the utility to execute

set UTILITY = ../../UTIL/2D2nc

# set number of parallel channels. Use batch script (like run.script.2D)
#set NC  = 1  # number channels
set NC  = 24  # number channels (cores)

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


set f = `echo $filename*_000$M$N1.2D`


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



