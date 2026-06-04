#! /bin/csh -f
#set filename = ../../OUT_3D/SOCRATES_9216x4608x74_M2005_021600-022100_9216_000
set filename = ../..//OUT_3D/SOCRATES_576x288x74_SAM1MOM_radearth16_nonudge3D_fix_tau10s_576_2018-02-17-00:00:00_000
#set filename = ../..//OUT_3D/SOCRATES_576x288x74_SAM1MOM_radearth16_nonudge3D_fix_576_000
set N1 = 8640
set N2 = 86400
#set N1 = 8640
#set N2 = 46800
set NN = 360
set NC  = 1  # number channels
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

set f = $filename$M${N1}_TABS

   if($CH == $NC) then
     ps -o pid,psr,args
     
     set string = `echo filename=\"$f\"`
     ncl $string avg_tropics_profile.ncl # >& out$CH
#     ../../UTIL/com3D2nc_sepfields_latlon $f >& out$CH
    wait
   else
   endif
  @ N1 = $N1 + $NN


  end

end



