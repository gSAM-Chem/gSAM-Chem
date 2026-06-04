#!/bin/csh

set echo

set N = 0
set NMAX = 70

while ($N <= $NMAX)

echo $N > number

if($N < 10) then
  set M = "000"
else if($N < 100) then
  set M = "00"
else if($N < 1000) then
  set M = "0"
else
  set M = ""
endif


ncl  plot_frame_pw.ncl
ctrans -d sun -res 512x512 gsnapp.ncgm > out.sun
convert out.sun FRAMES/image$M$N.jpg

@ N = $N + 1

end
