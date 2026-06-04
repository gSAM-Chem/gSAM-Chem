#!/bin/csh

foreach file (*.f90 *.inc)
    sed -i 's/DOUBLE PRECISION/REAL(8)/g' $file
end

