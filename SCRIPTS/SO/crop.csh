#!/bin/csh

foreach file (plot*.png)
    convert "$file" -crop 500x700+100+0 "${file:r}_cropped.png"
end

