nlats = 128
a = gaus(nlats/2)
asciiwrite("grid_out.txt", str_concat(sprintf ("%12.5f,", a(:,0))))
asciiwrite("grid_out", sprintf ("%12.5f,", a(:,0)))
