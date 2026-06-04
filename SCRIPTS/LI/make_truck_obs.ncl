obsfile = "20250903_Time_lat_lon.csv"
  ncol    = 3
  obs     = readAsciiTable(obsfile, ncol, "double", 0)   ; [nrec x 3]
  t_hrs   = obs(:,0)
  lat_o   = todouble(obs(:,1))
  lon_o   = todouble(obs(:,2) + 360.)
  nrec    = dimsizes(t_hrs)

  lonmin = min(lon_o)-0.1
  lonmax = max(lon_o)+0.1
  latmin = min(lat_o)-0.1
  latmax = max(lat_o)+0.1

  z_max = 3000.

  files = systemfunc("ls -1 ../../OUT_3D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_2025-09*.3D_atm.nc")
  nfiles = dimsizes(files)
  print(""+nfiles)
  f = addfiles(files,"r")
  files1 = systemfunc("ls -1 ../../OUT_3D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_2025-09*.3D_atm.nc")
  nfiles = dimsizes(files1)
  print(""+nfiles)
  g = addfiles(files1,"r")

  time0 = 14. ; staring file
  dtime = 1./3. ; tie period between the files (in hours)
  time_max = 22.

  nrec_max = 0
  do i=1,nrec-1
     if(t_hrs(i).gt.time_max) then
       break
     end if
     nrec_max = nrec_max +1
  end do
  print("nrec_max = "+nrec_max)


  z = f[0]->z({0.:z_max})
  nz = dimsizes(z)
  print("nz ="+nz)

  www = new((/nz,nrec_max/),float)
  qqq = new((/nz,nrec_max/),float)

  k0 = -1
  do i=0, nrec_max-1

     k = toint(floor((t_hrs(i)-time0)/dtime))
     print(k+"   "+t_hrs(i)+"  "+lat_o(i)+"  "+lon_o(i))
     if(k.ne.k0) then
       w = f[k]->W(0,{0.:z_max},{latmin:latmax},{lonmin:lonmax})
       qv = g[k]->QV(0,{0.:z_max},{latmin:latmax},{lonmin:lonmax})
     end if
     k0 = k
     www(:,i) = w(:,{lat_o(i)},{lon_o(i)})
     qqq(:,i) = qv(:,{lat_o(i)},{lon_o(i)})
     print(z+"  "+www(:,i)+"   "+qqq(:,i))
     

  end do

  www!1 = "time"
  www!0 = "z"
  www&time = t_hrs(0:nrec_max-1)
  www&z = z(0:nz-1)
  qqq!1 = "time"
  qqq!0 = "z"
  qqq&time = t_hrs(0:nrec_max-1)
  qqq&z = z(0:nz-1)

  system("rm -f out.nc")
  gg = addfile("out.nc","c")
  gg->W = www
  gg->QV = qqq

