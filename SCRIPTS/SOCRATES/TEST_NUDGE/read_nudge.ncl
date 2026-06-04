filenames = asciiread("filenames.txt",241,"string")
dir = "../../../GLOBAL_DATA/NUDGE/"
ntimes = dimsizes(filenames)
print("ntime="+ntimes)
do i=0,ntimes-1
 in = 0
 dims = fbinrecread(dir+filenames(i),in,3,"integer")
 print("dims:"+dims)
 nx = dims(0)
 ny = dims(1)
 nz = dims(2)
 if(i.eq.0) then
   fld = new((/nz,ny,nx/),float)
   time = new((/ntimes/),float)
   field = new((/ntimes,ny,nx/),float)
 end if
 in=in+1
 x = fbinrecread(dir+filenames(i),in,nx,"float")  
 in=in+1
 y = fbinrecread(dir+filenames(i),in,ny,"float")  
 in=in+1
 z = fbinrecread(dir+filenames(i),in,nz,"float")  
 in=in+1
 p = fbinrecread(dir+filenames(i),in,nz,"float")  
 time(i) = i
 do k=0,nz-1
  in=in+1
  fld(k,:,:) = fbinrecread(dir+filenames(i),in,(/ny,nx/),"float")
 end do
 in=in+1
 x = fbinrecread(dir+filenames(i),in,nx,"float")  
 in=in+1
 y = fbinrecread(dir+filenames(i),in,ny,"float")  
 in=in+1
 z = fbinrecread(dir+filenames(i),in,nz,"float")  
 in=in+1
 p = fbinrecread(dir+filenames(i),in,nz,"float")  
 do k=0,nz-1
  in=in+1
  fld(k,:,:) = fbinrecread(dir+filenames(i),in,(/ny,nx/),"float")
 end do
 in=in+1
 x = fbinrecread(dir+filenames(i),in,nx,"float")  
 in=in+1
 y = fbinrecread(dir+filenames(i),in,ny,"float")  
 in=in+1
 z = fbinrecread(dir+filenames(i),in,nz,"float")  
 in=in+1
 p = fbinrecread(dir+filenames(i),in,nz,"float")  
 do k=0,nz-1
  in=in+1
  fld(k,:,:) = fbinrecread(dir+filenames(i),in,(/ny,nx/),"float")
 end do
 field(i,:,:) = fld(1,:,:)
 print(i+" "+min(field(i,:,:))+" "+max(field(i,:,:)))
end do

system("rm ../../../OUT_MOMENTS/t_nudge.nc")
g=addfile("../../../OUT_MOMENTS/t_nudge.nc","c")

field!0="time"
field!1="lat"
field!2="lon"
field&time = time
field&lat = y
field&lon = x
g->T = field
