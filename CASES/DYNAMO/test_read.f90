open(1,file="snd",form="formatted")
read(1,*)
do while (.true.)
 read(1,*) day,n,pres
 do k=1,n
  read(1,*) 
 end do
 print*,day
end do

end
