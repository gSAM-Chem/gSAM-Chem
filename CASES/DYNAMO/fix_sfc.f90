character(70) s
read(*,'(a70)') s
print*,s
do while(.true.)
 read*,a,b,c,d,e
 write(*,'(6g15.8)')a,b,c,d,e,0.
end do

end
