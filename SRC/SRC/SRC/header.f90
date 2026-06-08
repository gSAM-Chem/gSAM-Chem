subroutine header()
use grid, only: version, version_date
integer          :: values(8)
character        :: date*8, time*10, zone*5
character(len=8) :: cdate          ! System date
character(len=8) :: ctime          ! System time


call date_and_time (date, time, zone, values)
  cdate(1:2) = date(5:6)
  cdate(3:3) = '/'
  cdate(4:5) = date(7:8)
  cdate(6:6) = '/'
  cdate(7:8) = date(3:4)
  ctime(1:2) = time(1:2)
  ctime(3:3) = ':'
  ctime(4:5) = time(3:4)
  ctime(6:6) = ':'
  ctime(7:8) = time(5:6)


print*,"   _____ _       _           _                                  "
print*,"  / ____| |     | |         | |                                 "
print*," | |  __| | ___ | |__   __ _| |                                 "
print*," | | |_ | |/ _ \| '_ \ / _` | |                                 "
print*," | |__| | | (_) | |_) | (_| | |                                 "
print*,"  \_____|_|\___/|_.__/ \__,_|_|         __                      "
print*,"  / ____|         | |                  / _|                     "
print*," | (___  _   _ ___| |_ ___ _ __ ___   | |_ ___  _ __            "
print*,"  \___ \| | | / __| __/ _ \ '_ ` _ \  |  _/ _ \| '__|           "
print*,"  ____) | |_| \__ \ ||  __/ | | | | | | || (_) | |              "
print*," |_____/ \__, |___/\__\___|_| |_| |_| |_| \___/|_|              "
print*,"          __/ |                                                 "
print*,"         |___/                          _               _       "
print*,"     /\  | |                           | |             (_)      "
print*,"    /  \ | |_ _ __ ___   ___  ___ _ __ | |__   ___ _ __ _  ___  "
print*,"   / /\ \| __| '_ ` _ \ / _ \/ __| '_ \| '_ \ / _ \ '__| |/ __| "
print*,"  / ____ \ |_| | | | | | (_) \__ \ |_) | | | |  __/ |  | | (__  "
print*," /_/  __\_\__|_| |_| |_|\___/|___/ .__/|_| |_|\___|_|  |_|\___| "
print*," |  \/  |         | |    | (_)   | |                            "
print*," | \  / | ___   __| | ___| |_ _ _|_| __ _                       "
print*," | |\/| |/ _ \ / _` |/ _ \ | | '_ \ / _` |                      "
print*," | |  | | (_) | (_| |  __/ | | | | | (_| |                      "
print*," |_|  |_|\___/ \__,_|\___|_|_|_| |_|\__, |                      "
print*,"                                     __/ |                      "
print*,"                                    |___/                       "

write(*,*) '    *****************************************************'
write(*,*) '    *****************************************************'
write(*,*) '    ***                   Global                      ***'
write(*,*) '    ***       System for Atmospheric Modeling         ***'
write(*,*) '    ***                   gSAM                        ***'
write(*,*) '              Version '//version//' ('//version_date//')  '  
write(*,*) '    **************************************************i***'
write(*,*) '    ***     (C) 2016- Marat Khairoutdinov             ***'
write(*,*) '    ***   School of Marine and Atmospheric Sciences   ***'
write(*,*) '    *** Institute for Advanced Computational Science  ***'
write(*,*) '    ***          Stony Brook University               ***'
write(*,*) '    ***          Stony Brook, New York, USA           ***'
write(*,*) '    *****************************************************'
write(*,*) '    *****************************************************'
write(*,*) '           DATE '//cdate//' TIME '//ctime
write(*,*) '    *****************************************************'
write(*,*)

return
end
