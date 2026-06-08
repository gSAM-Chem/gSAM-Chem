!=====================================================================
! compute solar unit vector (toward the sun) in model coordinate system,
! that is, x axis in to the east, y axis to the nofth, z axis - up.

subroutine sun_vector(year,day,lon,lat,sun)
! MK 2023

 use mpi_stuff, only: masterproc
 implicit none
 integer, intent(in) :: year   ! 4 digit year
 real(8), intent(in) :: day   ! fractional day of year
 real, intent(in) :: lat ! degrees north
 real, intent(in) :: lon ! degrees east (0 - 360.)
 real, intent(out) :: sun(3) ! unit vector
!local:
 real gmtime ! UTC time in decimal hour
 real xlat, xlon
 integer month, iday
 character(14) date
 real solarz, solara,pi

 call date_from_dayofyear(day,year,date)
 !print*,date
 read(date(5:6),   '(i2)') month
 read(date(7:8),   '(i2)') iday
 !print*,year, month, day
 gmtime = (day-int(day))*24.
 xlat = lat
 if(lon.gt.180.) then
    xlon = lon-360.
 else
    xlon = lon
 end if
 call sunpos_ultimate_azi_atan2(year, month, iday, gmtime, xlat, xlon, solarz, solara)
 pi = acos(-1.)
 sun(1) = -sin(solarz*pi/180.)*sin(solara*pi/180.)
 sun(2) = -sin(solarz*pi/180.)*cos(solara*pi/180.)
 sun(3) = cos(solarz*pi/180.)
 if(masterproc) print*,'sun:',sun


contains
!======================================================================
!======================================================================
!======================================================================
!======================================================================
! 2019-12-26 !Solar Geometry using subsolar point and atan2.
! by Taiping Zhang.
! Input variables:
! inyear: 4-digit year, e.g., 1998, 2020;
! inmon: month, in the range of 1 - 12;
! inday: day, in the range 1 - 28/29/30/31;
! gmtime: GMT in decimal hour, e.g., 15.2167;
! xlat: latitude in decimal degree, positive in Northern Hemisphere;
! xlon: longitude in decimal degree, positive for East longitude.
!
! Output variables:
! solarz: solar zenith angle in deg;
! azi: solar azimuth in deg the range -180 to 180, South-Clockwise
! Convention.
!
! Note: The user may modify the code to output other variables.
Subroutine sunpos_ultimate_azi_atan2(inyear, inmon, inday, gmtime, &
 xlat, xlon, solarz, azi)
 implicit none
 integer:: inyear, inmon, inday, nday(12), julday(0:12), xleap, i, &
 dyear, dayofyr
 real:: gmtime, xlat, xlon
 real:: n, L, g, lambda, epsilon, alpha, delta, R, EoT
 real:: solarz, azi, toadwn
 real, parameter:: rpd=acos(-1.0)/180
 real:: sunlat, sunlon, PHIo, PHIs, LAMo, LAMs, Sx, Sy, Sz
 data nday / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
 if((mod(inyear, 100) .ne. 0 .and. mod(inyear, 4) == 0) .or. &
 (mod(inyear, 100)==0.and.mod(inyear, 400)==0)) then
 nday(2)=29
 else
 nday(2)=28
 endif
 julday(0)=0
 do i=1, 12
 julday(i)=julday(i-1)+nday(i)
 enddo
! Note: julday(12) is equal to either 365 or 366.
 dyear=inyear-2000
 dayofyr=julday(inmon-1)+inday
 if(dyear<=0) then
 xleap=int(real(dyear)/4) !xleap has the SAME SIGN as dyear. !!!
 else
 xleap=int(real(dyear)/4)+1 !Since 2000 is a leap year.
 endif
! --- Astronomical Almanac for the Year 2019, Page C5 ---
 n=-1.5+dyear*365.0+xleap*1.0+dayofyr+gmtime/24
 L=modulo(280.460+0.9856474*n, 360.0)
 g=modulo(357.528+0.9856003*n, 360.0)
 lambda=modulo(L+1.915*sin(g*rpd)+0.020*sin(2*g*rpd), 360.0)
 epsilon=23.439-0.0000004*n
 alpha=modulo(atan2(cos(epsilon*rpd)*sin(lambda*rpd), &
 cos(lambda*rpd))/rpd, 360.0) !alpha in the same quadrant as lambda.
 delta=asin(sin(epsilon*rpd)*sin(lambda*rpd))/rpd
 R=1.00014-0.01671*cos(g*rpd)-0.00014*cos(2*g*rpd)
 EoT=modulo((L-alpha)+180.0, 360.0)-180.0 !In deg.
! --- Solar geometry ---
 sunlat=delta !In deg.
 sunlon=-15.0*(gmtime-12.0+EoT*4/60)
 PHIo=xlat*rpd
 PHIs=sunlat*rpd
 LAMo=xlon*rpd
 LAMs=sunlon*rpd
 Sx=cos(PHIs)*sin(LAMs-LAMo)
 Sy=cos(PHIo)*sin(PHIs)-sin(PHIo)*cos(PHIs)*cos(LAMs-LAMo)
 Sz=sin(PHIo)*sin(PHIs)+cos(PHIo)*cos(PHIs)*cos(LAMs-LAMo)

 solarz=acos(Sz)/rpd !In deg.
 azi=atan2(-Sx, -Sy)/rpd !In deg. South-Clockwise Convention.

end subroutine sunpos_ultimate_azi_atan2

end subroutine sun_vector
