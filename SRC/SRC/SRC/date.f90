! various date/time related utilities
! MK March 2020

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------

logical function leapyear(year)
implicit none
integer, intent(in) :: year
if(year.eq.0) then
 print*,'year 0 does not exist'
 stop
end if
if(mod(year,4).ne.0) then
 leapyear = .false.
elseif(mod(year,100).ne.0) then
 leapyear = .true.
elseif(mod(year,400).ne.0) then
 leapyear = .false.
else
 leapyear = .true.
end if
end

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------

subroutine date_int2char(date, datechar)
implicit none
integer(8), intent(in) ::  date
character(14), intent(out) :: datechar
write(datechar,'(i14.14)') date
end

!-----------------------------------------------------
!-----------------------------------------------------
!-----------------------------------------------------
! calculte date given reference date and seconds passed since that date
! Date format should be YYYYMMDDHHMMSS

subroutine secs_to_date (date_init, secs_since, date_out)

implicit none

character(14), intent(in) :: date_init  ! reference date
character(14), intent(out) :: date_out  ! output date
real(8), intent(in) :: secs_since ! seconds since

real(8) :: secs_left
integer              :: year, month, day, hour, minute, second
integer              :: days
integer:: days_in_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
logical, external :: leapyear

! Time parsing
read(date_init(1:4),   '(i4)') year
read(date_init(5:6),   '(i2)') month
read(date_init(7:8),   '(i2)') day
read(date_init(9:10),  '(i2)') hour
read(date_init(11:12), '(i2)') minute
read(date_init(13:14), '(i2)') second

! Constructing date format
if(year.eq.0) then
 print*,'year 0 does not exist'
 stop
end if

! compute number of days in secs_since

days = secs_since/86400
secs_left = secs_since - days*86400.

if(leapyear(year)) then
  days_in_month(2) = 29
else
  days_in_month(2) = 28
end if
do while(days.gt.0)
  if(day.eq.days_in_month(month)) then
   day = 1
   if(month.eq.12) then
      month = 1
      year = year+1
      if(leapyear(year)) then
        days_in_month(2) = 29
      else
        days_in_month(2) = 28
      end if
   else
      month = month+1
   end if
  else
      day = day+1
  end if
  days = days-1
end do
hour = secs_left/3600.
secs_left = secs_left - hour*3600.
minute = secs_left/60.
second = secs_left - minute*60.


! Time parsing
write(date_out(1:4),   '(i4.4)') year
write(date_out(5:6),   '(i2.2)') month
write(date_out(7:8),   '(i2.2)') day
write(date_out(9:10),  '(i2.2)') hour
write(date_out(11:12), '(i2.2)') minute
write(date_out(13:14), '(i2.2)') second

end


!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
! calculte fractional day-of-year given date 
! Date format should be YYYYMMDDHHMMSS

subroutine dayofyear_from_date (date, day_of_year)

implicit none

character(14), intent(in) :: date  ! reference date
real(8), intent(out) :: day_of_year 

integer:: days_in_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer              :: year, month, day, hour, minute, second
real(8) secs
logical, external :: leapyear

! Time parsing
read(date(1:4),   '(i4)') year
read(date(5:6),   '(i2)') month
read(date(7:8),   '(i2)') day
read(date(9:10),  '(i2)') hour
read(date(11:12), '(i2)') minute
read(date(13:14), '(i2)') second

if(year.eq.0) then
 print*,'year 0 does not exist'
 stop
end if

! Converting Fortran datetime format into seconds
if(leapyear(year)) then
  days_in_month(2) = 29
else
  days_in_month(2) = 28
end if

secs = 0.
secs = secs + sum(days_in_month(1:month-1))*86400._8
secs = secs + day*86400._8
secs = secs + hour*3600._8
secs = secs + minute*60._8
secs = secs + second

day_of_year = secs/86400._8

end

!------------------------------------------------------
! print date

subroutine print_date (date)

implicit none

character(14), intent(in) :: date  ! reference date

integer              :: year, month, day, hour, minute, second
character(3), parameter :: name(12) =  &
      (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)


! Time parsing
read(date(1:4),   '(i4)') year
read(date(5:6),   '(i2)') month
read(date(7:8),   '(i2)') day
read(date(9:10),  '(i2)') hour
read(date(11:12), '(i2)') minute
read(date(13:14), '(i2)') second

write(*,'(1x,i2.2,a1,i2.2,a1,i2.2,a1,i2,a1,a3,a1,i4.4)')  &
         hour,':',minute,':',second,' ',day,' ',name(month),' ',year

end


!----------------------------------------------------------------------------
! derive date from fractional day of the year and year

subroutine date_from_dayofyear(dayofyear,year,date)
implicit none
real(8), intent(in) :: dayofyear
integer, intent(in) :: year
character(14), intent(out) ::  date
character(14) date_init
real(8) secs 
secs = (dayofyear-1)*86400._8
write(date_init,'(i14.14)') year*10000000000_8+0101000000
call secs_to_date (date_init, secs, date)
end


!----------------------------------------------------------------------------------
! compute number of seconds from given data since Jan 1 1900, 00:00:00
! MK 2023

real(8) function utTimeSeconds(year, month, day, hour, minute, second)
implicit none
integer, intent(in) :: year, month, day, hour, minute, second
real(8) :: julian_date
utTimeSeconds = (julian_date(year,month,day,hour,minute,second) - &
                 julian_date(1900, 1, 1, 0, 0, 0)) * 86400._8
end

real(8) function julian_date(y, m, d, h, min, s)
implicit none
integer, intent(in) :: y, m, d, h, min, s
integer :: JDN, yr, mo
! Use local variables for calculations
yr = y
mo = m
if (mo <= 2) then
   yr = yr - 1
   mo = mo + 12
end if
JDN = d + (2 - (yr / 100) + (yr / 400)) + int(365.25 * yr) + int(30.6001 * (mo + 1))
julian_date = real(JDN,8) + (h / 24.0_8) + (min / 1440.0_8) + (s / 86400.0_8)
end function julian_date


