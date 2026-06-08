! compute heat index from temperature and relative humidity
! using method used by the US Weather Service

function heat_index(tabs, qv, pres) result (HI)

implicit none
real, intent(in) :: tabs ! temperature in K
real, intent(in) :: qv   ! vapor mixing ration in kg/kg
real, intent(in) :: pres ! pressure in mb
real HI  ! heat index in C 
real(8) T, RH
real, external :: qsatw
! first convert from K to F:

T = (tabs-273.15)*9./5.+32.
RH = qv/qsatw(tabs,pres)*100.
! first use a simple formula:
HI = 0.5 * (T + 61.0 + ((T-68.0)*1.2) + (RH*0.094)) 

! see if the average between HI and T is avove 80F. If true, use a more allborate method.

if(0.5*(HI+T).gt.80.) then

    HI = -42.379 + 2.04901523*T + 10.14333127*RH - .22475541*T*RH &
       - .00683783*T*T - .05481717*RH*RH + .00122874*T*T*RH + &
         .00085282*T*RH*RH - .00000199*T*T*RH*RH
    if(RH.lt.13..and.T.gt.80..and.T.lt.112.) then
        HI = HI - (13.-RH)/4.*SQRT((17-ABS(T-95.))/17.)
    end if
    if(RH.gt.85..and.T.gt.80..and.T.lt.87.) then
        HI = HI + (RH-85.)/10. * (87.-T)/5.
    end if

end if

end function heat_index
