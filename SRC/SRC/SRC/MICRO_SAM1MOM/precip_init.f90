  
subroutine precip_init

! Initialize precipitation related stuff

use vars
use microphysics
use micro_params
use params

implicit none

real pratio, coef1, coef2,estw,esti,rrr1,rrr2,nzeror1
integer i,j,k 
real, external :: esatw,esati

if(donomicro) return

gam3 = 3. 
gamr1 = 3.+b_rain
gamr2 = (5.+b_rain)/2.
gamr3 = 4.+b_rain
gams1 = 3.+b_snow
gams2 = (5.+b_snow)/2.
gams3 = 4.+b_snow
gamg1 = 3.+b_grau
gamg2 = (5.+b_grau)/2.
gamg3 = 4.+b_grau
gam3 = gamma(gam3) 
gamr1 = gamma(gamr1)
gamr2 = gamma(gamr2)
gamr3 = gamma(gamr3)
gams1 = gamma(gams1)
gams2 = gamma(gams2)
gams3 = gamma(gams3)
gamg1 = gamma(gamg1)
gamg2 = gamma(gamg2)
gamg3 = gamma(gamg3)
!if(masterproc) then
! print*,'gam3=',gam3
! print*,'gamr1,gamr2,gamr3:',gamr1,gamr2,gamr3
! print*,'gams1,gams2,gams3:',gams1,gams2,gams3
! print*,'gamg1,gamg2,gamg3:',gamg1,gamg2,gamg3
!endif
if(nint(gam3).ne.2) then 
   if(masterproc)print*,'cannot compute gamma-function in precip_init. Exiting...'
   call task_abort
end if

do k=1,nzm

! pratio = (1000. / pres(k)) ** 0.4
  pratio = sqrt(1.29 / rho(k))

! accretion by snow:

  coef1 = 0.25 * pi * nzeros * a_snow * gams1 * pratio/ &
            (pi * rhos * nzeros/rho(k) ) ** ((3+b_snow)/4.)
  accrsi(k) =  coef1 * esicoef
  accrsc(k) =  coef1 * esccoef 
          
! accretion by graupel:

  coef1 = 0.25*pi*nzerog*a_grau*gamg1*pratio/&
          (pi*rhog*nzerog/rho(k))**((3+b_grau)/4.)
  accrgi(k) =  coef1 * egicoef
  accrgc(k) =  coef1 * egccoef 
          
! accretion by rain:

  accrrc(k)=  0.25 * pi * nzeror * a_rain * gamr1 * pratio/ &
              (pi * rhor * nzeror / rho(k)) ** ((3+b_rain)/4.)* erccoef   

  do j=1,ny
   do i=1,nx

! adjust intersept parameter of rain - MK 2025
  
    nzeror1 = nzeror * nzeror_factor(i,j,k)

    rrr1=393./(tabs(i,j,k)+120.)*(tabs(i,j,k)/273.)**1.5
    rrr2=(tabs(i,j,k)/273.)**1.94*(1000./pp(i,j,k))
!!    rrr2=(tabs(i,j,k)/273.)**1.94*(1000./pres(k))

    estw = 100.*esatw(tabs(i,j,k))
    esti = 100.*esati(tabs(i,j,k))

! evaporation of snow:
 
    coef1  =(lsub/(tabs(i,j,k)*rv)-1.)*lsub/(therco*rrr1*tabs(i,j,k))
    coef2  = rv*tabs(i,j,k)/(diffelq*rrr2*esti)
    evaps1(i,j,k)  =  0.65*4.*nzeros/sqrt(pi*rhos*nzeros*rho(k))/(coef1+coef2) 
    evaps2(i,j,k)  =  0.49*4.*nzeros*gams2*sqrt(a_snow/(muelq*rrr1))/ & 
       (pi*rhos*nzeros)**((5+b_snow)/8.) / (coef1+coef2) &
               * rho(k)**((1+b_snow)/8.)*sqrt(pratio)  

! evaporation of graupel:
 
    coef1  =(lsub/(tabs(i,j,k)*rv)-1.)*lsub/(therco*rrr1*tabs(i,j,k))
    coef2  = rv*tabs(i,j,k)/(diffelq*rrr2*esti)
    evapg1(i,j,k)  = 0.65*4.*nzerog/sqrt(pi*rhog*nzerog*rho(k))/(coef1+coef2) 
    evapg2(i,j,k)  = 0.49*4.*nzerog*gamg2*sqrt(a_grau/(muelq*rrr1))/ &
        (pi * rhog * nzerog)**((5+b_grau)/8.) / (coef1+coef2) &
               * rho(k)**((1+b_grau)/8.)*sqrt(pratio)  

 
! evaporation of rain:

    coef1  =(lcond/(tabs(i,j,k)*rv)-1.)*lcond/(therco*rrr1*tabs(i,j,k))
    coef2  = rv*tabs(i,j,k)/(diffelq * rrr2 * estw)
    evapr1(i,j,k)  =  0.78 * 2. * pi * nzeror1 / &
        sqrt(pi * rhor * nzeror1*rho(k)) / (coef1+coef2)  
    evapr2(i,j,k)  =  0.31 * 2. * pi  * nzeror1 * gamr2 * &
         0.89 * sqrt(a_rain/(muelq*rrr1))/ &
        (pi * rhor * nzeror1)**((5+b_rain)/8.) / (coef1+coef2) & 
             * rho(k)**((1+b_rain)/8.)*sqrt(pratio) 

  end do  
 end do
end do

           
end subroutine precip_init


