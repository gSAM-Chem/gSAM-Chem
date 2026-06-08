  
subroutine precip_init

! Initialize precipitation related stuff

use vars
use microphysics
use micro_params
use params

implicit none

real pratio, coef1, coef2, estw, esti, rrr1, rrr2
real, external :: esatw,esati
integer k 

if(nint(gamma(3.)).ne.2) then 
   if(masterproc)print*,'cannot compute gamma-function in precip_init. Exiting...'
   call task_abort
end if

do k=1,nzm
	
! pratio = (1000. / pres(k)) ** 0.4
  pratio = sqrt(1.29 / rho(k))	

  rrr1=393./(tabs0(k)+120.)*(tabs0(k)/273.)**1.5
  rrr2=(tabs0(k)/273.)**1.94*(1000./pres(k))

  estw = 100.*esatw(tabs0(k))
  esti = 100.*esati(tabs0(k))

! ice crystal concentration as the function of T (#/L):
! factor of 2 is to account for concentration of crystals smaller than 50 mkm

  niz(k) = max(1.,min(54.,2.*3.3036*exp(-0.046*(tabs0(k)-273.15))))*1000.  

! gamma parameter:

  gami(k) = max(0.,-16.-0.27*(tabs0(k)-273.15))
  qi_min(k) = a_ice*niz(k)*gamma(b_ice+gami(k)+1.)/(rho(k)*gamma(gami(k)+1)) &
             *(gamma(gami(k)+1)/gamma(gami(k)+2)*d_ice_min*1.e-6)**b_ice
       

! ice/snow terminal velocity coefficients:

  bv_ice(k) = 1.17*(1000./pres(k))**0.0783
  av_ice(k) = 0.000603*(pres(k)/1000.)**0.315*10.**(6.*bv_ice(k))

!  if(nstep.eq.0.and.masterproc) then
!      if(k.eq.1) print*,'>>> k   tabs0(C)   niz   gami   qi_min   av_ice   bv_ice vt(50u) <<<<'
!      write(*,'(i4,7g15.5)') k,tabs0(k)-273.,niz(k),gami(k),qi_min(k),av_ice(k),bv_ice(k),  &
!      av_ice(k)*50.e-6**bv_ice(k)
!  end if

! riming of ice/snow :

  accric(k) = 0.25*pi*av_ice(k)*eiccoef*niz(k) &
       *gamma(bv_ice(k)+gami(k)+3.)/gamma(gami(k)+1) &
       *(rho(k)*gamma(gami(k)+1)/(a_ice*niz(k)*gamma(b_ice+gami(k)+1.)))**((2.+bv_ice(k))/b_ice)
          
! deposition/sublimation of ice/snow:
 
  coef1  =(lsub/(tabs0(k)*rv)-1.)*lsub/(therco*rrr1*tabs0(k))
  coef2  = rv*tabs0(k)/(diffelq*rrr2*esti)
  evapi1(k)  =  0.65*4.*niz(k)/rho(k)/(coef1+coef2) &
       *gamma(gami(k)+2)/gamma(gami(k)+1) &
       *(rho(k)*gamma(gami(k)+1)/(a_ice*niz(k)*gamma(b_ice+gami(k)+1.)))**(1./b_ice)
  evapi2(k)  =  0.44*4.*niz(k)/rho(k)/(coef1+coef2)*sqrt(av_ice(k)*rho(k)/(muelq*rrr1)) &
       *gamma(gami(k)+0.5*(5.+bv_ice(k)))/gamma(gami(k)+1) &
       *(rho(k)*gamma(gami(k)+1)/(a_ice*niz(k)*gamma(b_ice+gami(k)+1.)))**(0.5*(3.+bv_ice(k))/b_ice)

! accretion by graupel:

  coef1 = 0.25*pi*nzerog*a_grau*gamma(3.+b_grau)*pratio/&
          (pi*rhog*nzerog/rho(k))**((3+b_grau)/4.)
  accrgc(k) =  coef1 * egccoef 
          
! evaporation of graupel:
 
  coef1  =(lsub/(tabs0(k)*rv)-1.)*lsub/(therco*rrr1*tabs0(k))
  coef2  = rv*tabs0(k)/(diffelq*rrr2*esti)
  evapg1(k)  = 0.78*4.*nzerog/sqrt(pi*rhog*nzerog)/(coef1+coef2)/sqrt(rho(k)) 
  evapg2(k)  = 0.31*4.*nzerog*gamma((5.+b_grau)/2.)*sqrt(a_grau/(muelq*rrr1))/ &
        (pi * rhog * nzerog)**((5+b_grau)/8.) / (coef1+coef2) &
               * rho(k)**((1+b_grau)/8.)*sqrt(pratio)  

 
! accretion by rain:

  accrrc(k)=  0.25 * pi * nzeror * a_rain * gamma(3.+b_rain) * pratio/ &
              (pi * rhor * nzeror / rho(k)) ** ((3+b_rain)/4.)* erccoef   

! evaporation of rain:

  coef1  =(lcond/(tabs0(k)*rv)-1.)*lcond/(therco*rrr1*tabs0(k))
  coef2  = rv*tabs0(k)/(diffelq * rrr2 * estw)
  evapr1(k)  =  0.78 * 2. * pi * nzeror / &
        sqrt(pi * rhor * nzeror) / (coef1+coef2) / sqrt(rho(k)) 
  evapr2(k)  =  0.31 * 2. * pi  * nzeror * gamma((5.+b_rain)/2.) * &
		0.89 * sqrt(a_rain/(muelq*rrr1))/ &
        (pi * rhor * nzeror)**((5+b_rain)/8.) / (coef1+coef2) & 
             * rho(k)**((1+b_rain)/8.)*sqrt(pratio) 

end do

           
end subroutine precip_init


