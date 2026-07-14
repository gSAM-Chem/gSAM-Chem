module het_chem

   use grid, only: nx, ny, nzm, nz, dimx1_s, dimx2_s, dimy1_s, dimy2_s
   use vars, only: rho, dtn, tabs0, pres, qcl

   use chem_aqueous, only: naqchem_fields, molwt, iIEPOX, iTETROL, iIEPOX_SO4, &
                           iIP1NIT, iIPDINIT, aq_species_names, aq_gasprod_species_names, &
                           flag_aqchemvar_out3D, flag_aqchemgasvar_out3D, &
                           iepox_aqueous_tendencies, isop1nit_aqueous_tendencies

   use chem_aerosol, only: iTETROLr, iIEPOX_SO4r, iepox_aero_transfer_rate, narchem_fields, molwt_ar, ar_species_names
   use microphysics, only: micro_field, incl     ! , iqad, inad
   use micro_params, only : dopredictNc, Nc0
   use chemistry_params, only: p0, rhol, rho_aerosol, sigma_accum, do_iepox_droplet_chem, do_iepox_aero_chem
   use cloudchem_Parameters, only: ind_IEPOX, ind_ISOP1Nit, ind_ISOPDiNit, NVAR

   implicit none

   real, parameter :: pi = 3.1415927

   real, allocatable, dimension(:) :: Haq, NO3aq, SO4aq, HSO4aq        ! Constant (for now) aqueous concentrations
   real, allocatable, dimension(:) :: Haero, SO4aero, HSO4aero         ! Constant (for now) aerosol concentrations
   real actHaero                                                      ! H activity in aerosols
   real OrgMF                                                          ! Organic mass fraction of aerosol
   real FracTETROL                                                     ! Fraction of IEPOXg to convert to TETROL on aerosol
   real FracIEPOX_SO4                                                  ! Fraction of IEPOXg to convert to SO4 on aerosol

   logical :: isallocatedHetChem = .false.                             ! Is heterogeneous chemistry allocated already?

CONTAINS

   subroutine het_chem_initialize(hi_org, pHdrop, pHaero)
      implicit none

      logical, intent(in) :: hi_org               ! High organic switch
      real, intent(in) :: pHdrop, pHaero          ! pH of droplets and aerosols

      if (.not. isallocatedHetChem) then
         allocate (Haq(nzm), NO3aq(nzm), SO4aq(nzm), HSO4aq(nzm))
         allocate (Haero(nzm), SO4aero(nzm), HSO4aero(nzm))
         isallocatedHetChem = .true.
      end if

      ! For now set these aqueous concentrations as constants, could be put into namelist input in future
      Haq = 10**(-pHdrop)             ! Haq = 1.e-5 ! M H+ (Pye et al)
      NO3aq = 1.e-7                   ! M NO3-
      SO4aq = 1.e-7                   ! M SO4-2
      HSO4aq = 1.e-10                 ! M HSO4-

      ! Aerosol concentrations

      Haero = 10**(-pHaero)           ! Haero = 0.000038 ! M of H+ in aerosols based on RH = 0.8
      actHaero = 1.39                 !  based on RH=0.8
      SO4aero = 3.8e-5                ! M of nucleophile SO4-2
      HSO4aero = 6.3e-5               ! M of general acid HSO4-

      if (hi_org) then
         OrgMF = 0.85
      else
         OrgMF = 0.2
      end if

      FracTETROL = 0.85           ! portion of IEPOXg to convert to TETROL on aerosol
      FracIEPOX_SO4 = 0.15        ! portion of IEPOXg to convert to SO4 on aerosol
   end subroutine het_chem_initialize

    subroutine het_chem_driver(gchem_field, aqchem_field, aqchem_gasprod_field, aqchem_horiz_mean_tend, aqchem_gasprod_horiz_mean_tend, archem_field, archem_horiz_mean_tend)
      real, intent(inout) :: gchem_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, NVAR)                          ! in ppv air
      real, intent(inout) :: aqchem_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, naqchem_fields)               ! in kg/kg
      real, intent(inout) :: aqchem_gasprod_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, naqchem_fields)       ! in ppv air
      real, intent(inout) :: archem_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, narchem_fields)               ! in kg/kg

      real, intent(inout) :: aqchem_horiz_mean_tend(nzm, naqchem_fields)
      real, intent(inout) :: aqchem_gasprod_horiz_mean_tend(nzm, naqchem_fields)
      real, intent(inout) :: archem_horiz_mean_tend(nzm, narchem_fields)

      real :: aq_adjusted_tendency(naqchem_fields)
      real, parameter :: min_aerosol_radius = 1.e-8               ! Prevents division by zero

      real :: aq_tend(nzm, naqchem_fields)             ! mean tendency profiles of aq species
      real :: aq_gasprod_tend(nzm, naqchem_fields)     ! mean tendency profiles of gaseous products
      real :: ar_tend(nzm, narchem_fields)             ! mean tendency profiles of aerosol chemistry products

      real :: aqchem_conc(nzm, naqchem_fields)         ! In mol L-1 (M)
      real :: aqgas_conc(nzm, naqchem_fields)
      real :: archem_conc(nzm, narchem_fields)         ! In mol L-1 (M)

      real :: rho_tot_aerosol(nzm)        ! mean aerosol density (org+inorg component)
      real :: aero_transfer_rate(nzm)     !  fractional rate of transfer of IEPOXg to aerosol surface (/s)
      real :: aero_radius(nzm)            ! interstitial accumulation mode area weighted radius
      real :: rho_org_aerosol(nzm)
      real :: num_conc(nzm)
      real :: Rdrop(nzm)                  ! cloud droplet radius in m

      ! Temporary variables
      real :: qcloud(nzm)               ! temporary cloud water array
      real :: water_vol_frac(nzm)       ! temporary water volume array
      real :: pressure_atm(nzm)         ! temporary pressure array [atm]

      logical :: override_gamma = .true.
      logical :: do_debug_output = .false.
      real :: cloud_threshold = 5e-5   ! in kg/kg

      real :: IEPOX_transfer_ppv
      integer :: i, j, k, ispecies

      ! compute aqueous chemistry tendencies and apply them to aqueous fields
      aqchem_horiz_mean_tend(:, :) = 0.
      aqchem_gasprod_horiz_mean_tend(:, :) = 0.

      ! copy IEPOX into IEPOXg before calling aqueous
      aqchem_gasprod_field(:, :, :, iIEPOX) = gchem_field(:, :, :, ind_IEPOX)       ! in ppv air
      aqchem_gasprod_field(:, :, :, iIP1NIT) = gchem_field(:, :, :, ind_ISOP1Nit)   ! in ppv air
      aqchem_gasprod_field(:, :, :, iIPDINIT) = gchem_field(:, :, :, ind_ISOPDiNit) ! in ppv air

      pressure_atm = pres/p0                ! pres is in mbar; p0 is in mbar --> cancels out to atm

      do j = 1, ny
         do i = 1, nx
            qcloud = qcl(i,j,:)                        ! Cloud mixing ratio [kg/kg], used to be micro_field(i,j,:,iqcl)
            water_vol_frac = qcloud * (rho / rhol)         ! Converts kg water / kg air to m3 water / m3 air (rho, rhol in kg/m3)

            aq_tend(:, :) = 0.
            aq_gasprod_tend(:, :) = 0.

            ! Compute droplet radius - for now assume monodisperse
            do k = 1, nzm
               ! Number concentration is 1, or number concentration of liquid cloud droplets\

               if ( dopredictNc ) then
                  num_conc(k) = micro_field(i, j, k, incl)! max(1., micro_field(i, j, k, incl))                       ! Talk to Peter; number concentration is in number / kg
               else 
                  num_conc(k) = Nc0 * 100**3 / rho(k)      ! #/cm3 --> #/kg air
               endif

               if ( qcloud(k) .gt. cloud_threshold ) then
                  Rdrop(k) = ( (1.0 / rhol)*(qcloud(k) / num_conc(k))*(0.75/pi) )**(1./3.)   ! Calculates radius

                  if (Rdrop(k) .lt. min_aerosol_radius) then
                     Rdrop(k) = min_aerosol_radius 
                  endif
               else
                  Rdrop(k) = min_aerosol_radius                   ! avoid division by zero in aqueous subroutines
               endif
               ! Rdrop(k) = ( (rho(k) / rhol)*(qcloud(k) / num_conc(k))*(0.75/pi) )**(1./3.)   ! Calculates radius

               ! if (Rdrop(k) .lt. min_aerosol_radius) then
               !    Rdrop(k) = min_aerosol_radius                   ! avoid division by zero in aqueous subroutines
               ! end if
            end do

            if (do_debug_output) then
               do k = 1, nzm
                  if ( qcloud(k) > 0 ) then
                     write (*, *) 'k, Rdrop, num_conc, qcl', k, Rdrop(k), num_conc(k), qcloud(k)
                  endif
               end do
            end if

            ! Convert aqueous inputs kg/kg to M  (mol/L)
            ! if QC is zero but AQ field is nonzero, need to do something else
            do ispecies = 1, naqchem_fields
               aqchem_conc(:, ispecies) = aqchem_field(i, j, :, ispecies) * rhol / molwt(ispecies)       ! rhol in kg m-3; molwt in g/mol (1000s cancel out)
               ! aqchem_conc(:, ispecies) = aqchem_conc(:, ispecies)/dummy_water

               do k = 1, nzm
                  if (qcloud(k) .gt. cloud_threshold) then
                     aqchem_conc(k, ispecies) = aqchem_conc(k, ispecies)/qcloud(k)
                  else
                     aqchem_conc(k, ispecies) = 0. ! could convert this straggling aq to gas?
                  end if
               end do

               aqgas_conc(:, ispecies) = aqchem_gasprod_field(i, j, :, ispecies) * pressure_atm(:)    
            end do

            if (do_iepox_droplet_chem) then
               call iepox_aqueous_tendencies(nzm, tabs0, pressure_atm, Rdrop, &
                                             water_vol_frac, Haq, NO3aq, SO4aq, HSO4aq, &
                                             aqgas_conc, aqchem_conc, &  ! input conc fields
                                             aq_gasprod_tend, aq_tend, do_debug_output)   ! output tend fields

               if (maxval(abs(aq_gasprod_tend)) > 1.e3) then
                  print *, "aq_gasprod_tend exploded AFTER aqueous!"
                  print *, maxval(abs(aq_gasprod_tend))
                  stop
               endif

               call isop1nit_aqueous_tendencies(nzm, tabs0, pressure_atm, Rdrop, &
                                                water_vol_frac, aqgas_conc, aqchem_conc, &  ! input conc fields
                                                aq_gasprod_tend, aq_tend, do_debug_output)   ! output tend fields
            end if

            do ispecies = 1, naqchem_fields
               do k = 1, nzm
                  if (qcloud(k) .le. cloud_threshold) then
                     aq_tend(k, ispecies) = 0.
                     aq_gasprod_tend(k, ispecies) = 0.
                  end if
               end do
            end do

            !        if (do_debug_output) then
            !           write (*,*) 'preconv:  k, iepoxg, iepoxg_tend, iepoxa, iepoxa_tend'
            !           do k=1,3
            !              write(*,*) k, aqgas_conc(k, iIEPOX), aq_gasprod_tend(k, iIEPOX), aqchem_conc(k,iIEPOX), aq_tend(k, iIEPOX)
            !           end do
            !        end if

            ! Convert aqueous output tendencies back to model dims
            do ispecies = 1, naqchem_fields
               aq_tend(:, ispecies) = aq_tend(:, ispecies)*qcloud(:)*molwt(ispecies)/rhol      ! Reverses prev. calculations back to kg/kg
               aq_gasprod_tend(:, ispecies) = aq_gasprod_tend(:, ispecies)*(p0/pres(:))           ! Reverses prev. calculations back to mol/mol (ppv)
            end do

            !        if (do_debug_output) then
            !           write (*,*) 'postconv:  k,  iepoxg, iepoxg_tend, iepoxa, iepoxa_tend'
            !           do k=1,3
            !              write(*,*) k, aqgas_conc(k, iIEPOX), aq_gasprod_tend(k, iIEPOX), aqchem_conc(k,iIEPOX), aq_tend(k, iIEPOX)
            !           end do
            !        end if

            ! Prevent fields from going below 0 (adjusted tendency)
            do k = 1, nzm
               aq_adjusted_tendency = aq_tend(k, :)

               where (aqchem_field(i, j, k, :) + dtn * aq_tend(k, :) < 0.)
                  aq_adjusted_tendency = -aqchem_field(i, j, k, :)/dtn
               end where

               aqchem_field(i, j, k, :) = aqchem_field(i, j, k, :) + dtn * aq_adjusted_tendency
               aqchem_horiz_mean_tend(k, :) = aqchem_horiz_mean_tend(k, :) + aq_adjusted_tendency

               aq_adjusted_tendency = aq_gasprod_tend(k, :)

               where (aqchem_gasprod_field(i, j, k, :) + dtn * aq_gasprod_tend(k, :) < 0.)
                  aq_adjusted_tendency = -aqchem_gasprod_field(i, j, k, :)/dtn
               end where

               aqchem_gasprod_field(i, j, k, :) = aqchem_gasprod_field(i, j, k, :) + dtn * aq_adjusted_tendency
               aqchem_gasprod_horiz_mean_tend(k, :) = aqchem_gasprod_horiz_mean_tend(k, :) + aq_adjusted_tendency
               ! aq_adj_tendency = aq_gasprod_tend(k, iIEPOX)

               if ( ( gchem_field(i, j, k, ind_IEPOX) + dtn * aq_gasprod_tend(k, iIEPOX) ) .lt. 0. ) then
                  aq_adjusted_tendency(iIEPOX) = -gchem_field(i, j, k, ind_IEPOX) / dtn
               end if

               if ( ( gchem_field(i, j, k, ind_ISOP1Nit) + dtn * aq_gasprod_tend(k, iIP1NIT) ) .lt. 0. ) then
                  aq_adjusted_tendency(iIP1NIT) = -gchem_field(i, j, k, ind_ISOP1Nit) / dtn
               end if

               if ( ( gchem_field(i, j, k, ind_ISOPDiNit) + dtn * aq_gasprod_tend(k, iIPDINIT) ) .lt. 0. ) then
                  aq_adjusted_tendency(iIPDINIT) = -gchem_field(i, j, k, ind_ISOPDiNit) / dtn
               end if

               if (do_iepox_droplet_chem) then
                  gchem_field(i, j, k, ind_IEPOX) = gchem_field(i, j, k, ind_IEPOX) + dtn*aq_adjusted_tendency(iIEPOX)
                  gchem_field(i, j, k, ind_ISOP1Nit) = gchem_field(i, j, k, ind_ISOP1Nit) + dtn*aq_adjusted_tendency(iIP1NIT)
                  gchem_field(i, j, k, ind_ISOPDiNit) = gchem_field(i, j, k, ind_ISOPDiNit) + dtn*aq_adjusted_tendency(iIPDINIT)
               end if
            end do
         end do
      end do

      archem_horiz_mean_tend(:, :) = 0.

      ! if (do_iepox_aero_chem) then
      !     do j = 1,ny
      !        do i = 1,nx
      !        do_debug_output = (j.eq.1.and.i.eq.1)
      !        do_debug_output = .false.
      !        ar_tend(:,:) = 0.
      !        aero_radius = min_aerosol_radius
      !       do k = 1,nzm
      !            if (micro_field(i,j,k,inad).gt.1) then
      !                aero_radius(k) = 0.5 * ((1/rho_aerosol) * (micro_field(i,j,k,iqad)/micro_field(i,j,k,inad)) * (.75/pi)) **(1./3.)* &
      !                    EXP(3*LOG(sigma_accum)**2)
      !            end if
      !            if (do_debug_output) then
      !                write(*,*) 'k, qad, nad, sigma, radius=' , k, micro_field(i,j,k,iqad), micro_field(i,j,k,inad), sigma_accum, aero_radius(k)
      !            end if
      !       end do
      !
      !       where(aero_radius.lt.min_aerosol_radius)
      !           aero_radius = min_aerosol_radius
      !       end where
      !
      !       call iepox_aero_transfer_rate(nzm, tabs0, pressure_atm, aero_radius, Haero, &
      !               actHaero, SO4aero, HSO4aero, OrgMF, override_gamma, aero_transfer_rate, rho_org_aerosol, do_debug_output)
      !
      !       if (do_debug_output) then
      !          do k = 1,nzm
      !               write(*,*) 'k, radius, aero_transfer_rate= ', k, aero_radius(k), aero_transfer_rate(k)
      !           end do
      !      end if

      !      do k = 1,nzm
      !          ! Multiply aero_transfer_rate by aerosol number concentration
      !          aero_transfer_rate(k) = aero_transfer_rate(k) * micro_field(i,j,k,inad)
      !          ! limit IEPOXg loss
      !          if (aero_transfer_rate(k) * dtn.ge.1.) then
      !              aero_transfer_rate(k) = 1./dtn
      !          end if
      !      end do

      !     do k = 1,nzm
      !          if (aero_radius(k).ge.1.e-12) then  ! avoid division by zero for 0 size aerosol
      !             ! apply IEPOXg loss
      !             IEPOX_transfer_ppv = dtn * aero_transfer_rate(k)*gchem_field(i,j,k,ind_IEPOX)
      !              gchem_field(i,j,k, ind_IEPOX) = gchem_field(i,j,k, ind_IEPOX) - IEPOX_transfer_ppv       ! *(p0/pres(k)) ! convert to ppv
      !             ! distribute aerosol mass gain (converting to kg aerosol/kg air)
      !              archem_field(i,j,k, iTETROLr) = archem_field(i,j,k, iTETROLr) + &
      !              FracTETROL*IEPOX_transfer_ppv * molwt_ar(iTETROLr)/28.96  !  WHAT IS AIR MW CONSTANT CALLED - REPLACE HERE and next line
      !              archem_field(i,j,k, iIEPOX_SO4r) = archem_field(i,j,k,iIEPOX_SO4r) + &
      !              FracIEPOX_SO4 * IEPOX_transfer_ppv * molwt(iIEPOX_SO4r)/28.96
      !
      !              archem_horiz_mean_tend(k,iTETROLr) =  archem_horiz_mean_tend(k,iTETROLr) + FracTETROL * IEPOX_transfer_ppv * molwt_ar(iTETROLr)/28.96
      !              archem_horiz_mean_tend(k,iIEPOX_SO4r) =  archem_horiz_mean_tend(k,iIEPOX_SO4r) + FracIEPOX_SO4 * IEPOX_transfer_ppv * molwt_ar(iIEPOX_SO4r)/28.96

      !          end if
      !      end do
      !      end do
      !  end do
      !end if

   end subroutine het_chem_driver

   subroutine het_chem_finalize()
      implicit none

      if (isallocatedHetChem) then
         deallocate (Haq, NO3aq, SO4aq, HSO4aq)
         deallocate (Haero, SO4aero, HSO4aero)
      end if
   end subroutine het_chem_finalize

end module het_chem
