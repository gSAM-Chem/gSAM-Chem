module chemistry

  use grid, only: nx, ny, nzm, nz, &
       dimx1_s, dimx2_s, dimy1_s, dimy2_s, &
     z, zi, pres, adz, dz, dx, nx_gl, &
     time, dt, nstep, ncycle, nstat, nstatfrq, nrestart, day, &
     rank, dompi, masterproc, nsubdomains, save3Dbin, &
     case, caseid

  use rad, only : swDown, swUp
  use vars, only: rho, dtn, qcl, tabs0, pres, qv0, tabs, qpl, interactive_soil_wetness
  use cloudchem_Parameters, only: NVAR, NSPEC, NFIX, ind_IEPOX, ind_ISOP1Nit, ind_ISOPOOH, ind_ISOPDiNit, ind_OH, ind_O3, ind_CO, NREACT  ! NSPEC=NVAR+NFIX
  use params, only : dochem
  
  !! From KPP-generated files
  use cloudchem_Monitor, only: SPC_NAMES
  use cloudchem_Function, only: Fun
  use cloudchem_Global, only: RCONST, TEMP, C, ATOL, RTOL, SUN, C_M, C_H2O
  use cloudchem_Rates, only: Update_RCONST, Update_PHOTO
  use cloudchem_Integrator, only: Integrate

  use chemistry_params, only: p0, rhol, do_only_tropospheric_chemistry, do_transport_loss, &
                           do_OH_diurnal, OH_night, OH_day_peak, do_rainout, do_washout, do_convective_scavenging, wet_deposition_time_step, wet_deposition_species, k0_constants, cr_constants, dry_deposition_species, dry_deposition_velocities, &                           
                           do_megan_isoprene, do_surface_Isoprene_diurnal, do_bdsnp_no, &
                           do_IC_lightning, do_CTG_lightning, IC_decaria, CTG_decaria_reflectivity, CTG_price_and_rind, &
                           do_iepox_droplet_chem, do_iepox_aero_chem, hi_org, pHdrop, pHaero, & 
                           gas_init_name, gas_init_value, gas_out3D_name, flag_gchemvar_out3D, &
                           MW_air, avgd, rho_aerosol, sigma_accum, soil_wetness, tropopause_index, minimum_tropopause_height
  
  use chem_aqueous, only : aq_species_names, naqchem_fields, flag_aqchemvar_out3D, flag_aqchemgasvar_out3D, aq_gasprod_species_names
  use chem_aerosol, only : ar_species_names, narchem_fields, flag_archemvar_out3D
  use emissions, only : surface_emission_flux_driver, lightning_decaria_ic, lightning_decaria_ctg
  use deposition, only : dry_deposition_driver
  use wet_deposition, only : wet_deposition_driver
  use het_chem, only : het_chem_driver, het_chem_initialize, het_chem_finalize
  
  implicit none

  integer ngchem_fields, ngchem_fixed, ngchem_spec   ! equal to NVAR, NFIX, NSPEC  respectively
  integer nchem_fields_3Dsave 
  logical :: isallocatedCHEM = .false.

  real, allocatable, dimension(:,:,:,:) :: gchem_field            ! in ppv air
  real, allocatable, dimension(:,:,:,:) :: aqchem_field           ! in kg/kg
  real, allocatable, dimension(:,:,:,:) :: aqchem_gasprod_field   ! in kg/kg
  real, allocatable, dimension(:,:,:,:) :: archem_field           ! in kg/kg
  real, allocatable, dimension(:,:) :: gchem_profile_fixed        ! in ppv 
  real, allocatable, dimension(:) :: M_profile                    !  air density in molec/cm3
  real, allocatable, dimension(:,:) :: rate_const                 ! array of gas reaction rate constants

  real, allocatable, dimension(:,:) :: gchem_horiz_mean_tend      ! in ppv/s
  real, allocatable, dimension(:,:) :: aqchem_horiz_mean_tend     ! in kg/kg/s
  real, allocatable, dimension(:,:) :: aqchem_gasprod_horiz_mean_tend ! in kg/kg/s
  real, allocatable, dimension(:,:) :: archem_horiz_mean_tend     ! in kg/kg/s
  real, allocatable, dimension(:) :: g_depos_horiz_mean_tend_IEPOX ! in ppv/s
  real, allocatable, dimension(:) :: g_depos_horiz_mean_tend_ISOPOOH ! in ppv/s
  real, allocatable, dimension(:) :: soil_NO_emission_flux 
  real, allocatable, dimension(:) :: isop_emission_flux

   real, allocatable, dimension(:, :) :: tropopause_temp


  ! For wet deposition
  real, allocatable, dimension(:,:,:) :: previous_qcl  !  air density in molec/cm3
  real, allocatable, dimension(:,:,:) :: change_in_qcl  !  air density in molec/cm3

  real, allocatable, dimension(:,:,:) :: previous_qpl  ! array of gas reaction rate constants
  real, allocatable, dimension(:,:,:) :: change_in_qpl  ! array of gas reaction rate constants
  
  real, allocatable, dimension(:,:) :: & ! statistical arrays
       gchwle, &  ! resolved vertical flux
       gchadv, &  ! tendency due to vertical advection
       gchdiff, & ! tend. vertical diffusion
       gchwsb  ! SGS vertical flux

  real, allocatable, dimension(:,:) :: & ! statistical arrays
       aqchwle, &  ! resolved vertical flux
       aqchadv, &  ! tendency due to vertical advection
       aqchdiff, & ! tend. vertical diffusion
       aqchwsb  ! SGS vertical flux

  real, allocatable, dimension(:,:) :: & ! statistical arrays
       archwle, &  ! resolved vertical flux
       archadv, &  ! tendency due to vertical advection
       archdiff, & ! tend. vertical diffusion
       archwsb  ! SGS vertical flux

  real, allocatable, dimension(:,:,:) :: fluxbch, fluxtch  ! surface/top fluxes

  real gas_output_scale  ! convert all gas chem output to ppb
  real gas_input_scale  ! convert all gas chem input to parts per unit air
 
CONTAINS

  subroutine chem_setparm()
    implicit none
   
    integer ierr, ios, ios_missing_namelist, place_holder
    ngchem_fields = NVAR    ! number of advected che fields
    ngchem_fixed = NFIX     ! number of fixed chem profiles
    ngchem_spec = NSPEC  ! = NVAR + NFIX
    
    gas_output_scale = 1.e9  ! convert all gas chem output to ppbv
    gas_input_scale = 1./gas_output_scale  ! convert all input from ppb to ppunit
    
NAMELIST /CHEMISTRY/   do_only_tropospheric_chemistry, do_transport_loss, &
                           do_OH_diurnal, OH_night, OH_day_peak, do_rainout, do_washout, do_convective_scavenging, wet_deposition_species, wet_deposition_time_step, k0_constants, cr_constants, &
                           dry_deposition_species, dry_deposition_velocities, &
                           do_megan_isoprene, do_surface_Isoprene_diurnal, do_bdsnp_no, &
                           do_IC_lightning, do_CTG_lightning, IC_decaria, CTG_decaria_reflectivity, CTG_price_and_rind, &
                           do_iepox_droplet_chem, do_iepox_aero_chem, hi_org, pHdrop, pHaero, & 
                           gas_init_name, gas_init_value, gas_out3D_name, soil_wetness        
    ! read in namelist
    NAMELIST /BNCUIODSBJCB/ place_holder
    open(55,file='./CASES/'//trim(case)//'/prm', status='old',form='formatted') 
    read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
    rewind(55) !note that one must rewind before searching for new namelists
    read (55,CHEMISTRY,IOSTAT=ios)

    if (ios.ne.0) then
       if(ios.ne.ios_missing_namelist) then
           write(*,*) '****** ERROR: bad specification in CHEMISTRY namelist'
           rewind(55)
           read (55,CHEMISTRY) ! this should give a useful error message
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '****** No CHEMISTRY namelist in prm file *********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

    ! output namelist for documentation  
   if(masterproc) then
      open(unit=55,file='./OUT_STAT/'//trim(case)//'_'//trim(caseid)//'.nml',&
            form='formatted')      
      write (unit=55,nml=CHEMISTRY,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if
    
    ! allocate advection fields
   if ( dochem .and. ( .not.isallocatedCHEM) ) then
       ! allocate isoprene gas chemistry fields and related variables
       allocate(gchem_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,NVAR))
       allocate(gchem_profile_fixed(nzm, ngchem_fixed), M_profile(nzm))
       allocate(gchem_horiz_mean_tend(nzm, NVAR))
       allocate(rate_const(nzm, NREACT))
       allocate(g_depos_horiz_mean_tend_IEPOX(nzm))
       allocate(g_depos_horiz_mean_tend_ISOPOOH(nzm))
       allocate(soil_NO_emission_flux(nzm))
       allocate(isop_emission_flux(nzm))
       
       ! allocate aqueous IEPOX fields
       allocate(aqchem_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,naqchem_fields))
       allocate(aqchem_horiz_mean_tend(nzm, naqchem_fields))
       
       ! allocate gas product IEPOX fields
       allocate(aqchem_gasprod_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,naqchem_fields))
       allocate(aqchem_gasprod_horiz_mean_tend(nzm, naqchem_fields))

       ! allocate aerosol product fields
       allocate(archem_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,narchem_fields))
       allocate(archem_horiz_mean_tend(nzm, narchem_fields))
       
       allocate(gchwle(nz, ngchem_fields),gchadv(nz,ngchem_fields), &
            gchdiff(nz,ngchem_fields),gchwsb(nz,ngchem_fields))
       allocate(aqchwle(nz, naqchem_fields),aqchadv(nz,naqchem_fields), &
            aqchdiff(nz,naqchem_fields),aqchwsb(nz,naqchem_fields))
       ! Should have one here for gas products of aq chem FIX
       allocate(archwle(nz, narchem_fields),archadv(nz,narchem_fields), &
            archdiff(nz,narchem_fields),archwsb(nz,narchem_fields))

       allocate(fluxbch(nx, ny, ngchem_fields), fluxtch(nx, ny, ngchem_fields))

       ! Wet deposition
      allocate(previous_qcl(nx, ny, nzm), previous_qpl(nx, ny, nzm))
      allocate(change_in_qcl(nx, ny, nzm), change_in_qpl(nx, ny, nzm))

      allocate(tropopause_temp(nx, ny))
      
      previous_qcl(:,:,:) = -9999
      previous_qpl(:,:,:) = -9999

      change_in_qcl(:,:,:) = 0.
      change_in_qpl(:,:,:) = 0.

      interactive_soil_wetness(:, :) = soil_wetness

      isallocatedCHEM = .true.

    end if   

    if ( dochem .and. ( do_iepox_aero_chem .or. do_iepox_droplet_chem ) ) then 
      call het_chem_initialize(hi_org, pHdrop, pHaero)
   end if 

end subroutine chem_setparm
  
subroutine chem_init()
  ! called at start of run or restart
  implicit none
  integer i,j,k
  integer v_selected ! index of namelist input variables
  integer v ! index of kpp variable
  logical match
  
  rate_const = 0.
  M_profile = 0.001 * RHO * avgd / MW_air
  gchem_profile_fixed = 0. 
  
  do v_selected = 1,ngchem_spec
     match = .false.
     do v=ngchem_fields+1, ngchem_spec  ! search only over fixed-variable names
        if(gas_init_name(v_selected)==trim(SPC_NAMES(v))) then
           match=.true.
           exit
        end if
     end do   

     if(match) then
        gchem_profile_fixed(:, v-ngchem_fields) = gas_init_value(v_selected)* &
             gas_input_scale
        if (masterproc) then
           write(*,*) 'SET FIXED CHEM PROFILE: ', gas_init_name(v_selected),                 gas_init_value(v_selected)
        end if   
     end if
  end do   

  if(nrestart.eq.0) then
     ! initialize gas chem fields
     gchem_field = 0.
     ! compute conversion profile 
     do v_selected = 1,ngchem_spec
        match = .false.   
        do v = 1,ngchem_fields   ! search only over variable names  
           if(gas_init_name(v_selected)==trim(SPC_NAMES(v))) then
              match=.true.
              exit
           end if
        end do  
      
        if(match) then
           do i = 1, nx
              do j = 1, ny
                 gchem_field(i,j, :, v) = gas_init_value(v_selected) * gas_input_scale
              end do
           end do   
        end if
      
     end do
     ! initialize aqchem and archem fields as zero
     aqchem_field = 0.
     aqchem_gasprod_field = 0.
     archem_field = 0.
  end if  ! restart       
  
  ! set flags for 3d output based on namelist input 
  flag_gchemvar_out3D(:)=.false.
  nchem_fields_3Dsave = 0.

  do v = 1,ngchem_fields
     if(any(gas_out3D_name==trim(SPC_NAMES(v)))) then
        flag_gchemvar_out3D(v) = .true.
        nchem_fields_3Dsave = nchem_fields_3Dsave + 1
        if (masterproc) write(*,*) &
             'Chem 3d output field added: ', trim(SPC_NAMES(v))
     end if
  end do
  
   if ( do_iepox_droplet_chem ) then 
      flag_aqchemvar_out3D(:) = .true.
      flag_aqchemgasvar_out3D(:) = .true.
      nchem_fields_3Dsave = nchem_fields_3Dsave + naqchem_fields*2
   end if 

   if ( do_iepox_aero_chem ) then 
      flag_archemvar_out3D(:) = .true.
      nchem_fields_3Dsave = nchem_fields_3Dsave + narchem_fields
   end if
  
  ! initialize some statistics profiles, not output yet
  gchwle = 0.
  gchadv = 0.
  gchdiff = 0.
  gchwsb = 0.

  aqchwle = 0.
  aqchadv = 0.
  aqchdiff = 0.
  aqchwsb = 0.

  archwle = 0.
  archadv = 0.
  archdiff = 0.
  archwsb = 0.

  soil_NO_emission_flux(:) = 0.
  isop_emission_flux(:) = 0.

  ! top and bottom fluxes of fields
  fluxbch = 0.
  fluxtch = 0.

end subroutine chem_init 

subroutine chem_proc()
  implicit none

  integer :: i,j,k,n, ispecies
  real, dimension(nzm, NVAR) :: var_profile

  real, dimension(nzm, NVAR) :: gas_column_tend_profile ! in molecules/cm3/s
  real, dimension(NVAR) :: adjusted_tendency
  real, dimension(nzm,NFIX) :: fixed_profile

  real :: OH_conc ! OH gas concentration in molecules/cm3

  ! For KPP !
   INTEGER                :: IERR                      ! KPP success or failure flag
   INTEGER                :: ICNTRL (20)
   INTEGER                :: ISTATUS(20)
   REAL(8)                :: RCNTRL (20)
   REAL(8)                :: RSTATE (20)
  
  
  
  IERR      = 0                        ! KPP success or failure flag
  ISTATUS   = 0.0                   ! Rosenbrock output
  ICNTRL    = 0                        ! Rosenbrock input (integer)
  RCNTRL    = 0.0                   ! Rosenbrock input (real)
  RSTATE    = 0.0                   ! Rosenbrock output
  
  ICNTRL(1) = 1
  ICNTRL(2) = 1
  ICNTRL(4) = 40
  ICNTRL(15) = -1
  ICNTRL(16) = 1

  ATOL = 1.0e-2
  RTOL = 0.5e-2

  call t_startf ('chem_proc')

  gchem_horiz_mean_tend(:,:) = 0.
  g_depos_horiz_mean_tend_IEPOX(:) = 0.
  g_depos_horiz_mean_tend_ISOPOOH(:) = 0.
  M_profile = 0.001 * RHO * avgd / MW_air

  if ( do_only_tropospheric_chemistry .and. ( z(nzm) .gt. minimum_tropopause_height ) ) then 
      tropopause_index(:, :) = nzm
      call find_tropopause(tropopause_index, tropopause_temp)
   endif

  do j = 1, ny
     do i = 1, nx
         gas_column_tend_profile(:,:) = 0.

         do k = 1, nzm 
            C = 0.0
            TEMP = 0.0
            SUN = 0.0
            C_M = 0.0
            C_H2O = 0.0

            C_M = M_profile(k)                 ! Air density [molec/cm3]
            C_H2O = C_M * qv0(k) * 28.97 / 18.02 ! qv0 is in kg kg-1 INTERNALLY, g kg-1 in output!
            !!!! gchem_field(i, j, k, ind_H2O) = qv0(k) / 1000. * 28.97 / 18.02         ! Update water vapor mixing ratio (in mol mol-1 w.r.t. dry air)
            
            var_profile(k,:) = gchem_field(i,j,k,:) * M_profile(k)               ! Convert from mole fraction to [molec cm-3]
            fixed_profile(k,:) = gchem_profile_fixed(k,:) * M_profile(k)         ! Convert from mole fraction to [molec cm-3]

            TEMP = tabs(i, j, k)

            ! Set concentrations in KPP to the current box's concentration [molec cm-3]
            C(1:NVAR) = var_profile(k, :)
            !!!! C(NVAR+1:NSPEC) = fixed_profile(k, :)

            call Update_RCONST()

            SUN = MIN( ( swDown(i,k) + swUp(i, k) ) / 1300., 1.)
            call Update_PHOTO()

            rate_const(k, :) = RCONST                                            ! rate_const is the local version of RCONST

            if ( do_only_tropospheric_chemistry ) then
               if ( ( z(tropopause_index(i, j)) .gt. minimum_tropopause_height ) .and. ( k .gt. tropopause_index(i, j) + 5 ) ) then
                  rate_const(k, :) = 0.
               endif
            endif

            call Fun(var_profile(k, :), fixed_profile(k, :), rate_const(k, :), gas_column_tend_profile(k, :))     ! Calculate derivatives
            CALL Integrate(0.0, dtn, ICNTRL, RCNTRL, ISTATUS, RSTATE, IERR)       ! Use RCONST and Func to integrate the concentrations a timestep

            IF ( IERR < 0 ) THEN
               write(*,*) 'KPP failed here: ', i, j, k, IERR
               STOP
            ENDIF

            adjusted_tendency = ( ( C(1:NVAR)  / M_profile(k) ) - gchem_field(i,j,k,:) ) / dtn    ! Calculate the change in the species from one timestep to another
            gchem_field(i,j,k,:) = C(1:NVAR) / M_profile(k)
            gchem_horiz_mean_tend(k,:) = gchem_horiz_mean_tend(k,:) + adjusted_tendency ! / M_profile(k)
        enddo
     end do
  end do
  
  call check_nan

  gchem_horiz_mean_tend = gchem_horiz_mean_tend / (nx * ny)  

  ! Do dry deposition
  call dry_deposition_driver(gchem_field, g_depos_horiz_mean_tend_ISOPOOH, g_depos_horiz_mean_tend_IEPOX)

  ! Do wet deposition
  if ( do_convective_scavenging .or. do_rainout  .or. do_washout ) then
   if ( mod(nstep, wet_deposition_time_step) .eq. 0 ) then
      if ( ( all(previous_qcl .eq. -9999) ) .or. ( all(previous_qpl .eq. -9999) ) ) then
         previous_qcl = qcl
         previous_qpl = qpl
      else
         change_in_qcl = qcl - previous_qcl
         change_in_qpl = qpl - previous_qpl

         call wet_deposition_driver(gchem_field, M_profile, change_in_qcl, change_in_qpl)

         previous_qcl = qcl
         previous_qpl = qpl
      endif
    endif
  endif

   if ( do_iepox_aero_chem .or. do_iepox_droplet_chem ) then 
      call het_chem_driver(gchem_field, aqchem_field, aqchem_gasprod_field, aqchem_horiz_mean_tend, aqchem_gasprod_horiz_mean_tend, archem_field, archem_horiz_mean_tend)
   end if

  ! normalize tendencies
  aqchem_horiz_mean_tend = aqchem_horiz_mean_tend/(nx*ny)
  aqchem_gasprod_horiz_mean_tend = aqchem_gasprod_horiz_mean_tend/(nx*ny)
  archem_horiz_mean_tend = archem_horiz_mean_tend/(nx*ny)

  g_depos_horiz_mean_tend_ISOPOOH = g_depos_horiz_mean_tend_ISOPOOH/(nx*ny)
  g_depos_horiz_mean_tend_IEPOX = g_depos_horiz_mean_tend_IEPOX/(nx*ny)
  call t_stopf ('chem_proc')
end subroutine chem_proc

subroutine find_tropopause(tropopause_index, tropopause_temp)
   use grid, only : nx, ny, z
   use vars, only : tabs

   implicit none
   integer :: i,j,k
   
   integer, intent(out) :: tropopause_index(nx, ny)
   real, intent(out) :: tropopause_temp(nx, ny)

   do i = 1,nx
      do j = 1,ny
         tropopause_temp(i, j) = minval(tabs(i,j,:), dim = 1)
         tropopause_index(i, j) = minloc(tabs(i,j,:), dim = 1)
      enddo
   enddo
end subroutine find_tropopause

subroutine check_nan
   implicit none

   integer :: i, j, k, n

   do i = 1, nx
      do j = 1, ny
         do k = 1, nzm
            do n = 1, ngchem_fields
               if ( ( isnan(gchem_field(i,j,k,n)) ) .or. ( gchem_field(i,j,k,n) .gt. 1 ) ) then
                  print*, "******* There is a NaN or a high value in gchem_field. Check if your system is unstable! ****** ", gchem_field(i,j,k,n)
                  print*, "******* Value of M profile: ", M_profile(k)
                  print*, "i = ", i 
                  print*, "j = ", j 
                  print*, "k = ", k 
                  print*, "n = ", n 

                  write(*,*) 'TEMP=', tabs(i,j,k)
                  STOP
               endif
            enddo
         enddo
      enddo
   enddo
end subroutine check_nan

subroutine chem_hbuf_init(namelist, deflist, unitlist, status, average_type, count, chemcount)
   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*), average_type(*), count, chemcount, n

   character*16 name
   character*16 tend_name
   character*80 longname
   character*10 units

   chemcount = 0

   do n = 1, ngchem_fields
      name = trim( SPC_NAMES(n) )
      longname = trim( SPC_NAMES(n) )
      units = 'ppbv'
      call add_to_namelist(count, chemcount, name, longname, units, 0)

      tend_name = trim(SPC_NAMES(n)) // '+'
      longname = trim(SPC_NAMES(n)) // ' tendency due to gas reaction'
      units = 'ppbv/s'
      call add_to_namelist(count, chemcount, tend_name, longname, units, 0)
   end do

   if ( ngchem_fixed .gt. 1)  then
      do n = 1, ngchem_fixed
         name = trim(SPC_NAMES(n + ngchem_fields))
         longname = name // ' Fixed Species'
         units = 'ppbv'
         call add_to_namelist(count, chemcount, name, longname, units, 0)
      end do
   endif

   if ( masterproc ) then
      write(*,*) 'Added ', chemcount, ' arrays to statistics for gaseous chemical species!'
   end if

   if ( do_iepox_droplet_chem ) then
      do n = 1,naqchem_fields
         name = trim(aq_species_names(n))
         longname = trim(aq_species_names(n))
         units = 'kg/kg'
         call add_to_namelist(count,chemcount,name,longname,units,0)

         longname = trim(aq_species_names(n)) // ' tendency due to aqueous reaction'
         units = 'kg/kg/s'
         call add_to_namelist(count,chemcount,trim(aq_species_names(n))//'+',longname,units,0)
      end do

      do n = 1,naqchem_fields
         name = trim(aq_gasprod_species_names(n))
         longname = trim(aq_gasprod_species_names(n))
         units = 'ppbv'
         call add_to_namelist(count,chemcount,name,longname,units,0)

         longname = trim(aq_gasprod_species_names(n)) // 'gas tendency due to aqueous reaction'
         units = 'ppbv/s'
         call add_to_namelist(count,chemcount,trim(aq_gasprod_species_names(n))//'+',longname,units,0)
      end do
   end if
   
   if ( do_iepox_aero_chem ) then 
      do n = 1,narchem_fields
         name = trim(ar_species_names(n))
         longname = trim(ar_species_names(n))
         units = 'kg/kg'
         call add_to_namelist(count,chemcount,name,longname,units,0)

         longname = trim(ar_species_names(n)) // ' tendency due to aerosol chemistry'
         units = 'kg/kg/s'
         call add_to_namelist(count,chemcount,trim(ar_species_names(n))//'+',longname,units,0)
      end do
   end if

   name = 'IPOOHd+'
   longname = 'IPOOH gas tendency due to dry deposition'
   units = 'ppbv/s'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   name = 'IEPOXd+'
   longname = 'IEPOX gas tendence due to dry deposition'
   units = 'ppbv/s'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   name = 'SOILNO_E'
   longname = 'Soil NOx Emission Flux'
   units = 'ppbv m/s'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   name = 'ISOP_E'
   longname = 'Isoprene Emission Flux'
   units = 'ppbv m/s'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   name = 'SOIL_WET'
   longname = 'Soil wetness (interactive)'
   units = 'unitless'
   call add_to_namelist(count, chemcount, name, longname, units, 0)

   if ( masterproc ) then
      write(*,*) 'Added ', chemcount, ' arrays to statistics for gas and aqueous chemical species!'
   end if

end subroutine chem_hbuf_init

subroutine chem_finalize()
  ! deallocate
  implicit none
  integer :: ierr

  if(isallocatedCHEM) then
     deallocate(gchem_field, STAT=ierr)
     deallocate(aqchem_field, STAT=ierr)
     deallocate(aqchem_gasprod_field, STAT=ierr)
     deallocate(archem_field, STAT=ierr)
     deallocate(gchem_profile_fixed, gchwle, gchadv, gchdiff, gchwsb,M_profile,STAT=ierr)
     deallocate(aqchwle, aqchadv, aqchdiff, aqchwsb,STAT=ierr)
     deallocate(archwle, archadv, archdiff, archwsb,STAT=ierr)
     deallocate(gchem_horiz_mean_tend)
     deallocate(aqchem_horiz_mean_tend)
     deallocate(aqchem_gasprod_horiz_mean_tend)
     deallocate(archem_horiz_mean_tend)
     deallocate(g_depos_horiz_mean_tend_IEPOX)
     deallocate(g_depos_horiz_mean_tend_ISOPOOH)
     deallocate(rate_const)
     deallocate(fluxbch, fluxtch)
     deallocate(tropopause_temp)

     if ( do_iepox_aero_chem .or. do_iepox_droplet_chem ) then 
      call het_chem_finalize()
     endif 

     if(ierr.ne.0) then
        write(*,*) 'Failed to deallocated chem arrays on proc ', rank
     end if
  end if
  
end subroutine chem_finalize

! chem_flux is a wrapper function to call emissions module
subroutine chem_flux()
   soil_NO_emission_flux(:) = 0.
   isop_emission_flux(:) = 0.

   call surface_emission_flux_driver(fluxbch, M_profile, isop_emission_flux, soil_NO_emission_flux)

   if ( do_CTG_lightning ) then
      call lightning_decaria_ctg(gchem_field)
   endif

   if ( do_IC_lightning ) then 
      call lightning_decaria_ic(gchem_field)
   endif

end subroutine chem_flux

subroutine chem_print()
  implicit none
end subroutine chem_print  

subroutine chem_statistics()
  use hbuffer, only: hbuf_put
  implicit none
  ! average fields in space for .stat file

  real, dimension(nzm) :: tr0, tendency
  real, dimension(nzm) :: zeros
  real, allocatable, dimension(:) :: soil_wetness_array

  real factor_xy
  integer i,j,k,m, n, ii, jj, nn, ncond

  character*16 name
  character*16 tend_name
  factor_xy = 1./float(nx*ny)

  call t_startf ('chem_statistics')

  zeros(:) = 0.
  do n = 1,ngchem_fields
    ! compute horizontal mean of all gas chem fields  
    do k = 1,nzm
       tr0(k) = SUM(gchem_field(1:nx,1:ny,k,n))
    end do

    !if(n.eq.1) write(*,*) 'IP1O2 mean, (1,1,1) = ', tr0(1), gchem_field(1,1,1,1)*1.e9
    
    call hbuf_put(trim(SPC_NAMES(n)), tr0, gas_output_scale*factor_xy) ! factor is 1/(nx * ny)
    call hbuf_put(trim(SPC_NAMES(n))//'+', gchem_horiz_mean_tend(:, n), gas_output_scale)
  end do  

  if ( ngchem_fixed .gt. 1)  then
        do n = 1,ngchem_fixed
                name = trim(SPC_NAMES(n+ngchem_fields))
                call hbuf_put(name, gchem_profile_fixed(:, n), gas_output_scale)
        end do
  endif

   if ( do_iepox_droplet_chem ) then
      do n = 1,naqchem_fields
         ! compute horizontal mean of all aq chem fields  
         do k = 1,nzm
            tr0(k) = SUM(aqchem_field(1:nx,1:ny,k,n))
         end do

         call hbuf_put(trim(aq_species_names(n)), tr0, factor_xy) ! factor is 1/(nx * ny)
         call hbuf_put(trim(aq_species_names(n))//'+', aqchem_horiz_mean_tend(:, n), 1.)
      end do

      do n = 1,naqchem_fields
         do k = 1,nzm
            tr0(k) = SUM(aqchem_gasprod_field(1:nx,1:ny,k,n))
         end do
         
         call hbuf_put(trim(aq_gasprod_species_names(n)), tr0, factor_xy*gas_output_scale) ! factor is 1/(nx * ny)
         call hbuf_put(trim(aq_gasprod_species_names(n))//'+', aqchem_gasprod_horiz_mean_tend(:, n), gas_output_scale)   
      end do  
   end if

   if ( do_iepox_aero_chem ) then
      do n = 1,narchem_fields
         do k = 1,nzm
            tr0(k) = SUM(archem_field(1:nx,1:ny,k,n))
         end do
         
         call hbuf_put(trim(ar_species_names(n)), tr0, factor_xy) ! factor is 1/(nx * ny)
         call hbuf_put(trim(ar_species_names(n))//'+', archem_horiz_mean_tend(:, n), 1.)
      end do
   end if

  call hbuf_put('IPOOHd+', g_depos_horiz_mean_tend_ISOPOOH, gas_output_scale)
  call hbuf_put('IEPOXd+', g_depos_horiz_mean_tend_IEPOX, gas_output_scale)

  call hbuf_put('SOILNO_E', soil_NO_emission_flux, gas_output_scale * factor_xy)
  call hbuf_put('ISOP_E', isop_emission_flux, gas_output_scale * factor_xy)
  
  allocate(soil_wetness_array(nzm))
  soil_wetness_array = 0.                         ! To avoid uninitialized values in the output file
  soil_wetness_array(1) = SUM(interactive_soil_wetness(:, :))
  call hbuf_put('SOIL_WET', soil_wetness_array, factor_xy)
  deallocate(soil_wetness_array)

  call t_stopf ('chem_statistics')

end subroutine chem_statistics  

subroutine chem_write_fields3D(nfields1)
  use grid, only: nsubdomains, save3Dbin, nfiles3D
  implicit none
  
  integer, intent(inout) :: nfields1
  character *80 long_name
  character *16 name
  character *10 units
  integer :: i, j, k, f 
  real(4), dimension(nx,ny,nzm) :: tmp

  do f = 1,ngchem_fields
     if(flag_gchemvar_out3D(f)) then  
        nfields1=nfields1+1
        do k=1,nzm
           do j=1,ny
              do i=1,nx
                 tmp(i,j,k)=gchem_field(i,j,k,f)*gas_output_scale
              end do
           end do
        end do
        name=TRIM(SPC_NAMES(f))
        long_name=TRIM(SPC_NAMES(f))
        units='ppbv'
        call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains, nfiles3D, 1)
     end if
  end do

   if ( do_iepox_droplet_chem ) then
      do f = 1,naqchem_fields
         if(flag_aqchemvar_out3D(f)) then  
            nfields1=nfields1+1
            tmp = aqchem_field(1:nx,1:ny,:,f)
            name=TRIM(aq_species_names(f))
            long_name=TRIM(aq_species_names(f))
            units='kg/kg'
            call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
               save3Dbin,dompi,rank,nsubdomains, nfiles3D, 1)
         end if
      end do
   
      do f = 1,naqchem_fields
         if(flag_aqchemgasvar_out3D(f)) then  
            nfields1=nfields1+1
            tmp = aqchem_gasprod_field(1:nx, 1:ny,:,f) ! account for ghost cells - specify 1:nx, 1,ny
            name=TRIM(aq_gasprod_species_names(f))
            long_name=TRIM(aq_gasprod_species_names(f))
            units='kg/kg'
            call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
               save3Dbin,dompi,rank,nsubdomains, nfiles3D, 1)
         end if
      end do
   end if

   if ( do_iepox_aero_chem ) then
      do f = 1,narchem_fields
         if(flag_archemvar_out3D(f)) then  
            nfields1=nfields1+1
            tmp = archem_field(1:nx,1:ny,:,f)
            name=TRIM(ar_species_names(f))
            long_name=TRIM(ar_species_names(f))
            units='kg/kg'
            call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
               save3Dbin,dompi,rank,nsubdomains, nfiles3D, 1)
         end if
      end do
   end if
  
end subroutine chem_write_fields3D
end module chemistry
