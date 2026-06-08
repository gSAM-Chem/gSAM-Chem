module radar
  use micro_params, only: nzeror, nzeros, nzerog, rhor, rhos, rhog
  implicit none
  private
  public :: radar_reflectivity_3d

  ! Kind
  integer, parameter :: dp = 4

  ! Physical constants
  real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: K2w  = 0.93_dp    ! |K|^2 for liquid water (S-band)
  real(dp), parameter :: K2i  = 0.197_dp   ! |K|^2 for (dry) ice/snow/graupel (S-band)
  real(dp), parameter :: mm6m3_from_m = 1.0e18_dp   ! convert m^6 to mm^6

  ! Melt band for dielectric blending (linear):
  real(dp), parameter :: T_melt0 = 273.15_dp - 2.0_dp  ! start blending at -2 C
  real(dp), parameter :: T_melt1 = 273.15_dp  ! fully wet at 0 C

contains

  !---------------------------------------------------------------------------
  ! Compute 3-D equivalent reflectivity (dBZ) and composite (column max).
  !
  ! Inputs:
  !   qr, qs, qg : mixing ratios [kg/kg], dimensions (nz, ny, nx)
  !   rho        : air density [kg/m^3], dimension (nz)
  !
  ! Outputs:
  !   Ze_dBZ     : 3-D water-equivalent reflectivity [dBZ], (nz, ny, nx)
  !   Zcomp_dBZ  : 2-D composite (column max) [dBZ], (ny, nx)
  !   Marat Khairoutdinov (2025)
  !---------------------------------------------------------------------------
  subroutine radar_reflectivity_3d(qr, qs, qg, tabs, rho, nx, ny, nz, Ze_dBZ, Zcomp_dBZ)
    real, intent(in)  :: qr(nx,ny,nz), qs(nx,ny,nz), qg(nx,ny,nz), tabs(nx,ny,nz)
    real, intent(in)  :: rho(nz)                ! size nz
    real, intent(out) :: Ze_dBZ(nx,ny,nz)
    real, intent(out) :: Zcomp_dBZ(nx,ny)
    integer, intent(in) :: nx, ny, nz

    integer :: k, j, i
    real(dp) :: lwc_r, lwc_s, lwc_g                 ! kg/m^3
    real(dp) :: lam_r, lam_s, lam_g                 ! m^-1
    real(dp) :: Zr, Zs, Zg, Ze_lin                  ! mm^6 m^-3 (linear)
    real(dp) :: K2_ice_eff
    real(dp), parameter :: tiny_lwc  = 1.0e-8_dp    ! kg/m^3 guard
    real(dp), parameter :: Z_floor   = 1.0e-6_dp    ! mm^6 m^-3 floor to avoid log(0)
    real(dp), parameter :: dBZ_min   = -30.0_dp     ! low-end floor for display
    real(dp), parameter :: dBZ_max   = 75.0_dp      ! optional cap


    Zcomp_dBZ(:,:) = dBZ_min

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          ! LWC/IWC per species
          lwc_r = rho(k) * max(qr(i,j,k), 0.0_dp)
          lwc_s = rho(k) * max(qs(i,j,k), 0.0_dp)
          lwc_g = rho(k) * max(qg(i,j,k), 0.0_dp)

          ! If negligible content, skip species
          if (lwc_r > tiny_lwc) then
            lam_r = ( (pi * rhor * Nzeror) / lwc_r )**(0.25_dp)       ! m^-1
            Zr    = 720.0_dp * Nzeror / (lam_r**7) * mm6m3_from_m       ! mm^6 m^-3
          else
            Zr = 0.0_dp
          end if

          if (lwc_s > tiny_lwc) then
            lam_s = ( (pi * rhos * Nzeros) / lwc_s )**(0.25_dp)
            Zs    = 720.0_dp * Nzeros / (lam_s**7) * mm6m3_from_m
          else
            Zs = 0.0_dp
          end if

          if (lwc_g > tiny_lwc) then
            lam_g = ( (pi * rhog * Nzerog) / lwc_g )**(0.25_dp)
            Zg    = 720.0_dp * Nzerog / (lam_g**7) * mm6m3_from_m
          else
            Zg = 0.0_dp
          end if
        ! --- Wet-ice dielectric handling - ice particle coated with liquid water.
        ! Linear blend of |K|^2 between ice and water across [T_melt0, T_melt1]
          if (tabs(i,j,k) >= T_melt1) then
            K2_ice_eff = K2w
          else if (tabs(i,j,k) <= T_melt0) then
            K2_ice_eff = K2i
          else
            K2_ice_eff = K2i + (K2w - K2i) * ( (tabs(i,j,k) - T_melt0) / (T_melt1 - T_melt0) )
          end if

          ! Water-equivalent Ze (linear), applying dielectric factors
          Ze_lin = ( Zr * K2w + (Zs + Zg) * K2_ice_eff ) / K2w
          Ze_lin = max(Ze_lin, Z_floor)

          Ze_dBZ(i,j,k) = 10.0_dp * log10(Ze_lin)
          Ze_dBZ(i,j,k) = min( max(Ze_dBZ(i,j,k), dBZ_min), dBZ_max )

          ! Update composite (column max)
          if (Ze_dBZ(i,j,k) > Zcomp_dBZ(i,j)) Zcomp_dBZ(i,j) = Ze_dBZ(i,j,k)

        end do
      end do
    end do

  end subroutine radar_reflectivity_3d

end module radar

