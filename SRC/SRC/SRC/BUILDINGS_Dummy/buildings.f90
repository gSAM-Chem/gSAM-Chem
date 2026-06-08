module buildings
use grid, only: nx, ny
real:: as_wall = 0.
real:: al_wall = 0.
integer:: shadow_mask(nx,ny)=0.
real(4), allocatable :: shadow_height(:,:)
real n_face_sunny, n_face_shady ! number of sunny and shady faces
real t_face_sunny, t_face_shady ! average temperature of sunny and shady walls, deg C
real lw_face, sw_face, shf_face
real lw_face_sfc, lw_face_face, lw_sfc_face
real skt_str, tair_str, tair_str_max, tair_str_min
real skt_out, tair_out, tair_out_max, tair_out_min
logical:: doshadows = .false.
logical:: doirgroundfromwalls=.false.  
logical:: doirwallsfromground=.false.  
logical:: doirwallsfromwalls=.false.
integer:: k_face_max = 0 ! maximum vertical index for roof tops
integer, allocatable :: bld_mask(:,:)
integer(2), allocatable :: facemask(:,:) 
integer(2), allocatable :: cell_wall(:,:,:)
real, allocatable :: normal_face1(:,:,:)
real, allocatable :: normal_face2(:,:,:)
integer(2), allocatable :: cell_sunny(:,:,:)
integer(2), allocatable :: sunnycells(:,:)
integer(2), allocatable :: k_base(:,:)
integer, allocatable :: bld_index(:,:)
real, allocatable :: flx_wall(:,:,:)
real, allocatable :: flx_face(:,:,:)
real, allocatable :: flx_bld(:,:,:)
real, allocatable :: flx_lw(:,:,:)
real, allocatable :: flx_sw(:,:,:)
real, allocatable :: flx_sfc_wall(:,:,:)
real, allocatable :: flx_wall_wall(:,:,:)
real(8), allocatable :: building_cooling_energy(:)
real(8), allocatable :: building_heating_energy(:)
real, allocatable :: t_face(:,:,:)
real, allocatable :: t_wall(:,:,:)
real, allocatable :: t_bld(:)
integer, allocatable :: ind_wall_sfc(:,:)
real, allocatable :: flx_wall_sfc(:,:) 
real, allocatable :: sangle_wall_sfc(:,:)
logical:: doAC = .false.
real:: t_AC
logical:: doHEAT = .false.
real:: t_HEAT





contains
subroutine buildings_stat()
end subroutine buildings_stat
subroutine buildings_stepout()
end subroutine buildings_stepout



end module buildings
