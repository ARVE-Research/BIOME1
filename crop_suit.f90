module crop_suitability

use common_vars, only : sp,dp

implicit none

!subroutines in this module
public :: ph_soil, c_density, crop_suit

contains

!---------------------------------------------------------------------
subroutine ph_soil(phaq,ph_avg)
implicit none

!arguments
real(sp), dimension(:), intent(in)  :: phaq
real(dp),               intent(out) :: ph_avg

!local variables
integer :: layers

if (phaq(2) > 0.) then                          ! this subroutine takes the weighted average of the ph in top 30
  ph_avg = (2. * phaq(1) + phaq(2)) / 3.        ! cm of soil... if there's only one layer then it just uses the 
else                                            ! ph of the top layer of soil
  ph_avg = phaq(1)
end if

end subroutine ph_soil

!--------------------------------------------------------------------
subroutine c_density(totc,bulk,c_dens)
implicit none

!arguments
real(sp), dimension(:), intent(in) :: totc    !(g C kg-1)
real(sp), dimension(:), intent(in) :: bulk    !(kg dm-3)
real(dp), 		intent(out):: c_dens  !(kg m-2)

!local variables
real(dp) :: tmp
real(dp) :: depth  !(m)


tmp =  totc(1) * bulk(1) + 0.5d0 * (totc(2) * bulk(2))    

depth = 0.3d0

c_dens = tmp * depth 			!check units, output should be kg m-2


end subroutine c_density
!---------------------------------------------------------------

subroutine crop_suit(gdd5,aalpha,ph_avg,c_dens,S_index,P_index)

implicit none

!arguments
real (sp), intent(in)	:: gdd5 
real (sp), intent(in)	:: aalpha !mean annual alpha
real (dp), intent(in) 	:: ph_avg
real (dp), intent(in)	:: c_dens
real (dp), intent(out) 	:: S_index
real (dp), intent(out) 	:: P_index

!local variables
real(dp) :: f_gdd5 		!function for growing degree days for crop suitability
real(dp) :: f_alpha 		!function for annual mean of alpha for crop suitability
real(dp) :: f_ph_soil		!function for soil ph for crop suitability
real(dp) :: f_c_density		!function for carbon density in soil for crop suitability

f_gdd5 = 1.d0 / (1.d0 + exp(0.0052d0 *(1334.d0 - gdd5)))

f_alpha = 1.d0 / (1.d0 + exp(14.705d0 *(0.3295d0 - aalpha)))

if  (ph_avg <= 6.5d0) then
	f_ph_soil = -2.085d0 + 0.475d0 * ph_avg
else if (ph_avg > 6.5d0 .AND. ph_avg < 8.d0) then
	f_ph_soil = 1.d0
else 
	f_ph_soil = 1.d0 - 2.d0 * (ph_avg - 8.d0)
end if

f_ph_soil = max(f_ph_soil,0.d0)

f_c_density = (3.9157d0 / (1.d0 + exp(1.3766d0 * (3.468d0 - c_dens)))) * &
              (3.9157d0 / (1.d0 + exp(-0.0791d0*(-27.33d0 - c_dens))))

S_index = f_gdd5 * f_alpha * f_ph_soil * f_c_density ! crop suitabilty index calculation

P_index = f_gdd5 * f_alpha  !pasture suitability index calculation

!if (s_index < 0.1d0) then
  !write(0,'(a10,2f10.3)')'gdd',gdd5,f_gdd5
  !write(0,'(a10,2f10.3)')'alpha',aalpha,f_alpha
  !write(0,'(a10,2f10.3)')'soil C',c_dens,f_c_density
  !write(0,'(a10,2f10.3)')'soil pH',ph_avg,f_ph_soil
  !write(0,'(a10,f10.3)')'crop suit',S_index
  !read(*,*)
!end if

end subroutine crop_suit

end module crop_suitability
