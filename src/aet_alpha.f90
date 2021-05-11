module aetmod
  
use common_vars, only : dp

implicit none

!global variables

real(dp) :: m       		!instantaneous soil moisture (mm) 

!subroutines in this module

public :: aet_alpha

contains

!-----------------------------------------

subroutine aet_alpha(awc,dayl,prec,pet,aet,alpha)

use common_vars, only : dp

implicit none


!arguments

real(dp), intent(in)  :: awc   !total water holding capacity in the top meter of soil (mm); actually, total whc (mm3/mm3) multiplied by soil depth (mm) PMC 25 jan '10
real(dp), intent(in)  :: dayl  !day length (hr)
real(dp), intent(in)  :: prec  !daily precipitation (mm)
real(dp), intent(in)  :: pet   !daily potential evapotranspiration (mm)
real(dp), intent(out) :: aet   !daily actual evapotranspiration (mm)
real(dp), intent(out) :: alpha !ratio of AET/PET (fraction)

!parameters

real(dp), parameter :: etmax = 1.d0  !maximum equilibrium evapotranspiration rate (mm h-1)

!local variables

real(dp) :: Emd     !Emax: maximum evapotranspiration rate from saturated soils(mm d-1)
real(dp) :: supply  !water supply (mm) 
real(dp) :: demand  !water demand (mm)

!-------------------

!demand = PET (simple)

demand = pet

!first calculate supply based on previous day's soil moisture

Emd = etmax * dayl !first find total etmax for the whole day

supply = min(Emd,m) !Emd * m / awc

!then calculate aet for this day using supply

aet = min(supply,demand)

if (pet > 0.d0) then
  alpha = aet / pet
else
  alpha = 1.d0
end if

!then update soil moisture to reflect current day's aet

m=min(m+(prec-aet),awc) 

!write(*,'(7f10.3)')prec,emd,m,supply,pet,aet,alpha		

end subroutine aet_alpha

!------------------------

end module aetmod
