module radpet

public :: rad_and_pet  

contains

!-------------------------------

subroutine rad_and_pet(lat,elev,temp,prec,sunp,orbitvars,mdayl,mrad,ddayl,dpet,mpet,rad0,fdiff)

use common_vars,      only : sp,dp
use daily_insolation, only : dayins,phi,delta,dayl,ww
use airmass,          only : get_airmass,mc,ml,mbar,mo
use surfrad,          only : tcm,surface_rad,direct,diffuse
use netrad_pet,       only : pet,lw_rad,rl_pet

implicit none

!arguments

real(sp), intent(in) :: lat
real(sp), intent(in) :: elev

real(dp), dimension(365), intent(in) :: temp
real(dp), dimension(365), intent(in) :: prec
real(dp), dimension(365), intent(in) :: sunp

real(dp), dimension(3), intent(in) :: orbitvars

real(sp), intent(out) :: rad0

real(sp), dimension(12), intent(out) :: mdayl
real(sp), dimension(12), intent(out) :: mrad
real(sp), dimension(12), intent(out) :: mpet
real(sp), dimension(12), intent(out) :: fdiff  !diffuse fraction

real(dp), dimension(365), intent(out) :: ddayl
real(dp), dimension(365), intent(out) :: dpet

!variables

integer :: i
integer :: mon   !the month
integer :: dmo   !the day of the month
integer :: dyr   !the day of the year

integer, dimension(1) :: wm
integer, dimension(1) :: cm

real(sp) :: p_wm  !precipitation of the warmest month
real(sp) :: p_cm  !precipitation of the coldest month
real(sp) :: p     !relative atmospheric pressure (sea level = 1.0)

real(dp) :: xob
real(dp) :: perh
real(dp) :: ecc
real(dp) :: Pjj
real(dp) :: sw_rad
real(dp) :: oldpet

real(dp), dimension(365) :: dsun
real(dp), dimension(365) :: dfd  !daily diffuse fraction

real(dp) :: dayls
real(dp), parameter :: hsec = 3600.d0

!parameters

integer, dimension(12), parameter :: ndm = (/31,28,31,30,31,30,31,31,30,31,30,31/)

integer, dimension(12), parameter :: midday = (/16,44,75,105,136,166,197,228,258,289,319,350/)

!-------------------------------------------------------------------
!This routine needs as input:
!daily temp, prec, and % sun; elev, and lat
!output are: midmonth surface sw and day length
!            daily daylength and pet
!            total annual surface sw 
!-------------------------------------------------------------------

xob  = orbitvars(1)
perh = orbitvars(2)
ecc  = orbitvars(3)

wm   = maxloc(temp)
cm   = minloc(temp)
p_wm = prec(wm(1))
p_cm = prec(cm(1))

if (p_wm + p_cm > 0.) then
  Pjj = 2.d0 * (p_wm - p_cm) / (p_wm + p_cm)
  Pjj = max(Pjj,0.d0)
else
  Pjj = 0.d0
end if

!Pjj = max(Pjj,0.d0) !FLAG double check if this is allowed to be negative

tcm = minval(temp)

p = exp(-elev / 8000.)

!-----start a daily loop here

dyr=1
mpet = 0.

do mon=1,12           !month loop
  do dmo=1,ndm(mon)   !day of month loop
    
    !-------------------------------------------
    !calculate the daily top-of-atmosphere insolation

    phi = lat !phi (latitude in degrees) equals lat
    
    !write(*,*)'latitude',lat

    call dayins(ecc,perh,xob,phi,dyr,ww,dayl,delta) !the subroutine in daily_isolation.f90
    
    dayls = dayl * hsec * 0.001

    !---------------------------------------------
    !calculate elevation-corrected daily mean and instantaneous air mass

    call get_airmass(phi,delta,dayl,elev,mbar,mo,mc,ml)

    !---------------------------------------------
    !iterate to find surface shortwave, longwave, and PET

    pet = 0. !prec(dyr)
    oldpet = pet

    do i=1,10
      
      !write(0,*)pet

      !---------------------------------------------
      !calculate the daily surface shortwave insolation

      call surface_rad(ww,sunp(dyr)*0.01,dayl,mbar,mo,mc,ml,prec(dyr),Pjj,tcm,p,pet,direct,diffuse)

      sw_rad = direct + diffuse

      !---------------------------------------------
      !calculate the daily longwave, net radiation, and PET

      call rl_pet(temp(dyr),sunp(dyr)*0.01,dayl,sw_rad,lw_rad,pet) !subroutine in netrad_pet.f90
    
      pet = max(pet,0.d0)

      !---------------------------------------------
      !iterate to calculate surface radiation

      if (abs(pet - oldpet) < 0.1) exit

      oldpet = pet
      
      !if (direct > ww .or. diffuse < 0.d0 .or. pet < 0.d0) then
        !write(0,'(f7.3,3i5,8f9.3)')lat,mon,dmo,i,dayl,temp(dyr),prec(dyr),sunp(dyr),ww/dayls,direct/dayls,diffuse/dayls,pet
        !read(*,*)
      !end if

    end do

    ddayl(dyr) = dayl
    dsun(dyr)  = sw_rad
    dpet(dyr)  = pet

    if (sw_rad > 0.d0) then
      dfd(dyr) = diffuse / sw_rad
    else
      dfd(dyr) = 1.d0
    end if
          
    mpet(mon) = mpet(mon) + dpet(dyr)

    dyr = dyr + 1

  end do   !day of month
end do     !months in year

do mon=1,12
  mdayl(mon) = 3600.d0 * ddayl(midday(mon))  !convert to seconds

  if (mdayl(mon) > 0.d0) then
    mrad(mon) = 1000.d0 * dsun(midday(mon)) / mdayl(mon)
    fdiff(mon) = dfd(midday(mon))
  else
    mrad(mon) = 0.d0
  end if
end do

rad0=sum(dsun)

end subroutine rad_and_pet

end module radpet
