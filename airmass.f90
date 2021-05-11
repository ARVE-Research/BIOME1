module airmass

use common_vars, only : sp,dp,pi,pir

!module shared variables

real(dp) :: mbar
real(dp) :: mo
real(dp) :: mc
real(dp) :: ml

!module parameters

real(dp), parameter :: w   = 15.0d0  !solar angular velocity (degrees hr-1)
real(dp), parameter :: rw  = pir * w !solar angular velocity (radians hr-1)

!module subroutines and functions

public  :: get_airmass
private :: m
private :: F
private :: elev_corr

contains

!-------------------------------------------------------------------

subroutine get_airmass(lat,delta,dayl,elevation,mbar,mo,mc,ml)

!This code is based on the paper:
!X. Yin (1997) Optical air mass: Daily integration and its applications, Meteorol. Atmos. Phys. 63, 227-233
!Jed Kaplan, EPFL, 2008

implicit none

!arguments

real(dp), intent(in) :: lat          !latitude (degrees)
real(dp), intent(in) :: delta        !solar declination (degrees)
real(dp), intent(in) :: dayl         !day length (hours)
real(sp), intent(in) :: elevation    !elevation above sea level (m)

real(dp), intent(out) :: mbar        !daytime mean optical air mass (unitless, 1 at equatorial noon)
real(dp), intent(out) :: mo          !air mass at cosine zenith angle maximum
real(dp), intent(out) :: mc          !air mass at cosine zenith angle medium
real(dp), intent(out) :: ml          !air mass at cosine zenith angle bottom quarter range point

!parameters

real(dp), parameter :: m0  =  1.0d0  !air mass at 0 degree solar zenith angle
real(dp), parameter :: m80 =  5.6d0  !air mass at 80 degree solar zenith angle
real(dp), parameter :: m90 = 39.7d0  !air mass at 90 degree solar zenith angle

!real(dp), parameter :: cos80 = 0.173648177666930442d0 !cos(80) (degrees)

!local variables

real(dp) :: rlat      !latitude (radians)
real(dp) :: rdelta    !solar declination (radians)

real(dp) :: t1        !number of hours between sunrise/sunset and solar noon (hr)
real(dp) :: t80       !solar hour corresponding to the 80 degree zenith angle

real(dp) :: cos80
real(dp), dimension(3), target :: c00 != (/ 0.008307d0, 1.021d0, -0.01259d0 /) !air mass coefficients for solar zenith angle <=80 degrees
real(dp), dimension(3), target :: c80 != (/ 0.03716d0,  1.538d0, -1.689d0   /) !air mass coefficients for solar zenith angle  >80 degrees

real(dp) :: t    !solar hour (hr)

!real(dp) :: m    !instantaneous air mass (unitless)

real(dp) :: Z    !solar zenith angle (degrees)
real(dp) :: Zn   !lesser of solar zenith angle at sunset or at midnight (degrees)
real(dp) :: Z0   !zenith angle at solar noon (degrees)
real(dp) :: cosZ !cosine solar zenith angle (fraction), used in calculation of instantaneous air mass

real(dp) :: l

integer :: steps  !integer number of time steps
integer :: i      !counter

real(dp), allocatable, dimension(:) :: mvect

real(dp) :: sinlat
real(dp) :: coslat
real(dp) :: sindel
real(dp) :: cosdel

real(dp)                        :: a   !values in equation 2.6b
real(dp)                        :: b
real(dp), pointer, dimension(:) :: c

real(dp) :: tmp1
real(dp) :: tmp2
real(dp) :: tmp3

real(dp) :: tinv

real(dp) :: rZ0
real(dp) :: rZn

!-------------------------------------

cos80 = cos(80.d0 * pir)

c00(1) = 0.008307d0
c00(2) = (m0 - m80) * (c00(1) + 1.0d0) * (c00(1) + cos80) / (cos80 - 1.0d0)
c00(3) = m0 - c00(2) / (c00(1) + 1.0d0)

c80(1) = 0.037160d0
c80(2) = (m90 - m80) * c80(1) * (c80(1) + cos80) / cos80
c80(3) = m90 - c80(2) / c80(1)

!write(0,*)cos80
!write(0,*)c00
!write(0,*)c80

!write(*,*)'ld',lat,delta,dayl,elevation

!-------------------------------------
!calculate daily mean air mass (mbar)

if (dayl == 0.0d0) then

  mbar = m90
  mc   = m90
  ml   = m90

else
  
  !basic setup

  rlat   = pir * lat
  rdelta = pir * delta
  
  sinlat = sin(rlat)
  sindel = sin(rdelta)
  coslat = cos(rlat)
  cosdel = cos(rdelta)

  !------

  !Eqn. 2.5 
  if (abs(lat - delta) < 90.0d0 .and. abs(lat + delta) >= 90.0d0) then
   t1 = 12.0d0
  else
   t1 = (12.0d0 / pi) * acos(-tan(rlat) * tan(rdelta))
  end if
  
  tinv = 1.d0 / t1

  !Eqn. 2.9
  if (abs(lat + delta) >= 90.0d0) then
    Zn = acos(sinlat * sindel - coslat * cosdel) / pir
  else
    Zn = 90.0d0
  end if

  !Eqn. 2.10
  if (abs(lat - delta) >= 90.0d0) then
    Z0 = 90.0d0
  else
    Z0 = lat - delta
  end if
  
  rZ0 = Z0 * pir  !conver to radians
  rZn = Zn * pir
  
  !--------------------------

  b = coslat * cosdel
  
  !write(*,*)'c2',t1,Zn,Z0
  
  if (t1 == 0.d0) then

    mbar = m90
    
  else if (abs(Zn) <= 80.d0) then

    c => c00
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else if (abs(Z0) >= 80.d0) then

    c => c80
    a = c(1) + sinlat * sindel
    mbar = tinv * F(t1,a,b,c)

  else
    
    t80 = 1.d0 / w * acos((cos80 - sinlat * sindel) / (coslat * cosdel)) / pir  !Eqn. 2.8

    c => c00
    a = c(1) + sinlat * sindel
    
    !write(*,*)'crash',t80,a,b,c
    
    tmp1 = F(t80,a,b,c)

    c => c80
    a = c(1) + sinlat * sindel
    tmp2 = F(t1,a,b,c)

    c => c80
    a = c(1) + sinlat * sindel
    tmp3 = F(t80,a,b,c)
    
    mbar = tinv * (tmp1 + tmp2 - tmp3)
    
    !write(0,*)tinv,tmp1,tmp2,tmp3

  end if

  !--------------------------
  !calculate instantaneous air mass at max, mid, and bottom quarter solar zenith angle (m0, mc, ml)
  
  Z = Z0

  cosZ = cos(Z * pir)

  if (Z <= 80.d0) then
    c => c00
  else
    c => c80
  end if
  
  mo = m(cosZ,c)

  !--

  Z = (Z0 + Zn) / 2.d0
  
  cosz = (cos(rZ0) + cos(rZn)) / 2.d0

  if (Z <= 80.d0) then
    c => c00
  else
    c => c80
  end if

  mc = m(cosZ,c)

  !--

  Z = (Z0 + 3.d0 * Zn) / 4.d0
  
  cosz = (cos(rZ0) + 3.d0 * cos(rZn)) / 4.d0

  if (Z <= 80.d0) then
    c => c00
  else
    c => c80
  end if

  ml = m(cosZ,c)

  !--
  
  !write(0,'(a,2f8.3)')'lat,delta    ',lat,delta
  !write(0,'(a,2f8.3)')'t1,t80       ',t1,t80
  !write(0,'(a,2f8.3)')'Zn,Z0        ',Zn,Z0
  !write(0,'(a,4f8.3)')'mbar,mo,mc,ml',mbar,mo,mc,ml
  !read(*,*)

end if

!--------------------------
!correct calculated air mass for elevation

mbar = elev_corr(mbar,elevation)
mo = elev_corr(mo,elevation)
mc = elev_corr(mc,elevation)
ml = elev_corr(ml,elevation)

!------------------------------

end subroutine get_airmass

!------------------------------

real(dp) function m(cosZ,c)

!Instantaneous air mass m, equation 2.1 in Yin, 1997

implicit none

real(dp),               intent(in) :: cosZ
real(dp), dimension(:), intent(in) :: c

m = c(2) / (c(1) + cosZ) + c(3)

end function m

!------------------------------

real(dp) function F(t1,a,b,c)

!integral air mass function F, equation 2.6b in Yin, 1997
!section inside curly braces only - multiply result by 1/t1 to get mbar

implicit none

real(dp),               intent(in) :: t1
real(dp),               intent(in) :: a
real(dp),               intent(in) :: b
real(dp), dimension(:), intent(in) :: c

real(dp) :: wt1
real(dp) :: wpi

real(dp) :: e1
real(dp) :: e2

wpi  = 180.d0 / (pi * w)
wt1  = rw * t1

if (a > b) then
  
  F = wpi * c(2) / sqrt(a**2 - b**2) * acos((b + a * cos(wt1)) / (a + b * cos(wt1))) + c(3) * t1

else if (a < b) then
  
  e1 = sqrt((b + a) * (1.d0 + cos(wt1))) + sqrt((b - a) * (1.d0 - cos(wt1)))
  e2 = sqrt((b + a) * (1.d0 + cos(wt1))) - sqrt((b - a) * (1.d0 - cos(wt1)))
  
  F = wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2) + c(3) * t1

else

  F = wpi * c(2) / a * tan(wt1 / 2.d0) + c(3) * t1

end if

!write(0,*)'F ab ',a,b
!write(0,*)'F X  ',wpi * c(2) / sqrt(b**2 - a**2) * log(e1 / e2)
!write(0,*)'Fc3t1',c(3) * t1
  
end function F

!------------------------------ 

real(dp) function elev_corr(m,elevation)
  
implicit none

real(dp), intent(in) :: m
real(sp), intent(in) :: elevation

elev_corr = m * exp(-elevation / 8000.0d0)

end function elev_corr

!------------------------------

end module airmass
