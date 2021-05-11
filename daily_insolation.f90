module daily_insolation

use common_vars, only : dp

integer :: nd      !day number of the year

real(dp) :: phi     !latitude in degrees

real(dp) :: ww      !top of the atmosphere insolation
real(dp) :: dayl    !day length in hours
real(dp) :: delta   !declination in degrees

contains  

subroutine dayins(ecc,perh,xob,phi,nd,ww,dayl,delta)

!this subroutine calculates top-of-the-atmosphere insolation given orbital
!parameters, day of the year, and latitude
!output : ww=ly/day  or  kj m(-2) day(-1)  dayl=length of day (hours)

use common_vars, only : dp,pi,pir

implicit none

!arguments

integer, intent(in) :: nd

real(dp), intent(in) :: ecc
real(dp), intent(in) :: perh
real(dp), intent(in) :: xob
real(dp), intent(in) :: phi

real(dp), intent(out) :: ww      !top of the atmosphere insolation
real(dp), intent(out) :: dayl    !day length in hours
real(dp), intent(out) :: delta   !declination in degrees

!parameters

real(dp), parameter :: ss   = 1353.0d0
real(dp), parameter :: tau  =   86.4d0
real(dp), parameter :: test =    1.0e-8
real(dp), parameter :: step =  360.0d0/365.25d0

!variables

real(dp) :: xec
real(dp) :: dlam

real(dp) :: sf
real(dp) :: so
real(dp) :: xl
real(dp) :: xllp
real(dp) :: xee
real(dp) :: xse
real(dp) :: xlam
real(dp) :: dlamm
real(dp) :: anm
real(dp) :: ranm
real(dp) :: ranv
real(dp) :: anv
real(dp) :: tls

real(dp) :: rphi
real(dp) :: rau
real(dp) :: s
real(dp) :: rlam
real(dp) :: sd
real(dp) :: cd
real(dp) :: rdelta  !declination in radians
real(dp) :: spa
real(dp) :: cp
real(dp) :: aphi    !absolute value of the latitude
real(dp) :: adelta  !absolute value of the solar zenith angle
real(dp) :: tt

real(dp) :: at
real(dp) :: spd
real(dp) :: tp
real(dp) :: stp
real(dp) :: rdayl

!----------------------------------------------------------------------

sf  = tau*ss/pi
so  = sin(xob*pir)
xl  = perh+180.0d0

xllp = xl*pir
xee  = ecc*ecc
xse  = sqrt(1.0d0-xee)
xlam = (ecc/2.0d0+ecc*xee/8.0d0)*(1.0d0+xse)*sin(xllp)-xee/4.0d0*(0.5d0+xse) &
        *sin(2.0d0*xllp)+ecc*xee/8.0d0*(1.0d0/3.0d0+xse)*sin(3.0d0*xllp)

xlam = 2.0d0*xlam/pir
dlamm= xlam+(nd-80)*step
anm  = dlamm-xl

ranm = anm*pir
xec  = xee*ecc
ranv = ranm+(2.0d0*ecc-xec/4.0d0)*sin(ranm)+5.0d0/4.0d0*ecc*ecc*   &
             sin(2.0d0*ranm)+13.0d0/12.0d0*xec*sin(3.0d0*ranm)

anv  = ranv/pir
tls  = anv+xl

dlam = tls

!-----------------------------------------------

!----------------------------------------------------------------------------------------------------

rphi    =  phi * pir
ranv    =  (dlam - xl) * pir
rau     =  (1.0d0 - ecc * ecc) / (1.0d0 + ecc * cos(ranv))

s       =  sf / rau / rau
rlam    =  dlam * pir
sd      =  so * sin(rlam)
cd      =  sqrt(1.0d0 - sd * sd)

rdelta  =  atan(sd / cd)
delta   =  rdelta / pir
spa     =  sd * sin(rphi)

cp      =  cd * cos(rphi)
aphi    =  abs(phi)
adelta  =  abs(delta)

!singularity for aphi = 90 and delta = 0
!particular cases for lat = 0 or delta = 0

tt = abs(aphi - 90.0d0)

if (tt <= test .and. adelta <= test) then

  dayl = 0.00d0
  ww = 0.00d0

else if (adelta <= test) then

  dayl = 12.0d0
  ww = s * cos(rphi)

else if (aphi <= test) then

  dayl = 12.0d0
  ww = s * cos(rdelta)

else

  at = 90.0d0 - adelta
  spd = phi * delta

  if (aphi <= at) then

    tp = -spa / cp
    stp = sqrt(1.0d0 - tp * tp)
    rdayl = acos(tp)
    dayl = 24.0d0 * rdayl / pi
    ww = s * (rdayl * spa + cp * stp)

  else if (spd > 0.0d0) then

    dayl = 24.00d0
    ww = s * spa * pi

  else if (spd < 0.0d0) then

    dayl = 0.00d0
    ww = 0.00d0

  else

    tp =  - spa / cp
    stp = sqrt(1.0d0 - tp * tp)
    rdayl = acos(tp)
    dayl = 24.0d0 * rdayl / pi
    ww = s * (rdayl * spa + cp * stp)

  end if
end if

!------------------------------------

end subroutine dayins

end module daily_insolation
