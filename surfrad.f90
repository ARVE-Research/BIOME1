module surfrad

use common_vars, only : dp,sp,pi

!module shared variables

real(dp) :: direct
real(dp) :: diffuse
real(sp) :: tcm

!module subroutines

public :: surface_rad

contains

subroutine surface_rad(r0,sun,dayl,mbar,mo,mc,ml,prec,Pjj,tcm,p,pet,direct,diffuse)

!This code is based on the paper:
!X. Yin (1998) Temporally-aggregated atmospheric optical properties as a function of common climatic information:
!Systems development and application, Meteorol. Atmos. Phys. 68, 99-113
!Jed Kaplan, EPFL, 2008

implicit none

!-----------------------------------
!arguments

real(dp), intent(in)  :: Pjj  !precipitation equitability index
real(dp), intent(in)  :: dayl !daylength (hr)
real(dp), intent(in)  :: mbar !mean daily air mass
real(dp), intent(in)  :: mc   !airmass at medium cosine Z angle
real(dp), intent(in)  :: ml   !airmass at bottom-quarter cos Z angle
real(dp), intent(in)  :: mo   !airmass at max cosine zenith angle
real(dp), intent(in)  :: pet  !potential evapotranspiration mm/day
real(dp), intent(in)  :: r0   !top-of-atmospere insolation (kJ m-2 d-1)
real(sp), intent(in)  :: p    !relative atmospheric pressure (1=sea level)
real(dp), intent(in)  :: prec !precipitation mm/day
real(dp), intent(in)  :: sun  !bright sunshine duration fraction, n/N (fraction)
real(sp), intent(in)  :: tcm  !mean annual temperature (used as tropics indicator)

real(dp), intent(out) :: direct  !direct-beam downwelling shortwave (kJ m-2 d-1)
real(dp), intent(out) :: diffuse !diffuse downwelling shortwave (kJ m-2 d-1)

!-----------------------------------
!variables

real(dp) :: tau   !direct insolation atmospheric turbidity factor
real(dp) :: zeta0 !diffuse insolation atmospheric turbidity factor
real(dp) :: x     !tropics indicator (tropical = 1, else 0)
real(dp) :: fm    !atmospheric transmittance function

real(dp) :: ag = 0.17d0   !Surface shortwave albedo (average=0.17)

real(dp) :: j2w
real(dp) :: fdif
real(dp) :: stmp

!-----------------------------------
!parameters

real(dp), parameter :: kp  = 0.500d0 !links absorption coeff. to trans. coeff.
real(dp), parameter :: kag = 3.300d0
real(dp), parameter :: kan = 2.320d0
real(dp), parameter :: kn  = 0.686d0 !cloud parameter

!----------------------------------------------------------------------------

if (tcm < 10.0d0) then
  x = 0.0d0
else if (tcm > 20.0d0) then
  x = 1.0d0
else
  x = sin(pi / 2.0d0 * (tcm / 10.0d0 - 1.0d0))
end if

tau = exp(-0.115d0 * p * ((2.15d0 - 0.713d0 * x + exp(-6.74d0 / (prec + 1.0d0))) * exp(0.0971d0 * pet) - 0.650d0 * (1.0d0 - x) * Pjj))  !Eqn. 4.1

fm = 0.01452d0 * (mbar + ml) * exp(1.403d0 * tau) - 0.1528d0 * mo + mc + 0.48700d0 * (mc - ml) + 0.2323d0   !Eqn. 2.4 2nd term

direct = sun * tau**kp * r0 * tau**fm   !Eqn. 2.4

zeta0 = 0.503d0 * exp(-1.20d0 * p * exp(0.633d0 / (prec + 1.0d0) - 0.226d0 * pet)) * 3.300d0**ag * 2.32d0**(1.0d0 - sun) * (1.0d0 - 0.686d0 * (1.0d0 - sun))  !Eqn. 4.2

diffuse = zeta0 * kag**ag * kan**(1.0d0 - sun) * (1 - kn * (1.0d0 - sun)) * (tau**kp * r0 - direct)   !Eqn. 2.5

!----
!alternative calculation here as a test

!fdif = diffuse / (diffuse + direct)

!stmp = r0*(0.25d0+0.5d0*sun)

!direct  = stmp * (1.d0 - fdif)
!diffuse = stmp * fdif

!----

!if ((direct + diffuse) > 5000.) then

  !j2w = 1000.d0 / (3600.d0 * dayl)

  !write(0,*)r0*j2w,sun,dayl
  !write(0,*)mbar,mo,mc,ml
  !write(0,*)prec,Pjj,pet,p
  !write(0,*)tau,zeta0
  !write(0,*)direct*j2w,diffuse*j2w,j2w*(direct+diffuse)
  !write(0,*)'surface sw by Angstrom-Prescott method:',j2w*r0*(0.25d0+0.5d0*sun)
  !read(*,*)

!end if

end subroutine surface_rad

end module surfrad
