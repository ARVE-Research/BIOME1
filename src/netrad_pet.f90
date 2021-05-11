module netrad_pet

use common_vars, only : sp,dp
use esatdesdT

real(dp) :: lw_rad
real(dp) :: pet

public :: rl_pet

contains

subroutine rl_pet(temp,sun,dayl,sw_rad,lw_rad,pet)

!This code is based on the paper:
!A. Haxeltine and Prentice, I.C., BIOME3..., Glob. Biogeochem. Cycles, 10, 693-709
!With a new calculations of:
!downwelling longwave   (Josey et al., 2003. J. Geophys. Res., 108(C4), 3108, doi:10.1029/2002JC001418)
!dEsat/dT               (Oleson et al., 2004, CLM 3.0 technical note)
!lvap                   (Henderson-Sellers, 1984. Quart. J. R. Met. Soc., 110, 1186-1190)
!Other references:
!Linacre (1968) Agr. Meteorol., 5, 49-63
!Prentice et al. (1993) Ecological Modelling, 65, 51-70.
!Jed Kaplan, EPFL, 2008

implicit none

!arguments

real(dp), intent(in) :: temp      !surface air (2m) temperature (C)
real(dp), intent(in) :: sun       !bright sunshine duration (fraction)
real(dp), intent(in) :: dayl      !daylength (h)
real(dp), intent(in) :: sw_rad    !downwelling shortwave radiation (kJ m-2 d-1)

real(dp), intent(out) :: lw_rad   !longwave radiation (kJ m-2 d-1)
real(dp), intent(out) :: pet      !potential evapotranspiration (mm d-1)

!local variables

real(dp) :: Tk       !surface air temperature (K)
real(dp) :: albedo   !surface albedo (fraction)

real(dp) :: lvap     !Latent heat of vaporization of water (kJ kg-1)
real(dp) :: gamma    !psychrometer constant (Pa K-1)
real(dp) :: ss       !rate of increase of saturated vapor pressure with temperature (Pa K-1)
real(dp) :: f        !Linacre parameter (function of sunshine fraction)

real(dp) :: netrad   !net radiation (kJ m-2 d-1)

real(dp) :: Ql     !net longwave radiation (W m-2)
real(dp) :: Ql_up  !upwelling longwave radiation (W m-2)
real(dp) :: Ql_dn  !downwelling longwave radiation (W m-2)

real(dp) :: Ts     !surface temperature (K)

real(dp) :: n  !cloud fraction
real(dp) :: Tdew  !dew point temperature (K)
real(dp) :: D  !dew point depression (K)
real(dp) :: es !saturation vapor pressure

!parameters

real(dp), parameter :: sb = 5.6704e-8  !Stefan-Bolzmann constant (W m-2 K-4)
real(dp), parameter :: e  = 0.98d0     !emissivity ()
real(dp), parameter :: al = 0.045d0    !longwave reflectivity (albedo), Josey et al., pg 5-9

real(dp), parameter :: a  =  10.77d0   !parameters in Josey et al.
real(dp), parameter :: b  =   2.34d0
real(dp), parameter :: c  = -18.44d0

real(dp), parameter :: cs = 1.5d0 !shape parameter for the curve relating fractional cloud cover to fractional sunshine duration

!-------------------------------------------------

Tk = 273.0d0 + temp

albedo = 0.17d0

!calculate gamma, lvap

gamma = 65.05d0 + temp * 0.064d0  !psychrometer constant

lvap = 0.001 * 1.91846e6 * (Tk / (Tk - 33.91d0))**2  !(kJ kg-1) Eqn. from Henderson-Sellers (1984)

ss = desdT(Tk)

f = 0.2d0 + 0.8d0 * sun  !Linacre Eqn. 7

n = (cs - cs * sun) / (cs + sun)  !equation based on simple nonlinear fit to cloud vs. sunp (and vv). data from TMY3 database

!-------------------------------------------------
!calculate longwave radiation

!Ql = 697.8d0 * (0.245 - 0.158e-10 * Tk**4) * f         !Linacre Eqn. 23, converted to SI units (W m-2)

!write(*,'(a10,f9.3)')'linacre',Ql

!Ql = (0.2d0 + (1.d0 - 0.2d0) * sun) * (107.d0 - temp)  !Prentice et al. Eqn. 11 (W m-2)

!write(*,'(a10,f9.3)')'prentice',Ql

!----

Ts = Tk !approximation that mean daily surface temperature equals air temp.

Ql_up = e * sb * Ts**4                                          !black body upwelling longwave (W m-2)  !various sources e.g., Oleson et al.

!--

!Josey formulation for downwelling longwave

!Ql_dn = sb * (Tk + a*n**2 + b*n + c)**4                        !downwelling longwave (simplified) (W m-2) Josey et al. Eqn. 9,J1

!--

es = 0.01 * esat(Tk)   !saturation vapor pressure (mbar)

Tdew = 34.07d0 + 4157.d0 / log(2.1718e8 / es)  !Josey et al., Eqn. 10

D = Tdew - Tk

Ql_dn = sb * (Tk + a*n**2 + b*n + c + 0.84d0 * (D + 4.01d0))**4  !downwelling longwave (W m-2) Josey et al. Eqn. 14,J2

Ql = Ql_up - (1.d0 - al) * Ql_dn   !Josey et al., Eqn 1

!write(*,'(a10,3f9.3)')'josey',Ql,Ql_up,Ql_dn

!----

lw_rad = 0.001 * 3600.d0 * dayl * Ql  !daytime net longwave (kJ m-2 d-1)

!----

netrad = (1. - albedo) * sw_rad - lw_rad  !(kJ m-2 d-1)

!----
  
pet = max((ss / (ss + gamma)) * netrad / lvap, 0.d0)

!if (pet > 5.d0) then
!  write(0,*)'sun, cloud ',sun,n
!  write(0,*)'SW rad     ',sw_rad,dayl
!  write(0,*)'LW rad     ',lw_rad,Ql,Ql_up,Ql_dn
!  write(0,*)'gamma,ss   ',gamma,ss
!  write(0,*)'netrad, hs (W m-2)',netrad,1000.d0 * netrad / (3600.d0 * dayl)
!  write(0,*)'lvap, pet',lvap,pet
!  read(*,*)
!end if

!----

end subroutine rl_pet

!---------------------------

end module netrad_pet
