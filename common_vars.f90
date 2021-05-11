module common_vars

implicit none

!type

integer, parameter :: sp = selected_real_kind(4)
integer, parameter :: dp = selected_real_kind(13)

!basic constants

real(dp), parameter :: pi  = 3.14159265358979323846d0 !26433 83279 50288 41971 69399 37510 (unitless)
real(dp), parameter :: pir = pi/180.0d0
real(dp), parameter :: Tfreeze = 273.15d0
real(dp), parameter :: e = 2.718281828459045d0

!shared variables

real(dp), allocatable, dimension(:) :: time

real(dp), allocatable, target, dimension(:)     :: lon
real(dp), allocatable, target, dimension(:)     :: lat

real(sp), allocatable, target, dimension(:,:)   :: elev

real(sp), allocatable, target, dimension(:,:,:) :: temp
real(sp), allocatable, target, dimension(:,:,:) :: prec
real(sp), allocatable, target, dimension(:,:,:) :: sunp

integer(2), allocatable, target, dimension(:,:)   :: maxd   !Reynolds FAO soil depth
real(sp),   allocatable, target, dimension(:,:,:) :: tawc   !water holding capacity from FAO-WISE (cm m-1)
real(sp),   allocatable, target, dimension(:,:,:) :: totc
real(sp),   allocatable, target, dimension(:,:,:) :: phaq
real(sp),   allocatable, target, dimension(:,:,:) :: bulk

real, dimension(12), parameter :: mdays = (/ 31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)

character(100) :: climatefile
character(100) :: soilfile
character(100) :: depthfile

namelist /files/ climatefile,soilfile,depthfile

end module common_vars
