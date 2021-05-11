subroutine daily(mval,dval,means)

!linear interpolation of monthly to pseudo-daily values
use common_vars, only : dp

implicit none

integer, parameter, dimension(12) :: ndaymo = (/  31,28,31,30,31,30,31,31,30,31,30,31 /)
integer, parameter, dimension(14) :: midday = (/ -15,16,44,75,105,136,166,197,228,258,289,319,350,381 /) !middle day of each month

real,    intent(in),  dimension(12)  :: mval
logical, intent(in) :: means

real(dp), intent(out), dimension(365) :: dval

real(dp), dimension(14) :: emval
real(dp), dimension(13) :: slope

integer :: day
integer :: d
integer :: m
integer :: i

!-------------------------------------------------------
!interpolate the monthly values to daily ones in a cyclical format
!copy last month's data to first and vice versa

do i = 2,13
 emval(i) = mval(i-1)
end do
emval(1) = mval(12)
emval(14) = mval(1)

if (means) then

  !calculate slopes

  do m = 1,13 
    slope(m) = (emval(m+1) - emval(m)) / (midday(m+1) - midday(m))
  end do

  !calculate daily values based on monthly means

  m = 1
  do day = 1,365
    if (day > midday(m+1)) m = m + 1
    d = day - midday(m+1)
    dval(day) = slope(m) * d + emval(m+1)

  end do

else
  
  !distribute the total evenly among the days of the month
  
  day = 1
  do m = 1,12
    do d = 1,ndaymo(m)
      dval(day) = mval(m) / ndaymo(m)
      day = day + 1
    end do
  end do
  
end if

end subroutine daily
