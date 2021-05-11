module netcdf_error

integer :: status

contains

!-------------------------------------

subroutine handle_err(status)

! NOTE always add this subroutine to programs
! that use the NetCDF libraries. It handles the
! errors and you will be much happier you did.      

use netcdf

implicit none
      
integer, intent(in) :: status
      
if (status/=nf90_noerr) then
  write(0,*)nf90_strerror(status)
  stop
end if

end subroutine handle_err

end module netcdf_error
