module coordsmod
  
implicit none

public :: parsecoords
public :: calcpixels

integer, parameter :: dp = selected_real_kind(13)

contains

!-----------------------------------------------------

subroutine parsecoords(coordstring,val)

!subroutine to parse a coordinate string

implicit none

character(*),           intent(in)  :: coordstring
real(dp), dimension(4), intent(out) :: val

character(10), dimension(4) :: cval = '0'

integer :: i
integer :: lasti = 1
integer :: part  = 1

do i=1,len_trim(coordstring)
  if (coordstring(i:i) == '/') then
    cval(part) = coordstring(lasti:i-1)
    lasti=i+1
    part = part + 1
  end if
end do

cval(part) = coordstring(lasti:i-1)

read(cval,*)val

if (part < 4) then
  val(3)=val(2)
  val(4)=val(3)
  val(2)=val(1)
end if

end subroutine parsecoords

!-----------------------------------------------------

subroutine calcpixels(infile,bounds,gridref)

use typesizes
use netcdf
use netcdf_error

implicit none

!given infile and lon and lat boundaries, calculate srtx, cntx, srty, and cnty to extract data

!arguments

character(*),           intent(in)  :: infile
real(dp), dimension(4), intent(in)  :: bounds
integer,  dimension(4), intent(out) :: gridref

!local variables

integer :: ncid
integer :: dimid
integer :: varid

integer,  dimension(2) :: xpos
integer,  dimension(2) :: ypos

real(dp), allocatable, dimension(:) :: lon
real(dp), allocatable, dimension(:) :: lat

integer :: xlen
integer :: ylen

real(dp) :: xres
real(dp) :: yres

real(dp), dimension(2) :: lonrange
real(dp), dimension(2) :: latrange


integer, dimension(1) :: pos

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty

integer(2), allocatable, dimension(:) :: ivar

real :: scale_factor
real :: add_offset

integer :: tlen

integer :: i
integer :: j

!-----------------------------------------------------
!open input grid file and retrieve lon and lat vectors

status = nf90_open(infile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)
  
status = nf90_inq_dimid(ncid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

allocate(lon(xlen))

status = nf90_inq_dimid(ncid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

allocate(lat(ylen))

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)
  
status = nf90_get_att(ncid,varid,'actual_range',lonrange)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(ncid,varid,'actual_range',latrange)
if (status /= nf90_noerr) call handle_err(status)

xres = (lonrange(2) - lonrange(1)) / xlen

yres = (latrange(2) - latrange(1)) / ylen

!write(0,*)lonrange
!write(0,*)latrange
write(0,'(a,f6.2,a,f6.2,a)')'grid resolution: ',xres*60.,' x ',yres*60.,' minutes'

!-----------------------------------------------------
!find nearest pixel to minimum coordinates (single point or window)
 
pos = minloc(abs(lon - bounds(1)))
xpos(1) = pos(1)

pos = minloc(abs(lon - bounds(2)))
xpos(2) = pos(1)

pos = minloc(abs(lat - bounds(3)))
ypos(1) = pos(1)

pos = minloc(abs(lat - bounds(4)))
ypos(2) = pos(1)

!write(0,*)'nearest longitude in database:',lon(xpos(1))
!write(0,*)'nearest latitude in database: ',lat(ypos(1))

!-----------------------------------------------------
!calculate the indices of the starting point for the desired pixel, and the number of pixels to be retrieved 

srtx = minval(xpos)
if (lon(srtx) < bounds(1) .and. bounds(2) /= bounds(1)) srtx = srtx + 1
cntx = 1 + abs(maxval(xpos) - srtx)

srty = minval(ypos)
if (lat(srty) < bounds(3) .and. bounds(4) /= bounds(3)) srty = srty + 1
cnty = 1 + abs(maxval(ypos) - srty)

!write(0,*)srtx,srty
!write(0,*)cntx,cnty
 
gridref(1) = srtx
gridref(2) = cntx
gridref(3) = srty
gridref(4) = cnty

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

end subroutine calcpixels

!-----------------------------------------------------

end module coordsmod
