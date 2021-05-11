module read_data

use typesizes
use netcdf
use netcdf_error

implicit none

!shared variables

integer :: ncid
integer :: dimid
integer :: varid

!module subroutines

public :: read_climate
public :: read_soil

contains

!-------------------------------------------

subroutine read_climate(gridref)

use common_vars, only : climatefile,lon,lat,time,elev,temp,prec,sunp,mdays

implicit none

!arguments

integer, dimension(4), intent(in) :: gridref

!parameters

character(10), dimension(4), parameter :: varname = [ 'elv ', 'tmn ', 'pre ', 'sunp' ]

!local variables

integer :: x
integer :: y
integer :: i

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty
integer :: tlen

integer(2), allocatable, dimension(:,:,:) :: var_in

real :: missing
real :: scale_factor
real :: add_offset

!-------------------------------------------

srtx = gridref(1)
cntx = gridref(2)
srty = gridref(3)
cnty = gridref(4)

!open the input file

!this program needs an input file with monthly temp, prec, sun and elevation (m)
!calculated ancillary variables are warmest month and coldest month,
!amount of precipitation in these months,
!and Pjj = 2.0*(p_twm-p_tcm)/(p_twm+p_tcm)

status = nf90_open(climatefile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

write(0,*)'--reading climate-- '!,trim(climatefile),cntx,cnty,tlen

allocate(lon(cntx))
allocate(lat(cnty))
allocate(time(tlen))
allocate(elev(cntx,cnty))
allocate(temp(cntx,cnty,tlen))
allocate(prec(cntx,cnty,tlen))
allocate(sunp(cntx,cnty,tlen))

!-------------------------------------------
!get data

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lon,start=[srtx],count=[cntx])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ncid,varid,lat,start=[srty],count=[cnty])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)
  
status = nf90_get_var(ncid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

allocate(var_in(cntx,cnty,tlen))

do i=1,4

  write(0,*)varname(i)

  status = nf90_inq_varid(ncid,trim(varname(i)),varid)
  if (status /= nf90_noerr) call handle_err(status)

  status = nf90_get_att(ncid,varid,'missing_value',missing)
  if (status /= nf90_noerr) missing = -9999.   !call handle_err(status)

  status = nf90_get_att(ncid,varid,'scale_factor',scale_factor)
  if (status /= nf90_noerr) scale_factor = 1.  !call handle_err(status)

  status = nf90_get_att(ncid,varid,'add_offset',add_offset)
  if (status /= nf90_noerr) add_offset = 0.    !call handle_err(status)

  status = nf90_get_var(ncid,varid,var_in,start=[srtx,srty,1],count=[cntx,cnty,tlen])
  if (status /= nf90_noerr) call handle_err(status)
  
  select case (i)
  case(1)
    elev = real(var_in(:,:,1)) * scale_factor + add_offset
  case(2)
    temp = real(var_in) * scale_factor + add_offset
  case(3)
    prec = real(var_in) * scale_factor + add_offset
  case(4)
    sunp = real(var_in) * scale_factor + add_offset
  end select

end do

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)
  
!-------------------------------------------
!convert precipitation to mm/day

write(0,*)'convert precip to means'

forall (y = 1:cnty)
  forall (x = 1:cntx)
    prec(x,y,:) = prec(x,y,:) / mdays
  end forall
end forall

end subroutine read_climate

!------------------------------------------------------------------------------------------------------------------------

subroutine read_soil(gridref)

use common_vars, only : soilfile,depthfile,tawc,totc,phaq,bulk,maxd

implicit none

!arguments

integer, dimension(4), intent(in) :: gridref

!local variables

integer :: depth

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty

!---------------------------------

srtx = gridref(1)
cntx = gridref(2)
srty = gridref(3)
cnty = gridref(4)

write(0,*)'--reading soils-- '!,trim(soilfile),gridref

status = nf90_open(soilfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'depth',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=depth)
if (status /= nf90_noerr) call handle_err(status)

allocate(maxd(cntx,cnty))
allocate(tawc(cntx,cnty,depth))
allocate(totc(cntx,cnty,depth))
allocate(phaq(cntx,cnty,depth))
allocate(bulk(cntx,cnty,depth))

status = nf90_inq_varid(ncid,'tawc',varid)
if (status /= nf90_noerr) call handle_err(status)
  
write(0,*)'awc'
status = nf90_get_var(ncid,varid,tawc,start=[srtx,srty,1],count=[cntx,cnty,depth])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'totc',varid)
if (status /= nf90_noerr) call handle_err(status)
  
write(0,*)'totc'
status = nf90_get_var(ncid,varid,totc,start=[srtx,srty,1],count=[cntx,cnty,depth])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'phaq',varid)
if (status /= nf90_noerr) call handle_err(status)
  
write(0,*)'phaq'
status = nf90_get_var(ncid,varid,phaq,start=[srtx,srty,1],count=[cntx,cnty,depth])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'bulk',varid)
if (status /= nf90_noerr) call handle_err(status)
  
write(0,*)'bulk'
status = nf90_get_var(ncid,varid,bulk,start=[srtx,srty,1],count=[cntx,cnty,depth])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)
  
!-----

status = nf90_open(depthfile,nf90_nowrite,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ncid,'z',varid)
if (status /= nf90_noerr) call handle_err(status)
  
write(0,*)'depth'
status = nf90_get_var(ncid,varid,maxd,start=[srtx,srty],count=[cntx,cnty])
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)
  
end subroutine read_soil

!---------------------------------

end module read_data
