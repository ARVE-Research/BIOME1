program biome1

!Jed Kaplan, EPFL, 2009

use common_vars,      only : sp,dp,files,lon,lat,time,elev,temp,prec,sunp,tawc,phaq,bulk,mdays,maxd
use typesizes
use netcdf
use netcdf_error
use orbit
use radpet,           only : rad_and_pet
use aetmod,           only : m,aet_alpha
use read_data,        only : read_climate,read_soil
use crop_suitability, only : ph_soil, c_density, crop_suit
use biome1mod,        only : pftpar,pftpars,calcbiome,npfts
use coordsmod,        only : parsecoords,calcpixels

implicit none

!-----------------------------------------------------------
!variables

character(100) :: outfile

real(dp), dimension(3) :: orbitvars

integer :: x,y

integer, parameter :: tlen = 12

real(sp) :: rad0

real(sp), dimension(12) :: mdayl   !monthly mean day length (h)
real(sp), dimension(12) :: mrad    !monthly mean surface sw (W m-2)
real(sp), dimension(12) :: mpet    !monthly total PET (mm)
real(sp), dimension(12) :: fdiff   !diffuse fraction
real(sp), dimension(12) :: malpha  !monthly mean alpha

real(sp) :: drought                !integrated drought scalar, same units as alpha
real(sp) :: aalpha                 !mean annual alpha
real(dp) :: ph_avg                 !average pH over the soil column
real(dp) :: c_dens                 !carbon density of soil
real(dp) :: S_index                !index of land suitability to cultivation
real(dp) :: P_index                !index of land suitability to pasture

integer :: i,j
integer :: dyr  !day of year
integer :: mon  !month counter
integer :: dmo  !day of month counter
integer :: yr   !year counter for AET loop

real(dp) :: a1  !previous year's value of dalpha(1) to control iteration loop
real(dp) :: m1

real,       dimension(2) :: range
integer(2), dimension(2) :: irange

integer :: ncid
integer :: varid

integer :: srtx
integer :: srty
integer :: cntx
integer :: cnty

integer, dimension(4) :: gridref

real(dp) :: awc  !total available water content

character(80) :: argfile

!--------------------------------------------------------------------
!daily vectors

real(dp), dimension(365) :: dtemp  !daily temperature (C)
real(dp), dimension(365) :: dprec  !daily precipitation (mm)
real(dp), dimension(365) :: dsunp  !daily bright sunshine (fraction)
real(dp), dimension(365) :: ddayl  !daily dacntygth (h)
real(dp), dimension(365) :: dpet   !daily potential evapotranspiration (mm)
real(dp), dimension(365) :: daet   !daily actual evapotranspiration (mm)
real(dp), dimension(365) :: dalpha !daily alpha ratio

!output variables

real(sp), allocatable, dimension(:,:)   :: gdd5
real(sp), allocatable, dimension(:,:)   :: gridcrop  !crop suitability index
real(sp), allocatable, dimension(:,:)   :: gridpasture
real(sp), allocatable, dimension(:,:)   :: gridph
real(sp), allocatable, dimension(:,:)   :: gridcarbon
real(sp), allocatable, dimension(:,:)   :: g_annual_alpha
real(sp), allocatable, dimension(:,:,:) :: gridalpha

integer(2), allocatable, dimension(:,:) :: biome

real(sp) :: tcm
real(sp) :: twm
real(sp) :: gdd0

real, dimension(2) :: sdepth

character(80) :: coordstring
real(dp), dimension(4) :: bounds

character(80) :: note

logical, dimension(10000) :: biomenum
logical, dimension(npfts) :: present
integer :: val
integer :: pft

real :: aprec
real :: apet

!-------------------------------------------------------------------
!calculate orbital parameters for this year

call orbital_parameters(0.0d0,ecc,pre,perh,xob)

!write(0,*)'***************************************'
!write(0,'(A,F8.4,A)')'time         = ',ka,' ka'
!write(0,'(A,F8.4)')  'eccentricity = ',ecc
!write(0,'(A,F8.4)')  'prec. param. = ',pre
!write(0,'(A,F8.4)')  'long. perh.  = ',perh
!write(0,'(A,F8.4)')  'obliquity    = ',xob
!write(0,*)'***************************************'

orbitvars = (/xob,perh,ecc/)

!-------------------------------------------
!get input file names

call getarg(1,argfile)

open(10,file=argfile)

read(10,nml=files)

close(10)

!-------------------------------------------
!get bounding coordinates for run

call getarg(2,coordstring)

call parsecoords(coordstring,bounds)

call calcpixels(climatefile,bounds,gridref)

!write(0,*)bounds,gridref

srtx = gridref(1)
cntx = gridref(2)
srty = gridref(3)
cnty = gridref(4)

!-------------------------------------------
!create output file

call getarg(3,outfile)

write(0,'(2a,i5,a,i5)')'generating ',trim(outfile),cntx,' x ',cnty

! NB a new routine for generating an appropriate output file needs to be created
! call genoutputfile(trim(outfile),cntx,cnty)

!-------------------------------------------
!get climate data

call read_climate(gridref)

!-------------------------------------------
!get soils data 

call read_soil(gridref)

!-------------------------------------------
!initialize pft parameters

open(99,file='pftpars.namelist')
read(99,nml=pftpars)

!write(0,*)pftpar

!-------------------------------------------
!allocate output grids and initialize to null value

allocate(gdd5(cntx,cnty))
allocate(gridcrop(cntx,cnty))
allocate(gridpasture(cntx,cnty))
allocate(gridcarbon(cntx,cnty))
allocate(gridph(cntx,cnty))
allocate(g_annual_alpha(cntx,cnty))
allocate(gridalpha(cntx,cnty,tlen))
allocate(biome(cntx,cnty))

gdd5            = -9999.
gridpasture     = -9999.
gridcarbon      = -9999.
gridph          = -9999.
g_annual_alpha  = -9999.
gridalpha       = -9999.
gridcrop        = -9999.
biome           = -9999.

biomenum = .false.
!-------------------------------------------
!calculate the climate & soil variables

!start of grid cell loop
do y = 1,cnty

  write(note,'(a,i5,a,i5)')' calculating row', y,' of ',cnty
  call overprint(note)

  do x = 1,cntx
  
    !skip to the next gridcell if the input data is invalid (typically ocean gridpoints)

    if (elev(x,y) == -9999. .or. maxval(tawc(x,y,:)) <= 0.) then
      cycle
    else if (sunp(x,y,1) < 0. .or. prec(x,y,1) < 0.) then
      cycle
    end if
    
    !interpolate monthly vector to daily values
    
    call daily(temp(x,y,:),dtemp,.true.)
    call daily(prec(x,y,:),dprec,.true.)
    call daily(sunp(x,y,:),dsunp,.true.)
    
    !calculate growing degree days
    gdd5(x,y) = sum(dtemp - 5.d0,mask=dtemp > 5.d0)

    !write(0,*)x,y
    !write(0,*)lat(y),elev(x,y)
    !write(0,*)prec(x,y,:)

    call rad_and_pet(real(lat(y)),elev(x,y),dtemp,dprec,dsunp,orbitvars,mdayl,mrad,ddayl,dpet,mpet,rad0,fdiff)
    
    !----

    !calculate total available soil water content for calculation of alpha
    !awc is total water holding capacity of the soil
        
    !convert soil dataset layer-specific cm m-1 to total mm
    
    !awc = sum(tawc(x,y,:) * 2.,mask = tawc(x,y,:) > 0.)  !NB this statement assumes two layers of 20cm each
    
!    write(0,*)x,y,maxd(x,y)
    
    if (maxd(x,y) > 30) then
      sdepth(1) = 0.3
      sdepth(2) = 0.01 * (real(maxd(x,y)) - 30.)
    else if (maxd(x,y) > 0) then
      sdepth(1) = 0.01 * real(maxd(x,y))
      sdepth(2) = 0.
    else
      sdepth(1) = 0.3
      sdepth(2) = 0.7
    end if
    
    if (tawc(x,y,1) > 0.) then
      awc = 10. * sum(tawc(x,y,:) * sdepth)
    else
      awc = 0.
    end if
    
    if (awc > 0.d0) then
    
      m = awc  !for first day of first simulation year, set m (soil moisture) to awc
      m1 = m
      a1 = 0.d0 !a1 = previous year's value of alpha

      do !yr = 1,25        !run aet_alpha for up to 25 years to get equilibrium soil moisture conditions

        do dyr = 1,365
        
          call aet_alpha(awc,ddayl(dyr),dprec(dyr),dpet(dyr),daet(dyr),dalpha(dyr))
          !write(0,*)dyr,m
        end do

        !write(0,*)'test',yr,m,m1

        if (abs(m - m1) <= 1.) exit

        m1 = m
        a1 = dalpha(1)

      end do
      
      !create monthly means of alpha
      dyr = 1
      do mon=1,12             !month loop
        malpha(mon) = 0.d0
        do dmo = 1,int(mdays(mon))   !day of month loop

          malpha(mon) = malpha(mon) + dalpha(dyr) / mdays(mon)

          dyr = dyr + 1

        end do
      end do

    else
      dalpha = -9999.d0
      malpha = -9999.
    end if
    
    !create mean annual alpha
!    if (count(dtemp > 0.) > 0) then
!      aalpha = sum(dalpha,mask=dtemp > 0.) / count(dtemp > 0.)
!    else
      aalpha = sum(dalpha) / size(dalpha)
!    end if
    
    !alternative: calculate aalpha as integrated drought scalar only for days with alpha < 1
    
    if (count(dalpha < 1.) > 0.) then
      drought = 1. - sum(1. - dalpha) / count(dalpha < 1.)
    else
      drought = 1.
    end if
    
    !aalpha = drought

    !write(0,*)aalpha,drought
    !write(0,*)sum(dprec),sum(dpet)
    !read(*,*)

    !----

    g_annual_alpha(x,y)  = aalpha
    gridalpha(x,y,:) = malpha
    
    !soil variables

    call ph_soil(phaq(x,y,:),ph_avg)
    
    !gridph(x,y) = ph_avg
    gridph(x,y) = sum(dpet)  !store PET in gridpH temporarily

    call c_density(totc(x,y,:),bulk(x,y,:),c_dens)
    
    !gridcarbon(x,y) = c_dens

    !-----------------------------------------------------------------------
    !get index of land suitability for cultivation, S_index

    call crop_suit(gdd5(x,y),aalpha,ph_avg,c_dens,S_index,P_index)
    
    gridcrop(x,y) = S_index
    gridpasture(x,y) = P_index
    
    tcm = minval(temp(x,y,:))
    twm = maxval(temp(x,y,:))
    gdd0 = sum(dtemp,mask = dtemp > 0.)
    
    aprec = sum(dprec)
    apet = sum(dpet)
    
    call calcbiome(tcm,gdd5(x,y),gdd0,twm,aalpha,aprec,apet,biome(x,y))
    
    biomenum(biome(x,y)) = .true.
    
    !write(0,*)x,y,biome(x,y)
    !read(*,*)
    
    
   gridcarbon(x,y) = aprec - apet

    
  end do
end do

write(0,*)

write(0,*)
write(0,*)'number of unique biomes: ',count(biomenum)

j = 1

do i = 1,size(biomenum)
  if (biomenum(i)) then
    
    val = i
    
    do pft = npfts,1,-1
  
      if (mod(val,2) == 0) then
        present(pft) = .false.
      else
        present(pft) = .true.
      end if
  
      val = val / 2
  
    end do
    
    write(0,'(2i5,13l3)')i,j,present

    !where (biome == i) biome = j
    j = j + 1

  end if
end do
write(0,*)

!end of grid loop
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!write output

status = nf90_open(outfile,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_inq_varid(ncid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)
  
status = nf90_put_var(ncid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(lon)
range(2) = maxval(lon)

status = nf90_put_att(ncid,varid,'actual_range',range)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_varid(ncid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(lat)
range(2) = maxval(lat)

status = nf90_put_att(ncid,varid,'actual_range',range)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_varid(ncid,'times',varid)
if (status /= nf90_noerr) call handle_err(status)
  
status = nf90_put_var(ncid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_inq_varid(ncid,'gdd5',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,gdd5)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(gdd5,mask=gdd5/=-9999.)
range(2) = maxval(gdd5,mask=gdd5/=-9999.)

write(0,*)'gdd5',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)

!----
status = nf90_inq_varid(ncid,'crops',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,gridcrop)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(gridcrop,mask=gridcrop/=-9999.)
range(2) = maxval(gridcrop,mask=gridcrop/=-9999.)

write(0,*)'gridcrop',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_inq_varid(ncid,'pasture',varid)
if (status /= nf90_noerr) call handle_err(status)

gridpasture = minval(temp,dim=3)

where (gridpasture == -999.9) gridpasture = -9999.

status = nf90_put_var(ncid,varid,gridpasture)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(gridpasture,mask=gridpasture/=-9999)
range(2) = maxval(gridpasture,mask=gridpasture/=-9999)

write(0,*)'tcm',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_inq_varid(ncid,'soil_carbon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,gridcarbon)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(gridcarbon,mask=gridcarbon/=-9999.)
range(2) = maxval(gridcarbon,mask=gridcarbon/=-9999.)

write(0,*)'carbon',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_inq_varid(ncid,'alpha',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,gridalpha)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(gridalpha,mask=gridalpha/=-9999.)
range(2) = maxval(gridalpha,mask=gridalpha/=-9999.)

write(0,*)'alpha',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)

!----
status = nf90_inq_varid(ncid,'annual_alpha',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,g_annual_alpha)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(g_annual_alpha,mask=g_annual_alpha/=-9999.)
range(2) = maxval(g_annual_alpha,mask=g_annual_alpha/=-9999.)

write(0,*)'ann alpha',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)
  
!----

status = nf90_inq_varid(ncid,'soil_ph',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,gridph)
if (status /= nf90_noerr) call handle_err(status)

range(1) = minval(gridph,mask=gridph/=-9999.)
range(2) = maxval(gridph,mask=gridph/=-9999.)

write(0,*)'pH',range

status = nf90_put_att(ncid,varid,'valid_range',range)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_inq_varid(ncid,'biome',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,biome)
if (status /= nf90_noerr) call handle_err(status)

irange(1) = minval(biome,mask=biome/=-9999)
irange(2) = maxval(biome,mask=biome/=-9999)

write(0,*)'biome',irange

status = nf90_put_att(ncid,varid,'valid_range',irange)
if (status /= nf90_noerr) call handle_err(status)

!----

write(0,*)'closing output file'

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

!---------------------------------------------------------

end program biome1
