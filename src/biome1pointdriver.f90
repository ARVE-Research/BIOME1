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

integer :: x
integer :: y

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

real :: tmin
real :: val1
real :: val2

character(2) :: sitenum

character(100) :: soilpath = '/Volumes/arve/shared/med_soilvegchange/awc_vcf_orgc_out/site' 



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
!initialize pft parameters

open(99,file='pftpars.namelist')
read(99,nml=pftpars)

!write(0,*)pftpar

!-------------------------------------------
!get input file names

call getarg(1,climatefile)

!-------------------------------------------
!open input files

open(10,file=climatefile,status='old')

!-------------------------------------------------------------------
!allocate arrays

allocate(lon(1))
allocate(lat(1))
allocate(elev(1,1))
allocate(temp(1,1,12))
allocate(prec(1,1,12))
allocate(sunp(1,1,12))

!-------------------------------------------
!get climate and standard soils data

do

  read(10,*,end=99)sitenum,lon,lat,elev,tmin,temp,prec,sunp
  
  prec(1,1,:) = prec(1,1,:) / mdays

  !figure out name of scenario file

  soilfile=trim(soilpath)//sitenum//'.dat'
    
  open(20,file=soilfile,status='old')

  !-------------------------------------------
  !scan scenario soils dataset once to count number of scenarios

  cntx = 1
  do
    read(20,*,end=5)
    cntx = cntx + 1
  end do

  5 continue

  rewind(20)
  
  write(0,*)'working on site',sitenum
  
  !open output file for this site
  
  outfile='cforg/site'//sitenum//'.dat'		
  
  open(30,file=outfile,status='unknown')
  
  !-------------------------------------------
  !allocate output grids and initialize to null value

  cnty = 1

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
  !start the calculation loop for this site
  
  x = 1
  y = 1
  
  do

    read(20,*,end=98)awc,val1,val2	
    
    !write(0,*)awc,temp
    !read(*,*)

    !interpolate monthly vector to daily values

    call daily(temp(1,1,:),dtemp,.true.)
    call daily(prec(1,1,:),dprec,.true.)
    call daily(sunp(1,1,:),dsunp,.true.)

    !calculate growing degree days
    gdd5(1,1) = sum(dtemp - 5.d0,mask=dtemp > 5.d0)

    !write(0,*)1,1
    !write(0,*)lat(y),elev(1,1)
    !write(0,*)prec(1,1,:)

    call rad_and_pet(real(lat(y)),elev(1,1),dtemp,dprec,dsunp,orbitvars,mdayl,mrad,ddayl,dpet,mpet,rad0,fdiff)

    !----

    m = awc  !for first day of first simulation year, set m (soil moisture) to awc
    m1 = m
    a1 = 0.d0 !a1 = previous year's value of alpha

    do !yr = 1,25        !run aet_alpha for up to 25 years to get equilibrium soil moisture conditions

      do dyr = 1,365

        call aet_alpha(awc,ddayl(dyr),dprec(dyr),dpet(dyr),daet(dyr),dalpha(dyr))
        !write(0,*)dyr,m
      end do

      !write(0,*)'test',m,m1

      if (abs(m - m1) <= 1.) exit

      m1 = m
      a1 = dalpha(1)

    end do

    !create monthly means of alpha
    dyr = 1
    do mon=1,12             !month loop
      malpha(mon) = 0.d0
      do dmo=1,mdays(mon)   !day of month loop

        malpha(mon) = malpha(mon) + dalpha(dyr) / mdays(mon)

        dyr = dyr + 1

      end do
    end do

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

    !write(0,*)'hak',aalpha,drought
    !write(0,*)sum(dprec),sum(dpet)
    !read(*,*)

    !----

    g_annual_alpha(1,1)  = aalpha
    gridalpha(1,1,:) = malpha

    tcm = minval(temp(1,1,:))
    twm = maxval(temp(1,1,:))
    gdd0 = sum(dtemp,mask = dtemp > 0.)

    aprec = sum(dprec)
    apet = sum(dpet)
    
    !write(0,*)tcm,gdd5(1,1),gdd0,twm,aalpha,aprec,apet

    call calcbiome(tcm,gdd5(1,1),gdd0,twm,aalpha,aprec,apet,biome(1,1))

    biomenum(biome(1,1)) = .true.

    write(30,'(3f10.3,i7)'),awc,val1,val2,biome(1,1)	
    
    write(*,*)
    
    x = x + 1

  end do
  
  98 continue

  deallocate(gdd5)
  deallocate(gridcrop)
  deallocate(gridpasture)
  deallocate(gridcarbon)
  deallocate(gridph)
  deallocate(g_annual_alpha)
  deallocate(gridalpha)
  deallocate(biome)
  
  close(30)

  write(*,*)

end do

99 continue

!----------------------------------------------------------------------

end program biome1
