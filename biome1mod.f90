module biome1mod

implicit none

public :: calcbiome

integer, parameter :: npfts = 13

type pftparameters
  character(80) :: name
  real :: tcmax
  real :: tcmin
  real :: gdd5min
  real :: gdd0min
  real :: twmin
  real :: almin
  real :: almax
end type pftparameters

type(pftparameters), dimension(npfts) :: pftpar

namelist /pftpars/ pftpar

contains

!---------------------------------------------------

subroutine calcbiome(tcm,gdd5,gdd0,twm,alpha,aprec,apet,biome)

implicit none

!arguments

real, intent(in) :: tcm
real, intent(in) :: gdd5
real, intent(in) :: gdd0
real, intent(in) :: twm
real, intent(in) :: alpha  !annual mean
real, intent(in) :: aprec  !annual sum
real, intent(in) :: apet  !annual sum

integer(2), intent(out) :: biome

!parameters

real, parameter :: ud = -9999.

integer, parameter :: tests = 7

!local variables

integer :: pft

integer :: bit
integer :: bval

logical, dimension(npfts) :: present
logical, dimension(tests) :: passed

!---------

biome = -9999

!calculate present pfts

present = .true.

do pft = 1,npfts

  passed = .true.

  !tcmin
  if (pftpar(pft)%tcmin /= ud .and. tcm < pftpar(pft)%tcmin)      passed(1) = .false.

  !tcmax
  if (pftpar(pft)%tcmax /= ud .and. tcm > pftpar(pft)%tcmax)      passed(2) = .false.

  !GDD5
  if (pftpar(pft)%gdd5min /= ud .and. gdd5 < pftpar(pft)%gdd5min) passed(3) = .false.

  !GDD0
  if (pftpar(pft)%gdd0min /= ud .and. gdd0 < pftpar(pft)%gdd0min) passed(4) = .false.

  !twmin
  if (pftpar(pft)%twmin /= ud .and. twm < pftpar(pft)%twmin)      passed(5) = .false.

  !alpha min
  if (pftpar(pft)%almin /= ud .and. alpha < pftpar(pft)%almin)    passed(6) = .false.

  !alpha max
  if (pftpar(pft)%almax /= ud .and. alpha > pftpar(pft)%almax)    passed(7) = .false.
  
  if (all(passed)) then 
    present(pft) = .true.
  else
    present(pft) = .false.
  end if

end do

if (aprec - apet < -100.) then
  present(1:7) = .false.
end if

!choose biome based on present pfts

!calculate an unique integer value that represents the PFTs present

bval = 0

do pft = 1,npfts

  if (present(pft)) then
    bit = 1
  else
    bit = 0
  end if
  
  bval = bval + bit * 2**(npfts-pft)
  
end do

do

  if (present(1)) then
    if (present(2)) then
      biome = 2  !tropical raingreen
      exit
    else
      biome = 1  !tropical evergreen
      exit
    end if
  end if

  if (present(2)) then
    biome = 3   !tropical deciduous
    exit
  end if

  if (present(3)) then
    if (present(4)) then
      biome = 5   !temperate deciduous because temperate deciduous will always dominate in situations where it is present
      exit
    else
      biome = 4   !temperate evergreen (warm mixed)
      exit
    end if
  end if

  if (present(4)) then
    if (present(5) .and. present(6) .and. present(7)) then
      biome = 6  !cool mixed
      exit
    !else if (present(5) .and. present(7)) then
    else if (present(5) .and. tcm < 1.) then
      biome = 6
      exit
    else
      biome = 5  !temperate deciduous
      exit
    end if
  end if

  if (present(5) .and. present(6) .and. present (7)) then
    biome = 7   !cool conifer
    exit
  end if

  if (present(6) .and. present(7)) then
    biome = 8  !cold evergreen
    exit
  end if

  if (present(5) .and. present(7)) then
    biome = 9 !cold mixed
    exit
  end if

  if (present(7)) then
    biome = 10  !cold deciduous
    exit
  end if

  if (present(8)) then
    biome = 11  !xerophytic (sclerophyll)
    exit
  end if

  if (present(9)) then
    biome = 12  !warm grass
    exit
  end if

  if (present(10) .and. present(11)) then
    biome = 13  !cool grass
    exit
  else if (present(11)) then
    biome = 14  !tundra
    exit
  end if

  if (present(12)) then
    biome = 15  !hot desert
    exit
  end if

  if (present(13)) then
    biome = 16  !cool desert
    exit
  end if

  biome = 17  !barren

end do

!---

biome = bval

!write(0,*)tcm,gdd5,gdd0,twm,alpha
!write(0,*)present
!write(0,*)biome

end subroutine calcbiome

end module biome1mod
