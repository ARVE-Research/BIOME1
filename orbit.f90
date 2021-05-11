module orbit

use common_vars

real(dp) :: ka

real(dp) :: ecc
real(dp) :: pre
real(dp) :: perh
real(dp) :: xob

contains

subroutine orbital_parameters(ka,ecc,pre,perh,xob)

!-----------------------------------------------------------------------------
!This routine uses the orbital solutions of Berger (1978) and is valid only
!for calculations within +- 1.000.000 yr centered on 1950 AD. 
!For longer periods the Berger (1990) solution should be used.
!(Contact Berger for this 1990 solution).
!
!Recoded by J.O. Kaplan in 2002.
!
!Please refer to :
!  Berger A., 1978. A simple algorithm to compute long term
!                   variations of daily or monthly insolation.
!                   Contr. 18  Inst. of Astronomy and Geophysics,
!                   Universite Catholique de Louvain,
!                   Louvain-la-Neuve, Belgium
!
!  Berger A., 1978. Long term variations of daily insolation and
!                   Quaternary climatic changes.
!                   J. of Atmospheric Sciences 35, 2362-2367
!
!The function value returned by atan is assumed to be a real(dp)
!ranging from -pi/2 to pi/2
!
!Input parameters for the orbital solution are provided in a separate file.
!
!The read and write statements might have to be changed.
!-----------------------------------------------------------------------------

implicit none

!arguments

real(dp), intent(in) :: ka

real(dp), intent(out) :: ecc
real(dp), intent(out) :: pre
real(dp), intent(out) :: perh
real(dp), intent(out) :: xob

!parameters

integer, parameter :: nef = 19
integer, parameter :: nob = 47
integer, parameter :: nop = 78

real(dp), parameter :: pirr = pir/3600.0d0 
real(dp), parameter :: step = 360.0d0/365.25d0

!variables

integer :: i
integer :: term

integer :: neff
integer :: nobb
integer :: nopp

real(dp) :: y
real(dp) :: z
real(dp) :: xod
real(dp) :: xop
real(dp) :: prm                                         

real(dp) :: t
real(dp) :: xes
real(dp) :: xec
real(dp) :: arg
real(dp) :: tra
real(dp) :: rp
real(dp) :: prg

real(dp), dimension(19) :: ae
real(dp), dimension(19) :: be
real(dp), dimension(19) :: ce
real(dp), dimension(47) :: aob
real(dp), dimension(47) :: bob
real(dp), dimension(47) :: cob
real(dp), dimension(78) :: aop
real(dp), dimension(78) :: bop
real(dp), dimension(78) :: cop

!*******************************************                            
!   daily insolation - long term variation *                            
!*******************************************                            
!
!
!This program computes the total daily irradiation received at the top
!of the atmosphere for a given latitude and time in the ka (in kj m-2).

open(unit=8,status='old',file='INSOL.IN')

!                                                                       
!   1.earth orbital elements : eccentricity           ecc   table 1
!***************************   precessional parameter pre
!                              obliquity              xob   table 2     
!                              general precession     prg
!                              longitude perihelion   perh  table 3


!Read the header in the ASCII data file containing the parameters for the
!orbital solution

do i=1,6
  read(8,*)
end do

!read amplitude a  mean rate b  phase c                        
!they are immediately converted in radians                
!                                                                       
!nef  nob  nop  may be reduced to  19  18  9                   
!but the input data must be changed accordingly                

!-------------------------
!eccentricity                                                        

do i=1,nef
  read(8,*)term,ae(i),y,z
  be(i)=y*pirr
  ce(i)=z*pir
end do

!-------------------------
!obliquity                                                           

xod = 23.320556d0

do i=1,nob
  read(8,*)term,aob(i),y,z
  bob(i)=y*pirr
  cob(i)=z*pir
end do

!-------------------------
!general precession in longitude                                     

xop =  3.392506d0
prm = 50.439273d0

do i=1,nop
  read(8,*)term,aop(i),y,z
  bop(i)=y*pirr
  cop(i)=z*pir
end do

!-------------------------

neff=nef
nobb=nob
nopp=nop


!   3.numerical value for ecc pre xob
!************************************
!   t is negative for the past

t=ka*1000.0d0                                                      
xes=0.0d0                                                         
xec=0.0d0                                                         

do i=1,neff
  arg=be(i)*t+ce(i)
  xes=xes+ae(i)*sin(arg)
  xec=xec+ae(i)*cos(arg)
end do

ecc=sqrt(xes*xes+xec*xec)                                        
tra=abs(xec)                                                     

if(tra>1.0d-08) then
  rp=atan(xes/xec)
  if(xec>0.0d0) then !line 12
    if (xes>0.0d0) then !line 13
      perh=rp/pir
    else if (xes<0.0d0) then !line 14
      rp=rp+2.0d0*pi
      perh=rp/pir
    else !line 13
      perh=rp/pir
    end if
  elseif (xec<0.0d0) then !line 11
    rp=rp+pi
    perh=rp/pir
  else !line 10
    if (xes>0.0d0) then !line 17
      rp=pi/2.0d0
      perh=rp/pir
    else if (xes<0.0d0) then !line 15
      rp=1.5d0*pi
      perh=rp/pir
    else !line 16
      rp=0.0d0
      perh=rp/pir
    end if
  end if
else
  if (xes>0.0d0) then !line 17
    rp=pi/2.0d0
    perh=rp/pir
  else if (xes<0.0d0) then !line 15
    rp=1.5d0*pi
    perh=rp/pir
  else !line 16
    rp=0.0d0
    perh=rp/pir
  end if
end if

prg=prm*t

do i=1,nop
  arg=bop(i)*t+cop(i)
  prg=prg+aop(i)*sin(arg)
end do

prg=prg/3600.0d0+xop
perh=perh+prg

if (perh>0.0d0) then !line 53
  if(perh>360.0d0) then
    perh=perh-360.0d0
  end if
else if (perh<0.0d0) then
  perh=perh+360.0d0
end if

pre=ecc*sin(perh*pir)

xob=xod

do i=1,nobb
  arg=bob(i)*t+cob(i)
  xob=xob+aob(i)/3600.0d0*cos(arg)
end do

end subroutine orbital_parameters

end module orbit
