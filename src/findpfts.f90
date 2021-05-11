program findpfts

integer, parameter :: npfts = 13

integer :: inval
integer :: val

integer, dimension(npfts) :: binary
logical, dimension(npfts) :: present

integer :: pft,i

character(10) :: cval

call getarg(1,cval)

read(cval,*)inval

!inval = 278

val = inval

do i = 1,npfts
  
  pft = i !npfts - i + 1
  if (mod(val,2) == 0) then
    present(pft) = .false.
  else
    present(pft) = .true.
  end if
  
  !write(*,*)pft,inval,present(pft)

  val = val / 2
  
end do

write(*,*)inval,present

end program findpfts
