!figure out the length of the time dimension

status = nf90_inq_dimid(ncid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)


write(0,*)

write(0,*)
write(0,*)'number of unique biomes: ',count(biomenum)

stop

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

    where (biome == i) biome = j
    j = j + 1

  end if
end do
write(0,*)

!end of grid loop

!----------------------------------------------------------------------
!write output

