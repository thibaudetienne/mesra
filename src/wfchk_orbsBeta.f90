subroutine wfchk_orbsBeta(singular_values,lcao_Omatrix,lcao_Vmatrix,nfchk)

use declare

real*8 :: singular_values(norb),lcao_Omatrix(nbs,norb),lcao_Vmatrix(nbs,norb)
real*8,allocatable :: temporary(:)
character*(*) :: nfchk

open(72,file=fchk,status='old')
open(73,file=nfchk,form='formatted')

do
 read(72,'(A100)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:21) .eq. 'Beta Orbital Energies' ) exit
enddo

if (norb .lt. 2*neb) then
allocate(temporary(2*neb))
else
allocate(temporary(norb))
endif

temporary = 0.0D0
k = 0

do i =1,neb
 k = k + 1
 temporary(k) = singular_values(i)
 k = k + 1
 temporary(k) = singular_values(i)
enddo

write(73,'(a80)') line
write(73,'(5e16.8)') (temporary(i), i=1,norb)

deallocate(temporary)
allocate(temporary(2*norb*nbs))

do
 read(72,'(A80)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:20) .eq. 'Beta MO coefficients') exit
enddo

temporary = 0.0d0
k = 0

do i = 1,norb
 do j = 1, nbs 
  k = k + 1
  temporary(k) = lcao_Omatrix(j,i)
 enddo
 do j = 1, nbs
   k = k + 1
   temporary(k) = lcao_Vmatrix(j,i)
 enddo
enddo

write(73,'(a80)') line
write(73,'(5e16.8)') (temporary(i), i = 1,norb*nbs)

read(72,'(5e16.8)') (temporary(i), i = 1,norb*nbs)

deallocate(temporary)

!write(73,'(5e16.8)') ((lcao_Omatrix(j,i), j=1,nbs),(lcao_Vmatrix(j,i), j=1,nbs),i=1,nea)

do
 read(72,'(a80)',iostat=error) line
  if (error .ne. 0) exit
 write(73,'(a80)') line
enddo

close(72)
close(73)

end
