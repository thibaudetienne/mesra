! MESRA software
! Molecular Electronic Structure Reorganization: Analysis
! Copyright (C) 2019 Thibaud Etienne
! More information at mesrasoftware.wordpress.com
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License v2
! as published by the Free Software Foundation.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to
! 
! Free Software Foundation, Inc. 
! 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

subroutine wfchk_orbsAlpha(singular_values,lcao_Omatrix,lcao_Vmatrix,nfchk)

use declare

real*8 :: singular_values(nea),lcao_Omatrix(nbs,nea),lcao_Vmatrix(nbs,norb-nea)
real*8,allocatable :: temporary(:)
character*(*) :: nfchk

open(72,file=fchk,status='old')
open(73,file=nfchk,form='formatted')

do
 read(72,'(A100)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:22) .eq. 'Alpha Orbital Energies' ) exit
enddo

if (norb .lt. 2*nea) then
allocate(temporary(2*nea))
else
allocate(temporary(norb))
endif

temporary = 0.0D0
k = 0

do i =1,nea
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
  if (line(1:21) .eq. 'Alpha MO coefficients') exit
enddo

temporary = 0.0d0
k = 0

do i = 1,nea
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

if (unr) then
continue
else
do
 read(72,'(a80)',iostat=error) line
  if (error .ne. 0) exit
 write(73,'(a80)') line
enddo
endif

close(72)
close(73)

end
