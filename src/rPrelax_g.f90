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

subroutine rPrelax_g

! Reads the relaxed excited state density matrix.

use declare

! Opens the .fchk file where the excited state density matrix is stored.

open(10,file=fchk,status='old')

! Allocates and initializes the working arrays.

allocate(pxKrelaxed(nbs,nbs))
allocate(triangle(ntr))

! Locates the starting point where the matrix is written, and reads it and first
! stores it into a vector.

do
 read(10,'(a80)',iostat=error) line
  if (line(1:16) .eq. 'Total CI Density') then
   read(10,'(5e16.8)') (triangle(i), i=1,ntr)
   exit
  endif
enddo

! Recasts the vector components into a symmetric matrix.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pxKrelaxed(i,j) = triangle(k)
  pxKrelaxed(j,i) = pxKrelaxed(i,j)
 enddo
enddo

! Contracts the relaxed AO-space excited state density matrix with the 
! basis functions overlap matrix.

allocate(pxKrelaxedS(nbs,nbs))

pxKrelaxedS = matmul(pxKrelaxed,S)

! If the adiabatic connection of the Z-vector has not been required,
! computes the trace of the contracted matrix, deallocates the working 
! vector, and closes the .fchk file.

if (adiab) then
 continue
else
 call trace_mat(pxKrelaxedS,'pxKrelaxedS',nbs)
endif

deallocate(triangle)

close(10)

end
