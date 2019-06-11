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

subroutine rPxspin_g

! Reads the relaxed excited state density matrix.

use declare

! Opens the .fchk file where the excited state density matrix is stored.

open(10,file=fchk,status='old')

! Allocates and initializes the working arrays.

allocate(triangle(ntr))

! Locates the starting point where the matrix is written, and reads it and first
! stores it into a vector.

if (relax) then
 if (countrlx .eq. 2) then
  do
   read(10,'(a80)',iostat=error) line
   if (line(1:15) .eq. 'Spin CI Density') then
    read(10,'(5e16.8)') (triangle(i), i=1,ntr)
    exit
   endif
  enddo
 else if (countrlx .eq. 1) then
  do
   read(10,'(a80)',iostat=error) line
   if (line(1:22) .eq. 'Spin CI Rho(1) Density') then
    read(10,'(5e16.8)') (triangle(i), i=1,ntr)
    exit
   endif
  enddo
 endif
else
 do
  read(10,'(a80)',iostat=error) line
  if (line(1:22) .eq. 'Spin CI Rho(1) Density') then
   read(10,'(5e16.8)') (triangle(i), i=1,ntr)
   exit
  endif
 enddo
endif

! Recasts the vector components into a symmetric matrix.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pxKspin(i,j) = triangle(k)
  pxKspin(j,i) = pxKspin(i,j)
 enddo
enddo

deallocate(triangle)

close(10)

end
