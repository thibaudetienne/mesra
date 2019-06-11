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

subroutine rS

! Reads the basis functions overlap matrix S.

use declare

! Opens the "overlap" file.

open(10,file=fov,status='old')

! Allocates and initializes the "triangle" array for reading the S elements before dispatching them into an nbs*nbs matrix.

allocate(triangle(ntr))
triangle = 0.0d0

! Reads the S elements and places them into a vector.

do
 read(10,'(a80)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:19) .eq. ' Dump of file   514') then
   read(10,'(1x,5e20.8)') (triangle(i),i=1,ntr)
   exit
  endif
enddo

! Allocates and initializes S.

allocate(S(nbs,nbs))
S = 0.0d0

! Dispatches the "triangle" vector components into the square S matrix.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  S(i,j) = triangle(k)
  S(j,i) = S(i,j)
 enddo
enddo

! Deallocates triangle and closes the overlap file.

deallocate(triangle)

close(10)

end
