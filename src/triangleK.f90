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

subroutine triangleK(matrix,matriangle,dtmp)

! Converts a square matrix into its upper triangle, in the form of a vector.

use declare

integer :: dtmp
real*8 :: matrix(dtmp,dtmp),matriangle(dtmp*(dtmp+1)/2)

! Initializes the counter for the components of the vector.

k = 0

! Transforms the square matrix elements into vector components.

do i=1,dtmp
 do j=1,i
  k = k + 1
  matriangle(k) = matrix(i,j)
 enddo
enddo

end
