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

subroutine trace_mat(mat,whichmat,matsize)

! Computes and outputs the trace of a matrix.

use declare

integer :: matsize
real*8 :: mat(matsize,matsize)
character*(*) :: whichmat

! Initializes the counter x and computes the trace of the matrix mat.

x = 0.0d0

do i=1,matsize
 x = x + mat(i,i)
enddo

! Outputs the result, and reinitializes the counter x.

write(50,*) 'Trace of ',whichmat

if (dabs(x) .lt. 0.00001d0) then
 write(50,'(f12.5)') 0.00000
else
 write(50,'(f12.5)') x
endif

x = 0.0d0

end
