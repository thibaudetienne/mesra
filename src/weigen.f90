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

subroutine weigen(eigenfile,eigenvalues,eigendim)

! Writes the eigenvalues into a file.

use declare

integer :: eigendim
character*(*) :: eigenfile
real*8 :: eigenvalues(eigendim)

! Opens the file.

open(10,file=eigenfile,form='formatted')

! Writes the eigenvalues.

do i=1,eigendim
 write(10,'(i5,f10.6)') i,eigenvalues(i)
enddo

! Closes the file.

close(10)

end
