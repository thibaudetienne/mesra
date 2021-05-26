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

subroutine wfchk_dens(matrix,nfchk,dtmp)

! Outputs the elements of a density matrix into
! an .fchk file.

use declare

integer :: dtmp
real*8 :: matrix(dtmp),dumb(dtmp)
character*(*) :: nfchk

! Opens the old .fchk file and reads/rewrites the lines
! preceding the density matrix elements.

open(72,file=fchk,status='old')
open(73,file=nfchk,form='formatted')

do
 read(72,'(a80)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:17) .eq. 'Total SCF Density' ) exit
! write(73,'(a80)') line
enddo

! Writes the density matrix elements in the .fchk file.

write(73,'(a80)') line
write(73,'(5e16.8)') (matrix(i), i=1,dtmp)

! Skips the density matrix in the original file when reading it
! to continue to copy the remaining lines until the end of the file.

read(72,'(5e16.8)') (dumb(i), i=1,dtmp)

do
 read(72,'(a80)',iostat=error) line
  if (error .ne. 0) exit
 write(73,'(a80)') line
enddo

! Closes the old and new files.

close(72)
close(73)

end
