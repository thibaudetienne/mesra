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

subroutine rfchk

! Reads the .fchk file

use declare

! Initializes the parameters and quantities

nea = 0
neb = 0
nbs = 0
norb = 0
ntr = 0

unr = .false.

! Opens the .fchk file and reads nea, neb, nbs, norb.

open(10,file=fchk,status='old')

do
read(10,'(a80)',iostat=error) line
if (error .ne. 0) exit

if (line(1:25) .eq. 'Number of alpha electrons') then
read(line,'(56x,i5)') nea
endif

if (line(1:24) .eq. 'Number of beta electrons') then
read(line,'(56x,i5)') neb
endif

if (line(1:25) .eq. 'Number of basis functions') then
read(line,'(55x,i6)') nbs
endif

if (line(1:31) .eq. 'Number of independent functions') then
read(line,'(55x,i6)') norb
endif

enddo

! Outputs some information in the log file and on the user screen.

write(50,*) 'Number of alpha and beta electrons'
write(50,*) nea, neb
write(50,*) 'Number of atomic and molecular (spin)orbitals'
write(50,*) nbs, norb
write(6,*) 'Number of alpha and beta electrons'
write(6,*) nea, neb
write(6,*) 'Number of atomic and molecular (spin)orbitals'
write(6,*) nbs, norb

! If the number of alpha and beta electrons is different, we speak in terms of unrestricted SCF calculation.

if (nea .ne. neb) unr = .true.

! Rewinds and closes the file.

close(10)

! Calculates the number of elements in the upper triangle of an nbs*nbs square matrix.

ntr = (nbs*(nbs+1))/2

end
