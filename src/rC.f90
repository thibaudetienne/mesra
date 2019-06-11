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

subroutine rC

! Reads the LCAO coefficients matrix.

use declare

! Opens the file.

open(10,file=fchk,status='old')

! Allocates the vector (for reading the elements) and the matrix (to dispatch the elements after their first reading).

allocate(lvec(norb*nbs))
lvec = 0.0d0

allocate(C(nbs,norb))
C = 0.0d0

! Reads the fchk file until the orbital energies.

do
 read(10,'(A80)', iostat=error) line
  if (error.ne.0) exit
  if (line(1:22) .eq. 'Alpha Orbital Energies') then
   marker_Ea = line
   exit
  endif
enddo

if (unr) then
 do
  read(10,'(A80)', iostat=error) line
   if (error.ne.0) exit
   if (line(1:21) .eq. 'Beta Orbital Energies') then
    marker_Eb = line
    exit
   endif
 enddo
endif

! Reads the elements in the order given by Gaussian, and stores them into a vector.

do
 read(10,'(a80)',iostat=error) line
  if (error .ne. 0) exit
   if (line(1:21) .eq. 'Alpha MO coefficients') then
    marker_Ca = line
    read(10,'(5e16.8)') (lvec(i), i=1,norb*nbs)
    exit
   endif
enddo

! Dispatches the elements into an nbs*norb matrix.

k = 0

do i = 1,norb
 do j = 1,nbs
  k = k + 1
  C(j,i) = lvec(k)
 enddo
enddo

! Computes the adjoint of C, Cdagger.

write(50,*)

if (unr) then
call orth(C,'Calpha',nbs,norb)
else
call orth(C,'C',nbs,norb)
endif

allocate(Cdagger(norb,nbs))

Cdagger = transpose(C)

! If the number of alpha and beta electrons is different, reads also the beta C matrix, Cb.

if (unr) then

! By default, Ca = C when the coefficients are read.

! Allocation of the Ca array.

allocate(Ca(nbs,norb))

Ca = C

! Allocation and initialization of Cb, re-intialization of lvec to read the Cb elements and store them
! in the form of a vector first.

allocate(Cb(nbs,norb))

Cb = 0.0d0
lvec = 0.0d0

! Reads the elements in the order given by Gaussian, and stores them into a vector.

do
 read(10,'(a80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:20) .eq. 'Beta MO coefficients') then
  marker_Cb = line
  read(10,'(5e16.8)') (lvec(i), i=1,norb*nbs)
  exit
 endif
enddo

! Dispatches the Cb elements from the lvec vector into an nbs*norb Cb matrix.

k = 0

do i = 1,norb
 do j = 1,nbs
  k = k + 1
  Cb(j,i) = lvec(k)
 enddo
enddo

! Allocates and generates the adjoint of Ca and Cb.

allocate(Cadagger(norb,nbs))
allocate(Cbdagger(norb,nbs))

Cadagger = transpose(Ca)
Cbdagger = transpose(Cb)

call orth(Cb,'Cbeta',nbs,norb)

write(50,*)

endif

! Deallocates the intermediate lvec vector, and closes the .fchk file.

deallocate(lvec)

close(10)

end
