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

subroutine rP_g

! Reads the ground state density matrix from a Gaussian file.

use declare

! Opens the fchk file, allocates and initializes the triangle vector.
! (triangle stands for upper triangle of the symmetric matrix to read)

open(10,file=fchk,status='old')

allocate(triangle(ntr))

! Reads the AO-space ground state density matrix and stores its components
! into the triangle vector.

do
 read(10,'(a80)',iostat=error) line
  if (line(1:17) .eq. 'Total SCF Density') then
   read(10,'(5e16.8)') (triangle(i), i=1,ntr)
   exit
  endif
enddo

! Brings the triangle vector components into matrix elements.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pK(i,j) = triangle(k)
  pK(j,i) = pK(i,j)
 enddo
enddo

if (unr) then

allocate(trianglespin(ntr))
allocate(pKspin(nbs,nbs))

do
 read(10,'(a80)',iostat=error) line
  if (line(1:16) .eq. 'Spin SCF Density') then
   read(10,'(5e16.8)') (trianglespin(i), i=1,ntr)
   exit
  endif
enddo

! Brings the triangle vector components into matrix elements.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pKspin(i,j) = trianglespin(k)
  pKspin(j,i) = pKspin(i,j)
 enddo
enddo

pKalpha = (pK + pKspin)/2.0d0
pKbeta = (pK - pKspin)/2.0d0

endif

! Deallocates the triangle vectors and closes the fchk file.

deallocate(triangle)

if (unr) deallocate(trianglespin,pKspin)

close(10)

end

!! AM: 2020-12-31
subroutine rP_g_2

! Reads the ground state density matrix from a Gaussian file.
! For the second state.

use declare

! Opens the fchk file, allocates and initializes the triangle vector.
! (triangle stands for upper triangle of the symmetric matrix to read)

open(10,file=fchk2,status='old')

allocate(triangle(ntr))

! Reads the AO-space ground state density matrix and stores its components
! into the triangle vector.

do
 read(10,'(a80)',iostat=error) line
  if (line(1:17) .eq. 'Total SCF Density') then
   read(10,'(5e16.8)') (triangle(i), i=1,ntr)
   exit
  endif
enddo

! Brings the triangle vector components into matrix elements.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pxK(i,j) = triangle(k)
  pxK(j,i) = pxK(i,j)
 enddo
enddo

if (unr) then

allocate(trianglespin(ntr))
allocate(pKspin(nbs,nbs))

do
 read(10,'(a80)',iostat=error) line
  if (line(1:16) .eq. 'Spin SCF Density') then
   read(10,'(5e16.8)') (trianglespin(i), i=1,ntr)
   exit
  endif
enddo

! Brings the triangle vector components into matrix elements.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pKspin(i,j) = trianglespin(k)
  pKspin(j,i) = pKspin(i,j)
 enddo
enddo

!! FIXME AM should we declare pKxalpha and pKxbeta ?
pKalpha = (pxK + pKspin)/2.0d0
pKbeta = (pxK - pKspin)/2.0d0

endif

! Deallocates the triangle vectors and closes the fchk file.

deallocate(triangle)

if (unr) deallocate(trianglespin,pKspin)

close(10)

end
