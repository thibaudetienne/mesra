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

subroutine common_XY

! Transforms the transition density matrix into the MO space, extracts X and Y.

use declare

! Allocates the appropriate matrices.

if (unr) then
 if (countunr .eq. 1) then
  allocate(t(norb,norb))
  allocate(ttdagger(norb,norb))
  allocate(ttrsp(norb,norb))
 else if (countunr .eq. 2) then
  if (jobtype == 'daXY' .or. jobtype == 'rlxy_PA') then
   continue
  else
   allocate(t(norb,norb))
  endif
  allocate(ttdagger(norb,norb))
  allocate(ttrsp(norb,norb))
 endif
else
 allocate(t(norb,norb))
 allocate(ttdagger(norb,norb))
 allocate(ttrsp(norb,norb))
endif

! Shifts the values of C if countunr = 2 (beta electrons).

if (countunr .eq. 2) C = Cb
if (countunr .eq. 2) Cdagger = transpose(Cb)

call ao_to_mo(tK,t)

! If countunr = 2 (beta electrons),
! back-shifts the values of C into the alpha C matrix for further operations 

if (countunr .eq. 2) C = Ca
if (countunr .eq. 2) Cdagger = transpose(C)

! If the number of alpha and beta electrons is different, a normalization factor applies.

if (unr) t = t*dsqrt(2.0d0)/2.0d0

! Performs some tests on the transition density matrix.

call trace_mat(t,'t',norb)

ttrsp = transpose(t)

ttdagger = matmul(t,ttrsp)

call trace_mat(ttdagger,'t(t^dagger)',norb)

! Computes a transition orbitals normalization factor.

xy_norm = 0.0d0

do i=1,norb
 xy_norm = xy_norm + ttdagger(i,i)
enddo

! For the EOM formulation, outputs this normalization factor.

if (jobtype .eq. 'pNTOs' .or. jobtype .eq. 'orbsXY') then
write(6,*) '(x^dagger)x + (y^dagger)y'
write(6,'(f12.5)') xy_norm
endif

! Extracts X and Y from the transition density matrix (and their transpose).

allocate(xy_X(nel,norb-nel))
allocate(xy_Y(nel,norb-nel))
allocate(xy_Xt(norb-nel,nel))
allocate(xy_Yt(norb-nel,nel))

do i=1,nel
 do j=nel+1,norb
k = j-nel
xy_X(i,k) = t(i,j)
 enddo
enddo

xy_Xt = transpose(xy_X)

do i=nel+1,norb
 do j=1,nel
k = i-nel
xy_Yt(k,j) = t(i,j)
 enddo
enddo

xy_Y = transpose(xy_Yt)

deallocate(ttrsp,ttdagger)

end
