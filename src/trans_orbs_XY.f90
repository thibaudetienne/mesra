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

subroutine trans_orbs_XY

! Constructs the three types of transition matrices (t_1, t_2, and t_3)
! in order to extract their natural orbitals afterward.

! t_1 => pNTOs
! t_2 => CTOs
! t_3 => aNTOs

use declare

allocate(t_0p(nel,norb-nel))
allocate(t_0m(nel,norb-nel))
allocate(t_02m(nel,norb-nel))
allocate(t_02p(nel,norb-nel))
allocate(t_0pt(norb-nel,nel))
allocate(t_0pt_0pt(nel,nel))
allocate(t_1(nel,norb-nel))
allocate(t_2(nel,norb-nel))
allocate(t_3(nel,norb-nel))

! Prepares some preliminary matrices, that will be used to construct t_1, t_2, and t_3.

do i=1,nel
 do j=1,norb-nel
t_0p(i,j) = xy_X(i,j) + xy_Y(i,j)
t_0m(i,j) = xy_X(i,j) - xy_Y(i,j)
t_02m(i,j) = (xy_X(i,j))**2 - (xy_Y(i,j))**2
t_02p(i,j) = (xy_X(i,j))**2 + (xy_Y(i,j))**2
 enddo
enddo

! Computes norms that will be used for constructing t_1 and t_3.

t_0pt = transpose(t_0p)
t_0pt_0pt = matmul(t_0p,t_0pt)

t3_norm = 0.0d0

do i=1,nel
 t3_norm = t3_norm + t_0pt_0pt(i,i)
enddo

if (jobtype .eq. 'aNTOs' .or. jobtype .eq. 'orbsXY') then

write(6,*) '((x+y)^dagger)(x+y), and its square root'
write(6,'(2f12.5)') t3_norm,dsqrt(t3_norm)

endif

x2y2_norm = 0.0d0

do i=1,nel
 do j=1,norb-nel
x2y2_norm = x2y2_norm + t_02m(i,j)
 enddo
enddo

write(6,*) '(x^dagger)x - (y^dagger)y'
write(6,'(f12.5)') x2y2_norm

xy_residue = 0.0d0

do i=1,nel
 do j=1,norb-nel
xy_residue = xy_residue + dabs(t_02m(i,j)) - t_02m(i,j)
 enddo
enddo

xy_residue = x2y2_norm + xy_residue

if (jobtype .eq. 'CTOs' .or. jobtype .eq. 'orbsXY') then

write(6,*) 'CTOs residue'
write(6,'(f12.5)') xy_residue - x2y2_norm

endif

! Constructs t_1, t_2, and t_3.

do i=1,nel
 do j=1,norb-nel
t_1(i,j) = ((dabs(t_0p(i,j))/t_0p(i,j))*dsqrt(t_02p(i,j)))

t_2(i,j) = ((dabs(t_0m(i,j))/t_0m(i,j)))*&
((dabs(t_02m(i,j)))/t_02m(i,j))*&
dsqrt(dabs(t_02m(i,j)))/(dsqrt(xy_residue))

t_3(i,j) = t_0p(i,j)/dsqrt(t3_norm)
 enddo
enddo

deallocate(xy_X,xy_Y,xy_Xt,xy_Yt)
deallocate(t_0p,t_0m,t_02m,t_02p,t_0pt,t_0pt_0pt)

end
