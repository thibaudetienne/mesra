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

subroutine daz_main

! Computes the detachment/attachment density matrices
! related to the Z-vector relaxation process.

use declare

! Computes the trace of the Z-vector matrix in the spin-orbital space.

if (countunr .eq. 1) call trace_mat(zvec,'zvecAlpha',norb)
if (countunr .eq. 2) call trace_mat(zvec,'zvecBeta',norb)
if (countunr .eq. 3) call trace_mat(zvec,'zvec',norb)

write(50,*)

!! Computes the Z(Z^dagger) \oplus (Z^dagger)Z matrix (Hereafter named zzd_zdz).
!
!zzd_zdz = matmul(zvec,zvec)
!
!zzd_zdz = 0.5d0*zzd_zdz
!
!! Computes the trace of Z(Z^dagger).
!
!if (countunr .eq. 1) call trace_mat(zzd_zdz,'Z(Z^dagger) Alpha',norb)
!if (countunr .eq. 2) call trace_mat(zzd_zdz,'Z(Z^dagger) Beta',norb)
!if (countunr .eq. 3) call trace_mat(zzd_zdz,'Z(Z^dagger)',norb)

! Allocates the eigenvectors matrix, and the eigenvalues array.

allocate(Uvec(norb,norb))
allocate(lvec(norb))

! Performs the detachment/attachment procedure to the Z-vector matrix.

if (countunr .eq. 1) call det_at(zvec,'Zalpha',norb,Uvec,lvec)
if (countunr .eq. 2) call det_at(zvec,'Zbeta',norb,Uvec,lvec)
if (countunr .eq. 3) call det_at(zvec,'Z',norb,Uvec,lvec)

! Deallocates the allocated arrays.

deallocate(zvec,Uvec,lvec)

end
