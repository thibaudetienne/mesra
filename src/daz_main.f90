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

! Computes the Z(Z^dagger) \oplus (Z^dagger) \oplus Z matrix (Hereafter named zzd_zdz).

zzd_zdz = matmul(zvec,zvec)

! Computes the trace of zzd_zdz.

if (countunr .eq. 1) call trace_mat(zzd_zdz,'zzd_zdzAlpha',norb)
if (countunr .eq. 2) call trace_mat(zzd_zdz,'zzd_zdzBeta',norb)
if (countunr .eq. 3) call trace_mat(zzd_zdz,'zzd_zdz',norb)

! Allocates the eigenvectors matrix, and the eigenvalues array.

allocate(Uvec(norb,norb))
allocate(lvec(norb))

! Performs the detachment/attachment procedure to the Z-vector matrix.

if (countunr .eq. 1) call det_at(zvec,'Zalpha',norb,Uvec,lvec)
if (countunr .eq. 2) call det_at(zvec,'Zbeta',norb,Uvec,lvec)
if (countunr .eq. 3) call det_at(zvec,'Z',norb,Uvec,lvec)

! In case the population analysis-based derivation of the relaxed density-based descriptors
! was required, computes the theta_Z factor.

if (subjobtype .eq. 'rlxy_PA') then
 x = 0.0d0
 y = 0.0d0
 do i=1,norb
  if (lvec(i) .lt. 0.0d0) x = x - lvec(i)
  if (lvec(i) .gt. 0.0d0) y = y + lvec(i)
 enddo
 thetaZ = 0.5d0*(x+y)

write(6,*)
write(50,*)
write(6,*) 'Integral of relaxation-related detachment/attachment density'
write(6,'(f13.4)') thetaZ
write(50,*) 'Integral of relaxation-related detachment/attachment density'
write(50,'(f13.5)') thetaZ

if (countunr .eq. 1)  write(6,*) 
if (countunr .eq. 1)  write(50,*) 
if (countunr .eq. 1) thetaZalpha = thetaZ
if (countunr .eq. 2) thetaZbeta = thetaZ
endif

! Deallocates the allocated arrays.

deallocate(zvec,zzd_zdz,Uvec,lvec)

end
