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

subroutine common_dar_daz

! Subroutine common to the computation of the detachment/attachment density matrices from the 
! Z-vector or from the fully relaxed difference density matrices.

use declare
! If the unrelaxed quantum-topological metrics were already computed (if the jobtype is rlxy_LA), 
! then it is not necessary to read again the ground and unrelaxed excited states density matrices.

if (countunr .eq. 1 .or. countunr .eq. 3) call rP

! If the adiabatic connection of the Z-vector was not required, computes the trace
! of the relaxed excited state density matrix, for diagnosis purpose.

if (adiab) then
continue
else
if (countunr .eq. 1) call trace_mat(pxalphaRelaxed,'pxalphaRelaxed',norb)
if (countunr .eq. 2) call trace_mat(pxbetaRelaxed,'pxbetaRelaxed',norb)
if (countunr .eq. 3) call trace_mat(pxrelaxed,'pxrelaxed',norb)
endif

! Writes the Z-vector.

if (countunr .eq. 1 .or. countunr .eq. 3) allocate(zvec(norb,norb))

if (countunr .eq. 2) then
 if (jobtype == 'daz') allocate(zvec(norb,norb))
endif

if (countunr .eq. 1) zvec = pxalphaRelaxed - pxalpha
if (countunr .eq. 2) zvec = pxbetaRelaxed - pxbeta
if (countunr .eq. 3) zvec = pxrelaxed - px

! If the adiabatic connection of the Z-vector was not required, prints
! the Z-vector in an .fchk file, and allocates a zzd_zdz matrix for further operations.

if (adiab) then
 continue
else
if (countunr .eq. 1) call print_mat_mo_to_ao_fchk(zvec,'DeltaZalpha.fchk')
if (countunr .eq. 2) then
 call print_mat_mo_to_ao_fchk(zvec,'DeltaZbeta.fchk')
 zvec = zvec + pxalphaRelaxed - pxalpha
 call print_mat_mo_to_ao_fchk(zvec,'DeltaZ.fchk')
 zvec = zvec - pxalphaRelaxed + pxalpha
endif
if (countunr .eq. 3) call print_mat_mo_to_ao_fchk(zvec,'DeltaZ.fchk')
endif

if (countunr .eq. 1 .or. countunr .eq. 3) allocate(zzd_zdz(norb,norb))

if (countunr .eq. 2) then
 if (jobtype == 'daz') allocate(zzd_zdz(norb,norb))
endif

end
