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

subroutine no_adiab_z

! Actual computation of the detachment/attachment density matrices from the fully relaxed
! difference density matrix.

use declare

allocate(gD(norb,norb))
allocate(Uvec(norb,norb))
allocate(lvec(norb))

! Constructs the relaxed difference density matrix.

if (countunr .eq. 1) gD = pxalpharelaxed - palpha
if (countunr .eq. 2) gD = pxbetarelaxed - pbeta
if (countunr .eq. 3) gD = pxrelaxed - p

! Constructs the detachment/attachment density matrices from the relaxed difference density matrix.

if (countunr .eq. 1) call det_at(gD,'Ralpha',norb,Uvec,lvec)
if (countunr .eq. 2) call det_at(gD,'Rbeta',norb,Uvec,lvec)
if (countunr .eq. 3) call det_at(gD,'R',norb,Uvec,lvec)

! Prints the relaxed difference density matrix in an .fchk file.

if (countunr .eq. 1) call print_mat_mo_to_ao_fchk(gD,'DeltaRalpha.fchk')
if (countunr .eq. 2) call print_mat_mo_to_ao_fchk(gD,'DeltaRbeta.fchk')
if (countunr .eq. 3) call print_mat_mo_to_ao_fchk(gD,'DeltaR.fchk')

! Deallocates the working arrays.

deallocate(gD,Uvec,lvec)

end
