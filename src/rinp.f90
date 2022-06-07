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

subroutine rinp

! Reads the input.

use declare

! Finds the name of the input file, opens it.

call getarg(1,input)
input=trim(input)

open(10,file=input,status='old')

! Reads the job type. According to the job type, calls the appropriate subroutines to find complementary information.

read(10,*) jobtype

! "da" = detachment/attachment

if (jobtype(1:2) .eq. 'da') then

!! Antoine Marion: 2020-12-31
!!  - add option to perform DA with two densities comming from two different calculations
if (jobtype .eq. 'dau2') then
call rgen_info_2
call rPA_status

else ! regular
call rgen_info
call rPA_status

endif

! "adiabz" = adiabatic connection of the Z-vector to the difference density matrix.

else if (jobtype .eq. 'adiabZ') then

adiab = .true.
call rgen_info
call rPA_status

! "orbsXY" = transition orbitals, computed from a XY-type method.
! Those methods include: CIS, TDA and RPA, TDHF, TDDFT, BSE. They have the particularity to have a
! block-diagonal difference density matrix.

else if (jobtype .eq. 'orbsXY' .or. jobtype .eq. 'pNTOs' .or. jobtype .eq. 'CTOs' .or. &
jobtype .eq. 'aNTOs') then

call rgen_info
call rPA_status

! "daXY" = detachment/attachment, computed from a XY-type method.
! Those methods include: CIS, TDA and RPA, TDHF, TDDFT, BSE. They have the particularity to have a
! block-diagonal difference density matrix.

else if (jobtype .eq. 'daXY') then

call rgen_info
call rPA_status

! "rlxy_PA" = relaxation-adapted descriptors, computed with a population analysis only,
! and from a XY-type method.
! Those methods include: CIS, TDA and RPA, TDHF, TDDFT, BSE. They have the particularity to have a
! block-diagonal difference density matrix.

else if (jobtype .eq. 'rlxy_PA') then

call rgen_info
call rPA_status

! "qm_NI" Computes the density-based descriptors by Numerical Integration (NI). The grids are read from cube files.

else if (jobtype .eq. 'qm_NI') then

call rqm_NI_info

! "split" Separates grid points according to their sign and integrates.

else if (jobtype .eq. 'split') then

read(10,*) cubefile

! "cubeop" Performs operations on cubes: cubefile3 = op1*cubefile1 + op2*cubefile2

else if (jobtype .eq. 'cubeop') then

read(10,*) cubefile1
read(10,*) cubefile2
read(10,*) cubefile3
read(10,*) op1
read(10,*) op2

! "alphaddagger" Performs the population analysis, based on matrices already computed in a previous operation.
! NB: rdaggerf reads the name of the matrices .fchk files.

else if (jobtype .eq. 'alphaddagger') then

call rgen_info
call rPA_status
call rdaggerf

! "qmnirlx" Computes the relaxed density-based descriptors, based on calculations already performed previously.

else if (jobtype .eq. 'qmnirlx') then

read(10,*) theta_unrelaxed
read(10,*) phiS_unrelaxed
read(10,*) phi_unrelaxed
read(10,*) eta

! "incrcube" outputs the header to provide to the cubegen utility for getting cubes of increased size.

else if (jobtype .eq. 'incrcube') then

read(10,*) cubefile1
read(10,*) x

endif

! Closes the input file.

close(10)
end
