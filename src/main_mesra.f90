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

program main_mesra

use declare

write(6,*)
write(6,*) '           ********************'
write(6,*) '           *                  *'
write(6,*) '           *  MESRA software  *'
write(6,*) '           *                  *' 
write(6,*) '           ********************'
write(6,*)
! Reads the input.

call rinp

write(6,*) 'Jobtype: ', jobtype

! If the jobtype is 'incrcube', simply calls the appropriate subroutine.

if (jobtype .eq. 'incrcube') then
 call increase_cube_size(cubefile1,x)
 write(6,*)
 stop
endif

! Opens the mesra log file.

open(50,file='mesra-Job-'//trim(jobtype)//'.log',form='formatted')

write(50,*)
write(50,*) '********************'
write(50,*) '*                  *'
write(50,*) '*  MESRA software  *'
write(50,*) '*                  *' 
write(50,*) '********************'
write(50,*)
write(50,*) 'Jobtype: ', jobtype

! Outputs data related to the calculation and starts reading files.
! These operations are restricted to some types of mesra computations.

if (jobtype .eq. 'qm_NI') then
 continue
else if (jobtype .eq. 'split') then
 continue
else if (jobtype .eq. 'cubeop') then
 continue
else if (jobtype .eq. 'alphaddagger') then
 continue
else if (jobtype .eq. 'qmnidagger') then
continue
else

write(6,*) 'Number of excited states computed, and target transition'
write(6,*) nt,ns

write(50,*) 'Number of excited states computed, and target transition'
write(50,*) nt,ns

! Reads the fchk file.

 call rfchk

! Reads the AO overlap matrix.

 call rS

! Reads the LCAO coefficients matrix.

 call rC
endif

! Depending on the job type, calls the appropriate subroutine.

write(6,*)

if (jobtype .eq. 'dau') call dau_main
if (jobtype .eq. 'daz' .or. jobtype .eq. 'dar' .or. jobtype .eq. 'adiabZ') call zvec_main

if (jobtype .eq. 'orbsXY') call orbs_XY
if (jobtype .eq. 'pNTOs') call orbs_XY
if (jobtype .eq. 'CTOs') call orbs_XY
if (jobtype .eq. 'aNTOs') call orbs_XY

if (jobtype .eq. 'daXY') call daXY

if (jobtype .eq. 'rlxy_PA') call rlxy_PA

if (jobtype .eq. 'qm_NI') call qm_NI

if (jobtype .eq. 'split') call split(cubefile)

if (jobtype .eq. 'cubeop') call cubeop(cubefile1,cubefile2,cubefile3,op1,op2)

if (jobtype .eq. 'alphaddagger') call alphaddagger

if (jobtype .eq. 'qmnidagger') call qmNIdagger

! Deallocates the matrices.

if (jobtype .eq. 'qm_NI') then
 continue
else if (jobtype .eq. 'split') then
 continue
else if (jobtype .eq. 'cubeop') then
 continue
else if (jobtype .eq. 'alphaddagger') then
continue
else if (jobtype .eq. 'qmnidagger') then
continue
else
 call deal
endif

close(50)

end
