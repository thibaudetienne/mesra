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

subroutine daXY

! Creates the detachment/attachment density matrices for XY-type excited-state calculations.

use declare

! Allocates and reads the AO-space transition density matrix.

allocate(tK(nbs,nbs))
call rTK

write(50,*)

if (unr) then

write(6,*) '# Part A - Alpha density matrices'
write(50,*) '# Part A - Alpha density matrices'
write(50,*)

if (nea .eq. 0) then
continue
else

! nel and countunr are key numbers allowing us to use a single subroutine for both alpha, beta
! electrons of open-shell molecules, as well as for closed-shell molecules.

countunr = 1
nel = nea

! Calls the actual detachment/attachment subroutine.

call da_XY

endif

! For open-shell molecules, repeats the operations for the beta electrons.

if (neb .eq. 0) then
continue
else

countunr = 2

write(6,*) '# Part B - Beta density matrices'
write(50,*) '# Part B - Beta density matrices'
write(50,*)

nel = neb
tK = tKb

call da_XY

endif

else

! For closed-shell molecules, simply launches the detachment/attachment subroutine.

nel = nea
call da_XY

endif

end
