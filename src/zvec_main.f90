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

subroutine zvec_main

! Here, countunr can take three values:

! countunr = 1 : open-shell molecules, alpha electrons
! countunr = 2 : open-shell molecules, beta electrons
! countunr = 3 : closed-shell molecules

! This routine calls the appropriate subroutines according to the jobtype, and combines
! these calls according to the open/closed-shell type of the molecule.

use declare

relax = .true.


if (unr) then

if (jobtype == 'adiabZ') then
write(6,*) '# Part A - Alpha density matrices'
write(6,*)
else
write(6,*) '# Part A - Alpha density matrices'
write(50,*) '# Part A - Alpha density matrices'
write(50,*)
endif

countunr = 1

! First calls a subroutine common to this computation and the one of the
! detachment/attachment density matrices corresponding to the relaxed
! difference density matrix.

call common_dar_daz

if (jobtype == 'daz') call daz_main
if (jobtype == 'dar') call no_adiab_z
if (jobtype == 'adiabZ') call adiab_z

if (neb .eq. 0) then
continue
else

if (jobtype == 'adiabZ') then
write(6,*)
 write(6,*) '# Part B - Beta density matrices'
write(6,*)
else
 write(6,*) '# Part B - Beta density matrices'
 write(50,*) '# Part B - Beta density matrices'
write(50,*)
endif

 countunr = 2
 
 call common_dar_daz

if (jobtype == 'daz') call daz_main
if (jobtype == 'dar') call no_adiab_z
if (jobtype == 'adiabZ') call adiab_z
 
endif

else

countunr = 3

call common_dar_daz

!!!Gabriel Breuil 11/12-04-2019

allocate(newzzd_zdz(norb,norb))

newzzd_zdz=matmul(zvec,zvec)
zcoef=0.0d0

do i=1,norb
 zcoef=zcoef+newzzd_zdz(i,i)
enddo

zcoef=0.5d0*zcoef

write(6,*)
write(50,*)
write(6,*) 'Integral of ZZ^dagger'
write(6,'(f13.5)') zcoef
write(50,*) 'Integral of ZZ^dagger'
write(50,'(f13.5)') zcoef

deallocate(newzzd_zdz)

!!! END Gabriel Breuil 11-04-2019

if (jobtype == 'daz') call daz_main
if (jobtype == 'dar') call no_adiab_z
if (jobtype == 'adiabZ') call adiab_z

endif

end
