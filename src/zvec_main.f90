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

! Computes alpha zcoef (eta)

allocate(newzzd_zdz(norb,norb))

newzzd_zdz=matmul(zvec,zvec)
zcoefalpha=0.0d0

do i=1,norb
 zcoefalpha=zcoefalpha+newzzd_zdz(i,i)
enddo

zcoefalpha=0.5d0*zcoefalpha

write(6,*) 'Trace of alpha ZZ^dagger (eta)'
write(6,'(f13.5)') zcoefalpha
write(50,*) 'Trace of alpha ZZ^dagger (eta)'
write(50,'(f13.5)') zcoefalpha

if (jobtype == 'daz') write(6,*)
if (jobtype == 'daz') write(50,*)

deallocate(newzzd_zdz)

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

! Computes beta zcoef (eta)

allocate(newzzd_zdz(norb,norb))

newzzd_zdz=matmul(zvec,zvec)
zcoefbeta=0.0d0

do i=1,norb
 zcoefbeta=zcoefbeta+newzzd_zdz(i,i)
enddo

zcoefbeta=0.5d0*zcoefbeta

write(6,*) 'Trace of beta ZZ^dagger (eta)'
write(6,'(f13.5)') zcoefbeta
write(50,*) 'Trace of beta ZZ^dagger (eta)'
write(50,'(f13.5)') zcoefbeta

deallocate(newzzd_zdz)

if (jobtype == 'daz') call daz_main
if (jobtype == 'dar') call no_adiab_z
if (jobtype == 'adiabZ') call adiab_z
 
endif

else

countunr = 3

call common_dar_daz

! Implemented by GB

allocate(newzzd_zdz(norb,norb))

newzzd_zdz=matmul(zvec,zvec)
zcoef=0.0d0

do i=1,norb
 zcoef=zcoef+newzzd_zdz(i,i)
enddo

zcoef=0.5d0*zcoef

write(6,*) 'Trace of ZZ^dagger (eta)'
write(6,'(f13.5)') zcoef
write(50,*) 'Trace of ZZ^dagger (eta)'
write(50,'(f13.5)') zcoef

deallocate(newzzd_zdz)

! END of GB implementation

if (jobtype == 'daz') call daz_main
if (jobtype == 'dar') call no_adiab_z
if (jobtype == 'adiabZ') call adiab_z

endif

end
