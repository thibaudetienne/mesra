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

subroutine consistency(snatom,snatom1,sn1,sn11,sn2,sn21,sn3,sn31, &
& sdx,sdx1,sdy,sdy1,sdz,sdz1,sx0,sx01, &
& sy0,sy01,sz0,sz01,sTHR)

! Checks the consistency of the cube files.

use declare

integer :: snatom,snatom1,sn1,sn2,sn3,sn11,sn21,sn31
real*8 :: sdx,sdx1,sdy,sdy1,sdz,sdz1,sx0,sx01,sy0,sy01, &
& sz0,sz01,sTHR

! The characteristics extracted by the rdcubehdr subroutine are compared
! for the two cubes, and if they are not strictly equivalent, the user is alerted.

if (snatom .ne. snatom1) then
   write(6,*) 'Different number of atoms between the cubes, exit!', &
& snatom, snatom1
   stop
endif

if (sn1 .ne. sn11) then
   write(6,*) 'Different x-length between the cubes, exit!',sn1,sn11 
   stop
endif

if (sn2 .ne. sn21) then
   write(6,*) 'Different y-length between the cubes, exit!',sn2,sn21 
   stop
endif

if (sn3 .ne. sn31) then
   write(6,*) 'Different z-length between the cubes, exit!',sn3,sn31 
   stop
endif

if (dabs(sdx-sdx1) .gt. sTHR .or. dabs(sdy-sdy1) .gt. sTHR .or. &
& dabs(sdz-sdz1) .gt. sTHR) then
    write(6,*) 'Different distance increments between the cubes, exit!'
    stop
endif 

if (dabs(sx0-sx01) .gt. sTHR .or. dabs(sy0-sy01) .gt. sTHR .or. &
& dabs(sz0-sz01) .gt. sTHR) then
    write(6,*) 'The two cubes have a different origin. Exit!'
    stop
endif 
end
