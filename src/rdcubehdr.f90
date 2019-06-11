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

subroutine rdcubehdr(sfile,snatom,sx0,sy0,sz0,sn1,sn2,sn3,&
&                    sdx,sdy,sdz,sdummy)

! Reads the header of a cube file.

use declare

implicit none
integer :: sn1,sn2,sn3,snatom
real*8 :: sx0,sy0,sz0,sdx,sdy,sdz,sdummy
character*128 :: sfile

! Opens the cube file, and reads the characteristics of the cube.

open(10,file=sfile,form='formatted')

! snatom = number of atoms
! sx0, sy0 and sz0 = origin of the referential
! sn1, sn2, sn3 = number of steps from the origin in the x, y and z directions
! sdx, sdy, sdz = steps width

read(10,*)
read(10,*)
read(10,*) snatom,sx0,sy0,sz0
read(10,*) sn1,sdx,sdummy,sdummy
read(10,*) sn2,sdummy,sdy,sdummy
read(10,*) sn3,sdummy,sdummy,sdz

close(10)

return

end 
