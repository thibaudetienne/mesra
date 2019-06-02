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
