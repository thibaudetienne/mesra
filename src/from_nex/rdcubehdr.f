      subroutine rdcubehdr(sfile,snatom,sx0,sy0,sz0,sn1,sn2,sn3,
     $                    sdx,sdy,sdz,sdummy)

C       Produced at Théorie-Modélisation-Simulation'
C       SRSMC - Université de Lorraine'
C       Copyright 2014 by:'
C       X. Assfeld, A. Monari, T. Very, T. Etienne'


C   This file is part of NANCY_EX.

C   NANCY_EX is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.

C   NANCY_EX is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C   GNU General Public License for more details.

C   You should have received a copy of the GNU General Public License
C   along with NANCY_EX. If not, see <http://www.gnu.org/licenses/>.

C     Reads the header of a cube file

      use declare

      implicit none
      integer :: sn1,sn2,sn3,snatom
      real*8 :: sx0,sy0,sz0,sdx,sdy,sdz,sdummy
      character*128 :: sfile

      open(10,file=sfile,form='formatted')

      read(10,*)
      read(10,*)
      read(10,*) snatom,sx0,sy0,sz0
      read(10,*) sn1,sdx,sdummy,sdummy
      read(10,*) sn2,sdummy,sdy,sdummy
      read(10,*) sn3,sdummy,sdummy,sdz

      close(10)

      return
      end 
