      subroutine consistency(snatom,snatom1,sn1,sn11,sn2,sn21,sn3,sn31,
     $ sdx,sdx1,sdy,sdy1,sdz,sdz1,sx0,sx01,
     $ sy0,sy01,sz0,sz01,sTHR)

C       Produced at Théorie-Modélisation-Simulation'
C       SRSMC - Université de Lorraine'
C       Copyright 2014 by:'
C       X. Assfeld, A. Monari, T. Very, T. Etienne'


C   This file is part of NANCY_EX.
C
C   NANCY_EX is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.
C
C   NANCY_EX is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C   GNU General Public License for more details.

C   You should have received a copy of the GNU General Public License
C   along with NANCY_EX. If not, see <http://www.gnu.org/licenses/>.

      use declare

      implicit none
      integer :: snatom,snatom1,sn1,sn2,sn3,sn11,sn21,sn31
      real*8 :: sdx,sdx1,sdy,sdy1,sdz,sdz1,sx0,sx01,sy0,sy01,
     $ sz0,sz01,sTHR

C      Checks the consistency of the cube files

      if (snatom .ne. snatom1) then
         write(6,*) 'Different Number of atoms between cubes', 
     $               snatom, snatom1
         stop
      endif

      if (sn1 .ne. sn11) then
         write(6,*) 'Different x size of cube',sn1,sn11 
         stop
      endif

      if (sn2 .ne. sn21) then
         write(6,*) 'Different y size of cube',sn2,sn21 
         stop
      endif

      if (sn3 .ne. sn31) then
         write(6,*) 'Different z size of cube',sn3,sn31 
         stop
      endif

      if (dabs(sdx-sdx1) .gt. sTHR .or. dabs(sdy-sdy1) .gt. sTHR .or. 
     $     dabs(sdz-sdz1) .gt. sTHR) then
          write(6,*) 'Different increments'
          stop
      endif 
      if (dabs(sx0-sx01) .gt. sTHR .or. dabs(sy0-sy01) .gt. sTHR .or. 
     $     dabs(sz0-sz01) .gt. sTHR) then
          write(6,*) 'Different origins for the Cubes'
          stop

      endif 
      end
