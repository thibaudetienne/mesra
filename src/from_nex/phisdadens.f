      subroutine phisdadens

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

C   Computes the direct Detachment/Attachment-derived phi_S index

      use declare

      implicit none

      THR = 1.0D-9

C     Reads the input files

      read(5,*) input2
      file1=trim(input2)

C     Reads the cube 1 file header

      call rdcubehdr(file1,natom,x0,y0,z0,n1,n2,n3,
     $     dx,dy,dz,dummy)

      read(5,*) input3
      file2=trim(input3)

C     Reads the cube 2 file header

      call rdcubehdr(file2,natom1,x01,y01,z01,n11,n21,n31
     $    ,dx1,dy1,dz1,dummy1)

C     Checks the consistency of the two cube files

      call consistency(natom,natom1,n1,n11,n2,n21,n3,n31,
     $ dx,dx1,dy,dy1,dz,dz1,x0,x01,
     $ y0,y01,z0,z01,THR)

C     Reads cube files grid and performs numerical integrations

      call rddensities(x0,y0,z0,file1,file2,natom,n3,n2,n1,
     $           dx,dy,dz,integ1,integ2,
     $           integp, integm, xrp, yrp, zrp,
     $           xrm, yrm, zrm, x1c, y1c, z1c,
     $           x2c, y2c, z2c, rct, r12ct, phis, 1)

C     Outputs the results

      call printingdadens(n1,integ1, integ2, integm, integp,
     $ xrm, yrm, zrm, xrp, yrp, zrp,
     $ x1c, y1c, z1c, x2c, y2c, z2c,
     $ rct, r12ct, phis)

       end
