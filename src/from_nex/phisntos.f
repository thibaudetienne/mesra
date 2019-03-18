      subroutine phisntos

C       Produced at Théorie-Modélisation-Simulation'
C       SRSMC - Université de Lorraine
C       Copyright 2014 by:
C       X. Assfeld, A. Monari, T. Very, T. Etienne

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

      use declare

      implicit none

C     Reads the input

      read(5,*) input2
      file1=trim(input2)

C     Reads the header of the cube file

      call rdntoshdr(file1,n,x0,y0,z0,n1,n2,n3,
     $ dx,dy,dz,dummy)

C    Reads the NTOs cube file grid and performs integrations

      call rdntos(file1,n,x0,y0,z0,n1,n2,n3,dx,dy,dz,integ1,
     $ integ2, integm, integp, xrm, yrm, zrm, xrp, yrp, zrp,
     $ x1c, x2c, y1c, y2c, z1c, z2c, rct, r12ct, phis)

C     Outputs the results 

      call printingntos(n1,integ1,integ2,integm,integp,
     $      xrm,yrm,zrm,xrp,yrp,zrp,
     $      x1c,x2c,y1c,y2c,z1c,z2c,rct,r12ct,phis)

       end
