      program highcube

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

      implicit none
      integer :: n1,n2,n3,natom
      real*8 :: x0,y0,z0,dx,dy,dz,dummy,s, s1,s2,s3
      character*128 :: filename1,input1
      character*128 :: filename2,input2
      real*8 :: l1,l2,l3,l12,l22,l32
      integer :: n12, n22, n32
      integer, parameter :: out_unit=11, out_unit2=12
      real*8 :: x02,y02,z02

c     reads the input parameters; opens and reads the input files. For
c     more details, see the Cubegen utility manual

      call getarg(1,input1)
      filename1=trim(input1)

      open(10,file=filename1,form='formatted')

      read(10,*)
      read(10,*)
      read(10,*) natom,x0,y0,z0
      read(10,*) n1,dx,dummy,dummy
      read(10,*) n2,dummy,dy,dummy
      read(10,*) n3,dummy,dummy,dz

      close(10)

      dummy=0.000000

c     opens and reads the input file in which the percentage of
c     increase of each coordinate is written

      call getarg(2,input2)
      filename2=trim(input2)

      open(10,file=filename2,form='formatted')

      read(10,*) s

      close(10)

      s=ceiling(s)
      n12=ceiling(n1+(s*n1/100))

c     outputs the percentage of increase, original and final number of
c     points per coordinate

      write(6,'(a35,f5.1,2i5)') "s, n1, n12",s,n1,n12
      n22=ceiling(n2+(s*n2/100))
      write(6,'(a35,f5.1,2i5)') "s, n2, n22",s,n2,n22
      n32=ceiling(n3+(s*n3/100))
      write(6,'(a35,f5.1,2i5)') "s, n3, n32",s,n3,n32

c    setting the new origin of the enlarged cube

      x02 = x0 - ceiling(0.5*s*n1/100)*dx
      y02 = y0 - ceiling(0.5*s*n2/100)*dy
      z02 = z0 - ceiling(0.5*s*n3/100)*dz

      write(6,'(a35,3f14.6)') "x0, y0, z0",x0,y0,z0
      write(6,'(a35,3f14.6)') "x02,y02,z02",x02,y02,z02

c    computes the original and final lengths (for coordinate 1, 2 and 3)


      l1 = x0 + (n1-1)*dx
      l2 = y0 + (n2-1)*dy
      l3 = z0 + (n3-1)*dz 

      l12 = x02 + (n12-1)*dx
      l22 = y02 + (n22-1)*dy
      l32 = z02 + (n32-1)*dz

c     computes the length increase

      s1 = (l12 - x02) - (l1 - x0)
      s2 = (l22 - y02) - (l2 - y0)
      s3 = (l32 - z02) - (l3 - z0)

c     outputs informations about the new header construction

      write(6,'(a35,3f14.6)') "l1,l2,l3",l1,l2,l3
      write(6,'(a35,3f14.6)') "l12,l22,l32",l12,l22,l32
      write(6,'(a35,4f14.6)') "X in/out size + gain and ratio", l1-x0,
     $ l12-x02, s1, s1/(l1-x0)
      write(6,'(a35,4f14.6)') "Y in/out size + gain and ratio", l2-y0,
     $ l22-y02, s2, s2/(l2-y0)
      write(6,'(a35,4f14.6)') "Z in/out size + gain and ratio", l3-z0,
     $ l32-z02, s3, s3/(l3-z0)

c     prepares the new header to use with cubegen

      open(unit=out_unit,file='newheader-to-use',action='write',
     $status='replace',form='formatted')

      write(out_unit,'(i5,3f14.6)') -natom,x02/1.889726,y02/1.889726,
     $z02/1.889726

      write(out_unit,'(i5,3f14.6)') n12, dx/1.889726, dummy, dummy
      write(out_unit,'(i5,3f14.6)') n22, dummy, dy/1.889726, dummy
      write(out_unit,'(i5,3f14.6)') n32, dummy, dummy, dz/1.889726

      close(out_unit)

c     outputs the expected header to be obtained from cubegen use

      open(unit=out_unit2,file='expected-newheader',action='write',
     $status='replace',form='formatted')

      write(out_unit2,'(i5,3f14.6)') -natom,x02,y02,
     $z02
      write(out_unit2,'(i5,3f14.6)') n12, dx, dummy, dummy
      write(out_unit2,'(i5,3f14.6)') n22, dummy, dy, dummy
      write(out_unit2,'(i5,3f14.6)') n32, dummy, dummy, dz

      close(out_unit2)
      
      end 
