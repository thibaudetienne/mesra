       subroutine rddensities(sx0,sy0,sz0,sfile1,sfile2,snatom,sn3,sn2,
     $            sn1, sdx,sdy,sdz,sinteg1,sinteg2,sintegp,sintegm,
     $            sxrp, syrp, szrp, sxrm, syrm, szrm, sx1c, sy1c, sz1c,
     $            sx2c, sy2c, sz2c, srct,sr12ct,sphis, ss1)

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

C   Reads the cube files containing the densities and performs
C   numerical integrations

       use declare

       implicit none

       character*128 :: sfile1,sfile2
       integer :: snatom, sn3, sn2, sn1 , ss1, sn1p
       real*8 :: sdx, sdy, sdz, sinteg1, sinteg2, sintegp, sintegm,
     $           sxrp, syrp, szrp, sxrm, syrm, szrm, 
     $           sx1c, sy1c, sz1c, sx2c, sy2c, sz2c,
     $           srct,sr12ct, sphis, sx0, sy0, sz0
       
       if (sn1 .lt. 0) then
        sn1p = (-1)*sn1
       else
        sn1p = sn1
       endif

       open(10,file=sfile1,form='formatted')
       open(11,file=sfile2,form='formatted')

       do i=1,snatom+6
        read(10,*)
        read(11,*)
       enddo

       sinteg1 = 0.0D1
       sinteg2 = 0.0D1

       allocate (p1(sn3,sn2,sn1p))
       allocate (p2(sn3,sn2,sn1p))

c      Reads densities from cubes 1 and 2

       do i1 = 1,sn1p
         do i2 = 1,sn2
          read(10,'(6E13.5)') (p1(i3,i2,i1) , i3=1,sn3)
          read(11,'(6E13.5)') (p2(i3,i2,i1) , i3=1,sn3)
         enddo
       enddo

c      Performs the numerical integration of cubes 1 and 2 densities
c      over all the space span by the cube

       do i1=1,sn1p
         do i2 = 1,sn2
           do i3 = 1,sn3
           sinteg1 = sinteg1 + dabs(p1(i3,i2,i1))*sdx*sdy*sdz
           sinteg2 = sinteg2 + dabs(p2(i3,i2,i1))*sdx*sdy*sdz
           enddo
         enddo
       enddo

       close(10)
       close(11)

      allocate (diff(sn3,sn2,sn1p))
      allocate (pp(sn3,sn2,sn1p))
      allocate (pm(sn3,sn2,sn1p))

      pp = 0.0d01
      pm = 0.0d01

c     Computes the difference between the densities of cubes 1 and 2
c     and performs a segregation over the sign of the resulting entries

      do i1=1,sn1p
       do i2=1,sn2
        do i3=1,sn3
      diff(i3,i2,i1) = p2(i3,i2,i1) - p1(i3,i2,i1)
         if (diff(i3,i2,i1) .gt. 0.0d1) then
          pp(i3,i2,i1) = diff(i3,i2,i1)
          pm(i3,i2,i1) = 0.0d1
         else
          pm(i3,i2,i1) = diff(i3,i2,i1)
          pp(i3,i2,i1) = 0.0d1
         endif
        enddo
       enddo
      enddo

C     Starts the integration of the densities, the computation of the
C     charge-transfer integral, the densities centroids and
C     inter-centroid distances together with the phi_S computation

      allocate (x(sn3,sn2,sn1p))
      allocate (y(sn3,sn2,sn1p))
      allocate (z(sn3,sn2,sn1p))

      sphis = 0.0d1

      sintegp = 0.0d1
      sintegm = 0.0d1

      sxrp = 0.0d1
      syrp = 0.0d1
      szrp = 0.0d1
      sxrm = 0.0d1
      syrm = 0.0d1
      szrm = 0.0d1
      sx1c = 0.0d1
      sy1c = 0.0d1
      sz1c = 0.0d1
      sx2c = 0.0d1
      sy2c = 0.0d1
      sz2c = 0.0d1

      do i1=1,sn1p
       do i2=1,sn2
        do i3=1,sn3

      x(i3,i2,i1)=sx0+(i1-1)*sdx 
      y(i3,i2,i1)=sy0+(i2-1)*sdy
      z(i3,i2,i1)=sz0+(i3-1)*sdz

            sintegp = sintegp + pp(i3,i2,i1)*sdx*sdy*sdz
            sintegm = sintegm + pm(i3,i2,i1)*sdx*sdy*sdz

            sxrp = sxrp + 
     $(x(i3,i2,i1)*pp(i3,i2,i1)*sdx*sdy*sdz)
            sxrm = sxrm + 
     $(x(i3,i2,i1)*pm(i3,i2,i1)*sdx*sdy*sdz)

            syrp = syrp + 
     $(y(i3,i2,i1)*pp(i3,i2,i1)*sdx*sdy*sdz)
            syrm = syrm + 
     $(y(i3,i2,i1)*pm(i3,i2,i1)*sdx*sdy*sdz)

            szrp = szrp + 
     $(z(i3,i2,i1)*pp(i3,i2,i1)*sdx*sdy*sdz)
            szrm = szrm + 
     $(z(i3,i2,i1)*pm(i3,i2,i1)*sdx*sdy*sdz)

          enddo 
        enddo
      enddo

      sxrm = sxrm/sintegm
      syrm = syrm/sintegm
      szrm = szrm/sintegm
      sxrp = sxrp/sintegp
      syrp = syrp/sintegp
      szrp = szrp/sintegp

      srct = dsqrt((sxrp - sxrm)**2 + (syrp - syrm)**2
     $ + (szrp - szrm)**2)

      if (ss1 .eq. 1) then
      
      do i1 = 1,sn1p
        do i2 = 1,sn2
          do i3 = 1,sn3

            sx1c = sx1c + 
     $(x(i3,i2,i1)*p1(i3,i2,i1)*sdx*sdy*sdz)
            sx2c = sx2c + 
     $(x(i3,i2,i1)*p2(i3,i2,i1)*sdx*sdy*sdz)

            sy1c = sy1c + 
     $(y(i3,i2,i1)*p1(i3,i2,i1)*sdx*sdy*sdz)
            sy2c = sy2c + 
     $(y(i3,i2,i1)*p2(i3,i2,i1)*sdx*sdy*sdz)

            sz1c = sz1c + 
     $(z(i3,i2,i1)*p1(i3,i2,i1)*sdx*sdy*sdz)
            sz2c = sz2c + 
     $(z(i3,i2,i1)*p2(i3,i2,i1)*sdx*sdy*sdz)

      sphis = sphis +
     $dsqrt(dabs(p1(i3,i2,i1))*dabs(p2(i3,i2,i1)))
     $*sdx*sdy*sdz/((sinteg1 + sinteg2)/2)
c
          enddo
        enddo
      enddo

      sx1c = sx1c/sinteg1
      sy1c = sy1c/sinteg1
      sz1c = sz1c/sinteg1

      sx2c = sx2c/sinteg2
      sy2c = sy2c/sinteg2
      sz2c = sz2c/sinteg2

      sr12ct = dsqrt((sx1c - sx2c)**2 + (sy1c - sy2c)**2
     $ + (sz1c - sz2c)**2)

      deallocate(p1)
      deallocate(p2)
      deallocate(diff)
      deallocate(pp)
      deallocate(pm)
      deallocate(x)
      deallocate(y)
      deallocate(z)

      else
   
      deallocate(p2)
      deallocate(diff)
      deallocate(pp)
      deallocate(pm)
      deallocate(x)
      deallocate(y)
      deallocate(z)
      
      endif
      
      return 
      end 
