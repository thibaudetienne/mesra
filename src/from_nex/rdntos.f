      subroutine rdntos(sfile1,sn,sx0,sy0,sz0,sn1,sn2,sn3,sdx,sdy,sdz,
     $sinteg1, sinteg2, sintegm, sintegp, sxrm, syrm, szrm, sxrp, syrp, 
     $szrp, sx1c, sx2c, sy1c, sy2c, sz1c, sz2c, srct, sr12ct, sphis)

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

      use declare

      implicit none
      character*128 :: sfile1
      integer :: sn1, sn2, sn3, sn, sn1p
      real*8 :: sx0,sy0,sz0,sdx,sdy,sdz,sinteg1,sinteg2, sintegp,
     $ sintegm, sxrm, syrm, szrm, sxrp, syrp, szrp,
     $ srct, sx1c, sx2c, sy1c, sy2c, sz1c, sz2c, sr12ct, sphis

      if (sn1 .lt. 0) then
       sn1p = (-1)*sn1
      else
       sn1p = sn1
      endif

      allocate ( nto1(sn3,sn2,sn1p) )
      allocate ( nto2(sn3,sn2,sn1p) )
      allocate ( p1(sn3,sn2,sn1p) )
      allocate ( p2(sn3,sn2,sn1p) )

      open(10,file=sfile1,form='formatted')

      do i=1,sn+7
       read(10,*)
      enddo

      sphis = 0.0d1
      occ = 0.0d1
      virt = 0.0d1

      do i1=1,sn1p
       do i2=1,sn2
        read(10,'(6E13.6)') (nto1(i3,i2,i1), nto2(i3,i2,i1), i3=1,sn3)
         do i3=1,n3
          p1(i3,i2,i1) = (nto1(i3,i2,i1)**2)
          p2(i3,i2,i1) = (nto2(i3,i2,i1)**2)
          sinteg1 = sinteg1 + p1(i3,i2,i1)*sdx*sdy*sdz
          sinteg2 = sinteg2 + p2(i3,i2,i1)*sdx*sdy*sdz
         enddo
       enddo
      enddo

C     Takes the difference between the OCC/VIRT densities
C     and performs a segregation over the sign of the resulting entries

      allocate (diff(sn3,sn2,sn1p))
      allocate (pp(sn3,sn2,sn1p))
      allocate (pm(sn3,sn2,sn1p))

      pp = 0.0d01
      pm = 0.0d01

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

             sphis = sphis + dsqrt(((nto1(i3,i2,i1))**2)*((nto2(i3,i2,i
     $1))**2))*sdx*sdy*sdz

          enddo
        enddo
      enddo
     
      sx1c = sx1c/sinteg1
      sy1c = sy1c/sinteg1
      sz1c = sz1c/sinteg1

      sx2c = sx2c/sinteg2
      sy2c = sy2c/sinteg2
      sz2c = sz2c/sinteg2

      sxrm = sxrm/sintegm
      syrm = syrm/sintegm
      szrm = szrm/sintegm

      sxrp = sxrp/sintegp
      syrp = syrp/sintegp
      szrp = szrp/sintegp

      srct = dsqrt( (sxrp - sxrm)**2 + (syrp - syrm)**2
     $ + (szrp - szrm)**2)

      sr12ct = dsqrt((sx1c - sx2c)**2 + (sy1c - sy2c)**2
     $ + (sz1c - sz2c)**2)
     
      close(10)

      return
      end 
