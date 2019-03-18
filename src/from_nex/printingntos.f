      subroutine printingntos(sn1,sinteg1, sinteg2, sintegm, sintegp,
     $ sxrm, syrm, szrm, sxrp, syrp, szrp, sx1c, sx2c, sy1c, sy2c, sz1c,
     $ sz2c, srct, sr12ct, sphis)

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
      integer :: sn1
      real*8 :: sinteg1, sinteg2, sintegm, sintegp, sxrm, syrm, szrm,
     $ sxrp, syrp, szrp, sx1c, sy1c, sz1c,
     $ sx2c, sy2c, sz2c, srct, sr12ct, sphis, ata, sss12
      character*128 :: string,string1,string2

      ata=0.52917725

c     This subroutine prints the elements computed previously  

      sss12=0.5d0*(sinteg1+sinteg2)

      string='1. Squared NTOs integration'
      write(6,'(a60)') string
      write(6,*) " " 
      string='Squared Occupied NTO integral value is'
      write(6,'(a60,f10.2)') string,
     $                      sinteg1
      write(6,*) " " 
      string='Squared Virtual NTO integral value is '
      write(6,'(a60,f10.2)') string,
     $           sinteg2
      write(6,*) " " 
      string='2. Charge from squared Occ/Virt densities'
      write(6,'(a60)') string
      write(6,*) " " 
      string='chi- from negative density variation is '
      write(6,'(a60,f10.2)') string,
     $           (-1)*sintegm
      write(6,*) " " 
      string='chi+ from positive density variation is '
      write(6,'(a60,f10.2)') string,
     $           sintegp
      write(6,*) " " 
      string='chi and phi indices values'
      write(6,'(a60,f10.2,f13.2)') string,
     $           0.5d0*(sintegp-sintegm), 0.5d0*(sintegp-sintegm)/sss12
      write(6,*) " " 
      string1='3. Centroid coordinates (bohr) from p-/p+ densities, '
      string2='zeta(+/-)'
      write(6,'(a53,a9)') string1,string2
      write(6,*) " " 
      string='Negative centroid coordinates'
      write(6,'(a60,3f13.5)') string,
     $ sxrm, syrm, szrm
      write(6,*) " " 
      string='Positive centroid coordinates'
      write(6,'(a60,3f13.5)') string,
     $                      sxrp, syrp, szrp
      write(6,*) " " 
      string="CT distance from p-/p+ in Bohr and Angstrom"
      write(6,'(a60,f13.5,f13.5)') string,
     $srct,srct*ata
      write(6,*) " " 
      string='4. Centroid coordinates (bohr) from squared NTOs, zeta'
      write(6,'(a60)') string
      write(6,*) " " 
      string='Squared Occ NTO centroid coordinates'
      write(6,'(a60,3f13.5)') string,
     $       sx1c, sy1c, sz1c
      write(6,*) " " 
      string='Squared Virt NTO centroid coordinates'
      write(6,'(a60,3f13.5)') string,
     $        sx2c, sy2c, sz2c
      write(6,*) " " 
      string='CT distance from Occ/Virt in Bohr and Angstrom'
      write(6,'(a60,2f13.5)') string,
     $sr12ct,sr12ct*ata
      write(6,*) " " 
      string='5. phi_S (Occ/Virt overlap) and psi indices calculation'
      write(6,'(a60)') string
      write(6,*) " " 
      string='phi_S index value'
      write(6,'(a60,f10.2)') string,
     $ sphis
      write(*,*) " " 
      string='psi index value'
      write(6,'(a52,f18.2)') string,
     $(2.0d0/dacos(-1.0d0))*atan(sss12*sphis/(0.5d0*(sintegp-sintegm)))
      write(6,*) " " 
      string1='NB: Centroid coordinates are given with respect '
      string2='to the cube used for the numerical integration'
      write(6,'(a48,a56)') string1,string2
      write(*,*) " " 
      end subroutine

