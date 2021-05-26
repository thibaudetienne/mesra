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

subroutine qmNIrlx

! Computes and outputs the relaxed descriptors, when the unrelaxed ones are read
! from the log file of a numerical integration calculation.

use declare

pi = dacos(-1.0d0)

phiS_relaxed=phiS_unrelaxed+(1.0d0-phiS_unrelaxed)*((eta/(theta_unrelaxed+eta)))
!phi_relaxed=phi_unrelaxed-phi_unrelaxed*((eta/(theta_unrelaxed+eta)))
!psi_relaxed=2.0d0*(datan(phiS_relaxed/phi_relaxed))/(pi)
Rcoef=(eta/(theta_unrelaxed+eta))

write(6,*) 'Relaxed descriptors, from numerical integration'
write(50,*) 'Relaxed descriptors, from numerical integration'

write(6,*)
write(50,*)

write(6,*) 'eta and R'
write(6,'(2f10.4)') eta,Rcoef

write(6,*) 'Relaxed phiS'
write(6,'(f10.4)') phiS_relaxed

!write(6,*) 'phi_relaxed'
!write(6,'(f10.4)') phi_relaxed

!write(6,*) 'psi_relaxed'
!write(6,'(f10.4)') psi_relaxed

write(50,*) 'eta and R'
write(50,'(2f10.4)') eta, Rcoef
write(50,*) 'Relaxed phiS'
write(50,'(f10.4)') phiS_relaxed
!write(50,*) 'phi_relaxed'
!write(50,'(f10.4)') phi_relaxed
!write(50,*) 'psi_relaxed'
!write(50,'(f10.4)') psi_relaxed

end
