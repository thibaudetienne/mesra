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

subroutine rlxy_PAops

! 1) performs an unrelaxed calculation of the descriptors;
! 2) extracts an information from the relaxation transformation;
! 3) gathers the info into the relaxed descriptors.

! NB: alpha and beta calculations are separated for open-shell molecules.

use declare

write(6,*) '1) Unrelaxed picture of the transition'
write(6,*)
write(50,*) '1) Unrelaxed picture of the transition'
write(50,*)

 call daXY

if (unr) then
 continue
else
 thetaU = theta
 phiSPAU = phiSPA
 phiPAU = phiPA
 psiPAU = psiPA
endif

 PA = .false.
 scanPA = .false.

write(6,*) '2) Relaxation'
write(6,*) 
write(50,*) '2) Relaxation'
write(50,*) 

 jobtype = 'daz'
 subjobtype = 'rlxy_PA'
 call zvec_main 

if (unr) then

!write(6,*)
!write(6,*) '3) Alpha relaxed descriptors'
!write(6,*)
!write(50,*)
!write(50,*) '3) Alpha relaxed descriptors'
!write(50,*)
!
! phiSPAdaggeralpha =  phiSPAUalpha + (1-phiSPAUalpha)*((thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha))
! phiPAdaggeralpha = phiPAUalpha - phiPAUalpha*((thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha))
! psiPAdaggeralpha = 2.0d0*(datan(phiSPAdaggeralpha/phiPAdaggeralpha))/(pi)
! 
! write(6,*) 'Relative amplitude of the relaxation'
! write(6,'(f13.4)') (thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha)
! write(6,*) 'Generalized hole/particle overlap (phiS, from population analysis)'
! write(6,'(f13.4)') phiSPAdaggeralpha
! write(6,*) 'Generalized effectively displaced charge (phi, from population analysis)'
! write(6,'(f13.4)') phiPAdaggeralpha
! write(6,*) 'Generalized psi metric (from population analysis)'
! write(6,'(f13.4)') psiPAdaggeralpha
! write(50,*) 'Relative amplitude of the relaxation'
! write(50,'(f13.4)') (thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha)
! write(50,*) 'Generalized hole/particle overlap (phiS, from population analysis)'
! write(50,'(f13.4)') phiSPAdaggeralpha
! write(50,*) 'Generalized effectively displaced charge (phi, from population analysis)'
! write(50,'(f13.4)') phiPAdaggeralpha
! write(50,*) 'Generalized psi metric (from population analysis)'
! write(50,'(f13.4)') psiPAdaggeralpha
!
!write(6,*)
!write(6,*) '4) Beta relaxed descriptors'
!write(6,*)
!write(50,*)
!write(50,*) '4) Beta relaxed descriptors'
!write(50,*)
!
! phiSPAdaggerbeta =  phiSPAUbeta + (1-phiSPAUbeta)*((thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta))
! phiPAdaggerbeta = phiPAUbeta - phiPAUbeta*((thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta))
! psiPAdaggerbeta = 2.0d0*(datan(phiSPAdaggerbeta/phiPAdaggerbeta))/(pi)
! 
! write(6,*) 'Relative amplitude of the relaxation'
! write(6,'(f13.4)') (thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta)
! write(6,*) 'Generalized hole/particle overlap (phiS, from population analysis)'
! write(6,'(f13.4)') phiSPAdaggerbeta
! write(6,*) 'Generalized effectively displaced charge (phi, from population analysis)'
! write(6,'(f13.4)') phiPAdaggerbeta
! write(6,*) 'Generalized psi metric (from population analysis)'
! write(6,'(f13.4)') psiPAdaggerbeta
! write(50,*) 'Relative amplitude of the relaxation'
! write(50,'(f13.4)') (thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta)
! write(50,*) 'Generalized hole/particle overlap (phiS, from population analysis)'
! write(50,'(f13.4)') phiSPAdaggerbeta
! write(50,*) 'Generalized effectively displaced charge (phi, from population analysis)'
! write(50,'(f13.4)') phiPAdaggerbeta
! write(50,*) 'Generalized psi metric (from population analysis)'
! write(50,'(f13.4)') psiPAdaggerbeta
!
!else
!
!write(6,*)
!write(6,*) '3) Relaxed descriptors with old method'
!write(6,*)
!write(50,*)
!write(50,*) '3) Relaxed descriptors with old method'
!write(50,*)
!
! phiSPAdagger =  phiSPAU + (1-phiSPAU)*((thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ))
! phiPAdagger = phiPAU - phiPAU*((thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ))
! psiPAdagger = 2.0d0*(datan(phiSPAdagger/phiPAdagger))/(pi)
! 
! write(6,*) 'Relative amplitude of the relaxation'
! write(6,'(f13.4)') (thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ)
! write(6,*) 'Generalized hole/particle overlap (phiS, from population analysis)'
! write(6,'(f13.4)') phiSPAdagger
! write(6,*) 'Generalized effectively displaced charge (phi, from population analysis)'
! write(6,'(f13.4)') phiPAdagger
! write(6,*) 'Generalized psi metric (from population analysis)'
! write(6,'(f13.4)') psiPAdagger
! write(50,*) 'Relative amplitude of the relaxation'
! write(50,'(f13.4)') (thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ)
! write(50,*) 'Generalized hole/particle overlap (phiS, from population analysis)'
! write(50,'(f13.4)') phiSPAdagger
! write(50,*) 'Generalized effectively displaced charge (phi, from population analysis)'
! write(50,'(f13.4)') phiPAdagger
! write(50,*) 'Generalized psi metric (from population analysis)'
! write(50,'(f13.4)') psiPAdagger


!!!Gabriel Breuil 11-04-2019


write(6,*)
write(6,*) '3) Relaxed density-based descriptors'
write(6,*)
write(50,*)
write(50,*) '3) Relaxed density-based descriptors'
write(50,*)

 phiSPAnewdagger =  phiSPAU + (1-phiSPAU)*((zcoef/(thetaU+zcoef)))
 phiPAnewdagger = phiPAU - phiPAU*((zcoef/(thetaU+zcoef)))
 psiPAnewdagger = 2.0d0*(datan(phiSPAnewdagger/phiPAnewdagger))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
 write(6,*) 'Relaxed hole/particle overlap (relaxed phiS, from population analysis)'
 write(6,'(f13.4)') phiSPAnewdagger
 write(6,*) 'Relaxed charge displacement within D/A (relaxed phi, from population analysis)'
 write(6,'(f13.4)') phiPAnewdagger
 write(6,*) 'Relaxed psi (relaxed psi, from population analysis)'
 write(6,'(f13.4)') psiPAnewdagger
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
 write(50,*) 'Relaxed hole/particle overlap (relaxed phiS, from population analysis)'
 write(50,'(f13.4)') phiSPAnewdagger
 write(50,*) 'Relaxed D/A contribution to the net displaced charge (relaxed phi, from population analysis)'
 write(50,'(f13.4)') phiPAnewdagger
 write(50,*) 'Relaxed psi (from population analysis)'
 write(50,'(f13.4)') psiPAnewdagger

!!!END Gabriel Breuil 04-11-2019

endif

write(6,*)

end
