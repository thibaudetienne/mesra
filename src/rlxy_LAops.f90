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

subroutine rlxy_LAops

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
 phiSLAU = phiSLA
 phiLAU = phiLA
 psiLAU = psiLA
endif

 LA = .false.
 scanLA = .false.

write(6,*) '2) Relaxation'
write(6,*) 
write(50,*) '2) Relaxation'
write(50,*) 

 jobtype = 'daz'
 subjobtype = 'rlxy_LA'
 call zvec_main 

if (unr) then

write(6,*)
write(6,*) '3) Alpha relaxed descriptors'
write(6,*)
write(50,*)
write(50,*) '3) Alpha relaxed descriptors'
write(50,*)

 phiSLAdaggeralpha =  phiSLAUalpha + (1-phiSLAUalpha)*((thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha))
 phiLAdaggeralpha = phiLAUalpha - phiLAUalpha*((thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha))
 psiLAdaggeralpha = 2.0d0*(datan(phiSLAdaggeralpha/phiLAdaggeralpha))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha)
 write(6,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAdaggeralpha
 write(6,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(6,'(f13.4)') phiLAdaggeralpha
 write(6,*) 'Generalized psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAdaggeralpha
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha)
 write(50,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAdaggeralpha
 write(50,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(50,'(f13.4)') phiLAdaggeralpha
 write(50,*) 'Generalized psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAdaggeralpha

write(6,*)
write(6,*) '4) Beta relaxed descriptors'
write(6,*)
write(50,*)
write(50,*) '4) Beta relaxed descriptors'
write(50,*)

 phiSLAdaggerbeta =  phiSLAUbeta + (1-phiSLAUbeta)*((thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta))
 phiLAdaggerbeta = phiLAUbeta - phiLAUbeta*((thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta))
 psiLAdaggerbeta = 2.0d0*(datan(phiSLAdaggerbeta/phiLAdaggerbeta))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta)
 write(6,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAdaggerbeta
 write(6,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(6,'(f13.4)') phiLAdaggerbeta
 write(6,*) 'Generalized psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAdaggerbeta
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta)
 write(50,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAdaggerbeta
 write(50,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(50,'(f13.4)') phiLAdaggerbeta
 write(50,*) 'Generalized psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAdaggerbeta

else

write(6,*)
write(6,*) '3) Relaxed descriptors with old method'
write(6,*)
write(50,*)
write(50,*) '3) Relaxed descriptors with old method'
write(50,*)

 phiSLAdagger =  phiSLAU + (1-phiSLAU)*((thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ))
 phiLAdagger = phiLAU - phiLAU*((thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ))
 psiLAdagger = 2.0d0*(datan(phiSLAdagger/phiLAdagger))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ)
 write(6,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAdagger
 write(6,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(6,'(f13.4)') phiLAdagger
 write(6,*) 'Generalized psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAdagger
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ)
 write(50,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAdagger
 write(50,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(50,'(f13.4)') phiLAdagger
 write(50,*) 'Generalized psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAdagger


!!!Gabriel Breuil 11-04-2019


write(6,*)
write(6,*) '3) Relaxed descriptors with new method'
write(6,*)
write(50,*)
write(50,*) '3) Relaxed descriptors with new method'
write(50,*)

 phiSLAnewdagger =  phiSLAU + (1-phiSLAU)*((zcoef/(thetaU+zcoef))**(1.0d0+zcoef))
 phiLAnewdagger = phiLAU - phiLAU*((zcoef/(thetaU+zcoef))**(1.0d0+zcoef))
 psiLAnewdagger = 2.0d0*(datan(phiSLAnewdagger/phiLAnewdagger))/(pi)
 
 write(6,*) 'Relative amplitude of the new relaxation'
 write(6,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
 write(6,*) 'Generalized hole/particle overlap (new phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAnewdagger
 write(6,*) 'Generalized effectively displaced charge (new phi, from linear algebra)'
 write(6,'(f13.4)') phiLAnewdagger
 write(6,*) 'Generalized new psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAnewdagger
 write(50,*) 'Relative amplitude of the new relaxation'
 write(50,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
 write(50,*) 'Generalized hole/particle overlap (new phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAnewdagger
 write(50,*) 'Generalized effectively displaced charge (new phi, from linear algebra)'
 write(50,'(f13.4)') phiLAnewdagger
 write(50,*) 'Generalized new psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAnewdagger

!!!END Gabriel Breuil 04-11-2019

endif

write(6,*)

end
