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

 phiSLAdagalpha =  phiSLAUalpha + (1-phiSLAUalpha)*((thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha))
 phiLAdagalpha = phiLAUalpha - phiLAUalpha*((thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha))
 psiLAdagalpha = 2.0d0*(datan(phiSLAdagalpha/phiLAdagalpha))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha)
 write(6,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAdagalpha
 write(6,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(6,'(f13.4)') phiLAdagalpha
 write(6,*) 'Generalized psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAdagalpha
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (thetaZalpha/(thetaUalpha+thetaZalpha))**(1.0d0+thetaZalpha)
 write(50,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAdagalpha
 write(50,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(50,'(f13.4)') phiLAdagalpha
 write(50,*) 'Generalized psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAdagalpha

write(6,*)
write(6,*) '4) Beta relaxed descriptors'
write(6,*)
write(50,*)
write(50,*) '4) Beta relaxed descriptors'
write(50,*)

 phiSLAdagbeta =  phiSLAUbeta + (1-phiSLAUbeta)*((thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta))
 phiLAdagbeta = phiLAUbeta - phiLAUbeta*((thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta))
 psiLAdagbeta = 2.0d0*(datan(phiSLAdagbeta/phiLAdagbeta))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta)
 write(6,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAdagbeta
 write(6,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(6,'(f13.4)') phiLAdagbeta
 write(6,*) 'Generalized psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAdagbeta
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (thetaZbeta/(thetaUbeta+thetaZbeta))**(1.0d0+thetaZbeta)
 write(50,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAdagbeta
 write(50,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(50,'(f13.4)') phiLAdagbeta
 write(50,*) 'Generalized psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAdagbeta

else

write(6,*)
write(6,*) '3) Relaxed descriptors with old method'
write(6,*)
write(50,*)
write(50,*) '3) Relaxed descriptors with old method'
write(50,*)

 phiSLAdag =  phiSLAU + (1-phiSLAU)*((thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ))
 phiLAdag = phiLAU - phiLAU*((thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ))
 psiLAdag = 2.0d0*(datan(phiSLAdag/phiLAdag))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ)
 write(6,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAdag
 write(6,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(6,'(f13.4)') phiLAdag
 write(6,*) 'Generalized psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAdag
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (thetaZ/(thetaU+thetaZ))**(1.0d0+thetaZ)
 write(50,*) 'Generalized hole/particle overlap (phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAdag
 write(50,*) 'Generalized effectively displaced charge (phi, from linear algebra)'
 write(50,'(f13.4)') phiLAdag
 write(50,*) 'Generalized psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAdag


!!!Gabriel Breuil 11-04-2019


write(6,*)
write(6,*) '3) Relaxed descriptors with new method'
write(6,*)
write(50,*)
write(50,*) '3) Relaxed descriptors with new method'
write(50,*)

 phiSLAnewdag =  phiSLAU + (1-phiSLAU)*((zcoef/(thetaU+zcoef))**(1.0d0+zcoef))
 phiLAnewdag = phiLAU - phiLAU*((zcoef/(thetaU+zcoef))**(1.0d0+zcoef))
 psiLAnewdag = 2.0d0*(datan(phiSLAnewdag/phiLAnewdag))/(pi)
 
 write(6,*) 'Relative amplitude of the new relaxation'
 write(6,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
 write(6,*) 'Generalized hole/particle overlap (new phiS, from linear algebra)'
 write(6,'(f13.4)') phiSLAnewdag
 write(6,*) 'Generalized effectively displaced charge (new phi, from linear algebra)'
 write(6,'(f13.4)') phiLAnewdag
 write(6,*) 'Generalized new psi metric (from linear algebra)'
 write(6,'(f13.4)') psiLAnewdag
 write(50,*) 'Relative amplitude of the new relaxation'
 write(50,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
 write(50,*) 'Generalized hole/particle overlap (new phiS, from linear algebra)'
 write(50,'(f13.4)') phiSLAnewdag
 write(50,*) 'Generalized effectively displaced charge (new phi, from linear algebra)'
 write(50,'(f13.4)') phiLAnewdag
 write(50,*) 'Generalized new psi metric (from linear algebra)'
 write(50,'(f13.4)') psiLAnewdag

!!!END Gabriel Breuil 04-11-2019

endif

write(6,*)

end
