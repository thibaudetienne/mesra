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

 call zvec_main 

if (unr) then

write(6,*)
write(6,*) '3) Alpha relaxed hole/particle overlap integral'
write(6,*)
write(50,*)
write(50,*) '3) Alpha relaxed hole/particle overlap integral'
write(50,*)

 phiSPArlxalpha =  phiSPAUalpha + (1.0d0-phiSPAUalpha)*(zcoefalpha/(thetaUalpha+zcoefalpha))
! phiPArlxalpha = phiPAUalpha - phiPAUalpha*(zcoefalpha/(thetaUalpha+zcoefalpha))
! psiPArlxalpha = 2.0d0*(datan(phiSPArlxalpha/phiPArlxalpha))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation',thetaUalpha
 write(6,'(f13.4)') (zcoefalpha/(thetaUalpha+zcoefalpha))
 write(6,*) 'Relaxed/particle overlap (phiS, from PA)'
 write(6,'(f13.4)') phiSPArlxalpha
! write(6,*) 'Z-vector-adapted D/A contribution to the net charge displacement (phi, from PA)'
! write(6,'(f13.4)') phiPArlxalpha
! write(6,*) 'Generalized psi (from PA)'
! write(6,'(f13.4)') psiPArlxalpha
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (zcoefalpha/(thetaUalpha+zcoefalpha))
 write(50,*) 'Relaxed/particle overlap (phiS, from PA)'
 write(50,'(f13.4)') phiSPArlxalpha
! write(50,*) 'Z-vector-adapted D/A contribution to the net charge displacement (phi, from PA)'
! write(50,'(f13.4)') phiPArlxalpha
! write(50,*) 'Generalized psi (from PA)'
! write(50,'(f13.4)') psiPArlxalpha

write(6,*)
write(6,*) '4) Beta relaxed hole/particle overlap integral'
write(6,*)
write(50,*)
write(50,*) '4) Beta relaxed hole/particle overlap integral'
write(50,*)

 phiSPArlxbeta =  phiSPAUbeta + (1.0d0-phiSPAUbeta)*((zcoefbeta/(thetaUbeta+zcoefbeta)))
! phiPArlxbeta = phiPAUbeta - phiPAUbeta*((zcoefbeta/(thetaUbeta+zcoefbeta)))
! psiPArlxbeta = 2.0d0*(datan(phiSPArlxbeta/phiPArlxbeta))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (zcoefbeta/(thetaUbeta+zcoefbeta))
 write(6,*) 'Relaxed/particle overlap (phiS, from PA)'
 write(6,'(f13.4)') phiSPArlxbeta
! write(6,*) 'Z-vector-adapted D/A contribution to the net charge displacement (phi, from PA)'
! write(6,'(f13.4)') phiPArlxbeta
! write(6,*) 'Relaxed psi descriptor (from PA)'
! write(6,'(f13.4)') psiPArlxbeta
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (zcoefbeta/(thetaUbeta+zcoefbeta))
 write(50,*) 'Relaxed/particle overlap (phiS, from PA)'
 write(50,'(f13.4)') phiSPArlxbeta
! write(50,*) 'Z-vector-adapted D/A contribution to the net charge displacement (phi, from PA)'
! write(50,'(f13.4)') phiPArlxbeta
! write(50,*) 'Relaxed psi descriptor (from PA)'
! write(50,'(f13.4)') psiPArlxbeta

else

write(6,*)
write(6,*) '3) Relaxed hole/particle overlap integral'
write(6,*)
write(50,*)
write(50,*) '3) Relaxed hole/particle overlap integral'
write(50,*)

 phiSPArlx =  phiSPAU + (1.0d0-phiSPAU)*(zcoef/(thetaU+zcoef))
! phiPArlx = phiPAU - phiPAU*(zcoef/(thetaU+zcoef))
! psiPArlx = 2.0d0*(datan(phiSPArlx/phiPArlx))/(pi)
 
 write(6,*) 'Relative amplitude of the relaxation'
 write(6,'(f13.4)') (zcoef/(thetaU+zcoef))
 write(6,*) 'Relaxed hole/particle overlap (phiS, from PA)'
 write(6,'(f13.4)') phiSPArlx
! write(6,*) 'Z-vector-adapted D/A contribution to the net charge displacement (phi, from PA)'
! write(6,'(f13.4)') phiPArlx
! write(6,*) 'Relaxed psi descriptor (from PA)'
! write(6,'(f13.4)') psiPArlx
 write(50,*) 'Relative amplitude of the relaxation'
 write(50,'(f13.4)') (zcoef/(thetaU+zcoef))
 write(50,*) 'Relaxed/particle overlap (phiS, from PA)'
 write(50,'(f13.4)') phiSPArlx
! write(50,*) 'Z-vector-adapted D/A contribution to the net charge displacement (phi, from PA)'
! write(50,'(f13.4)') phiPArlx
! write(50,*) 'Relaxed psi descriptor (from PA)'
! write(50,'(f13.4)') psiPArlx


!!!!Gabriel Breuil 11-04-2019
!
!
!write(6,*)
!write(6,*) '3) Relaxed density-based descriptors'
!write(6,*)
!write(50,*)
!write(50,*) '3) Relaxed density-based descriptors'
!write(50,*)
!
! phiSPAnewdagger =  phiSPAU + (1-phiSPAU)*((zcoef/(thetaU+zcoef)))
! phiPAnewdagger = phiPAU - phiPAU*((zcoef/(thetaU+zcoef)))
! psiPAnewdagger = 2.0d0*(datan(phiSPAnewdagger/phiPAnewdagger))/(pi)
! 
! write(6,*) 'Relative amplitude of the relaxation'
! write(6,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
! write(6,*) 'Relaxed hole/particle overlap (relaxed phiS, from PA)'
! write(6,'(f13.4)') phiSPAnewdagger
! write(6,*) 'Relaxed charge displacement within D/A (relaxed phi, from PA)'
! write(6,'(f13.4)') phiPAnewdagger
! write(6,*) 'Relaxed psi (relaxed psi, from PA)'
! write(6,'(f13.4)') psiPAnewdagger
! write(50,*) 'Relative amplitude of the relaxation'
! write(50,'(f13.4)') (zcoef/(thetaU+zcoef))**(1.0d0+zcoef)
! write(50,*) 'Relaxed hole/particle overlap (relaxed phiS, from PA)'
! write(50,'(f13.4)') phiSPAnewdagger
! write(50,*) 'Relaxed D/A contribution to the net displaced charge (relaxed phi, from PA)'
! write(50,'(f13.4)') phiPAnewdagger
! write(50,*) 'Relaxed psi (from PA)'
! write(50,'(f13.4)') psiPAnewdagger
!
!!!!END Gabriel Breuil 04-11-2019

endif

write(6,*)

end
