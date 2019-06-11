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

subroutine qmNIdagger

! Computes and outputs the relaxed descriptors, when the unrelaxed ones are read
! from the log file of a numerical integration calculation.

use declare

pi = dacos(-1.0d0)
lambda_dagger = 0.0d0

phiS_relaxed =  phiS_unrelaxed + (1.0d0-phiS_unrelaxed)*((theta_z/(theta_unrelaxed+theta_z))**(1.0d0+theta_z))
phi_relaxed = phi_unrelaxed - phi_unrelaxed*((theta_z/(theta_unrelaxed+theta_z))**(1.0d0+theta_z))
psi_relaxed = 2.0d0*(datan(phiS_relaxed/phi_relaxed))/(pi)
lambda_dagger = (theta_z/(theta_unrelaxed+theta_z))**(1.0d0+theta_z)
eta = eta*0.5d0

!!!Gabriel Breuil 12-04-2019
write(6,*) 'Relaxed descriptors obtained with an old qmnidagger calculation'
write(50,*) 'Relaxed descriptors obtained with an old qmnidagger calculation'
!!!End Gabriel Breuil

write(6,*) 'theta_z, eta, (lambda_dagger)**(1+thetaZ), phiSdagger, phidagger, psidagger'
write(6,'(6f10.4)') theta_z,eta,lambda_dagger, phiS_relaxed,phi_relaxed,psi_relaxed
write(50,*) 'theta_z, eta, (lambda_dagger)**(1+thetaZ), phiSdagger, phidagger, psidagger'
write(50,'(6f10.4)') theta_z,eta,lambda_dagger, phiS_relaxed,phi_relaxed,psi_relaxed

!!!Gabriel Breuil 12-04-2019

newphiS_relaxed=phiS_unrelaxed+(1.0d0-phiS_unrelaxed)*((zcoef/(theta_unrelaxed+zcoef))**(1.0d0+zcoef))
newphi_relaxed=phi_unrelaxed+(1.0d0-phi_unrelaxed)*((zcoef/(theta_unrelaxed+zcoef))**(1.0d0+zcoef))
newpsi_relaxed=2.0d0*(datan(newphiS_relaxed/newphi_relaxed))/(pi)
newlambda_dagger=(zcoef/(theta_unrelaxed+zcoef))**(1.0d0+zcoef)

write(6,*) 'Relaxed descriptors obtained with a new qmnidagger calculation'
write(50,*) 'Relaxed descriptors obtained with a new qmnidagger calculation'
write(6,*) 'zcoef, newlambda_dagger, newphiSdagger, newphidagger, newpsidagger'
write(6,'(6f10.4)') zcoef,newlambda_dagger,newphiS_relaxed,newphi_relaxed,newpsi_relaxed
write(50,*) 'zcoef, newlambda_dagger, newphiSdagger, newphidagger, newpsidagger'
write(50,'(6f10.4)') zcoef,newlambda_dagger,newphiS_relaxed,newphi_relaxed,newpsi_relaxed

!!!End Gabriel Breuil

end
