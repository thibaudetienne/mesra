subroutine qmNIdag

! Computes and outputs the relaxed descriptors, when the unrelaxed ones are read
! from the log file of a numerical integration calculation.

use declare

pi = dacos(-1.0d0)
lambda_dag = 0.0d0

phiS_relaxed =  phiS_unrelaxed + (1.0d0-phiS_unrelaxed)*((theta_z/(theta_unrelaxed+theta_z))**(1.0d0+theta_z))
phi_relaxed = phi_unrelaxed - phi_unrelaxed*((theta_z/(theta_unrelaxed+theta_z))**(1.0d0+theta_z))
psi_relaxed = 2.0d0*(datan(phiS_relaxed/phi_relaxed))/(pi)
lambda_dag = (theta_z/(theta_unrelaxed+theta_z))**(1.0d0+theta_z)
eta = eta*0.5d0

!!!Gabriel Breuil 12-04-2019
write(6,*) 'Relaxed descriptors obtained with an old qmnidag calculation'
write(50,*) 'Relaxed descriptors obtained with an old qmnidag calculation'
!!!End Gabriel Breuil

write(6,*) 'theta_z, eta, (lambda_dag)**(1+thetaZ), phiSdag, phidag, psidag'
write(6,'(6f10.4)') theta_z,eta,lambda_dag, phiS_relaxed,phi_relaxed,psi_relaxed
write(50,*) 'theta_z, eta, (lambda_dag)**(1+thetaZ), phiSdag, phidag, psidag'
write(50,'(6f10.4)') theta_z,eta,lambda_dag, phiS_relaxed,phi_relaxed,psi_relaxed

!!!Gabriel Breuil 12-04-2019

newphiS_relaxed=phiS_unrelaxed+(1.0d0-phiS_unrelaxed)*((zcoef/(theta_unrelaxed+zcoef))**(1.0d0+zcoef))
newphi_relaxed=phi_unrelaxed+(1.0d0-phi_unrelaxed)*((zcoef/(theta_unrelaxed+zcoef))**(1.0d0+zcoef))
newpsi_relaxed=2.0d0*(datan(newphiS_relaxed/newphi_relaxed))/(pi)
newlambda_dag=(zcoef/(theta_unrelaxed+zcoef))**(1.0d0+zcoef)

write(6,*) 'Relaxed descriptors obtained with a new qmnidag calculation'
write(50,*) 'Relaxed descriptors obtained with a new qmnidag calculation'
write(6,*) 'zcoef, newlambda_dag, newphiSdag, newphidag, newpsidag'
write(6,'(6f10.4)') zcoef,newlambda_dag,newphiS_relaxed,newphi_relaxed,newpsi_relaxed
write(50,*) 'zcoef, newlambda_dag, newphiSdag, newphidag, newpsidag'
write(50,'(6f10.4)') zcoef,newlambda_dag,newphiS_relaxed,newphi_relaxed,newpsi_relaxed

!!!End Gabriel Breuil

end
