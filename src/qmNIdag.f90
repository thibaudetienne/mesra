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

write(6,*) 'theta_z, eta, (lambda_dag)**(1+thetaZ), phiSdag, phidag, psidag'
write(6,'(6f10.4)') theta_z,eta,lambda_dag, phiS_relaxed,phi_relaxed,psi_relaxed
write(50,*) 'theta_z, eta, (lambda_dag)**(1+thetaZ), phiSdag, phidag, psidag'
write(50,'(6f10.4)') theta_z,eta,lambda_dag, phiS_relaxed,phi_relaxed,psi_relaxed

end
