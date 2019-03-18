subroutine QuantumMetricsLA(d_xy,a_xy)

! Computes the topological metrics from a linear algebraic analysis of the detachment/attachment
! through population analysis of (S^x)P(S^y) (P = D,A), written d_xy and a_xy

use declare

real*8 :: d_xy(nbs,nbs), a_xy(nbs,nbs)

! Allocates the arrays and initializes the quantities and arrays

allocate(fd(nbs))
allocate(fa(nbs))
allocate(fdta(nbs))
allocate(famd(nbs))
allocate(famdm(nbs))
allocate(famdp(nbs))

fd = 0.0d0 
fa = 0.0d0
fdta = 0.0d0
famd = 0.0d0
famdm = 0.0d0
famdp = 0.0d0
chimLA = 0.0d0
chipLA = 0.0d0
phiSLA = 0.0d0
phiLA = 0.0d0
psiLA = 0.0d0
theta = 0.0d0
pi = 0.0d0
trD = 0.0d0
trA = 0.0d0

! Computes the trace of detachment/attachment density matrices.

do i=1,nbs
trD = trD + d_xy(i,i)
trA = trA + a_xy(i,i)
enddo

! Computes their average.

theta = 0.5d0*(trD + trA)

! fd(i) and fa(i) are the diagonal entries of d_xy and a_xy.
! They are used to produce:
! fdta, with "t" standing for "times", so we have the ith diagonal entry of d_xy times the ith diagonal entry of a_xy;
! famd, with "m" standing for "minus", so we have the ith diagonal entry of a_xy minus the ith diagonal entry of d_xy.
! The famd results are split into:
! famdm (the last "m" standing for "minus", i.e., negative);
! famdp (the "p" standing for "plus", i.e., positive.
! The chim (-) and chip (+) are then computed from famdm and famdp.

do i=1,nbs
 fd(i) = d_xy(i,i)
 fa(i) = a_xy(i,i)
 fdta(i) = fd(i)*fa(i)
 if (fd(i) .gt. 0.0d0 .and. fa(i) .gt. 0.0d0) phiSLA = phiSLA + (dsqrt(fdta(i)))/theta
 if (fd(i) .lt. 0.0d0 .or. fa(i) .lt. 0.0d0) then
  fd(i) = 0.0d0
  fa(i) = 0.0d0
 endif
 famd(i) = fa(i) - fd(i)
 if (famd(i) .lt. 0.0d0) then
  famdm(i) = -1.0d0*famd(i)
 else
  famdp(i) = famd(i)
 endif
 chimLA = chimLA + famdm(i)
 chipLA = chipLA + famdp(i)
enddo

! The phi metric is approximated.

phiLA = 0.5d0*(chimLA + chipLA)/theta

pi = dacos(-1.0d0)

! The psiLA metric is approximated.

psiLA = 2.0d0*(datan(phiSLA/phiLA))/(pi)

if (adiab) then
 continue
else

! The quantum metrics evaluated from the D/A population analysis are then printed.

write(6,*)
write(50,*)

write(6,*) 'Integral of detachment/attachment density'
write(6,'(f11.4)') theta
write(6,*) 'Hole/particle overlap (phiS, from linear algebra)'
write(6,'(f11.4)') phiSLA
write(6,*) 'Effectively displaced charge (phi, from linear algebra)'
write(6,'(f11.4)') phiLA
write(6,*) 'psi metric (from linear algebra)'
write(6,'(f11.4)') psiLA
write(50,*) 'Integral of detachment/attachment density'
write(50,'(f11.4)') theta
write(50,*) 'Hole/particle overlap (phiS, from linear algebra)'
write(50,'(f11.4)') phiSLA
write(50,*) 'Effectively displaced charge (phi, from linear algebra)'
write(50,'(f11.4)') phiLA
write(50,*) 'psi metric (from linear algebra)'
write(50,'(f11.4)') psiLA
write(6,*)
write(50,*)
endif

if (jobtype == 'rlxy_LA') then
 if (unr) then
  if (countunr .eq. 1) then
   thetaUalpha = theta
   phisLAUalpha = phiSLA
   phiLAUalpha = phiLA
  else if (countunr .eq. 2) then
   thetaUbeta = theta
   phisLAUbeta = phiSLA
   phiLAUbeta = phiLA
  endif
 endif
endif

! The working arrays are finally deallocated.

deallocate(fd,fa,fdta,famd,famdm,famdp)

end
