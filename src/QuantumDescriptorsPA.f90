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

subroutine QuantumDescriptorsPA(d_xy,a_xy)

! Computes the density-based descriptors from a population analysis of the detachment/attachment
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
chimPA = 0.0d0
chipPA = 0.0d0
phiSPA = 0.0d0
phiPA = 0.0d0
psiPA = 0.0d0
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
 if (fd(i) .gt. 0.0d0 .and. fa(i) .gt. 0.0d0) phiSPA = phiSPA + (dsqrt(fdta(i)))/theta
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
 chimPA = chimPA + famdm(i)
 chipPA = chipPA + famdp(i)
enddo

! The phi descriptor is approximated.

phiPA = 0.5d0*(chimPA + chipPA)/theta

pi = dacos(-1.0d0)

! The psi descriptor is approximated.

psiPA = 2.0d0*(datan(phiSPA/phiPA))/(pi)

if (adiab) then
 continue
else

! The density-based descriptors evaluated from the D/A population analysis are then printed.

write(6,*) 'Integral of detachment/attachment density'
write(6,'(f11.4)') theta
write(6,*) 'Hole/particle overlap (phiS, from population analysis)'
write(6,'(f11.4)') phiSPA
write(6,*) 'Fraction of D/A contributing to the net displaced charge (phi, from population analysis)'
write(6,'(f11.4)') phiPA
write(6,*) 'psi (from population analysis)'
write(6,'(f11.4)') psiPA
write(50,*) 'Integral of detachment/attachment density'
write(50,'(f11.4)') theta
write(50,*) 'Hole/particle overlap (phiS, from population analysis)'
write(50,'(f11.4)') phiSPA
write(50,*) 'Fraction of D/A contributing to the net displaced charge (phi, from population analysis)'
write(50,'(f11.4)') phiPA
write(50,*) 'psi (from population analysis)'
write(50,'(f11.4)') psiPA
write(6,*)
write(50,*)
endif

if (jobtype == 'rlxy_PA') then
 if (unr) then
  if (countunr .eq. 1) then
   thetaUalpha = theta
   phisPAUalpha = phiSPA
   phiPAUalpha = phiPA
  else if (countunr .eq. 2) then
   thetaUbeta = theta
   phisPAUbeta = phiSPA
   phiPAUbeta = phiPA
  endif
 endif
endif

! The working arrays are finally deallocated.

deallocate(fd,fa,fdta,famd,famdm,famdp)

end
