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

subroutine adiab_z

use declare

if (countunr .eq. 1) then
open(10,file='adiabZalpha.ac',form='formatted')
write(10,*) '# xi, Trace of ZZ^dagger, R, Relaxed phiS (Lowdin)'
write(10,*)
write(6,*)
write(6,*) '# xi, Trace of ZZ^dagger, R, Relaxed phiS (Lowdin)'
open(11,file='acVec_Zalpha.ac',form='formatted')
write(11,*) '# xi, i, Auto-correlation vector'
write(11,*)
endif

if (countunr .eq. 2) then
write(6,*)
write(6,*) '# xi, Trace of ZZ^dagger, R, Relaxed phiS (Lowdin)'
open(10,file='adiabZbeta.ac',form='formatted')
write(10,*) '# xi, Trace of ZZ^dagger, R, Relaxed phiS (Lowdin)'
write(10,*)
open(11,file='acVec_Zbeta.ac',form='formatted')
write(11,*) '# xi, i, Auto-correlation vector'
write(11,*)
endif

if (countunr .eq. 3) then
open(10,file='adiabZ.ac',form='formatted')
write(10,*) '# xi, Trace of ZZ^dagger, R, Relaxed phiS (Lowdin)'
write(10,*)
write(6,*)
write(6,*) '# xi, Trace of ZZ^dagger, R, Relaxed phiS (Lowdin)'
open(11,file='acVec_Z.ac',form='formatted')
write(11,*) '# xi, i, Auto-correlation vector'
write(11,*)
endif

allocate(zvecxi(norb,norb))
allocate(U0(norb,norb))
allocate(U0t(norb,norb))
allocate(U0tU(norb,norb))
allocate(Uvec(norb,norb))
allocate(gD(norb,norb))
allocate(lvec(norb))

PA = .true.

xi = 0.0d0

write(6,*)

do iteration=0,100

xi = iteration*0.01d0
zvecxi = xi*zvec

if (countunr .eq. 1) pxalpharelaxed = pxalpha + zvecxi
if (countunr .eq. 2) pxbetarelaxed = pxbeta + zvecxi
if (countunr .eq. 3) pxrelaxed = px + zvecxi

if (allocated(zzd_zdz)) deallocate(zzd_zdz)
allocate(zzd_zdz(norb,norb))

zzd_zdz = matmul(zvecxi,zvecxi)

eta = 0.0d0

do i=1,norb
eta = eta + zzd_zdz(i,i)
enddo

eta = 0.5d0*eta

if (countunr .eq. 1) gD = pxalpharelaxed - palpha
if (countunr .eq. 2) gD = pxbetarelaxed - pbeta
if (countunr .eq. 3) gD = pxrelaxed - p

if (countunr .eq. 2) C = Cb
if (countunr .eq. 2) Cdagger = Cbdagger

if (countunr .eq. 1) call det_at(gD,'adiabAlpha',norb,Uvec,lvec)
if (countunr .eq. 2) call det_at(gD,'adiabBeta',norb,Uvec,lvec)
if (countunr .eq. 3) call det_at(gD,'adiab',norb,Uvec,lvec)

if (iteration .eq. 0) phiSPA0 = phiSPA
if (iteration .eq. 0) phiPA0 = phiPA
if (iteration .eq. 0) theta0 = theta
if (iteration .eq. 0) U0 = Uvec

! Implemented by GB
if (iteration .eq. 0) then
phiS_relaxedPA0=phiSPA0+(1-phiSPA0)*((eta/(theta0+eta)))
!phi_relaxedPA0=phiPA0+(1-phiPA0)*((eta/(theta0+eta)))
endif
! End GB

U0t = transpose(U0)
U0tU = matmul(U0t,Uvec)

do i=1,norb
write(11,'(f10.2,i8,f12.8)') xi,i,dabs(U0tU(i,i))
enddo

if (PA) t_PA = .true.
if (scanPA) t_scanPA = .true.

!t_PA = .false.
!t_scanPA = .false.
!if (PA) t_PA = .true.
!if (scanPA) t_scanPA = .true.
!
!PA = .false.
!scanPA = .false.
!call det_at(zvecxi,'adiab',norb,Uvec,lvec)
!if (t_PA) PA = .true.
!if (t_scanPA) scanPA = .true.
!
!trDZ = 0.0d0
!trAZ = 0.0d0
!thetaZ = 0.0d0
!
!do i=1,norb
! if (lvec(i) .lt. 0.0d0) trDZ = trDZ - lvec(i)
! if (lvec(i) .gt. 0.0d0) trAZ = trAZ + lvec(i)
!enddo
!
!thetaZ = 0.5d0*(trDZ + trAZ)

pi = dacos(-1.0d0)

!phiSPAdagger =  phiSPA0 + (1-phiSPA0)*((thetaZ/(theta0+thetaZ))**(1+thetaZ))
!phiPAdagger =  phiPA0 - phiPA0*((thetaZ/(theta0+thetaZ))**(1+thetaZ))
!psiPAdagger = 2.0d0*datan(phiSPAdagger/phiPAdagger)/(pi)

! Implemented by GB
phiS_relaxedPA=phiS_relaxedPA0+(1-phiS_relaxedPA0)*((eta/(theta0+eta)))
!phi_relaxedPA=phi_relaxedPA0+(1-phi_relaxedPA0)*((eta/(theta0+eta)))
!psi_relaxedPA=2.0d0*datan(phiS_relaxedPA/phi_relaxedPA)/(pi)

!write(6,*) 'alter_phiSdagger', phiSPA0 + (1-phiSPA0)*((eta/(theta0+eta))**(1+eta))
!phiSdagger =  phiSPA0 + (1-phiSPA0)*((eta/(theta0+eta))**(1+eta))
!write(6,*) 'alter_phiSdagger', phiSPA0 + (1-phiSPA0)*eta/thetaZ
write(10,'(f8.2,3f8.4)') xi, eta, (eta/(theta0+eta)),phiS_relaxedPA
write(6,'(f8.2,3f8.4)') xi, eta, (eta/(theta0+eta)),phiS_relaxedPA
! End of GB implementation
enddo

deallocate(zvecxi)
deallocate(U0,U0t,Uvec,U0tU)
deallocate(gD)
deallocate(lvec)

close(10)
close(11)

end
