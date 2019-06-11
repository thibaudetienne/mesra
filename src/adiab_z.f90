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
open(10,file='adiabZalpha.log',form='formatted')
write(10,*) '# xi, Trace of ZZdagger, ThetaZ, phiSLA(Lowdin), phiSLAdagger, phiLA(Lowdin), phiLAdagger, psiLA(Lowdin), psiLAdagger'
write(10,*)
open(11,file='acVec_Zbeta.log',form='formatted')
write(11,*) '# xi, i, Auto-correlation vector'
write(11,*)
endif

if (countunr .eq. 2) then
open(10,file='adiabZbeta.log',form='formatted')
write(10,*) '# xi, Trace of ZZdagger, ThetaZ, phiSLA(Lowdin), phiSLAdagger, phiLA(Lowdin), phiLAdagger, psiLA(Lowdin), psiLAdagger'
write(10,*)
open(11,file='acVec_Zbeta.log',form='formatted')
write(11,*) '# xi, i, Auto-correlation vector'
write(11,*)
endif

if (countunr .eq. 3) then
open(10,file='adiabZ.log',form='formatted')
write(10,*) '# xi, Trace of ZZdagger, ThetaZ, phiSLA(Lowdin), phiSLAdagger, phiSLAnewdagger, phiLA(Lowdin), &
               & phiLAdagger, phiLAnewdagger, psiLA(Lowdin), psiLAdagger, psiLAnewdagger'
write(10,*)
open(11,file='acVec_Z.log',form='formatted')
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

LA = .true.

xi = 0.0d0

do iteration=0,100

xi = iteration*0.01d0
zvecxi = xi*zvec
write(6,'(a3,f10.2)') 'xi ',xi

if (countunr .eq. 1) pxalpharelaxed = pxalpha + zvecxi
if (countunr .eq. 2) pxbetarelaxed = pxbeta + zvecxi
if (countunr .eq. 3) pxrelaxed = px + zvecxi

zzd_zdz = matmul(zvecxi,zvecxi)

alter_thetaZ = 0.0d0

do i=1,norb
alter_thetaZ = alter_thetaZ + zzd_zdz(i,i)
enddo

alter_thetaZ = 0.5d0*alter_thetaZ

if (countunr .eq. 1) gD = pxalpharelaxed - palpha
if (countunr .eq. 2) gD = pxbetarelaxed - pbeta
if (countunr .eq. 3) gD = pxrelaxed - p

if (countunr .eq. 2) C = Cb
if (countunr .eq. 2) Cdagger = Cbdagger

if (countunr .eq. 1) call det_at(gD,'adiabAlpha',norb,Uvec,lvec)
if (countunr .eq. 2) call det_at(gD,'adiabBeta',norb,Uvec,lvec)
if (countunr .eq. 3) call det_at(gD,'adiab',norb,Uvec,lvec)

if (iteration .eq. 0) phiSLA0 = phiSLA
if (iteration .eq. 0) phiLA0 = phiLA
if (iteration .eq. 0) theta0 = theta
if (iteration .eq. 0) U0 = Uvec
!!!Gabriel Breuil 12-04-2019
if (iteration .eq. 0) then
phiSLAnew0=phiSLA0+(1-phiSLA0)*((alter_thetaZ/(theta0+alter_thetaZ))**(1+alter_thetaZ))
phiLAnew0=phiLA0+(1-phiLA0)*((alter_thetaZ/(theta0+alter_thetaZ))**(1+alter_thetaZ))
endif
!!!End Gabriel Breuil

U0t = transpose(U0)
U0tU = matmul(U0t,Uvec)

do i=1,norb
write(11,'(f10.2,i8,f12.8)') xi,i,dabs(U0tU(i,i))
enddo

t_LA = .false.
t_scanLA = .false.
if (LA) t_LA = .true.
if (scanLA) t_scanLA = .true.

LA = .false.
scanLA = .false.
call det_at(zvecxi,'adiab',norb,Uvec,lvec)
if (t_LA) LA = .true.
if (t_scanLA) scanLA = .true.

trDZ = 0.0d0
trAZ = 0.0d0
thetaZ = 0.0d0

do i=1,norb
 if (lvec(i) .lt. 0.0d0) trDZ = trDZ - lvec(i)
 if (lvec(i) .gt. 0.0d0) trAZ = trAZ + lvec(i)
enddo

thetaZ = 0.5d0*(trDZ + trAZ)

pi = dacos(-1.0d0)

!write(6,*) 'thetaZ ',thetaZ
phiSLAdagger =  phiSLA0 + (1-phiSLA0)*((thetaZ/(theta0+thetaZ))**(1+thetaZ))
phiLAdagger =  phiLA0 - phiLA0*((thetaZ/(theta0+thetaZ))**(1+thetaZ))
psiLAdagger = 2.0d0*datan(phiSLAdagger/phiLAdagger)/(pi)
!!!Gabriel Breuil 12-04-2019
phiSLAnewdagger=phiSLAnew0+(1-phiSLAnew0)*((alter_thetaZ/(theta0+alter_thetaZ))**(1+alter_thetaZ))
phiLAnewdagger=phiLAnew0+(1-phiLAnew0)*((alter_thetaZ/(theta0+alter_thetaZ))**(1+alter_thetaZ))
psiLAnewdagger=2.0d0*datan(phiSLAnewdagger/phiLAnewdagger)/(pi)


!write(6,*) 'alter_phiSdagger', phiSLA0 + (1-phiSLA0)*((alter_thetaZ/(theta0+alter_thetaZ))**(1+alter_thetaZ))
!phiSdagger =  phiSLA0 + (1-phiSLA0)*((alter_thetaZ/(theta0+alter_thetaZ))**(1+alter_thetaZ))
!write(6,*) 'alter_phiSdagger', phiSLA0 + (1-phiSLA0)*alter_thetaZ/thetaZ

write(10,'(f10.2,11f10.4)') xi, alter_thetaZ,thetaZ,phiSLA,phiSLAdagger,phiSLAnewdagger,phiLA,&
        & phiLAdagger,phiLAnewdagger,psiLA,psiLAdagger,psiLAnewdagger
!!!End Gabriel Breuil
enddo

deallocate(zvecxi)
deallocate(U0,U0t,Uvec,U0tU)
deallocate(gD)
deallocate(lvec)

close(10)
close(11)

end
