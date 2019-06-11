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

subroutine launch_svd_and_rotate

use declare

real*8 :: lambda_lambdadagger_trace
real*8 :: lambda_percent

allocate(lambda(nel))
allocate(left_eig(nel,nel))
allocate(right_eig_t(norb-nel,norb-nel))
allocate(right_eig(norb-nel,norb-nel))

call svd(mat_to_svd,nel,norb-nel,lambda,left_eig,right_eig_t)

if (unr) then
 if (countunr .eq. 1) then
  eigname = 'svd_'//trim(orbs_filename)//'A.eig'
  call weigen(eigname,lambda,nel)
 endif
 if (countunr .eq. 2) then 
  eigname = 'svd_'//trim(orbs_filename)//'B.eig'
  call weigen(eigname,lambda,nel)
 endif
else
 eigname = 'svd_'//trim(orbs_filename)//'.eig'
 call weigen(eigname,lambda,nel)
endif

do i=1,nel
lambda_lambdadagger_trace = lambda_lambdadagger_trace + (lambda(i))**2.0d0
enddo

write(50,*) 'Trace of lambda(lambda^dagger)'
write(50,'(f12.5)') lambda_lambdadagger_trace
write(50,*)

write(6,*) 'If the square of a singular value contributes'
write(6,*) 'by more than ten percent to the sum of squared' 
write(6,*) 'singular values, it is reported below'
write(6,*)
write(6,*) 'Singular value, squared singular value, and percentage of its contribution'
write(6,*)

j = 0

do i=1,nel
lambda_percent = (lambda(i)**2.0d0)/lambda_lambdadagger_trace
if (lambda_percent .gt. 0.1d0) then
j = j + 1
write(6,'(2f12.5,f12.1)') lambda(i),lambda(i)**2.0d0,lambda_percent*100.0d0
endif
enddo

if (j .eq. 0) then
 write(6,*) 'There is no singular value superior to 0.1'
 write(50,*) 'There is no singular value superior to 0.1'
endif

! test the SVD

right_eig = transpose(right_eig_t)

allocate(rotated_O_lcao(nbs,nel))
allocate(rotated_V_lcao(nbs,norb-nel))

rotated_O_lcao = matmul(O_lcao,left_eig)
rotated_V_lcao = matmul(V_lcao,right_eig)

allocate(tLK(nel,nbs))
allocate(OdaggerSO(nel,nel))

tLK = transpose(rotated_O_lcao)
tLK = matmul(tLK,S)
OdaggerSO = matmul(tLK,rotated_O_lcao)

call trace_mat(OdaggerSO,'(O^dagger)SO',nel)

deallocate(tLK,OdaggerSO)

allocate(tLK(norb-nel,nbs))
allocate(VdaggerSV(norb-nel,norb-nel))

tLK = transpose(rotated_V_lcao)
tLK = matmul(tLK,S)
VdaggerSV = matmul(tLK,rotated_V_lcao)

call trace_mat(VdaggerSV,'(V^dagger)SV',norb-nel)

deallocate(tLK,VdaggerSV)

deallocate(left_eig,right_eig_t,right_eig)

end
