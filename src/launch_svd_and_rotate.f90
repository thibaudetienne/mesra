subroutine launch_svd_and_rotate

use declare

real*8 :: lambda_lambdadag_trace
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
lambda_lambdadag_trace = lambda_lambdadag_trace + (lambda(i))**2.0d0
enddo

write(50,*) 'Trace of lambda(lambda^dag)'
write(50,'(f12.5)') lambda_lambdadag_trace
write(50,*)

write(6,*) 'If the square of a singular value contributes'
write(6,*) 'by more than ten percent to the sum of squared' 
write(6,*) 'singular values, it is reported below'
write(6,*)
write(6,*) 'Singular value, squared singular value, and percentage of its contribution'
write(6,*)

j = 0

do i=1,nel
lambda_percent = (lambda(i)**2.0d0)/lambda_lambdadag_trace
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
allocate(OdagSO(nel,nel))

tLK = transpose(rotated_O_lcao)
tLK = matmul(tLK,S)
OdagSO = matmul(tLK,rotated_O_lcao)

call trace_mat(OdagSO,'(O^dag)SO',nel)

deallocate(tLK,OdagSO)

allocate(tLK(norb-nel,nbs))
allocate(VdagSV(norb-nel,norb-nel))

tLK = transpose(rotated_V_lcao)
tLK = matmul(tLK,S)
VdagSV = matmul(tLK,rotated_V_lcao)

call trace_mat(VdagSV,'(V^dag)SV',norb-nel)

deallocate(tLK,VdagSV)

deallocate(left_eig,right_eig_t,right_eig)

end
