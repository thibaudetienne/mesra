subroutine da(Umatrix,mvector,matrix0,matrix1,matrix2,matrix_dimension)

! Diagonalizes a matrix, separates its eigenvalues according to their sign,
! backtransforms the result into detachment/attachment matrices.

use declare

integer :: matrix_dimension
real*8 :: matrix0(matrix_dimension,matrix_dimension)
real*8 :: matrix1(matrix_dimension,matrix_dimension)
real*8 :: matrix2(matrix_dimension,matrix_dimension)
real*8 :: Umatrix(matrix_dimension,matrix_dimension)
real*8 :: mvector(matrix_dimension)
real*8 :: Mm(matrix_dimension,matrix_dimension)
real*8 :: m(matrix_dimension)
real*8 :: kp(matrix_dimension,matrix_dimension)
real*8 :: km(matrix_dimension,matrix_dimension)
real*8 :: Mmdag(matrix_dimension,matrix_dimension)

! Initializes the arrays.

Mm = 0.0d0
m = 0.0d0
kp = 0.0d0
km = 0.0d0
matrix1 = 0.0d0
matrix2 = 0.0d0
Mmtranspose = 0.0d0

! Diagonalizes the matrix.

call diag(matrix0,matrix_dimension,m,Mm)

! Splits the eigenvalues according to their sign.

do i=1,matrix_dimension
 if (m(i) .lt. 0.0d0) then
  km(i,i) = -1.0d0*m(i)
 else
  kp(i,i) = m(i)
 endif
enddo

! Backtransforms the resulting matrices.

Mmdag = transpose(Mm)

call double_mat_prod(Mm,km,Mmdag,matrix1,norb)
call double_mat_prod(Mm,kp,Mmdag,matrix2,norb)

! Computes the trace of the resulting two matrices.

x = 0.0d0
y = 0.0d0

do i=1,matrix_dimension
 x = x + matrix1(i,i)
 y = y + matrix2(i,i)
enddo

trmat1 = x
trmat2 = y

! (Re-initializes x and y.)

x = 0.0d0
y = 0.0d0

if (adiab) then
continue
else

! Prints the largest eigenvalues.

write(6,*)
write(6,*) 'Largest eigenvalues (> 0.1)'
write(6,*)
write(50,*)
write(50,*) 'Largest eigenvalues (> 0.1)'
write(50,*)

j = 0

do i=1,matrix_dimension
 if (m(i) .gt. 0.1d0)  then
  j = j + 1
  write(6,'(A6,f9.5)') 'eig+ ', m(i)
  write(50,'(A6,f9.5)') 'eig+ ', m(i)
 endif
 if (m(i) .lt. -0.1d0) then
  j = j + 1
  write(6,'(A6,f9.5)') 'eig- ',m(i)
  write(50,'(A6,f9.5)') 'eig- ',m(i)
 endif
enddo
  write(50,*)

if (j .eq. 0) then
 write(6,*) 'No eigenvalue superior to 0.1'
 write(50,*) 'No eigenvalue superior to 0.1'
endif

endif

! Stores the eigenvectors and eigenvalues.

Umatrix = Mm
mvector = m

end
