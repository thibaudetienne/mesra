subroutine double_mat_prod(mat1,mat2,mat3,mat4,dimmat)

! Performs a double matrix product with one subroutine.

integer :: dimmat
real*8 :: mat1(dimmat,dimmat),mat2(dimmat,dimmat),mat3(dimmat,dimmat)
real*8 :: mat4(dimmat,dimmat),temp_mat(dimmat,dimmat)

temp_mat = matmul(mat1,mat2)
mat4 = matmul(temp_mat,mat3)


end
