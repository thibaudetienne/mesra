subroutine LA_analysis(matrix_dLA,matrix_aLA,x_PA,y_PA)

! Performs the actual D/A population analysis.

use declare

real*8 :: x_PA,y_PA
real*8 :: matrix_dLA(nbs,nbs),matrix_aLA(nbs,nbs)

! Allocates the (S^x)P(S^y) arrays (P = D,A) 

allocate(SxDSy(nbs,nbs))
allocate(SxASy(nbs,nbs))

! Creates S^x, S^y and (S^x)P(S^y) (P = D,A).

call LA_mat(matrix_dLA,matrix_aLA,SxDSy,SxASy,x_PA,y_PA)

! Computes the descriptors from (S^x)P(S^y) (P = D,A).

call QuantumMetricsLA(SxDSy,SxASy)

! Deallocates the (S^x)P(S^y) matrices (P = D,A).

deallocate(SxDSy,SxASy)

end
