subroutine LinearAlgebra(matrix_d,matrix_a,LA_x)

! D/A population analysis

use declare

real*8 :: LA_x,LA_y
real*8 :: matrix_d(nbs,nbs), matrix_a(nbs,nbs)

! If the "scanpa" option has been selected in the input,
! scans the x and y exponents in S^x P S^y (with P = D,A)

 if (scanLA) then

! Defines the x and y exponents as decimal numbers

  do iteration=0,100
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
   LA_x = iteration*0.01d0
   LA_y = 1.0d0 - LA_x 

! Calls the population analysis subroutine

   call LA_analysis(matrix_d,matrix_a,LA_x,LA_y)

  enddo

 else

! If the "scanpa" option has not been selected, a single-point D/A
! population analysis is performed with the x and y exponents given
! in the input (actually, only x is given, since y = 1 - x)

  LA_x = LA_x*0.01d0
  LA_y = 1.0d0 - LA_x

! Calls the population analysis subroutine

  call LA_analysis(matrix_d,matrix_a,LA_x,LA_y)

 endif

end
