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

subroutine PopulationAnalysis(matrix_d,matrix_a,PA_x)

! D/A population analysis

use declare

real*8 :: PA_x,PA_y
real*8 :: matrix_d(nbs,nbs), matrix_a(nbs,nbs)

! If the "scanpa" option has been selected in the input,
! scans the x and y exponents in S^x P S^y (with P = D,A)

 if (scanPA) then

! Defines the x and y exponents as decimal numbers

  do iteration=0,100
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
   PA_x = iteration*0.01d0
   PA_y = 1.0d0 - PA_x 

! Calls the population analysis subroutine

   call PA_analysis(matrix_d,matrix_a,PA_x,PA_y)

  enddo

 else

! If the "scanpa" option has not been selected, a single-point D/A
! population analysis is performed with the x and y exponents given
! in the input (actually, only x is given, since y = 1 - x)

  PA_x = PA_x*0.01d0
  PA_y = 1.0d0 - PA_x

! Calls the population analysis subroutine

  call PA_analysis(matrix_d,matrix_a,PA_x,PA_y)

 endif

end
