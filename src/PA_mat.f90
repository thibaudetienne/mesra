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

subroutine PA_mat(d_matao,a_matao,xy_d,xy_a,x_PA,y_PA)

! Computes S^x and S^y, creates the actual (S^x)P(S^y) (P = D,A) matrices, and computes their trace.

use declare

real*8 :: x_PA,y_PA
real*8 :: d_matao(nbs,nbs),a_matao(nbs,nbs),xy_d(nbs,nbs),xy_a(nbs,nbs)
real*8 :: uS(nbs,nbs),lS(nbs),Sx(nbs,nbs),Sy(nbs,nbs)
real*8 :: uSt(nbs,nbs),lSx(nbs),lSy(nbs),diag_Sx(nbs,nbs),diag_Sy(nbs,nbs)

! Initializes all the quantities and arrays for this subroutine.

uS = 0.0d0
lS = 0.0d0
uSt = 0.0d0
lSx = 0.0d0     
lSy = 0.0d0
diag_Sx = 0.0d0
diag_Sy = 0.0d0
Sx = 0.0d0
Sy = 0.0d0
xy_d = 0.0d0
xy_a = 0.0d0

! Diagonalizes S to create S^x and S^y.

call diag(S,nbs,lS,uS)

uSt = transpose(uS)

do i=1,nbs
 if (lS(i) .lt. 0.0d0) lS(i) = 0.0d0
 lSx(i) = (lS(i))**x_PA
 lSy(i) = (lS(i))**y_PA
 diag_Sx(i,i) = lSx(i)
 diag_Sy(i,i) = lSy(i)
enddo

call double_mat_prod(uS,diag_Sx,uSt,Sx,nbs)
call double_mat_prod(uS,diag_Sy,uSt,Sy,nbs)

! Creates (S^x)P(S^y) (P = D,A).

call double_mat_prod(Sx,d_matao,Sy,xy_d,nbs)
call double_mat_prod(Sx,a_matao,Sy,xy_a,nbs)

if (adiab) then
 continue
else

! Computes their trace.
! (x and y here have nothing to do with the x and y exponents used before. Here, they are simple real numbers.)

call trace_mat(xy_d,'(S^x)DS^y',nbs)
call trace_mat(xy_a,'(S^x)AS^y',nbs)

endif

end
