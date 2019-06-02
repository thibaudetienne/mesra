subroutine LA_mat(d_matao,a_matao,xy_d,xy_a,x_LA,y_LA)

! Computes S^x and S^y, creates the actual (S^x)P(S^y) (P = D,A) matrices, and computes their trace.

use declare

real*8 :: x_LA,y_LA
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
 lSx(i) = (lS(i))**x_LA
 lSy(i) = (lS(i))**y_LA
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

call trace_mat(xy_d,'SxDSy',nbs)
call trace_mat(xy_a,'SxASy',nbs)

x = 0.0d0
y = 0.0d0

do i=1,nbs
 x = x + xy_d(i,i)
 y = y + xy_a(i,i)
enddo

write(50,*) 'Tr SxDSy'
write(50,'(f12.5)') x
write(50,*) 'Tr SxASy'
write(50,'(f12.5)')  y

endif

end
