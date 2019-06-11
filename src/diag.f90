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

subroutine diag(a,n,wr,vr)

! Diagonalizes a square matrix a to produce its eigenvectors vr and eigenvalues wr.

implicit none
integer :: n,q,i,error,lwork,j,IL,IU
real*8 :: a(n,n), ABSTOL, DLAMCH,VL,VU
real*8, allocatable :: lap_wrk(:)  
real*8 :: wr(n), vr(n,n), tvr(n,n), vtvr(n,n)
integer, allocatable :: jfail(:), iwork(:)
real*8 :: tp(n,n)

! Initial set of parameters, and primary allocation of arrays.

ABSTOL = 2*DLAMCH('S')
allocate(lap_wrk(1))
allocate(jfail (n))
allocate(iwork (5*n))

tp = a
lwork = -1

! Performs a first call to set optimal space before the actual use of DSYEV.

call DSYEVX('V','A','U',n,tp,n,VL,VU,IL,IU,ABSTOL,q,&
 wr, vr, n,lap_wrk,lwork,iwork,jfail,error) 
 lwork = lap_wrk(1)

deallocate(lap_wrk)
allocate(lap_wrk(lwork))

! Actual use of DSYEV to diagonalize a.

call DSYEVX('V','A','U',n,tp,n,VL,VU,IL,IU,ABSTOL,q,&
 wr, vr, n,lap_wrk,lwork,iwork,jfail,error) 

! Deallocates the arrays.

deallocate(lap_wrk)
deallocate(iwork)
deallocate(jfail)

RETURN
END
