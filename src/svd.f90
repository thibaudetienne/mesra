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

SUBROUTINE SVD(a,n,m,s,u,vt)

implicit none
integer :: n,m,n1,error,lwork,i,j
integer, allocatable :: iwork(:)
real*8 :: a (n,m)
real*8 :: s(n), U(n,n), VT(m,m)
real*8, allocatable :: work(:)

allocate(work(1))
allocate(iwork(8*n))

S = 0.0D0
U = 0.0D0
VT = 0.0D0

!     First call to DGESDD to set the proper parameters for the SVD      

call DGESDD( 'A', n, m, A, n, S, U, n, VT, m,&
                  WORK, -1, iwork, error )

lwork = work(1)   

deallocate(work)
allocate (work (lwork) )

!     SVD

call DGESDD( 'A', n, m, A, n, S, U, n, VT, m,&
                  WORK, LWORK, iwork, error )

deallocate(work)
deallocate(iwork) 

RETURN
END
