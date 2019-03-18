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
