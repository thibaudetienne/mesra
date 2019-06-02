subroutine det_at(central_mat,fchk_out,sqmatdim,U,mvec)

! Actual computation of the detachment/attachment density matrices
! in the canonical space.

use declare

integer :: sqmatdim
real*8 :: central_mat(sqmatdim,sqmatdim)
real*8 :: U(sqmatdim,sqmatdim)
character*(*) :: fchk_out
real*8 :: mvec(sqmatdim)

! Allocates the MO-space detachment (gamma_d) and MO-space attachment (gamma_a) density matrices.

allocate(gamma_d(sqmatdim,sqmatdim))
allocate(gamma_a(sqmatdim,sqmatdim))

! Computes them.

call da(U,mvec,central_mat,gamma_d,gamma_a,sqmatdim)

! If the difference density matrix is not built from an adiabatic connection of the
! Z-vector, then the eigenvalues are stored into an .eig file (weigen = "Write EIGENvalues").

if (adiab) then
 continue
else
 
 eigname = 'da'//trim(fchk_out)//'.eig'

 call weigen(eigname,mvec,sqmatdim)

! The consistency of the trace of the detachment/attachment density matrices is checked.
 
 trD = trmat1
 trA = trmat2

 call trace_mat(gamma_d,'gamma_d',sqmatdim)
 call trace_mat(gamma_a,'gamma_a',sqmatdim)
endif

! The detachment/attachment density matrices are transformed into the AO space.

allocate(gamma_d_ao(nbs,nbs))
allocate(gamma_a_ao(nbs,nbs))

call mo_to_ao(gamma_d,gamma_d_ao)
call mo_to_ao(gamma_a,gamma_a_ao)

if (adiab) then
 continue
else

! They are then converted into vector components, for storing them as the upper triangle of a density matrix.

allocate(triangled(ntr))
allocate(trianglea(ntr))

call triangleK(gamma_d_ao,triangled,nbs)
call triangleK(gamma_a_ao,trianglea,nbs)

! They are written in detachment.fchk and attachment.fchk.

call wfchk_dens(triangled,'detachment'//trim(fchk_out)//'.fchk',ntr)
call wfchk_dens(trianglea,'attachment'//trim(fchk_out)//'.fchk',ntr) 

! Construction of the DS and AS matrices, and verification of their trace.
! (tp stands for "TemPorary".)

allocate(tp(nbs,nbs))

tp = matmul(gamma_d_ao,S)
call trace_mat(tp,'DS',nbs)

tp = matmul(gamma_a_ao,S)
call trace_mat(tp,'AS',nbs)

deallocate(tp)

endif

! D/A population analysis, if requested.

LinearX = xLA
if (LA) call LinearAlgebra(gamma_d_ao,gamma_a_ao,LinearX)

! Deallocation of the arrays.

deallocate(gamma_d,gamma_a,gamma_d_ao,gamma_a_ao)

if (adiab) then
 continue
else
 deallocate(triangled,trianglea)
endif

end

