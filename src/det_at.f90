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

PopX = xPA
if (PA) call PopulationAnalysis(gamma_d_ao,gamma_a_ao,PopX)

! Deallocation of the arrays.

deallocate(gamma_d,gamma_a,gamma_d_ao,gamma_a_ao)

if (adiab) then
 continue
else
 deallocate(triangled,trianglea)
endif

end

