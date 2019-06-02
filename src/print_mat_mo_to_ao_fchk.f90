subroutine print_mat_mo_to_ao_fchk(mat_to_print,mat_name)

! Converts a square matrix from MO to AO space, re-write it as the upper triangle of a square matrix,
! and writes the result in an fchk file. 

use declare

character*(*) :: mat_name
real*8 :: mat_to_print(norb,norb)
real*8 :: mat_to_print_ao(nbs,nbs)
real*8 :: triangle_mat_to_print_ao(ntr)

! Brings the matrix in the AO space.

call mo_to_ao(mat_to_print,mat_to_print_ao)

! Brings it in the triangle shape (symmetric matrices).

call triangleK(mat_to_print_ao,triangle_mat_to_print_ao,nbs)

! Writes the result in an .fchk file (wfchk = "Write in an .FCHK file")

call wfchk_dens(triangle_mat_to_print_ao,mat_name,ntr)

end
