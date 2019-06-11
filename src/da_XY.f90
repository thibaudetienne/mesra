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

subroutine da_XY

use declare

! Extracts X and Y.

call common_XY

! Since tK is allocated when read by tK in daXY, we should deallocate it between
! the calculations related to the alpha and beta electrons.

if (unr) then
if (countunr .eq. 2) deallocate(tK,t)
else
deallocate(tK,t)
endif

! Constructs the difference density matrix, and the detachment/attachment density matrices.

allocate(xy_XXt(nel,nel))
allocate(xy_XtX(norb-nel,norb-nel))
allocate(xy_YYt(nel,nel))
allocate(xy_YtY(norb-nel,norb-nel))

xy_XXt = matmul(xy_X,xy_Xt)
xy_XtX = matmul(xy_Xt,xy_X)

xy_YYt = matmul(xy_Y,xy_Yt)
xy_YtY = matmul(xy_Yt,xy_Y)

allocate(gamma_d_XY(norb,norb))
allocate(gamma_a_XY(norb,norb))

gamma_d_XY = 0.0d0
gamma_a_XY = 0.0d0

allocate(diff_mat_unrelaxed(norb,norb))

diff_mat_unrelaxed = 0.0d0

do i=1,nel
 do j=1,nel
diff_mat_unrelaxed(i,j) = -1.0d0*(xy_XXt(i,j) + xy_YYt(i,j))
gamma_d_XY(i,j) = -1.0d0*diff_mat_unrelaxed(i,j)
 enddo
enddo

call trace_mat(gamma_d_XY,'gamma_d_XY',norb)

do i=nel+1,norb
 do j=nel+1,norb
diff_mat_unrelaxed(i,j) = (xy_XtX(i-nel,j-nel) + xy_YtY(i-nel,j-nel))
gamma_a_XY(i,j) = diff_mat_unrelaxed(i,j)
 enddo
enddo

call trace_mat(gamma_a_XY,'gamma_a_XY',norb)

call trace_mat(diff_mat_unrelaxed,'diff_mat_unrelaxed',norb)

! Prints the matrices in fchk files.

if (unr) then
 if (countunr .eq. 1) then
call print_mat_mo_to_ao_fchk(diff_mat_unrelaxed,'DeltaUAlpha_XY.fchk')
call print_mat_mo_to_ao_fchk(gamma_d_XY,'detachmentAlpha_XY.fchk')
call print_mat_mo_to_ao_fchk(gamma_a_XY,'attachmentAlpha_XY.fchk')
 else if (countunr .eq. 2) then
call print_mat_mo_to_ao_fchk(diff_mat_unrelaxed,'DeltaUBeta_XY.fchk')
call print_mat_mo_to_ao_fchk(gamma_d_XY,'detachmentBeta_XY.fchk')
call print_mat_mo_to_ao_fchk(gamma_a_XY,'attachmentBeta_XY.fchk')
 endif
else
call print_mat_mo_to_ao_fchk(diff_mat_unrelaxed,'DeltaU_XY.fchk')
call print_mat_mo_to_ao_fchk(gamma_d_XY,'detachment_XY.fchk')
call print_mat_mo_to_ao_fchk(gamma_a_XY,'attachment_XY.fchk')
endif

! Computes the descriptors from population analysis (PA = Population Analysis), if requested.

if (PA) then

 allocate(gamma_d_XY_ao(nbs,nbs))
 allocate(gamma_a_XY_ao(nbs,nbs))

 if (countunr .eq. 2) then
  C = Cb
  Cdagger = transpose(Cb)
 endif

 call mo_to_ao(gamma_d_XY,gamma_d_XY_ao)
 call mo_to_ao(gamma_a_XY,gamma_a_XY_ao)

! xPA is the x exponent in the generalized population analysis.

 PopX = xPA
 
 call PopulationAnalysis(gamma_d_XY_ao,gamma_a_XY_ao,PopX)

 deallocate(gamma_d_XY_ao,gamma_a_XY_ao)

endif

deallocate(xy_XXt,xy_XtX,xy_YYt,xy_YtY)
deallocate(diff_mat_unrelaxed)
deallocate(gamma_d_XY,gamma_a_XY)
deallocate(xy_X,xy_Y,xy_Xt,xy_Yt)

end
