subroutine svd_and_rotation_XY

! Constructs the appropriate transition matrices, performs their singular value decomposition,
! and rotates the atomic-orbital space to produce transition orbitals.

use declare

! Calls a subroutines that read X and Y to create the appropriate transition matrices.

call common_XY
call trans_orbs_XY

! Prepares the rotation.

allocate(O_lcao(nbs,nel))
allocate(V_lcao(nbs,norb-nel))
allocate(mat_to_svd(nel,norb-nel))

do i=1,nel
 do j=1,nbs
O_lcao(j,i) = lcao_coeff_mat(j,i)
 enddo
enddo

do i=1,norb-nel
 do j=1,nbs
V_lcao(j,i) = lcao_coeff_mat(j,i+nel)
 enddo
enddo

! Starts with the pNTOs.

if (jobtype .eq. 'pNTOs' .or. jobtype .eq. 'orbsXY') then
 write(6,*)
 write(6,*) '# pNTOs computation'
 write(6,*)
 if (0.5d0*(xy_norm-1.0d0) .ge. 0.0d0) then
  write(6,*) 'Amplitude of the particle-hole correlation'
  write(6,'(f12.5)') 0.5d0*(xy_norm-1.0d0)
 endif

! mat_to_svd is a key array that allows to use the launch_svd_and_rotate subroutine
! for any matrix to svd. orbs_filename is used for creating the fchk files.

 mat_to_svd = t_1
 orbs_filename = 'pNTOs'

! Launches the SVD and rotation of orbitals.

 call launch_svd_and_rotate

! Checks the orthonormality conservation.

 call orth(rotated_O_lcao,'pNTOsO',nbs,nel)
 call orth(rotated_V_lcao,'pNTOsV',nbs,norb-nel)

! Uses the appropriate subroutine for writing the fchk files. For closed-shell molecules,
! the orbitals to be output are usually the alpha orbitals (arbitrarily).

 if (unr) then
  if (countunr .eq. 1) then
   call wfchk_orbsAlpha(lambda,rotated_O_lcao,rotated_V_lcao,'pNTOsA.fchk')
  else if (countunr .eq. 2) then
   call wfchk_orbsBeta(lambda,rotated_O_lcao,rotated_V_lcao,'pNTOsB.fchk')
  endif
 else
  call wfchk_orbsAlpha(lambda,rotated_O_lcao,rotated_V_lcao,'pNTOs.fchk')
 endif

 deallocate(lambda,rotated_O_lcao,rotated_V_lcao)
endif

! The procedure is similar for CTOs.

if (jobtype .eq. 'CTOs' .or. jobtype .eq. 'orbsXY') then
 write(6,*)
 write(6,*) '# CTOs computation'
 write(6,*)

 mat_to_svd = t_2
 orbs_filename = 'CTOs'

 call launch_svd_and_rotate
 call orth(rotated_O_lcao,'CTOsO',nbs,nel)
 call orth(rotated_V_lcao,'CTOsV',nbs,norb-nel)

 if (unr) then
  if (countunr .eq. 1) then
   call wfchk_orbsAlpha(lambda,rotated_O_lcao,rotated_V_lcao,'CTOsA.fchk')
  else if (countunr .eq. 2) then
   call wfchk_orbsBeta(lambda,rotated_O_lcao,rotated_V_lcao,'CTOsB.fchk')
  endif
 else
  call wfchk_orbsAlpha(lambda,rotated_O_lcao,rotated_V_lcao,'CTOs.fchk')
 endif
 
 deallocate(lambda,rotated_O_lcao,rotated_V_lcao)
endif

! The procedure is similar for aNTOs.

if (jobtype .eq. 'aNTOs' .or. jobtype .eq. 'orbsXY') then
 write(6,*)
 write(6,*) '# aNTOs computation'
 write(6,*)

 mat_to_svd = t_3
 orbs_filename = 'aNTOs'

 call launch_svd_and_rotate
 call orth(rotated_O_lcao,'aNTOsO',nbs,nel)
 call orth(rotated_V_lcao,'aNTOsV',nbs,norb-nel)

 if (unr) then
  if (countunr .eq. 1) then
   call wfchk_orbsAlpha(lambda,rotated_O_lcao,rotated_V_lcao,'aNTOsA.fchk')
  else if (countunr .eq. 2) then
   call wfchk_orbsBeta(lambda,rotated_O_lcao,rotated_V_lcao,'aNTOsB.fchk')
  endif
 else
  call wfchk_orbsAlpha(lambda,rotated_O_lcao,rotated_V_lcao,'aNTOs.fchk')
 endif

 deallocate(lambda,rotated_O_lcao,rotated_V_lcao)
endif

deallocate(mat_to_svd,O_lcao,V_lcao)

! Manages the outputs.

if (countunr .eq. 1) then
continue
else
write(6,*)
endif

if (unr) then
continue
else
write(6,*)
endif

end
