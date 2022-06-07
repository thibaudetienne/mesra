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

subroutine rdaggerf

! Reads the name of the files containing the unrelaxed detachment/attachment density matrices
! and the detachment/attachment density matrix constructed by diagonalizing the Z-vector.

use declare

! shell_statement = 0 (1) means (un)restricted calculation.

read(10,*) shell_statement

if (shell_statement .eq. 0) then
read(10,*) detachmentU
read(10,*) attachmentU
read(10,*) detachmentR
read(10,*) attachmentR
else if (shell_statement .eq. 1) then
read(10,*) detachmentUalpha
read(10,*) attachmentUalpha
read(10,*) detachmentRalpha
read(10,*) attachmentRalpha
read(10,*) detachmentUbeta
read(10,*) attachmentUbeta
read(10,*) detachmentRbeta
read(10,*) attachmentRbeta
endif

end
