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

subroutine rPA_status

! Reads the status to be decided relatively to the population analysis to perform.

use declare

! Reads the status from the input.

read(10,*) PA_status

! PA = population analysis => which exponents (x)?
! noPA = no population analysis;
! scan = a scan of the x exponents from 0 to 1.

 if (PA_status .eq. 'PA') then
  PA = .true.
  scanPA = .false.
  read(10,*) xPA
 else if (PA_status .eq. 'noPA') then
  PA = .false.
  scanPA = .false.
 else if (PA_status .eq. 'scan') then
  PA = .true.
  scanPA = .true.
 endif

end
