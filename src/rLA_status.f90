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

subroutine rLA_status

! Reads the status to be decided relatively to the population analysis to perform.

use declare

! Reads the status from the input.

read(10,*) LA_status

! LA = population analysis (Linear Algebra) => which exponents (x)?
! noLA = no population analysis;
! scan = a scan of the x exponents from 0 to 1.

 if (LA_status .eq. 'LA') then
  LA = .true.
  scanLA = .false.
  read(10,*) xLA
 else if (LA_status .eq. 'noLA') then
  LA = .false.
  scanLA = .false.
 else if (LA_status .eq. 'scan') then
  LA = .true.
  scanLA = .true.
 endif

end
