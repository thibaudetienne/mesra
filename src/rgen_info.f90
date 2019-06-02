subroutine rgen_info

! Reads the general information.

use declare

read(10,*) soft 	! Software
read(10,*) fchk		! .fchk filename
read(10,*) fov		! overlap matrix file
read(10,*) fdens	! density matrix file
read(10,*) nt		! total number of transitions computed
read(10,*) ns		! the number of the transition of interest

end
