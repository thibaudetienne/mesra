subroutine relax_D

use declare

call rPrelax_g

allocate(pxrelaxed(norb,norb))

call ao_to_mo(pKrelaxed,pxrelaxed)

call trace_mat(pxrelaxed,'pxrelaxed',norb)

allocate(zvec(norb,norb))

zvec = pxrelaxed - px

allocate(zzd_zdz(norb,norb))

if (adiab) then

call adiab_z

else

call no_adiab_z

endif

end

