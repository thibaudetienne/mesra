subroutine deal_dens

use declare

if (jobtype .eq. 'orbsXY' .or. jobtype .eq. 'pNTOs' .or. jobtype .eq. 'CTOs' .or. &
jobtype .eq. 'aNTOs') then
 continue
else if (jobtype .eq. 'daXY') then
 continue
else if (jobtype .eq. 'rlxy_LA') then
 continue
else
 deallocate(p,px)
 deallocate(pK,pxK)
 deallocate(pKS,pxKS)
 if (jobtype .eq. 'dar') then
  if (unr) then
   continue
  else
  deallocate(pxKrelaxed,pxKrelaxedS,pxrelaxed,zvec,zzd_zdz)
  endif
 endif
endif

if (unr) then
 if (jobtype .eq. 'dau' .or. jobtype .eq. 'daz' .or. jobtype .eq. 'dar') then
  deallocate(palpha,pbeta,pKalpha,pKalphaS,pKbeta,pKbetaS,pxalpha,pxbeta,pxKalpha,pxKalphaS,pxKbeta,pxKbetaS)
 endif
if (relax) then
 deallocate(pxKspinRelaxed,pxalpharelaxed,pxbetarelaxed)
endif
endif

end
