module declare

implicit none

! General

integer :: error,i,j,k,nLA,nAdiab,jobN,countunr,countrlx,shell_statement,iteration
real*8 :: x,y,z,xi,xLA,yLA,trmat1,trmat2,LinearX
character*128 :: LA_status,eigname

! Input reading

character*128 :: jobname,jobtype,qmLA_obj,subjobtype
character*128 :: soft,fchk,fov,fdens,input,fdensx
integer :: nt,ns
logical :: relax,adiab,LA,LAr,xy,scanLA,qmLA,t_LA,t_scanLA
character*128,allocatable :: dummy_array(:)

! Output
character*128 :: logfile

! Data related to the calculation

character*128 :: line
integer :: nea,neb,nel,nbs,norb,ntr
integer :: n,l,nL
logical :: unr

! LCAO matrices

real*8,allocatable :: triangle(:),trianglespin(:),lvec(:),triangled(:),trianglea(:)
real*8,allocatable :: C(:,:),Ca(:,:),Cb(:,:),Cdag(:,:),Cadag(:,:),Cbdag(:,:)
character*128 :: marker_Ca, marker_Cb, marker_Ea, marker_Eb

! State density matrices

real*8,allocatable :: p(:,:),pK(:,:),pKS(:,:)
real*8,allocatable :: px(:,:),pxK(:,:),pxKS(:,:)
integer :: n0
real*8,allocatable :: pKspin(:,:),pKalpha(:,:),pKbeta(:,:)
real*8,allocatable :: pKalphaS(:,:),pKbetaS(:,:)
real*8,allocatable :: palpha(:,:),pbeta(:,:)
real*8,allocatable :: pxalpha(:,:),pxKalpha(:,:),pxbeta(:,:),pxKbeta(:,:)
real*8,allocatable :: pxKspin(:,:),pxKalphaS(:,:),pxKbetaS(:,:)
real*8,allocatable :: pxKspinRelaxed(:,:),pxKalphaRelaxed(:,:),pxKbetaRelaxed(:,:),pxalphaRelaxed(:,:),pxbetaRelaxed(:,:)
real*8,allocatable :: pxKalphaRelaxedS(:,:),pxKbetaRelaxedS(:,:)

! Overlap matrix

real*8,allocatable :: S(:,:)

! Temporary matrices

real*8,allocatable :: tp(:,:),tLK(:,:),tKL(:,:),tLL(:,:),tKK(:,:)
real*8,allocatable :: tKK1(:,:),tKK2(:,:)

! Difference density matrix, detachment/attachment

real*8,allocatable :: gD(:,:),gamma_d(:,:),gamma_a(:,:)
real*8,allocatable :: gamma_d1(:,:),gamma_a1(:,:),gamma_d_ao(:,:),gamma_a_ao(:,:)
real*8,allocatable :: triangle_gamma_d_ao(:),triangle_gamma_a_ao(:)
real*8,allocatable :: Mmg(:,:),Mmd(:,:)
real*8,allocatable :: SxDSy(:,:),SxASy(:,:),Uvec(:,:)
real*8 :: trD,trA,trDZ,trAZ,phiSLA,chimLA,chipLA,phiLA,psiLA,theta,theta0,thetaZ,alter_thetaZ
real*8 :: pi,phiSLA0,phiLA0
real*8 :: phiSdag,phidag,psidag
real*8,allocatable :: fd(:),fa(:),fdta(:),famd(:),famdp(:),famdm(:)
real*8,allocatable :: gamma_d_XY_ao(:,:),gamma_a_XY_ao(:,:)

! Relaxed difference density matrix

real*8,allocatable :: pxKrelaxed(:,:),pxKrelaxedS(:,:),pxrelaxed(:,:)
real*8,allocatable :: zvec(:,:),zzd_zdz(:,:),zvecxi(:,:)
real*8,allocatable :: U0(:,:),U0t(:,:),U0tU(:,:)

! Orbitals

real*8,allocatable :: tK(:,:),t(:,:),TS(:,:),diff_mat_unrelaxed(:,:)
real*8,allocatable :: tKb(:,:),tb(:,:),TSb(:,:)
real*8,allocatable :: ttrsp(:,:),ttdag(:,:)
real*8,allocatable :: xy_X(:,:),xy_Y(:,:),xy_Xt(:,:),xy_Yt(:,:),XYt(:,:),YXt(:,:)
real*8,allocatable :: xy_XXt(:,:),xy_XtX(:,:),xy_YYt(:,:),xy_YtY(:,:)
real*8,allocatable :: t_0p(:,:),t_0m(:,:),t_02m(:,:),t_02p(:,:),t_0pt(:,:),t_0pt_0pt(:,:)
real*8,allocatable :: t_02mt(:,:),t_02mt_02mt(:,:)
real*8,allocatable :: t_1(:,:),t_2(:,:),t_3(:,:)
real*8 :: xy_norm,xy_residue,t3_norm,x2y2_norm
real*8,allocatable :: gamma_d_XY(:,:),gamma_a_XY(:,:)

character*128 :: orbs_filename

! SVD

real*8,allocatable :: lcao_coeff_mat(:,:)
real*8,allocatable :: lambda(:),left_eig(:,:),right_eig_t(:,:),mat_to_svd(:,:),right_eig(:,:)
real*8,allocatable :: O_lcao(:,:),V_lcao(:,:),rotated_O_lcao(:,:),rotated_V_lcao(:,:)
real*8,allocatable :: OdagSO(:,:),VdagSV(:,:)

! AO matrices

real*8,allocatable :: gD_ao(:,:),T_ao(:,:),zvec_ao(:,:)
real*8,allocatable :: triangle_gD_ao(:),triangle_T_ao(:),triangle_zvec_ao(:)

! Relaxed metrics (rlxy_LA)


real*8 :: thetaUalpha, thetaUbeta
real*8 :: phiSLAUalpha, phiSLAUbeta
real*8 :: phiLAUalpha, phiLAUbeta
real*8 :: thetaZalpha, thetaZbeta

real*8 :: phiSLAU,phiLAU,psiLAU,thetaU,phiSLAdag,phiLAdag,psiLAdag
real*8 :: phiSLAdagalpha, phiSLAdagbeta
real*8 :: phiLAdagalpha, phiLAdagbeta
real*8 :: psiLAdagalpha, psiLAdagbeta

! Numerical Integration

character*128 :: file1,file2
character*128 :: input1,input2,input3,filename
integer :: natom,n1,n2,n3,i1,i2,i3
integer :: natom1,n11,n21,n31,nci
real*8 :: dx,dy,dz,int1,integ1,integ2,phis,dummy,dummy1
real*8 :: integp, integm, integrp, integrm
real*8 ::dx1,dy1,dz1,THR,occ,virt
real*8 :: integ1r,integ2r,x1c,y1c,z1c,x2c,y2c,z2c,r12ct
real*8 :: xrp, yrp, zrp, xrm, yrm, zrm, rct
real*8, allocatable :: xNI(:,:,:), yNI(:,:,:), zNI(:,:,:)
real*8, allocatable :: nto1(:,:,:), nto2(:,:,:)
real*8 :: x0,y0,z0,x01,y01,z01
real*8, allocatable :: p1(:,:,:),p2(:,:,:),r(:,:,:)
real*8, allocatable :: pp(:,:,:),pm(:,:,:),diff(:,:,:)

! Split and cubeop

character*128 :: cubefile,cubefile1,cubefile2,cubefile3,op1,op2

! dag

character*128 :: detachmentU,attachmentU,detachmentR,attachmentR
character*128 :: detachmentUalpha,attachmentUalpha,detachmentRalpha,attachmentRalpha
character*128 :: detachmentUbeta,attachmentUbeta,detachmentRbeta,attachmentRbeta

! qmNIr

real*8 :: theta_unrelaxed,phiS_unrelaxed,phi_unrelaxed,lambda_dag,theta_z
real*8 :: phiS_relaxed,phi_relaxed,psi_relaxed,eta
end
