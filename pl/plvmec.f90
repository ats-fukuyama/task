! plvmec.f90

MODULE plvmec
  USE bpsd

  TYPE(bpsd_device_type),PRIVATE,SAVE  :: device
  TYPE(bpsd_equ1D_type),PRIVATE,SAVE   :: equ1D
  TYPE(bpsd_metric1D_type),PRIVATE,SAVE:: metric1D

!... data from wout_file 

  CHARACTER(LEN=4)         :: version        ! version of vmec2000
                                             !   accept '6.90' or '8.40'
  INTEGER                  :: ns             ! number of radial grids
  INTEGER                  :: nsp1           ! = ns + 1
  INTEGER                  :: nsm1           ! = ns - 1
  INTEGER                  :: nfp            ! number of toroidal field periods
  INTEGER                  :: mpol           ! number of poloidal modes ;
                                             !   0 <= m <= mpol-1
  INTEGER                  :: ntor           ! number of toroidal modes ;
                                             !   -ntor <= n <= ntor
  INTEGER                  :: mnmax          ! number of Fourier modes
                                             !   for R, Z, and lambda
  INTEGER                  :: mnmax2         ! number of Fourier modes
                                             !   for B and gsqrt
                                             !   (=mnmax for version 6.9)
  INTEGER                  :: isgn           ! sign of  gsqrth in vmec2000
                                             !   (left-hand coord. -> isgn=-1)
  REAL(rkind), ALLOCATABLE :: xm(:)          ! poloidal mode number
                                             !   for R, Z and lambda; j=1:mnmax
  REAL(rkind), ALLOCATABLE :: xn(:)          ! toroidal mode number
                                             !   for R, Z and lambda; j=1:mnmax
  REAL(rkind), ALLOCATABLE :: xm2(:)         ! poloidal mode number
                                             !   for B and gsqrt;    j=1:mnmax2
  REAL(rkind), ALLOCATABLE :: xn2(:)         ! toroidal mode number
                                             !   for B and gsqrt;    j=1:mnmax2
  REAL(rkind), ALLOCATABLE :: rmnc(:,:)      ! Rmn on full grids;
                                             !   cos,  (1:mnmax,0:ns)
  REAL(rkind), ALLOCATABLE :: zmns(:,:)      ! Zmn on full grids;
                                             !   sin,  (1:mnmax,0:ns)
  REAL(rkind), ALLOCATABLE :: lmnsh(:,:)     ! lambda_mn on half grids;
                                             !   sin,  (1:mnmax,0:ns)
  REAL(rkind), ALLOCATABLE :: bmnh(:,:)      ! Bmn on half grids;
                                             !   cos,  (1:mnmax2,0:ns)
  REAL(rkind), ALLOCATABLE :: gmnh(:,:)      ! gsqrt_mn on half grids;
                                             ! cos, (1:mnmax2,0:ns)
  REAL(rkind), ALLOCATABLE :: bsubumnh(:,:)  ! bsubhalf mesh
  REAL(rkind), ALLOCATABLE :: bsubvmnh(:,:)  ! half mesh
  REAL(rkind), ALLOCATABLE :: bsubsmn(:,:)   ! full mesh
  REAL(rkind), ALLOCATABLE :: bsupumnh(:,:)  ! half mesh
  REAL(rkind), ALLOCATABLE :: bsupvmnh(:,:)  ! half mesh
  REAL(rkind), ALLOCATABLE :: currvmn(:,:)   ! full mesh
  REAL(rkind), ALLOCATABLE :: aiotah(:)      ! rotational transform
                                             !   on half grids ; j=0,ns+1
  REAL(rkind), ALLOCATABLE :: presh(:)       ! pressure on half grids [Pa];
                                             !   j=0,ns+1
  REAL(rkind), ALLOCATABLE :: Itorh(:)       ! twopi*isgn*Itorh/mu0
                                             !   = net toroidal current [A]
                                             !   on half grids ; j=0,ns+1
  REAL(rkind), ALLOCATABLE :: Ipolh(:)       ! twopi*isgn*Ipolh/mu0
                                             !   = net poloidal current[A]
                                             !   on half grids ; j=0,ns+1
  REAL(rkind), ALLOCATABLE :: vprimh(:)      ! (dV/ds)/twopi**2 on half grids,
                                             !   here V is volume [m3];
                                             !   j=0,ns+1
  REAL(rkind)              :: phiedge        ! (total toroidal flux at the
                                             !   edge)/(isgn*twopi) 

!... equilibrium quantities from vmec

  REAL(rkind), ALLOCATABLE :: iota_vmec(:)   ! rotational transform
                                             !   from vmec2000 ; j=0,ns
  REAL(rkind), ALLOCATABLE :: pprim(:)       ! dp/ds from vmec2000,
                                             !   where p[Pa] ; j=0,ns
  REAL(rkind), ALLOCATABLE :: S11(:)         ! S11 on a full mesh ; j=0,ns
  REAL(rkind), ALLOCATABLE :: S12(:)         ! S12 on a full mesh ; j=0,ns
  REAL(rkind), ALLOCATABLE :: Bsqav(:)       ! <B^2> [T2] on a full mesh ;
                                             !   j=0,ns
  REAL(rkind), ALLOCATABLE :: vprim(:)       ! dV/ds [m3] on a full mesh ;
                                             !   j=0,ns
  REAL(rkind), ALLOCATABLE :: grdssq_av(:)   ! <|grad(s)|^2> on a full mesh ;
                                             !   j=0,ns
  REAL(rkind), ALLOCATABLE :: Itorf(:)       ! net toroidal current [A]
                                             !   on a full mesh ; j=0,ns
  REAL(rkind), ALLOCATABLE :: Ipolf(:)       ! net poloidal current [A]
                                             !   on a full mesh ; j=0,ns
  REAL(rkind), ALLOCATABLE :: s(:)           ! radial grid points (full mesh)
                                             !   s is normalized tor. flux ;
                                             !    j=0:ns

  PRIVATE
  PUBLIC pl_vmec
      
CONTAINS

  SUBROUTINE pl_vmec(file_name,ierr)

    USE task_kinds
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT):: file_name
    INTEGER,INTENT(OUT):: ierr

    CALL read_vmec_equil(file_name,ierr)
    IF(ierr.NE.0) RETURN
    CALL put_vmec_bpsd
    RETURN

  END SUBROUTINE pl_vmec

  SUBROUTINE read_vmec_equil(file_name,ierr)

!  read wout-file by vmec2000 and calculate S11, S12, <B2>, p', V',
!    and so on.
!                                                  (July 2007, Y.N)

    USE task_kinds
    USE task_constants
    USE libfio
    IMPLICIT NONE

    CHARACTER(LEN=*),INTENT(INOUT):: file_name
    INTEGER         ,INTENT(OUT):: ierr
    CHARACTER(LEN=80) :: cdummy
    REAL(rkind)       :: dummy
    INTEGER           :: idmy, iasym, nbsets, ntheta2, nzeta
    INTEGER           :: j, mn, mxm, nxn, l, k, iarg
    REAL(rkind)       :: hs, dtheta, dzeta, arg, dnorm, ohs
    REAL(rkind)       :: R, Z, R_t, Z_t, R_z, Z_z, l_th, l_zh
    REAL(rkind)       :: g_tt, g_tz, g_zz, l_t, l_z, gsqrth, Bh, gsqrt, Bf
    REAL(rkind), ALLOCATABLE :: xtheta(:), xzeta(:)
    REAL(rkind), ALLOCATABLE :: tsin(:,:,:), tcos(:,:,:)
    REAL(rkind), ALLOCATABLE :: tsin2(:,:,:), tcos2(:,:,:)
    REAL(rkind), ALLOCATABLE :: l_th_save(:,:),l_zh_save(:,:)
    REAL(rkind), ALLOCATABLE :: gsqrth_save(:,:),Bh_save(:,:)

!......
!    open(8,file=trim(file_name),form='formatted')
    
    CALL fropen(8,TRIM(file_name),1,0,'pl_vmec',ierr)
    IF(ierr.NE.0) RETURN

    READ(8,'(a15,a4)') cdummy,version

    IF(version(1:4) == '6.90') THEN
       READ(8,*) (dummy,j=1,7)
       READ(8,*) nfp,nsp1,mpol,ntor,mnmax,idmy,idmy,iasym,idmy,idmy
       mnmax2 = mnmax
       IF(iasym == 1) THEN
          PRINT*,'asymmetric equilibrium' ; STOP
       ENDIF
       READ(8,*) idmy,idmy,idmy,idmy,idmy,idmy
       READ(8,'(a)') cdummy
    ELSE IF(version(1:4) == '8.40') THEN
       READ(8,*) (dummy,j=1,7)
       READ(8,*) nfp,nsp1,mpol,ntor,mnmax,mnmax2,idmy,idmy,iasym,idmy,idmy
       IF(iasym == 1) THEN
          PRINT*,'asymmetric equilibrium' ; STOP
       ENDIF
       READ(8,*) idmy,idmy,nbsets,idmy,idmy,idmy
       IF(nbsets > 0) READ(8,*) (dummy,j=1,nbsets)
       READ(8,'(a)') cdummy
    ELSE
       PRINT*,'unsupported version : version =',version
       STOP
    ENDIF

    ns      = nsp1 - 1
    nsm1    = ns - 1
    ntheta2 = 4*mpol+4
    nzeta   = 4*ntor+4

!... allocation of arrays
    ALLOCATE(xm(1:mnmax), xn(1:mnmax), xm2(1:mnmax2), xn2(1:mnmax2))
    ALLOCATE(rmnc(1:mnmax,0:ns), zmns(1:mnmax,0:ns))
    ALLOCATE(lmnsh(1:mnmax,0:ns), bmnh(1:mnmax2,0:ns), gmnh(1:mnmax2,0:ns))
    ALLOCATE(aiotah(0:nsp1), presh(0:nsp1), vprimh(0:nsp1), Itorh(0:nsp1), &
         Itorf(0:ns))
    ALLOCATE(Ipolh(0:nsp1), Ipolf(0:ns))
    ALLOCATE(iota_vmec(0:ns))
    ALLOCATE(xtheta(1:ntheta2), xzeta(1:nzeta))
    ALLOCATE(tsin(1:mnmax,1:nzeta,1:ntheta2))
    ALLOCATE(tcos(1:mnmax,1:nzeta,1:ntheta2))
    ALLOCATE(tsin2(1:mnmax2,1:nzeta,1:ntheta2))
    ALLOCATE(tcos2(1:mnmax2,1:nzeta,1:ntheta2))
    ALLOCATE(pprim(0:ns))
    ALLOCATE(l_th_save(1:nzeta,1:ntheta2), l_zh_save(1:nzeta,1:ntheta2))
    ALLOCATE(gsqrth_save(1:nzeta,1:ntheta2), Bh_save(1:nzeta,1:ntheta2))
    ALLOCATE(S11(0:ns), S12(0:ns), Bsqav(0:ns), vprim(0:ns), grdssq_av(0:ns))
    ALLOCATE(s(0:ns))

    DO j=0,ns
       IF(version == '6.90') THEN
          DO mn=1,mnmax                                           
             IF(j==0) THEN
                READ(8,*) mxm,nxn
                xm(mn)     = mxm ; xn(mn)     = nxn
                xm2(mn)    = mxm ; xn2(mn)    = nxn
             ENDIF
             READ(8,*) rmnc(mn,j),zmns(mn,j), &             ! full mesh
                       lmnsh(mn,j),bmnh(mn,j),gmnh(mn,j), & ! half mesh
                       bsubumnh(mn,j),bsubvmnh(mn,j), &     ! half mesh
                       bsubsmn(mn,j), &                     ! full mesh    
                       bsupumnh(mn,j),bsupvmnh(mn,j), &     ! half mesh
                       currvmn(mn,j)                        ! full mesh
          ENDDO
       ELSE IF(version == '8.40') THEN
          DO mn=1,mnmax                                           
             IF(j==0) THEN
                READ(8,*) mxm,nxn
                xm(mn)     = mxm ; xn(mn)     = nxn
             ENDIF
             READ(8,*) rmnc(mn,j),zmns(mn,j), &             ! full mesh
                       lmnsh(mn,j)                          ! half mesh
          ENDDO
          DO mn=1,mnmax2                                          
             IF(j==0) THEN
                READ(8,*) mxm,nxn
                xm2(mn) = mxm ; xn2(mn) = nxn
             ENDIF
             READ(8,*) bmnh(mn,j),gmnh(mn,j), &             ! half mesh
                       dummy,dummy, &                       ! half mesh
                       dummy, &                             ! full mesh    
                       dummy,dummy                          ! half mesh
          ENDDO
       ENDIF
    ENDDO

    IF(version == '6.90') THEN
       READ(8,*) (aiotah(j),dummy,presh(j),dummy,phiedge, &
            Itorh(j),Ipolh(j),dummy,vprimh(j), &
            dummy,dummy,dummy,dummy,j=1,ns)
       READ(8,*) (dummy,j=1,6)
       READ(8,*) isgn
       DO j=1,nsm1
          iota_vmec(j) = HALF*(aiotah(j)+aiotah(j+1))
       ENDDO
       aiotah(0)    = ONEHALF*aiotah(1) - HALF*aiotah(2)
       aiotah(nsp1) = ONEHALF*aiotah(ns) - HALF*aiotah(nsm1)
       iota_vmec(0)    = aiotah(0)
       iota_vmec(ns)   = aiotah(nsp1)

    ELSE IF(version == '8.40') THEN
       READ(8,*) (iota_vmec(j),dummy,dummy,dummy,dummy,dummy, j=0,ns) ! full
       READ(8,*) (aiotah(j),dummy,presh(j),dummy,phiedge, &
            Itorh(j),Ipolh(j),vprimh(j),dummy,dummy,j=1,ns)           ! half
       READ(8,*) (dummy,j=1,6)
       READ(8,*) isgn
       aiotah(0)    = iota_vmec(0)
       aiotah(nsp1) = iota_vmec(ns)
    ENDIF

    Itorh(0)     = ZERO
    Itorh(nsp1)  = ONEHALF*Itorh(ns)  - HALF*Itorh(nsm1)
    Ipolh(nsp1)  = ONEHALF*Ipolh(ns)  - HALF*Ipolh(nsm1)

    DO j=1,nsm1
       Itorf(j) = HALF*(Itorh(j)+Itorh(j+1))
       Ipolf(j) = HALF*(Ipolh(j)+Ipolh(j+1))
    ENDDO
    Itorf(0)  = Itorh(0)
    Itorf(ns) = Itorh(nsp1)
    Ipolf(0)  = Ipolh(1)
    Ipolf(ns) = Ipolh(nsp1)

    DO mn=1,mnmax                                                 
       IF(xm(mn) /= ZERO) THEN                                           
          rmnc(mn,0) = ZERO                                              
          zmns(mn,0) = ZERO                                              
       ENDIF
    ENDDO

!... set up equilibrium parameters                                                   
    hs     = ONE/ns
    dtheta  = HALF/(ntheta2-1)                                            
    dzeta   = ONE/(nzeta*nfp)                                               

    DO l = 1, ntheta2
       xtheta(l) = (l-1)*dtheta                                         
    ENDDO

    DO k = 1, nzeta                                               
       xzeta(k) = dzeta*(k-1)                                          
    ENDDO

    DO l = 1, ntheta2
       DO k = 1, nzeta
          DO mn = 1, mnmax
             arg  = xm(mn)*xtheta(l)-xn(mn)*xzeta(k)
             iarg = INT(arg)
             arg  = twopi * ( arg - iarg )
             tsin(mn,k,l) = sin(arg)
             tcos(mn,k,l) = cos(arg)
          ENDDO
       ENDDO
    ENDDO

    DO l = 1, ntheta2
       DO k = 1, nzeta
          DO mn = 1, mnmax2
             arg  = xm2(mn)*xtheta(l)-xn2(mn)*xzeta(k)
             iarg = INT(arg)
             arg  = twopi * ( arg - iarg )
             tsin2(mn,k,l) = sin(arg)
             tcos2(mn,k,l) = cos(arg)
          ENDDO
       ENDDO
    ENDDO
 
!... get p' on full grid points
     
    DO j=1,nsm1
       pprim(j) = (presh(j+1)-presh(j))*ns      ! dp/ds [Pa]
    ENDDO
    pprim(0)  = pprim(1)
    pprim(ns) = pprim(nsm1)    !  pprim(ns) = zero
                               ! here we assume p'=0 at plasma boundary

    dnorm = ONE/(nzeta*(ntheta2-1))
    ohs   = ns

    DO l=1,ntheta2
       DO k=1,nzeta
          l_th_save(k,l) = ZERO
          l_zh_save(k,l) = ZERO
          DO mn = 1, mnmax
             l_th_save(k,l) = l_th_save(k,l) + xm(mn)*lmnsh(mn,1)*tcos(mn,k,l)
             l_zh_save(k,l) = l_zh_save(k,l) - xn(mn)*lmnsh(mn,1)*tcos(mn,k,l)
          ENDDO
          gsqrth_save(k,l) = ZERO
          Bh_save(k,l)     = ZERO
          DO mn = 1, mnmax2
             gsqrth_save(k,l) = gsqrth_save(k,l) + gmnh(mn,1)*tcos2(mn,k,l)
             Bh_save(k,l)     = Bh_save(k,l)      + bmnh(mn,1)*tcos2(mn,k,l)
          ENDDO
       ENDDO
    ENDDO

    DO j=1,nsm1
       S11(j) = ZERO ; S12(j) = ZERO ; Bsqav(j) = ZERO
       vprim(j) = ZERO ; grdssq_av(j) = ZERO
       DO l=1,ntheta2
          DO k=1,nzeta
             R = ZERO ; Z = ZERO ; R_t = ZERO ; Z_t = ZERO
             R_z = ZERO ; Z_z = ZERO ; l_th = ZERO ; l_zh = ZERO
             DO mn = 1, mnmax
                R    = R    + rmnc(mn,j)*tcos(mn,k,l)
                Z    = Z    + zmns(mn,j)*tcos(mn,k,l)
                R_t  = R_t  - xm(mn)*rmnc(mn,j)*tsin(mn,k,l)
                R_z  = R_z  + xn(mn)*rmnc(mn,j)*tsin(mn,k,l)           
                Z_t  = Z_t  + xm(mn)*zmns(mn,j)*tcos(mn,k,l)           
                Z_z  = Z_z  - xn(mn)*zmns(mn,j)*tcos(mn,k,l) 
                l_th = l_th + xm(mn)*lmnsh(mn,j+1)*tcos(mn,k,l)           
                l_zh = l_zh - xn(mn)*lmnsh(mn,j+1)*tcos(mn,k,l) 
             ENDDO
             g_tt = R_t*R_t + Z_t*Z_t
             g_tz = R_t*R_z + Z_t*Z_z
             g_zz = R_z*R_z + Z_z*Z_z + R*R
             l_t  = HALF*(l_th_save(k,l)+l_th)
             l_z  = HALF*(l_zh_save(k,l)+l_zh)
             l_th_save(k,l) = l_th
             l_zh_save(k,l) = l_zh
                       
             gsqrth = ZERO ; Bh = ZERO
             DO mn = 1, mnmax2
                gsqrth = gsqrth + gmnh(mn,j+1)*tcos2(mn,k,l)
                Bh     = Bh     + bmnh(mn,j+1)*tcos2(mn,k,l)
             ENDDO
             gsqrt = HALF*(gsqrth_save(k,l)+gsqrth)
             Bf    = HALF*(Bh_save(k,l)+Bh)
             gsqrth_save(k,l) = gsqrth
             Bh_save(k,l)     = Bh

             IF(l==1.OR.l==ntheta2) THEN
                S11(j)    = S11(j) + HALF*g_tt/gsqrt
                S12(j)    = S12(j) + HALF*(g_tz*(1+l_t)-g_tt*l_z)/gsqrt
                Bsqav(j)  = Bsqav(j) + HALF*Bf*Bf*gsqrt
                vprim(j)  = vprim(j) + HALF*gsqrt
                grdssq_av(j) = grdssq_av(j) + HALF*(g_tt*g_zz-g_tz**2)/gsqrt
             ELSE
                S11(j)    = S11(j) + g_tt/gsqrt
                S12(j)    = S12(j) + (g_tz*(1+l_t)-g_tt*l_z)/gsqrt
                Bsqav(j)  = Bsqav(j) + Bf*Bf*gsqrt
                vprim(j)  = vprim(j) + gsqrt
                grdssq_av(j) = grdssq_av(j) + (g_tt*g_zz-g_tz**2)/gsqrt
             ENDIF
          ENDDO
       ENDDO
       S11(j)       = dnorm*S11(j)
       S12(j)       = dnorm*S12(j)
       Bsqav(j)     = Bsqav(j)/vprim(j)
       grdssq_av(j) = grdssq_av(j)/vprim(j)
       vprim(j)     = isgn*dnorm*vprim(j)
 
!<<check>>
!    print*,'vprim(j), (vprimh(j)+vprimh(j+1))/2 = ',
!    vprim(j),HALF*(vprimh(j)+vprimh(j+1))
 
    ENDDO

!... calculation for the plasma edge

    S11(ns) = ZERO
    S12(ns) = ZERO
    Bsqav(ns) = ZERO
    vprim(ns) = ZERO
    grdssq_av(ns) = ZERO
    DO l=1,ntheta2
       DO k=1,nzeta
          R = ZERO ; Z = ZERO ; R_t = ZERO ; Z_t = ZERO
          R_z = ZERO ; Z_z = ZERO ; l_th = ZERO ; l_zh = ZERO
          DO mn = 1, mnmax
             R    = R    + rmnc(mn,ns)*tcos(mn,k,l)
             Z    = Z    + zmns(mn,ns)*tcos(mn,k,l)
             R_t  = R_t  - xm(mn)*rmnc(mn,ns)*tsin(mn,k,l)
             R_z  = R_z  + xn(mn)*rmnc(mn,ns)*tsin(mn,k,l)           
             Z_t  = Z_t  + xm(mn)*zmns(mn,ns)*tcos(mn,k,l)           
             Z_z  = Z_z  - xn(mn)*zmns(mn,ns)*tcos(mn,k,l) 
             l_th = l_th + xm(mn)*lmnsh(mn,ns-1)*tcos(mn,k,l)           
             l_zh = l_zh - xn(mn)*lmnsh(mn,ns-1)*tcos(mn,k,l) 
          ENDDO
          g_tt = R_t*R_t + Z_t*Z_t
          g_tz = R_t*R_z + Z_t*Z_z
          g_zz = R_z*R_z + Z_z*Z_z + R*R
          l_t  = ONEHALF*l_th_save(k,l) - HALF*l_th
          l_z  = ONEHALF*l_zh_save(k,l) - HALF*l_zh
        
          gsqrth = ZERO ; Bh = ZERO
          DO mn = 1, mnmax2
             gsqrth = gsqrth + gmnh(mn,ns-1)*tcos2(mn,k,l)
             Bh     = Bh     + bmnh(mn,ns-1)*tcos2(mn,k,l)
          ENDDO
          gsqrt = ONEHALF*gsqrth_save(k,l) - HALF*gsqrth
          Bf    = ONEHALF*Bh_save(k,l)     - HALF*Bh
          IF(l==1.OR.l==ntheta2) THEN
             S11(ns)    = S11(ns) + HALF*g_tt/gsqrt
             S12(ns)    = S12(ns) + HALF*(g_tz*(1+l_t)-g_tt*l_z)/gsqrt
             Bsqav(ns)  = Bsqav(ns) + HALF*Bf*Bf*gsqrt
             vprim(ns)  = vprim(ns) + HALF*gsqrt
             grdssq_av(ns) = grdssq_av(ns) + HALF*(g_tt*g_zz-g_tz*g_tz)/gsqrt
          ELSE
             S11(ns)    = S11(ns) + g_tt/gsqrt
             S12(ns)    = S12(ns) + (g_tz*(1+l_t)-g_tt*l_z)/gsqrt
             Bsqav(ns)  = Bsqav(ns) + Bf*Bf*gsqrt
             vprim(ns)  = vprim(ns) + gsqrt
             grdssq_av(ns) = grdssq_av(ns) + (g_tt*g_zz-g_tz*g_tz)/gsqrt
          ENDIF
       ENDDO
    ENDDO
    S11(ns)   = dnorm*S11(ns)
    S12(ns)   = dnorm*S12(ns)
    Bsqav(ns) = Bsqav(ns)/vprim(ns)
    grdssq_av(ns) = grdssq_av(ns)/vprim(ns)
    vprim(ns) = isgn*dnorm*vprim(ns)

    S11(0)       = ZERO
    S12(0)       = ZERO
    Bsqav(0)     = two*Bsqav(1) - Bsqav(2)
    vprim(0)     = two*vprim(1) - vprim(2)
    grdssq_av(0) = ZERO

    DO j = 1, ns
       s(j)     = j*hs                                          
    ENDDO
    s(0)      = ZERO  ; s(ns)        = ONE                                                 
!... force iota_vmec to be (Itor(vmec)/phiedge-S12)/S11

    DO j = 1, ns
       iota_vmec(j) = (Itorf(j)/phiedge-S12(j))/S11(j)
    ENDDO
    iota_vmec(0) = 2*iota_vmec(1) - iota_vmec(2)

!... check consistency
!      N.B. Unit of Itorf & Ipolf is not [A] and
!      vprim=d(V/twopi**2)/ds at this point.

    PRINT*, '* Phi_tor_edge = twopi*isgn*phiedge [Wb] =', twopi*isgn*phiedge
    PRINT*, '* I_tor_edge (in vmec)  [kA]             =', &
         0.001_dp*twopi*isgn/RMU0*Itorf(ns)
    PRINT*
    PRINT*,'      s           S11          S12           -S12/S11      iota(vmec)       Itor[kA]        Itor[kA]'
    PRINT*,'                                                                       phip*(S11*iota+S12)   (vmec)'
    DO j = 1,ns
       WRITE(6,'(f10.4,2(e15.5),4(f15.8))') &
            s(j),S11(j),S12(j),-S12(j)/S11(j),iota_vmec(j), &
            0.001_dp*twopi*isgn/RMU0*phiedge*(S11(j)*iota_vmec(j)+S12(j)), &
            0.001_dp*twopi*isgn/RMU0*Itorf(j)
    ENDDO
    PRINT*

    PRINT*,'      s        pprim        pprim_check       <B^2>       <B^2>_check      vprim    <|grad(s)|^2>'
    DO j = 1,nsm1
       WRITE(6,'(f10.4,6(e15.5))') s(j),pprim(j),&
            -isgn*phiedge*((Ipolh(j+1)-Ipolh(j)) &
                          +(Itorh(j+1)-Itorh(j))*iota_vmec(j)) &
                         /hs/RMU0/vprim(j),&
            Bsqav(j), &
            isgn*(Itorf(j)*iota_vmec(j)+Ipolf(j))*phiedge/vprim(j), &
            vprim(j),grdssq_av(j)
    ENDDO
    PRINT*
    
!... physical value

    DO j=0,ns
       Itorf(j) = twopi*isgn/RMU0*Itorf(j)  ! net toroidal current [A]
       Ipolf(j) = twopi*isgn/RMU0*Ipolf(j)  ! net poloidal current [A]
       vprim(j) = twopi**2*vprim(j)         ! dV/ds [m3]
    ENDDO

    DEALLOCATE(xtheta, xzeta)
    DEALLOCATE(tsin, tcos)
    DEALLOCATE(tsin2, tcos2)
    DEALLOCATE(l_th_save, l_zh_save)
    DEALLOCATE(gsqrth_save, Bh_save)

    CLOSE(8)

  END SUBROUTINE read_vmec_equil

!=======================================================================

  SUBROUTINE put_vmec_bpsd

    USE task_kinds
    USE task_constants
    IMPLICIT NONE
    INTEGER       :: ierr, nr, mn
    LOGICAL, SAVE :: init_flag

    IF(init_flag) THEN
       equ1D%nrmax=0
       metric1D%nrmax=0
       init_flag=.FALSE.
    ENDIF

    DO mn=1,mnmax
       IF(xm(mn).EQ.0.d0.AND.xn(mn).EQ.0.d0) THEN
          device%rr=rmnc(mn,0)
          device%zz=zmns(mn,0)
       ENDIF
       IF(xm(mn).EQ.1.d0.AND.xn(mn).EQ.0.d0) THEN
          device%ra=rmnc(mn,ns)
          device%rb=device%ra*1.2d0
          device%elip=rmnc(mn,ns)/zmns(mn,ns)
       ENDIF
    ENDDO

    device%bb=sqrt(Bsqav(0))
    device%ip=Itorf(ns)
    device%trig=0.d0
    bpsd_debug_flag=.true.
    CALL bpsd_put_data(device,ierr)
    bpsd_debug_flag=.false.

    equ1D%time=0.D0
    IF(equ1D%nrmax.NE.ns+1) THEN
       IF(ALLOCATED(equ1D%rho)) DEALLOCATE(equ1D%rho)
       IF(ALLOCATED(equ1D%data)) DEALLOCATE(equ1D%data)
       equ1D%nrmax=ns+1
       ALLOCATE(equ1D%rho(ns+1))
       ALLOCATE(equ1D%data(ns+1))
    ENDIF

    metric1D%time=0.D0
    IF(metric1D%nrmax.NE.ns+1) THEN
       IF(ALLOCATED(metric1D%rho)) DEALLOCATE(metric1d%rho)
       IF(ALLOCATED(metric1D%data)) DEALLOCATE(metric1d%data)
       metric1D%nrmax=ns+1
       ALLOCATE(metric1D%rho(ns+1))
       ALLOCATE(metric1D%data(ns+1))
    ENDIF

    DO nr=0,ns
       equ1D%rho(nr+1)=SQRT(s(nr))
       equ1D%data(nr+1)%psit=s(nr)
       equ1D%data(nr+1)%psip=0.d0
       equ1D%data(nr+1)%ppp=pprim(nr)
       equ1D%data(nr+1)%piq=iota_vmec(nr)
       equ1D%data(nr+1)%pip=Ipolf(nr)
       equ1D%data(nr+1)%pit=Itorf(nr)
    ENDDO

    CALL bpsd_put_data(equ1D,ierr)

    DO nr=0,ns
       metric1D%rho(nr+1)=SQRT(s(nr))
       metric1D%data(nr+1)%pvol=vprim(nr)
       metric1D%data(nr+1)%psur=vprim(nr)/(2.d0*pi*device%rr)
       metric1D%data(nr+1)%dvpsit=vprim(nr)
       metric1D%data(nr+1)%dvpsip=vprim(nr)/iota_vmec(nr)
       metric1D%data(nr+1)%aver2= device%rr**2
       metric1D%data(nr+1)%aver2i=1.d0/device%rr**2
       metric1D%data(nr+1)%aveb2= Bsqav(nr)
       metric1D%data(nr+1)%aveb2i=1.d0/Bsqav(nr)
       metric1D%data(nr+1)%avegv2=vprim(nr)**2*grdssq_av(nr)
       metric1D%data(nr+1)%avegvr2=vprim(nr)**2*grdssq_av(nr) &
                                    /device%rr**2
       metric1D%data(nr+1)%avegpp2=0.d0
       metric1D%data(nr+1)%rr=device%rr
       metric1D%data(nr+1)%rs=device%ra*sqrt(s(nr))
       metric1D%data(nr+1)%elip=device%elip
       metric1D%data(nr+1)%trig=0.d0
    ENDDO
    CALL bpsd_put_data(metric1D,ierr)
  END SUBROUTINE put_vmec_bpsd

END MODULE plvmec
