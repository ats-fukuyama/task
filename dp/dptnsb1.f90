! dprnsb1.f90

MODULE dptnsb1

!     ***************************************
!         COMPONENTS OF DIELECTRIC TENSOR
!             MAGNETIC FIELD     (0,   0,   B)
!             WAVE NUMBER VECTOR (k_x, 0, k_z)

!             CLDISP(1)=EPS_XX
!             CLDISP(2)=EPS_ZZ - EPS_XX
!             CLDISP(3)=EPS_YY - EPS_XX
!             CLDISP(4)=EPS_ZX
!             CLDISP(5)=EPS_XY
!             CLDISP(6)=EPS_YZ

  PRIVATE
  PUBLIC dp_tnsb1,test_dpbes

CONTAINS

!     ****** Cold beam model with FLR ******

  SUBROUTINE dp_tnsb1(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

    USE libbes,ONLY: BESEINX
    USE dpcomm
    USE plprof
    USE plprofw
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
    INTEGER,INTENT(IN):: NS
    TYPE(pl_mag_type),INTENT(IN):: mag
    TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
    COMPLEX(rkind),INTENT(OUT):: CLDISP(6)

    REAL(rkind),ALLOCATABLE:: BESJN(:),DBESJN(:),DDBESJN(:)
    COMPLEX(rkind):: cwp2,cwc,crkprprww,crkppprwc,crkprrkpp2
    REAL(rkind):: omegac,rkpr,rkpp,zeta
    REAL(rkind):: bdbz,bdbx
    COMPLEX(rkind):: cfact1,cfact2,cxx,cxy,cxz,cyy,cyz,czz
    INTEGER:: nhmax,nh,ierr
      
    cwp2=plfw(ns)%rn*1.D20*pz(ns)*pz(ns)*AEE*AEE/(EPS0*AMP*pa(ns)*cw*cw)
    omegac=mag%babs*pz(ns)*AEE/(AMP*pa(ns))
    cwc=omegac/cw
    
    rkpr=REAL(ckpr)
    rkpp=REAL(ckpp)
    zeta=rkpp*plfw(ns)%rupp/omegac
    crkprprww=rkpr*plfw(ns)%rupr/cw
    crkppprwc=rkpp*plfw(ns)%rupr/omegac
    IF(ABS(rkpp).LT.1.D-80) THEN
       crkprrkpp2=(rkpr*cwc/1.D-4)**2
    ELSE
       crkprrkpp2=(rkpr*cwc/rkpp)**2
    END IF

    nhmax=MAX(ABS(ncmin(ns)),ABS(ncmax(ns)))
    ALLOCATE(BESJN(-nhmax:nhmax))
    ALLOCATE(DBESJN(-nhmax:nhmax))
    ALLOCATE(DDBESJN(-nhmax:nhmax))
    CALL DPBESJ(nhmax,zeta,BESJN,DBESJN,DDBESJN,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I6,ES12.4,I6)') 'nh,zeta,ierr=',nh,zeta,ierr
       STOP
    END IF
    
    cxx=1.D0
    cxy=0.D0
    cxz=0.D0
    cyy=1.D0
    cyz=0.D0
    czz=1.D0
    DO nh=-nhmax,nhmax
       cfact1=nh*cwc/(1.D0-nh*cwc-crkprprww)
       cfact2=crkprrkpp2/(1.D0-nh*cwc-crkprprww)**2
       IF(ABS(zeta).LT.1.D-80) THEN
          SELECT CASE(nh)
          CASE(0)
             bdbz=-0.5D0
          CASE(-1,1)
             bdbz= 0.25D0
          CASE DEFAULT
             bdbz= 0.D0
          END SELECT
       ELSE
          bdbz=BESJN(nh)*DBESJN(nh)/zeta
       END IF
       bdbx=bdbz+BESJN(nh)*DDBESJN(nh)+DBESJN(nh)**2

       cxx=cxx+cfact1*2.D0*nh**2*bdbz &
              +cfact2*nh**2*BESJN(nh)**2
       cxy=cxy+cfact1*CI*nh*bdbx &
              +cfact2*CI*nh*zeta*BESJN(nh)*DBESJN(nh)
       cxz=cxz+cfact1*2.D0*nh*crkppprwc*bdbz &
              +cfact2     *nh*crkppprwc*BESJN(nh)**2
       cyy=cyy+cfact1*2.D0 &
                     *(zeta*DBESJN(nh)*DDBESJN(nh)+DBESJN(nh)**2) &
              +cfact2*zeta**2*DBESJN(nh)**2
       cyz=cyz-cfact1*CI*crkppprwc*bdbx &
              -cfact2*CI*crkppprwc*zeta*BESJN(nh)*DBESJN(nh)
       czz=czz+cfact1*2.D0*crkppprwc**2*bdbz &
              +cfact2     *crkppprwc**2*BESJN(nh)**2
    END DO
    CLDISP(1)=-cwp2*cxx
    CLDISP(2)=-cwp2*(czz-cxx)
    CLDISP(3)=-cwp2*(cyy-cxx)
    CLDISP(4)=-cwp2*cxz
    CLDISP(5)=-cwp2*cxy
    CLDISP(6)=-cwp2*cyz

!    WRITE(21,'(3ES12.4)') ww,rkpr,rkpp
!    WRITE(21,'(3ES12.4)') wp2,omegac,zeta
!    WRITE(21,'(6ES12.4)') cldisp(1),cldisp(2),cldisp(3)
!    WRITE(21,'(6ES12.4)') cldisp(4),cldisp(5),cldisp(6)
    RETURN
  END SUBROUTINE dp_tnsb1

  SUBROUTINE test_dpbes
    USE dpcomm,ONLY: rkind
    IMPLICIT NONE
    
    REAL(rkind),ALLOCATABLE:: BESJN(:),DBESJN(:),DDBESJN(:)
    INTEGER:: nhmax=5
    REAL(rkind):: x=10.D0
    INTEGER:: nh,ierr

1   CONTINUE
    WRITE(6,'(A,I6,ES12.4)') '## input nhmax and x:',nhmax,x
    READ(5,*,ERR=1,END=9000) nhmax,x
    
    ALLOCATE(BESJN(-nhmax:nhmax))
    ALLOCATE(DBESJN(-nhmax:nhmax))
    ALLOCATE(DDBESJN(-nhmax:nhmax))

    CALL DPBESJ(nhmax,x,BESJN,DBESJN,DDBESJN,ierr)

    DO nh=-nhmax,nhmax
       WRITE(6,'(I6,4ES12.4)') nh,BESJN(nh),DBESJN(nh),DDBESJN(nh), &
            BESJN(nh)*DBESJN(nh)/x
    END DO

    DEALLOCATE(BESJN,DBESJN,DDBESJN)
    
    GO TO 1
    
9000 CONTINUE
    RETURN
  END SUBROUTINE test_dpbes

  SUBROUTINE DPBESJ(nhmax,zeta,BESJN_,DBESJN_,DDBESJN_,ierr)

    USE dpcomm,ONLY: rkind
    USE dpsub
    USE libbes
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nhmax
    REAL(rkind),INTENT(IN):: zeta
    REAL(rkind),INTENT(OUT):: BESJN_(-nhmax:nhmax)
    REAL(rkind),INTENT(OUT):: DBESJN_(-nhmax:nhmax)
    REAL(rkind),INTENT(OUT):: DDBESJN_(-nhmax:nhmax)
    INTEGER,INTENT(OUT):: ierr
    REAL(rkind),ALLOCATABLE:: BESJN(:),DBESJN(:),DDBESJN(:),temp(:)
    INTEGER:: nh,nh80,nh32
    REAL(rkind):: alzeta

    ALLOCATE(BESJN(-nhmax-2:nhmax+2))
    ALLOCATE(DBESJN(-nhmax-1:nhmax+1))
    ALLOCATE(DDBESJN(-nhmax:nhmax))
    ALLOCATE(temp(0:nhmax+2))

    ! --- avoid too small value:  0 for for LOG10(zeta**nh) < -80
    !                             power expression for LOG10(zeta**nh) < -32
    
    IF(zeta.LT.1.D-80) THEN
       temp(0)=1.D0
       DO nh=1,nhmax+2
          temp(nh)=0.D0
       END DO
       ierr=0
    ELSE IF(zeta.GT.0.1D0) THEN
       CALL BESJNV(nhmax+2,zeta,temp,ierr)
       IF(ierr.NE.0) RETURN
    ELSE 
       alzeta=LOG10(zeta)
       nh80=NINT(-80.D0/alzeta)
       nh32=NINT(-32.D0/alzeta)
       IF(nhmax+2.LE.nh32) THEN
          CALL BESJNV(nhmax+2,zeta,temp,ierr)
          IF(ierr.NE.0) RETURN
       ELSE IF(nhmax+2.LE.nh80) THEN
          CALL BESJNV(nh32,zeta,temp,ierr)
          IF(ierr.NE.0) RETURN
          DO nh=nh32+1,nhmax+2
             temp(nh)=(0.5D0*zeta)**nh/dkaijou(nh)
          END DO
       ELSE
          CALL BESJNV(nh32,zeta,temp,ierr)
          IF(ierr.NE.0) RETURN
          DO nh=nh32+1,nh80
             temp(nh)=(0.5D0*zeta)**nh/dkaijou(nh)
          END DO
          DO nh=nh80+1,nhmax+2
             temp(nh)=0.D0
          END DO
       END IF
    END IF
    
    DO nh=0,nhmax+2
       BESJN(nh)=temp(nh)
    END DO
    DO nh=1,nhmax+2
       BESJN(-nh)=(-1)**(nh)*temp(nh)
    END DO
    DO nh=-nhmax-1,nhmax+1
       DBESJN(nh)=0.5D0*(BESJN(nh-1)-BESJN(nh+1))
    END DO
    DO nh=-nhmax,nhmax
       DDBESJN(nh)=0.5D0*(DBESJN(nh-1)-DBESJN(nh+1))
    END DO

    DO nh=-nhmax,nhmax
       BESJN_(nh)=BESJN(nh)
       DBESJN_(nh)=DBESJN(nh)
       DDBESJN_(nh)=DDBESJN(nh)
    END DO
  END SUBROUTINE DPBESJ
END MODULE dptnsb1

