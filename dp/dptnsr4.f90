MODULE DPTNSR4

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
  PUBLIC DPTNTW

CONTAINS

!     ****** KINETIC MODEL WITH FLR and Drift (by T. Watanabe) ******

  SUBROUTINE DPTNTW(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFNV
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
      COMPLEX(rkind),ALLOCATABLE:: CALAM(:)

      COMPLEX(rkind):: CWP,CBETA,CGZ,CZ,CDZ,CDDZ,CDDDZ,CKPRL
      COMPLEX(rkind):: CWN,CFN,CGN,CHIXX,CHIYY,CHIZZ,CHIXY,CHIYZ,CHIZX
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM
      REAL(rkind):: UPR,WC,WTPR,VTPR,WTPP
      INTEGER:: I,NHMAX,NC,NH

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)))+1
      ALLOCATE(CALAM(0:NHMAX))

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      WC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))     ! WC is not devided by CW
      UPR=plfw(NS)%RUPR                           ! u_para

      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))  ! v_para^2
      VTPR=SQRT(2.D0*WTPR)                     ! v_para
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))  ! v_perp^2

      CBETA=CKPP*CKPP*WTPP/(WC*WC)           ! l.c. lambda
      DO NH=0,NHMAX
         CALAM(NH)=BESEINX(NH,DBLE(CBETA))     ! I_n(cbeta) EXP(-cbeta)
      END DO

      DO NC=NCMIN(NS),NCMAX(NS)

         IF(ABS(CKPR).LE.0.D0) THEN
            CKPRL=1.D-4
         ELSE
            CKPRL=CKPR
         END IF
         CGZ= (CW+NC*WC-CKPR*UPR)/(ABS(CKPRL)*VTPR)
         CALL DSPFNV(CGZ,CZ,CDZ,CDDZ,CDDDZ)

         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM      ! d Lambda/dx
         CBLAM=0.5D0*(CLAMM-CLAMP)           ! n Lambda/x

         CWN=(CW-CKPR*UPR)/(ABS(CKPRL)*VTPR)*CZ+0.5D0*(1.D0-WTPP/WTPR)*CDZ
         CFN=NC*CKPR*UPR/(ABS(CKPRL)*VTPR)*CZ+0.5D0*((CW/WC)*(WTPP/WTPR)+NC*(WTPP/WTPR-1.D0))*CDZ
         CGN=(CW+NC*WC*(1.D0-WTPR/WTPP-UPR*UPR/WTPP))/(ABS(CKPR)*VTPR)*CZ &
            -CKPR*UPR*(CW+NC*WC*(1.D0-2.D0*WTPR/WTPP))/(ABS(CKPR)*VTPR)**2*CDZ &
            +0.5D0*(CW+NC*WC*(1.D0-2.D0*WTPR/WTPP))/(ABS(CKPR)*VTPR)*CDDZ

         CHIXX=NC*CBLAM*CWN
         CHIYY=(NC*CBLAM-2.D0*CKPP**2*WTPP/WC**2*CDLAM)*CWN
         CHIZZ=CLAM*CGN
         CHIXY=-CI*NC*CDLAM*CWN
         CHIYZ= CI*CKPP/CKPR*CDLAM*CFN
         CHIZX= CKPP/CKPR*CBLAM*CFN

         CLDISP(1)=CLDISP(1)+CWP*CHIXX
         CLDISP(2)=CLDISP(2)+CWP*(CHIZZ-CHIXX)
         CLDISP(3)=CLDISP(3)+CWP*(CHIYY-CHIZZ)
         CLDISP(4)=CLDISP(4)+CWP*CHIZX
         CLDISP(5)=CLDISP(5)+CWP*CHIXY
         CLDISP(6)=CLDISP(6)+CWP*CHIYZ

      ENDDO
      RETURN
    END SUBROUTINE DPTNTW
  END MODULE DPTNSR4
