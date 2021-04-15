MODULE DPTNSR1

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

  PUBLIC DPTNCL,DPTNCC,DPTNIM,DPTNRM,DPTNWP,DPTNUP,DPTNHP
  PUBLIC DPTNKL,DPTNKS,DPTNKP,DPTNKR

CONTAINS

!     ****** COLLISIONLESS COLD MODEL ******

  SUBROUTINE DPTNCL(CW,NS,mag,plfw,CLDISP)

      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CYY,CZZ

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
      CYY=1.D0/(1.D0-CWC*CWC)

!      CLDISP(1)=-CWP/(1.D0-CWC*CWC)   XXX Nag Fortran compiler bug XXX
      CLDISP(1)=-CWP*CYY
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
!      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC*CYY
      CLDISP(6)= 0.D0

      RETURN
  END SUBROUTINE DPTNCL

!     ****** COLLISIONAL COLD MODEL ******

  SUBROUTINE DPTNCC(CW,NS,mag,plfw,CLDISP)

      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CWNU

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF
      CWP=CWP/CWNU
      CWC=CWC/CWNU

      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)= 0.D0
      RETURN
  END SUBROUTINE DPTNCC

!     ****** IDEAL MHD MODEL ******

  SUBROUTINE DPTNIM(CKPR,NS,mag,plfw,CLDISP)

      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CKPR
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CPARA
      REAL(rkind):: RKPR,VTE

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS))
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))

      RKPR=DBLE(CKPR)
      IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
      VTE=SQRT(plfw(NS)%RTPR*AEE*1.D3/(PA(NS)*AMP))
      CPARA=plfw(NS)%RN*1.D20*AEE*AEE &
           /(EPS0*AMP*PA(NS)*RKPR**2*VTE**2)

      CLDISP(1)= CWP/(CWC*CWC)
      CLDISP(2)= CPARA-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)= 0.D0
      CLDISP(6)= 0.D0
      RETURN
  END SUBROUTINE DPTNIM

!     ****** RESISTIVE MHD MODEL ******

  SUBROUTINE DPTNRM(CW,NS,mag,plfw,CLDISP)

      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CWNU

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF
      CWP=CWP/CWNU
      CWC=CWC/CWNU

      CLDISP(1)= CWP/(CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)= 0.D0
      CLDISP(6)= 0.D0
      RETURN
  END SUBROUTINE DPTNRM

!     ****** WARM PLAMSA MODEL ******

  SUBROUTINE DPTNWP(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CWNU,CTPR,CTPP
      REAL(rkind):: WTPR,WTPP,WTPX

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF
      CWP=CWP/CWNU
      CWC=CWC/CWNU

      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)

      CTPR=1.D0/(1.D0-CKPR**2*WTPR/(CW*CW*CWNU))
      CTPP=WTPP/((1.D0-CWC*CWC)*CW*CW*CWNU)

      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP*CTPR*(1.D0-CKPP**2*CTPP)&
                -CLDISP(1)
      CLDISP(3)= CWP*CKPP**2*CTPP*CTPR
      CLDISP(4)=-CKPP*CKPR*CTPP*CTPR*WTPX
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)=(0.D0, 1.D0)*CWC*CLDISP(4)
      RETURN
  END SUBROUTINE DPTNWP

!     ******  UNMAGNETIZE KINETIC DISPERSION ******

  SUBROUTINE DPTNUP(CW,CKPR,CKPP,NS,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CK2,CGZ,CZ,CDZ,CWP,CFX,CFZ,CLP,CLT,CWNU
      REAL(rkind):: RT,VT2

      RT=(plfw(NS)%RTPR+2.D0*plfw(NS)%RTPP)/3.D0
      CK2=SQRT(CKPR**2+CKPP**2)
      VT2=RT*AEE*1.D3/(AMP*PA(NS))
      CGZ=CW/SQRT(2.D0*CK2*VT2)
      CALL DSPFN(CGZ,CZ,CDZ)

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF
      CWP=CWP/CWNU

      CFX=CKPP/CK2
      CFZ=CKPR/CK2
      CLP=-CWP*CGZ**2*CDZ
      CLT= CWP*CGZ*CZ

      CLDISP(1)=  CFZ**2*CLT+CFX**2*CLP
      CLDISP(2)= (CFZ**2-CFX**2)*(CLP-CLT)
      CLDISP(3)= -CFX**2*(CLP-CLT)
      CLDISP(4)= CFX*CFZ*(CLP-CLT)
      CLDISP(5)= 0.D0
      CLDISP(6)= 0.D0
      RETURN
  END SUBROUTINE DPTNUP

!     ****** KINETIC MODEL WITHOUT FLR ******

  SUBROUTINE DPTNHP(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CALAM(0:2)
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NC

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF

      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=0.D0
      CALAM(0)=1.D0
      CALAM(1)=0.D0
      CALAM(2)=0.D0

      DO NC=-1,1
         CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
         CGZ= (CWNU-NC*CWC)*CPR
         CALL DSPFN(CGZ,CZ,CDZ)

         CK=CKPP/CKPR

         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-NC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-NC*(1.D0-1.D0/WTPX))*CGZ*CDZ

         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)

         CLDISP(1)=CLDISP(1)+CWP*   NC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-NC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
      ENDDO
      RETURN
  END SUBROUTINE DPTNHP

!     ****** KINETIC MODEL WITH LOWEST FLR ******

  SUBROUTINE DPTNKL(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
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
      COMPLEX(rkind):: CALAM(0:3)
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NC,NHMAX

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)))+1
      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF

      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      DO NC=0,NHMAX
         CALAM(NC)=BESEINX(NC,REAL(CBETA))
      END DO

      DO NC=NCMIN(NS),NCMAX(NS)
         CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
         CGZ= (CWNU-NC*CWC)*CPR
         CALL DSPFN(CGZ,CZ,CDZ)

         CK=CKPP/CKPR

         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-NC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-NC*(1.D0-1.D0/WTPX))*CGZ*CDZ

         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)

         CLDISP(1)=CLDISP(1)+CWP*   NC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-NC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
      ENDDO
      RETURN
  END SUBROUTINE DPTNKL

!     ****** KINETIC MODEL WITH FLR (SYMMETRIC) ******

  SUBROUTINE DPTNKS(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
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
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NHMAX,NC,NH
      COMPLEX(rkind):: C1,C2,C3

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)))+1
      ALLOCATE(CALAM(0:NHMAX))

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF

      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      C1=CKPP*CKPP*WTPP
      C2=CWC*CWC*CW*CW
      C3=1.D0/C2
      CBETA=C1*C3
!      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      DO NH=0,NHMAX
         CALAM(NH)=BESEINX(NH,REAL(CBETA))
      END DO

      DO NC=NCMIN(NS),NCMAX(NS)
         IF(ABS(CKPR).LE.0.D0) THEN
            CPR=CW/SQRT(2.D0*1.D-4**2*WTPR)
            CK=CKPP/1.D-4
         ELSE
            C1=2.D0*CKPR**2*WTPR
            C2=SQRT(C1)
            C3=1.D0/C2
            CPR=CW*C3
!            CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
            C3=1.D0/CKPR
            CK=CKPP*C3
!            CK=CKPP/CKPR
        ENDIF
         CGZ= (CWNU-NC*CWC)*CPR
         CALL DSPFN(CGZ,CZ,CDZ)

         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-NC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-NC*(1.D0-1.D0/WTPX))*CGZ*CDZ

         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)

         CLDISP(1)=CLDISP(1)+CWP*   NC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-NC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
!         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(4)=0.D0
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
!         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
         CLDISP(6)=0.D0
      ENDDO
      RETURN
    END SUBROUTINE DPTNKS

!     ****** KINETIC MODEL WITH FLR ******

    SUBROUTINE DPTNKP(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
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
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NHMAX,NC,NH

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NHMAX=MAX(ABS(NCMIN(NS)),ABS(NCMAX(NS)))+1
      ALLOCATE(CALAM(0:NHMAX))

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF
      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      DO NH=0,NHMAX
         CALAM(NH)=BESEINX(NH,REAL(CBETA))
      END DO
!      CALL LAMBDA(NHMAX,CBETA,CALAM,IERR)
!      IF(IERR.EQ.1) WRITE(6,*) 'XX LAMBDA: N out of range'
!      IF(IERR.EQ.2) WRITE(6,*) 'XX LAMBDA: CBETA out of range'

      DO NC=NCMIN(NS),NCMAX(NS)

         IF(ABS(CKPR).LE.0.D0) THEN
            CPR=CW/SQRT(2.D0*1.D-4**2*WTPR)
            CK=CKPP/1.D-4
         ELSE
            CPR=CW/SQRT(2.D0*CKPR**2*WTPR)
            CK=CKPP/CKPR
         ENDIF
         CGZ= (CWNU-NC*CWC)*CPR
         CALL DSPFN(CGZ,CZ,CDZ)

         CFW=CPR*CZ+0.5D0*(1.D0-WTPX)*CDZ
         CFF=0.5D0*(WTPX/CWC-NC*(WTPX-1.D0))*CDZ
         CFG=-CPR*(1.D0-NC*(1.D0-1.D0/WTPX))*CGZ*CDZ

         CLAM=CALAM(ABS(NC))
         CLAMM=CALAM(ABS(NC-1))
         CLAMP=CALAM(ABS(NC+1))
         CDLAM=0.5D0*(CLAMM+CLAMP)-CLAM
         CBLAM=0.5D0*(CLAMM-CLAMP)

         CLDISP(1)=CLDISP(1)+CWP*   NC*CBLAM*CFW
         CLDISP(2)=CLDISP(2)+CWP*     (CLAM*CFG-NC*CBLAM*CFW)
         CLDISP(3)=CLDISP(3)-CWP* 2.D0*CBETA*CDLAM*CFW
         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
      ENDDO
      RETURN
  END SUBROUTINE DPTNKP

!     ****** CALCULATE DETERMINANT OF DISPERSION TENSOR ******
!     ******                    IN                      ******
!     ******     WEAKLY RELATIVISTIC THERMAL PLASMA     ******

  SUBROUTINE DPTNKR(CW,CKPR,CKPP,NS,mag,plfw,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE dpsub
      USE plprof
      USE plprofw
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plfw_type),DIMENSION(nsmax),INTENT(IN):: plfw
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
!      COMPLEX(rkind),ALLOCATABLE:: CALAM(:)
      COMPLEX(rkind):: CWP,CWC,CBETA,CNPP,CNPR,CZ,CWNU
      COMPLEX(rkind):: CSUM1,CSUM2,CSUM3,CSUM4,CSUM5
      COMPLEX(rkind):: CN1,CN2,CN3,CF1,CF2,CF3,CF4,CF5,CF6,CF7,CPART
      COMPLEX(rkind):: CE11,CE12,CE22,CE23,CE31,CE33
      REAL(rkind):: DELZ,WTPR,WTPP,RMU,DKAI
      INTEGER:: NC,N,NSIG

      DELZ=1.D-6

      CWP=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=ABS(mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW))

      IF(plfw(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plfw(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plfw(NS)%RNUC/CW
      END IF

      WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plfw(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)

      RMU=VC**2/WTPR
      CNPP=CKPP*VC/CW
      CNPR=CKPR*VC/CW

      CSUM1=(0.D0,0.D0)
      CSUM2=(0.D0,0.D0)
      CSUM3=(0.D0,0.D0)
      CSUM4=(0.D0,0.D0)
      CSUM5=(0.D0,0.D0)

      DO NC=NCMIN(NS),NCMAX(NS)
         IF(NC.NE.0)THEN
            N=ABS(NC)
            NSIG=NC/N
            CZ=NC*RMU*CWC
            DKAI=DKAIJOU(N)
            CN1=N**2*CBETA**(N-1)/(2**N*DKAI)
            CN2=N   *CBETA**(N-1)/(2**N*DKAI)
            CN3=     CBETA**N    /(2**N*DKAI)

            CF1=CFQ(N+3.D0/2.D0,CZ     ,CNPR     ,RMU)
            CF2=CFQ(N+5.D0/2.D0,CZ+DELZ,CNPR     ,RMU)
            CF3=CFQ(N+5.D0/2.D0,CZ-DELZ,CNPR     ,RMU)
            CF4=CFQ(N+5.D0/2.D0,CZ     ,CNPR+DELZ,RMU)
            CF5=CFQ(N+5.D0/2.D0,CZ     ,CNPR-DELZ,RMU)

            CSUM1=CSUM1+     CN1*CF1
            CSUM2=CSUM2+NSIG*CN1*CF1
            CSUM3=CSUM3+     CN2*(CF2-CF3)/(DELZ*2.D0)
            CSUM4=CSUM4+NSIG*CN2*(CF2-CF3)/(DELZ*2.D0)
            CSUM5=CSUM5+     CN3*(CF4*(CNPR+DELZ)-CF5*(CNPR-DELZ)) &
                             /(DELZ*2.D0)
         ENDIF
      ENDDO
      CF6=CFQ(  5.D0/2.D0,(0.D0,0.D0),CNPR+DELZ,RMU)
      CF7=CFQ(  5.D0/2.D0,(0.D0,0.D0),CNPR-DELZ,RMU)
      CPART=(CF6*(CNPR+DELZ)-CF7*(CNPR-DELZ))/(DELZ*2.D0)

      CE11=   -CWP*RMU          *CSUM1
      CE12= CI*CWP*RMU          *CSUM2
      CE22=   -CWP*RMU          *CSUM1
      CE23=-CI*CWP/CWC*RMU*CNPP*CNPR*CSUM4
      CE31=   -CWP/CWC*RMU*CNPP*CNPR*CSUM3
      CE33=   -CWP*RMU*(CPART+CSUM5)   

      CLDISP(1)=CE11
      CLDISP(2)=CE33-CE11
      CLDISP(3)=CE22-CE11
      CLDISP(4)=CE31
      CLDISP(5)=CE12
      CLDISP(6)=CE23

      RETURN
  END SUBROUTINE DPTNKR

END MODULE DPTNSR1
