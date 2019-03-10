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

  SUBROUTINE DPTNCL(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      CLDISP(1)=-CWP/(1.D0-CWC*CWC)
      CLDISP(2)=-CWP-CLDISP(1)
      CLDISP(3)= 0.D0
      CLDISP(4)= 0.D0
      CLDISP(5)=(0.D0,-1.D0)*CWP*CWC/(1.D0-CWC*CWC)
      CLDISP(6)= 0.D0

      RETURN
  END SUBROUTINE DPTNCL

!     ****** COLLISIONAL COLD MODEL ******

  SUBROUTINE DPTNCC(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CWNU

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
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

  SUBROUTINE DPTNIM(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CPARA
      REAL(rkind):: RKPR,VTE

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS))
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))

      RKPR=DBLE(CKPR)
      IF(ABS(RKPR).LT.1.D-8) RKPR=1.D-8
      VTE=SQRT(plf(NS)%RTPR*AEE*1.D3/(PA(NS)*AMP))
      CPARA=plf(NS)%RN*1.D20*AEE*AEE &
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

  SUBROUTINE DPTNRM(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CWNU,CPARA

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
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

  SUBROUTINE DPTNWP(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CWP,CWC,CWNU,CPARA,CTPR,CTPP
      REAL(rkind):: WTPR,WTPP,WTPX

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF
      CWP=CWP/CWNU
      CWC=CWC/CWNU

      WTPR=plf(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plf(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
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

  SUBROUTINE DPTNUP(CW,CKPR,CKPP,NS,plf,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CK2,CGZ,CZ,CDZ,CWP,CFX,CFZ,CLP,CLT,CWNU
      REAL(rkind):: RT,VT2

      RT=(plf(NS)%RTPR+2.D0*plf(NS)%RTPP)/3.D0
      CK2=SQRT(CKPR**2+CKPP**2)
      VT2=RT*AEE*1.D3/(AMP*PA(NS))
      CGZ=CW/SQRT(2.D0*CK2*VT2)
      CALL DSPFN(CGZ,CZ,CDZ)

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
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

  SUBROUTINE DPTNHP(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CALAM(0:2)
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NC

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF

      WTPR=plf(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plf(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
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

  SUBROUTINE DPTNKL(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind):: CALAM(0:3)
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NC

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF

      WTPR=plf(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plf(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      CALAM(0)=1.D0-      CBETA+0.750D0*CBETA*CBETA
      CALAM(1)=     0.5D0*CBETA-0.500D0*CBETA*CBETA
      CALAM(2)=                 0.125D0*CBETA*CBETA
      CALAM(3)=0.D0

      DO NC=-2,2
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

  SUBROUTINE DPTNKS(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE libbes,ONLY: BESEIN
      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind),ALLOCATABLE:: CALAM(:)
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NCMIN,NCMAX,NHMAX,NC,IERR,NH

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NCMIN=NDISP1(NS)
      NCMAX=NDISP2(NS)
      NHMAX=MAX(ABS(NCMIN),ABS(NCMAX))+1
      ALLOCATE(CALAM(0:NHMAX))

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF

      WTPR=plf(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plf(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      DO NH=0,NHMAX
         CALAM(NH)=BESEIN(NH,REAL(CBETA))
      END DO

!      CALL LAMBDA(NHMAX,CBETA,CALAM,IERR)
!      IF(IERR.EQ.1) WRITE(6,*) 'XX LAMBDA: N out of range'
!      IF(IERR.EQ.2) WRITE(6,*) 'XX LAMBDA: CBETA out of range'

      DO NC=NCMIN,NCMAX
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
!         CLDISP(4)=CLDISP(4)-CWP*   CK*CBLAM*CFF
         CLDISP(4)=0.D0
         CLDISP(5)=CLDISP(5)+CWP*CI*NC*CDLAM*CFW
!         CLDISP(6)=CLDISP(6)+CWP*CI*CK*CDLAM*CFF
         CLDISP(6)=0.D0
      ENDDO
      RETURN
    END SUBROUTINE DPTNKS

!     ****** KINETIC MODEL WITH FLR ******

    SUBROUTINE DPTNKP(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE libbes,ONLY: BESEIN
      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind),ALLOCATABLE:: CALAM(:)
      COMPLEX(rkind):: CWP,CWC,CBETA,CPR,CGZ,CZ,CDZ,CK,CFW,CFF,CFG
      COMPLEX(rkind):: CLAM,CLAMM,CLAMP,CDLAM,CBLAM,CWNU
      REAL(rkind):: WTPR,WTPP,WTPX
      INTEGER:: I,NCMIN,NCMAX,NHMAX,NC,IERR,NH

      DO I=1,6
         CLDISP(I)=0.D0
      ENDDO

      NCMIN=NDISP1(NS)
      NCMAX=NDISP2(NS)
      NHMAX=MAX(ABS(NCMIN),ABS(NCMAX))+1
      ALLOCATE(CALAM(0:NHMAX))

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF

      WTPR=plf(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plf(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      WTPX=SQRT(WTPR/WTPP)
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)
      DO NH=0,NHMAX
         CALAM(NH)=BESEIN(NH,REAL(CBETA))
      END DO
!      CALL LAMBDA(NHMAX,CBETA,CALAM,IERR)
!      IF(IERR.EQ.1) WRITE(6,*) 'XX LAMBDA: N out of range'
!      IF(IERR.EQ.2) WRITE(6,*) 'XX LAMBDA: CBETA out of range'

      DO NC=NCMIN,NCMAX

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

  SUBROUTINE DPTNKR(CW,CKPR,CKPP,NS,mag,plf,CLDISP)

      USE libdsp,ONLY: DSPFN
      USE dpcomm
      USE plprof
      IMPLICIT NONE
      COMPLEX(rkind),INTENT(IN):: CW,CKPR,CKPP
      INTEGER,INTENT(IN):: NS
      TYPE(pl_mag_type),INTENT(IN):: mag
      TYPE(pl_plf_type),DIMENSION(nsmax),INTENT(IN):: plf
      COMPLEX(rkind),INTENT(OUT):: CLDISP(6)
      COMPLEX(rkind),ALLOCATABLE:: CALAM(:)
      COMPLEX(rkind):: CWP,CWC,CBETA,CNPP,CNPR,CZ,CWNU
      COMPLEX(rkind):: CSUM1,CSUM2,CSUM3,CSUM4,CSUM5
      COMPLEX(rkind):: CN1,CN2,CN3,CF1,CF2,CF3,CF4,CF5,CF6,CF7,CPART
      COMPLEX(rkind):: CE11,CE12,CE22,CE23,CE31,CE33
      REAL(rkind):: DELZ,WTPR,WTPP,RMU,DKAI
      INTEGER:: NCMIN,NCMAX,NC,N,NSIG

      DELZ=1.D-6

      NCMIN=NDISP1(NS)
      NCMAX=NDISP2(NS)

      CWP=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE/(EPS0*AMP*PA(NS)*CW*CW)
      CWC=ABS(mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW))

      IF(plf(NS)%RNUC.EQ.0.D0) THEN
         CWNU=DCMPLX(1.D0,plf(NS)%RZCL)
      ELSE
         CWNU=1.D0+CI*plf(NS)%RNUC/CW
      END IF

      WTPR=plf(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
      WTPP=plf(NS)%RTPP*1.D3*AEE/(AMP*PA(NS))
      CBETA=CKPP*CKPP*WTPP/(CWC*CWC*CW*CW)

      RMU=VC**2/WTPR
      CNPP=CKPP*VC/CW
      CNPR=CKPR*VC/CW

      CSUM1=(0.D0,0.D0)
      CSUM2=(0.D0,0.D0)
      CSUM3=(0.D0,0.D0)
      CSUM4=(0.D0,0.D0)
      CSUM5=(0.D0,0.D0)

      DO NC=NCMIN,NCMAX
         IF(NC.NE.0)THEN
            N=ABS(NC)
            NSIG=NC/N
            CZ=NC*RMU*CWC
            DKAI=DKAIJO(N)
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

!     ****** CALCULATE ! ******

  FUNCTION DKAIJO(N)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: CI,PI
      IMPLICIT NONE
      INTEGER,INTENT(IN):: N
      REAL(rkind):: DKAIJO
      INTEGER,PARAMETER::  &
          K(0:10)=(/1,1,2,6,24,120,720,5040,40320,362880,3628800/)
      REAL(rkind):: D
      INTEGER:: I

      IF(N.LT.0) THEN
         WRITE(6,*) 'XX DKAIJO: WRONG ARGUMENT: ',N
      ELSEIF(N.LE.10) THEN
         DKAIJO=DBLE(K(N))
      ELSE
         D=DBLE(K(10))
         DO I=11,N
            D=D*DBLE(I)
         ENDDO
         DKAIJO=D
      ENDIF
      RETURN
  END FUNCTION DKAIJO

!     ****** CALCULATE F ******

  FUNCTION CFQ(Q,CZ,CNPR,RMU)

      USE bpsd_kinds,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q,RMU
      COMPLEX(rkind),INTENT(IN):: CZ,CNPR
      COMPLEX(rkind):: CFQ
      COMPLEX(rkind):: CFQ0

      CFQ0=CFQZ(Q,RMU-CZ)
      CFQ=CFQ0 &
         +RMU*CNPR**2/2.D0*(   CFQZ(Q-1,RMU-CZ) &
                            -2*CFQ0 &
                            +  CFQZ(Q+1,RMU-CZ))
      RETURN
  END FUNCTION CFQ
!     
!     ****** CALCULATE SHKAROFSKY ******
!           *** WITH Z-FUNCTION ***

      FUNCTION CFQZ(Q,CZ)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: CI,PI
      USE libdsp,ONLY: DSPFN
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q
      COMPLEX(rkind),INTENT(IN):: CZ
      COMPLEX(rkind):: CFQZ,CZ2,CDZ2
      COMPLEX(rkind):: CSUM,CTERM,CGZ
      INTEGER:: NUMAX,NU,TEMP

      IF(ABS(CZ).GT.15.D0) THEN
         NUMAX=20
         CSUM=0.D0
         DO NU=0,NUMAX
            CTERM=-DGAMM(Q+NU)/(-CZ)**(NU+1)
            CSUM=CSUM+CTERM
            IF(ABS(CTERM).LE.1.D-12) GOTO 100
         ENDDO
  100    CONTINUE
         TEMP=DBLE(SQRT(CZ))
         IF(ABS(TEMP).LT.1.D-12) THEN
            CSUM=CSUM-  CI*PI*(-CZ)**(Q-1)*EXP(CZ)
         ELSEIF(TEMP.LT.0.D0) THEN
            CSUM=CSUM-2*CI*PI*(-CZ)**(Q-1)*EXP(CZ)
         ENDIF
         CFQZ=CSUM/DGAMM(Q)
      ELSE
         NUMAX=NINT(Q-3.D0/2.D0)
         CSUM=(0.D0,0.D0)
         DO NU=0,NUMAX 
            CTERM=(-CZ)**NU*DGAMM(Q-1-NU)
            CSUM=CSUM+CTERM
         ENDDO
         CGZ=CI*SQRT(CZ)
         CALL DSPFN(CGZ,CZ2,CDZ2)
         CSUM=CSUM+SQRT(PI)*(-CZ)**NUMAX*(CI*SQRT(CZ)*CZ2)
         CFQZ=CSUM/DGAMM(Q)
      ENDIF
      RETURN
  END FUNCTION CFQZ
!     
!     ****** CALCULATE SHKAROFSKY ******
!           *** WITH Z-FUNCTION ***

  FUNCTION CFQZ_Z(Q,CZ)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: CI,PI
      USE libdsp
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q
      COMPLEX(rkind),INTENT(IN):: CZ
      COMPLEX(rkind):: CFQZ_Z
      COMPLEX(rkind):: CSUM,CTERM,CGZ,CZ2,CDZ2
      INTEGER:: NUMAX,NU

      NUMAX=NINT(Q-3.D0/2.D0)
      CSUM=(0.D0,0.D0)
      DO NU=0,NUMAX 
         CTERM=(-CZ)**NU*DGAMM(Q-1-NU)
         CSUM=CSUM+CTERM
      ENDDO
      CGZ=CI*SQRT(CZ)
      CALL DSPFN(CGZ,CZ2,CDZ2)
      CSUM=CSUM+SQRT(PI)*(-CZ)**NUMAX*(CI*SQRT(CZ)*CZ2)
      CFQZ_Z=CSUM/DGAMM(Q)
      RETURN
  END FUNCTION CFQZ_Z
!     
!     ****** CALCULATE SHKAROFSKY ******
!     *** WITH ASYMPTOTIC EXPANSION ***

  FUNCTION CFQZ_EXP(Q,CZ)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: PI
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q
      COMPLEX(rkind),INTENT(IN):: CZ
      COMPLEX(rkind):: CFQZ_EXP
      COMPLEX(rkind):: CFQZ,CSUM
      REAL(rkind):: RREZ,RIMZ,SIG
      INTEGER:: NU

      CFQZ=(0.D0,0.D0)

         RREZ=DBLE(CZ)
         IF(RREZ.GE.-(11.D0+(Q-5.D0/2.D0)*1.35D0).AND.RREZ.LE.11.D0)THEN
            DO NU=0,50
               CSUM=(-CZ)**NU*DGAMM(Q-1-NU)
               CFQZ=CFQZ+CSUM
            ENDDO
            CFQZ=(CFQZ-PI*(-CZ)**(Q-3.D0/2.D0)*SQRT(CZ)*EXP(CZ)) &
                 /DGAMM(Q)
         ELSE
            RIMZ=DIMAG(CZ)
            IF(RIMZ.EQ.0.D0.AND.RREZ.LT.0.D0)THEN
               SIG=1.D0
            ELSE
               SIG=0.D0
            ENDIF
            DO NU=0,10
               CSUM=DGAMM(Q+NU)*((-CZ)**(-1-NU))
               CFQZ=CFQZ+CSUM
            ENDDO
            CFQZ=(CFQZ+(CFQZ+DGAMM(Q+(NU+1))*(-CZ)**(-1-NU-1)))/2.D0
            CFQZ_EXP &
                 =-(CFQZ-PI*SIG*(-CZ)**(Q-3.D0/2.D0)*SQRT(CZ)*EXP(CZ)) &
                 /DGAMM(Q)
         ENDIF
      RETURN
  END FUNCTION CFQZ_EXP

!****************************************************
!                   gamma function                  *
!                in double precision                *
!      COPYRIGHT : M.Mori  JUNE 30 1989  V.1        *
!****************************************************

  FUNCTION DGAMM(X)

      USE bpsd_kinds,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X
      REAL(rkind):: DGAMM
      INTEGER,PARAMETER:: IN=19
      REAL(rkind),PARAMETER:: &
           C(0:19)=(/  1.0D0, &
                      -0.4227843350984671D0, &
                      -0.2330937364217867D0, &
                       0.1910911013876915D0, &
                      -0.2455249000540002D-1, &
                      -0.1764524455014432D-1, &
                       0.8023273022267347D-2, &
                      -0.8043297756042470D-3, &
                      -0.3608378162548D-3, &
                       0.1455961421399D-3, &
                      -0.175458597517D-4, &
                      -0.25889950224D-5, &
                       0.13385015466D-5, &
                      -0.2054743152D-6, &
                      -0.1595268D-9, &
                       0.62756218D-8, &
                      -0.12736143D-8, &
                       0.923397D-10, &
                       0.120028D-10, &
                      -0.42202D-11 /)
      REAL(rkind):: XX,A,FCTR,Z,Y
      INTEGER:: M,MG,I

      IF (X .GT. 57.0D0) GO TO 901

      XX = X
      IF (XX .LE. 1.5D0) THEN
        IF (XX .GE. 0.5D0) THEN
          A = XX - 1.0D0
          FCTR = 1.0D0
        ELSE
          M = INT(XX)
          A = XX - M
          IF (A .EQ. 0.0D0) THEN
            GO TO 902
          ELSE IF (A .GE. -0.5D0) THEN
            MG = IABS(M) + 1
          ELSE
            MG = IABS(M) + 2
            A = A + 1.0D0
          END IF
          Z = 1.0D0
          DO I = 1, MG
            Z = Z * XX
            XX = XX + 1.0D0
          ENDDO
          FCTR = 1.0D0 / Z
        END IF

      ELSE
        M = INT (XX)
        A = XX - M
        IF (A .LE. 0.5D0) THEN
          MG = M - 1
        ELSE
          MG = M
          A = A - 1.0D0
        END IF
        Z = 1.0D0
        DO I = 1, MG
          Z = Z * (XX - 1.0D0)
          XX = XX - 1.0D0
       ENDDO
        FCTR = Z
      END IF

      Y = C(IN)
      DO I = IN - 1, 0, -1
        Y = C(I) + A * Y
      ENDDO

      DGAMM = FCTR / ((1.0D0 + A) * Y)
      RETURN

  901 CONTINUE
      WRITE (6,2001) X
 2001 FORMAT (' (FUNC.DGAMM) X(=',D23.16,')', &
              ' must be smaller than 57.0')
      DGAMM = 1.0D75
      RETURN

  902 CONTINUE
      WRITE (6,2002) X
 2002 FORMAT (' (FUNC.DGAMM) invalid argument', &
              ' X =',D23.16)
      DGAMM = 1.0D75
      RETURN

  END FUNCTION DGAMM
END MODULE DPTNSR1
