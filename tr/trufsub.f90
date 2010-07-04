!     $Id$
!     ***********************************************************

!           1D UFILE READER

!     ***********************************************************

!     input:

!     KFID        : Variable Name
!     DT          : Time Step Width
!     NTLMAX      : Maximum Time Step for TASK/TR
!     TLMAX       : Maximum Time for TASK/TR
!     ICK         : Check Indicator of UFILE Consistency
!     MDLXP       : Select UFILE or MDSplus

!     output:

!     TL(NTUM)    : Total Time Data (The Number of DT)
!     F1(NTUM)    : Functional Values
!     IERR        : Error Indicator

!     *****************************************************

      SUBROUTINE UF1D(KFID,KUFDEV,KUFDCG,DT,TL,F1,NTLMAX,NTXMAX,TLMAX,ICK,MDLXP,IERR)

      USE TRCOMM, ONLY : NTUM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      CHARACTER(LEN=10),INTENT(IN)   :: KFID
      CHARACTER(LEN=80),INTENT(INOUT):: KUFDEV,KUFDCG
      REAL(8)          ,INTENT(IN)   :: DT
      REAL(8)          ,INTENT(INOUT):: TLMAX
      REAL(8),DIMENSION(NTUM),INTENT(OUT):: TL, F1
      INTEGER(4)       ,INTENT(IN)   :: NTXMAX, MDLXP
      INTEGER(4)       ,INTENT(INOUT):: NTLMAX, ICK
      INTEGER(4)       ,INTENT(OUT)  :: IERR
      INTEGER(4):: MDCHK


      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TL,F1,NTXMAX,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TL,F1,NTXMAX,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,'(A13,A10,A14)') '## UF1D: NO "',KFID,'" FILE EXISTS.'
         RETURN
      ENDIF
!
!     Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
!
      TL(2:NTXMAX)=(TL(2:NTXMAX)-TL(1))
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(ICK.EQ.2.AND.TLMAX.EQ.0.D0) GOTO 1000
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX UF1D:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
            STOP
         ENDIF
      ENDIF
 1000 CONTINUE
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
      IF(ICK.NE.2) ICK=1

      RETURN
      END SUBROUTINE UF1D

!     *** ROUTINE FOR TEXT ***

      SUBROUTINE UF1DT(KFID,KUFDEV,KUFDCG,NTS,FOUT,AMP,MDLXP,IERR)

      USE TRCOMM, ONLY : NTUM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      CHARACTER(LEN=10),INTENT(IN)   :: KFID
      CHARACTER(LEN=80),INTENT(INOUT):: KUFDEV,KUFDCG
      INTEGER(4)       ,INTENT(IN)   :: NTS, MDLXP
      INTEGER(4)       ,INTENT(OUT)  :: IERR
      REAL(8)          ,INTENT(IN)   :: AMP
      REAL(8)          ,INTENT(OUT)  :: FOUT
      INTEGER(4)             :: NTXMAX, MDCHK
      REAL(8),DIMENSION(NTUM):: F1, TL


      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TL,F1,NTXMAX,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TL,F1,NTXMAX,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,'(A14,A10,A14)') '## UF1DG: NO "',KFID,'" FILE EXISTS.'
         RETURN
      ENDIF

      IF(IERR.EQ.1) THEN
         FOUT=0.D0
      ELSE
         FOUT=F1(NTS)*AMP
      ENDIF

      RETURN
      END SUBROUTINE UF1DT

!     *** ROUTINE FOR GRAPHIC ***

      SUBROUTINE UF1DG(KFID,KUFDEV,KUFDCG,GTL,TL,FOUT,AMP,NINMAX,MDLXP,IERR)

      USE TRCOMM, ONLY : NTUM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      CHARACTER(LEN=10),INTENT(IN)   :: KFID
      CHARACTER(LEN=80),INTENT(INOUT):: KUFDEV, KUFDCG
      REAL(8),          INTENT(IN)   :: AMP
      REAL(8),DIMENSION(NTUM),INTENT(INOUT):: TL, FOUT
      REAL(4),DIMENSION(NTUM),INTENT(IN)   :: GTL
      INTEGER(4),             INTENT(IN)   :: NINMAX, MDLXP
      INTEGER(4),             INTENT(OUT)  :: IERR
      INTEGER(4):: NTXMAX, MDCHK, NIN, IERRL
      REAL(8)   :: TLN, F0
      REAL(8),DIMENSION(NTUM):: F1(NTUM)


      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TL,F1,NTXMAX,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TL,F1,NTXMAX,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,'(A14,A10,A14)') '## UF1DG: NO "',KFID,'" FILE EXISTS.'
         RETURN
      ENDIF
!
!     Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
!
      TL(2:NTXMAX)=(TL(2:NTXMAX)-TL(1))
      TL(1)=0.D0
      DO NIN=1,NINMAX
         TLN=DBLE(GTL(NIN))
         CALL TIMESPL(TLN,F0,TL,F1,NTXMAX,NTUM,IERRL)
!         IF(IERRL.NE.0) WRITE(6,600) "XX UF1DG: TIMESPL ",KFID,": IERR=",IERR
         FOUT(NIN)=F0
      ENDDO
!
      IF(IERR.EQ.1) THEN
         FOUT(1:NINMAX)=0.D0
      ELSE
         FOUT(1:NINMAX)=FOUT(1:NINMAX)*AMP
      ENDIF

 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END SUBROUTINE UF1DG

!     ***********************************************************

!           2D UFILE READER

!     ***********************************************************

!     *** FOR STEADY STATE SIMULATION ***

!     input:

!     KFID        : Variable Name
!     DR          : Radial Step Width
!     AMP         : Amplitude Factor of FOUT
!     NRMAX       : Radial Node Number
!     NSW         : Mesh Selector (0:RM, 1:RG)
!     ID          : Boundary Condition in the center and/or edge for Spline
!     MDLXP       : Select UFILE or MDSplus

!     output:

!     TL(NTUM)    : Total Time Data (The Number of DT)
!     FOUT(NRMP)  : Functional Values
!     IERR        : Error Indicator

!     *****************************************************

      SUBROUTINE UF2DS(KFID,KUFDEV,KUFDCG,DR,TL,FOUT,AMP,NRMAX,NSW,ID,MDLXP,IERR)

      USE TRCOMM, ONLY : NRMP, NRUM, NTUM
      USE TRCOM1, ONLY : KDIRX, NTXMAX
      IMPLICIT NONE
      CHARACTER(LEN=10),INTENT(IN):: KFID
      CHARACTER(LEN=80),INTENT(INOUT):: KUFDEV, KUFDCG
      REAL(8),DIMENSION(NTUM),INTENT(OUT):: TL
      REAL(8),DIMENSION(NRMP),INTENT(OUT):: FOUT
      REAL(8),INTENT(IN) :: DR,  AMP
      INTEGER(4),INTENT(IN):: NRMAX, NSW, ID, MDLXP
      INTEGER(4),INTENT(OUT):: IERR

      INTEGER(4):: IERRP, IERRS, MDCHK, NRL, NRLMAX, NTX
      REAL(8)   :: F0, FCTR, RSL
      REAL(8),DIMENSION(NRUM):: RL
      REAL(8),DIMENSION(4,NRUM):: U
      REAL(8),DIMENSION(NTUM,NRUM):: F2


      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,NRLMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RL,TL,F2,NRLMAX,NTXMAX,IERR)
      ENDIF

      IF(KUFDEV.EQ.'jet'.AND.(KUFDCG.EQ.'35156'.OR.KUFDCG.EQ.'35171') &
     &     .AND.KFID.EQ.'GRHO1') THEN
         DO NTX=1,NTXMAX
            F2(NTX,2)=FCTR(RL(3),RL(4),F2(NTX,3),F2(NTX,4))
            F2(NTX,1)=FCTR(RL(2),RL(3),F2(NTX,2),F2(NTX,3))
         ENDDO
      ENDIF
      IERRP=IERR
      CALL PRETREAT0(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,ID,IERRP)
!
!     Error check in case of a single radial point
!
      IF(IERRP.NE.0.AND.NRLMAX.LE.1) THEN
         FOUT(1:NRMAX)=0.D0
         IERR=1
         IF(IERRP.LT.0) IERR=IERRP
         RETURN
      ENDIF
!
!     Calculate values suitable for arbitrary radial mesh using spline
!
      DO NRL=1,NRMAX
         IF(NSW.EQ.0) THEN
            RSL=(DBLE(NRL)-0.5D0)*DR
         ELSEIF(NSW.EQ.1) THEN
            RSL= DBLE(NRL)       *DR
         ENDIF
         CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERRS)
         IF(IERRS.NE.0) WRITE(6,600) "XX UF2DS: SPL1DF ",KFID,": IERRS=",IERRS
         FOUT(NRL)=F0*AMP
      ENDDO

 600  FORMAT(' ',A17,A10,A7,I2)
      RETURN
      END SUBROUTINE UF2DS

!     *****************************************************

!     input:

!     KFID        : Variable name
!     DR          : Radial Step Width
!     AMP         : Amplitude Factor of FOUT
!     RHOA        : Normalized Radius at Arbitrary Surface
!     NRAMAX      : Radial Node Number Corresponding to RHOA
!     NRMAX       : Radial Node Number
!     MDLXP       : Select UFILE or MDSplus

!     output:

!     PV          : Surface(Peripheral) Functional Value
!     PVS         : Surface(Peripheral) Functional Value
!                   Corresponding to RHOA
!     FOUT(NRMP)  : Functional Values
!     IERR        : Error Indicator

!     *****************************************************

      SUBROUTINE UF2DSP(KFID,KUFDEV,KUFDCG,DR,PV,PVA,FOUT,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)

      USE TRCOMM, ONLY : NRMP, NRUM, NTUM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      CHARACTER(LEN=10), INTENT(IN):: KFID
      CHARACTER(LEN=80), INTENT(INOUT):: KUFDEV, KUFDCG
      REAL(8) ,          INTENT(IN)  :: DR, AMP, RHOA
      REAL(8) ,          INTENT(OUT)  :: PV, PVA
      REAL(8),DIMENSION(NRMP), INTENT(OUT)::FOUT
      INTEGER(4),        INTENT(IN):: NRAMAX, NRMAX, MDLXP
      INTEGER(4),        INTENT(OUT):: IERR
      INTEGER(4):: IERRP, MDCHK, NRL, NRLMAX, NTXMAX
      REAL(8)   :: F0, RGN, RMN
      REAL(8),DIMENSION(NRUM):: RL
      REAL(8),DIMENSION(NTUM):: TL
      REAL(8),DIMENSION(NTUM,NRUM):: F2
      REAL(8),DIMENSION(4,NRUM):: U


      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,NRLMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RL,TL,F2,NRLMAX,NTXMAX,IERR)
      ENDIF
      IERRP=IERR
      CALL PRETREAT0(KFID,RL,TL,F2,U,NRLMAX,NTXMAX,1,IERRP)
!
!     Error check in case of a single radial point
!
      IF(IERRP.NE.0.OR.NRLMAX.LE.1) THEN
         FOUT(1:NRMAX)=0.D0
         RETURN
!200607   IERR=1
!200607   IF(IERR.LT.0) IERR=IERRP
      ENDIF
!
!     Calculate values suitable for arbitrary radial mesh using spline
!
      DO NRL=1,NRMAX
         RMN=(DBLE(NRL)-0.5D0)*DR
         CALL SPL1DF(RMN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0) WRITE(6,600) "XX UF2DSP: SPL1DF ",KFID,": IERR=",IERR
         FOUT(NRL)=F0*AMP
      ENDDO
!
!     Edge values
!
      RGN=DBLE(NRMAX)*DR
      CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
      IF(IERR.NE.0) WRITE(6,600) "XX UF2DSP: SPL1DF ",KFID,": IERR=",IERR
      PV=F0*AMP
      IF(RHOA.NE.1.D0) THEN
         RGN=DBLE(NRAMAX)*DR
         CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0) WRITE(6,600) "XX UF2DSP: SPL1DF ",KFID,": IERR=",IERR
         PVA=F0*AMP
      ENDIF

 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END SUBROUTINE UF2DSP

!     *** FOR TIME EVOLUTION SIMULATION ***

!     input:

!     KFID        : Variable name
!     DR          : Radial Step Width
!     DT          : Time Step Width
!     AMP         : Amplitude Factor of FOUT
!     NRMAX       : Radial Node Number
!     TLMAX       : Maximum Time
!     NSW         : Mesh Selector (0:RM, 1:RG)
!     ID          : Boundary Condition in the center and/or edge for Spline
!     ICK         : Check Indicator
!     MDLXP       : Select UFILE or MDSplus

!     output:

!     TL(NTUM)    : Total Time Data (The Number of DT)
!     FOUT(NRMP)  : Functional Values
!     NTAMAX      : Maximum Time Step Number Corresponding to RHOA
!     IERR        : Error Indicator

!     *****************************************************

      SUBROUTINE UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TL,FOUT,AMP,NTLMAX,NTXMAX,NRMAX,TLMAX,NSW,ID,ICK,MDLXP,IERR)

      USE TRCOMM, ONLY : NRMP, NRUM, NTUM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      CHARACTER(LEN=10),           INTENT(IN):: KFID
      CHARACTER(LEN=80),           INTENT(INOUT):: KUFDEV, KUFDCG
      REAL(8),DIMENSION(NTUM),     INTENT(OUT):: TL
      REAL(8),DIMENSION(NTUM,NRMP),INTENT(OUT):: FOUT
      REAL(8),                     INTENT(IN):: DR, DT, AMP
      REAL(8),                     INTENT(INOUT):: TLMAX
      INTEGER(4),                  INTENT(IN):: NTXMAX, NRMAX, NSW, ID, MDLXP
      INTEGER(4),                  INTENT(INOUT):: ICK
      INTEGER(4),                  INTENT(OUT):: NTLMAX, IERR
      INTEGER(4):: MDCHK, NRL, NRLMAX, NTX
      REAL(8)   :: F0, RSL
      REAL(8),DIMENSION(NRUM):: RL
      REAL(8),DIMENSION(NTUM,NRUM):: F2
      REAL(8),DIMENSION(4,NRUM):: U
      REAL(8),DIMENSION(NRUM):: DERIV, TMP0

      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,NRLMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RL,TL,F2,NRLMAX,NTXMAX,IERR)
      ENDIF
!
!     Error check in case of a single radial point
!
      IF(IERR.NE.0.AND.NRLMAX.LE.1) THEN
         FOUT(1:NTUM,1:NRMP)=0.D0
         IF(NRLMAX.LE.1) IERR=1
         RETURN
      ENDIF
      DERIV(1)=0.D0
      DERIV(NRLMAX)=0.D0
!
!     Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
!
      TL(2:NTXMAX)=(TL(2:NTXMAX)-TL(1))
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(ICK.EQ.2.AND.TLMAX.EQ.0.D0) GOTO 1000
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX UF2DT:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
            STOP
         ENDIF
      ENDIF
 1000 CONTINUE
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
      IF(ICK.NE.2) ICK=1
!
!     Calculate values suitable for arbitrary radial mesh using spline
!
      DO NTX=1,NTXMAX
         TMP0(1:NRLMAX)=F2(NTX,1:NRLMAX)
         CALL SPL1D(RL,TMP0,DERIV,U,NRLMAX,ID,IERR)
         DO NRL=1,NRMAX
            IF(NSW.EQ.0) THEN
               RSL=(DBLE(NRL)-0.5D0)*DR
            ELSEIF(NSW.EQ.1) THEN
               RSL= DBLE(NRL)       *DR
            ENDIF
            CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERR)
            IF(IERR.NE.0) WRITE(6,600) "XX TRFILE: SPL1DF ",KFID,": IERR=",IERR
            FOUT(NTX,NRL)=F0*AMP
         ENDDO
      ENDDO

 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END SUBROUTINE UF2DT

!     *** FOR THE VARIABLE WE'D LIKE TO USE PERIPHERAL VALUES ***

      SUBROUTINE UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TL,FOUT,AMP, &
     &                  NTLMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TLMAX,NSW,ICK,MDLXP,IERR)

      USE TRCOMM, ONLY : NRMP, NRUM, NTUM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      CHARACTER(LEN=10),           INTENT(IN)   :: KFID
      CHARACTER(LEN=80),           INTENT(INOUT):: KUFDEV,KUFDCG
      REAL(8),                     INTENT(IN)   :: DR, DT, AMP, RHOA
      REAL(8),                     INTENT(INOUT):: TLMAX
      REAL(8),DIMENSION(NTUM),     INTENT(OUT)  :: PV,PVA,TL
      REAL(8),DIMENSION(NTUM,NRMP),INTENT(OUT)  :: FOUT
      INTEGER(4),                  INTENT(IN)   :: NTXMAX, NRAMAX, NRMAX, NSW, MDLXP
      INTEGER(4),                  INTENT(INOUT):: ICK
      INTEGER(4),                  INTENT(OUT)  :: NTLMAX, IERR
      INTEGER(4):: MDCHK, NRL, NRLMAX, NTX
      REAL(8)   :: RSL, F0, RGN
      REAL(8),DIMENSION(NRUM)     :: DERIV, RL, TMP0
      REAL(8),DIMENSION(4,NRUM)   :: U
      REAL(8),DIMENSION(NTUM,NRUM):: F2


      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,NRLMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RL,TL,F2,NRLMAX,NTXMAX,IERR)
      ENDIF
!
!     All the variables which this subroutine should handle
!        have a derivative of zero on-axis.
!
      DERIV(1)=0.D0
!
!     Error check in case of a single radial point
!
      IF(IERR.NE.0.AND.NRLMAX.LE.1) THEN
         FOUT(1:NTUM,1:NRMP)=0.D0
         IF(NRLMAX.LE.1) IERR=1
         RETURN
      ENDIF
!
!     Time mesh normalization (from t=a to t=b -> t=0 to t=b-a)
!
      TL(2:NTXMAX)=(TL(2:NTXMAX)-TL(1))
      TL(1)=0.D0
      IF(ICK.NE.0) THEN
         IF(ICK.EQ.2.AND.TLMAX.EQ.0.D0) GOTO 1000
         IF(TLMAX.NE.TL(NTXMAX)) THEN
            WRITE(6,*) 'XX UF2DTP:',KFID,'UFILE HAS AN ERROR!'
            WRITE(6,*) TLMAX,TL(NTXMAX)
!            STOP
         ENDIF
      ENDIF
 1000 CONTINUE
      TLMAX=TL(NTXMAX)
      NTLMAX=INT(DINT(TL(NTXMAX)*1.D2)*1.D-2/DT)
      IF(ICK.NE.2) ICK=1
!
!     Calculate values suitable for arbitrary radial mesh using spline
!
      DO NTX=1,NTXMAX
         TMP0(1:NRLMAX)=F2(NTX,1:NRLMAX)
         CALL SPL1D(RL,TMP0,DERIV,U,NRLMAX,1,IERR)
         DO NRL=1,NRMAX
            IF(NSW.EQ.0) THEN
               RSL=(DBLE(NRL)-0.5D0)*DR
            ELSEIF(NSW.EQ.1) THEN
               RSL= DBLE(NRL)       *DR
            ENDIF
            CALL SPL1DF(RSL,F0,RL,U,NRLMAX,IERR)
            IF(IERR.NE.0) WRITE(6,600) "XX UF2DTP: SPL1DF ",KFID,": IERR=",IERR
            FOUT(NTX,NRL)=F0*AMP
         ENDDO
!
!     Edge values
!
         RGN=DBLE(NRMAX)*DR
         CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
         IF(IERR.NE.0) WRITE(6,600) "XX UF2DTP: SPL1DF ",KFID,": IERR=",IERR
         PV(NTX)=F0*AMP
         IF(RHOA.NE.1.D0) THEN
            RGN=DBLE(NRAMAX)*DR
            CALL SPL1DF(RGN,F0,RL,U,NRLMAX,IERR)
            IF(IERR.NE.0) WRITE(6,600) "XX UF2DTP: SPL1DF ",KFID,": IERR=",IERR
            PVA(NTX)=F0*AMP
         ENDIF
      ENDDO

 600  FORMAT(' ',A18,A10,A7,I2)
      RETURN
      END SUBROUTINE UF2DTP
!     ***********************************************************

!           PRETREATMENT SUBROUTINE FOR UFILE INTERFACE

!     ***********************************************************

!     This subroutine is only used if MDLUF=2 and any subroutines
!     do not call this except UF2DS and UF2DSP.
!     In this subroutine, one can choose arbitrary slice time if one
!     did not choose it ahead of time. If one choose the time,
!     this checks its value consistent with the range of time.
!     One of the main function in this section is to make "Spline Array"
!     in order to interpolate various profiles radially.

!     IERR: negative value means "no data file".

      SUBROUTINE PRETREAT0(KFID,RL,TL,F2,U,NRFMAX,NTXMAX,ID,IERR)

      USE TRCOMM, ONLY : MDLJQ, NRUM, NTUM, TIME_INT
      IMPLICIT NONE
      CHARACTER(LEN=10)           ,INTENT(IN):: KFID
      REAL(8),DIMENSION(NRUM)     ,INTENT(IN) :: RL
      REAL(8),DIMENSION(NTUM)     ,INTENT(IN) :: TL
      REAL(8),DIMENSION(NTUM,NRUM),INTENT(IN) :: F2
      REAL(8),DIMENSION(4,NRUM)   ,INTENT(OUT):: U
      INTEGER(4),                  INTENT(IN) :: NRFMAX, NTXMAX, ID
      INTEGER(4),                  INTENT(INOUT):: IERR
      INTEGER(4):: NTSL, NTX, NTX_MIN
      REAL(8)   :: TL_MIN, TL_MIN_OLD
      REAL(8),DIMENSION(NRUM):: DERIV, TMP

      IF(IERR.EQ.1) THEN
         U(1:4,1:NRUM)=0.D0
         IF(KFID.EQ.'CURTOT'.AND.IERR.NE.0.AND.MDLJQ.EQ.0) MDLJQ=1
         IERR=-1
         RETURN
      ENDIF

      DERIV(1)=0.D0
      DERIV(NRFMAX)=0.D0

      NTSL=1
      IF(NTXMAX.NE.1) THEN
         IF(TIME_INT.LE.0.D0) THEN
 100        WRITE(6,500) 'INPUT ARBITRARY TIME:',TL(1),' -',TL(NTXMAX)
            READ(5,*,ERR=100) TIME_INT
            IF(TIME_INT.LT.TL(1).OR.TIME_INT.GT.TL(NTXMAX)) GOTO 100
            DO NTX=1,NTXMAX
               IF(ABS(TL(NTX)-TIME_INT).LE.1.D-5) THEN
                  NTSL=NTX
                  GOTO 1000
               ENDIF
            ENDDO

            TL_MIN=TL(NTXMAX)
            DO NTX=1,NTXMAX
               TL_MIN_OLD=TL_MIN
               TL_MIN=MIN(ABS(TL(NTX)-TIME_INT),TL_MIN)
               IF(TL_MIN_OLD.EQ.TL_MIN) THEN
                  NTX_MIN=NTX-1
                  GOTO 200
               ENDIF
            ENDDO
         ELSE
            IF(TIME_INT.GE.TL(1).AND.TIME_INT.LE.TL(NTXMAX)) THEN
               DO NTX=1,NTXMAX
                  IF(ABS(TL(NTX)-TIME_INT).LE.1.D-5) THEN
                     NTSL=NTX
                     GOTO 1000
                  ENDIF
               ENDDO

               TL_MIN=TL(NTXMAX)
               DO NTX=1,NTXMAX
                  TL_MIN_OLD=TL_MIN
                  TL_MIN=MIN(ABS(TL(NTX)-TIME_INT),TL_MIN)
                  IF(TL_MIN_OLD.EQ.TL_MIN) THEN
                     NTX_MIN=NTX-1
                     GOTO 200
                  ENDIF
               ENDDO
            ENDIF
            WRITE(6,*) '## DESIGNATED TIME_INT IS NOT IN THE RANGE.'
            WRITE(6,*) 'TIME_INT=',SNGL(TIME_INT),'RANGE=',SNGL(TL(1)),'-',SNGL(TL(NTXMAX))
         ENDIF

 200     WRITE(6,500) 'TIME_INT=',TIME_INT,' HAS BEEN REPLACED BY',TL(NTX_MIN)
         TIME_INT=TL(NTX_MIN)
         NTSL=NTX_MIN
 1000    CONTINUE
      ENDIF

      TMP(1:NRFMAX)=F2(NTSL,1:NRFMAX)
      CALL SPL1D(RL,TMP,DERIV,U,NRFMAX,ID,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX PRETREAT0: SPL1D',KFID,': IERR=',IERR

      RETURN
 500  FORMAT(' ',A,F9.5,A,F9.5)
      END SUBROUTINE PRETREAT0
