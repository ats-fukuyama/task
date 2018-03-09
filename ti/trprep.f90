! MODULE tiprep

MODULE tiprep

CONTAINS

  SUBROUTINE ti_prep

    USE ticomm
    IMPLICIT NONE
    INTEGER:: NS,NSA,NZ,NEQ,IERR

!   *** count NSAMAX: number of active particles ***

    NSA=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(-1,0,1,5,6)
          NSA=NSA+1
       CASE(10)
          DO NZ=NZMIN_NS(NS),NZMAX_NS(NS)
             NSA=NSA+1
          END DO
       CASE DEFAULT
          WRITE(6,'(A,2I5)') &
               'XX ti_prep: Undefined ID_NS(NS): ID_NS,NS:', &
               ID_NS(NS),NS
          STOP
       END SELECT
    END DO
    NSMAX=NSA

!   *** ALLOCATE array for NSAMAX and NRMAX ***

    CALL allocate_ticomm(IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX tiprep: allocate_ticomm ERROR: IERR=',IERR
       STOP
    END IF

!   *** Define NSA variables ***

    NSA=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(-1,0,1,5,6)
          NSA=NSA+1
          PMA(NSA)=PM(NS)
          PZA(NSA)=PZ(NS)   ! for ID_NS=5,6, PZ will be ealuated later
          NS_NSA(NSA)=NS
       CASE(10)
          DO NZ=NZMIN_NS(NS),NZMAX_NS(NS)
             NSA=NSA+1
             PMA(NSA)=PM(NS)
             PZA(NSA)=DBLE(NZ)
             NS_NSA(NSA)=NS
          END DO
       END SELECT
    END DO
    IF(NSA.NE.NSAMAX) THEN
       WRITE(6,*) 'XX ti_prep: INCONSISTENT NSAMAX'
       STOP
    END IF

!   *** Display NSA variables ***

    WRITE(6,*) 'NSA  ','NS   ','ID   ','PMA            ','PZA'
    DO NSA=1,NSAMAX
       WRITE(6,'(3I5,1P2E12.4)') &
            NSA,NS_NSA(NSA),ID_NS(NS_NSA(NSA)),PMA(NSA),PZA(NSA)
    END DO

!   *** Count NEQMAX: Size of equation ***

    NEQ=0
    IF(MODEL_EQB.EQ.1) NEQ=NEQ+1
    DO NSA=1,NSAMAX
       NS=NS_NSA(NSA)
       IF(ID_NS(NS).GE.0) THEN
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+1
       END IF
    END DO
    NEQMAX=NEQ
    WRITE(6,'(A,I5)') 'NEQMAX=',NEQMAX

!   *** ALLOCATE array for NEQMAX and NRMAX ***

    CALL allocate_neqmax(IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX tiprep: allocate_neqmax ERROR: IERR=',IERR
       STOP
    END IF

!   *** Define NEQ variables ***

    NEQ=0
    IF(MODEL_EQB.EQ.1) THEN
       NEQ=NEQ+1
       NSA_NEQ(NEQ)=0
       NV_NEQ(NEQ)=0
    END IF
    DO NSA=1,NSAMAX
       NS=NS_NSA(NSA)
       IF(ID_NS(NS).GE.0) THEN
          IF(MODEL_EQN.EQ.1) THEN
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=1
          END IF
          IF(MODEL_EQU.EQ.1) THEN
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=2
          END IF
          IF(MODEL_EQT.EQ.1) THEN
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=3
          END IF
       END IF
    END DO
    IF(NEQ.NE.NEQMAX) THEN
       WRITE(6,*) 'XX ti_prep: INCONSISTENT NEQMAX'
       STOP
    END IF

!   *** Display NEQ variables ***

    WRITE(6,*) 'NEQ  ','NS   ','NSA   ','NV  '
    DO NEQ=1,NEQMAX
       WRITE(6,'(4I5)') &
            NEQ,NS_NSA(NSA_NEQ(NEQ))),NSA_NEQ(NEQ),NV_NEQ(NEQ)
    END DO

    RETURN
  END SUBROUTINE ti_prep
END MODULE tiprep






    NSCL=0
    DO NSA=1,NSAMAX
       IF(PZA(NSA).NE.0.D0) THEN
          NSCL=NSCL+1
          NSCL_NSA(NSA)=NSCL
          NSA_NSCL(NSCL)=NSA
       END IF
    END DO
    NSCLMAX


!     ====== INTIALIZE ======

    NT=0
    TIME_INT=0.D0
    CNB=1.D0
    SUMPBM=0.D0
    IREAD=0

    RETURN    
  END SUBROUTINE trprep

!     ***********************************************************

!           MODEL SELECTOR

!     ***********************************************************

      SUBROUTINE TR_EQS_SELECT(INIT)

      USE TRCOMM, ONLY : AMP, AMZ, MDDIAG, MDLEOI, MDLEQ0, MDLEQB, MDLEQE, MDLEQN, MDLEQT, MDLEQU, MDLEQZ, MDLKAI, MDLUF, &
                         MDLWLD, MDNCLS, NEA, NEQM, NEQMAX, NEQMAXM, NNS, NREDGE, NRMAX, NSCMAX, NSLMAX, NSM, NSMAX,      &
                         NSNMAX, NSS, NST, NSTM, NSTMAX, NSV, NSZMAX, PA, PZ, RGFLS, RQFLS
      USE TRCOM1, ONLY : INS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INIT
      INTEGER(4) :: IND, INDH, INDHD, MDANOM, MDSLCT, NEQ, NEQ1, NEQI, NEQRMAX, NEQS, NEQT, NNSC, NNSMAX, NNSN, &
                    NS, NSSN, NSVN
      INTEGER(4), SAVE :: NSSMAX


!     If INS is zero, all particles designated by NSMAX are fully
!     calculated. If INS=1, we handle only electrons and bulk ions
!     IF INS=2, three species are employed in the simulation but
!     just one of them is calculated.

      NSCMAX=NSMAX+NSZMAX          ! the number of charged particles
      NSTMAX=NSMAX+NSZMAX+NSNMAX   ! the number of all particles

      IF(NSMAX.EQ.1.AND.MDLEOI.EQ.0) MDLEOI=1
      IF(INIT.EQ.0) THEN
         INS=0
      ELSE
         IF(NSSMAX.NE.NSMAX) THEN
            INS=0
         ENDIF
      ENDIF
      NSSMAX=NSMAX
      IF(MDLUF.NE.0) THEN
         CALL CHECK_IMPURITY(MDSLCT)
         IF(MDSLCT.EQ.0) THEN
            IF(NSMAX.EQ.1) THEN
               INS=1
               NSMAX=2
            ENDIF
         ELSE
            IF(NSMAX.EQ.1) INS=2
            NSMAX=3
!            PA(3)=12.D0
!            PZ(3)=6.D0
         ENDIF
      ENDIF
!     *** for NCLASS ***
      IF(NSMAX.EQ.1) THEN
         NSLMAX=2
      ELSE
         NSLMAX=NSMAX
      ENDIF
!     ***

!      IF(MDLEQT.EQ.0) THEN
!         NEQMAX=MDLEQB+(MDLEQN+MDLEQT+MDLEQU)*2+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
!      ELSEIF(MDLEQT.EQ.1) THEN
!         NEQMAX=MDLEQB+2+(MDLEQT+MDLEQU)*2+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
!      ENDIF
      IF(MDLEQT.EQ.0) THEN
         NEQMAX=MDLEQB+(MDLEQN+MDLEQT+MDLEQU)*NSMAX+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ELSEIF(MDLEQT.EQ.1) THEN
         NEQMAX=MDLEQB+NSMAX+(MDLEQT+MDLEQU)*NSMAX+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ENDIF

      IF(MDLEQN.EQ.0.AND.(MDLEQE.EQ.1.OR.MDLEQE.EQ.2)) THEN
         WRITE(6,*)  'XX TR_EQS_SELECT : MDLEQE can be 1 or 2 when MDLEQN is 1.'
         STOP
      ENDIF

      NSS(1:NEQMAXM)=-1  ! 0: EM, particle species
      NSV(1:NEQMAXM)=-1  ! 1: density, 2: temperature 3: toroidal flow
      NNS(1:NEQMAXM)=0   ! 
      NST(1:NEQMAXM)=0
      NEQ=0
      IF(MDLEQB.EQ.1) THEN
         NEQ=NEQ+1
         NSS(NEQ)=0
         NSV(NEQ)=0
      ENDIF
      IND=0
      INDH=0
      INDHD=0
      DO NS=1,NSM
         IF(MDLEQN.EQ.1.OR.MDLEQT.EQ.1) THEN
            IF(MDLEQN.EQ.1) IND=1
            CALL TR_TABLE(NS,NEQ,1,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
         IF(MDLEQT.EQ.1) THEN
            CALL TR_TABLE(NS,NEQ,2,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
         IF(MDLEQU.EQ.1) THEN
            CALL TR_TABLE(NS,NEQ,3,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
      ENDDO
      IF(MDLEQZ.EQ.1) THEN
         IF(NSZMAX.EQ.1) THEN
            NEQ=NEQ+1
            NSS(NEQ  )=5
            NSV(NEQ  )=1
         ELSEIF(NSZMAX.EQ.2) THEN
            NEQ=NEQ+2
            NSS(NEQ-1)=5
            NSV(NEQ-1)=1
            NSS(NEQ  )=6
            NSV(NEQ  )=1
         ENDIF
      ENDIF
      IF(MDLEQ0.EQ.1) THEN
         NEQ=NEQ+2
         NSS(NEQ-1)=7
         NSV(NEQ-1)=1
         NSS(NEQ  )=8
         NSV(NEQ  )=1
      ENDIF

      IF(INS.NE.0) THEN
         NEQI=0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.NE.MDLEOI.OR.NSVN.NE.2) THEN
               NEQI=NEQI+1
               NNS(NEQI)=NEQ
            ENDIF
         ENDDO
      ENDIF

      NNSC=0
      DO NEQ=1,NEQMAX
         NNSN=NNS(NEQ)
         IF(NNSN.NE.0) THEN
            NNSC=NNSC+1
         ELSE
            GOTO 1000
         ENDIF
      ENDDO
 1000 CONTINUE
      NNSMAX=NNSC
      NEQRMAX=NEQMAX-NNSMAX

      NEQS=1
      NEQT=1
      DO NEQ=1,NEQMAX
         NNSN=NNS(NEQ)
         DO NEQ1=1,NEQMAX
            IF(NEQ1.GE.NEQS) THEN
               IF(NEQ1.LT.NNSN.OR.NNSN.EQ.0) THEN
                  NST(NEQ1)=NEQT
                  NEQT=NEQT+1
                  IF(NEQT.GT.NEQRMAX) GOTO 1200
               ELSEIF(NNSN.EQ.NEQ1) THEN
                  NST(NEQ1)=0
                  NEQS=NNSN+1
                  GOTO 1100
               ENDIF
            ENDIF
         ENDDO
 1100    CONTINUE
      ENDDO
 1200 CONTINUE

      DO NS=1,NSMAX
         IF(PZ(NS).NE.0.D0) THEN
            AMZ(NS)=PA(NS)*AMP/PZ(NS)**2
         ELSE
            AMZ(NS)=0.D0
         ENDIF
      ENDDO

      WRITE(6,600) 'NEQ','NSS','NSV','NNS','NST'
      DO NEQ=1,NEQMAX
         WRITE(6,610) NEQ,NSS(NEQ),NSV(NEQ),NNS(NEQ),NST(NEQ)
      ENDDO
 600  FORMAT(' ',5(' ',A))
 610  FORMAT(' ',5I4)

!     *** EQUATION SELECTOR ***

!     Format : NEA(species,equation) for all equations

      NEA(0:NSTM,0:3)=0
      DO NEQ=1,NEQMAX
         NEA(NSS(NEQ),NSV(NEQ))=NEQ
      ENDDO

!     CHECK WHETHER TURBULENT TRANSPORT MODEL HAS OFF-DIAGONAL PARTS

      IF(MDLKAI.EQ.61.OR.(MDLKAI.EQ.63.AND.MDLWLD.EQ.1)) THEN
         MDANOM=1
      ELSE
         MDANOM=0
      ENDIF

!     |-----------------------|
!     |MDDIAG |NCLASS |MDANOM |
!     |-------|-------|-------|
!     |   0   |   *   |   *   |
!     |   1   |   o   |   *   |
!     |   2   |   *   |   o   |
!     |   3   |   o   |   o   |
!     |-----------------------|

      IF(MDNCLS.EQ.0) THEN
         IF(MDANOM.EQ.1) THEN
            MDDIAG=2
         ELSE
            MDDIAG=0
         ENDIF
         RGFLS(1:NRMAX,1:5,1:NSMAX)=0.D0
         RQFLS(1:NRMAX,1:5,1:NSMAX)=0.D0
      ELSE
         IF(MDANOM.EQ.1) THEN
            MDDIAG=3
         ELSE
            MDDIAG=1
         ENDIF
      ENDIF

!     *** GRID POINT OF EDGE REGION ***

      NREDGE=NINT(0.93*NRMAX)

      RETURN
      END SUBROUTINE TR_EQS_SELECT

!     ***********************************************************

!           SORTER AS MAIN PART OF MODEL SELECTOR

!     ***********************************************************

      SUBROUTINE TR_TABLE(NS,NEQ,NSW,IND,INDH,INDHD)

      USE TRCOMM, ONLY : AME, AMP, MDLEQE, NEQMAX, NNS, NSMAX, NSS, NSV, PA
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NS, NSW
      INTEGER(4),INTENT(INOUT):: NEQ, IND, INDH, INDHD
      INTEGER(4) :: NEQI, NEQII, NNSN, NSVN
      REAL(8)    :: REM

      REM=AME/AMP
      IF(NS.LE.NSMAX) THEN
         IF(ABS(PA(NS)-REM).LE.1.D-10) THEN
!     electron
            NEQ=NEQ+1
            NSS(NEQ)=1
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1.AND.MDLEQE.EQ.0) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 100
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 100           CONTINUE
            ELSEIF(NSW.EQ.2.AND.IND.EQ.1) THEN
               IF(MDLEQE.EQ.0) THEN
                  IF(NSS(1).EQ.0) THEN
                     NNS(1)=2
                  ELSEIF(NSS(1).EQ.1) THEN
                     NNS(1)=1
                  ENDIF
               ELSEIF(MDLEQE.EQ.2) THEN
                  NNS(1)=3
               ENDIF
            ENDIF
         ELSEIF(PA(NS).EQ.1.D0.OR.ABS(PA(NS)-2.D0).LT.0.5D0) THEN
!     regard the particle whose mass is 1.0 as HYDROGEN
!     regard the particle whose mass is between 1.5 and 2.5 as DEUTERIUM
!     If bulk particle is hydrogen, INDH=1
            IF(PA(NS).EQ.1.D0) INDH=1
            NEQ=NEQ+1
!     If plasma is composed of hydrogen and deuterium, INDHD=1
            IF(INDH.EQ.1.AND.ABS(PA(NS)-2.D0).LT.0.5D0) THEN
               NSS(NEQ)=3
               INDHD=1
            ELSE
               NSS(NEQ)=2
            ENDIF
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 200
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 200           CONTINUE
            ENDIF
         ELSEIF(ABS(PA(NS)-3.D0).LT.0.5D0) THEN
!     regard the particle whose mass is between 2.5 and 3.5 as TRITIUM
            NEQ=NEQ+1
            IF(INDHD.EQ.1) THEN
               NSS(NEQ)=4
            ELSE
               NSS(NEQ)=3
            ENDIF
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 300
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 300           CONTINUE
            ENDIF
         ELSEIF(ABS(PA(NS)-4.D0).LT.0.5D0) THEN
!     regard the particle whose mass is between 3.5 and 4.5 as HELIUM
            NEQ=NEQ+1
            NSS(NEQ)=4
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 400
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 400           CONTINUE
            ENDIF
         ELSEIF(ABS(PA(NS)-12.D0).LT.3.D0.AND.NSMAX.EQ.3) THEN
!     regard the particle whose mass is between 9.0 and 15.0 as CARBON
            NEQ=NEQ+1
            NSS(NEQ)=3
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 500
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 500           CONTINUE
            ENDIF
         ELSEIF(PA(NS).EQ.0.D0) THEN
            IND=-1
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE TR_TABLE
