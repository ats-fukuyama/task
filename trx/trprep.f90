! trprep.f90

MODULE trprep

  PRIVATE
  PUBLIC tr_prep

CONTAINS

  SUBROUTINE tr_prep(ierr)

    USE trcomm
    USE trcom0, ONLY : NSTM
    USE trcom1, ONLY : NTAMAX,KDIRX
    USE trprof
    USE trbpsd
    USE trmetric
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

      CALL ALLOCATE_TRCOMM(IERR)
      IF(IERR.NE.0) RETURN

      CALL TR_EQS_SELECT(0)

      NRAMAX=INT(RHOA*NRMAX)
      DR = 1.D0/DBLE(NRMAX)

      NROMAX=NRMAX
      NRMAX=NRAMAX

      NT    = 0
      T     = 0.D0
      TPRE  = 0.D0
      TST   = 0.D0
      VSEC  = 0.D0
      NGR   = 0
      NGT   = 0
      NGST  = 0
      RIP   = RIPS

      CALL tr_prof

!     *** Initialize bpsd data ***

      CALL tr_bpsd_init
      CALL tr_bpsd_put(ierr)
      IF(ierr.NE.0) THEN
         write(6,'(A,I5)') 'XX tr_bpsd_put in tr_prof: ierr=',ierr
         STOP
      END IF

!     *** CALCULATE METRIC FACTOR ***

      CALL tr_set_metric(ierr)
      IF(ierr.NE.0) THEN
         write(6,'(A,I5)') 'XX tr_set_metric in tr_metric: ierr=',ierr
         STOP
      END IF
      
!     *** CALCULATE ANEAVE and ANC, ANFE ***

      CALL tr_prof_impurity

!     *** CALCULATE AJ, QP, BP, EZ ***

      IF(MODELG.EQ.2) CALL tr_prof_current

!     *** Initialize bpsd data ***

      CALL tr_bpsd_put(ierr)
      IF(ierr.NE.0) THEN
         write(6,'(A,I5)') 'XX tr_bpsd_put in tr_prof: ierr=',ierr
         STOP
      END IF

!     *** initilize graphic data ***

      GRG(1)=0.0
      GRM(1:NRMAX)  =SNGL(RM(1:NRMAX))
      GRG(2:NRMAX+1)=SNGL(RG(1:NRMAX))

      RETURN

    END SUBROUTINE tr_prep
  END MODULE trprep


!     ***********************************************************

!           MODEL SELECTOR

!     ***********************************************************

  SUBROUTINE TR_EQS_SELECT(INIT)

    USE TRCOMM, ONLY : &
         AMM, AMZ, MDDIAG, MDLEOI, MDLEQ0, MDLEQB, MDLEQE, MDLEQN, &
         MDLEQT, MDLEQU, MDLEQZ, model_chi_tb, &
         MDLWLD, MDNCLS, NEA, NEQM, NEQMAX, NEQMAXM, NNS, NREDGE, &
         NRMAX, NSCMAX, NSLMAX, NSM, NSMAX,      &
         NSNMAX, NSS, NST, NSTM, NSTMAX, NSV, NSZMAX, PA, PZ, RGFLS, RQFLS
    USE TRCOM1, ONLY : INS
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: INIT
    INTEGER:: IND, INDH, INDHD, MDANOM, MDSLCT, NEQ, NEQ1, NEQI
    INTEGER:: NEQRMAX, NEQS, NEQT, NNSC, NNSMAX, NNSN, NS, NSSN, NSVN
    INTEGER, SAVE :: NSSMAX


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
            AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
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

      IF(model_chi_tb.EQ.61.OR.(model_chi_tb.EQ.63.AND.MDLWLD.EQ.1)) THEN
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

      USE TRCOMM, ONLY : AME, AMM, MDLEQE, NEQMAX, NNS, NSMAX, NSS, NSV, PA, rkind
      IMPLICIT NONE
      INTEGER,INTENT(IN)   :: NS, NSW
      INTEGER,INTENT(INOUT):: NEQ, IND, INDH, INDHD
      INTEGER :: NEQI, NEQII, NNSN, NSVN
      REAL(rkind)    :: REM

      REM=AME/AMM
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
