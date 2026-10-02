! wfwaveK.f90

!     ********** WF WAVE SOLVER **********

MODULE wfwaveK

  PRIVATE
  PUBLIC wf_wave_kinetic
  PUBLIC wf_cmcalc_kinetic
  PUBLIC wf_wpre_kinetic

CONTAINS

  SUBROUTINE wf_wave_kinetic

    USE wfcomm
    IMPLICIT NONE
    INTEGER :: IERR
    REAL :: GTMAIN,GTSOLV,GCPUT0,GCPUT1,GCPUT2,GCPUT3

    GTMAIN=0.0
    GTSOLV=0.0
  
    CALL GUTIME(GCPUT0)
  
    IF (nrank.EQ.0) WRITE(6,*) '--- WFWPRE start ---'
    CALL wf_wpre_kinetic(IERR)
    IF(IERR.NE.0) goto 9000

    IF (nrank.EQ.0) WRITE(6,*) '--- CVCALC start ---'
    CALL wf_cvcalc
  
    CALL GUTIME(GCPUT1)
  
    IF (nrank.EQ.0) WRITE(6,*) '--- CVSOLV start ---'
    CALL wf_cvsolv
  
    CALL GUTIME(GCPUT2)

    IF (nrank.EQ.0) WRITE(6,*) '--- CALFLD start ---'
    CALL wf_calfld
    CALL wf_pwrabs
    CALL wf_pwrrad
    CALL wf_lpefld

    CALL GUTIME(GCPUT3)
    GTSOLV=GTSOLV+GCPUT2-GCPUT1
    GTMAIN=GTMAIN+GCPUT3-GCPUT2+GCPUT1-GCPUT0

!  IF (nrank.EQ.0) WRITE (6,'(A/5F12.3)') &
!       "GCPUT0,GCPUT1,GCPUT2,GCPUT3,GTSOLV=", &
!        GCPUT0,GCPUT1,GCPUT2,GCPUT3,GTSOLV
!  IF (nrank.EQ.0) CALL LPEFLD

    IF (nrank.EQ.0) WRITE(6,100) GTMAIN,GTSOLV
100 FORMAT(' ','****** CPU TIME : MAIN = ',F10.3,' SEC',5X,&
                               ': SOLV = ',F10.3,' SEC ******')
  
9000 CONTINUE

    RETURN
  END SUBROUTINE wf_wave_kinetic

!     ********** WF WAVE PREPARATION **********

  SUBROUTINE wf_wpre_kinetic(IERR)

    USE wfcomm
    USE wfparm
    USE plload,ONLY: pl_load
    USE femmeshprep
    USE feminterpolate
    USE libbes
    IMPLICIT NONE
    REAL(rkind):: RGAMMA
    INTEGER,INTENT(OUT) :: IERR

    IERR=0

    CALL fem_meshprep

    CALL fem_setup_zone

    CALL pl_load(ierr)
    IF(IERR.NE.0) return

    IF(MODELG.EQ.11) THEN
       RGAMMA=ABS(Hpitch1*RRCH)
       IF(RGAMMA.LT.1.D-5) THEN
          HA1=0.D0
       ELSE
          HA1=RGAMMA*BESKNX(1,2.D0*RGAMMA)+BESKNX(2,2.D0*RGAMMA)
       ENDIF
       RKAP=SQRT((1.D0+2.D0*HA1)/(1.D0-2.D0*HA1))
       WRITE(6,'(A,1P2E12.4)') 'HA1,RKAP=',HA1,RKAP
    ENDIF
  

    SELECT CASE(MODELG)
    CASE(0,12)
       CALL wf_bpsi(RA,0.D0,PSIA)
    CASE(1:10,13)
       CALL wf_bpsi(RR+RA,0.D0,PSIA)
    END SELECT

    CALL wf_lpelmt

    IF (nrank.EQ.0) WRITE(6,*) '----- wf_setbdy start ---'
    CALL wf_setbdy(IERR)
    IF(IERR.NE.0) return

    IF (nrank.EQ.0) WRITE(6,*) '----- wf_setlsd start ---'
    CALL wf_setlsd
  
    IF (nrank.EQ.0) WRITE(6,*) '----- wf_modant start ---'
    CALL wf_modant(IERR)
    IF(IERR.NE.0) return
  
    IF (nrank.EQ.0) WRITE(6,*) '----- wf_setewg start ---'
    CALL wf_setewg
    IF(IERR.NE.0) return
  
    IF (nrank.EQ.0) WRITE(6,*) '----- wf_defmlen start ---'
    CALL wf_defmlen
  
    CALL wffld_allocate

    IF (nrank.EQ.0) CALL wf_view

    RETURN
  END SUBROUTINE wf_wpre_kinetic

!     ****** DIELECTRIC TENSOR ******

  SUBROUTINE wf_dtensr(NE,DTENS)

    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NE
    INTEGER :: IN,I,J,ID
    COMPLEX(rkind),INTENT(OUT):: DTENS(NSM,3,3,3)
    INTEGER    :: NS,NN
    REAL(rkind)    :: R,Z,WW,WP(NSM),WC(NSM),BABS,AL(3),RN(NSM),RTPR(NSM)
    REAL(rkind)    :: RTPP(NSM),RZCL(NSM),FP,FR,FZ,DR,DZ,F
    COMPLEX(rkind) :: CWP,CWC,CDT0,CDX0,CDP0,CDT,CDP,CDX,CDAMP
    COMPLEX(rkind) :: CRR,CRP,CRZ,CPR,CPP,CPZ,CZR,CZP,CZZ

    ! ----- initialize -----  

    DO J=1,3
       DO I=1,3
          DO IN=1,3
             DO NS=1,NSMAX
                DTENS(NS,IN,I,J)=(0.d0,0.d0)
             END DO
          END DO
       END DO
    END DO

    WW=2.D0*PI*RF*1.D6

    DO NS=1,NSMAX
       WP(NS)=PZ(NS)*PZ(NS)*AEE*AEE*1.D20/(PA(NS)*AMP*EPS0*WW*WW)
       WC(NS)=PZ(NS)*AEE/(PA(NS)*AMP*WW)
    ENDDO

    ! ----- collisional cold plasma model -----

    DO IN=1,3
     
       NN=NDELM(IN,NE)
       R=RNODE(NN)
       Z=ZNODE(NN)
     
       CALL WFSMAG(R,Z,BABS,AL)
       FR=AL(1)
       FP=AL(2)
       FZ=AL(3)
     
       CALL WFSDEN(R,Z,RN,RTPR,RTPP,RZCL)

       DO NS=1,NSMAX
        
          CWP = WP(NS)*RN(NS)/(1.D0+CII*RZCL(NS))
          CWC = WC(NS)*BABS  /(1.D0+CII*RZCL(NS))
          CDT0= CWP/(1.D0-CWC**2)
          CDX0= CII*CWP*CWC/(1.D0-CWC**2)
          CDP0= CWP
        
          CDT=CDT0       
          CDP=CDP0-CDT0
          CDX=CDX0

          CRR= CDT   +CDP*FR*FR
          CRP= CDX*FZ+CDP*FR*FP
          CRZ=-CDX*FP+CDP*FR*FZ
          CPR=-CDX*FZ+CDP*FP*FR
          CPP= CDT   +CDP*FP*FP
          CPZ= CDX*FR+CDP*FP*FZ
          CZR= CDX*FP+CDP*FZ*FR
          CZP=-CDX*FR+CDP*FZ*FP
          CZZ= CDT   +CDP*FZ*FZ
        
          DTENS(NS,IN,1,1)=DTENS(NS,IN,1,1)-CRR
          DTENS(NS,IN,1,2)=DTENS(NS,IN,1,2)-CRP
          DTENS(NS,IN,1,3)=DTENS(NS,IN,1,3)-CRZ
          DTENS(NS,IN,2,1)=DTENS(NS,IN,2,1)-CPR
          DTENS(NS,IN,2,2)=DTENS(NS,IN,2,2)-CPP
          DTENS(NS,IN,2,3)=DTENS(NS,IN,2,3)-CPZ
          DTENS(NS,IN,3,1)=DTENS(NS,IN,3,1)-CZR
          DTENS(NS,IN,3,2)=DTENS(NS,IN,3,2)-CZP
          DTENS(NS,IN,3,3)=DTENS(NS,IN,3,3)-CZZ

       END DO


       IF(WDAMP.GT.0.D0) THEN
          CDAMP=CII*PZCL(NSMAX)
          F=FDAMP
          IF(R-BDRMIN.LT.WDAMP) THEN
             ID=1
             IF(MDAMP.EQ.1.AND. &
                  Z.GT.ZDAMP_MIN.AND.Z.LT.ZDAMP_MAX) ID=0
             IF(ID.EQ.1) THEN
                DR=R-BDRMIN
!               DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DR)/(DR-CDAMP)
                DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DR)/(DR-CDAMP)
                DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DR)/(DR-CDAMP)
             END IF
          END IF
          IF(BDRMAX-R.LT.WDAMP) THEN
             ID=1
             IF(MDAMP.EQ.2.AND. &
                    Z.GT.ZDAMP_MIN.AND.Z.LT.ZDAMP_MAX) ID=0
             IF(ID.EQ.1) THEN
                DR=BDRMAX-R
!               DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DR)/(DR-CDAMP)
                DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DR)/(DR-CDAMP)
                DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DR)/(DR-CDAMP)
             END IF
          END IF
          IF(Z-BDZMIN.LT.WDAMP) THEN
             ID=1
             IF(MDAMP.EQ.3.AND. &
                  R.GT.RDAMP_MIN.AND.R.LT.RDAMP_MAX) ID=0
             IF(ID.EQ.1) THEN
                DZ=Z-BDZMIN
                DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DZ)/(DZ-CDAMP)
                DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DZ)/(DZ-CDAMP)
!               DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DZ)/(DZ-CDAMP)
             END IF
          END IF
          IF(BDZMAX-Z.LT.WDAMP) THEN
             ID=1
             IF(MDAMP.EQ.4.AND. &
                  R.GT.RDAMP_MIN.AND.R.LT.RDAMP_MAX) ID=0
             IF(ID.EQ.1) THEN
                DZ=BDZMAX-Z
                DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DZ)/(DZ-CDAMP)
                DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DZ)/(DZ-CDAMP)
!               DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DZ)/(DZ-CDAMP)
             END IF
          END IF
       END IF
    END DO

    RETURN
  END SUBROUTINE wf_dtensr

!     ***** rotation tensor *****

  SUBROUTINE wf_mutensr(NE,MU)

    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NE
    INTEGER :: ISD,NSD,I,J,K
    REAL(rkind),INTENT(OUT):: MU(3,3,6)
    REAL(rkind) :: A(3),B(3),C(3),L(3)

    DO ISD=1,3
       NSD=ABS(NSDELM(ISD,NE))
       IF(MODELWF.EQ.0) THEN
          L(ISD)=LSID(NSD)
       ELSE
          IF(NSDELM(ISD,NE).GT.0.D0) THEN
             L(ISD)=LSID(NSD)
          ELSE
             L(ISD)=-LSID(NSD)
          END IF
       END IF
    END DO

    CALL WFABC(NE,A,B,C)

    DO K=1,6
       DO J=1,3
          DO I=1,3
             MU(I,J,K)=0.d0
          END DO
       END DO
    END DO

    MU(1,1,1)= L(1)*B(2)
    MU(1,3,1)= L(1)*C(2)
    MU(1,2,4)= 1.D0
    MU(1,1,3)=-L(3)*B(3)
    MU(1,3,3)=-L(3)*C(3)

    MU(2,1,1)=-L(1)*B(1)
    MU(2,3,1)=-L(1)*C(1)
    MU(2,2,5)= 1.D0
    MU(2,1,2)= L(2)*B(3)
    MU(2,3,2)= L(2)*C(3)

    MU(3,1,2)=-L(2)*B(2)
    MU(3,3,2)=-L(2)*C(2)
    MU(3,2,6)= 1.D0
    MU(3,1,3)= L(3)*B(1)
    MU(3,3,3)= L(3)*C(1)

    RETURN
  END SUBROUTINE wf_mutensr

!     ****** current COEFFICIENT VECTOR CALCULATION ******
!     LIF: Line Integral of interpolation Function

  SUBROUTINE wf_cvcalc

    USE wfcomm
    IMPLICIT NONE
    INTEGER        :: NE,NA,IJ,IV,I
    REAL(rkind)    :: RW,PHASE,MU(3,3,6),A(3),B(3),C(3)
    REAL(rkind)    :: R1,Z1,R2,Z2,LIF(3),R21,Z21
    COMPLEX(rkind) :: CJ(3),CVJ

    RW=2.D0*PI*RF*1.D6

    DO NE=1,NEMAX
       DO IV=1,6
          CVTOT(IV,NE)=(0.d0,0.d0)
       ENDDO
    ENDDO
  
    DO NA=1,NAMAX
       PHASE =APH(NA)*PI/180.D0
       CVJ=CII*RW*RMU0*AJ(NA)*EXP(CII*(PHASE))
       IF(JNUM(NA).EQ.1) THEN
          NE=JELMT(1,NA)
          CALL WFABC(NE,A,B,C)
          R1=RJ(1,NA)
          Z1=ZJ(1,NA)
          CJ(1)=0.D0
          CJ(2)=CVJ*RR
          CJ(3)=0.D0
          SELECT CASE(MODELG)
          CASE(0,12)
             DO I=1,3
                LIF(I)=A(I)*RR &
                      +B(I)*R1*RR &
                      +C(I)*Z1*RR
             END DO
          CASE(1:10,13)
             DO I=1,3
                LIF(I)=A(I)*R1 &
                      +B(I)*R1*R1 &
                      +C(I)*Z1*R1
             END DO
          END SELECT
          CALL wf_mutensr(NE,MU)
          DO I=1,3
             DO IV=1,6
                CVTOT(IV,NE)= CVTOT(IV,NE) &
                             +LIF(I)*( MU(I,1,IV)*CJ(1) &
                                      +MU(I,2,IV)*CJ(2) &
                                      +MU(I,3,IV)*CJ(3))
             END DO
          END DO
       ELSE
          DO IJ=2,JNUM(NA)
             NE=JELMT(IJ,NA)
             CALL WFABC(NE,A,B,C)
             R1=RJ(IJ-1,NA)
             Z1=ZJ(IJ-1,NA)
             R2=RJ(IJ,NA)
             Z2=ZJ(IJ,NA)
             R21=R2-R1
             Z21=Z2-Z1

             CJ(1)=CVJ*R21  
             CJ(2)=(0.d0,0.d0) 
             CJ(3)=CVJ*Z21  

             SELECT CASE(MODELG)
             CASE(0,12)
                DO I=1,3
                   LIF(I)=A(I)*RR &
                         +B(I)*(R1+R2)*RR/2.D0 &
                         +C(I)*(Z1+Z2)*RR/2.D0
                END DO
             CASE(1:10,13)
                DO I=1,3
                   LIF(I)=A(I)*(R1+R2)/2.d0 &
                         +B(I)*(R2**2+R1*R2+R1**2)/3.d0 &
                         +C(I)*(R2*Z1+Z2*R1+2.d0*R2*Z2+2.d0*R1*Z1)/6.d0
                END DO
             END SELECT
             CALL wf_mutensr(NE,MU)

!        WRITE(16,*) NE
!        DO I=1,3
!           DO J=1,3
!              WRITE(16,'(2I3,1P6E12.4)') I,J,MU(I,J,1:6)
!           END DO
!        END DO

             DO I=1,3
                DO IV=1,6
                   CVTOT(IV,NE)= CVTOT(IV,NE)&
                                +LIF(I)*( MU(I,1,IV)*CJ(1)&
                                         +MU(I,2,IV)*CJ(2)&
                                         +MU(I,3,IV)*CJ(3))
!                         TEMP=LIF(I)*( MU(I,1,IV)*CJ(1)&
!                                    +MU(I,2,IV)*CJ(2)&
!                                    +MU(I,3,IV)*CJ(3))
!                         IF(ABS(TEMP).NE.0.D0) THEN
!                           WRITE(16,'(2I6,1P3E12.4)') I,IV,LIF(I),TEMP
!                           WRITE(16,'(1P3E12.4)') MU(I,1,IV),CJ(1)
!                           WRITE(16,'(1P3E12.4)') MU(I,2,IV),CJ(2)
!                           WRITE(16,'(1P3E12.4)') MU(I,3,IV),CJ(3)
!                        END IF
                END DO
             END DO
          END DO
       END IF
    END DO

!  DO NE=1,NEMAX
!     DO IV=1,6
!        IF(nrank.EQ.0.and.CVTOT(IV,NE).NE.(0.d0,0.d0)) &
!                                   & WRITE(16,*) NE,IV,CVTOT(IV,NE)
!     END DO
!  END DO

    RETURN
  END SUBROUTINE wf_cvcalc

!     ****** local ELEMENT MATRIX  ******

  SUBROUTINE wf_cmcalc_kinetic(NE)

    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NE
    INTEGER:: I,J

    ! --- initialize CM ---

    DO J=1,6
       DO I=1,6
          CM(I,J)=(0.D0,0.D0)
       END DO
    END DO

    CALL wf_cmcalcv(NE)

    CALL wf_cmcalcs(NE)

    CALL wf_cmcalcp(NE)

    RETURN
  END SUBROUTINE wf_cmcalc_kinetic


!     ****** LOCAL ELEMENT MATRIX : volume integral ******

  SUBROUTINE wf_cmcalcv(NE)

    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(in) :: NE 
    INTEGER :: I,J,K,M,N,ISD,NSD
    COMPLEX(rkind) :: CM1(3,3)
    REAL(rkind) :: S,L(3)
    REAL(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
    REAL(rkind) :: R(3),Z(3)

    ! --- initialize ---

    S=SELM(NE)
    CALL WFABC(NE,A,B,C)
    CALL WFNODE(NE,R,Z)

    DO ISD=1,3
       NSD=ABS(NSDELM(ISD,NE))
       IF(MODELWF.EQ.0) THEN
          L(ISD)=LSID(NSD)
       ELSE
          IF(NSDELM(ISD,NE).GT.0.D0) THEN
             L(ISD)=LSID(NSD)
          ELSE
             L(ISD)=-LSID(NSD)
          END IF
       END IF
    END DO

    DO ISD=1,3
       M=ISD
       N=ISD+1
       IF(N.gt.3) N=N-3
       AW(ISD)=L(ISD)*(A(M)*B(N)-A(N)*B(M))
       BW(ISD)=L(ISD)*(B(M)*C(N)-B(N)*C(M))
       CW(ISD)=L(ISD)*(C(M)*A(N)-C(N)*A(M))
    END DO

  ! ----- rotErotF term -----

  ! --- E1F1 ---

    DO J=1,3
       DO I=1,3
          CM1(I,J)=(0.d0,0.d0)
       END DO
    END DO

    SELECT CASE(MODELG)
    CASE(0,12)
       DO K=1,3 
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
                        +(RKZ**2)*RR &
                         *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*Z(K) &
                           +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*R(K) &
                           +BW(I)*BW(J)*(R(K)**2+Z(K)**2)) &
                         *S*AIF1(K) &
                        +4.d0*BW(I)*BW(J)*RR*S*AIF1(K)
             END DO
          END DO
       END DO
    CASE(1:10,13)
       DO K=1,3 
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
                        +(real(NPH)**2)/R(K) &
                         *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*Z(K) &
                           +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*R(K) &
                           +BW(I)*BW(J)*(R(K)**2+Z(K)**2)) &
                         *S*AIF1(K) &
                        +4.d0*BW(I)*BW(J)*R(K)*S*AIF1(K)
             END DO
          END DO
       END DO
    END SELECT
    DO J=1,3
       DO I=1,3
          CM(I,J)=CM(I,J)+CM1(I,J)
       END DO
    END DO

    ! --- E1F2 --- 

    DO J=1,3
       DO I=1,3
          CM1(I,J)=(0.d0,0.d0)
       END DO
    END DO

    SELECT CASE(MODELG)
    CASE(0,12)
       DO K=1,3
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
                        +(CII*RKZ*RR) &
                         *(-B(I) &
                            *(AW(J)-BW(J)*Z(K)) &
                           +C(I) &
                            *(CW(J)-BW(J)*R(K))) &
                         *S*AIF1(K)
             END DO
          END DO
       END DO
    CASE(1:10,13)
       DO K=1,3
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
!
                        +(CII*real(NPH)) &
                         *(-B(I) &
                          *(AW(J)-BW(J)*Z(K)) &
                          +C(I) &
                          *(CW(J)-BW(J)*R(K))) &
                         *S*AIF1(K)&
!
                        +(CII*real(NPH))&
                         *(-(AW(J)-BW(J)*Z(K))/R(K))&
                         *S*AIF2(I,K)
             END DO
          END DO
       END DO
    END SELECT

    DO J=1,3
       DO I=1,3
          CM(I+3,J)=CM(I+3,J)+CM1(I,J)
       END DO
    END DO

    ! --- E2F1 ---

    DO J=1,3
       DO I=1,3
          CM1(I,J)=(0.d0,0.d0)
       END DO
    END DO

    SELECT CASE(MODELG)
    CASE(0,12)
       DO K=1,3
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
                        -(CII*RKZ*RR) &
                         *(-B(J) &
                            *(AW(I)-BW(I)*Z(K)) &
                           +C(J) &
                            *(CW(I)-BW(I)*R(K))) &
                         *S*AIF1(K)
             END DO
          END DO
       END DO
    CASE(1:10,13)
       DO K=1,3
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
!
                       -(CII*real(NPH)) &
                        *(-B(J) &
                           *(AW(I)-BW(I)*Z(K)) &
                          +C(J)&
                           *(CW(I)-BW(I)*R(K))) &
                        *S*AIF1(K) &
!
                       -(CII*real(NPH)) &
                        *(-(AW(I)-BW(I)*Z(K))/R(K)) &
                        *S*AIF2(J,K)
             END DO
          END DO
       END DO
    END SELECT
    DO J=1,3
       DO I=1,3
          CM(I,J+3)=CM(I,J+3)+CM1(I,J)
       END DO
    END DO

  ! --- E2F2 ---

    DO J=1,3
       DO I=1,3
          CM1(I,J)=(0.d0,0.d0)
       END DO
    END DO

    SELECT CASE(MODELG)
    CASE(0,12)
       DO K=1,3
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
                        +(B(I)*B(J)+C(I)*C(J))*RR*S*AIF1(K)
             END DO
          END DO
       END DO
    CASE(1:10,13)
       DO K=1,3
          DO J=1,3
             DO I=1,3
                CM1(I,J)=CM1(I,J) &
                        +(B(I)*B(J)+C(I)*C(J))*R(K)*S*AIF1(K) &
                        +B(J)*S*AIF2(I,K) &
                        +B(I)*S*AIF2(J,K) &
                        +1.D0/(R(K))*S*AIF3(I,J,K)
             END DO
          END DO
       END DO
    END SELECT
    
    DO J=1,3
       DO I=1,3
          CM(I+3,J+3)=CM(I+3,J+3)+CM1(I,J)
       END DO
    END DO

    RETURN
  END SUBROUTINE wf_cmcalcv

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

  SUBROUTINE wf_cmcalcs(NE)

    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NE 
    INTEGER :: I,J,K,I1,J1,K1,ID,ISD,NSD,IND(2),ND1,ND2
    REAL(rkind) :: S,L(3),RN,ZN
    REAL(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
    REAL(rkind) :: R(3),Z(3)
    COMPLEX(rkind):: CTEMP

    ! --- if no boundary side, return ---

    ID=0
    DO ISD=1,3
       NSD=ABS(NSDELM(ISD,NE))
       IF(KASID(NSD).EQ.1) ID=1
    END DO
    IF(ID.EQ.0) RETURN

    ! --- initialize ---

    S=SELM(NE)
    CALL WFABC(NE,A,B,C)
    CALL WFNODE(NE,R,Z)

    DO ISD=1,3
       NSD=ABS(NSDELM(ISD,NE))
       IF(KASID(NSD).EQ.1) THEN
          IF(MODELWF.EQ.0) THEN
             L(ISD)=LSID(NSD)
          ELSE
             IF(NSDELM(ISD,NE).GT.0.D0) THEN
                L(ISD)=LSID(NSD)
             ELSE
                L(ISD)=-LSID(NSD)
             END IF
          END IF

          ND1=ISD
          ND2=ISD+1
          IF(ND2.GT.3) ND2=ND2-3
          AW(ISD)=L(ISD)*(A(ND1)*B(ND2)-A(ND2)*B(ND1))
          BW(ISD)=L(ISD)*(B(ND1)*C(ND2)-B(ND2)*C(ND1))
          CW(ISD)=L(ISD)*(C(ND1)*A(ND2)-C(ND2)*A(ND1))
       END IF
    END DO

    DO ISD=1,3
       NSD=ABS(NSDELM(ISD,NE))
       IF(KASID(NSD).EQ.1) THEN

          IND(1)=ISD
          IND(2)=ISD+1
          IF(IND(2).GT.3) IND(2)=IND(2)-3
          RN= (Z(IND(2))-Z(IND(1)))/L(ISD)
          ZN=-(R(IND(2))-R(IND(1)))/L(ISD)

          ! ----- rotErotF term -----

          ! --- E1F1 ---

          I=ISD
          J=ISD
          SELECT CASE(MODELG)
          CASE(0,12)
             DO K1=1,2
                K=IND(K1)
                CTEMP=CM(I,J)
                CM(I,J)=CM(I,J) &
                       -2.D0*BW(J)*(ZN*(AW(I)-BW(I)*Z(K)) &
                                   +RN*(CW(I)-BW(I)*R(K)))*RR*L(ISD)*AIE1(K)
!                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
!                         'CM:',NE,I,J,K,CM(I,J),CM(I,J)-CTEMP
!                    WRITE(21,'(5X,1P5E12.4)') &
!                         BW(J),ZN,AW(I),BW(I),Z(K)
!                    WRITE(21,'(5X,1P5E12.4)') &
!                         RN,CW(I),R(K),RR,L(ISD),AIE1(K)
             END DO
          CASE(1:10,13)
             DO K1=1,2
                K=IND(K1)
                CM(I,J)=CM(I,J) &
                       -2.D0*BW(J)*(ZN*(AW(I)-BW(I)*Z(K)) &
                                   +RN*(CW(I)-BW(I)*R(K)))*R(K)*L(ISD)*AIE1(K)
             END DO
          END SELECT

          ! --- E2F1 --- 

          I=ISD
          DO J1=1,2
             J=IND(J1)
             SELECT CASE(MODELG)
             CASE(0,12)
                DO K1=1,2
                   K=IND(K1)
                   CTEMP=CM(I,J+3)
                   CM(I,J+3)=CM(I,J+3) &
                            +CII*RKZ*RR &
                             *(-ZN*(-CW(J)+BW(J)*R(K)) &
                               -RN*( AW(J)-BW(J)*Z(K)))*L(ISD)*AIE2(J,K)
                END DO
             CASE(1:10,13)
                DO K1=1,2
                   K=IND(K1)
                   CM(I,J+3)=CM(I,J+3) &
                            +CII*NPH &
                             *(-ZN*(-CW(J)+BW(J)*R(K)) &
                               -RN*( AW(J)-BW(J)*Z(K)))*L(ISD)*AIE2(J,K)
                END DO
             END SELECT
          END DO

          ! --- E2F2 ---

          DO J1=1,2
             J=IND(J1)
             DO I1=1,2
                I=IND(I1)
                SELECT CASE(MODELG)
                CASE(0,12)
                   DO K1=1,2
                      K=IND(K1)
                      CTEMP=CM(I+3,J+3)
                      CM(I+3,J+3)=CM(I+3,J+3) &
                                 +(RN*B(J)+ZN*C(J))*RR &
                                  *L(ISD)*AIE3(I,J,K)
                   END DO
                CASE(1:10,13)
                   DO K1=1,2
                      K=IND(K1)
                      CM(I+3,J+3)=CM(I+3,J+3) &
                                 +((RN*B(J)+ZN*C(J))*R(K) &
                                   +RN*(A(J)+B(J)*R(K)+C(J)*Z(K))) &
                                  *L(ISD)*AIE3(I,J,K)
                   END DO
                END SELECT
             END DO
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE wf_cmcalcs

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

  SUBROUTINE wf_cmcalcp(NE)

    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(in) :: NE 
    INTEGER :: I,J,K,NS,IN,JJ,II
    REAL(rkind) :: RW,WC,WC2
    REAL(rkind) :: S
    REAL(rkind) :: A(3),B(3),C(3)
    REAL(rkind) :: R(3),Z(3),MU(3,3,6)
    COMPLEX(rkind):: CM2(6,6)
    COMPLEX(rkind) :: DTENS(NSM,3,3,3)
    COMPLEX(rkind) :: DTENST(3,3,3)

    ! --- initialize ---

    RW=2.D0*PI*RF*1.D6
    WC=RW/VC
    WC2=WC**2

    S=SELM(NE)  
    CALL WFABC(NE,A,B,C)
    CALL WFNODE(NE,R,Z)

    ! ----- dielectric tensor term -----

    CALL wf_dtensr(NE,DTENS)
    CALL wf_mutensr(NE,MU)

    DO J=1,6
       DO I=1,6
          CM2(I,J)=(0.D0,0.D0)
       END DO
    END DO

    DO J=1,3
       DO I=1,3
          DO IN=1,3
             IF(I.EQ.J) THEN
                DTENST(IN,I,J)=(1.d0,0.d0)
             ELSE
                DTENST(IN,I,J)=(0.d0,0.d0)
             END IF
          END DO
       END DO
    END DO

    ! --- assemble dielectric tensor ---

    DO J=1,3
       DO I=1,3
          DO IN=1,3
             DO NS=1,NSMAX
                DTENST(IN,I,J)=DTENST(IN,I,J)+DTENS(NS,IN,I,J)
             END DO
          END DO
       END DO
    END DO

    SELECT CASE(MODELG)
    CASE(0,12)
       DO JJ=1,6
          DO II=1,6
             DO K=1,3
                DO J=1,3
                   DO I=1,3
                      CM2(II,JJ)= CM2(II,JJ)&
                                +((MU(I,1,II)*DTENST(J,1,1)&
                                  +MU(I,2,II)*DTENST(J,2,1)&
                                  +MU(I,3,II)*DTENST(J,3,1))&
                                  *MU(K,1,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,2)&
                                  +MU(I,2,II)*DTENST(J,2,2)&
                                  +MU(I,3,II)*DTENST(J,3,2))&
                                  *MU(K,2,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,3)&
                                  +MU(I,2,II)*DTENST(J,2,3)&
                                  +MU(I,3,II)*DTENST(J,3,3))&
                                  *MU(K,3,JJ))&
                                 *RR*S*AIF3(I,J,K)
                   END DO
                END DO
             END DO
          END DO
       END DO
    CASE(1:10,13)
       DO JJ=1,6
          DO II=1,6
             DO K=1,3
                DO J=1,3
                   DO I=1,3
                      CM2(II,JJ)= CM2(II,JJ)&
                                +((MU(I,1,II)*DTENST(J,1,1)&
                                  +MU(I,2,II)*DTENST(J,2,1)&
                                  +MU(I,3,II)*DTENST(J,3,1))&
                                  *MU(K,1,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,2)&
                                  +MU(I,2,II)*DTENST(J,2,2)&
                                  +MU(I,3,II)*DTENST(J,3,2))&
                                  *MU(K,2,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,3)&
                                  +MU(I,2,II)*DTENST(J,2,3)&
                                  +MU(I,3,II)*DTENST(J,3,3))&
                                  *MU(K,3,JJ))&
                                 *R(J)*S*AIF3(I,J,K)
                   END DO
                END DO
             END DO
          END DO
       END DO
    END SELECT

    DO J=1,6
       DO I=1,6
          CM(I,J)=CM(I,J)-WC2*CM2(I,J)
       END DO
    END DO

!  WRITE(*,*) 'NE=',NE
!  DO I=1,6
!     WRITE(6,'(6(A1,E9.3,A1,E9.3,A1))')&
!          "(",real(CM(I,1)),',',aimag(CM(I,1)),")",&
!          "(",real(CM(I,2)),',',aimag(CM(I,2)),")",&
!          "(",real(CM(I,3)),',',aimag(CM(I,3)),")",&
!          "(",real(CM(I,4)),',',aimag(CM(I,4)),")",&
!          "(",real(CM(I,5)),',',aimag(CM(I,5)),")",&
!          "(",real(CM(I,6)),',',aimag(CM(I,6)),")"
!  END DO

    RETURN
  END SUBROUTINE wf_cmcalcp

!     ****** SOLV MATRIX EQUATION *****

SUBROUTINE wf_cvsolv

  use wfcomm
  use libmpi
  use libmtx
  use libqsort
  implicit none
  integer :: ISD,NSD,nnd,nnd1,nnd2
  integer :: NE,NN,nv,nvmax
  integer :: I,J,KK,LL
  integer :: JNSD,JNN,INSD,INN
  integer :: IN,INV,JNV,KB
  integer :: itype
  integer :: its
  integer :: JMIN,JMAX,MILEN,MJLEN
  integer :: NNZ,NNZMAX,NNZME      !Number of Non-Zero Matrix Element
  integer,dimension(:),ALLOCATABLE :: NEFLAG
  integer :: ORIENTJ,ORIENTI
  real :: cputime1,cputime2
  real(rkind) :: x,y,val
  complex(rkind):: CEB
  complex(rkind),dimension(:),ALLOCATABLE :: CRVP,CEQP
  integer(long),dimension(:),ALLOCATABLE :: NSEQ
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: VAL_SORT
  INTEGER(long),DIMENSION(:),ALLOCATABLE:: NV_SORT
  INTEGER,DIMENSION(:),ALLOCATABLE:: ntyp_nv,nnsd_nv
  INTEGER(long):: IX,IY,nv_old

  ! ----- initialize ------
  
  do NV=1,MLEN
     CSV(NV) =(0.d0,0.d0)
  end do

  ! --- decide istart,iend ---

  IF(nrank.EQ.0) write(6,*) 'MLEN=',MLEN
  call mtxc_setup(MLEN,istart,iend,0)
  call mtxc_cleanup

  ! ----- set NV_NSD & NV_NN -----

  ALLOCATE(ntyp_nv(nnmax+nsdmax),nnsd_nv(nnmax+nsdmax))

  NV=0
  do NSD=1,NSDMAX
     if(KASID(NSD).eq.1) then
        NVNSD(NSD)=0
     else
        NV=NV+1
        NVNSD(NSD)=NV
        NTYP_NV(NV)=2
        NNSD_NV(NV)=NSD
     end if
  end do
  do NN=1,NNMAX
     if(KANOD(NN).eq.1) then
        NVNN(NN)=0
     else
        NV=NV+1
        NVNN(NN)=NV
        NTYP_NV(NV)=1
        NNSD_NV(NV)=NN
     end if
  end do
  NVMAX=NV

  ALLOCATE(NV_SORT(NVMAX),VAL_SORT(NVMAX))

  DO NV=1,NVMAX
     IF(NTYP_NV(NV).EQ.1) THEN
        NN=NNSD_NV(NV)
        X=RNODE(NN)
        Y=ZNODE(NN)
        VAL=sort_weight_x*X+sort_weight_y*Y
     ELSE
        NSD=NNSD_NV(NV)
        NND1=NDSID(1,NSD)
        NND2=NDSID(2,NSD)
        X=0.5D0*(RNODE(NND1)+RNODE(NND2))
        Y=0.5D0*(ZNODE(NND1)+ZNODE(NND2))
     END IF
     VAL=sort_weight_x*X+sort_weight_y*Y
     nv_sort(NV)=NV
     val_sort(NV)=VAL
  END DO

  CALL qsort_dl(val_sort,nv_sort)

  DO nv=1,nvmax
     nv_old=nv_sort(nv)
     IF(ntyp_nv(nv_old).EQ.1) THEN
        nnd=nnsd_nv(nv_old)
        nvnn(nnd)=nv
!        X=RNODE(nnd)
!        Y=ZNODE(nnd)
!        VAL=sort_weight_x+sort_weight_y*Y
!        write(21,'(I10,A,I10,1P3E12.4)') nv,' node ',nnd,X,Y,VAL
     ELSE
        nsd=nnsd_nv(nv_old)
        nvnsd(nsd)=nv
!        NND1=NDSID(1,NSD)
!        NND2=NDSID(2,NSD)
!        X=0.5D0*(RNODE(NND1)+RNODE(NND2))
!        Y=0.5D0*(ZNODE(NND1)+ZNODE(NND2))
!        VAL=sort_weight_x*X+sort_weight_y*Y
!        write(21,'(I10,A,I10,1P3E12.4)') nv,' side ',nsd,X,Y,VAL
     END IF
  END DO
     
  DEALLOCATE(NV_SORT,VAL_SORT,ntyp_nv,nnsd_nv)

!  do NSD=1,NSDMAX
!     write(*,*) NSD,KASID(NSD),NV_NSD(NSD)
!  end do
!  do NN=1,NNMAX
!     write(*,*) NN,KANOD(NN),NV_NN(NN)
!  end do

  ! ----- set NEFLAG ------

  allocate(NEFLAG(NEMAX))
  do NE=1,NEMAX
     NEFLAG(NE)=0
  end do

  do NE=1,NEMAX
     do ISD=1,3
        NSD=ABS(NSDELM(ISD,NE))
        NV=NVNSD(NSD)
        if(NV.ge.istart.and.&
           NV.le.iend ) then
           NEFLAG(NE)=1
        end if
     end do

     if(NEFLAG(NE).eq.0) then
        do IN=1,3
           NN=NDELM(IN,NE)
           NV=NVNN(NN)
           if(NV.ge.istart.and.&
              NV.le.iend) then
              NEFLAG(NE)=1
           end if
        end do
     end if

  end do

  ! ----- set MJLEN -----

  JMIN=MLEN
  JMAX=0
  do NE=1,NEMAX

     if(NEFLAG(NE).eq.0) goto 8100

     LL=0
     DO J=1,6
        if(J.ge.1.and.J.le.3) then
           JNSD=ABS(NSDELM(J,NE))
           JNV =NVNSD(JNSD)
        else
           JNN=NDELM(J-3,NE)
           JNV=NVNN(JNN)
        end if
        if (JNV.eq.0) goto 8110
        LL=JNV

        KK=0
        DO I=1,6
           if(I.ge.1.and.I.le.3) then
              INSD=ABS(NSDELM(I,NE))
              INV=NVNSD(INSD)
           else
              INN=NDELM(I-3,NE)
              INV=NVNN(INN)
           end if
           if(INV.eq.0) goto 8120
           KK=INV

           if((KK.ge.istart).and.&
              (KK.le.iend  )) then
              if(LL.lt.JMIN) JMIN=LL
              if(LL.gt.JMAX) JMAX=LL
           end if

8120       continue
        ENDDO
8110    continue
     ENDDO
     
8100 continue
  end do

! ------ Count non-zero component -------

  NNZ=0

  DO NE=1,NEMAX
     IF(NEFLAG(NE).NE.0) THEN

        LL=0
        DO J=1,6
           IF(J.ge.1.and.J.le.3) then
              JNSD=NSDELM(J,NE)
              if(JNSD.lt.0) then
                 JNSD=-JNSD
              end if
              JNV =NVNSD(JNSD)
           else
              JNN=NDELM(J-3,NE)
              JNV=NVNN(JNN)
           END IF
           LL=JNV

           if ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

              KK=0
              DO I=1,6
                 if(I.ge.1.and.I.le.3) then
                    INSD=NSDELM(I,NE)
                    if(INSD.lt.0) then
                       INSD=-INSD
                    end if
                    INV =NVNSD(INSD)
                 else
                    INN=NDELM(I-3,NE)
                    INV=NVNN(INN)
                 end if
                 KK=INV

                 if((KK.ge.istart).and.&
                      (KK.le.iend  )) then
                    NNZ=NNZ+1
                 end if
              END DO
           END if
        ENDDO
     END IF
  END DO

  NNZMAX=NNZ
  MILEN=iend-istart+1
  MJLEN=JMAX-JMIN+1

  ! ----- set CEQP,CRVP -----

  allocate(CEQP(NNZMAX),NSEQ(NNZMAX),CRVP(MILEN))
  DO NNZ=1,NNZMAX
     CEQP(NNZ)=(0.d0,0.d0)
     NSEQ(NNZ)=0  ! (i-istart)*MJLEN+j-jmin
  END DO
  DO I=1,MILEN
     CRVP(I)=(0.d0,0.d0)
  END DO

! ------ set grobal matrix -------

  NNZ=0

  do NE=1,NEMAX
     if(NEFLAG(NE).eq.0) goto 8000
     CALL wf_cmcalc_kinetic(NE)
     
!    === ASSEMBLY ===
!    If KK (or LL) is out of assigned range,  
!      CVSOLV do not save the matrix element. 

!    --- inside of the boundary ---
     LL=0
     DO J=1,6
        ORIENTJ=1
        if(J.ge.1.and.J.le.3) then
           JNSD=NSDELM(J,NE)
           if(JNSD.lt.0) then
              JNSD=-JNSD
              ORIENTJ=-1
           end if
           JNV =NVNSD(JNSD)
        else
           JNN=NDELM(J-3,NE)
           JNV=NVNN(JNN)
        end if
        LL=JNV

        if ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

           KK=0
           DO I=1,6
              ORIENTI=1
              if(I.ge.1.and.I.le.3) then
                 INSD=NSDELM(I,NE)
                 if(INSD.lt.0) then
                    INSD=-INSD
                    ORIENTI=-1
                 end if
                 INV =NVNSD(INSD)
              else
                 INN=NDELM(I-3,NE)
                 INV=NVNN(INN)
              end if
              KK=INV
              if((KK.ge.istart).and.&
                   (KK.le.iend  )) then
                 if(abs(CM(I,J)).ne.0.d0) THEN 
                    NNZ=NNZ+1
                    IX=KK-istart
                    IY=MJLEN
                    NSEQ(NNZ)=IX*IY+LL-jmin
                    CEQP(NNZ)=ORIENTJ*ORIENTI*CM(I,J)
                 END if
              end if
           END DO
        END if
     ENDDO

!    --- Contribution from the boundary electric field ---

     LL=0
     DO J=1,6
        ORIENTJ=1
        if(J.ge.1.and.J.le.3) then
           JNSD=NSDELM(J,NE)
           if(JNSD.lt.0) then
              JNSD=-JNSD
              ORIENTJ=-1
           end if
           KB =KBSID(JNSD)
           IF(KB.NE.0) CEB=CEBSD(KB)
        else
           JNN=NDELM(J-3,NE)
           KB=KBNOD(JNN)
           IF(KB.NE.0) CEB=CEBND(KB)
        end if

        IF(KB.NE.0) THEN
           IF(ABS(CEB).GT.0.D0) THEN
              KK=0
              DO I=1,6
                 ORIENTI=1
                 if(I.ge.1.and.I.le.3) then
                    INSD=NSDELM(I,NE)
                    if(INSD.lt.0) then
                       INSD=-INSD
                       ORIENTI=-1
                    end if
                    INV =NVNSD(INSD)
                 else
                    INN=NDELM(I-3,NE)
                    INV=NVNN(INN)
                 end if
                 KK=INV
                 if((KK.ge.istart).and.&
                      (KK.le.iend  )) then
                    CRVP(KK-istart+1)=CRVP(KK-istart+1) &
                         -ORIENTI*ORIENTJ*CM(I,J)*CEB
                 end if
              END DO
           END IF
        END IF
     ENDDO

!    --- Contribution from the antenna current ---

     KK=0
     DO I=1,6
        ORIENTI=1
        if(I.ge.1.and.I.le.3) then
           INSD=NSDELM(I,NE)
           if(INSD.lt.0) then
              INSD=-INSD
              ORIENTI=-1
           end if
           INV =NVNSD(INSD)
        else
           INN=NDELM(I-3,NE)
           INV=NVNN(INN)
        end if
        KK=INV
        if((KK.ge.istart).and.&
           (KK.le.iend  )) then
           CRVP(KK-istart+1)=CRVP(KK-istart+1)+ORIENTI*CVTOT(I,NE)
        end if
     ENDDO

8000 continue
  end do

  IF(nrank.EQ.0) write(6,*) 'wf_cvsolv: sort started (kinetic)'
  CALL qsort_lc(NSEQ,CEQP)
  IF(nrank.EQ.0) write(6,*) 'wf_cvsolv: reduction started (kinetic)'
  NNZME=1
  DO NNZ=2,NNZMAX
     IF(NSEQ(NNZ).EQ.NSEQ(NNZME)) THEN
        CEQP(NNZME)=CEQP(NNZME)+CEQP(NNZ)
     ELSE
        NNZME=NNZME+1
        NSEQ(NNZME)=NSEQ(NNZ)
        CEQP(NNZME)=CEQP(NNZ)
     END IF
  END DO

  if(nrank.eq.0) write(6,'(A77)') &
  '      nrank     istart       iend      MILEN      MJLEN     NNZMAX      NNZME'
  call mtx_barrier
  write(6,'(7I11)') nrank,istart,iend,iend-istart+1, &
                    JMAX-JMIN+1,NNZMAX,NNZME

  ! ----- initialize for parallel computing -----

  itype = 0
  call mtxc_setup(MLEN,istart,iend,nzmax=NNZME)

  do NNZ=1,NNZME
     IF(ABS(CEQP(NNZ)).GT.0.D0) THEN
        i=NSEQ(NNZ)/MJLEN
        j=NSEQ(NNZ)-i*MJLEN
!        if(i+istart.lt.0) then
        if(i.LT.0.OR.i.GT.iend-istart) then
           WRITE(6,'(A/6I12)') 'NNZ,NSEQ(NNZ),i,istart,i+istart,iend=', &
                               NNZ,NSEQ(NNZ),i,istart,i+istart,iend
           STOP
        END if
        call mtxc_set_matrix(i+istart,j+jmin,CEQP(NNZ))
     END IF
  end do

  do i=istart,iend
     call mtxc_set_source(i,CRVP(i-istart+1))
  end do

  call GUTIME(cputime1)

  call mtxc_solve(itype,tolerance,its)
  !zmumps always return "its = 0"
  if(nrank.eq.0) write(6,*) 'Iteration Number=',its

  call GUTIME(cputime2)

  call mtxc_gather_vector(CSV)

  deallocate(CEQP,NSEQ,CRVP)
  deallocate(NEFLAG)
  call mtxc_cleanup
  RETURN
END SUBROUTINE wf_cvsolv

!     ******* ELECTRIC AND MAGNETIC FIELD CALCULATION *******

  SUBROUTINE wf_calfld

    USE wfcomm
    IMPLICIT NONE
    INTEGER :: NN,NSD,NV,NE,NSIDE,ND1,ND2,NBSD,I
    REAL(rkind):: RW,DX,DY,PF1,PF2,PFLUX,S
    INTEGER:: NSDA(3)
    REAL(rkind):: LSD(3),A(3),B(3),C(3)
    COMPLEX(rkind):: CESUM,CEN(3)

    RW=2.D0*PI*RF*1.D6

    DO NSD=1,NSDMAX ! initialise side E field
       CESD(NSD)=(0.d0,0.d0)
    ENDDO
    DO NN=1,NNMAX   ! initialize node E filed
       CEND(NN) =(0.d0,0.d0)
    END DO

    DO NSD=1,NSDMAX
       NV=NVNSD(NSD)
       IF (NV.EQ.0) THEN
          IF(KBSID(NSD).NE.0) THEN
             CESD(NSD)=CEBSD(KBSID(NSD))
          ELSE
             CESD(NSD)=(0.d0,0.d0)
          END IF
       ELSE
          CESD(NSD)=CSV(NV)
       END IF
       !     if(nrank.EQ.0) WRITE(6,*) NSD,CESD(NSD),KASID(NSD)
    END DO
    
    DO NN=1,NNMAX
       NV=NVNN(NN)
       IF (NV.EQ.0) THEN
          IF(KBNOD(NN).NE.0) THEN
             CEND(NN)=CEBND(KBNOD(NN))
          ELSE
             CEND(NN)=(0.d0,0.d0)
          END IF
       ELSE
          CEND(NN)=CSV(NV)
       END IF
    END DO

    ! --- calculate B-field in a element ---

    DO NE=1,NEMAX
       CESUM=(0.D0,0.D0)
       DO NSIDE=1,3
          NSD=ABS(NSDELM(NSIDE,NE))
          IF(NSDELM(NSIDE,NE).GT.0) THEN
             CESUM=CESUM+CESD(NSD)*LSID(NSD)
          ELSE
             CESUM=CESUM-CESD(NSD)*LSID(NSD)
          END IF
       END DO
       CBELM(NE)=CESUM/(-CI*RW*SELM(NE))
    END DO

  ! --- calculate B-field parallel to a side ---
  
    DO NE=1,NEMAX
       CALL WFABC(NE,A,B,C)
       S=SELM(NE)
       DO I=1,3
          NSDA(I)=ABS(NSDELM(I,NE))
          LSD(I)=LSID(NSDA(I))
          CEN(I)=CEND(NDELM(I,NE))
       END DO
       DO I=1,3
          CBSD(NSDA(I))=-2.D0*S/(CI*RW*LSD(I)) &
               *(CEN(1)*(B(1)*B(I)+C(1)*C(I)) &
                +CEN(2)*(B(2)*B(I)+C(2)*C(I)) &
                +CEN(3)*(B(3)*B(I)+C(3)*C(I)))
          IF(NSDELM(I,NE).NE.0) CBSD(NSDA(I))=-CBSD(NSDA(I))
       END DO
    END DO

    ! --- calculate power flux across a side (CESD*CBELM) ---
  
    PFLUXSD(1:NSDMAX)=0.D0
    PFLUXSDX(1:NSDMAX)=0.D0
    PFLUXSDY(1:NSDMAX)=0.D0
    DO NE=1,NEMAX
       DO NSIDE=1,3
          NSD=ABS(NSDELM(NSIDE,NE))
          ND1=NDSID(1,NSD)
          ND2=NDSID(2,NSD)
          DX=RNODE(ND2)-RNODE(ND1)
          DY=ZNODE(ND2)-ZNODE(ND1)
          DX=DX/SQRT(DX**2+DY**2)
          DY=DY/SQRT(DX**2+DY**2)
          IF(NSDELM(NSIDE,NE).GT.0) THEN ! positive direction
             IF(KASID(NSD).EQ.1) THEN ! boundary side
                PFLUX=REAL(DCONJG(CESD(NSD))*CBELM(NE))*LSID(NSD)/RMU0
             ELSE
                PFLUX=0.5D0*REAL(DCONJG(CESD(NSD))*CBELM(NE))*LSID(NSD)/RMU0
             END IF
          ELSE                           ! negative direction
             IF(KASID(NSD).EQ.1) THEN ! boundary side
                PFLUX=-REAL(DCONJG(CESD(NSD))*CBELM(NE))*LSID(NSD)/RMU0
             ELSE
                PFLUX=-0.5D0*REAL(DCONJG(CESD(NSD))*CBELM(NE))*LSID(NSD)/RMU0
             END IF
          END IF
          PFLUXSDX(NSD)=PFLUXSDX(NSD)-PFLUX*DY
          PFLUXSDY(NSD)=PFLUXSDX(NSD)+PFLUX*DX
          PFLUXSD(NSD)=PFLUXSD(NSD)+PFLUX
          IF(KASID(NSD).NE.0) PFLUXBDY(KBSID(NSD))=PFLUXBDY(KBSID(NSD))+PFLUX
       END DO
    END DO

    ! --- calculate power flux across a side (CEND*CBSD) ---

    PFLUXND(1:NSDMAX)=0.D0
    PFLUXNDX(1:NNMAX)=0.D0
    PFLUXNDY(1:NNMAX)=0.D0
    DO NE=1,NEMAX
       DO NSIDE=1,3
          NSD=ABS(NSDELM(NSIDE,NE))
          ND1=NDSID(1,NSD)
          ND2=NDSID(2,NSD)
          DX=RNODE(ND2)-RNODE(ND1)
          DY=ZNODE(ND2)-ZNODE(ND1)
          DX=DX/SQRT(DX**2+DY**2)
          DY=DY/SQRT(DX**2+DY**2)
          IF(NSDELM(NSIDE,NE).GT.0) THEN ! positive direction
             PF1=0.5D0*REAL(DCONJG(CEND(ND1))*CBSD(NSD))*LSID(NSD)/RMU0
             PF2=0.5D0*REAL(DCONJG(CEND(ND2))*CBSD(NSD))*LSID(NSD)/RMU0
          ELSE
             PF1=-0.5D0*REAL(DCONJG(CEND(ND1))*CBSD(NSD))*LSID(NSD)/RMU0
             PF2=-0.5D0*REAL(DCONJG(CEND(ND2))*CBSD(NSD))*LSID(NSD)/RMU0
          END IF
          PFLUXNDX(ND1)=PFLUXNDX(ND1)-PF1*DY
          PFLUXNDY(ND1)=PFLUXNDY(ND1)+PF1*DX
          PF2=0.5D0*REAL(DCONJG(CEND(ND2))*CBSD(NSD))*LSID(NSD)/RMU0
          PFLUXNDX(ND2)=PFLUXNDX(ND2)-PF2*DY
          PFLUXNDY(ND2)=PFLUXNDY(ND2)+PF2*DX
          PFLUXND(NSD)=PFLUXND(NSD)+PF1+PF2
       END DO
    END DO

    PFLUXX(1:NNMAX)=PFLUXNDX(1:NNMAX)
    PFLUXY(1:NNMAX)=PFLUXNDY(1:NNMAX)
    DO NSD=1,NSDMAX
       ND1=NDSID(1,NSD)
       ND2=NDSID(2,NSD)
       PFLUXX(ND1)=PFLUXX(ND1)+PFLUXSDX(NSD)
       PFLUXY(ND1)=PFLUXY(ND1)+PFLUXSDY(NSD)
    END DO

    WRITE(6,'(A)') 'PFLUXBDY:'

    PFLUXBDY(1:NBSID)=0.D0
    PFLUX=0.D0
    DO NSD=1,NSDMAX
       IF(KASID(NSD).NE.0) THEN
          NBSD=KBSID(NSD)
          PFLUXBDY(NBSD)=PFLUXBDY(NBSD) &
               +PFLUXSD(NSD)+PFLUXND(NSD)
          PFLUX=PFLUX+PFLUXBDY(NBSD)
          IF(PFLUXBDY(NBSD).NE.0.D0) THEN
             WRITE(6,'(A,5ES12.4)') 'PFLUX XY:', &
                  0.5D0*(RNODE(NDSID(1,NSD))+RNODE(NDSID(2,NSD))), &
                  0.5D0*(ZNODE(NDSID(1,NSD))+ZNODE(NDSID(2,NSD))), &
                  PFLUXBDY(NBSD),PFLUXSD(NSD),PFLUXND(NSD)
          END IF
       END IF
    END DO
     
    WRITE(6,'(A,1ES12.4)') 'PFLUX_IN=',PFLUX

    RETURN
  END SUBROUTINE wf_calfld

!     ******* POWER ABSORPTION *******

  SUBROUTINE wf_pwrabs

    USE wfcomm
    USE libmpi
    IMPLICIT NONE

    INTEGER    :: NE,IN,NN,NSD,NS
    INTEGER    :: I,J,K,II,JJ
    REAL(rkind)    :: RW,S,MU(3,3,6),R(3),Z(3)
    COMPLEX(rkind) :: DTENS(NSM,3,3,3),CTENS(NSM,3,3,3)
    COMPLEX(rkind) :: CIWE,CINT(NSM,6,6),CE(6)
    INTEGER,ALLOCATABLE:: nelm_len_nrank(:),nelm_pos_nrank(:)
    REAL(rkind),ALLOCATABLE:: rdata(:),rdata_tot(:)
    INTEGER:: ipos,n,nblk,ndata,nelm1,nelm2,nsize_high,nsize_low

    ALLOCATE(nelm_len_nrank(0:nsize-1),nelm_pos_nrank(0:nsize-1))
    nblk=nemax/nsize
    nsize_high=nemax-nblk*nsize
    nsize_low=nsize-nsize_high
    ipos=0
    DO n=0,nsize_high-1
       nelm_len_nrank(n)=nblk+1
       nelm_pos_nrank(n)=ipos
       ipos=ipos+nelm_len_nrank(n)
    END DO
    DO n=nsize_high,nsize-1
       nelm_len_nrank(n)=nblk
       nelm_pos_nrank(n)=ipos
       ipos=ipos+nelm_len_nrank(n)
    END DO
    IF(ipos.NE.nemax) THEN
       WRITE(6,'(A,2I8)') 'XX ne parallel error: nemax,ipos=',nemax,ipos
       STOP
    END IF

    ! --- initialize ---
  
    RW=2.D0*PI*RF*1.D6
    CIWE=CII*RW*EPS0

    DO NE=nelm_pos_nrank(nrank)+1,nelm_pos_nrank(nrank)+nelm_len_nrank(nrank)
       S=SELM(NE)

       ! --- calculate conductivity tensor ---

       CALL wf_dtensr(NE,DTENS)
       DO NS=1,NSMAX
          DO IN=1,3
             DO J=1,3
                DO I=1,3
                   CTENS(NS,IN,I,J)=-CIWE*DTENS(NS,IN,I,J)
                END DO
             END DO
          END DO
       END DO

       CALL WFNODE(NE,R,Z)
       CALL wf_mutensr(NE,MU)
     
       CINT=0.d0

       DO NS=1,NSMAX
          SELECT CASE(MODELG)
          CASE(0,12)
             DO JJ=1,6
                DO II=1,6
                   DO K=1,3
                      DO J=1,3
                         DO I=1,3
                            CINT(NS,II,JJ)= CINT(NS,II,JJ) &
                                          +((MU(I,1,II)*CTENS(NS,J,1,1) &
                                            +MU(I,2,II)*CTENS(NS,J,2,1) &
                                            +MU(I,3,II)*CTENS(NS,J,3,1)) &
                                            *MU(K,1,JJ) &
                                           +(MU(I,1,II)*CTENS(NS,J,1,2) &
                                            +MU(I,2,II)*CTENS(NS,J,2,2) &
                                            +MU(I,3,II)*CTENS(NS,J,3,2)) &
                                            *MU(K,2,JJ) &
                                           +(MU(I,1,II)*CTENS(NS,J,1,3) &
                                            +MU(I,2,II)*CTENS(NS,J,2,3) &
                                            +MU(I,3,II)*CTENS(NS,J,3,3)) &
                                            *MU(K,3,JJ)) &
                                           *RR*S*AIF3(I,J,K)
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          CASE(1:10,13)
             DO JJ=1,6
                DO II=1,6
                   DO K=1,3
                      DO J=1,3
                         DO I=1,3
                            CINT(NS,II,JJ)= CINT(NS,II,JJ) &
                                          +((MU(I,1,II)*CTENS(NS,J,1,1) &
                                            +MU(I,2,II)*CTENS(NS,J,2,1) &
                                            +MU(I,3,II)*CTENS(NS,J,3,1)) &
                                            *MU(K,1,JJ) &
                                           +(MU(I,1,II)*CTENS(NS,J,1,2) &
                                            +MU(I,2,II)*CTENS(NS,J,2,2) &
                                            +MU(I,3,II)*CTENS(NS,J,3,2)) &
                                            *MU(K,2,JJ) &
                                           +(MU(I,1,II)*CTENS(NS,J,1,3) &
                                            +MU(I,2,II)*CTENS(NS,J,2,3) &
                                            +MU(I,3,II)*CTENS(NS,J,3,3)) &
                                            *MU(K,3,JJ)) &
                                           *R(J)*S*AIF3(I,J,K)
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END SELECT
       END DO

       DO I=1,3
          NSD=NSDELM(I,NE)
          IF(NSD.LT.0) THEN
             NSD=-NSD
             CE(I)=-CESD(NSD)
          ELSE
             CE(I)=CESD(NSD)
          END IF
       END DO
       DO I=1,3
          NN=NDELM(I,NE)
          CE(I+3)=CEND(NN)
       END DO

       DO NS=1,NSMAX
          PABS_ns_nelm(NS,NE)=0.d0
          DO JJ=1,6
             DO II=1,6
                PABS_ns_nelm(NS,NE)=PABS_ns_nelm(NS,NE) &
                            +0.5d0*real(CONJG(CE(II))*CINT(NS,II,JJ)*CE(JJ))
             END DO
          END DO
       END DO

    END DO

    nelm1=nelm_pos_nrank(nrank)+1
    nelm2=nelm_pos_nrank(nrank)+nelm_len_nrank(nrank)
    ndata=nelm_len_nrank(nrank)
    ALLOCATE(rdata(ndata),rdata_tot(nemax))
    DO ns=1,nsmax
       rdata(1:ndata)=pabs_ns_nelm(ns,nelm1:nelm2)
       CALL mtx_allgatherv_real8(rdata,ndata,rdata_tot,nemax, &
            nelm_len_nrank,nelm_pos_nrank)
       pabs_ns_nelm(ns,1:nemax)=rdata_tot(1:nemax)
    END DO

    DO NS=1,NSMAX
       PABST(NS)=0.d0
       DO NE=1,NEMAX
          PABST(NS)=PABST(NS)+PABS_ns_nelm(NS,NE)
       END DO
    END DO

    PABSTT=0.D0
    DO NS=1,NSMAX
       PABSTT=PABSTT+PABST(NS)
    END DO

    RETURN
  END SUBROUTINE wf_pwrabs

  !     ******* POWER RADIATION *******

  SUBROUTINE wf_pwrrad

    USE wfcomm
    IMPLICIT NONE

    INTEGER    :: NE,NA,I,NN,IV
    INTEGER    :: IJ,NSD
    REAL(rkind)    :: PHASE,RW,LIF(3),A(3),B(3),C(3)
    REAL(rkind)    :: R1,R2,Z1,Z2,R21,Z21,MU(3,3,6)
    COMPLEX(rkind) :: CE(6),CJ(3),CVJ

    ! --- initialize ---

    RW=2.D0*PI*RF*1.D6

    DO NA=1,NAMAX
       PHASE =APH(NA)*PI/180.D0
       CVJ=AJ(NA)*EXP(CII*(PHASE))
       CIMP(NA)=(0.d0,0.d0)
       DO IJ=2,JNUM(NA)
          NE=JELMT(IJ,NA)
          R1=RJ(IJ-1,NA)
          Z1=ZJ(IJ-1,NA)
          R2=RJ(IJ,NA)
          Z2=ZJ(IJ,NA)
          R21=R2-R1
          Z21=Z2-Z1
          CJ(1)=CVJ*R21
          CJ(2)=(0.d0,0.d0)
          CJ(3)=CVJ*Z21

          CALL WFABC(NE,A,B,C)

          SELECT CASE(MODELG)
          CASE(0,12)
             DO I=1,3
                LIF(I)= A(I)*RR &
                       +B(I)*RR*(R2+R1)/2.d0 &
                       +C(I)*RR*(Z1+Z2)/2.d0           
             END DO
          CASE(1:10,13)
             DO I=1,3
                LIF(I)= A(I)*(R1+R2)/2.d0 &
                       +B(I)*(R2**2+R1*R2+R1**2)/3.d0 &
                       +C(I)*(R2*Z1+Z2*R1+2.d0*R2*Z2+2.d0*R1*Z1)/6.d0
             END DO
          END SELECT
          DO I=1,3
             NSD=NSDELM(I,NE)
             IF(NSD.LT.0) THEN
                NSD=-NSD
                CE(I)=-CESD(NSD)
             ELSE
                CE(I)=CESD(NSD)
             END IF
          END DO
          DO I=1,3
             NN=NDELM(I,NE)
             CE(I+3)=CEND(NN)
          END DO

          CALL wf_mutensr(NE,MU)

          DO I=1,3
             DO IV=1,6
                CIMP(NA)=CIMP(NA) &
                         -0.5d0*LIF(I)*CONJG(CE(IV)) &
                                      *( MU(I,1,IV)*CJ(1) &
                                        +MU(I,2,IV)*CJ(2) &
                                        +MU(I,3,IV)*CJ(3))
             END DO
          END DO
       END DO
    END DO

    CTIMP=(0.d0,0.d0)

    DO NA=1,NAMAX
       CTIMP=CTIMP+CIMP(NA)
    END DO

    RETURN
  END SUBROUTINE wf_pwrrad

  !     ******* OUTPUT FIELD DATA *******

  SUBROUTINE wf_lpefld

    USE wfcomm
    IMPLICIT NONE
    INTEGER:: NS,NA

    IF(NPRINT.LT.1) RETURN
    IF(nrank.NE.0) RETURN
    
    WRITE(6,120) DBLE(CTIMP),PABSTT
120 FORMAT(1H ,'RADIATED POWER =',1PE12.4/ &
           1H ,'ABSORBED POWER =',1PE12.4)

    DO NS=1,NSMAX
       WRITE(6,126) NS,PABST(NS)
126    FORMAT(1H ,'      PABS(',I2,') =',1PE12.4)
    END DO

    WRITE(6,130)
130 FORMAT(1H ,' I JNUM', '  AJ(I)','  APH(I)','  AWD(I)', &
            ' APOS(I)',' XJ(I)','  YJ(I)', &
            8X,'LOADING IMP.[ohm]')
    DO NA=1,NAMAX
       WRITE(6,140) NA,JNUM(NA),AJ(NA),APH(NA),AWD(NA),APOS(NA), &
                               RJ(1,NA),ZJ(1,NA),CIMP(NA)
140    FORMAT(1H ,I2,I3,0PF8.2,F8.2,1X,4F7.4,2X,'(',1P2E12.4,')')
    END DO

    IF(NPRINT.LT.2) RETURN

    ! field output

    RETURN
  END SUBROUTINE wf_lpefld
  
  !     ******* OUTPUT ELEMENT DATA *******

  SUBROUTINE wf_lpelmt

    USE wfcomm
    IMPLICIT NONE

    INTEGER :: I,J,NA

    IF(NPRINT.LT.3) RETURN
    IF(nrank.NE.0) RETURN
     
    WRITE(6,110) NNMAX
110 FORMAT(/' ','NODE DATA     : #### NNMAX =',I5,' ####'/ &
                 ' ',2('  NNMAX',' KANOD', &
                 9X,'R',14X,'Z',9X))
    WRITE(6,115) (I,KANOD(I),RNODE(I),ZNODE(I), &
                  I=1,NNMAX)
115 FORMAT((' ',2(2I6,2X,1P2E15.7,2X)))
  
    WRITE(6,120) NEMAX,(I,(NDELM(J,I),J=1,3),I=1,NEMAX)
120 FORMAT(/' ','ELEMENT DATA  : #### NEMAX =',I5,' ####'/ &
                (' ',4(I6,'(',3I5,')',2X)))
  
    WRITE(6,125) NEMAX,(I,(NSDELM(J,I),J=1,3),I=1,NEMAX)
125 FORMAT(/' ','SIDE    DATA  : #### NEMAX =',I5,' ####'/ &
                (' ',2(I8,'(',3I8,')',2X)))
  
    DO NA=1,NAMAX
       WRITE(6,140) NA,JNUM0(NA)
140    FORMAT(/' ','ORIGINAL ANTENNA DATA : NA =',I5,' JNUM0 =',I5/ &
                    ' ',2('  NO.',13X,' RJ0',11X,' ZJ0',6X))
       WRITE(6,150) (I,RJ0(I,NA),ZJ0(I,NA),I=1,JNUM0(NA))
150    FORMAT((' ',2(I5,8X,1P2E15.7)))
     
       WRITE(6,154) NA,JNUM(NA)
154    FORMAT(/' ','MODIFIED ANTENNA DATA : NA =',I5,' JNUM  =',I5/ &
                     ' ',2('  NO.',' JELM',8X,' JR ',11X,' JZ ',6X))
       WRITE(6,156) (I,JELMT(I,NA),RJ(I,NA),ZJ(I,NA),I=1,JNUM(NA))
156    FORMAT((' ',2(2I5,3X,1P2E15.7)))
    ENDDO
  
    RETURN
  END SUBROUTINE wf_lpelmt
END MODULE wfwaveK
