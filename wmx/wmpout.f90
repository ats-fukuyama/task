! wmpout.f90

MODULE wmpout

  PRIVATE
  PUBLIC wm_pout,wm_pout_init,wm_pout_sum

CONTAINS

!     ****** DISPLAY OUTPUT DATA ******

  SUBROUTINE wm_pout

    USE wmcomm
    IMPLICIT NONE
    INTEGER:: NS,ND,NDX,NN,MD,MDX,MM,NR

    IF(NPRINT.LT.1) RETURN

    IF(NRANK.EQ.0) THEN
       WRITE(6,601) CPRAD,PABSTT,PCURT
       WRITE(6,602) (PABST(NS),NS=1,NSMAX)
    ENDIF

    IF(NPRINT.LT.2) RETURN

    IF(NRANK.EQ.0) THEN
       WRITE(6,*) '   MD   IMP                        PABSKT'
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          NN=NPH0+NHC*ND
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             MM=NTH0+MD
             WRITE(6,603) NN,MM,CPRADK(MDX,NDX), &
                          (PABSKT(MDX,NDX,NS),NS=1,NSMAX)
          ENDDO
       ENDDO
    ENDIF

    IF(NPRINT.LT.4) RETURN

    IF(NRANK.EQ.0) THEN
       DO NR=1,NRMAX+1
          IF(NRANK.EQ.0) &
               WRITE(6,'(A,I3,6ES10.2)') 'NR,E1,E2,E3=',NR, &
               CEFLD(1,1,1,NR),CEFLD(2,1,1,NR),CEFLD(3,1,1,NR)
       ENDDO
    END IF

    RETURN

601 FORMAT(' ',5X,3X,'RANT=',ES12.4,10X,'LANT=',ES12.4/ &
           ' ',5X,3X,'PABS=',ES12.4,10X,'IDRV=',ES12.4)
602 FORMAT(' ',5X,3X,27X,'PABS=',3ES12.4)
603 FORMAT(' ',I5,I5,3X,2ES12.4,3X,3ES12.4)
604 FORMAT(' ',8X,3F8.4,24X,2F8.4)
605 FORMAT(' ',9F8.4)

  END SUBROUTINE wm_pout

! --- initial setup to sum up em field and pabs ---

  SUBROUTINE wm_pout_init

    USE wmcomm
    IMPLICIT NONE
    INTEGER:: NS,NR,NHH,NTH,ND,NDX,MD,MDX

    DO NS=1,NSMAX
       DO NR=1,NRMAX+1
          DO ND=NDMIN_F,NDMAX_F
             NDX=ND-NDMIN_F+1
             DO MD=MDMIN_F,MDMAX_F
                MDX=MD-MDMIN_F+1
                PABSK(MDX,NDX,NR,NS)=0d0
                CPABSK(MDX,NDX,NR,NS)=0d0
             END DO
          END DO
       END DO
    END DO
    
    DO NS=1,NSMAX
       DO NR=1,NRMAX
          DO NHH=1,NHHMAX_F
             DO NTH=1,NTHMAX_F
                PABS(NTH,NHH,NR,NS)=0d0
                CPABS(NTH,NHH,NR,NS)=0d0
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO NR=1,NRMAX
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             PCUR(NTH,NHH,NR)=0d0
          ENDDO
       ENDDO
    ENDDO
    
    DO NR=1,NRMAX+1
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             CEFLDK(1,MDX,NDX,NR)=0d0
             CEFLDK(2,MDX,NDX,NR)=0d0
             CEFLDK(3,MDX,NDX,NR)=0d0
             CBFLDK(1,MDX,NDX,NR)=0d0
             CBFLDK(2,MDX,NDX,NR)=0d0
             CBFLDK(3,MDX,NDX,NR)=0d0
          ENDDO
       ENDDO
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             CEFLD(1,NTH,NHH,NR)=0d0
             CEFLD(2,NTH,NHH,NR)=0d0
             CEFLD(3,NTH,NHH,NR)=0d0
             CBFLD(1,NTH,NHH,NR)=0d0
             CBFLD(2,NTH,NHH,NR)=0d0
             CBFLD(3,NTH,NHH,NR)=0d0
             CEN(1,NTH,NHH,NR)  =0d0
             CEN(2,NTH,NHH,NR)  =0d0
             CEN(3,NTH,NHH,NR)  =0d0
             CEP(1,NTH,NHH,NR)  =0d0
             CEP(2,NTH,NHH,NR)  =0d0
             CEP(3,NTH,NHH,NR)  =0d0
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE wm_pout_init

! --- post processing to sum up pabs ---
  
  SUBROUTINE wm_pout_sum

    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),ALLOCATABLE:: DS(:),DSS(:,:,:)
    COMPLEX(rkind):: cw
    REAL(rkind):: DTH,DHH,DRHO,DPSIPDRHO,DSSS,FACT,FACTSQ
    INTEGER:: NR,NHH,NTH,NS,NHHF,NTHF,ND,NDX,MD,MDX
    
    ALLOCATE(DS(nrmax),DSS(nthmax,nhhmax,nrmax))
    
    CW=2.D0*PI*CMPLX(RF,RFI)*1.D6
    DTH=2.D0*PI/DBLE(NTHMAX)
    DHH=2.D0*PI/DBLE(NHHMAX*NHC)

    DO NR=1,NRMAX
       PCURR(NR)=0.D0
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH*DHH
          ENDDO
       ENDDO
    ENDDO

    DO NS=1,NSMAX
       DO NR=1,NRMAX
          PABSR(NR,NS)=0.D0
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                PABSR(NR,NS)=PABSR(NR,NS)+PABS(NTH,NHH,NR,NS)*DTH*DHH
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!    DO NR=1,NRMAX
!       WRITE(6,'(A,I6,64ES12.4)') 'P:',NR,PABSR(NR,1),PABSR(NR,2), &
!            PABS(1,1,NR,1),PABS(1,1,NR,2),PABSK(1,1,NR,1),PABSK(1,1,NR,1)
!       WRITE(6,'(A,I6,4ES12.4)') 'CPK:',NR, &
!            CPABSK(1,1,NR,1),CPABSK(1,1,NR,2)
!       WRITE(6,'(A,I6,4ES12.4)') 'CPR:',NR, &
!            CPABS(1,1,NR,1),CPABS(1,1,NR,2)
!    END DO

    PCURT=0.D0
    DO NR=1,NRMAX
       PCURT=PCURT+PCURR(NR)
    ENDDO

    DO NS=1,NSMAX
       PABST(NS)=0.D0
       DO NR=1,NRMAX
          PABST(NS)=PABST(NS)+PABSR(NR,NS)
       ENDDO
    ENDDO

    PABSTT=0.D0
    DO NS=1,NSMAX
       PABSTT=PABSTT+PABST(NS)
    ENDDO

! --- evaluate factor for given RF power PRFIN ---
    
    FACT=1.D0
    IF(PRFIN.GT.0.D0.AND.PABSTT.GT.0.D0) THEN
       FACT=PRFIN/PABSTT
    ENDIF
    FACTSQ=SQRT(FACT)

! --- calculate section area ---
    
    NR=1
       DS(NR)=0.D0
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             DSS(NTH,NHH,NR)=0.D0
          ENDDO
       ENDDO
    DO NR=2,NRMAX
       DS(NR)=0.D0
       DRHO=0.5D0*(XRHO(NR+1)-XRHO(NR-1))
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             NHHF=(NHH-1)*FACTOR_NHH+1
             NTHF=(NTH-1)*FACTOR_NTH+1
             IF(MODELG.EQ.3) THEN
                DPSIPDRHO=2.D0*PSITA*XRHO(NR)/QPS(NR)
             ELSE
                DPSIPDRHO=2.D0*PSIPA*XRHO(NR)
             ENDIF
             DSSS=DPSIPDRHO*DRHO*RJ(NTHF,NHHF,NR)
             DSS(NTH,NHH,NR)=1.D0/DSSS
             DS(NR)=DS(NR)+DSSS*DTH*DHH
          ENDDO
       ENDDO
       DS(NR)=1.D0/DS(NR)
    ENDDO

! --- normalize absorbed power ---
    
    PABSTT=FACT*PABSTT
    DO NS=1,NSMAX
       PABST(NS)=FACT*PABST(NS)
       DO NR=1,NRMAX
          PABSR(NR,NS)=FACT*PABSR(NR,NS)*DS(NR)
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                PABS(NTH,NHH,NR,NS)=FACT*PABS(NTH,NHH,NR,NS) &
                                   *DSS(NTH,NHH,NR)
                CPABS(NTH,NHH,NR,NS)=FACT*CPABS(NTH,NHH,NR,NS) &
                                   *DSS(NTH,NHH,NR)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

! --- normalize electromagnetic field ---
    
    DO NR=1,NRMAX+1
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             CEFLDK(1,MDX,NDX,NR)=FACTSQ*CEFLDK(1,MDX,NDX,NR)
             CEFLDK(2,MDX,NDX,NR)=FACTSQ*CEFLDK(2,MDX,NDX,NR)
             CEFLDK(3,MDX,NDX,NR)=FACTSQ*CEFLDK(3,MDX,NDX,NR)
             CBFLDK(1,MDX,NDX,NR)=FACTSQ*CBFLDK(1,MDX,NDX,NR)
             CBFLDK(2,MDX,NDX,NR)=FACTSQ*CBFLDK(2,MDX,NDX,NR)
             CBFLDK(3,MDX,NDX,NR)=FACTSQ*CBFLDK(3,MDX,NDX,NR)
          ENDDO
       ENDDO
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             CEFLD(1,NTH,NHH,NR)=FACTSQ*CEFLD(1,NTH,NHH,NR)
             CEFLD(2,NTH,NHH,NR)=FACTSQ*CEFLD(2,NTH,NHH,NR)
             CEFLD(3,NTH,NHH,NR)=FACTSQ*CEFLD(3,NTH,NHH,NR)
             CBFLD(1,NTH,NHH,NR)=FACTSQ*CBFLD(1,NTH,NHH,NR)
             CBFLD(2,NTH,NHH,NR)=FACTSQ*CBFLD(2,NTH,NHH,NR)
             CBFLD(3,NTH,NHH,NR)=FACTSQ*CBFLD(3,NTH,NHH,NR)
             CEN(1,NTH,NHH,NR)  =FACTSQ*CEN(1,NTH,NHH,NR)
             CEN(2,NTH,NHH,NR)  =FACTSQ*CEN(2,NTH,NHH,NR)
             CEN(3,NTH,NHH,NR)  =FACTSQ*CEN(3,NTH,NHH,NR)
             CEP(1,NTH,NHH,NR)  =FACTSQ*CEP(1,NTH,NHH,NR)
             CEP(2,NTH,NHH,NR)  =FACTSQ*CEP(2,NTH,NHH,NR)
             CEP(3,NTH,NHH,NR)  =FACTSQ*CEP(3,NTH,NHH,NR)
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE wm_pout_sum
END MODULE wmpout

