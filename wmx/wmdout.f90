! wmdout.f90

MODULE wmdout

  PRIVATE
  PUBLIC wm_dout

CONTAINS

!     ****** OUTPUT CALCULATION DATA AND PARAMETER ******

  SUBROUTINE wm_dout(IERR)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    LOGICAL:: LEX
    CHARACTER(LEN=6):: KNAMEO,KNAMES,KNAMET
    CHARACTER(LEN=13):: KNAMEP
    CHARACTER(LEN=10):: KNAMER
    INTEGER:: IST,NS,NA,NR,NHH,NTH,NSU,NSW

!     ###### WRITTING THE DATA OF ABSORPED POWER ######

7000 CONTINUE
    WRITE(6,*) 'INPUT DATA NAME OF PARAMETER : A6'
    READ(5,'(A6)',ERR=7000,END=9000) KNAMES
8000 CONTINUE
    WRITE(6,*) 'INPUT DATA NAME OF ABSORPED POWER : A6'
    READ(5,'(A6)',ERR=8000,END=9000) KNAMET

    KNAMER='gout'//KNAMES

    INQUIRE(FILE=KNAMER,EXIST=LEX)
    IF(LEX) THEN
       WRITE(6,*) 'XX FILE(',KNAMER,' ) HAS ALREADY EXISTED'
       GO TO 300
    ELSE
       OPEN(11,FILE=KNAMER,IOSTAT=IST,STATUS='NEW',ERR=100)
       GO TO 200
100    WRITE(6,*) 'XX FILE OPEN ERROR gout :  IOSTAT=',IST
       GO TO 9000
    ENDIF

200 CONTINUE
    WRITE(11,*) NSMAX,NAMAX,NRMAX,NTHMAX,NHHMAX
    WRITE(11,*) NHC
    WRITE(11,*) MODELG,MODELJ,MODELA,MODELK,MODELM
    WRITE(11,*) BB,RR,RA,RB,Q0,QA,RKAP,RDLT
    WRITE(11,*) PROFN1,PROFN2,PROFT1,PROFT2
    WRITE(11,*) ZEFF
    WRITE(11,*) PNA,PNAL,PTA
    WRITE(11,*) RF,RFI
    WRITE(11,*) RD
    DO NS=1,NSMAX
       WRITE(11,*) PA(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(11,*) PZ(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(11,*) PN(NS),PNS(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(11,*) PZCL(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(11,*) PTPR(NS),PTPP(NS),PTS(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(11,*) MODELP(NS)
    ENDDO
    DO NA=1,NAMAX
       WRITE(11,*) AJ(NA),APH(NA)
    ENDDO
    DO NA=1,NAMAX
       WRITE(11,*) THJ1(NA),THJ2(NA)
    ENDDO
    DO NA=1,NAMAX
       WRITE(11,*) PHJ1(NA),PHJ2(NA)
    ENDDO
    DO NA=1,NAMAX
       WRITE(11,*) BETAJ(NA)
    ENDDO
    WRITE(11,*) RGMAX,RGMIN
    WRITE(11,*) ZGMAX,ZGMIN
    WRITE(11,*) NSUMAX,NSWMAX

    DO NR=1,NRMAX+1
       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             WRITE(11,*) RPST(NTH,NHH,NR),ZPST(NTH,NHH,NR),BPST(NTH,NHH,NR)
          ENDDO
       ENDDO
    ENDDO

    DO NSU=1,NSUMAX
       DO NHH=1,NHHMAX
          WRITE(11,*) RSU(NSU,NHH),ZSU(NSU,NHH)
       ENDDO
    ENDDO

    DO NSW=1,NSWMAX
       DO NHH=1,NHHMAX
          WRITE(11,*) RSW(NSW,NHH),ZSW(NSW,NHH)
       ENDDO
    ENDDO

    CLOSE(11)
!         
300 CONTINUE
    CALL NAMEK(KNAMEO,NPH0)
    KNAMEP=KNAMEO//KNAMET

    INQUIRE(FILE=KNAMEP,EXIST=LEX)
    IF(LEX) THEN
       WRITE(6,*) 'XX FILE (',KNAMEP,') HAS ALREADY EXISTED'
       GO TO 9000
    ELSE
       OPEN(10,FILE=KNAMEP,IOSTAT=IST,STATUS='NEW',ERR=400)
       GO TO 500
400    WRITE(6,*) 'XX FILE OPEN ERROR pdata : IOSTAT=',IST
       GO TO 9000
    ENDIF

500 CONTINUE
    WRITE(10,*) NSMAX,NAMAX,NRMAX,NTHMAX,NHHMAX
    WRITE(10,*) NHC
    WRITE(10,*) MODELG,MODELJ,MODELA,MODELK,MODELM
    WRITE(10,*) BB,RR,RA,RB,Q0,QA,RKAP,RDLT
    WRITE(10,*) PROFN1,PROFN2,PROFT1,PROFT2
    WRITE(10,*) ZEFF
    WRITE(10,*) PNA,PNAL,PTA
    WRITE(10,*) RF,RFI
    WRITE(10,*) RD

    DO NS=1,NSMAX
       WRITE(10,*) PA(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(10,*) PZ(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(10,*) PN(NS),PNS(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(10,*) PZCL(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(10,*) PTPR(NS),PTPP(NS),PTS(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(10,*) MODELP(NS)
    ENDDO
    DO NA=1,NAMAX
       WRITE(10,*) AJ(NA),APH(NA)
    ENDDO
    DO NA=1,NAMAX
       WRITE(10,*) THJ1(NA),THJ2(NA)
    ENDDO
    DO NA=1,NAMAX
       WRITE(10,*) PHJ1(NA),PHJ2(NA)
    ENDDO
    DO NA=1,NAMAX
       WRITE(10,*) BETAJ(NA)
    ENDDO
!     
    DO NS=1,NSMAX
       DO NR=1,NRMAX
          DO NHH=1,NHHMAX
             DO NTH=1,NTHMAX
                WRITE(10,*) PABS(NTH,NHH,NR,NS)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    CLOSE(10)

    RETURN

9000 CONTINUE
    IERR=1
    RETURN

  END SUBROUTINE wm_dout

!     ****** NAME DATA FILE OF ABSORBED POWER ******

  SUBROUTINE NAMEK(KNAMEP,NPH0)
    IMPLICIT NONE
    CHARACTER(LEN=6),INTENT(OUT):: KNAMEP
    INTEGER,INTENT(IN):: NPH0
    
    IF(NPH0.GE.0) THEN
       KNAMEP='pd+000'
       SELECT CASE(NPH0)
       CASE(0:9)
          WRITE(KNAMEP(6:6),'(I1)') NPH0
       CASE(10:99)
          WRITE(KNAMEP(5:6),'(I2)') NPH0
       CASE(100:999)
          WRITE(KNAMEP(4:6),'(I3)') NPH0
       END SELECT
    ELSE
       KNAMEP='pd-000'
       SELECT CASE(-NPH0)
       CASE(1:9)
          WRITE(KNAMEP(6:6),'(I1)') -NPH0
       CASE(10:99)
          WRITE(KNAMEP(5:6),'(I2)') -NPH0
       CASE(100:999)
          WRITE(KNAMEP(4:6),'(I3)') -NPH0
       END SELECT
    ENDIF

    RETURN
  END SUBROUTINE NAMEK
END MODULE wmdout
