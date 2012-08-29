MODULE wmtest
  USE wmcont

CONTAINS

  SUBROUTINE wm_test(CEFLD3D,MDM,NPHM,NRM,NTHMAX_IN,NPHMAX_IN,NRMAX_IN, &
                     RGMIN,RGMAX,RAXIS,NCONT)

    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM
    REAL(8),INTENT(IN):: RGMIN,RGMAX,RAXIS
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NCONT
    INTEGER:: MODE,NPH

1   WRITE(6,'(A)') '## INPUT: 1-6:3D, 7:2D, 0:EXIT'
    READ(5,*,ERR=1,END=9000) MODE

    SELECT CASE(MODE)
    CASE(1:6)
       CALL wmg3d2(CEFLD3D,MDM,NPHM,NRM,NTHMAX_IN,NPHMAX_IN,NRMAX_IN, &
                   RGMIN,RGMAX,RAXIS,MODE,NCONT)
    CASE(7)
2      WRITE(6,'(A,I,A)') '## INPUT: NPH=1..',NPHMAX_IN,' 0 for exit'
       READ(5,*,ERR=2,END=1) NPH
       IF(NPH.EQ.0) GO TO 1
3      WRITE(6,'(A)') '## INPUT: 1-6:2D, 0:EXIT'
       READ(5,*,ERR=3,END=2) MODE
       IF(MODE.EQ.0) GO TO 2
       CALL wmg3d1(CEFLD3D,MDM,NPHM,NRM,NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NPH,MODE)
       GO TO 3
    CASE(0)
       GOTO 9000
    END SELECT
    GO TO 1
9000    RETURN
  END SUBROUTINE wm_test

  SUBROUTINE wmg3d2(CEFLD3D,MDM,NPHM,NRM,NTHMAX_IN,NPHMAX_IN,NRMAX_IN, &
                    RGMIN,RGMAX,RAXIS,MODE,NCONT)

    USE bpsd_constants,ONLY:PI,CI
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM
    REAL(8),INTENT(IN):: RGMIN,RGMAX,RAXIS
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,MODE,NCONT
    REAL(4),DIMENSION(:,:),ALLOCATABLE:: Z
    REAL(4),DIMENSION(:),ALLOCATABLE:: R
    REAL(4),DIMENSION(:),ALLOCATABLE:: PH
    REAL(4),DIMENSION(:),ALLOCATABLE:: ZL,WLN
    INTEGER,DIMENSION(:),ALLOCATABLE:: ILN
    REAL(4),DIMENSION(:,:),ALLOCATABLE:: RGB
    INTEGER:: nrmax,nphmax,nr,nth,nph,nn,nglmax,ngl,ispl
    INTEGER:: nr_in,nth_in,nph_in
    REAL(4):: rmin,rmax,rr,rb,phmin,phmax,dph,zmin,zmax
    REAL(4):: rmin1,rmax1,rstep,zmin1,zmax1,zstep,fsign
    COMPLEX(8):: value
    INTERFACE 
       FUNCTION NGULEN(R)
         REAL(4),INTENT(IN):: R
         INTEGER:: NGULEN
       END FUNCTION NGULEN
       FUNCTION GUCLIP(D)
         REAL(8),INTENT(IN):: D
         REAL(4):: GUCLIP
       END FUNCTION GUCLIP
    END INTERFACE

    NPHMAX=128
    NGLMAX=NCONT
    NRMAX=2*(NRMAX_IN+1)

    ALLOCATE(Z(NRMAX,NPHMAX),R(NRMAX),PH(NPHMAX))
    ALLOCATE(ZL(NGLMAX),WLN(NGLMAX),ILN(NGLMAX),RGB(3,NGLMAX))

    RMIN=SNGL(RGMIN)
    RMAX=SNGL(RGMAX)
    RR=SNGL(RAXIS)
    RB=RR-RMIN
    DO NR=1,NRMAX_IN+1
       R(NR)=RR-RB*REAL(NRMAX_IN+1-NR)/REAL(NRMAX_IN)
    END DO
    RB=RMAX-RR
    DO NR=NRMAX_IN+2,NRMAX
       R(NR)=RR+RB*REAL(NR-NRMAX_IN-2)/REAL(NRMAX_IN)
    END DO
    PHMIN=0.D0
    PHMAX=2.D0*PI
    DPH=(PHMAX-PHMIN)/(NPHMAX-1)
    DO NPH=1,NPHMAX
       PH(NPH)=PHMIN+DPH*(NPH-1)
    END DO

    DO NPH=1,NPHMAX
       DO NR=1,NRMAX
          IF(NR.LE.NRMAX_IN+1) THEN
             NR_IN=NRMAX_IN-NR+2
             NTH_IN=NTHMAX_IN/2+1
             fsign=-1.0
          ELSE
             NR_IN=NR-NRMAX_IN-1
             NTH_IN=1
             fsign=1.0
          ENDIF
          VALUE=(0.D0,0.D0)
          SELECT CASE(MODE)
          CASE(1,2)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(1,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(3,4)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+fsign*CEFLD3D(2,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(5,6)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(3,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          END SELECT
          SELECT CASE(MODE)
          CASE(1,3,5)
             Z(NR,NPH)=GUCLIP(REAL(VALUE))
          CASE(2,4,6)
             Z(NR,NPH)=GUCLIP(IMAG(VALUE))
          END SELECT
       END DO
    END DO

!    DO NR=1,NRMAX
!       WRITE(6,'(I5,1P2E12.4)') NR,R(NR),Z(NR,1)
!    END DO


    CALL pages
      CALL SETLNW(0.07)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
!
      CALL GMNMX2(Z,NRMAX,1,NRMAX,1,1,NPHMAX,1,ZMIN,ZMAX)
      CALL GQSCAL(ZMIN,ZMAX,ZMIN1,ZMAX1,ZSTEP)
      CALL GQSCAL(-RMAX,RMAX,RMIN1,RMAX1,RSTEP)
!
      zmax1=MAX(ABS(zmin1),ABS(zmax1))
      DO ngl=1,nglmax
         ZL(NGL)=-zmax1+2*zmax1*REAL(ngl-1)/REAL(NGLMAX-1)
         WLN(NGL)=0.0
         ILN(NGL)=0
         IF(ZL(NGL)>=0.D0) THEN
            RGB(1,NGL)=1.0
            RGB(2,NGL)=0.0
            RGB(3,NGL)=0.0
         ELSE
            RGB(1,NGL)=0.0
            RGB(2,NGL)=0.0
            RGB(3,NGL)=1.0
         ENDIF
      END DO
      ISPL=0

      CALL GDEFIN(2.,17.,2.,17.,-RMAX,RMAX,-RMAX,RMAX)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,RSTEP,0.,RSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,RSTEP*2,0.,RSTEP*2,NGULEN(2*RSTEP))

      CALL CONTG4X(Z,R,PH,NRMAX,NRMAX,NPHMAX,ZL,RGB,ILN,WLN,NGLMAX,ISPL)
    CALL pagee

    DEALLOCATE(Z,R,PH)
    DEALLOCATE(ZL,WLN,ILN,RGB)

    RETURN
  END SUBROUTINE wmg3d2

  SUBROUTINE wmg3d1(CEFLD3D,MDM,NPHM,NRM,NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NPH,MODE)
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NPH,MODE
    REAL(4),DIMENSION(NRM,MDM):: Z
    INTEGER:: NR,NTH
    CHARACTER(LEN=1):: K2,K3,K4
    INTERFACE 
       FUNCTION GUCLIP(D)
         REAL(8),INTENT(IN):: D
         REAL(4):: GUCLIP
       END FUNCTION GUCLIP
    END INTERFACE

    K2='E'
    SELECT CASE(MODE)
    CASE(1)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CEFLD3D(1,NTH,NPH,NR)))
          END DO
       END DO
       K3='R'
       K4='R'
    CASE(2)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(IMAG(CEFLD3D(1,NTH,NPH,NR)))
          END DO
       END DO
       K3='R'
       K4='I'
    CASE(3)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CEFLD3D(2,NTH,NPH,NR)))
          END DO
       END DO
       K3='T'
       K4='R'
    CASE(4)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(IMAG(CEFLD3D(2,NTH,NPH,NR)))
          END DO
       END DO
       K3='T'
       K4='I'
    CASE(5)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CEFLD3D(3,NTH,NPH,NR)))
          END DO
       END DO
       K3='Z'
       K4='R'
    CASE(6)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(IMAG(CEFLD3D(3,NTH,NPH,NR)))
          END DO
       END DO
       K3='Z'
       K4='I'
    END SELECT

    CALL pages

      CALL SETCHS(0.3,0.0)
      CALL WMGXEQ(Z,1,K2,K3)

      CALL MOVE(20.0,17.5)
      IF(K2.EQ.'E') CALL TEXT('E ',2)
      IF(K2.EQ.'B') CALL TEXT('B ',2)
      IF(K2.EQ.'P') CALL TEXT('P ',2)
      IF(K2.EQ.'J') CALL TEXT('J ',2)
      CALL MOVE(20.0,17.1)
      IF(K2.NE.'P') THEN
         IF(K3.EQ.'R') CALL TEXT('r ',2)
         IF(K3.EQ.'T') CALL TEXT('theta ',6)
         IF(K3.EQ.'Z') CALL TEXT('phi ',4)
         IF(K3.EQ.'S') CALL TEXT('s ',2)
         IF(K3.EQ.'H') CALL TEXT('h ',2)
         IF(K3.EQ.'B') CALL TEXT('b ',2)
         IF(K3.EQ.'+') CALL TEXT('+ ',2)
         IF(K3.EQ.'-') CALL TEXT('- ',2)
         IF(K3.EQ.'P') CALL TEXT('P ',2)
      ELSE
         IF(K3.EQ.'E') CALL TEXT('e ',2)
         IF(K3.EQ.'D') CALL TEXT('D ',2)
         IF(K3.EQ.'T') CALL TEXT('T ',2)
         IF(K3.EQ.'1') CALL TEXT('1 ',2)
         IF(K3.EQ.'2') CALL TEXT('2 ',2)
         IF(K3.EQ.'3') CALL TEXT('3 ',2)
         IF(K3.EQ.'4') CALL TEXT('4 ',2)
         IF(K3.EQ.'5') CALL TEXT('5 ',2)
         IF(K3.EQ.'6') CALL TEXT('6 ',2)
      ENDIF
      CALL MOVE(20.0,16.7)
      IF(K4.EQ.'R') CALL TEXT('Real ',5)
      IF(K4.EQ.'I') CALL TEXT('Imag ',5)
      IF(K4.EQ.'A') CALL TEXT('Abs ',4)
!
      CALL WMGPRM('C',K3,0,0,0,0)
    CALL pagee
    RETURN
  END SUBROUTINE wmg3d1
END MODULE wmtest
