MODULE wmtest
  USE wmcont

CONTAINS

  SUBROUTINE wm_test(CEFLD3D,CPABS3D,PABST3D,PABSTT3D, &
                     MDM,NPHM,NRM,NSM, &
                     NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX, &
                     RGMIN,RGMAX,RAXIS,NCONT,NGRAPH)
    USE grf1d_mod
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    COMPLEX(8),DIMENSION(MDM,NPHM,NRM,NSM),INTENT(IN):: CPABS3D
    REAL(8),DIMENSION(NPHM,NSM),INTENT(IN):: PABST3D
    REAL(8),DIMENSION(NPHM),INTENT(IN):: PABSTT3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM,NSM
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX
    REAL(8),INTENT(IN):: RGMIN,RGMAX,RAXIS
    INTEGER,INTENT(IN):: NCONT,NGRAPH
    INTEGER:: MODE0,MODE1,NN,NPH
    REAL(4),DIMENSION(NPHMAX_IN):: GX
    REAL(4),DIMENSION(NPHMAX_IN,1):: GY
    INTERFACE 
       FUNCTION GUCLIP(D)
         REAL(8),INTENT(IN):: D
         REAL(4):: GUCLIP
       END FUNCTION GUCLIP
    END INTERFACE

1   WRITE(6,'(A)') '## INPUT: 1:3D, 2:2D, 3:pabs, 4:Ewg, 0:EXIT'
    READ(5,*,ERR=1,END=9000) MODE0

    SELECT CASE(MODE0)
    CASE(1)
2      WRITE(6,'(A)') '## INPUT: 1-9:E, 11..:PABS 0:EXIT'
       READ(5,*,ERR=2,END=1) MODE1
       IF(MODE1.EQ.0) GO TO 1
       CALL wmg3d1(CEFLD3D,CPABS3D,MDM,NPHM,NRM,NSM, &
                   NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX, &
                   RGMIN,RGMAX,RAXIS,MODE1,NCONT,NGRAPH)
       GO TO 2
    CASE(2)
3      WRITE(6,'(A,I4,A,I4)') '## INPUT: MODE:1-9:E 11..:PABS 0:EXIT, NN:', &
            -NPHMAX_IN/2,'..',NPHMAX_IN/2-1
       READ(5,*,ERR=3,END=1) MODE1,NN
       IF(MODE1.EQ.0) GO TO 1
       CALL wmg3d2(CEFLD3D,CPABS3D,MDM,NPHM,NRM,NSM, &
                   NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX,NN,MODE1,NGRAPH)
       GO TO 3
    CASE(3)
       DO NPH=1,NPHMAX_IN
          GX(NPH)=REAL(NPH-NPHMAX_IN/2-1)
          GY(NPH,1)=GUCLIP(PABSTT3D(NPH))
       END DO
       CALL pages
       CALL grf1d(0,GX,GY,NPHMAX_IN,NPHMAX_IN,1,'@Pabs vs n_phi@',0)
       CALL pagee
    CASE(4)
4      WRITE(6,'(A)') '## INPUT: 1-9:E, 11..:PABS 0:EXIT'
       READ(5,*,ERR=4,END=1) MODE1
       IF(MODE1.EQ.0) GO TO 1
       CALL wmg3d4(CEFLD3D,CPABS3D,MDM,NPHM,NRM,NSM, &
                   NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX, &
                   RGMIN,RGMAX,RAXIS,MODE1,NCONT)
       GO TO 4
    CASE(0)
       GOTO 9000
    END SELECT
    GO TO 1
9000    RETURN
  END SUBROUTINE wm_test

  SUBROUTINE wmg3d1(CEFLD3D,CPABS3D,MDM,NPHM,NRM,NSM, &
                    NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX, &
                    RGMIN,RGMAX,RAXIS,MODE,NCONT,NGRAPH)

    USE bpsd_constants,ONLY:PI,CI
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    COMPLEX(8),DIMENSION(MDM,NPHM,NRM,NSM),INTENT(IN):: CPABS3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM,NSM
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX
    REAL(8),INTENT(IN):: RGMIN,RGMAX,RAXIS
    INTEGER,INTENT(IN):: MODE,NCONT,NGRAPH
    REAL(4),DIMENSION(:,:),ALLOCATABLE:: Z
    REAL(4),DIMENSION(:),ALLOCATABLE:: R
    REAL(4),DIMENSION(:),ALLOCATABLE:: PH
    REAL(4),DIMENSION(:),ALLOCATABLE:: ZL,WLN
    INTEGER,DIMENSION(:),ALLOCATABLE:: ILN
    REAL(4),DIMENSION(:,:),ALLOCATABLE:: RGB
    INTEGER,DIMENSION(:,:),ALLOCATABLE:: KA(:,:,:)
    INTEGER:: nrmax,nphmax,nr,nth,nph,nn,nglmax,ngl,ispl
    INTEGER:: nr_in,nth_in,nph_in
    REAL(4):: rmin,rmax,rr,rb,phmin,phmax,dph,zmin,zmax
    REAL(4):: rmin1,rmax1,rstep,zmin1,zmax1,zstep,fsign,phl
    COMPLEX(8):: value
    INTEGER,parameter:: NRGBA=5
    REAL(4):: GRGBA(3,NRGBA) = reshape( &
                               (/ 0.0,0.0,1.0, &
                                  0.0,1.0,1.0, &
                                  1.0,1.0,1.0, &
                                  1.0,1.0,0.0, &
                                  1.0,0.0,0.0/), &
                                  shape(GRGBA))
    REAL(4):: GLA(NRGBA) = (/0.0,0.40,0.5,0.60,1.0/)
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
    NRMAX=2*NRMAX_IN+1

    ALLOCATE(Z(NRMAX,NPHMAX),R(NRMAX),PH(NPHMAX))
    ALLOCATE(ZL(NGLMAX),WLN(NGLMAX),ILN(NGLMAX),RGB(3,NGLMAX+1))
    ALLOCATE(KA(8,NRMAX,NPHMAX))

    RMIN=SNGL(RGMIN)
    RMAX=SNGL(RGMAX)
    RR=SNGL(RAXIS)
    RB=RR-RMIN
    DO NR=1,NRMAX_IN+1
       R(NR)=RR-RB*REAL(NRMAX_IN+1-NR)/REAL(NRMAX_IN)
    END DO
    RB=RMAX-RR
    DO NR=NRMAX_IN+2,NRMAX
       R(NR)=RR+RB*REAL(NR-NRMAX_IN-1)/REAL(NRMAX_IN)
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
             NR_IN=NR-NRMAX_IN
             NTH_IN=1
             fsign=1.0
          ENDIF
          VALUE=(0.D0,0.D0)
          SELECT CASE(MODE)
          CASE(1,2,3)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(1,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(4,5,6)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+fsign*CEFLD3D(2,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(7,8,9)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(3,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(11:16)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CPABS3D(NTH_IN,NPH_IN,NR_IN,MODE-10)
             END DO
          END SELECT
          SELECT CASE(MODE)
          CASE(1,4,7)
             Z(NR,NPH)=GUCLIP(REAL(VALUE))
          CASE(2,5,8)
             Z(NR,NPH)=GUCLIP(IMAG(VALUE))
          CASE(3,6,9)
             Z(NR,NPH)=GUCLIP(ABS(VALUE))
          CASE(11:16)
             Z(NR,NPH)=GUCLIP(REAL(VALUE))
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
      zstep=zmax1/REAL(NGLMAX-1)
      DO ngl=1,nglmax
         ZL(NGL)=-zmax1+2*zmax1*REAL(ngl-1)/REAL(NGLMAX-1)
         WLN(NGL)=0.0
         ILN(NGL)=0
         SELECT CASE(NGRAPH)
         CASE(1)
            IF(ZL(NGL)>=0.D0) THEN
               RGB(1,NGL)=1.0
               RGB(2,NGL)=0.0
               RGB(3,NGL)=0.0
            ELSE
               RGB(1,NGL)=0.0
               RGB(2,NGL)=0.0
               RGB(3,NGL)=1.0
            ENDIF
         CASE(2)
               CALL GUSRGB(REAL(NGL-1.5)/REAL(NGLMAX-1), &
                           RGB(1:3,NGL),NRGBA,GLA,GRGBA)
         END SELECT
      END DO
      IF(NGRAPH==2) THEN
         CALL GUSRGB(REAL(NGLMAX+0.5)/REAL(NGLMAX-1), &
                     RGB(1:3,NGLMAX+1),NRGBA,GLA,GRGBA)
      END IF
      ISPL=0

      CALL GDEFIN(2.,17.,2.,17.,-RMAX,RMAX,-RMAX,RMAX)
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,RSTEP,0.,RSTEP,1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,RSTEP*2,0.,RSTEP*2,NGULEN(2*RSTEP))
      CALL SETLNW(0.035)

      SELECT CASE(NGRAPH)
      CASE(1)
         IF(DBLE(ZMIN)*DBLE(ZMAX).GE.0.d0) THEN ! To avoid floating overflow
            zstep=(zmax1-zmin1)/REAL(NGLMAX-1)
            CALL SETRGB(1.0,0.0,0.0)
            CALL CONTQ4(Z,R,PH,NRMAX,NRMAX,NPHMAX,ZMIN1,ZSTEP,NGLMAX, &
                        0,KA)
         ELSE
            zstep=zmax1/REAL((NGLMAX-1)/2)
            CALL SETRGB(1.0,0.0,0.0)
            CALL CONTQ4(Z,R,PH,NRMAX,NRMAX,NPHMAX,0.5*ZSTEP,ZSTEP,NGLMAX/2, &
                        0,KA)
            CALL SETRGB(0.0,0.0,1.0)
            CALL CONTQ4(Z,R,PH,NRMAX,NRMAX,NPHMAX,-0.5*ZSTEP,-ZSTEP,NGLMAX/2, &
                        0,KA)
         ENDIF
      CASE(2)
         CALL CONTF4(Z,R,PH,NRMAX,NRMAX,NPHMAX,ZL,RGB,NGLMAX)
      END SELECT

      DPH=2.D0*PI/NPHMAX
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,1.0,0.0)
      CALL MOVE2D(GUCLIP(RGMIN),0.0)
      DO NPH=1,NPHMAX
         PHL=DPH*NPH
         CALL DRAW2D(GUCLIP(RGMIN*COS(PHL)),GUCLIP(RGMIN*SIN(PHL)))
      END DO
      CALL MOVE2D(GUCLIP(RAXIS),0.0)
      DO NPH=1,NPHMAX
         PHL=DPH*NPH
         CALL DRAW2D(GUCLIP(RAXIS*COS(PHL)),GUCLIP(RAXIS*SIN(PHL)))
      END DO
      CALL MOVE2D(GUCLIP(RGMAX),0.0)
      DO NPH=1,NPHMAX
         PHL=DPH*NPH
         CALL DRAW2D(GUCLIP(RGMAX*COS(PHL)),GUCLIP(RGMAX*SIN(PHL)))
      END DO

      CALL move(2.0,17.1)
      SELECT CASE(MODE)
      CASE(1)
         CALL text('Er real',7)
      CASE(2)
         CALL text('Er imag',7)
      CASE(3)
         CALL text('Er abs ',7)
      CASE(4)
         CALL text('Eth real',8)
      CASE(5)
         CALL text('Eth imag',8)
      CASE(6)
         CALL text('Eth abs ',8)
      CASE(7)
         CALL text('Eph real',8)
      CASE(8)
         CALL text('Eph imag',8)
      CASE(9)
         CALL text('Eph abs ',8)
      CASE(11:16)
         CALL text('Pabs ',5)
         CALL numbi(MODE-6,'(I1)',4)
      END SELECT

      CALL WMGPRM('C','R',0,0,0,0)
      
    CALL pagee

    DEALLOCATE(Z,R,PH)
    DEALLOCATE(ZL,WLN,ILN,RGB)

    RETURN
  END SUBROUTINE wmg3d1

  SUBROUTINE wmg3d2(CEFLD3D,CPABS3D,MDM,NPHM,NRM,NSM, &
                    NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX, &
                    NN,MODE,NGRAPH)
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    COMPLEX(8),DIMENSION(MDM,NPHM,NRM,NSM),INTENT(IN):: CPABS3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM,NSM
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX,NN,MODE,NGRAPH
    REAL(4),DIMENSION(NRM,MDM):: Z
    INTEGER:: NR,NTH,NPH
    CHARACTER(LEN=1):: K2,K3,K4
    INTERFACE 
       FUNCTION GUCLIP(D)
         REAL(8),INTENT(IN):: D
         REAL(4):: GUCLIP
       END FUNCTION GUCLIP
    END INTERFACE

    NPH=NN+NPHMAX_IN/2+1
!         write(6,'(A,1P2E12.4)') 'CEFLD3D(1,5,NPH,5)=',CEFLD3D(1,5,NPH,5)
    SELECT CASE(MODE)
    CASE(1)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CEFLD3D(1,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='R'
       K4='R'
    CASE(2)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(IMAG(CEFLD3D(1,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='R'
       K4='I'
    CASE(3)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(ABS(CEFLD3D(1,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='R'
       K4='A'
    CASE(4)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CEFLD3D(2,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='T'
       K4='R'
    CASE(5)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(IMAG(CEFLD3D(2,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='T'
       K4='I'
    CASE(6)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(ABS(CEFLD3D(2,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='T'
       K4='A'
    CASE(7)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CEFLD3D(3,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='Z'
       K4='R'
    CASE(8)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(IMAG(CEFLD3D(3,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='Z'
       K4='I'
    CASE(9)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(ABS(CEFLD3D(3,NTH,NPH,NR)))
          END DO
       END DO
       K2='E'
       K3='Z'
       K4='A'
    CASE(11:16)
       DO NR=1,NRMAX_IN+1
          DO NTH=1,NTHMAX_IN
             Z(NR,NTH)=GUCLIP(REAL(CPABS3D(NTH,NPH,NR,MODE-10)))
          END DO
       END DO
       K2='P'
       K3=' '
       K4=' '
    END SELECT

    IF(NGRAPH.EQ.0) THEN
       CALL WMGFWR(Z,1,K2,K3,K4)
       GOTO 9000
    ELSEIF(NGRAPH.EQ.1) THEN
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL WMGXEQ(Z,1,K2,K3)
    ELSE IF(NGRAPH.EQ.2) THEN
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL WMGXEQP(Z,1,K2,K3)
    ELSE IF(NGRAPH.EQ.3) THEN
       IF(NTHMAX_IN.LE.2) THEN
          WRITE(6,*) 'XX NTHMAX:',NTHMAX_IN,'  CONDITION:NTHMAX >= 4'
          GO TO 9000
       ENDIF
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL WMG3DA(Z,1,K2,K3,0)
    ELSE IF(NGRAPH.EQ.4) THEN
       IF(NTHMAX_IN.LE.2) THEN
          WRITE(6,*) 'XX NTHMAX:',NTHMAX_IN,'  CONDITION:NTHMAX >= 4'
          GO TO 9000
       ENDIF
       CALL PAGES
       CALL SETCHS(0.3,0.0)
       CALL WMG3DA(Z,1,K2,K3,4)
    ELSE
       WRITE(6,*) 'XX WMGREQ: UNDEFINED NGRAPH: NGRAPH=',NGRAPH
       GO TO 9000
    ENDIF

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
9000 CONTINUE
    RETURN
  END SUBROUTINE wmg3d2

  SUBROUTINE wmg3d4(CEFLD3D,CPABS3D,MDM,NPHM,NRM,NSM, &
                    NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX, &
                    RGMIN,RGMAX,RAXIS,MODE,NCONT)

    USE bpsd_constants,ONLY:PI,CI
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,MDM,NPHM,NRM),INTENT(IN):: CEFLD3D
    COMPLEX(8),DIMENSION(MDM,NPHM,NRM,NSM),INTENT(IN):: CPABS3D
    INTEGER,INTENT(IN):: MDM,NPHM,NRM,NSM
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,NSMAX
    REAL(8),INTENT(IN):: RGMIN,RGMAX,RAXIS
    INTEGER,INTENT(IN):: MODE,NCONT
    REAL(4),DIMENSION(:,:),ALLOCATABLE:: Z
    INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: KA
    
    REAL(4),DIMENSION(:),ALLOCATABLE:: TH
    REAL(4),DIMENSION(:),ALLOCATABLE:: PH
    REAL(4),DIMENSION(:),ALLOCATABLE:: ZL,WLN
    INTEGER,DIMENSION(:),ALLOCATABLE:: ILN
    REAL(4),DIMENSION(:,:),ALLOCATABLE:: RGB
    INTEGER:: nr,nth,nph,nn,nglmax,ngl,ispl,iprd
    INTEGER:: nr_in,nth_in,nph_in
    REAL(4):: rmin,rmax,rr,rb,phmin,phmax,dph,zmin,zmax,thmin,thmax,dth
    REAL(4):: rmin1,rmax1,rstep,zmin1,zmax1,zstep,fsign,phl,zstep1
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

    NGLMAX=NCONT

    ALLOCATE(Z(NPHMAX_IN+1,NTHMAX_IN+1),PH(NPHMAX_IN+1),TH(NTHMAX_IN+1))
    ALLOCATE(KA(2,NPHMAX_IN+1,NTHMAX_IN+1))
    ALLOCATE(ZL(NGLMAX),WLN(NGLMAX),ILN(NGLMAX),RGB(3,NGLMAX))

    PHMIN=-PI
    PHMAX= PI
    DPH=(PHMAX-PHMIN)/NPHMAX_IN
    DO NPH=1,NPHMAX_IN+1
       PH(NPH)=PHMIN+DPH*(NPH-1)
    END DO

    THMIN=-PI
    THMAX= PI
    DTH=(THMAX-THMIN)/NTHMAX_IN
    DO NTH=1,NTHMAX_IN+1
       TH(NTH)=THMIN+DTH*(NTH-1)
    END DO

    NR=NRMAX_IN+1
    DO NPH=1,NPHMAX_IN+1
       DO NTH=1,NTHMAX_IN+1
          IF(NTH<=NTHMAX_IN/2) THEN
             NTH_IN=NTH+NTHMAX_IN/2
          ELSE
             NTH_IN=NTH-NTHMAX_IN/2
          ENDIF
          VALUE=(0.D0,0.D0)
          SELECT CASE(MODE)
          CASE(1,2,3)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(1,NTH_IN,NPH_IN,NR) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(4,5,6)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(2,NTH_IN,NPH_IN,NR) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(7,8,9)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(3,NTH_IN,NPH_IN,NR) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(11:16)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CPABS3D(NTH_IN,NPH_IN,NR,MODE-10)
             END DO
          END SELECT
          SELECT CASE(MODE)
          CASE(1,4,7)
             Z(NPH,NTH)=GUCLIP(REAL(VALUE))
          CASE(2,5,8)
             Z(NPH,NTH)=GUCLIP(IMAG(VALUE))
          CASE(3,6,9)
             Z(NPH,NTH)=GUCLIP(ABS(VALUE))
          CASE(11:16)
             Z(NPH,NTH)=GUCLIP(REAL(VALUE))
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
      CALL GMNMX2(Z,NPHMAX_IN+1,1,NPHMAX_IN+1,1,1,NTHMAX_IN+1,1,ZMIN,ZMAX)
      CALL GQSCAL(ZMIN,ZMAX,ZMIN1,ZMAX1,ZSTEP)
!
      zmax1=MAX(ABS(zmin1),ABS(zmax1))
      zstep1=2*zmax1/REAL(NGLMAX-1)
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
      IPRD=0

      CALL GDEFIN(2.,17.,2.,17.,-GUCLIP(PI),GUCLIP(PI), &
                                -GUCLIP(PI),GUCLIP(PI))
      CALL GFRAME
      CALL SETLNW(0.035)
      CALL GSCALE(0.,0.25*GUCLIP(PI),0.,0.25*GUCLIP(PI),1.0,0)
      CALL SETLNW(0.07)
      CALL GVALUE(0.,0.50*GUCLIP(PI),0.,0.50*GUCLIP(PI),2)

!      CALL CONTG2(Z,PH,TH,NPHMAX_IN+1,NPHMAX_IN+1,NTHMAX_IN+1,ZL,RGB,ILN,WLN, &
!                  NGLMAX,ISPL,IPRD,KA)
      CALL SETLNW(0.035)
      CALL SETRGB(1.0,0.0,0.0)
      CALL CONTR2(Z,PH,TH,NPHMAX_IN+1,NPHMAX_IN+1,NTHMAX_IN+1, &
                  0.5*zstep1,zstep1,NGLMAX/2,IPRD)
      CALL SETRGB(0.0,0.0,1.0)
      CALL CONTR2(Z,PH,TH,NPHMAX_IN+1,NPHMAX_IN+1,NTHMAX_IN+1, &
                  -0.5*zstep1,-zstep1,NGLMAX/2,IPRD)

      CALL move(2.0,17.1)
      SELECT CASE(MODE)
      CASE(1)
         CALL text('Er real',7)
      CASE(2)
         CALL text('Er imag',7)
      CASE(3)
         CALL text('Er abs ',7)
      CASE(4)
         CALL text('Eth real',8)
      CASE(5)
         CALL text('Eth imag',8)
      CASE(6)
         CALL text('Eth abs ',8)
      CASE(7)
         CALL text('Eph real',8)
      CASE(8)
         CALL text('Eph imag',8)
      CASE(9)
         CALL text('Eph abs ',8)
      CASE(11:16)
         CALL text('Pabs ',5)
         CALL numbi(MODE-6,'(I1)',4)
      END SELECT

      CALL WMGPRM('C','R',0,0,0,0)
      
    CALL pagee

    DEALLOCATE(Z,TH,PH,KA)
    DEALLOCATE(ZL,WLN,ILN,RGB)

    RETURN
  END SUBROUTINE wmg3d4

END MODULE wmtest
