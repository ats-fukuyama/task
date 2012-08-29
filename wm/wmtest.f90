MODULE wmtest
  USE wmcont

CONTAINS

  SUBROUTINE wm_test(CEFLD3D,NTHMAX_IN,NPHMAX_IN,NRMAX_IN,RR,RA,MODE)

    USE bpsd_constants,ONLY:PI,CI
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(3,NTHMAX_IN,NPHMAX_IN,NRMAX_IN),INTENT(IN):: CEFLD3D
    REAL(8),INTENT(IN):: RR,RA
    INTEGER,INTENT(IN):: NTHMAX_IN,NPHMAX_IN,NRMAX_IN,MODE
    INTEGER,PARAMETER:: NRM=400
    INTEGER,PARAMETER:: NPHM=128
    REAL(4),DIMENSION(NRM,NPHM):: Z
    REAL(4),DIMENSION(NRM):: R
    REAL(4),DIMENSION(NPHM):: PH
    INTEGER,PARAMETER:: NGLM=30  
    REAL(4),DIMENSION(NGLM):: ZL,WLN
    INTEGER,DIMENSION(NGLM):: ILN
    REAL(4),DIMENSION(3,NGLM):: RGB
    INTEGER:: nrmax,nphmax,nr,nth,nph,nn,nglmax,ngl,ispl
    INTEGER:: nr_in,nth_in,nph_in
    REAL(4):: rmin,rmax,dr,phmin,phmax,dph,zmin,zmax
    REAL(4):: rmin1,rmax1,rstep,zmin1,zmax1,zstep
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
    NRMAX=2*NRMAX_IN
    RMIN=SNGL(RR-RA)
    RMAX=SNGL(RR+RA)
    DR=(RMAX-RMIN)/(NRMAX-1)
    DO NR=1,NRMAX
       R(NR)=RMIN+DR*(NR-1)
    END DO
    PHMIN=0.D0
    PHMAX=2.D0*PI
    DPH=(PHMAX-PHMIN)/(NPHMAX-1)
    DO NPH=1,NPHMAX
       PH(NPH)=PHMIN+DPH*(NPH-1)
    END DO

    DO NPH=1,NPHMAX
       DO NR=1,NRMAX
          IF(NR.LE.NRMAX_IN) THEN
             NR_IN=NRMAX_IN-NR+1
             NTH_IN=NTHMAX_IN/2+1
          ELSE
             NR_IN=NR-NRMAX_IN
             NTH_IN=1
          ENDIF
          VALUE=(0.D0,0.D0)
          SELECT CASE(MODE)
          CASE(0,1)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(1,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(2,3)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(2,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          CASE(4,5)
             DO NPH_IN=1,NPHMAX_IN
                NN=NPH_IN-NPHMAX_IN/2-1
                VALUE=VALUE+CEFLD3D(3,NTH_IN,NPH_IN,NR_IN) &
                           *EXP(CI*NN*PH(NPH))
             END DO
          END SELECT
          SELECT CASE(MODE)
          CASE(0,2,4)
             Z(NR,NPH)=GUCLIP(REAL(VALUE))
          CASE(1,3,5)
             Z(NR,NPH)=GUCLIP(IMAG(VALUE))
          END SELECT
       END DO
    END DO


    CALL pages
      CALL SETLNW(0.07)
      CALL SETCHS(.3,0.)
      CALL SETFNT(32)
!
      CALL GMNMX2(Z,NRM,1,NPHMAX,1,1,NPHMAX,1,ZMIN,ZMAX)
      CALL GQSCAL(ZMIN,ZMAX,ZMIN1,ZMAX1,ZSTEP)
      CALL GQSCAL(-RMAX,RMAX,RMIN1,RMAX1,RSTEP)
!
      nglmax=30
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

      CALL CONTG4X(Z,R,PH,NRM,NRMAX,NPHMAX,ZL,RGB,ILN,WLN,NGLMAX,ISPL)
    CALL pagee

    RETURN
  END SUBROUTINE wm_test
END MODULE wmtest
