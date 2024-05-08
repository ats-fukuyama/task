!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE FOR COMPARISON

!     ***********************************************************

      SUBROUTINE TRCOMP(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: K2
      INTEGER,       INTENT(IN):: INQ

      IF(RHOA.NE.1.D0) NRMAX=NROMAX

      IF(K2.EQ.'1') CALL TRCMP1(INQ)
      IF(K2.EQ.'7') CALL TRCMP7(INQ)

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRCOMP

!     ***********************************

!        COMPARE WITH DIFFERENT MODELS

!     ***********************************

      SUBROUTINE TRCMP1(INQ)

      USE TRCOMM, ONLY : ADNC, AJBS, AKNC, AKNCP, AKNCT, DR, ETA,&
           & ETANC, GAD, GAK, GET, GJB, GRG, GRM, MDLETA, MDLKNC,&
           & MDLNCL, NRMAX, NRMP, NSM, NSMAX, PNSS, PTS, RGFLS, RHOG,&
           & RHOM, RJCB, RKEV, RN, RQFLS, RT, rkind
      USE libitp
      IMPLICIT NONE
      INTEGER,INTENT(IN):: INQ
      INTEGER                  :: IERR, MDLETASTCK, MDLKNCSTCK, MDLNCLSTCK, MODE, NR, NS
      REAL(rkind), DIMENSION(NRMAX)     :: AJBSSTCK, ETASTCK
      REAL(rkind), DIMENSION(NRMAX,NSMAX) :: ADNCL, AKNCL, DNN, DTN, RNN
      REAL(rkind), DIMENSION(NRMP,NSMAX):: RGFLSUM, RQFLSUM
      REAL :: GLOG


      ETASTCK(1:NRMAX)=ETA(1:NRMAX)
      AJBSSTCK(1:NRMAX)=AJBS(1:NRMAX)

!     *** Bootstrap Current and Neoclassical Resistivity ***

      CALL TRAJBS
      GJB(1:NRMAX,1)=SNGL(AJBS(1:NRMAX))*1.E-6
      CALL TRAJBSNEW
      GJB(1:NRMAX,2)=SNGL(AJBS(1:NRMAX))*1.E-6
      CALL TRAJBSSAUTER
      GJB(1:NRMAX,3)=SNGL(AJBS(1:NRMAX))*1.E-6

      MDLETASTCK=MDLETA
      MDLNCLSTCK=MDLNCL
      IF(MDLNCL.EQ.1) MDLNCL=0
      DO MDLETA=1,4
         CALL TRCFET
         DO NR=1,NRMAX
            GET(NR,MDLETA)=GLOG(ETA(NR),1.D-10,1.D0)
         ENDDO
      ENDDO
      MDLETA=MDLETASTCK
      MDLNCL=MDLNCLSTCK

      IF(MDLNCL.EQ.0) MDLNCL=1
      CALL TR_NCLASS(IERR)
      CALL TRAJBS_NCLASS
      GJB(1:NRMAX,4)=SNGL(AJBS(1:NRMAX))*1.E-6
      DO NR=1,NRMAX
         GET(NR,5)=GLOG(ETANC(NR),1.D-10,1.D0)
      ENDDO
      MDLNCL=MDLNCLSTCK
      ETA(1:NRMAX)=ETASTCK(1:NRMAX)
      AJBS(1:NRMAX)=AJBSSTCK(1:NRMAX)

!     *** Neoclassical Particle and Heat Flux Diffusivity ***

      MODE=0
      IF(MODE.EQ.0) THEN
         DO NS=1,2
            AKNCL(1:NRMAX,NS)=AKNCP(1:NRMAX,NS,NS)+AKNCT(1:NRMAX,NS,NS)
         ENDDO
         ADNCL(1:NRMAX,1:2)=ADNC(1:NRMAX,1:2)
      ELSE
         DO NR=1,NRMAX
            DO NS=1,NSMAX
               RGFLSUM(NR,NS)=SUM(RGFLS(NR,1:5,NS))
               RQFLSUM(NR,NS)=SUM(RQFLS(NR,1:5,NS))
            ENDDO
         ENDDO
         DO NS=1,NSMAX
            RNN(1:NRMAX-1,NS)=(RN(1:NRMAX-1+1,NS)+RN(1:NRMAX-1,NS))*0.5D0
            DNN(1:NRMAX-1,NS)=(RN(1:NRMAX-1+1,NS)-RN(1:NRMAX-1,NS))     *RJCB(1:NRMAX-1)/DR
            DTN(1:NRMAX-1,NS)=(RT(1:NRMAX-1+1,NS)-RT(1:NRMAX-1,NS))*RKEV*RJCB(1:NRMAX-1)/DR
         ENDDO
         NR=NRMAX
         DO NS=1,NSMAX
            RNN(NR,NS)=PNSS(NS)
            DNN(NR,NS)=DERIV3P(PNSS(NS),RN(NR,NS),RN(NR-1,NS),RHOG(NR),RHOM(NR),RHOM(NR-1))
            DTN(NR,NS)=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS), RHOG(NR),RHOM(NR),RHOM(NR-1))*RKEV
         ENDDO
         ADNCL(1:NRMAX,1:NSMAX)=-RGFLSUM(1:NRMAX,1:NSMAX)/DNN(1:NRMAX,1:NSMAX)
         AKNCL(1:NRMAX,1:NSMAX)=-RQFLSUM(1:NRMAX,1:NSMAX)/(RNN(1:NRMAX,1:NSMAX)*DTN(1:NRMAX,1:NSMAX))
      ENDIF

      MDLKNCSTCK=MDLKNC
      DO MDLKNC=1,3,2
         CALL TRCFNC
         GAK(2:NRMAX+1,MDLKNC+2) = SNGL(AKNC(1:NRMAX,1))
         GAK(2:NRMAX+1,MDLKNC+3) = SNGL(AKNC(1:NRMAX,2))
         GAK(1,MDLKNC+2) = SNGL(AKNC(1,1))
         GAK(1,MDLKNC+3) = SNGL(AKNC(1,2))
      ENDDO
      MDLKNC=MDLKNCSTCK
      GAD(2:NRMAX+1,1) = SNGL(ADNCL(1:NRMAX,1))
      GAD(2:NRMAX+1,2) = SNGL(ADNCL(1:NRMAX,2))
      GAK(2:NRMAX+1,1) = SNGL(AKNCL(1:NRMAX,1))
      GAK(2:NRMAX+1,2) = SNGL(AKNCL(1:NRMAX,2))
      GAD(2:NRMAX+1,3) = SNGL(ADNC(1:NRMAX,1))
      GAD(2:NRMAX+1,4) = SNGL(ADNC(1:NRMAX,2))
      GAD(1,1) = SNGL(ADNCL(1,1))
      GAD(1,2) = SNGL(ADNCL(1,2))
      GAK(1,1) = SNGL(AKNCL(1,1))
      GAK(1,2) = SNGL(AKNCL(1,2))
      GAD(1,3) = SNGL(ADNC(1,1))
      GAD(1,4) = SNGL(ADNC(1,2))

!     *** Graphic Routine ***

      CALL PAGES
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GJB,NRMP,NRMAX,4,'@JBS [MA/m$+2$=]  vs r@',2+INQ)
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GET,NRMP,NRMAX,5,'@LOG:ETA  vs r @',11+INQ)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GAD,NRMP,NRMAX+1,4,'@ADNCE, ADNCD [m$+2$=/s]  vs r@',2+INQ)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GAK,NRMP,NRMAX+1,6,'@AKNCE, AKNCD [m$+2$=/s]  vs r @',2+INQ)
      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRCMP1


!     ****************************************************

!        COMPARE WITH RADIAL ELECTRIC FIELD MODELS

!     ****************************************************

      SUBROUTINE TRCMP7(INQ)

      USE TRCOMM
      IMPLICIT NONE
      INTEGER, INTENT(IN):: INQ
      INTEGER:: MDLER_ORG
      REAL(rkind),DIMENSION(NRMAX):: ER_ORG

      CALL PAGES

      MDLER_ORG=MDLER
      ER_ORG(1:NRMAX)=ER(1:NRMAX)
      DO MDLER=1,4
         CALL TRERAD
         GYR(1:NRMAX,MDLER)=SNGL(ER(1:NRMAX))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX,4,'@ER  vs r@',2+INQ)

      MDLER=MDLER_ORG
      ER(1:NRMAX)=ER_ORG(1:NRMAX)

      CALL PAGEE
      CALL TRGRTM

      RETURN
      END SUBROUTINE TRCMP7
