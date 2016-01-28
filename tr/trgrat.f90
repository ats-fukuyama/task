!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRT0(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA, NT, NGTSTP
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: K2
      INTEGER(4),       INTENT(IN):: INQ

      IF(NT.LT.NGTSTP) RETURN
      IF(RHOA.NE.1.D0) NRMAX=NROMAX

      IF(K2.EQ.'1') CALL TRGRT1(INQ)
      IF(K2.EQ.'2') CALL TRGRT2(INQ)
      IF(K2.EQ.'5') CALL TRGRT5(INQ)
      IF(K2.EQ.'6') CALL TRGRT6(INQ)
      IF(K2.EQ.'7') CALL TRGRT7(INQ)
      IF(K2.EQ.'8') CALL TRGRT8(INQ)

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRGRT0

!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRX0(K2,INQ)

      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: K2
      INTEGER(4),       INTENT(IN):: INQ

      IF(K2.EQ.'1') CALL TRGRX1(INQ)
      RETURN
      END SUBROUTINE TRGRX0

!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRA0(K2,INQ)

      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: K2
      INTEGER(4),       INTENT(IN):: INQ

      IF(K2.EQ.'1') THEN
         CALL TRGRT6(INQ)
         CALL TRGRT7(INQ)
         CALL TRGRX1(INQ)
         CALL TRGRR1(INQ)
         CALL TRGRR2(INQ)
         CALL TRGRR4(INQ)
         CALL TRGRR7(INQ)
      ELSEIF(K2.EQ.'2') THEN
         CALL TRGRT6(INQ)
         CALL TRGRT7(INQ)
         CALL TRGRX1(INQ)
         CALL TRGRG1(INQ)
         CALL TRGRG3(INQ)
         CALL TRGRR4(INQ)
         CALL TRGRG5(INQ)
      ENDIF
      RETURN
      END SUBROUTINE TRGRA0

!     ***********************************************************

!           GRAPHIC : TEMPORAL CONTOUR
!                     TE, TD, NE, AJ, TE(XP), TD(XP), AJBS, AJCD

!     ***********************************************************

      SUBROUTINE TRGRT1(INQ)

      USE TRCOMM, ONLY : GT, GVRT, GRM, RTU, NGT, NTM, NRMAX, MDLUF
      USE TRCOM1, ONLY : TMU, TMU1, NTXMAX, NTXMAX1, NTUM
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4) :: I, NR, IERR
      REAL(4),DIMENSION(NTM,NRMAX) :: GYL
      REAL(8) :: TSL, RTEL, RTDL
      REAL(4) :: GUCLIP
      CALL PAGES

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I, 1)
      ENDDO
      CALL TRGR1DC( 3.0,12.0,14.0,17.0,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@TE [keV]  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I, 2)
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@TD [keV]  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I, 5)
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 5.4, 8.4,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@NE [10$+20$=/m$+3$=]  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I, 9) * 1.E-6
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 1.1, 4.1,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@AJ [MA]  vs t@')

      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         DO I=1,NGT
            TSL = DBLE(GT(I))
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,RTEL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
               GYL(I,NR) = GUCLIP(RTEL)
            ENDDO
         ENDDO
         CALL TRGR1DC(15.0,24.0,14.0,17.0,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@TE [keV] (XP) vs t@')

         DO I=1,NGT
            TSL = DBLE(GT(I))
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,RTDL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
               GYL(I,NR) = GUCLIP(RTDL)
            ENDDO
         ENDDO
         CALL TRGR1DC(15.0,24.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@TD [keV] (XP) vs t@')
      ENDIF

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I,13) * 1.E-6
      ENDDO
      CALL TRGR1DC(15.0,24.0, 5.4, 8.4,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@AJBS [MA]  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = (GVRT(1:NRMAX,I,11) + GVRT(1:NRMAX,I,12)) * 1.E-6
      ENDDO
      CALL TRGR1DC(15.0,24.0, 1.1, 4.1,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@AJCD [MA]  vs t@')

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRT1

!     ***********************************************************

!           GRAPHIC : TEMPORAL CONTOUR
!                     QP, VTOR, VLOOP, ER

!     ***********************************************************

      SUBROUTINE TRGRT2(INQ)

      USE TRCOMM, ONLY : GT, GRG, GRM, GVRT, NGT, NTM, NRMAX
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4) :: I
      REAL(4),DIMENSION(NTM,NRMAX) :: GYL

      CALL PAGES

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I,27)
      ENDDO
      CALL TRGR1DC( 3.0,12.0,14.0,17.0,GT,GRG,GYL,NTM,NGT,NRMAX,NRMAX,'@QP  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I,31)
      ENDDO
      CALL TRGR1DC(15.0,24.0,14.0,17.0,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@VLOOP [V]  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I,60)
      ENDDO
      CALL TRGR1DC( 3.0,12.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@VTOR [m/s]  vs t@')

      DO I=1,NGT
         GYL(I,1:NRMAX) = GVRT(1:NRMAX,I,63)
      ENDDO
      CALL TRGR1DC(15.0,24.0, 9.7,12.7,GT,GRM,GYL,NTM,NGT,NRMAX,NRMAX,'@ER [V/m] vs t@')

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRT2

!     ***********************************************************

!           GRAPHIC : TIME DEPENDENCE : RR,RA,RKAP,IP,BB,H98Y2
!                     TIME DEPENDENCE : TRANSIENT

!     ***********************************************************

      SUBROUTINE TRGRT5(INQ)

      USE TRCOMM, ONLY : DR, GT, GTS, GVT, GYT, IZERO, NGPST, NGST, NGT, NRMAX, NTM
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4) :: I, IX, K, M, N, NGPSTH
      REAL(4)    :: GD, GW

      IF(NGST.EQ.0) THEN

      CALL PAGES

      DO I=1,NGT
         GYT(I,1)=GVT(I,97)
         GYT(I,2)=GVT(I,98)
      ENDDO
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,2,'@RR, RA [M]  vs t@',2+INQ)

      DO I=1,NGT
         GYT(I,1)=GVT(I,99)
         GYT(I,2)=GVT(I,101)!GVT(I,34)
      ENDDO
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,2,'@BT, IP  vs t@',2+INQ)

      DO I=1,NGT
         GYT(I,1)=GVT(I,100)
      ENDDO
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,'@RKAP  vs t@',2+INQ)

      DO I=1,NGT
         GYT(I,1)=GVT(I,81)
         GYT(I,2)=GVT(I,102)
         GYT(I,3)=GVT(I,103)
      ENDDO
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,3,'@TAUE89,TAUE98 [s], H98Y2  vs t@',2+INQ)

      CALL PAGEE

      ELSE

    1 WRITE(6,*) ' CHOSE ONE (TE=1,TD=2,TT=3,TA=4) '
      READ(5,*,ERR=1,END=900) M
!
      CALL PAGES

      IX = INT((NRMAX+1-IZERO)/(NGPST-1))
      GW = 150.0/(15.0*NGPST/2.0-5.0)
      GD = 0.5*GW
!      IR = INT((TSST/DT)/NGTSTP)

      DO N=0,NGPST/2-1
         K=4*N

         GYT(1:NGST,1) = GVT(1:NGST,K+89)
         GYT(1:NGST,2) = GVT(1:NGST,K+90)
         GYT(1:NGST,3) = GVT(1:NGST,K+91)
         GYT(1:NGST,4) = GVT(1:NGST,K+92)

         CALL TRGR1D(3.0,12.0,17-(N+1)*GW-N*GD,17-N*(GW+GD), &
     &            GTS,GYT,NTM,NGST,M,'@TE,TD,TT,TA [keV]@',1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(10.0,17.1-N*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)

      ENDDO

      DO N=NGPST/2,NGPST-1
         K=4*N

         GYT(1:NGST,1) = GVT(1:NGST,K+89)
         GYT(1:NGST,2) = GVT(1:NGST,K+90)
         GYT(1:NGST,3) = GVT(1:NGST,K+91)
         GYT(1:NGST,4) = GVT(1:NGST,K+92)

         NGPSTH=N-NGPST/2
         CALL TRGR1D(15.5,24.5,17-(NGPSTH+1)*GW-NGPSTH*GD,17-NGPSTH*(GW+GD), &
     &               GTS,GYT,NTM,NGST,M,'@TE,TD,TT,TA [keV]  vs t@',1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(22.5,17.1-NGPSTH*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)

      ENDDO

      CALL PAGEE

      ENDIF

  900 RETURN
      END SUBROUTINE TRGRT5

!     ***********************************************************

!           GRAPHIC : TIME DEPENDENCE : TEMPERATURE,CURRENT,
!                                       INPUT POWER,OUTPUT POWER
!cpub This subroutine will be called if choosing (T6)
!     ***********************************************************

      SUBROUTINE TRGRT6(INQ)

      USE TRCOMM, ONLY : GT, GVT, GYT, MDLNF, NGT, NTM
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4) :: I

      CALL PAGES

      GYT(1:NGT,1)=GVT(1:NGT, 9)
      GYT(1:NGT,2)=GVT(1:NGT,10)
      GYT(1:NGT,3)=GVT(1:NGT,13)
      GYT(1:NGT,4)=GVT(1:NGT,14)
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,'@TE,TD,<TE>,<TD> [keV]  vs t@',2+INQ)



!      GYT(1:NGT,1)=GVT(1:NGT,34)
      GYT(1:NGT,1)=GVT(1:NGT,101)	! AJTTOR
      GYT(1:NGT,2)=GVT(1:NGT,38)		! AJBST
!      GYT(1:NGT,2)=GVT(1:NGT,35)	! AJOHT
!      GYT(1:NGT,3)=GVT(1:NGT,36)	! AJNBT
!      GYT(1:NGT,4)=GVT(1:NGT,37)	! AJRFT
!      GYT(1:NGT,5)=GVT(1:NGT,38)	! AJBST
!      GYT(1:NGT,6)=GVT(1:NGT,34)	! AJT
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,2,'@IP,IBS [MA]  vs t@',2+INQ)
!cpub      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,1,'@IP,IOH,INB,IRF [MA]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,39)
      GYT(1:NGT,2)=GVT(1:NGT,40)
      GYT(1:NGT,3)=GVT(1:NGT,41)
      GYT(1:NGT,4)=GVT(1:NGT,42)+GVT(1:NGT,43)+GVT(1:NGT,44)+GVT(1:NGT,45)
      GYT(1:NGT,5)=GVT(1:NGT,46)
      GYT(1:NGT,6)=GVT(1:NGT,89)+GVT(1:NGT,90)+GVT(1:NGT,91)+GVT(1:NGT,92)
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,6,'@PIN,POH,PNB,PRF,PNF,PEX [MW]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,57)
      GYT(1:NGT,2)=GVT(1:NGT,58)
      GYT(1:NGT,3)=GVT(1:NGT,59)
      GYT(1:NGT,4)=GVT(1:NGT,60)
      GYT(1:NGT,5)=GVT(1:NGT,61)+GVT(1:NGT,62)+GVT(1:NGT,63)+GVT(1:NGT,64)
!     &           +GVT(1:NGT,89)+GVT(1:NGT,90)
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,5,'@POUT,PCX,PIE,PRL,PCON [MW]  vs t@',2+INQ)

! figure(1,2)
      GYT(1:NGT,1)=GVT(1:NGT,77)
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,'@q0  vs t@',2+INQ)

! figure(1,2)
!      GYT(1:NGT,1)=GVT(1:NGT,87)
!      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,'@QF  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,13)
      GYT(1:NGT,2)=GVT(1:NGT,14)
      GYT(1:NGT,3)=GVT(1:NGT,15)
      GYT(1:NGT,4)=GVT(1:NGT,16)
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,4,'@TEAV,TDAV,TTAV,TAAV [keV]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT, 9)
      GYT(1:NGT,2)=GVT(1:NGT,10)
      GYT(1:NGT,3)=GVT(1:NGT,11)
      GYT(1:NGT,4)=GVT(1:NGT,12)
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,4,'@TE0,TD0,TT0,TA0 [keV]  vs t@',2+INQ)

      IF(MDLNF.EQ.0) THEN
         GYT(1:NGT,1)=GVT(1:NGT,1)
         GYT(1:NGT,2)=GVT(1:NGT,5)
         GYT(1:NGT,3)=GVT(1:NGT,104)
         CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,3,'@NE0,<NE>,<NEL> [10$+20$=/m$+3$=]  vs t@',2+INQ)
      ELSE
         GYT(1:NGT,1)=GVT(1:NGT,1)
         GYT(1:NGT,2)=GVT(1:NGT,2)
         GYT(1:NGT,3)=GVT(1:NGT,3)
         GYT(1:NGT,4)=GVT(1:NGT,4)
         CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,3,'@NE0,ND0,NT0,NA0 [10$+20$=/m$+3$=]  vs t@',2+INQ)
      ENDIF

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRT6

!     ***********************************************************

!           GRAPHIC : TIME DEPENDENCE : TAUE,QF,BETAP,BETA

!     ***********************************************************

      SUBROUTINE TRGRT7(INQ)

      USE TRCOMM, ONLY : GT, GVT, GYT, NGT, NTM, RA, BB
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4) :: I

      CALL PAGES

      GYT(1:NGT,1)=GVT(1:NGT,33)
      GYT(1:NGT,2)=GVT(1:NGT,31)+GVT(1:NGT,29)
      GYT(1:NGT,3)=GVT(1:NGT,31)
      GYT(1:NGT,4)=GVT(1:NGT,17)
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,'@WF,WB,WI,WE [MJ]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,85)*100.0
      GYT(1:NGT,2)=GVT(1:NGT,84)*100.0
      GYT(1:NGT,3)=GVT(1:NGT,85)*100.0/(GVT(1:NGT,101)/SNGL(RA*BB))!(GVT(1:NGT,34)/SNGL(RA*BB))
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,3,'@BETAa,BETA0,[%],BETAN  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,79)
      GYT(1:NGT,2)=GVT(1:NGT,80)
      GYT(1:NGT,3)=GVT(1:NGT,81)
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,3,'@TAUE1,TAUE2,TAUE89  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,83)
      GYT(1:NGT,2)=GVT(1:NGT,82)
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,2,'@BETAPa,BETAP0  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,74)
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,1,'@VLOOP [V]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,86)
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,1,'@ZEFF0  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,103)
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,1,'@H98Y2  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,77)
      CALL TRGR1D(15.0,24.0, 1.1, 4.1,GT,GYT,NTM,NGT,1,'@Q0  vs t@',2+INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRT7

!     ***********************************************************

!           GRAPHIC : TIME DEPENDENCE : TRANSIENT

!     ***********************************************************

      SUBROUTINE TRGRT8(INQ)

      USE TRCOMM, ONLY : DR, GRM, GT, GTS, GYT, GVT, MDLUF, NGT, NGST, NGPST, NRMAX, NTM, IZERO
!!!      USE TRCOM1, ONLY : NTXMAX, TMU
      IMPLICIT NONE
      INTEGER(4),INTENT(IN):: INQ
      INTEGER(4),DIMENSION(8,NTM,NRMAX):: KATR
      INTEGER(4) :: I, IERR, NR, M, IX, N, K, NGPSTH
      REAL(8)    :: RTDL, RTEL, TSL
      REAL(4),DIMENSION(NTM,NRMAX)  :: GYL
      REAL(4)    :: GW, GD

      IF(NGST.EQ.0) RETURN

    1 WRITE(6,*) ' CHOSE ONE (TE=1,TD=2,TT=3,TA=4) '
      READ(5,*,ERR=1,END=900) M

      CALL PAGES

      IX = INT((NRMAX+1-IZERO)/(NGPST-1))
      GW = 150.0/(15.0*NGPST/2.0-5.0)
      GD = 0.5*GW
!      IR = INT((TSST/DT)/NGTSTP)

      DO N=0,NGPST/2-1
         K=4*N

         GYT(1:NGST,1) = GVT(1:NGST,K+89)
         GYT(1:NGST,2) = GVT(1:NGST,K+90)
         GYT(1:NGST,3) = GVT(1:NGST,K+91)
         GYT(1:NGST,4) = GVT(1:NGST,K+92)

         CALL TRGR1D(3.0,12.0,17-(N+1)*GW-N*GD,17-N*(GW+GD), &
                     GTS,GYT,NTM,NGST,M,'@TE,TD,TT,TA [keV]@',1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(10.0,17.1-N*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)
        
      ENDDO

      DO N=NGPST/2,NGPST-1
         K=4*N

         GYT(1:NGST,1) = GVT(1:NGST,K+89)
         GYT(1:NGST,2) = GVT(1:NGST,K+90)
         GYT(1:NGST,3) = GVT(1:NGST,K+91)
         GYT(1:NGST,4) = GVT(1:NGST,K+92)

         NGPSTH=N-NGPST/2
         CALL TRGR1D(15.5,24.5,17-(NGPSTH+1)*GW-NGPSTH*GD,17-NGPSTH*(GW+GD), &
                     GTS,GYT,NTM,NGST,M,'@TE,TD,TT,TA [keV]  vs t@',1+INQ)
         CALL SETCHS(0.3,0.0)
         CALL SETFNT(32)
         CALL SETLIN(-1,-1,7)
         CALL MOVE(22.5,17.1-NGPSTH*(GW+GD))
         CALL NUMBR(REAL((IZERO-1+N*IX)*DR),'(F7.3)',7)
         CALL TEXT('[m]',3)
        
      ENDDO

      CALL PAGEE

900   RETURN
      END SUBROUTINE TRGRT8

!     ***********************************************************

!           GRAPHIC : TIME LISSAGE : tau vs eps*betap,li,

!     ***********************************************************

      SUBROUTINE TRGRX1(INQ)

      USE TRCOMM, ONLY : GVT, GYT, NGT, NTM, RA, RR
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4) :: I

      CALL PAGES

      GYT(1:NGT,1)=SNGL(RA/RR)*GVT(1:NGT,83)
      GYT(1:NGT,2)=GVT(1:NGT,80)/GVT(1:NGT,81)
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GYT(1,1),GYT(1,2),NTM,NGT,1,'@tauE/tauE89 vs eps*betap@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,75)
      GYT(1:NGT,2)=GVT(1:NGT,80)/GVT(1:NGT,81)
      CALL TRGR1D(15.5,24.5,11.0,17.0,GYT(1,1),GYT(1,2),NTM,NGT,1,'@tauE/tauE89 vs li@',2+INQ)

      GYT(1:NGT,1)=SNGL(RA/RR)*GVT(1:NGT,83)
      GYT(1:NGT,2)=GVT(1:NGT,38)/GVT(1:NGT,101)!GVT(1:NGT,34)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GYT(1,1),GYT(1,2),NTM,NGT,1,'@Ibs/Ip vs eps*betap@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,77)
      GYT(1:NGT,2)=GVT(1:NGT,80)/GVT(1:NGT,81)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GYT(1,1),GYT(1,2),NTM,NGT,1,'@tauE/tauE89 vs q0@',2+INQ)

      CALL PAGEE

      RETURN
      END SUBROUTINE TRGRX1
