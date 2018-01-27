!
!     ***********************************************************

!           CALCULATE Z-EFFF

!     ***********************************************************

      SUBROUTINE TRZEFF

      USE TRCOMM, ONLY : ANC,ANFE,MDLEQN,MDLUF,NRMAX,PZ,PZC,PZFE,RN,RNF,RT,ZEFF
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8)   :: TE, TRZEC, TRZEFE

      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            TE =RT(NR,1)
            ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +PZ(4)  *PZ(4)  *RN(NR,4) &
     &                +TRZEC(TE)**2   *ANC (NR) +TRZEFE(TE)**2  *ANFE(NR))/RN(NR,1)
         ENDDO
      ELSE
         IF(MDLEQN.EQ.0) THEN ! fixed density
            DO NR=1,NRMAX
               ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +PZ(2)  *PZ(2)  *RNF(NR,1))/RN(NR,1)
            ENDDO
         ENDIF
      ENDIF

      DO NR=1,NRMAX
         TE=RT(NR,1)
         PZC(NR)=TRZEC(TE)
         PZFE(NR)=TRZEFE(TE)
      ENDDO
!
      RETURN
      END SUBROUTINE TRZEFF

!     ***********************************************************

!           CALCULATE Z-C

!     **********************************************************

      REAL(8) FUNCTION TRZEC(TE)

      IMPLICIT NONE
      REAL(8):: TE,TEL
      REAL(8):: BC00= 2.093736D+03, BC01= 5.153766D+03, BC02= 5.042105D+03, &
                BC03= 2.445345D+03, BC04= 5.879207D+02, BC05= 5.609128D+01
      REAL(8):: BC10= 7.286051D+00, BC11= 4.506650D+01, BC12= 1.595483D+02, &
                BC13= 2.112702D+02, BC14= 1.186473D+02, BC15= 2.400040D+01
      REAL(8):: BC20= 5.998782D+00, BC21= 6.808206D-03, BC22=-3.096242D-02, &
                BC23=-4.794194D-02, BC24= 1.374732D-01, BC25= 5.157571D-01
      REAL(8):: BC30= 5.988153D+00, BC31= 7.335329D-02, BC32=-1.754858D-01, &
                BC33= 2.034126D-01, BC34=-1.145930D-01, BC35= 2.517150D-02
      REAL(8):: BC40=-3.412166D+01, BC41= 1.250454D+02, BC42=-1.550822D+02, &
                BC43= 9.568297D+01, BC44=-2.937297D+01, BC45= 3.589667D+00

      TEL = LOG10(TE)
      IF(TE.LE.3.D-3) THEN
         TRZEC=0.D0
      ELSEIF(TE.LE.2.D-2) THEN
         TRZEC= BC00+(BC01*TEL)+(BC02*TEL**2)+(BC03*TEL**3) &
                    +(BC04*TEL**4)+(BC05*TEL**5)
      ELSEIF(TE.LE.0.2D0) THEN
         TRZEC= BC10+(BC11*TEL)+(BC12*TEL**2)+(BC13*TEL**3) &
                    +(BC14*TEL**4)+(BC15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         TRZEC= BC20+(BC21*TEL)+(BC22*TEL**2)+(BC23*TEL**3) &
                    +(BC24*TEL**4)+(BC25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         TRZEC= BC30+(BC31*TEL)+(BC32*TEL**2)+(BC33*TEL**3) &
                    +(BC34*TEL**4)+(BC35*TEL**5)
      ELSEIF(TE.LE.100.D0) THEN
         TRZEC= BC40+(BC41*TEL)+(BC42*TEL**2)+(BC43*TEL**3) &
                    +(BC44*TEL**4)+(BC45*TEL**5)
      ELSE
         TRZEC= 6.D0
      ENDIF

      RETURN
      END FUNCTION TRZEC

!     ***********************************************************

!           CALCULATE Z-FE

!     ***********************************************************

      REAL(8) FUNCTION TRZEFE(TE)

      IMPLICIT NONE
      REAL(8):: TE, TEL
      REAL(8):: BF10= 8.778318D+00, BF11=-5.581412D+01, BF12=-1.225124D+02, &
                BF13=-1.013985D+02, BF14=-3.914244D+01, BF15=-5.851445D+00
      REAL(8):: BF20= 1.959496D+01, BF21= 1.997852D+01, BF22= 9.593373D+00, &
                BF23=-7.609962D+01, BF24=-1.007190D+02, BF25=-1.144046D+01
      REAL(8):: BF30= 1.478150D+01, BF31= 6.945311D+01, BF32=-1.991202D+02, &
                BF33= 2.653804D+02, BF34=-1.631423D+02, BF35= 3.770026D+01
      REAL(8):: BF40= 2.122081D+01, BF41= 6.891607D+00, BF42=-4.076853D+00, &
                BF43= 1.577171D+00, BF44=-5.139468D-01, BF45= 8.934176D-02

      TEL = LOG10(TE)

      IF(TE.LE.3.D-3) THEN
         TRZEFE=0.D0
      ELSEIF(TE.LE.0.2D0) THEN
         TRZEFE = BF10+(BF11*TEL)+(BF12*TEL**2)+(BF13*TEL**3) &
                      +(BF14*TEL**4)+(BF15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         TRZEFE = BF20+(BF21*TEL)+(BF22*TEL**2)+(BF23*TEL**3) &
                      +(BF24*TEL**4)+(BF25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         TRZEFE = BF30+(BF31*TEL)+(BF32*TEL**2)+(BF33*TEL**3) &
                      +(BF34*TEL**4)+(BF35*TEL**5)
      ELSEIF(TE.LE.100.D0) THEN
         TRZEFE = BF40+(BF41*TEL)+(BF42*TEL**2)+(BF43*TEL**3) &
                      +(BF44*TEL**4)+(BF45*TEL**5)
      ELSE
         TRZEFE = 26.D0
      ENDIF

      RETURN
      END FUNCTION TRZEFE

!     ***********************************************************

!           RADIATION POWER - C

!     ***********************************************************

      REAL(8) FUNCTION TRRPC(TE)

      IMPLICIT NONE
      REAL(8):: TE,TEL,ARG
      REAL(8):: AC00=1.965300D+03, AC01=4.572039D+03, AC02=4.159590D+03, AC03=1.871560D+03, AC04=4.173889D+02, AC05=3.699382D+01
      REAL(8):: AC10=7.467599D+01, AC11=4.549038D+02, AC12=8.372937D+02, AC13=7.402515D+02, AC14=3.147607D+02, AC15=5.164578D+01
      REAL(8):: AC20=-2.120151D+01,AC21=-3.668933D-01,AC22=7.295099D-01, AC23=-1.944827D-01,AC24=-1.263576D-01,AC25=-1.491027D-01
      REAL(8):: AC30=-2.121979D+01,AC31=-2.346986D-01,AC32=4.093794D-01, AC33=7.874548D-02, AC34=-1.841379D-01,AC35=5.590744D-02
      REAL(8):: AC40=-2.476796D+01,AC41=9.408181D+00, AC42=-9.657446D+00,AC43=4.999161D+00, AC44=-1.237382D+00,AC45=1.160610D-01

      TEL = LOG10(TE)
      IF(TE.LE.3.D-3) THEN
         TEL=LOG10(3.D-3)
         ARG=AC00+(AC01*TEL)+(AC02*TEL**2)+(AC03*TEL**3)+(AC04*TEL**4)+(AC05*TEL**5)
      ELSEIF(TE.LE.2.D-2) THEN
         ARG=AC00+(AC01*TEL)+(AC02*TEL**2)+(AC03*TEL**3)+(AC04*TEL**4)+(AC05*TEL**5)
      ELSEIF(TE.LE.0.2D0) THEN
         ARG=AC10+(AC11*TEL)+(AC12*TEL**2)+(AC13*TEL**3)+(AC14*TEL**4)+(AC15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         ARG=AC20+(AC21*TEL)+(AC22*TEL**2)+(AC23*TEL**3)+(AC24*TEL**4)+(AC25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         ARG=AC30+(AC31*TEL)+(AC32*TEL**2)+(AC33*TEL**3)+(AC34*TEL**4)+(AC35*TEL**5)
      ELSEIF(TE.LE.1.D2) THEN
         ARG=AC40+(AC41*TEL)+(AC42*TEL**2)+(AC43*TEL**3)+(AC44*TEL**4)+(AC45*TEL**5)
      ENDIF
      TRRPC = 10.D0**(ARG-13.D0)

      RETURN
      END FUNCTION TRRPC

!     ***********************************************************

!           RADIATION POWER - FE

!     ***********************************************************

      REAL(8) FUNCTION TRRPFE(TE)

      IMPLICIT NONE
      REAL(8)::TE,TEL,ARG
      REAL(8)::AF10=-2.752599D+01,AF11=-3.908228D+01,AF12=-6.469423D+01,AF13=-5.555048D+01,AF14=-2.405568D+01,AF15=-4.093160D+00
      REAL(8)::AF20=-1.834973D+01,AF21=-1.252028D+00,AF22=-7.533115D+00,AF23=-3.289693D+00,AF24=2.866739D+01, AF25=2.830249D+01
      REAL(8)::AF30=-1.671042D+01,AF31=-1.646143D+01,AF32=3.766238D+01, AF33=-3.944080D+01,AF34=1.918529D+01, AF35=-3.509238D+00
      REAL(8)::AF40=-2.453957D+01,AF41=1.795222D+01, AF42=-2.356360D+01,AF43=1.484503D+01, AF44=-4.542323D+00,AF45=5.477462D-01


      TEL = LOG10(TE)
      IF(TE.LE.3.D-3) THEN
         TEL=LOG10(3.D-3)
         ARG=AF10+(AF11*TEL)+(AF12*TEL**2)+(AF13*TEL**3)+(AF14*TEL**4)+(AF15*TEL**5)
      ELSEIF(TE.LE.0.2D0) THEN
         ARG=AF10+(AF11*TEL)+(AF12*TEL**2)+(AF13*TEL**3)+(AF14*TEL**4)+(AF15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         ARG=AF20+(AF21*TEL)+(AF22*TEL**2)+(AF23*TEL**3)+(AF24*TEL**4)+(AF25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         ARG=AF30+(AF31*TEL)+(AF32*TEL**2)+(AF33*TEL**3)+(AF34*TEL**4)+(AF35*TEL**5)
      ELSEIF(TE.LE.1.D2) THEN
         ARG=AF40+(AF41*TEL)+(AF42*TEL**2)+(AF43*TEL**3)+(AF44*TEL**4)+(AF45*TEL**5)
      ENDIF
      TRRPFE = 10.D0**(ARG-13.D0)

      RETURN
      END FUNCTION TRRPFE

!     ***********************************************************

!           POWER LOSS

!     ***********************************************************

      SUBROUTINE TRLOSS

      USE TRCOMM, ONLY : AEE, ANC, ANFE, ANNU, DT, KUFDEV, MDLUF, MDLPR, &
           NRMAX, NT, NTUM, PCX, PIE, PRB, PRC, PRSUM, PRL, PRLU, RKEV, RN, &
           RT, SCX, SIE, TSCX, TSIE
      USE TRCOM1, ONLY : NTAMAX, NTXMAX, PNBI, TMU
      USE tr_cytran_mod, ONLY: tr_cytran
      IMPLICIT NONE
      INTEGER(4):: IERR, NR
      REAL(8):: ANDX, ANE, ANHE, ANT, EION, PLC, PLD, PLFE, PLHE, PLTT, PRLL, SCH, SION, TD, TE, TN, TNU, TRRPC, TRRPFE, TSL


      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            DO NR=1,NRMAX
               PRSUM(NR)=PRLU(1,NR)
            ENDDO
         ELSE
            TSL=DT*DBLE(NT)
            IF(KUFDEV.EQ.'X') THEN
               DO NR=1,NRMAX
                  CALL TIMESPL(TSL,PRLL,TMU,PRLU(1,NR),NTXMAX,NTUM,IERR)
                  PRSUM(NR)=PRLL
               ENDDO
            ELSE
               DO NR=1,NRMAX
                  IF(PNBI.LT.12.D6) THEN
                     PRSUM(NR)=PRLU(1,NR)
                  ELSE
                     PRSUM(NR)=PRLU(2,NR)
                     IF(NT.EQ.NTAMAX) PRSUM(NR)=PRLU(3,NR)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSEIF(MDLUF.EQ.2) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            DO NR=1,NRMAX
               PRSUM(NR)=PRLU(1,NR)
            ENDDO
         ELSE
            TSL=DT*DBLE(NT)
            DO NR=1,NRMAX
               PRSUM(NR)=PRLU(1,NR)
            ENDDO
         ENDIF
      ELSE
         DO NR=1,NRMAX
            ANE =RN(NR,1)
            ANDX=RN(NR,2)
            ANT =RN(NR,3)
            ANHE=RN(NR,4)
            TE  =RT(NR,1)
!     Radiation loss caused by impurities
            PLFE  = ANE*ANFE(NR)*TRRPFE(TE)*1.D40
            PLC   = ANE*ANC (NR)*TRRPC (TE)*1.D40
            PRL(NR)=PLFE+PLC
!     Bremsstrahlung
            PLD   = ANE*ANDX*5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
            PLTT  = ANE*ANT *5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
            PLHE  = ANE*ANHE*5.35D-37*2.D0**2*SQRT(ABS(TE))*1.D40
            PRB(NR)= PLD+PLTT+PLHE
!     Sumup
            SELECT CASE(MDLPR)
            CASE(0)
               PRSUM(NR)=PRL(NR)+PRB(NR)
            CASE(1,2)
               PRSUM(NR)=PRL(NR)+PRB(NR)+PRC(NR)
            END SELECT
         ENDDO
      ENDIF

!     ****** IONIZATION LOSS ******

      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         EION  = 13.64D0
!         SIONS  = 1.76D-13/SQRT(TE*1.D3)*FKN(EION/(TE*1.D3))
         TN    = MAX(TE*1.D3/EION,1.D-2)
         SION  = 1.D-11*SQRT(TN)*EXP(-1.D0/TN)/(EION**1.5D0*(6.D0+TN))
         PIE(NR) = ANE*ANNU(NR)*SION*1.D40*EION*AEE
         SIE(NR) = ANE*ANNU(NR)*SION*1.D20
         TSIE(NR)= ANE         *SION*1.D20
      ENDDO

!     ****** CHARGE EXCHANGE LOSS ******

      DO NR=1,NRMAX
         ANE =RN(NR,1)
         ANDX=RN(NR,2)
         TE  =RT(NR,1)
         TD  =RT(NR,2)
         TNU = 0.D0

         TN    = LOG10(MAX(TD*1.D3,50.D0))
         SCH= 1.57D-16*SQRT(ABS(TD)*1.D3)*(TN*TN-14.63D0*TN+53.65D0)
         SCX(NR) =ANDX*ANNU(NR)*SCH*1.D20
         TSCX(NR)=ANDX         *SCH*1.D20

         PCX(NR)=1.5D0*ANDX*ANNU(NR)*SCH*(TD-TNU)*RKEV*1.D40
      ENDDO

      RETURN
      END SUBROUTINE TRLOSS

      REAL(8) FUNCTION FKN(X)

      IMPLICIT NONE
      REAL(8) X,GAMMA

      GAMMA=0.577215664901532D0
      FKN=-(LOG10(X)+GAMMA-X+X**2/4.D0-X**3/1.8D1)

      RETURN
      END FUNCTION FKN
