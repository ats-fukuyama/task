C     $Id$
C
C *****************
C     MAIN LOOP
C *****************
C
      SUBROUTINE FPLOOP
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION RJNS(NRM,NSAM),RJN(NRM),RJ3(NRM),E3(NRM),DELE(NRM)
      
C
      IF(MODELE.NE.0) CALL FPNEWE

C     +++++ Time loop +++++

      DO NT=1,NTMAX

C     +++++ Iteration loop for toroidal electric field +++++

         L=0
C
         IF(MODELE.NE.0) THEN
            DO NR=1,NRMAX
               E3(NR)=0.D0
               RJ3(NR)=0.D0
            ENDDO
         ENDIF
C
    1    L=L+1

C     +++++ NSA loop +++++


c         IF(NT.eq.NTMAX)THEN
c         DO NSA=1,2
c         DO NSB=1,2
c          Do NP=1,NPMAX
c         NSA=1
c         NSB=1
c         NP=2
c         write(*,*) " "
c            write(*,1600) Nt,NSa,PG(NP),
c     &           DPP(2,NP,1,NSA),DCTT2(2,NP,1,NSB,NSA),
c     &           DCPT2(2,NP,1,NSB,NSA),FCPP2(2,NP,1,NSB,NSA),
c     &           FCTH2(2,NP,1,NSB,NSA)
c          END DO
c         write(*,*) " "
c         write(7,*) " "
c         END DO
c         END DO
c         END IF
c
c 1600    FORMAT(2I2,6E14.6)

         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     F(NTH,NP,NR)=FNS(NTH,NP,NR,NS)
                  END DO
               END DO
            END DO

            CALL FPEXEC(NSA,IERR)
            IF(IERR.NE.0) GOTO 250

            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS1(NTH,NP,NR,NS)=F1(NTH,NP,NR) 
                  ENDDO
               ENDDO
            ENDDO
         ENDDO


C     +++++ update velocity distribution function +++++

         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            FNS(NTH,NP,NR,NS)=FNS1(NTH,NP,NR,NS)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

C     +++++ end of NSA loop +++++

         IF (MOD(NT,NTSTPC).EQ.0) CALL FPCOEF

C     ----- calculation of current density -----

         IF(MODELE.NE.0) THEN
            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               DO NR=2,NRMAX
                  RSUM=0.D0
                  DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM=RSUM+VOL(NTH,NP)*FNS1(NTH,NP,NR,NS)*PM(NP)
                  ENDDO
                  ENDDO
                  RJNS(NR,NSA)=AEFP(NSA)*RNFP0(NSA)*1.D20
     &                    *PTFP0(NSA)*DELP*RSUM/(AMFP(NSA)*RM(NR)*RA)
               ENDDO
               RJNS(1,NSA)=(4.D0*RJNS(2,NSA)-RJNS(3,NSA))/3.D0
            ENDDO

C     ----- calculation of toroidal electric field -----

            DELEM=0.D0
            DO NR=1,NRMAX
               RJNL=0.D0
               DO NSA=1,NSAMAX
                  RJNL=RJNL+RJNS(NR,NSA)
               END DO
               RJN(NNR)=RJNL
               IF(ABS(RJNL-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJNL-RJ2(NR))*(E2(NR)-E3(NR))
     &                    /(RJNL-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJNL
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
            ENDDO

            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 1
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'

            DO NR=1,NRMAX
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
            ENDDO
            CALL FPNEWE
         ENDIF
C     +++++ end of toroidal electric field loop +++++
C
  250    CONTINUE
C
C         CALL FPGRAC('F -2',F,4)
C         CALL FPGRAC('F1-2',F1,4)


C     +++++ calculate and save global data +++++

         TIMEFP=TIMEFP+DELT

         ISAVE=0
         IF (MOD(NT,NTSTP2).EQ.0) THEN
            CALL FPSGLB
            CALL FPWRT2
         ENDIF
         IF (MOD(NT,NTSTP1).EQ.0) THEN
            CALL FPSPRF
            CALL FPWRT1
         ENDIF

      call FPWRT3


      IF(NT.eq.NTMAX)THEN
         open(8,file='radial_profile.dat')
         DO NR=1,NRMAX
            WRITE(8,645) RPCT2(NR,1,1,NTG1),RPCT2(NR,2,1,NTG1),
     &         RPCT2(NR,3,1,NTG1),RPCT2(NR,1,2,NTG1),RPCT2(NR,2,2,NTG1),
     &         RPCT2(NR,3,2,NTG1),RPCT2(NR,1,3,NTG1),RPCT2(NR,2,3,NTG1),
     &         RPCT2(NR,3,3,NTG1)
         END DO
c         write(8,*)" " 
c         write(8,*)" " 
      END IF
 645  FORMAT(9E14.6)
         IF(IERR.NE.0) RETURN
      ENDDO
C     +++++ end of time loop +++++
C
      RETURN
      END
C
C ************************************
C     PREDICTION OF ELECTRIC FIELD
C ************************************
C

      SUBROUTINE FPNEWE
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=2,NRMAX
         BP(NR)=BP(NR)+(E1(NR)-E1(NR-1))*DELT/(RA*DELR)
      ENDDO
C
      DO NR=1,NRMAX
         RJ2(NR)=(RG(NR+1)*BP(NR+1)-RG(NR)*BP(NR))
     &           /(RMU0*RM(NR)*DELR*RA)
      ENDDO
C
      DO NR=1,NRMAX
         E2(NR)=RJ2(NR)*E1(NR)/RJ1(NR)
      ENDDO
C
      RETURN
      END
