
c     $Id$
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
     &     ,RSUMF(NSAM),RSUMF0(NSAM),DEPS2(NSAM)
C
      IF(MODELE.NE.0) CALL FPNEWE

C     +++++ Time loop +++++
c      open(8,file='dfdt_n10b_comp06.dat')
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

         NCHECK=0
         DEPS=1.D0
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS2(NTH,NP,NR,NS)=FNS(NTH,NP,NR,NS)
                  END DO
               END DO
            END DO
         END DO

         DO NSA=1,NSAMAX
            DEPS2(NSA)=1.D0
         END DO
         DO WHILE(DEPS.gt.EPSFP.and.NCHECK.le.LMAXFP)
            NCHECK=NCHECK+1
            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               DO NR=1,NRMAX
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNS(NTH,NP,NR,NS)
               END DO
               END DO
               END DO

               IF(DEPS2(NSA).ge.EPSFP)THEN
c                  write(*,*)NSA,EPSFP,DEPS2(NSA)
                  CALL FPEXEC(NSA,IERR)
                  IF(IERR.NE.0) GOTO 250
               ELSE
                  DO NR=1,NRMAX
                  DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     F1(NTH,NP,NR)=FNS2(NTH,NP,NR,NS)
                  END DO
                  END DO
                  END DO
               END IF

               RSUMF(NSA)=0.D0
               RSUMF0(NSA)=0.D0
               DO NR=1,NRMAX
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNS1(NTH,NP,NR,NS)=F1(NTH,NP,NR)
                  RSUMF(NSA)=RSUMF(NSA)
     &                   +ABS(FNS1(NTH,NP,NR,NS)-FNS2(NTH,NP,NR,NS))**2
                  RSUMF0(NSA)=RSUMF0(NSA)
     &                   +ABS(FNS2(NTH,NP,NR,NS))**2
               ENDDO
               ENDDO
               ENDDO
            ENDDO
            DEPS=0.D0
            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               DEPS1=RSUMF(NSA)/RSUMF0(NSA)
               DEPS2(NSA)=RSUMF(NSA)/RSUMF0(NSA)
               IF(DEPS1.ge.DEPS) NSTEST=NS
               DEPS=MAX(DEPS,DEPS1)
            END DO

      write(6,1274)DEPS,NCHECK,NT,NSTEST,(RSUMF(NSA)/RSUMF0(NSA),
     &             NSA=1,NSAMAX)

C     +++++ update velocity distribution function +++++

            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               DO NR=1,NRMAX
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNS2(NTH,NP,NR,NS)=FNS(NTH,NP,NR,NS)
                  FNS(NTH,NP,NR,NS)=FNS1(NTH,NP,NR,NS)
               ENDDO
               ENDDO
               ENDDO
            ENDDO

C     +++++ end of NSA loop +++++

            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               IF(DEPS2(NSA).ne.0.D0)THEN
                  IF (MOD(NT,NTSTPC).EQ.0) CALL FPCOEF(NSA)
               END IF
            END DO

            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               DO NR=1,NRMAX
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNS(NTH,NP,NR,NS)=FNS2(NTH,NP,NR,NS)
                  FNS2(NTH,NP,NR,NS)=FNS1(NTH,NP,NR,NS)
               ENDDO
               ENDDO
               ENDDO
            ENDDO
         
         END DO
C END OF DOWHILE
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

 1274    FORMAT("DEPS",E14.6,3I3,3E14.6)

c         NP=2
c            Write(8,646) NT
c     & ,DCPP(2,NP,1,1),DCPP(2,NP,1,2),DCPP(2,NP,1,3)
c     & ,DCPT(2,NP,1,1),DCPT(2,NP,1,2),DCPT(2,NP,1,3)
c     & ,DCTT(2,NP,1,1),DCTT(2,NP,1,2),DCTT(2,NP,1,3)
c     & ,FCPP(2,NP,1,1),FCPP(2,NP,1,2),FCPP(2,NP,1,3)


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
c         open(8,file='radial_profile.dat')
c         open(8,file='coef.dat')
         open(9,file='power_test.dat')
c         Do NP=1,NPMAX
c            Write(8,646) NP
c     & ,FNS(14,NP,1,3),FNS(15,NP,1,3),FNS(16,NP,1,3)
c     & ,FNS(17,NP,1,3),FNS(18,NP,1,3),FNS(19,NP,1,3)
c     & ,FNS(20,NP,1,3),FNS(21,NP,1,3),FNS(22,NP,1,3)
c     & ,FNS(23,NP,1,3),FNS(24,NP,1,3),FNS(25,NP,1,3)

c     & ,DCPP(2,NP,1,1),DCPP(2,NP,1,2),DCPP(2,NP,1,3)
c     & ,DCPT(2,NP,1,1),DCPT(2,NP,1,2),DCPT(2,NP,1,3)
c     & ,DCTT(2,NP,1,1),DCTT(2,NP,1,2),DCTT(2,NP,1,3)
c     & ,FCPP(2,NP,1,1),FCPP(2,NP,1,2),FCPP(2,NP,1,3)
c         END DO

c         DO NR=1,NRMAX
c            WRITE(8,645) RM(NR),RPCT2(NR,1,1,NTG1),RPCT2(NR,2,1,NTG1),
c     &         RPCT2(NR,3,1,NTG1),RPCT2(NR,1,2,NTG1),RPCT2(NR,2,2,NTG1),
c     &         RPCT2(NR,3,2,NTG1),RPCT2(NR,1,3,NTG1),RPCT2(NR,2,3,NTG1),
c     &         RPCT2(NR,3,3,NTG1),
c     &         RPWT(NR,1,NTG1),RPWT(NR,2,NTG1),RPWT(NR,3,NTG1)
c         END DO

c         close(8)
c         open(8,file='coef_theta.dat')
c         Do NTH=1,NTHMAX
c            Write(8,646) NTH
c     & ,DCPP(NTH,2,1,1),DCPP(NTH,2,1,2),DCPP(NTH,2,1,3)
c     & ,DCPT(NTH,2,1,1),DCPT(NTH,2,1,2),DCPT(NTH,2,1,3)
c     & ,DCTT(NTH,2,1,1),DCTT(NTH,2,1,2),DCTT(NTH,2,1,3)
c     & ,FCPP(NTH,2,1,1),FCPP(NTH,2,1,2),FCPP(NTH,2,1,3)
c         END DO
c         close(8)

         DO NTI=1,NTMAX
            dw=PPCT(3,NTI)+PPWT(3,NTI)+PPET(3,NTI)
            WRITE(9,645) PTG(NTI)*1000
     &           ,PPCT2(1,1,NTI),PPCT2(2,1,NTI),PPCT2(3,1,NTI)
     &           ,PPCT2(1,2,NTI),PPCT2(2,2,NTI),PPCT2(3,2,NTI)
     &           ,PPCT2(1,3,NTI),PPCT2(2,3,NTI),PPCT2(3,3,NTI)
     &           ,PPWT(1,NTI),PPWT(2,NTI),PPWT(3,NTI)
     &           ,PWT(1,NTI),PWT(2,NTI),PWT(3,NTI)
c<<<<<<< fploop.f
cccccccccccccc
c     &           ,PNT(3,NTI),PNT(1,NTI),PNT(2,NTI)
c     &           ,dw/1000.D0,PWT(3,NTI)-PWT(3,NTI-1)
c     &           ,(PWT(3,NTI+1)-PWT(3,NTI-1))/2.D0 
c     &           ,PWT(3,NTI+1)-PWT(3,NTI)
c     &    , (PWT(3,NTI+1)-PWT(3,NTI-1)-(PWT(3,NTI+2)-PWT(3,NTI-2))/8.D0)
c     &           /1.5D0

c密度のズレ補正
c     &           ,(PWT(3,NTI)*PNT(3,NTI)-PWT(3,NTI-1)*PNT(3,NTI-1))
c     &           /(PPWT(3,NTI)*PNT(3,NTI)+PPCT(3,NTI)*PNT(3,NTI)+
c     &           PPET(3,NTI)*PNT(3,NTI))

     &           ,(PWT(3,NTI)-PWT(3,NTI-1))
     &           /(PPWT(3,NTI)+PPCT(3,NTI)+PPET(3,NTI))

c     &           ,(PWT(3,NTI+1)*PNT(3,NTI+1)-PWT(3,NTI)*PNT(3,NTI))
c     &           /(PPWT(3,NTI)*PNT(3,NTI)+PPCT(3,NTI)*PNT(3,NTI)+
c     &           PPET(3,NTI)*PNT(3,NTI))

c     &           ,( PWT(3,NTI+1)*PNT(3,NTI+1)-PWT(3,NTI-1)*PNT(3,NTI-1)- 
c     &      (PWT(3,NTI+2)*PNT(3,NTI+2)-PWT(3,NTI-2)*PNT(3,NTI-2))/8.D0)
c     &           /(PPWT(3,NTI)*PNT(3,NTI)+PPCT(3,NTI)*PNT(3,NTI)+
c     &           PPET(3,NTI)*PNT(3,NTI))/1.5D0
c=======
c>>>>>>> 1.52

         END DO
         close(9)
      END IF
 645  FORMAT(17E14.6)
 646  FORMAT(I3,17E14.6)
         IF(IERR.NE.0) RETURN
      ENDDO
      close(8)
C     +++++ end of time loop +++++
C
      IF(MODELA.eq.1)THEN
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS(NTH,NP,NR,NS) = FNS(NTH,NP,NR,NS)/RCOEF(NR)
                  END DO
               END DO
            END DO
         END DO
      END IF
      
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
