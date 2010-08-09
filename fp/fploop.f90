!     $Id$

! *****************
!     MAIN LOOP
! *****************

      MODULE fploop

      USE fpcomm
      use fpexec
      use fpcoef
      use fpsave
      contains

!-----------------------------

      SUBROUTINE FP_LOOP

      USE libmtx
      IMPLICIT NONE
      real(8),dimension(NRSTART:NREND,NSAMAX):: RJNS
      real(8),dimension(NRSTART:NREND):: RJN,RJ3,E3,DELE
      real(8),dimension(NSAMAX)::RSUMF,RSUMF0!,DEPS2

      integer:: NT, NR, NP, NTH, NSA, NTI, NSBA
      integer:: L, IERR, NSTEST, I!, NCHECK
      real(8):: DEPS, DEPS1, RSUM, DELEM, RJNL, dw, RSUM1, RSUM2
      real(4):: gut, gut2, gut3, gut4
      real(8),DIMENSION(nprocs):: RSUMA

      IF(MODELE.NE.0) CALL FPNEWE

!     +++++ Time loop +++++
!      CALL GUTIME(gut3)

      DO NT=1,NTMAX
!         CALL GUTIME(gut)
!         write(*,*)"time1", gut, NSA
         
!     +++++ Iteration loop for toroidal electric field +++++

         L=0

         IF(MODELE.NE.0) THEN
            DO NR=NRSTART,NREND
               E3(NR)=0.D0
               RJ3(NR)=0.D0
            ENDDO
         ENDIF

    1    L=L+1

         NCHECK=0
         DEPS=1.D0
         DO NSA=1,NSAMAX
            NSBA=NSB_NSA(NSA)
            DO NR=NRSTART,NREND
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS2(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSBA)
                     FNS1(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSBA)
                  END DO
               END DO
            END DO
            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS22(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSBA) ! before
                  END DO
               END DO
            END DO
!            DEPS2(NSA)=1.D2
         END DO

         DO WHILE(DEPS.gt.EPSFP.and.NCHECK.le.LMAXFP) ! start do while
            NCHECK=NCHECK+1
            DO NSA=1,NSAMAX 
               NSBA=NSB_NSA(NSA)
               DO NR=NRSTART,NREND
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNS2(NTH,NP,NR,NSBA)
               END DO
               END DO
               END DO

               IF(DEPS.ge.EPSFP)THEN
                  CALL fp_exec(NSA,IERR) ! F1 and FNS are changed
                  IF(IERR.NE.0) GOTO 250
               ELSE
                  DO NR=NRSTART,NREND
                  DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     F1(NTH,NP,NR)=FNS1(NTH,NP,NR,NSA)
                     FNS(NTH,NP,NR,NSA)=FNS1(NTH,NP,NR,NSA)
                  END DO
                  END DO
                  END DO
               END IF

!test
!            NR=1
!            RSUM1=0.D0
!            RSUM2=0.D0
!            DO NP=1,NPMAX
!              DO NTH=1,NTHMAX
!                  RSUM1=RSUM1+VOLP(NTH,NP,NSA)*FNS(NTH,NP,NR,NSA)
!                  RSUM2=RSUM2+VOLP(NTH,NP,NSA)*FNS2(NTH,NP,NR,NSA)
!               END DO
!            END DO
!            DO NP=1,NPMAX
!               DO NTH=1,NTHMAX
!                  FNS(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSA)*RSUM2/RSUM1
!                  F1(NTH,NP,NR)=F1(NTH,NP,NR)*RSUM2/RSUM1
!!                  write(*,*) NP, NTH, RSUM2/RSUM1
!               END DO
!            END DO
!!end of test

               RSUMF(NSA)=0.D0
               RSUMF0(NSA)=0.D0
               DO NR=NRSTART,NREND
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUMF(NSA)=RSUMF(NSA) &
                         +ABS(FNS1(NTH,NP,NR,NSA)-FNS(NTH,NP,NR,NSA))**2
                  RSUMF0(NSA)=RSUMF0(NSA) &
                         +ABS(FNS2(NTH,NP,NR,NSA))**2
                  FNS1(NTH,NP,NR,NSA)=F1(NTH,NP,NR)
               ENDDO
               ENDDO
               ENDDO
            ENDDO ! END OF NSA

            DO NSA=1,NSAMAX
               CALL mtx_gather_real8(RSUMF(NSA),RSUMA)
               RSUMF(NSA)=0.D0
               DO i=1,nprocs
                  RSUMF(NSA)=RSUMF(NSA)+RSUMA(i)
               ENDDO
               RSUMA(1)=RSUMF(NSA)
               CALL mtx_broadcast_real8(RSUMA,1)
               RSUMF(NSA)=RSUMA(1)

               CALL mtx_gather_real8(RSUMF0(NSA),RSUMA)
               RSUMF0(NSA)=0.D0
               DO i=1,nprocs
                  RSUMF0(NSA)=RSUMF0(NSA)+RSUMA(i)
               ENDDO
               RSUMA(1)=RSUMF0(NSA)
               CALL mtx_broadcast_real8(RSUMA,1)
               RSUMF0(NSA)=RSUMA(1)
            ENDDO

            DEPS=0.D0
            NSTEST=0
            DO NSA=1,NSAMAX
               DEPS1=RSUMF(NSA)/RSUMF0(NSA)
               IF(DEPS1.ge.DEPS) NSTEST=NSA 
               DEPS=MAX(DEPS,DEPS1)
            END DO
            IF(DEPS.le.EPSFP)THEN ! exit do while
               NCHECK=NCHECK+LMAXFP
            END IF
            IF(nrank.eq.0) THEN
               write(6,'(A,1PE12.4,3I4,1P4E12.4)') 'DEPS',&
                    DEPS,NCHECK,NREND,NSTEST,(RSUMF(NSA)/RSUMF0(NSA), &
                    NSA=1,NSAMAX)
            ENDIF


            DO NSA=1,NSAMAX
!            CALL GUTIME(gut)
!            write(*,*)"time1", gut, NSA
                  IF (MOD(NT,NTCLSTEP).EQ.0) CALL FP_COEF(NSA)
!            CALL GUTIME(gut)
!            write(*,*)"time2", gut, NSA
            END DO

         END DO ! END OF DOWHILE

!         open(8,file='dfdt_r0c1.dat')
         DO NSA=1,NSAMAX
            NSBA=NSB_NSA(NSA)
            DO NR=NRSTART,NREND
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               FNS(NTH,NP,NR,NSBA)=FNS1(NTH,NP,NR,NSA)
!               IF(NR.eq.1.and.NSA.eq.1)THEN
!                  WRITE(8,'(1P14E14.6)') PM(NP,NSBA)*COSM(NTH), PM(NP,NSBA)*SINM(NTH), &
!                       ( FNS(NTH,NP,NR,NSA)-FNS22(NTH,NP,NR,NSA) )/FNS22(NTH,NP,NR,NSA)
!               END IF
            ENDDO
!            WRITE(8,*) " "
!            WRITE(8,*) " "
            ENDDO
            ENDDO
         ENDDO
!         close(8)
!     ----- calculation of current density -----

         IF(MODELE.NE.0) THEN
            DO NSA=1,NSAMAX
            NSBA=NSB_NSA(NSA)
               DO NR=2,NRMAX
                  RSUM=0.D0
                  DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM=RSUM+VOLP(NTH,NP,NSBA)*FNS1(NTH,NP,NR,NSA)*PM(NP,NSBA)
                  ENDDO
                  ENDDO
                  RJNS(NR,NSA)=AEFP(NSA)*RNFP0(NSA)*1.D20 &
                          *PTFP0(NSA)*DELP(NSBA)*RSUM/(AMFP(NSA)*RM(NR)*RA)
               ENDDO
               RJNS(1,NSA)=(4.D0*RJNS(2,NSA)-RJNS(3,NSA))/3.D0
            ENDDO

!     ----- calculation of toroidal electric field -----

            DELEM=0.D0
            DO NR=NRSTART,NREND
               RJNL=0.D0
               DO NSA=1,NSAMAX
                  RJNL=RJNL+RJNS(NR,NSA)
               END DO
               RJN(NR)=RJNL
               IF(ABS(RJNL-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJNL-RJ2(NR))*(E2(NR)-E3(NR)) &
                          /(RJNL-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJNL
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
            ENDDO

            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 1
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'

            DO NR=NRSTART,NREND
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
            ENDDO
            CALL FPNEWE
         ENDIF
!     +++++ end of toroidal electric field loop +++++
!
  250    CONTINUE
!
!         CALL FPGRAC('F -2',F,4)
!         CALL FPGRAC('F1-2',F1,4)


!     +++++ calculate and save global data +++++

         TIMEFP=TIMEFP+DELT

         ISAVE=0
         IF (MOD(NT,NTG1STEP).EQ.0) THEN
            CALL FPSSUB
            IF(NRANK.EQ.0) THEN
               CALL FPSGLB
               CALL FPWRTGLB
            ENDIF
         ENDIF
         IF (MOD(NT,NTG2STEP).EQ.0) THEN
            CALL FPSSUB
            IF(NRANK.EQ.0) THEN
               CALL FPSPRF
               CALL FPWRTPRF
            ENDIF
         ENDIF
         IF(NRANK.EQ.0.AND.NTG1.GT.0) call FPWRTSNAP

         IF(NT.eq.NTMAX.or.NTMAX.eq.0)THEN
            open(9,file='power_SNA_1s_D6.dat')

!       ,DCPP(2,NP,1,1),DCPP(2,NP,1,2),DCPP(2,NP,1,3),DCPP(2,NP,1,4) &
!       ,DCPT(2,NP,1,1),DCPT(2,NP,1,2),DCPT(2,NP,1,3),DCPT(2,NP,1,4) &
!       ,DCTT(2,NP,1,1),DCTT(2,NP,1,2),DCTT(2,NP,1,3),DCTT(2,NP,1,4) &
!       ,FCPP(2,NP,1,1),FCPP(2,NP,1,2),FCPP(2,NP,1,3),FCPP(2,NP,1,4)
!       ,DCPP2(2,NP,1,1,1),DCTT2(2,NP,1,1,1),DCTP2(2,NP,1,1,1),DCPT2(2,NP,1,1,1) &
!       ,FCPP2(2,NP,1,1,1) &
!       ,DCPP2(2,NP,1,2,1),DCTT2(2,NP,1,2,1),DCTP2(2,NP,1,2,1),DCPT2(2,NP,1,2,1) &
!       ,FCPP2(2,NP,1,2,1) 
!         END DO

            DO NTI=1,NTG1
               WRITE(9,645) NTI, PTG(NTI)*1000 &
                    ,PPCT2(1,1,NTI),PPCT2(2,1,NTI),PPCT2(3,1,NTI),PPCT2(4,1,NTI) &
                    ,PPCT2(1,2,NTI),PPCT2(2,2,NTI),PPCT2(3,2,NTI),PPCT2(4,2,NTI) &
                    ,PPCT2(1,3,NTI),PPCT2(2,3,NTI),PPCT2(3,3,NTI),PPCT2(4,3,NTI) &
                    ,PPCT2(1,4,NTI),PPCT2(2,4,NTI),PPCT2(3,4,NTI),PPCT2(4,4,NTI) &
                    ,PPWT(1,NTI),PPWT(2,NTI),PPWT(3,NTI),PPWT(4,NTI)             &
                    ,PWT(1,NTI),PWT(2,NTI),PWT(3,NTI),PWT(4,NTI)                 &
                    ,PNT(1,NTI),PNT(2,NTI),PNT(3,NTI),PNT(4,NTI)                 &
                    ,PTT2(1,NTI),PTT2(2,NTI),PTT2(3,NTI),PTT2(4,NTI)                 &
                    ,PSPBT(2,NTI),PSPFT(4,NTI)
            END DO
            close(9)

!         open(8,file='deriv_W_r1c4_12.dat')
!         DO NTI=2,NTG1
!            Write(8,645) NTI, PTG(NTI)*1000 &
!                 ,( PWT(1,NTI)-PWT(1,NTI-1) )/( PTG(NTI)-PTG(NTI-1) ) &
!                 ,( PWT(2,NTI)-PWT(2,NTI-1) )/( PTG(NTI)-PTG(NTI-1) ) &
!                 ,( PWT(3,NTI)-PWT(3,NTI-1) )/( PTG(NTI)-PTG(NTI-1) ) &
!                 ,PPCT(1,NTI)+PPWT(1,NTI)+PSPBT(1,NTI)+PSPFT(1,NTI) 
!                 ,PPCT(2,NTI)+PPWT(2,NTI)+PSPBT(2,NTI)+PSPFT(2,NTI) &
!                 ,PPCT(3,NTI)+PPWT(3,NTI)+PSPBT(3,NTI)+PSPFT(3,NTI) &
!                 ,PWTD(1,NTI),PWTD(2,NTI),PWTD(3,NTI) &
!                 ,( PPCT(1,NTI)+PPWT(1,NTI)+PSPBT(1,NTI)+PSPFT(1,NTI) &
!                 +PPCT(1,NTI-1)+PPWT(1,NTI-1)+PSPBT(1,NTI-1)+PSPFT(1,NTI-1) )*0.5D0 &
!                 ,( PPCT(2,NTI)+PPWT(2,NTI)+PSPBT(2,NTI)+PSPFT(2,NTI) &
!                 +PPCT(2,NTI-1)+PPWT(2,NTI-1)+PSPBT(2,NTI-1)+PSPFT(2,NTI-1) )*0.5D0 &
!                 ,( PPCT(3,NTI)+PPWT(3,NTI)+PSPBT(3,NTI)+PSPFT(3,NTI) &
!                 +PPCT(3,NTI-1)+PPWT(3,NTI-1)+PSPBT(3,NTI-1)+PSPFT(3,NTI-1) )*0.5D0 
!         END DO
!         close(8)
!            open(8,file='PNT_r1c4_2.dat')
!            DO NTI=1,NTMAX
!               WRITE(*,*) NTI, PNT(1,NTI)/PNT(1,1)
!            END DO
!            close(8)

         END IF

 645  FORMAT(I3,40E14.6)
 646  FORMAT(I3,17E14.6)
 647  FORMAT(12E14.6) 
         IF(IERR.NE.0) RETURN


!         CALL GUTIME(gut2)
!         write(*,*)"1 loop time", nrank, gut2-gut

      ENDDO ! END OF NT LOOP

!      CALL GUTIME(gut4)
!      write(*,*)"total loop time", nrank, gut4-gut3
         
!     +++++ end of time loop +++++
      
      RETURN
      END SUBROUTINE FP_LOOP

! ************************************
!     PREDICTION OF ELECTRIC FIELD
! ************************************

      SUBROUTINE FPNEWE

      IMPLICIT NONE
      integer:: NR

      DO NR=2,NRMAX
         BP(NR)=BP(NR)+(E1(NR)-E1(NR-1))*DELT/(RA*DELR)
      ENDDO

      DO NR=NRSTART,NREND
         RJ2(NR)=(RG(NR+1)*BP(NR+1)-RG(NR)*BP(NR)) &
                 /(RMU0*RM(NR)*DELR*RA)
      ENDDO

      DO NR=NRSTART,NREND
         E2(NR)=RJ2(NR)*E1(NR)/RJ1(NR)
      ENDDO

      RETURN
      END SUBROUTINE FPNEWE

      END MODULE fploop
