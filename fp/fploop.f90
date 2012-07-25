!     $Id$

! *****************
!     MAIN LOOP
! *****************

      MODULE fploop

      USE fpcomm
      use fpexec
      use fpdrexec
      use fpcoef
      use fpsave
      contains

!-----------------------------

      SUBROUTINE FP_LOOP

      USE libmtx
      IMPLICIT NONE
      real(kind8),dimension(NRSTART:NREND,NSAMAX):: RJNS
      real(kind8),dimension(NRSTART:NREND):: RJN,RJ3,E3,DELE
      real(kind8),dimension(NSAMAX)::RSUMF,RSUMF0,RSUM_SS!,DEPS2

      integer:: NT, NR, NP, NTH, NSA, NTI, NSBA
      integer:: L, IERR, I!, N_IMPL
      real(kind8):: RSUM, DELEM, RJNL, dw, RSUM1, RSUM2
      real(4):: gut, gut1, gut2, gut3, gut4, gut5, gut6
      real(4):: gut_ex, gut_calc, gut_1step
      real(kind8),DIMENSION(nprocs):: RSUMA
      integer,dimension(NSAMAX)::NTCLSTEP2
      real(kind8):: DEPS_MAX, DEPS, DEPS1
      real(kind8),dimension(NSAMAX):: DEPS_MAXV, DEPSV
      integer,dimension(NSAMAX):: ILOC
      integer:: NSTEST
      real(kind8):: temp_send, temp_recv
      character:: fmt*40
      integer:: modeld_temp, modela_temp

      IF(MODELE.NE.0) CALL FPNEWE

!     +++++ Time loop +++++

!      IF(NRANK.eq.29)THEN
!         DO NP=1,NPMAX
!            DO NTH=1,NTHMAX
!               diff_f(NTH,NP)=FNS(NTH,NP,30,1)
!            END DO
!         END DO
!      END IF


      DO NT=1,NTMAX
         
!     +++++ Iteration loop for toroidal electric field +++++

         L=0

         IF(MODELE.NE.0) THEN
            DO NR=NRSTART,NREND
               E3(NR)=0.D0
               RJ3(NR)=0.D0
            ENDDO
         ENDIF

    1    L=L+1

         N_IMPL=0
         DEPS=1.D0
         DO NSA=1,NSAMAX
            NSBA=NSB_NSA(NSA)
            DO NR=NRSTART,NREND ! local
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS2(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSBA) ! old step: invariant during N_IMPL 
                     FNS1(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSBA) ! new step: variant each N_IMPL
                  END DO
               END DO
            END DO
            DO NR=1,NRMAX ! global
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS22(NTH,NP,NR,NSA)=FNS(NTH,NP,NR,NSBA) ! old step
                  END DO
               END DO
            END DO
         END DO

         gut_EX = 0.D0
         gut_CALC= 0.D0
         gut_1step= 0.D0
         CALL GUTIME(gut5)

!         DO WHILE(DEPS.gt.EPSFP.and.N_IMPL.le.LMAXFP) ! start do while
         DO WHILE(N_IMPL.le.LMAXFP) ! start do while
            N_IMPL=N_IMPL+1
            DO NSA=1,NSAMAX 
               NSBA=NSB_NSA(NSA)
               DO NR=NRSTART,NREND
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNS2(NTH,NP,NR,NSBA)
               END DO
               END DO
               END DO

               CALL GUTIME(gut1)
               modeld_temp=modeld
               modeld=0
               CALL fp_exec(NSA,IERR) ! F1 and FNS are changed
               modeld=modeld_temp
               IF(MODELD.ge.1)THEN
                  CALL fp_drexec(NSA,IERR)
               END IF

               IF(IERR.NE.0) GOTO 250
               CALL GUTIME(gut2)
               GUT_EX = GUT_EX + (gut2-gut1)

               RSUMF(NSA)=0.D0
               RSUMF0(NSA)=0.D0
               RSUM_SS(NSA)=0.D0
               DO NR=NRSTART,NREND
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUMF(NSA)=RSUMF(NSA) &
                         +ABS(FNS1(NTH,NP,NR,NSA)-F1(NTH,NP,NR))**2
                  RSUMF0(NSA)=RSUMF0(NSA) &
                         +ABS(FNS2(NTH,NP,NR,NSA))**2
                  RSUM_SS(NSA)=RSUM_SS(NSA) &
                         +ABS(FNS2(NTH,NP,NR,NSA)-F1(NTH,NP,NR))**2
                  FNS1(NTH,NP,NR,NSA)=F1(NTH,NP,NR)
               ENDDO
               ENDDO
               ENDDO
            ENDDO ! END OF NSA

            DEPS=0.D0
            DO NSA=1,NSAMAX
               DEPSV(NSA) = RSUMF(NSA)/RSUMF0(NSA)
               DEPS1 = DEPSV(NSA)
               DEPS=MAX(DEPS,DEPS1)
               DEPS_SS(NSA)=RSUM_SS(NSA)/RSUMF0(NSA)/DELT ! steady state
            END DO

            DEPS_MAX=0.D0
            CALL mtx_reduce_real8(DEPS,1,DEPS_MAX)
            DEPS = DEPS_MAX
            IF(DEPS.le.EPSFP)THEN
               N_IMPL=1+LMAXFP ! exit dowhile
            ENDIF

            CALL mtx_maxloc_real8(DEPSV,NSAMAX,DEPS_MAXV,ILOC)

            IF(nrank.eq.0) THEN
               WRITE(fmt,'(a16,I1,a6,I1,a3)') &
                    '(A,1PE12.4,I4,1P',NSAMAX,'E12.4,',NSAMAX,'I4)'
               WRITE(6,fmt) 'DEPS',&
                    DEPS,N_IMPL,(DEPS_MAXV(NSA),NSA=1,NSAMAX) &
                    ,(ILOC(NSA),NSA=1,NSAMAX)
            ENDIF

            CALL GUTIME(gut3)
            DO NSA=1,NSAMAX
                  IF (MOD(NT,NTCLSTEP).EQ.0) CALL FP_COEF(NSA)
            END DO
            CALL GUTIME(gut4)
            GUT_CALC = GUT_CALC + (gut4-gut3)

         END DO ! END OF DOWHILE
         CALL GUTIME(gut6)
         GUT_1step = gut6-gut5

         IF(NRANK.eq.0) &
              WRITE(*,'(" GUT_EX = ", E14.6, " GUT_CALC = ", E14.6, " GUT_1step = ", E14.6)') &
              GUT_EX, GUT_CALC, GUT_1step

!     ----- calculation of current density -----

         IF(MODELE.NE.0) THEN
            DO NSA=1,NSAMAX
            NSBA=NSB_NSA(NSA)
               DO NR=2,NRMAX
                  RSUM=0.D0
                  DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM=RSUM+VOLP(NTH,NP,NSBA)*FNS(NTH,NP,NR,NSA)*PM(NP,NSBA)
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

!         IF(NT.eq.NTMAX.or.NTMAX.ne.0)THEN
!!            open(9,file='power_D_5s_D0_taul1000_2kev_NB.dat')
!!            DO NTI=1,NTG1
!!               WRITE(9,645) NTI, PTG(NTI)*1000 &
!!                    ,PPCT2(1,1,NTI),PPCT2(2,1,NTI),PPCT2(3,1,NTI),PPCT2(4,1,NTI),PPCT(1,NTI) &
!!                    ,PPCT2(1,2,NTI),PPCT2(2,2,NTI),PPCT2(3,2,NTI),PPCT2(4,2,NTI),PPCT(2,NTI) &
!!                    ,PPCT2(1,3,NTI),PPCT2(2,3,NTI),PPCT2(3,3,NTI),PPCT2(4,3,NTI),PPCT(3,NTI) &
!!                    ,PPCT2(1,4,NTI),PPCT2(2,4,NTI),PPCT2(3,4,NTI),PPCT2(4,4,NTI),PPCT(4,NTI) &
!!                    ,PPWT(1,NTI),PPWT(2,NTI),PPWT(3,NTI),PPWT(4,NTI) &
!!                    ,PDR(1,NTI),PDR(2,NTI),PDR(3,NTI),PDR(4,NTI) &
!!                    ,PWT(1,NTI),PWT(2,NTI),PWT(3,NTI),PWT(4,NTI) &
!!                    ,PNT(1,NTI),PNT(2,NTI),PNT(3,NTI),PNT(4,NTI) &
!!                    ,PTT2(1,NTI),PTT2(2,NTI),PTT2(3,NTI),PTT2(4,NTI) &
!!                    ,PTT_BULK(1,NTI),PTT_BULK(2,NTI),PTT_BULK(3,NTI),PTT_BULK(4,NTI) &
!!                    ,PSPBT(2,NTI),PSPFT(2,NTI),PSPFT(3,NTI),PSPFT(4,NTI) &
!!                   ,PECT(1,NTI)
!!            END DO
!!            close(9)
!         END IF

         IF(IERR.NE.0) RETURN


!         CALL GUTIME(gut2)
!         write(*,*)"1 loop time", nrank, gut2-gut

      ENDDO ! END OF NT LOOP

      IF(NRANK.eq.0)THEN
            open(9,file='fns_e_a1w4s2_pinch.dat')
            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     WRITE(9,'(3I4,3E16.8)') NR, NP, NTH, FNS(NTH,NP,NR,1), PM(NP,1)*COSM(NTH), PM(NP,1)*SINM(NTH)
                  END DO
               END DO
               WRITE(9,*) " "
               WRITE(9,*) " "
            END DO
            close(9)
!            open(9,file='fns_D.dat')
!            DO NR=1,NRMAX
!               DO NP=1,NPMAX
!                  DO NTH=1,NTHMAX
!                     WRITE(9,'(3I4,E16.8)') NR, NP, NTH, FNS(NTH,NP,NR,2)
!                  END DO
!               END DO
!               WRITE(9,*) " "
!               WRITE(9,*) " "
!            END DO
!            close(9)

!            open(9,file='fns_T.dat')
!            DO NR=1,NRMAX
!               DO NP=1,NPMAX
!                  DO NTH=1,NTHMAX
!                     WRITE(9,'(3I4,E16.8)') NR, NP, NTH, FNS(NTH,NP,NR,3)
!                  END DO
!               END DO
!               WRITE(9,*) " "
!               WRITE(9,*) " "
!            END DO
!            close(9)
!            open(9,file='fns_He.dat')
!            DO NR=1,NRMAX
!               DO NP=1,NPMAX
!                  DO NTH=1,NTHMAX
!                     WRITE(9,'(3I4,E16.8)') NR, NP, NTH, FNS(NTH,NP,NR,4)
!                  END DO
!               END DO
!               WRITE(9,*) " "
!               WRITE(9,*) " "
!            END DO
!            close(9)
      END IF

!      IF(NRANK.eq.29)THEN
!      open(8,file='diff_FNS.dat')
!      DO NP=1,NPMAX
!         DO NTH=1,NTHMAX
!            diff_f(NTH,NP)=FNS(NTH,NP,30,1)-diff_f(NTH,NP)
!!            WRITE(8,*) NTH, NP, diff_f(NTH,NP)
!            WRITE(8,*) PM(NP,1)*COSM(NTH), PM(NP,1)*SINM(NTH), diff_f(NTH,NP)
!         END DO 
!         WRITE(8,*) " "
!         WRITE(8,*) " "
!      END DO
!      CLOSE(8)
!      END IF


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
