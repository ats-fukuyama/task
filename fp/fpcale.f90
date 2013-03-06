!
! *************************
!     SAVE DATA ROUTINE
! *************************
!
      MODULE fpcaleind

      USE fpcomm
      USE libmpi
      USE libmtx
      USE plprof
      USE equnit_mod

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE UPDATE_FEPP

      USE fpcoef, only: FP_CALE
      IMPLICIT NONE
      INTEGER:: NSA, NR

!      IF(NRANK.eq.0)THEN
!         WRITE(6,'(A,32E10.2)') "     EM=",(EM(NR),NR=1,32)
!         WRITE(6,'(A,32E10.2)') "     EP=",(EP(NR),NR=1,32)
!      END IF

      DO NSA=NSASTART,NSAEND
         CALL FP_CALE(NSA)
      END DO

      END SUBROUTINE UPDATE_FEPP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE E_IND_EVOL

      IMPLICIT NONE
      INTEGER:: NSA,NR,NP,NTH,IERR
      double precision:: coef0,coefm,coefp,RHS,E_EDGE
      integer:: imtxstart1,imtxend1,its,j

      EPM(:)=EP(:)
      CALL mtx_set_communicator(comm_nr) !3D

      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)

!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
!         coef0=2.D0/(RA*DELR)**2 + RMU0*SIGP(NR)/DELT
         coef0=2.D0 + RMU0*SIGP(NR)/DELT*(RA*DELR)**2
         CALL mtx_set_matrix(NR,NR,coef0)
         CALL mtx_set_vector(NR,EM(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         coefm=-RG(NR)/( RM(NR) )
         coefp=-RG(NR+1)/( RM(NR) )
         IF(NR.ne.1)THEN
            CALL mtx_set_matrix(NR,NR-1,coefm)
         END IF
         IF(NR.ne.NRMAX)THEN
            CALL mtx_set_matrix(NR,NR+1,coefp)
         END IF
      END DO
!---- RIGHT HAND SIDE
      DO NR=NRSTART,NREND
         IF(NR.ne.NRMAX)THEN
            RHS=RMU0/DELT*( RJ_M(NR)-RJ_P(NR)+SIGM(NR)*EPM(NR) )*(RA*DELR)**2
            CALL mtx_set_source(NR,RHS)
         ELSE
            CALL INDUCTANCE_EDGE(E_EDGE)
            RHS=RMU0/DELT*( RJ_M(NR)-RJ_P(NR)+SIGM(NR)*EPM(NR) )*(RA*DELR)**2 &
                 + RG(NRMAX+1)/RM(NRMAX)*E_EDGE
            CALL mtx_set_source(NR,RHS)
            WRITE(*,*) N_IMPL,"E_EDGE=",E_EDGE
         END IF
      END DO

!---- SOLVE

      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 

      CALL mtx_gather_vector(EP)

      CALL mtx_cleanup
      CALL mtx_reset_communicator
!      IF(NRANK.eq.0) WRITE(*,'(I3,A,7E14.6)') N_IMPL," RHS=", &
!           (RMU0/DELT*( RJ_M(NR)-RJ_P(NR)+SIGP(NR)*EP(NR) )*(RA*DELR)**2,NR=1,4),SIGP(1),EPM(1),coef0

      END SUBROUTINE E_IND_EVOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE conductivity_sigma_ind

      IMPLICIT NONE
      INTEGER:: NSA,NR,NP,NTH,IERR,its
      double precision:: DELE
      double precision,parameter:: FACT=1.01D0, DELE2=1.D-5
      double precision,dimension(NRMAX):: RJ_PO, RJ_PP
      double precision, &
         dimension(NTHMAX,NPSTART:NPENDWG,NRSTART:NRENDWM,NSAMAX)::FEPPS
      double precision, &
         dimension(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX)::FETHS

      IF(N_IMPL.ne.0)THEN
         SIGM(:)=SIGP(:)
         CALL FPCURRENT(RJ_P)
         CALL j_to_i(RJ_P,RI_P)
      END IF

      FNS0(:,:,:,:) =FNSM(:,:,:,:)
      FEPPS(:,:,:,:)=FEPP(:,:,:,:)
      FETHS(:,:,:,:)=FETH(:,:,:,:)

      IF(N_IMPL.eq.0)THEN
         CALL FPCURRENT(RJ_M)
         CALL j_to_i(RJ_M,RI_M)
      END IF

      DO NSA=NSASTART,NSAEND
         CALL fp_exec_sigma(NSA,IERR,its) ! partial j/ partial t
      END DO
      CALL FPCURRENT(RJ_PP)

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  FEPP(NTH,NP,NR,NSA)= AEFP(NSA)*(EP(NR)+DELE2)/PTFP0(NSA)*COSM(NTH)
                  FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
               END DO
            END DO
         END DO
      END DO
      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  FETH(NTH,NP,NR,NSA)=-AEFP(NSA)*(EP(NR)+DELE2)/PTFP0(NSA)*SING(NTH) 
                  FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA) 
               END DO
            END DO
         END DO
      END DO

      DO NSA=NSASTART,NSAEND
         CALL fp_exec_sigma(NSA,IERR,its) ! dj/ dt
      END DO
      CALL FPCURRENT(RJ_PO)

      DO NR=1,NRMAX
         DELE=(FACT-1.D0)*EP(NR)
!         SIGP(NR)=(RJ_PO(NR) - RJ_PP(NR) )/DELE
         SIGP(NR)=(RJ_PO(NR) - RJ_PP(NR) )/DELE2*10
      END DO

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  FEPP(NTH,NP,NR,NSA)= FEPPS(NTH,NP,NR,NSA)
                  FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
               END DO
            END DO
         END DO
      END DO
      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  FETH(NTH,NP,NR,NSA)=FETHS(NTH,NP,NR,NSA)
                  FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA) 
               END DO
            END DO
         END DO
      END DO

      IF(N_IMPL.eq.0) SIGM(:)=SIGP(:)
      IF(NRANK.eq.0)THEN
         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      SIGMA= ",(SIGP(NR),NR=1,8)
         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      RJ_M = ",(RJ_M(NR),NR=1,8)
         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      RJ_P = ",(RJ_P(NR),NR=1,8)
         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      RJ_PO= ",(RJ_PO(NR),NR=1,8)
         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      RJ_PP= ",(RJ_PP(NR),NR=1,8)
      END IF
      END SUBROUTINE conductivity_sigma_ind
!-------------------------------------------------
      SUBROUTINE fp_exec_sigma(NSA,IERR,its)

      USE libmtx
      USE fpexec
      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NSBA
      integer:: NTHS, NLL
      integer:: IERR,its,i,j,ll1
      integer:: imtxstart1,imtxend1

      NSBA=NSB_NSA(NSA)

      CALL mtx_set_communicator(comm_nrnp) !3D

!     ----- Set up matrix solver -----
      CALL mtx_setup(imtxsize,imtxstart1,imtxend1,imtxwidth)
      IF(imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend) THEN
         WRITE(6,*) 'XX fp_exec: '
         WRITE(6,*) '   imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend'
         WRITE(6,*) '   imtxstart1,imtxstart = ',imtxstart1,imtxstart
         WRITE(6,*) '   imtxend1,imtxend     = ',imtxend1,imtxend
         STOP
      ENDIF

!     ----- Set up weight array -----

      CALL FPWEIGHT(NSA,IERR)

!     ----- Set up index array NMA -----
!               NM: line number of the coefficient matrix
!               NL: 

      CALL SET_FM_NMA(NSA,FNSM)

      DO NM=NMSTART,NMEND
         NLMAX(NM)=0
         BM(NM)=0.D0
         DO NL=1,NLMAXM
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
         ENDDO
      ENDDO

!     ----- Calculate matrix coefficients in a row -----

      DO NR=NRSTART,NREND
         IF(MODELA.EQ.0) THEN
            DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               CALL FPSETM(NTH,NP,NR,NSA,NLMAX(NM))
            ENDDO
            ENDDO
         ELSE
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX/2
                  NM=NMA(NTH,NP,NR)
                  CALL FPSETM(NTH,NP,NR,NSA,NLMAX(NM))
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  NM=NMA(NTH,NP,NR)
                  CALL FPSETM(NTH,NP,NR,NSA,NLMAX(NM))
               ENDDO
            ENDDO
         ENDIF
      ENDDO

!     ----- Diagonal term -----

      DO NR=NRSTART,NREND ! LHS
         IF(MODELA.EQ.0) THEN
            DO NP=NPSTART,NPEND
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               BM(NM)=(1.D0+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                     +DELT*SPP(NTH,NP,NR,NSA)
               IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
                  CALL mtx_set_matrix(nm,nm,1.D0-rimpl*delt*dl(nm))
                  CALL mtx_set_vector(nm,FM(NM))
               ENDIF
            ENDDO
            ENDDO
         ELSE
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX/2
                  NM=NMA(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                        +DELT*SPP(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)
                  IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
                     CALL mtx_set_matrix(nm,nm, &
                                         RLAMDA(NTH,NR)-RIMPL*DELT*DL(NM))
                     CALL mtx_set_vector(nm,FM(NM))
                  ENDIF
               ENDDO
               DO NTH=NTHMAX/2+1,ITU(NR)
                  NTHS=NTHMAX+1-NTH
                  NM=NMA(NTH,NP,NR)
                  BM(NM)=0.D0
                  IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
                     CALL mtx_set_matrix(nm,nm,1.d0)
                     CALL mtx_set_matrix(nm,nm+NTHS-NTH,-1.d0)
                     CALL mtx_set_vector(nm,FM(NM))
                  ENDIF
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  NM=NMA(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                        +DELT*SPP(NTH,NP,NR,NSA)*RLAMDA(NTH,NR)
                  IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
                     CALL mtx_set_matrix(nm,nm, &
                                         RLAMDA(NTH,NR)-RIMPL*DELT*DL(NM))
                     CALL mtx_set_vector(nm,FM(NM))
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO

!     ----- Off diagonal term -----

      DO NM=NMSTART,NMEND ! LHS
         IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
            DO NL=1,NLMAX(NM)
               IF(LL(NM,NL).NE.0) THEN
                  CALL mtx_set_matrix(nm,LL(NM,NL),-RIMPL*DELT*AL(NM,NL))
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DO NR=NRSTART,NREND
      DO NP=NPSTART,NPEND
      DO NTH=NTHMAX/2+1,ITU(NR)
         NM=NMA(NTH,NP,NR)
         IF(LL(NM,NL).NE.0) WRITE(6,'(A,5I5,1PE12.4)') &
              'NR,NP,NTH,NM.NL,AL=',NR,NP,NTH,NM,NL,AL(NM,NL)
      ENDDO
      ENDDO
      ENDDO

!     ----- Source vector: contribution from off-diagonal term -----

      DO NM=NMSTART,NMEND ! RHS
         DO NL=1,NLMAX(NM)
            IF(LL(NM,NL).NE.0) THEN
               BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM(LL(NM,NL))
            ENDIF
         ENDDO
         IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
            CALL mtx_set_source(nm,BM(NM))
         ENDIF
      ENDDO

!     ----- Solve matrix equation -----

      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) ! ncom is nessesary for MUMPS not PETSc
      ierr=0

!     ----- Get solution vector -----

      CALL mtx_gather_vector(BMTOT)
      
      DO NR=NRSTARTW,NRENDWM
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FNS0(NTH,NP,NR,NSBA)=BMTOT(NM)
            ENDDO
         ENDDO
      ENDDO

!     ----- Clean up matrix solver -----
      CALL mtx_cleanup

      CALL mtx_reset_communicator

      RETURN
      END SUBROUTINE FP_EXEC_SIGMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INDUCTANCE_EDGE(E_EDGE)

      IMPLICIT NONE
      INTEGER:: NR
      double precision:: SUM, SUM2, E_EDGE, L_EDGE, CAP_L

      SUM=0.D0
      SUM2=0.D0
      DO NR=1,NRMAX
         SUM=SUM+RMU0/(2.D0*PI*RM(NR))*RI_P(NR)*DELR
      END DO
      DO NR=1,NRMAX
         SUM2=SUM2+RMU0/(2.D0*PI*RM(NR))*RI_M(NR)*DELR
      END DO
      L_EDGE=(SUM-SUM2)/(RI_P(NRMAX)-RI_M(NRMAX))

!      E_EDGE=-L_EDGE*( RI_P(NRMAX)-RI_M(NRMAX) )/DELT

      CAP_L=L_EDGE*AEE**2*RNS(NRMAX,1)/AMFP(1)

      E_EDGE=-( L_EDGE/DELT*(RJ_P(NRMAX)-RJ_M(NRMAX))+CAP_L*EP(NRMAX) )/(1.D0+CAP_L)

      END SUBROUTINE INDUCTANCE_EDGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE j_to_i(RJ,RI)

      IMPLICIT NONE
      INTEGER:: NR 
      double precision,dimension(NRMAX),INTENT(IN)::RJ
      double precision,dimension(NRMAX),INTENT(OUT)::RI

      RI(:)=0.D0
      DO NR=1,NRMAX
         RI(NR)=RI(NR)+RJ(NR)*VOLR(NR)/(2.D0*PI*RR)
      END DO

      END SUBROUTINE j_to_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FPCURRENT(RJ_T)
!
      USE fpmpi
      IMPLICIT NONE
      integer:: NR, NSA, NSB, NSBA, NP, NTH, NS, NSW, N
      integer:: IERR
      real(8):: RSUM2, FACT, PV
      double precision,dimension(NRSTART:NREND,NSAMAX)::RJ_L
      double precision,dimension(NRMAX,NSAMAX)::RJ_G
      double precision,dimension(NRMAX),INTENT(OUT)::RJ_T

      CALL mtx_set_communicator(comm_np) 
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         NSBA=NSB_NSA(NSA)
         DO NR=NRSTART,NRENDX
            RSUM2=0.D0
            IF(MODELA.eq.0) THEN
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNS0(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                       &
                             +VOLP(NTH,NP,NSBA)*FNS0(NTH,NP,NR,NSBA) &
                             *PM(NP,NSBA)*COSM(NTH)/PV
                     END DO
                  END DO
               ENDIF
            ELSE
               IF(MODELR.EQ.0) THEN
                  DO NP=NPSTART,NPEND
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNS0(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)*RLAMDA(NTH,NR)
                     END DO
                  ENDDO
               ELSE
                  DO NP=NPSTART,NPEND
                     PV=SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
                     DO NTH=1,NTHMAX
                        RSUM2 = RSUM2                        &
                             +VOLP(NTH,NP,NSBA)*FNS0(NTH,NP,NR,NSBA)  &
                             *PM(NP,NSBA)*COSM(NTH)/PV*RLAMDA(NTH,NR)
                     END DO
                  END DO
               ENDIF               
            END IF

            CALL p_theta_integration(RSUM2) 

            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)*RCOEFNG(NR)
            RJ_L(NR,NSA) = RSUM2*FACT*AEFP(NSA)*PTFP0(NSA) &
                           /AMFP(NSA)
         ENDDO ! NR
      ENDDO ! NSA
   
      CALL mtx_set_communicator(comm_nsanr)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL fp_gatherv_real8_sav(RJ_L,SAVLEN(NRANK+1),RJ_G,N,NSA)
      END DO

      RJ_T(:)=0.D0
      DO NSA=1,NSAMAX
         DO NR=1,NRMAX
            RJ_T(NR)=RJ_T(NR)+RJ_G(NR,NSA)
         END DO
      END DO
      CALL mtx_reset_communicator 
      CALL mtx_broadcast_real8(RJ_T,NRMAX)
      
      RETURN
      END SUBROUTINE FPCURRENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

    end MODULE fpcaleind
