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

      USE fpcoef, only: FP_CALE, FP_CALE_IND
      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH

      DO NSA=NSASTART,NSAEND
         CALL FP_CALE_IND(NSA)
      END DO

!     UPDATE FPP, FTH
      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA) &
                 +FEPP_IND(NTH,NP,NR,NSA)
         END DO
         END DO
!
         DO NP=NPSTARTW,NPENDWM 
         DO NTH=1,NTHMAX+1 
            FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA) &
                 +FETH_IND(NTH,NP,NR,NSA)
         END DO
         END DO
         END DO
      END DO

      IF(NPENDWG.eq.NPMAX+1)THEN
      DO NSA=NSASTART,NSAEND
      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            FPP(NTH,NPMAX+1,NR,NSA)=max(0.D0,FPP(NTH,NPMAX+1,NR,NSA))
         END DO
      END DO
      END DO
      END IF

      END SUBROUTINE UPDATE_FEPP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE E_IND_EVOL!(L_EDGE)

      IMPLICIT NONE
!      double precision,intent(in):: L_EDGE!, l_j!,E_EDGE
      double precision:: E_EDGE, L_EDGE
      INTEGER:: NSA,NR,NP,NTH,IERR
      double precision:: coef0,coefm,coefp,RHS,tauIp,time_now
      integer:: imtxstart1,imtxend1,its,I_IMPL,NITE
      double precision:: EPS_EIND
      double precision,dimension(NRMAX):: EPS_EP
      real:: gut1, gut2

      NITE=1
      CALL GUTIME(gut1)
      I_IMPL=0
!      DO WHILE(I_IMPL.le.NITE)
         I_IMPL=I_IMPL+1

      EPM(:)=EP(:)
      CALL mtx_set_communicator(comm_nr) !3D

      CALL mtx_setup(NRMAX,imtxstart1,imtxend1)
!---- DIAGONAL TERM
      DO NR=NRSTART,NREND
         coef0=2.D0/(RA*DELR)**2 + RMU0*SIGP(NR)/DELT 
!         coef0=2.D0/(RA*DELR)**2 
         CALL mtx_set_matrix(NR,NR,coef0)
         CALL mtx_set_vector(NR,EPM(NR))
      END DO
!---- OFF DIAGONAL
      DO NR=NRSTART,NREND
         coefm=-RG(NR)/( RM(NR)*(RA*DELR)**2 )
         coefp=-RG(NR+1)/( RM(NR)*(RA*DELR)**2 )
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
            RHS=RMU0/DELT*( RJ_M(NR)-RJ_P(NR)+SIGP(NR)*EP(NR) )
         ELSEIF(NR.eq.NRMAX)THEN ! boundary condition E_ind at RM(NRMAX)+DELR
!            E_EDGE=0.D0 ! 
            CALL EDGE_FIELD(E_EDGE,L_EDGE)
            RHS=RMU0/DELT*( RJ_M(NR)-RJ_P(NR)+SIGP(NR)*EP(NR) ) &
                 + RG(NR+1)/( RM(NR)*(RA*DELR)**2 )*E_EDGE
         END IF
         CALL mtx_set_source(NR,RHS)
      END DO

!---- SOLVE
      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) 
!      IF(NRANK.eq.0) write(6,*) 'E_IND_EVOL, Number of iterations    =',its
      IF(NRSTART.eq.NRMAX.and.NPSTART.eq.1.and.NSASTART.eq.1.and.I_IMPL.eq.1)  &
           WRITE(*,'(A,4E16.8)') "L_EDGE,E_EDGE,RI_P,RI_M=", L_EDGE, E_EDGE, RI_P(NRMAX), RI_M(NRMAX)
      CALL mtx_gather_vector(EP)
      CALL mtx_cleanup
      CALL mtx_reset_communicator

!         EPS_EIND=0.D0
!         DO NR=1,NRMAX
!            EPS_EIND = EPS_EIND + (EP(NR)-EPM(NR))**2/(EP(NR))**2
!            EPS_EP(NR) = (EP(NR)-EPM(NR))**2/(EP(NR))**2
!         END DO
!         EPS_EIND=SQRT(EPS_EIND)
         IF(NRSTART.eq.NRMAX.and.NPSTART.eq.1.and.NSASTART.eq.1.and.MOD(I_IMPL,1).eq.0)THEN
            WRITE(*,'(A,I5,2E14.6)') "I_IMPL,E_EDGE",I_IMPL, E_EDGE, SIGP_E
!            WRITE(*,'(16E10.2)') (EPS_EP(NR),NR=1,NRMAX)
!            WRITE(*,'(16E10.2)') (EP(NR),NR=1,NRMAX)
         END IF
!      END DO

      CALL GUTIME(gut2)
!      IF(NRANK.eq.0) WRITE(*,*) "GUT_EVOL_ITERATION = ",gut2-gut1

      END SUBROUTINE E_IND_EVOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE EDGE_FIELD(E_EDGE,L_EDGE)

      IMPLICIT NONE
      DOUBLE PRECISION,intent(out):: E_EDGE, L_EDGE
      DOUBLE PRECISION:: l_j
      INTEGER:: ISW_L

      ISW_L=3

      IF(ISW_L.eq.0)THEN
         L_EDGE=RMU0*RI_P(NRMAX)/(4*PI**2*RA**2*RG(NRMAX+1)**2*RJ_P(NRMAX))
         E_EDGE=-L_EDGE/(2*PI*RR*(1.D0+EPSRM2(NRMAX)))*(RI_P(NRMAX)-RI_M(NRMAX))/DELT
      ELSEIF(ISW_L.eq.1)THEN
         l_j=RMU0*RI_P(NRMAX)/(2*PI*RA*RG(NRMAX+1) )*(RJ_P(NRMAX)-RJ_P(NRMAX-1))/DELR
         L_EDGE=l_j
         E_EDGE=-L_EDGE/(2*PI*RR*(1.D0+EPSRM2(NRMAX)))*(RJ_P(NRMAX)-RJ_M(NRMAX))/DELT
      ELSEIF(ISW_L.eq.2)THEN
         l_j=RMU0*RI_P(NRMAX)/(2*PI*RA*RG(NRMAX+1) )*(RJ_P(NRMAX)-RJ_P(NRMAX-1))/DELR
         L_EDGE=l_j
         E_EDGE=L_EDGE/(2*PI*RR*(1.D0+EPSRM2(NRMAX))*DELT+L_EDGE*SIGP(NRMAX) ) &
              *(RJ_M(NRMAX)-RJ_P(NRMAX)+SIGP(NRMAX)*EP(NRMAX) )
      ELSE
         l_j=RMU0*RI_P(NRMAX)/(2*PI*RA*RG(NRMAX+1) )*(0.D0-RJ_P(NRMAX) )/DELR
         L_EDGE=l_j
         E_EDGE=L_EDGE/(2*PI*RR*(1.D0+EPSRM2(NRMAX)+DELR)*DELT+L_EDGE*SIGP_E ) &
              *(RJ_M(NRMAX)-RJ_P(NRMAX)+SIGP_E*E_EDGEM )
         E_EDGEM=E_EDGE
      END IF

      END SUBROUTINE EDGE_FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE conductivity_sigma_ind

      USE fpcoef, only: FP_CALE_LAV
      IMPLICIT NONE
      INTEGER:: NSA,NR,NP,NTH,IERR,its, isw_test
      double precision:: DELE, DELT2
      double precision,parameter:: FACT=1.01D0, DELE2=1.D-5
      double precision,dimension(NRMAX):: RJ_PO, RJ_PP
      double precision, &
         dimension(NTHMAX,NPSTART:NPENDWG,NRSTART:NRENDWM,NSAMAX)::FEPPS
      double precision, &
         dimension(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX)::FETHS
      double precision:: fact2, taue_col, sigma_sp, E_dri

      IF(N_IMPL.eq.0)THEN
         CALL FPCURRENT(RJ_M)
         CALL j_to_i(RJ_M,RI_M)
      END IF

      ISW_TEST=1  ! 0: pd, 1: sp, 2: J/E

      IF(ISW_TEST.eq.0)THEN
         FNS0(:,:,:,:) =FNSM(:,:,:,:)
         FEPPS(:,:,:,:)=FEPP(:,:,:,:)
         FETHS(:,:,:,:)=FETH(:,:,:,:)

         DELT2=DELT
         DELT=1.D-4
         DO NSA=NSASTART,NSAEND
            CALL fp_exec_sigma(NSA,IERR,its) ! partial j/ partial t
         END DO
         CALL FPCURRENT(RJ_PP)
         
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPENDWG
                  DO NTH=1,NTHMAX
                     FEPP(NTH,NP,NR,NSA)= AEFP(NSA)*(E1(NR)+DELE2)/PTFP0(NSA)*COSM(NTH)
                  END DO
               END DO
            END DO
         END DO
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX+1
                     FETH(NTH,NP,NR,NSA)=-AEFP(NSA)*(E1(NR)+DELE2)/PTFP0(NSA)*SING(NTH) 
                  END DO
               END DO
            END DO
         END DO
         IF(MODELA.eq.1)THEN
            DO NSA=NSASTART,NSAEND
               DO NR=NRSTART,NREND
                  CALL FP_CALE_LAV(NR,NSA)
               END DO
            END DO
         END IF
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPENDWG
                  DO NTH=1,NTHMAX
                     FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA) &
                          +FEPP_IND(NTH,NP,NR,NSA)
                  END DO
               END DO
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX+1
                     FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA) &
                          +FETH_IND(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
         END DO
         
         DO NSA=NSASTART,NSAEND
            CALL fp_exec_sigma(NSA,IERR,its) ! dj/ dt
         END DO
         CALL FPCURRENT(RJ_PO)
         DELT=DELT2
         
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPENDWG
                  DO NTH=1,NTHMAX
                     FEPP(NTH,NP,NR,NSA)= FEPPS(NTH,NP,NR,NSA)
                     FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA) &
                          +FEPP_IND(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
         END DO
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTART,NREND
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX+1
                     FETH(NTH,NP,NR,NSA)=FETHS(NTH,NP,NR,NSA)
                     FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA) &
                          +FETH_IND(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
         END DO
         DO NR=1,NRMAX
            SIGP(NR)=(RJ_PO(NR) - RJ_PP(NR) )/DELE2
         END DO
      ELSEIF(ISW_TEST.eq.1)THEN
         IF(NTG1.ge.2)THEN
            DO NR=NRSTART,NREND
               FACT2=AEFP(1)**2*AEFD(2)**2*LNLAM(NR,2,1)/(4.D0*PI*EPS0**2)
               taue_col=3.D0*SQRT(0.5D0*PI)/FACT2*SQRT(AMFP(1)*(AEE*RT_IMPL(NR,1)*1.D3)**3)/RN_IMPL(NR,2)*1.D-20! wesson P. 69
               E_dri=SQRT(RT_IMPL(1,1)*1.D3*AEE*AMFP(1))/(AEFP(1)*taue_col)
               sigma_sp=1.96D0*RNS(NR,1)*1.D20*AEFP(1)**2*taue_col/AMFP(1) ! P. 174
            END DO
            CALL mtx_set_communicator(comm_nr)
            CALL mtx_allgather1_real8(sigma_sp,SIGP)
            CALL mtx_reset_communicator

            taue_col=3.D0*SQRT(0.5D0*PI)/FACT2*SQRT(AMFP(1)*(AEE*RT_E*1.D3)**3)/RN_E*1.D-20
            SIGP_E=1.96D0*RN_E*1.D20*AEFP(1)**2*taue_col/AMFP(1)
         ELSE
            DO NR=NRSTART,NREND
               SIGP(NR)=0.D0
            END DO
         END IF
      ELSEIF(ISW_TEST.eq.2)THEN
         DO NR=1,NRMAX
            IF(N_IMPL.eq.0)THEN
               SIGP(NR)=(RJ_M(NR))/EP(NR)
            ELSE
               SIGP(NR)=(RJ_P(NR))/EP(NR)
            END IF
         END DO
      END IF

      IF(N_IMPL.eq.0) SIGM(:)=SIGP(:)
      IF(NRANK.eq.0)THEN
         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      SIGMA= ",(SIGP(NR),NR=NRSTART,NREND)
!         IF(NTG1.ge.2) WRITE(*,'(A,E16.8)') "E_dri=", E_dri
!         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      RJ_M = ",(RJ_M(NR),NR=1,8)
!         WRITE(6,'(I3,A,8E14.6)') N_IMPL,"      RJ_P = ",(RJ_P(NR),NR=1,8)
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
      SUBROUTINE INDUCTIVE_FIELD

      IMPLICIT NONE
      double precision:: E_IND, L_IND, l_j, SUM
      double precision,dimension(NRMAX):: SUM2
      integer:: NR

!      DO NR=NRSTART,NREND
!         CALL POLOIDAL_FLUX(NR,L_IND)
!      END DO
!      E_EDGEM2=EP(NRMAX)
      CALL E_IND_EVOL!(L_IND)

!      SUM=0.D0
!      DO NR=1,NRMAX
!         SUM=SUM+SQRT( (EP(NR)-EPM(NR) )**2/(EPM(NR))**2 )
!         SUM2(NR)=SQRT( (EP(NR)-EPM(NR) )**2/(EPM(NR))**2 )
!      END DO
!      IF(NPSTART.eq.1.and.NRSTART.eq.1.and.NSASTART.eq.1) &
!           WRITE(*,'(A,32E10.2)') "EPS_EP=", (SUM2(NR),NR=1,NRMAX)
      CALL mtx_broadcast_real8(EP,NRMAX)

      CALL mtx_set_communicator(comm_nr)
!      CALL mtx_allgather1_real8(E_IND,EP)
      CALL mtx_allgather_real8(PSIPM_P,NREND-NRSTART+1,PSIPG_P)
      CALL mtx_allgather_real8(PSIPM_M,NREND-NRSTART+1,PSIPG_M)
      CALL mtx_reset_communicator

      END SUBROUTINE INDUCTIVE_FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INDUCTIVE_FIELD_A1

      IMPLICIT NONE
!      INTEGER,intent(in):: NR
      double precision:: E_IND, DELH, SUM, EM, ETAL, E_NEW
      INTEGER:: NG, NTH, NR

      DO NR=NRSTART,NREND

      DO NG=1,NAVMAX
         DO NTH=1,NTHMAX
            EM_PHIM(NG,NTH,NR)=EP_PHIM(NG,NTH,NR)
         END DO
         DO NTH=1,NTHMAX+1
            EM_PHIG(NG,NTH,NR)=EP_PHIG(NG,NTH,NR)
         END DO
      END DO
      EPM(NR)=EP(NR)

!      DO NTH=1,NTHMAX
      DO NTH=1,ITL(NR)-1
         DELH=2.D0*ETAM(NTH,NR)/NAVMAX
         SUM=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            EM=EM_PHIM(NG,NTH,NR)
            CALL POLOIDAL_FLUX_PHI(NR,NTH,NG,ETAL,EM,E_IND)
            EP_PHIM(NG,NTH,NR)=E_IND
            SUM=SUM+E_IND*DELH
         END DO
         ETHM(NTH,NR)=SUM * RCOEFNG(NR)
      END DO
      DO NTH=ITL(NR)+1,ITU(NR)-1
         ETHM(NTH,NR)=0.D0
      END DO
      DO NTH=ITU(NR)+1,NTHMAX
         DELH=2.D0*ETAM(NTH,NR)/NAVMAX
         SUM=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            EM=EM_PHIM(NG,NTH,NR)
            CALL POLOIDAL_FLUX_PHI(NR,NTH,NG,ETAL,EM,E_IND)
            EP_PHIM(NG,NTH,NR)=E_IND
            SUM=SUM+E_IND*DELH
         END DO
         ETHM(NTH,NR)=SUM * RCOEFNG(NR)
      END DO
      ETHM(ITL(NR),NR) = RLAMDA(ITL(NR),NR)/4.D0      &
           *( ETHM(ITL(NR)-1,NR)/RLAMDA(ITL(NR)-1,NR) &
             +ETHM(ITU(NR)+1,NR)/RLAMDA(ITU(NR)+1,NR) )
      ETHM(ITU(NR),NR) = ETHM(ITL(NR),NR)

!      DO NTH=1,NTHMAX+1
      DO NTH=1,ITL(NR)
         DELH=2.D0*ETAG(NTH,NR)/NAVMAX
         SUM=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            EM=EM_PHIG(NG,NTH,NR)
            CALL POLOIDAL_FLUX_PHI(NR,NTH,NG,ETAL,EM,E_IND)
            EP_PHIG(NG,NTH,NR)=E_IND
            SUM=SUM+E_IND*DELH
         END DO
         ETHG(NTH,NR)=SUM * RCOEFNG(NR)
      END DO
      DO NTH=ITL(NR)+1,ITU(NR)
         ETHG(NTH,NR) = 0.D0
      END DO
      DO NTH=ITU(NR)+1,NTHMAX
         DELH=2.D0*ETAG(NTH,NR)/NAVMAX
         SUM=0.D0
         DO NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            EM=EM_PHIG(NG,NTH,NR)
            CALL POLOIDAL_FLUX_PHI(NR,NTH,NG,ETAL,EM,E_IND)
            EP_PHIG(NG,NTH,NR)=E_IND
            SUM=SUM+E_IND*DELH
         END DO
         ETHG(NTH,NR)=SUM * RCOEFNG(NR)
      END DO

      E_NEW=EP_PHIG(1,1,NR) ! it means inductive field on equator
!     unfortunately EP_PHIG(1,NTHMAX/2+1,NR) has no value
      CALL mtx_set_communicator(comm_nr)
      CALL mtx_allgather1_real8(E_NEW,EP)
      CALL mtx_reset_communicator

      END DO

      END SUBROUTINE INDUCTIVE_FIELD_A1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE POLOIDAL_FLUX(NR,L_IND)

      IMPLICIT NONE
      INTEGER,intent(in):: NR
      double precision:: E_IND ! inductive field
      double precision,intent(out):: L_IND ! inductance of whole plasma
      double precision:: l_j, metric
      double precision:: SUM, SUM2, PSIP, PSIM, rv, delrho, EPS_IND,EPM2
      INTEGER:: NR2, NSWE, N_IMPL2

      NSWE=0

      IF(NSWE.eq.0)THEN ! cylindrical
         SUM=0.D0
         SUM2=0.D0
!         DO NR2=1,NR
!            SUM=SUM+RMU0/(2.D0*PI*RM(NR2))*RI_P(NR2)*DELR
!         END DO
!         DO NR2=1,NR
!            SUM2=SUM2+RMU0/(2.D0*PI*RM(NR2))*RI_M(NR2)*DELR
!         END DO
!         PSIPM_P(NR)=SUM
!         PSIPM_M(NR)=SUM2
!         L_IND=(SUM-SUM2)/(RI_P(NR)-RI_M(NR)) ! inductance 0<RM<RM(NR)
         L_IND=RMU0*RI_P(NRMAX)/(4*PI**2*RA**2*RG(NRMAX+1)**2*RJ_P(NRMAX))
         metric=1.D0
      ELSEIF(NSWE.eq.1)THEN ! torus
         SUM=0.D0
         SUM2=0.D0
         DO NR2=1,NRMAX
            rv = EPSRM2(NR2)*RR
!            DELrho = (EPSRM2(NR2+1)-EPSRM2(NR2) ) * RR
            DELrho = (EPSRG2(NR2+1)-EPSRG2(NR2) ) * RR
            SUM = SUM + RMU0*RI_P(NR2)*(RR/rv - 1.D0) * DELrho
            SUM2 =SUM2+ RMU0*RI_M(NR2)*(RR/rv - 1.D0) * DELrho
         END DO
         PSIP = SUM + &
              RMU0*RI_P(NRMAX)*( RR*LOG(RR/RA) - (RR-RA) )
         PSIM = SUM2 + &
              RMU0*RI_M(NRMAX)*( RR*LOG(RR/RA) - (RR-RA) )

         SUM=0.D0
         SUM2=0.D0
         DO NR2=1,NR
            rv = EPSRM2(NR2)*RR
!            DELrho = (EPSRM2(NR2+1)-EPSRM2(NR2) ) * RR
            DELrho = (EPSRG2(NR2+1)-EPSRG2(NR2) ) * RR
            SUM = SUM + RMU0*RI_P(NR2)*(RR/rv + 1.D0) * DELrho
            SUM2= SUM2+ RMU0*RI_M(NR2)*(RR/rv + 1.D0) * DELrho
         END DO
         PSIP = PSIP - SUM
         PSIM = PSIM - SUM2
         PSIPM_P(NR)=PSIP
         PSIPM_M(NR)=PSIM
         L_IND=(PSIP-PSIM)/(RI_P(NR)-RI_M(NR)) ! inductance 0<RM<RM(NR)
!         metric=2.D0*PI*(RR+RA*RM(NR))
         metric=2.D0*PI*RR*(1.D0+EPSRM2(NR))
      END IF

!      l_j = L_IND*( RI_P(NR)-RI_M(NR) )/( RJ_P(NR)-RJ_M(NR) )

!      E_IND = l_j/(l_j*SIGP(NR)+DELT*metric) &
!           *(RJ_M(NR)-RJ_P(NR)+SIGP(NR)*EPM(NR)) ! implicit 

!      IF(NSASTART.eq.1.and.NPSTART.eq.1.and.NRSTART.eq.1) THEN
!         WRITE(*,'(I3,E16.8,A,3E16.8)') N_IMPL, SQRT(EPS_IND)," l_j, E, Vloop= ",l_j, E_IND, E_IND*2*PI*RR
!      END IF

      END SUBROUTINE POLOIDAL_FLUX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE POLOIDAL_FLUX_PHI(NR,NTH,NG,ETAL,EM,E_IND)

      IMPLICIT NONE
      INTEGER,intent(in):: NR, NTH, NG
      double precision,intent(in):: ETAL,EM
      double precision,intent(out):: E_IND 
      double precision:: L_IND, l_j ! inductance
      double precision:: SUM, SUM2, PSIP, PSIM, rv, delrho
      double precision:: xmax, zz, DELxx, xx, PIP, MIP,metric
      INTEGER:: NR2, nx
      INTEGER,parameter:: NXMAX=100

      IF(ETAL.le.0.5D0*PI)THEN
         xmax=0.D0
      ELSE
         xmax=EPSRM2(NR)*RR*ABS(COS(ETAL))
      END IF
      zz = EPSRM2(NR)*RR*SIN(ETAL)
      DELxx=(RR-xmax)/NXMAX

      SUM=0.D0
      SUM2=0.D0
      DO nx=1,nxmax
         xx = xmax + DELxx * (nx-0.5D0)
         rv = SQRT(xx**2 + zz*2)
         CALL INTERPOLATION(rv,RI_P,PIP)
         SUM = SUM + RMU0*PIP*xx/(rv**2)*(RR-xx)*DELxx
         CALL INTERPOLATION(rv,RI_M,MIP)
         SUM2= SUM2+ RMU0*MIP*xx/(rv**2)*(RR-xx)*DELxx
      END DO
      PSIP = SUM
      PSIM = SUM2

      SUM=0.D0
      SUM2=0.D0
      IF(ETAL.le.0.5D0*PI)THEN
         xmax=EPSRM2(NR)*RR*ABS(COS(ETAL))
         DELxx=(RR-xmax)/NXMAX
         DO nx=1,nxmax
            xx = xmax + DELxx * (nx-0.5D0)
            rv = SQRT(xx**2 + zz*2)
            CALL INTERPOLATION(rv,RI_P,PIP)
            SUM = SUM + RMU0*PIP*xx/(rv**2)*(RR+xx)*DELxx
            CALL INTERPOLATION(rv,RI_M,MIP)
            SUM2= SUM2+ RMU0*MIP*xx/(rv**2)*(RR+xx)*DELxx
         END DO
      ELSE
         SUM=0.D0
         SUM2=0.D0
      END IF
      
      PSIP = PSIP - SUM
      PSIM = PSIM - SUM2
      L_IND=(PSIP-PSIM)/(RI_P(NR)-RI_M(NR)) ! inductance 0<RM<RM(NR)

      metric=2.D0*PI*RR*(1.D0+EPSRM2(NR)*COS(ETAL))
      l_j = L_IND*( RI_P(NR)-RI_M(NR) )/( RJ_P(NR)-RJ_M(NR) )

      E_IND = l_j/(l_j*SIGP(NR)+DELT*metric) &
           *(RJ_M(NR)-RJ_P(NR)+SIGP(NR)*EM) ! implicit 


      END SUBROUTINE POLOIDAL_FLUX_PHI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE INTERPOLATION(rv,IP,XIP)

      IMPLICIT NONE
      double precision,intent(in):: rv
      double precision,dimension(NRMAX),intent(in):: IP
      double precision,intent(out):: XIP
      INTEGER:: NR, NRL, NRU
      double precision:: rvmax, rv1, rv2

      NRL=0
      NRU=NRMAX
      rvmax = EPSRM2(NRMAX)*RR
      IF(rv.ge.rvmax)THEN
         XIP=IP(NRMAX)
      ELSE
         DO NR=1,NRMAX
            rvmax = EPSRM2(NR)*RR
            IF(rv.ge.rvmax)THEN
               NRL=NR
            END IF
         END DO
         DO NR=NRMAX,1,-1
            rvmax = EPSRM2(NR)*RR
            IF(rv.le.rvmax)THEN
               NRU=NR
            END IF
         END DO
         IF(NRL.eq.0)THEN
            rv1=0.D0
            rv2=EPSRM2(NRU)*RR
            XIP=( IP(NRU) )/(rv2-rv1)*(rv-rv1)
         ELSE
            rv1=EPSRM2(NRL)*RR
            rv2=EPSRM2(NRU)*RR
            XIP=IP(NRL)+( IP(NRU)-IP(NRL) )/(rv2-rv1)*(rv-rv1)            
         END IF
      END IF

      END SUBROUTINE INTERPOLATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE j_to_i(RJ,RI)

      IMPLICIT NONE
      INTEGER:: NR 
      double precision,dimension(NRMAX),INTENT(IN)::RJ
      double precision,dimension(NRMAX),INTENT(OUT)::RI

      RI(:)=0.D0
      RI(1)=RJ(1)*VOLR(1)/(2.D0*PI*RR)
      DO NR=2,NRMAX
         RI(NR)=RI(NR-1)+RJ(NR)*VOLR(NR)/(2.D0*PI*RR)
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

            FACT=RNFP0(NSA)*1.D20/RFSADG(NR)!*RCOEFNG(NR)
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
