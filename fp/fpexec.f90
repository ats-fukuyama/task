!     $Id$

!     **************************
!        EXECUTE TIME ADVANCE
!     **************************

      MODULE fpexec

      use fpcomm

      contains

      SUBROUTINE fp_exec(NSA,IERR,its)

      USE libmpi
      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NSBA, NN
      integer:: NTHS, NLL
      integer:: IERR,its,i,j,ll1
      integer:: imtxstart1,imtxend1
!      integer,optional:: methodKSP, methodPC
      real(8),dimension(nmend-nmstart+1):: BM_L
      real(8),dimension(nthmax):: sendbuf_p, recvbuf_p
      real(8),dimension(nthmax*(npend-npstart+1)):: sendbuf_r, recvbuf_r

      NSBA=NSB_NSA(NSA)

!      CALL mtx_set_communicator(comm_nr) !2D
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

      DO NR=NRSTART,NREND ! LHS: set_matrix, RHS: BM
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

!      DO NR=NRSTART,NREND
!      DO NP=NPSTART,NPEND
!      DO NTH=NTHMAX/2+1,ITU(NR)
!         NM=NMA(NTH,NP,NR)
!         IF(LL(NM,NL).NE.0) WRITE(6,'(A,5I5,1PE12.4)') &
!              'NR,NP,NTH,NM.NL,AL=',NR,NP,NTH,NM,NL,AL(NM,NL)
!      ENDDO
!      ENDDO
!      ENDDO

!     ----- Source vector: contribution from off-diagonal term -----

      DO NM=NMSTART,NMEND ! RHS
         DO NL=1,NLMAX(NM)
            NN=LL(NM,NL)
            IF(NN.NE.0) THEN
               IF(NN.ge.NMSTART-NTHMAX.and.NN.le.NMEND+NTHMAX)THEN
                  BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM(NN)
               ELSEIF(NN.lt.NMSTART-NTHMAX)THEN
                  BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM_shadow_m(NN)
               ELSE
                  BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM_shadow_p(NN)
               END IF
            ENDIF
         ENDDO
         IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
            CALL mtx_set_source(nm,BM(NM))
         ENDIF
      ENDDO

!      DO NM=NMSTART,NMEND ! RHS
!         DO NL=1,NLMAX(NM)
!            IF(LL(NM,NL).NE.0) THEN
!               BM(NM)=BM(NM)+(1.D0-RIMPL)*DELT*AL(NM,NL)*FM(LL(NM,NL))
!            ENDIF
!         ENDDO
!         IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
!            CALL mtx_set_source(nm,BM(NM))
!         ENDIF
!      ENDDO

!     ----- Solve matrix equation -----

      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) ! ncom is nessesary for MUMPS not PETSc
!      IF(MODELD_temp.eq.0)THEN
!         if(nrank.eq.0) then
!            write(6,*) 'Number of iterations, NSA    =',its,NSA
!         endif
!      END IF
      ierr=0

!     ----- Get solution vector -----

      CALL mtx_vector(BM_L)
      DO NR=NRSTART, NREND
         DO NP=NPSTART, NPEND
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FNS0(NTH,NP,NR,NSBA)=BM_L(NM-NMSTART+1)
            END DO
         END DO
      END DO
!     shadow requires to communicate
      CALL mtx_set_communicator(comm_np)
      DO NR=NRSTART, NREND
         CALL shadow_comm_np(NR,NSBA)
      END DO
      CALL mtx_set_communicator(comm_nr)
      CALL shadow_comm_nr(NSBA)
      CALL mtx_set_communicator(comm_nrnp) !3D


!      CALL mtx_gather_vector(BMTOT)
!      DO NR=NRSTARTW,NRENDWM
!         DO NP=NPSTARTW,NPENDWM
!            DO NTH=1,NTHMAX
!               NM=NMA(NTH,NP,NR)
!               FNS0(NTH,NP,NR,NSBA)=BMTOT(NM)
!            ENDDO
!         ENDDO
!      ENDDO

!     ----- Clean up matrix solver -----
      CALL mtx_cleanup

      CALL mtx_reset_communicator

      RETURN
      END SUBROUTINE FP_EXEC

!     ---------------------------------

      SUBROUTINE SET_FM_NMA(NSA,func_in)

      IMPLICIT NONE
      integer:: NTH, NP, NR, NSA, NSBA, NM, NRS, NPS
      double precision,dimension(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND), &
           intent(IN):: func_in

      IF(NRSTART.eq.1)THEN
         NRS=1
      ELSE
         NRS=NRSTART-1
      END IF

      DO NR=NRSTARTW,NRENDWM
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NTH+NTHMAX*(NP-1)+NPMAX*NTHMAX*(NR-1)
               NMA(NTH,NP,NR)=NM
            END DO
         END DO
      END DO

      NSBA=NSB_NSA(NSA)
      DO NR=NRSTART,NREND
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FM(NM)=func_in(NTH,NP,NR,NSBA)
            ENDDO
         ENDDO
      ENDDO
      NR=NRSTARTW
      IF(NR.ne.NRSTART)THEN
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FM_shadow_m(NM)=func_in(NTH,NP,NR,NSBA)
            ENDDO
         ENDDO
      END IF
      NR=NRENDWM
      IF(NR.ne.NREND)THEN
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FM_shadow_p(NM)=func_in(NTH,NP,NR,NSBA)
            ENDDO
         ENDDO
      END IF

!      DO NR=NRSTARTW,NRENDWM
!         DO NP=NPSTARTW,NPENDWM
!            DO NTH=1,NTHMAX
!               NM=NMA(NTH,NP,NR)
!               FM(NM)=func_in(NTH,NP,NR,NSA)
!            ENDDO
!         ENDDO
!      ENDDO
      
      END SUBROUTINE SET_FM_NMA

!
!     ***************************
!        CALCULATION OF WEIGHT
!     ***************************
!
      SUBROUTINE FPWEIGHT(NSA,IERR) ! proposed by Chang and Cooper [30] in Karney

      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NTHA, NTHB, NTB, NSBA
      integer:: IERR
      real(8):: DFDTH, FVEL, DVEL, DFDP, DFDB

!     +++++ calculation of weigthing (including off-diagonal terms) +++++

      real(8)::EPSWT=1.D-70

      NSBA=NSB_NSA(NSA)
      DO NR=NRSTART,NREND
!      DO NP=1,NPMAX+1
      DO NP=NPSTART,NPENDWG
      DO NTH=1,NTHMAX
         DFDTH=0.D0
         IF(NP.NE.1) THEN
            NTHA=MIN(NTH+1,NTHMAX)
            NTHB=MAX(NTH-1,1)
            IF(NP.EQ.NPMAX+1) THEN
               IF(ABS(F(NTH,NP-1,NR)).GT.EPSWT) THEN
                  DFDTH= (F(NTHA,NP-1,NR)-F(NTHB,NP-1,NR)) &
                        /(2.D0*PG(NP,NSBA)*DELTH*F(NTH,NP-1,NR))
               ENDIF
            ELSE
               IF(ABS(F(NTH,NP-1,NR)).GT.EPSWT) THEN
                  IF(ABS(F(NTH,NP,NR)).GT.EPSWT) THEN
                     DFDTH= (F(NTHA,NP-1,NR)-F(NTHB,NP-1,NR)) &
                           /(4.D0*PG(NP,NSBA)*DELTH*F(NTH,NP-1,NR)) &
                          + (F(NTHA,NP  ,NR)-F(NTHB,NP  ,NR)) &
                           /(4.D0*PG(NP,NSBA)*DELTH*F(NTH,NP  ,NR))
                  ELSE
                     DFDTH= (F(NTHA,NP-1,NR)-F(NTHB,NP-1,NR)) &
                           /(2.D0*PG(NP,NSBA)*DELTH*F(NTH,NP-1,NR)) 
                  ENDIF
               ELSE
                  IF(ABS(F(NTH,NP,NR)).GT.EPSWT) THEN
                     DFDTH= (F(NTHA,NP  ,NR)-F(NTHB,NP  ,NR)) &
                           /(2.D0*PG(NP,NSBA)*DELTH*F(NTH,NP  ,NR))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         FVEL=FPP(NTH,NP,NR,NSA)-DPT(NTH,NP,NR,NSA)*DFDTH
         WEIGHP(NTH,NP,NR,NSA)=FPWEGH(-DELP(NSBA)*FVEL,DPP(NTH,NP,NR,NSA))
!         FVEL=FCPP(NTH,NP,NR,NSA)-DCPT(NTH,NP,NR,NSA)*DFDTH
!         WEIGHP(NTH,NP,NR,NSA)=FPWEGH(-DELP(NSBA)*FVEL,DCPP(NTH,NP,NR,NSA))
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NREND
!      DO NP=1,NPMAX
      DO NP=NPSTART,NPEND
         DFDP=-PM(NP,NSBA)*RTFP0(NSA)/RTFP(NR,NSA)
         DFDB=DFDP
      DO NTH=1,NTHMAX+1
        IF(NTH.EQ.1) THEN
            IF(F(NTH,NP,NR).GT.EPSWT) THEN
               IF(NP.EQ.1) THEN
                  DFDP= (F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                       /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
               ELSEIF(NP.EQ.NPMAX) THEN
                  DFDP= (F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                       /(     DELP(NSBA)*F(NTH  ,NP,NR))
               ELSE
                  DFDP= (F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                       /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
               ENDIF
            ENDIF
         ELSEIF(NTH.EQ.NTHMAX+1) THEN
            IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
               IF(NP.EQ.1) THEN
                  DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                       /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR))
               ELSEIF(NP.EQ.NPMAX) THEN
                  DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                       /(     DELP(NSBA)*F(NTH-1,NP,NR))
               ELSE
                  DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                       /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR))
               ENDIF
            ENDIF
         ELSE IF(NTH.EQ.ITL(NR)) THEN
            NTB=ITU(NR)
            IF(NP.EQ.1) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSEIF(NP.EQ.NPMAX) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(ABS(F(NTH,NP,NR)).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(     DELP(NSBA)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(     DELP(NSBA)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSE
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ENDIF
            IF(NP.EQ.1) THEN
               IF(F(NTB-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTHMAX-NTB+2,1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTB-1,NP,NR)) &
                          +(F(NTB  ,NP+1,NR)-F(NTHMAX-NTB+1,1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTB  ,NP,NR))
                  ELSE
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTHMAX-NTB+2,1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTB-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB  ,NP+1,NR)-F(NTHMAX-NTB+1,1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTB  ,NP,NR))
                  ENDIF
               ENDIF
            ELSEIF(NP.EQ.NPMAX) THEN
               IF(F(NTB-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB-1,NP,NR)-F(NTB-1,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTB-1,NP,NR)) &
                          +(F(NTB  ,NP,NR)-F(NTB  ,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTB  ,NP,NR))
                  ELSE
                     DFDB= (F(NTB-1,NP,NR)-F(NTB-1,NP-1,NR)) &
                          /(     DELP(NSBA)*F(NTB-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB  ,NP,NR)-F(NTB  ,NP-1,NR)) &
                          /(     DELP(NSBA)*F(NTB  ,NP,NR))
                  ENDIF
               ENDIF
            ELSE
               IF(F(NTB-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTB-1,NP-1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTB-1,NP,NR)) &
                          +(F(NTB  ,NP+1,NR)-F(NTB  ,NP-1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTB  ,NP,NR))
                  ELSE
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTB-1,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTB-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB  ,NP+1,NR)-F(NTB  ,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTB  ,NP,NR))
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF(NP.EQ.1) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSEIF(NP.EQ.NPMAX) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(     DELP(NSBA)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(     DELP(NSBA)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSE
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(4.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NSBA)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF(NTH.EQ.ITL(NR)) THEN
!            FVEL=FTH(NTH,NP,NR,NSA)-DTP(NTH,NP,NR,NSA)*DFDP &
!                -FTH(NTB,NP,NR,NSA)+DTP(NTB,NP,NR,NSA)*DFDB
!            DVEL=DTT(NTH,NP,NR,NSA)+DTT(NTB,NP,NR,NSA)
            FVEL=FCTH(NTH,NP,NR,NSA)-DCTP(NTH,NP,NR,NSA)*DFDP &
                -FCTH(NTB,NP,NR,NSA)+DCTP(NTB,NP,NR,NSA)*DFDB
            DVEL=DCTT(NTH,NP,NR,NSA)+DCTT(NTB,NP,NR,NSA)
         ELSE
!            FVEL=FTH(NTH,NP,NR,NSA)-DTP(NTH,NP,NR,NSA)*DFDP
!            DVEL=DTT(NTH,NP,NR,NSA) 
            FVEL=FCTH(NTH,NP,NR,NSA)-DCTP(NTH,NP,NR,NSA)*DFDP
            DVEL=DCTT(NTH,NP,NR,NSA)
         ENDIF
         WEIGHT(NTH,NP,NR,NSA) &
                 =FPWEGH(-DELTH*PM(NP,NSBA)*FVEL,DVEL)
!         WEIGHT(NTH,NP,NR,NSA) &
!                 =0.5D0
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NRENDWG
      DO NP=NPSTART,NPEND
      DO NTH=1,NTHMAX
         FVEL=FRR(NTH,NP,NR,NSA)
         DVEL=DRR(NTH,NP,NR,NSA)
         WEIGHR(NTH,NP,NR,NSA)=FPWEGH(-DELR*FVEL,DVEL)
!         WEIGHR(NTH,NP,NR,NSA)=0.5D0
!         IF(NR.ne.1)THEN
!            WEIGHR(NTH,NP,NR,NSA)=(4.D0*RG(NR)+DELR)/(8.D0*RG(NR))
!         ELSE
!            WEIGHR(NTH,NP,NR,NSA)=0.5D0
!         ENDIF
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FPWEIGHT

! ************************************************
!     WEIGHTING FUNCTION FOR CONVECTION EFFECT
! ************************************************

      FUNCTION FPWEGH(X,Y)

      IMPLICIT NONE
      real(8):: X, Y, Z
      real(8):: FPWEGH

      IF(ABS(Y).LT.1.D-70) THEN
         IF(X.GT.0.D0) THEN
            FPWEGH=0.D0
         ELSEIF(X.LT.0.D0) THEN
            FPWEGH=1.D0
         ELSE
            FPWEGH=0.5D0
         ENDIF
      ELSE
         Z=X/Y
         IF(ABS(Z).LT.1.D-5)THEN
            FPWEGH=0.5D0-Z/12.D0+Z**3/720.D0
         ELSE IF(Z.GE.100.D0)THEN
            FPWEGH=1.D0/Z
         ELSE IF(Z.LE.-100.D0)THEN
            FPWEGH=1.D0/Z+1.D0
         ELSE
            FPWEGH=1.D0/Z-1.D0/(EXP(Z)-1.D0)
         END IF
      ENDIF
      RETURN
      END FUNCTION FPWEGH

! ******************************************
!     Calculation of matrix coefficients
! ******************************************

      SUBROUTINE FPSETM(NTH,NP,NR,NSA,NL)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTH,NP,NR,NSA
      INTEGER,INTENT(OUT):: NL
      INTEGER:: NM, NTHA, NTHB
                     ! if NTH=ITL(NR)-1 (NTH:untrapped, NTH+1:trapped)
                     !       flux to NTH-1 as usual
                     !       flux to NTH+1 as usual
                     ! if NTH=ITL(NR) (NTH: trapped, NTH-1:untrapped)
                     !       flux to NTH-1 and ITU(NR)+1
                     !       flux to NTH+1 as usual
                     ! if NTH=ITU(NR)+1 (NTH:untrapped, NTH-1:trapped)
                     !       flux to ITL(NR)
                     !       flux to NTH+1 as usual

                     ! --- we assume monotonic trapping boundary 
                     ! if NTH<=ITL(NR)-1 
                     !    if NTH<=ITL(NR+1)-1 (NR:untrapped, NR+1:untrapped)
                     !       flux to NR+1 as usual
                     !    if NTH>=ITL(NR+1)   (NR:untrapped, NR+1:trapped)
                     !       flux to NR+1 as usual
                     ! if NTH>=ITL(NR)
                     !    if NTH<=ITL(NR-1)-1 (NR:trapped, NR-1:untrapped)
                     !       flux(NR-1) to NTH and ITU(NR-1)+ITL(NR-1)-NTH
                     !    if NTH>=ITL(NR-1)   (NR:trapped, NR-1:trapped)
                     !       flux(NR-1) to NTH as usual
                     ! if NTH>=ITU(NR)+1
                     !    if NTH<=ITU(NR+1) (NR:untrapped, NR+1:trapped)
                     !       flux(NR+1) to ITL(NR+1)+ITU(NR+1)-NTH
                     !    if NTH>=ITU(NR+1)+1 (NR:untrapped, NR+1:untrapped)
                     !       flux(NR+1) to NTH as usual

      INTEGER:: NTB  ! if NTH=ITL(NR)
                     !    then NTB=ITU(NR)+1, else NTB=0
      INTEGER:: NTBM ! if NTH>=ITL(NR) and NTH<=ITL(NR-1)-1, 
                     !    then NTBM=ITU(NR-1)+ITL(NR-1)-NTH, else NTBM=0
      INTEGER:: NTBP ! if NTH>=ITU(NR)+1 and NTH<=ITU(NR+1)
                     !    then NTBP=ITL(NR+1)+ITU(NR+1)-NTH, else NTBP=0

      integer:: IERR, NSBA
      real(8):: DFDTH, FVEL, DVEL, DFDP, DFDB
      real(8):: DPPM, DPPP, DPTM, DPTP, SL, DTPM, DTPP, DTTM, DTTP
      real(8):: WPM, WPP, WTM, WTP, VPM, VPP, VTM, VTP
      real(8):: WTB, VTB, WRBM, VRBM, WRBP, VRBP
      real(8):: DIVDPP, DIVDPT, DIVDTP, DIVDTT, DIVFPP, DIVFTH
      real(8):: RL,DRRM,DRRP,WRM,WRP,VRM,VRP,DIVDRR,DIVFRR
      real(8):: PL
      DOUBLE PRECISION:: WPBM, VPBM, WPBP, VPBP

      NSBA=NSB_NSA(NSA)

      NTB=0
      NTBM=0
      NTBP=0
      IF(MODELA.ne.0)THEN
!         IF(NTH.EQ.ITL(NR)+1) THEN
!            NTB=ITU(NR)+1
         IF(NTH.EQ.ITL(NR)) THEN
            NTB=ITU(NR)+1
         ENDIF
         IF(NR-1.GE.1) THEN
!            IF(NTH.GE.ITL(NR)+1.AND.NTH.LE.ITL(NR-1)) THEN 
            IF(NTH.GE.ITL(NR).AND.NTH.LT.ITL(NR-1)) THEN 
!               NTBM=ITU(NR-1)+ITL(NR-1)-NTH+1
!               NTBM=ITU(NR)
               NTBM=NTHMAX+1-NTH
            ENDIF
         ENDIF
      END IF ! MODELA

      NL=0
      NM=NMA(NTH,NP,NR)
      PL=PM(NP,NSBA)
      DPPM=PG(NP,NSBA  )**2
      DPPP=PG(NP+1,NSBA)**2
      DPTM=PG(NP,NSBA  )
      DPTP=PG(NP+1,NSBA)
      SL=SINM(NTH)
      DTPM=SING(NTH  )
      DTPP=SING(NTH+1)
      DTTM=SING(NTH  )/PL
      DTTP=SING(NTH+1)/PL
      IF(NTB.NE.0) THEN
         DPTM=0.5D0*DPTM
         DPTP=0.5D0*DPTP
         DTPM=0.5D0*DTPM
         DTTM=0.5D0*DTTM
      ENDIF
      RL=RM(NR)
      DRRM=RG(NR  )
      DRRP=RG(NR+1)
      IF(NTBM.NE.0) THEN
         DRRM=0.5D0*DRRM
      ENDIF
!     delta
      WPM=WEIGHP(NTH  ,NP  ,NR  ,NSA)
      WPP=WEIGHP(NTH  ,NP+1,NR  ,NSA)
      WTM=WEIGHT(NTH  ,NP  ,NR  ,NSA)
      WTP=WEIGHT(NTH+1,NP  ,NR  ,NSA)
      WRM=WEIGHR(NTH  ,NP  ,NR  ,NSA)
      WRP=WEIGHR(NTH  ,NP  ,NR+1,NSA)
!     epsilon
      VPM=1.D0-WPM
      VPP=1.D0-WPP
      VTM=1.D0-WTM
      VTP=1.D0-WTP
      VRM=1.D0-WRM
      VRP=1.D0-WRP
      IF(NTB.NE.0) THEN
         WTB=WEIGHT( NTB, NP  ,NR  ,NSA)
         VTB=1.D0-WTB
         WPBM=WEIGHP( NTB-1, NP  ,NR  ,NSA)
         VPBM=1.D0-WPBM
         WPBP=WEIGHP( NTB-1, NP+1,NR  ,NSA)
         VPBP=1.D0-WPBP
      ENDIF
      IF(NTBM.NE.0) THEN
         WRBM=WEIGHR(NTBM,NP  ,NR,NSA)
         VRBM=1.D0-WRBM
      ENDIF
      DIVDPP=1.D0/(     PL*PL*DELP(NSBA) *DELP(NSBA))
      DIVDPT=1.D0/(2.D0*PL*PL*DELP(NSBA) *DELTH)
      DIVDTP=1.D0/(2.D0*PL*SL*DELTH*DELP(NSBA))
      DIVDTT=1.D0/(     PL*SL*DELTH*DELTH)
      DIVFPP=1.D0/(     PL*PL*DELP(NSBA))
      DIVFTH=1.D0/(     PL*SL*DELTH)
!      DIVDRR=1.D0/(     RL   *DELR *DELR)
!      DIVFRR=1.D0/(     RL   *DELR)
      DIVDRR=1.D0/(     RL   *DELR *DELR)*RFSADG(NR)
      DIVFRR=1.D0/(     RL   *DELR)*RFSADG(NR)
!      IF(NP.EQ.NPMAX) THEN
!         DIVDTP=2.D0*DIVDTP
!      ENDIF

      IF(MODELD.GT.0) THEN
!         IF(NR.NE.1.AND.NP.NE.NPMAX) THEN
         IF(NR.NE.1) THEN
            NL=NL+1
            LL(NM,NL)=NMA(NTH,NP,NR-1)
            AL(NM,NL)=DRR(NTH  ,NP  ,NR,NSA)    *DIVDRR*DRRM &
!                         *RLAMDAG(NTH,NR-1)/RFSADG(NR-1) &
                     +FRR(NTH  ,NP  ,NR,NSA)*WRM*DIVFRR*DRRM 
!                         *RLAMDAG(NTH,NR-1)/RFSADG(NR-1)
            IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
               LL(NM,NL)=0
               AL(NM,NL)=0.D0
               NL=NL-1
            ENDIF
            IF(NTBM.NE.0) THEN
               NL=NL+1
               LL(NM,NL)=NMA(NTBM,NP,NR-1)
               AL(NM,NL)=DRR(NTBM ,NP  ,NR,NSA)     *DIVDRR*DRRM &
!                            *RLAMDAG(NTBM,NR-1)/RFSADG(NR-1) &
                        +FRR(NTBM ,NP  ,NR,NSA)*WRBM*DIVFRR*DRRM 
!                            *RLAMDAG(NTBM,NR-1)/RFSADG(NR-1)
               IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
                  LL(NM,NL)=0
                  AL(NM,NL)=0.D0
                  NL=NL-1
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF(NP.NE.1.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP-1,NR)
         AL(NM,NL)=+DPT(NTH  ,NP  ,NR,NSA)*WPM*DIVDPT*DPTM &
                   +DTP(NTH  ,NP  ,NR,NSA)*WTM*DIVDTP*DTPM

         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(NP.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP-1,NR)
         AL(NM,NL)=+DPP(NTH  ,NP  ,NR,NSA)    *DIVDPP*DPPM &
                   +FPP(NTH  ,NP  ,NR,NSA)*WPM*DIVFPP*DPPM &
                   -DTP(NTH+1,NP  ,NR,NSA)*WTP*DIVDTP*DTPP &
                   +DTP(NTH  ,NP  ,NR,NSA)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            AL(NM,NL)=AL(NM,NL) &
                   -DTP(NTB  ,NP  ,NR,NSA)*WTB*DIVDTP*DTPM 
         ENDIF
      ENDIF

      IF(NP.NE.1.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP-1,NR)
         AL(NM,NL)=-DPT(NTH  ,NP  ,NR,NSA)*WPM*DIVDPT*DPTM &
                   -DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         IF(NTB.ne.0)THEN
            AL(NM,NL)=AL(NM,NL) &
                 -DPT(NTH  ,NP  ,NR,NSA)*WPM*DIVDPT*DPTM
         END IF
      ENDIF

      IF(NP.NE.1.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB,NP-1,NR)
         AL(NM,NL)= &
              +DPT(NTH, NP  ,NR,NSA)*WPBM*DIVDPT*DPTM &
              -DTP(NTB, NP  ,NR,NSA)*VTB *DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP,NR)
         AL(NM,NL)=-DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP &
                   +DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM &
                   +DTT(NTH  ,NP  ,NR,NSA)    *DIVDTT*DTTM &
                   +FTH(NTH  ,NP  ,NR,NSA)*WTM*DIVFTH*DTPM 
!         IF(NP.EQ.NPMAX) THEN !
!            AL(NM,NL)=AL(NM,NL) &
!                   -DTP(NTH  ,NP  ,NR,NSA)*VTM*DIVDTP*DTPM
!         ENDIF
      ENDIF

      IF(NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP,NR)
         AL(NM,NL)=+DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP &
                   -DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM &
                   +DTT(NTH+1,NP  ,NR,NSA)    *DIVDTT*DTTP &
                   -FTH(NTH+1,NP  ,NR,NSA)*VTP*DIVFTH*DTPP 
         IF(NTB.ne.0)THEN
            AL(NM,NL)=AL(NM,NL) &
                   +DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP &
                   -DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM 
         END IF
!         IF(NP.EQ.NPMAX) THEN !
!            AL(NM,NL)=AL(NM,NL) &
!                   +DTP(NTH+1,NP  ,NR,NSA)*WTP*DIVDTP*DTPP
!         ENDIF
      ENDIF

      IF(NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB,NP,NR)
         AL(NM,NL)=-DPT(NTH,NP+1,NR,NSA)*WPBP*DIVDPT*DPTP &
                   +DPT(NTH,NP  ,NR,NSA)*VPBM*DIVDPT*DPTM &
                   +DTT(NTB,NP  ,NR,NSA)    *DIVDTT*DTTM &
                   -FTH(NTB,NP  ,NR,NSA)*VTB*DIVFTH*DTPM
!         IF(NP.EQ.NPMAX) THEN
!            AL(NM,NL)=AL(NM,NL) &
!                   +DTP(NTB,NP  ,NR,NSA)*VTB*DIVDTP*DTPM
!         ENDIF
      ENDIF

      IF(NP.NE.NPMAX.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP+1,NR)
         AL(NM,NL)=-DPT(NTH  ,NP+1,NR,NSA)*VPP*DIVDPT*DPTP &
                   -DTP(NTH  ,NP  ,NR,NSA)*WTM*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(NP.NE.NPMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP+1,NR)
         AL(NM,NL)=+DPP(NTH  ,NP+1,NR,NSA)    *DIVDPP*DPPP &
                   -FPP(NTH  ,NP+1,NR,NSA)*VPP*DIVFPP*DPPP &
                   +DTP(NTH+1,NP  ,NR,NSA)*WTP*DIVDTP*DTPP &
                   -DTP(NTH  ,NP  ,NR,NSA)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            AL(NM,NL)=AL(NM,NL) &
                   +DTP(NTB  ,NP  ,NR,NSA)*WTB*DIVDTP*DTPM 
         ENDIF
      ENDIF

      IF(NP.NE.NPMAX.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP+1,NR)
         AL(NM,NL)=+DPT(NTH  ,NP+1,NR,NSA)*VPP*DIVDPT*DPTP &
                   +DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         IF(NTB.ne.0)THEN
            AL(NM,NL)=AL(NM,NL) &
                   +DPT(NTH  ,NP+1,NR,NSA)*VPP*DIVDPT*DPTP
         END IF
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(NP.NE.NPMAX.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB,NP+1,NR)
         AL(NM,NL)= &
              -DPT(NTH,NP+1,NR,NSA)*VPBP*DIVDPT*DPTP &
              +DTP(NTB,NP  ,NR,NSA)*VTB *DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(MODELD.GT.0) THEN!
!         IF(NR.NE.NRMAX.AND.NP.NE.NPMAX) THEN
         IF(NR.NE.NRMAX) THEN
            NL=NL+1
            LL(NM,NL)=NMA(NTH,NP,NR+1)
            AL(NM,NL)=DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
!                         *RLAMDAG(NTH,NR+1)/RFSADG(NR+1) &
                     -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR*DRRP 
!                         *RLAMDAG(NTH,NR+1)/RFSADG(NR+1)
            IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
               LL(NM,NL)=0
               AL(NM,NL)=0.D0
               NL=NL-1
            ENDIF
!            IF(NTBP.NE.0) THEN
!               NL=NL+1
!               LL(NM,NL)=NMA(NTBP,NP,NR+1)
!               AL(NM,NL)=DRR(NTBP ,NP  ,NR+1,NSA)     *DIVDRR*DRRP &
!                            *RLAMDAG(NTBP,NR+1)/RFSADG(NR+1) &
!                        -FRR(NTBP ,NP  ,NR+1,NSA)*VRBP*DIVFRR*DRRP &
!                            *RLAMDAG(NTBP,NR+1)/RFSADG(NR+1)
!               IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
!                  LL(NM,NL)=0
!                  AL(NM,NL)=0.D0
!                  NL=NL-1
!               ENDIF
!            ENDIF
         ENDIF
      ENDIF

      DL(NM)=-DPP(NTH  ,NP+1,NR  ,NSA)    *DIVDPP*DPPP &
             -FPP(NTH  ,NP+1,NR  ,NSA)*WPP*DIVFPP*DPPP &
             -DPP(NTH  ,NP  ,NR  ,NSA)    *DIVDPP*DPPM &
             +FPP(NTH  ,NP  ,NR  ,NSA)*VPM*DIVFPP*DPPM &
             -DTT(NTH+1,NP  ,NR  ,NSA)    *DIVDTT*DTTP &
             -FTH(NTH+1,NP  ,NR  ,NSA)*WTP*DIVFTH*DTPP &
!
             -DTT(NTH  ,NP  ,NR  ,NSA)    *DIVDTT*DTTM &
             +FTH(NTH  ,NP  ,NR  ,NSA)*VTM*DIVFTH*DTPM &
             +PPL(NTH  ,NP  ,NR  ,NSA)

      IF(NTB.NE.0) THEN
         DL(NM)=DL(NM) &
             -DTT(NTB,NP  ,NR  ,NSA)    *DIVDTT*DTTM &
             -FTH(NTB,NP  ,NR  ,NSA)*VTB*DIVFTH*DTPM
      ENDIF

!      IF(NP.EQ.NPMAX) THEN !
!         DL(NM)=DL(NM) &
!             +DTP(NTH+1,NP  ,NR  ,NSA)*WTP*DIVDTP*DTPP &
!             -DTP(NTH  ,NP  ,NR  ,NSA)*VTM*DIVDTP*DTPM
!         IF(NTB.NE.0) THEN
!            DL(NM)=DL(NM) &
!                +DTP(NTB,NP  ,NR  ,NSA)*VTB*DIVDTP*DTPM
!         ENDIF
!      ENDIF

      IF(MODELD.GT.0) THEN
         DL(NM)= DL(NM) &
              -DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) &
              -FRR(NTH  ,NP  ,NR+1,NSA)*WRP*DIVFRR*DRRP &
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) &
              -DRR(NTH  ,NP  ,NR  ,NSA)    *DIVDRR*DRRM &
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) &
              +FRR(NTH  ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM 
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) 
         IF(NTBM.ne.0)THEN
            DL(NM)= DL(NM) &
                 -DRR(NTBM ,NP  ,NR  ,NSA)    *DIVDRR*DRRM &
                               *RLAMDAG(NTBM,NR)/RFSADG(NR) &
                 +FRR(NTBM ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM &
                               *RLAMDAG(NTBM,NR)/RFSADG(NR)
         END IF
!            IF(NR.NE.NRMAX) THEN
!               DL(NM)= DL(NM) &
!                 -DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
!                     *RLAMDAG(NTH,NR)/RFSADG(NR) &
!                 -FRR(NTH  ,NP  ,NR+1,NSA)*WRP*DIVFRR*DRRP &
!                     *RLAMDAG(NTH,NR)/RFSADG(NR)
!               IF(NTBP.NE.0) THEN
!                  DL(NM)= DL(NM) &
!                    -DRR(NTBP ,NP  ,NR+1,NSA)     *DIVDRR*DRRP &
!                        *RLAMDAG(NTBP,NR)/RFSADG(NR) &
!                    -FRR(NTBP ,NP  ,NR+1,NSA)*WRBP*DIVFRR*DRRP &
!                        *RLAMDAG(NTBP,NR)/RFSADG(NR) 
!               ENDIF
!            ENDIF
!            IF(NR.NE.1) THEN
!               DL(NM)= DL(NM) &
!                 -DRR(NTH  ,NP  ,NR  ,NSA)    *DIVDRR*DRRM &
!                     *RLAMDAG(NTH,NR)/RFSADG(NR) &
!                 +FRR(NTH  ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM &
!                     *RLAMDAG(NTH,NR)/RFSADG(NR) 
!               IF(NTBM.NE.0) THEN
!                  DL(NM)= DL(NM) &
!                    -DRR(NTBM ,NP  ,NR  ,NSA)     *DIVDRR*DRRM &
!                        *RLAMDAG(NTBM,NR)/RFSADG(NR) &
!                    +FRR(NTBM ,NP  ,NR  ,NSA)*VRBM*DIVFRR*DRRM &
!                        *RLAMDAG(NTBM,NR)/RFSADG(NR) 
!               ENDIF
!            ENDIF
      ENDIF

      SPP(NTH,NP,NR,NSA) &
              =( SPPB(NTH,NP,NR,NSA) &
                +SPPF(NTH,NP,NR,NSA) &
                +SPPS(NTH,NP,NR,NSA) &
                +SPPI(NTH,NP,NR,NSA) )

      IF(MODELD.GT.0.AND.NR.EQ.NRMAX) THEN
         SPPD(NTH,NP,NSA)= FS2(NTH,NP,NSA) &
              *(DRR(NTH,NP,NR+1,NSA)    *DIVDRR &
               -FRR(NTH,NP,NR+1,NSA)*VRP*DIVFRR)*DRRP
         SPP(NTH,NP,NR,NSA) = SPP(NTH,NP,NR,NSA) &
              + SPPD(NTH,NP,NSA)
      ENDIF

!         IF(NP.NE.NPMAX) THEN
!            SPP(NTH,NP,NR,NSA)=SPP(NTH,NP,NR,NSA) &
!                 +FS2(NTH,NP,NSA)*(DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR &
!                                  -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR) &
!                      *DRRP
!            IF(NTBP.NE.0) THEN
!               SPP(NTH,NP,NR,NSA)=SPP(NTH,NP,NR,NSA) &
!                    +FS2(NTH,NP,NSA)*(DRR(NTBP ,NP  ,NR+1,NSA)     *DIVDRR &
!                                     -FRR(NTBP ,NP  ,NR+1,NSA)*VRBP*DIVFRR) &
!                      *DRRP
!            ENDIF
!         ENDIF
      RETURN
      END SUBROUTINE FPSETM

      END MODULE fpexec
