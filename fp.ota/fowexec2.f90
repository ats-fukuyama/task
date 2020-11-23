!     $Id: fpexec.f90,v 1.29 2013/02/08 07:36:24 nuga Exp $

!     **************************
!        EXECUTE TIME ADVANCE
!     **************************

      MODULE fpexec

      use fpcomm
      use fowcomm

      contains

      SUBROUTINE fp_exec(NSA,IERR,its)

      USE libmpi
      USE libmtx
      USE fpmpi
      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NN, NS
      integer:: NTHS, NLL
      integer:: IERR,its,i,j,ll1
      integer:: imtxstart1,imtxend1
      real(8),dimension(nmend-nmstart+1):: BM_L
      real(8),dimension(nthmax):: sendbuf_p, recvbuf_p
      real(8),dimension(nthmax*(npend-npstart+1)):: sendbuf_r, recvbuf_r

      NS=NS_NSA(NSA)

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


!     ----- Solve matrix equation -----

      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC) ! ncom is nessesary for MUMPS not PETSc
      ierr=0

!     ----- Get solution vector -----

      CALL mtx_get_vector(BM_L)

      DO NR=NRSTART, NREND
         DO NP=NPSTART, NPEND
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               IF(ABS(BM_L(NM-NMSTART+1)).LT.1.D-100) THEN
                  FNS0(NTH,NP,NR,NSA)=0.D0
               ELSE
                  FNS0(NTH,NP,NR,NSA)=BM_L(NM-NMSTART+1)
               END IF
            END DO
         END DO
      END DO
!     shadow requires to communicate
      CALL mtx_set_communicator(comm_np)
      DO NR=NRSTART, NREND
         CALL shadow_comm_np(NR,NSA)
      END DO
      CALL mtx_set_communicator(comm_nr)
      CALL shadow_comm_nr(NSA)
      CALL mtx_set_communicator(comm_nrnp) !3D


!     ----- Clean up matrix solver -----
      CALL mtx_cleanup

      CALL mtx_reset_communicator

      RETURN
      END SUBROUTINE FP_EXEC

!     ---------------------------------

      SUBROUTINE SET_FM_NMA(NSA,func_in)

      IMPLICIT NONE
      integer:: NTH, NP, NR, NSA, NM, NRS, NPS
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

      DO NR=NRSTART,NREND
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FM(NM)=func_in(NTH,NP,NR,NSA)
            ENDDO
         ENDDO
      ENDDO
      NR=NRSTARTW
      IF(NR.ne.NRSTART)THEN
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FM_shadow_m(NM)=func_in(NTH,NP,NR,NSA)
            ENDDO
         ENDDO
      END IF
      NR=NRENDWM
      IF(NR.ne.NREND)THEN
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FM_shadow_p(NM)=func_in(NTH,NP,NR,NSA)
            ENDDO
         ENDDO
      END IF

      
      END SUBROUTINE SET_FM_NMA

!
!     ***************************
!        CALCULATION OF WEIGHT
!     ***************************
!
      SUBROUTINE FPWEIGHT(NSA,IERR) ! proposed by Chang and Cooper [30] in Karney

      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NTHA, NTHB, NTB, NS
      integer:: IERR
      real(8):: DFDTH, FVEL, DVEL, DFDP, DFDB

!     +++++ calculation of weigthing (including off-diagonal terms) +++++

      real(8)::EPSWT=1.D-70

      NS=NS_NSA(NSA)
      DO NR=NRSTART,NREND
      DO NP=NPSTART,NPENDWG
      DO NTH=1,NTHMAX
         DFDTH=0.D0
         IF(NP.NE.1) THEN
            NTHA=MIN(NTH+1,NTHMAX)
            NTHB=MAX(NTH-1,1)
            IF(NP.EQ.NPMAX+1) THEN
               IF(ABS(F(NTH,NP-1,NR)).GT.EPSWT) THEN
                  DFDTH= (F(NTHA,NP-1,NR)-F(NTHB,NP-1,NR)) &
                        /(2.D0*PG(NP,NS)*DELTH*F(NTH,NP-1,NR))
               ENDIF
            ELSE
               IF(ABS(F(NTH,NP-1,NR)).GT.EPSWT) THEN
                  IF(ABS(F(NTH,NP,NR)).GT.EPSWT) THEN
                     DFDTH= (F(NTHA,NP-1,NR)-F(NTHB,NP-1,NR)) &
                           /(4.D0*PG(NP,NS)*DELTH*F(NTH,NP-1,NR)) &
                          + (F(NTHA,NP  ,NR)-F(NTHB,NP  ,NR)) &
                           /(4.D0*PG(NP,NS)*DELTH*F(NTH,NP  ,NR))
                  ELSE
                     DFDTH= (F(NTHA,NP-1,NR)-F(NTHB,NP-1,NR)) &
                           /(2.D0*PG(NP,NS)*DELTH*F(NTH,NP-1,NR)) 
                  ENDIF
               ELSE
                  IF(ABS(F(NTH,NP,NR)).GT.EPSWT) THEN
                     DFDTH= (F(NTHA,NP  ,NR)-F(NTHB,NP  ,NR)) &
                           /(2.D0*PG(NP,NS)*DELTH*F(NTH,NP  ,NR))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         FVEL=FPP(NTH,NP,NR,NSA)-DPT(NTH,NP,NR,NSA)*DFDTH
         WEIGHP(NTH,NP,NR,NSA)=FPWEGH(-DELP(NS)*FVEL,DPP(NTH,NP,NR,NSA))
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NREND
!      DO NP=1,NPMAX
      DO NP=NPSTART,NPEND
         DFDP=-PM(NP,NS)*RTFP0(NSA)/RTFP(NR,NSA)
         DFDB=DFDP
      DO NTH=1,NTHMAX+1
        IF(NTH.EQ.1) THEN
            IF(F(NTH,NP,NR).GT.EPSWT) THEN
               IF(NP.EQ.1) THEN
                  DFDP= (F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                       /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
               ELSEIF(NP.EQ.NPMAX) THEN
                  DFDP= (F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                       /(     DELP(NS)*F(NTH  ,NP,NR))
               ELSE
                  DFDP= (F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                       /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
               ENDIF
            ENDIF
         ELSEIF(NTH.EQ.NTHMAX+1) THEN
            IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
               IF(NP.EQ.1) THEN
                  DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                       /(2.D0*DELP(NS)*F(NTH-1,NP,NR))
               ELSEIF(NP.EQ.NPMAX) THEN
                  DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                       /(     DELP(NS)*F(NTH-1,NP,NR))
               ELSE
                  DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                       /(2.D0*DELP(NS)*F(NTH-1,NP,NR))
               ENDIF
            ENDIF
         ELSE IF(NTH.EQ.ITL(NR)) THEN
            NTB=ITU(NR)
            IF(NP.EQ.1) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSEIF(NP.EQ.NPMAX) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(ABS(F(NTH,NP,NR)).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(     DELP(NS)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(     DELP(NS)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSE
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ENDIF
            IF(NP.EQ.1) THEN
               IF(F(NTB-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTHMAX-NTB+2,1,NR)) &
                          /(4.D0*DELP(NS)*F(NTB-1,NP,NR)) &
                          +(F(NTB  ,NP+1,NR)-F(NTHMAX-NTB+1,1,NR)) &
                          /(4.D0*DELP(NS)*F(NTB  ,NP,NR))
                  ELSE
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTHMAX-NTB+2,1,NR)) &
                          /(2.D0*DELP(NS)*F(NTB-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB  ,NP+1,NR)-F(NTHMAX-NTB+1,1,NR)) &
                          /(2.D0*DELP(NS)*F(NTB  ,NP,NR))
                  ENDIF
               ENDIF
            ELSEIF(NP.EQ.NPMAX) THEN
               IF(F(NTB-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB-1,NP,NR)-F(NTB-1,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTB-1,NP,NR)) &
                          +(F(NTB  ,NP,NR)-F(NTB  ,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTB  ,NP,NR))
                  ELSE
                     DFDB= (F(NTB-1,NP,NR)-F(NTB-1,NP-1,NR)) &
                          /(     DELP(NS)*F(NTB-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB  ,NP,NR)-F(NTB  ,NP-1,NR)) &
                          /(     DELP(NS)*F(NTB  ,NP,NR))
                  ENDIF
               ENDIF
            ELSE
               IF(F(NTB-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTB-1,NP-1,NR)) &
                          /(4.D0*DELP(NS)*F(NTB-1,NP,NR)) &
                          +(F(NTB  ,NP+1,NR)-F(NTB  ,NP-1,NR)) &
                          /(4.D0*DELP(NS)*F(NTB  ,NP,NR))
                  ELSE
                     DFDB= (F(NTB-1,NP+1,NR)-F(NTB-1,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTB-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTB,NP,NR).GT.EPSWT) THEN
                     DFDB= (F(NTB  ,NP+1,NR)-F(NTB  ,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTB  ,NP,NR))
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF(NP.EQ.1) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTHMAX-NTH+2,1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTHMAX-NTH+1,1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSEIF(NP.EQ.NPMAX) THEN
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP,NR)-F(NTH-1,NP-1,NR)) &
                          /(     DELP(NS)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP,NR)-F(NTH  ,NP-1,NR)) &
                          /(     DELP(NS)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ELSE
               IF(F(NTH-1,NP,NR).GT.EPSWT) THEN
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH-1,NP,NR)) &
                          +(F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(4.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ELSE
                     DFDP= (F(NTH-1,NP+1,NR)-F(NTH-1,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH-1,NP,NR))
                  ENDIF
               ELSE
                  IF(F(NTH,NP,NR).GT.EPSWT) THEN
                     DFDP= (F(NTH  ,NP+1,NR)-F(NTH  ,NP-1,NR)) &
                          /(2.D0*DELP(NS)*F(NTH  ,NP,NR))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF(NTH.EQ.ITL(NR)) THEN
            FVEL=FCTH(NTH,NP,NR,NSA)-DCTP(NTH,NP,NR,NSA)*DFDP &
                -FCTH(NTB,NP,NR,NSA)+DCTP(NTB,NP,NR,NSA)*DFDB
            DVEL=DCTT(NTH,NP,NR,NSA)+DCTT(NTB,NP,NR,NSA)
         ELSE
            FVEL=FCTH(NTH,NP,NR,NSA)-DCTP(NTH,NP,NR,NSA)*DFDP
            DVEL=DCTT(NTH,NP,NR,NSA)
         ENDIF
         WEIGHT(NTH,NP,NR,NSA) &
                 =FPWEGH(-DELTH*PM(NP,NS)*FVEL,DVEL)
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NRENDWG
      DO NP=NPSTART,NPEND
      DO NTH=1,NTHMAX
         FVEL=FRR(NTH,NP,NR,NSA)
         DVEL=DRR(NTH,NP,NR,NSA)
         WEIGHR(NTH,NP,NR,NSA)=FPWEGH(-DELR*FVEL,DVEL)
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

      integer:: IERR, NS
      real(8):: DFDTH, FVEL, DVEL, DFDP, DFDB
      real(8):: DPPM, DPPP, DPTM, DPTP, SL, DTPM, DTPP, DTTM, DTTP
      ! real(8):: WPM, WPP, WTM, WTP, VPM, VPP, VTM, VTP
      real(8):: WTB, VTB, WRBM, VRBM, WRBP, VRBP
      real(8):: DIVDPP, DIVDPT, DIVDTP, DIVDTT, DIVFPP, DIVFTH
      real(8):: RL,DRRM,DRRP,WRM,WRP,VRM,VRP,DIVDRR,DIVFRR
      real(8):: PL
      DOUBLE PRECISION:: WPBM, VPBM, WPBP, VPBP
      ! fow extension
      double precision :: DPRM, DPRP, DTRM, DTRP, DRPM, DRPP, DRTM, DRTP
      double precision :: DIVDPR, DIVDTR, DIVDRP, DIVDRT
      double precision :: deltap, deltath, deltaps
      double precision :: wp(-1:1,0:1,-1:1), wt(0:1,-1:1,-1:1), wr(-1:1,-1:1,0:1)
      double precision :: vp(-1:1,0:1,-1:1), vt(0:1,-1:1,-1:1), vr(-1:1,-1:1,0:1)
      double precision :: DIJP(3,3), DIJM(3,3), DIVDIJ(3,3), del(3), W(3,-1:1,-1:1,0:1), V(3,-1:1,-1:1,0:1)
      integer :: i,j,k

      NS=NS_NSA(NSA)

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
            IF(NTH.GE.ITL(NR).AND.NTH.LT.ITL(NR-1)) THEN 
               NTBM=NTHMAX+1-NTH
            ENDIF
         ENDIF
      END IF ! MODELA

      NL=0
      NM=NMA(NTH,NP,NR)

      ! normalize coefficients
      PL   = PM(NP,NS)
      DPPM = PG(NP,NS  )**2
      DPPP = PG(NP+1,NS)**2
      DPTM = PG(NP,NS  )
      DPTP = PG(NP+1,NS)
      DPRM = PG(NP,NS  )**2
      DPRP = PG(NP+1,NS)**2

      SL   = sin(thetam(nth,np,nr,nsa))
      DTPM = sin(thetamg(nth,np,nr,nsa))
      DTPP = sin(thetamg(nth+1,np,nr,nsa))
      DTTM = sin(thetamg(nth,np,nr,nsa))/pl
      DTTP = sin(thetamg(nth+1,np,nr,nsa))pl
      DTRM = sin(thetamg(nth,np,nr,nsa))
      DTRP = sin(thetamg(nth+1,np,nr,nsa))

      RL   = psim(nr)
      DRPM = psimg(nr  )
      DRPP = psimg(nr+1)
      DRTM = psimg(nr  )/pl
      DRTP = psimg(nr+1)/pl
      DRRM = psimg(nr  )
      DRRP = psimg(nr+1)

      ! weight function, delta, epsilon
      do i = -1, 1
         do j = -1, 1
            do k = -1, 1
               if ( i /= -1 ) then
                  wt(i,j,k) = WEIGHT(NTH+i,NP+j,NR+k,NSA)
                  vt(i,j,k) = 1.d0 - wt(i,j,k)   
               end if
               
               if ( j /= -1 ) then
                  wp(i,j,k) = WEIGHP(NTH+i,NP+j,NR+k,NSA)
                  vp(i,j,k) = 1.d0 - wp(i,j,k)   
               end if

               if ( k /= -1 ) then
                  wr(i,j,k) = WEIGHR(NTH+i,NP+j,NR+k,NSA)
                  vr(i,j,k) = 1.d0 - wr(i,j,k)   
               end if

            end do
         end do
      end do

      ! discretized (div(d/dX))_Y
      DIVDPP = 1.d0/(pl**2 * deltap**2)
      DIVDPT = 1.d0/(pl**2 * deltap * deltath * 2.d0)
      DIVDPR = 1.d0/(pl**2 * deltap * deltaps * 2.d0)
      DIVDTP = 1.d0/(pl * sl * deltap * deltath * 2.d0)
      DIVDTT = 1.d0/(pl * sl * deltath**2 )
      DIVDTR = 1.d0/(pl * sl * deltath * deltaps * 2.d0)
      DIVDRP = 1.d0/(rl * deltap  * deltaps * 2.d0)
      DIVDRT = 1.d0/(rl * deltath * deltaps * 2.d0)
      DIVDRR = 1.d0/(rl * deltaps**2)
      DIVFPP = 1.d0/(pl**2 * deltap)
      DIVFTH = 1.d0/(pl * sl * deltath)
      DIVFRR = 1.d0/(rl * deltaps)

      ! f00-
      IF( nr /= 1 .and. np /= npmax .and. nth /= nthmax ) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH,NP,NR-1)
         AL(NM,NL)=DPR(NTH  ,NP+1,NR  ,NSA)*DIVDPR*DPRP*wp( 0, 1,-1)&
                  +DPR(NTH  ,NP  ,NR  ,NSA)*DIVDPR*DPRM*vp( 0, 0,-1)&
                  -DTR(NTH+1,NP  ,NR  ,NSA)*DIVDTR*DTRP*wt( 1, 0,-1)&
                  +DTR(NTH  ,NP  ,NR  ,NSA)*DIVDTR*DTRM*vt( 0, 0,-1)&
                  +DRR(NTH  ,NP  ,NR  ,NSA)*DIVDRR*DRRM             &
                  +FRR(NTH  ,NP  ,NR  ,NSA)*DIVFRR*DRRM*wr( 0, 0, 0)
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      ! f++0
      IF(NP.NE.1.AND.NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP+1,NR)
         AL(NM,NL)=DPT(nth  ,np+1,nr  ,nsa)*divdpt*dpt*vp(1,1,0)&
                  +DTP(nth+1,np  ,nr  ,nsa)*divdtp*dtp*vp(1,1,0)
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      do i = -1, 1
         do j = -1, 1
            do k = -1, 1
               if ( i*j*k /= 0 ) then                          ! (p,thetam,psim)->(     integer,     integer,     integer)
                  cycle
               else if ( j*k /= 0 .and. i = 0 ) then           ! (p,thetam,psim)->(half integer,     integer,     integer)
               else if ( k*i /= 0 .and. j = 0 ) then           ! (p,thetam,psim)->(     integer,half integer,     integer)
               else if ( i*j /= 0 .and. k = 0 ) then           ! (p,thetam,psim)->(     integer,     integer,half integer)
               else if ( i = 0 .and. j = 0 .and. k /= 0 ) then ! (p,thetam,psim)->(half integer,half integer,     integer)
               else if ( j = 0 .and. k = 0 .and. i /= 0 ) then ! (p,thetam,psim)->(     integer,half integer,half integer)
               else if ( k = 0 .and. i = 0 .and. j /= 0 ) then ! (p,thetam,psim)->(half integer,     integer,half integer)
               else if ( i = 0 .and. j = 0 .and. k = 0 ) then  ! (p,thetam,psim)->(half integer,half integer,half integer)
               end if
            end do
         end do
      end do


      SPP(NTH,NP,NR,NSA) &
              =( SPPB(NTH,NP,NR,NSA) &
                +SPPF(NTH,NP,NR,NSA) &
                +SPPS(NTH,NP,NR,NSA) &
                +SPPL(NTH,NP,NR,NSA) &
                +SPPI(NTH,NP,NR,NSA) &
                +SPPL_CX(NTH,NP,NR,NSA) )

      IF(MODELD.GT.0.AND.NR.EQ.NRMAX) THEN
         SPPD(NTH,NP,NSA)= FS2(NTH,NP,NSA) &
              *(DRR(NTH,NP,NR+1,NSA)    *DIVDRR &
               -FRR(NTH,NP,NR+1,NSA)*VRP*DIVFRR)*DRRP/RLAMDA_RG(NTH,NRMAX+1)
         SPP(NTH,NP,NR,NSA) = SPP(NTH,NP,NR,NSA) &
              + SPPD(NTH,NP,NSA) 
      ENDIF

      RETURN
      END SUBROUTINE FPSETM

      END MODULE fpexec
