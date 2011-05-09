!     $Id$

!     **************************
!        EXECUTE TIME ADVANCE
!     **************************

      MODULE fpexec

      use fpcomm

      contains

      SUBROUTINE fp_exec(NSA,IERR)

      USE libmtx
      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NSBA
      integer:: NTHS, NLL
      integer:: IERR,its,i,j,ll1
      integer:: imtxstart1,imtxend1

      NSBA=NSB_NSA(NSA)

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

!     ----- Set up index array -----

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NM=NTH+NTHMAX*(NP-1)+NPMAX*NTHMAX*(NR-1)
               NMA(NTH,NP,NR)=NM
               FM(NM)=FNS22(NTH,NP,NR,NSBA)
            ENDDO
         ENDDO
      ENDDO

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
            DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               CALL FPSETM(NTH,NP,NR,NSA,NLMAX(NM))
            ENDDO
            ENDDO
         ELSE
            DO NP=1,NPMAX
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
!         write(6,*) 'NLMAX=',NLMAX
      ENDDO

!     ----- Diagonal term -----

      DO NR=NRSTART,NREND
         IF(MODELA.EQ.0) THEN
            DO NP=1,NPMAX
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
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX/2
                  NM=NMA(NTH,NP,NR)
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                        +DELT*SPP(NTH,NP,NR,NSA)
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
                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM)) &
                         *FM(NM) &
                        +DELT*SPP(NTH,NP,NR,NSA)
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

      DO NM=NMSTART,NMEND
         IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
            DO NL=1,NLMAX(NM)
               IF(LL(NM,NL).NE.0) THEN
                  CALL mtx_set_matrix(nm,LL(NM,NL),-RIMPL*DELT*AL(NM,NL))
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DO NR=NRSTART,NREND
      DO NP=1,NPMAX
      DO NTH=NTHMAX/2+1,ITU(NR)
         NM=NMA(NTH,NP,NR)
         IF(LL(NM,NL).NE.0) WRITE(6,'(A,5I5,1PE12.4)') &
              'NR,NP,NTH,NM.NL,AL=',NR,NP,NTH,NM,NL,AL(NM,NL)
      ENDDO
      ENDDO
      ENDDO

!     ----- Source vector: contribution from off-diagonal term -----

      DO NM=NMSTART,NMEND
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

      CALL mtx_solve(imtx,epsm,its)
      if(nrank.eq.0) then
         write(6,*) 'Number of iterations =',its
      endif
      ierr=0

!     ----- Get solution vector -----

      CALL mtx_gather_vector(BMTOT)
      
      DO NR=NRSTART,NREND
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               F1(NTH,NP,NR)=BMTOT(NM)
            ENDDO
         ENDDO
      ENDDO

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               FNS(NTH,NP,NR,NSBA)=BMTOT(NM)
            ENDDO
         ENDDO
      ENDDO

!     ----- Clean up matrix solver -----

      CALL mtx_cleanup

      RETURN
      END SUBROUTINE FP_EXEC

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
      DO NP=1,NPMAX+1
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
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NREND
      DO NP=1,NPMAX
         DFDP=-PM(NP,NSBA)*RTFP0(NSA)/RTFP(NR,NSA)
!         DFDP=-PM(NP,NSBA)*RTFP0(NSA)/RTFP(NR,NSA)/SQRT(1.D0+THETA0(NSA)*PM(NP,NSBA)**2)
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
            FVEL=FTH(NTH,NP,NR,NSA)-DTP(NTH,NP,NR,NSA)*DFDP &
                -FTH(NTB,NP,NR,NSA)+DTP(NTB,NP,NR,NSA)*DFDB
            DVEL=DTT(NTH,NP,NR,NSA)+DTT(NTB,NP,NR,NSA)
         ELSE
            FVEL=FTH(NTH,NP,NR,NSA)-DTP(NTH,NP,NR,NSA)*DFDP
            DVEL=DTT(NTH,NP,NR,NSA)
         ENDIF
         WEIGHT(NTH,NP,NR,NSA) &
                 =FPWEGH(-DELTH*PM(NP,NSBA)*FVEL,DVEL)
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NREND+1
      DO NP=1,NPMAX
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
      integer:: NSA, NP, NTH, NR, NL, NM, NTHA, NTHB, NTB, NTBM, NTBP
      integer:: IERR, NSBA
      real(8):: DFDTH, FVEL, DVEL, DFDP, DFDB
      real(8):: DPPM, DPPP, DPTM, DPTP, SL, DTPM, DTPP, DTTM, DTTP
      real(8):: WPM, WPP, WTM, WTP, VPM, VPP, VTM, VTP
      real(8):: WTB, VTB, WRBM, VRBM, WRBP, VRBP
      real(8):: DIVDPP, DIVDPT, DIVDTP, DIVDTT, DIVFPP, DIVFTH
      real(8):: RL,DRRM,DRRP,WRM,WRP,VRM,VRP,DIVDRR,DIVFRR
      real(8):: PL

      NSBA=NSB_NSA(NSA)

      NTB=0
      NTBM=0
      NTBP=0
      IF(NTH.EQ.ITL(NR)) THEN
         NTB=ITU(NR)
      ENDIF
      IF(NR-1.GE.1) THEN
         IF(NTH.GE.ITL(NR).AND.NTH.LT.ITL(NR-1)) THEN
            NTBM=NTHMAX-NTH+1
         ENDIF
      ENDIF
      IF(NR+1.LE.NRMAX) THEN
         IF(NTH.LE.ITL(NR).AND.NTH.GT.ITL(NR+1)) THEN
            NTBP=NTHMAX-NTH+1
         ENDIF
      ENDIF

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
         DTPM=0.5D0*DTPM
         DTTM=0.5D0*DTTM
      ENDIF
      RL=RM(NR)
      DRRM=RG(NR  )
      DRRP=RG(NR+1)
      IF(NTBM.NE.0) THEN
         DRRM=0.5D0*DRRM
      ENDIF
      IF(NTBP.NE.0) THEN
         DRRP=0.5D0*DRRP
      ENDIF
      WPM=WEIGHP(NTH  ,NP  ,NR  ,NSA)
      WPP=WEIGHP(NTH  ,NP+1,NR  ,NSA)
      WTM=WEIGHT(NTH  ,NP  ,NR  ,NSA)
      WTP=WEIGHT(NTH+1,NP  ,NR  ,NSA)
      WRM=WEIGHR(NTH  ,NP  ,NR  ,NSA)
      WRP=WEIGHR(NTH  ,NP  ,NR+1,NSA)
      VPM=1.D0-WPM
      VPP=1.D0-WPP
      VTM=1.D0-WTM
      VTP=1.D0-WTP
      VRM=1.D0-WRM
      VRP=1.D0-WRP
      IF(NTB.NE.0) THEN
         WTB=WEIGHT(NTB+1,NP  ,NR  ,NSA)
         VTB=1.D0-WTB
      ENDIF
      IF(NTBM.NE.0) THEN
         WRBM=WEIGHR(NTBM,NP  ,NR-1  ,NSA)
         VRBM=1.D0-WRBM
      ENDIF
      IF(NTBP.NE.0) THEN
         WRBP=WEIGHR(NTBP,NP  ,NR+1  ,NSA)
         VRBP=1.D0-WRBP
      ENDIF
      DIVDPP=1.D0/(     PL*PL*DELP(NSBA) *DELP(NSBA))
      DIVDPT=1.D0/(2.D0*PL*PL*DELP(NSBA) *DELTH)
      DIVDTP=1.D0/(2.D0*PL*SL*DELTH*DELP(NSBA))
      DIVDTT=1.D0/(     PL*SL*DELTH*DELTH)
      DIVFPP=1.D0/(     PL*PL*DELP(NSBA))
      DIVFTH=1.D0/(     PL*SL*DELTH)
      DIVDRR=1.D0/(     RL   *DELR *DELR)
      DIVFRR=1.D0/(     RL   *DELR)
      IF(NP.EQ.NPMAX) THEN
         DIVDTP=2.D0*DIVDTP
      ENDIF

      IF(MODELD.GT.0) THEN
         IF(NR.NE.1.AND.NP.NE.NPMAX) THEN
            NL=NL+1
            LL(NM,NL)=NMA(NTH,NP,NR-1)
            AL(NM,NL)=DRR(NTH  ,NP  ,NR,NSA)    *DIVDRR*DRRM &
                     +FRR(NTH  ,NP  ,NR,NSA)*WRM*DIVFRR*DRRM
            IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
               LL(NM,NL)=0
               AL(NM,NL)=0.D0
               NL=NL-1
            ENDIF
            IF(NTBM.NE.0) THEN
               NL=NL+1
               LL(NM,NL)=NMA(NTBM,NP,NR-1)
               AL(NM,NL)=DRR(NTBM ,NP  ,NR,NSA)     *DIVDRR*DRRM &
                        +FRR(NTBM ,NP  ,NR,NSA)*WRBM*DIVFRR*DRRM
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
                   -DTP(NTB+1,NP  ,NR,NSA)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF

      IF(NP.NE.1.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP-1,NR)
         AL(NM,NL)=-DPT(NTH  ,NP  ,NR,NSA)*WPM*DIVDPT*DPTM &
                   -DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
      ENDIF

      IF(NP.NE.1.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP-1,NR)
         AL(NM,NL)=-DTP(NTB+1,NP  ,NR,NSA)*WTB*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(NTH.NE.1) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH-1,NP,NR)
         AL(NM,NL)=+DTT(NTH  ,NP  ,NR,NSA)    *DIVDTT*DTTM &
                   +FTH(NTH  ,NP  ,NR,NSA)*WTM*DIVFTH*DTPM &
                   -DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP &
                   +DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL) &
                   -DTP(NTH  ,NP  ,NR,NSA)*WTM*DIVDTP*DTPM
         ENDIF
      ENDIF

      IF(NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP,NR)
         AL(NM,NL)=+DTT(NTH+1,NP  ,NR,NSA)    *DIVDTT*DTTP &
                   -FTH(NTH+1,NP  ,NR,NSA)*VTP*DIVFTH*DTPP &
                   +DPT(NTH  ,NP+1,NR,NSA)*WPP*DIVDPT*DPTP &
                   -DPT(NTH  ,NP  ,NR,NSA)*VPM*DIVDPT*DPTM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL) &
                   +DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         ENDIF
      ENDIF

      IF(NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP,NR)
         AL(NM,NL)= DTT(NTB+1,NP  ,NR,NSA)    *DIVDTT*DTTM &
                   -FTH(NTB+1,NP  ,NR,NSA)*WTB*DIVFTH*DTPM
         IF(NP.EQ.NPMAX) THEN
            AL(NM,NL)=AL(NM,NL) &
                   +DTP(NTB+1,NP  ,NR,NSA)*WTB*DIVDTP*DTPM
         ENDIF
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
                   +DTP(NTB+1,NP  ,NR,NSA)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF

      IF(NP.NE.NPMAX.AND.NTH.NE.NTHMAX) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTH+1,NP+1,NR)
         AL(NM,NL)=+DPT(NTH  ,NP+1,NR,NSA)*VPP*DIVDPT*DPTP &
                   +DTP(NTH+1,NP  ,NR,NSA)*VTP*DIVDTP*DTPP
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(NP.NE.NPMAX.AND.NTB.NE.0) THEN
         NL=NL+1
         LL(NM,NL)=NMA(NTB+1,NP+1,NR)
         AL(NM,NL)=+DTP(NTB+1,NP  ,NR,NSA)*WTB*DIVDTP*DTPM
         IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
            LL(NM,NL)=0
            AL(NM,NL)=0.D0
            NL=NL-1
         ENDIF
      ENDIF

      IF(MODELD.GT.0) THEN!
         IF(NR.NE.NRMAX.AND.NP.NE.NPMAX) THEN
            NL=NL+1
            LL(NM,NL)=NMA(NTH,NP,NR+1)
            AL(NM,NL)=DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
                     -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR*DRRP
            IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
               LL(NM,NL)=0
               AL(NM,NL)=0.D0
               NL=NL-1
            ENDIF
            IF(NTBP.NE.0) THEN
               NL=NL+1
               LL(NM,NL)=NMA(NTBP,NP,NR+1)
               AL(NM,NL)=DRR(NTBP ,NP  ,NR+1,NSA)     *DIVDRR*DRRP &
                        -FRR(NTBP ,NP  ,NR+1,NSA)*VRBP*DIVFRR*DRRP
               IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
                  LL(NM,NL)=0
                  AL(NM,NL)=0.D0
                  NL=NL-1
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      DL(NM)=-DPP(NTH  ,NP+1,NR  ,NSA)    *DIVDPP*DPPP &
             -FPP(NTH  ,NP+1,NR  ,NSA)*WPP*DIVFPP*DPPP &
             -DPP(NTH  ,NP  ,NR  ,NSA)    *DIVDPP*DPPM &
             +FPP(NTH  ,NP  ,NR  ,NSA)*VPM*DIVFPP*DPPM &
             -DTT(NTH+1,NP  ,NR  ,NSA)    *DIVDTT*DTTP &
             -FTH(NTH+1,NP  ,NR  ,NSA)*WTP*DIVFTH*DTPP &
             -DTT(NTH  ,NP  ,NR  ,NSA)    *DIVDTT*DTTM &
             +FTH(NTH  ,NP  ,NR  ,NSA)*VTM*DIVFTH*DTPM &
             +PPL(NTH  ,NP  ,NR  ,NSA)                 

      IF(NTB.NE.0) THEN
         DL(NM)=DL(NM) &
             -DTT(NTB+1,NP  ,NR  ,NSA)    *DIVDTT*DTTM &
             -FTH(NTB+1,NP  ,NR  ,NSA)*VTB*DIVFTH*DTPM
      ENDIF

      IF(NP.EQ.NPMAX) THEN
         DL(NM)=DL(NM) &
             +DTP(NTH+1,NP  ,NR  ,NSA)*WTP*DIVDTP*DTPP &
             -DTP(NTH  ,NP  ,NR  ,NSA)*VTM*DIVDTP*DTPM
         IF(NTB.NE.0) THEN
            DL(NM)=DL(NM) &
                +DTP(NTB+1,NP  ,NR  ,NSA)*VTB*DIVDTP*DTPM
         ENDIF
      ENDIF

      IF(MODELD.GT.0) THEN
         IF(NP.NE.NPMAX) THEN
            IF(NR.NE.NRMAX) THEN
               DL(NM)= DL(NM) &
                 -DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
                 -FRR(NTH  ,NP  ,NR+1,NSA)*WRP*DIVFRR*DRRP
               IF(NTBP.NE.0) THEN
                  DL(NM)= DL(NM) &
                    -DRR(NTBP ,NP  ,NR+1,NSA)     *DIVDRR*DRRP &
                    -FRR(NTBP ,NP  ,NR+1,NSA)*WRBP*DIVFRR*DRRP 
               ENDIF
            ENDIF
            IF(NR.NE.1) THEN
               DL(NM)= DL(NM) &
                 -DRR(NTH  ,NP  ,NR  ,NSA)    *DIVDRR*DRRM &
                 +FRR(NTH  ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM 
               IF(NTBM.NE.0) THEN
                  DL(NM)= DL(NM) &
                    -DRR(NTBM ,NP  ,NR  ,NSA)     *DIVDRR*DRRM &
                    +FRR(NTBM ,NP  ,NR  ,NSA)*VRBM*DIVFRR*DRRM 
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      SPP(NTH,NP,NR,NSA) &
              =( SPPB(NTH,NP,NR,NSA) &
                +SPPF(NTH,NP,NR,NSA) &
                +SPPS(NTH,NP,NR,NSA) )*RLAMDAG(NTH,NR)
      IF(MODELD.GT.0.AND.NR.EQ.NRMAX) THEN
         IF(NP.NE.NPMAX) THEN
            SPP(NTH,NP,NR,NSA)=SPP(NTH,NP,NR,NSA) &
                 +FS2(NTH,NP,NSA)*(DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR &
                                  -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR) &
                      *DRRP
            IF(NTBP.NE.0) THEN
               SPP(NTH,NP,NR,NSA)=SPP(NTH,NP,NR,NSA) &
                    +FS2(NTH,NP,NSA)*(DRR(NTBP ,NP  ,NR+1,NSA)     *DIVDRR &
                                     -FRR(NTBP ,NP  ,NR+1,NSA)*VRBP*DIVFRR) &
                      *DRRP
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE FPSETM

      END MODULE fpexec
