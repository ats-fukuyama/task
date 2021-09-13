!     $Id: fpdrexec.f90,v 1.7 2013/01/22 16:21:46 fukuyama Exp $

!     **************************
!        EXECUTE TIME ADVANCE
!     **************************

      MODULE fpdrexec

      use fpcomm
      use fpexec

      contains

      SUBROUTINE fp_drexec(NSA,IERR,its1)

      USE libmpi
      USE libmtx
      IMPLICIT NONE
      integer:: NSA, NP, NTH, NR, NL, NM, NSBA
      integer:: NTHS, NLL
      integer:: IERR,its,i,j,ll1,its1
      integer:: imtxstart1,imtxend1
!      integer,optional:: methodKSP, methodPC

      NSBA=NSB_NSA(NSA)

      CALL mtx_set_communicator(comm_nrnp) !3D

!     ----- Set up matrix solver -----
      CALL mtx_setup(imtxsize,imtxstart1,imtxend1,imtxwidth)
      IF(imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend) THEN
         WRITE(6,*) 'XX fp_drexec: '
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

      CALL SET_FM_NMA(NSA,FNS0)
!      CALL SET_FM_NMA(NSA,FNSP)

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
               CALL FPSETMDR(NTH,NP,NR,NSA,NLMAX(NM))
            ENDDO
            ENDDO
         ELSE
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX/2
                  NM=NMA(NTH,NP,NR)
                  CALL FPSETMDR(NTH,NP,NR,NSA,NLMAX(NM))
               ENDDO
               DO NTH=ITU(NR)+1,NTHMAX
                  NM=NMA(NTH,NP,NR)
                  CALL FPSETMDR(NTH,NP,NR,NSA,NLMAX(NM))
               ENDDO
            ENDDO
         ENDIF
!         write(6,*) 'NLMAX=',NLMAX
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
!                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                  BM(NM)=(1.D0+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                        +DELT*SPP(NTH,NP,NR,NSA)
                  IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
                     CALL mtx_set_matrix(nm,nm, &
!                                         RLAMDA(NTH,NR)-RIMPL*DELT*DL(NM))
                                         1.D0-RIMPL*DELT*DL(NM))
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
!                  BM(NM)=(RLAMDA(NTH,NR)+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                  BM(NM)=(1.D0+(1.D0-RIMPL)*DELT*DL(NM))*FM(NM) &
                        +DELT*SPP(NTH,NP,NR,NSA)
                  IF(nm.GE.imtxstart.AND.nm.LE.imtxend) THEN
                     CALL mtx_set_matrix(nm,nm, &
!                                         RLAMDA(NTH,NR)-RIMPL*DELT*DL(NM))
                                         1.D0-RIMPL*DELT*DL(NM))
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

      CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC)! ncom is necessary for MUMPS
      if(nrank.eq.0) then
         write(6,'(A,3I4)') 'Number of iterations p, r, NSA =',its1, its, NSA
      endif
      ierr=0

!     ----- Get solution vector -----

      CALL mtx_gather_vector(BMTOT)
      
      DO NR=NRSTART,NREND
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX
               NM=NMA(NTH,NP,NR)
               F1(NTH,NP,NR)=BMTOT(NM)
            ENDDO
         ENDDO
      ENDDO

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
      END SUBROUTINE FP_DREXEC

! ******************************************
!     Calculation of matrix coefficients
! ******************************************

      SUBROUTINE FPSETMDR(NTH,NP,NR,NSA,NL)

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
      REAL(rkind):: DFDTH, FVEL, DVEL, DFDP, DFDB
      REAL(rkind):: DPPM, DPPP, DPTM, DPTP, SL, DTPM, DTPP, DTTM, DTTP
      REAL(rkind):: WPM, WPP, WTM, WTP, VPM, VPP, VTM, VTP
      REAL(rkind):: WTB, VTB, WRBM, VRBM, WRBP, VRBP
      REAL(rkind):: DIVDPP, DIVDPT, DIVDTP, DIVDTT, DIVFPP, DIVFTH
      REAL(rkind):: RL,DRRM,DRRP,WRM,WRP,VRM,VRP,DIVDRR,DIVFRR
      REAL(rkind):: PL
      REAL(rkind):: WPBM, VPBM, WPBP, VPBP

      NSBA=NSB_NSA(NSA)

      NTBM=0
      NTBP=0
      IF(MODELA.ne.0)THEN
         IF(NR-1.GE.1) THEN
            IF(NTH.GE.ITL(NR).AND.NTH.LT.ITL(NR-1)) THEN 
               NTBM=NTHMAX+1-NTH
            ENDIF
         ENDIF
      END IF ! MODELA

      NL=0
      NM=NMA(NTH,NP,NR)
      RL=RM(NR)
      DRRM=RG(NR  )
      DRRP=RG(NR+1)
      IF(NTBM.NE.0) THEN
         DRRM=0.5D0*DRRM
      ENDIF

!     delta
      WRM=WEIGHR(NTH  ,NP  ,NR  ,NSA)
      WRP=WEIGHR(NTH  ,NP  ,NR+1,NSA)
!     epsilon
      VRM=1.D0-WRM
      VRP=1.D0-WRP
      IF(NTBM.NE.0) THEN
         WRBM=WEIGHR(NTBM,NP  ,NR,NSA)
         VRBM=1.D0-WRBM
      ENDIF
      DIVDRR=1.D0/(     RL   *DELR *DELR)
      DIVFRR=1.D0/(     RL   *DELR)
!      DIVDRR=1.D0/(     RL   *DELR *DELR)*RFSADG(NR)
!      DIVFRR=1.D0/(     RL   *DELR)*RFSADG(NR)
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
                     +FRR(NTH  ,NP  ,NR,NSA)*WRM*DIVFRR*DRRM !&
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
                        +FRR(NTBM ,NP  ,NR,NSA)*WRBM*DIVFRR*DRRM !&
!                            *RLAMDAG(NTBM,NR-1)/RFSADG(NR-1)
               IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
                  LL(NM,NL)=0
                  AL(NM,NL)=0.D0
                  NL=NL-1
               ENDIF
            ENDIF
         ENDIF

!         IF(NR.NE.NRMAX.AND.NP.NE.NPMAX) THEN
         IF(NR.NE.NRMAX) THEN
            NL=NL+1
            LL(NM,NL)=NMA(NTH,NP,NR+1)
            AL(NM,NL)=DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
!                         *RLAMDAG(NTH,NR+1)/RFSADG(NR+1) &
                     -FRR(NTH  ,NP  ,NR+1,NSA)*VRP*DIVFRR*DRRP !&
!                         *RLAMDAG(NTH,NR+1)/RFSADG(NR+1)
            IF(ABS(AL(NM,NL)).LT.1.D-70) THEN
               LL(NM,NL)=0
               AL(NM,NL)=0.D0
               NL=NL-1
            ENDIF
         ENDIF

         DL(NM)=  &
              -DRR(NTH  ,NP  ,NR+1,NSA)    *DIVDRR*DRRP &
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) &
              -FRR(NTH  ,NP  ,NR+1,NSA)*WRP*DIVFRR*DRRP &
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) &
              -DRR(NTH  ,NP  ,NR  ,NSA)    *DIVDRR*DRRM &
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) &
              +FRR(NTH  ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM !&
!                            *RLAMDAG(NTH,NR)/RFSADG(NR) 
         IF(NTBM.ne.0)THEN
            DL(NM)= DL(NM) &
                 -DRR(NTBM ,NP  ,NR  ,NSA)    *DIVDRR*DRRM &
!                               *RLAMDAG(NTBM,NR)/RFSADG(NR) &
                 +FRR(NTBM ,NP  ,NR  ,NSA)*VRM*DIVFRR*DRRM !&
!                               *RLAMDAG(NTBM,NR)/RFSADG(NR)
         END IF
      ENDIF

!      SPP(NTH,NP,NR,NSA) &
!              =( SPPB(NTH,NP,NR,NSA) &
!                +SPPF(NTH,NP,NR,NSA) &
!                +SPPS(NTH,NP,NR,NSA) )

      IF(MODELD.GT.0.AND.NR.EQ.NRMAX) THEN
         SPPD(NTH,NP,NSA)= &!SPPD(NTH,NP,NSA)+ &
              FS2(NTH,NP,NSA) &
!              *RLAMDAG(NTH,NRMAX+1)/RFSADG(NRMAX+1) &
              *(DRR(NTH,NP,NR+1,NSA)    *DIVDRR &
               -FRR(NTH,NP,NR+1,NSA)*VRP*DIVFRR)*DRRP

         SPP(NTH,NP,NR,NSA) = &! SPP(NTH,NP,NR,NSA) +&
                SPPD(NTH,NP,NSA)
      ENDIF


      RETURN
      END SUBROUTINE FPSETMDR
      
      END MODULE fpdrexec
