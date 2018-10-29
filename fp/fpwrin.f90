!     $Id: fpwrin.f90,v 1.6 2013/01/20 23:24:03 fukuyama Exp $
!
!***********************************************************************
!     load ray data and calculate spline coefficient
!***********************************************************************
!
      MODULE fpwrin

      USE fpcomm

      INTEGER:: NRAYMAX,NITMAXM
      INTEGER,PARAMETER:: NCRMAXM=5
      INTEGER,DIMENSION(:),POINTER:: NITMAX  !(NRAYMAX)

      REAL(rkind),DIMENSION(:,:),POINTER:: RAYIN !(8,NRAYMAX)
      COMPLEX(rkind),DIMENSION(:,:),POINTER:: &  !(0:NITMAXM,NRAYMAX)
           CEXS,CEYS,CEZS
      REAL(rkind),DIMENSION(:,:),POINTER:: &  !(0:NITMAXM,NRAYMAX)
           RKXS,RKYS,RKZS,RXS,RYS,RZS,RAYRB1,RAYRB2, &
           BNXS,BNYS,BNZS,BABSS
      REAL(rkind),DIMENSION(:,:,:),POINTER:: &  !(0:8,0:NITMAXM,NRAYMAX)
           RAYS

      REAL(rkind),DIMENSION (:,:),POINTER:: PSIX,SI      !(0:NITM,NRAYM)
      COMPLEX(rkind),DIMENSION (:,:,:),POINTER:: CU1,CU2,CU3
                                                         !(4,0:NITM,NRAYM)
      REAL(rkind),DIMENSION (:,:,:),POINTER:: U4,U5,U6,U7,U8,U9,U10
                                                         !(4,0:NITM,NRAYM)

      INTEGER,DIMENSION(:,:),POINTER:: NCRMAX            !(NRM,NRAYM)
      COMPLEX(rkind),DIMENSION(:,:,:,:),POINTER:: CECR   !(3,NCRM,NRM,NRAYM)
      REAL(rkind),DIMENSION(:,:,:,:),POINTER:: RKCR,RCR  !(3,NCRM,NRM,NRAYM)

      REAL(rkind),DIMENSION(:),POINTER:: RRMIN,RRMAX     !(NRM)

      REAL(rkind),DIMENSION(:,:,:,:),POINTER:: ARGB      !(NRM,NTHM,NAVM,NRAYM)
      REAL(rkind),DIMENSION(:,:,:,:),POINTER:: RBCR      !(2,NCRM,NRM,NRAYM)
      COMPLEX(rkind),DIMENSION(:,:,:,:,:),POINTER:: CEB 
                                                       !(3,NRM,NTHM,NAVM,NRAYM)
      REAL(rkind),DIMENSION(:,:,:,:,:),POINTER:: RKB,RBB 
                                                       !(3,NRM,NTHM,NAVM,NRAYM)

      CONTAINS

!     ----- allocate most of arrays -----

      SUBROUTINE fp_wr_allocate

      IMPLICIT NONE

      ALLOCATE(RAYIN(8,NRAYMAX))
      ALLOCATE(CEXS(0:NITMAXM,NRAYMAX))
      ALLOCATE(CEYS(0:NITMAXM,NRAYMAX))
      ALLOCATE(CEZS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RKXS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RKYS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RKZS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RXS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RYS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RZS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RAYRB1(0:NITMAXM,NRAYMAX))
      ALLOCATE(RAYRB2(0:NITMAXM,NRAYMAX))
      ALLOCATE(BNXS(0:NITMAXM,NRAYMAX))
      ALLOCATE(BNYS(0:NITMAXM,NRAYMAX))
      ALLOCATE(BNZS(0:NITMAXM,NRAYMAX))
      ALLOCATE(BABSS(0:NITMAXM,NRAYMAX))
      ALLOCATE(RAYS(0:8,0:NITMAXM,NRAYMAX))

      ALLOCATE(PSIX(0:NITMAXM,NRAYMAX))
      ALLOCATE(SI(0:NITMAXM,NRAYMAX))
      ALLOCATE(CU1(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(CU2(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(CU3(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U4(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U5(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U6(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U7(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U8(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U9(4,0:NITMAXM,NRAYMAX))
      ALLOCATE(U10(4,0:NITMAXM,NRAYMAX))

      ALLOCATE(NCRMAX(NRMAX,NRAYMAX))
      ALLOCATE(CECR(3,NCRMAXM,NRMAX,NRAYMAX))
      ALLOCATE(RKCR(3,NCRMAXM,NRMAX,NRAYMAX))
      ALLOCATE(RCR(3,NCRMAXM,NRMAX,NRAYMAX))
      ALLOCATE(RRMIN(NRMAX))
      ALLOCATE(RRMAX(NRMAX))

      ALLOCATE(ARGB(NRMAX,NTHMAX+1,NAVMAX,NRAYMAX))
      ALLOCATE(RBCR(2,NCRMAXM,NRMAX,NRAYMAX))
      ALLOCATE(CEB(3,NRMAX,NTHMAX,NAVMAX,NRAYMAX))
      ALLOCATE(RKB(3,NRMAX,NTHMAX,NAVMAX,NRAYMAX))
      ALLOCATE(RBB(3,NRMAX,NTHMAX,NAVMAX,NRAYMAX))

      RETURN
      END SUBROUTINE fp_wr_allocate

!     ----- deallocate arrays -----

      SUBROUTINE fp_wr_deallocate

      IMPLICIT NONE

      IF(ASSOCIATED(NITMAX)) DEALLOCATE(NITMAX)

      IF(ASSOCIATED(RAYIN)) DEALLOCATE(RAYIN)
      IF(ASSOCIATED(CEXS)) DEALLOCATE(CEXS)
      IF(ASSOCIATED(CEYS)) DEALLOCATE(CEYS)
      IF(ASSOCIATED(CEZS)) DEALLOCATE(CEZS)
      IF(ASSOCIATED(RKXS)) DEALLOCATE(RKXS)
      IF(ASSOCIATED(RKYS)) DEALLOCATE(RKYS)
      IF(ASSOCIATED(RKZS)) DEALLOCATE(RKZS)
      IF(ASSOCIATED(RXS)) DEALLOCATE(RXS)
      IF(ASSOCIATED(RYS)) DEALLOCATE(RYS)
      IF(ASSOCIATED(RZS)) DEALLOCATE(RZS)
      IF(ASSOCIATED(RAYRB1)) DEALLOCATE(RAYRB1)
      IF(ASSOCIATED(RAYRB2)) DEALLOCATE(RAYRB2)
      IF(ASSOCIATED(BNXS)) DEALLOCATE(BNXS)
      IF(ASSOCIATED(BNYS)) DEALLOCATE(BNYS)
      IF(ASSOCIATED(BNZS)) DEALLOCATE(BNZS)
      IF(ASSOCIATED(RAYS)) DEALLOCATE(RAYS)

      IF(ASSOCIATED(PSIX))  DEALLOCATE(PSIX)
      IF(ASSOCIATED(SI))  DEALLOCATE(SI)
      IF(ASSOCIATED(CU1))  DEALLOCATE(CU1)
      IF(ASSOCIATED(CU2))  DEALLOCATE(CU2)
      IF(ASSOCIATED(CU3))  DEALLOCATE(CU3)
      IF(ASSOCIATED(U4))  DEALLOCATE(U4)
      IF(ASSOCIATED(U5))  DEALLOCATE(U5)
      IF(ASSOCIATED(U6))  DEALLOCATE(U6)
      IF(ASSOCIATED(U7))  DEALLOCATE(U7)
      IF(ASSOCIATED(U8))  DEALLOCATE(U8)
      IF(ASSOCIATED(U9))  DEALLOCATE(U9)
      IF(ASSOCIATED(U10))  DEALLOCATE(U10)

      IF(ASSOCIATED(NCRMAX))  DEALLOCATE(NCRMAX)
      IF(ASSOCIATED(CECR))  DEALLOCATE(CECR)
      IF(ASSOCIATED(RKCR))  DEALLOCATE(RKCR)
      IF(ASSOCIATED(RCR))  DEALLOCATE(RCR)
      IF(ASSOCIATED(RRMIN))  DEALLOCATE(RRMIN)
      IF(ASSOCIATED(RRMAX))  DEALLOCATE(RRMAX)

      IF(ASSOCIATED(ARGB))  DEALLOCATE(ARGB)
      IF(ASSOCIATED(RBCR))  DEALLOCATE(RBCR)
      IF(ASSOCIATED(CEB))  DEALLOCATE(CEB)
      IF(ASSOCIATED(RKB))  DEALLOCATE(RKB)
      IF(ASSOCIATED(RBB))  DEALLOCATE(RBB)

      RETURN
      END SUBROUTINE fp_wr_deallocate

      SUBROUTINE fp_wr_read(IERR)

      USE plprof
      USE libfio
      USE libmpi
      USE libmtx
      IMPLICIT NONE

      INTEGER,INtent(OUT):: IERR
      COMPLEX(rkind),DIMENSION (:),POINTER:: CFD
      REAL(rkind),DIMENSION (:),POINTER:: FD
      INTEGER,DIMENSION(2):: idata
      REAL(rkind),DIMENSION(:),POINTER:: rdata
      tYPE(pl_mag_type):: MAG
      INTEGER:: NRAY,NIT,I,NITMX,NR,NCR
      REAL(rkind):: RHON,RHOL,PSICR,PSIPRE,PSIL,SICR
      COMPLEX(rkind):: CEX,CEY,CEZ
      REAL(rkind):: RKX,RKY,RKZ,RX,RY,RZ

      CALL fp_wr_deallocate

      IF(nrank.EQ.0) THEN
         CALL FROPEN(21,KNAMWR,0,MODEFR,'WR',IERR)
         idata(1)=IERR
      ENDIF
      CALL mtx_broadcast_integer(idata,1)
      IERR=idata(1)
      IF(IERR.NE.0) THEN
         IF(nrank.EQ.0) WRITE(6,*) 'XX FPWRIN: FROPEN: IERR=',IERR
         RETURN
      ENDIF

      IF(nrank.EQ.0) THEN
         REWIND(21)
         READ(21) NRAYMAX
         ALLOCATE(NITMAX(NRAYMAX))

         NITMAXM=1
         DO NRAY=1,NRAYMAX
            READ(21) NITMAX(NRAY)
            NITMAXM=MAX(NITMAX(NRAY),NITMAXM)
         ENDDO
         idata(1)=NRAYMAX
         idata(2)=NITMAXM
      ENDIF
      CALL mtx_broadcast_integer(idata,2)
      NRAYMAX=idata(1)
      NITMAXM=idata(2)
      IF(nrank.NE.0) ALLOCATE(NITMAX(NRAYMAX))
      CALL mtx_broadcast_integer(NITMAX,NRAYMAX)

      CALL fp_wr_allocate
      ALLOCATE(rdata(NITMAXM+1))

      IF(nrank.EQ.0) THEN
         DO NRAY=1,NRAYMAX
            READ(21) (RAYIN(I,NRAY),I=1,8)
            READ(21,ERR=9101,END=9001) (CEXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21,ERR=9102,END=9002) (CEYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21,ERR=9103,END=9003) (CEZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21) (RKXS(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21) (RKYS(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21) (RKZS(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21) (RXS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
            READ(21) (RYS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
            READ(21) (RZS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
            READ(21) (RAYRB1(NIT,NRAY),NIT=0,NITMAX(NRAY))
            READ(21) (RAYRB2(NIT,NRAY),NIT=0,NITMAX(NRAY))
            DO I=0,8
               READ(21) (RAYS(I,NIT,NRAY),NIT=0,NITMAX(NRAY))
            ENDDO
            READ(21) (BNXS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
            READ(21) (BNYS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
            READ(21) (BNZS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
            READ(21) (BABSS(NIT,NRAY) ,NIT=0,NITMAX(NRAY))
         ENDDO
         CLOSE(21)
         WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
         WRITE(6,'(A,I5)') 'NRAYMAX=',NRAYMAX
      ENDIF

      DO NRAY=1,NRAYMAX
         CALL mtx_broadcast_real8(RAYIN(1:8,NRAY),8)
         CALL mtx_broadcast_complex8(CEXS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_complex8(CEYS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_complex8(CEZS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RKXS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RKYS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RKZS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RXS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RYS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RZS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RAYRB1(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(RAYRB2(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(BNXS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(BNYS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(BNZS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         CALL mtx_broadcast_real8(BABSS(0:NITMAX(NRAY),NRAY),NITMAX(NRAY)+1)
         DO I=0,8
            IF(nrank.EQ.0) THEN
               DO NIT=0,NITMAX(NRAY)
                  rdata(NIT+1)=RAYS(I,NIT,NRAY)
               ENDDO
            ENDIF
            CALL mtx_broadcast_real8(rdata,NITMAX(NRAY)+1)
            IF(nrank.NE.0) THEN
               DO NIT=0,NITMAX(NRAY)
                  RAYS(I,NIT,NRAY)=rdata(NIT+1)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      DEALLOCATE(rdata)

!      NRAY=1
!      DO NIT=NITMAX(NRAY)-10,NITMAX(NRAY)
!         CALL pl_mag(RXS(NIT,NRAY),RYS(NIT,NRAY),RZS(NIT,NRAY),RHON,MAG)
!         WRITE(6,'(A,I5,1P5E12.4)') 'NIT,X,Y,Z,B,B=',NIT, &
!              RXS(NIT,NRAY),RYS(NIT,NRAY),RZS(NIT,NRAY), &
!              MAG.BABS,BABSS(NIT,NRAY)
!      END DO
!     ----- Calculate spline coefficients -----
!
      ALLOCATE(CFD(0:NITMAXM))
      ALLOCATE(FD(0:NITMAXM))
      DO NRAY=1,NRAYMAX
         NITMX=NITMAX(NRAY)
         DO NIT=0,NITMX
            CALL pl_mag(RXS(NIT,NRAY),RYS(NIT,NRAY),RZS(NIT,NRAY),RHON,MAG)
            PSIX(NIT,NRAY)=RHON**2
            SI(NIT,NRAY)=RAYS(0,NIT,NRAY)
         ENDDO

         CALL CSPL1D(SI(0,NRAY),CEXS(0,NRAY),CFD,CU1(1,0,NRAY),NITMX+1,0,IERR)
         CALL CSPL1D(SI(0,NRAY),CEYS(0,NRAY),CFD,CU2(1,0,NRAY),NITMX+1,0,IERR)
         CALL CSPL1D(SI(0,NRAY),CEZS(0,NRAY),CFD,CU3(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RKXS(0,NRAY),FD,U4(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RKYS(0,NRAY),FD,U5(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RKZS(0,NRAY),FD,U6(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RXS(0,NRAY),FD,U7(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RYS(0,NRAY),FD,U8(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),RZS(0,NRAY),FD,U9(1,0,NRAY),NITMX+1,0,IERR)
         CALL SPL1D (SI(0,NRAY),PSIX(0,NRAY),FD,U10(1,0,NRAY),NITMX+1,0,IERR)
      ENDDO
      DEALLOCATE(cfd,fd)
!
!     ----- Find crossing point -----
!
      DO NRAY=1,NRAYMAX
         NITMX=NITMAX(NRAY)
         WRITE(6,'(A,I5,1PE12.4)') 'NRAY,RF_WR=',NRAY,RAYIN(1,NRAY)
         DO NR=1,NRMAX
            PSICR =RM(NR)**2
            PSIPRE=PSIX(0,NRAY)
            NCR=0
            DO NIT=1,NITMX
               PSIL=PSIX(NIT,NRAY)
               IF((PSIPRE-PSICR)*(PSIL-PSICR).LT.0.D0.OR. &
                   PSIL-PSICR.EQ.0.D0) THEN
                  CALL FPCROS(PSICR,NIT,NRAY,SICR)
                  WRITE(6,'(A,3I6,1P2E12.4)') '# NR,NRAY,NIT,PSICR,SICR=', &
                                                 NR,NRAY,NIT,PSICR,SICR
                  CALL FPCREK(SICR,NRAY,CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ)
                  WRITE(6,'(A,1P6E12.4)')     '# CE=',CEX,CEY,CEZ
                  WRITE(6,'(A,1P6E12.4)')     '# KR=',RKX,RKY,RKZ,RX,RY,RZ
                  
                  IF(NCR.LT.NCRMAXM) THEN
                     NCR=NCR+1
                     CECR(1,NCR,NR,NRAY)=FACT_WR*CEX
                     CECR(2,NCR,NR,NRAY)=FACT_WR*CEY
                     CECR(3,NCR,NR,NRAY)=FACT_WR*CEZ
                     RKCR(1,NCR,NR,NRAY)=RKX
                     RKCR(2,NCR,NR,NRAY)=RKY
                     RKCR(3,NCR,NR,NRAY)=RKZ
                     RCR(1,NCR,NR,NRAY)=RX
                     RCR(2,NCR,NR,NRAY)=RY
                     RCR(3,NCR,NR,NRAY)=RZ
                  ENDIF
               ENDIF
               PSIPRE=PSIL
            ENDDO
            NCRMAX(NR,NRAY)=NCR
         ENDDO
      ENDDO

      DO NR=1,NRMAX
         RHOL =RM(NR)
!         CALL pl_rrmx(RHOL,RRMIN(NR),RRMAX(NR))
         CALL pl_rrminmax(RHOL,RRMIN(NR),RRMAX(NR)) 
!         WRITE(6,602) NR,RHOL,RRMIN(NR),RRMAX(NR)
!  602    FORMAT('# NR,RHOL,RRMIN,RRMAX=',I3,1P3E15.7)
      ENDDO

  900 IERR=0
      RETURN
9001  write(6,*) '9001:NIT=',NIT
      stop
9002  write(6,*) '9002:NIT=',NIT
      stop
9003  write(6,*) '9003:NIT=',NIT
      stop
9101  write(6,*) '9101:NIT=',NIT
      stop
9102  write(6,*) '9102:NIT=',NIT
      stop
9103  write(6,*) '9103:NIT=',NIT
      stop
      RETURN
      END SUBROUTINE fp_wr_read
!
!***********************************************************************
!     Calculate crossing point by Newton's method
!***********************************************************************
!
      SUBROUTINE FPCROS(PSICR,NIT,NRAY,SICR)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: PSICR
      INTEGER,INTENT(IN):: NIT,NRAY
      REAL(rkind),INTENT(OUT):: SICR
      REAL(rkind):: Y,YNEW,X,DX,DYDX,PSIL
      INTEGER:: ICOUNT,IERR,NITMX

      Y   =PSIX(NIT-1,NRAY)-PSICR
      YNEW=PSIX(NIT  ,NRAY)-PSICR
      X   =SI(NIT-1,NRAY)
      DX  =SI(NIT,NRAY)-SI(NIT-1,NRAY)
      IF(DX.EQ.0.D0) GOTO 8000
      DYDX=(YNEW-Y)/DX
      Y   =YNEW

      ICOUNT=0

  100 CONTINUE
         ICOUNT=ICOUNT+1

!         WRITE(6,'(A,I3,1P4E15.7)') 'FPCROS:',ICOUNT,X,DX,Y,YNEW
         IF(ABS(Y).LE.EPS_WR) GOTO 7000
         IF(ICOUNT.GT.LMAX_WR) GOTO 8100
         IF(DYDX.EQ.0.D0) GOTO 8200
         DX=-Y/DYDX
         IF(DX.EQ.0.D0) GOTO 8300
         X=X+DX
         CALL SPL1DF(X,PSIL,SI(0,NRAY),U10(1,0,NRAY),NITMAX(NRAY),IERR)
         YNEW=PSIL-PSICR
         DYDX=(YNEW-Y)/DX
         Y=YNEW
      GOTO 100

 7000 SICR=X
      RETURN

 8000 WRITE(6,*) 'XX FPCROS: INVALID INITIAL VALUES: X0=X1'
      GOTO 7000
 8100 WRITE(6,*) 'XX FPCROS: ICOUNT EXCEEDS LMAX_WR AT SI =',X
      GOTO 7000
 8200 WRITE(6,*) 'XX FPCROS: DYDX BECOMES 0 AT SI =',X
      GOTO 7000
 8300 WRITE(6,*) 'XX FPCROS: DX BECOMES 0 AT SI = ',X
      GOTO 7000
      END SUBROUTINE FPCROS
!
!***********************************************************************
!     Calculate wave field, wave number and position at crossing point
!***********************************************************************
!
      SUBROUTINE FPCREK(SICR,NRAY,CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: SICR
      INTEGER,INTENT(IN):: NRAY
      COMPLEX(rkind),INTENT(OUT):: CEX,CEY,CEZ
      REAL(rkind),INTENT(OUT):: RKX,RKY,RKZ,RX,RY,RZ
      INTEGER:: NITMX,IERR

      NITMX=NITMAX(NRAY)
      CALL CSPL1DF(SICR,CEX,SI(0,NRAY),CU1(1,0,NRAY),NITMX,IERR)
      CALL CSPL1DF(SICR,CEY,SI(0,NRAY),CU2(1,0,NRAY),NITMX,IERR)
      CALL CSPL1DF(SICR,CEZ,SI(0,NRAY),CU3(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RKX,SI(0,NRAY),U4(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RKY,SI(0,NRAY),U5(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RKZ,SI(0,NRAY),U6(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RX,SI(0,NRAY),U7(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RY,SI(0,NRAY),U8(1,0,NRAY),NITMX,IERR)
      CALL SPL1DF(SICR,RZ,SI(0,NRAY),U9(1,0,NRAY),NITMX,IERR)

      RETURN
      END SUBROUTINE FPCREK

      END MODULE fpwrin
