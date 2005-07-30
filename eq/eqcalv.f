C     $Id$
C
C     ***** CALCULATE FLUX AVERAGE FOR EQCALC *****
C
      SUBROUTINE EQCALV(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      DIMENSION XA(NTVM),YA(2,NTVM)
      DIMENSION DERIV(NRVM)
C
      IERR=0
C
      CALL EQCALX(IERR)
      IF(IERR.NE.0) RETURN
C
      DR=(REDGE-RAXIS)/(NRVMAX-1)
      PSIPV(1)=0.D0
C
      DO NRV=2,NRVMAX
         RINIT=RAXIS+DR*(NRV-1)
         ZINIT=ZAXIS
         PSIL=PSIG(RINIT,ZINIT)
         PSIPL=PSIL-PSI0
         PSIPN=PSIPL/PSIPA
         PSITN=EQPSITN(PSIPN)
         PSITL=PSITA*PSITN
         RMINL=SQRT(PSITL/(BB*PI))
         IF(MDLEQA.EQ.0) THEN
            PSIN=PSIPN
         ELSE
            PSIN=PSITN
         ENDIF
         CALL EQFPSI(PSIN,FIPL,DFIPL)
         PSIPV(NRV)=PSIPL
C
         CALL EQMAGS(RINIT,ZINIT,NTVMAX,XA,YA,NA,IERR)
C
         SUMV=0.D0
         SUMQP=0.D0
         SUMAVBB2=0.D0
         SUMAVIR2=0.D0
         SUMAVRR2=0.D0
         DO N=2,NA
            CALL EQPSID(YA(1,N),YA(2,N),DPSIDR,DPSIDZ)
            R=YA(1,N)
            H=XA(N)-XA(N-1)
            BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
            BTL=FIPL/(2.D0*PI*R)
            B2L=BPL**2+BTL**2
C
            SUMV    =SUMV    +H/BPL
            SUMQP   =SUMQP   +H/(BPL*R*R)
            SUMAVBB2=SUMAVBB2+H*B2L/BPL
            SUMAVIR2=SUMAVIR2+H/(BPL*R*R)
            SUMAVRR2=SUMAVRR2+H*BPL
         ENDDO
         QPV(NRV)=FIPL*SUMQP/(4.D0*PI**2)
         AVBB2(NRV)=SUMAVBB2/(BB**2*SUMV)
         AVIR2(NRV)=SUMAVIR2*RR**2/SUMV
         AVRR2(NRV)=SUMAVRR2*QPV(NRV)**2/(RMINL**2*BB**2)
      ENDDO
C
      QPV(1)  =(4*QPV(2)  -QPV(3)  )/3.D0
      AVBB2(1)=(4*AVBB2(2)-AVBB2(3))/3.D0
      AVIR2(1)=(4*AVIR2(2)-AVIR2(3))/3.D0
      AVRR2(1)=(4*AVRR2(2)-AVRR2(3))/3.D0
C
C     ----- CALCULATE TOROIDAL FLUX -----
C
      PSITV(1)=0.D0
      DO NRV=2,NRVMAX
         PSITV(NRV)=PSITV(NRV-1)
     &        +0.5D0*(QPV(NRV)+QPV(NRV-1))*(PSIPV(NRV)-PSIPV(NRV-1))
      ENDDO
      PSITA=PSITV(NRVMAX)
C
      DO NRV=1,NRVMAX
         PSIPNV(NRV)=PSIPV(NRV)/PSIPA
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,AVBB2,DERIV,UAVBB2,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVBB2: IERR=',IERR
      CALL SPL1D(PSIPNV,AVIR2,DERIV,UAVIR2,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVIR2: IERR=',IERR
      CALL SPL1D(PSIPNV,AVRR2,DERIV,UAVRR2,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVRR2: IERR=',IERR
C
C     ----- calculate PSIITB -----
C
      CALL EQCNVA(RHOITB**2,PSIITB)
      RETURN
      END
C
C     *************************************
C        Convert Psipn to Psitn 
C     *************************************
C
      FUNCTION EQPSITN(PSIPNL)
C
      INCLUDE 'eqcomc.inc'
C
      CALL SPL1DF(PSIPNL,PSITL,PSIPNV,UPSITV,NRVMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQPSITN: SPL1DF: PSITL: IERR=',IERR
      EQPSITN=PSITL/PSITA
      RETURN
      END
C
C     *************************************
C        Calculate QPV 
C     *************************************
C
      FUNCTION EQQPV(PSIPNL)
C
      INCLUDE 'eqcomc.inc'
C
      CALL SPL1DF(PSIPNL,QPVL,PSIPNV,UQPV,NRVMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQQPV: SPL1DF: QPVL: IERR=',IERR
      EQQPV=QPVL
      RETURN
      END
C
C     **************************************
C        Calculate bounce-averaged metric
C     **************************************
C
      SUBROUTINE EQAVMT(PSIPNL,AVBB2L,AVIR2L,AVRR2L)
C
      INCLUDE 'eqcomc.inc'
C
      CALL SPL1DF(PSIPNL,AVBB2L,PSIPNV,UAVBB2,NRVMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQAVMT: SPL1DF: AVBB2: IERR=',IERR
      CALL SPL1DF(PSIPNL,AVIR2L,PSIPNV,UAVIR2,NRVMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQAVMT: SPL1DF: AVIR2: IERR=',IERR
      CALL SPL1DF(PSIPNL,AVRR2L,PSIPNV,UAVRR2,NRVMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQAVMT: SPL1DF: AVRR2: IERR=',IERR
      RETURN
      END
C
C     ****************************
C        Calculate IPOL from JP
C     ****************************
C
      SUBROUTINE EQIPJP
C
      INCLUDE 'eqcomc.inc'
      DIMENSION DERIV(NRVM)
C
      FIPV(NRVMAX)=2.D0*PI*BB*RR
C
      DO NRV=NRVMAX,2,-1
         PSIPL=PSIPV(NRV)
         PSIPNL=PSIPL/PSIPA
         CALL EQPPSI(PSIPNL,PPSI,DPPSI)
         CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
         PSIPLP=PSIPL
         ALP=RMU0*DPPSI/(BB**2*AVBB2(NRV))
         BLP=-RMU0*HJPSI/(BB*AVBB2(NRV))
C
         PSIPL=PSIPV(NRV-1)
         PSIPNL=PSIPL/PSIPA
         CALL EQPPSI(PSIPNL,PPSI,DPPSI)
         CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
         PSIPLM=PSIPL
         ALM=RMU0*DPPSI/(BB**2*AVBB2(NRV-1))
         BLM=-RMU0*HJPSI/(BB*AVBB2(NRV-1))
C
         AL=0.5D0*(ALP+ALM)*0.5D0*(PSIPLP-PSIPLM)
         BL=0.5D0*(BLP+BLM)      *(PSIPLP-PSIPLM)
         FIPV(NRV-1)=((1.D0+AL)*FIPV(NRV)-BL)/(1.D0-AL)
      ENDDO
C
C      DO NRV=1,NRVMAX
C         WRITE(6,'(A,I5,1P2E12.4)') 'NRV:',NRV,PSIPV(NRV),FIPV(NRV)
C      ENDDO
         WRITE(6,'(A,1P2E12.4)') 'FIPV:',FIPV(1),FIPV(NRVMAX)
C
      CALL SPL1D(PSIPNV,FIPV,DERIV,UFIPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVRR2: IERR=',IERR
      END
C
C     ****************************
C        Calculate IPOL from QP
C     ****************************
C
      SUBROUTINE EQIPQP
C
      INCLUDE 'eqcomc.inc'
C
      RETURN
      END
C
C     **************************************
C        Calculate bounce-averaged metric
C     **************************************
C
      SUBROUTINE EQFIPV(PSIPNL,FIPL,DFIPL)
C
      INCLUDE 'eqcomc.inc'
C
      CALL SPL1DD(PSIPNL,FIPL,DFIPL,PSIPNV,UFIPV,NRVMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQFIPV: SPL1DD: FIPL: IERR=',IERR
      DFIPL=DFIPL/PSIPA
      RETURN
      END
