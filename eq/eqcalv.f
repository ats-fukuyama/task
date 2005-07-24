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
         PSIPN=-PSIPL/PSI0
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
      FUNCTION EQPSITN(PSIN)
C
      INCLUDE 'eqcomc.inc'
C
      CALL SPL1DF(PSIN,PSITL,PSIPNV,UPSITV,NRVMAX,IERR)
      EQPSITN=PSITL/PSITA
C      WRITE(6,'(A,1P6E12.4)') 'EQPSITN:',PSIN,PSIPL,PSIPA,
C     &              EQPSITN,PSITL,PSITA
      RETURN
      END
C
C     *************************************
C        Calculate QPV 
C     *************************************
C
      FUNCTION EQQPV(PSIN)
C
      INCLUDE 'eqcomc.inc'
C
      CALL SPL1DF(PSIN,QPVL,PSIPNV,UQPV,NRVMAX,IERR)
      EQQPV=QPVL
      RETURN
      END
