C     $Id$
C
C     ****** FREE BOUNDARY EQUILIBRIUM SOLVER ******
C
      SUBROUTINE EQCALX(ID,IERR)
C
      INCLUDE '../eq/eqcomx.inc'
      DIMENSION PSIBL(4),RSLT(4),PSIBLX(4)
      EXTERNAL EQCALX_SUB,EQCALX_LOOP
C
      IERR=0
C
      DRG=(RGMAX-RGMIN)/(NRGMAX-1)
      DO NRG=1,NRGMAX
         RG(NRG)=RGMIN+DRG*(NRG-1)
      ENDDO
C
      DZG=(ZGMAX-ZGMIN)/(NZGMAX-1)
      DO NZG=1,NZGMAX
         ZG(NZG)=ZGMIN+DZG*(NZG-1)
      ENDDO
C
      MWMAX=8*(NRGMAX+2)-1
      MLMAX=4*NRGMAX*NZGMAX
      NBND=4*(NRGMAX+2)
C
      IF(ID.EQ.0) CALL EQPSIX_INIT
      CALL EQCALFMA
C
      DO I=1,4
         PSIBL(I)=PSIB(I-1)
      ENDDO
C
      IF(MDLEQX.EQ.0) THEN
C
         CALL EQCALX_LOOP(PSIBL,RSLT,4,IERR)
C
      ELSEIF(MDLEQX.EQ.1) THEN
C
         CALL NEWTON4(4,PSIBL,PSIBLX,EQCALX_LOOP,IERR,ICOUNT,
     &                DELNW,EPSNW,NLPNW)
C
      ELSEIF(MDLEQX.EQ.2) THEN
C
         CALL NEWTON4(4,PSIBL,PSIBLX,EQCALX_SUB,IERR,ICOUNT,
     &                DELNW,EPSNW,NLPNW)
C
      ELSE
         WRITE(6,*) 'XX EQCALX: UNKNOWN MDLEQX:', MDLEQX
      ENDIF
C
      DPS=1.D0/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPNL=DPS*(NPS-1)
         PSIPS(NPS)=PSIPA*PSIPNL
         CALL EQPPSI(PSIPNL,PPSI,DPPSI)
         CALL EQFPSI(PSIPNL,FPSI,DFPSI)
         PPPS(NPS)=PPSI
         TTPS(NPS)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
C         WRITE(6,'(I5,1P5E12.4)')
C     &        NPS,PSIPNL,PPSI,FPSI,PPPS(NPS),TTPS(NPS)
      ENDDO
C
C      CALL EQGRAX
C
      RETURN
      END
C
C     ****** SOLVER LOOP ******
C
      SUBROUTINE EQCALX_LOOP(PSIBL,RSLT,NEQ,IERR)
C
      INCLUDE '../eq/eqcomx.inc'
      DIMENSION PSIBL(NEQ),RSLT(NEQ)
C
      DO NLOOP=1,NLPMAX
C
         CALL EQCALX_SUB(PSIBL,RSLT,NEQ,IERR)
         IF(IERR.NE.0) RETURN
C
         SUM0=0.D0
         SUM1=0.D0
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            SUM0=SUM0+PSIX(NRG,NZG)**2
            SUM1=SUM1+DELPSIX(NRG,NZG)**2
         ENDDO
         ENDDO
         SUM2=SQRT(SUM1/SUM0)
C
         IF(SUM2.LT.EPSEQ) THEN
            IF(NPRINT.GE.1) THEN
               WRITE(6,'(A,1P4E14.6,I5)')
     &           'SUM2,R/ZAXIS,PSI0=',SUM2,RAXIS,ZAXIS,PSI0,NLOOP
            ENDIF
            GOTO 1000
         ELSE
            IF((NPRINT.EQ.1.AND.NLOOP.EQ.1).OR.
     &          NPRINT.GE.2) THEN
               WRITE(6,'(A,1P4E14.6)')
     &           'SUM2,R/ZAXIS,PSI0=',SUM2,RAXIS,ZAXIS,PSI0
            ENDIF
         ENDIF
      ENDDO
C
      IF(NPRINT.GE.1) THEN
         WRITE(6,'(A,1P4E14.6,I5)')
     &           'SUM2,R/ZAXIS,PSI0=',SUM2,RAXIS,ZAXIS,PSI0,NLOOP
      ENDIF
      WRITE(6,*) 'XX EQCALX: NLOOP exceeds NLPMAX'
      IERR=10
C
 1000 CONTINUE
C
      RETURN
      END
C
C     ****** CORE SOLVER ******
C
      SUBROUTINE EQCALX_SUB(PSIBL,RSLT,NEQ,IERR)
C
      USE libbnd
      INCLUDE '../eq/eqcomx.inc'
      DIMENSION PSIBL(NEQ),RSLT(NEQ)
      DIMENSION FMAX(MWM,MLM)
C
      IERR=0
C
      DO I=1,NEQ
         PSIB(I-1)=PSIBL(I)
      ENDDO
      WRITE(6,'(4(A,1PE12.4))') 'PSIB=',PSIB(0),'      =',PSIB(1),
     &                        '      =',PSIB(2),'      =',PSIB(3)
C
      DO ML=1,MLMAX
      DO MW=1,MWMAX
         FMAX(MW,ML)=FMA(MW,ML)
      ENDDO
      ENDDO
C
      CALL EQCALFVB(IERR)
      IF(IERR.NE.0) RETURN
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         PSIXOLD(NRG,NZG)=PSIX(NRG,NZG)
      ENDDO
      ENDDO
C
      CALL BANDRD(FMAX,FVB,MLMAX,MWMAX,MWM,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'BANDRD ERROR: IERR =',IERR
         RETURN
      ENDIF
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         L=4*((NZG-1)*NRGMAX+NRG-1)+1
         PSIX(NRG,NZG)=FVB(L)
         UPSIX(1,NRG,NZG)=FVB(L)
         UPSIX(2,NRG,NZG)=FVB(L+1)
         UPSIX(3,NRG,NZG)=FVB(L+2)
         UPSIX(4,NRG,NZG)=FVB(L+3)
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         DELPSIX(NRG,NZG)=PSIX(NRG,NZG)-PSIXOLD(NRG,NZG)
         HJTRZ(NRG,NZG)=FJRZ(NRG,NZG)
      ENDDO
      ENDDO
C
      CALL EQAXIS(IERR)
C
      IF(PSI0*RIP.GE.0.D0) THEN
         WRITE(6,*) 
     &        'XX EQCALX_SUB: Invaild sign of PSI0: no plasma region'
         WRITE(6,'(A,1P4E12.4)') 'RIP,PSI0,AXIS=',RIP,PSI0,RAXIS,ZAXIS
         IERR=100
      ENDIF
C
      IF(IDEBUG.EQ.1) THEN
         CALL EQGX1D
         CALL EQGX2D
      ENDIF
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         PSIRZ(NRG,NZG)=PSIX(NRG,NZG)-PSI0
         HJTRZ(NRG,NZG)=FJRZ(NRG,NZG)
      ENDDO
      ENDDO
C
      CALL EQCALR(RRR,RRA,RRK,RRD)
C
      RSLT(1)=RRR-RR
      RSLT(2)=RRA-RA
      RSLT(3)=RRK-RKAP
      RSLT(4)=RRD-RDLT
C
      RETURN
      END
C
C     ******************************************
C     *     Set FEM element integral table     *
C     ******************************************
C
      SUBROUTINE EQTABL
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION K1AL(4,4),M1AL(4,4)
C
      DATA K1AL/1,2,1,2,
     &          3,4,3,4,
     &          1,2,1,2,
     &          3,4,3,4/
      DATA M1AL/1,1,2,2,
     &          1,1,2,2,
     &          3,3,4,4,
     &          3,3,4,4/
C
      DO J=1,4
      DO I=1,4
         K1A(I,J)=K1AL(I,J)
         M1A(I,J)=M1AL(I,J)
      ENDDO
      ENDDO
C
      RK(1,1,1,1)= 3.D0/5.D0
      RK(1,1,1,2)= 3.D0/5.D0
      RK(1,1,2,1)=-3.D0/5.D0
      RK(1,1,2,2)=-3.D0/5.D0
      RK(1,2,1,1)=-3.D0/5.D0
      RK(1,2,1,2)=-3.D0/5.D0
      RK(1,2,2,1)= 3.D0/5.D0
      RK(1,2,2,2)= 3.D0/5.D0
c
      RK(2,1,1,1)=-1.D0/70.D0
      RK(2,1,1,2)= 4.D0/35.D0
      RK(2,1,2,1)= 1.D0/70.D0
      RK(2,1,2,2)=-4.D0/35.D0
      RK(2,2,1,1)= 4.D0/35.D0
      RK(2,2,1,2)=-1.D0/70.D0
      RK(2,2,2,1)=-4.D0/35.D0
      RK(2,2,2,2)= 1.D0/70.D0
c
      RK(3,1,1,1)=-1.D0/70.D0
      RK(3,1,1,2)= 4.D0/35.D0
      RK(3,1,2,1)= 4.D0/35.D0
      RK(3,1,2,2)=-1.D0/70.D0
      RK(3,2,1,1)= 1.D0/70.D0
      RK(3,2,1,2)=-4.D0/35.D0
      RK(3,2,2,1)=-4.D0/35.D0
      RK(3,2,2,2)= 1.D0/70.D0
c
      RK(4,1,1,1)= 43.D0/420.D0
      RK(4,1,1,2)= 13.D0/420.D0
      RK(4,1,2,1)=-1.D0/60.D0
      RK(4,1,2,2)=-1.D0/60.D0
      RK(4,2,1,1)=-1.D0/60.D0
      RK(4,2,1,2)=-1.D0/60.D0
      RK(4,2,2,1)= 13.D0/420.D0
      RK(4,2,2,2)= 43.D0/420.D0
C
      RKK(1,1,1,1)= 9.D0/70.D0
      RKK(1,1,1,2)=-9.D0/70.D0
      RKK(1,1,2,1)=-9.D0/70.D0
      RKK(1,1,2,2)= 9.D0/70.D0
      RKK(1,2,1,1)=-9.D0/70.D0
      RKK(1,2,1,2)= 9.D0/70.D0
      RKK(1,2,2,1)= 9.D0/70.D0
      RKK(1,2,2,2)=-9.D0/70.D0
c
      RKK(2,1,1,1)= 1.D0/140.D0
      RKK(2,1,1,2)=-3.D0/140.D0
      RKK(2,1,2,1)=-1.D0/140.D0
      RKK(2,1,2,2)= 3.D0/140.D0
      RKK(2,2,1,1)= 3.D0/140.D0
      RKK(2,2,1,2)=-1.D0/140.D0
      RKK(2,2,2,1)=-3.D0/140.D0
      RKK(2,2,2,2)= 1.D0/140.D0
c
      RKK(3,1,1,1)= 1.D0/140.D0
      RKK(3,1,1,2)=-3.D0/140.D0
      RKK(3,1,2,1)= 3.D0/140.D0
      RKK(3,1,2,2)=-1.D0/140.D0
      RKK(3,2,1,1)=-1.D0/140.D0
      RKK(3,2,1,2)= 3.D0/140.D0
      RKK(3,2,2,1)=-3.D0/140.D0
      RKK(3,2,2,2)= 1.D0/140.D0
c
      RKK(4,1,1,1)= 1.D0/120.D0
      RKK(4,1,1,2)=-1.D0/168.D0
      RKK(4,1,2,1)=-1.D0/840.D0
      RKK(4,1,2,2)= 1.D0/840.D0
      RKK(4,2,1,1)=-1.D0/840.D0
      RKK(4,2,1,2)= 1.D0/840.D0
      RKK(4,2,2,1)= 1.D0/168.D0
      RKK(4,2,2,2)=-1.D0/120.D0
C
      RL(1,1,1,1)= 43.D0/140.D0
      RL(1,1,1,2)=  9.D0/140.D0
      RL(1,1,2,1)=  9.D0/140.D0
      RL(1,1,2,2)=  9.D0/140.D0
      RL(1,2,1,1)=  9.D0/140.D0
      RL(1,2,1,2)=  9.D0/140.D0
      RL(1,2,2,1)=  9.D0/140.D0
      RL(1,2,2,2)= 43.D0/140.D0
c
      RL(2,1,1,1)= 97.D0/2520.D0
      RL(2,1,1,2)=  1.D0/72.D0
      RL(2,1,2,1)=  1.D0/72.D0
      RL(2,1,2,2)= 43.D0/2520.D0
      RL(2,2,1,1)=-43.D0/2520.D0
      RL(2,2,1,2)=- 1.D0/72.D0
      RL(2,2,2,1)=- 1.D0/72.D0
      RL(2,2,2,2)=-97.D0/2520.D0
c
      RL(3,1,1,1)= 97.D0/2520.D0
      RL(3,1,1,2)=  1.D0/72.D0
      RL(3,1,2,1)=-43.D0/2520.D0
      RL(3,1,2,2)=- 1.D0/72.D0
      RL(3,2,1,1)=  1.D0/72.D0
      RL(3,2,1,2)= 43.D0/2520.D0
      RL(3,2,2,1)=- 1.D0/72.D0
      RL(3,2,2,2)=-97.D0/2520.D0
C
      RL(4,1,1,1)= 2.D0/315.D0
      RL(4,1,1,2)= 1.D0/315.D0
      RL(4,1,2,1)=-1.D0/280.D0
      RL(4,1,2,2)=-1.D0/280.D0
      RL(4,2,1,1)=-1.D0/280.D0
      RL(4,2,1,2)=-1.D0/280.D0
      RL(4,2,2,1)= 1.D0/315.D0
      RL(4,2,2,2)= 2.D0/315.D0
C
      RLL(1,1,1,1)= 97.D0/2520.D0
      RLL(1,1,1,2)=-43.D0/2520.D0
      RLL(1,1,2,1)=  1.D0/72.D0
      RLL(1,1,2,2)=- 1.D0/72.D0
      RLL(1,2,1,1)=  1.D0/72.D0
      RLL(1,2,1,2)=- 1.D0/72.D0
      RLL(1,2,2,1)= 43.D0/2520.D0
      RLL(1,2,2,2)=-97.D0/2520.D0
C
      RLL(2,1,1,1)= 2.D0/315.D0
      RLL(2,1,1,2)=-1.D0/280.D0
      RLL(2,1,2,1)= 1.D0/315.D0
      RLL(2,1,2,2)=-1.D0/280.D0
      RLL(2,2,1,1)=-1.D0/280.D0
      RLL(2,2,1,2)= 1.D0/315.D0
      RLL(2,2,2,1)=-1.D0/280.D0
      RLL(2,2,2,2)= 2.D0/315.D0
C
      RLL(3,1,1,1)= 2.D0/315.D0
      RLL(3,1,1,2)=-1.D0/280.D0
      RLL(3,1,2,1)=-1.D0/280.D0
      RLL(3,1,2,2)= 1.D0/315.D0
      RLL(3,2,1,1)= 1.D0/315.D0
      RLL(3,2,1,2)=-1.D0/280.D0
      RLL(3,2,2,1)=-1.D0/280.D0
      RLL(3,2,2,2)= 2.D0/315.D0
C
      RLL(4,1,1,1)= 1.D0/840.D0
      RLL(4,1,1,2)=-1.D0/1260.D0
      RLL(4,1,2,1)=-1.D0/1260.D0
      RLL(4,1,2,2)= 1.D0/1260.D0
      RLL(4,2,1,1)=-1.D0/1260.D0
      RLL(4,2,1,2)= 1.D0/1260.D0
      RLL(4,2,2,1)= 1.D0/1260.D0
      RLL(4,2,2,2)=-1.D0/840.D0
C
      RM(1,1,1)= 13.D0/ 35.D0
      RM(1,1,2)=  9.D0/ 70.D0
      RM(1,2,1)=  9.D0/ 70.D0
      RM(1,2,2)= 13.D0/ 35.D0
      RM(2,1,1)= 11.D0/210.D0
      RM(2,1,2)= 13.D0/420.D0
      RM(2,2,1)=-13.D0/420.D0
      RM(2,2,2)=-11.D0/210.D0
      RM(3,1,1)= 11.D0/210.D0
      RM(3,1,2)=-13.D0/420.D0
      RM(3,2,1)= 13.D0/420.D0
      RM(3,2,2)=-11.D0/210.D0
      RM(4,1,1)=  1.D0/105.D0
      RM(4,1,2)=- 1.D0/140.D0
      RM(4,2,1)=- 1.D0/140.D0
      RM(4,2,2)=  1.D0/105.D0
C
      RN(1,1,1)= 6.D0/ 5.D0
      RN(1,1,2)=-6.D0/ 5.D0
      RN(1,2,1)=-6.D0/ 5.D0
      RN(1,2,2)= 6.D0/ 5.D0
      RN(2,1,1)= 1.D0/10.D0
      RN(2,1,2)=-1.D0/10.D0
      RN(2,2,1)= 1.D0/10.D0
      RN(2,2,2)=-1.D0/10.D0
      RN(3,1,1)= 1.D0/10.D0
      RN(3,1,2)= 1.D0/10.D0
      RN(3,2,1)=-1.D0/10.D0
      RN(3,2,2)=-1.D0/10.D0
      RN(4,1,1)= 2.D0/15.D0
      RN(4,1,2)=-1.D0/30.D0
      RN(4,2,1)=-1.D0/30.D0
      RN(4,2,2)= 2.D0/15.D0
C
      RETURN
      END
C
C   ************************************************ 
C   **              Initial Psi                 **
C   ************************************************
C
      SUBROUTINE EQPSIX_INIT
C
      USE libspl1d
      INCLUDE '../eq/eqcomx.inc'
      DIMENSION DERIV(NRVM)
C
      RAXIS=RR
      ZAXIS=0.0D0
C
C     --- assuming elliptic crosssection, flat current profile ---
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MOD(MDLEQFL,5).EQ.4) THEN
         PSITA=PI*RKAP*RA**2*BB
         PSIPA=PSITA/QA
         PSI0=-PSIPA
      ELSE
         PSI0=-0.5D0*RMU0*RIP*1.D6*RR
         PSIPA=-PSI0
         PSITA=PI*RKAP*RA**2*BB
      ENDIF
C
      DRHO=1.D0/(NRVMAX-1)
      DO NRV=1,NRVMAX
         RHOL=(NRV-1)*DRHO
         PSIPNV(NRV)=RHOL*RHOL
         PSIPV(NRV)=PSIPA*RHOL*RHOL
         PSITV(NRV)=PSITA*RHOL*RHOL
         QPV(NRV)=PSITA/PSIPA
         TTV(NRV)=2.D0*PI*RR*BB
      ENDDO
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         RHON2=((RG(NRG)-RR)**2+(ZG(NZG)/RKAP)**2)/RA**2
         PSIX(NRG,NZG)=PSI0*(1-RHON2)
         HJTRZ(NRG,NZG)=0.D0
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         UPSIX(1,NRG,NZG)=PSIX(NRG,NZG)
         IF(NRG.EQ.1) THEN
            NRGP=NRG+1
            NRGN=NRG
         ELSEIF(NRG.EQ.NRGMAX) THEN
            NRGP=NRG
            NRGN=NRG-1
         ELSE
            NRGP=NRG+1
            NRGN=NRG-1
         ENDIF
         IF(NZG.EQ.1) THEN
            NZGP=NZG+1
            NZGN=NZG
         ELSEIF(NZG.EQ.NZGMAX) THEN
            NZGP=NZG
            NZGN=NZG-1
         ELSE
            NZGP=NZG+1
            NZGN=NZG-1
         ENDIF
         UPSIX(2,NRG,NZG)=(PSIX(NRGP,NZG)-PSIX(NRGN,NZG))
     &                   /(RG(NRGP)-RG(NRGN))
         UPSIX(3,NRG,NZG)=(PSIX(NRG,NZGP)-PSIX(NRG,NZGN))
     &                   /(ZG(NZGP)-ZG(NZGN))
         UPSIX(4,NRG,NZG)=(PSIX(NRGP,NZGP)-PSIX(NRGP,NZGN)
     &                    -PSIX(NRGN,NZGP)+PSIX(NRGN,NZGN))
     &                   /(RG(NRGP)-RG(NRGN))/(ZG(NZGP)-ZG(NZGN))
      ENDDO
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
C      WRITE(6,'(A,I5,1P4E12.4)')
C     &     ('NRV:',NRV,PSIPNV(NRV),PSITV(NRV),QPV(NRV),TTV(NRV),
C     &      NRV=1,NRVMAX)
C
      RETURN
      END
C
C     ********************************************
C     *     Calculate FEM coefficient matrix     *
C     ********************************************
C
      SUBROUTINE EQCALFMA
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION FACTK(4),FACTL(4),FACTM(4),FACTN(4)
      DIMENSION RG_SAVE(NRGM),ZG_SAVE(NRGM)
      SAVE NRGMAX_SAVE,NZGMAX_SAVE,RG_SAVE,ZG_SAVE,INIT
      DATA INIT/0/
C
      IF(INIT.EQ.1) THEN
         IF(NRGMAX.NE.NRGMAX_SAVE) INIT=0
         IF(NZGMAX.NE.NZGMAX_SAVE) INIT=0
      ENDIF
      IF(INIT.EQ.1) THEN
         DO NRG=1,NRGMAX
            IF(RG(NRG).NE.RG_SAVE(NRG)) INIT=0
         ENDDO
         DO NZG=1,NZGMAX
            IF(ZG(NZG).NE.ZG_SAVE(NZG)) INIT=0
         ENDDO
      ENDIF
      IF(INIT.EQ.1) RETURN
C
      IF(INIT.EQ.0) THEN
         NRGMAX_SAVE=NRGMAX
         NZGMAX_SAVE=NZGMAX
         DO NRG=1,NRGMAX
            RG_SAVE(NRG)=RG(NRG)
         ENDDO
         DO NZG=1,NZGMAX
            ZG_SAVE(NZG)=ZG(NZG)
         ENDDO
         INIT=1
      ENDIF
C
      DO N=1,MLMAX
      DO M=1,MWMAX
         FMA(M,N)=0
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX-1
      DO NRG=1,NRGMAX-1
         DRG=RG(NRG+1)-RG(NRG)
         FACTK(1)=1.D0/DRG
         FACTK(2)=1.D0
         FACTK(3)=1.D0
         FACTK(4)=DRG
         FACTL(1)=DRG
         FACTL(2)=DRG**2
         FACTL(3)=DRG**2
         FACTL(4)=DRG**3
C
         DZG=ZG(NZG+1)-ZG(NZG)
         FACTM(1)=DZG
         FACTM(2)=DZG**2
         FACTM(3)=DZG**2
         FACTM(4)=DZG**3
         FACTN(1)=1.D0/DZG
         FACTN(2)=1.D0
         FACTN(3)=1.D0
         FACTN(4)=DZG
C
         DO K=1,4
         DO L=1,4
            K1=K1A(K,L)
            M1=M1A(K,L)
            L1=K1A(K,L)
            N1=M1A(K,L)
         DO I1=1,2
         DO I2=1,2
         DO J1=1,2
         DO J2=1,2
            N=4*((NZG-1)*NRGMAX+NRG-1)+K+4*(I1-1)+4*NRGMAX*(J1-1)
            M=4*((NZG-1)*NRGMAX+NRG-1)+L+4*(I2-1)+4*NRGMAX*(J2-1)
         DO I3=1,2
            FMA(NBND+M-N,N)=FMA(NBND+M-N,N)
     &                     -(RK(K1,I1,I2,I3)/RG(NRG+I3-1)
     &                      -RKK(K1,I1,I2,I3)*DRG/RG(NRG+I3-1)**2)
     &                                   *FACTK(K1)
     &                      *RM(M1,J1,J2)*FACTM(M1)
     &                     -(RL(L1,I1,I2,I3)/RG(NRG+I3-1)
     &                      -RLL(L1,I1,I2,I3)*DRG/RG(NRG+I3-1)**2)
     &                                   *FACTL(L1)
     &                      *RN(N1,J1,J2)*FACTN(N1) 
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- SET BOUNDARY ARRAY -----
C
      NBA=0
      DO NRG=1,NRGMAX-1
         NBA=NBA+1
         NRGB(NBA)=NRG
         NZGB(NBA)=1
      ENDDO
      DO NZG=1,NZGMAX-1
         NBA=NBA+1
         NRGB(NBA)=NRGMAX
         NZGB(NBA)=NZG
      ENDDO
      DO NRG=NRGMAX,2,-1
         NBA=NBA+1
         NRGB(NBA)=NRG
         NZGB(NBA)=NZGMAX
      ENDDO
      DO NZG=NZGMAX,2,-1
         NBA=NBA+1
         NRGB(NBA)=1
         NZGB(NBA)=NZG
      ENDDO
      NBAMAX=NBA
C
C     ----- SET BOUNDARY CONDITION AND PLASMA CURRENT FACTOR -----
C
      DR=1.D-6
      DZ=1.D-6
      DO NBA=1,NBAMAX
         NRG=NRGB(NBA)
         NZG=NZGB(NBA)
         IF(NRG.EQ.1.OR.NRG.EQ.NRGMAX) THEN
         DO K=1,3,2
            N=4*((NZG-1)*NRGMAX+NRG-1)+K
            DO MW=1,MWMAX
               FMA(MW,N)=0.D0
            ENDDO
            FMA(NBND,N)=1.D0
         ENDDO
         ELSE
         DO K=1,2
            N=4*((NZG-1)*NRGMAX+NRG-1)+K
            DO MW=1,MWMAX
               FMA(MW,N)=0.D0
            ENDDO
            FMA(NBND,N)=1.D0
         ENDDO
         ENDIF
C
         R1=RG(NRG)
         Z1=ZG(NZG)
         DO NZG2=2,NZGMAX-1
         DO NRG2=2,NRGMAX-1
            R2=RG(NRG2)
            Z2=ZG(NZG2)
            DRG=0.5D0*(RG(NRG2+1)-RG(NRG2-1))
            DZG=0.5D0*(ZG(NZG2+1)-ZG(NZG2-1))
            F00=RGSGRF(R1,Z1,R2,Z2)
            FPP=RGSGRF(R1+DR,Z1+DZ,R2,Z2)
            FMP=RGSGRF(R1-DR,Z1+DZ,R2,Z2)
            FPM=RGSGRF(R1+DR,Z1-DZ,R2,Z2)
            FMM=RGSGRF(R1-DR,Z1-DZ,R2,Z2)
            FDR=((FPP+FPM)-(FMP+FMM))/(4.D0*DR)
            FDZ=((FPP+FMP)-(FPM+FMM))/(4.D0*DZ)
            FRZ=((FPP+FMM)-(FPM+FMP))/(4.D0*DR*DZ)
            PSIBGR(NRG2,NZG2,NBA,1)=DRG*DZG*F00
            PSIBGR(NRG2,NZG2,NBA,2)=DRG*DZG*FDR
            PSIBGR(NRG2,NZG2,NBA,3)=DRG*DZG*FDZ
            PSIBGR(NRG2,NZG2,NBA,4)=DRG*DZG*FRZ
         ENDDO
         ENDDO
      ENDDO
C
C      OPEN(16,FILE='eqxdata',FORM='FORMATTED')
C      DO 9999 ML=1,MLMAX
C         WRITE(16,*) 'ML=',ML
C         WRITE(16,'(10F7.2)') (FMA(MW,ML),MW=1,MWMAX)
C 9999 CONTINUE
C      CLOSE(16)
C
      RETURN
      END
C
C     **********************************
C     *     Calulate FEM RHS vector    *
C     **********************************
C
      SUBROUTINE EQCALFVB(IERR)
C
      USE libspl1d
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION FACTM(4),FACTH(4)
      DIMENSION RJ(4),PSIBRZ(4),SUML(4)
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: HJ0,HJ1,HJ2

      ALLOCATE(HJ0(NRGM,NZGM),HJ1(NRGM,NZGM),HJ2(NRGM,NZGM))

      IMDLEQF=MOD(MDLEQF,5)
C
C     ----- Plasma current -----
C
      IF(IMDLEQF.EQ.0) THEN
         RRC=RAXIS
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=(1.D0-RRC**2/RG(NRG)**2)*HJ0(NRG,NZG)
               HJ2(NRG,NZG)=-(RRC/RG(NRG))*HJPSI
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               HJ2(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C        ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         FJ2=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
            FJ2=FJ2+HJ2(NRG,NZG)*DVOL
         ENDDO
         ENDDO
         TJ=(RIP*1.D6-FJ1)/FJ2
         RIPX=RIP
C
         DO NRG=1,NRGMAX
         DO NZG=1,NZGMAX
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ1(NRG,NZG)+TJ*HJ2(NRG,NZG))
            FJRZ(NRG,NZG)=FJPRZ(NRG,NZG)
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
               TTRZ(NRG,NZG)=SQRT((2.D0*PI*BB*RR)**2
     &                           +2.D0*RMU0*RRC
     &                           *(2.D0*PI*TJ*HJPSID-RRC*PPSI))
               PPRZ(NRG,NZG)=PPSI
            ELSE
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C     ----- Given pressure and poloidal current profiles -----
C
      ELSEIF(IMDLEQF.EQ.1) THEN

         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=-2.D0*PI*BB*RR*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
               HJ2(NRG,NZG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               HJ2(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C        ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         FJ2=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
            FJ2=FJ2+HJ2(NRG,NZG)*DVOL
         ENDDO
         ENDDO
C
C        ----- Adjust plasma current -----
C
         IF(FJ1.GT.0.D0) THEN
            TJ=(-FJ1+SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ELSE
            TJ=(-FJ1-SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ENDIF
         RIPX=RIP
         IERR=0
C
C         WRITE(6,'(A,1P4E12.4)') 'FJ0,FJ1,FJ2,TJ=',FJ0,FJ1,FJ2,TJ
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)
     &                                  +TJ*HJ1(NRG,NZG)
     &                                  +TJ*TJ*HJ2(NRG,NZG))
            FJRZ(NRG,NZG)=FJPRZ(NRG,NZG)
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               PPRZ(NRG,NZG)=PPSI
            ELSE
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C     ----- Given pressure and parallel current profiles -----
C
      ELSEIF(IMDLEQF.EQ.2) THEN
         CALL EQIPJP
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFIPV(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=-2.D0*PI*BB*RR       *DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
               HJ2(NRG,NZG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               HJ2(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C        ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         FJ2=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
            FJ2=FJ2+HJ2(NRG,NZG)*DVOL
         ENDDO
         ENDDO
C
         IF(FJ1.GT.0.D0) THEN
            TJ=(-FJ1+SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ELSE
            TJ=(-FJ1-SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ENDIF
         RIPX=RIP
C
         DO NRG=1,NRGMAX
         DO NZG=1,NZGMAX
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)
     &                                  +TJ*HJ1(NRG,NZG)
     &                                  +TJ*TJ*HJ2(NRG,NZG))
            FJRZ(NRG,NZG)=FJPRZ(NRG,NZG)
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
               PPRZ(NRG,NZG)=PPSI
            ELSE
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C     ----- Given pressure and parallel current profiles -----
C     ----- Given pressure and safety factor profiles ------
C
      ELSEIF(IMDLEQF.EQ.3.OR.IMDLEQF.EQ.4) THEN
         IF(IMDLEQF.EQ.3) THEN
            CALL EQIPJP
         ELSE
            CALL EQIPQP
         ENDIF
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFIPV(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=-FPSI*DFPSI/(2.D0*PI*RMU0*RG(NRG))
               TTRZ(NRG,NZG)=FPSI
               PPRZ(NRG,NZG)=PPSI
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)+HJ1(NRG,NZG))
            FJRZ(NRG,NZG)=FJPRZ(NRG,NZG)
         ENDDO
         ENDDO
C
C           ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
         ENDDO
         ENDDO
         RIPX=(FJ0+FJ1)*1.D-6
C
C         WRITE(6,'(A,1P3E12.4)') 'RIPX=',RIPX,FJP*1.D-6,FJT*1.D-6
C
      ENDIF
C
C     ----- CALCULATE POLOIDAL CURRENT -----
C
      IMDLEQF=MOD(MDLEQF,5)
      DO NRV=1,NRVMAX
         PSIPNL=PSIPNV(NRV)
         IF (IMDLEQF.EQ.0) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            CALL EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
            TTVL=SQRT((2.D0*PI*BB*RR)**2
     &                  +2.D0*RMU0*RRC
     &                  *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                  *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
         ELSEIF (IMDLEQF.EQ.1) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            TTVL=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSEIF (IMDLEQF.EQ.2) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSEIF (IMDLEQF.EQ.3) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=FPSI
         ELSEIF (IMDLEQF.EQ.4) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=FPSI
         ENDIF
         TTV(NRV)=TTVL
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
C     ----- Poloidal Field Coil current -----
C
      DO NPFC=1,NPFCMAX
         RIPFCL=RIPFC(NPFC)
         RPFCL=RPFC(NPFC)
         ZPFCL=ZPFC(NPFC)
         CALL EQNRZX(RPFCL,ZPFCL,NRC,NZC,IERR)
         IF(NZC.GE.1.AND.NZC.LT.NZGMAX.AND.
     &      NRC.GE.1.AND.NRC.LT.NRGMAX) THEN
            DRG=RG(NRC+1)-RG(NRC)
            DZG=ZG(NZC+1)-ZG(NZC)
            FACTR=(RPFCL-RG(NRC))/DRG
            FACTZ=(ZPFCL-ZG(NZC))/DZG
            FACTRC=1.D0-FACTR
            FACTZC=1.D0-FACTZ
C
            FJRZ(NRC  ,NZC  )=FJRZ(NRC  ,NZC  )
     &                       +2.D0*PI*RMU0*FACTRC*FACTZC
     &                        *RIPFCL*1.D6/(DRG*DZG)
            FJRZ(NRC+1,NZC  )=FJRZ(NRC+1,NZC  )
     &                       +2.D0*PI*RMU0*FACTR *FACTZC
     &                        *RIPFCL*1.D6/(DRG*DZG)
            FJRZ(NRC  ,NZC+1)=FJRZ(NRC  ,NZC+1)
     &                       +2.D0*PI*RMU0*FACTRC*FACTZ 
     &                        *RIPFCL*1.D6/(DRG*DZG)
            FJRZ(NRC+1,NZC+1)=FJRZ(NRC+1,NZC+1)
     &                       +2.D0*PI*RMU0*FACTR *FACTZ 
     &                        *RIPFCL*1.D6/(DRG*DZG)
         ELSE
            WRITE(6,*) 'XX (PFC) OUT OF REGION: ZPFC=',ZPFCL
            WRITE(6,*) 'XX (PFC) OUT OF REGION: RPFC=',RPFCL
         ENDIF
      ENDDO
C
C     ----- Set RHS vector FVB -----
C
      DO N=1,MLMAX
         FVB(N)=0
      ENDDO
C
      DO NZG=1,NZGMAX-1
      DO NRG=1,NRGMAX-1
         DRG=RG(NRG+1)-RG(NRG)
         FACTH(1)=DRG
         FACTH(2)=DRG**2
         FACTH(3)=DRG**2
         FACTH(4)=DRG**3
C
         DZG=ZG(NZG+1)-ZG(NZG)
         FACTM(1)=DZG
         FACTM(2)=DZG**2
         FACTM(3)=DZG**2
         FACTM(4)=DZG**3
C
         DO I1=1,2
         DO I2=1,2
         DO J1=1,2
         DO J2=1,2
            RJ(1)=FJRZ(NRG+I2-1,NZG+J2-1)
            IF(NRG+I2-1.GT.1.AND.
     &         NRG+I2-1.LT.NRGMAX) THEN
               RJ(2)=(FJRZ(NRG+I2,NZG+J2-1)-FJRZ(NRG+I2-2,NZG+J2-1))
     &              /(2.D0*DRG)
            ELSE  
               RJ(2)=0
            ENDIF
            IF(NZG+J2-1.GT.1.AND.
     &         NZG+J2-1.LT.NZGMAX) THEN
               RJ(3)=(FJRZ(NRG+I2-1,NZG+J2)-FJRZ(NRG+I2-1,NZG+J2-2))
     &              /(2.D0*DZG)
            ELSE  
               RJ(3)=0
            ENDIF
            IF(NRG+I2-1.GT.1.AND.
     &         NRG+I2-1.LT.NRGMAX.AND.
     &         NZG+J2-1.GT.1.AND.
     &         NZG+J2-1.LT.NZGMAX) THEN
               RJ(4)=(FJRZ(NRG+I2,NZG+J2)-FJRZ(NRG+I2,NZG+J2-2)
     &               -FJRZ(NRG+I2-2,NZG+J2)+FJRZ(NRG+I2-2,NZG+J2-2))
     &              /(4.D0*DZG*DRG)
            ELSE
               RJ(4)=0
            ENDIF
C
            DO K=1,4
            DO L=1,4
               N=4*((NZG-1)*NRGMAX+NRG-1)+K+4*(I1-1)+4*NRGMAX*(J1-1)
               K1=K1A(K,L)
               M1=M1A(K,L)
               FVB(N)=FVB(N)+RM(K1,I1,I2)*FACTH(K1)
     &                      *RM(M1,J1,J2)*FACTM(M1)*RMU0*1.D6*RJ(L)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C     ----- SET BOUNDARY CONDITION -----
C
      DO NBA=1,NBAMAX
         NRG=NRGB(NBA)
         NZG=NZGB(NBA)
         R1=RG(NRG)
         Z1=ZG(NZG)
         CALL EQPSIB(R1,Z1,PSIBRZ)
         IF(NRG.EQ.1.OR.NRG.EQ.NRGMAX) THEN
         DO K=1,3,2
            SUML(K)=PSIBRZ(K)
            DO NZG2=2,NZGMAX-1
            DO NRG2=2,NRGMAX-1
               SUML(K)=SUML(K)-PSIBGR(NRG2,NZG2,NBA,K)*FJRZ(NRG2,NZG2)
            ENDDO
            ENDDO
            N=4*((NZG-1)*NRGMAX+NRG-1)+K
            FVB(N)=SUML(K)
         ENDDO
         ELSE
         DO K=1,2
            SUML(K)=PSIBRZ(K)
            DO NZG2=2,NZGMAX-1
            DO NRG2=2,NRGMAX-1
               SUML(K)=SUML(K)-PSIBGR(NRG2,NZG2,NBA,K)*FJRZ(NRG2,NZG2)
            ENDDO
            ENDDO
            N=4*((NZG-1)*NRGMAX+NRG-1)+K
            FVB(N)=SUML(K)
         ENDDO
         ENDIF
      ENDDO

      DEALLOCATE(HJ0,HJ1,HJ2)
      RETURN
      END
C
C     **********************************
C     *     Calulate plasma current    *
C     **********************************
C
      SUBROUTINE EQRJPX(IERR)
C
      USE libspl1d
      INCLUDE '../eq/eqcomx.inc'
C
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: HJ0,HJ1,HJ2

      ALLOCATE(HJ0(NRGM,NZGM),HJ1(NRGM,NZGM),HJ2(NRGM,NZGM))
C
C     ----- calculate flux average from PSIRZ(R,Z) -----
C
      IERR=0
      CALL EQCALV(IERR)
      IF(IERR.NE.0) RETURN
C
      IMDLEQF=MOD(MDLEQF,5)
C
      IF(IMDLEQF.EQ.0) THEN
         RRC=RAXIS
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=(1.D0-RRC**2/RG(NRG)**2)*HJ0(NRG,NZG)
               HJ2(NRG,NZG)=-(RRC/RG(NRG))*HJPSI
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               HJ2(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C        ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         FJ2=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
            FJ2=FJ2+HJ2(NRG,NZG)*DVOL
         ENDDO
         ENDDO
         TJ=(RIP*1.D6-FJ1)/FJ2
         RIPX=RIP
C
         DO NRG=1,NRGMAX
         DO NZG=1,NZGMAX
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ1(NRG,NZG)+TJ*HJ2(NRG,NZG))
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
               TTRZ(NRG,NZG)=SQRT((2.D0*PI*BB*RR)**2
     &                           +2.D0*RMU0*RRC
     &                           *(2.D0*PI*TJ*HJPSID-RRC*PPSI))
               PPRZ(NRG,NZG)=PPSI
            ELSE
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C     ----- Given pressure and poloidal current profiles -----
C
      ELSEIF(IMDLEQF.EQ.1) THEN

         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=-2.D0*PI*BB*RR*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
               HJ2(NRG,NZG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               HJ2(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C        ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         FJ2=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
            FJ2=FJ2+HJ2(NRG,NZG)*DVOL
         ENDDO
         ENDDO
C
C        ----- Adjust plasma current -----
C
         IF(FJ1.GT.0.D0) THEN
            TJ=(-FJ1+SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ELSE
            TJ=(-FJ1-SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ENDIF
         RIPX=RIP
C
C         WRITE(6,'(A,1P4E12.4)') 'FJ0,FJ1,FJ2,TJ=',FJ0,FJ1,FJ2,TJ
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)
     &                                  +TJ*HJ1(NRG,NZG)
     &                                  +TJ*TJ*HJ2(NRG,NZG))
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               PPRZ(NRG,NZG)=PPSI
            ELSE
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C     ----- Given pressure and parallel current profiles -----
C
      ELSEIF(IMDLEQF.EQ.2) THEN
         CALL EQIPJP
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFIPV(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=-2.D0*PI*BB*RR       *DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
               HJ2(NRG,NZG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                              /(2.D0*PI*RMU0*RG(NRG))
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               HJ2(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C        ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         FJ2=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
            FJ2=FJ2+HJ2(NRG,NZG)*DVOL
         ENDDO
         ENDDO
C
         IF(FJ1.GT.0.D0) THEN
            TJ=(-FJ1+SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ELSE
            TJ=(-FJ1-SQRT(FJ1**2+4.D0*FJ2*(RIP*1.D6-FJ0)))
     &         /(2.D0*FJ2)
         ENDIF
         RIPX=RIP
C
         DO NRG=1,NRGMAX
         DO NZG=1,NZGMAX
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)
     &                                  +TJ*HJ1(NRG,NZG)
     &                                  +TJ*TJ*HJ2(NRG,NZG))
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
               PPRZ(NRG,NZG)=PPSI
            ELSE
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
C     ----- Given pressure and parallel current profiles -----
C     ----- Given pressure and safety factor profiles ------
C
      ELSEIF(IMDLEQF.EQ.3.OR.IMDLEQF.EQ.4) THEN
         IF(IMDLEQF.EQ.3) THEN
            CALL EQIPJP
         ELSE
            CALL EQIPQP
         ENDIF
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIX(NRG,NZG)*PSI0.GT.0.D0.AND.
     &           ZG(NZG).LE.ZLIMP.AND.
     &           ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=1.D0-PSIX(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFIPV(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)=-2.D0*PI*RG(NRG)*DPPSI
               HJ1(NRG,NZG)=-FPSI*DFPSI/(2.D0*PI*RMU0*RG(NRG))
               TTRZ(NRG,NZG)=FPSI
               PPRZ(NRG,NZG)=PPSI
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
               TTRZ(NRG,NZG)=2.D0*PI*BB*RR
               PPRZ(NRG,NZG)=0.D0
            ENDIF
            FJPRZ(NRG,NZG)=2.D0*PI*RMU0*(HJ0(NRG,NZG)+HJ1(NRG,NZG))
         ENDDO
         ENDDO
C
C           ----- Integrate plasma current -----
C
         FJ0=0.D0
         FJ1=0.D0
         DO NZG=2,NZGMAX-1
         DO NRG=2,NRGMAX-1
            DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
            DZG=0.5D0*(ZG(NRG+1)-ZG(NRG-1))
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
         ENDDO
         ENDDO
         RIPX=(FJ0+FJ1)*1.D-6
C
C         WRITE(6,'(A,1P3E12.4)') 'RIPX=',RIPX,FJP*1.D-6,FJT*1.D-6
C
      ENDIF
C
C     ----- CALCULATE POLOIDAL CURRENT -----
C
      IMDLEQF=MOD(MDLEQF,5)
      DO NRV=1,NRVMAX
         PSIPNL=PSIPNV(NRV)
         IF (IMDLEQF.EQ.0) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            CALL EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
            TTVL=SQRT((2.D0*PI*BB*RR)**2
     &                  +2.D0*RMU0*RRC
     &                  *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                  *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
         ELSEIF (IMDLEQF.EQ.1) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            TTVL=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSEIF (IMDLEQF.EQ.2) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSEIF (IMDLEQF.EQ.3) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=FPSI
         ELSEIF (IMDLEQF.EQ.4) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=FPSI
         ENDIF
         TTV(NRV)=TTVL
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
      DEALLOCATE(HJ0,HJ1,HJ2)
      RETURN
      END
C
C     ****** BOUNDARY PSI ON WALL ******
C
      SUBROUTINE EQPSIB(R,Z,PSIBRZ)
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION PSIBRZ(4),PSIV(4,0:5)
C
      R01=RR
      R02=RR**2
      R03=RR**3
      R04=RR**4
      R05=RR**5
C
      PSIV(1,0)= 1.D0
      PSIV(2,0)= 0.D0
      PSIV(3,0)= 0.D0
      PSIV(4,0)= 0.D0
C
      PSIV(1,1)=(R**2-R02)/(2.D0*R01)
      PSIV(2,1)= R/R01
      PSIV(3,1)= 0.D0
      PSIV(4,1)= 0.D0
C
      PSIV(1,2)=-(R**2*Z**2
     &           -(R**2-R02)**2/4.D0)/R02
      PSIV(2,2)=-(2.D0*R*Z**2
     &           -R*(R**2-R02))/R02
      PSIV(3,2)=-2.D0*R**2*Z/R02
      PSIV(4,2)=-4.D0*R*Z/R02
C
      PSIV(1,3)=(R**2*Z**4
     &          -1.5D0*(R**2-R02)*R**2*Z**2
     &          +0.125D0*(R**2-R02)**3)/R03
      PSIV(2,3)=(2.D0*R*Z**4
     &          -3.D0*R**3.D0*Z**2
     &          -3.D0*R*(R**2-R02)*Z**2
     &          +0.75D0*R*(R**2-R02)**2)/R03
      PSIV(3,3)=(4.D0*R**2*Z**3
     &          -3.D0*(R**2-R03)*R**2*Z)/R03
      PSIV(4,3)=(8.D0*R*Z**3
     &          -6.D0*R**3*Z
     &          -6.D0*(R**2-R02)*R*Z)/R03
C
      PSIV(1,4)=-(0.2D0*R**2*Z**6
     &           -0.25D0*(3.D0*R**2-2.D0*R02)*R**2*Z**4
     &           +0.375D0*(R**2-R02)**2*R**2*Z**2
     &           -0.015625D0*(R**2-R02)**4)/R04
      PSIV(2,4)=-(0.4D0*R*Z**6-1.5D0*R**3*Z**4
     &           -0.5D0*R*(3.D0*R**2-2.D0*R02)*Z**4
     &           +1.5D0*R**3*(R**2-R02)*Z**2
     &           +0.75D0*R*(R**2-R02)**2*Z**2
     &           -0.125D0*R*(R**2-R02)**3)/R04
      PSIV(3,4)=-(1.2D0*R**2*Z**5
     &           -(3.D0*R**2-2.D0*R02)*R**2*Z**3
     &           +0.75D0*(R**2-R02)**2*R**2*Z)/R04
      PSIV(4,4)= (2.4D0*R*Z**5
     &           -6.D0*R**3*Z**3
     &           -2.D0*(3.D0*R**2-2.D0*R02)*Z**3
     &           +3.D0*R**3*(R**2-R02)*Z
     &           +1.5D0*R*(R**2-R02)*Z)/R04
C
      PSIV(1,5)= (4.D0/7.D0*R**2*Z**8
     &           -2.D0*(2.D0*R**2-R02)*R**2*Z**6
     &           +2.5D0*(2.D0*R**2-R02)*(R**2-R02)*R**2*Z**4
     &           -1.25D0*(R**2-R02)**3*R**2*Z**2
     &           +0.03125D0*(R**2-R02)**5)/R05
      PSIV(2,5)= (8.D0/7.D0*R*Z**8
     &           -8*R**3*Z**6
     &           -4.D0*R*(2.D0*R**2-R02)*Z**6
     &           +5.D0*R**3*(4.D0*R**2-3.D0*R02)*Z**4
     &           +5.D0*R*(2.D0*R**2-R02)*(R**2-R02)*Z**4
     &           -7.5D0*R**3*(R**2-R02)**2*Z**2
     &           -2.5D0*R*(R**2-R02)**3*Z**2
     &           +0.3125D0*R*(R**2-R02)**4)/R05
      PSIV(3,5)= (32.D0/7.D0*R**2*Z**7
     &           -12.D0*(2.D0*R**2-R02)*R**2*Z**5
     &           +10.D0*(2.D0*R**2-R02)*(R**2-R02)*R**2*Z**3
     &           -2.5D0*(R**2-R02)**3*R**2*Z)/R05
      PSIV(4,5)= (64.D0/7.D0*R*Z**7
     &           -48.D0*R**3*Z**5
     &           -24.D0*R*(2.D0*R**2-R02)*Z**5
     &           +20.D0*R**3*(4.D0*R**2-3.D0*R02)*Z**3
     &           +20.D0*R*(2.D0*R**2-R02)*(R**2-R02)*Z**3
     &           -15.D0*R**3*(R**2-R02)**2*Z
     &           -5.D0*R*(R**2-R02)**3*Z)/R05
C
      PSIBRZ(1)=0.D0
      PSIBRZ(2)=0.D0
      PSIBRZ(3)=0.D0
      PSIBRZ(4)=0.D0
      DO I=0,5
         PSIBRZ(1)=PSIBRZ(1)+2.D0*PI*PSIV(1,I)*PSIB(I)
         PSIBRZ(2)=PSIBRZ(2)+2.D0*PI*PSIV(2,I)*PSIB(I)
         PSIBRZ(3)=PSIBRZ(3)+2.D0*PI*PSIV(3,I)*PSIB(I)
         PSIBRZ(4)=PSIBRZ(4)+2.D0*PI*PSIV(4,I)*PSIB(I)
      ENDDO
C
      RETURN
      END
C
C  ************** RMAX,RMIN at PSI=0,RR,RDLT,RKAP **********************
C
      SUBROUTINE EQCALR(RRR,RRA,RRK,RRD)
C
      INCLUDE '../eq/eqcomq.inc'
      DIMENSION XA(NTVM),YA(2,NTVM)
C
      RINIT=REDGE
      ZINIT=ZAXIS
C
      NUMT=0
      NUMB=0
      RMIN=RINIT
      RMAX=RINIT
      ZMIN=ZINIT
      ZMAX=ZINIT
C
      CALL EQMAGS(RINIT,ZINIT,NTVMAX,XA,YA,NA,IERR)
C
      DO N=2,NA
         R=YA(1,N)
         Z=YA(2,N)
C
         IF(Z.GE.ZMAX) NUMT=N
         IF(Z.LE.ZMIN) NUMB=N
         RMIN=MIN(RMIN,R)
         RMAX=MAX(RMAX,R)
         ZMIN=MIN(ZMIN,Z)
         ZMAX=MAX(ZMAX,Z)
      ENDDO
      RRR=0.5D0*(RMIN+RMAX)
      RRA=0.5D0*(RMAX-RMIN)            
      RRK=(ZMAX-ZMIN)/(RMAX-RMIN)
      IF(NUMT.EQ.0.OR.NUMB.EQ.0) THEN
         RRD=0.D0
      ELSE
         RRD=(RRR-0.5D0*(YA(1,NUMT)+YA(1,NUMB)))/RRA
      ENDIF
C
C      WRITE(6,'(4(A,1PE12.4))') 'RMIN=',RMIN,'  RMAX=',RMAX,
C     &                        '  ZMIN=',ZMIN,'  ZMAX=',ZMAX
      WRITE(6,'(4(A,1PE12.4))') 'RRR =',RRR, '  RRA =',RRA,
     &                        '  RRK =',RRK, '  RRD =',RRD
C
      RETURN
      END
C
C     ---- calculate Psi(R,Z=0) -----
C
      FUNCTION PSIXZ0(R)
      INCLUDE '../eq/eqcomx.inc'
      PSIXZ0=PSIXF(R,ZAXIS)
      RETURN
      END
C
C     ---- calculate Psi(R,Z) -----
C
      FUNCTION PSIXF(R,Z)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
      PSIXF=0.D0
C
      CALL EQNRZX(R,Z,NRG,NZG,IERR)
      IF(IERR.NE.0.AND.IERR.NE.2.AND.IERR.NE.4) RETURN
C
      DRG=RG(NRG+1)-RG(NRG)
      DZG=ZG(NZG+1)-ZG(NZG)
      VRG=(R-RG(NRG))/DRG
      VZG=(Z-ZG(NZG))/DZG
C
      HR0= HRMT0(1.D0-VRG)
      HZ0= HRMT0(1.D0-VZG)
      HR1=-HRMT1(1.D0-VRG)*DRG
      HZ1=-HRMT1(1.D0-VZG)*DZG
C
      HR0C= HRMT0(VRG)
      HZ0C= HRMT0(VZG)
      HR1C= HRMT1(VRG)*DRG
      HZ1C= HRMT1(VZG)*DZG
C
      PSIXF=UPSIX(1,NRG  ,NZG  )*HR0 *HZ0
     &     +UPSIX(2,NRG  ,NZG  )*HR1 *HZ0
     &     +UPSIX(3,NRG  ,NZG  )*HR0 *HZ1
     &     +UPSIX(4,NRG  ,NZG  )*HR1 *HZ1
     &     +UPSIX(1,NRG+1,NZG  )*HR0C*HZ0
     &     +UPSIX(2,NRG+1,NZG  )*HR1C*HZ0
     &     +UPSIX(3,NRG+1,NZG  )*HR0C*HZ1
     &     +UPSIX(4,NRG+1,NZG  )*HR1C*HZ1
     &     +UPSIX(1,NRG  ,NZG+1)*HR0 *HZ0C
     &     +UPSIX(2,NRG  ,NZG+1)*HR1 *HZ0C
     &     +UPSIX(3,NRG  ,NZG+1)*HR0 *HZ1C
     &     +UPSIX(4,NRG  ,NZG+1)*HR1 *HZ1C
     &     +UPSIX(1,NRG+1,NZG+1)*HR0C*HZ0C
     &     +UPSIX(2,NRG+1,NZG+1)*HR1C*HZ0C
     &     +UPSIX(3,NRG+1,NZG+1)*HR0C*HZ1C
     &     +UPSIX(4,NRG+1,NZG+1)*HR1C*HZ1C
C
      RETURN
      END
C
C     ---- calculate Psi(R,Z) and its derivatives -----
C
      SUBROUTINE EQPSIX(R,Z,DPSIDR,DPSIDZ,IERR)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
C
      CALL EQNRZX(R,Z,NRG,NZG,IERR)
      IF(IERR.GT.4) RETURN
C
      DRG=RG(NRG+1)-RG(NRG  )
      DZG=ZG(NZG+1)-ZG(NZG  )
      VRG=(R-RG(NRG  ))/DRG
      VZG=(Z-ZG(NZG  ))/DZG
C
C      WRITE(6,'(A,1P2E12.4)') 'DRG,DZG=',DRG,DZG
C
      HR0= HRMT0(1.D0-VRG)
      HZ0= HRMT0(1.D0-VZG)
      HR1=-HRMT1(1.D0-VRG)*DRG
      HZ1=-HRMT1(1.D0-VZG)*DZG
C
      HR0C= HRMT0(VRG)
      HZ0C= HRMT0(VZG)
      HR1C= HRMT1(VRG)*DRG
      HZ1C= HRMT1(VZG)*DZG
C
      DHR0=-DHRMT0(1.D0-VRG)/DRG
      DHZ0=-DHRMT0(1.D0-VZG)/DZG
      DHR1= DHRMT1(1.D0-VRG)
      DHZ1= DHRMT1(1.D0-VZG)
C
      DHR0C= DHRMT0(VRG)/DRG
      DHZ0C= DHRMT0(VZG)/DZG
      DHR1C= DHRMT1(VRG)
      DHZ1C= DHRMT1(VZG)
C
C      WRITE(6,'(A:1P4E12.4)')
C     &     'V:',VZG,DHRMT0(1.D0-VRG),DHRMT0(VRG)
C      WRITE(6,'(A:1P4E12.4)')
C     &     'H:',HR0,DHZ0,HR0C,DHZ0C
C      WRITE(6,'(A:1P4E12.4)')
C     &     'U:',UPSIX(1,NRG  ,NZG  ),UPSIX(1,NRG+1,NZG  ),
C     &          UPSIX(1,NRG  ,NZG+1),UPSIX(1,NRG+1,NZG+1)
C
      DPSIDR=UPSIX(1,NRG  ,NZG  )*DHR0 *HZ0
     &      +UPSIX(2,NRG  ,NZG  )*DHR1 *HZ0
     &      +UPSIX(3,NRG  ,NZG  )*DHR0 *HZ1
     &      +UPSIX(4,NRG  ,NZG  )*DHR1 *HZ1
     &      +UPSIX(1,NRG+1,NZG  )*DHR0C*HZ0
     &      +UPSIX(2,NRG+1,NZG  )*DHR1C*HZ0
     &      +UPSIX(3,NRG+1,NZG  )*DHR0C*HZ1
     &      +UPSIX(4,NRG+1,NZG  )*DHR1C*HZ1
     &      +UPSIX(1,NRG  ,NZG+1)*DHR0 *HZ0C
     &      +UPSIX(2,NRG  ,NZG+1)*DHR1 *HZ0C
     &      +UPSIX(3,NRG  ,NZG+1)*DHR0 *HZ1C
     &      +UPSIX(4,NRG  ,NZG+1)*DHR1 *HZ1C
     &      +UPSIX(1,NRG+1,NZG+1)*DHR0C*HZ0C
     &      +UPSIX(2,NRG+1,NZG+1)*DHR1C*HZ0C
     &      +UPSIX(3,NRG+1,NZG+1)*DHR0C*HZ1C
     &      +UPSIX(4,NRG+1,NZG+1)*DHR1C*HZ1C
C
      DPSIDZ=UPSIX(1,NRG  ,NZG  )*HR0 *DHZ0
     &      +UPSIX(2,NRG  ,NZG  )*HR1 *DHZ0
     &      +UPSIX(3,NRG  ,NZG  )*HR0 *DHZ1
     &      +UPSIX(4,NRG  ,NZG  )*HR1 *DHZ1
     &      +UPSIX(1,NRG+1,NZG  )*HR0C*DHZ0
     &      +UPSIX(2,NRG+1,NZG  )*HR1C*DHZ0
     &      +UPSIX(3,NRG+1,NZG  )*HR0C*DHZ1
     &      +UPSIX(4,NRG+1,NZG  )*HR1C*DHZ1
     &      +UPSIX(1,NRG  ,NZG+1)*HR0 *DHZ0C
     &      +UPSIX(2,NRG  ,NZG+1)*HR1 *DHZ0C
     &      +UPSIX(3,NRG  ,NZG+1)*HR0 *DHZ1C
     &      +UPSIX(4,NRG  ,NZG+1)*HR1 *DHZ1C
     &      +UPSIX(1,NRG+1,NZG+1)*HR0C*DHZ0C
     &      +UPSIX(2,NRG+1,NZG+1)*HR1C*DHZ0C
     &      +UPSIX(3,NRG+1,NZG+1)*HR0C*DHZ1C
     &      +UPSIX(4,NRG+1,NZG+1)*HR1C*DHZ1C
C
      IERR=0
C
      RETURN
      END
C
C     ---- calculate NRZ and NGZ for (R,Z) -----
C
      SUBROUTINE EQNRZX(R,Z,NRG,NZG,IERR)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
C
      IF(RG(NRGMAX).EQ.RG(1)) THEN
         IERR=9
         RETURN
      ENDIF
      IF(ZG(NZGMAX).EQ.ZG(1)) THEN
         IERR=8
         RETURN
      ENDIF
C
      FRG=1.D0/(RG(NRGMAX)-RG(1))
      NRG=NINT((R-RG(1))*FRG*(NRGMAX-1))+1
      IF(NRG.LT.1) THEN
         IERR=1
         NRG=1
      ENDIF
      IF(NRG.GT.NRGMAX-1) THEN
         IERR=2
         NRG=NRGMAX-1
      ENDIF
      FZG=1.D0/(ZG(NZGMAX)-ZG(1))
      NZG=NINT((Z-ZG(1))*FZG*(NZGMAX-1))+1
      IF(NZG.LT.1) THEN
         IERR=3
         NZG=1
      ENDIF
      IF(NZG.GT.NZGMAX-1) THEN
         IERR=4
         NZG=NZGMAX-1
      ENDIF
C
 5001 IF(NRG.GE.NRGMAX-1) GOTO 5002
      IF((R-RG(NRG+1))*FRG.LE.0.D0) GOTO 5002
         NRG=NRG+1
         GOTO 5001
 5002 CONTINUE
 5003 IF(NRG.LE.1) GOTO 5004
      IF((R-RG(NRG  ))*FRG.GE.0.D0) GOTO 5004
         NRG=NRG-1
         GOTO 5003
 5004 CONTINUE
      IF(NRG.LT.1)     NRG=1
C
 5005 IF(NZG.GE.NZGMAX-1) GOTO 5006
      IF((Z-ZG(NZG+1))*FZG.LE.0.D0) GOTO 5006
         NZG=NZG+1
         GOTO 5005
 5006 CONTINUE
 5007 IF(NZG.LE.1) GOTO 5008
      IF((Z-ZG(NZG  ))*FZG.GE.0.D0) GOTO 5008
         NZG=NZG-1
         GOTO 5007
 5008 CONTINUE
      IF(NZG.LT.1)     NZG=1
C
C      WRITE(6,'(A,1P,E12.4,I5,2E12.4)') 'R:',R,NRG,RG(NRG-1),RG(NRG)
C      WRITE(6,'(A,1P,E12.4,I5,2E12.4)') 'Z:',Z,NZG,ZG(NZG-1),ZG(NZG)
      RETURN
      END
C
      FUNCTION HRMT0(X)
      IMPLICIT NONE
      REAL*8 HRMT0,X
      HRMT0=X*X*(3.D0-2.D0*X)
      RETURN
      END
C
      FUNCTION HRMT1(X)
      IMPLICIT NONE
      REAL*8 HRMT1,X
      HRMT1=X*X*(X-1.D0)
      RETURN
      END
C
      FUNCTION DHRMT0(X)
      IMPLICIT NONE
      REAL*8 DHRMT0,X
      DHRMT0=6.D0*X*(1.D0-X)
      RETURN
      END
C
      FUNCTION DHRMT1(X)
      IMPLICIT NONE
      REAL*8 DHRMT1,X
      DHRMT1=X*(3.D0*X-2.D0)
      RETURN
      END
C
C     ****** PSI generated by current FJRZ ******
C
      FUNCTION EQPSIV(R1,Z1)
C
      INCLUDE '../eq/eqcomx.inc'
C
      SUML=0.D0
      DO NZG=2,NZGMAX-1
      DO NRG=2,NRGMAX-1
         R2=RG(NRG)
         Z2=ZG(NZG)
         DRG=0.5D0*(RG(NRG+1)-RG(NRG-1))
         DZG=0.5D0*(ZG(NZG+1)-ZG(NZG-1))
         SUML=SUML
     &       -DRG*DZG*FJRZ(NRG,NZG)*RGSGRF(R1,Z1,R2,Z2)
      ENDDO
      ENDDO
      EQPSIV=SUML
      RETURN
      END
C
C     ***** Green function of Grad-Shafranov equation *****
C
      FUNCTION RGSGRF(R1,Z1,R2,Z2)
C
      IMPLICIT NONE
      REAL*8 RGSGRF,R1,Z1,R2,Z2
      REAL*8 PI,RK,ELLFC,ELLEC
      INTEGER IERR1,IERR2
      DATA PI/3.14159265358979D0/
C
      RK=4.D0*R1*R2/((R1+R2)**2+(Z1-Z2)**2)
C
      RGSGRF=SQRT(R1*R2)/(2.D0*PI*RK)
     &      *((2.D0-RK*RK)*ELLFC(RK,IERR1)-2.D0*ELLEC(RK,IERR2))
      IF(IERR1.NE.0) WRITE(6,*) 'XX RGSGRF: ELLFC: IERR1=',IERR1
      IF(IERR2.NE.0) WRITE(6,*) 'XX RGSGRF: ELLFC: IERR2=',IERR1
      RETURN
      END
