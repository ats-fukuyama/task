C     $Id$
C
C     ***** Calculated Flux Functions from PSIRZ *****
C
      SUBROUTINE EQCALQ(NRMAX1,NTHMAX1,NSUMAX1,IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      IERR=0
C
C     ----- CHECK NRMAX,NTHMAX,NSUMAX are not greater than *M -----
C
      IF(NRMAX1.GT.NRM) THEN
         WRITE(6,'(A,2I5)') 
     &        'NRMAX1.GT.NRM: NRMAX1,NRM=',NRMAX1,NRM
         IERR=IERR+1
      ENDIF
      IF(NTHMAX1.GT.NTHM) THEN
         WRITE(6,'(A,2I5)') 
     &        'NTHMAX1.GT.NTHM: NTHMAX1,NTHM=',NTHMAX1,NTHM
         IERR=IERR+2
      ENDIF
      IF(NSUMAX1.GT.NSUM) THEN
         WRITE(6,'(A,2I5)') 
     &        'NSUMAX1.GT.NSUM: NSUMAX1,NSUM=',NSUMAX1,NSUM
         IERR=IERR+4
      ENDIF
      IF(IERR.NE.0) RETURN
C
      NRMAX=NRMAX1
      NTHMAX=NTHMAX1
      NSUMAX=NSUMAX1
C
      CALL EQSETP
C
      IF(NSUMAX.LE.0) THEN
         CALL EQCALQP(IERR)
      ELSE
         CALL EQCALQP(IERR)
         CALL EQCALQV(IERR)
      ENDIF
C
      CALL EQSETS(IERR)
      RETURN
      END
C
C     ***** SETUP DATA (spline PSIRZ and find axis) *****
C
      SUBROUTINE EQSETP
C
      INCLUDE '../eq/eqcomq.inc'
C
      DIMENSION PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM)
      DIMENSION DERIV(NPSM)
      EXTERNAL EQPSID
C
      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,URZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for PSIRZ: IERR=',IERR
      CALL SPL1D(PSIPS,PPPS,  DERIV,UPPPS, NPSMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PPPS: IERR=',IERR
      CALL SPL1D(PSIPS,TTPS,  DERIV,UTTPS, NPSMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTPS: IERR=',IERR
C
      DELT=1.D-8
      EPS=1.D-8
      ILMAX=20
      LIST=0
      RINIT=RR
      ZINIT=0.D0
      CALL NEWTN(EQPSID,RINIT,ZINIT,RAXIS,ZAXIS,
     &            DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) WRITE(6,*) 'XX EQSETP: NEWTN ERROR: IER=',IER
      SAXIS=PSIG(RAXIS,ZAXIS)
C
      RETURN
      END
C
C     ***** CALCULATE FLUX VARIABLES IN PLASMA *****
C
      SUBROUTINE EQCALQP(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      EXTERNAL PSIAX,EQDERV
      DIMENSION XA(NNM),YA(2,NNM)
      DIMENSION XCHI0(NNM),XCHI1(NNM),RCHI(NNM),ZCHI(NNM),DXCHI(NNM)
      DIMENSION URCHI(4,NNM),UZCHI(4,NNM)
C
      IERR=0
C
C     ----- SET REDGE, DR, DTH -----
C
      REDGE=ZBRENT(PSIAX,RR,RR+RB,1.D-8)
C
C      WRITE(6,*) REDGE,RAXIS,RB,RA
C
      IF(NSUMAX.EQ.0) THEN
         DR=(REDGE-RAXIS)/(NRMAX-1)
         NRPMAX=NRMAX
      ELSE
         DR=RB/(NRMAX-1)
         NRPMAX=NINT(RA/DR)+1
         DR=(REDGE-RAXIS)/(NRPMAX-1)
      ENDIF
      DTH=2*PI/NTHMAX
C
C     ----- SET NUMBER OF DIVISION for integration -----
C
      NMAX=200
      IF(NMAX.GT.NNM) NMAX=NNM
C
C     ----- CALCULATE PSI,PPS,TTS,FTS and RPS, ZPS on magnetic surfaces -----
C
      NR=1
      DO NTH=1,NTHMAX+1
         RPS(NTH,NR)=RAXIS
         ZPS(NTH,NR)=ZAXIS
      ENDDO
      PSS(1)=PSIG(RAXIS,ZAXIS)
      PPS(1)=PPFUNC(PSS(1))
      TTS(1)=TTFUNC(PSS(1))
C
      DO NR=2,NRPMAX
         RINIT=RAXIS+DR*(NR-1)
         ZINIT=ZAXIS
         PSS(NR)=PSIG(RINIT,ZINIT)
         PPS(NR)=PPFUNC(PSS(NR))
         TTS(NR)=TTFUNC(PSS(NR))
C
         CALL EQMAGS(RINIT,ZINIT,NMAX,XA,YA,NA,IERR)
C
         SUMS=0.D0
         SUMV=0.D0
         SUMQ=0.D0
         RMIN=RINIT
         RMAX=RINIT
         BMIN=ABS(2.D0*BB)
         BMAX=0.D0
         SUMAVBR=0.D0
         SUMAVRR=0.D0
         SUMAVR1=0.D0
         SUMAVR2=0.D0
         SUML=0.D0
         XCHI0(1)=0.D0
         XCHI1(1)=0.D0
         RCHI(1)=RINIT
         ZCHI(1)=ZINIT
         DO N=2,NA
            H=XA(N)-XA(N-1)
            R=0.5D0*(YA(1,N-1)+YA(1,N))
            Z=0.5D0*(YA(2,N-1)+YA(2,N))
            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
            BPRL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI)
            B=SQRT((TTS(NR)/(2.D0*PI*R))**2+(BPRL/R)**2)
C
            SUMS=SUMS+H/BPRL
            SUMV=SUMV+H*R/BPRL
            SUMQ=SUMQ+H/(R*BPRL)
C
            SUMAVBR=SUMAVBR+H*BPRL/R
            SUMAVRR=SUMAVRR+H/(R*BPRL)
            SUMAVR1=SUMAVR1+H*R
            SUMAVR2=SUMAVR2+H*R*BPRL
            SUML   =SUML   +H*R/BPRL
C
            XCHI1(N)=SUMQ
            RCHI(N)=YA(1,N)
            ZCHI(N)=YA(2,N)
C
            R=YA(1,N)
            Z=YA(2,N)
            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
            BPRL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI)
            B=SQRT((TTS(NR)/(2.D0*PI*R))**2+(BPRL/R)**2)
C
            RMIN=MIN(RMIN,R)
            RMAX=MAX(RMAX,Z)
            BMIN=MIN(BMIN,B)
            BMAX=MAX(BMAX,B)
         ENDDO
C
         SPS(NR)=SUMS
         VPS(NR)=SUMV
         QPS(NR)=SUMQ*TTS(NR)/(2.D0*PI)**2
         RLEN(NR)=XA(NA)
         RRMIN(NR)=RMIN
         RRMAX(NR)=RMAX
         BBMIN(NR)=BMIN
         BBMAX(NR)=BMAX
         AVBR(NR)=SUMAVBR/SUML
         AVRR(NR)=SUMAVRR/SUML
         AVR1(NR)=SUMAVR1/SUML
         AVR2(NR)=SUMAVR2/SUML
C
C        ----- CALCULATE POLOIDAL COORDINATES -----
C
         IF(NTHMAX.GT.0) THEN
            DO N=1,NA
               XCHI0(N)=2.D0*PI*XA(N)/XA(NA)
               XCHI1(N)=2.D0*PI*XCHI1(N)/XCHI1(NA)
            ENDDO
            RCHI(NA)=RCHI(1)
            ZCHI(NA)=ZCHI(1)
C     
            IF(MDLEQC.EQ.0) THEN
               CALL SPL1D(XCHI0,RCHI,DXCHI,URCHI,NA,4,IERR)
               IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RCHI: IERR=',IERR
               CALL SPL1D(XCHI0,ZCHI,DXCHI,UZCHI,NA,4,IERR)
               IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for ZCHI: IERR=',IERR
            ELSE
               CALL SPL1D(XCHI1,RCHI,DXCHI,URCHI,NA,4,IERR)
               IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RCHI: IERR=',IERR
               CALL SPL1D(XCHI1,ZCHI,DXCHI,UZCHI,NA,4,IERR)
               IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for ZCHI: IERR=',IERR
            ENDIF
C
            RPS(1,NR)=YA(1,1)
            ZPS(1,NR)=YA(2,1)
            DO NTH=2,NTHMAX+1
               TH=DTH*(NTH-1)
               IF(MDLEQC.EQ.0) THEN
                  CALL SPL1DF(TH,RPS(NTH,NR),XCHI0,URCHI,NA,IERR)
                  CALL SPL1DF(TH,ZPS(NTH,NR),XCHI0,UZCHI,NA,IERR)
               ELSE
                  CALL SPL1DF(TH,RPS(NTH,NR),XCHI1,URCHI,NA,IERR)
                  CALL SPL1DF(TH,ZPS(NTH,NR),XCHI1,UZCHI,NA,IERR)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C
C     +++++ SETUP AXIS DATA +++++
C
      NR=1
      SPS(NR)=(4*SPS(2)-SPS(3))/3.D0
      VPS(NR)=(4*VPS(2)-VPS(3))/3.D0
      QPS(NR)=(4*QPS(2)-QPS(3))/3.D0
      RLEN(NR)=0.D0
      RRMIN(NR)=RAXIS
      RRMAX(NR)=RAXIS
      BBMIN(NR)=ABS(TTS(NR)/(2.D0*PI*RAXIS))
      BBMAX(NR)=ABS(TTS(NR)/(2.D0*PI*RAXIS))
C
C     ----- CALCULATE TOROIDAL FLUX -----
C
      FTS(1)=0.D0
      DO NR=2,NRPMAX
         FTS(NR)=FTS(NR-1)
     &           +0.5D0*(QPS(NR)+QPS(NR-1))*(PSS(NR)-PSS(NR-1))
      ENDDO
C
C     ----- CALCULATE EDGE VALUE -----
C
      RINIT=REDGE
      TTSA=TTFUNC(0.D0)
      CALL EQMAGS(RINIT,ZINIT,NMAX,XA,YA,NA,IERR)
C
      SUMS=0.D0
      SUMV=0.D0
      SUMQ=0.D0
      DO N=2,NA
         H=XA(N)-XA(N-1)
         R=0.5D0*(YA(1,N-1)+YA(1,N))
         Z=0.5D0*(YA(2,N-1)+YA(2,N))
         CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
         BPRL=SQRT(DPSIDR**2+DPSIDZ**2)
         BPL=BPRL/(2.D0*PI*R)
         BTL=TTSA/(2.D0*PI*R)
         B=SQRT(BPL**2+BTL**2)
C
         SUMS=SUMS+H/BPRL
         SUMV=SUMV+H*R/BPRL
         SUMQ=SUMQ+H/(R*BPRL)
      ENDDO
C
      SPSA=SUMS
      VPSA=SUMV*2.D0*PI
      QPSA=SUMQ*TTSA/(2.D0*PI)
      FTSA=FTS(NRPMAX)
C
      DO NR=2,NRPMAX
         AVBR(NR)=AVBR(NR)*QPS(NR)**2/(4.D0*FTS(NR)*FTSA)
         AVR1(NR)=AVR1(NR)*QPS(NR)/(2.D0*SQRT(FTS(NR)*FTSA))
         AVR2(NR)=AVR2(NR)*QPS(NR)**2/(4.D0*FTS(NR)*FTSA)
      ENDDO
      AVBR(1)=(4.D0*AVBR(2)-AVBR(3))/3.D0
      AVRR(1)=(4.D0*AVRR(2)-AVRR(3))/3.D0
      AVR1(1)=(4.D0*AVR1(2)-AVR1(3))/3.D0
      AVR2(1)=(4.D0*AVR2(2)-AVR2(3))/3.D0
C
C     ----- CALCULATE PLASMA SURFACE -----
C
      CALL EQCALF(REDGE,ZAXIS,NTHMAX,RSU,ZSU,IERR)
C
      RGMIN=RSU(1)
      RGMAX=RSU(1)
      ZGMIN=ZSU(1)
      ZGMAX=ZSU(1)
      DO NTH=2,NTHMAX
         RGMIN=MIN(RGMIN,RSU(NTH))
         RGMAX=MAX(RGMAX,RSU(NTH))
         ZGMIN=MIN(ZGMIN,ZSU(NTH))
         ZGMAX=MAX(ZGMAX,ZSU(NTH))
      ENDDO
      DO NTH=1,NTHMAX
         RSW(NTH)=RSU(NTH)
         ZSW(NTH)=ZSU(NTH)
      ENDDO
C
      RETURN
      END
C
C     ***** CALCULATE FLUX VARIABLES IN VACUUM *****
C
      SUBROUTINE EQCALQV(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      IERR=0
C
C     +++++ SETUP VACUUM DATA +++++
C
      DR=RB/(NRMAX-1)
      NRPMAX=NINT(RA/DR)+1
      DR=(RR+RB-REDGE)/(NRMAX-NRPMAX)
      DO NR=NRPMAX+1,NRMAX
         RL=REDGE+DR*(NR-NRPMAX)
         ZL=ZAXIS
         PSS(NR)=PSIG(RL,ZL)
         PPS(NR)=0.D0
         TTS(NR)=2.D0*PI*BB*RR
C
C         FACTOR=(RL-RPS(1,NRPMAX))/(RPS(1,NRPMAX)-RPS(1,NRPMAX-1))
C         WRITE(6,'(A,I5,1P3E12.4)') 'NR,RL,FACTOR,PSS(NR)=',
C     &                               NR,RL,FACTOR,PSS(NR)
C         SPS(NR)=SPS(NRPMAX)+(SPS(NRPMAX)-SPS(NRPMAX-1))*FACTOR
C         VPS(NR)=VPS(NRPMAX)+(VPS(NRPMAX)-VPS(NRPMAX-1))*FACTOR
C         QPS(NR)=QPS(NRPMAX)+(QPS(NRPMAX)-QPS(NRPMAX-1))*FACTOR
C         RLEN(NR)=RLEN(NRPMAX)+(RLEN(NRPMAX)-RLEN(NRPMAX-1))*FACTOR
C         FTS(NR)=FTS(NRPMAX)+(FTS(NRPMAX)-FTS(NRPMAX-1))*FACTOR
C         AVBR(NR)=AVBR(NRPMAX)+(AVBR(NRPMAX)-AVBR(NRPMAX-1))*FACTOR
C         AVRR(NR)=AVRR(NRPMAX)+(AVRR(NRPMAX)-AVRR(NRPMAX-1))*FACTOR
C         AVR1(NR)=AVR1(NRPMAX)+(AVR1(NRPMAX)-AVR1(NRPMAX-1))*FACTOR
C         AVR2(NR)=AVR2(NRPMAX)+(AVR2(NRPMAX)-AVR2(NRPMAX-1))*FACTOR
C
         FACTOR=(RL-RAXIS)/(REDGE-RAXIS)
         SPS(NR)= SPS(NRPMAX)*FACTOR**2
         VPS(NR)= VPS(NRPMAX)*FACTOR**2
         QPS(NR)= QPS(NRPMAX)*FACTOR**2
         RLEN(NR)=RLEN(NRPMAX)*FACTOR
         FTS(NR)= FTS(NRPMAX)*FACTOR**2
         AVBR(NR)=AVBR(NRPMAX)
         AVRR(NR)=AVRR(NRPMAX)
         AVR1(NR)=AVR1(NRPMAX)
         AVR2(NR)=AVR2(NRPMAX)
C
         RMIN=RR
         RMAX=RR
         DO NTH=1,NTHMAX+1
C            DELR=RPS(NTH,NRPMAX)-RPS(NTH,NRPMAX-1)
C            DELZ=ZPS(NTH,NRPMAX)-ZPS(NTH,NRPMAX-1)
C            RPS(NTH,NR)=RPS(NTH,NRPMAX)+FACTOR*DELR
C            ZPS(NTH,NR)=ZPS(NTH,NRPMAX)+FACTOR*DELZ
            RPS(NTH,NR)=RAXIS+(RPS(NTH,NRPMAX)-RAXIS)*FACTOR
            ZPS(NTH,NR)=       ZPS(NTH,NRPMAX)       *FACTOR
            RMIN=MIN(RMIN,RPS(NTH,NR))
            RMAX=MAX(RMAX,RPS(NTH,NR))
         ENDDO
         RRMIN(NR)=RMIN
         RRMAX(NR)=RMAX
         BBMIN(NR)=BBMIN(NRPMAX)/FACTOR
         BBMAX(NR)=BBMAX(NRPMAX)*FACTOR
      ENDDO
      FTSB=FTS(NRMAX)
C
C      ----- CALCULATE PLASMA SURFACE -----
C
      CALL EQCALF(REDGE,ZAXIS,NSUMAX,RSU,ZSU,IERR)
C
C      +++++ CALCULATE WALL DATA +++++
C
      FACTOR=(RB+RR-RAXIS)/(REDGE-RAXIS)
      DO NSU=1,NSUMAX+1
         RSW(NSU)=RAXIS+(RSU(NSU)-RAXIS)*FACTOR
         ZSW(NSU)=       ZSU(NSU)       *FACTOR
      ENDDO
C
      RGMIN=RSW(1)
      RGMAX=RSW(1)
      ZGMIN=ZSW(1)
      ZGMAX=ZSW(1)
      DO NSU=2,NSUMAX
         RGMIN=MIN(RGMIN,RSW(NSU))
         RGMAX=MAX(RGMAX,RSW(NSU))
         ZGMIN=MIN(ZGMIN,ZSW(NSU))
         ZGMAX=MAX(ZGMAX,ZSW(NSU))
      ENDDO
C
      RETURN
      END
C
C     ***** CALCULATE SPLINES AND INTEGRAL QUANTITIES *****
C
      SUBROUTINE EQSETS(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      DIMENSION DERIV(NPSM)
      DIMENSION D01(NTHMP,NRM),D10(NTHMP,NRM),D11(NTHMP,NRM)
C
      IERR=0
      DTH=2*PI/NTHMAX
C
      DO NR=1,NRMAX
         RHOT(NR)=SQRT(FTS(NR)/FTSA)
      ENDDO
C
      CALL SPL1D(PSS,FTS,DERIV,UFTS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for FTS: IERR=',IERR
      CALL SPL1D(FTS,PSS,DERIV,UFTT,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSS: IERR=',IERR
C
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,PPS,DERIV,UPPS,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PPS: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,TTS,DERIV,UTTS,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTS: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,QPS,DERIV,UQPS,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPS: IERR=',IERR
C
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,VPS,DERIV,UVPS,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for VPS: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,SPS,DERIV,USPS,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for SPS: IERR=',IERR
      CALL SPL1D(RHOT,RLEN,DERIV,URLEN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RLEN: IERR=',IERR
C
      CALL SPL1D(RHOT,RRMIN,DERIV,URRMIN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RRMIN: IERR=',IERR
      CALL SPL1D(RHOT,RRMAX,DERIV,URRMAX,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RRMAX: IERR=',IERR
      CALL SPL1D(RHOT,BBMIN,DERIV,UBBMIN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for BBMIN: IERR=',IERR
      CALL SPL1D(RHOT,BBMAX,DERIV,UBBMAX,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for BBMAX: IERR=',IERR
C
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVBR,DERIV,UAVBR,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVBR: IERR=',IERR
      CALL SPL1D(RHOT,AVRR,DERIV,UAVRR,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVRR: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVR1,DERIV,UAVR1,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVR1: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVR2,DERIV,UAVR2,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVR2: IERR=',IERR
C
C        +++++ CALCULATE DERIVATIVES +++++
C     
      DTH=2.D0*PI/NTHMAX
      DO NTH=1,NTHMAX+1
         THIT(NTH)=DTH*(NTH-1)
      ENDDO
C
      CALL SPL2D(THIT,RHOT,RPS,D10,D01,D11,URPS,
     &           NTHMP,NTHMAX+1,NRMAX,4,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for RPS: IERR=',IERR
C
      CALL SPL2D(THIT,RHOT,ZPS,D10,D01,D11,UZPS,
     &           NTHMP,NTHMAX+1,NRMAX,4,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for ZPS: IERR=',IERR
C
      DO NTH=1,NTHMAX+1
         THITL=THIT(NTH)
         DO NR=1,NRMAX
            RHOTL=RHOT(NR)
            CALL SPL2DD(THITL,RHOTL,RPSL,DRCHIL,DRRHOL,
     &                  THIT,RHOT,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
            CALL SPL2DD(THITL,RHOTL,ZPSL,DZCHIL,DZRHOL,
     &                  THIT,RHOT,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
            DRPSI(NTH,NR)=DRRHOL/(2.D0*FTSA)
            DZPSI(NTH,NR)=DZRHOL/(2.D0*FTSA)
            DRCHI(NTH,NR)=DRCHIL
            DZCHI(NTH,NR)=DZCHIL
C
            CALL EQPSID(RPSL,ZPSL,DPSIDR,DPSIDZ)
            BPR(NTH,NR)= DPSIDZ/(2.D0*PI*RPSL)
            BPZ(NTH,NR)=-DPSIDR/(2.D0*PI*RPSL)
            BPT(NTH,NR)=SQRT(BPR(NTH,NR)**2+BPZ(NTH,NR)**2)
            BTP(NTH,NR)= TTS(NR)/(2.D0*PI*RPSL)
C            WRITE(6,'(2I3,1P6E12.4)') 
C     &           NTH,NR,RPSL,ZPSL,
C     &           BPR(NTH,NR),BPZ(NTH,NR),BPT(NTH,NR),BTP(NTH,NR)
         ENDDO
      ENDDO
C
C        +++++ CALCULATE MAGNETIC FIELD +++++
C
C      DO NR=1,NRMAX,5
C         DO NTH=1,NTHMAX
C            WRITE(6,'(2I5,1P3E12.4)') 
C     &           NR,NTH,RPS(NTH,NR),DRPSI(NTH,NR),DRCHI(NTH,NR)
C         ENDDO
C      ENDDO
C
C        +++++ CALCULATE INTEGRATED QUANTITIES +++++
C
      NDPMAX=100
      DELPS=-SAXIS/NDPMAX
      PPSL=FNPPS(0.D0)
      VPSL=FNVPS(0.D0)
      SPSL=FNSPS(0.D0)
      SUMV =0.5D0*VPSL*DELPS
      SUMS =0.5D0*SPSL*DELPS
      SUMPV=0.5D0*PPSL*VPSL*DELPS
      SUMPS=0.5D0*PPSL*SPSL*DELPS
      DO NDP=1,NDPMAX-1
         PSIL=SAXIS+DELPS*NDP
         PSIN=1.D0-PSIL/SAXIS
         PPSL=FNPPS(PSIN)
         VPSL=FNVPS(PSIN)
         SPSL=FNSPS(PSIN)
         SUMV =SUMV +VPSL*DELPS
         SUMS =SUMS +SPSL*DELPS
         SUMPV=SUMPV+PPSL*VPSL*DELPS
         SUMPS=SUMPS+PPSL*SPSL*DELPS
      ENDDO
      PPSL=FNPPS(1.D0)
      VPSL=FNVPS(1.D0)
      SPSL=FNSPS(1.D0)
      SUMV =SUMV +0.5D0*VPSL*DELPS
      SUMS =SUMS +0.5D0*SPSL*DELPS
      SUMPV=SUMPV+0.5D0*PPSL*VPSL*DELPS
      SUMPS=SUMPS+0.5D0*PPSL*SPSL*DELPS
      PVOL=SUMV
      PAREA=SUMS
      RAAVE=SQRT(PAREA/PI)
      PVAVE=SUMPV/SUMV
      PSAVE=SUMPS/SUMS
      BPA=RMU0*RIP*1.D6/FNRLEN(1.D0)
      BETAT=PVAVE/(BB**2/(2.D0*RMU0))
      BETAP=PSAVE/(BPA**2/(2.D0*RMU0))
      QAXIS=FNQPS(0.D0)
      QSURF=FNQPS(1.D0)
C
      IF(NPRINT.GE.2) THEN
         WRITE(6,'(A,1P4E12.4)') 
     &        'PVOL,PAREA,PVAVE,PSAVE  =',PVOL,PAREA,PVAVE,PSAVE
         WRITE(6,'(A,1P4E12.4)') 
     &        'BETAT,BETAP,QAXIS,QSURF =',BETAT,BETAP,QAXIS,QSURF
      ENDIF
C
      RETURN
      END
C
C     ***** CALCULATE FLUX SURFACE *****
C
      SUBROUTINE EQCALF(RINIT,ZINIT,NTHUMAX,RU,ZU,IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      DIMENSION RU(NTHUMAX+1),ZU(NTHUMAX+1)
      DIMENSION XA(NNM),YA(2,NNM)
      DIMENSION RCHI(NNM),ZCHI(NNM),DXCHI(NNM)
      DIMENSION URCHI(4,NNM),UZCHI(4,NNM)
C
      IERR=0
      DTH=2*PI/NTHUMAX
C
C     ----- SET NUMBER OF DIVISION for integration -----
C

      NMAX=200
      IF(NMAX.GT.NNM) NMAX=NNM
C
C     ----- CALCULATE PSI,PPS,TTS,FTS and RPS, ZPS on magnetic surfaces -----
C
      CALL EQMAGS(RINIT,ZINIT,NMAX,XA,YA,NA,IERR)
C
      FACTOR=2.D0*PI/XA(NA)
      RCHI(1)=RINIT
      ZCHI(1)=ZINIT
      XA(1)=FACTOR*XA(1)
      DO N=2,NA-1
         RCHI(N)=YA(1,N)
         ZCHI(N)=YA(2,N)
         XA(N)=FACTOR*XA(N)
      ENDDO
      RCHI(NA)=RCHI(1)
      ZCHI(NA)=ZCHI(1)
      XA(NA)=FACTOR*XA(NA)
C
      CALL SPL1D(XA,RCHI,DXCHI,URCHI,NA,4,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RCHI: IERR=',IERR
      CALL SPL1D(XA,ZCHI,DXCHI,UZCHI,NA,4,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for ZCHI: IERR=',IERR
C
      RU(1)=RINIT
      ZU(1)=ZINIT
      DO NTH=2,NTHUMAX
         TH=DTH*(NTH-1)
         CALL SPL1DF(TH,RU(NTH),XA,URCHI,NA,IERR)
         CALL SPL1DF(TH,ZU(NTH),XA,UZCHI,NA,IERR)
      ENDDO
      RU(NTHUMAX+1)=RINIT
      ZU(NTHUMAX+1)=ZINIT
C
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PP(PSI) *****
C
      FUNCTION PPFUNC(PSIL)
C
      INCLUDE '../eq/eqcomq.inc'
C
C      write(6,*) PSIL
      CALL SPL1DF(PSIL,PPL,PSIPS,UPPPS,NPSMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX PPFUNC: SPL1DF ERROR : IERR=',IERR
C         WRITE(6,*) PSIL,PSIPS(1),PSIPS(NPSMAX)
      ENDIF
      PPFUNC=PPL
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF TT(PSI) *****
C
      FUNCTION TTFUNC(PSIL)
C
      INCLUDE '../eq/eqcomq.inc'
C
C      write(6,*) PSIL
      CALL SPL1DF(PSIL,TTL,PSIPS,UTTPS,NPSMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TTFUNC: SPL1DF ERROR : IERR=',IERR
      TTFUNC=TTL
      RETURN
      END
