C     $Id$
C
C     ***** SETUP DATA *****
C
      SUBROUTINE EQSETP
C
      INCLUDE '../eq/eqcomq.inc'
C
      DIMENSION PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM)
      DIMENSION DERIV(NPSM)
C      DIMENSION GX(101),GY(101,2)
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
C      DR1=0.005D0
C      WRITE(6,'(1P3E14.6)') RAXIS,ZAXIS,SAXIS
C      DO N=1,101
C         R=RR+DR1*(N-51)
C         pp=PSIG(R,ZAXIS)
C         WRITE(6,'(1P3E14.6)') R,ZAXIS,pp
C         GX(N)=GUCLIP(R)
C         GY(N,1)=GUCLIP(pp)
C      ENDDO
C      CALL GSOPEN
C      CALL PAGES
C      CALL EQGR1D(3.0,23.0,2.0,17.0,GX,GY,101,101,1,'@psi-R@',1)
C      CALL PAGEE
C
C      WRITE(6,*) RAXIS,ZAXIS,SAXIS
C
      RETURN
      END
C
C     ***** 
C
      SUBROUTINE EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      EXTERNAL PSIAX,EQDERV
      DIMENSION Y(2),DYDX(2),YOUT(2)
      DIMENSION XA(NNM),YA(2,NNM)
      DIMENSION XCHI0(NNM),XCHI1(NNM),RCHI(NNM),ZCHI(NNM),DXCHI(NNM)
      DIMENSION URCHI(4,NNM),UZCHI(4,NNM)
      DIMENSION DERIV(NPSM)
C
      IERR=0
C
C     ----- CHECK NRMAX,NTHMAX,NSUMAX are not greater than *M -----
C
      NRMAX=NRMAX1
      NTHMAX=NTHMAX1
      NSUMAX=NSUMAX1
C      WRITE(6,*) 'NRMAX,NTHMAX,NSUMAX=',NRMAX,NTHMAX,NSUMAX
      IF(NRMAX.GT.NRM) THEN
         WRITE(6,'(A,2I5)') 'NRMAX.GT.NRM: NRMAX,NRM=',NRMAX,NRM
         IERR=IERR+1
      ENDIF
      IF(NTHMAX.GT.NTHM) THEN
         WRITE(6,'(A,2I5)') 'NTHMAX.GT.NTHM: NTHMAX,NTHM=',NTHMAX,NTHM
         IERR=IERR+2
      ENDIF
      IF(NSUMAX.GT.NSUM) THEN
         WRITE(6,'(A,2I5)') 'NSUMAX.GT.NSUM: NSUMAX,NSUM=',NSUMAX,NSUM
         IERR=IERR+4
      ENDIF
      IF(IERR.NE.0) RETURN
C
C     ----- SET DR, DTH -----
C
      RLIM=RG(NRGMAX)
C
      REDGE=ZBRENT(PSIAX,RR,RLIM,1.D-8)
C
C      SEDGE=PSIG(REDGE,ZAXIS)
C      SLIM =PSIG(RLIM, ZAXIS)
C      WRITE(6,*) REDGE,RLIM,SEDGE,SLIM
C
      RAEQ=(REDGE-RAXIS)
      DR=(RB-RA+REDGE-RAXIS)/(NRMAX-1)
      NRPMAX=INT((REDGE-RAXIS)/DR)+1
      DTH=2*PI/(NTHMAX-1)
C
C     ----- SET NUMBER OF DIVISION for integration -----
C

      NMAX=200
      IF(NMAX.GT.NNM) NMAX=NNM
C
C     ----- CALCULATE PSI,PPS,TTS,FTS and RPS, ZPS on magnetic surfaces -----
C
      NR=1
      DO NTH=1,NTHMAX
         RPS(NTH,NR)=RAXIS
         ZPS(NTH,NR)=ZAXIS
      ENDDO
      PSS(1)=PSIG(RAXIS,ZAXIS)
      PPS(1)=PPFUNC(PSS(1))
      TTS(1)=TTFUNC(PSS(1))
      FTS(1)=0.D0
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
         BMIN=2.D0*BB
         BMAX=0.D0
C         BMIN=ABS(TTS(NR)/RAXIS)
C         BMAX=ABS(TTS(NR)/RAXIS)
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
            BPR=SQRT(DPSIDR**2+DPSIDZ**2)
            B=SQRT((TTS(NR)/R)**2+(BPR/R)**2)
C
            SUMS=SUMS+H/BPR
            SUMV=SUMV+H*R/BPR
            SUMQ=SUMQ+H/(R*BPR)
            XCHI1(N)=SUMQ
            RCHI(N)=YA(1,N)
            ZCHI(N)=YA(2,N)
C
            SUMAVBR=SUMAVBR+H*BPR/R
            SUMAVRR=SUMAVRR+H/(R*BPR)
            SUMAVR1=SUMAVR1+H*R
            SUMAVR2=SUMAVR2+H*R*BPR
            SUML   =SUML   +H*R/BPR
C
            R=YA(1,N)
            Z=YA(2,N)
            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
            BPR=SQRT(DPSIDR**2+DPSIDZ**2)
            B=SQRT((TTS(NR)/R)**2+(BPR/R)**2)
C
            RMIN=MIN(RMIN,R)
            RMAX=MAX(RMAX,Z)
            BMIN=MIN(BMIN,B)
            BMAX=MAX(BMAX,B)
C            IF(NR.EQ.NRPMAX) THEN
C               WRITE(6,'(I5,1P5E12.4)') 
C     &              N,XA(N),YA(1,N),YA(2,N),SUMQ
C            ENDIF
         ENDDO
C
         SPS(NR)=SUMS
         VPS(NR)=SUMV*2.D0*PI
         QPS(NR)=SUMQ*TTS(NR)/(2.D0*PI)
         DO N=1,NA
            XCHI0(N)=2.D0*PI*XA(N)/XA(NA)
            XCHI1(N)=2.D0*PI*XCHI1(N)/SUMQ
         ENDDO
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
C            WRITE(6,'(I5,1P,6E12.4)') 
C     &      NR,SUMAVR1,SUMAVR2,SUML,SUMAVR1/SUML,SUMAVR2/SUML
C         WRITE(6,'(A,I5,1P3E12.4)') 'NR,PSS,AVBR,AVRR=', 
C     &                            NR,PSS(NR)-PSS(1),AVBR(NR),AVRR(NR)
C
C         WRITE(6,'(I5,1P6E12.4/5X,1P5E12.4)') 
C     &        NR,PSS(NR),PPS(NR),TTS(NR),SPS(NR),VPS(NR),QPS(NR),
C     &        RLEN(NR),RRMIN(NR),RRMAX(NR),BBMIN(NR),BBMAX(NR)
C
C        +++++ RECALCULATE MAGSURFACE FOR MESH DATA +++++
C
C         DO N=1,NA
C            WRITE(6,'(A,I5,1P4E12.4)') 'N,X1,X2,R,Z=',
C     &            N,XCHI0(N),XCHI1(N),RCHI(N),ZCHI(N)
C         ENDDO
C         STOP
C
         RCHI(NA)=RCHI(1)
         ZCHI(NA)=ZCHI(1)
C         CALL SPL1D(XCHI0,RCHI,DXCHI,URCHI,NA,4,IERR)
         CALL SPL1D(XCHI1,RCHI,DXCHI,URCHI,NA,4,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RCHI: IERR=',IERR
C         CALL SPL1D(XCHI0,ZCHI,DXCHI,UZCHI,NA,4,IERR)
         CALL SPL1D(XCHI1,ZCHI,DXCHI,UZCHI,NA,4,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for ZCHI: IERR=',IERR
C
         RPS(1,NR)=YA(1,1)
         ZPS(1,NR)=YA(2,1)
         DO NTH=2,NTHMAX
            TH=DTH*(NTH-1)
C            CALL SPL1DF(TH,RPS(NTH,NR),XCHI0,URCHI,NA,IERR)
C            CALL SPL1DF(TH,ZPS(NTH,NR),XCHI0,UZCHI,NA,IERR)
            CALL SPL1DF(TH,RPS(NTH,NR),XCHI1,URCHI,NA,IERR)
            CALL SPL1DF(TH,ZPS(NTH,NR),XCHI1,UZCHI,NA,IERR)
         ENDDO
C         IF(NR.EQ.20) THEN
C            DO NTH=1,NTHMAX
C               TH=DTH*(NTH-1)
C               WRITE(6,'(A,I5,1P3E12.4)') 'NTH,TH,R,Z=',
C     &              NTH,TH,RPS(NTH,NR),ZPS(NTH,NR)
C            ENDDO
C            STOP
C         ENDIF
C         WRITE(6,'(A,I5,1P4E12.4)') 
C     &        'NR,RPS,ZPS,RPS,ZPS=',
C     &        NR,RPS(1,NR),ZPS(1,NR),RPS(NTHMAX,NR),ZPS(NTHMAX,NR)
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
      BBMIN(NR)=ABS(TTS(NR)/RAXIS)
      BBMAX(NR)=ABS(TTS(NR)/RAXIS)
      AVBR(NR)=0.D0
      AVRR(NR)=1.D0/RAXIS**2
      AVR1(NR)=0.D0
      AVR2(NR)=0.D0
C
      FTS(1)=0.D0
      DO NR=2,NRPMAX
         IF(NR.EQ.2) THEN
            FTS(NR)=FTS(NR-1)
     &             +0.5D0*(QPS(NR)+QPS(NR-1))*(PSS(NR)-PSS(NR-1))
C     &             +QPS(NR)*(PSS(NR)-PSS(NR-1))
         ELSE
            FTS(NR)=FTS(NR-1)
     &             +0.5D0*(QPS(NR)+QPS(NR-1))*(PSS(NR)-PSS(NR-1))
         ENDIF
      ENDDO
      FTSA=FTS(NRPMAX)
C
C     +++++ SETUP VACUUM DATA +++++
C
      DO NR=NRPMAX+1,NRMAX
         Y(1)=RAXIS+DR*(NR-1)
         Y(2)=ZAXIS
         PSS(NR)=PSIG(Y(1),Y(2))
         PPS(NR)=0.D0
         TTS(NR)=BB*RR
C
         FACTOR=(Y(1)-RPS(1,NRPMAX))/(RPS(1,NRPMAX)-RPS(1,NRPMAX-1))
         SPS(NR)=SPS(NRPMAX)+(SPS(NRPMAX)-SPS(NRPMAX-1))*FACTOR
         VPS(NR)=VPS(NRPMAX)+(VPS(NRPMAX)-VPS(NRPMAX-1))*FACTOR
         QPS(NR)=QPS(NRPMAX)+(QPS(NRPMAX)-QPS(NRPMAX-1))*FACTOR
         RLEN(NR)=RLEN(NRPMAX)+(RLEN(NRPMAX)-RLEN(NRPMAX-1))*FACTOR
         FTS(NR)=FTS(NRPMAX)+(FTS(NRPMAX)-FTS(NRPMAX-1))*FACTOR
         AVBR(NR)=AVBR(NRPMAX)+(AVBR(NRPMAX)-AVBR(NRPMAX-1))*FACTOR
         AVRR(NR)=AVRR(NRPMAX)+(AVRR(NRPMAX)-AVRR(NRPMAX-1))*FACTOR
         AVR1(NR)=AVR1(NRPMAX)+(AVR1(NRPMAX)-AVR1(NRPMAX-1))*FACTOR
         AVR2(NR)=AVR2(NRPMAX)+(AVR2(NRPMAX)-AVR2(NRPMAX-1))*FACTOR
C
C         WRITE(6,'(I5,1P6E12.4)') 
C     &        NR,PSS(NR),PPS(NR),TTS(NR),SPS(NR),VPS(NR),QPS(NR)
C
         RMIN=RR
         RMAX=RR
         DO NTH=1,NTHMAX
            DELR=RPS(NTH,NRPMAX)-RPS(NTH,NRPMAX-1)
            DELZ=ZPS(NTH,NRPMAX)-ZPS(NTH,NRPMAX-1)
            RPS(NTH,NR)=RPS(NTH,NRPMAX)+FACTOR*DELR
            ZPS(NTH,NR)=ZPS(NTH,NRPMAX)+FACTOR*DELZ
            RMIN=MIN(RMIN,RPS(NTH,NR))
            RMAX=MAX(RMAX,RPS(NTH,NR))
         ENDDO
         RRMIN(NR)=RMIN
         RRMAX(NR)=RMAX
         BBMIN(NR)=BBMIN(NRPMAX)
     &        +(BBMIN(NRPMAX)-BBMIN(NRPMAX-1))*FACTOR
         BBMAX(NR)=BBMAX(NRPMAX)
     &        +(BBMAX(NRPMAX)-BBMAX(NRPMAX-1))*FACTOR
      ENDDO
      FTSB=FTS(NRMAX)
C
C      DO NR=1,NRMAX
C         WRITE(6,'(I5,1P4E12.4)') 
C     &        NR,SQRT(ABS(FTS(NR)/FTSA)),FTS(NR)/FTSA,PSS(NR),FTS(NR)
C      ENDDO
C      PAUSE
C
C        +++++ CALCULATE PLASMA SURFACE DATA +++++
C
      RINIT=REDGE
      ZINIT=ZAXIS
      CALL EQMAGS(RINIT,ZINIT,NMAX,XA,YA,NA,IERR)
C
      DX=XA(NA)/(NSUMAX-1)
      NEQ=2
      NDIV=10
      H=DX/NDIV
      X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
      RSU(1)=Y(1)
      ZSU(1)=Y(2)
      DO NSU=2,NSUMAX-1
         DO J=1,NDIV
            CALL EQDERV(X,Y,DYDX)
            CALL RK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
            X=X+H
            Y(1)=YOUT(1)
            Y(2)=YOUT(2)
         ENDDO
         RSU(NSU)=Y(1)
         ZSU(NSU)=Y(2)
C         WRITE(6,*) X,Y(1),Y(2),PSIG(Y(1),Y(2))
      ENDDO
      RSU(NSUMAX)=RSU(1)
      ZSU(NSUMAX)=ZSU(1)
C
C        +++++ CALCULATE WALL SURFACE DATA +++++
C
      RINIT=REDGE-0.01D0*RR
      ZINIT=ZAXIS
      CALL EQMAGS(RINIT,ZINIT,NMAX,XA,YA,NA,IERR)
C
      DX=XA(NA)/(NSUMAX-1)
      NDIV=10
      H=DX/NDIV
      X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
      RSW(1)=Y(1)
      ZSW(1)=Y(2)
      DO NSU=2,NSUMAX-1
         DO J=1,NDIV
            CALL EQDERV(X,Y,DYDX)
            CALL RK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
            X=X+H
            Y(1)=YOUT(1)
            Y(2)=YOUT(2)
         ENDDO
         RSW(NSU)=Y(1)
         ZSW(NSU)=Y(2)
      ENDDO
      RSW(NSUMAX)=RSW(1)
      ZSW(NSUMAX)=ZSW(1)
C
C        +++++ CALCULATE WALL DATA +++++
C
      FACTOR=(RPS(1,NRMAX)-RPS(1,NRPMAX))/(RSU(1)-RSW(1))
      DO NSU=1,NSUMAX
         DELR=RSU(NSU)-RSW(NSU)
         DELZ=ZSU(NSU)-ZSW(NSU)
         RSW(NSU)=RSU(NSU)+FACTOR*DELR
         ZSW(NSU)=ZSU(NSU)+FACTOR*DELZ
      ENDDO
C
      RGMIN=RSW(1)
      RGMAX=RSW(1)
      DO NSU=2,NSUMAX
         RGMIN=MIN(RGMIN,RSW(NSU))
         RGMAX=MAX(RGMAX,RSW(NSU))
      ENDDO
      ZGMIN=ZSW(1)
      ZGMAX=ZSW(1)
      DO NSU=2,NSUMAX
         ZGMIN=MIN(ZGMIN,ZSW(NSU))
         ZGMAX=MAX(ZGMAX,ZSW(NSU))
      ENDDO
C
C        +++++ CALCULATE SPLINE COEFFICIENTS +++++
C
C      WRITE(6,'(A,1P,2E12.4)') 'RA,RAEQ=',RA,RAEQ
      DO NR=2,NRMAX
         IF(NR.EQ.NRMAX) THEN
            DPSIDRHO=(PSS(NR  )-PSS(NR-1))/(DR/RAEQ)
         ELSE
            DPSIDRHO=(PSS(NR+1)-PSS(NR-1))/(2*DR/RAEQ)
         ENDIF
         AVBR(NR)=AVBR(NR)/DPSIDRHO**2
         AVR1(NR)=AVR1(NR)/DPSIDRHO
         AVR2(NR)=AVR2(NR)/DPSIDRHO**2
      ENDDO
      AVBR(1)=(4.D0*AVBR(2)-AVBR(3))/3.D0
      AVR1(1)=(4.D0*AVR1(2)-AVR1(3))/3.D0
      AVR2(1)=(4.D0*AVR2(2)-AVR2(3))/3.D0
C
      IF(NPRINT.GE.3) THEN
         DO NR=1,NRMAX
            WRITE(6,'(I5,1P,5E12.4)') 
     &      NR,PSS(NR)-PSS(1),RLEN(NR),AVR1(NR),AVR2(NR),
     &           AVRR(NR)
         ENDDO
      ENDIF
C
      CALL SPL1D(PSS,PPS,DERIV,UPPS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PPS: IERR=',IERR
      CALL SPL1D(PSS,TTS,DERIV,UTTS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTS: IERR=',IERR
      CALL SPL1D(PSS,QPS,DERIV,UQPS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPS: IERR=',IERR
      CALL SPL1D(PSS,FTS,DERIV,UFTS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for FTS: IERR=',IERR
      CALL SPL1D(FTS,PSS,DERIV,UFTT,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSS: IERR=',IERR
      CALL SPL1D(PSS,VPS,DERIV,UVPS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for VPS: IERR=',IERR
      CALL SPL1D(PSS,SPS,DERIV,USPS,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for SPS: IERR=',IERR
      CALL SPL1D(PSS,RLEN,DERIV,URLEN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RLEN: IERR=',IERR
      CALL SPL1D(PSS,RRMIN,DERIV,URRMIN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RRMIN: IERR=',IERR
      CALL SPL1D(PSS,RRMAX,DERIV,URRMAX,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for RRMAX: IERR=',IERR
      CALL SPL1D(PSS,BBMIN,DERIV,UBBMIN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for BBMIN: IERR=',IERR
      CALL SPL1D(PSS,BBMAX,DERIV,UBBMAX,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for BBMAX: IERR=',IERR
      CALL SPL1D(PSS,AVBR,DERIV,UAVBR,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVBR: IERR=',IERR
      CALL SPL1D(PSS,AVRR,DERIV,UAVRR,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVRR: IERR=',IERR
      CALL SPL1D(PSS,AVR1,DERIV,UAVR1,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVR1: IERR=',IERR
      CALL SPL1D(PSS,AVR2,DERIV,UAVR2,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVR2: IERR=',IERR
C
C        +++++ CALCULATE DERIVATIVES +++++
C     
      DO NTH=1,NTHMAX
         NR=1
            DRPSI(NTH,NR)=0.D0
            DZPSI(NTH,NR)=0.D0
         DO NR=2,NRMAX-1
            DPSI=PSS(NR+1)-PSS(NR-1)
            DRPSI(NTH,NR)=(RPS(NTH,NR+1)-RPS(NTH,NR-1))/DPSI
            DZPSI(NTH,NR)=(ZPS(NTH,NR+1)-ZPS(NTH,NR-1))/DPSI
         ENDDO
         NR=NRMAX
            DPSI=PSS(NR)-PSS(NR-1)
            DRPSI(NTH,NR)=(RPS(NTH,NR)-RPS(NTH,NR-1))/DPSI
            DZPSI(NTH,NR)=(ZPS(NTH,NR)-ZPS(NTH,NR-1))/DPSI
      ENDDO
C
      DO NR=1,NRMAX
         PSIN=1.D0-PSS(NR)/SAXIS
         IF(PSIN.GT.0.D0) THEN
            FLTT=FNTTS(PSIN)
            FLQP=FNQPS(PSIN)
            FLRL=FNRLEN(PSIN)
            FACTOR=FLTT*FLRL/(2.D0*PI*FLQP)
         ENDIF
C         WRITE(6,'(1P5E12.4)') PSIN,FLTT,FLQP,FLRL,FACTOR
C
         NTH=1
            IF(PSIN.GT.0.D0) THEN
               CALL EQPSID(RPS(NTH,NR),ZPS(NTH,NR),DPSIDR,DPSIDZ)
               DPSI=SQRT(DPSIDR**2+DPSIDZ**2)
               DCHI=FACTOR*DTH/(RPS(NTH,NR)*DPSI)
            ELSE
               DCHI=DTH
            ENDIF
            DRCHI(NTH,NR)=(RPS(NTH+1,NR)-RPS(NTHMAX,NR))/(2*DCHI)
            DZCHI(NTH,NR)=(ZPS(NTH+1,NR)-ZPS(NTHMAX,NR))/(2*DCHI)
         DO NTH=2,NTHMAX-1
            IF(PSIN.GT.0.D0) THEN
               CALL EQPSID(RPS(NTH,NR),ZPS(NTH,NR),DPSIDR,DPSIDZ)
               DPSI=SQRT(DPSIDR**2+DPSIDZ**2)
               DCHI=FACTOR*DTH/(RPS(NTH,NR)*DPSI)
            ELSE
               DCHI=DTH
            ENDIF
            DRCHI(NTH,NR)=(RPS(NTH+1,NR)-RPS(NTH-1,NR))/(2*DCHI)
            DZCHI(NTH,NR)=(ZPS(NTH+1,NR)-ZPS(NTH-1,NR))/(2*DCHI)
         ENDDO
         NTH=NTHMAX
            IF(PSIN.GT.0.D0) THEN
               CALL EQPSID(RPS(NTH,NR),ZPS(NTH,NR),DPSIDR,DPSIDZ)
               DPSI=SQRT(DPSIDR**2+DPSIDZ**2)
               DCHI=FACTOR*DTH/(RPS(NTH,NR)*DPSI)
            ELSE
               DCHI=DTH
            ENDIF
            DRCHI(NTH,NR)=(RPS(1,NR)-RPS(NTH-1,NR))/(2*DCHI)
            DZCHI(NTH,NR)=(ZPS(1,NR)-ZPS(NTH-1,NR))/(2*DCHI)
      ENDDO
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
C     ***** 
C
      SUBROUTINE EQMAGS(RINIT,ZINIT,NMAX,XA,YA,N,IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      EXTERNAL EQDERV
      DIMENSION Y(2),DYDX(2),YOUT(2)
      DIMENSION XA(NMAX),YA(2,NMAX)
C
      X=0.D0
      Y(1)=RINIT
      Y(2)=ZINIT
      N=1
      XA(N)=X
      YA(1,N)=Y(1)
      YA(2,N)=Y(2)
      NEQ=2
C      H=4.D0*PI*(Y(1)-RAXIS)/NMAX
      H=1.2D0*2.D0*PI*RKAP*(Y(1)-RAXIS)/NMAX
      IMODE=0
      DO I=2,NMAX
         CALL EQDERV(X,Y,DYDX)
         CALL RK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
         IF(IMODE.EQ.0) THEN
            IF(YOUT(2).LT.ZINIT) IMODE=1
         ELSE
            IF(YOUT(2).GT.ZINIT) GOTO 1000
         ENDIF
         X=X+H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
C         WRITE(6,'(I5,1P3E12.4)') N,X,Y(1),Y(2)
         N=N+1
         XA(N)=X
         YA(1,N)=Y(1)
         YA(2,N)=Y(2)
      ENDDO
      WRITE(6,*) 'XX RK4: NOT ENOUGH N'
      IERR=1
      RETURN
C
 1000 CONTINUE
      H=0.1D0*H
      DO I=1,11
         CALL EQDERV(X,Y,DYDX)
         CALL RK4(X,Y,DYDX,YOUT,H,NEQ,EQDERV)
         IF(YOUT(2).GT.ZINIT) GOTO 2000
         X=X+H
         Y(1)=YOUT(1)
         Y(2)=YOUT(2)
      ENDDO
      WRITE(6,*) 'XX RK4: UNEXPECTED BEHAVIOR'
      IERR=2
      RETURN
C
 2000 CONTINUE
      DEL=(ZINIT-Y(2))/(YOUT(2)-Y(2))
      X=X+H*DEL
      Y(1)=Y(1)+(YOUT(1)-Y(1))*DEL
      Y(2)=Y(2)+(YOUT(2)-Y(2))*DEL
      N=N+1
      XA(N)=X
      YA(1,N)=Y(1)
      YA(2,N)=Y(2)
      IERR=0
C
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTIONS *****
C
      FUNCTION FNPSIN(RHON)
C
      INCLUDE '../eq/eqcomq.inc'
C
      FTL=FTSA*RHON*RHON
      CALL SPL1DF(FTL,PSIL,FTS,UFTT,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPSIN: SPL1DF ERROR : IERR=',IERR
      FNPSIN=1.D0-PSIL/SAXIS
      RETURN
      END
C
      FUNCTION FNFTS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,FTL,PSS,UFTS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNFTS: SPL1DF ERROR : IERR=',IERR
      FNFTS=FTL
      RETURN
      END
C
      FUNCTION FNPPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,PPL,PSS,UPPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNPPS: SPL1DF ERROR : IERR=',IERR
      FNPPS=PPL
      RETURN
      END
C
      FUNCTION FNTTS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      IF(PSIN.LT.1.D0) THEN
         PSIL=SAXIS*(1.D0-PSIN)
C         write(6,*) PSIL
         CALL SPL1DF(PSIL,TTL,PSS,UTTS,NRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 
     &        'XX FNTTS: SPL1DF ERROR : IERR=',IERR
      ELSE
         TTL=RR*BB
      ENDIF
      FNTTS=TTL
      RETURN
      END
C
      FUNCTION FNQPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
C      write(6,*) PSIL
      CALL SPL1DF(PSIL,QPL,PSS,UQPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNQPS: SPL1DF ERROR : IERR=',IERR
      FNQPS=QPL
      RETURN
      END
C
      FUNCTION FNVPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,VPL,PSS,UVPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNVPS: SPL1DF ERROR : IERR=',IERR
      FNVPS=VPL
      RETURN
      END
C
      FUNCTION FNSPS(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,SPL,PSS,USPS,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNSPS: SPL1DF ERROR : IERR=',IERR
      FNSPS=SPL
      RETURN
      END
C
      FUNCTION FNRLEN(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
C      write(6,*) PSIL
      CALL SPL1DF(PSIL,RLENL,PSS,URLEN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRLEN: SPL1DF ERROR : IERR=',IERR
      FNRLEN=RLENL
      RETURN
      END
C
      FUNCTION FNRRMN(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,RRMINL,PSS,URRMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMN: SPL1DF ERROR : IERR=',IERR
      FNRRMN=RRMINL
      RETURN
      END
C
      FUNCTION FNRRMX(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,RRMAXL,PSS,URRMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNRRMX: SPL1DF ERROR : IERR=',IERR
      FNRRMX=RRMAXL
      RETURN
      END
C
      FUNCTION FNBBMN(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,BBMINL,PSS,UBBMIN,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNBBMN: SPL1DF ERROR : IERR=',IERR
      FNBBMN=BBMINL
      RETURN
      END
C
      FUNCTION FNBBMX(PSIN)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIL=SAXIS*(1.D0-PSIN)
      CALL SPL1DF(PSIL,BBMAXL,PSS,UBBMAX,NRMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX FNBBMX: SPL1DF ERROR : IERR=',IERR
      FNBBMX=BBMAXL
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PSI(R,Z) *****
C
      FUNCTION PSIG(R,Z)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL2DF(R,Z,PSIL,RG,ZG,URZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX PSIG: SPL2DF ERROR : IERR=',IERR
      PSIG=PSIL
      RETURN
      END
C
C     ***** INTERPOLATE SUBROUTINE PSI,DPSIDR,DPSIDZ(R,Z) *****
C
      SUBROUTINE  EQPSID(R,Z,DPSIDR,DPSIDZ)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL2DD(R,Z,PSI,DPSIDR,DPSIDZ,
     &            RG,ZG,URZ,NRGM,NRGMAX,NZGMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQPSID: SPL2DD ERROR : IERR=',IERR
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PSI(R,ZAXIS) *****
C
      FUNCTION PSIAX(R)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PSIAX=PSIG(R,ZAXIS)
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
         WRITE(6,*) PSIL,PSIPS(1),PSIPS(NPSMAX)
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
C
C     ***** DERIVATIVES *****
C
      SUBROUTINE EQDERV(X,Y,DYDX)
C
      INCLUDE '../eq/eqcomq.inc'
      DIMENSION Y(2),DYDX(2)
C
      CALL EQPSID(Y(1),Y(2),PSIRL,PSIZL)
C
      PSID=SQRT(PSIRL**2+PSIZL**2)
C
      DYDX(1)=-PSIZL/PSID
      DYDX(2)= PSIRL/PSID
C      WRITE(6,'(1P5E12.4)') X,Y(1),Y(2),DYDX(1),DYDX(2)
      RETURN
      END
C
C     ****** SIMPLE RUNGE-KUTTA METHOD ******
C
      SUBROUTINE RK4(X,Y,DYDX,YOUT,H,N,DERIVS)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N,NMAX
      DIMENSION DYDX(N),Y(N),YOUT(N)
      EXTERNAL DERIVS
      PARAMETER (NMAX=50)
      INTEGER I
      REAL*8 H6,HH,XH,DYM(NMAX),DYT(NMAX),YT(NMAX)
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      DO 11 I=1,N
        YT(I)=Y(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
        YT(I)=Y(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(X+H,YT,DYT)
      DO 14 I=1,N
        YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.D0*DYM(I))
14    CONTINUE
      RETURN
      END
C
C     ****** TWO-DIMENSIONAL NEWTON METHOD ******
C
      SUBROUTINE NEWTN(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,IER)
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SUB
C
      IER=0
      ITER=0
      IF(ABS(X).GT.1.D0) THEN
         HX=DELT*X
      ELSE
         HX=DELT
      ENDIF
      IF(ABS(Y).GT.1.D0) THEN
         HY=DELT*Y
      ELSE
         HY=DELT
      ENDIF
      CALL SUB(X,   Y,   S0,G0)
      IF(LIST.GT.0) WRITE(6,600) X,Y,S0,G0
C      WRITE(6,*) X,Y,S0,G0
      IF(LIST.GT.0) WRITE(6,*) 'HX=',HX,'HY=',HY
      IF(LIST.GT.0) CALL GUFLSH
    1 CALL SUB(X+HX,Y,   SX,GX)
      CALL SUB(X,   Y+HY,SY,GY)
      FXX =(SX-S0)/HX
      FXY1=(GX-G0)/HX
      FYY =(GY-G0)/HY
      FXY2=(SY-S0)/HY
      FXY=0.5D0*(FXY1+FXY2)
      IF(LIST.GT.1) WRITE(6,601) SX,SY,GX,GY
      IF(LIST.GT.1) WRITE(6,601) FXX,FYY,FXY1,FXY2
      IF(LIST.GT.1) CALL GUFLSH
      DF=SQRT(S0*S0+G0*G0)
      DET=FXX*FYY-FXY*FXY
      H11= FYY/DET
      H12=-FXY/DET
      H21=-FXY/DET
      H22= FXX/DET
C
      DX=-(H11*S0+H12*G0)
      DY=-(H21*S0+H22*G0)
      TT=1.0D0
    2 X=X+TT*DX
      Y=Y+TT*DY
      CALL SUB(X,   Y,   SN,GN)
      IF(LIST.GT.0) WRITE(6,600) DX,DY
      IF(LIST.GT.0) WRITE(6,600) X,Y,SN,GN
      IF(LIST.GT.0) CALL GUFLSH
      DFN=SQRT(SN*SN+GN*GN)
      IF(DFN.GT.DF) THEN
         X=X-TT*DX
         Y=Y-TT*DY
         TT=0.5D0*TT
         ITER=ITER+1
         IF(TT.LE.1.D-3) GOTO 8000
         IF(ITER.LE.ILMAX) GOTO 2
      ELSE
         S0=SN
         G0=GN
         DF=DFN
         ITER=ITER+1
         IF(DF.LE.EPS) GO TO 9000
         IF(ITER.LE.ILMAX) GO TO 1
      ENDIF
C
      IER=2
      IF(LIST.GT.0)
     &WRITE(6,*) 'XX NEWTN: LOOP COUNT EXCEEDS UPPER BOUND.'
      GOTO 9000
C
 8000 IER=1
      IF(LIST.GT.0)
     &WRITE(6,*) 'XX NEWTN: DOES NOT CONVERGE.'
      GOTO 9000
C
 9000 XX=X
      YY=Y
      RETURN
  600 FORMAT(" ",6X,'X,Y,FX,FY = ',1P4E15.7)
  601 FORMAT(" ",6X,'FXX,YY,XY = ',1P4E15.7)
      END
