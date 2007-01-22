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
      CALL EQCALQP(IERR)
C
      IF(NSUMAX.GT.0) CALL EQCALQV(IERR)
C
      CALL EQSETS(IERR)
C
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
      DIMENSION HJTRG(NRGM,NZGM),HJTZG(NRGM,NZGM),HJTRZG(NRGM,NZGM)
      DIMENSION DERIV(NPSM)
      EXTERNAL EQPSID
C
      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,UPSIRZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for PSIRZ: IERR=',IERR
C
      CALL SPL2D(RG,ZG,HJTRZ,HJTRG,HJTZG,HJTRZG,UHJTRZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for HJTRZ: IERR=',IERR
C
C      DO NPS=1,NPSMAX
C         WRITE(6,'(A,I5,1P3E12.4)') 'NPS:',NPS,PSIPS(NPS),
C     &                             PPPS(NPS),TTPS(NPS)
C      ENDDO
C
      CALL SPL1D(PSIPS,PPPS,  DERIV,UPPPS, NPSMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PPPS: IERR=',IERR
      CALL SPL1D(PSIPS,TTPS,  DERIV,UTTPS, NPSMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTPS: IERR=',IERR
C
      RAXIS=RR
      ZAXIS=0.D0
      PSI0=PSIG(RAXIS,ZAXIS)
      PSIPA=-PSI0
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
      EXTERNAL EQDERV
      DIMENSION XA(NTVM),YA(2,NTVM)
      DIMENSION XCHI0(NTVM),XCHI1(NTVM)
      DIMENSION RCHI(NTVM),ZCHI(NTVM),DXCHI(NTVM)
      DIMENSION URCHI(4,NTVM),UZCHI(4,NTVM)
C
      IERR=0
C
      CALL EQAXIS(IERR)
      IF(IERR.NE.0) RETURN
C
C     ----- SET DR, DTH -----
C
      IF(NSUMAX.EQ.0) THEN
         NRPMAX=NRMAX
      ELSE
         DR=(RB-RA+REDGE-RAXIS)/(NRMAX-1)
         NRPMAX=NINT((REDGE-RAXIS)/DR)+1
      ENDIF
      DR=(REDGE-RAXIS)/(NRPMAX-1)
      DTH=2*PI/NTHMAX
C
C     ----- SET NUMBER OF DIVISION for integration -----
C
      IF(NTVMAX.GT.NTVM) NTVMAX=NTVM
C
C     ----- CALCULATE PSI,PPS,TTS,PSIT and RPS, ZPS on magnetic surfaces -----
C
      NR=1
      DO NTH=1,NTHMAX+1
         RPS(NTH,NR)=RAXIS
         ZPS(NTH,NR)=ZAXIS
      ENDDO
      PSIP(1)=PSIG(RAXIS,ZAXIS)-PSI0
      PPS(1)=PPFUNC(PSIP(1))
      TTS(1)=TTFUNC(PSIP(1))
C
C         WRITE(6,'(A,I5,1P3E12.4)') 'NR:',NR,
C     &        PSIP(NR),PPS(NR),TTS(NR)
C
      DO NR=2,NRPMAX
         RINIT=RAXIS+DR*(NR-1)
         ZINIT=ZAXIS
         PSIP(NR)=PSIG(RINIT,ZINIT)-PSI0
         PPS(NR)=PPFUNC(PSIP(NR))
         TTS(NR)=TTFUNC(PSIP(NR))
C
C         WRITE(6,'(A,I5,1P3E12.4)') 'NR:',NR,
C     &        PSIP(NR),PPS(NR),TTS(NR)
C
         CALL EQMAGS(RINIT,ZINIT,NTVMAX,XA,YA,NA,IERR)
C
         SUMS=0.D0
         SUMV=0.D0
         SUMQ=0.D0
         SUMAVRHR=0.D0
         SUMAVIR2=0.D0
         SUMAVRH1=0.D0
         SUMAVRH2=0.D0
         SUMAVBB2=0.D0
         SUMAVIB2=0.D0
         SUMAVRHB=0.D0
C
         RMIN=RAXIS
         RMAX=RAXIS
         ZMIN=ZAXIS
         ZMAX=ZAXIS
         BMIN=ABS(2.D0*BB)
         BMAX=0.D0
C
         XCHI0(1)=0.D0
         XCHI1(1)=0.D0
         RCHI(1)=RINIT
         ZCHI(1)=ZINIT
C
         DO N=2,NA
            H=XA(N)-XA(N-1)
            R=0.5D0*(YA(1,N-1)+YA(1,N))
            Z=0.5D0*(YA(2,N-1)+YA(2,N))
            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
            BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
            BTL=TTS(NR)/(2.D0*PI*R)
            B2L=BTL**2+BPL**2
C
            SUMS=SUMS+H/(BPL*R)
            SUMV=SUMV+H/BPL
            SUMQ=SUMQ+H/(BPL*R*R)
C
            SUMAVRHR=SUMAVRHR+H*BPL
            SUMAVIR2=SUMAVIR2+H/(BPL*R*R)
            SUMAVRH1=SUMAVRH1+H*R
            SUMAVRH2=SUMAVRH2+H*R*R*BPL
            SUMAVBB2=SUMAVBB2+H*B2L/BPL
            SUMAVIB2=SUMAVIB2+H/(B2L*BPL)
            SUMAVRHB=SUMAVRHB+H*R*R*BPL/B2L
C
            XCHI1(N)=SUMAVIR2
            RCHI(N)=YA(1,N)
            ZCHI(N)=YA(2,N)
C
            R=YA(1,N)
            Z=YA(2,N)
            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
            BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
            BTL=TTS(NR)/(2.D0*PI*R)
            B2L=BTL**2+BPL**2
            B=SQRT(B2L)
C
            RMIN=MIN(RMIN,R)
            RMAX=MAX(RMAX,R)
            ZMIN=MIN(ZMIN,Z)
            ZMAX=MAX(ZMAX,Z)
            BMIN=MIN(BMIN,B)
            BMAX=MAX(BMAX,B)
         ENDDO
C
         SPS(NR)=SUMS/(2.D0*PI)
         VPS(NR)=SUMV
         QPS(NR)=SUMAVIR2*TTS(NR)/(4.D0*PI**2)
         RLEN(NR)=XA(NA)
         RRMIN(NR)=RMIN
         RRMAX(NR)=RMAX
         ZZMIN(NR)=ZMIN
         ZZMAX(NR)=ZMAX
         BBMIN(NR)=BMIN
         BBMAX(NR)=BMAX
         AVERHR(NR)=SUMAVRHR/SUMV
         AVEIR2(NR)=SUMAVIR2/SUMV
         AVERH1(NR)=SUMAVRH1/SUMV
         AVERH2(NR)=SUMAVRH2/SUMV
         AVEBB2(NR)=SUMAVBB2/SUMV
         AVEIB2(NR)=SUMAVIB2/SUMV
         AVERHB(NR)=SUMAVRHB/SUMV
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
      ZZMIN(NR)=ZAXIS
      ZZMAX(NR)=ZAXIS
      BBMIN(NR)=ABS(TTS(NR)/(2.D0*PI*RAXIS))
      BBMAX(NR)=ABS(TTS(NR)/(2.D0*PI*RAXIS))
C
C     ----- CALCULATE TOROIDAL FLUX -----
C
      PSIT(1)=0.D0
      DO NR=2,NRPMAX
         PSIT(NR)=PSIT(NR-1)
     &           +2.0D0*QPS(NR)*QPS(NR-1)/(QPS(NR)+QPS(NR-1))
     &                 *(PSIP(NR)-PSIP(NR-1))
      ENDDO
      PSITA=PSIT(NRPMAX)
      PSIPA=PSIP(NRPMAX)
C
      RST(1)=0.D0
      DO NR=2,NRPMAX
         RST(NR)=SQRT(PSIT(NR)/(PI*BB))
      ENDDO
      RSTA=RST(NRPMAX)
C
      RHOT(1)=0.D0
      DO NR=2,NRPMAX
         RHOT(NR)=SQRT(PSIT(NR)/PSITA)
      ENDDO
C
C     ----- CALCULATE EDGE VALUE -----
C
      RINIT=REDGE
      ZINIT=ZAXIS
      TTSA=TTFUNC(0.D0)
      CALL EQMAGS(RINIT,ZINIT,NTVMAX,XA,YA,NA,IERR)
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
C
      DO NR=2,NRPMAX
         AVERHR(NR)=AVERHR(NR)*QPS(NR)**2/(RST(NR)**2*BB**2)/RSTA**2
         AVEIR2(NR)=AVEIR2(NR)*RR**2
         AVERH1(NR)=AVERH1(NR)*QPS(NR)/(RST(NR)*BB)/RSTA
         AVERH2(NR)=AVERH2(NR)*QPS(NR)**2/(RST(NR)**2*BB**2)/RSTA**2
         AVEBB2(NR)=AVEBB2(NR)/BB**2
         AVEIB2(NR)=AVEIB2(NR)*BB**2
         AVERHB(NR)=AVERHB(NR)*QPS(NR)**2/(RST(NR)**2*RR**2)/RSTA**2
         DPPSL=DPPFUNC(PSIP(NR))
         DTTSL=DTTFUNC(PSIP(NR))
         TTSL= TTFUNC(PSIP(NR))
         AVEJPR(NR)=-(TTSL*DPPSL/BB+AVEBB2(NR)*BB*DTTSL/RMU0)
         AVEJTR(NR)=-(2.D0*PI*RR*DPPSL
     &               +AVEIR2(NR)*TTSL*DTTSL/(2.D0*PI*RMU0*RR))
C         WRITE(6,'(A,1P6E12.4)') 'AVEJ=',
C     &        TTSL*DPPSL/BB,AVEBB2(NR)*BB*DTTSL/RMU0,
C     &        2.D0*PI*RR*DPPSL,AVEIR2(NR)*TTSL*DTTSL/(2.D0*PI*RMU0*RR),
C     &        AVEJPR(NR),AVEJTR(NR)
      ENDDO
      AVERHR(1)=(4.D0*AVERHR(2)-AVERHR(3))/3.D0
      AVEIR2(1)=(4.D0*AVEIR2(2)-AVEIR2(3))/3.D0
      AVERH1(1)=(4.D0*AVERH1(2)-AVERH1(3))/3.D0
      AVERH2(1)=(4.D0*AVERH2(2)-AVERH2(3))/3.D0
      AVEBB2(1)=(4.D0*AVEBB2(2)-AVEBB2(3))/3.D0
      AVEIB2(1)=(4.D0*AVEIB2(2)-AVEIB2(3))/3.D0
      AVERHB(1)=(4.D0*AVERHB(2)-AVERHB(3))/3.D0
      AVEJPR(1)=(4.D0*AVEJPR(2)-AVEJPR(3))/3.D0
      AVEJTR(1)=(4.D0*AVEJTR(2)-AVEJTR(3))/3.D0
C
C     ----- CALCULATE PLASMA SURFACE -----
C
      CALL EQCALF(REDGE,ZAXIS,NTHMAX,RSU,ZSU,IERR)
C
      IF(MDLEQF.LT.10) THEN
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
      ENDIF
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
      EXTERNAL EQDERV
      integer np,i
      parameter (np=3)
      real*8 nra,rpspa,zpspa,dy,spspa,vpspa,qpspa,rlenpa,psitpa
      real*8 averhrpa,aveir2pa,averh1pa,averh2pa,avebb2pa,aveib2pa
      real*8 averhbpa
      DIMENSION XA(NTVM),YA(2,NTVM)
      DIMENSION XCHI0(NTVM),XCHI1(NTVM)
      DIMENSION RCHI(NTVM),ZCHI(NTVM),DXCHI(NTVM)
      DIMENSION URCHI(4,NTVM),UZCHI(4,NTVM)
      dimension nra(np),rpspa(np),zpspa(np),spspa(np),vpspa(np)
      dimension qpspa(np),rlenpa(np),psitpa(np),averhrpa(np)
      dimension aveir2pa(np),averh1pa(np),averh2pa(np),avebb2pa(np)
      dimension aveib2pa(np),averhbpa(np)
      DIMENSION THW(NTHMP)
      DIMENSION RPSW(NTHMP),DRPSW(NTHMP),URPSW(4,NTHMP)
      DIMENSION ZPSW(NTHMP),DZPSW(NTHMP),UZPSW(4,NTHMP)
C
      IERR=0
C
C     +++++ SETUP VACUUM DATA +++++
C
      DR=(REDGE-RAXIS)/(NRPMAX-1)
      REDGEM=REDGE-DR
      DR=(RB-RA+REDGE-RAXIS)/(NRMAX-1)
      NRPMAX=NINT((REDGE-RAXIS)/DR)+1
      DR=(RB-RA)/(NRMAX-NRPMAX)
      IF(MDLEQF.LT.10) THEN
C         REDGEM=RA-RA/(NRPMAX-1)
         DO NR=NRPMAX+1,NRMAX
            RL=REDGE+DR*(NR-NRPMAX)
            ZL=ZAXIS
            PSIP(NR)=PSIG(RL,ZL)-PSI0
            PPS(NR)=0.D0
            TTS(NR)=2.D0*PI*BB*RR
C
C            FACTOR=(RL-REDGE)/(REDGE-REDGEM)
            FACTOR=(RL-RAXIS)/(REDGE-RAXIS)
c            SPS(NR)= SPS(NRPMAX)*FACTOR**2
            do i=1,np
               nra(i)=nr-np-1+i
               spspa(i)=sps(nr-np-1+i)
            enddo
c     
            call polint(nra,spspa,np,nr,sps(nr),dy) 
c
c

c            VPS(NR)= VPS(NRPMAX)*FACTOR**2
               do i=1,np
                  nra(i)=nr-np-1+i
                  vpspa(i)=vps(nr-np-1+i)
               enddo
c
                  call polint(nra,vpspa,np,nr,vps(nr),dy)
c
c
c            QPS(NR)= QPS(NRPMAX)*FACTOR**2
               do i=1,np
                  nra(i)=nr-np-1+i
                  qpspa(i)=qps(nr-np-1+i)
               enddo
c
                  call polint(nra,qpspa,np,nr,qps(nr),dy)
c
c
c            RLEN(NR)=RLEN(NRPMAX)*FACTOR
               do i=1,np
                  nra(i)=nr-np-1+i
                  rlenpa(i)=rlen(nr-np-1+i)
               enddo
c
                  call polint(nra,rlenpa,np,nr,rlen(nr),dy)
c
c
c            PSIT(NR)= PSIT(NRPMAX)*FACTOR**2
               do i=1,np
                  nra(i)=nr-np-1+i
                  psitpa(i)=psit(nr-np-1+i)
               enddo
c
                  call polint(nra,psitpa,np,nr,psit(nr),dy)
c
c
c            AVERHR(NR)=AVERHR(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  averhrpa(i)=averhr(nr-np-1+i)
               enddo
c
                  call polint(nra,averhrpa,np,nr,averhr(nr),dy)
c
c
c            AVEIR2(NR)=AVEIR2(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  aveir2pa(i)=aveir2(nr-np-1+i)
               enddo
c
                  call polint(nra,aveir2pa,np,nr,aveir2(nr),dy)
c
c
c            AVERH1(NR)=AVERH1(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  averh1pa(i)=averh1(nr-np-1+i)
               enddo
c
                  call polint(nra,averh1pa,np,nr,averh1(nr),dy)
c
c
c            AVERH2(NR)=AVERH2(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  averh2pa(i)=averh2(nr-np-1+i)
               enddo
c
                  call polint(nra,averh2pa,np,nr,averh2(nr),dy)
c
c
c            AVEBB2(NR)=AVEBB2(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  avebb2pa(i)=avebb2(nr-np-1+i)
               enddo
c
                  call polint(nra,avebb2pa,np,nr,avebb2(nr),dy)
c
c
c            AVEIB2(NR)=AVEIB2(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  aveib2pa(i)=aveib2(nr-np-1+i)
               enddo
c
                  call polint(nra,aveib2pa,np,nr,aveib2(nr),dy)
c
c
c            AVERHB(NR)=AVERHB(NRPMAX)
               do i=1,np
                  nra(i)=nr-np-1+i
                  averhbpa(i)=averhb(nr-np-1+i)
               enddo
c
                  call polint(nra,averhbpa,np,nr,averhb(nr),dy)
c
c
C
            RMIN=RAXIS
            RMAX=RAXIS
            ZMIN=ZAXIS
            ZMAX=ZAXIS
            DO NTH=1,NTHMAX+1
C               RPS(NTH,NR)=RPS(NTH,NRPMAX)
C     &                    +(RPS(NTH,NRPMAX)-RPS(NTH,NRPMAX-1))*FACTOR
C               ZPS(NTH,NR)=ZPS(NTH,NRPMAX)       *FACTOR
C     &                    +(ZPS(NTH,NRPMAX)-ZPS(NTH,NRPMAX-1))*FACTOR
C               RPS(NTH,NR)=RAXIS+(RPS(NTH,NRPMAX)-RAXIS)*FACTOR !
C               ZPS(NTH,NR)=       ZPS(NTH,NRPMAX)       *FACTOR !
C
C
c
c
               do i=1,np
                  nra(i)=nr-np-1+i
                  rpspa(i)=rps(nth,nr-np-1+i)
               enddo
c
                  call polint(nra,rpspa,np,nr,rps(nth,nr),dy)
c
                  do i=1,np
                     nra(i)=nr-np-1+i
                     zpspa(i)=zps(nth,nr-np-1+i)
                  enddo
               call polint(nra,zpspa,np,nr,zps(nth,nr),dy)
c
c
C
               RMIN=MIN(RMIN,RPS(NTH,NR))
               RMAX=MAX(RMAX,RPS(NTH,NR))
               ZMIN=MIN(ZMIN,ZPS(NTH,NR))
               ZMAX=MAX(ZMAX,ZPS(NTH,NR))
            ENDDO
            RRMIN(NR)=RMIN
            RRMAX(NR)=RMAX
            ZZMIN(NR)=ZMIN
            ZZMAX(NR)=ZMAX
            BBMIN(NR)=BBMIN(NRPMAX)/FACTOR
            BBMAX(NR)=BBMAX(NRPMAX)*FACTOR
         ENDDO
c
c
cc         do nr=1,nrmax
cc            do nth=1,nthmax
cc               print*, nr,nth,rps(nth,nr)
cc           enddo
cc         enddo
C
         DO NR=NRPMAX+1,NRMAX
            RST(NR)=SQRT(PSIT(NR)/(PI*BB))
         ENDDO
      ELSE
C
C     ****** Free boundary without X point ******
C
         DR=(RB-RA)/(NRMAX-NRPMAX)
         DTH=2*PI/NTHMAX
C
C     ----- CALCULATE PSI,PPS,TTS,PSIT and RPS, ZPS on mag surfaces -----
C
         DO NR=NRPMAX+1,NRMAX
            PSIP(NR)=PSIP(NRPMAX)
            PPS(NR)=0.D0
            TTS(NR)=2.D0*PI*BB*RR
            SPS(NR)=SPS(NRPMAX)
            VPS(NR)=VPS(NRPMAX)
            QPS(NR)=QPS(NRPMAX)
            RLEN(NR)=RLEN(NRPMAX)
            RRMIN(NR)=RRMIN(NRPMAX)
            RRMAX(NR)=RRMAX(NRPMAX)
            ZZMIN(NR)=ZZMIN(NRPMAX)
            ZZMAX(NR)=ZZMAX(NRPMAX)
            BBMIN(NR)=BBMIN(NRPMAX)
            BBMAX(NR)=BBMAX(NRPMAX)
            AVERHR(NR)=AVERHR(NRPMAX)
            AVEIR2(NR)=AVEIR2(NRPMAX)
            AVERH1(NR)=AVERH1(NRPMAX)
            AVERH2(NR)=AVERH2(NRPMAX)
            AVEBB2(NR)=AVEBB2(NRPMAX)
            AVEIB2(NR)=AVEIB2(NRPMAX)
            AVERHB(NR)=AVERHB(NRPMAX)
            DO NTH=1,NTHMAX+1
               RPS(NTH,NR)=RPS(NTH,NRPMAX)
               ZPS(NTH,NR)=ZPS(NTH,NRPMAX)
            ENDDO
         ENDDO
C
         DO NR=NRPMAX+1,NRMAX
            RINIT=REDGE+DR*(NR-NRPMAX)
            ZINIT=ZAXIS
            PSIP(NR)=PSIG(RINIT,ZINIT)-PSI0
            PPS(NR)=0.D0
            TTS(NR)=2.D0*PI*BB*RR
C
C         WRITE(6,'(A,I5,1P3E12.4)') 'NR:',NR,
C     &        PSIP(NR),PPS(NR),TTS(NR)
C
            CALL EQMAGS(RINIT,ZINIT,NTVMAX,XA,YA,NA,IERR)
            IF(IERR.NE.0) GOTO 1000
C
            SUMS=0.D0
            SUMV=0.D0
            SUMQ=0.D0
            SUMAVRHR=0.D0
            SUMAVIR2=0.D0
            SUMAVRH1=0.D0
            SUMAVRH2=0.D0
            SUMAVBB2=0.D0
            SUMAVIB2=0.D0
            SUMAVRHB=0.D0
C
            RMIN=RAXIS
            RMAX=RAXIS
            ZMIN=ZAXIS
            ZMAX=ZAXIS
            BMIN=ABS(2.D0*BB)
            BMAX=0.D0
C
            XCHI0(1)=0.D0
            XCHI1(1)=0.D0
            RCHI(1)=RINIT
            ZCHI(1)=ZINIT
C
            DO N=2,NA
               H=XA(N)-XA(N-1)
               R=0.5D0*(YA(1,N-1)+YA(1,N))
               Z=0.5D0*(YA(2,N-1)+YA(2,N))
               CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
               BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
               BTL=TTS(NR)/(2.D0*PI*R)
               B2L=BTL**2+BPL**2
C
               SUMS=SUMS+H/(BPL*R)
               SUMV=SUMV+H/BPL
               SUMQ=SUMQ+H/(BPL*R*R)
C
               SUMAVRHR=SUMAVRHR+H*BPL
               SUMAVIR2=SUMAVIR2+H/(BPL*R*R)
               SUMAVRH1=SUMAVRH1+H*R
               SUMAVRH2=SUMAVRH2+H*R*R*BPL
               SUMAVBB2=SUMAVBB2+H*B2L/BPL
               SUMAVIB2=SUMAVIB2+H/(B2L*BPL)
               SUMAVRHB=SUMAVRHB+H*R*R*BPL/B2L
C
               XCHI1(N)=SUMAVIR2
               RCHI(N)=YA(1,N)
               ZCHI(N)=YA(2,N)
C
               R=YA(1,N)
               Z=YA(2,N)
               CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
               BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
               BTL=TTS(NR)/(2.D0*PI*R)
               B2L=BTL**2+BPL**2
               B=SQRT(B2L)
C
               RMIN=MIN(RMIN,R)
               RMAX=MAX(RMAX,R)
               ZMIN=MIN(ZMIN,Z)
               ZMAX=MAX(ZMAX,Z)
               BMIN=MIN(BMIN,B)
               BMAX=MAX(BMAX,B)
            ENDDO
C
            SPS(NR)=SUMS/(2.D0*PI)
            VPS(NR)=SUMV
            QPS(NR)=SUMAVIR2*TTS(NR)/(4.D0*PI**2)
            RLEN(NR)=XA(NA)
            RRMIN(NR)=RMIN
            RRMAX(NR)=RMAX
            ZZMIN(NR)=ZMIN
            ZZMAX(NR)=ZMAX
            BBMIN(NR)=BMIN
            BBMAX(NR)=BMAX
            AVERHR(NR)=SUMAVRHR/SUMV
            AVEIR2(NR)=SUMAVIR2/SUMV
            AVERH1(NR)=SUMAVRH1/SUMV
            AVERH2(NR)=SUMAVRH2/SUMV
            AVEBB2(NR)=SUMAVBB2/SUMV
            AVEIB2(NR)=SUMAVIB2/SUMV
            AVERHB(NR)=SUMAVRHB/SUMV
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
                  IF(IERR.NE.0) 
     &                 WRITE(6,*) 'XX SPL1D for RCHI: IERR=',IERR
                  CALL SPL1D(XCHI0,ZCHI,DXCHI,UZCHI,NA,4,IERR)
                  IF(IERR.NE.0) 
     &                 WRITE(6,*) 'XX SPL1D for ZCHI: IERR=',IERR
               ELSE
                  CALL SPL1D(XCHI1,RCHI,DXCHI,URCHI,NA,4,IERR)
                  IF(IERR.NE.0) 
     &                 WRITE(6,*) 'XX SPL1D for RCHI: IERR=',IERR
                  CALL SPL1D(XCHI1,ZCHI,DXCHI,UZCHI,NA,4,IERR)
                  IF(IERR.NE.0) 
     &                 WRITE(6,*) 'XX SPL1D for ZCHI: IERR=',IERR
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
 1000    CONTINUE
C
C     ----- CALCULATE TOROIDAL FLUX -----
C
         DO NR=NRPMAX+1,NRMAX
            PSIT(NR)=PSIT(NR-1)
     &              +2.0D0*QPS(NR)*QPS(NR-1)/(QPS(NR)+QPS(NR-1))
     &                    *(PSIP(NR)-PSIP(NR-1))
         ENDDO
C
         DO NR=NRPMAX+1,NRMAX
            RST(NR)=SQRT(PSIT(NR)/(PI*BB))
         ENDDO
C
         DO NR=NRPMAX+1,NRMAX
            AVERHR(NR)=AVERHR(NR)*QPS(NR)**2/(RST(NR)**2*BB**2)/RSTA**2
            AVEIR2(NR)=AVEIR2(NR)*RR**2
            AVERH1(NR)=AVERH1(NR)*QPS(NR)/(RST(NR)*BB)/RSTA
            AVERH2(NR)=AVERH2(NR)*QPS(NR)**2/(RST(NR)**2*BB**2)/RSTA**2
            AVEBB2(NR)=AVEBB2(NR)/BB**2
            AVEIB2(NR)=AVEIB2(NR)*BB**2
            AVERHB(NR)=AVERHB(NR)*QPS(NR)**2/(RST(NR)**2*RR**2)/RSTA**2
            AVEJPR(NR)=0.D0
            AVEJTR(NR)=0.D0
         ENDDO
      ENDIF
C
      PSITB=PSIT(NRMAX)
      PSIPB=PSIP(NRMAX)
C
      DO NR=NRPMAX+1,NRMAX
         RST(NR)=SQRT(PSIT(NR)/(PI*BB))
      ENDDO
      RSTB=RST(NRMAX)
C
      DO NR=NRPMAX+1,NRMAX
         RHOT(NR)=SQRT(PSIT(NR)/PSITA)
      ENDDO
C
C      ----- CALCULATE PLASMA SURFACE -----
C
      CALL EQCALF(REDGE,ZAXIS,NSUMAX,RSU,ZSU,IERR)
C
C      +++++ CALCULATE WALL DATA +++++
C
C     FACTOR=(RB+RR-RAXIS)/(REDGE-RAXIS)
C
      DO NTH=1,NTHMAX
         THW(NTH)=(NTH-1)*2*PI/NTHMAX
         RPSW(NTH)=RPS(NTH,NRMAX)
         ZPSW(NTH)=ZPS(NTH,NRMAX)
C         WRITE(6,'(A,I5,1P3E12.4)') 
C     &        'NTH: ',NTH,THW(NTH),RPSW(NTH),ZPSW(NTH)
      ENDDO
      NTH=NTHMAX+1
      THW(NTH)=2*PI
      RPSW(NTH)=RPS(1,NRMAX)
      ZPSW(NTH)=ZPS(1,NRMAX)
C         WRITE(6,'(A,I5,1P3E12.4)') 
C     &        'NTH: ',NTH,THW(NTH),RPSW(NTH),ZPSW(NTH)
      CALL SPL1D(THW,RPSW,DRPSW,URPSW,NTHMAX+1,4,IERR)
      CALL SPL1D(THW,ZPSW,DZPSW,UZPSW,NTHMAX+1,4,IERR)
      DO NSU=1,NSUMAX+1
         THWL=(NSU-1)*2*PI/NSUMAX
         CALL SPL1DF(THWL,RSW(NSU),THW,URPSW,NTHMAX+1,IERR)
         CALL SPL1DF(THWL,ZSW(NSU),THW,UZPSW,NTHMAX+1,IERR)
C         WRITE(6,'(A,I5,1P3E12.4)') 
C     &        'NSU: ',NSU,THWL,RSW(NSU),ZSW(NSU)
      ENDDO
C
      IF(MDLEQF.LT.10) THEN
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
      ENDIF
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
      DIMENSION DERIV(NRM)
      DIMENSION D01(NTHMP,NRM),D10(NTHMP,NRM),D11(NTHMP,NRM)
C
C      DO NR=1,NRMAX
C         WRITE(6,'(A,I5,1P5E12.4)')
C     &        'NR:',NR,PSIP(NR),PSIT(NR),QPS(NR),TTS(NR),AVEJPR(NR)
C      ENDDO
C
      IERR=0
      DTH=2*PI/NTHMAX
C
      CALL SPL1D(PSIP,PSIT,DERIV,UPSIT,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSIT: IERR=',IERR
      CALL SPL1D(PSIT,PSIP,DERIV,UPSIP,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSIP: IERR=',IERR
C
      DO NR=1,NRMAX
         RHOT(NR)=SQRT(PSIT(NR)/PSITA)
      ENDDO
C
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,PSIP,DERIV,UPSIPS,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSIP: IERR=',IERR
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
      CALL SPL1D(RHOT,VPS,DERIV,UVPS,NRMAX,0,IERR)
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
      CALL SPL1D(RHOT,ZZMIN,DERIV,UZZMIN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for ZZMIN: IERR=',IERR
      CALL SPL1D(RHOT,ZZMAX,DERIV,UZZMAX,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for ZZMAX: IERR=',IERR
      CALL SPL1D(RHOT,BBMIN,DERIV,UBBMIN,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for BBMIN: IERR=',IERR
      CALL SPL1D(RHOT,BBMAX,DERIV,UBBMAX,NRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for BBMAX: IERR=',IERR
C
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVERHR,DERIV,UAVERHR,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVERHR: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVEIR2,DERIV,UAVEIR2,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVEIR2: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVERH1,DERIV,UAVERH1,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVERH1: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVERH2,DERIV,UAVERH2,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVERH2: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVEBB2,DERIV,UAVEBB2,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVEBB2: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVEIB2,DERIV,UAVEIB2,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVEIB2: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVERHB,DERIV,UAVERHB,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVERHB: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVEJPR,DERIV,UAVEJPR,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVEJPR: IERR=',IERR
      DERIV(1)=0.D0
      CALL SPL1D(RHOT,AVEJTR,DERIV,UAVEJTR,NRMAX,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for AVEJTR: IERR=',IERR
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
      DO NR=1,NRMAX
         RHOTL=RHOT(NR)
         CALL SPL1DF(RHOTL,QPSL,RHOT,UQPS,NRMAX,IERR)
         DO NTH=1,NTHMAX+1
            THITL=THIT(NTH)
            CALL SPL2DD(THITL,RHOTL,RPSL,DRCHIL,DRRHOL,
     &                  THIT,RHOT,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
            CALL SPL2DD(THITL,RHOTL,ZPSL,DZCHIL,DZRHOL,
     &                  THIT,RHOT,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
            DRPSI(NTH,NR)=DRRHOL*QPSL/(2.D0*PSITA)
            DZPSI(NTH,NR)=DZRHOL*QPSL/(2.D0*PSITA)
            DRCHI(NTH,NR)=DRCHIL
            DZCHI(NTH,NR)=DZCHIL
C
C            IF(NTH.EQ.1) WRITE(6,'(I5,1P4E12.4)')
C     &           NR,DRRHOL,DRPSI(NTH,NR),QPSL,PSITA
C
            CALL EQPSID(RPSL,ZPSL,DPSIDR,DPSIDZ)
            BPR(NTH,NR)= DPSIDZ/(2.D0*PI*RPSL)
            BPZ(NTH,NR)=-DPSIDR/(2.D0*PI*RPSL)
            BPT(NTH,NR)=SQRT(BPR(NTH,NR)**2+BPZ(NTH,NR)**2)
            BTP(NTH,NR)= TTS(NR)/(2.D0*PI*RPSL)
C
C            IF(NTH.EQ.1) WRITE(6,'(2I3,1P6E12.4)') 
C     &           NTH,NR,RPSL,ZPSL,
C     &           BPR(NTH,NR),BPZ(NTH,NR),BPT(NTH,NR),BTP(NTH,NR)
C
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
      DELPS=-PSI0/NDPMAX
      PPSL=FNPPS(0.D0)
      VPSL=FNVPS(0.D0)
      SPSL=FNSPS(0.D0)
      SUMV =0.5D0*VPSL*DELPS
      SUMS =0.5D0*SPSL*DELPS
      SUMPV=0.5D0*PPSL*VPSL*DELPS
      SUMPS=0.5D0*PPSL*SPSL*DELPS
      DO NDP=1,NDPMAX-1
         PSIL=PSI0+DELPS*NDP
         PSIN=1.D0-PSIL/PSI0
         RHON=FNRHON(PSIN)
         PPSL=FNPPS(RHON)
         VPSL=FNVPS(RHON)
         SPSL=FNSPS(RHON)
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
      DIMENSION XA(NTVM),YA(2,NTVM)
      DIMENSION RCHI(NTVM),ZCHI(NTVM),DXCHI(NTVM)
      DIMENSION URCHI(4,NTVM),UZCHI(4,NTVM)
C
      IERR=0
      DTH=2*PI/NTHUMAX
C
C     ----- SET NUMBER OF DIVISION for integration -----
C

      NMAX=200
      IF(NMAX.GT.NTVM) NMAX=NTVM
C
C     ----- CALCULATE PSIP, PSIT, PPS, TTS, RPS and ZPS -----
C     -----              on magnetic surfaces           -----
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
C     ***** INTERPOLATE FUNCTION OF PP(PSIP) *****
C
      FUNCTION PPFUNC(PSIPL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(PSIPL,PPL,PSIPS,UPPPS,NPSMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX PPFUNC: SPL1DF ERROR : IERR=',IERR
         WRITE(6,*) PSIPL,PSIPS(1),PSIPS(NPSMAX)
      ENDIF
      PPFUNC=PPL
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF TT(PSIP) *****
C
      FUNCTION TTFUNC(PSIPL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DF(PSIPL,TTL,PSIPS,UTTPS,NPSMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TTFUNC: SPL1DF ERROR : IERR=',IERR
      TTFUNC=TTL
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF PP(PSIP) *****
C
      FUNCTION DPPFUNC(PSIPL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DD(PSIPL,PPL,DPPL,PSIPS,UPPPS,NPSMAX,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX PPFUNC: SPL1DF ERROR : IERR=',IERR
         WRITE(6,*) PSIPL,PSIPS(1),PSIPS(NPSMAX)
      ENDIF
      DPPFUNC=DPPL
      RETURN
      END
C
C     ***** INTERPOLATE FUNCTION OF TT(PSIP) *****
C
      FUNCTION DTTFUNC(PSIPL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL1DD(PSIPL,TTL,DTTL,PSIPS,UTTPS,NPSMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX TTFUNC: SPL1DF ERROR : IERR=',IERR
      DTTFUNC=DTTL
      RETURN
      END
C
C
C **********************************************
      SUBROUTINE polint(nra,psa,n,nr,ps,dy)
      implicitnone
      INTEGER n,NMAX,nr
      REAL*8 dy,nra(n),ps,psa(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
c
c      do i=1,n
c         print*, i,nra(i),psa(i)
c      enddo
c
      ns=1
      dif=abs(nr-nra(1))
      do 11 i=1,n
        dift=abs(nr-nra(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=psa(i)
        d(i)=psa(i)
11    continue
      ps=psa(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=nra(i)-nr
          hp=nra(i+m)-nr
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
        endif
        ps=ps+dy
13    continue
      return
      END
c
c
c ***************************************
c      SUBROUTINE ratint(xa,ya,n,x,y,dy)
c      INTEGER n,NMAX
c      REAL*8 dy,x,y,xa(n),ya(n),TINY
c      PARAMETER (NMAX=10,TINY=1.e-25)
c      INTEGER i,m,ns
c      REAL dd,h,hh,t,w,c(NMAX),d(NMAX)
c
c      do i=1,n
c         print*,i,xa(i),ya(i)
c      enddo
c
c      ns=1
c      hh=abs(x-xa(1))
c      do 11 i=1,n
c        h=abs(x-xa(i))
c        if (h.eq.0.)then
c          y=ya(i)
c          dy=0.0
c          return
c        else if (h.lt.hh) then
c          ns=i
c          hh=h
c        endif
c        c(i)=ya(i)
c        d(i)=ya(i)+TINY
c11    continue
c      y=ya(ns)
c      ns=ns-1
c      do 13 m=1,n-1
c        do 12 i=1,n-m
c          w=c(i+1)-d(i)
c          h=xa(i+m)-x
c          t=(xa(i)-x)*d(i)/h
c          dd=t-c(i+1)
c          if(dd.eq.0.)pause 'failure in ratint'
c          dd=w/dd
c          d(i)=c(i+1)*dd
c          c(i)=t*dd
c12      continue
c        if (2*ns.lt.n-m)then
c          dy=c(ns+1)
c        else
c          dy=d(ns)
c          ns=ns-1
c        endif
c        y=y+dy
c13    continue
c      return
c      END
