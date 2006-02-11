C
C********************************************
C*         graphic  OUTPUT             *
C********************************************
C
      SUBROUTINE EQGX2D
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION RSPL(4*NRGM),ZSPL(4*NZGM)
      DIMENSION GR(4*NRGM),GZ(4*NZGM),GF(4*NRGM,4*NZGM)
      DIMENSION KA(8,4*NRGM,4*NZGM)
C
      GX1=2.0
      GX2=17.0
      GY1=2.0
      GY2=17.0
C
      DRSPL=(RGMAX-RGMIN)/(4*NRGMAX-1)
      DO NR=1,4*NRGMAX
         RSPL(NR)=RGMIN+DRSPL*(NR-1)
      ENDDO
      DZSPL=(ZGMAX-ZGMIN)/(4*NZGMAX-1)
      DO NZ=1,4*NZGMAX
         ZSPL(NZ)=ZGMIN+DZSPL*(NZ-1)
      ENDDO
C
      DO NR=1,4*NRGMAX
         GR(NR)=GUCLIP(RSPL(NR))
      ENDDO
      DO NZ=1,4*NZGMAX
         GZ(NZ)=GUCLIP(ZSPL(NZ))
      ENDDO
C
      DO NZ=1,4*NZGMAX
      DO NR=1,4*NRGMAX
         GF(NR,NZ)=PSIXF(RSPL(NR),ZSPL(NZ))
      ENDDO
      ENDDO
C
      CALL PAGES
      CALL GSUB2D(GX1,GX2,GY1,GY2,GR,GZ,GF,
     &            4*NRGM,4*NRGMAX,4*NZGMAX,KA,'/PSIXF(R,Z)/')
      CALL PAGEE
C
      DO NR=1,NRGMAX
         GR(NR)=GUCLIP(RG(NR))
      ENDDO
      DO NZ=1,NZGMAX
         GZ(NZ)=GUCLIP(ZG(NZ))
      ENDDO
C
      DO NZ=1,NZGMAX
      DO NR=1,NRGMAX
         GF(NR,NZ)=GUCLIP(HJTRZ(NR,NZ))
      ENDDO
      ENDDO
C
      CALL PAGES
      CALL GSUB2D(GX1,GX2,GY1,GY2,GR,GZ,GF,
     &            4*NRGM,NRGMAX,NZGMAX,KA,'/FJTRZ(R,Z)/')
      CALL PAGEE
C
      RETURN
      END
C
C     ****** CONTOUR PLOT ******
C
      SUBROUTINE GSUB2D(GX1,GX2,GY1,GY2,GX,GY,GF,
     &                  NXM,NXMAX,NYMAX,KA,KTITLE)
C
      DIMENSION GX(NXM),GY(NYMAX),GF(NXM,NYMAX)
      DIMENSION KA(8,NXM,NYMAX)
      CHARACTER KTITLE*(*)
C
      CALL GMNMX2(GF,NXM,1,NXMAX,1,1,NYMAX,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GFPMIN,GFPMAX,GFSTEP)
      CALL GQSCAL(GX(1),GX(NXMAX),GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(GY(1),GY(NYMAX),GYMIN,GYMAX,GYSTEP)
      IF(GX(1)*GX(NXMAX).LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GX(1)
      ENDIF
      IF(GY(1)*GY(NYMAX).LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GX(1)
      ENDIF
      GFSTEP=0.5*GFSTEP
C
      GXP=GX2-GX1
      GYP=GY2-GY1
      GXLEN=GX(NXMAX)-GX(1)
      GYLEN=GY(NYMAX)-GY(1)
      IF(GXLEN*GYP.GT.GYLEN*GXP) THEN
         GYP=GXP*GYLEN/GXLEN
      ELSE
         GXP=GYP*GXLEN/GYLEN
      ENDIF
C
      CALL SETLIN(0,2,7)
      CALL SETCHS(0.3,0.0)
      CALL MOVE(GX1,GY2+0.2)
      CALL TEXTX(KTITLE)
C
      CALL GDEFIN(GX1,GX1+GXP,GY1,GY1+GYP,
     &            GX(1),GX(NXMAX),GY(1),GY(NYMAX))
      CALL GFRAME
      CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.2,9)
      CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
C
      CALL SETLIN(0,-1,7)
      IF(GFMIN*GFMAX.LE.0.0) THEN
      CALL CONTQ2(GF,GX,GY,NXM,NXMAX,NYMAX,0.0,GFMAX-GFMIN,1,
     &            0,4,KA)
      CALL CONTQ2(GF,GX,GY,NXM,NXMAX,NYMAX, GFSTEP, GFSTEP,20,
     &            0,0,KA)
      CALL CONTQ2(GF,GX,GY,NXM,NXMAX,NYMAX,-GFSTEP,-GFSTEP,20,
     &            0,2,KA)
      ELSEIF(GFMIN.GT.0.0) THEN
      CALL CONTQ2(GF,GX,GY,NXM,NXMAX,NYMAX, GFPMIN, GFSTEP,20,
     &            0,0,KA)
      ELSE
      CALL CONTQ2(GF,GX,GY,NXM,NXMAX,NYMAX, GFPMIN, GFSTEP,20,
     &            0,2,KA)
      ENDIF
C
      CALL EQGPRX(GX2,GY2,GFMIN,GFMAX,GFSTEP)
C
      RETURN
      END
C
C     ****** DRAW PARM ******
C
      SUBROUTINE EQGPRX(GX2,GY2,GFMIN,GFMAX,GFSTEP)
C
      INCLUDE '../eq/eqcomx.inc'
C
      CALL SETLIN(-1,-1,7)
      CALL MOVE(GX2+0.5,GY2-0.3)
      CALL TEXTX('/MAX =/')
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      CALL MOVE(GX2+0.5,GY2-0.8)
      CALL TEXTX('/MIN =/')
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(GX2+0.5,GY2-1.3)
      CALL TEXTX('/STEP=/')
      CALL NUMBR(GFSTEP,'(1PE12.4)',12)
C
      CALL MOVE(17.5,15.0)
      CALL TEXT('PSIB(0)=',8)
      CALL NUMBD(PSIB(0),  '(1PE11.3)',11)
      CALL MOVE(17.5,14.5)
      CALL TEXT('PSIN(1)=',8)
      CALL NUMBD(PSIB(1),  '(1PE11.3)',11)
      CALL MOVE(17.5,14.0)
      CALL TEXT('PSIB(2)=',8)
      CALL NUMBD(PSIB(2),  '(1PE11.3)',11)
      CALL MOVE(17.5,13.5)
      CALL TEXT('PSIB(3)=',8)
      CALL NUMBD(PSIB(3),  '(1PE11.3)',11)
      CALL MOVE(17.5,13.0)
      CALL TEXT('PSIB(4)=',8)
      CALL NUMBD(PSIB(4),  '(1PE11.3)',11)
      CALL MOVE(17.5,12.5)
      CALL TEXT('PSIB(5)=',8)
      CALL NUMBD(PSIB(5),  '(1PE11.3)',11)
C
      CALL MOVE(17.5,12.0)
      CALL TEXT('NRGMAX = ',9)
      CALL NUMBI(NRGMAX,    '(I2)',2)
      CALL MOVE(17.5,11.5)
      CALL TEXT('NZGMAX = ',9)
      CALL NUMBI(NZGMAX,    '(I2)',2)
      CALL MOVE(17.5,11.0)
      CALL TEXT('RR     =',8)
      CALL NUMBD(RR,       '(1PE11.3)',11)
      CALL MOVE(17.5,10.5)
      CALL TEXT('RA     =',8)
      CALL NUMBD(RA,       '(1PE11.3)',11)
      CALL MOVE(17.5,10.0)
      CALL TEXT('RKAP   =',8)
      CALL NUMBD(RKAP,     '(1PE11.3)',11)
      CALL MOVE(17.5,9.5)
      CALL TEXT('RDLT   =',8)
      CALL NUMBD(RDLT,     '(1PE11.3)',11)
      CALL MOVE(17.5,9.0)
      CALL TEXT('RIP    =',8)
      CALL NUMBD(RIP,      '(1PE11.3)',11)
      CALL MOVE(17.5,8.5)
      CALL TEXT('BB     =',8)
      CALL NUMBD(BB,       '(1PE11.3)',11)
      CALL MOVE(17.5,8.0)
      CALL TEXT('P0     =',8)
      CALL NUMBD(P0,       '(1PE11.3)',11)
      CALL MOVE(17.5,7.5)
      CALL TEXT('P1     =',8)
      CALL NUMBD(P1,       '(1PE11.3)',11)
      CALL MOVE(17.5,7.0)
      CALL TEXT('PROFP0 =',8)
      CALL NUMBD(PROFP0,   '(1PE11.3)',11)
      CALL MOVE(17.5,6.5)
      CALL TEXT('PROFP1 =',8)
      CALL NUMBD(PROFP1,   '(1PE11.3)',11)
      CALL MOVE(17.5,6.0)
      CALL TEXT('PROFT  =',8)
      CALL NUMBD(PROFT,    '(1PE11.3)',11)
      CALL MOVE(17.5,5.5)
      CALL TEXT('RMIN   =',8)
      CALL NUMBD(RMIN,     '(1PE11.3)',11)
      CALL MOVE(17.5,5.0)
      CALL TEXT('RMAX   =',8)
      CALL NUMBD(RMAX,     '(1PE11.3)',11)
      CALL MOVE(17.5,4.5)
      CALL TEXT('ZMIN   =',8)
      CALL NUMBD(ZMIN,     '(1PE11.3)',11)
      CALL MOVE(17.5,4.0)
      CALL TEXT('ZMAX   =',8)
      CALL NUMBD(ZMAX,     '(1PE11.3)',11)
      CALL MOVE(17.5,3.5)
      CALL TEXT('EPS    =',8)
      CALL NUMBD(EPS,      '(1PE11.3)',11)
C
      RETURN
      END
C
C     ****** CORE SOLVER ******
C
      SUBROUTINE EQGX1D
C
      INCLUDE '../eq/eqcomx.inc'
      PARAMETER (NIM=301)
      DIMENSION GX(NIM),GY(NIM,3)
C
      NIMAX=301
      DRG=(RGMAX-RGMIN)/(NIMAX-1)
      DO I=1,NIMAX
         R=RGMIN+DRG*(I-1)
         Z=0.D0
         PSIL=PSIXF(R,Z)
         CALL EQPSIX(R,Z,DPSIDR,DPSIDZ,IERR)
         GX(I)=GUCLIP(R)
         GY(I,1)=GUCLIP(PSIL)
         GY(I,2)=GUCLIP(DPSIDR)
         GY(I,3)=GUCLIP(DPSIDZ)
      ENDDO
C
      CALL PAGES
      CALL GRF1D(1,GX,GY(1,1),NIM,NIMAX,1,'@PSI vs R@',0)
      CALL GRF1D(2,GX,GY(1,2),NIM,NIMAX,1,'@DPSIDR vs R@',0)
      CALL GRF1D(3,GX,GY(1,3),NIM,NIMAX,1,'@DPSIDZ vs R@',0)
      CALL PAGEE
C
      DZG=(ZGMAX-ZGMIN)/(NIMAX-1)
      DO I=1,NIMAX
         R=RR
         Z=ZGMIN+DZG*(I-1)
         PSIL=PSIXF(R,Z)
         CALL EQPSIX(R,Z,DPSIDR,DPSIDZ,IERR)
         GX(I)=GUCLIP(Z)
         GY(I,1)=GUCLIP(PSIL)
         GY(I,2)=GUCLIP(DPSIDR)
         GY(I,3)=GUCLIP(DPSIDZ)
      ENDDO
C     
      CALL PAGES
      CALL GRF1D(1,GX,GY(1,1),NIM,NIMAX,1,'@PSI vs Z@',0)
      CALL GRF1D(2,GX,GY(1,2),NIM,NIMAX,1,'@DPSIDR vs Z@',0)
      CALL GRF1D(3,GX,GY(1,3),NIM,NIMAX,1,'@DPSIDZ vs Z@',0)
      CALL PAGEE
C
      RETURN
      END
