C     $Id$
C
C     ***** LOAD JFT2M EQ DATA *****
C
      SUBROUTINE EQRJFT(KNAMEQ)
C
      INCLUDE 'eqrcom.inc'
C
      DIMENSION PSI(NRGM,NZGM),RCU(NRGM,NZGM)
      DIMENSION PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM)
      DIMENSION READP(NGRZM)
      CHARACTER*72 KNAMEQ,KNAM0D,KNAM1D
      CHARACTER*80 DREAD
C
C     ----- R and Z grid and PSI data read -----
C
      DO IK=72,1,-1
         IKMAX=IK
         IF(KNAMEQ(IK:IK).NE.' ') GOTO 100
      ENDDO
  100 KNAM0D=KNAMEQ(1:IKMAX)//'.0d'
      KNAM1D=KNAMEQ(1:IKMAX)//'.1d'
C
      OPEN(7,FILE=KNAM1D,FORM='FORMATTED',STATUS='OLD')
C
      WRITE(6,*) 'Now reading psi data from :',KNAM1D
C
      READ(7,700) NSHOT,TIME
      READ(7,710) NRGMAX
      READ(7,710) NZGMAX
      READ(7,710) NSUMAX
      WRITE(6,600) NSHOT,TIME
      WRITE(6,611) NRGMAX
      WRITE(6,612) NZGMAX
      WRITE(6,613) NSUMAX
  700 FORMAT(9X,I8,13X,E11.3)
  710 FORMAT(9X,I8)
  600 FORMAT(1H ,'NSHOT  =',I8,'   TIME(MS)=',1PE11.3)
  611 FORMAT(1H ,'NRGMAX =',I8)
  612 FORMAT(1H ,'NZGMAX =',I8)
  613 FORMAT(1H ,'NSUMAX =',I8) 
C
C     ----- R grid data read -----
C
      READ(7,715) DREAD
      READ(7,715) DREAD
  715 FORMAT(A80)
C
      CALL FMREAD(7,7,NRGMAX,RG)
C      WRITE(6,615) (RG(NRG),NRG=1,NRGMAX)
  615 FORMAT(7F10.5)
C
C     ----- Z grid data read -----
C
      READ(7,715) DREAD
      READ(7,715) DREAD
C
      CALL FMREAD(7,7,NZGMAX,ZG)
C      WRITE(6,615) (ZG(NZG),NZG=1,NZGMAX)
C
C     ----- PSI data read -----
C
      READ(7,715) DREAD
      READ(7,715) DREAD
C
      NPSI=NRGMAX*NZGMAX
      CALL FMREAD(7,7,NPSI,READP)
      DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            J=NRGMAX*(NZG-1)+NRG
            PSIRZ(NRG,NZG)=READP(J)
         ENDDO
C         WRITE(6,615) (PSIRZ(NRG,NZG),NRG=1,NRGMAX)
      ENDDO
C
C      CALL EQGPRZ(RG,ZG,PSIRZ,NRGM,NRGMAX,NZGMAX)
C
      READ(7,715) DREAD
      READ(7,715) DREAD
C
      NPSI=NRGMAX*NZGMAX
      CALL FMREAD(7,7,NPSI,READP)
C
      READ(7,715) DREAD
      READ(7,715) DREAD
C
      NPSI=NRGMAX*NZGMAX
      CALL FMREAD(7,7,NPSI,READP)
C
C     ----- Contour data read -----
C
      READ(7,715) DREAD
      READ(7,715) DREAD
C
      CALL FMREAD(7,7,NSUMAX,RSU)
C
      READ(7,715) DREAD
      READ(7,715) DREAD
C
      CALL FMREAD(7,7,NSUMAX,ZSU)
C
C      WRITE(6,615) (RSU(NSU),NSU=1,NSUMAX)
C      WRITE(6,615) (ZSU(NSU),NSU=1,NSUMAX)
C
      READ(7,716) RSEPU,ZSEPU
      READ(7,716) RSEPL,ZSEPL
  716 FORMAT(11X,F10.0,14X,F10.0)
C      WRITE(6,*) RSEPU,ZSEPU,RSEPL,ZSEPL
C
      CLOSE(7)
C
      OPEN(8,FILE=KNAM0D,FORM='FORMATTED',STATUS='OLD')
      NDATA=39
C
      WRITE(6,*) 'Now reading psi0d data from :',KNAM0D
C
      DO I=1,6
         READ(8,715) DREAD
      ENDDO
C
      CALL FMREAD(6,8,NDATA,READP)
C      WRITE(6,620) (READP(J),J=1,NDATA)
C
      CLOSE(8)
C
      RR    = READP( 3)
      RAXIS = READP( 7)
      ZAXIS = READP( 8)
      RKAP  = READP( 9)
      RDLT  = READP(10)
      RA    = READP(17)
      RIP   = READP(20)
      BB    = READP(21)
      QA    = READP(11)
C
      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,U,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
C
      CALL SPL2DF(RAXIS,ZAXIS,SAXIS,RG,ZG,U,NRGM,NRGMAX,NZGMAX,IERR)
C
      WRITE(6,*) 'RAXIS,ZAXIS,SAXIS=',RAXIS,ZAXIS,SAXIS
C
      RRPL=RR
      RAPL=RA
      RKAPPL=RKAP
      RDLTPL=RDLT
      BBPL=BB
      RIPPL=RIP
      Q0PL=QPS(NPSMAX)
      QAPL=QSU
C
      RETURN
  620 FORMAT(7(1X,1PE10.3))
      END
C
C     ***** READ FORMATED FILE *****
C
      SUBROUTINE FMREAD(NI,IFN,NDATA,DATA)
C
      IMPLICIT REAL*8(A-F,H,O-Z)
      DIMENSION DATA(NDATA)
C
      IDATA=(NDATA-1)/NI+1
      DO I=1,IDATA
         NST=1+NI*(I-1)
         NED=  NI*I
         IF(I.EQ.IDATA) NED=NDATA
         READ(IFN,*,END=9000,ERR=8000) (DATA(J),J=NST,NED)
      ENDDO
      RETURN
C
 8000 WRITE(6,*) 'XX FMREAD: FILE READ ERROR'
      RETURN
 9000 WRITE(6,*) 'XX FMREAD: END OF FILE'
      RETURN
      END
