C
C     ****** TASK/W1 ******
C
C     ******* DIVISION IN X DIRECTION *******
C
      SUBROUTINE W1SETX(IERR)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER:: NXVH,NXP1,NSACM,NS,NX,NX1,NX2,I,IERR
      REAL(rkind):: RW,FWC,WCH,WCL,XACM,DX,DXACM,FVT,WT,X,FACT
      REAL(rkind):: X1,DX1,DX2
C
      NXT=NXP+2*NXV
         IF(NXP.GT.NXPM) GO TO 9100
         IF(NXV.GT.NXVM) GO TO 9100
      NXVH=NXV/2
C
      IF(MOD(NMODEL,2).EQ.1) THEN
         NXP1=NXP
      ELSE
         NXP1=NXP-1
      ENDIF
C
      RW=2.D6*PI*RF
      NSACM=0
      DO 10  NS=1,NSMAX
         FWC=AEE*PZ(NS)*BB/(AMP*PA(NS))
         WCH=FWC/(1.D0+RA/RR)
         WCL=FWC/(1.D0-RA/RR)
         IF(ABS(WCH).LT.RW.AND.ABS(WCL).GT.RW) THEN
            XACM=RR*(FWC/RW-1.D0)
            NSACM=NS
         ENDIF
   10 CONTINUE
C
      IF(NSACM.EQ.0.OR.DXFACT.LT.1.D-6) THEN
         DX=2.D0*RA/DBLE(NXP1)
         DO 20 NX=1,NXP1+1
            XA(NX)=DBLE(NX-1)*DX-RA
   20    CONTINUE
      ELSE
         FVT=AEE*1.D3/(AMP*PA(NSACM))
         FWC=AEE*PZ(NSACM)*BB/(AMP*PA(NSACM))
         DXACM=DXWDTH*SQRT(FVT*PTPP(NSACM))/FWC
         WT=2.D0*RA+DXFACT*DXACM
     &             *(ATAN((RA-XACM)/DXACM)-ATAN((-RA-XACM)/DXACM))
         X=-RA
         DO 25 NX=1,NXP1
            XA(NX)=X
            FACT=1.D0+DXFACT/(1.D0+((X-XACM)/DXACM)**2)
            X1=X+0.5D0*WT/(FACT*(NXP1))
            FACT=1.D0+DXFACT/(1.D0+((X1-XACM)/DXACM)**2)
            X=X+WT/(FACT*(NXP1))
   25    CONTINUE
         WRITE(6,601) NSACM,XACM,DXACM,RA-X
         XA(NXP1+1)=RA
      ENDIF
C
      IF(MOD(NMODEL,2).EQ.1) THEN
         DO 27 NX=1,NXP
            XAM(NX)=0.5D0*(XA(NX)+XA(NX+1))
   27    CONTINUE
      ELSE
         DO 28 NX=1,NXP
            XAM(NX)=XA(NX)
   28    CONTINUE
      ENDIF
C
      DX1=(RB-RD)/DBLE(NXVH-1)
      DX2=(RD-RA)/DBLE(NXVH-1)
      DO 30 I=1,NXV
         IF(I.LE.NXVH) THEN
            X=DX1*DBLE(I-1)-RB
         ELSE
            X=DX2*DBLE(I-NXVH-1)-RD
         ENDIF
         NX1=NXP+NXV+I
         XAM(NX1)=X
         NX2=NXP+NXV+1-I
         XAM(NX2)=-X
   30 CONTINUE
      IERR=0
      RETURN
C
 9100 WRITE(6,691) NXP,NXPM,NXV,NXVM
      IERR=10000
      RETURN
  601 FORMAT(1H ,'** DX-ACCUM. : NS,X,XWIDTH,XERR = ',
     &           I3,1P2E12.4,1PE12.4)
  691 FORMAT(1H ,'!! ERROR IN W1SETX   : DIMENSION OVER'/
     &       1H ,'NXP=',I5,5X,'NXPM=',I5/
     &       1H ,'NXV=',I5,5X,'NXVM=',I5/)
      END
C
C     ******* DIVISION IN Z DIRECTION *******
C
      SUBROUTINE W1SETZ(IERR)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER:: NZH,NZ,IERR
      REAL(rkind):: AKZ0
C
         IF(NZP.GT.NZPM) GO TO 9100
      NZH=NZP/2
C
      DZ=RZ/DBLE(NZP)
      IF(NZP.EQ.1) THEN
         ZA(1)=0.D0
         AKZ(1)=RKZ
         IF(ABS(RKZ).LT.1.D-32) RKZ=0.001D0
      ELSE
         DO 40 NZ=1,NZP
            ZA(NZ)=DZ*DBLE(NZ-1) - 0.5D0*RZ
   40    CONTINUE
         AKZ0=2.D0*PI/RZ
         AKZ(1)   =0.001D0
         AKZ(NZH+1)=DBLE(NZH)*AKZ0
         DO 50 NZ=2,NZH
            AKZ(NZ)      =DBLE(NZ-1)*AKZ0
            AKZ(NZP-NZ+2)=-AKZ(NZ)
   50    CONTINUE
      ENDIF
C
      IERR=0
      RETURN
C
 9100 WRITE(6,691) NZP,NZPM
      IERR=10000
      RETURN
  691 FORMAT(1H ,'!! ERROR IN W1SETZ   : DIMENSION OVER'/
     &       1H ,'NZP=',I5,5X,'NZPM=',I5/)
      END
C
C     ******* SETING ANTENNA CURRENT DENSITY *******
C
      SUBROUTINE W1ANTS
C
      USE w1comm
      IMPLICIT NONE
      INTEGER:: NZ,NA
      REAL(rkind):: ANTPOS,PHASE
      COMPLEX(rkind):: CPHASE
C
      DO 10 NZ = 1 , NZP
         CJ1( NZ ) = ( 0.D0 , 0.D0 )
         CJ2( NZ ) = ( 0.D0 , 0.D0 )
         CJ3( NZ ) = ( 0.D0 , 0.D0 )
         CJ4( NZ ) = ( 0.D0 , 0.D0 )
   10 CONTINUE
C
      IF(NZP.EQ.1) THEN
         NZANT1(1)=1
         CJ1(1)   =AJYH(1)/DZ
         NZANT2(1)=1
         CJ2(1)   =AJYL(1)/DZ
         CJ3(1)   =AJZH(1)
         CJ4(1)   =AJZL(1)
      ELSE
      DO 100 NA = 1 , NAMAX
         ANTPOS            = RZ*ALYH(NA)/360.D0 + (RZ + DZ)*.5D0
         NZANT1(NA)        = INT( ANTPOS/DZ ) + 1
         PHASE             = APYH(NA)*PI/180.D0
         CPHASE            = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
         CJ1( NZANT1(NA) ) = AJYH( NA ) / DZ * CPHASE
     &                       + CJ1( NZANT1(NA) )
         CJ3( NZANT1(NA) ) = AJZH( NA ) * CPHASE
     &                       + CJ3( NZANT1(NA) )
C
         ANTPOS            = RZ*ALYL(NA)/360.D0 + (RZ + DZ)*.5D0
         NZANT2(NA)        = INT( ANTPOS/DZ ) + 1
         PHASE             = APYL(NA)*PI/180.D0
         CPHASE            = CMPLX( COS( PHASE ) , SIN( PHASE ) )
         CJ2( NZANT2(NA) ) = AJYL( NA ) / DZ * CPHASE
     &                       + CJ2( NZANT2(NA) )
         CJ4( NZANT2(NA) ) = AJZL( NA ) * CPHASE
     &                       + CJ4( NZANT2(NA) )
  100 CONTINUE
      ENDIF
C
      RETURN
      END
C
C     ****** INTERFACE FOR FFT ******
C
      SUBROUTINE W1FFTL(CAL,N,KEY)
C
      USE w1comm
      IMPLICIT NONE
C
      INTEGER,INTENT(IN):: N,KEY
      COMPLEX(rkind),INTENT(INOUT):: CAL(N)

      INTEGER,SAVE:: N_SAVE=0
      INTEGER:: LP,IND,I
C
      IF(N.GT.1) THEN
         LP=LOG(DBLE(N))/LOG(2.D0)+0.5D0
         IF(N.EQ.N_SAVE) THEN
            IND=0
         ELSE
            IND=1
            N_SAVE=N
         ENDIF
         DO 10 I=1,N
           CT(I)=CAL(I)
   10    CONTINUE
         CALL FFT2L(CT,CAL,CFFT,LFFT,N/2,IND,KEY+1,LP)
         IF(KEY.EQ.0) CAL(N/2+1)=0.D0
      ENDIF
      RETURN
      END
C
C     ****** SUBSTITUTION BY VIRUE OF SYMMETRY IN Z-DIRECTION ******
C
      SUBROUTINE W1SYMS(NZ,NSYML)
C
      USE w1comm
C
      NZ0 = NZP +2 -NZ
C
      IF(NSYML.EQ.1) THEN
         DO 70 NX = 1, NXT
            CE2DA(NZ,NX,1) =  CE2DA(NZ0,NX,1)
            CE2DA(NZ,NX,2) =  CE2DA(NZ0,NX,2)
            CE2DA(NZ,NX,3) = -CE2DA(NZ0,NX,3)
   70    CONTINUE
      ELSEIF(NSYML.EQ.-1) THEN
         DO 80 NX = 1, NXT
            CE2DA(NZ,NX,1) = -CE2DA(NZ0,NX,1)
            CE2DA(NZ,NX,2) = -CE2DA(NZ0,NX,2)
            CE2DA(NZ,NX,3) =  CE2DA(NZ0,NX,3)
   80    CONTINUE
      ENDIF
      PAKT(NZ,1) =  PAKT(NZ0,1)
      PAKT(NZ,2) =  PAKT(NZ0,2)
      PAKT(NZ,3) =  PAKT(NZ0,3)
      PANTK(NZ)  =  PANTK(NZ0)
      DO 90 NS = 1, NSMAX
         PAK(NZ,NS) =  PAK(NZ0,NS)
   90 CONTINUE
C
      RETURN
      END
