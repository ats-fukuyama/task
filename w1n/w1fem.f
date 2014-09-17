C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 2) ******
C
C     ******* BAND MATRIX COEFFICIENT *******
C
      SUBROUTINE W1BNDA(IERR)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: I,J,K,NSF,N,KML,L,NX,N1,N2,M,IND
      REAL(rkind):: RKV,DX
      REAL(rkind):: DS01,DS02,DS11,DS12,DS11A,DS12A,DS21,DS22
      REAL(rkind):: DT11,DT12,DT11A,DT12A
      REAL(rkind),DIMENSION(2,2,2)::  DS0,DS1,DS2,DT0,DT1,DT2,DU0
C
      RKV=2.D6*PI*RF/VC
C
      DO 10 I=1,2
      DO 10 J=1,2
      DO 10 K=1,2
         IF(I.EQ.J.AND.I.EQ.K) THEN
            DS0(I,J,K)=0.25D0
         ELSE
            DS0(I,J,K)=1.D0/12.D0
         ENDIF
         IF(J.EQ.K) THEN
            DS1(I,J,K)=1.D0/3.D0
         ELSE
            DS1(I,J,K)=1.D0/6.D0
         ENDIF
         IF(I.EQ.1) THEN
            DS1(I,J,K)=-DS1(I,J,K)
         ENDIF
         IF(I.EQ.J) THEN
            DS2(I,J,K)= 0.5D0
         ELSE
            DS2(I,J,K)=-0.5D0
         ENDIF
         IF(I.EQ.1) THEN
            IF(J.EQ.K) THEN
               DT0(I,J,K)=1.D0/3.D0
            ELSE
               DT0(I,J,K)=1.D0/6.D0
            ENDIF
         ELSE
            DT0(I,J,K)=0.D0
         ENDIF
         IF(I.EQ.1) THEN
            IF(J.EQ.1) THEN
               DT1(I,J,K)=-0.5D0
            ELSE
               DT1(I,J,K)= 0.5D0
            ENDIF
         ELSE
            DT1(I,J,K)=0.D0
         ENDIF
         IF(I.EQ.1.AND.J.EQ.1) THEN
            DU0(I,J,K)=0.5D0
         ELSE
            DU0(I,J,K)=0.D0
         ENDIF
   10 CONTINUE
C
      NSF=3*NXP+4
      DO 20 I=1,11
      DO 20 N=1,NSF
         CF(I,N)=(0.D0,0.D0)
   20 CONTINUE
      DO 30 N=1,NSF
         CA(N)=(0.D0,0.0D0)
   30 CONTINUE
C
      KML=0
      CF(KML+6,1)=CGIN(1,1)
      CF(KML+7,1)=CGIN(2,1)
      CF(KML+9,1)=(-1.D0,0.D0)
      CF(KML+5,2)=CGIN(1,3)
      CF(KML+6,2)=CGIN(2,3)
      CF(KML+9,2)=(-1.D0,0.D0)
      CF(KML+3,4)=CGIN(1,2)
      CF(KML+4,4)=CGIN(2,2)
      CF(KML+2,5)=CGIN(1,4)
      CF(KML+3,5)=CGIN(2,4)
C
      CF(KML+6,NSF-1)=CGOT(1,1)
      CF(KML+7,NSF-1)=CGOT(2,1)
      CF(KML+4,NSF-1)=(-1.D0,0.D0)
      CF(KML+5,NSF  )=CGOT(1,3)
      CF(KML+6,NSF  )=CGOT(2,3)
      CF(KML+4,NSF  )=(-1.D0,0.D0)
      CF(KML+8,NSF-3)=CGOT(1,2)
      CF(KML+9,NSF-3)=CGOT(2,2)
      CF(KML+7,NSF-2)=CGOT(1,4)
      CF(KML+8,NSF-2)=CGOT(2,4)
      CF(KML+6,NSF-4)=1.D0
C
      CA(1    )=-CGIN(3,1)
      CA(2    )=-CGIN(3,3)
      CA(4    )=-CGIN(3,2)
      CA(5    )=-CGIN(3,4)
      CA(NSF-1)=-CGOT(3,1)
      CA(NSF  )=-CGOT(3,3)
      CA(NSF-3)=-CGOT(3,2)
      CA(NSF-2)=-CGOT(3,4)
C
      DO 6000 I=1,2
      DO 6000 J=1,2
         L=3*(J-I+2)-1
         DS01=DS0(I,J,1)
         DS02=DS0(I,J,2)
         DS11=DS1(I,J,1)+0.5D0*DS1(1,I,J)
         DS12=DS1(I,J,2)+0.5D0*DS1(2,I,J)
         DS11A=DS1(J,I,1)+0.5D0*DS1(1,J,I)
         DS12A=DS1(J,I,2)+0.5D0*DS1(2,J,I)
         DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
         DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
         DT11=DT1(I,J,1)+0.5D0*DT1(I,1,J)
         DT12=DT1(I,J,2)+0.5D0*DT1(I,2,J)
         DT11A=DT1(J,I,1)+0.5D0*DT1(J,1,I)
         DT12A=DT1(J,I,2)+0.5D0*DT1(J,2,I)
C
         DO 5000 NX=1,NXP-1
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            M=3*(NX+I-1)-1
            CF(L+1,M+1)=CF(L+1,M+1)
     &                 +CD0(1,N1)*DU0(I,J,1)*DX
     &                 +CD0(1,N2)*DU0(I,J,2)*DX
            CF(L+2,M+1)=CF(L+2,M+1)
     &                 +CD0(2,N1)*DT0(I,J,1)*DX
     &                 +CD0(2,N2)*DT0(I,J,2)*DX
            CF(L,  M+2)=CF(L,  M+2)
     &                 -CD0(2,N1)*DT0(J,I,1)*DX
     &                 -CD0(2,N2)*DT0(J,I,2)*DX
            CF(L+1,M+2)=CF(L+1,M+2)
     &                 +CD0(3,N1)*DS01*DX
     &                 +CD0(3,N2)*DS02*DX
     &                 +CD2(3,N1)*DS21/DX
     &                 +CD2(3,N2)*DS22/DX
            CF(L+1,M+3)=CF(L+1,M+3)
     &                 +CD0(4,N1)*DS01*DX
     &                 +CD0(4,N2)*DS02*DX
     &                 +CD2(4,N1)*DS21/DX
     &                 +CD2(4,N2)*DS22/DX
            CF(L+3,M+1)=CF(L+3,M+1)
     &                 -CI*CD1(1,N1)*DT11
     &                 -CI*CD1(1,N2)*DT12
            CF(L+2,M+2)=CF(L+2,M+2)
     &                 +CI*CD1(2,N1)*DS11
     &                 +CI*CD1(2,N2)*DS12
            CF(L-1,M+3)=CF(L-1,M+3)
     &                 +CI*CD1(1,N1)*DT11A
     &                 +CI*CD1(1,N2)*DT12A
            CF(L,  M+3)=CF(L,  M+3)
     &                 +CI*CD1(2,N1)*DS11A
     &                 +CI*CD1(2,N2)*DS12A
 5000 CONTINUE
 6000 CONTINUE
C
C     WRITE(6,999) ((L,M,CF(L,M),L=30,41),M=195,197)
C 999 FORMAT((1H ,3(3X,I3,I4,1P2E14.6)))
C
      CALL BANDCD(CF,CA,NSF,11,6*MATLM+5,IND)
         IF(IND.NE.0) WRITE(6,601) IND
      IERR=0
      RETURN
C
  601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
      END
C
C     ******* ELECTROMAGNETIC FIELD IN PLASMA *******
C
      SUBROUTINE W1EPWA(NZ)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NZ
      INTEGER:: NS,NX,I,J,K,N1,N2,M,L
      REAL(rkind):: RKV,RCE,DX,PABSL,PFLXL
      REAL(rkind):: DS01,DS02,DS11,DS12,DS11A,DS12A,DS21,DS22
      REAL(rkind):: DT11,DT12,DT11A,DT12A
      REAL(rkind),DIMENSION(2,2,2)::  DS0,DS1,DS2,DT0,DT1,DT2,DU0
      COMPLEX(rkind):: CM11,CM12,CM13,CM21,CM22,CM23,CM31,CM32,CM33
      COMPLEX(rkind):: CD11,CD12,CD13,CD21,CD22,CD23,CD31,CD32,CD33
      COMPLEX(rkind):: CABSL,CFLXL
C
      RKV=2.D6*PI*RF/VC
C
      RCE=VC*EPS0
C
      DO 100 NS=1,NSMAX
      DO 100 NX=1,NXP
         PABS(NX,NS)=0.D0
  100 CONTINUE
      DO 110 NX=1,NXP
         FLUX(NX)=0.D0
  110 CONTINUE
C
      DO 200 NX=1,NXP
         CE2DA(NZ,NX,1)=CA(3*NX)
         CE2DA(NZ,NX,2)=CA(3*NX+1)
         CE2DA(NZ,NX,3)=CA(3*NX+2)
  200 CONTINUE
C
      DO 10 I=1,2
      DO 10 J=1,2
      DO 10 K=1,2
         IF(I.EQ.J.AND.I.EQ.K) THEN
            DS0(I,J,K)=0.25D0
         ELSE
            DS0(I,J,K)=1.D0/12.D0
         ENDIF
         IF(J.EQ.K) THEN
            DS1(I,J,K)=1.D0/3.D0
         ELSE
            DS1(I,J,K)=1.D0/6.D0
         ENDIF
         IF(I.EQ.1) THEN
            DS1(I,J,K)=-DS1(I,J,K)
         ENDIF
         IF(I.EQ.J) THEN
            DS2(I,J,K)= 0.5D0
         ELSE
            DS2(I,J,K)=-0.5D0
         ENDIF
         IF(I.EQ.1) THEN
            IF(J.EQ.K) THEN
               DT0(I,J,K)=1.D0/3.D0
            ELSE
               DT0(I,J,K)=1.D0/6.D0
            ENDIF
         ELSE
            DT0(I,J,K)=0.D0
         ENDIF
         IF(I.EQ.1) THEN
            IF(J.EQ.1) THEN
               DT1(I,J,K)=-0.5D0
            ELSE
               DT1(I,J,K)= 0.5D0
            ENDIF
         ELSE
            DT1(I,J,K)=0.D0
         ENDIF
         IF(I.EQ.1.AND.J.EQ.1) THEN
            DU0(I,J,K)=0.5D0
         ELSE
            DU0(I,J,K)=0.D0
         ENDIF
   10 CONTINUE
C
      DO 5000 I=1,2
      DO 5000 J=1,2
         DS01=DS0(I,J,1)
         DS02=DS0(I,J,2)
         DS11=DS1(I,J,1)+0.5D0*DS1(1,I,J)
         DS12=DS1(I,J,2)+0.5D0*DS1(2,I,J)
         DS11A=DS1(J,I,1)+0.5D0*DS1(1,J,I)
         DS12A=DS1(J,I,2)+0.5D0*DS1(2,J,I)
         DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
         DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
         DT11=DT1(I,J,1)+0.5D0*DT1(I,1,J)
         DT12=DT1(I,J,2)+0.5D0*DT1(I,2,J)
         DT11A=DT1(J,I,1)+0.5D0*DT1(J,1,I)
         DT12A=DT1(J,I,2)+0.5D0*DT1(J,2,I)
C
         DO 5000 NS=1,NSMAX
         DO 5000 NX=1,NXP-1
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            M=3*(NX+I-1)-1
            L=3*(NX+J-1)-1
            CM11       =+CM0(1,N1,NS)*DU0(I,J,1)*DX
     &                  +CM0(1,N2,NS)*DU0(I,J,2)*DX
            CM12       =+CM0(2,N1,NS)*DT0(I,J,1)*DX
     &                  +CM0(2,N2,NS)*DT0(I,J,2)*DX
            CM21       =-CM0(2,N1,NS)*DT0(J,I,1)*DX
     &                  -CM0(2,N2,NS)*DT0(J,I,2)*DX
            CM22       =+CM0(3,N1,NS)*DS01*DX
     &                  +CM0(3,N2,NS)*DS02*DX
     &                  +CM2(3,N1,NS)*DS21/DX
     &                  +CM2(3,N2,NS)*DS22/DX
            CM33       =+CM0(4,N1,NS)*DS01*DX
     &                  +CM0(4,N2,NS)*DS02*DX
     &                  +CM2(4,N1,NS)*DS21/DX
     &                  +CM2(4,N2,NS)*DS22/DX
            CM13       =-CI*CM1(1,N1,NS)*DT11
     &                  -CI*CM1(1,N2,NS)*DT12
            CM23       =+CI*CM1(2,N1,NS)*DS11
     &                  +CI*CM1(2,N2,NS)*DS12
            CM31       =+CI*CM1(1,N1,NS)*DT11A
     &                  +CI*CM1(1,N2,NS)*DT12A
            CM32       =+CI*CM1(2,N1,NS)*DS11A
     &                  +CI*CM1(2,N2,NS)*DS12A
C
         CABSL=DCONJG(CA(M+1))*(CM11*CA(L+1)+CM12*CA(L+2)+CM13*CA(L+3))
     &        +DCONJG(CA(M+2))*(CM21*CA(L+1)+CM22*CA(L+2)+CM23*CA(L+3))
     &        +DCONJG(CA(M+3))*(CM31*CA(L+1)+CM32*CA(L+2)+CM33*CA(L+3))
         PABSL=-CI*RCE*CABSL
         IF(I.EQ.1.AND.J.EQ.1) THEN
            PABS(NX  ,NS)=PABS(NX  ,NS)+PABSL
         ELSEIF(I.EQ.2.AND.J.EQ.2) THEN
            PABS(NX+1,NS)=PABS(NX+1,NS)+PABSL
         ELSE
            PABS(NX  ,NS)=PABS(NX  ,NS)+0.5D0*PABSL
            PABS(NX+1,NS)=PABS(NX+1,NS)+0.5D0*PABSL
         ENDIF
 5000 CONTINUE
C
      DO 8000 I=1,2
      DO 8000 J=1,2
         DS21=DS1(J,I,1)+0.25D0*DS1(1,I,J)
         DS22=DS1(J,I,2)+0.25D0*DS1(2,I,J)
         DS11=DS0(J,I,1)
         DS12=DS0(J,I,2)
C
         DO 7000 NX=1,NXP-1
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            L=3*(NX+I-2)+2
            M=3*(NX+J-2)+2
            CD11=+CD2(1,N1)*DS21
     &           +CD2(1,N2)*DS22
            CD12=+CD2(2,N1)*DS21
     &           +CD2(2,N2)*DS22
            CD21=-CD2(2,N1)*DS21
     &           -CD2(2,N2)*DS22
            CD22=+CD2(3,N1)*DS21
     &           +CD2(3,N2)*DS22
            CD33=+CD2(4,N1)*DS21
     &           +CD2(4,N2)*DS22
            CD13=0.D0
            CD23=+CI*CD1(2,N1)*DS11*DX
     &           +CI*CD1(2,N2)*DS12*DX
            CD31=+CI*CD1(1,N1)*DS11*DX
     &           +CI*CD1(1,N2)*DS12*DX
            CD32=0.D0
C
         CFLXL=DCONJG(CA(L+1))*(CD11*CA(M+1)+CD12*CA(M+2)+CD13*CA(M+3))
     &        +DCONJG(CA(L+2))*(CD21*CA(M+1)+CD22*CA(M+2)+CD23*CA(M+3))
     &        +DCONJG(CA(L+3))*(CD31*CA(M+1)+CD32*CA(M+2)+CD33*CA(M+3))
         PFLXL=CI*RCE*CFLXL/DX
         FLUX(NX)     =FLUX(NX)     +PFLXL
 7000 CONTINUE
 8000 CONTINUE
      FLUX(NXP)=FLUX(NXP-1)
      RETURN
      END
C
C     ******* BAND MATRIX COEFFICIENT *******
C
      SUBROUTINE W1BNDC(IERR)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: I,J,K,NSF,N,KML,L,NX,N1,N2,M,IND
      REAL(rkind):: RKV,DX
      REAL(rkind):: DS01,DS02,DS11,DS12,DS21,DS22
      REAL(rkind):: DT11,DT12
      REAL(rkind),DIMENSION(2,2,2)::  DS0,DS1,DS2
C
      RKV=2.D6*PI*RF/VC
C
      DO 10 I=1,2
      DO 10 J=1,2
      DO 10 K=1,2
         IF(I.EQ.J.AND.I.EQ.K) THEN
            DS0(I,J,K)=0.25D0
         ELSE
            DS0(I,J,K)=1.D0/12.D0
         ENDIF
         IF(J.EQ.K) THEN
            DS1(I,J,K)=1.D0/3.D0
         ELSE
            DS1(I,J,K)=1.D0/6.D0
         ENDIF
         IF(I.EQ.1) THEN
            DS1(I,J,K)=-DS1(I,J,K)
         ENDIF
         IF(I.EQ.J) THEN
            DS2(I,J,K)= 0.5D0
         ELSE
            DS2(I,J,K)=-0.5D0
         ENDIF
   10 CONTINUE
C
      NSF=3*NXP+4
      DO 20 I=1,11
      DO 20 N=1,NSF
         CF(I,N)=(0.D0,0.D0)
   20 CONTINUE
      DO 30 N=1,NSF
         CA(N)=(0.D0,0.0D0)
   30 CONTINUE
C
      KML=0
      CF(KML+6,1)=CGIN(1,1)
      CF(KML+7,1)=CGIN(2,1)
      CF(KML+9,1)=(-1.D0,0.D0)
      CF(KML+5,2)=CGIN(1,3)
      CF(KML+6,2)=CGIN(2,3)
      CF(KML+9,2)=(-1.D0,0.D0)
      CF(KML+3,4)=CGIN(1,2)
      CF(KML+4,4)=CGIN(2,2)
      CF(KML+2,5)=CGIN(1,4)
      CF(KML+3,5)=CGIN(2,4)
C
      CF(KML+6,NSF-1)=CGOT(1,1)
      CF(KML+7,NSF-1)=CGOT(2,1)
      CF(KML+4,NSF-1)=(-1.D0,0.D0)
      CF(KML+5,NSF  )=CGOT(1,3)
      CF(KML+6,NSF  )=CGOT(2,3)
      CF(KML+4,NSF  )=(-1.D0,0.D0)
      CF(KML+8,NSF-3)=CGOT(1,2)
      CF(KML+9,NSF-3)=CGOT(2,2)
      CF(KML+7,NSF-2)=CGOT(1,4)
      CF(KML+8,NSF-2)=CGOT(2,4)
C
      CA(1    )=-CGIN(3,1)
      CA(2    )=-CGIN(3,3)
      CA(4    )=-CGIN(3,2)
      CA(5    )=-CGIN(3,4)
      CA(NSF-1)=-CGOT(3,1)
      CA(NSF  )=-CGOT(3,3)
      CA(NSF-3)=-CGOT(3,2)
      CA(NSF-2)=-CGOT(3,4)
C
      DO 6000 I=1,2
      DO 6000 J=1,2
         L=3*(J-I+2)-1
         DS01=DS0(I,J,1)
         DS02=DS0(I,J,2)
         DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
         DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
         DS11=(DS1(J,I,1)+0.5D0*DS1(1,I,J))
         DS12=(DS1(J,I,2)+0.5D0*DS1(2,I,J))
         DT11=(DS1(I,J,1)+0.5D0*DS1(1,I,J))
         DT12=(DS1(I,J,2)+0.5D0*DS1(2,I,J))
C
         DO 5000 NX=1,NXP-1
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            M=3*(NX+I-1)-1
            CF(L+1,M+1)=CF(L+1,M+1)
     &                 +CD0(1,N1)*DS01*DX
     &                 +CD0(1,N2)*DS02*DX
     &                 +CD2(1,N1)*DS21/DX
     &                 +CD2(1,N2)*DS22/DX
            CF(L+2,M+1)=CF(L+2,M+1)
     &                 +CD0(2,N1)*DS01*DX
     &                 +CD0(2,N2)*DS02*DX
     &                 +CD2(2,N1)*DS21/DX
     &                 +CD2(2,N2)*DS22/DX
            CF(L,  M+2)=CF(L,  M+2)
     &                 -CD0(2,N1)*DS01*DX
     &                 -CD0(2,N2)*DS02*DX
     &                 -CD2(2,N1)*DS21/DX
     &                 -CD2(2,N2)*DS22/DX
            CF(L+1,M+2)=CF(L+1,M+2)
     &                 +CD0(3,N1)*DS01*DX
     &                 +CD0(3,N2)*DS02*DX
     &                 +CD2(3,N1)*DS21/DX
     &                 +CD2(3,N2)*DS22/DX
            CF(L+1,M+3)=CF(L+1,M+3)
     &                 +CD0(4,N1)*DS01*DX
     &                 +CD0(4,N2)*DS02*DX
     &                 +CD2(4,N1)*DS21/DX
     &                 +CD2(4,N2)*DS22/DX
            CF(L+3,M+1)=CF(L+3,M+1)
     &                 -CI*CD1(1,N1)*DS11
     &                 -CI*CD1(1,N2)*DS12
            CF(L+2,M+2)=CF(L+2,M+2)
     &                 +CI*CD1(2,N1)*DT11
     &                 +CI*CD1(2,N2)*DT12
            CF(L-1,M+3)=CF(L-1,M+3)
     &                 +CI*CD1(1,N1)*DT11
     &                 +CI*CD1(1,N2)*DT12
            CF(L,  M+3)=CF(L,  M+3)
     &                 +CI*CD1(2,N1)*DS11
     &                 +CI*CD1(2,N2)*DS12
 5000 CONTINUE
 6000 CONTINUE
C
C     WRITE(6,999) ((L,M,CF(L,M),L=30,41),M=195,197)
C 999 FORMAT((1H ,3(3X,I3,I4,1P2E14.6)))
C
      CALL BANDCD(CF,CA,NSF,11,6*MATLM+5,IND)
         IF(IND.NE.0) WRITE(6,601) IND
      IERR=0
      RETURN
C
  601 FORMAT(1H ,'!! ERROR IN BANDCD : IND = ',I5)
      END
C
C     ******* ELECTROMAGNETIC FIELD IN PLASMA *******
C
      SUBROUTINE W1EPWC(NZ)
C
      USE w1comm
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NZ
      INTEGER:: NS,NX,I,J,K,N1,N2,M,L
      REAL(rkind):: RKV,RCE,DX,PABSL,PFLXL
      REAL(rkind):: DS01,DS02,DS11,DS12,DS21,DS22
      REAL(rkind):: DT11,DT12
      REAL(rkind),DIMENSION(2,2,2)::  DS0,DS1,DS2
      COMPLEX(rkind):: CM11,CM12,CM13,CM21,CM22,CM23,CM31,CM32,CM33
      COMPLEX(rkind):: CD11,CD12,CD13,CD21,CD22,CD23,CD31,CD32,CD33
      COMPLEX(rkind):: CABSL,CFLXL
C
      RKV=2.D6*PI*RF/VC
C
      RCE=VC*EPS0
C
      DO 100 NS=1,NSMAX
      DO 100 NX=1,NXP
         PABS(NX,NS)=0.D0
  100 CONTINUE
      DO 110 NX=1,NXP
         FLUX(NX)=0.D0
  110 CONTINUE
C
      DO 200 NX=1,NXP
         CE2DA(NZ,NX,1)=CA(3*NX)
         CE2DA(NZ,NX,2)=CA(3*NX+1)
         CE2DA(NZ,NX,3)=CA(3*NX+2)
  200 CONTINUE
C
      DO 10 I=1,2
      DO 10 J=1,2
      DO 10 K=1,2
         IF(I.EQ.J.AND.I.EQ.K) THEN
            DS0(I,J,K)=0.25D0
         ELSE
            DS0(I,J,K)=1.D0/12.D0
         ENDIF
         IF(J.EQ.K) THEN
            DS1(I,J,K)=1.D0/3.D0
         ELSE
            DS1(I,J,K)=1.D0/6.D0
         ENDIF
         IF(I.EQ.1) THEN
            DS1(I,J,K)=-DS1(I,J,K)
         ENDIF
         IF(I.EQ.J) THEN
            DS2(I,J,K)= 0.5D0
         ELSE
            DS2(I,J,K)=-0.5D0
         ENDIF
   10 CONTINUE
C
      DO 5000 I=1,2
      DO 5000 J=1,2
         DS01=DS0(I,J,1)
         DS02=DS0(I,J,2)
         DS21=DS2(I,J,1)+0.25D0*(DS2(I,1,J)+DS2(J,1,I))
         DS22=DS2(I,J,2)+0.25D0*(DS2(I,2,J)+DS2(J,2,I))
         DS11=(DS1(J,I,1)+0.5D0*DS1(1,I,J))
         DS12=(DS1(J,I,2)+0.5D0*DS1(2,I,J))
         DT11=(DS1(I,J,1)+0.5D0*DS1(1,I,J))
         DT12=(DS1(I,J,2)+0.5D0*DS1(2,I,J))
C
         DO 5000 NS=1,NSMAX
         DO 5000 NX=1,NXP-1
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            M=3*(NX+I-1)-1
            L=3*(NX+J-1)-1
            CM11=+CM0(1,N1,NS)*DS01*DX
     &           +CM0(1,N2,NS)*DS02*DX
     &           +CM2(1,N1,NS)*DS21/DX
     &           +CM2(1,N2,NS)*DS22/DX
            CM12=+CM0(2,N1,NS)*DS01*DX
     &           +CM0(2,N2,NS)*DS02*DX
     &           +CM2(2,N1,NS)*DS21/DX
     &           +CM2(2,N2,NS)*DS22/DX
            CM21=-CM0(2,N1,NS)*DS01*DX
     &           -CM0(2,N2,NS)*DS02*DX
     &           -CM2(2,N1,NS)*DS21/DX
     &           -CM2(2,N2,NS)*DS22/DX
            CM22=+CM0(3,N1,NS)*DS01*DX
     &           +CM0(3,N2,NS)*DS02*DX
     &           +CM2(3,N1,NS)*DS21/DX
     &           +CM2(3,N2,NS)*DS22/DX
            CM33=+CM0(4,N1,NS)*DS01*DX
     &           +CM0(4,N2,NS)*DS02*DX
     &           +CM2(4,N1,NS)*DS21/DX
     &           +CM2(4,N2,NS)*DS22/DX
            CM13=-CI*CM1(1,N1,NS)*DS11
     &           -CI*CM1(1,N2,NS)*DS12
            CM23=+CI*CM1(2,N1,NS)*DT11
     &           +CI*CM1(2,N2,NS)*DT12
            CM31=+CI*CM1(1,N1,NS)*DT11
     &           +CI*CM1(1,N2,NS)*DT12
            CM32=+CI*CM1(2,N1,NS)*DS11
     &           +CI*CM1(2,N2,NS)*DS12
C
         CABSL=DCONJG(CA(M+1))*(CM11*CA(L+1)+CM12*CA(L+2)+CM13*CA(L+3))
     &        +DCONJG(CA(M+2))*(CM21*CA(L+1)+CM22*CA(L+2)+CM23*CA(L+3))
     &        +DCONJG(CA(M+3))*(CM31*CA(L+1)+CM32*CA(L+2)+CM33*CA(L+3))
         PABSL=-CI*RCE*CABSL
         IF(I.EQ.1.AND.J.EQ.1) THEN
            PABS(NX  ,NS)=PABS(NX  ,NS)+PABSL
         ELSEIF(I.EQ.2.AND.J.EQ.2) THEN
            PABS(NX+1,NS)=PABS(NX+1,NS)+PABSL
         ELSE
            PABS(NX  ,NS)=PABS(NX  ,NS)+0.5D0*PABSL
            PABS(NX+1,NS)=PABS(NX+1,NS)+0.5D0*PABSL
         ENDIF
 5000 CONTINUE
C
      DO 8000 I=1,2
      DO 8000 J=1,2
         DS21=DS1(J,I,1)+0.25D0*DS1(1,I,J)
         DS22=DS1(J,I,2)+0.25D0*DS1(2,I,J)
         DS11=DS0(J,I,1)
         DS12=DS0(J,I,2)
C
         DO 7000 NX=1,NXP-1
            N1=NX
            N2=NX+1
            DX=RKV*(XA(N2)-XA(N1))
            L=3*(NX+I-2)+2
            M=3*(NX+J-2)+2
            CD11=+CD2(1,N1)*DS21
     &           +CD2(1,N2)*DS22
            CD12=+CD2(2,N1)*DS21
     &           +CD2(2,N2)*DS22
            CD21=-CD2(2,N1)*DS21
     &           -CD2(2,N2)*DS22
            CD22=+CD2(3,N1)*DS21
     &           +CD2(3,N2)*DS22
            CD33=+CD2(4,N1)*DS21
     &           +CD2(4,N2)*DS22
            CD13=0.D0
            CD23=+CI*CD1(2,N1)*DS11*DX
     &           +CI*CD1(2,N2)*DS12*DX
            CD31=+CI*CD1(1,N1)*DS11*DX
     &           +CI*CD1(1,N2)*DS12*DX
            CD32=0.D0
C
         CFLXL=DCONJG(CA(L+1))*(CD11*CA(M+1)+CD12*CA(M+2)+CD13*CA(M+3))
     &        +DCONJG(CA(L+2))*(CD21*CA(M+1)+CD22*CA(M+2)+CD23*CA(M+3))
     &        +DCONJG(CA(L+3))*(CD31*CA(M+1)+CD32*CA(M+2)+CD33*CA(M+3))
         PFLXL=CI*RCE*CFLXL/DX
         FLUX(NX)     =FLUX(NX)     +PFLXL
 7000 CONTINUE
 8000 CONTINUE
      FLUX(NXP)=FLUX(NXP-1)
      RETURN
      END
