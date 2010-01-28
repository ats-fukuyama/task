C
********************************************************************
C                          DPDKDTF
C dielectric tensor with Drift Kinetic effects by Fourier expansion
C
C                         2009/12/15    
C                                     programed by T.OKAMOTO
C*******************************************************************
C
      SUBROUTINE DPDKDTF(CW,CKPR,NS,NR,MM,CLDISP1,CLDISP2,CLDISP3)
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION CLDISP1(9,2*MM+1)
      DIMENSION CLDISP2(3),CLDISP3(3)
C
      CALL DPDKDTFC(CW, CKPR, NS,NR, MM,MATRICSQR, MATRICSQC, MATRICSCP)
      CALL DPDKDTFS(CW, CKPR, NS,NR, MM, CLDISP1,CLDISP3,CLDISP3)
C     
      RETURN
      END SUBROUTINE

C
C     *************************************
C             SET MATRIX 
C     **************************************
C
      SUBROUTINE DPDKDTFC(CW, CKPR, NS, NR, MM, MATRICSQR,
     &     MATRICSQC, MATRICSCP)
C
C------------------------------------------------------------------!
C     CW : FREQUENCY
C     CKPR: WAVE NUMBER
C     NS : PARTICLE SPECIES
C     MM : POLOIDAL NUMBER
C------------------------------------------------------------------!
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION DRFP(NTHMAX,NPMAX)
      DIMENSION CWAST(NTHMAX,NPMAX)
      DIMENSION DPFP(NTHMAX,NPMAX)
      DIMENSION MATRIXA(2*MM+1,2*MM+1)
      DIMENSION MATRIXD(2*MM+1,2*MM+1)
      DIMENSION MATRIXF(2*MM+1,2*MM+1)
      DIMENSION MATRIXG(2*MM+1,2*MM+1)
      DIMENSION MATRIXQR(2*MM+1,2*MM+1)
      DIMENSION MATRIXQC(2*MM+1,2*MM+1)
      DIMENSION MATRIXQP(2*MM+1,2*MM+1)
      DIMENSION TCSM2(NTHMAX)
      DIMENSION TSNM2(NTHMAX)
      DIMENSION PV(NPMAX)
C
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      PV(NP)=PG(NP)/AM
      CKPRX=CKPR*PV(NP)*TCSM(NTH)
      DKPRX=DBLE(CKPRX)
      TCSM2(NTH)=TCSM(NTH)**2
      TSNM2(NTH)=TSNM(NTH)**2
C
C     ****** derivation in the radial direction ******
C
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP(NTH,NP) = (FP(NTH,NP,NR+1)
     &        -FP(NTH,NP,NR-1))/(2*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/2*DELP)/FP(NTH,NP,NR)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
      ENDDO
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/2*DELP)/FP(NTH,NP,NR)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
      ENDDO
C
C
C     ***** SETING MATRIX ******
C
C     ----- Coefficient matrixA -----
      DO i=1,2*MM+1
      DO j=1,2*MM+1
         IF (j.EQ.i-1) THEN
            MATRIXA(i,j)=CI*AM*PV(NP)**2*(TSNM2(NTH)+2*TCSM2(NTH))
     &           *(MM-1)/(4*AE*BABS*RSL)
         ELSE IF (j.EQ.i) THEN
            MATRIXA(i,j)=-CI*(CW-DKPRX)
         ELSE IF (j.EQ.i+1) THEN
            MATRIXA(i,j)=CI*AM*PV(NP)**2*(TSNM2(NTH)+2*TCSM2(NTH))
     &           *(MM+1)/4*AE*BABS*RSL
         ELSE
            MATRIXA(i,j)=0
         END IF
      ENDDO
      ENDDO
C    ------- Get Ainv -----
C
      CALL DGETRI(2*MM+1,MATRIXA, LDA, IPIV, WORK, LWORK, INFO)
C
C    ----- Coefficient matrixD -------
C
      DO i=1,2*MM+1
      DO j=1,2*MM+1
         IF (j.EQ.i-1) THEN
            MATRIXD(i,j)=-AE*FP(NTH,NP,NR)*(DPFP(NTH,NP)
     &           -CWAST(NTH,NP)/CW)*PV(NP)*TCSM(NTH)
         ELSE IF (j.EQ.i+1) THEN
            MATRIXD(i,j)=-AE*FP(NTH,NP,NR)*(DPFP(NTH,NP)
     &           -CWAST(NTH,NP)/CW)*PV(NP)*TCSM(NTH)
         ELSE
            MATRIXD(i,j)=0
         END IF
      ENDDO
      ENDDO
C
C     ----- Coefficient matrixF -----
C 
      DO i=1,2*MM+1
      DO j=1,2*MM+1
         IF (j.EQ.i-1) THEN
            MATRIXF(i,j)=-AE*FP(NTH,NP,NR)*(DPFP(NTH,NP)
     &           -CWAST(NTH,NP)/CW)*PV(NP)**2
     &           *(TSNM2(NTH)+2*TCSM2(NTH))/4
         ELSE IF (j.EQ.i+1) THEN
            MATRIXF(i,j)=-AE*FP(NTH,NP,NR)*(DPFP(NTH,NP)
     &           -CWAST(NTH,NP)/CW)*PV(NP)**2
     &           *(TSNM2(NTH)+2*TCSM2(NTH))/4
         ELSE
            MATRIXF(i,j)=0
         END IF
      ENDDO
      ENDDO
C
C     ------ Coefficient matrixG -----
      DO i=1,2*MM+1
      Do j=1,2*MM+1
         IF (j.EQ.i) THEN
            MATRIXG(i,j)=-AE*FP(NTH,NP,NR)*(DPFP(NTH,NP)-
     &           CWAST(NTH,NP)/CW)*PV(NP)**2
     &           *(TSNM2(NTH)+2*TCSM2(NTH))/4
         ELSE
            MATRIXG(i,j)=-AE*FP(NTH,NP,NR)*(DPFP(NTH,NP)-
     &           CWAST(NTH,NP)/CW)*PV(NP)**2
     &           *(TSNM2(NTH)+2*TCSM2(NTH))/4
         ENDIF
      ENDDO
      ENDDO
C
C    --------  Get QR,QC,QP : Ainv*(D,F,G)=QR,QC,QP-------  
C
      CALL DGEMM ( TRANSA, TRANSM, 2*MM+1, 2*MM+1, 2*MM+1, ALPHA,
     & MATRIXA, LDA, MATRIXD, LDB, BETA, MATRICSQR)
      CALL DGEMM ( TRANSA, TRANSM, 2*MM+1, 2*MM+1, 2*MM+1, ALPHA,
     & MATRIXA, LDA, MATRIXF, LDB, BETA, MATRICSQC)
      CALL DGEMM ( TRANSA, TRANSM, 2*MM+1, 2*MM+1, 2*MM+1, ALPHA,
     & MATRIXA, LDA, MATRIXG, LDB, BETA, MATRICSQP)
C
C
      RETURN
      END SUBROUTINE DPDKDTFC
C
C
C     ***********************************
C         VELOCITY INTEGRATION 
C     ***********************************
C
      SUBROUTINE DPDKDTFS(CW,CKPR,NS,NR,MM,CLDISP1,CLDISP2M,CLDISP3M)
C
C------------------------------------------------------------------!
C     CW : FREQUENCY
C     CKPR: WAVE NUMBER
C     NS : PARTICLE SPECIES
C     MM : POLOIDAL NUMBER
C------------------------------------------------------------------!
C
      INCLUDE 'dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'

      DIMENSION CLDISP1(9,2*MM+1)
      DIMENSION CLDISPRCP(9,2*MM+1)
      DIMENSION CLDISP2M(3)
      DIMENSION CLDISP3M(3)
C
      DIMENSION DRFP(NTHMAX,NPMAX)
      DIMENSION CWAST(NTHMAX,NPMAX)
      DIMENSION DPFP(NTHMAX,NPMAX)
      DIMENSION MATRIXQR(2*MM+1,2*MM+1)
      DIMENSION MATRIXQC(2*MM+1,2*MM+1)
      DIMENSION MATRIXQP(2*MM+1,2*MM+1)
      DIMENSION TCSM2(NTHMAX)
      DIMENSION TSNM2(NTHMAX)
      DIMENSION PV(NPMAX)
C
C
C     ****** derivation in the radial direction ******
C
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
         DRFP(NTH,NP) = (FP(NTH,NP,NR+1)
     &        -FP(NTH,NP,NR-1))/(2*DELR)
         CWAST(NTH,NP)=MM*(DRFP(NTH,NP)/FP(NTH,NP,NR))*BABS/(AA*RSL)
      ENDDO
      ENDDO
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/2*DELP)/FP(NTH,NP,NR)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
      ENDDO
C
      DO NTH=1,NTHMAX
         IF (NP == 1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &           -FP(NTH,NP,NR))/DELP)/FP(NTH,NP,NR)
         ELSE IF (NP.GE.2 .AND. NP.LE.NPMAX-2) THEN
            DO NP=2,NPMAX-2
               DPFP(NTH,NP) = ((FP(NTH,NP+1,NR)
     &              - FP(NTH,NP-1,NR))/2*DELP)/FP(NTH,NP,NR)
            ENDDO
         ELSE IF (NP == NPMAX-1) THEN
            DPFP(NTH,NP) = ((FP(NTH,NP,NR)
     &           -FP(NTH,NP-1,NR))/DELP)/FP(NTH,NP,NR)
         END IF
      ENDDO
C
C
      CLDISPRCP(1,:) = (0.D0,0.D0)
      CLDISPRCP(2,:) = (0.D0,0.D0)
      CLDISPRCP(3,:) = (0.D0,0.D0)
      CLDISPRCP(4,:) = (0.D0,0.D0)
      CLDISPRCP(5,:) = (0.D0,0.D0)
      CLDISPRCP(6,:) = (0.D0,0.D0)
      CLDISPRCP(7,:) = (0.D0,0.D0)
      CLDISPRCP(8,:) = (0.D0,0.D0)
      CLDISPRCP(9,:) = (0.D0,0.D0)
C
      CLDISPRCP(1,:) = MATRIXQR(MM+1,:)+MATRIXQR(MM-1,:)
      CLDISPRCP(2,:) = MATRIXQC(MM+1,:)+MATRIXQC(MM-1,:)
      CLDISPRCP(3,:) = MATRIXQP(MM+1,:)+MATRIXQP(MM-1,:)
      CLDISPRCP(4,:) = CLDISPRCP(1,:)
      CLDISPRCP(5,:) = CLDISPRCP(2,:)
      CLDISPRCP(6,:) = CLDISPRCP(3,:)
      CLDISPRCP(7,:) = MATRIXQR(MM,:)
      CLDISPRCP(8,:) = MATRIXQC(MM,:)
      CLDISPRCP(9,:) = MATRIXQP(MM,:)
C
C
      CLDISP2M(1) = 0.D0
      CLDISP2M(2) = 0.D0
      CLDISP2M(3) = 0.D0
C
C
      CLDISP3M(1) = 0.D0
      CLDISP3M(2) = 0.D0
      CLDISP3M(3) = 0.D0
C
C
      DO NP=1,NPMAX-1
      DO NTH=1,NTHMAX
C
         AM=PA(NS)*AMP
         AE=PZ(NS)*AEE
         PV(NP)=PG(NP)/AM
         CKPRX=CKPR*PV(NP)*TCSM(NTH)
         DKPRX=DBLE(CKPRX)
         TCSM2(NTH)=TCSM(NTH)**2
         TSNM2(NTH)=TSNM(NTH)**2
C
         COE1=PV(NP)**4*(TSNM2(NTH)+2*TCSM2(NTH))*TSNM(NTH)*DELP*DELTH
         COE2=PV(NP)**3*TCSM(NTH)*TSNM(NTH)*DELP*DELTH
         COE3=PV(NP)**2*TSNM(NTH)*DELP*DELTH
C
         CLDISP1(1,:) = CLDISP1(1,:) + COE1*CLDISPRCP(1,:)
         CLDISP1(2,:) = CLDISP1(2,:) + COE1*CLDISPRCP(2,:)
         CLDISP1(3,:) = CLDISP1(3,:) + COE1*CLDISPRCP(3,:) 
         CLDISP1(4,:) = CLDISP1(4,:) + COE1*CLDISPRCP(4,:)
         CLDISP1(5,:) = CLDISP1(5,:) + COE1*CLDISPRCP(5,:)
         CLDISP1(6,:) = CLDISP1(6,:) + COE1*CLDISPRCP(6,:)
         CLDISP1(7,:) = CLDISP1(7,:) + COE2*CLDISPRCP(7,:)
         CLDISP1(8,:) = CLDISP1(8,:) + COE2*CLDISPRCP(8,:)
         CLDISP1(9,:) = CLDISP1(9,:) + COE2*CLDISPRCP(9,:)
C
         CLDISP2M(1) = CLDISP2M(1) + DRFP(NTH,NP)*COE1
         CLDISP2M(2) = CLDISP2M(2) + DRFP(NTH,NP)*COE1
         CLDISP2M(3) = CLDISP2M(3) + DRFP(NTH,NP)*COE2
C
         CLDISP3M(1) = CLDISP3M(1) + FP(NTH,NP,NR)*COE3
         CLDISP3M(2) = CLDISP3M(2) + FP(NTH,NP,NR)*COE3
         CLDISP3M(3) = CLDISP3M(3) + FP(NTH,NP,NR)*COE3
C
      ENDDO
      ENDDO
C
      COE4=PI*AM/(2*BARS*RR)
      COE5=2*PI*AE
C
      CLDISP1(1,:) = COE4*CLDISP1(1,:)
      CLDISP1(2,:) = COE4*CLDISP1(2,:)
      CLDISP1(3,:) = COE4*CLDISP1(3,:)
      CLDISP1(4,:) = -CI*CLDISP1(1,:)
      CLDISP1(5,:) = -CI*CLDISP1(2,:)
      CLDISP1(6,:) = -CI*CLDISP1(3,:)
      CLDISP1(7,:) = COE5*CLDISP1(7,:)
      CLDISP1(8,:) = COE5*CLDISP1(8,:)
      CLDISP1(9,:) = COE5*CLDISP1(9,:)
C
      COE6=PI*AM/(2*CW*RR)
      COE7=2*PI*AE*BARS/CW
C
      CLDISP2M(1) = -COE6*CLDISP2M(1)
      CLDISP2M(2) = COE6*CLDISP2M(2)
      CLDISP2M(3) = -CI*COE7*CLDISP2M(3)
C
      COE8=2*PI*AE*BARS
C
      CLDISP3M(1) = -CI*COE8*CLDISP3M(1)
      CLDISP3M(2) = COE8*CLDISP3M(2)
      CLDISP3M(3) = COE8*CLDISP3M(3)
C
C
      RETURN
      END SUBROUTINE DPDKDTFS
C
C
