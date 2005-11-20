C     $Id$
C
C     ****** FREE BOUNDARY EQUILIBRIUM SOLVER ******
C
      SUBROUTINE EQCALX(ID,IERR)
C
      INCLUDE '../eq/eqcomx.inc'
C
      IERR=0
C
      MWMAX=8*(NRGMAX+2)-1
      MLMAX=4*NRGMAX*NZGMAX
      NBND=4*(NRGMAX+2)
C
      DO NLOOP=1,NLPMAX
         CALL EQCALFMA
         CALL EQCALFVB(ID)
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZOLD(NRG,NZG)=PSIRZ(NRG,NZG)
         ENDDO
         ENDDO
C
         CALL BANDRD(FMA,FVB,MLMAX,MWMAX,MWM,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'BANDRD ERROR: IERR =',IERR
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            N=4*((NZG-1)*NRGMAX+NRG-1)+1
            PSIRZ(NRG,NZG)=FVB(N)
         ENDDO
         ENDDO
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            DELPSIRZ(NRG,NZG)=PSIRZ(NRG,NZG)-PSIRZOLD(NRG,NZG)
            HJTRZ(NRG,NZG)=FJRZ(NRG,NZG)
         ENDDO
         ENDDO
C
         PSI0=PSIRZ(1,1)
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSI0=MIN(PSI0,PSIRZ(NRG,NZG))
         ENDDO
         ENDDO
         PSIPA=-PSI0
C
         SUM=0.D0
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            SUM=SUM+DELPSIRZ(NRG,NZG)**2
         ENDDO
         ENDDO
         WRITE(6,*) 'SUM=',SUM
         IF(SUM.LT.EPSEQ) GOTO 1000
      ENDDO
      IF(ID.EQ.0) GOTO 9000
      WRITE(6,*) 'XX NO CONVERGENCE IN EQCALX'
 1000 CONTINUE
C
      DPS=1.D0/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPNL=DPS*(NPS-1)
         PSIPS(NPS)=PSIPA*PSIPNL
         CALL EQPPSI(PSIPNL,PPSI,DPPSI)
         CALL EQFPSI(PSIPNL,FPSI,DFPSI)
         PPPS(NPS)=PPSI
         TTPS(NPS)=SQRT(BB**2*RR**2+TJ*FPSI*DFPSI)
         WRITE(6,'(I5,1P3E12.4)')
     &        NPS,PSIPNL,PPSI,FPSI
      ENDDO

C
 9000 RETURN
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
      DIMENSION K1AL(4,4),M1AL(4,4),L1AL(4,4),N1AL(4,4)
C
      DATA K1AL/1,1,2,2,
     &          1,1,2,2,
     &          3,3,4,4,
     &          3,3,4,4/
      DATA M1AL/1,2,1,2,
     &          3,4,3,4,
     &          1,2,1,2,
     &          3,4,3,4/
      DATA L1AL/1,1,2,2,
     &          1,1,2,2,
     &          3,3,4,4,
     &          3,3,4,4/
      DATA N1AL/1,2,1,2,
     &          3,4,3,4,
     &          1,2,1,2,
     &          3,4,3,4/
C
      DO J=1,4
      DO I=1,4
         K1A(I,J)=K1AL(I,J)
         L1A(I,J)=L1AL(I,J)
         M1A(I,J)=M1AL(I,J)
         N1A(I,J)=N1AL(I,J)
      ENDDO
      ENDDO

C
      RK(1,1,1)= 6.D0/ 5.D0
      RK(1,1,2)=-6.D0/ 5.D0
      RK(1,2,1)=-6.D0/ 5.D0
      RK(1,2,2)= 6.D0/ 5.D0
      RK(2,1,1)=-1.D0/10.D0
      RK(2,1,2)=-1.D0/10.D0
      RK(2,2,1)= 1.D0/10.D0
      RK(2,2,2)= 1.D0/10.D0
      RK(3,1,1)=-1.D0/10.D0
      RK(3,1,2)= 1.D0/10.D0
      RK(3,2,1)=-1.D0/10.D0
      RK(3,2,2)= 1.D0/10.D0
      RK(4,1,1)= 2.D0/15.D0
      RK(4,1,2)=-1.D0/30.D0
      RK(4,2,1)=-1.D0/30.D0
      RK(4,2,2)= 2.D0/15.D0
C
      RL(1,1,1)= 13.D0/ 35.D0
      RL(1,1,2)=  9.D0/ 70.D0
      RL(1,2,1)=  9.D0/ 70.D0
      RL(1,2,2)= 13.D0/ 35.D0
      RL(2,1,1)=-11.D0/210.D0
      RL(2,1,2)= 13.D0/420.D0
      RL(2,2,1)=-13.D0/420.D0
      RL(2,2,2)= 11.D0/210.D0
      RL(3,1,1)=-11.D0/210.D0
      RL(3,1,2)=-13.D0/420.D0
      RL(3,2,1)= 13.D0/420.D0
      RL(3,2,2)= 11.D0/210.D0
      RL(4,1,1)=  1.D0/105.D0
      RL(4,1,2)=- 1.D0/140.D0
      RL(4,2,1)=- 1.D0/140.D0
      RL(4,2,2)=  1.D0/105.D0
C
      DO J=1,2
      DO I=1,2
      DO K=1,4
         RH(K,I,J)= RL(K,I,J)
         RM(K,I,J)= RL(K,I,J)
         RN(K,I,J)= RK(K,I,J)
      ENDDO
      ENDDO
      ENDDO
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
      DO N=1,MLMAX
      DO M=1,MWMAX
         FMA(M,N)=0
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX-1
      DO NRG=1,NRGMAX-1
         RMID=0.5D0*(RG(NRG)+RG(NRG+1))
         DRG=RG(NRG)-RG(NRG+1)
         FACTK(1)=1.D0/DRG**2
         FACTK(2)=1.D0/DRG
         FACTK(3)=1.D0/DRG
         FACTK(4)=1.D0
         FACTL(1)=1.D0
         FACTL(2)=1.D0*DRG
         FACTL(3)=1.D0*DRG
         FACTL(4)=1.D0*DRG**2
C
         DZG=ZG(NZG)-ZG(NZG+1)
         FACTM(1)=1.D0
         FACTM(2)=1.D0*DZG
         FACTM(3)=1.D0*DZG
         FACTM(4)=1.D0*DZG**2
         FACTN(1)=1.D0/DZG**2
         FACTN(2)=1.D0/DZG
         FACTN(3)=1.D0/DZG
         FACTN(4)=1.D0
C
         DO I1=1,2
         DO I2=1,2
         DO J1=1,2
         DO J2=1,2
         DO K=1,4
         DO L=1,4
            N=4*((NZG-1)*NRGMAX+NRG-1)+K+4*(I1-1)+4*NRGMAX*(J1-1)
            M=4*((NZG-1)*NRGMAX+NRG-1)+L+4*(I2-1)+4*NRGMAX*(J2-1)
C
            K1=K1A(K,L)
            M1=M1A(K,L)
            L1=L1A(K,L)
            N1=N1A(K,L)
            FMA(NBND+M-N,N)=FMA(NBND+M-N,N)+RK(K1,I1,I2)/RMID*FACTK(K1)
     &                                     *RM(M1,J1,J2)     *FACTM(M1)
     &                                     +RL(L1,I1,I2)/RMID*FACTL(L1)
     &                                     *RN(N1,J1,J2)     *FACTN(N1) 
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      DO NZG=1,NZGMAX,NZGMAX-1
      DO NRG=1,NRGMAX
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         DO MW=1,MWMAX
            FMA(MW,N)=0.D0
         ENDDO
         FMA(NBND,N)=1.D0
      ENDDO
      ENDDO
C
      DO NRG=1,NRGMAX,NRGMAX-1
      DO NZG=2,NZGMAX-1
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         DO MW=1,MWMAX
            FMA(MW,N)=0.D0
         ENDDO
         FMA(NBND,N)=1.D0
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
C
C     **********************************
C     *     Calulate FEM RHS vector    *
C     **********************************
C
      SUBROUTINE EQCALFVB(ID)
C
      INCLUDE '../eq/eqcomx.inc'
C
      DIMENSION FACTM(4),FACTH(4)
      DIMENSION RJ(4),PSIBRZ(4)
      DIMENSION HJ0(NRGM,NZGM),HJ1(NRGM,NZGM)
C
      DRG=(RGMAX-RGMIN)/(NRGMAX-1)
      DZG=(ZGMAX-ZGMIN)/(NZGMAX-1)
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
        FJRZ(NRG,NZG)=0.D0
      ENDDO
      ENDDO
C
      IF(ID.EQ.0) THEN
         FJRZ((NRGMAX+1)/2,(NZGMAX+1)/2)=-RIP/(DRG*DZG)
         ID=1
      ELSEIF(ID.EQ.1) THEN
         PSI0=PSIRZ(1,1)
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSI0=MIN(PSI0,PSIRZ(NRG,NZG))
         ENDDO
         ENDDO
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIRZ(NRG,NZG)*PSI0.LT.0.D0) THEN
               PSIPNL=PSIRZ(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               PPRZ(NRG,NZG)=PPSI
            ELSE
               PPRZ(NRG,NZG)=0.D0
            ENDIF
         ENDDO
         ENDDO
C
         FJ0=0.D0
         FJ1=0.D0
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            IF(PSIRZ(NRG,NZG)*PSI0.GT.0.D0.AND.
     &         ZG(NZG).LE.ZLIMP.AND.
     &         ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=PSIRZ(NRG,NZG)/PSI0
               CALL EQPPSI(PSIPNL,PPSI,DPPSI)
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               HJ0(NRG,NZG)= DPPSI
               HJ1(NRG,NZG)= FPSI*DFPSI
            ELSE
               HJ0(NRG,NZG)=0.D0
               HJ1(NRG,NZG)=0.D0
            ENDIF
C
            DVOL=DRG*DZG
            FJ0=FJ0+HJ0(NRG,NZG)*DVOL
            FJ1=FJ1+HJ1(NRG,NZG)*DVOL
         ENDDO
         ENDDO
C
         TJ=(RIP-FJ0)/FJ1
C
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            FJRZ(NRG,NZG)=HJ0(NRG,NZG)+TJ*HJ1(NRG,NZG)
            IF(PSIRZ(NRG,NZG).GT.0.D0.AND.
     &         ZG(NZG).LE.ZLIMP.AND.
     &         ZG(NZG).GE.ZLIMM) THEN
               PSIPNL=PSIRZ(NRG,NZG)/PSI0
               CALL EQFPSI(PSIPNL,FPSI,DFPSI)
               TTRZ(NRG,NZG)=SQRT(BB**2*RR**2+TJ*FPSI*DFPSI)
            ELSE
               TTRZ(NRG,NZG)=BB*RR
            ENDIF
         ENDDO
         ENDDO
      ENDIF
C
      DO NPFC=1,NPFCMAX
         NZC=INT((ZPFC(NPFC)-ZGMIN)/DZG)+1
         IF(NZC.GE.1.AND.NZC.LT.NZGMAX) THEN
C            DO NRG=1,NRGMAX
C            DO NZG=1,NZC+2
C               FJRZ(NRG,NZG)=0.D0
C            ENDDO
C            ENDDO
            DO I=1,3
               IF(I.EQ.1) THEN
                  RDCL=RPFC(NPFC)
                  RIDCL=RIPFC(NPFC)
               ELSEIF(I.EQ.2) THEN
                  RDCL=RPFC(NPFC)-WPFC(NPFC)
                  RIDCL=-0.5D0*RIPFC(NPFC)
               ELSEIF(I.EQ.3) THEN
                  RDCL=RPFC(NPFC)+WPFC(NPFC)
                  RIDCL=-0.5D0*RIPFC(NPFC)
               ENDIF
               NRC=INT((RDCL-RGMIN)/DRG)+1
               IF(NRC.GE.1.AND.NRC.LT.NRGMAX) THEN
                  FACTR=(RDCL-RG(NRC))/DRG
                  FACTZ=(ZPFC(NPFC)-ZG(NZC))/DZG
                  FACTRC=1.D0-FACTR
                  FACTZC=1.D0-FACTZ
C
                  FJRZ(NRC  ,NZC  )=FJRZ(NRC  ,NZC  )
     &                             +FACTRC*FACTZC*RIPFC(NPFC)/(DRG*DZG)
                  FJRZ(NRC+1,NZC  )=FJRZ(NRC+1,NZC  )
     &                             +FACTR *FACTZC*RIPFC(NPFC)/(DRG*DZG)
                  FJRZ(NRC  ,NZC+1)=FJRZ(NRC  ,NZC+1)
     &                             +FACTRC*FACTZ *RIPFC(NPFC)/(DRG*DZG)
                  FJRZ(NRC+1,NZC+1)=FJRZ(NRC+1,NZC+1)
     &                             +FACTR *FACTZ *RIPFC(NPFC)/(DRG*DZG)
               ELSE
                  WRITE(6,*) 'XX (RDC) OUT OF REGION'
               ENDIF
            ENDDO
         ELSE
            WRITE(6,*) 'XX (ZPFC) OUT OF REGION: ZPFC=',ZPFC(NPFC)
         ENDIF
      ENDDO
C
      DO N=1,MLMAX
         FVB(N)=0
      ENDDO
C
      DO NZG=1,NZGMAX-1
      DO NRG=1,NRGMAX-1
         DRG=RG(NRG)-RG(NRG+1)
         FACTH(1)=1.D0
         FACTH(2)=1.D0*DRG
         FACTH(3)=1.D0*DRG
         FACTH(4)=1.D0*DRG**2
C
         DZG=ZG(NZG)-ZG(NZG+1)
         FACTM(1)=1.D0
         FACTM(2)=1.D0*DZG
         FACTM(3)=1.D0*DZG
         FACTM(4)=1.D0*DZG**2
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
            RJ(4)=0
            DO K=1,4
            DO L=1,4
               N=4*((NZG-1)*NRGMAX+NRG-1)+K+4*(I1-1)+4*NRGMAX*(J1-1)
C
               K1=K1A(K,L)
               M1=M1A(K,L)
               FVB(N)=FVB(N)+RH(K1,I1,I2)*FACTH(K1)
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
C      WRITE(6,*) (FVB(ML),ML=1,MLMAX)
C
      DO NZG=1,NZGMAX,NZGMAX-1
      DO NRG=1,NRGMAX
         R1=RG(NRG)
         Z1=ZG(NZG)
         CALL EQPSIB(R1,Z1,PSIBRZ)
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         FVB(N)=PSIBRZ(K)
      ENDDO
      ENDDO
C
      DO  NRG=1,NRGMAX,NRGMAX-1
      DO NZG=2,NZGMAX-1
         R1=RG(NRG)
         Z1=ZG(NZG)
         CALL EQPSIB(R1,Z1,PSIBRZ)
         K=1
         N=4*((NZG-1)*NRGMAX+NRG-1)+K
         FVB(N)=PSIBRZ(K)
      ENDDO
      ENDDO
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
      PSIV(1,5)= (64.D0/7.D0*R*Z**7
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
         PSIBRZ(1)=PSIBRZ(1)+PSIV(1,I)*PSIB(I)
         PSIBRZ(2)=PSIBRZ(2)+PSIV(2,I)*PSIB(I)
         PSIBRZ(3)=PSIBRZ(3)+PSIV(3,I)*PSIB(I)
         PSIBRZ(4)=PSIBRZ(4)+PSIV(4,I)*PSIB(I)
      ENDDO
C
      RETURN
      END
