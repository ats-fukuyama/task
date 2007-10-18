C     $Id$
C*******************************
C*** DIS FUNC INT TEST VOL.1 ***
C*******************************
C
      SUBROUTINE FPCALC_NL(NR,NS)
C
      INCLUDE 'fpcomm.inc'
C
      PARAMETER (N=NPM+2,M=NTHM+2,LNM=5)
      DIMENSION PLM(-1:LNM,M),PLG(-1:LNM,M),FPL(-1:LNM,N)
      DIMENSION D1PLM(-1:LNM,M)
      DIMENSION D1PLG(-1:LNM,M),D2PLG(-1:LNM,M)
      DIMENSION RM1M(-1:LNM,N),RM3M(-1:LNM,N)
      DIMENSION RM2M(-1:LNM,N),RM4M(-1:LNM,N)
      DIMENSION RM1G(-1:LNM,N),RM3G(-1:LNM,N)
      DIMENSION RM2G(-1:LNM,N),RM4G(-1:LNM,N)
C
      DIMENSION TX(NTHM+2),TY(NTHM+2),DF(NTHM+2)
      DIMENSION UTY(4,NTHM+2),UTY0(NTHM+2)
      DIMENSION TX1(NPM+2),TY1(NPM+2),DF1(NPM+2)
      DIMENSION UTY1(4,NPM+2),UTY10(NPM+2)
      DIMENSION PHYM(-1:LNM,N),PSYM(-1:LNM,N)
      DIMENSION D1PSYM(-1:LNM,N)
      DIMENSION PSYG(-1:LNM,N)
      DIMENSION D1PHYG(-1:LNM,N),D1PSYG(-1:LNM,N),D2PSYG(-1:LNM,N)
C
C     additional definition
      AMFD=PA(NS)*AMP
      AEFD=PZ(NS)*AEE
      PTFPL=PTFP(NR)
      PTFDL=PTFD(NR,NS)
      VTFPL=VTFP(NR)
      VTFDL=VTFD(NR,NS)
      RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)

      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
      RGAMH=RNUD(NR,NS)*SQRT(2.D0)*VTFD(NR,NS)*AMFP/(RNFP0*PTFP0*1.D20)
      IF(MODELR.eq.0)THEN
         TMC2FD0=0.D0
         TMC2FP0=0.D0
      ELSE IF(MODELR.eq.1)THEN
         TMC2FD0=(PTFD0/(AMFD*VC))**2
         TMC2FP0=(PTFP0/(AMFP*VC))**2
      END IF
C     end of additional definition 
      LLMIN=0

      DO NTH=1,NTHMAX
         CALL DPLEG(COSM(NTH),LLMAX,PLM(0,NTH),IER)
         PLM(-1,NTH)=0.D0
      END DO

      DO NTH=1,NTHMAX+1
         CALL DPLEG(COSG(NTH),LLMAX,PLG(0,NTH),IER)
         PLG(-1,NTH)=0.D0
      END DO

      DO L=LLMIN,LLMAX
         DO NTH=1,NTHMAX
            D1PLM(L,NTH)=L/SINM(NTH)*(COSM(NTH)*PLM(L,NTH)-PLM(L-1,NTH))
C     D2PLM(L,NTH)=-(L/(SINM(NTH)**2)+L**2)*PLM(L,NTH)
C     &              +L*COSM(NTH)/(SINM(NTH)**2)*PLM(L-1,NTH) 
         END DO
      END DO

      DO L=LLMIN,LLMAX
         NTH=1 
         D1PLG(L,NTH)=0.D0
         D2PLG(L,NTH)=0.D0
      END DO
      DO L=LLMIN,LLMAX
         NTH=NTHMAX+1 
         D1PLG(L,NTH)=0.D0
         D2PLG(L,NTH)=0.D0
      END DO

      DO L=LLMIN,LLMAX
         DO NTH=2,NTHMAX
            D1PLG(L,NTH)=L/SING(NTH)*(COSG(NTH)*PLG(L,NTH)-PLG(L-1,NTH))
            D2PLG(L,NTH)=-(L/(SING(NTH)**2)+L**2)*PLG(L,NTH)
     &           +L*COSG(NTH)/(SING(NTH)**2)*PLG(L-1,NTH) 
         END DO
      END DO

      DO NP=1,NPMAX
         DO L=LLMIN,LLMAX
            TX(1)=0.D0
            TY(1)=0.D0
            DO NTH=1,NTHMAX
               TX(NTH+1)=THM(NTH)
               TY(NTH+1)=FNS(NTH,NP,NR,NS)*PLM(L,NTH)*SINM(NTH)
            END DO
            TX(NTHMAX+2)=PI
            TY(NTHMAX+2)=0.D0
            DF(1)       = 0.D0
            DF(NTHMAX+2)= 0.D0
C     CALL SPLC(TX,NTHMAX+2,TY,DF,IOPT,TEMP,NTHM+2,IER)
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
     
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPL(L,NP)=0.5D0*(2*L+1.D0)*SUM1
         END DO
      END DO
C     \hat{M}_l calculation 
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(1-L))*RGAMB**(1+L)
         END DO
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM1M(L,NP)=PSUM-SUM2
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM1G(L,NPG)=PSUM-SUM3
         END DO
      END DO
C     End of \hat{M}_l calc
C     calc of N_l  
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(2+L))*RGAMB**(-L)
         END DO
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM2M(L,NP)=SUM4
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM2G(L,NPG)=SUM5
         END DO
      END DO
C     end of N_l calc  
C     calc of hat{M}_l^+ 
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(3-L))*RGAMB**(L-1)
         END DO
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            CALL SPL1DI(PCRIT,SUM6,TX1,UTY1,UTY10,NPMAX+2,IER)
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            RM3M(L,NP)=PSUM-SUM6
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM7,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM3G(L,NPG)=PSUM-SUM7
         END DO
      END DO
C     end of hat{M}_l^+
C     calc of hat{N}_l^+    
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(4+L))*RGAMB**(-L-2)
         END DO
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM8,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM4M(L,NP)=SUM8
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM9,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM4G(L,NPG)=SUM9
         END DO
      END DO
C     end of hat{N}_l^+
      DO L=LLMIN,LLMAX
         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            pabar=PTFP0/(AMFP*RGAMA)
            pbbar=PTFD0/AMFD
            PHYM(L,NP)=-1.D0/(2*L+1)*((PM(NP)**(-L-1))*RM2M(L,NP)
     &           *pabar**(-L-1)*pbbar**L
     &           +(PM(NP)**L*RM1M(L,NP))*pabar**L*pbbar**(-L-1) )

            PSYM(L,NP)=-0.5D0/(2*L+1)*(
     &              1.D0/(2*L+3)*((PM(NP)**(-L-1))*RM4M(L,NP)
     &           *pabar**(-L-1)*pbbar**(L+2)
     &                           +(PM(NP)**( L+2))*RM1M(L,NP)
     &           *pabar**(L+2)*pbbar**(-L-1))
     &             -1.D0/(2*L-1)*((PM(NP)**(-L+1))*RM2M(L,NP)
     &           *pabar**(-L+1)*pbbar**(L)
     &                              +(PM(NP)**L  )*RM3M(L,NP)
     &           *pabar**(L)*pbbar**(-L+1)))

            D1PSYM(L,NP)=-0.5D0/(2*L+1)*(
     &          1.D0/(2*L+3)*( (L+2)*(PM(NP)**(L+1))*RM1M(L,NP)
     &           *pabar**(L+1)*pbbar**(-L-1)
     &                       -(L+1)*(PM(NP)**(-L-2))*RM4M(L,NP)
     &           *pabar**(-L-2)*pbbar**(L+2))
     &         -1.D0/(2*L-1)*( L   *(PM(NP)**( L-1))*RM3M(L,NP)
     &           *pabar**(L-1)*pbbar**(-L+1)
     &                       -(L-1)*(PM(NP)**(-L  ))*RM2M(L,NP)
     &           *pabar**(-L)*pbbar**(L)))*pabar
     &           *RGAMA**(-2)
         END DO
      END DO

C
      DO 182 L=LLMIN,LLMAX
        NP=1
C        PHYG(L,NP)=0.D0
        PSYG(L,NP)=0.D0
        D1PHYG(L,NP)=0.D0
        D1PSYG(L,NP)=0.D0
        D2PSYG(L,NP)=0.D0
  182 CONTINUE
C
      DO L=LLMIN,LLMAX
         DO NP=2,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
            pabar=PTFP0/(AMFP*RGAMA)
            pbbar=PTFD0/AMFD
            PSYG(L,NP)=-0.5D0/(2*L+1)*(
     &           1.D0/(2*L+3)*((PG(NP)**(-L-1))*RM4G(L,NP)
     &           *pabar**(-L-1)*pbbar**(L+2)
     &           +(PG(NP)**( L+2))*RM1G(L,NP)
     &           *pabar**(L+2)*pbbar**(-L-1))
     &           -1.D0/(2*L-1)*( (PG(NP)**(-L+1))*RM2G(L,NP)
     &           *pabar**(-L+1)*pbbar**(L)
     &           +(PG(NP)**  L   )*RM3G(L,NP)
     &           *pabar**(L)*pbbar**(-L+1) ) )

            D1PHYG(L,NP)=-1.D0/(2*L+1)*(L   *(PG(NP)**(L-1))*RM1G(L,NP)
     &           *pabar**(L-1)*pbbar**(-L-1)
     &           -(L+1)*(PG(NP)**(-L-2))*RM2G(L,NP)
     &           *pabar**(-L-2)*pbbar**(L) )*pabar
     &          /RGAMA**2

            D1PSYG(L,NP)=-0.5D0/(2*L+1)*(
     &           1.D0/(2*L+3)*( (L+2)*(PG(NP)**( L+1))*RM1G(L,NP)
     &           *pabar**(L+1)*pbbar**(-L-1)
     &           -(L+1)*(PG(NP)**(-L-2))*RM4G(L,NP)
     &           *pabar**(-L-2)*pbbar**(L+2) )
     &           -1.D0/(2*L-1)*( L   *(PG(NP)**( L-1))*RM3G(L,NP)
     &           *pabar**(L-1)*pbbar**(-L+1)
     &           -(L-1)*(PG(NP)**(-L  ))*RM2G(L,NP)
     &           *pabar**(-L)*pbbar**(L) ) )*pabar
     &           *RGAMA**(-2)

            D2PSYG(L,NP)=-0.5D0/(2*L+1)*(
     &           DBLE(L+1)*(L+2)/(2*L+3)*( (PG(NP)**(-L-3))*RM4G(L,NP)
     &           *pabar**(-L-3)*pbbar**(L+2)
     &           +(PG(NP)**  L   )*RM1G(L,NP)
     &           *pabar**(L)*pbbar**(-L-1) )
     &           -DBLE(L  )*(L-1)/(2*L-1)*( (PG(NP)**(-L-1))*RM2G(L,NP)
     &           *pabar**(-L-1)*pbbar**(L)
     &           +(PG(NP)**( L-2))*RM3G(L,NP)
     &           *pabar**(L-2)*pbbar**(-L+1) ) )*(pabar)**2
     &           /RGAMA**4
         END DO
      END DO

C   
C*************************************************
C*****        KAKUSAN KEISU NO KEISAN        *****
C*************************************************
C
      FACT=-4.D0*PI*RGAMH*1.D20
      FACT2=-4.D0*PI*RGAMH*1.D20
      L0MIN=0
      L0MAX=1

      DO NP=1,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RGAMA2=SQRT(1.D0-(PG(NP)/RGAMA)**2*TMC2FP0)
         DO NTH=1,NTHMAX
            WA=0 
            DO L=L0MIN,L0MAX
               WA=WA+D2PSYG(L,NP)*PLM(L,NTH) 
            END DO
            DCPP(NTH,NP,NR)=FACT*WA*AMFP/PTFP0/(RGAMA2)**6
         END DO
      END DO

C      RNUL=RNUD(NR,NS)
      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RGAMA2=SQRT(1.D0-(PM(NP)/RGAMA)**2*TMC2FP0)
         DO NTH=1,NTHMAX+1
            WB=0 
            DO L=L0MIN,L0MAX
               WB=WB+( 1.D0/PM(NP)*D1PSYM(L,NP)*PLG(L,NTH)
     &              +1.D0/(PM(NP)**2)*PSYM(L,NP)*D2PLG(L,NTH) ) 
            END DO
            DCTT(NTH,NP,NR)=FACT*WB*AMFP/PTFP0/(RGAMA2)**6
C            DCTT(NTH,NP,NR)=FACT*WB 
C     &           +0.5D0*RNUL*ZEFF/PV
         END DO
      END DO
C
      DO NP=1,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RGAMA2=SQRT(1.D0-(PG(NP)/RGAMA)**2*TMC2FP0)
         DO NTH=1,NTHMAX
            WC=0
            WCTEST=0
            DO L=L0MIN,L0MAX
               WC=WC+D1PHYG(L,NP)*PLM(L,NTH)
            END DO
            FCPP(NTH,NP,NR)=FACT2*WC*PTFP0/AMFD/(RGAMA2)**3
         END DO
      END DO

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RGAMA2=SQRT(1.D0-(PM(NP)/RGAMA)**2*TMC2FP0)
         DO NTH=1,NTHMAX+1
            WD=0
            DO L=L0MIN,L0MAX
               WD=WD+1.D0/PM(NP)*PHYM(L,NP)*D1PLG(L,NTH)
            END DO
            FCTH(NTH,NP,NR)=FACT2*WD*PTFP0/AMFD/(RGAMA2)**3
         END DO
      END DO

      NP=1
      DO NTH=1,NTHMAX
         DCPT(NTH,NP,NR)=0.D0
      END DO

      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RGAMA2=SQRT(1.D0-(PG(NP)/RGAMA)**2*TMC2FP0)
         DO NTH=1,NTHMAX
            WE=0
            DO L=L0MIN,L0MAX
               WE=WE+( 1.D0/PG(NP)*D1PSYG(L,NP)*D1PLM(L,NTH)
     &              -1.D0/(PG(NP)**2)*PSYG(L,NP)*D1PLM(L,NTH) )
            END DO
            DCPT(NTH,NP,NR)=FACT*WE*AMFP/PTFP0/(RGAMA2)**6
         END DO
      END DO

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RGAMA2=SQRT(1.D0-(PM(NP)/RGAMA)**2*TMC2FP0)
         DO NTH=1,NTHMAX+1
            WF=0
            DO L=L0MIN,L0MAX
               WF=WF+( 1.D0/PM(NP)*D1PSYM(L,NP)*D1PLG(L,NTH)
     &              -1.D0/(PM(NP)**2)*PSYM(L,NP)*D1PLG(L,NTH) )
            END DO
            DCTP(NTH,NP,NR)=FACT*WF*AMFP/PTFP0/(RGAMA2)**6
         END DO
      END DO

      DO NP=2, NPMAX+1
         PNFP=PG(NP)
         DCPPLL=DCPP(1,NP,NR)
         FCPPLL=FCPP(1,NP,NR)
         RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
         vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
         ptatb=PG(NP)/RGAMA
         PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
         WRITE(6,'(I5,1P5E12.4)') NP,
     &      DCPPLL, FCPPLL,
     &      (FCPPLL/DCPPLL/PNFP*RGAMA),PCRIT
C     &      PNFP/RGAMA,
C     &      (D1PHYG(0,NP)+D1PHYG(0,NP))/AMFD/
C     &        (D2PSYG(0,NP)+D2PSYG(0,NP))/PNFP*RGAMA**(-2)
C     &        RM1G(0,NP),RM2G(0,NP),RM3G(0,NP),RM4G(0,NP),PCRIT
C     &        RM2G(0,NP)/RM4G(0,NP)*AMFD/PTFD0
      END DO

      IF (MODELA.EQ.1) THEN
C
         DO 1200 NTH=1,NTHMAX
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
         DO 1200 NP=1,NPMAX+1
            SUM1=0.D0
            SUM2=0.D0
            SUM3=0.D0
C
         DO 1100 NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            X=EPSR(NR)*COS(ETAL)*RR
            PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
            IF (COSM(NTH).GE.0.D0) THEN
               PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
            ELSE
               PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
            ENDIF
C
            SUM1=SUM1+DCPP(NTH,NP,NR)*COSM(NTH)/PCOS
            SUM2=SUM2+FCPP(NTH,NP,NR)*COSM(NTH)/PCOS
            SUM3=SUM3+DCPT(NTH,NP,NR)/SQRT(PSIB)
 1100 CONTINUE
         DCPP(NTH,NP,NR)=SUM1*DELH/PI
         FCPP(NTH,NP,NR)=SUM2*DELH/PI
         DCPT(NTH,NP,NR)=SUM3*DELH/PI
 1200 CONTINUE
C
         DO 1400 NTH=1,NTHMAX+1
            DELH=2.D0*ETAG(NTH,NR)/NAVMAX
         DO 1400 NP=1,NPMAX
            SUM4=0.D0
            SUM5=0.D0
            SUM6=0.D0
C
         DO 1300 NG=1,NAVMAX
            ETAL=DELH*(NG-0.5D0)
            X=EPSR(NR)*COS(ETAL)*RR
            PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
            IF(NTH.NE.NTHMAX/2) THEN
               ARG=1.D0-PSIB*SING(NTH)**2
               IF(ARG.GT.0.D0) THEN
                  IF (COSG(NTH).GE.0.D0) THEN
                     PCOS= SQRT(ARG)
                  ELSE
                     PCOS=-SQRT(ARG)
                  ENDIF
               ELSE
                  PCOS=0.D0
               ENDIF
               SUM4=SUM4+DCTT(NTH,NP,NR)*PCOS/(PSIB*COSG(NTH))
            ENDIF
            SUM5=SUM5+FCTH(NTH,NP,NR)/SQRT(PSIB)
            SUM6=SUM6+DCTP(NTH,NP,NR)/SQRT(PSIB)
 1300 CONTINUE
         DCTT(NTH,NP,NR)=SUM4*DELH/PI
         FCTH(NTH,NP,NR)=SUM5*DELH/PI
         DCTP(NTH,NP,NR)=SUM6*DELH/PI
 1400 CONTINUE      
C
      DO 2200 NP=1,NPMAX+1
         DO 2100 NTH=ITL(NR)+1,NTHMAX/2
            DCPP(NTH,NP,NR)=(DCPP(NTH,NP,NR)
     &                      +DCPP(NTHMAX-NTH+1,NP,NR))/2.D0
            FCPP(NTH,NP,NR)=(FCPP(NTH,NP,NR)
     &                      +FCPP(NTHMAX-NTH+1,NP,NR))/2.D0
            DCPT(NTH,NP,NR)=(DCPT(NTH,NP,NR)
     &                      +DCPT(NTHMAX-NTH+1,NP,NR))/2.D0
            DCPP(NTHMAX-NTH+1,NP,NR)=DCPP(NTH,NP,NR)
            FCPP(NTHMAX-NTH+1,NP,NR)=FCPP(NTH,NP,NR)
            DCPT(NTHMAX-NTH+1,NP,NR)=DCPT(NTH,NP,NR)
 2100    CONTINUE
         DCPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                    *( DCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +DCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +DCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +DCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         FCPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                    *( FCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +FCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +FCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +FCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DCPT(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                    *( DCPT(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +DCPT(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +DCPT(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +DCPT(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DCPP(ITU(NR),NP,NR)=DCPP(ITL(NR),NP,NR)
         FCPP(ITU(NR),NP,NR)=FCPP(ITL(NR),NP,NR)
         DCPT(ITU(NR),NP,NR)=DCPT(ITL(NR),NP,NR)
 2200 CONTINUE
      DO 2300 NP=1,NPMAX      
      DO 2300 NTH=ITL(NR)+1,NTHMAX/2
            DCTT(NTH,NP,NR)=(DCTT(NTH,NP,NR)
     &                      +DCTT(NTHMAX-NTH+2,NP,NR))/2.D0
            FCTH(NTH,NP,NR)=(FCTH(NTH,NP,NR)
     &                      +FCTH(NTHMAX-NTH+2,NP,NR))/2.D0
            DCTP(NTH,NP,NR)=(DCTP(NTH,NP,NR)
     &                      +DCTP(NTHMAX-NTH+2,NP,NR))/2.D0
            DCTT(NTHMAX-NTH+2,NP,NR)=DCTT(NTH,NP,NR)
            FCTH(NTHMAX-NTH+2,NP,NR)=FCTH(NTH,NP,NR)
            DCTP(NTHMAX-NTH+2,NP,NR)=DCTP(NTH,NP,NR)
 2300 CONTINUE
 3000 CONTINUE
      ENDIF   
      RETURN
      END
