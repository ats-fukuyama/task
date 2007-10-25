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
      DIMENSION PLM(M,-1:LNM),PLG(M,-1:LNM)
      DIMENSION D1PLM(M,-1:LNM)
      DIMENSION D1PLG(M,-1:LNM),D2PLG(M,-1:LNM)
      DIMENSION PLTEMP(0:LNM)
      DIMENSION FPL(N,-1:LNM)
      DIMENSION RM1M(N,-1:LNM),RM3M(N,-1:LNM)
      DIMENSION RM2M(N,-1:LNM),RM4M(N,-1:LNM)
      DIMENSION RM1G(N,-1:LNM),RM3G(N,-1:LNM)
      DIMENSION RM2G(N,-1:LNM),RM4G(N,-1:LNM)
C
      DIMENSION TX(NTHM+2),TY(NTHM+2),DF(NTHM+2)
      DIMENSION UTY(4,NTHM+2),UTY0(NTHM+2)
      DIMENSION TX1(NPM+2),TY1(NPM+2),DF1(NPM+2)
      DIMENSION UTY1(4,NPM+2),UTY10(NPM+2)
      DIMENSION PHYM(N,-1:LNM),PSYM(N,-1:LNM)
      DIMENSION D1PSYM(N,-1:LNM)
      DIMENSION PSYG(N,-1:LNM)
      DIMENSION D1PHYG(N,-1:LNM),D1PSYG(N,-1:LNM),D2PSYG(N,-1:LNM)
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
         CALL DPLEG(COSM(NTH),LLMAX,PLTEMP,IER)
         DO L=0,LLMAX
            PLM(NTH,L)=PLTEMP(L)
         ENDDO
         PLM(NTH,-1)=0.D0
      END DO

      DO NTH=1,NTHMAX+1
         CALL DPLEG(COSG(NTH),LLMAX,PLTEMP,IER)
         DO L=0,LLMAX
            PLG(NTH,L)=PLTEMP(L)
         ENDDO
         PLG(NTH,-1)=0.D0
      END DO

      DO L=LLMIN,LLMAX
         DO NTH=1,NTHMAX
            D1PLM(NTH,L)=L/SINM(NTH)*(COSM(NTH)*PLM(NTH,L)-PLM(NTH,L-1))
C     D2PLM(NTH,L)=-(L/(SINM(NTH)**2)+L**2)*PLM(NTH,L)
C     &              +L*COSM(NTH)/(SINM(NTH)**2)*PLM(NTH,L-1) 
         END DO
      END DO

      DO L=LLMIN,LLMAX
         NTH=1 
         D1PLG(NTH,L)=0.D0
         D2PLG(NTH,L)=-0.5D0*L*(L+1)
      END DO
      DO L=LLMIN,LLMAX
         NTH=NTHMAX+1 
         D1PLG(NTH,L)=0.D0
         D2PLG(NTH,L)=-0.5D0*L*(L+1)*(-1)**L
      END DO

      DO L=LLMIN,LLMAX
         DO NTH=2,NTHMAX
            D1PLG(NTH,L)=L/SING(NTH)*(COSG(NTH)*PLG(NTH,L)-PLG(NTH,L-1))
            D2PLG(NTH,L)=-(L/(SING(NTH)**2)+L**2)*PLG(NTH,L)
     &           +L*COSG(NTH)/(SING(NTH)**2)*PLG(NTH,L-1) 
         END DO
      END DO

c$$$!     plot of Legendre polynomials and their derivatives
c$$$      CALL PAGES
c$$$      CALL GRD1D(0,thm,plm,M,NTHMAX,LLMAX+2,'@PLM:@',0)
c$$$      CALL PAGEE
c$$$      CALL PAGES
c$$$      CALL GRD1D(0,thm,d1plm,M,NTHMAX,LLMAX+2,'@D1PLM:@',0)
c$$$      CALL PAGEE
c$$$C
c$$$      CALL PAGES
c$$$      CALL GRD1D(0,thg,plg,M,NTHMAX+1,LLMAX+2,'@PLG:@',0)
c$$$      CALL PAGEE
c$$$      CALL PAGES
c$$$      CALL GRD1D(0,thg,d1plg,M,NTHMAX+1,LLMAX+2,'@D1PLG:@',0)
c$$$      CALL PAGEE
c$$$      CALL PAGES
c$$$      CALL GRD1D(0,thg,d2plg,M,NTHMAX+1,LLMAX+2,'@D2PLG:@',0)
c$$$      CALL PAGEE

      DO L=LLMIN,LLMAX
         DO NP=1,NPMAX
            TX(1)=0.D0
            TY(1)=0.D0
            DO NTH=1,NTHMAX
               TX(NTH+1)=THM(NTH)
               TY(NTH+1)=FNS(NTH,NP,NR,NS)*PLM(NTH,L)*SINM(NTH)
            END DO
            TX(NTHMAX+2)=PI
            TY(NTHMAX+2)=0.D0
            DF(1)       = 0.D0
            DF(NTHMAX+2)= 0.D0
C     CALL SPLC(TX,NTHMAX+2,TY,DF,IOPT,TEMP,NTHM+2,IER)
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
     
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPL(NP,L)=0.5D0*(2*L+1.D0)*SUM1
         END DO
      END DO

      CALL PAGES
      CALL GRD1D(0,pm,fpl,N,NPMAX,LLMAX+2,'@FPL:@',0)
      CALL PAGEE

C     \hat{M}_l calculation 
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP)**(1-L))*RGAMB**(1+L)
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
            RM1M(NP,L)=PSUM-SUM2
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM1G(NPG,L)=PSUM-SUM3
         END DO
      END DO
C     End of \hat{M}_l calc

      CALL PAGES
      CALL GRD1D(0,pg,rm1g,N,NPMAX+1,LLMAX+2,'@RM1G:@',0)
      CALL PAGEE

C     calc of N_l  
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP)**(2+L))*RGAMB**(-L)
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
            RM2M(NP,L)=SUM4
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM2G(NPG,L)=SUM5
         END DO
      END DO
C     end of N_l calc  

      CALL PAGES
      CALL GRD1D(0,pg,rm2g,N,NPMAX+1,LLMAX+2,'@RM1G:@',0)
      CALL PAGEE

C     calc of hat{M}_l^+ 
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP)**(3-L))*RGAMB**(L-1)
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
            RM3M(NP,L)=PSUM-SUM6
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM7,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM3G(NPG,L)=PSUM-SUM7
         END DO
      END DO
C     end of hat{M}_l^+

      CALL PAGES
      CALL GRD1D(0,pg,rm3g,N,NPMAX+1,LLMAX+2,'@RM3G:@',0)
      CALL PAGEE

C     calc of hat{N}_l^+    
      DO L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP)**(4+L))*RGAMB**(-L-2)
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
            RM4M(NP,L)=SUM8
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.gt.PMAX) PCRIT=PMAX
            CALL SPL1DI(PCRIT,SUM9,TX1,UTY1,UTY10,NPMAX+2,IER)
            RM4G(NPG,L)=SUM9
         END DO
      END DO
C     end of hat{N}_l^+

      CALL PAGES
      CALL GRD1D(0,pg,rm4g,N,NPMAX+1,LLMAX+2,'@RM4G:@',0)
      CALL PAGEE

      DO L=LLMIN,LLMAX
         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            pabar=PTFP0/(AMFP*RGAMA)
            pbbar=PTFD0/AMFD
            pmbar=PM(NP)*pabar/pbbar
            PHYM(NP,L)=-1.D0/(2*L+1)
     &               *((pmbar**(-L-1))*RM2M(NP,L)
     &                +(pmbar**L      *RM1M(NP,L)))/pbbar

            PSYM(NP,L)=-0.5D0/(2*L+1)
     &               *(1.D0/(2*L+3)*((pmbar**(-L-1))*RM4M(NP,L)
     &                              +(pmbar**( L+2))*RM1M(NP,L))
     &                -1.D0/(2*L-1)*((pmbar**(-L+1))*RM2M(NP,L)
     &                              +(pmbar**L     )*RM3M(NP,L)))
     &               *pbbar

            D1PSYM(NP,L)=-0.5D0/(2*L+1)
     &               *(1.D0/(2*L+3)*( (L+2)*(pmbar**( L+1))*RM1M(NP,L)
     &                               -(L+1)*(pmbar**(-L-2))*RM4M(NP,L))
     &                -1.D0/(2*L-1)*(  L   *(pmbar**( L-1))*RM3M(NP,L)
     &                               -(L-1)*(pmbar**(-L  ))*RM2M(NP,L)))
     &               *pabar/RGAMA**2
         END DO
      END DO

      CALL PAGES
      CALL GRD1D(0,pm,phym,N,NPMAX,LLMAX+2,'@PHYM:@',0)
      CALL PAGEE
      CALL PAGES
      CALL GRD1D(0,pm,psym,N,NPMAX,LLMAX+2,'@PSYM:@',0)
      CALL PAGEE
      CALL PAGES
      CALL GRD1D(0,pm,d1psym,N,NPMAX,LLMAX+2,'@D1PSYM:@',0)
      CALL PAGEE
C
      DO 182 L=LLMIN,LLMAX
        NP=1
C        PHYG(NP,L)=0.D0
        PSYG(NP,L)=0.D0
        D1PHYG(NP,L)=0.D0
        D1PSYG(NP,L)=0.D0
        D2PSYG(NP,L)=0.D0
  182 CONTINUE
C
      DO L=LLMIN,LLMAX
         DO NP=2,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
            pabar=PTFP0/(AMFP*RGAMA)
            pbbar=PTFD0/AMFD
            pgbar=PG(NP)*pabar/pbbar
            PSYG(NP,L)=-0.5D0/(2*L+1)
     &           *(1.D0/(2*L+3)*((pgbar**(-L-1))*RM4G(NP,L)
     &                          +(pgbar**( L+2))*RM1G(NP,L))
     &            -1.D0/(2*L-1)*((pgbar**(-L+1))*RM2G(NP,L)
     &                          +(pgbar**  L   )*RM3G(NP,L)))
     &           *pbbar

            D1PHYG(NP,L)=-1.D0/(2*L+1)
     &           *(  L   *(pgbar**( L-1))*RM1G(NP,L)
     &             -(L+1)*(pgbar**(-L-2))*RM2G(NP,L))
     &           *pabar/(pbbar**2*RGAMA**2)

            D1PSYG(NP,L)=-0.5D0/(2*L+1)
     &           *(1.D0/(2*L+3)*( (L+2)*(pgbar**( L+1))*RM1G(NP,L)
     &                           -(L+1)*(pgbar**(-L-2))*RM4G(NP,L))
     &            -1.D0/(2*L-1)*(  L   *(pgbar**( L-1))*RM3G(NP,L)
     &                           -(L-1)*(pgbar**(-L  ))*RM2G(NP,L)))
     &           *pabar/RGAMA**2

            D2PSYG(NP,L)=-0.5D0/(2*L+1)
     &           *(DBLE(L+1)*(L+2)/(2*L+3)*((pgbar**(-L-3))*RM4G(NP,L)
     &                                     +(pgbar**  L   )*RM1G(NP,L))
     &            -DBLE(L  )*(L-1)/(2*L-1)*((pgbar**(-L-1))*RM2G(NP,L)
     &                                     +(pgbar**( L-2))*RM3G(NP,L)))
     &           *(pabar)**2/(pbbar*RGAMA**4)
         END DO
      END DO

      CALL PAGES
      CALL GRD1D(0,pg,psyg,N,NPMAX+1,LLMAX+2,'@PSYG:@',0)
      CALL PAGEE
      CALL PAGES
      CALL GRD1D(0,pg,d1phyg,N,NPMAX+1,LLMAX+2,'@D1PHYG:@',0)
      CALL PAGEE
      CALL PAGES
      CALL GRD1D(0,pg,d1psyg,N,NPMAX+1,LLMAX+2,'@D1PSYG:@',0)
      CALL PAGEE
      CALL PAGES
      CALL GRD1D(0,pg,d2psyg,N,NPMAX+1,LLMAX+2,'@D2PSYG:@',0)
      CALL PAGEE
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
               WA=WA+D2PSYG(NP,L)*PLM(NTH,L) 
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
               WB=WB+( 1.D0/ PM(NP)    *D1PSYM(NP,L)*PLG(NTH,L)
     &              +  1.D0/(PM(NP)**2)*PSYM(NP,L)*D2PLG(NTH,L) ) 
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
               WC=WC+D1PHYG(NP,L)*PLM(NTH,L)
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
               WD=WD+1.D0/PM(NP)*PHYM(NP,L)*D1PLG(NTH,L)
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
               WE=WE+( 1.D0/ PG(NP)    *D1PSYG(NP,L)*D1PLM(NTH,L)
     &                -1.D0/(PG(NP)**2)*PSYG(NP,L)  *D1PLM(NTH,L) )
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
               WF=WF+( 1.D0/ PM(NP)    *D1PSYM(NP,L)*D1PLG(NTH,L)
     &                -1.D0/(PM(NP)**2)*PSYM(NP,L)  *D1PLG(NTH,L) )
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
C     &      (D1PHYG(NP,0)+D1PHYG(NP,0))/AMFD/
C     &        (D2PSYG(NP,0)+D2PSYG(NP,0))/PNFP*RGAMA**(-2)
C     &        RM1G(NP,0),RM2G(NP,0),RM3G(NP,0),RM4G(NP,0),PCRIT
C     &        RM2G(NP,0)/RM4G(NP,0)*AMFD/PTFD0
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
