C     $Id$
C
C ************************************************************
C
C      CALCULATION OF NONLINEAR COLLISIONAL OPERATOR
C
C ************************************************************
C
      SUBROUTINE FPCALC_NL(NR,NSB,NSA)
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
C     ----- definition of local quantities -----
C
      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA)
     &     /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
      vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
      pabbar=(PTFP0(NSA)*AMFD(NSB))/(PTFD0(NSB)*AMFP(NSA))

      IF(MODELR.eq.0)THEN
         TMC2FD0=0.D0
         TMC2FP0=0.D0
      ELSE IF(MODELR.eq.1)THEN
         TMC2FD0=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
         TMC2FP0=(PTFP0(NSA)/(AMFP(NSA)*VC))**2
      END IF
C
C     ----- calculation of Legendre Polynomials -----
C
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

      IF(MOD(IDBGFP,2).EQ.1) THEN
C
C     +++ plot of Legendre polynomials and their derivatives +++
C
         CALL PAGES
         CALL GRD1D(1,thm,plm,M,NTHMAX,LLMAX+2,'@PLM:@',0)
         CALL GRD1D(2,thm,d1plm,M,NTHMAX,LLMAX+2,'@D1PLM:@',0)
         CALL PAGEE
C
         CALL PAGES
         CALL GRD1D(1,thg,plg,M,NTHMAX+1,LLMAX+2,'@PLG:@',0)
         CALL GRD1D(2,thg,d1plg,M,NTHMAX+1,LLMAX+2,'@D1PLG:@',0)
         CALL GRD1D(3,thg,d2plg,M,NTHMAX+1,LLMAX+2,'@D2PLG:@',0)
         CALL PAGEE
      ENDIF
C
C     ----- Legendre expansion of distribution funstion FNS -----
C
      NS=NS_NSA(NSA)
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
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPL(NP,L)=0.5D0*(2*L+1.D0)*SUM1
         END DO
      END DO

C
C     ----- calculation of \hat{M}_l -----
C
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
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM1M(NP,L)=PSUM-SUM2
            ELSE
               RM1M(NP,L)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM1G(NPG,L)=PSUM-SUM3
            ELSE
               RM1G(NPG,L)=0.D0
            ENDIF
         END DO
      END DO

C
C     ----- calculation of \hat{N}_l -----
C
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
         CALL SPL1DI(PMAX,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) then
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM2M(NP,L)=SUM4
            ELSE
               RM2M(NP,L)=PSUM
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) then
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM2G(NPG,L)=SUM5
            ELSE
               RM2G(NPG,L)=PSUM
            ENDIF
         END DO
      END DO

C
C     ----- calculation of \hat{M}_l^+ -----
C
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
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) then
               CALL SPL1DI(PCRIT,SUM6,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM3M(NP,L)=PSUM-SUM6
            else
               RM3M(NP,L)=0.d0
            endif
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) then
               CALL SPL1DI(PCRIT,SUM7,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM3G(NPG,L)=PSUM-SUM7
            else
               RM3G(NPG,L)=0.d0
            endif
         END DO
      END DO

C
C     ----- calculation of \hat{N}_l^+ -----
C
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
         CALL SPL1DI(PMAX,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) then
               CALL SPL1DI(PCRIT,SUM8,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM4M(NP,L)=SUM8
            else
               RM4M(NP,L)=PSUM
            endif
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG)**2*TMC2FP0)
            ptatb=PG(NPG)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX) then
               CALL SPL1DI(PCRIT,SUM9,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM4G(NPG,L)=SUM9
            else
               RM4G(NPG,L)=PSUM
            endif
         END DO
      END DO

      IF(MOD(IDBGFP/2,2).EQ.1) THEN
C
C     +++ plot of Legendre expansion and M_l, N_l +++
C
         CALL PAGES
         CALL GRD1D(0,pm,fpl,N,NPMAX,LLMAX+2,'@FPL:@',0)
         CALL PAGEE
         CALL PAGES
         CALL GRD1D(1,pg,rm1g,N,NPMAX+1,LLMAX+2,'@RM1G:@',0)
         CALL GRD1D(2,pg,rm1m,N,NPMAX,  LLMAX+2,'@RM1M:@',0)
         CALL GRD1D(3,pg,rm3g,N,NPMAX+1,LLMAX+2,'@RM3G:@',0)
         CALL GRD1D(4,pg,rm3m,N,NPMAX,  LLMAX+2,'@RM3M:@',0)
         CALL PAGEE
         CALL PAGES
         CALL GRD1D(1,pg,rm2g,N,NPMAX+1,LLMAX+2,'@RM2G:@',0)
         CALL GRD1D(2,pg,rm2m,N,NPMAX,  LLMAX+2,'@RM2M:@',0)
         CALL GRD1D(3,pg,rm4g,N,NPMAX+1,LLMAX+2,'@RM4G:@',0)
         CALL GRD1D(4,pg,rm4m,N,NPMAX,  LLMAX+2,'@RM4M:@',0)
         CALL PAGEE
      ENDIF
C
C     ----- calculation of phi_l, psi_l and their derivatives -----
C
      DO L=LLMIN,LLMAX
         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            pabar=PTFP0(NSA)/(AMFP(NSA)*RGAMA)
            pbbar=PTFD0(NSB)/AMFD(NSB)
            pmbar=PM(NP)*pabar/pbbar
            PHYM(NP,L)=-1.D0/(2*L+1)
     &               *((pmbar**(-L-1))*RM2M(NP,L)
     &                +(pmbar**L      *RM1M(NP,L)))
     &               *pabbar

            PSYM(NP,L)=-0.5D0/(2*L+1)
     &               *(1.D0/(2*L+3)*((pmbar**(-L-1))*RM4M(NP,L)
     &                              +(pmbar**( L+2))*RM1M(NP,L))
     &                -1.D0/(2*L-1)*((pmbar**(-L+1))*RM2M(NP,L)
     &                              +(pmbar**L     )*RM3M(NP,L)))
     &               /pabbar

            D1PSYM(NP,L)=-0.5D0/(2*L+1)
     &               *(1.D0/(2*L+3)*( (L+2)*(pmbar**( L+1))*RM1M(NP,L)
     &                               -(L+1)*(pmbar**(-L-2))*RM4M(NP,L))
     &                -1.D0/(2*L-1)*(  L   *(pmbar**( L-1))*RM3M(NP,L)
     &                               -(L-1)*(pmbar**(-L  ))*RM2M(NP,L)))
     &               /RGAMA**3
         END DO
      END DO

C
      DO 182 L=LLMIN,LLMAX
        NP=1
        PSYG(NP,L)=0.D0
        D1PHYG(NP,L)=0.D0
        D1PSYG(NP,L)=0.D0
        D2PSYG(NP,L)=0.D0
  182 CONTINUE
C
      DO L=LLMIN,LLMAX
C         DO NP=2,NPMAX+1
         DO NP=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
            pabar=PTFP0(NSA)/(AMFP(NSA)*RGAMA)
            pbbar=PTFD0(NSB)/AMFD(NSB)
            IF(NP.EQ.1) THEN
               pgbar=0.001D0*PG(2)*pabar/pbbar
            ELSE
               pgbar=PG(NP)*pabar/pbbar
            ENDIF
            PSYG(NP,L)=-0.5D0/(2*L+1)
     &           *(1.D0/(2*L+3)*((pgbar**(-L-1))*RM4G(NP,L)
     &                          +(pgbar**( L+2))*RM1G(NP,L))
     &            -1.D0/(2*L-1)*((pgbar**(-L+1))*RM2G(NP,L)
     &                          +(pgbar**  L   )*RM3G(NP,L)))
     &           /pabbar

            D1PHYG(NP,L)=-1.D0/(2*L+1)
     &           *(  L   *(pgbar**( L-1))*RM1G(NP,L)
     &             -(L+1)*(pgbar**(-L-2))*RM2G(NP,L))
     &           *pabbar**2/RGAMA**3

            D1PSYG(NP,L)=-0.5D0/(2*L+1)
     &           *(1.D0/(2*L+3)*( (L+2)*(pgbar**( L+1))*RM1G(NP,L)
     &                           -(L+1)*(pgbar**(-L-2))*RM4G(NP,L))
     &            -1.D0/(2*L-1)*(  L   *(pgbar**( L-1))*RM3G(NP,L)
     &                           -(L-1)*(pgbar**(-L  ))*RM2G(NP,L)))

            D2PSYG(NP,L)=-0.5D0/(2*L+1)
     &           *(DBLE(L+1)*(L+2)/(2*L+3)*((pgbar**(-L-3))*RM4G(NP,L)
     &                                     +(pgbar**  L   )*RM1G(NP,L))
     &            -DBLE(L  )*(L-1)/(2*L-1)*((pgbar**(-L-1))*RM2G(NP,L)
     &                                     +(pgbar**( L-2))*RM3G(NP,L)))
     &           *pabbar/RGAMA**6
         END DO
      END DO

      IF(MOD(IDBGFP/4,2).EQ.1) THEN
C
C     +++ plot of Phi, Psi and their derivatives +++
C
         CALL PAGES
         CALL GRD1D(1,pm,psym,N,NPMAX,LLMAX+2,'@PSYM:@',0)
         CALL GRD1D(2,pg,psyg,N,NPMAX+1,LLMAX+2,'@PSYG:@',0)
         CALL GRD1D(3,pm,d1psym,N,NPMAX,LLMAX+2,'@D1PSYM:@',0)
         CALL GRD1D(4,pg,d1psyg,N,NPMAX+1,LLMAX+2,'@D1PSYG:@',0)
         CALL PAGEE
         CALL PAGES
         CALL GRD1D(1,pg,d2psyg,N,NPMAX+1,LLMAX+2,'@D2PSYG:@',0)
         CALL GRD1D(3,pm,phym,N,NPMAX,LLMAX+2,'@PHYM:@',0)
         CALL GRD1D(4,pg,d1phyg,N,NPMAX+1,LLMAX+2,'@D1PHYG:@',0)
         CALL PAGEE
      ENDIF
C
C     ----- calculation of local diffusion coefficienst -----
C   
      FACT=-4.D0*PI*RGAMH*1.D20

C      L0MIN=0

      DO NP=1,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         DO NTH=1,NTHMAX
            WA=0 
            WC=0
            DO L=LLMIN,LLMAX
               WA=WA+D2PSYG(NP,L)*PLM(NTH,L) 
               WC=WC+D1PHYG(NP,L)*PLM(NTH,L)
            END DO
            DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)
     &                      +FACT*WA*RGAMA**6
            FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)
     &                      +FACT*WC*RGAMA**3*AMFP(NSA)/AMFD(NSB)
         END DO
      END DO

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         DO NTH=1,NTHMAX+1
            WB=0 
            WD=0
            DO L=LLMIN,LLMAX
               WB=WB+( 1.D0/ PM(NP)    *D1PSYM(NP,L)*PLG(NTH,L)*RGAMA**4
     &          +  1.D0/(PM(NP)**2)*PSYM(NP,L)  *D2PLG(NTH,L) *RGAMA**2) 
               WD=WD+  1.D0/ PM(NP)    *PHYM(NP,L)  *D1PLG(NTH,L)
            END DO
            DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)
     &                     +FACT*WB
            FCTH2(NTH,NP,NR,NSB,NSA)=FCTH2(NTH,NP,NR,NSB,NSA)
     &                     +FACT*WD*RGAMA*AMFP(NSA)/AMFD(NSB)
         END DO
      END DO
C
      NP=1
      DO NTH=1,NTHMAX
         DCPT2(NTH,NP,NR,NSB,NSA)=0.D0
      END DO

      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         DO NTH=1,NTHMAX
            WE=0
            DO L=LLMIN,LLMAX
               WE=WE+( 1.D0/ PG(NP) *D1PSYG(NP,L)*D1PLM(NTH,L)*RGAMA**4
     &           -1.D0/(PG(NP)**2)*PSYG(NP,L)  *D1PLM(NTH,L) )*RGAMA**2
            END DO
            DCPT2(NTH,NP,NR,NSB,NSA)=DCPT2(NTH,NP,NR,NSB,NSA)
     &                     +FACT*WE
         END DO
      END DO

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         DO NTH=1,NTHMAX+1
            WF=0
            DO L=LLMIN,LLMAX
               WF=WF+( 1.D0/ PM(NP) *D1PSYM(NP,L)*D1PLG(NTH,L)*RGAMA**4
     &           -1.D0/(PM(NP)**2)*PSYM(NP,L)  *D1PLG(NTH,L) )*RGAMA**2
            END DO
            DCTP2(NTH,NP,NR,NSB,NSA)=DCTP2(NTH,NP,NR,NSB,NSA)
     &                     +FACT*WF
         END DO
      END DO
c
c      DO NP=2, NPMAX+1
c         PNFP=PG(NP)
c         RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
c         ptatb=PG(NP)/RGAMA
c         PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))*ptatb
c         ratio=AMFP/AMFD
c         if(ratio.eq.1.D0)WRITE(6,'(I5,1P5E14.5)') NP,
c     &      D1PHYG(NP,0)/D2PSYG(NP,0)/PNFP,
c     &      (FCPP(1,NP,NR)/DCPP(1,NP,NR)/PNFP*RGAMA),
c     &        PLM(1,0)
c      END DO
C     
      RETURN
      END
C
      SUBROUTINE FPCALC_NLAV(NR,NSA)
C
      INCLUDE 'fpcomm.inc'
      DIMENSION sum11(NSBM),sum12(NSBM),sum13(NSBM)
     &         ,sum14(NSBM),sum15(NSBM),sum16(NSBM)
C     
         DO NTH=1,NTHMAX
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            DO NP=1,NPMAX+1
               DO NSB=1,NSBMAX
                  sum11(NSB)=0.D0
                  sum12(NSB)=0.D0
                  sum13(NSB)=0.D0
               END DO
C     
               DO NG=1,NAVMAX
                  ETAL=DELH*(NG-0.5D0)
                  X=EPSR(NR)*COS(ETAL)*RR
                  PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
                  IF (COSM(NTH).GE.0.D0) THEN
                     PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
                  ELSE
                     PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
                  ENDIF
C
                  DO NSB=1,NSBMAX
                     sum11(NSB)=sum11(NSB)
     &                    +DCPP2(NTH,NP,NR,NSB,NSA)*COSM(NTH)/PCOS
                     sum12(NSB)=sum12(NSB)
     &                    +FCPP2(NTH,NP,NR,NSB,NSA)*COSM(NTH)/PCOS
                     sum13(NSB)=sum13(NSB)
     &                    +DCPT2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                  END DO
               END DO
               DO NSB=1,NSBMAX
                  DCPP2(NTH,NP,NR,NSB,NSA)=SUM11(NSB)*DELH/PI
                  FCPP2(NTH,NP,NR,NSB,NSA)=SUM12(NSB)*DELH/PI
                  DCPT2(NTH,NP,NR,NSB,NSA)=SUM13(NSB)*DELH/PI                  
               END DO
            END DO
         END DO
C         
         DO NTH=1,NTHMAX+1
            DELH=2.D0*ETAG(NTH,NR)/NAVMAX
            DO NP=1,NPMAX
               DO NSB=1,NSBMAX
                  sum14(NSB)=0.D0
                  sum15(NSB)=0.D0
                  sum16(NSB)=0.D0
               END DO
C     
               DO NG=1,NAVMAX
                  ETAL=DELH*(NG-0.5D0)
                  X=EPSR(NR)*COS(ETAL)*RR
                  PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
                  IF(NTH.NE.NTHMAX/2+1) THEN
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
                     DO NSB=1,NSBMAX
                        SUM14(NSB)=SUM14(NSB)
     &                       +DCTT2(NTH,NP,NR,NSB,NSA)*PCOS
     &                       /(PSIB*COSG(NTH))
                     END DO
                  ENDIF
                  DO NSB=1,NSBMAX
                     SUM15(NSB)=SUM15(NSB)
     &                    +FCTH2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                     SUM16(NSB)=SUM16(NSB)
     &                    +DCTP2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                  END DO

               END DO
               Do NSB=1,NSBMAX
                  DCTT2(NTH,NP,NR,NSB,NSA)=SUM14(NSB)*DELH/PI
                  FCTH2(NTH,NP,NR,NSB,NSA)=SUM15(NSB)*DELH/PI
                  DCTP2(NTH,NP,NR,NSB,NSA)=SUM16(NSB)*DELH/PI                  
               END DO
            END DO
         END DO
C
         DO NP=1,NPMAX+1
            DO NTH=ITL(NR)+1,NTHMAX/2
               DO NSB=1,NSBMAX
                  DCPP2(NTH,NP,NR,NSB,NSA)
     &                 =(DCPP2(NTH,NP,NR,NSB,NSA)
     &                  +DCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA))/2.D0
                  FCPP2(NTH,NP,NR,NSB,NSA)
     &                 =(FCPP2(NTH,NP,NR,NSB,NSA)
     &                  +FCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA))/2.D0
                  DCPT2(NTH,NP,NR,NSB,NSA)
     &                 =(DCPT2(NTH,NP,NR,NSB,NSA)
     &                  +DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA))/2.D0
                  DCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA)
     &                 =DCPP2(NTH,NP,NR,NSB,NSA)
                  FCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA)
     &                 =FCPP2(NTH,NP,NR,NSB,NSA)
                  DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA)
     &                 =DCPT2(NTH,NP,NR,NSB,NSA)
               END DO
            END DO
            DO NSB=1,NSBMAX
               DCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &           *( DCPP2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR)
     &             +DCPP2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR)
     &             +DCPP2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR)
     &             +DCPP2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               FCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &           *( FCPP2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR)
     &             +FCPP2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR)
     &             +FCPP2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR)
     &             +FCPP2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               DCPT2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &           *( DCPT2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR)
     &             +DCPT2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR)
     &             +DCPT2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR)
     &             +DCPT2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
               FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
               DCPT2(ITU(NR),NP,NR,NSB,NSA)=DCPT2(ITL(NR),NP,NR,NSB,NSA)
            END DO
         END DO

         DO NP=1,NPMAX
            DO NTH=ITL(NR)+1,NTHMAX/2
               DO NSB=1,NSBMAX
                  DCTT2(NTH,NP,NR,NSB,NSA)
     &                 =(DCTT2(NTH,NP,NR,NSB,NSA)
     &                  +DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA))/2.D0
                  FCTH2(NTH,NP,NR,NSB,NSA)
     &                 =(FCTH2(NTH,NP,NR,NSB,NSA)
     &                  +FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA))/2.D0
                  DCTP2(NTH,NP,NR,NSB,NSA)
     &                 =(DCTP2(NTH,NP,NR,NSB,NSA)
     &                  +DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA))/2.D0
                  DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA)
     &                 =DCTT2(NTH,NP,NR,NSB,NSA)
                  FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA)
     &                 =FCTH2(NTH,NP,NR,NSB,NSA)
                  DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA)
     &                 =DCTP2(NTH,NP,NR,NSB,NSA)
               END DO
            END DO
         END DO

      RETURN
      END
