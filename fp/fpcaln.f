C     $Id$
C*******************************
C*** DIS FUNC INT TEST VOL.1 ***
C*******************************
C
      SUBROUTINE FPCALN
C
      INCLUDE 'fpcomm.inc'
C
      PARAMETER (N=NPM+2,M=NTHM+2,LNM=5)
      DIMENSION PLM(-1:LNM,M),PLG(-1:LNM,M),FPL(-1:LNM,N)
C      DIMENSION D2PLM(-1:LNM,M)
      DIMENSION D1PLM(-1:LNM,M)
      DIMENSION D1PLG(-1:LNM,M),D2PLG(-1:LNM,M)
      DIMENSION RM1M(-1:LNM,N),RM3M(-1:LNM,N)
      DIMENSION RM2M(-1:LNM,N),RM4M(-1:LNM,N)
      DIMENSION RM1G(-1:LNM,N),RM3G(-1:LNM,N)
      DIMENSION RM2G(-1:LNM,N),RM4G(-1:LNM,N)
C
      DIMENSION TX(NTHM+2),TY(NTHM+2),TEMP(NTHM+2,3)
      DIMENSION TX1(NPM+2),TY1(NPM+2),TEMP1(NPM+2,3)
      DIMENSION IOPT(2),DF(2)
      DIMENSION IOPT1(2),DF1(2)
      DIMENSION PHYM(-1:LNM,N),PSYM(-1:LNM,N)
      DIMENSION D1PSYM(-1:LNM,N)
C      DIMENSION D1PHYM(-1:LNM,N),D2PSYM(-1:LNM,N)
C      DIMENSION PHYG(-1:LNM,N)
      DIMENSION PSYG(-1:LNM,N)
      DIMENSION D1PHYG(-1:LNM,N),D1PSYG(-1:LNM,N),D2PSYG(-1:LNM,N)
C      DIMENSION SUM1(NTHM,NPMP,NRM),SUM2(NTHM,NPMP,NRM)
C      DIMENSION SUM3(NTHM,NPMP,NRM),SUM4(NTHMP,NPM,NRM)
C      DIMENSION SUM5(NTHM,NPMP,NRM),SUM6(NTHM,NPMP,NRM)
C      DIMENSION DELH(NTHMP,NRM),ETAL(NTHMP,NRM),PCOS(NTHMP,NRM)
C      DIMENSION PSI(NTHMP,NRM),Y(NTHMP,NRM),X(NTHMP,NRM)
C
      LLMIN=0
C      
      DO 60 NTH=1,NTHMAX
       CALL DPLEG(COSM(NTH),LLMAX,PLM(0,NTH),IER)
       PLM(-1,NTH)=0.D0
   60 CONTINUE
C
      DO 62 NTH=1,NTHMAX+1
       CALL DPLEG(COSG(NTH),LLMAX,PLG(0,NTH),IER)
       PLG(-1,NTH)=0.D0
   62 CONTINUE
C
      DO 65 L=LLMIN,LLMAX
      DO 65 NTH=1,NTHMAX
       D1PLM(L,NTH)=L/SINM(NTH)*(COSM(NTH)*PLM(L,NTH)-PLM(L-1,NTH))
C       D2PLM(L,NTH)=-(L/(SINM(NTH)**2)+L**2)*PLM(L,NTH)
C     &              +L*COSM(NTH)/(SINM(NTH)**2)*PLM(L-1,NTH) 
   65 CONTINUE 
C
      DO 66 L=LLMIN,LLMAX
       NTH=1 
       D1PLG(L,NTH)=0.D0
       D2PLG(L,NTH)=0.D0
   66 CONTINUE 
      DO 67 L=LLMIN,LLMAX
       NTH=NTHMAX+1 
       D1PLG(L,NTH)=0.D0
       D2PLG(L,NTH)=0.D0
   67 CONTINUE 

      DO 68 L=LLMIN,LLMAX
      DO 68 NTH=2,NTHMAX
       D1PLG(L,NTH)=L/SING(NTH)*(COSG(NTH)*PLG(L,NTH)-PLG(L-1,NTH))
       D2PLG(L,NTH)=-(L/(SING(NTH)**2)+L**2)*PLG(L,NTH)
     &              +L*COSG(NTH)/(SING(NTH)**2)*PLG(L-1,NTH) 
   68 CONTINUE 
C
      DO 1000 NR=1,NRMAX
C
      DO 90 NP=1,NPMAX
      DO 90 L=LLMIN,LLMAX
         TX(1)=0.D0
         TY(1)=0.D0
         DO 80 NTH=1,NTHMAX
            TX(NTH+1)=THM(NTH)
            TY(NTH+1)=F(NTH,NP,NR)*PLM(L,NTH)*SINM(NTH)
   80    CONTINUE
         TX(NTHMAX+2)=PI
         TY(NTHMAX+2)=0.D0
         IOPT(1) = 1
         DF(1)   = 0.D0
         IOPT(2) = 1
         DF(2)   = 0.D0
         CALL SPLC(TX,NTHMAX+2,TY,DF,IOPT,TEMP,NTHM+2,IER)
C
         CALL SPLQ(TX,NTHMAX+2,TY,TEMP,NTHM+2,0.D0,PI,SUM,IER)
         FPL(L,NP)=0.5D0*(2*L+1.D0)*SUM
   90 CONTINUE
C  
      DO 110 L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO 100 NNP=1,NPMAX
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(1-L))
  100    CONTINUE          
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         IOPT1(1) = 1
         DF1(1)   = 0.D0
         IOPT1(2) = 1
         DF1(2)   = 0.D0
         CALL SPLC(TX1,NPMAX+2,TY1,DF1,IOPT1,TEMP1,NPM+2,IER)
C
         DO 115 NP=1,NPMAX
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,PM(NP),PMAX,
     &                                      RM1M(L,NP),IER)
  115 CONTINUE      
         DO 116 NPG=1,NPMAX+1
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,PG(NPG),PMAX,
     &                                      RM1G(L,NPG),IER)
  116 CONTINUE      
  110 CONTINUE
C
      DO 130 L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO 120 NNP=1,NPMAX
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(2+L))
  120    CONTINUE          
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         IOPT1(1) = 1
         DF1(1)   = 0.D0
         IOPT1(2) = 1
         DF1(2)   = 0.D0
         CALL SPLC(TX1,NPMAX+2,TY1,DF1,IOPT1,TEMP1,NPM+2,IER)
C
         DO 135 NP=1,NPMAX
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,0.D0,PM(NP),
     &                                      RM2M(L,NP),IER)
  135 CONTINUE      
         DO 136 NPG=1,NPMAX+1
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,0.D0,PG(NPG),
     &                                      RM2G(L,NPG),IER)
  136 CONTINUE      
  130 CONTINUE
C
      DO 150 L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO 140 NNP=1,NPMAX
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(3-L))
  140    CONTINUE          
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         IOPT1(1) = 1
         DF1(1)   = 0.D0
         IOPT1(2) = 1
         DF1(2)   = 0.D0
         CALL SPLC(TX1,NPMAX+2,TY1,DF1,IOPT1,TEMP1,NPM+2,IER)
C
         DO 155 NP=1,NPMAX
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,PM(NP),PMAX,
     &                                      RM3M(L,NP),IER)
  155 CONTINUE      
         DO 156 NPG=1,NPMAX+1
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,PG(NPG),PMAX,
     &                                      RM3G(L,NPG),IER)
  156 CONTINUE      
  150 CONTINUE
C
      DO 170 L=LLMIN,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO 160 NNP=1,NPMAX
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(L,NNP)*(PM(NNP)**(4+L))
  160    CONTINUE          
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         IOPT1(1) = 1
         DF1(1)   = 0.D0
         IOPT1(2) = 1
         DF1(2)   = 0.D0
         CALL SPLC(TX1,NPMAX+2,TY1,DF1,IOPT1,TEMP1,NPM+2,IER)
C
         DO 175 NP=1,NPMAX
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,0.D0,PM(NP),
     &                                      RM4M(L,NP),IER)
  175 CONTINUE      
         DO 176 NPG=1,NPMAX+1
         CALL SPLQ(TX1,NPMAX+2,TY1,TEMP1,NPM+2,0.D0,PG(NPG),
     &                                      RM4G(L,NPG),IER)
  176 CONTINUE      
  170 CONTINUE
C
      DO 180 L=LLMIN,LLMAX
      DO 180 NP=1,NPMAX
        PHYM(L,NP)=-1.D0/(2*L+1)*((PM(NP)**(-L-1))*RM2M(L,NP)
     &                           +(PM(NP)**L     )*RM1M(L,NP))
C    
        PSYM(L,NP)=-0.5D0/(2*L+1)*(
     &              1.D0/(2*L+3)*((PM(NP)**(-L-1))*RM4M(L,NP)
     &                           +(PM(NP)**( L+2))*RM1M(L,NP))
     &             -1.D0/(2*L-1)*((PM(NP)**(-L+1))*RM2M(L,NP)  
     &                              +(PM(NP)**L  )*RM3M(L,NP)) )
C     
C        D1PHYM(L,NP)=-1.D0/(2*L+1)*( L   *(PM(NP)**( L-1))*RM1M(L,NP)
C     &                             -(L+1)*(PM(NP)**(-L-2))*RM2M(L,NP))
C       
        D1PSYM(L,NP)=-0.5D0/(2*L+1)*(
     &          1.D0/(2*L+3)*((L+2)*(PM(NP)**( L+1))*RM1M(L,NP)
     &                       -(L+1)*(PM(NP)**(-L-2))*RM4M(L,NP))
     &         -1.D0/(2*L-1)*( L   *(PM(NP)**( L-1))*RM3M(L,NP)  
     &                       -(L-1)*(PM(NP)**(-L  ))*RM2M(L,NP))) 
C       
C        D2PSYM(L,NP)=-0.5D0/(2*L+1)*(
C     &     DBLE(L+1)*(L+2)/(2*L+3)*((PM(NP)**(-L-3))*RM4M(L,NP)
C     &                             +(PM(NP)** L    )*RM1M(L,NP))
C     &    -DBLE(L  )*(L-1)/(2*L-1)*((PM(NP)**(-L-1))*RM2M(L,NP)
C     &                             +(PM(NP)**( L-2))*RM3M(L,NP))) 
  180 CONTINUE
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
      DO 184 L=LLMIN,LLMAX
      DO 184 NP=2,NPMAX+1
C        PHYG(L,NP)=-1.D0/(2*L+1)*((PG(NP)**(-L-1))*RM2G(L,NP)
C     &                           +(PG(NP)** L    )*RM1G(L,NP))
C    
        PSYG(L,NP)=-0.5D0/(2*L+1)*(
     &              1.D0/(2*L+3)*((PG(NP)**(-L-1))*RM4G(L,NP)
     &                           +(PG(NP)**( L+2))*RM1G(L,NP))
     &             -1.D0/(2*L-1)*((PG(NP)**(-L+1))*RM2G(L,NP)  
     &                           +(PG(NP)**  L   )*RM3G(L,NP)) )
C     
        D1PHYG(L,NP)=-1.D0/(2*L+1)*(L   *(PG(NP)**( L-1))*RM1G(L,NP)
     &                            -(L+1)*(PG(NP)**(-L-2))*RM2G(L,NP))
C       
        D1PSYG(L,NP)=-0.5D0/(2*L+1)*(
     &          1.D0/(2*L+3)*((L+2)*(PG(NP)**( L+1))*RM1G(L,NP)
     &                       -(L+1)*(PG(NP)**(-L-2))*RM4G(L,NP))
     &         -1.D0/(2*L-1)*( L   *(PG(NP)**( L-1))*RM3G(L,NP)  
     &                       -(L-1)*(PG(NP)**(-L  ))*RM2G(L,NP))) 
C       
        D2PSYG(L,NP)=-0.5D0/(2*L+1)*(
     &      DBLE(L+1)*(L+2)/(2*L+3)*((PG(NP)**(-L-3))*RM4G(L,NP)
     &                              +(PG(NP)**  L   )*RM1G(L,NP))
     &     -DBLE(L  )*(L-1)/(2*L-1)*((PG(NP)**(-L-1))*RM2G(L,NP)
     &                              +(PG(NP)**( L-2))*RM3G(L,NP))) 
C     IF(L.EQ.0) WRITE(6,*) NP,D2PSYG(L,NP),RM1G(L,NP),RM4G(L,NP)
  184 CONTINUE
C   
C*************************************************
C*****        KAKUSAN KEISU NO KEISAN        *****
C*************************************************
C
      FACT=-4.D0*PI*RNU0 
      AMASS=1.D0
      BMASS=1.D0
      RMASS=AMASS/BMASS
      FACT2=FACT*RMASS
      L0MIN=0
      L0MAX=1
C
      DO 200 NP=1,NPMAX+1
      DO 200 NTH=1,NTHMAX
         WA=0 
         DO 190 L=L0MIN,L0MAX
            WA=WA+D2PSYG(L,NP)*PLM(L,NTH) 
  190    CONTINUE       
         DCPP(NTH,NP,NR)=FACT*WA
  200 CONTINUE               
C
         RNUL=RNU0*RNFP(NR)
      DO 220 NP=1,NPMAX
         PX=PM(NP)
         PV=PX/SQRT(1.D0+THETA0*PX*PX)
      DO 220 NTH=1,NTHMAX+1
         WB=0 
         DO 210 L=L0MIN,L0MAX
            WB=WB+( 1.D0/PM(NP)*D1PSYM(L,NP)*PLG(L,NTH)
     &             +1.D0/(PM(NP)**2)*PSYM(L,NP)*D2PLG(L,NTH) ) 
  210    CONTINUE       
         DCTT(NTH,NP,NR)=FACT*WB 
     &                           +0.5D0*RNUL*ZEFF/PV
  220 CONTINUE               
C
      DO 240 NP=1,NPMAX+1
      DO 240 NTH=1,NTHMAX
         WC=0 
         DO 230 L=L0MIN,L0MAX
            WC=WC+D1PHYG(L,NP)*PLM(L,NTH) 
  230    CONTINUE       
         FCPP(NTH,NP,NR)=FACT2*WC
  240 CONTINUE               
C
      DO 260 NP=1,NPMAX
      DO 260 NTH=1,NTHMAX+1
         WD=0 
         DO 250 L=L0MIN,L0MAX
            WD=WD+1.D0/PM(NP)*PHYM(L,NP)*D1PLG(L,NTH) 
  250    CONTINUE       
         FCTH(NTH,NP,NR)=FACT2*WD
  260 CONTINUE               
C
      NP=1
      DO 270 NTH=1,NTHMAX
         DCPT(NTH,NP,NR)=0.D0
  270 CONTINUE               
      DO 280 NP=2,NPMAX+1
      DO 280 NTH=1,NTHMAX
         WE=0
         DO 275 L=L0MIN,L0MAX
            WE=WE+( 1.D0/PG(NP)*D1PSYG(L,NP)*D1PLM(L,NTH)
     &             -1.D0/(PG(NP)**2)*PSYG(L,NP)*D1PLM(L,NTH) ) 
  275    CONTINUE       
         DCPT(NTH,NP,NR)=FACT*WE
  280 CONTINUE               
C
      DO 300 NP=1,NPMAX
      DO 300 NTH=1,NTHMAX+1
         WF=0
         DO 290 L=L0MIN,L0MAX
            WF=WF+( 1.D0/PM(NP)*D1PSYM(L,NP)*D1PLG(L,NTH)
     &             -1.D0/(PM(NP)**2)*PSYM(L,NP)*D1PLG(L,NTH) ) 
  290    CONTINUE
         DCTP(NTH,NP,NR)=FACT*WF
  300 CONTINUE    
C
 1000 CONTINUE
C
      IF (MODELA.EQ.1) THEN
C
      DO 2000 NR=1,NRMAX
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
 2000 CONTINUE
C
      DO 3000 NR=1,NRMAX
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
