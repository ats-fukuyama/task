C
C     *** NUMERICAL vs. ANALYTIC SOLUTIONS ***
C
C     HISTORY
C
C     2002/01/15
C     2002/02/08
C     2002/02/11
C     2002/02/12
C     2002/02/18
C     2002/02/19
C     2002/03/01
C     2002/03/05
C     2002/03/06
C     2002/03/09
C     2002/03/13
C     2002/03/19
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
      PARAMETER (NXM=64*128,NZM=64*128)
      DIMENSION GX(NXM,2),GY(NXM,2),GJY(NXM,2)
      DIMENSION RX(NXM),ZX(NZM),PSI(NXM,NZM),PSIN(NXM),HJ(NXM)
      DIMENSION PSIZAX(NXM),RNUM(NXM),HJT(NXM)
      CHARACTER TXT*80
C
      CALL GSOPEN
C
C     MUL : Only a multiple of 64
C
      MUL   = 64
      NXMAX = 64*MUL
      NZMAX = 64*MUL
C
C     *********************************************
C     *            Given Values                   *
C     *********************************************
C
C     *** CONSTANTS ***
C
C        PI    : Pi
C        RMU0  : Permeability of free space
C        AMP   : Proton mass
C        AEE   : Electron charge
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      AMP    = 1.6726231D-27
      AEE    = 1.60217733D-19
C
C     *** CONFIGURATION PARAMETERS ***
C
C        RR    : Plasma major radius                             (m)
C        RA    : Plasma minor radius                             (m)
C        ZE    : Plasma range of z axis                          (m)
C        EPSA  : Cross-section parameter ; -1 < EPSA < 1
C        VPHI  : Troidal rotation velocity                     (m/s)
C        TMP   : Plasma temperature                            (keV)
C        PINT  : Plasma pressure                               (MPa)
C
      RR     = 3.D0
      RA     = 1.D0
      ZE     = RA
      EPSA   = 0.D0
      VPHI   = 5.D4
      TMP    = 6.D0*1.D3
      PINT   = 1.D-1*1.D6
      RC     = RR-RA
C
C     *** CONSTANTS OBTAINED BY ANALYTIC CALCULATION ***
C
      RAXIS  = 3.16275D0
      ZAXIS  = 0.D0
C
C     *** CHANGE PARAMETER ***
C
      NCHG   = 1
      NLRC   = 0
C
C     *** ARBITRARY CONSTANT ***
C
C        HJ0   : Current density at R=RC                   (MA/m**2)
C
      HJ0    = 1.335477D-1*1.D6
C
C     *********************************************
C     *            Calculate constant             *
C     *********************************************
C
      OTC = (VPHI**2*AMP)/(RR**2*(TMP*AEE))
C
C     *********************************************
C     *            ANALYTIC PSI,J                 *
C     *********************************************
C
      CC1=((EPSA-1.D0)*RAXIS**2/(8.D0*RR**2))+(1.D0/(2.D0*RR**2*OTC))
     &     *(EXP(OTC*RAXIS**2/2.D0)-1.D0)
      A1=(CC1/RR**2)*(RAXIS**2-(RR-RA)**2)
      A1R=(CC1/RR**2)*(RAXIS**2-(RR+RA)**2)
      A2=-((EPSA-1.D0)*RAXIS**4)/(1.6D1*RR**4)+(1/(RR**4*OTC**2))
     &     *(1.D0+RAXIS**2*OTC/2.D0-EXP(RAXIS**2*OTC/2.D0))
      A3=((EPSA-1.D0)*(RR-RA)**4)/(1.6D1*RR**4)-(1/(RR**4*OTC**2))
     &     *(1.D0+(RR-RA)**2*OTC/2.D0-EXP((RR-RA)**2*OTC/2.D0))
      A3R=((EPSA-1.D0)*(RR+RA)**4)/(1.6D1*RR**4)-(1/(RR**4*OTC**2))
     &     *(1.D0+(RR+RA)**2*OTC/2.D0-EXP((RR+RA)**2*OTC/2.D0))
      IF (A1+A2+A3.LT.0.D0) THEN
         AA=-(A1+A2+A3)
      ELSEIF (A1+A2+A3.EQ.0.D0) THEN
         WRITE(6,*) "ERROR-LEFT"
         STOP
      ELSEIF (A1+A2+A3.GT.0.D0) THEN
         AA=A1+A2+A3
      ENDIF
C
      IF (A1R+A2+A3R.LT.0.D0) THEN
         AAR=-(A1R+A2+A3R)
      ELSEIF (A1R+A2+A3R.EQ.0.D0) THEN
         WRITE(6,*) "ERROR-RIGHT"
         STOP
      ELSEIF (A1R+A2+A3R.GT.0.D0) THEN
         AAR=A1R+A2+A3R
      ENDIF
      PP   = SQRT(RMU0*RR**4*PINT/AA)
      PPR  = SQRT(RMU0*RR**4*PINT/AAR)
      HMT = 2.D0*RMU0*PINT*RC**2
     &     *((RR**4/PP*RC)*RMU0*HJ0-EXP(RC**2*OTC/2.D0))
      HMTR = 2.D0*RMU0*PINT*RC**2
     &     *((RR**4/PPR*RC)*RMU0*HJ0-EXP(RC**2*OTC/2.D0))
      HM   = (HMT*PP)/(2.D0*RMU0*PINT*RR**2)
      HMR  = (HMTR*PPR)/(2.D0*RMU0*PINT*RR**2)
      CC2  = -CC1*PP*(RR-RA)**2/RR**2+PP*A3
      CC2R = -CC1*PPR*(RR+RA)**2/RR**2+PPR*A3R
C
C     *** CALCULATE RKAP ***
C
C        RKAP  : Plasma shape elongation
C
      HNUM = EXP(OTC*RAXIS**2/2.D0)-(1.D0-EPSA)/2.D0
      IF (NLRC.EQ.0) THEN
         DEN = ((RR**2*HM)/(RAXIS**2*PP))+(1.D0-EPSA)/2.D0
      ELSE
         DEN = ((RR**2*HMR)/(RAXIS**2*PPR))+(1.D0-EPSA)/2.D0
      ENDIF
      RKAPRA = SQRT(HNUM/DEN)
C
C     ***
C
      ZDLTMX=0.D0
      RRMIN=3.D0
      RRMAX=3.D0
      ZZMIN=0.D0
      ZZMAX=0.D0
      PSIM=1.D0
      PSIMAX=0.D0
      DX = (2.D0*RA)/(NXMAX-1)
      DZ = (2.D0*ZE)/(NZMAX-1)
C
      DO NX=1,NXMAX
      DO NZ=1,NZMAX
C
         IF (NX.LE.NXMAX/2) THEN
            RX(NX)=(RR-RA)+DX*DBLE(NX-1)
         ELSE
            RX(NX)=(RR+DBLE(DX/2))+DX*DBLE(NX-NXMAX/2-1)
         ENDIF
         IF (NX.LE.NZMAX/2) THEN
            ZX(NZ)=-DZ*(DBLE(NZMAX/2-1)+0.5D0)+DZ*DBLE(NZ-1)
         ELSE
            ZX(NZ)=DBLE(DZ/2)+DZ*DBLE(NZ-NZMAX/2-1)
         ENDIF
C
         IF (NLRC.EQ.0) THEN
         PSI(NX,NZ)=CC1*PP*RX(NX)**2/RR**2
     &        -(HM*ZX(NZ)**2/(2.D0*RR**2))
     &        +PP*((EPSA-1.D0)/4.D0)*(ZX(NZ)**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PP/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2
         ELSE
         PSI(NX,NZ)=CC1*PPR*RX(NX)**2/RR**2
     &        -(HM*ZX(NZ)**2/(2.D0*RR**2))
     &        +PPR*((EPSA-1.D0)/4.D0)*(ZX(NZ)**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PPR/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2R
         ENDIF
         IF (PSI(NX,NZ).GT.0.D0) THEN
            IF (PSI(NX,NZ).GT.PSIMAX) PSIMAX=PSI(NX,NZ)
            IF (ZX(NZ).GT.ZDLTMX) THEN
               ZDLTMX=ZX(NZ)
               RDLTMX=RX(NX)
            ENDIF
            IF (RX(NX).LT.RRMIN) RRMIN=RX(NX)
            IF (RX(NX).GT.RRMAX) RRMAX=RX(NX)
            IF (ZX(NZ).LT.ZZMIN) ZZMIN=ZX(NZ)
            IF (ZX(NZ).GT.ZZMAX) ZZMAX=ZX(NZ)
         ENDIF
         IF (ABS(PSI(NX,NZ)).LT.PSIM) THEN
            PSIM=ABS(PSI(NX,NZ))
            PSIMIN=PSI(NX,NZ)
         ENDIF
C
      ENDDO
C      write(6,*) rx(nx)
C
         IF (NLRC.EQ.0) THEN
         PSIZAX(NX)=CC1*PP*RX(NX)**2/RR**2
     &        -(HM*ZAXIS**2/(2.D0*RR**2))
     &        +PP*((EPSA-1.D0)/4.D0)*(ZAXIS**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PP/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2
         HJ(NX)=(HM/(RMU0*RR**2*RX(NX))+(RX(NX)*PP/(RMU0*RR**4))
     &           *EXP(RX(NX)**2*OTC/2.D0))*1.D-6
         ELSE
         PSIZAX(NX)=CC1*PPR*RX(NX)**2/RR**2
     &        -(HMR*ZAXIS**2/(2.D0*RR**2))
     &        +PPR*((EPSA-1.D0)/4.D0)*(ZAXIS**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*(RX(NX)**2/RR**2)
     &        +(PPR/(RR**4*OTC**2))*(1.D0+(OTC*RX(NX)**2/2.D0)
     &        -EXP(OTC*RX(NX)**2/2.D0))+CC2R
         HJ(NX)=(HMR/(RMU0*RR**2*RX(NX))+(RX(NX)*PPR/(RMU0*RR**4))
     &           *EXP(RX(NX)**2*OTC/2.D0))*1.D-6
         ENDIF
         IF (PSIZAX(NX).GT.PSIMAX) PSIMAX=PSIZAX(NX)
C
      ENDDO
C      write(6,*) rrmin,rrmax,zzmin,zzmax,psimin,psimax
C
C     *********************************************
C     *          Calculate RKAP,RDLT              *
C     *********************************************
C
      A=ABS(RRMAX-RRMIN)
      B=ABS(ZZMAX-ZZMIN)
      RKAP=B/A
      C=0.5D0*(RRMAX-RRMIN)
      D=ABS((RRMAX-C)-RDLTMX)
      RDLT=D/C
C
C     *********************************************
C     *            NUMERICAL PSI                  *
C     *********************************************
C
      PSIN(1)  =  -1.1196792401011622D-2
      PSIN(2)  =  -3.3449933606798572D-2
      PSIN(3)  =  -5.5540849469167464D-2
      PSIN(4)  =  -7.7447405882131033D-2
      PSIN(5)  =  -9.9147117624733099D-2
      PSIN(6)  =  -0.1206171497340683D0
      PSIN(7)  =  -0.1418343190170958D0
      PSIN(8)  =  -0.1627750957135530D0
      PSIN(9)  =  -0.1834156053255615D0
      PSIN(10) =  -0.2037316306335839D0
      PSIN(11) =  -0.2236986139234690D0
      PSIN(12) =  -0.2432916594557582D0
      PSIN(13) =  -0.2624855362166691D0
      PSIN(14) =  -0.2812546810008735D0
      PSIN(15) =  -0.2995732018902857D0
      PSIN(16) =  -0.3174148822119035D0
      PSIN(17) =  -0.3347531850833447D0
      PSIN(18) =  -0.3515612586901082D0
      PSIN(19) =  -0.3678119424885558D0
      PSIN(20) =  -0.3834777746007136D0
      PSIN(21) =  -0.3985310007737615D0
      PSIN(22) =  -0.4129435854398460D0
      PSIN(23) =  -0.4266872256684879D0
      PSIN(24) =  -0.4397333692248303D0
      PSIN(25) =  -0.4520532386707346D0
      PSIN(26) =  -0.4636178647630797D0
      PSIN(27) =  -0.4743981349771284D0
      PSIN(28) =  -0.4843648684874897D0
      PSIN(29) =  -0.4934889422355945D0
      PSIN(30) =  -0.5017415309619615D0
      PSIN(31) =  -0.5090946700659170D0
      PSIN(32) =  -0.5155233473299582D0
      PSIN(33) =  -0.5210095125263042D0
      PSIN(34) =  -0.5254860209572469D0
      PSIN(35) =  -0.5289464711556562D0
      PSIN(36) =  -0.5313552798785677D0
      PSIN(37) =  -0.5326780154369875D0
      PSIN(38) =  -0.5328800691538271D0
      PSIN(39) =  -0.5319265334337927D0
      PSIN(40) =  -0.5297821628760251D0
      PSIN(41) =  -0.5264113518940063D0
      PSIN(42) =  -0.5217781187802599D0
      PSIN(43) =  -0.5158460932801249D0
      PSIN(44) =  -0.5085785063488917D0
      PSIN(45) =  -0.4999381813447626D0
      PSIN(46) =  -0.4898875261894374D0
      PSIN(47) =  -0.4783885261860112D0
      PSIN(48) =  -0.4654027372806923D0
      PSIN(49) =  -0.4508912796171261D0
      PSIN(50) =  -0.4348148312734998D0
      PSIN(51) =  -0.4171336221007859D0
      PSIN(52) =  -0.3978074276000611D0
      PSIN(53) =  -0.3767955627906666D0
      PSIN(54) =  -0.3540568760308886D0
      PSIN(55) =  -0.3295497427600113D0
      PSIN(56) =  -0.3032320591358474D0
      PSIN(57) =  -0.2750612355457167D0
      PSIN(58) =  -0.2449941899716957D0
      PSIN(59) =  -0.2129873411930705D0
      PSIN(60) =  -0.1789966018104755D0
      PSIN(61) =  -0.1429773710773267D0
      PSIN(62) =  -0.1048845275249571D0
      PSIN(63) =  -6.4672421368418131D-2
      PSIN(64) =  -2.2294866680252694D-2
C
      HJT(1)  =   0.5364402205789185D0
      HJT(2)  =   0.5399095599397759D0
      HJT(3)  =   0.5434503093957147D0
      HJT(4)  =   0.5470593529625080D0
      HJT(5)  =   0.5507337556295782D0
      HJT(6)  =   0.5544707504432505D0
      HJT(7)  =   0.5582677266807472D0
      HJT(8)  =   0.5621222190089521D0
      HJT(9)  =   0.5660318975335898D0
      HJT(10) =   0.5699945586546486D0
      HJT(11) =   0.5740081166528505D0
      HJT(12) =   0.5780705959398849D0
      HJT(13) =   0.5821801239121031D0
      HJT(14) =   0.5863349243535565D0
      HJT(15) =   0.5905333113397294D0
      HJT(16) =   0.5947736835981956D0
      HJT(17) =   0.5990545192867240D0
      HJT(18) =   0.6033743711532322D0
      HJT(19) =   0.6077318620453984D0
      HJT(20) =   0.6121256807408116D0
      HJT(21) =   0.6165545780712844D0
      HJT(22) =   0.6210173633173936D0
      HJT(23) =   0.6255129008515238D0
      HJT(24) =   0.6300401070096574D0
      HJT(25) =   0.6345979471739284D0
      HJT(26) =   0.6391854330495614D0
      HJT(27) =   0.6438016201212476D0
      HJT(28) =   0.6484456052753148D0
      HJT(29) =   0.6531165245752199D0
      HJT(30) =   0.6578135511789500D0
      HJT(31) =   0.6625358933878878D0
      HJT(32) =   0.6672827928175509D0
      HJT(33) =   0.6720485506057368D0
      HJT(34) =   0.6768323989277270D0
      HJT(35) =   0.6816386212535697D0
      HJT(36) =   0.6864665801909384D0
      HJT(37) =   0.6913156638531058D0
      HJT(38) =   0.6961852846118922D0
      HJT(39) =   0.7010748779230777D0
      HJT(40) =   0.7059839012194182D0
      HJT(41) =   0.7109118328667627D0
      HJT(42) =   0.7158581711791204D0
      HJT(43) =   0.7208224334888248D0
      HJT(44) =   0.7258041552682392D0
      HJT(45) =   0.7308028892996928D0
      HJT(46) =   0.7358182048905916D0
      HJT(47) =   0.7408496871308603D0
      HJT(48) =   0.7458969361900741D0
      HJT(49) =   0.7509595666518276D0
      HJT(50) =   0.7560372068830610D0
      HJT(51) =   0.7611294984362190D0
      HJT(52) =   0.7662360954822658D0
      HJT(53) =   0.7713566642727137D0
      HJT(54) =   0.7764908826289486D0
      HJT(55) =   0.7816384394572510D0
      HJT(56) =   0.7867990342880129D0
      HJT(57) =   0.7919723768377592D0
      HJT(58) =   0.7971581865926640D0
      HJT(59) =   0.8023561924123415D0
      HJT(60) =   0.8075661321527717D0
      HJT(61) =   0.8127877523072879D0
      HJT(62) =   0.8180208076646272D0
      HJT(63) =   0.8232650609831067D0
      HJT(64) =   0.8285202826800393D0
C
      RNUM(1)  = 2.015625000000000D0
      RNUM(2)  = 2.046875000000000D0
      RNUM(3)  = 2.078125000000000D0
      RNUM(4)  = 2.109375000000000D0
      RNUM(5)  = 2.140625000000000D0
      RNUM(6)  = 2.171875000000000D0
      RNUM(7)  = 2.203125000000000D0
      RNUM(8)  = 2.234375000000000D0
      RNUM(9)  = 2.265625000000000D0
      RNUM(10) = 2.296875000000000D0
      RNUM(11) = 2.328125000000000D0
      RNUM(12) = 2.359375000000000D0
      RNUM(13) = 2.390625000000000D0
      RNUM(14) = 2.421875000000000D0
      RNUM(15) = 2.453125000000000D0
      RNUM(16) = 2.484375000000000D0
      RNUM(17) = 2.515625000000000D0
      RNUM(18) = 2.546875000000000D0
      RNUM(19) = 2.578125000000000D0
      RNUM(20) = 2.609375000000000D0
      RNUM(21) = 2.640625000000000D0
      RNUM(22) = 2.671875000000000D0
      RNUM(23) = 2.703125000000000D0
      RNUM(24) = 2.734375000000000D0
      RNUM(25) = 2.765625000000000D0
      RNUM(26) = 2.796875000000000D0
      RNUM(27) = 2.828125000000000D0
      RNUM(28) = 2.859375000000000D0
      RNUM(29) = 2.890625000000000D0
      RNUM(30) = 2.921875000000000D0
      RNUM(31) = 2.953125000000000D0
      RNUM(32) = 2.984375000000000D0
      RNUM(33) = 3.015625000000000D0
      RNUM(34) = 3.046875000000000D0
      RNUM(35) = 3.078125000000000D0
      RNUM(36) = 3.109375000000000D0
      RNUM(37) = 3.140625000000000D0
      RNUM(38) = 3.171875000000000D0
      RNUM(39) = 3.203125000000000D0
      RNUM(40) = 3.234375000000000D0
      RNUM(41) = 3.265625000000000D0
      RNUM(42) = 3.296875000000000D0
      RNUM(43) = 3.328125000000000D0
      RNUM(44) = 3.359375000000000D0
      RNUM(45) = 3.390625000000000D0
      RNUM(46) = 3.421875000000000D0
      RNUM(47) = 3.453125000000000D0
      RNUM(48) = 3.484375000000000D0
      RNUM(49) = 3.515625000000000D0
      RNUM(50) = 3.546875000000000D0
      RNUM(51) = 3.578125000000000D0
      RNUM(52) = 3.609375000000000D0
      RNUM(53) = 3.640625000000000D0
      RNUM(54) = 3.671875000000000D0
      RNUM(55) = 3.703125000000000D0
      RNUM(56) = 3.734375000000000D0
      RNUM(57) = 3.765625000000000D0
      RNUM(58) = 3.796875000000000D0
      RNUM(59) = 3.828125000000000D0
      RNUM(60) = 3.859375000000000D0
      RNUM(61) = 3.890625000000000D0
      RNUM(62) = 3.921875000000000D0
      RNUM(63) = 3.953125000000000D0
      RNUM(64) = 3.984375000000000D0
C
C     *********************************************
C     *           DISPLAY PARAMETERS              *
C     *********************************************
C
 200  FORMAT(3H  |,8X,'PP =',1X,1F18.15,3H  |)
 210  FORMAT(3H  |,8X,'HM =',1X,1F18.15,3H  |)
 220  FORMAT(3H  |,' RKAP_AXIS =',1X,1F18.15,3H  |)
 230  FORMAT(3H  |,6X,'RKAP =',1X,1F18.15,3H  |)
 240  FORMAT(3H  |,6X,'DELT =',1X,1F18.15,3H  |)
 250  FORMAT(2H  ,' ------ CALCULATION RESULTS ------ ')
 260  FORMAT(2H  ,' --------------------------------- ')
      WRITE(6,250)
      WRITE(6,200) PP
      WRITE(6,210) HM
      WRITE(6,220) RKAPRA
      WRITE(6,230) RKAP
      WRITE(6,240) RDLT
      WRITE(6,260)
C
C     *********************************************
C     *               DRAW GRAPH                  *
C     *********************************************
C
      IF (MUL.EQ.64) THEN
         NCNT=0
         DO NX=1,NXMAX
            IF (MOD(NX,MUL).EQ.1) THEN
               GX((NX-1)/MUL+1,1)   = GUCLIP(RX(NX+NCNT))
               IF (NCHG.EQ.0) THEN
                  GY((NX-1)/MUL+1,1) 
     &                 = GUCLIP(-PSI(NX+NCNT,NZMAX/(MUL*2)))
               ELSE
                  GY((NX-1)/MUL+1,1) = GUCLIP(-PSIZAX(NX+NCNT))
               ENDIF
               GJY((NX-1)/MUL+1,1)  = GUCLIP(HJ(NX+NCNT))
               NCNT=NCNT+(MUL/64)
            ENDIF
         ENDDO
      ELSEIF (MUL.EQ.1) THEN
         DO N=1,NXMAX
            GX(N,1) = GUCLIP(RX(N))
            IF (NCHG.EQ.0) THEN
               GY(N,1) = GUCLIP(-PSI(N,NZMAX/2))
            ELSEIF (NCHG.EQ.1) THEN
               GY(N,1)  = GUCLIP(-PSIZAX(N))
            ENDIF
            GJY(N,1)  = GUCLIP(HJ(N))
         ENDDO
      ELSE
         NCNT=0
         DO NX=1,NXMAX
            IF (MOD(NX,MUL).EQ.1) THEN
               IF ((NX-1)/MUL+1.EQ.64) NCNT=NCNT+1
               GX((NX-1)/MUL+1,1)   = GUCLIP(RX(NX+NCNT))
               write(6,*) NX+NCNT,NCNT,(NX-1)/MUL+1,GX((NX-1)/MUL+1,1)
               IF (NCHG.EQ.0) THEN
                  GY((NX-1)/MUL+1,1) 
     &                 = GUCLIP(-PSI(NX+NCNT,NZMAX/(MUL*2)))
               ELSE
                  GY((NX-1)/MUL+1,1) = GUCLIP(-PSIZAX(NX+NCNT))
               ENDIF
               GJY((NX-1)/MUL+1,1)  = GUCLIP(HJ(NX+NCNT))
               NCNT=NCNT+(MUL/64)
            ENDIF
         ENDDO
      ENDIF
C     
      DO 110 N=1,NXMAX
         GX(N,2)  = GUCLIP(RNUM(N))
         GY(N,2)  = GUCLIP(PSIN(N))
         GJY(N,2) = GUCLIP(HJT(N))
  110 CONTINUE    
C
      CALL PAGES
      TXT='/ANALYTIC AND NUMERICAL PSI/'
      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GY,NXM,NXMAX/MUL,2,TXT,7.5)
      CALL PAGEE
C
      CALL PAGES
      TXT='/ANALYTIC AND NUMERICAL J/'
      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GJY,NXM,NXMAX/MUL,2,TXT,7.5)
      CALL PAGEE
      GOTO 9000
C
 9000 CALL GSCLOS
      STOP
      END
C
C     *********************************************
C     *               SUBROUTINES                 *
C     *********************************************
C
      SUBROUTINE GRAPH1(PXMIN,PXMAX,PYMIN,PYMAX,GX,GY,NXM,
     &                   NXMAX,NGMAX,TXT,POS)
C
      DIMENSION GX(NXM,2),GY(NXM,NGMAX)
      CHARACTER*80 TXT
C
      CALL GMNMX1(GX(1,1),1,NXMAX,1,XMIN,XMAX)
      CALL GMNMX1(GX(1,2),1,NXMAX,1,XMIN1,XMAX1)
      XMIN=MIN(XMIN,XMIN1)
      XMAX=MAX(XMAX,XMAX1)
      CALL GMNMX1(GY(1,1),1,NXMAX,1,YMIN,YMAX)
      DO 100 NG=2,NGMAX
         CALL GMNMX1(GY(1,NG),1,NXMAX,1,YMIN1,YMAX1)
         YMIN=MIN(YMIN,YMIN1)
         YMAX=MAX(YMAX,YMAX1)
  100 CONTINUE
      CALL GQSCAL(XMIN,XMAX,GXMIN,GXMAX,GXSCAL)
      CALL GQSCAL(YMIN,YMAX,GYMIN,GYMAX,GYSCAL)
C
C     Origin should have a scale mark
C
      IF(GXMIN*GXMAX.GT.0.0) THEN
         GXORG=GXMIN
      ELSE
         GXORG=0.0
      ENDIF
      IF(GYMIN*GYMAX.GT.0.0) THEN
         GYORG=GYMIN
      ELSE
         GYORG=0.0
      ENDIF
C
C     symbol font has better minus sign
C
      CALL SETFNT(2)
      CALL SETCHS(0.35,0.0)
      CALL GDEFIN(PXMIN,PXMAX,PYMIN,PYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,0,7)
      CALL GFRAME
C
C     NGULEN choose appropriate format of scale values
C
      CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.3,9)
      CALL GVALUE(GXORG,2*GXSCAL,0.0,0.0,NGULEN(2*GXSCAL))
      CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.3,9)
      CALL GVALUE(0.0,0.0,GYORG,2*GYSCAL,NGULEN(2*GYSCAL))
C      DO 200 NG=1,NGMAX
C         CALL SETLIN(-1,-1,7-MOD(NG-1,5))
         CALL SETLIN(-1,-1,7)
         CALL GPLOTP(GX(1,1),GY(1,1),1,NXMAX,1,0,0,0)
         CALL SETLIN(-1,-1,6)
         CALL GPLOTP(GX(1,2),GY(1,2),1,NXMAX,1,0,0,0)
C         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,MOD(NG-1,8))
C  200 CONTINUE
      CALL SETLIN(-1,-1,7)
C
      CALL GTEXTX(POS,17.3,TXT,2)
C
      RETURN
      END
C
C
      SUBROUTINE GRAPH2(PXMIN,PXMAX,PYMIN,PYMAX,GX,GY,NXM,
     &                   NXMAX,NGMAX,TXT,POS)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX)
      CHARACTER*80 TXT
C
      CALL GMNMX1(GX,1,NXMAX,1,XMIN,XMAX)
      CALL GMNMX1(GY(1,1),1,NXMAX,1,YMIN,YMAX)
      DO 100 NG=2,NGMAX
         CALL GMNMX1(GY(1,NG),1,NXMAX,1,YMIN1,YMAX1)
         YMIN=MIN(YMIN,YMIN1)
         YMAX=MAX(YMAX,YMAX1)
  100 CONTINUE
      CALL GQSCAL(XMIN,XMAX,GXMIN,GXMAX,GXSCAL)
      CALL GQSCAL(YMIN,YMAX,GYMIN,GYMAX,GYSCAL)
C
C     Origin should have a scale mark
C
      IF(GXMIN*GXMAX.GT.0.0) THEN
         GXORG=GXMIN
      ELSE
         GXORG=0.0
      ENDIF
      IF(GYMIN*GYMAX.GT.0.0) THEN
         GYORG=GYMIN
      ELSE
         GYORG=0.0
      ENDIF
C
C     symbol font has better minus sign
C
      CALL SETFNT(2)
      CALL SETCHS(0.35,0.0)
      CALL GDEFIN(PXMIN,PXMAX,PYMIN,PYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,0,7)
      CALL GFRAME
C
C     NGULEN choose appropriate format of scale values
C
      CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.3,9)
      CALL GVALUE(GXORG,2*GXSCAL,0.0,0.0,NGULEN(2*GXSCAL))
      CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.3,9)
      CALL GVALUE(0.0,0.0,GYORG,2*GYSCAL,NGULEN(2*GYSCAL))
      DO 200 NG=1,NGMAX
         CALL SETLIN(-1,-1,7-MOD(NG-1,5))
         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,MOD(NG-1,8))
  200 CONTINUE
      CALL SETLIN(-1,-1,7)
C
      CALL GTEXTX(POS,17.3,TXT,2)
C
      RETURN
      END
