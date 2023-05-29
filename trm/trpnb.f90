!     ***********************************************************

!           NEUTRAL BEAM

!     ***********************************************************

      SUBROUTINE TRPWNB

      USE TRCOMM, ONLY : MDLNB, NRMAX, SNB, TAUB
      IMPLICIT NONE

      IF(MDLNB.EQ.0) THEN
         TAUB(1:NRMAX)=1.D0
      ELSEIF(MDLNB.EQ.1) THEN
         CALL TRNBIA
         SNB(1:NRMAX)=0.D0
         CALL TRAJNB
      ELSEIF(MDLNB.EQ.2) THEN
         CALL TRNBIA
         CALL TRAJNB
      ELSEIF(MDLNB.EQ.3) THEN
         CALL TRNBIB
         SNB(1:NRMAX)=0.D0
         CALL TRAJNB
      ELSEIF(MDLNB.EQ.4) THEN
         CALL TRNBIB
         CALL TRAJNB
      ENDIF
      RETURN
      END SUBROUTINE TRPWNB

!     ***********************************************************

!           NEUTRAL BEAM : (GAUSSIAN PROFILE)

!     ***********************************************************

      SUBROUTINE TRNBIA

      USE TRCOMM, ONLY : DR, DVRHO, NRMAX, PNB, PNBENG, PNBR0, PNBRW, PNBTOT, RA, &
           &             RKEV, RM, SNB, rkind
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind)   :: PNB0, SUM

      IF(PNBTOT.LE.0.D0) RETURN

      SUM = 0.D0
      DO NR=1,NRMAX
         SUM=SUM+DEXP(-((RA*RM(NR)-PNBR0)/PNBRW)**2)*DVRHO(NR)*DR
      ENDDO

      PNB0=PNBTOT*1.D6/SUM

      DO NR=1,NRMAX
         PNB(NR)=PNB0*DEXP(-((RA*RM(NR)-PNBR0)/PNBRW)**2)
         SNB(NR)=PNB(NR)/(PNBENG*RKEV)*1.D-20
      ENDDO

      RETURN
      END SUBROUTINE TRNBIA

!     ***********************************************************

!           NEUTRAL BEAM : (PENCIL BEAM)

!     ***********************************************************

      SUBROUTINE TRNBIB

      USE TRCOMM, ONLY : GPNB, NRMAX, PNB, PNBENG, PNBRTG, PNBRW, PNBTOT, PNBVW, PNBVY, &
           &             RKEV, RTG, SNB, NRNBMAX, rkind
      IMPLICIT NONE
      INTEGER:: I, J, NRNB
      INTEGER,SAVE:: NRNBMAX_save=0
      REAL(rkind)   ::DRTG, DVY, RDD, VY,sum,x
      REAL(rkind),ALLOCATABLE :: AR(:)

      IF(PNBTOT.LE.0.D0) RETURN

      ALLOCATE(AR(NRNBMAX))
      IF(NRNBMAX.NE.NRNBMAX_save) THEN
         IF(ALLOCATED(GPNB)) DEALLOCATE(GPNB)
         ALLOCATE(GPNB(4*NRMAX,NRNBMAX))
         sum=0.D0
         DO NRNB=1,NRNBMAX
            x=2.D0*(DBLE(NRNB)/DBLE(NRNBMAX+1)-0.5D0)  ! -1.0..1.0
            AR(NRNB)=EXP(-3.D0*x**2)
            sum=sum+AR(NRNB)
         END DO
         DO NRNB=1,NRNBMAX
            AR(NRNB)=AR(NRNB)/sum
         END DO
         NRNBMAX_save=NRNBMAX
      END IF
      SNB(1:NRMAX) = 0.D0
      GPNB(1:4*NRMAX,1:NRNBMAX) = 0.0

      DRTG=2.D0*PNBRW/NRNBMAX
      DVY =2.D0*PNBVW/NRNBMAX
      DO I=1,NRNBMAX
         DO J=1,NRNBMAX
            RDD=AR(J)*AR(I)
            RTG(J)=PNBRTG-PNBRW+0.5D0*DRTG+DRTG*(J-1)
            VY    =PNBVY -PNBVW+0.5D0*DVY +DVY *(I-1)
            CALL TRNBPB(J,RTG(J),VY,RDD)
         ENDDO
      ENDDO

      PNB(1:NRMAX) = SNB(1:NRMAX)*1.D20*PNBENG*RKEV
      DEALLOCATE(AR)

      RETURN
      END SUBROUTINE TRNBIB

!     ***********************************************************

!           NEUTRAL BEAM CALCULATION

!     ***********************************************************

      SUBROUTINE TRNBPB(J,R0,VY,RDD)

!     J   (in): J-th NBI
!     R0  (in): tangential radius of NBI beam (m)
!     VY  (in): vertical position of NBI (m)
!     RDD (in): NBI deposition rate

!     KL  : judging condition parameter
!           0 : beam energy has decreased at stopping condition
!           1 : beam energy has not decreased at stopping condition yet
!           2 : median center of beam line

      USE TRCOMM, ONLY : ANC, ANFE, DR, DVRHO, GBAN, GBL, GBP1, GBR, GBRH, GPNB, &
           &             NLMAX, NRMAX, PNBENG, PNBTOT, PZC, &
           &             PZFE, RA, RG, RKEV, RN, RR, RT, SNB, rkind
      IMPLICIT NONE
      integer, intent(in) :: J
      real(rkind), intent(in) :: R0, VY, RDD
      INTEGER:: I, IB, IDL, IM, KL
      REAL(rkind)   :: ANL, ANL0, COST, COSTV, DCX, DEX, DFX, DL, DOX, DRG, P1, P1SUM, &
           &       RADIUS1, RADIUS2, RADIUSG, &
     &             RAL, RG1, RG2, SGM, SUML, TEX, XL, ZCX, ZFX, ZOX
      REAL   :: GUCLIP   ! GSAF

      IF(PNBTOT <= 0.D0) RETURN

!  COSTV : cosine between midplane and vertical position of NB
      COSTV=SQRT(RA**2-VY**2)/RA
!  RAL : minor radius of NB injection point projected to the midplane
      RAL=RA*COSTV
!  XL : maximum distance from injection point to wall
      IF(R0 >= RR-RAL) THEN
         XL=2.D0*SQRT((RR+RAL)**2-R0**2)
      ELSE
         XL=SQRT((RR+RAL)**2-R0**2)-SQRT((RR-RAL)**2-R0**2)
      ENDIF
!  ANL : beam intensity [num/s]
!        (the number of beam particles per unit length
!         divided by thier velocity)
      ANL=RDD*PNBTOT*1.D6/(PNBENG*RKEV*1.D20)
!  IB : number of vertical grid point from midplane to NB injection point,
!       which defines innermost radial grid point of the heating region for each NB chord
      IB = INT(ABS(VY/(DR*RA)))+1
!  I : radial grid point of current NB deposition position
      I=NRMAX-IB+1

!  SUML : total distance from injection point in eye direction of NB
      SUML=0.D0
!  ANL0 : stored ANL for truncation of calculation
      ANL0=ANL
!  DL : arbitrary minute length in direction of NB
!      DL=XL/1000
      DL=XL/500
!  COST : cosine between midplane and NB chord
      COST=SQRT((RR+RAL)**2-R0**2)/(RR+RAL)
!  IDL : serial number of DL
      IDL=0

!  IM : radial grid point turned back at magnetic axis
   10 IF(I > 0) THEN
         IM=I+IB-1
      ELSE
         IF(ABS(I) <= NRMAX) THEN
            IM=IB+ABS(I)-1
         ELSEIF(ABS(I) <= 2*NRMAX) THEN
            IM=IB-ABS(I)+2*NRMAX
         ELSE
            IM=IB+ABS(I)-2*NRMAX-1
         ENDIF
      ENDIF

!      WRITE(6,'(3(1X,I4),3F15.7)') IB,I,IM,SUML,DL,XL

!      IF(I.NE.0.AND.(IM.GE.1.AND.IM.LE.NRMAX)) THEN
      IF(I /= 0 .AND. (IM >= IB .AND. IM <= NRMAX)) THEN

      P1SUM = 0.D0
      RG1=(RG(IM  )*RA)**2-VY**2
      RG2=(RG(IM-1)*RA)**2-VY**2
      IF(RG2 < 0.D0) RG2=0.D0 ! Turn-around point
      IF(I > 0) THEN
         RADIUS1 = RR+SQRT(RG1) ! First LFS
      ELSE
         IF(ABS(I) <= NRMAX) THEN
            RADIUS1 = RR-SQRT(RG2) ! First HFS
         ELSEIF(ABS(I) <= 2*NRMAX) THEN
            RADIUS1 = RR-SQRT(RG2) ! Second HFS after turn-around
         ELSE
            RADIUS1 = RR+SQRT(RG1) ! Second LFS after turn-around
         ENDIF
      ENDIF
      IF(I == -3*NRMAX) THEN ! Shine-through
         NLMAX(J)=IDL ! Maximum number of points in the eye direction of NB
         RETURN
      ENDIF

      ! Plasma parameters
      DEX = RN(IM,1)*1.D1
      TEX = RT(IM,1)
      DCX = ANC(IM)*1.D1
      DFX = ANFE(IM)*1.D1
      DOX = 0.D0
      ZCX = PZC(IM)
      ZFX = PZFE(IM)
      ZOX = 8.D0

!  Accumulate P1 inside one grid (P1 denotes the deposition power at each region)
 20   IF(SUML+DL >= XL) DL=XL-SUML ! Shine-through
      IDL=IDL+1
      GBL (IDL)  =GUCLIP(SUML)
      GBAN(IDL,J)=GUCLIP(ANL)

      ! Calculate beam stopping cross-section
      CALL TRBSCS(DEX,TEX,PNBENG,DCX,DFX,DOX,ZCX,ZFX,ZOX,SGM)

      P1=(DEX*0.1D0)*SGM*ANL*DL
      GBP1(IDL,J)=GUCLIP(P1)
      IF(P1 < 0.D0) P1=0.D0
      IF(ANL > ANL0*1.D-3) THEN ! Beam intensity still remains.
         KL=1
      ELSE ! All the power is absorbed in a plasma.
         P1=ANL
         KL=0
      ENDIF

      SUML=SUML+DL
!  inside the torus
      IF(SUML < XL) THEN
         IF(KL == 1) THEN
            DRG=SQRT(RG1)-SQRT(RG2)
            ! RADIUSG: Tangential NBI radius where NB existed at one step past.
            RADIUSG=SQRT( (SUML-DL)**2 +(RR+RAL)*(RR+RAL-2.D0*(SUML-DL)*COST))
            GBR (IDL,J)=GUCLIP(RADIUSG)
            GBRH(IDL,J)=GUCLIP(ABS((RADIUSG-RR)/RAL))
            ! RADIUS2: Tangential NBI radius where NB currently exists.
            RADIUS2=SQRT(SUML**2+(RR+RAL)*(RR+RAL-2.D0*SUML*COST))
!            write(6,'(2I4,5F13.7)') I,IM,ABS(RADIUS1-RADIUS2),DRG -ABS(RADIUS1-RADIUS2),RADIUS2-R0,RADIUS1,RADIUS2
!            write(6,'(2I4,5F13.7)') I,IM,DRG-ABS(RADIUS1-RADIUS2),RADIUS1,RADIUS2,RADIUS2-R0,SUML
!            write(6,'(2I4,5F13.7)') I,IM,ANL,P1,P1SUM,DRG-ABS(RADIUS1-RADIUS2),RADIUS2-R0
!            write(6,'(2I4,5F13.7)') I,IM,DRG-ABS(RADIUS1-RADIUS2),RADIUS1,RADIUS2,SUML,RADIUS2-R0
!            write(6,'(2I4,3F13.7)') I,IM,RG1,RG2,DRG
            IF(RADIUS2-R0 > 1.D-6) THEN
               IF(DRG-ABS(RADIUS1-RADIUS2) > 1.D-6) THEN
!  inside the grid
                  ANL=ANL-P1
                  P1SUM=P1SUM+P1
                  GOTO 20
               ELSE
!  run off the grid
                  P1=P1SUM
                  IDL=IDL-1
                  SUML=SUML-DL
               ENDIF
            ELSE
!  innermost grid
               ANL=ANL-P1
               P1=P1SUM
               KL=2
            ENDIF
         ENDIF
      ENDIF

!  NB source flux (Half mesh)
      SNB(IM) = SNB(IM)+P1/(DVRHO(IM)*DR)
!!$      IF(IM.GT.1) SNB(IM-1) = SNB(IM-1)+0.25D0*P1/(DVRHO(IM-1)*DR)
!!$      IF(IM.GT.1.AND.IM.LT.NRMAX) THEN
!!$         SNB(IM  ) = SNB(IM  )+0.5D0 *P1/(DVRHO(IM  )*DR)
!!$      ELSEIF(IM.EQ.1.OR.IM.EQ.NRMAX) THEN
!!$         SNB(IM  ) = SNB(IM  )+       P1/(DVRHO(IM  )*DR)
!!$      ENDIF
!!$      IF(IM.LT.NRMAX) SNB(IM+1) = SNB(IM+1)+0.25D0*P1/(DVRHO(IM+1)*DR)

!  for graphics
      IF(I > 0) THEN
         GPNB(IM        ,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ELSEIF(I < -1 .AND. I >= -NRMAX) THEN
         GPNB(IM+  NRMAX,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ELSEIF(I <= -NRMAX-1 .AND. I >= -2*NRMAX) THEN
         GPNB(IM+2*NRMAX,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ELSE
         GPNB(IM+3*NRMAX,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ENDIF

!      WRITE(6,'(4(1X,I4),4F15.7)') J,I,KL,IM,SUML,DL,ANL,P1
      IF(KL == 0) THEN !  below the lower limits of intensity
         NLMAX(J)=IDL-1
         RETURN
      ELSE IF(KL == 2) THEN !  inversion of the label number because of innermost grid
         IF(I > 0) THEN
            I=-2*NRMAX-ABS(I)
         ELSE
            I=-2*(NRMAX-ABS(I))+I
            IF(ABS(I) < NRMAX) I=I-NRMAX
         ENDIF
         GOTO 10
      ENDIF

      ENDIF

      IF(XL-(SUML+DL) > 1.D-6) THEN
!  inside the torus
         I=I-1
      ELSE
!  outside the torus
         NLMAX(J)=IDL-1
         RETURN
      ENDIF

      GOTO 10

      END SUBROUTINE TRNBPB

!     ***********************************************************

!           CALCULATE BEAM STOPPING CROSS-SECTIONS

!     ***********************************************************

      SUBROUTINE TRBSCS(ANEX,TEX,EB,ANCX,ANFEX,ANOX,PZCX,PZFEX,PZOX,SGM)

!     R.K. Janev, C.D. BOLEY and D.E. Post, Nucl. Fusion 29 (1989) 2125
!     < Input >
!        ANEX  : ELECTRON DENSITY     [10^19]
!        TEX   : ELECTRON TEMPERATURE [keV]
!        EB    : BEAM ENERGY          [keV/u]
!        ANCX  : CARBON DENSITY       [10^19]
!        ANFEX : IRON DENSITY         [10^19]
!        ANOX  : OXIGEN DENSITY       [10^19]
!        PZCX  : CARBON CHARGE NUMBER
!        PZFEX : IRON CHARGE NUMBER
!        PZOX  : OXIGEN CHARGE NUMBER
!     < Output >
!        SGM   : BEAM STOPPING CROSS-SECTION [10^-20 m^2]

!     NOTE: "C" DENOTES CARBON, "F" IRON AND "O" OXIGEN.

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind), INTENT(IN) :: ANEX, TEX, EB, ANCX, ANFEX, ANOX, PZCX, PZFEX, PZOX
      REAL(rkind), INTENT(OUT):: SGM
      INTEGER          ::  I, J, K
      REAL(rkind)             ::  S1, S2, S3, S4
      REAL(rkind),DIMENSION(2,3,2):: A
      REAL(rkind),DIMENSION(2,3,2):: C, F, O

!     DATA FORMAT FOR A           DATA FORMAT FOR B
!         111  211  121  221          111  211  311  121
!         131  231  112  212          221  321  112  212
!         122  222  132  232          312  122  222  322

      DATA A/ 4.40D+00, 2.30D-01, 7.46D-02,-2.55D-03, &
     &        3.16D-03, 1.32D-03,-2.49D-02,-1.15D-02, &
     &        2.27D-03,-6.20D-04,-2.78D-05, 3.38D-05/

      DATA C/-1.49D+00, 5.18D-01,-3.36D-02,-1.19D-01, &
     &        2.92D-02,-1.79D-03,-1.54D-02, 7.18D-03, &
     &        3.41D-04,-1.50D-02, 3.66D-03,-2.04D-04/

      DATA F/-1.03D+00, 3.22D-01,-1.87D-02,-5.58D-02, &
     &        1.24D-02,-7.43D-04, 1.06D-01,-3.75D-02, &
     &        3.53D-03,-3.72D-03, 8.61D-04,-5.12D-05/

      DATA O/-1.41D+00, 4.77D-01,-3.05D-02,-1.08D-01, &
     &        2.59D-02,-1.57D-03,-4.08D-04, 1.57D-03, &
     &        7.35D-04,-1.38D-02, 3.33D-03,-1.86D-04/

      S1=0.D0
      S2=0.D0
      S3=0.D0
      S4=0.D0
      DO I=1,2
      DO J=1,3
      DO K=1,2
         S1=S1+A(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1) &
     &                 *(LOG(TEX))**(K-1)
      ENDDO
      ENDDO
      ENDDO
      DO I=1,2
      DO J=1,3
      DO K=1,2
         S2=S2+C(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)*(LOG(TEX))**(K-1)
         S3=S3+F(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)*(LOG(TEX))**(K-1)
         S4=S4+O(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)*(LOG(TEX))**(K-1)
      ENDDO
      ENDDO
      ENDDO

      SGM=EXP(S1)/EB*(1.D0+(ANCX *PZCX *(PZCX -1.D0)*S2 &
     &                     +ANFEX*PZFEX*(PZFEX-1.D0)*S3 &
     &                     +ANOX *PZOX *(PZOX -1.D0)*S4)/ANEX)

      RETURN
      END SUBROUTINE TRBSCS

!     ***********************************************************

!          NEUTRAL BEAM DRIVEN CURRENT

!     ***********************************************************

      SUBROUTINE TRAJNB

        USE TRCOMM, ONLY : &
             AEE, AJNB, AME, AMM, ANC, ANFE, EPS0, EPSRHO, MDLUF, NRMAX, &
             PA, PBCL, PBIN, PI, PNBCD, PNBENG, &
             PZ, PZC, PZFE, RKEV, RN, RNF, RT, RTF, RW, TAUB, ZEFF, rkind
      IMPLICIT NONE
      REAL(rkind)    :: AMA, AMB, AMD, AMT, ANE, COULOG, EC, EPS, HY, HYB, &
           P2, P3, P4, PAB, PZB, TAUS, TAUS0, TE, VB, VC3,  &
           VCA3, VCD3, VCR, VCT3, VE, WB, XB, ZEFFM, ZN
      INTEGER :: NR


      AMD=AMM*PA(2)
      AMT=AMM*PA(3)
      AMA=AMM*PA(4)
      AMB=AMM*PA(2)
      PAB=PA(2)
      PZB=PZ(2)
      VB=SQRT(2.D0*PNBENG*RKEV/AMB)

      IF(MDLUF.NE.0) THEN
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         IF(ANE.EQ.0.D0) THEN
            P4=0.D0
            TAUS=0.D0
         ELSE
         P4 = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
         VCD3 = P4*RN(NR,2)*PZ(2)**2/AMD
         VCT3 = P4*RN(NR,3)*PZ(3)**2/AMT
         VCA3 = P4*RN(NR,4)*PZ(4)**2/AMA
         VC3  = VCD3+VCT3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYB  = HY(VB/VCR)
         TAUS = 0.2D0*PAB*ABS(TE)**1.5D0 /(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))
         TAUB(NR) = 0.5D0*TAUS*(1.D0-HYB)
         ENDIF
      ENDDO
      ELSE
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
!         WB =RW(NR,1)*1.5D0
         WB =RW(NR,1)
         IF(ANE.EQ.0.D0) THEN
            P4=0.D0
            TAUS=0.D0
         ELSE
         P4 = 3.D0*SQRT(0.5D0*PI)*AME/ANE*(ABS(TE)*RKEV/AME)**1.5D0
         VCD3 = P4*RN(NR,2)*PZ(2)**2/AMD
         VCT3 = P4*RN(NR,3)*PZ(3)**2/AMT
         VCA3 = P4*RN(NR,4)*PZ(4)**2/AMA
         VC3  = VCD3+VCT3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYB  = HY(VB/VCR)
         TAUS = 0.2D0*PAB*ABS(TE)**1.5D0/(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))
         TAUB(NR) = 0.5D0*TAUS*(1.D0-HYB)
         RNF(NR,1)= 2.D0*LOG(1.D0+(VB/VCR)**3)*WB/(3.D0*(1.D0-HYB)*PNBENG)
         ENDIF

         IF(RNF(NR,1).GT.0.D0) THEN
!            RTF(NR,1)= WB/(1.5D0*RNF(NR,1))
            RTF(NR,1)= WB/RNF(NR,1)
         ELSE
            RTF(NR,1)= 0.D0
         ENDIF
         PBIN(NR)   = WB*RKEV*1.D20/TAUB(NR)
         PBCL(NR,1) =   (1.D0-HYB)*PBIN(NR)
         PBCL(NR,2) = VCD3/VC3*HYB*PBIN(NR)
         PBCL(NR,3) = VCT3/VC3*HYB*PBIN(NR)
         PBCL(NR,4) = VCA3/VC3*HYB*PBIN(NR)
      ENDDO

      IF(PNBCD.LE.0.D0) RETURN

!     D. R .Mikkelsen and C. E. Singer, J. Plasma Phys. 4 237 (1983)
!        xi_0 corresponds to PNBCD
!        H(r)*P_b/V_p corresponds to PBIN(NR)
!
      TAUS0=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*AMB*AME &
     &     /(1.D20*AEE**4*PZB**2*COULOG(1,2,ANE,TE))
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         EPS = EPSRHO(NR)
         VE  = SQRT(ABS(TE)*RKEV/AME)
         IF(ANE.EQ.0.D0) THEN
            TAUS=0.D0
            ZEFFM=0.D0
            XB=0.D0
            AJNB(NR)=0.D0
         ELSE
         TAUS=TAUS0*VE**3/ANE
         ZEFFM = (PZ(2)  *PZ(2)  *RN(NR,2)/PA(2) &
     &           +PZ(3)  *PZ(3)  *RN(NR,3)/PA(3) &
     &           +PZ(4)  *PZ(4)  *RN(NR,3)/PA(4) &
     &           +PZC(NR) *PZC(NR) *ANC(NR)/12.D0 &
     &           +PZFE(NR)*PZFE(NR)*ANFE(NR)/52.D0)/ANE
         EC  = 14.8D0*TE*PAB*ZEFFM**(2.D0/3.D0)
         VCR = VB*SQRT(ABS(EC)/PNBENG)
         P2  = (1.55D0+0.85D0/ZEFF(NR))*SQRT(EPS) &
     &        -(0.2D0+1.55D0/ZEFF(NR))*EPS
         XB  = VB/VCR
         ZN  = 0.8D0*ZEFF(NR)/PAB
         P3  = XB*XB/(4.D0+3.D0*ZN+XB*XB*(XB+1.39D0+0.61D0*ZN**0.7D0))

         AJNB(NR) = PNBCD*2.D0*AEE*PZB*TAUS/(AMB*VCR) &
     &            *(1.D0-PZB*(1.D0-P2)/ZEFF(NR))*P3*PBIN(NR)
      ENDIF
      ENDDO

      ENDIF

      RETURN
      END SUBROUTINE TRAJNB

!     ***********************************************************

      FUNCTION HY(V)

      USE TRCOMM, ONLY : PI,rkind
      IMPLICIT NONE
      REAL(rkind), INTENT(IN) :: V
      REAL(rkind) :: HY

      HY = 2.D0*(LOG((V**3+1.D0)/(V+1.D0)**3)/6.D0 &
     &      +(ATAN((2.D0*V-1.D0)/SQRT(3.D0))+PI/6.D0)/SQRT(3.D0))/V**2
      RETURN
      END FUNCTION HY
