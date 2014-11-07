!     $Id$
!
! ************************************************************
!
!      CALCULATION OF NONLINEAR COLLISIONAL OPERATOR
!
! ************************************************************
!
      MODULE fpcalcn

      USE fpcomm

      contains

!----------------------------------

      SUBROUTINE FPCALC_NL(NR,NSB,NSA)
!
      USE libgrf,ONLY: grd1d
      USE libspf, ONLY: dpleg
      USE fpmpi
      IMPLICIT NONE
!
!      integer,parameter::N=NPM+2, M=NTHM+2, LNM=5
      integer,parameter::LNM=5
      real(8),dimension(NTHMAX+3,-1:LNM):: PLM,PLG
      real(8),DIMENSION(NTHMAX+3,-1:LNM):: D1PLM, D1PLG, D2PLG
      real(8),DIMENSION(0:LNM):: PLTEMP
      real(8),DIMENSION(NPSTART:NPEND):: FPLL 
      real(8),DIMENSION(NPMAX):: FPL_recv
      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL
      double precision,dimension(-1:LNM):: FPLS1
      double precision:: FPLS1_temp
!      real(8),DIMENSION(NPMAX+3,-1:LNM):: RM1M, RM2M, RM3M, RM4M
!      real(8),DIMENSION(NPMAX+3,-1:LNM):: RM1G, RM2G, RM3G, RM4G
      real(8),DIMENSION(NPSTARTW:NPENDWM,-1:LNM):: RM1M, RM2M, RM3M, RM4M
      real(8),DIMENSION(NPSTART :NPENDWG,-1:LNM):: RM1G, RM2G, RM3G, RM4G
      real(8),DIMENSION(NTHMAX+3):: TX,TY,DF
      real(8),dimension(4,NTHMAX+3):: UTY
      real(8),DIMENSION(NTHMAX+3):: UTY0
      real(8),DIMENSION(NPMAX+3):: TX1, TY1, DF1
      real(8),dimension(4,NPMAX+3):: UTY1
      real(8),DIMENSION(NPMAX+3):: UTY10
      real(8),DIMENSION(NPMAX+3,-1:LNM):: PHYM, PSYM, D1PSYM
      real(8),DIMENSION(NPMAX+3,-1:LNM):: PSYG, D1PHYG, D1PSYG, D2PSYG

      integer:: NR, NSB, NSA, NTH, L,  NP, NNP, NPG, NSBA, NPS
      integer:: N,M
      real(8):: RGAMH, FACT, WA, WC, WB, WD, WE, WF
      real(8):: SUM1, PSUM, SUM2, SUM3, SUM4, SUM5, SUM6, SUM7, SUM8, SUM9
      real(8):: RGAMA, RGAMB, vtatb, pabbar, ptatb, PCRIT, pabar, pbbar, pgbar, pmbar
      integer:: LLMIN, IER

!      N=NPMAX+2
!      M=NTHMAX+2
!      LNM=5

!
!     ----- definition of local quantities -----
!
      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
           /(RNFP0(NSA)*PTFP0(NSA)*1.D20)*RNFD0(NSB)
      vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
      pabbar=(PTFP0(NSA)*AMFD(NSB))/(PTFD0(NSB)*AMFP(NSA))
      NSBA=NSB_NSA(NSA)

!
!     ----- calculation of Legendre Polynomials -----
!
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
            D2PLG(NTH,L)=-(L/(SING(NTH)**2)+L**2)*PLG(NTH,L) &
                 +L*COSG(NTH)/(SING(NTH)**2)*PLG(NTH,L-1) 
         END DO
      END DO

      IF(MOD(IDBGFP,2).EQ.1) THEN
!
!     +++ plot of Legendre polynomials and their derivatives +++
!
         CALL PAGES
         CALL GRD1D(1,thm,plm,M,NTHMAX,LLMAX+2,'@PLM:@',0)
         CALL GRD1D(2,thm,d1plm,M,NTHMAX,LLMAX+2,'@D1PLM:@',0)
         CALL PAGEE
!
         CALL PAGES
         CALL GRD1D(1,thg,plg,M,NTHMAX+1,LLMAX+2,'@PLG:@',0)
         CALL GRD1D(2,thg,d1plg,M,NTHMAX+1,LLMAX+2,'@D1PLG:@',0)
         CALL GRD1D(3,thg,d2plg,M,NTHMAX+1,LLMAX+2,'@D2PLG:@',0)
         CALL PAGEE
      ENDIF
!
!     ----- Legendre expansion of distribution funstion FNS -----
!
      CALL mtx_set_communicator(comm_np)
      DO L=LLMIN,LLMAX
         DO NP=NPSTART,NPEND
            TX(1)=0.D0
            TY(1)=0.D0
            DO NTH=1,NTHMAX
               TX(NTH+1)=THM(NTH)
               TY(NTH+1)=FNSB(NTH,NP,NR,NSB)*PLM(NTH,L)*SINM(NTH)
            END DO
            TX(NTHMAX+2)=PI
            TY(NTHMAX+2)=0.D0
            DF(1)= FNSB(1,NP,NR,NSB)*PLM(1,L)!*COSM(1)
            DF(NTHMAX+2)= (-1)**(L+1)*FNSB(NTHMAX,NP,NR,NSB)
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPLL(NP)=0.5D0*(2*L+1.D0)*SUM1
         END DO
!      END DO
         CALL fpl_comm(FPLL,FPL_recv) 
         DO NP=1,NPMAX
            FPL(NP,L)=FPL_recv(NP)
         END DO

!         TX(1)=0.D0
!         TY(1)=0.D0
!         FPLS1_temp=FPMXWL_calcn(0.D0,NR,NSB)
!         DO NTH=1,NTHMAX
!            TX(NTH+1)=THM(NTH)
!            TY(NTH+1)=FPLS1_temp*PLM(NTH,L)*SINM(NTH)
!         END DO
!         TX(NTHMAX+2)=PI
!         TY(NTHMAX+2)=0.D0
!         DF(1)= FPLS1_temp
!         DF(NTHMAX+2)= (-1)**(L+1)*FPLS1_temp ! satisfy until l=2
!         CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
!         CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
!         CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
!         FPLS1(L)=0.5D0*(2.D0*L+1.D0)*SUM1
      END DO ! LLMAX
      CALL mtx_reset_communicator

!
!      CALL INTEGRATION_BACKGROUND_F(FPL,NR,NSB,NSA,RM1M,RM2M,RM3M,RM4M,RM1G,RM2G,RM3G,RM4G)
      CALL INTEGRATION_BACKGROUND_F_FINE(FPL,NR,NSB,NSA,RM1M,RM2M,RM3M,RM4M,RM1G,RM2G,RM3G,RM4G)
!

      IF(MOD(IDBGFP/2,2).EQ.1) THEN
!
!     +++ plot of Legendre expansion and M_l, N_l +++
!
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
!
!     ----- calculation of phi_l, psi_l and their derivatives -----
!
      DO L=LLMIN,LLMAX
!         DO NP=1,NPMAX
         DO NP=NPSTART,NPENDWM
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            pabar=PTFP0(NSA)/(AMFP(NSA)*RGAMA)
            pbbar=PTFD0(NSB)/AMFD(NSB)
            pmbar=PM(NP,NSBA)*pabar/pbbar
            PHYM(NP,L)=-1.D0/(2*L+1) &
                     *((pmbar**(-L-1))*RM2M(NP,L) &
                      +(pmbar**L      *RM1M(NP,L))) &
                     *pabbar

            PSYM(NP,L)=-0.5D0/(2*L+1) &
                     *(1.D0/(2*L+3)*((pmbar**(-L-1))*RM4M(NP,L) &
                                    +(pmbar**( L+2))*RM1M(NP,L)) &
                      -1.D0/(2*L-1)*((pmbar**(-L+1))*RM2M(NP,L) &
                                    +(pmbar**L     )*RM3M(NP,L))) &
                     /pabbar

            D1PSYM(NP,L)=-0.5D0/(2*L+1) &
                     *(1.D0/(2*L+3)*( (L+2)*(pmbar**( L+1))*RM1M(NP,L) &
                                     -(L+1)*(pmbar**(-L-2))*RM4M(NP,L)) &
                      -1.D0/(2*L-1)*(  L   *(pmbar**( L-1))*RM3M(NP,L) &
                                     -(L-1)*(pmbar**(-L  ))*RM2M(NP,L))) &
                     /RGAMA**3
         END DO
      END DO


!      DO 182 L=LLMIN,LLMAX
!        NP=1
!        PSYG(NP,L)=0.D0
!        D1PHYG(NP,L)=0.D0
!        D1PSYG(NP,L)=0.D0
!        D2PSYG(NP,L)=0.D0
!  182 CONTINUE

      IF(NPSTART.eq.1)THEN
         NPS=2
      ELSE
         NPS=NPSTART
      END IF
      DO L=LLMIN,LLMAX
!         DO NP=1,NPMAX+1
         DO NP=NPS,NPENDWG
            RGAMA=SQRT(1.D0+PG(NP,NSBA)**2*THETA0(NSA))
            pabar=PTFP0(NSA)/(AMFP(NSA)*RGAMA)
            pbbar=PTFD0(NSB)/AMFD(NSB)
!            IF(NP.EQ.1) THEN
!               pgbar=0.001D0*PG(2,NSBA)*pabar/pbbar
!            ELSE
               pgbar=PG(NP,NSBA)*pabar/pbbar
!            ENDIF
            PSYG(NP,L)=-0.5D0/(2*L+1) &
                 *(1.D0/(2*L+3)*((pgbar**(-L-1))*RM4G(NP,L) &
                                +(pgbar**( L+2))*RM1G(NP,L)) &
                  -1.D0/(2*L-1)*((pgbar**(-L+1))*RM2G(NP,L) &
                                +(pgbar**  L   )*RM3G(NP,L))) &
                 /pabbar

            D1PHYG(NP,L)=-1.D0/(2*L+1) &
                 *(  L   *(pgbar**( L-1))*RM1G(NP,L) &
                   -(L+1)*(pgbar**(-L-2))*RM2G(NP,L)) &
                 *pabbar**2

            D1PSYG(NP,L)=-0.5D0/(2*L+1) &
                 *(1.D0/(2*L+3)*( (L+2)*(pgbar**( L+1))*RM1G(NP,L) &
                                 -(L+1)*(pgbar**(-L-2))*RM4G(NP,L)) &
                  -1.D0/(2*L-1)*(  L   *(pgbar**( L-1))*RM3G(NP,L) &
                                 -(L-1)*(pgbar**(-L  ))*RM2G(NP,L)))

            D2PSYG(NP,L)=-0.5D0/(2*L+1) &
                 *(DBLE(L+1)*(L+2)/(2*L+3)*((pgbar**(-L-3))*RM4G(NP,L) &
                                           +(pgbar**  L   )*RM1G(NP,L)) &
                  -DBLE(L  )*(L-1)/(2*L-1)*((pgbar**(-L-1))*RM2G(NP,L) &
                                           +(pgbar**( L-2))*RM3G(NP,L))) &
                 *pabbar
         END DO
         IF(NPSTART.eq.1)THEN
            NP=1
            pabar=PTFP0(NSA)/AMFP(NSA)
            pbbar=PTFD0(NSB)/AMFD(NSB)
            pgbar=PG(1,NSBA)*pabar/pbbar
            PSYG(NP,L)=-0.5D0/(2*L+1) &
                 *(1.D0/(2*L+3)*( &
                                (pgbar**( L+2))*RM1G(NP,L)) &
                  -1.D0/(2*L-1)*( &
                                (pgbar**  L   )*RM3G(NP,L))) &
                 /pabbar

            IF(L.eq.0)THEN
               D1PSYG(NP,L)=-0.5D0/(2*L+1) &
                    *(1.D0/(2*L+3)*( (L+2)*(pgbar**( L+1))*RM1G(NP,L) &
                    ) &
                    -1.D0/(2*L-1))
               D1PHYG(NP,L)=-1.D0/(2*L+1)*pabbar**2
            ELSE
               D1PSYG(NP,L)=-0.5D0/(2*L+1) &
                    *(1.D0/(2*L+3)*( (L+2)*(pgbar**( L+1))*RM1G(NP,L) &
                    ) &
                    -1.D0/(2*L-1)*(  L   *(pgbar**( L-1))*RM3G(NP,L) &
                    ))
               D1PHYG(NP,L)=-1.D0/(2*L+1) &
                    *(  L   *(pgbar**( L-1))*RM1G(NP,L) &
                    ) &
                    *pabbar**2
            END IF
            IF(L.eq.0.or.L.eq.1)THEN
               D2PSYG(NP,L)=-0.5D0/(2*L+1) &
                    *(DBLE(L+1)*(L+2)/(2*L+3)*( &
                    (pgbar**  L   )*RM1G(NP,L)) &
                    ) &
                    *pabbar
            ELSE
               D2PSYG(NP,L)=-0.5D0/(2*L+1) &
                    *(DBLE(L+1)*(L+2)/(2*L+3)*( &
                    (pgbar**  L   )*RM1G(NP,L)) &
                    -DBLE(L  )*(L-1)/(2*L-1)*( &
                    (pgbar**( L-2))*RM3G(NP,L))) &
                    *pabbar
            END IF

         END IF
      END DO

      IF(MOD(IDBGFP/4,2).EQ.1) THEN
!
!     +++ plot of Phi, Psi and their derivatives +++
!
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
!
!     ----- calculation of local diffusion coefficienst -----
!   
      FACT=-4.D0*PI*RGAMH*1.D20

!      L0MIN=0

!      DO NP=1,NPMAX+1
      DO NP=NPSTART,NPENDWG
         RGAMA=SQRT(1.D0+PG(NP,NSBA)**2*THETA0(NSA))
         DO NTH=1,NTHMAX
            WA=0 
            WC=0
            DO L=LLMIN,LLMAX
               WA=WA+D2PSYG(NP,L)*PLM(NTH,L) 
               WC=WC+D1PHYG(NP,L)*PLM(NTH,L)
            END DO
            IF(NP.eq.1)THEN
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA) &
                    +FACT*WA*RGAMA**6
               FCPP2(NTH,NP,NR,NSB,NSA)=0.D0
            ELSE
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA) &
                    +FACT*WA*RGAMA**6
               FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA) &
                    +FACT*WC*RGAMA**3*AMFP(NSA)/AMFD(NSB)
            ENDIF
         END DO
      END DO

!      DO NP=1,NPMAX
      DO NP=NPSTART,NPENDWM
         RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
         DO NTH=1,NTHMAX+1
            WB=0 
            WD=0
            DO L=LLMIN,LLMAX
               WB=WB &
                 + 1.D0/ PM(NP,NSBA)    *D1PSYM(NP,L)*PLG(NTH,L)  *RGAMA**4 &
                 + 1.D0/(PM(NP,NSBA)**2)*PSYM(NP,L)  *D2PLG(NTH,L)*RGAMA**2
               WD=WD &
                 + 1.D0/ PM(NP,NSBA)    *PHYM(NP,L)  *D1PLG(NTH,L)
            END DO
            DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA) &
                           +FACT*WB
            FCTH2(NTH,NP,NR,NSB,NSA)=FCTH2(NTH,NP,NR,NSB,NSA) &
                           +FACT*WD*RGAMA*AMFP(NSA)/AMFD(NSB)
         END DO
      END DO
!
      IF(NPSTART.eq.1)THEN
         NP=1
         DO NTH=1,NTHMAX
            DCPT2(NTH,NP,NR,NSB,NSA)=0.D0
         END DO
      END IF

!      DO NP=2,NPMAX+1
      IF(NPSTART.eq.1)THEN
         NPS=2
      ELSE
         NPS=NPSTART
      END IF
      DO NP=NPS,NPENDWG
         RGAMA=SQRT(1.D0+PG(NP,NSBA)**2*THETA0(NSA))
         DO NTH=1,NTHMAX
            WE=0
            DO L=LLMIN,LLMAX
               WE=WE+( 1.D0/ PG(NP,NSBA) *D1PSYG(NP,L)*D1PLM(NTH,L)*RGAMA**4 &
                 -1.D0/(PG(NP,NSBA)**2)*PSYG(NP,L)  *D1PLM(NTH,L) )*RGAMA**2
            END DO
            DCPT2(NTH,NP,NR,NSB,NSA)=DCPT2(NTH,NP,NR,NSB,NSA) &
                           +FACT*WE
         END DO
      END DO

!      DO NP=1,NPMAX
      DO NP=NPSTART,NPENDWM
         RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
         DO NTH=1,NTHMAX+1
            WF=0
            DO L=LLMIN,LLMAX
               WF=WF+( 1.D0/ PM(NP,NSBA) *D1PSYM(NP,L)*D1PLG(NTH,L)*RGAMA**4 &
                 -1.D0/(PM(NP,NSBA)**2)*PSYM(NP,L)  *D1PLG(NTH,L) )*RGAMA**2
            END DO
            DCTP2(NTH,NP,NR,NSB,NSA)=DCTP2(NTH,NP,NR,NSB,NSA) &
                           +FACT*WF
         END DO
      END DO

      RETURN
      END SUBROUTINE FPCALC_NL
!-----------------------------------------
      SUBROUTINE FPCALC_NLAV(NR,NSA)
!
      IMPLICIT NONE
      integer,intent(in):: NR, NSA
      integer:: NSB, NTH, NP, NG
      real(8):: DELH, ETAL, X, PSIB, PCOS, ARG
      real(8):: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9
      real(8):: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
      INTEGER:: ISW_LAV, INTH
     
      ISW_LAV=0
! INTEGRATION OF BOUNCE AVERAGING
      DO NSB = 1,NSBMAX
!         DO NP=1,NPMAX+1
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               IF(NTH.ne.ITL(NR).and.NTH.ne.ITU(NR))THEN
                  DELH=2.D0*ETAM(NTH,NR)/NAVMAX
                  sum1=0.D0
                  sum2=0.D0
                  sum3=0.D0
                  DO NG=1,NAVMAX
                     ETAL=DELH*(NG-0.5D0)
                     X=EPSRM2(NR)*COS(ETAL)
                     PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
                     PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)

                     sum1=sum1 + DCPP2(NTH,NP,NR,NSB,NSA)*ABS(COSM(NTH))/PCOS
                     sum2=sum2 + DCPT2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                     sum3=sum3 + FCPP2(NTH,NP,NR,NSB,NSA)*ABS(COSM(NTH))/PCOS
                  END DO ! END NAVMAX 
                  
                  DCPP2(NTH,NP,NR,NSB,NSA)=SUM1*DELH * RCOEFNG(NR)
                  DCPT2(NTH,NP,NR,NSB,NSA)=SUM2*DELH * RCOEFNG(NR)
                  FCPP2(NTH,NP,NR,NSB,NSA)=SUM3*DELH * RCOEFNG(NR)
               END IF
            END DO ! END NTH
!            IF(ISW_LAV.eq.1)THEN
!               INTH=0
!               DO NTH=ITL(NR),ITL(NR)+1
!                  INTH=INTH+1
!                  DELH=2.D0*ETAM(NTH,NR)/NAVMAX
!                  sum7=0.D0
!                  sum8=0.D0
!                  sum9=0.D0
!                  
!                  temp7 = DCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp8 = FCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp9 = DCPT2B(INTH,NP,NR,NSB,NSA)
!                  
!                  IF (COSM(NTH).GE.0.D0) THEN
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX 
!                  ELSE ! SIGN OF PCOS
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX
!                  END IF
!                  
!                  DCPP2B(INTH,NP,NR,NSB,NSA)=SUM7*DELH
!                  FCPP2B(INTH,NP,NR,NSB,NSA)=SUM8*DELH
!                  DCPT2B(INTH,NP,NR,NSB,NSA)=SUM9*DELH
!               END DO ! END NTH
!               DO NTH=ITU(NR),ITU(NR)+1
!                  INTH=INTH+1
!                  DELH=2.D0*ETAM(NTH,NR)/NAVMAX
!                  sum7=0.D0
!                  sum8=0.D0
!                  sum9=0.D0
!                  
!                  temp7 = DCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp8 = FCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp9 = DCPT2B(INTH,NP,NR,NSB,NSA)
!                  
!                  IF (COSM(NTH).GE.0.D0) THEN
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX 
!                  ELSE ! SIGN OF PCOS
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX
!                  END IF
!
!                  DCPP2B(INTH,NP,NR,NSB,NSA)=SUM7*DELH
!                  FCPP2B(INTH,NP,NR,NSB,NSA)=SUM8*DELH
!                  DCPT2B(INTH,NP,NR,NSB,NSA)=SUM9*DELH
!               END DO ! END NTH
!            END IF ! NEW GRID ISW
         END DO ! END NP
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               IF(NTH.NE.NTHMAX/2+1) THEN
                  DELH=2.D0*ETAG(NTH,NR)/NAVMAX
                  sum4=0.D0
                  sum5=0.D0
                  sum6=0.D0
                  DO NG=1,NAVMAX
                     ETAL=DELH*(NG-0.5D0)
                     X=EPSRM2(NR)*COS(ETAL)
                     PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
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
                     sum4=sum4 + DCTP2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                     sum5=sum5 + DCTT2(NTH,NP,NR,NSB,NSA)*PCOS/(PSIB*COSG(NTH))
                     sum6=sum6 + FCTH2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                  END DO ! END NAVMAX
                  DCTP2(NTH,NP,NR,NSB,NSA)=sum4*DELH * RCOEFNG(NR)
                  DCTT2(NTH,NP,NR,NSB,NSA)=sum5*DELH * RCOEFNG(NR)
                  FCTH2(NTH,NP,NR,NSB,NSA)=sum6*DELH * RCOEFNG(NR)
               ELSE
                  DCTP2(NTH,NP,NR,NSB,NSA)=0.D0
                  DCTT2(NTH,NP,NR,NSB,NSA)=0.D0!?
                  FCTH2(NTH,NP,NR,NSB,NSA)=0.D0
               ENDIF ! END NTH!=pi/2
            END DO ! END NTH
         END DO ! END NP
      END DO ! END NSB
! END OF INTEGRATION

! BALANCE TRAPPED REGION for P direction
      DO NSB=1,NSBMAX
         DO NP=NPSTART,NPENDWG
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCPP2(NTH,NP,NR,NSB,NSA) &
                    =(DCPP2(NTH,NP,NR,NSB,NSA) &
                    +DCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0
               DCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA) &
                    =DCPP2(NTH,NP,NR,NSB,NSA)
            END DO ! END NTH
            IF(ISW_LAV.ne.1)THEN
               DCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
                    *( DCPP2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +DCPP2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
                    +DCPP2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
                    +DCPP2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
            ELSE
!               DCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0 &
!                    *( DCPP2B(1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
!                    +DCPP2B(2,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR)   &
!                    +DCPP2B(3,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR)   &
!                    +DCPP2B(4,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
            END IF
         END DO ! END NP
      END DO ! END NSB
      DO NSB=1,NSBMAX
!         DO NP=1,NPMAX+1
         DO NP=NPSTART,NPENDWG
            DO NTH=ITL(NR)+1,NTHMAX/2
               FCPP2(NTH,NP,NR,NSB,NSA) &
                    =(FCPP2(NTH,NP,NR,NSB,NSA) &
                    +FCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0
               FCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA) &
                    =FCPP2(NTH,NP,NR,NSB,NSA)
            END DO ! END NTH
            IF(ISW_LAV.ne.1)THEN
               FCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
                    *( FCPP2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +FCPP2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
                    +FCPP2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
                    +FCPP2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
            ELSE
!               FCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
!                    *( FCPP2B(1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
!                    +FCPP2B(2,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
!                    +FCPP2B(3,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
!                    +FCPP2B(4,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
            END IF
         END DO ! END NP
      END DO ! END NSB
      DO NSB=1,NSBMAX
!         DO NP=1,NPMAX+1
         DO NP=NPSTART,NPENDWG
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCPT2(NTH,NP,NR,NSB,NSA) &
                    =(DCPT2(NTH,NP,NR,NSB,NSA) &
!                     +DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0 ! symmetry
                     -DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0 ! assymmetry
               DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA) &
!                    =DCPT2(NTH,NP,NR,NSB,NSA) ! symmetry
                    =-DCPT2(NTH,NP,NR,NSB,NSA) ! assymmetry
            END DO ! END NTH
            IF(ISW_LAV.ne.1)THEN
               DCPT2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
                    *( DCPT2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +DCPT2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
!                    +DCPT2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) & ! symmetry
!                    +DCPT2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
                    -DCPT2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) & ! assymetry
                    -DCPT2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               DCPT2(ITU(NR),NP,NR,NSB,NSA)=DCPT2(ITL(NR),NP,NR,NSB,NSA) ! symmetry
               DCPT2(ITU(NR),NP,NR,NSB,NSA)=-DCPT2(ITL(NR),NP,NR,NSB,NSA) ! assymmetry
            ELSE
!               DCPT2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
!                    *( DCPT2B(1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
!                    +DCPT2B(2,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
!                    +DCPT2B(3,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
!                    +DCPT2B(4,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               DCPT2(ITU(NR),NP,NR,NSB,NSA)=-DCPT2(ITL(NR),NP,NR,NSB,NSA)
            END IF
         END DO ! END NP
      END DO ! END NSB
! END OF BALANCE TRAPPED REGION for P direction

! BALANCE TRAPPED REGION for THETA direction
      DO NSB=1,NSBMAX
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCTT2(NTH,NP,NR,NSB,NSA)        &
                    =(DCTT2(NTH,NP,NR,NSB,NSA) &
                    +DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
                    =DCTT2(NTH,NP,NR,NSB,NSA)
            END DO
         END DO
      END DO
      DO NSB=1,NSBMAX
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=ITL(NR)+1,NTHMAX/2
               FCTH2(NTH,NP,NR,NSB,NSA)        &
                    =(FCTH2(NTH,NP,NR,NSB,NSA) &
!                     +FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0 ! symmetry
                     -FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0 ! a
               FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
!                    =FCTH2(NTH,NP,NR,NSB,NSA) ! symmetry
                    =-FCTH2(NTH,NP,NR,NSB,NSA)
            END DO
         END DO
      END DO
      DO NSB=1,NSBMAX
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCTP2(NTH,NP,NR,NSB,NSA)        &
                    =(DCTP2(NTH,NP,NR,NSB,NSA) &
!                     +DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0 ! symmetry
                     -DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
!                    =DCTP2(NTH,NP,NR,NSB,NSA) ! symmetry
                    =-DCTP2(NTH,NP,NR,NSB,NSA)
            END DO
         END DO
      END DO
! END OF BALANCE TRAPPED REGION for THETA direction

      
      RETURN
      END SUBROUTINE FPCALC_NLAV

!-----------------------

      SUBROUTINE INTEGRATION_BACKGROUND_F(FPL,NR,NSB,NSA,RM1M,RM2M,RM3M,RM4M,RM1G,RM2G,RM3G,RM4G)

      IMPLICIT NONE

      integer,parameter::LNM=5
!      real(8),dimension(NTHMAX+3,-1:LNM):: PLM,PLG
!      real(8),DIMENSION(NTHMAX+3,-1:LNM):: D1PLM, D1PLG, D2PLG
!      real(8),DIMENSION(0:LNM):: PLTEMP
      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(IN):: FPL
      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(OUT):: RM1M, RM2M, RM3M, RM4M
      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(OUT):: RM1G, RM2G, RM3G, RM4G
      real(8),DIMENSION(NPMAX+3):: TX1, TY1, DF1, UTY10
      real(8),dimension(4,NPMAX+3):: UTY1

      integer:: NR, NSB, NSA, NTH, L,  NP, NNP, NPG, NSBA
      integer:: N,M
      real(8):: RGAMH, FACT, WA, WC, WB, WD, WE, WF
      real(8):: SUM1, PSUM, SUM2, SUM3, SUM4, SUM5, SUM6, SUM7, SUM8, SUM9
      real(8):: RGAMA, RGAMB, vtatb, pabbar, ptatb, PCRIT, pabar, pbbar, pgbar, pmbar
      integer:: LLMIN, IER

!
!     ----- definition of local quantities -----
!
      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
           /(RNFP0(NSA)*PTFP0(NSA)*1.D20)*RNFD0(NSB)
      vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
      pabbar=(PTFP0(NSA)*AMFD(NSB))/(PTFD0(NSB)*AMFP(NSA))
      NSBA=NSB_NSA(NSA)
!
!     ----- calculation of \hat{M}_l -----
!
      DO L=0,LLMAX
!         IF(L.lt.2)THEN
            TX1(1)=0.D0
!         ELSE
!            TX1(1)=PM(1,NSB)
!         ENDIF
         IF(L.eq.0)THEN
            TY1(1)=0.D0
         ELSEIF(L.eq.1)THEN
            TY1(1)=FPL(1,L)
         ELSE
            TY1(1)=FPL(1,L)*PM(1,NSB)**(1-L)
         ENDIF
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            TX1(NNP+1)=PM(NNP,NSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSB)**(1-L))*RGAMB**(1+L)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
         TY1(NPMAX+2)=FPL(NPMAX,L)*(PM(NPMAX,NSB)**(1-L))
!         TY1(NPMAX+2)=0.D0
         IF(L.eq.0)THEN
!            DF1(1) = FPL(1,L)
            DF1(1) = FPL(1,L) - (FPL(2,L)-FPL(1,L))/(PM(2,NSB)-PM(1,NSB))*PM(1,NSB)
         ELSEIF(L.eq.1)THEN
            DF1(1) = 0.D0
         ELSE
            DF1(1) = (1-L)*PM(1,NSB)**(-L)*FPL(1,NSB)
         ENDIF
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM1M(NP,L)=PSUM-SUM2
            ELSE
               RM1M(NP,L)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM1G(NPG,L)=PSUM-SUM3
            ELSE
               RM1G(NPG,L)=0.D0
            ENDIF
         END DO
      END DO

!
!     ----- calculation of \hat{N}_l -----
!
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            TX1(NNP+1)=PM(NNP,NSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSB)**(2+L))*RGAMB**(-L)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
!         TY1(NPMAX+2)=0.D0
         TY1(NPMAX+2)=FPL(NPMAX,L)*(PM(NPMAX,NSB)**(2+L))
         DF1(1)   = 0.D0
         DF1(NPMAX+2) = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM2M(NP,L)=SUM4
            ELSE
               RM2M(NP,L)=PSUM
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM2G(NPG,L)=SUM5
            ELSE
               RM2G(NPG,L)=PSUM
            ENDIF
         END DO
      END DO

!
!     ----- calculation of \hat{M}_l^+ -----
!
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            TX1(NNP+1)=PM(NNP,NSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSB)**(3-L))*RGAMB**(L-1)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
!         TY1(NPMAX+2)=0.D0
         TY1(NPMAX+2)=FPL(NPMAX,L)*(PM(NPMAX,NSB)**(3-L))
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM6,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM3M(NP,L)=PSUM-SUM6
            else
               RM3M(NP,L)=0.d0
            endif
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM7,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM3G(NPG,L)=PSUM-SUM7
            else
               RM3G(NPG,L)=0.d0
            endif
         END DO
      END DO

!
!     ----- calculation of \hat{N}_l^+ -----
!
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            TX1(NNP+1)=PM(NNP,NSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSB)**(4+L))*RGAMB**(-L-2)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
!         TY1(NPMAX+2)=0.D0
         TY1(NPMAX+2)=FPL(NPMAX,L)*(PM(NPMAX,NSB)**(4+L))
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM8,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM4M(NP,L)=SUM8
            else
               RM4M(NP,L)=PSUM
            endif
         END DO
         DO NPG=1,NPMAX+1
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM9,TX1,UTY1,UTY10,NPMAX+2,IER)
               RM4G(NPG,L)=SUM9
            else
               RM4G(NPG,L)=PSUM
            endif
         END DO
      END DO

      END SUBROUTINE INTEGRATION_BACKGROUND_F
!-------------------------------------------

      SUBROUTINE INTEGRATION_BACKGROUND_F_FINE(FPL,NR,NSB,NSA,RM1M,RM2M,RM3M,RM4M,RM1G,RM2G,RM3G,RM4G)

      IMPLICIT NONE

      integer,parameter::LNM=5
!      real(8),dimension(NTHMAX+3,-1:LNM):: PLM,PLG
!      real(8),DIMENSION(NTHMAX+3,-1:LNM):: D1PLM, D1PLG, D2PLG
!      real(8),DIMENSION(0:LNM):: PLTEMP
      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(IN):: FPL
      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL0
!      double precision,dimension(-1:LNM),intent(in):: FPLS1
!      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(OUT):: RM1M, RM2M, RM3M, RM4M
!      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(OUT):: RM1G, RM2G, RM3G, RM4G
      real(8),DIMENSION(NPSTARTW:NPENDWM,-1:LNM), INTENT(OUT):: RM1M, RM2M, RM3M, RM4M
      real(8),DIMENSION(NPSTART :NPENDWG,-1:LNM), INTENT(OUT):: RM1G, RM2G, RM3G, RM4G
      real(8),DIMENSION(2*NPMAX+3):: TX1, TY1, DF1, UTY10
      real(8),dimension(4,2*NPMAX+3):: UTY1

      integer:: NR, NSB, NSA, NTH, L,  NP, NNP, NPG, NSBA
      integer:: N,M
      real(8):: RGAMH, FACT, WA, WC, WB, WD, WE, WF
      real(8):: SUM1, PSUM, SUM2, SUM3, SUM4, SUM5, SUM6, SUM7, SUM8, SUM9
      real(8):: RGAMA, RGAMB, vtatb, pabbar, ptatb, PCRIT, pabar, pbbar, pgbar, pmbar
      real(8):: testP, testF
      integer:: LLMIN, IER, NPF

!
!     ----- definition of local quantities -----
!
      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
           /(RNFP0(NSA)*PTFP0(NSA)*1.D20)*RNFD0(NSB)
      vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
      pabbar=(PTFP0(NSA)*AMFD(NSB))/(PTFD0(NSB)*AMFP(NSA))
      NSBA=NSB_NSA(NSA)

!    SET FPL0
      DO L=0,LLMAX
         FPL0(1,L)=0.5D0*( FPL(1,L)+FPL(2,L) & 
              + (FPL(1,L)-FPL(2,L))/( DELP(NSB)*(PM(1,NSB)+PM(2,NSB)) ) &
              *(PM(1,NSB)**2+PM(2,NSB)**2) )
!         FPL0(1,L)=FPLS1(L)
         TX1(1)=0.D0
         TY1(1)=FPL0(1,L)
         DO NP=1,NPMAX
            TX1(NP+1)=PM(NP,NSB)
            TY1(NP+1)=FPL(NP,L)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
         TY1(NPMAX+2)=0.D0
         DF1(1)=0.D0
         DF1(NPMAX+2)=0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         DO NP=1,NPMAX-1
            testP=PM(NP,NSB)/PMAX(NSB)
            CALL SPL1DF(testP,testF,TX1,UTY1,NPMAX+2,IER)
            FPL0(NP+1,L)=testF
         END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     ----- calculation of \hat{M}_l -----
!
      DO L=0,LLMAX
         IF(L.lt.2)THEN
            TX1(1)=0.D0
!         ELSE
!            TX1(1)=PM(1,NSB)
         ENDIF
         IF(L.eq.0)THEN
            TY1(1)=0.D0
         ELSEIF(L.eq.1)THEN
            TY1(1)=FPL0(1,L)
         ELSE
            TY1(1)=FPL0(1,L)*PM(1,NSB)**(1-L)
         ENDIF
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(1-L))*RGAMB**(1+L)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            IF(PM(NNP,NSB).ge.1.D0)THEN
               RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
               NPF=NPF+1
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(1-L))*RGAMB**(1+L)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         IF(L.eq.0)THEN
            DF1(1) = FPL0(1,L)
         ELSEIF(L.eq.1)THEN
            DF1(1) = 0.D0
         ELSE
            DF1(1) = (1-L)*PM(1,NSB)**(-L)*FPL0(1,NSB)
         ENDIF
         DF1(NPF+1)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPF+1,IER)
               RM1M(NP,L)=PSUM-SUM2
            ELSE
               RM1M(NP,L)=0.D0
            ENDIF
         END DO
!         DO NPG=1,NPMAX+1
         DO NPG=NPSTART,NPENDWG
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPF+1,IER)
               RM1G(NPG,L)=PSUM-SUM3
            ELSE
               RM1G(NPG,L)=0.D0
            ENDIF
         END DO
      END DO

!
!     ----- calculation of \hat{N}_l -----
!
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+L))*RGAMB**(-L)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(2+L))*RGAMB**(-L)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         DF1(1)   = 0.D0
         DF1(NPF+1)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPF+1,IER)
               RM2M(NP,L)=SUM4
            ELSE
               RM2M(NP,L)=PSUM
            ENDIF
         END DO
!         DO NPG=1,NPMAX+1
         DO NPG=NPSTART,NPENDWG
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPF+1,IER)
               RM2G(NPG,L)=SUM5
            ELSE
               RM2G(NPG,L)=PSUM
            ENDIF
         END DO
      END DO

!
!     ----- calculation of \hat{M}_l^+ -----
!
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(3-L))*RGAMB**(L-1)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(3-L))*RGAMB**(L-1)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         DF1(1)   = 0.D0
         DF1(NPF+1)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM6,TX1,UTY1,UTY10,NPF+1,IER)
               RM3M(NP,L)=PSUM-SUM6
            else
               RM3M(NP,L)=0.d0
            endif
         END DO
!         DO NPG=1,NPMAX+1
         DO NPG=NPSTART,NPENDWG
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM7,TX1,UTY1,UTY10,NPF+1,IER)
               RM3G(NPG,L)=PSUM-SUM7
            else
               RM3G(NPG,L)=0.d0
            endif
         END DO
      END DO

!
!     ----- calculation of \hat{N}_l^+ -----
!
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(4+L))*RGAMB**(-L-2)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(4+L))*RGAMB**(-L-2)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         DF1(1)   = 0.D0
         DF1(NPF+1)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX(NSB),PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            ptatb=PM(NP,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM8,TX1,UTY1,UTY10,NPF+1,IER)
               RM4M(NP,L)=SUM8
            else
               RM4M(NP,L)=PSUM
            endif
         END DO
!         DO NPG=1,NPMAX+1
         DO NPG=NPSTART,NPENDWG
            RGAMA=SQRT(1.D0+PG(NPG,NSBA)**2*THETA0(NSA))
            ptatb=PG(NPG,NSBA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-THETA0(NSB)*vtatb**2*ptatb**2))*ptatb
            IF(PCRIT.le.PMAX(NSB)) then
               CALL SPL1DI(PCRIT,SUM9,TX1,UTY1,UTY10,NPF+1,IER)
               RM4G(NPG,L)=SUM9
            else
               RM4G(NPG,L)=PSUM
            endif
         END DO
      END DO

      END SUBROUTINE INTEGRATION_BACKGROUND_F_FINE
!-------------------------------------------
! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL_calcn(PML,NR,NS)

      USE plprof
      USE libbes,ONLY: besekn 
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_calcn

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.eq.0)THEN
         RL=0.D0
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NRSTART-1) THEN
         RL=RM(NRSTART)-DELR
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NREND+1.and.NR.ne.NRMAX+1) THEN
         RL=RM(NREND)+DELR
         RHON=MIN(RL,1.D0)
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NREND)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF
      CALL PL_PROF(RHON,PLF)
      IF(NT_init.eq.0)THEN
         RNFDL=PLF(NS)%RN/RNFD0L
         RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
      ELSE
         RNFDL=PLF(NS)%RN/RNFD0L
!         RNFDL=PNT(NR,NS,NTG2)
         RTFDL=RT_IMPL(NR,NS)
      END IF

      FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
      EX=PML**2/(2.D0*RTFDL/RTFD0L)
      FPMXWL_calcn=FACT*EXP(-EX)

      RETURN
      END FUNCTION FPMXWL_calcn

!      SUBROUTINE integration_background_f_trapezoid(FPL,NR,NSB,NSA,RM1M,RM2M,RM3M,RM4M,RM1G,RM2G,RM3G,RM4G)
!
!      IMPLICIT NONE
!
!      integer,parameter::LNM=5
!      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(IN):: FPL
!      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(OUT):: RM1M, RM2M, RM3M, RM4M
!      real(8),DIMENSION(NPMAX+3,-1:LNM), INTENT(OUT):: RM1G, RM2G, RM3G, RM4G
!
!      integer:: NR,NSA,NSB,NP,NSBA,NPB
!      real(8):: RGAMH, vtatb, pcrit
!!
!!     ----- definition of local quantities -----
!!
!      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
!           /(RNFP0(NSA)*PTFP0(NSA)*1.D20)*RNFD0(NSB)
!      vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
!      NSBA=NSB_NSA(NSA)
!      
!!
!!     ----- calculation of \hat{M}_l -----
!!
!      DO L=0,LLMAX
!         DO NP=1,NPMAX
!            pcrit=vtatb*PM(NP,NSBA)
!
!            NPB=1
!            DO WHILE(NPB.ne.0)
!               IF(pcrit.le.PM(NPB,NSB))THEN
!               END IF
!               NPB=NPB+1
!            END DO
!         END DO
!
!      END DO
!
!      END SUBROUTINE integration_background_f_trapezoid
!-----------------------------------------------------
      END MODULE fpcalcn



