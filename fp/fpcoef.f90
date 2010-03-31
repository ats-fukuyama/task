!     $Id$
!
! ************************************************************
!
!                      CALCULATION OF D AND F
!
! ************************************************************
!
      MODULE fpcoef

      USE fpcomm
      USE fpcalc
      USE fpcalw
      USE fpcalwm
      USE fpcalwr

      interface
         DOUBLE PRECISION FUNCTION BESEKN(N,X)
           real(8) :: X
           integer :: N
         end function BESEKN
      end interface 

      contains

      SUBROUTINE FP_COEF(NSA)

      IMPLICIT NONE
      integer:: NSA, NR, NTH, NP, NS
      real(8):: FPMAX

      ISAVE=0
      NS=NS_NSA(NSA)

      DO NR=NRSTART,NREND
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO

!
!     ----- Parallel electric field accleration term -----
!
      CALL FP_CALE(NSA)
!
!     ----- Quasi-linear wave-particle interaction term -----
!

      IF(MODELW(NS).EQ.0) THEN
         CALL FP_CALW(NSA)
      ELSEIF(MODELW(NS).EQ.1) THEN
         CALL FP_CALWR(NSA)
      ELSEIF(MODELW(NS).EQ.2) THEN
         CALL FP_CALWR(NSA)
      ELSEIF(MODELW(NS).EQ.3) THEN
         CALL FP_CALWM(NSA)
      ELSEIF(MODELW(NS).EQ.4) THEN
!         MODELA=0
         CALL FP_CALWM(NSA)
!         MODELA=1
      ELSE
         IF(nrank.eq.0) WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW(NS)
      ENDIF

!     ----- Collisional slowing down and diffusion term -----

      CALL FP_CALC(NSA)

!     ----- Sum up velocity diffusion terms -----

      DO NR=NRSTART,NREND
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO

         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=DCPP(NTH,NP,NR,NSA)+DWPP(NTH,NP,NR,NSA)
            DPT(NTH,NP,NR,NSA)=DCPT(NTH,NP,NR,NSA)+DWPT(NTH,NP,NR,NSA)
            FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
!
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=DCTP(NTH,NP,NR,NSA)+DWTP(NTH,NP,NR,NSA)
            DTT(NTH,NP,NR,NSA)=DCTT(NTH,NP,NR,NSA)+DWTT(NTH,NP,NR,NSA)
            FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
      ENDDO

!      DO NR=NRSTART, NREND
!      IF(NR.eq.8.and.NSA.eq.1)THEN
!         open(8,file='Dw_R8_N1_a1.dat')
!         DO NTH=1,NTHMAX+1
!            DO NP=1,NPMAX
!               WRITE(8,'(3I4,1P6E12.4)')NR,NP,NTH &
!                    ,PG(NP,NSA)*COSM(NTH), PG(NP,NSA)*SINM(NTH) &
!                    ,DWPP(NTH,NP,NR,NSA) &
!                    ,DWPT(NTH,NP,NR,NSA),DWTT(NTH,NP,NR,NSA),DWTP(NTH,NP,NR,NSA)
!            END DO
!         END DO
!         close(8)
!      END IF
!      END DO


!     ----- Radial diffusion term -----

      IF(MODELD.GT.0) CALL FP_CALR(NSA)
!      IF(NSA.eq.1.and.NR.eq.6)THEN
!         open(8,file='DRR_a1_R6.dat')
!         DO NTH=1,NTHMAX
!            DO NP=1,NPMAX
!               WRITE(8,'(3I4,1P6E12.4)')NR,NP,NTH &
!                    ,PG(NP,NSA)*COSM(NTH), PG(NP,NSA)*SINM(NTH) &
!                    ,DRR(NTH,NP,NR,NSA)
!            END DO
!            WRITE(8,*) " "
!            WRITE(8,*) " "
!         END DO
!         close(8)
!      END IF

!      IF(NSA.eq.1)THEN
!         DO NR=NRSTART,NREND+1
!            WRITE(*,'(I3,1P4E14.6)') NR, RG(NR),DRR(1,2,NR,NSA), RCOEF_GG(NR)
!         END DO
!      END IF

!     ----- Particle source term -----

      CALL FP_CALS(NSA)

!
!     ****************************
!     Boundary condition at p=pmax
!     ****************************
!
      DO NR=NRSTART,NREND
      DO NTH=1,NTHMAX
!         DPP(NTH,NPMAX,NR,NSA)=0.5D0*DPP(NTH,NPMAX,NR,NSA)
!         DPT(NTH,NPMAX,NR,NSA)=0.5D0*DPT(NTH,NPMAX,NR,NSA)
!         DPP(NTH,NPMAX+1,NR,NSA)=0.D0
!         DPT(NTH,NPMAX+1,NR,NSA)=0.D0
         FPMAX=FPP(NTH,NPMAX+1,NR,NSA)
         FPP(NTH,NPMAX+1,NR,NSA)=MAX(FPMAX,0.D0)
      ENDDO
      ENDDO
!
!      DO NR=NRSTART,NREND
!      DO NTH=1,NTHMAX+1
!         DTP(NTH,NPMAX,NR,NSA)=0.D0
!         DTT(NTH,NPMAX,NR,NSA)=0.D0
!         FTH(NTH,NPMAX,NR,NSA)=0.D0
!       ENDDO
!       ENDDO
!


      RETURN
      END SUBROUTINE fp_coef

! ****************************************
!     Parallel electric field
! ****************************************

      SUBROUTINE FP_CALE(NSA)

      IMPLICIT NONE
      integer:: NSA, NSB, NR, NTH, NP
      real(8):: PSP, SUML, ANGSP, SPL, FPMAX
      integer:: NG
      real(8):: FACT, DELH, sum11, ETAL, X, PSIB, PCOS, sum15, ARG

      DO NR=NRSTART,NREND
      DO NP=1,NPMAX+1
      DO NTH=1,NTHMAX
         FEPP(NTH,NP,NR,NSA)= AEFP(NSA)*E2(NR)/PTFP0(NSA)*COSM(NTH)
      ENDDO
      ENDDO
      ENDDO

      DO NR=NRSTART,NREND
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX+1
         FETH(NTH,NP,NR,NSA)=-AEFP(NSA)*E2(NR)/PTFP0(NSA)*SING(NTH)
      ENDDO
      ENDDO
      ENDDO

      IF (MODELA.EQ.0) RETURN

      DO NR=NRSTART,NREND
         IF(MODELA.eq.1)then
!         FACT=1.D0/SQRT(1.D0-EPSR(NR)**2)
            FACT=1.D0
            DO NP=1,NPMAX+1
               DO NTH=1,ITL(NR)-1
                  FEPP(NTH,NP,NR,NSA)= FACT*FEPP(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
            
            DO NP=1,NPMAX+1
               DO NTH=ITL(NR),ITU(NR)
                  FEPP(NTH,NP,NR,NSA)= 0.D0
               ENDDO
            ENDDO
            
            DO NP=1,NPMAX+1
               DO NTH=ITU(NR)+1,NTHMAX
                  FEPP(NTH,NP,NR,NSA)= FACT*FEPP(NTH,NP,NR,NSA)
               ENDDO
               FEPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
                    *( FEPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +FEPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                    +FEPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                    +FEPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR))
               
               FEPP(ITU(NR),NP,NR,NSA)=FEPP(ITL(NR),NP,NR,NSA)
            ENDDO
            
            DO NP=1,NPMAX
               DO NTH=1,ITL(NR)
                  FETH(NTH,NP,NR,NSA)=FACT*FETH(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
            
            DO NP=1,NPMAX
               DO NTH=ITL(NR)+1,ITU(NR)
                  FETH(NTH,NP,NR,NSA)= 0.D0
               ENDDO
            ENDDO
            
            DO NP=1,NPMAX
               DO NTH=ITU(NR)+1,NTHMAX+1
                  FETH(NTH,NP,NR,NSA)=FACT*FETH(NTH,NP,NR,NSA)
               ENDDO
            ENDDO
            

      END IF

      ENDDO

      RETURN
      END SUBROUTINE FP_CALE

! ****************************************
!     Radial transport
! ****************************************

      SUBROUTINE FP_CALR(NSA)

      IMPLICIT NONE
      integer:: NSA, NSBA, NS, NR, NTH, NP, NG
      real(8):: RHON, RTFPL, FACTR, FACTP, FACTRN, FACTRT, SV
      real(8):: PSIB, PCOS, X, ETAL, sumd, sumf, DELH
      real(8):: DNDR, NEDGE, FACT

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)
      DO NR=NRSTART,NREND+1
         RHON=RG(NR)
         IF(MODELD.EQ.2.OR.MODELD.EQ.4.OR.MODELD.eq.5) THEN
            SV=MAX(PNS(NS)/PN(NS),1.D-3)
            FACTRN=PROFN1*PROFN2*RHON**(PROFN1-1.D0)/((1-RHON**PROFN1)+SV)
            SV=MAX(PTS(NS)/PTPP(NS),1.D-3)
            FACTRT=PROFT1*PROFT2*RHON**(PROFT1-1.D0)/((1-RHON**PROFT1)+SV)

            IF(PROFN2.ge.1.D0)THEN
               NEDGE=(RNFP0(NSA)-RNFPS(NSA))*(1.D0-RG(NR)**PROFN1)**PROFN2+RNFPS(NSA)
               DNDR=-PROFN1*RHON**(PROFN1-1.D0)          &
                    *PROFN2*( RNFP0(NSA)-RNFPS(NSA) )    &
                    *( 1.D0-RHON**PROFN1 )**(PROFN2-1.D0)
            ELSE
               IF(NR.ne.NRMAX+1)THEN
                  NEDGE=(RNFP0(NSA)-RNFPS(NSA))*(1.D0-RG(NR)**PROFN1)**PROFN2+RNFPS(NSA)
                  DNDR=-PROFN1*RHON**(PROFN1-1.D0)          &
                       *PROFN2*( RNFP0(NSA)-RNFPS(NSA) )    &
                       *( 1.D0-RHON**PROFN1 )**(PROFN2-1.D0)
               ELSEIF(NR.eq.NRMAX+1)THEN
                  NEDGE=(RNFP0(NSA)-RNFPS(NSA))*(1.D0-RG(NRMAX)**PROFN1)**PROFN2+RNFPS(NSA)
                  DNDR=( -NEDGE+RNFPS(NSA))/DELR
               END IF
            END IF

!            IF(NSA.eq.1) write(*,*) NR, DNDR
         ENDIF
         IF(MODELD.EQ.2.OR.MODELD.EQ.3.OR.MODELD.EQ.4) THEN
            RTFPL=RTFP(NR,NSA)/RTFP0(NSA)
         ENDIF
      DO NP=1,NPMAX+1
         SELECT CASE(MODELD)
         CASE(1)
            FACTR=0.D0
            FACTP=1.D0
         CASE(2)
            FACTR=-FACTRN+(1.5D0-0.5D0*PM(NP,NSBA)**2/RTFPL)*FACTRT
            FACTP=1.D0
         CASE(3)
            FACTR=0.D0
            FACTP=1.D0/SQRT(1.D0+PM(NP,NSBA)**2/RTFPL)
         CASE(4)
            FACTR=-FACTRN+(1.5D0-0.5D0*PM(NP,NSBA)**2/RTFPL)*FACTRT
            FACTP=1.D0/SQRT(1.D0+PM(NP,NSBA)**2/RTFPL)
         CASE(5)
            FACTR=DNDR/NEDGE
            FACTP=1.D0
        END SELECT
      DO NTH=1,NTHMAX
         FACT= (DRR0-DRRS)*(1.D0-RHON**2)+DRRS 
!         FACT=1.D0
         DRR(NTH,NP,NR,NSA)= FACT &
              *FACTP      /(RA*RA)*RLAMDA_GG(NTH,NR)!/PTFP0(NSA)**2
         FRR(NTH,NP,NR,NSA)= FACT &
              *FACTP*FACTR/(RA*RA)*RLAMDA_GG(NTH,NR)!/PTFP0(NSA)
      ENDDO
      ENDDO

      IF(MODELA.eq.1)THEN! Bounce average for radial diffusion coef.
         DO NP=1,NPMAX+1
            DO NTH=ITLG(NR)+1,NTHMAX/2
               DRR(NTH,NP,NR,NSA) &
                    =(DRR(NTH,NP,NR,NSA) &
                    +DRR(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
               FRR(NTH,NP,NR,NSA) &
                    =(FRR(NTH,NP,NR,NSA) &
                    +FRR(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
               DRR(NTHMAX-NTH+1,NP,NR,NSA) &
                    =DRR(NTH,NP,NR,NSA)
               FRR(NTHMAX-NTH+1,NP,NR,NSA) &
                    =FRR(NTH,NP,NR,NSA)
            END DO
         END DO
         DO NP=1,NPMAX+1
            DRR(ITLG(NR),NP,NR,NSA) = RLAMDA_GG(ITLG(NR),NR)/4.D0 &
                 *( DRR(ITLG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)-1,NR) &
                 +DRR(ITLG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)+1,NR)   &
                 +DRR(ITUG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)-1,NR)   &
                 +DRR(ITUG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)+1,NR)  )
            FRR(ITLG(NR),NP,NR,NSA) = RLAMDA_GG(ITLG(NR),NR)/4.D0 &
                 *( FRR(ITLG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)-1,NR) &
                 +FRR(ITLG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITLG(NR)+1,NR)   &
                 +FRR(ITUG(NR)-1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)-1,NR)   &
                 +FRR(ITUG(NR)+1,NP,NR,NSA)/RLAMDA_GG(ITUG(NR)+1,NR)  )
            DRR(ITUG(NR),NP,NR,NSA)=DRR(ITLG(NR),NP,NR,NSA)
            FRR(ITUG(NR),NP,NR,NSA)=FRR(ITLG(NR),NP,NR,NSA)
         END DO
      END IF
      ENDDO ! end of bounce average


      IF(MODELA.eq.2)THEN 
      DO NR=NRSTART, NREND+1
         DO NTH=1,NTHMAX
            DELH=2.D0*ETAMG(NTH,NR)/NAVMAX
            DO NP=1,NPMAX+1
               sumd=0.D0
               sumf=0.D0
               DO NG=1,NAVMAX
                  ETAL=DELH*(NG-0.5D0)
                  X=EPSRG(NR)*COS(ETAL)*RR
                  PSIB=(1.D0+EPSRG(NR))/(1.D0+X/RR)
                  IF (COSM(NTH).GE.0.D0) THEN
                     PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
                  ELSE
                     PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
                  ENDIF
                  sumd=sumd &
                       +DRR(NTH,NP,NR,NSA)*COSM(NTH)/PCOS
                  sumf=sumf &
                       +FRR(NTH,NP,NR,NSA)*COSM(NTH)/PCOS
               END DO
               DRR(NTH,NP,NR,NSA)=SUMd*DELH/PI!*RCOEF(NR)
               FRR(NTH,NP,NR,NSA)=SUMf*DELH/PI!*RCOEF(NR)
            END DO
         END DO

         DO NP=1,NPMAX+1
            DO NTH=ITL(NR)+1,NTHMAX/2
               DRR(NTH,NP,NR,NSA) &
                    =(DRR(NTH,NP,NR,NSA) &
                    +DRR(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
               FRR(NTH,NP,NR,NSA) &
                    =(FRR(NTH,NP,NR,NSA) &
                    +FRR(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
               DRR(NTHMAX-NTH+1,NP,NR,NSA) &
                    =DRR(NTH,NP,NR,NSA)
               FRR(NTHMAX-NTH+1,NP,NR,NSA) &
                    =FRR(NTH,NP,NR,NSA)
            END DO
            DRR(ITL(NR),NP,NR,NSA)=RLAMDAG(ITL(NR),NR)/4.D0     &
                 *( DRR(ITL(NR)-1,NP,NR,NSA)/RLAMDAG(ITL(NR)-1,NR) &
                 +DRR(ITL(NR)+1,NP,NR,NSA)/RLAMDAG(ITL(NR)+1,NR) &
                 +DRR(ITU(NR)-1,NP,NR,NSA)/RLAMDAG(ITU(NR)-1,NR) &
                 +DRR(ITU(NR)+1,NP,NR,NSA)/RLAMDAG(ITU(NR)+1,NR))
            FRR(ITL(NR),NP,NR,NSA)=RLAMDAG(ITL(NR),NR)/4.D0     &
                 *( FRR(ITL(NR)-1,NP,NR,NSA)/RLAMDAG(ITL(NR)-1,NR) &
                 +FRR(ITL(NR)+1,NP,NR,NSA)/RLAMDAG(ITL(NR)+1,NR) &
                 +FRR(ITU(NR)-1,NP,NR,NSA)/RLAMDAG(ITU(NR)-1,NR) &
                 +FRR(ITU(NR)+1,NP,NR,NSA)/RLAMDAG(ITU(NR)+1,NR))
            DRR(ITU(NR),NP,NR,NSA)=DRR(ITL(NR),NP,NR,NSA)
            FRR(ITU(NR),NP,NR,NSA)=FRR(ITL(NR),NP,NR,NSA)
         END DO
      ENDDO
      END IF 

!      DO NR=NRSTART,NREND
!      IF(NR.eq.20.and.NSA.eq.1)THEN
!         open(8,file='drr_L20_r.dat')
!         DO NTH=1,NTHMAX
!            DO NP=1,NPMAX
!               WRITE(8,'(2I5,1P4E12.4)') NP,NTH,PG(NP,NSA)*COSM(NTH),PG(NP,NSA)*SINM(NTH) &
!                    ,DRR(NTH,NP,NR,NSA),RLAMDA(NTH,NR)
!            END DO
!            WRITE(8,*) " "
!            WRITE(8,*) " "
!         END DO
!         close(8)
!      END IF
!      END DO

      RETURN
      END SUBROUTINE FP_CALR

! ****************************************
!     Particle source and loss
! ****************************************

      SUBROUTINE FP_CALS(NSA)

      USE fpnfrr
      IMPLICIT NONE
      integer:: NSA, NSB, NSBA, NR, NTH, NP, NS, ID
      integer:: NBEAM, NSABEAM, NSAX
      real(8):: PSP, SUML, ANGSP, SPL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

!     ----- Particle source term -----

      DO NR=NRSTART,NREND
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               PPL(NTH,NP,NR,NSA)=0.D0
               SPPB(NTH,NP,NR,NSA)=0.D0
               SPPF(NTH,NP,NR,NSA)=0.D0
               SPPS(NTH,NP,NR,NSA)=0.D0
            ENDDO
         ENDDO
      ENDDO

!     ----- NBI source term -----

      DO NBEAM=1,NBEAMMAX
         IF(NSSPB(NBEAM).EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPBENG(NBEAM)*AEE)/PTFP0(NSA)
            ANGSP=PI*SPBANG(NBEAM)/180.D0
            SUML=0.D0
            DO NP=1,NPMAX-1
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NTH=1,NTHMAX
                     IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                        DO NR=1,NRMAX
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                           SUML=SUML &
                                +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)*RLAMDAG(NTH,NR)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
            DO NP=1,NPMAX-1
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NTH=1,NTHMAX
                     IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                        DO NR=NRSTART,NREND
                           SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                           SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                + SPBTOT(NBEAM)*SPL/SUML*RLAMDA(NTH,NR)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!  ----  FOR ELECTRON (NS=1)
      IF(NS_NSA(NSA).EQ.1) THEN
         DO NBEAM=1,NBEAMMAX
            NSABEAM=0
            DO NSAX=1,NSAMAX
               IF(NS_NSA(NSAX).EQ.NSSPB(NBEAM)) NSABEAM=NSAX
            ENDDO

            IF(NSABEAM.NE.0) THEN
               PSP=SQRT(2.D0*AMFP(NSA)**2*SPBENG(NBEAM)*AEE &
                    /AMFP(NSABEAM))/PTFP0(NSA)
               ANGSP=PI*SPBANG(NBEAM)/180.D0
               SUML=0.D0
               DO NP=0,NPMAX-1
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NTH=1,NTHMAX
                        IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                           DO NR=1,NRMAX
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SUML=SUML &
                                   +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)*RLAMDAG(NTH,NR)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               DO NP=0,NPMAX-1
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NTH=1,NTHMAX
                        IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                           DO NR=NRSTART,NREND
                              SPL=EXP(-(RM(NR)-SPBR0(NBEAM))**2/SPBRW(NBEAM)**2)
                              SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                                   +SPBTOT(NBEAM)*SPL/SUML*RLAMDA(NTH,NR)/RNFP0(NSABEAM)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

         END DO
      END IF

!     ----- Fixed fusion source term -----

      IF(MODELS.EQ.1) THEN
         IF(NSSPF.EQ.NS_NSA(NSA)) THEN
            PSP=SQRT(2.D0*AMFP(NSA)*SPFENG*AEE)/PTFP0(NSA)
            SUML=0.D0
            DO NP=1,NPMAX-1
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=1,NRMAX
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SUML=SUML &
                            +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)*RLAMDAG(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
            DO NP=1,NPMAX-1
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=NRSTART,NREND
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             + SPFTOT*SPL/SUML*RLAMDA(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         END IF

!     ----- electron -----

         IF(NS_NSA(NSA).EQ.1) THEN
            NSABEAM=0
            DO NSAX=1,NSAMAX
               IF(NS_NSA(NSAX).EQ.NSSPF) NSABEAM=NSAX
            ENDDO
            IF(NSABEAM.NE.0) THEN
            PSP=SQRT(2.D0*AMFP(NSA)**2*SPFENG*AEE &
                    /(AMFP(NSABEAM)))*PTFP0(NSA)
            SUML=0.D0
            DO NP=0,NPMAX-1
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=1,NRMAX
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SUML=SUML &
                            +SPL*VOLP(NTH,NP,NSBA)*VOLR(NR)*RLAMDAG(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            SUML=SUML*RNFP0(NSA)
            DO NP=0,NPMAX-1
               IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                  DO NR=NRSTART,NREND
                     SPL=EXP(-(RM(NR)-SPFR0)**2/SPFRW**2)
                     DO NTH=1,NTHMAX
                        SPPB(NTH,NP,NR,NSA)=SPPB(NTH,NP,NR,NSA) &
                             +SPFTOT*SPL/SUML*RLAMDA(NTH,NR)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            ENDIF
         END IF
      END IF

!     ----- Calcluated fusion source term -----

      IF(MODELS.EQ.2) THEN
         DO ID=1,6
!            IF(NCHECK.eq.0.or.NCHECK.gt.LMAXFP.and.NR.eq.1)THEN
!               WRITE(6,'(A,4I5,1PE12.4)') 'ID,NSA,NS,NSA1_NF,ENG1_NF=', &
!                 ID,NSA,NS,NSA1_NF(ID),ENG1_NF(ID)
!            END IF
            IF(NSA.EQ.NSA1_NF(ID)) THEN
               PSP=SQRT(2.D0*AMFP(NSA)*ENG1_NF(ID)*AEE)/PTFP0(NSA)
               DO NP=1,NPMAX-1
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
!                     IF(NCHECK.eq.0.or.NCHECK.gt.LMAXFP)THEN
!                        write(6,'(A,I5,1P3E12.4)') ' |-NP,PSP,PG=',&
!                          NP,PSP,PG(NP,NSBA),PG(NP+1,NSBA)
!                     END IF
                     DO NR=NRSTART,NREND
                        SUML=0.D0
                        DO NTH=1,NTHMAX
                           SUML=SUML+VOLP(NTH,NP,NSBA)*RLAMDA(NTH,NR)
                        ENDDO
                        CALL NF_REACTION_RATE(NR,ID)
                        DO NTH=1,NTHMAX
                           SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                               +RATE_NF(NR,ID)/SUML/RNFP0(NSA)*RLAMDA(NTH,NR)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               IF(PSP.ge.PG(NPMAX,NSBA))THEN
                  NP=NPMAX
                  IF(NCHECK.eq.0.or.NCHECK.gt.LMAXFP)THEN
                     write(6,'(A,I5,1P3E12.4)') '  |-NP,PSP,PG=',&
                       NP,PSP,PMAX(NSBA)
                     WRITE(6,*) ' |-  OUT OF RANGE PMAX'
                  END IF
                  DO NR=NRSTART,NREND
                     SUML=0.D0
                     DO NTH=1,NTHMAX
                        SUML=SUML+VOLP(NTH,NP,NSBA)*RLAMDA(NTH,NR)
                     ENDDO
                     CALL NF_REACTION_RATE(NR,ID)
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                             +RATE_NF(NR,ID)/SUML/RNFP0(NSA)*RLAMDA(NTH,NR)
                     ENDDO
                  ENDDO
               END IF
            ENDIF
            IF(NSA.EQ.NSA2_NF(ID)) THEN
               PSP=SQRT(2.D0*AMFP(NSA)*ENG2_NF(ID)*AEE)/PTFP0(NSA)
               DO NP=1,NPMAX-1
                  IF(PG(NP,NSBA).LE.PSP.AND.PG(NP+1,NSBA).GT.PSP) THEN
                     DO NR=NRSTART,NREND
                        SUML=0.D0
                        DO NTH=1,NTHMAX
                           SUML=SUML+VOLP(NTH,NP,NSBA)*RLAMDA(NTH,NR)
                        ENDDO
                        DO NTH=1,NTHMAX
                           SPPF(NTH,NP,NR,NSA)=SPPF(NTH,NP,NR,NSA) &
                                +RATE_NF(NR,ID)/SUML/RNFP0(NSA)*RLAMDA(NTH,NR)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
            IF(NSA.EQ.NSA1_NF(ID)) THEN
               DO NR=NRSTART,NREND
                  DO NP=1,NPMAX
                     DO NTH=1,NTHMAX
                        SPPF(NTH,NP,NR,NSB1_NF(ID))=                  &
                             SPPF(NTH,NP,NR,NSB1_NF(ID))              &
                             -RATE_NF_D1(NR,ID,NTH,NP)*RLAMDA(NTH,NR) &
                             /RNFP0(NSB1_NF(ID))
                        SPPF(NTH,NP,NR,NSB2_NF(ID))=                  &
                             SPPF(NTH,NP,NR,NSB2_NF(ID))              &
                             -RATE_NF_D2(NR,ID,NTH,NP)*RLAMDA(NTH,NR) &
                             /RNFP0(NSB2_NF(ID))
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!     ----- Particle loss and source terms -----
!

      IF(TLOSS(NS).EQ.0.D0) THEN
         DO NR=NRSTART,NREND
            DO NTH=1,NTHMAX
               DO NP=1,NPMAX-1
                  PPL(NTH,NP,NR,NSA)=0.D0
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO NR=NRSTART,NREND
            DO NTH=1,NTHMAX
               DO NP=1,NPMAX-1
                  PPL(NTH,NP,NR,NSA)=-RLAMDA(NTH,NR)/TLOSS(NS)
                  SPPS(NTH,NP,NR,NSA)=FPMXWL(PM(NP,NSBA),NR,NS) &
                                     /TLOSS(NS)*RLAMDA(NTH,NR)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE FP_CALS

! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      real(8):: FPMXWL

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.EQ.NRSTART-1) THEN
         RL=RM(NRSTART)-DELR
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NREND+1) THEN
         RL=RM(NREND)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF
      CALL PL_PROF(RHON,PLF)
      RNFDL=PLF(NS)%RN/RNFD0L
      RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0

      IF (MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(-EX)
         ENDIF

      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
            DKBSL=BESEKN(2,Z)
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
             *RTFD0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         IF(EX.LT.-100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(EX)
         ENDIF
      END IF

      RETURN
      END FUNCTION FPMXWL
!-------------------------------------------------------------
!
!-------------------------------------------------------------

      SUBROUTINE FPMXWL_EDGE(NP,NSA,FL)

      implicit none
      integer :: NP, NSA, NSBA, NS
      real(8):: FL1, FL2
      real(8),intent(out):: FL

      NS=NS_NSA(NSA)
      NSBA=NSB_NSA(NSA)

      FL1=FPMXWL(PM(NP,NSBA),NRMAX,NS)
      FL2=FPMXWL(PM(NP,NSBA),NRMAX+1,NS)

      FL=FL2+(FL2-FL1)

      RETURN
      END SUBROUTINE FPMXWL_EDGE


      END MODULE fpcoef
