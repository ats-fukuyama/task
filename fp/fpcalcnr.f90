!     $Id: fpcalcnr.f90,v 1.22 2013/02/08 07:36:24 nuga Exp $
!
! ************************************************************
!
!      CALCULATION OF NONLINEAR-RELAIVISTIC COLLISIONAL OPERATOR
!
! ************************************************************
!
      MODULE fpcalcnr

      USE fpcomm
      USE libspf, ONLY: dpleg
!      USE fpcoef, ONLY: FPMXWL
      real(rkind):: PNFP_NLR, THETA0L_NLR, THETAL_NLR


      contains

!-------------------------------------------------------------

      SUBROUTINE FPCALC_NLR(NR,NSB,NSA)

      USE libde,ONLY: DEHIFT
      USE libspl1d
      USE libgrf,ONLY: grd1d
      USE fpmpi
      Implicit none
!      PARAMETER (N=NPM+2,M=NTHM+2,LNM=5)
      integer,parameter::LNM=5
      real(rkind),DIMENSION(NTHMAX+3,-1:LNM):: PLM, PLG, D1PLM, D1PLG, D2PLG
      real(rkind),DIMENSION(0:LNM):: PLTEMP
      real(rkind),DIMENSION(NPSTART:NPEND):: FPLL
      real(rkind),DIMENSION(NPMAX):: FPL_recv
      real(rkind),DIMENSION(NPMAX,-1:LNM):: FPL
      double precision,dimension(-1:LNM):: FPLS1
      double precision:: FPLS1_temp
      integer:: NPS

      real(rkind),DIMENSION(NTHMAX+3):: TX,TY,DF
      real(rkind),DIMENSION(4,NTHMAX+3):: UTY
      real(rkind),dimension(NTHMAX+3)::UTY0
      real(rkind),DIMENSION(NPMAX+3):: TX1,TY1,DF1,UTY10
      real(rkind),DIMENSION(4,NPMAX+3):: UTY1

!!      real(rkind),DIMENSION(-2:LLMAX+2, -1:2):: FKLF_J,FKLF_Y
      real(rkind),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(rkind),DIMENSION(-2:LLMAX+2, 0:2):: DERJ, DERY

      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX, 0:2, -1:2):: RJABM
      real(rkind),DIMENSION(NPSTART:NPENDWG , 0:LLMAX, 0:2, -1:2):: RJABG 
      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX, 0:2, -1:2):: RYABM
      real(rkind),DIMENSION(NPSTART:NPENDWG , 0:LLMAX, 0:2, -1:2):: RYABG

      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX):: &
           DPSI02M,DPSI022M,PSI0M,PSI02M,PSI022M,PSI1M,PSI11M,DPSI11M

      real(rkind),DIMENSION(NPSTART:NPENDWG,  0:LLMAX):: &
           DPSI02G,DPSI022G,PSI0G,PSI02G,PSI022G,DPSI1G,DPSI11G


      integer:: NP, NTH, NSA, NSB, L, NR, LLMIN, NI, NA, NNP, NPG, NSSA, NSSB
      integer:: IER, LTEST, INTH
      real(rkind):: RGAMH, SUM1, SUM2, SUM3, SUM4, SUM5
      real(rkind):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2, RUFD
      real(rkind):: SUMA, SUMB, SUMC, SUMD, SUME, SUMF, SUMG, SUMH
      real(rkind):: vtatb, ptatb, PMAX2, RINT0, RINT2, ES0, ES2, testF, testP
      real(rkind):: DKBSL0, DKBSL1, DKBSL2, Z
!
!----- DEFINITION OF LOCAL QUANTITIES -------------
! 
      NSSA=NS_NSA(NSA)
      NSSB=NS_NSB(NSB)
      THETA0L_NLR=THETA0(NSSB)
      IF(MODEL_DISRUPT.eq.0)THEN 
         THETAL_NLR =(PTFD(NR,NSB)/(AMFD(NSB)*VC))**2 
      ELSE
         THETAL_NLR =THETA0(NSSB)*RT_quench(NR)/RTFD0(NSB)
      END IF
      FACT=AEFP(NSA)**2*AEFD(NSB)**2*LNLAM(NR,NSB,NSA)/(4.D0*PI*EPS0**2)
!      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
!           /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
      RGAMH=FACT*AMFP(NSA)/PTFP0(NSA)**3
!
!   --- calculation of Legendre Polynomials and its derivatives---
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
         CALL GRD1D(1,thm,plm,NTHMAX+3,NTHMAX,LLMAX+2,'@PLM:@',0)
         CALL GRD1D(2,thm,d1plm,NTHMAX+3,NTHMAX,LLMAX+2,'@D1PLM:@',0)
         CALL PAGEE

         CALL PAGES
         CALL GRD1D(1,thg,plg,NTHMAX+3,NTHMAX+1,LLMAX+2,'@PLG:@',0)
         CALL GRD1D(2,thg,d1plg,NTHMAX+3,NTHMAX+1,LLMAX+2,'@D1PLG:@',0)
         CALL GRD1D(3,thg,d2plg,NTHMAX+3,NTHMAX+1,LLMAX+2,'@D2PLG:@',0)
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
            DF(1)= FNSB(1,NP,NR,NSB)
            DF(NTHMAX+2)= (-1)**(L+1)*FNSB(NTHMAX,NP,NR,NSB)
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPLL(NP)=0.5D0*(2.D0*L+1.D0)*SUM1
!            WRITE(8,*) L, NP, FPLL(NP,L)
         END DO
         CALL fpl_comm(FPLL,FPL_recv)
         DO NP=1,NPMAX
            FPL(NP,L)=FPL_recv(NP)
         END DO
!
         TX(1)=0.D0
         TY(1)=0.D0
         FPLS1_temp=FPMXWL_calcnr(0.D0,NR,NSSB)
         DO NTH=1,NTHMAX
            TX(NTH+1)=THM(NTH)
            TY(NTH+1)=FPLS1_temp*PLM(NTH,L)*SINM(NTH)
         END DO
         TX(NTHMAX+2)=PI
         TY(NTHMAX+2)=0.D0
         DF(1)= FPLS1_temp
         DF(NTHMAX+2)= (-1)**(L+1)*FPLS1_temp ! satisfy until l=2
         CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
         CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
         CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
         FPLS1(L)=0.5D0*(2.D0*L+1.D0)*SUM1
!            WRITE(8,*) L, NP, FPL(NP,L)
      END DO ! LLMAX
      CALL mtx_reset_communicator
!      close(8)
!
!---- INTEGRAL ABBREVIATIONS
!
!      CALL INTEGRATION_RJAB_RYAB(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)
      CALL INTEGRATION_RJAB_RYAB_FINE(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM,FPLS1)

!---- END OF INTEGRALS
!
!
!---- PSI AND ITS DERIVATIONS
!
!
!---------- FOR MID OF GRID
!
!      IF(NSA.eq.1.and.NSB.eq.2) THEN
!         open(8,file='RJYM.dat')
!      END IF
      DO L = 0,LLMAX
         DO NP = NPSTARTW, NPENDWM
            RGAMA=SQRT(1.D0+PM(NP,NSSA)**2*THETA0(NSSA))
            RUFP = (PTFP0(NSA)*PM(NP,NSSA))/AMFP(NSA)

            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            CALL DERIV_JY(RUFP,RJ_1,RY_1,DERJ,DERY)

            DPSI02M(NP,L) = 0.5D0*(                               &
                 DERY(L,0)*RJABM(NP,L,1,1)                        &
                 -(RY_1(L-1,1)+RUFP*DERY(L-1,1) )*RJABM(NP,L,0,2) &
                 +(RJ_1(L+1,1)+RUFP*DERJ(L+1,1) )*RYABM(NP,L,0,0) &
                 -DERJ(L,2)*RYABM(NP,L,1,1) )

            DPSI022M(NP,L) = VC/8.D0*(                            &
                 DERY(L,0)*RJABM(NP,L,2,0)                        &
                 -2.D0*(RY_1(L-1,1)+RUFP*DERY(L-1,1) )            &
                 *(RJABM(NP,L,1,1)+RJABM(NP,L,2,0)/VC)            &
                 +(2.D0*RUFP*RY_1(L-2,0)+RUFP**2*DERY(L-2,0) )    &
                 *RJABM(NP,L,0,2)                                 &
                 +(2.D0*RUFP*RJ_1(L+2,0)+RUFP**2*DERJ(L+2,0) )    &
                 *RYABM(NP,L,0,0)                                 &
                 -2.D0*(RJ_1(L+1,1)+2.D0*RUFP/VC*RJ_1(L+2,0)      &
                 +RUFP*( DERJ(L+1,1)+RUFP/VC*DERJ(L+2,0) ) )      &
                 *RYABM(NP,L,1,1)                                 &
                 +DERJ(L,2)*RYABM(NP,L,2,0) )


            PSI02M(NP,L) = 0.5D0*(                       &
                   RY_1(L,0)*RJABM(NP,L,1,1)             &
                 - RUFP*RY_1(L-1,1)*RJABM(NP,L,0,2)      &
                 + RUFP*RJ_1(L+1,1)*RYABM(NP,L,0,0)      &
                 - RJ_1(L,2)*RYABM(NP,L,1,1)  )


            PSI022M(NP,L) = VC/8.D0*(                                 &
                 RY_1(L,0)*RJABM(NP,L,2,0)                            &
                 -2.D0*RUFP*RY_1(L-1,1)                               &
                 *(RJABM(NP,L,1,1)+RJABM(NP,L,2,0)/VC)                &
                 +RUFP**2*RY_1(L-2,0)*RJABM(NP,L,0,2)                 &
                 +RUFP**2*RJ_1(L+2,0)*RYABM(NP,L,0,0)                 &
                 -2.D0*RUFP                                           &
                 *(RJ_1(L+1,1)+RUFP/VC*RJ_1(L+2,0) )*RYABM(NP,L,1,1)  &
                 +RJ_1(L,2)*RYABM(NP,L,2,0) )

            PSI11M(NP,L) = 0.5D0 * (                &
                 RY_1(L,1)*RJABM(NP,L,1,0)          &
                 -RUFP*RY_1(L-1,0)*RJABM(NP,L,0,1)  &
                 +RUFP*RJ_1(L+1,0)*RYABM(NP,L,0,1)  &
                 -RJ_1(L,1)*RYABM(NP,L,1,0)         &
                 )

            PSI1M(NP,L) = ( RY_1(L,1)*RJABM(NP,L,0,1) &
                 +RJ_1(L,1)*RYABM(NP,L,0,1) )/VC

         END DO
      END DO

!
!----------- END OF MID 
!
!
!----------- ON GRID
!
      IF(NPSTART.eq.1)THEN
         NPS=2
      ELSE
         NPS=NPSTART
      END IF
      DO L=0,LLMAX
         DO NP = NPS, NPENDWG
            RGAMA=SQRT(1.D0+PG(NP,NSSA)**2*THETA0(NSSA))
            RUFP = (PTFP0(NSA)*PG(NP,NSSA))/AMFP(NSA)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            CALL DERIV_JY(RUFP,RJ_1,RY_1,DERJ,DERY)

            DPSI02G(NP,L) = 0.5D0*(                               &
                 DERY(L,0)*RJABG(NP,L,1,1)                        &
                 -(RY_1(L-1,1)+RUFP*DERY(L-1,1) )*RJABG(NP,L,0,2) &
                 +(RJ_1(L+1,1)+RUFP*DERJ(L+1,1) )*RYABG(NP,L,0,0) &
                 -DERJ(L,2)*RYABG(NP,L,1,1) )

            DPSI022G(NP,L) = VC/8.D0*(                            &
                 DERY(L,0)*RJABG(NP,L,2,0)                        &
                 -2.D0*(RY_1(L-1,1)+RUFP*DERY(L-1,1) )            &
                 *(RJABG(NP,L,1,1)+RJABG(NP,L,2,0)/VC)            &
                 +(2.D0*RUFP*RY_1(L-2,0)+RUFP**2*DERY(L-2,0) )    &
                 *RJABG(NP,L,0,2)                                 &
                 +(2.D0*RUFP*RJ_1(L+2,0)+RUFP**2*DERJ(L+2,0) )    &
                 *RYABG(NP,L,0,0)                                 &
                 -2.D0*(RJ_1(L+1,1)+2.D0*RUFP/VC*RJ_1(L+2,0)      &
                 +RUFP*( DERJ(L+1,1)+RUFP/VC*DERJ(L+2,0) ) )      &
                 *RYABG(NP,L,1,1)                                 &
                 +DERJ(L,2)*RYABG(NP,L,2,0) )

            PSI0G(NP,L) = ( RY_1(L,0)*RJABG(NP,L,0,0)  &
                 + RJ_1(L,0)*RYABG(NP,L,0,0) )/VC

            PSI02G(NP,L) = 0.5D0*(                     &
                   RY_1(L,0)*RJABG(NP,L,1,1)           &
                 - RUFP*RY_1(L-1,1)*RJABG(NP,L,0,2)    &
                 + RUFP*RJ_1(L+1,1)*RYABG(NP,L,0,0)    &
                 - RJ_1(L,2)*RYABG(NP,L,1,1)  )

            PSI022G(NP,L) = VC/8.D0*(                                &
                 RY_1(L,0)*RJABG(NP,L,2,0)                           &
                 -2.D0*RUFP*RY_1(L-1,1)                              &
                 *(RJABG(NP,L,1,1)+RJABG(NP,L,2,0)/VC)               &
                 +RUFP**2*RY_1(L-2,0)*RJABG(NP,L,0,2)                &
                 +RUFP**2*RJ_1(L+2,0)*RYABG(NP,L,0,0)                &
                 -2.D0*RUFP                                          &
                 *(RJ_1(L+1,1)+RUFP/VC*RJ_1(L+2,0) )*RYABG(NP,L,1,1) &
                 +RJ_1(L,2)*RYABG(NP,L,2,0) )

            DPSI1G(NP,L)=(DERY(L,1)*RJABG(NP,L,0,1)  &
                 + DERJ(L,1)*RYABG(NP,L,0,1) )/VC

            DPSI11G(NP,L) = 0.5D0*(                              &
                  DERY(L,1)*RJABG(NP,L,1,0)                      &
                 -(RY_1(L-1,0)+RUFP*DERY(L-1,0))*RJABG(NP,L,0,1) &
                 +(RJ_1(L+1,0)+RUFP*DERJ(L+1,0))*RYABG(NP,L,0,1) &
                 -DERJ(L,1)*RYABG(NP,L,1,0) )

         END DO
      END DO

!
!----------END OF ON GIRD
!---- END OF PSI AND ITS DERIVATIVES

!
!--- LOCAL DIFFUSION COEFFICIENTS
!
      FACT= 4.D0*PI*RGAMH*1.D20 &
           * (PTFP0(NSA) / AMFP(NSA))*RNFD0(NSB)
      FACT2=4.D0*PI*RGAMH*1.D20 &
           * (PTFP0(NSA) / AMFP(NSA))**2*RNFD0(NSB)
!-----DCPP & FCPP-----------------
!      DO NP=2,NPMAX+1
      DO NP=NPS,NPENDWG
         RGAMA=SQRT(1.D0+PG(NP,NSSA)**2*THETA0(NSSA))
         RUFP = (PTFP0(NSA)*PG(NP,NSSA))/AMFP(NSA)
         DO NTH=1,NTHMAX
            DO L=LLMAX,LLMIN,-1
               SUMA = DPSI02G(NP,L) *PLM(NTH,L)
               SUMB = DPSI022G(NP,L)*PLM(NTH,L)
               SUMC = PSI0G(NP,L)   *PLM(NTH,L)
               SUMD = PSI02G(NP,L)  *PLM(NTH,L)
               SUME = PSI022G(NP,L) *PLM(NTH,L)
               SUMF = DPSI1G(NP,L)  *PLM(NTH,L)
               SUMG = DPSI11G(NP,L) *PLM(NTH,L)

               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA) &
                 + FACT * RGAMA / RUFP *(                        &
                 2.D0*RGAMA**2* SUMA                             &
                 -8.D0*(RGAMA/VC)**2* SUMB                       &
                 -RUFP* SUMC - L*(L+1)/RUFP* SUMD                &
                 +(8.D0*(RUFP/VC)**2+4.D0*L*(L+1))/(RUFP*VC**2)  &
                 *SUME    )

               FCPP2(NTH,NP,NR,NSB,NSA) = FCPP2(NTH,NP,NR,NSB,NSA) &
                    +FACT2 * AMFP(NSA)/AMFD(NSB)*RGAMA             &
                    *( -SUMF + 2.D0/VC**2*SUMG )
            END DO
         END DO 
      END DO


!! p->0 limit
      IF(NPSTART.eq.1)THEN
         DO NTH=1,NTHMAX
            FCPP2(NTH,1,NR,NSB,NSA) = 0.D0
            
            PNFP_NLR=0.D0
            PCRIT=0.D0
            CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R2)
            PNFP_NLR=PCRIT
            CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R2)
            DCPP2(NTH,1,NR,NSB,NSA) = RGAMH/(3.D0*RINT0)*(  &
                 (AMFD(NSB)*PTFP0(NSA))                     &
                 /(AMFP(NSA)*PTFD0(NSB))*RINT2 )            &
                 *RNFD(NR,NSB)*1.D20
         END DO
      END IF
! FOR BAUNDARY
!      DCPP2B(1,1,NR,NSB,NSA) = &
!           (DCPP2(ITL(NR),1,NR,NSB,NSA)+DCPP2(ITL(NR)-1,1,NR,NSB,NSA))*0.5D0
!      DCPP2B(2,1,NR,NSB,NSA) = &
!           (DCPP2(ITL(NR),1,NR,NSB,NSA)+DCPP2(ITL(NR)+1,1,NR,NSB,NSA))*0.5D0
!      FCPP2B(1,1,NR,NSB,NSA) = 0.D0
!      FCPP2B(2,1,NR,NSB,NSA) = 0.D0
!
!! end of p->0 limit


!-----DCTT & FCTH--------------------

!      DO NP=1,NPMAX
      DO NP=NPSTARTW,NPENDWM
         RGAMA=SQRT(1.D0+PM(NP,NSSA)**2*THETA0(NSSA))
         RUFP = (PTFP0(NSA)*PM(NP,NSSA))/AMFP(NSA)
         DO NTH=1,NTHMAX+1
            SUMA = 0.D0
            SUMB = 0.D0
            SUMC = 0.D0
            SUMD = 0.D0
            SUME = 0.D0
            SUMF = 0.D0
            SUMG = 0.D0
            SUMH = 0.D0
            DO L=LLMAX,LLMIN,-1
               SUMA = SUMA + DPSI02M(NP,L) * PLG(NTH,L)
               SUMB = SUMB + PSI02M(NP,L) * PLG(NTH,L)
               SUMD = SUMD + DPSI022M(NP,L) * PLG(NTH,L)
               SUME = SUME + PSI022M(NP,L) * PLG(NTH,L)
            END DO
            DO L=LLMIN,LLMAX
               SUMC = SUMC + PSI02M(NP,L) * D2PLG(NTH,L)
               SUMF = SUMF + PSI022M(NP,L) * D2PLG(NTH,L)
               SUMG = SUMG + PSI1M(NP,L) * D1PLG(NTH,L)
               SUMH = SUMH + PSI11M(NP,L) * D1PLG(NTH,L)
            END DO

            DCTT2(NTH,NP,NR,NSB,NSA) = DCTT2(NTH,NP,NR,NSB,NSA)   &
                 +FACT/RGAMA/RUFP                                 &
                 *(-RGAMA**2*SUMA - RUFP/VC**2*SUMB - SUMC/RUFP   &
                 +4.D0*RGAMA**2/VC**2*SUMD                        &
                 - (4.D0*RUFP/VC**4*SUME - 4.D0/RUFP/VC**2*SUMF ) &
                 )

            FCTH2(NTH,NP,NR,NSB,NSA) = FCTH2(NTH,NP,NR,NSB,NSA) &
                 + FACT2 * AMFP(NSA)/AMFD(NSB)/RGAMA/RUFP       & 
                 *(- SUMG + 2.D0/VC**2*SUMH )

         END DO
      END DO

!-----DCPT & DCTP
      IF(NPSTARTW.eq.1)THEN
         DO NTH=1,NTHMAX
            DCPT2(NTH,1,NR,NSB,NSA)=0.D0
         END DO
      END IF

      DO NP=NPS,NPENDWG
         RGAMA=SQRT(1.D0+PG(NP,NSSA)**2*THETA0(NSSA))
         RUFP = (PTFP0(NSA)*PG(NP,NSSA))/AMFP(NSA)
         DO NTH=1,NTHMAX
            SUMA = 0.D0
            SUMB = 0.D0
            SUMC = 0.D0
            SUMD = 0.D0
            DO L=LLMIN,LLMAX
!            DO L=LLMAX,LLMIN,-1
               SUMA = SUMA + DPSI022G(NP,L) * D1PLM(NTH,L)
               SUMB = SUMB + DPSI02G(NP,L) * D1PLM(NTH,L)
               SUMC = SUMC + PSI022G(NP,L) * D1PLM(NTH,L)
               SUMD = SUMD + PSI02G(NP,L) * D1PLM(NTH,L)
            END DO
            DCPT2(NTH,NP,NR,NSB,NSA) = DCPT2(NTH,NP,NR,NSB,NSA)  &
                 +FACT*RGAMA/RUFP                                &
                 *( 4.D0/VC**2*SUMA - SUMB                       &
                 - 4.D0/(RUFP*VC**2)*SUMC + SUMD/RUFP )
         END DO
      END DO

      DO NP=NPSTARTW,NPENDWM
         RGAMA=SQRT(1.D0+PM(NP,NSSA)**2*THETA0(NSSA))
         RUFP = (PTFP0(NSA)*PM(NP,NSSA))/AMFP(NSA)
         DO NTH=1,NTHMAX+1
            SUMA = 0.D0
            SUMB = 0.D0
            SUMC = 0.D0
            SUMD = 0.D0
            DO L=LLMIN,LLMAX
!            DO L=LLMAX,LLMIN,-1
               SUMA = SUMA + DPSI022M(NP,L) * D1PLG(NTH,L)
               SUMB = SUMB + DPSI02M(NP,L) * D1PLG(NTH,L)
               SUMC = SUMC + PSI022M(NP,L) * D1PLG(NTH,L)
               SUMD = SUMD + PSI02M(NP,L) * D1PLG(NTH,L)
            END DO
            DCTP2(NTH,NP,NR,NSB,NSA) = DCTP2(NTH,NP,NR,NSB,NSA) &
                 +FACT*RGAMA/RUFP                               &
                 *( 4.D0/VC**2*SUMA - SUMB                      &
                 - 4.D0/(RUFP*VC**2)*SUMC + SUMD/RUFP )
         END DO
      END DO

      RETURN
      END SUBROUTINE FPCALC_NLR
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---- RECURRENCE EQUATION OF FIRST KIND LEGENDRE FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE FKLF_JY(RUFP,RJ_1,RY_1)

      IMPLICIT NONE

      real(rkind),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      integer:: L, NA
      real(rkind):: RUFP, RGAMA, RZ, RSIGMA, ra1, ra2

!     FKLF_J(L,NA) = P^{-L-1/2}_{A-1/2}
!     FKLF_Y(L,NA) = P^{L+1/2}_{A-1/2}
!     RJ_1(L,NA) = j_{L[1]NA}
!     RY_1(L,NA) = y_{L[1]NA}

!
!---- INITIALIZE
!
      Do L = -2, LLMAX+2
!         Do NA=-1,2
!            FKLF_J(L,NA)=0.D0
!            FKLF_Y(L,NA)=0.D0
!         END DO
         DO NA=0,2
            RJ_1(L,NA)=0.D0
            RY_1(L,NA)=0.D0
         END DO
      END DO
!
!---- END OF INITIALIZATION
!
      RGAMA= sqrt(1.D0+(RUFP/VC)**2)
      RZ = RUFP/VC
      IF(RZ.le.1.d-1)THEN
         RSIGMA =                 &
              35.D0/1152.D0*RZ**9 &
              -5.D0/112.D0*RZ**7  &
              +3.D0/40.D0*RZ**5   &
              -RZ**3/6.D0         &
              +RZ

      ELSE
         RSIGMA = LOG(RZ+SQRT(1.D0+RZ**2))
      END IF

!
!---- FIRST KIND LEGENDRE FUN!TION FOR J
!
!      FKLF_J(0,0) = sqrt(2.D0*VC/PI/RUFP)*RSIGMA
!      FKLF_J(0,1) = sqrt(2.D0*RUFP/PI/VC)
!      FKLF_J(0,2) = sqrt(2.D0*RUFP/PI/VC)*RGAMA
!      FKLF_J(0,-1) = FKLF_J(0,1)

!
!---- FIRST KIND LEGENDRE FUNCTION FOR Y
!
!      FKLF_Y(0,0) = sqrt(2.D0*VC/PI/RUFP)
!      FKLF_Y(0,1) = sqrt(2.D0*VC/PI/RUFP)*RGAMA
!      FKLF_Y(0,2) = sqrt(2.D0*RUFP/PI/VC)*(VC/RUFP+2.D0*RUFP/VC)
!      FKLF_Y(0,-1) = FKLF_Y(0,1)
!
!---- RECURRENCE EQUATION 
!
!      Do L = 0, LLMAX+1
!         Do NA = 0, 2
!            FKLF_J(L+1,NA) = (RGAMA*FKLF_J(L,NA)-FKLF_J(L,NA-1)) &
!                 *VC/RUFP/DBLE(NA+L+1)
!            FKLF_Y(L+1,NA) =                                           &
!                 ( (NA-L-1)*RGAMA*FKLF_Y(L,NA)-(NA+L)*FKLF_Y(L,NA-1) ) &
!                 *VC/RUFP
!         END DO
!         FKLF_J(L+1,-1)=FKLF_J(L+1,1)
!         FKLF_Y(L+1,-1)=FKLF_Y(L+1,1)
!      END DO

!
!--- ANALYZED J, Y
!

         RJ_1(0,0) = RSIGMA/RZ
         RJ_1(0,1) = 1.D0
         RJ_1(0,2) = RGAMA
         
         RJ_1(1,0) = (RGAMA*RSIGMA-RZ)/RZ**2
         RJ_1(1,1) = (RZ*RGAMA-RSIGMA)*0.5D0/RZ**2
         RJ_1(1,2) = RZ / 3.D0
         
         IF(RZ.le.2.D-4)THEN
            RJ_1(2,0) = ( (4.D0/15.D0)*RZ**2  &
                 -(15.D0/112.D0+3.D0/80.D0)*RZ**4 )/4.D0
            RJ_1(2,1) =(3.D0/12.D0*RZ**2 - 29.D0/80.D0*RZ**4 )/6.D0
            RJ_1(2,2) = (8.D0/5.D0*RZ**2 - 4.D0/7.D0*RZ**4 )/24.D0
         ELSE
            RJ_1(2,0) = ( (2.D0*RGAMA**2+1.D0)*RSIGMA-3.D0*RGAMA*RZ ) &
                 /4.D0/RZ**3
            RJ_1(2,1) =((RGAMA**2+2.D0)*RZ-3.D0*RGAMA*RSIGMA)/6.D0/RZ**3
            RJ_1(2,2) =  &
                (2.D0*RGAMA*RZ**3-3.D0*RZ*RGAMA+3.D0*RSIGMA)/24.D0/RZ**3
         END IF

         RY_1(0,0) = -1.D0/RZ
         RY_1(0,1) = -RGAMA/RZ
         RY_1(0,2) = -(1.D0+2.D0*RZ**2)/RZ

         RY_1(1,0) = -RGAMA/RZ**2
         RY_1(1,1) = -1.D0/RZ**2
         RY_1(1,2) = -(1.D0-2.D0*RZ**2)*RGAMA/RZ**2

         RY_1(2,0) = -(1.D0+2.D0*RGAMA**2)/RZ**3
         RY_1(2,1) = -3.D0*RGAMA/RZ**3
         RY_1(2,2) = -3.D0/RZ**3

         IF(LLMAX.ge.1)THEN
            IF(RZ.le.1.D-2)THEN
            RJ_1(3,0) =                                   &
                 ( (81.D0/80.D0-75.D0/112.D0)*RZ**3       &
                    -13.D0/320.D0*RZ**5-9.D0/160.D0*RZ**7 &
                    )/36.D0
            RJ_1(3,1) =                                   &
                 ( (11.D0/16.D0-0.9D0+75.D0/112.D0)*RZ**3 &
                 + (1.D0/8.D0+15.D0/28.D0)*RZ**5)/48.D0
            RJ_1(3,2) =                             &
                 ( (29.D0/240.D0-5.D0/112.D0)*RZ**3 &
                 -(11.D0/336.D0+3.D0/320.D0)*RZ**5 )*15.D0/120.D0

            ELSE
            RJ_1(3,0) =                                   &
                 ( (6.D0*RGAMA**2+9.D0)*RGAMA*RSIGMA      &
                 -(11.D0*RGAMA**2+4.D0)*RZ )/36.D0/RZ**4
            RJ_1(3,1) =                                   &
                 ((2.D0*RGAMA**2+13.D0)*RGAMA*RZ          &
                 -3.D0*(4.D0*RGAMA**2+1.D0)*RSIGMA)/48.D0/RZ**4
            RJ_1(3,2) =                                       &
                 (2.D0*RGAMA**2*RZ**3-(7.D0*RGAMA**2+8.D0)*RZ &
                 +15.D0*RGAMA*RSIGMA )/120.D0/RZ**4
            END IF
            RY_1(3,0) = -(3.D0*RGAMA*(3.D0+2.D0*RGAMA**2))/RZ**4
            RY_1(3,1) = -3.D0*(1.D0+4.D0*RGAMA**2)/RZ**4
            RY_1(3,2) = -15.D0*RGAMA/RZ**4

         END IF

         IF(LLMAX.ge.2)THEN
            IF(RZ.le.1.D-1)THEN
               RJ_1(4,0) = ( (125.D0/128.D0+3675.D0/1152.D0-249.D0/70.D0)*RZ**4 &
               )/576.D0
               ra2=35.D0/1152.D0-5.D0/224.D0-3.D0/320.D0-1.D0/96.D0-5.D0/128.D0
               RJ_1(4,1) = ( &
                    (-105.D0*ra2-480.D0/105.D0)*RZ**4  &
                    )/720.D0
               RJ_1(4,2) = (301.D0/128.D0+3675.D0/1152.D0-225.D0/56.D0) & 
                    *RZ**4/1440.D0
            ELSE
               ra1 = 105.D0*(RSIGMA-RGAMA*RZ)
               RJ_1(4,0) =                               &
                    (3.D0*(8.D0*RZ**4+4.D1*RZ**2)*RSIGMA &
                    -50.D0*RGAMA*RZ**3+RA1)/576.D0/RZ**5

!               RJ_1(4,1) = ((6.D0*rgama**4+83.D0*RGAMA**2+16.D0)*RZ &
!                    -(60.D0*RGAMA**2+45.D0)*RGAMA*RSIGMA )/720.D0/RZ**5
               RJ_1(4,1) = ((6.D0*RZ**4+95.D0*RZ**2+105.D0)*RZ &
                    -(60.D0*RZ**2+105.D0)*RGAMA*RSIGMA )/720.D0/RZ**5
               RJ_1(4,2) = ((90.D0*RGAMA**2+15.D0)*RSIGMA +              &
                    (4.D0*RGAMA**2*RZ**2-24.D0*RGAMA**2-81.D0)*RGAMA*RZ) &
                    /1440.D0/RZ**5               
            END IF
            RY_1(4,0) = &
                 -3.D0*(8.D0*RGAMA**4+24.D0*RGAMA**2+3.D0)/RZ**5
            RY_1(4,1) = &
                 -15.D0*RGAMA*(4.D0*RGAMA**2+3.D0)/RZ**5
            RY_1(4,2) = &
                 -15.D0*(6.D0*RGAMA**2+1.D0)/RZ**5
         END IF

      DO NA=0,2
!         Do L=0,LLMAX+2
!             RJ_1(L,NA)=SQRT(PI*VC/2.D0/RUFP)*FKLF_J(L,NA)
!             RY_1(L,NA)=SQRT(PI*VC/2.D0/RUFP)*FKLF_Y(L,NA)*(-1)**(-L-1)
!         END DO
         RJ_1(-1,NA) = - RY_1(0,NA)
         RY_1(-1,NA) = RJ_1(0,NA)
         RJ_1(-2,NA) = RY_1(1,NA)
         RY_1(-2,NA) = - RJ_1(1,NA)
      END DO


 935  FORMAT(4E12.4)
      Return
      END SUBROUTINE FKLF_JY

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     DERIVATION J Y
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DERIV_JY(RUFP,RJ_1,RY_1,DERJ,DERY)

      IMPLICIT NONE

!      real(rkind),DIMENSION(-2:LLMAX+2, -1:2):: FKLF_J,FKLF_Y
      real(rkind),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1, DERJ, DERY
      integer:: L, NA
      real(rkind):: RUFP, RGAMA, RZ, RSIGMA

!
      Do L = -2, LLMAX+2
         DO NA=0,2
            DERY(L,NA)=0.D0
            DERJ(L,NA)=0.D0
         END DO
      END DO
!
!---- END OF INITIALIZATION
!
      RGAMA= sqrt(1.D0+(RUFP/VC)**2)

      DO NA = 0, 2
         DO L= -1, LLMAX+2
           DERJ(L,NA) = RJ_1(L-1,NA)/(VC*RGAMA) - (L+1)/RUFP*RJ_1(L,NA)
         END DO
         DO L= -2, LLMAX+1
           DERY(L,NA) =-RY_1(L+1,NA)/(VC*RGAMA) + L/RUFP*RY_1(L,NA)
         END DO
      END DO

      RETURN
      END SUBROUTINE DERIV_JY

!----------------------------------------------------

      SUBROUTINE INTEGRATION_RJAB_RYAB(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)

      USE libspl1d
      IMPLICIT NONE

      integer,parameter::LNM=5

      real(rkind),DIMENSION(NPMAX,-1:LNM):: FPL

      real(rkind),DIMENSION(NTHMAX+3):: TX,TY,DF
      real(rkind),DIMENSION(4,NTHMAX+3):: UTY
      real(rkind),dimension(NTHMAX+3)::UTY0
      real(rkind),DIMENSION(NPMAX+3):: TX1,TY1,DF1,UTY10
      real(rkind),DIMENSION(4,NPMAX+3):: UTY1

      real(rkind),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX, 0:2, -1:2):: RJABM
      real(rkind),DIMENSION(NPSTART:NPENDWG , 0:LLMAX, 0:2, -1:2):: RJABG 
      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX, 0:2, -1:2):: RYABM
      real(rkind),DIMENSION(NPSTART:NPENDWG , 0:LLMAX, 0:2, -1:2):: RYABG

      integer:: NP, NTH, NSA, NSB, L, LLMIN, NI, NA, NNP, NPG, NSSA, NSSB
      integer:: IER, NS, NPS
      real(rkind):: SUM1, SUM2, SUM3, SUM4, SUM5
      real(rkind):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2
      real(rkind):: vtatb, pabbar, ptatb, PMAX2, testF

      NSSA=NS_NSA(NSA)
      NSSB=NS_NSB(NSB)
      IF(NPSTART.eq.1)THEN
         NPS=2
      ELSE
         NPS=NPSTART
      END IF

      DO NI = 0, 2
      DO NA = 0, 2

      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSSB)**2*THETA0(NSSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSSB))/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=PM(NNP,NSSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSSB)**(2+NI))/RGAMB &
                 *RJ_1(L+NI,NA)
         END DO
         TX1(NPMAX+2)=PMAX(NSSB)
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         PMAX2=PMAX(NSSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSSA)
            IF(PCRIT.le.PMAX(NSSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
!         DO NPG=2,NPMAX+1
         DO NPG=NPSTART,NPENDWG            
            IF(NPG.eq.1)THEN
               RJABG(1,L,NI,NA)=0.D0
            ELSE
               PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSSA)
               IF(PCRIT.le.PMAX(NSSB)) THEN
                  CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
                  RJABG(NPG,L,NI,NA)=SUM3*(PTFD0(NSB)/AMFD(NSB))**NI
               ELSE
                  RJABG(NPG,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
               ENDIF
            END IF
         END DO
      END DO
!-------
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSSB)**2*THETA0(NSSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSSB))/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=PM(NNP,NSSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSSB)**(2+NI))/RGAMB &
                 *RY_1(L-NI,NA)
         END DO
         TX1(NPMAX+2)=PMAX(NSSB)
         TY1(NPMAX+2)=0.D0
         IF(NI.ne.0)THEN
            DF1(1) = 0.D0
         ELSE
            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*FPL(1,L)
!            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*(3.D0*FPL(1,L)-FPL(2,L))/2.D0
         END IF
         DF1(NPMAX+2)   = 0.D0
         PMAX2=PMAX(NSSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSSA)
            IF(PCRIT.le.PMAX(NSSB)) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
!         DO NPG=1,NPMAX+1
         DO NPG=NPSTART,NPENDWG
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSSA)
            IF(PCRIT.le.PMAX(NSSB)) THEN
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABG(NPG,L,NI,NA)=(PSUM-SUM5)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABG(NPG,L,NI,NA)=0.D0
            ENDIF
         END DO
      END DO
      END DO
      END DO

      END SUBROUTINE INTEGRATION_RJAB_RYAB

!------
      REAL*8 FUNCTION FPFN0R2(X)
!
      IMPLICIT NONE
!      real(rkind)::FPFN0R2
      real(rkind),INTENT(IN):: X
      real(rkind)::PN
!
      PN=X
      FPFN0R2=PN**2*FPRMXW2(PN)
!
      RETURN
      END FUNCTION FPFN0R2
!
! ===============================================================
!
      REAL*8 FUNCTION FPFN2R2(X)

!      real(rkind):: FPFN2R2
      real(rkind),INTENT(IN):: X
      real(rkind):: A, PN, B

      A=1.D0
      PN=A*(X+PNFP_NLR)
      B=PN*SQRT(1.D0+PN**2*THETA0L_NLR)
      FPFN2R2=A*B*FPRMXW2(PN)

      RETURN
      END FUNCTION FPFN2R2

!
! ===============================================================
!
      FUNCTION FPRMXW2(PN)

      real(rkind):: FPRMXW2
      real(rkind),INTENT(IN):: PN
      real(rkind):: EX

      EX=(1.D0-SQRT(1.D0+PN**2*THETA0L_NLR))/THETAL_NLR
      IF (EX.LT.-100.D0)THEN
         FPRMXW2=0.D0
      ELSE
         FPRMXW2=EXP(EX)
      ENDIF

      RETURN
      END FUNCTION FPRMXW2

!-------------------------------------------

      SUBROUTINE INTEGRATION_RJAB_RYAB_FINE(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM,FPLS1)

      USE libspl1d
      IMPLICIT NONE

      integer,parameter::LNM=5

      integer,intent(IN):: NSA, NSB
      real(rkind),DIMENSION(NPMAX,-1:LNM),INTENT(IN):: FPL
      real(rkind),DIMENSION(NPMAX,-1:LNM):: FPL0
      double precision,dimension(-1:LNM),intent(in):: FPLS1

      real(rkind),DIMENSION(2*NPMAX+3):: TX1,TY1,DF1,UTY10
      real(rkind),DIMENSION(4,2*NPMAX+3):: UTY1

      real(rkind),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RJABM
      real(rkind),DIMENSION(NPSTART:NPENDWG , 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RJABG 
      real(rkind),DIMENSION(NPSTARTW:NPENDWM, 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RYABM
      real(rkind),DIMENSION(NPSTART:NPENDWG , 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RYABG

      integer:: NP, NTH, L, LLMIN, NI, NA, NNP, NPG, NSSA, NSSB
      integer:: IER, NS, NPF
      real(rkind):: SUM1, SUM2, SUM3, SUM4, SUM5
      real(rkind):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2
      real(rkind):: vtatb, pabbar, ptatb, PMAX2, testF, testP
      integer:: N_fine_range

      THETA0L_NLR=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
      NSSA=NS_NSA(NSA)
      NSSB=NS_NSB(NSB)
      N_fine_range=1

!      open(9,file='spl_j_fine.dat')
!      open(10,file='spl_y_fine.dat')

      DO L=0,LLMAX
         FPL0(1,L)=FPLS1(L) 
         TX1(1)=0.D0
         TY1(1)=FPL0(1,L)
         DO NP=1,NPMAX
            TX1(NP+1)=PM(NP,NSSB)
            TY1(NP+1)=FPL(NP,L)
         END DO
         TX1(NPMAX+2)=PMAX(NSSB)
         TY1(NPMAX+2)=0.D0
         DF1(1)=0.D0
         DF1(NPMAX+2)=0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         DO NP=1,NPMAX-1
            testP=PM(NP,NSSB)/PMAX(NSSB)*N_fine_range
            CALL SPL1DF(testP,testF,TX1,UTY1,NPMAX+2,IER)
            FPL0(NP+1,L)=testF
         END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO NI = 0, 2
      DO NA = 0, 2

      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSSB)/PMAX(NSSB)*N_fine_range
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSSB))
            RUFP = (PTFD0(NSB)*testP)/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+NI))/RGAMB &
                 *RJ_1(L+NI,NA)
!            IF(L.eq.1) TY1(NNP+1)=0.D0
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSSB)**2*THETA0(NSSB))
            RUFP = (PTFD0(NSSB)*PM(NNP,NSSB))/AMFD(NSB)
            IF(PM(NNP,NSSB).ge.1.D0*N_fine_range)THEN
               NPF=NPF+1
               CALL FKLF_JY(RUFP,RJ_1,RY_1)
               TX1(NPF)=PM(NNP,NSSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSSB)**(2+NI))/RGAMB &
                    *RJ_1(L+NI,NA)
!               IF(L.eq.1) TY1(NPF)=0.D0
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSSB)
         TY1(NPF+1)=0.D0
         DF1(1)   = 0.D0
         DF1(NPF+1)   = 0.D0
         PMAX2=PMAX(NSSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPF
!         IF(NSB.eq.1.and.nsa.eq.1) THEN
!            write(9,'(3I2,1P14E14.6)') L, NA,NI,TX1(NP),TY1(NP)
!         END IF
!         END DO

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSSA)
            IF(PCRIT.le.PMAX(NSSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPF+1,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
!         DO NPG=2,NPMAX+1
         DO NPG=NPSTART,NPENDWG
            IF(NPG.eq.1)THEN
               RJABG(1,L,NI,NA)=0.D0
            ELSE
               PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSSA)
               IF(PCRIT.le.PMAX(NSSB)) THEN
                  CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPF+1,IER)
                  RJABG(NPG,L,NI,NA)=SUM3*(PTFD0(NSB)/AMFD(NSB))**NI
               ELSE
                  RJABG(NPG,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
               ENDIF
            END IF
         END DO
!      WRITE(9,*) " "
!      WRITE(9,*) " "
!      WRITE(10,*) " "
!      WRITE(10,*) " "
      END DO
!-------
      DO L=0,LLMAX
         TX1(1)=0.D0
         IF(L.eq.1.and.NI.eq.0)THEN
            TY1(1)=-(AMFD(NSB)*VC/PTFD0(NSB))**2*FPL0(1,L)
         ELSEIF(L.eq.2.and.NI.eq.0)THEN
            TY1(1)=-3.D0*(AMFD(NSB)*VC/PTFD0(NSB))**3*FPL0(1,L)/(PM(1,NSSB)*VTFP0(NSA)/VTFD0(NSB)/PMAX(NSSB))
         ELSE
            TY1(1)=0.D0
         END IF
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSSB)/PMAX(NSSB)*N_fine_range
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSSB))
            RUFP = (PTFD0(NSB)*testP)/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+NI))/RGAMB &
                 *RY_1(L-NI,NA)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSSB)**2*THETA0(NSSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSSB))/AMFD(NSB)
            IF(PM(NNP,NSSB).ge.1.D0*N_fine_range)THEN
               NPF=NPF+1
               CALL FKLF_JY(RUFP,RJ_1,RY_1)
               TX1(NPF)=PM(NNP,NSSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSSB)**(2+NI))/RGAMB &
                    *RY_1(L-NI,NA)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSSB)
         TY1(NPF+1)=0.D0
         IF(L.eq.0.and.NI.ne.0)THEN
            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*FPL0(1,L)
!            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*(3.D0*FPL(1,L)-FPL(2,L))/2.D0
         ELSEIF(L.eq.2.and.NI.eq.1)THEN
            DF1(1) = - (AMFD(NSB)*VC/PTFD0(NSB))**2*FPL0(1,L)
         ELSEIF(L.eq.2.and.NI.eq.0)THEN
            DF1(1) = 3.D0*(AMFD(NSB)*VC/PTFD0(NSB))**3*FPL0(1,L)/(PM(1,NSSB)*VTFP0(NSA)/VTFD0(NSA)/PMAX(NSSB))**2
         ELSE
            DF1(1) = 0.D0
         END IF
         DF1(NPF+1)   = 0.D0
         PMAX2=PMAX(NSSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)
!         CALL SPL1DI_inv(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPF
!            IF(NSB.eq.1.and.nsa.eq.2) THEN
!               write(10,'(3I2,1P14E14.6)') L,NA,NI,TX1(NP),TY1(NP)
!            END IF
!         END DO

         DO NP=NPSTARTW,NPENDWM
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSSA)
            IF(PCRIT.le.PMAX(NSSB)) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPF+1,IER)
!               CALL SPL1DI_inv(PCRIT,SUM4,TX1,UTY1,UTY10,NPF+1,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0(NSB)/AMFD(NSB))**NI
!               IF(NSA.eq.1.and.NSB.eq.2)THEN
!                  WRITE(*,'(A,3I4,4E14.6)') "not zero", l, na, NP, PCRIT, RYABM(NP,L,NI,NA), PSUM, SUM4
!               END IF
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
         DO NPG=NPSTART,NPENDWG
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSSA)
            IF(PCRIT.le.PMAX(NSSB)) THEN
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPF+1,IER)
               RYABG(NPG,L,NI,NA)=(PSUM-SUM5)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABG(NPG,L,NI,NA)=0.D0
            ENDIF

         END DO
!      WRITE(9,*) " "
!      WRITE(9,*) " "
!      WRITE(10,*) " "
!      WRITE(10,*) " "
      END DO
!      WRITE(9,*) " "
!      WRITE(9,*) " "
!      WRITE(10,*) " "
!      WRITE(10,*) " "
      END DO
      END DO

      END SUBROUTINE INTEGRATION_RJAB_RYAB_FINE

! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL_calcnr(PML,NR,NS)

      USE plprof
      USE libbes,ONLY: beseknx
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_prf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_calcnr

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
         RNFDL=RN_TEMP(NR,NS)/RNFD0L
         RTFDL=RT_TEMP(NR,NS)
!            RTFDL=RT_BULK(NR,NS)
!         IF(NRANK.eq.0) WRITE(*,'(A,2I5,10E14.6)') "TEST_cnr ", NR, NS, RNFDL, RN_TEMP(NR,NS)/RNFD0L, RTFDL, RT_TEMP(NR,NS), RN_TEMP(NR,NS)
      END IF

      THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
      THETAL=THETA0L*RTFDL/RTFD0L
      Z=1.D0/THETAL
      DKBSL=BESEKNX(2,Z)
      FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
           *RTFD0L
!      EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
      EX=0.D0
      FPMXWL_calcnr=FACT*EXP(EX)

      RETURN
      END FUNCTION FPMXWL_calcnr


      END MODULE fpcalcnr

