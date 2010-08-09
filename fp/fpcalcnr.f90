!     $Id$
!
! ************************************************************
!
!      CALCULATION OF NONLINEAR-RELAIVISTIC COLLISIONAL OPERATOR
!
! ************************************************************
!
      MODULE fpcalcnr

      USE fpcomm
      real(8):: PNFP, TMC2FD0, TMC2FD


      contains

!-------------------------------------------------------------

      SUBROUTINE FPCALC_NLR(NR,NSB,NSA)

      Implicit none
!      PARAMETER (N=NPM+2,M=NTHM+2,LNM=5)
      integer,parameter::LNM=5
      real(8),DIMENSION(NTHMAX+3,-1:LNM):: PLM, PLG, D1PLM, D1PLG, D2PLG
      real(8),DIMENSION(0:LNM):: PLTEMP
      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL, FPL0

      real(8),DIMENSION(NTHMAX+3):: TX,TY,DF
      real(8),DIMENSION(4,NTHMAX+3):: UTY
      real(8),dimension(NTHMAX+3)::UTY0
      real(8),DIMENSION(NPMAX+3):: TX1,TY1,DF1,UTY10
      real(8),DIMENSION(4,NPMAX+3):: UTY1

!!      real(8),DIMENSION(-2:LLMAX+2, -1:2):: FKLF_J,FKLF_Y
      real(8),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(8),DIMENSION(-2:LLMAX+2, 0:2):: DERJ, DERY

      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2):: RJABM,RJABG 
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2):: RYABM,RYABG

      real(8),DIMENSION(NPMAX+1, 0:LLMAX):: DPSI02M,DPSI022M
      real(8),DIMENSION(NPMAX+1, 0:LLMAX):: PSI0M,PSI02M,PSI022M
      real(8),DIMENSION(NPMAX+1, 0:LLMAX):: PSI1M,PSI11M,DPSI11M

      real(8),DIMENSION(NPMAX+2, 0:LLMAX):: DPSI02G,DPSI022G
      real(8),DIMENSION(NPMAX+2, 0:LLMAX):: PSI0G,PSI02G,PSI022G
      real(8),DIMENSION(NPMAX+2, 0:LLMAX):: DPSI1G, DPSI11G


      integer:: NP, NTH, NSA, NSB, L, NR, LLMIN, NI, NA, NNP, NPG, NSBA
      integer:: IER
      real(8):: RGAMH, SUM1, SUM2, SUM3, SUM4, SUM5
      real(8):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2, RUFD
      real(8):: SUMA, SUMB, SUMC, SUMD, SUME, SUMF, SUMG, SUMH
      real(8):: vtatb, ptatb, PMAX2, RINT0, RINT2, ES0, ES2, testF, testP
      real(8):: DKBSL0, DKBSL1, DKBSL2, Z
      interface
         DOUBLE PRECISION FUNCTION BESEKN(N,X)
           real(8) :: X
           integer :: N
         end function BESEKN
      end interface
!
!----- DEFINITION OF LOCAL QUANTITIES -------------
!
      TMC2FD0=THETA0(NSB)
      TMC2FD =(PTFD(NR,NSB)/(AMFD(NSB)*VC))**2 
      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
           /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
      NSBA=NSB_NSA(NSA)
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
!      open(8,file='FPL_t1.dat')
!      NS=NS_NSB(NSB)
      DO L=LLMIN,LLMAX
         DO NP=1,NPMAX
            TX(1)=0.D0
            TY(1)=0.D0
            DO NTH=1,NTHMAX
               TX(NTH+1)=THM(NTH)
               TY(NTH+1)=FNS(NTH,NP,NR,NSB)*PLM(NTH,L)*SINM(NTH)
            END DO
            TX(NTHMAX+2)=PI
            TY(NTHMAX+2)=0.D0
            DF(1)= FNS(1,NP,NR,NSB)
            DF(NTHMAX+2)= (-1)**(L+1)*FNS(NTHMAX,NP,NR,NSB)
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPL(NP,L)=0.5D0*(2.D0*L+1.D0)*SUM1
!            WRITE(8,*) L, NP, FPL(NP,L)
         END DO
!         WRITE(8,*)" "
!         WRITE(8,*)" "
      END DO
!      close(8)
!
!---- INTEGRAL ABBREVIATIONS
!
!      CALL INTEGRATION_RJAB_RYAB(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)
      IF(NCALCNR.eq.0)THEN
         CALL INTEGRATION_RJAB_RYAB_FINE(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)
      ELSE
         CALL INTEGRATION_RJAB_RYAB_weighp(NR,NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)
      END IF
!
!---- END OF INTEGRALS
!

!
!---- PSI AND ITS DERIVATIONS
!
!
!---------- FOR MID OF GRID
!
      DO L = 0,LLMAX
         DO NP = 1, NPMAX
            RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
            RUFP = (PTFP0(NSA)*PM(NP,NSBA))/AMFP(NSA)

!            CALL FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
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
      DO L=0,LLMAX
         DO NP = 2, NPMAX+1
            RGAMA=SQRT(1.D0+PG(NP,NSBA)**2*THETA0(NSA))
            RUFP = (PTFP0(NSA)*PG(NP,NSBA))/AMFP(NSA)
!            CALL FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
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
      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP,NSBA)**2*THETA0(NSA))
         RUFP = (PTFP0(NSA)*PG(NP,NSBA))/AMFP(NSA)
         DO NTH=1,NTHMAX
            DO L=LLMIN,LLMAX
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
      DO NTH=1,NTHMAX
         FCPP2(NTH,1,NR,NSB,NSA) = 0.D0

         PNFP=0.D0
!         THETA0(NSB)=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
         PCRIT=0.D0
         CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R2)
         PNFP=PCRIT
         CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R2)
         DCPP2(NTH,1,NR,NSB,NSA) = RGAMH/(3.D0*RINT0)*(  &
              (AMFD(NSB)*PTFP0(NSA))                     &
              /(AMFP(NSA)*PTFD0(NSB))*RINT2 )            &
              *RNFD(NR,NSB)*1.D20

!         DCPP2(NTH,1,NR,NSB,NSA)                                          &
!              =RGAMH*RNFD(NR,NSB)*1.D20*(2.D0/(3.D0*SQRT(PI)))             &
!              *(PTFP0(NSA)/(SQRT(2.D0)*PTFD(NR,NSB)))*AMFD(NSB)/AMFP(NSA) &
!              +DCPP2(NTH,1,NR,NSB,NSA)
      END DO
!! end of p->0 limit


!-----DCTT & FCTH--------------------

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
         RUFP = (PTFP0(NSA)*PM(NP,NSBA))/AMFP(NSA)
         DO NTH=1,NTHMAX+1
            SUMA = 0.D0
            SUMB = 0.D0
            SUMC = 0.D0
            SUMD = 0.D0
            SUME = 0.D0
            SUMF = 0.D0
            SUMG = 0.D0
            SUMH = 0.D0
            DO L=LLMIN,LLMAX
               SUMA = SUMA + DPSI02M(NP,L) * PLG(NTH,L)
               SUMB = SUMB + PSI02M(NP,L) * PLG(NTH,L)
               SUMC = SUMC + PSI02M(NP,L) * D2PLG(NTH,L)
               SUMD = SUMD + DPSI022M(NP,L) * PLG(NTH,L)
               SUME = SUME + PSI022M(NP,L) * PLG(NTH,L)
               SUMF = SUMF + PSI022M(NP,L) * D2PLG(NTH,L)
               SUMG = SUMG + PSI1M(NP,L) * D1PLG(NTH,L)
               SUMH = SUMH + PSI11M(NP,L) * D1PLG(NTH,L)
            END DO
            DCTT2(NTH,NP,NR,NSB,NSA) = DCTT2(NTH,NP,NR,NSB,NSA)   &
                 +FACT/RGAMA/RUFP                                 &
                 *(-RGAMA**2*SUMA - RUFP/VC**2*SUMB + SUMC/RUFP   &
                 +4.D0*RGAMA**2/VC**2*SUMD                        &
                 - (4.D0*RUFP/VC**4*SUME - 4.D0/RUFP/VC**2*SUMF ) &
                 )

            FCTH2(NTH,NP,NR,NSB,NSA) = FCTH2(NTH,NP,NR,NSB,NSA) &
                 + FACT2 * AMFP(NSA)/AMFD(NSB)/RGAMA/RUFP       & 
                 *(- SUMG + 2.D0/VC**2*SUMH )

         END DO
      END DO

!-----DCPT & DCTP
      DO NTH=1,NTHMAX
         DCPT2(NTH,1,NR,NSB,NSA)=0.D0
      END DO

      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP,NSBA)**2*THETA0(NSA))
         RUFP = (PTFP0(NSA)*PG(NP,NSBA))/AMFP(NSA)
         DO NTH=1,NTHMAX
            SUMA = 0.D0
            SUMB = 0.D0
            SUMC = 0.D0
            SUMD = 0.D0
            DO L=LLMIN,LLMAX
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

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP,NSBA)**2*THETA0(NSA))
         RUFP = (PTFP0(NSA)*PM(NP,NSBA))/AMFP(NSA)
         DO NTH=1,NTHMAX+1
            SUMA = 0.D0
            SUMB = 0.D0
            SUMC = 0.D0
            SUMD = 0.D0
            DO L=LLMIN,LLMAX
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

!! p -> infty limit 
!      DO NTH=1,NTHMAX
!         RGAMA=SQRT(1.D0+THETA0(NSA)*PMAX(NSBA)**2)
!         RUFD=PTFD0(NSB)/AMFD(NSB)
!         RUFP=PTFP0(NSA)/AMFP(NSA)*PMAX(NSA)
!         Z=1.D0/THETA(NR,NSB)
!         DKBSL0=BESEKN(0,Z)
!         DKBSL1=BESEKN(1,Z)
!         DKBSL2=BESEKN(2,Z)

!         DCPP2(NTH,NPMAX+1,NR,NSB,NSA) = &
!              FACT / (4.D0*PI) &
!              * RUFD**2 &
!              /(RUFP/RGAMA)**3*DKBSL1/DKBSL2 &
!              *(1.D0- DKBSL0/DKBSL1*(RUFD/RGAMA/VC)**2 )

!         DCTT2(NTH,NPMAX,NR,NSB,NSA) = &
!              FACT / (4.D0*PI) &
!              *0.5D0/(RUFP/RGAMA) &
!              *(1.D0-DKBSL1/DKBSL2*( (RUFD/RUFP)**2+RUFD**2/(RGAMA*VC)**2) &
!              +DKBSL0/DKBSL2*(RUFD/RUFP)**2*(RUFD/RGAMA/VC)**2)

!         FCPP2(NTH,NPMAX+1,NR,NSB,NSA) = &
!              -FACT2 / (4.D0*PI) &
!              *AMFP(NSA)/AMFD(NSB) &
!              /(RUFP/RGAMA)**2*DKBSL1/DKBSL2 &
!              *(1.D0-DKBSL0/DKBSL1*(RUFD/RGAMA/VC)**2)

!      END DO
!      DCTT2(NTHMAX+1,NPMAX,NR,NSB,NSA) = DCTT2(NTHMAX,NPMAX,NR,NSB,NSA)

!! end of p -> infty limit

      RETURN
      END SUBROUTINE FPCALC_NLR
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---- RECURRENCE EQUATION OF FIRST KIND LEGENDRE FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      SUBROUTINE FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
      SUBROUTINE FKLF_JY(RUFP,RJ_1,RY_1)

      IMPLICIT NONE

!      real(8),DIMENSION(-2:LLMAX+2, -1:2):: FKLF_J,FKLF_Y
      real(8),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      integer:: L, NA
      real(8):: RUFP, RGAMA, RZ, RSIGMA, ra1

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
              5.D0/384.D0*RZ**9   &
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
            RSIGMA = LOG(RZ+SQRT(1.D0+RZ**2))
            IF(RZ.le.3.D-2)THEN
               RJ_1(4,0) = (9.D0/5.D0)*RZ**4/576.D0
               RJ_1(4,1) = ((357.D0/64.D0)*RZ**4  &
                    -(999.D0/224.D0+19.D0/16.D0)*RZ**6 )/720.D0
               RJ_1(4,2) = (525.D0/96.D0-267.D0/56.D0)*RZ**4/1440.D0
            ELSE
               ra1 = 105.D0*(RSIGMA-RGAMA*RZ)
               RJ_1(4,0) =                               &
                    (3.D0*(8.D0*RZ**4+4.D1*RZ**2)*RSIGMA &
                    -50.D0*RGAMA*RZ**3+RA1)/576.D0/RZ**5

               RJ_1(4,1) = ((6.D0*rgama**4+83.D0*RGAMA**2+16.D0)*RZ &
                    -(60.D0*RGAMA**2+45.D0)*RGAMA*RSIGMA )/720.D0/RZ**5
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

!      real(8),DIMENSION(-2:LLMAX+2, -1:2):: FKLF_J,FKLF_Y
      real(8),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1, DERJ, DERY
      integer:: L, NA
      real(8):: RUFP, RGAMA, RZ, RSIGMA

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

      IMPLICIT NONE

      integer,parameter::LNM=5

      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL

      real(8),DIMENSION(NTHMAX+3):: TX,TY,DF
      real(8),DIMENSION(4,NTHMAX+3):: UTY
      real(8),dimension(NTHMAX+3)::UTY0
      real(8),DIMENSION(NPMAX+3):: TX1,TY1,DF1,UTY10
      real(8),DIMENSION(4,NPMAX+3):: UTY1

      real(8),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2):: RJABM,RJABG 
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2):: RYABM,RYABG

      integer:: NP, NTH, NSA, NSB, L, LLMIN, NI, NA, NNP, NPG, NSBA
      integer:: IER, NS
      real(8):: SUM1, SUM2, SUM3, SUM4, SUM5
      real(8):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2
      real(8):: vtatb, pabbar, ptatb, PMAX2, testF

!      THETA0(NSB)=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
      NSBA=NSB_NSA(NSA)

      DO NI = 0, 2
      DO NA = 0, 2

      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSB))/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=PM(NNP,NSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSB)**(2+NI))/RGAMB &
                 *RJ_1(L+NI,NA)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         PMAX2=PMAX(NSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
         DO NPG=2,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABG(NPG,L,NI,NA)=SUM3*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
              RJABG(NPG,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
      RJABG(1,L,NI,NA)=0.D0
      END DO
!-------
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSB))/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=PM(NNP,NSB)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP,NSB)**(2+NI))/RGAMB &
                 *RY_1(L-NI,NA)
         END DO
         TX1(NPMAX+2)=PMAX(NSB)
         TY1(NPMAX+2)=0.D0
         IF(NI.ne.0)THEN
            DF1(1) = 0.D0
         ELSE
            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*FPL(1,L)
!            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*(3.D0*FPL(1,L)-FPL(2,L))/2.D0
         END IF
         DF1(NPMAX+2)   = 0.D0
         PMAX2=PMAX(NSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
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
      FUNCTION FPFN0R2(X)
!
      IMPLICIT NONE
      real(8)::FPFN0R2
      real(8)::PN, X
!
      PN=X
      FPFN0R2=PN**2*FPRMXW2(PN)
!
      RETURN
      END FUNCTION FPFN0R2
!
! ===============================================================
!
      FUNCTION FPFN2R2(X)

      real(8):: FPFN2R2
      real(8):: X, A, PN, B

      A=1.D0
      PN=A*(X+PNFP)
      B=PN*SQRT(1.D0+PN**2*TMC2FD0)
      FPFN2R2=A*B*FPRMXW2(PN)

      RETURN
      END FUNCTION FPFN2R2

!
! ===============================================================
!
      FUNCTION FPRMXW2(PN)

      real(8):: FPRMXW2
      real(8):: PN, EX

      EX=(1.D0-SQRT(1.D0+PN**2*TMC2FD0))/TMC2FD
      IF (EX.LT.-100.D0)THEN
         FPRMXW2=0.D0
      ELSE
         FPRMXW2=EXP(EX)
      ENDIF

      RETURN
      END FUNCTION FPRMXW2

!-------------------------------------------

      SUBROUTINE INTEGRATION_RJAB_RYAB_FINE(NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)

      IMPLICIT NONE

      integer,parameter::LNM=5

      real(8),DIMENSION(NPMAX+3,-1:LNM),INTENT(IN):: FPL
      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL0

      real(8),DIMENSION(2*NPMAX+3):: TX1,TY1,DF1,UTY10
      real(8),DIMENSION(4,2*NPMAX+3):: UTY1

      real(8),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RJABM,RJABG 
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RYABM,RYABG

      integer:: NP, NTH, NSA, NSB, L, LLMIN, NI, NA, NNP, NPG, NSBA
      integer:: IER, NS, NPF
      real(8):: SUM1, SUM2, SUM3, SUM4, SUM5
      real(8):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2
      real(8):: vtatb, pabbar, ptatb, PMAX2, testF, testP

      TMC2FD0=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
      NSBA=NSB_NSA(NSA)

!      open(9,file='spl_j_fine.dat')
!      open(10,file='spl_y_fine.dat')

      DO L=0,LLMAX
         FPL0(1,L)=0.5D0*( FPL(1,L)+FPL(2,L) & 
              + (FPL(1,L)-FPL(2,L))/( DELP(NSB)*(PM(1,NSB)+PM(2,NSB)) ) &
              *(PM(1,NSB)**2+PM(2,NSB)**2) )
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
      DO NI = 0, 2
      DO NA = 0, 2

      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*testP)/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+NI))/RGAMB &
                 *RJ_1(L+NI,NA)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSB))/AMFD(NSB)
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               CALL FKLF_JY(RUFP,RJ_1,RY_1)
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(2+NI))/RGAMB &
                    *RJ_1(L+NI,NA)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         DF1(1)   = 0.D0
         DF1(NPF+1)   = 0.D0
         PMAX2=PMAX(NSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPF
!         IF(NSB.eq.1.and.nsa.eq.2) THEN
!            write(9,'(3I2,1P14E14.6)') L, NA,NI,TX1(NP),TY1(NP)
!         END IF
!         END DO

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPF+1,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
         DO NPG=2,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPF+1,IER)
               RJABG(NPG,L,NI,NA)=SUM3*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
              RJABG(NPG,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
      RJABG(1,L,NI,NA)=0.D0
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
            TY1(1)=-3.D0*(AMFD(NSB)*VC/PTFD0(NSB))**3*FPL0(1,L)/(PM(1,NSB)*VTFP0(NSA)/VTFD0(NSB)/PMAX(NSB))
         ELSE
            TY1(1)=0.D0
         END IF
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*testP)/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+NI))/RGAMB &
                 *RY_1(L-NI,NA)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSB))/AMFD(NSB)
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               CALL FKLF_JY(RUFP,RJ_1,RY_1)
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(2+NI))/RGAMB &
                    *RY_1(L-NI,NA)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         IF(L.eq.0.and.NI.ne.0)THEN
            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*FPL0(1,L)
!            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*(3.D0*FPL(1,L)-FPL(2,L))/2.D0
         ELSEIF(L.eq.2.and.NI.eq.1)THEN
            DF1(1) = - (AMFD(NSB)*VC/PTFD0(NSB))**2*FPL0(1,L)
         ELSEIF(L.eq.2.and.NI.eq.0)THEN
            DF1(1) = 3.D0*(AMFD(NSB)*VC/PTFD0(NSB))**3*FPL0(1,L)/(PM(1,NSB)*VTFP0(NSA)/VTFD0(NSA)/PMAX(NSB))**2
         ELSE
            DF1(1) = 0.D0
         END IF
         DF1(NPF+1)   = 0.D0
         PMAX2=PMAX(NSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPF
!            IF(NSB.eq.1.and.nsa.eq.2) THEN
!               write(10,'(3I2,1P14E14.6)') L,NA,NI,TX1(NP),TY1(NP)
!            END IF
!         END DO

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPF+1,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
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

!------
!-------------------------------------------
      SUBROUTINE INTEGRATION_RJAB_RYAB_weighp(NR,NSB,NSA,FPL,RJABG,RJABM,RYABG,RYABM)

      IMPLICIT NONE

      integer,parameter::LNM=5

      real(8),DIMENSION(NTHMAX+3,-1:LNM):: PLM, PLG, D1PLM, D1PLG, D2PLG
      real(8),DIMENSION(0:LNM):: PLTEMP
      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL, FPL0

      real(8),DIMENSION(NTHMAX+3):: TX,TY,DF
      real(8),DIMENSION(4,NTHMAX+3):: UTY
      real(8),dimension(NTHMAX+3)::UTY0
!      real(8),DIMENSION(NPMAX+3,-1:LNM),INTENT(IN):: FPL
!      real(8),DIMENSION(NPMAX+3,-1:LNM):: FPL0

      real(8),DIMENSION(2*NPMAX+3):: TX1,TY1,DF1,UTY10
      real(8),DIMENSION(4,2*NPMAX+3):: UTY1

      real(8),DIMENSION(-2:LLMAX+2, 0:2):: RJ_1, RY_1
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RJABM,RJABG 
      real(8),DIMENSION(NPMAX+4, 0:LLMAX, 0:2, -1:2),INTENT(OUT):: RYABM,RYABG

      real(8),dimension(NTHMAX, 2*NPMAX-1):: FNSMG
      real(8),DIMENSION(2*NPMAX+3,-1:LNM):: FPL2

      integer:: NP, NTH, NSA, NSB, L, LLMIN, NI, NA, NNP, NPG, NSBA
      integer:: IER, NS, NPF, NR
      real(8):: SUM1, SUM2, SUM3, SUM4, SUM5
      real(8):: PSUM, PCRIT, RGAMA, RGAMB, RUFP, FACT, FACT2
      real(8):: vtatb, pabbar, ptatb, PMAX2, testF, testP, WPP

      NSBA=NSB_NSA(NSA)
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
! interpolation of FNS
!
      DO NTH=1,NTHMAX
         DO NP=1, NPMAX-1
            WPP=1.D0-WEIGHP(NTH,NP+1,NR,NSB)
            FNSMG(NTH,2*NP-1)=FNS(NTH,NP,NR,NSB)
            FNSMG(NTH,2*NP  )= &
                 (1.D0-WPP)*FNS(NTH,NP+1,NR,NSB)+WPP*FNS(NTH,NP,NR,NSB)

         END DO
         FNSMG(NTH,2*NPMAX-1)=FNS(NTH,NPMAX,NR,NSB)
      END DO

!      DO NP=1,NPMAX-1
!         WRITE(*,*) PM(NP,NSB), FNSMG(1,2*NP-1)!, WEIGHP(1,NP,NR,NSB)
!         WRITE(*,*) PG(NP+1,NSB)  , FNSMG(1,2*NP  ), FNSMG(1,2*NP)-FNSMG(1,2*NP-1)
!      END DO

!
!     ----- Legendre expansion of distribution funstion FNS -----
!
!      open(8,file='FPL_t1.dat')
!      NS=NS_NSB(NSB)
      DO L=LLMIN,LLMAX
         DO NP=1,2*NPMAX-1
            TX(1)=0.D0
            TY(1)=0.D0
            DO NTH=1,NTHMAX
               TX(NTH+1)=THM(NTH)
               TY(NTH+1)=FNSMG(NTH,NP)*PLM(NTH,L)*SINM(NTH)
            END DO
            TX(NTHMAX+2)=PI
            TY(NTHMAX+2)=0.D0
            DF(1)= FNS(1,NP,NR,NSB)
            DF(NTHMAX+2)= (-1)**(L+1)*FNS(NTHMAX,NP,NR,NSB)
            CALL SPL1D(TX,TY,DF,UTY,NTHMAX+2,3,IER)
            CALL SPL1DI0(TX,UTY,UTY0,NTHMAX+2,IER)
            CALL SPL1DI(PI,SUM1,TX,UTY,UTY0,NTHMAX+2,IER)
            FPL2(NP,L)=0.5D0*(2.D0*L+1.D0)*SUM1
!            WRITE(8,*) L, NP, FPL(NP,L)
         END DO
!         WRITE(8,*)" "
!         WRITE(8,*)" "
      END DO
!      close(8)


!      open(9,file='spl_j_fine.dat')
!      open(10,file='spl_y_fine.dat')

      DO L=0,LLMAX
         FPL0(1,L)=0.5D0*( FPL2(1,L)+FPL2(2,L) & 
              + (FPL2(1,L)-FPL2(2,L))/( 0.5D0*DELP(NSB)*(PM(1,NSB)+PG(2,NSB)) ) &
              *(PM(1,NSB)**2+PG(2,NSB)**2) )
         TX1(1)=0.D0
         TY1(1)=FPL0(1,L)
         DO NP=1,2*NPMAX-1
            TX1(NP+1)=PG(1,NSB)+0.5D0*DELP(NSB)*NP
            TY1(NP+1)=FPL2(NP,L)
         END DO
         TX1(2*NPMAX+1)=PMAX(NSB)
         TY1(2*NPMAX+1)=0.D0
         DF1(1)=0.D0
         DF1(2*NPMAX+1)=0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,2*NPMAX+1,3,IER)
         DO NP=1,NPMAX-1
            testP=PM(NP,NSB)/PMAX(NSB)
            CALL SPL1DF(testP,testF,TX1,UTY1,2*NPMAX+1,IER)
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
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*testP)/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+NI))/RGAMB &
                 *RJ_1(L+NI,NA)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSB))/AMFD(NSB)
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               CALL FKLF_JY(RUFP,RJ_1,RY_1)
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(2+NI))/RGAMB &
                    *RJ_1(L+NI,NA)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         DF1(1)   = 0.D0
         DF1(NPF+1)   = 0.D0
         PMAX2=PMAX(NSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPF
!         IF(NSB.eq.1.and.nsa.eq.2) THEN
!            write(9,'(3I2,1P14E14.6)') L, NA,NI,TX1(NP),TY1(NP)
!         END IF
!         END DO

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPF+1,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
         DO NPG=2,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPF+1,IER)
               RJABG(NPG,L,NI,NA)=SUM3*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
              RJABG(NPG,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
      RJABG(1,L,NI,NA)=0.D0
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
            TY1(1)=-3.D0*(AMFD(NSB)*VC/PTFD0(NSB))**3*FPL0(1,L)/(PM(1,NSB)*VTFP0(NSA)/VTFD0(NSB)/PMAX(NSB))
         ELSE
            TY1(1)=0.D0
         END IF
         DO NNP=1,NPMAX-1
            testP=PM(NNP,NSB)/PMAX(NSB)
            RGAMB=SQRT(1.D0+testP**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*testP)/AMFD(NSB)
            CALL FKLF_JY(RUFP,RJ_1,RY_1)
            TX1(NNP+1)=testP
            TY1(NNP+1)=FPL0(NNP+1,L)*(testP**(2+NI))/RGAMB &
                 *RY_1(L-NI,NA)
         END DO
         NPF=NPMAX
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP,NSB)**2*THETA0(NSB))
            RUFP = (PTFD0(NSB)*PM(NNP,NSB))/AMFD(NSB)
            IF(PM(NNP,NSB).ge.1.D0)THEN
               NPF=NPF+1
               CALL FKLF_JY(RUFP,RJ_1,RY_1)
               TX1(NPF)=PM(NNP,NSB)
               TY1(NPF)=FPL(NNP,L)*(PM(NNP,NSB)**(2+NI))/RGAMB &
                    *RY_1(L-NI,NA)
            END IF
         END DO
         TX1(NPF+1)=PMAX(NSB)
         TY1(NPF+1)=0.D0
         IF(L.eq.0.and.NI.ne.0)THEN
            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*FPL0(1,L)
!            DF1(1) = - AMFD(NSB)*VC/PTFD0(NSB)*(3.D0*FPL(1,L)-FPL(2,L))/2.D0
         ELSEIF(L.eq.2.and.NI.eq.1)THEN
            DF1(1) = - (AMFD(NSB)*VC/PTFD0(NSB))**2*FPL0(1,L)
         ELSEIF(L.eq.2.and.NI.eq.0)THEN
            DF1(1) = 3.D0*(AMFD(NSB)*VC/PTFD0(NSB))**3*FPL0(1,L)/(PM(1,NSB)*VTFP0(NSA)/VTFD0(NSA)/PMAX(NSB))**2
         ELSE
            DF1(1) = 0.D0
         END IF
         DF1(NPF+1)   = 0.D0
         PMAX2=PMAX(NSB)
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPF+1,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPF+1,IER)
         CALL SPL1DI(PMAX2,PSUM,TX1,UTY1,UTY10,NPF+1,IER)

!         DO NP=1,NPF
!            IF(NSB.eq.1.and.nsa.eq.2) THEN
!               write(10,'(3I2,1P14E14.6)') L,NA,NI,TX1(NP),TY1(NP)
!            END IF
!         END DO

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPF+1,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG,NSBA)
            IF(PCRIT.le.PMAX(NSB)) THEN
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

      END SUBROUTINE INTEGRATION_RJAB_RYAB_weighp

!------

      END MODULE fpcalcnr

