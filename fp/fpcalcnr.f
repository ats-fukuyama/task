C     $Id$
C
C ************************************************************
C
C      CALCULATION OF NONLINEAR-RELAIVISTIC COLLISIONAL OPERATOR
C
C ************************************************************
C
      SUBROUTINE FPCALC_NLR(NR,NSB,NSA)
C
      INCLUDE 'fpcomm.inc'

      PARAMETER (N=NPM+2,M=NTHM+2,LNM=5)
      DIMENSION PLM(M,-1:LNM),PLG(M,-1:LNM)
      DIMENSION D1PLM(M,-1:LNM)
      DIMENSION D1PLG(M,-1:LNM),D2PLG(M,-1:LNM)
      DIMENSION PLTEMP(0:LNM)
      DIMENSION FPL(N,-1:LNM)

      DIMENSION TX(NTHM+2),TY(NTHM+2),DF(NTHM+2)
      DIMENSION UTY(4,NTHM+2),UTY0(NTHM+2)
      DIMENSION TX1(NPM+2),TY1(NPM+2),DF1(NPM+2)
      DIMENSION UTY1(4,NPM+2),UTY10(NPM+2)

      DIMENSION FKLF_J(-2:LLMAX+2, -1:2),FKLF_Y(-2:LLMAX+2, -1:2)
      DIMENSION RJ_1(-2:LLMAX+2, 0:2), RY_1(-2:LLMAX+2, 0:2)
      DIMENSION RJABM(N+1,0:LLMAX,0:2,-1:2),RJABG(N+1,0:LLMAX,0:2,-1:2)
      DIMENSION RYABM(N+1,0:LLMAX,0:2,-1:2),RYABG(N+1,0:LLMAX,0:2,-1:2)

      DIMENSION DPSI02M(N,0:LLMAX),DPSI022M(N,0:LLMAX)
      DIMENSION PSI0M(N,0:LLMAX),PSI02M(N,0:LLMAX),PSI022M(N,0:LLMAX)
      DIMENSION PSI1M(N,0:LLMAX),PSI11M(N,0:LLMAX),DPSI11M(N,0:LLMAX)

      DIMENSION DPSI02G(N+1,0:LLMAX),DPSI022G(N+1,0:LLMAX)
      DIMENSION PSI0G(N+1,0:LLMAX),PSI02G(N+1,0:LLMAX)
      DIMENSION PSI022G(N+1,0:LLMAX)
      DIMENSION DPSI1G(N+1,0:LLMAX), DPSI11G(N+1,0:LLMAX)

      DIMENSION DERJ(-2:LLMAX+2, 0:2), DERY(-2:LLMAX+2, 0:2)

C
C----- DEFINITION OF LOCAL QUANTITIES -------------
C
      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA)
     &     /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
      TMC2FD0=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
      TMC2FP0=(PTFP0(NSA)/(AMFP(NSA)*VC))**2

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

         CALL PAGES
         CALL GRD1D(1,thg,plg,M,NTHMAX+1,LLMAX+2,'@PLG:@',0)
         CALL GRD1D(2,thg,d1plg,M,NTHMAX+1,LLMAX+2,'@D1PLG:@',0)
         CALL GRD1D(3,thg,d2plg,M,NTHMAX+1,LLMAX+2,'@D2PLG:@',0)
         CALL PAGEE
      ENDIF

C
C     ----- Legendre expansion of distribution funstion FNS -----
C
      NS=NS_NSB(NSB)
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
C---- INTEGRAL ABBREVIATIONS
C
      DO NI = 0, 2
      DO NA = 0, 2

      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         NNP=1
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            RUFP = (PTFD0(NSB)*PM(NNP))/AMFD(NSB)
            CALL FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP)**(2+NI))/RGAMB
     &           *RJ_1(L+NI,NA)

         END DO
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABG(NPG,L,NI,NA)=SUM3*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RJABG(NPG,L,NI,NA)=PSUM*(PTFD0(NSB)/AMFD(NSB))**NI
            ENDIF
         END DO
      END DO
C-------
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            RUFP = (PTFD0(NSB)*PM(NNP))/AMFD(NSB)
            CALL FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
            TX1(NNP+1)=PM(NNP)
            TY1(NNP+1)=FPL(NNP,L)*(PM(NNP)**(2+NI))/RGAMB
     &           *RY_1(L-NI,NA)
         END DO
         TX1(NPMAX+2)=PMAX
         TY1(NPMAX+2)=0.D0
         DF1(1)   = 0.D0
         DF1(NPMAX+2)   = 0.D0
         CALL SPL1D(TX1,TY1,DF1,UTY1,NPMAX+2,3,IER)
         CALL SPL1DI0(TX1,UTY1,UTY10,NPMAX+2,IER)
         CALL SPL1DI(PMAX,PSUM,TX1,UTY1,UTY10,NPMAX+2,IER)

         DO NP=1,NPMAX
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NPG)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABG(NPG,L,NI,NA)=(PSUM-SUM5)*(PTFD0(NSB)/AMFD(NSB))**NI
            ELSE
               RYABG(NPG,L,NI,NA)=0.D0
            ENDIF

         END DO
      END DO
      END DO
      END DO
C
C---- END OF INTEGRALS
C

C
C---- PSI AND IT'S DERIVATIONS
C
C
C---------- FOR MID OF GRID
C
      DO L = 0,LLMAX
         DO NP = 1, NPMAX
            RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
            RUFP = (PTFP0(NSA)*PM(NP))/AMFP(NSA)

            CALL FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
            CALL DERIV_JY(RUFP,RJ_1,RY_1,DERJ,DERY)

            DPSI02M(NP,L) = 0.5D0*( 
     &           DERY(L,0)*RJABM(NP,L,1,1)
     &           -(RY_1(L-1,1)+RUFP*DERY(L-1,1) )*RJABM(NP,L,0,2)
     &           +(RJ_1(L+1,1)+RUFP*DERJ(L+1,1) )*RYABM(NP,L,0,0)
     &           -DERJ(L,2)*RYABM(NP,L,1,1) )

            DPSI022M(NP,L) = VC/8.D0*(
     &           DERY(L,0)*RJABM(NP,L,2,0)
     &           -2.D0*(RY_1(L-1,1)+RUFP*DERY(L-1,1) )
     &           *(RJABM(NP,L,1,1)+RJABM(NP,L,2,0)/VC)
     &           +(2.D0*RUFP*RY_1(L-2,0)+RUFP**2*DERY(L-2,0) )
     &           *RJABM(NP,L,0,2)
     &           +(2.D0*RUFP*RJ_1(L+2,0)+RUFP**2*DERJ(L+2,0) )
     &           *RYABM(NP,L,0,0)
     &           -2.D0*(RJ_1(L+1,1)+2.D0*RUFP/VC*RJ_1(L+2,0)
     &           +RUFP*( DERJ(L+1,1)+RUFP/VC*DERJ(L+2,0) ) )
     &           *RYABM(NP,L,1,1)
     &           +DERJ(L,2)*RYABM(NP,L,2,0) )


            PSI02M(NP,L) = 0.5D0*(
     &             RY_1(L,0)*RJABM(NP,L,1,1) 
     &           - RUFP*RY_1(L-1,1)*RJABM(NP,L,0,2)
     &           + RUFP*RJ_1(L+1,1)*RYABM(NP,L,0,0)
     &           - RJ_1(L,2)*RYABM(NP,L,1,1)  )


            PSI022M(NP,L) = VC/8.D0*(
     &           RY_1(L,0)*RJABM(NP,L,2,0)
     &           -2.D0*RUFP*RY_1(L-1,1)
     &           *(RJABM(NP,L,1,1)+RJABM(NP,L,2,0)/VC)
     &           +RUFP**2*RY_1(L-2,0)*RJABM(NP,L,0,2)
     &           +RUFP**2*RJ_1(L+2,0)*RYABM(NP,L,0,0)
     &           -2.D0*RUFP
     &           *(RJ_1(L+1,1)+RUFP/VC*RJ_1(L+2,0) )*RYABM(NP,L,1,1)
     &           +RJ_1(L,2)*RYABM(NP,L,2,0) )

            PSI11M(NP,L) = 0.5D0 * (
     &           RY_1(L,1)*RJABM(NP,L,1,0)
     &           -RUFP*RY_1(L-1,0)*RJABM(NP,L,0,1)
     &           +RUFP*RJ_1(L+1,0)*RYABM(NP,L,0,1)
     &           -RJ_1(L,1)*RYABM(NP,L,1,0)
     &           )

            PSI1M(NP,L) = ( RY_1(L,1)*RJABM(NP,L,0,1)
     &           +RJ_1(L,1)*RYABM(NP,L,0,1) )/VC

         END DO
      END DO
C
C----------- END OF MID 
C
C
C----------- ON GRID
C
      DO L=0,LLMAX
         DO NP = 2, NPMAX+1
            RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
            RUFP = (PTFP0(NSA)*PG(NP))/AMFP(NSA)
            CALL FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)
            CALL DERIV_JY(RUFP,RJ_1,RY_1,DERJ,DERY)

            DPSI02G(NP,L) = 0.5D0*( 
     &           DERY(L,0)*RJABG(NP,L,1,1)
     &           -(RY_1(L-1,1)+RUFP*DERY(L-1,1) )*RJABG(NP,L,0,2)
     &           +(RJ_1(L+1,1)+RUFP*DERJ(L+1,1) )*RYABG(NP,L,0,0)
     &           -DERJ(L,2)*RYABG(NP,L,1,1) )

            DPSI022G(NP,L) = VC/8.D0*(
     &           DERY(L,0)*RJABG(NP,L,2,0)
     &           -2.D0*(RY_1(L-1,1)+RUFP*DERY(L-1,1) )
     &           *(RJABG(NP,L,1,1)+RJABG(NP,L,2,0)/VC)
     &           +(2.D0*RUFP*RY_1(L-2,0)+RUFP**2*DERY(L-2,0) )
     &           *RJABG(NP,L,0,2)
     &           +(2.D0*RUFP*RJ_1(L+2,0)+RUFP**2*DERJ(L+2,0) )
     &           *RYABG(NP,L,0,0)
     &           -2.D0*(RJ_1(L+1,1)+2.D0*RUFP/VC*RJ_1(L+2,0)
     &           +RUFP*( DERJ(L+1,1)+RUFP/VC*DERJ(L+2,0) ) )
     &           *RYABG(NP,L,1,1)
     &           +DERJ(L,2)*RYABG(NP,L,2,0) )

            PSI0G(NP,L) = ( RY_1(L,0)*RJABG(NP,L,0,0) 
     &           + RJ_1(L,0)*RYABG(NP,L,0,0) )/VC

            PSI02G(NP,L) = 0.5D0*(
     &             RY_1(L,0)*RJABG(NP,L,1,1) 
     &           - RUFP*RY_1(L-1,1)*RJABG(NP,L,0,2)
     &           + RUFP*RJ_1(L+1,1)*RYABG(NP,L,0,0)
     &           - RJ_1(L,2)*RYABG(NP,L,1,1)  )

            PSI022G(NP,L) = VC/8.D0*(
     &           RY_1(L,0)*RJABG(NP,L,2,0)
     &           -2.D0*RUFP*RY_1(L-1,1)
     &           *(RJABG(NP,L,1,1)+RJABG(NP,L,2,0)/VC)
     &           +RUFP**2*RY_1(L-2,0)*RJABG(NP,L,0,2)
     &           +RUFP**2*RJ_1(L+2,0)*RYABG(NP,L,0,0)
     &           -2.D0*RUFP
     &           *(RJ_1(L+1,1)+RUFP/VC*RJ_1(L+2,0) )*RYABG(NP,L,1,1)
     &           +RJ_1(L,2)*RYABG(NP,L,2,0) )

            DPSI1G(NP,L)=(DERY(L,1)*RJABG(NP,L,0,1)
     &           + DERJ(L,1)*RYABG(NP,L,0,1) )/VC

            DPSI11G(NP,L) = 0.5D0*(
     &            DERY(L,1)*RJABG(NP,L,1,0)
     &           -(RY_1(L-1,0)+RUFP*DERY(L-1,0))*RJABG(NP,L,0,1)
     &           +(RJ_1(L+1,0)+RUFP*DERJ(L+1,0))*RYABG(NP,L,0,1)
     &           -DERJ(L,1)*RYABG(NP,L,1,0) )




c            rz = RUFP/VC
c            RSIGMA=LOG(rz+SQRT(1.D0+rz**2))
c         RSIGMA = RZ - RZ**3/6.D0 +9.D0/120.D0*RZ**5
c            rgama=SQRT(1+rz**2)
c            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NP)


         END DO
      END DO

C
C----------END OF ON GIRD
C---- END OF PSI AND IT'S DERIVATIVES

C
C--- LOCAL DIFFUSION COEFFICIENTS
C
      FACT=4.D0*PI*RGAMH*1.D20
     &     * (PTFP0(NSA) / AMFP(NSA))*RNFP0(NSB) 
      FACT2=4.D0*PI*RGAMH*1.D20
     &     * (PTFP0(NSA) / AMFP(NSA))**2*RNFP0(NSB)
C-----DCPP & FCPP-----------------
      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RUFP = (PTFP0(NSA)*PG(NP))/AMFP(NSA)
         DO NTH=1,NTHMAX
            DO L=LLMIN,LLMAX
               SUMA = DPSI02G(NP,L)*PLM(NTH,L)
               SUMB = DPSI022G(NP,L)*PLM(NTH,L)
               SUMC = PSI0G(NP,L) * PLM(NTH,L)
               SUMD = PSI02G(NP,L) * PLM(NTH,L)
               SUME = PSI022G(NP,L) * PLM(NTH,L)
               SUMF = DPSI1G(NP,L) * PLM(NTH,L)
               SUMG = DPSI11G(NP,L) * PLM(NTH,L)

               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)
     &           + FACT * RGAMA / RUFP *(
     &           2.D0*RGAMA**2* SUMA
     &           -8.D0*(RGAMA/VC)**2* SUMB
     &           -RUFP* SUMC - L*(L+1)/RUFP* SUMD
     &           +(8.D0*(RUFP/VC)**2+4.D0*L*(L+1))/(RUFP*VC**2)
     &           *SUME    )

               FCPP2(NTH,NP,NR,NSB,NSA) = FCPP2(NTH,NP,NR,NSB,NSA)
     &              +FACT2 * AMFP(NSA)/AMFD(NSB)*RGAMA
     &              *( -SUMF + 2.D0/VC**2*SUMG )
            END DO
         END DO
      END DO


      DO NTH=1,NTHMAX
         DCPP2(NTH,1,NR,NSB,NSA)
     &        =RGAMH*RNFD(NR,NS)*1.D20*(2.D0/(3.D0*SQRT(PI)))
     &       *(PTFP0(NSA)/(SQRT(2.D0)*PTFD(NR,NSB)))*AMFD(NSB)/AMFP(NSA)
     &        +DCPP2(NTH,1,NR,NSB,NSA)
         FCPP2(NTH,1,NR,NSB,NSA)=0.D0
      END DO

C-----DCTT & FCTH--------------------

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RUFP = (PTFP0(NSA)*PM(NP))/AMFP(NSA)
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
            DCTT2(NTH,NP,NR,NSB,NSA) = DCTT2(NTH,NP,NR,NSB,NSA)
     &           +FACT/RGAMA/RUFP
     &           *(-RGAMA**2*SUMA - RUFP/VC**2*SUMB + SUMC/RUFP
     &           +4.D0*RGAMA**2/VC**2*SUMD 
     &           - (4.D0*RUFP/VC**4*SUME - 4.D0/RUFP/VC**2*SUMF )
     &           )

            FCTH2(NTH,NP,NR,NSB,NSA) = FCTH2(NTH,NP,NR,NSB,NSA)
     &           + FACT2 * AMFP(NSA)/AMFD(NSB)/RGAMA/RUFP
     &           *(- SUMG + 2.D0/VC**2*SUMH )

         END DO
      END DO

C-----DCPT & DCTP
      DO NTH=1,NTHMAX
         DCPT2(NTH,1,NR,NSB,NSA)=0.D0
      END DO

      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RUFP = (PTFP0(NSA)*PG(NP))/AMFP(NSA)
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
            DCPT2(NTH,NP,NR,NSB,NSA) = DCPT2(NTH,NP,NR,NSB,NSA)
     &           +FACT*RGAMA/RUFP
     &           *( 4.D0/VC**2*SUMA - SUMB 
     &           - 4.D0/(RUFP*VC**2)*SUMC + SUMD/RUFP )
         END DO
      END DO

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RUFP = (PTFP0(NSA)*PM(NP))/AMFP(NSA)
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
            DCTP2(NTH,NP,NR,NSB,NSA) = DCTP2(NTH,NP,NR,NSB,NSA)
     &           +FACT*RGAMA/RUFP
     &           *( 4.D0/VC**2*SUMA - SUMB 
     &           - 4.D0/(RUFP*VC**2)*SUMC + SUMD/RUFP )
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C---- RECURRENCE EQUATION OF FIRST KIND LEGENDRE FUNCTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FKLF_JY(RUFP,FKLF_J,FKLF_Y,RJ_1,RY_1)

      INCLUDE 'fpcomm.inc'
      DIMENSION FKLF_J(-2:LLMAX+2, -1:2),FKLF_Y(-2:LLMAX+2, -1:2)
      DIMENSION RJ_1(-2:LLMAX+2, 0:2), RY_1(-2:LLMAX+2, 0:2)


C     FKLF_J(L,NA) = P^{-L-1/2}_{A-1/2}
C     FKLF_Y(L,NA) = P^{L+1/2}_{A-1/2}
C     RJ_1(L,NA) = j_{L[1]NA}
C     RY_1(L,NA) = y_{L[1]NA}

C
C---- INITIALIZE
C
      Do L = -2, LLMAX+2
         Do NA=-1,2
            FKLF_J(L,NA)=0.D0
            FKLF_Y(L,NA)=0.D0
         END DO
         DO NA=0,2
            RJ_1(L,NA)=0.D0
            RY_1(L,NA)=0.D0
         END DO
      END DO
C
C---- END OF INITIALIZATION
C
      RGAMA= sqrt(1.D0+(RUFP/VC)**2)
      RZ = RUFP/VC
      IF(RZ.le.1.d-1)THEN
         RSIGMA = 
     &        5.D0/384.D0*RZ**9
     &        -5.D0/112.D0*RZ**7
     &        +3.D0/40.D0*RZ**5
     &        -RZ**3/6.D0 
     &        +RZ

      ELSE
         RSIGMA = LOG(RZ+SQRT(1.D0+RZ**2))
      END IF

C
C---- FIRST KIND LEGENDRE FUNCTION FOR J
C
      FKLF_J(0,0) = sqrt(2.D0*VC/PI/RUFP)
     &     *RSIGMA
      FKLF_J(0,1) = sqrt(2.D0*RUFP/PI/VC)
      FKLF_J(0,2) = sqrt(2.D0*RUFP/PI/VC)*RGAMA
      FKLF_J(0,-1) = FKLF_J(0,1)

C
C---- FIRST KIND LEGENDRE FUNCTION FOR Y
C
      FKLF_Y(0,0) = sqrt(2.D0*VC/PI/RUFP)
      FKLF_Y(0,1) = sqrt(2.D0*VC/PI/RUFP)*RGAMA
      FKLF_Y(0,2) = sqrt(2.D0*RUFP/PI/VC)
     &     *(VC/RUFP+2.D0*RUFP/VC)
      FKLF_Y(0,-1) = FKLF_Y(0,1)
C
C---- RECURRENCE EQUATION 
C
      Do L = 0, LLMAX+1
         Do NA = 0, 2
            FKLF_J(L+1,NA) = (RGAMA*FKLF_J(L,NA)-FKLF_J(L,NA-1))
     &           *VC/RUFP/DBLE(NA+L+1)
            FKLF_Y(L+1,NA) = 
     &           ( (NA-L-1)*RGAMA*FKLF_Y(L,NA)-(NA+L)*FKLF_Y(L,NA-1) )
     &           *VC/RUFP
         END DO

         FKLF_J(L+1,-1)=FKLF_J(L+1,1)
         FKLF_Y(L+1,-1)=FKLF_Y(L+1,1)
      END DO

C
C--- ANALYZED J, Y
C

         RJ_1(0,0) = RSIGMA/RZ
         RJ_1(0,1) = 1.D0
         RJ_1(0,2) = RGAMA
         
         RJ_1(1,0) = (RGAMA*RSIGMA-RZ)/RZ**2
         RJ_1(1,1) = (RZ*RGAMA-RSIGMA)*0.5D0/RZ**2
         RJ_1(1,2) = RZ / 3.D0
         
         IF(RZ.le.2.D-4)THEN
            RJ_1(2,0) = ( (4.D0/15.D0)*RZ**2 
     &           -(15.D0/112.D0+3.D0/80.D0)*RZ**4 )/4.D0
            RJ_1(2,1) =(3.D0/12.D0*RZ**2 - 29.D0/80.D0*RZ**4 )/6.D0
            RJ_1(2,2) = (8.D0/5.D0*RZ**2 - 4.D0/7.D0*RZ**4 )/24.D0
         ELSE
            RJ_1(2,0) = ( (2.D0*RGAMA**2+1.D0)*RSIGMA-3.D0*RGAMA*RZ )
     &           /4.D0/RZ**3
            RJ_1(2,1) =((RGAMA**2+2.D0)*RZ-3.D0*RGAMA*RSIGMA)/6.D0/RZ**3
            RJ_1(2,2) = 
     &          (2.D0*RGAMA*RZ**3-3.D0*RZ*RGAMA+3.D0*RSIGMA)/24.D0/RZ**3
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
            RJ_1(3,0) = 
     &           ( (81.D0/80.D0-75.D0/112.D0)*RZ**3
     &              -13.D0/320.D0*RZ**5-9.D0/160.D0*RZ**7
     &              )/36.D0
            RJ_1(3,1) =
     &           ( (11.D0/16.D0-0.9D0+75.D0/112.D0)*RZ**3
     &           + (1.D0/8.D0+15.D0/28.D0)*RZ**5)/48.D0
            RJ_1(3,2) =
     &           ( (29.D0/240.D0-5.D0/112.D0)*RZ**3
     &           -(11.D0/336.D0+3.D0/320.D0)*RZ**5 )*15.D0/120.D0

            ELSE
            RJ_1(3,0) = 
     &           ( (6.D0*RGAMA**2+9.D0)*RGAMA*RSIGMA
     &           -(11.D0*RGAMA**2+4.D0)*RZ )/36.D0/RZ**4
            RJ_1(3,1) = 
     &           ((2.D0*RGAMA**2+13.D0)*RGAMA*RZ
     &           -3.D0*(4.D0*RGAMA**2+1.D0)*RSIGMA)/48.D0/RZ**4
            RJ_1(3,2) = 
     &           (2.D0*RGAMA**2*RZ**3-(7.D0*RGAMA**2+8.D0)*RZ
     &           +15.D0*RGAMA*RSIGMA )/120.D0/RZ**4
            END IF
            RY_1(3,0) = -(3.D0*RGAMA*(3.D0+2.D0*RGAMA**2))/RZ**4
            RY_1(3,1) = -3.D0*(1.D0+4.D0*RGAMA**2)/RZ**4
            RY_1(3,2) = -15.D0*RGAMA/RZ**4

         END IF

         IF(LLMAX.ge.2)THEN
            RSIGMA = LOG(RZ+SQRT(1.D0+RZ**2))
            IF(RZ.le.3.D-2)THEN
               RJ_1(4,0) = (9.D0/5.D0)*RZ**4/576.D0
               RJ_1(4,1) = ((357.D0/64.D0)*RZ**4 
     &              -(999.D0/224.D0+19.D0/16.D0)*RZ**6 )/720.D0
               RJ_1(4,2) = (525.D0/96.D0-267.D0/56.D0)*RZ**4/1440.D0
            ELSE
               ra1 = 105.D0*(RSIGMA-RGAMA*RZ)
               RJ_1(4,0) = 
c((24.D0*RGAMA**4+72.D0*RGAMA**2+9.D0)*RSIGMA
c     &              -5.D0*(10.D0*RGAMA**2+11.D0)*RGAMA*RZ)/576.D0/RZ**5
     &              (3.D0*(8.D0*RZ**4+4.D1*RZ**2)*RSIGMA
     &              -50.D0*RGAMA*RZ**3+RA1)/576.D0/RZ**5

               RJ_1(4,1) = ((6.D0*rgama**4+83.D0*RGAMA**2+16.D0)*RZ
     &              -(60.D0*RGAMA**2+45.D0)*RGAMA*RSIGMA )/720.D0/RZ**5
               RJ_1(4,2) = ((90.D0*RGAMA**2+15.D0)*RSIGMA + 
     &              (4.D0*RGAMA**2*RZ**2-24.D0*RGAMA**2-81.D0)*RGAMA*RZ)
     &              /1440.D0/RZ**5               
            END IF
            RY_1(4,0) = 
     &           -3.D0*(8.D0*RGAMA**4+24.D0*RGAMA**2+3.D0)/RZ**5
            RY_1(4,1) = 
     &           -15.D0*RGAMA*(4.D0*RGAMA**2+3.D0)/RZ**5
            RY_1(4,2) = 
     &           -15.D0*(6.D0*RGAMA**2+1.D0)/RZ**5
         END IF

      DO NA=0,2
c         Do L=0,LLMAX+2
c             RJ_1(L,NA)=SQRT(PI*VC/2.D0/RUFP)*FKLF_J(L,NA)
c             RY_1(L,NA)=SQRT(PI*VC/2.D0/RUFP)*FKLF_Y(L,NA)*(-1)**(-L-1)
c         END DO
         RJ_1(-1,NA) = - RY_1(0,NA)
         RY_1(-1,NA) = RJ_1(0,NA)
         RJ_1(-2,NA) = RY_1(1,NA)
         RY_1(-2,NA) = - RJ_1(1,NA)
      END DO


 935  FORMAT(4E12.4)
      Return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DERIVATION J Y
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DERIV_JY(RUFP,RJ_1,RY_1,DERJ,DERY)

      INCLUDE 'fpcomm.inc'
      DIMENSION FKLF_J(-2:LLMAX+2, -1:2),FKLF_Y(-2:LLMAX+2, -1:2)
      DIMENSION RJ_1(-2:LLMAX+2, 0:2), RY_1(-2:LLMAX+2, 0:2)
      DIMENSION DERJ(-2:LLMAX+2, 0:2), DERY(-2:LLMAX+2, 0:2)

C
      Do L = -2, LLMAX+2
         DO NA=0,2
            DERY(L,NA)=0.D0
            DERJ(L,NA)=0.D0
         END DO
      END DO
C
C---- END OF INITIALIZATION
C
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
      END

