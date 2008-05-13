C     $Id$
C
C ************************************************************
C
C      CALCULATION OF NONLINEAR-RELAIVISTIC COLLISIONAL OPERATOR
C
C ************************************************************
C
      SUBROUTINE FPCALC_NLR(NR,NS)
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
      AMFD=PA(NS)*AMP
      AEFD=PZ(NS)*AEE
c      AMFP=PA(NSFP)*AMP
      PTFPL=PTFP(NR)
      PTFDL=PTFD(NR,NS)
      VTFPL=VTFP(NR)
      VTFDL=VTFD(NR,NS)
      RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)


      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
      RGAMH=RNUD(NR,NS)*SQRT(2.D0)*VTFD(NR,NS)*AMFP/(RNFP0*PTFP0*1.D20)
      TMC2FD0=(PTFD0/(AMFD*VC))**2
      TMC2FP0=(PTFP0/(AMFP*VC))**2

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
c      NI=1
c      NA=1

      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         NNP=1
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            RUFP = (PTFD0*PM(NNP))/AMFD
c            RUFP = (PTFP0*PM(NNP))/AMFP
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
            PCRIT=(AMFD*PTFP0)/(AMFP*PTFD0)*PM(NP)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM2,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABM(NP,L,NI,NA)=SUM2*(PTFD0/AMFD)**NI
            ELSE
               RJABM(NP,L,NI,NA)=PSUM*(PTFD0/AMFD)**NI
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD*PTFP0)/(AMFP*PTFD0)*PG(NPG)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM3,TX1,UTY1,UTY10,NPMAX+2,IER)
               RJABG(NPG,L,NI,NA)=SUM3*(PTFD0/AMFD)**NI
            ELSE
               RJABG(NPG,L,NI,NA)=PSUM*(PTFD0/AMFD)**NI
            ENDIF
         END DO
      END DO
C-------
      DO L=0,LLMAX
         TX1(1)=0.D0
         TY1(1)=0.D0
         DO NNP=1,NPMAX
            RGAMB=SQRT(1.D0+PM(NNP)**2*TMC2FD0)
            RUFP = (PTFD0*PM(NNP))/AMFD
c            RUFP = (PTFP0*PM(NNP))/AMFP
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
            PCRIT=(AMFD*PTFP0)/(AMFP*PTFD0)*PM(NP)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM4,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABM(NP,L,NI,NA)=(PSUM-SUM4)*(PTFD0/AMFD)**NI
            ELSE
               RYABM(NP,L,NI,NA)=0.D0
            ENDIF
         END DO
         DO NPG=1,NPMAX+1
            PCRIT=(AMFD*PTFP0)/(AMFP*PTFD0)*PG(NPG)
            IF(PCRIT.le.PMAX) THEN
               CALL SPL1DI(PCRIT,SUM5,TX1,UTY1,UTY10,NPMAX+2,IER)
               RYABG(NPG,L,NI,NA)=(PSUM-SUM5)*(PTFD0/AMFD)**NI
            ELSE
               RYABG(NPG,L,NI,NA)=0.D0
            ENDIF

c         IF(NI.eq.1.and.NA.eq.1.and.L.eq.0)
c     &           WRITE(*,*)NPG,L, RYABG(NPG,L,0,0), PCRIT, PTFD0

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
            RUFP = (PTFP0*PM(NP))/AMFP

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
            RUFP = (PTFP0*PG(NP))/AMFP
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




            rz = RUFP/VC
            sinh=LOG(rz+SQRT(1.D0+rz**2))
            rgama=SQRT(1+rz**2)

            IF(L.eq.LLMAX)then
               write(*,765) NP
c     &              ,DPSI022G(NP,L)/VC**2*4.D0
c     &              ,-DPSI02G(NP,L)
c     &              ,-PSI022G(NP,L)*4.D0/RUFP/VC**2
c     &              ,PSI02G(NP,L)/RUFP

c     &              ,(DPSI022G(NP,L)/VC**2*4.D0
c     &              -DPSI02G(NP,L)
c     &              -PSI022G(NP,L)*4.D0/RUFP/VC**2
c     &              +PSI02G(NP,L)/RUFP)
c*rgama/rufp

     &              ,RJABG(NP,L,1,1)*RY_1(L,0)
     &              ,-RJABG(NP,L,0,2)*RUFP*RY_1(L-1,1)
     &              ,RYABG(NP,L,0,0)*RUFP*RJ_1(L+1,1)
     &              ,-RYABG(NP,L,1,1)*RJ_1(L,2)
     &              ,RJABG(NP,L,1,1),RY_1(L,0)
           END IF

         END DO
      END DO
      

 765  FORMAT(I2, 6E14.6)
C
C----------END OF ON GIRD
C---- END OF PSI AND IT'S DERIVATIVES

C
C--- LOCAL DIFFUSION COEFFICIENTS
C
      FACT=4.D0*PI*RGAMH*1.D20
     &     * (PTFP0 / AMFP) 
      FACT2=4.D0*PI*RGAMH*1.D20
     &     * (PTFP0 / AMFP)**2
C-----DCPP & FCPP-----------------
      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RUFP = (PTFP0*PG(NP))/AMFP
         DO NTH=1,NTHMAX
            DO L=LLMIN,LLMAX
               SUMA = DPSI02G(NP,L)*PLM(NTH,L)
               SUMB = DPSI022G(NP,L)*PLM(NTH,L)
               SUMC = PSI0G(NP,L) * PLM(NTH,L)
               SUMD = PSI02G(NP,L) * PLM(NTH,L)
               SUME = PSI022G(NP,L) * PLM(NTH,L)
               SUMF = DPSI1G(NP,L) * PLM(NTH,L)
               SUMG = DPSI11G(NP,L) * PLM(NTH,L)

               DCPP(NTH,NP,NR)=DCPP(NTH,NP,NR)
     &           + FACT * RGAMA / RUFP *(
     &           2.D0*RGAMA**2* SUMA
     &           -8.D0*(RGAMA/VC)**2* SUMB
     &           -RUFP* SUMC - L*(L+1)/RUFP* SUMD
     &           +(8.D0*(RUFP/VC)**2+4.D0*L*(L+1))/(RUFP*VC**2)
     &           *SUME    )

               FCPP(NTH,NP,NR) = FCPP(NTH,NP,NR)
     &              +FACT2 * AMFP/AMFD*RGAMA
     &              *( -SUMF + 2.D0/VC**2*SUMG )

            END DO
         END DO
      END DO

      DO NTH=1,NTHMAX
         DCPP(NTH,1,NR)=RGAMH*RNFD(NR,NS)*1.D20*(2.D0/(3.D0*SQRT(PI)))
     &        *(PTFP0/(SQRT(2.D0)*PTFDL)) +DCPP(NTH,1,NR)
         FCPP(NTH,1,NR)=0.D0
      END DO
C-----DCTT & FCTH--------------------

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RUFP = (PTFP0*PM(NP))/AMFP
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
            DCTT(NTH,NP,NR) = DCTT(NTH,NP,NR)
     &           +FACT/RGAMA/RUFP
     &           *(-RGAMA**2*SUMA - RUFP/VC**2*SUMB + SUMC/RUFP
     &           +4.D0*RGAMA**2/VC**2*SUMD 
     &           - (4.D0*RUFP/VC**4*SUME - 4.D0/RUFP/VC**2*SUMF )
     &           )

            FCTH(NTH,NP,NR) = FCTH(NTH,NP,NR)
     &           + FACT2 * AMFP/AMFD/RGAMA/RUFP
     &           *(- SUMG + 2.D0/VC**2*SUMH )

         END DO
      END DO

C-----DCPT & DCTP
      DO NTH=1,NTHMAX
         DCPT(NTH,1,NR)=0.D0
      END DO

      DO NP=2,NPMAX+1
         RGAMA=SQRT(1.D0+PG(NP)**2*TMC2FP0)
         RUFP = (PTFP0*PG(NP))/AMFP
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
            DCPT(NTH,NP,NR) = DCPT(NTH,NP,NR)
     &           +FACT*RGAMA/RUFP
     &           *( 4.D0/VC**2*SUMA - SUMB 
     &           - 4.D0/(RUFP*VC**2)*SUMC + SUMD/RUFP )
         END DO
      END DO

      DO NP=1,NPMAX
         RGAMA=SQRT(1.D0+PM(NP)**2*TMC2FP0)
         RUFP = (PTFP0*PM(NP))/AMFP
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
            DCTP(NTH,NP,NR) = DCTP(NTH,NP,NR)
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
      RGAMMA= sqrt(1.D0+(RUFP/VC)**2)
C
C---- FIRST KIND LEGENDRE FUNCTION FOR J
C
      RZ=RUFP/VC
      RSIGMA=LOG(RZ+SQRT(1.D0+RZ**2))
      FKLF_J(0,0) = sqrt(2.D0*VC/PI/RUFP)
     &     *RSIGMA
      FKLF_J(0,1) = sqrt(2.D0*RUFP/PI/VC)
      FKLF_J(0,2) = sqrt(2.D0*RUFP/PI/VC)*RGAMMA
      FKLF_J(0,-1) = FKLF_J(0,1)
C
C---- FIRST KIND LEGENDRE FUNCTION FOR Y
C
      FKLF_Y(0,0) = sqrt(2.D0*VC/PI/RUFP)
      FKLF_Y(0,1) = sqrt(2.D0*VC/PI/RUFP)*RGAMMA
      FKLF_Y(0,2) = sqrt(2.D0*RUFP/PI/VC)
     &     *(VC/RUFP+2.D0*RUFP/VC)
      FKLF_Y(0,-1) = FKLF_Y(0,1)
C
C---- RECURRENCE EQUATION 
C
      Do L = 0, LLMAX+1
         Do NA = 0, 2
            FKLF_J(L+1,NA) = (RGAMMA*FKLF_J(L,NA)-FKLF_J(L,NA-1))
     &           *VC/RUFP/DBLE(NA+L+1)
            FKLF_Y(L+1,NA) = 
     &           ( (NA-L-1)*RGAMMA*FKLF_Y(L,NA)-(NA+L)*FKLF_Y(L,NA-1) )
     &           *VC/RUFP
         END DO
         FKLF_J(L+1,-1)=FKLF_J(L+1,1)
         FKLF_Y(L+1,-1)=FKLF_Y(L+1,1)
      END DO

C
C---- MINUS L
C
c      Do NA = -1,2
c         FKLF_J(-1,NA) = - FKLF_Y(0,NA)
c         FKLF_Y(-1,NA) = FKLF_J(0,NA)
c         FKLF_J(-2,NA) = FKLF_Y(1,NA)
c         FKLF_Y(-2,NA) = - FKLF_J(1,NA)
c      END DO


c      DO L=0,2
c         write(*,935) (FKLF_J(L,NA), NA=-1,2)
c      END DO
c      write(*,*) RUFP/VC
C
C--- TRANSFORM P TO J OR Y
C

      DO NA=0,2
         Do L=0,LLMAX+2
             RJ_1(L,NA)=SQRT(PI*VC/2.D0/RUFP)*FKLF_J(L,NA)
             RY_1(L,NA)=SQRT(PI*VC/2.D0/RUFP)*FKLF_Y(L,NA)*(-1)**(-L-1)
         END DO
         RJ_1(-1,NA) = - RY_1(0,NA)
         RY_1(-1,NA) = RJ_1(0,NA)
         RJ_1(-2,NA) = RY_1(1,NA)
         RY_1(-2,NA) = - RJ_1(1,NA)
      END DO

c      DO L=0,2
c         write(*,935) (RJ_1(L,NA), NA=0,2)
c      END DO
c      write(*,*) RUFP/VC

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
      RGAMMA= sqrt(1.D0+(RUFP/VC)**2)

      DO NA = 0, 2
         DO L= -1, LLMAX+2
           DERJ(L,NA) = RJ_1(L-1,NA)/(VC*RGAMMA) - (L+1)/RUFP*RJ_1(L,NA)
         END DO
         DO L= -2, LLMAX+1
           DERY(L,NA) =-RY_1(L+1,NA)/(VC*RGAMMA) + L/RUFP*RY_1(L,NA)
         END DO
      END DO


      RETURN
      END

