C     $Id$
C
C *****************************
C     PREPARATION OF FPLOOP
C *****************************

      SUBROUTINE FPPREP(IERR)

      INCLUDE 'fpcomm.inc'

      EXTERNAL FPFN0U, FPFN0T, FPFN1A, FPFN2A

C     ----- NS_NSA and NS_NSB array -----

      NSAMAX=NSFPMA-NSFPMI+1
      DO NSA=1,NSAMAX
         NS_NSA(NSA)=NSFPMI+NSA-1
      ENDDO
      NSBMAX=NSMAX
      DO NSB=1,NSBMAX
         NS_NSB(NSB)=NSB
      ENDDO
      DO NSA=1,NSAMAX
         WRITE(6,'(A,2I3)') 'NSA,NS=',NSA,NS_NSA(NSA)
      ENDDO
      DO NSB=1,NSBMAX
         WRITE(6,'(A,2I3)') 'NSB,NS=',NSB,NS_NSB(NSB)
      ENDDO

      CALL FPMESH(IERR)

C     ----- Initialize velocity distribution function of all species -----

      DO NS=1,NSMAX
         DO NR=1,NRMAX
            DO NP=1,NPMAX
               FL=FPMXWL(PM(NP),NR,NS)
               DO NTH=1,NTHMAX
                  FNS(NTH,NP,NR,NS)=FL
               END DO
            ENDDO
         END DO
      END DO

C     ----- set boundary distribution functions -----

      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         DO NP=1,NPMAX
            FL=FPMXWL(PM(NP),0,NS)
            DO NTH=1,NTHMAX
               FS1(NTH,NP,NSA)=FL
            ENDDO
         ENDDO
C
         DO NP=1,NPMAX
            FL=FPMXWL(PM(NP),NRMAX+1,NS)
            DO NTH=1,NTHMAX
               FS2(NTH,NP,NSA)=FL
            ENDDO
         ENDDO
      ENDDO

C     ----- set parameters for target species -----

      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         AEFP(NSA)=PZ(NS)*AEE
         AMFP(NSA)=PA(NS)*AMP
         RNFP0(NSA)=PN(NS)
         RNFPS(NSA)=PNS(NS)
         RTFP0(NSA)=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
         RTFPS(NSA)=PTS(NS)
C
         PTFP0(NSA)=SQRT(RTFP0(NSA)*1.D3*AEE*AMFP(NSA))
         VTFP0(NSA)=SQRT(RTFP0(NSA)*1.D3*AEE/AMFP(NSA))
      ENDDO

      DO NSB=1,NSBMAX
         NS=NS_NSB(NSB)
         AEFD(NSB)=PZ(NS)*AEE
         AMFD(NSB)=PA(NS)*AMP
         RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
C
         PTFD0(NSB)=SQRT(RTFD0L*1.D3*AEE*AMFD(NSB))
         VTFD0(NSB)=SQRT(RTFD0L*1.D3*AEE/AMFD(NSB))
      ENDDO
C
C     ----- set profile data -----
C
      DO NR=1,NRMAX
C
         RHON=RM(NR)
         CALL PLPROF(RHON)
C
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            RNFP(NR,NSA)=RN(NS)
            RTFP(NR,NSA)=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
            PTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE*AMFP(NSA))
            VTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE/AMFP(NSA))
         ENDDO

         RNE=RN(1)
         RTE=(RTPR(1)+2.D0*RTPP(1))/3.D0

         DO NSB=1,NSBMAX
            NS=NS_NSB(NSB)
            AEFD(NSB)=PZ(NS)*AEE
            AMFD(NSB)=PA(NS)*AMP
            RNFD(NR,NSB)=RN(NS)
            RTFD(NR,NSB)=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
            PTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE*AMFD(NSB))
            VTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE/AMFD(NSB))
         ENDDO

         DO NSA=1,NSAMAX
            NSFP=NS_NSA(NSA)
            DO NSB=1,NSBMAX
               NS=NS_NSB(NSB)

               IF(NSFP.EQ.1.AND.NS.EQ.1) THEN
                  RLNRL=14.9D0-0.5D0*LOG(RNE)+LOG(RTE)
               ELSEIF(NSFP.EQ.1.OR.NS.EQ.1) THEN
                  RLNRL=15.2D0-0.5D0*LOG(RNE)+LOG(RTE)
               ELSE
                  RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTFD(NR,NS))
               ENDIF
               FACT=AEFP(NSA)**2*AEFD(NSB)**2*RLNRL/(4.D0*PI*EPS0**2)
               RNUD(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20
     &                 /(SQRT(2.D0)*VTFD(NR,NSB)*PTFP0(NSA)**2)
               RNUF(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20
     &                 /(2*AMFD(NSB)*VTFD(NR,NSB)**2*PTFP0(NSA))
            ENDDO
         ENDDO
      ENDDO
c     ----- normalize f with respect to n -----

      NTEST=0
      IF(NTEST.eq.1)THEN
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            NR=1
            RSUM=0.D0
            DO  NP=1,NPMAX
               DO  NTH=1,NTHMAX
                  RSUM11=RSUM11
     &                 +VOL(NTH,NP)*RLAMDA(NTH,NR)*FNS(NTH,NP,NR,NS)
               END DO
            END DO

            WRITE(*,'(A,I5,1P2E12.4)') 'NS,RNFP,RSUM=',
     &                  NS,RNFP(1,NSA),RSUM11

            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS(NTH,NP,NR,NS)=FNS(NTH,NP,NR,NS)
     &                    *RNFP(NR,NSA)/RSUM11
                  END DO
               END DO
            END DO
         END DO
      END IF

C
C     ----- set relativistic parameters -----
C
      IF (MODELR.EQ.0) THEN
C
         DO NSA=1,NSAMAX
            THETA0(NSA)=0.D0
            DO NR=1,NRMAX
               THETA(NR,NSA)=0.D0
               DKBSR(NR,NSA)=0.D0
            ENDDO
         ENDDO
C
      ELSE
C
         DO NSA=1,NSAMAX
            THETA0(NSA)=RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)
            DO NR=1,NRMAX
               THETA(NR,NSA)=THETA0(NSA)*RTFP(NR,NSA)/RTFP0(NSA)
               Z=1.D0/THETA(NR,NSA)
               IF(Z.LE.100.D0) THEN
                  DKBSR(NR,NSA)=BESKN(2,Z)
               ELSE
                  DKBSR(NR,NSA)=SQRT(PI/(2.D0*Z))*EXP(-Z)
     &             *( 1.D0 + 15.D0/(8.D0*Z) + 105.D0/(128.D0*Z**2) )
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL FPCOEF
      CALL FPSGLB
      CALL FPWRT2
      CALL FPSPRF
      CALL FPWRT1

      IERR=0
      RETURN
      END

C     ***** create mesh quantities *****

      SUBROUTINE FPMESH(IERR)
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER*(80) LINE
      EXTERNAL FPFN0U,FPFN0T,FPFN1A,FPFN2A,FPFN2U,FPFN2T
C
C     ----- exec EQ -----
C
      IF(MODELG.EQ.3) THEN
         CALL EQLOAD(MODELG,KNAMEQ,IERR)
         IF(IERR.EQ.0) THEN
            write(LINE,'(A,I5)') 'nrmax=',51
            call eqparm(2,line,ierr)
            write(LINE,'(A,I5)') 'nthmax=',64
            call eqparm(2,line,ierr)
            write(LINE,'(A,I5)') 'nsumax=',64
            call eqparm(2,line,ierr)
            CALL EQCALQ(IERR)
            CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
         ELSE
            write(6,*) 'XX FPMESH:EQLOAD:IERR=',IERR
         ENDIF
      ENDIF
C
C     ----- set radial mesh -----
C
      IF(NRMAX.EQ.1) THEN
         DELR=DELR1
      ELSE
         DELR=(RMAX-RMIN)/NRMAX
      ENDIF
C
      IF(NRMAX.EQ.1) THEN
         RM(1)=R1
         RG(1)=R1-0.5D0*DELR
         RG(2)=R1+0.5D0*DELR
      ELSE
         DO 10 NR=1,NRMAX
            RM(NR)=RMIN+DELR*(NR-1)+0.5D0*DELR
            RG(NR)=RMIN+DELR*(NR-1)
   10    CONTINUE
         RG(NRMAX+1)=RMAX
      ENDIF
C
C     ----- load WR resluts -----
C
      IF(MODELW.EQ.1.OR.MODELW.EQ.2) THEN
         CALL FPLDWR(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF
      IF(MODELW.EQ.4) THEN
         CALL FPWMREAD(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF

C
C     ----- set poloidal magneticl field -----
C
      DO NR=1,NRMAX+1
         RHON=RG(NR)
         CALL FPSETB(RHON,0.5D0*PI,BT,BP(NR))
         EPSR(NR)=RSRHON(RHON)/RR
c         WRITE(*,*) BT, BP(NR)
      ENDDO
C
C     ----- set parallel current density -----
C
      DO NR=1,NRMAX
         RJ1(NR)=(RG(NR+1)*BP(NR+1)-RM(NR)*BP(NR))
     &          /(RMU0*RM(NR)*DELR)
      ENDDO
C
C     ----- set parallel electric field -----
C
      DO NR=1,NRMAX
         E1(NR)=E0
         E2(NR)=E0
         RJ2(NR)=RJ1(NR)
      ENDDO
C
C     ----- set momentum space mesh -----
C
      DELP =PMAX/NPMAX
      DELTH=PI/NTHMAX
C
      DO NP=1,NPMAX
        PG(NP)=DELP*(NP-1)
        PM(NP)=DELP*(NP-0.5D0)
      ENDDO
      PG(NPMAX+1)=PMAX
C
      DO NTH=1,NTHMAX
         THG(NTH)=DELTH*(NTH-1)
         THM(NTH)=DELTH*(NTH-0.5D0)
C
         SINM(NTH)=SIN(THM(NTH))
         COSM(NTH)=COS(THM(NTH))
         SING(NTH)=SIN(THG(NTH))
         COSG(NTH)=COS(THG(NTH))
      ENDDO
      THG(NTHMAX+1)=PI
      SING(NTHMAX+1)=0.D0
      COSG(NTHMAX+1)=-1.D0
C
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         VOL(NTH,NP)=2.D0*PI*SINM(NTH)*PM(NP)**2*DELP*DELTH
      ENDDO
      ENDDO

C
C     ----- set bounce-average parameters -----
C
      IF (MODELA.NE.0) THEN
         DO NR=1,NRMAX
            A1=ACOS(SQRT(2.D0*EPSR(NR)/(1.D0+EPSR(NR))))
            DO NTH=1,NTHMAX/2
               IF (THG(NTH).LE.A1.AND.THG(NTH+1).GE.A1) GOTO 201
            ENDDO
C
  201       CONTINUE
C
C            WRITE(6,*) 'NR,NTHC=',NR,NTH
            IF(NTH.EQ.NTHMAX/2) NTH=NTHMAX/2-1
C
            ITL(NR)=NTH
            ITU(NR)=NTHMAX-NTH+1
            EPSL=COSM(ITL(NR))**2/(2.D0-COSM(ITL(NR))**2)
C            WRITE(6,'(A,1P2E12.4)') 'EPSR(NR)=',EPSR(NR),EPSL
C            WRITE(6,'(A,2I8)') 'ITL,ITU=',ITL(NR),ITU(NR)
            EPSR(NR)=EPSL
            FACT=(1.D0+EPSL)/(2.D0*EPSL)
C
            DO NTH=1,ITL(NR)
               ETAM(NTH,NR)=PI/2.D0
            ENDDO
C
            DO NTH=ITL(NR)+1,ITU(NR)-1
               A1=FACT*COSM(NTH)**2
               ETAM(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
            ENDDO
C
            DO NTH=ITU(NR),NTHMAX
               ETAM(NTH,NR)=PI/2.D0
            ENDDO
C
            DO NTH=1,ITL(NR)
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
C
            DO NTH=ITL(NR)+1,ITU(NR)
               A1=FACT*COSG(NTH)**2
               ETAG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
            ENDDO
C
            DO NTH=ITU(NR)+1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
C
            NRX=NR
            DO NTH=1,NTHMAX/2
               NTHX=NTH
               IF(NTH.LT.ITL(NR)) THEN
                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0U)
C                  CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2U)
               ELSEIF(NTH.GT.ITL(NR)) THEN
                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0T)
C                  CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2T)
               ELSE
                  RINT0=0.D0
               ENDIF
               CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A)
               RLAMDA(NTH,NR)=RINT0*ABS(COSM(NTH))/PI
C               RLAMDC(NTH,NR)=RINT2/(PI*(1.D0+EPSR(NR))*ABS(COSG(NTH)))
               RLAMDC(NTH,NR)=RINT2/(PI*(1.D0+EPSR(NR))*(COSG(NTH)))
            ENDDO
            RLAMDA(ITL(NR),NR)=0.5D0*(RLAMDA(ITL(NR)-1,NR)
     &                               +RLAMDA(ITL(NR)+1,NR))
            DO NTH=1,NTHMAX/2
               RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
               RLAMDC(NTHMAX-NTH+2,NR)=RLAMDC(NTH,NR)
            ENDDO
            RLAMDC(NTHMAX/2+1,NR)=0.D0
         ENDDO
C
       ELSE
         DO NR=1,NRMAX
            ITL(NR)=0
            ITU(NR)=0
            DO NTH=1,NTHMAX
               ETAM(NTH,NR)=PI/2.D0
               RLAMDA(NTH,NR)=1.D0
            ENDDO
            DO NTH=1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
         ENDDO
      END IF

      IERR=0
      RETURN
      END
C
C ****************************************
C     MAXWELLIAN VELOCITY DISTRIBUTION
C ****************************************
C
      FUNCTION FPMXWL(PML,NR,NS)
C
      INCLUDE 'fpcomm.inc'
C
      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.EQ.0) THEN
         RL=RM(1)-DELR
         RHON=RL
         CALL PLPROF(RHON)
         RNFDL=RN(NS)/RNFD0L
         RT=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
         RTFDL=RT/RTFD0L
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=RL
         CALL PLPROF(RHON)
         RNFDL=RN(NS)/RNFD0L
         RT=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
         RTFDL=RT/RTFD0L
      ELSE
         RL=RM(NR)
         RHON=RL
         CALL PLPROF(RHON)
         RNFDL=RN(NS)/RNFD0L
         RTFDL=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
      ENDIF

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
         IF(Z.LE.150.D0) THEN
            DKBSL=BESKN(2,Z)
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL)
     &        *RTFD0L*EXP(-1.D0/THETAL)
C            FACT=RNFDL*SQRT(THETA0**3)/(4.D0*PI*THETAL*DKBSL)
C     &           *EXP(-1.D0/THETAL)
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         ELSE
            DAPPROX = SQRT( pi/(2.D0*Z) )*
     &           ( 1.D0+15.D0/(8.D0*Z) +15.D0*7.D0/(8.D0*Z)**2/2.D0 )
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DAPPROX)
     &        *RTFD0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL

         ENDIF
         IF(EX.LT.-100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(EX)
         ENDIF
      ENDIF

      RETURN
      END
c-----------------------
      FUNCTION FPBEAM(PML,NR,NS)
      INCLUDE 'fpcomm.inc'

      AMFD=PA(NS)*AMP
      AEFD=PZ(NS)*AEE
c      RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
c      PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)
c      RNFD0=PN(NS)
c      RNFDL=RNFD(NR,NS)
      rsigma=1.d-1

c      IF(NS.eq.NSBM)THEN
c      EX=-(PML-1.D0)**2/(2.D0*rsigma**2)
c      FACT=RNFDL*PTFD0**3/(SQRT(2.D0*PI)*rsigma)
c      FPBEAM=FACT*EXP(EX)
c      END IF
      FPBEAM=RNFDL/(4.D0*PI)

      RETURN
      END

C
C *************************
C     INITIAL DATA SAVE
C *************************
C
      SUBROUTINE FPCINI
C
      INCLUDE 'fpcomm.inc'
C
      ISAVE=0
C
      DO NSA=1,NSAMAX
      DO NR=1,NRMAX
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
      ENDDO

      RETURN
      END
C
C ***************************************************************
C
C                       SET OF INTEGRAND
C
C ***************************************************************
C
C ============================================================
C
      REAL*8 FUNCTION  FPFN0U(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=(1.D0+EPSR(NRX))*SINM(NTHX)**2
      FPFN0U=A0*SQRT(A1/(A1-A2))
C
      RETURN
      END
C
C ============================================================
C
      REAL*8 FUNCTION FPFN0T(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=2.D0*EPSR(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN0T=A0*SQRT(A1/A2)
C
      RETURN
      END
C
C ============================================================
C
      REAL*8 FUNCTION FPFN1A(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      FPFN1A=A0*SQRT(A1)
C
      RETURN
      END
C
C ============================================================
C
      REAL*8 FUNCTION FPFN2A(X,XM,XP)
C
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSR(NRX))*SING(NTHX)**2
      FPFN2A=A0*SQRT(A1*A2)
C
      RETURN
      END
C--------------------------------
      REAL*8 FUNCTION  FPFN2U(X,XM,XP)
C                           
      INCLUDE 'fpcomm.inc'
C
      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=(1.D0+EPSR(NRX))*SING(NTHX)**2
      FPFN2U=A0*SQRT(A1/(A1-A2))
C
      RETURN
      END
C------------------------------
      REAL*8 FUNCTION FPFN2T(X,XM,XP)
C                        
      INCLUDE 'fpcomm.inc'
C
      XX=X
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSR(NRX)*COS(A0*XP)
      A2=2.D0*EPSR(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN2T=A0*SQRT(A1/A2)
C
      RETURN
      END
C
C *****************
C     MAIN LOOP
C *****************
C
      SUBROUTINE FPLOOP
C
      INCLUDE 'fpcomm.inc'
C
      DIMENSION RJNS(NRM,NSAM),RJN(NRM),RJ3(NRM),E3(NRM),DELE(NRM)
      
C
      IF(MODELE.NE.0) CALL FPNEWE

C     +++++ Time loop +++++

      DO NT=1,NTMAX

C     +++++ Iteration loop for toroidal electric field +++++

         L=0
C
         IF(MODELE.NE.0) THEN
            DO NR=1,NRMAX
               E3(NR)=0.D0
               RJ3(NR)=0.D0
            ENDDO
         ENDIF
C
    1    L=L+1

C     +++++ NSA loop +++++

         IF (MOD(NT-1,NTSTPC).EQ.0.AND.NT.NE.1) CALL FPCOEF

         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     F(NTH,NP,NR)=FNS(NTH,NP,NR,NS)
                  END DO
               END DO
            END DO

            CALL FPEXEC(NSA,IERR)
            IF(IERR.NE.0) GOTO 250

            DO NR=1,NRMAX
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     FNS1(NTH,NP,NR,NS)=F1(NTH,NP,NR) 
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

C     +++++ end of NSA loop +++++

C     ----- calculation of current density -----

         IF(MODELE.NE.0) THEN
            DO NSA=1,NSAMAX
               NS=NS_NSA(NSA)
               DO NR=2,NRMAX
                  RSUM=0.D0
                  DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM=RSUM+VOL(NTH,NP)*FNS1(NTH,NP,NR,NS)*PM(NP)
                  ENDDO
                  ENDDO
                  RJNS(NR,NSA)=AEFP(NSA)*RNFP0(NSA)*1.D20
     &                    *PTFP0(NSA)*DELP*RSUM/(AMFP(NSA)*RM(NR)*RA)
               ENDDO
               RJNS(1,NSA)=(4.D0*RJNS(2,NSA)-RJNS(3,NSA))/3.D0
            ENDDO

C     ----- calculation of toroidal electric field -----

            DELEM=0.D0
            DO NR=1,NRMAX
               RJNL=0.D0
               DO NSA=1,NSAMAX
                  RJNL=RJNL+RJNS(NR,NSA)
               END DO
               RJN(NNR)=RJNL
               IF(ABS(RJNL-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJNL-RJ2(NR))*(E2(NR)-E3(NR))
     &                    /(RJNL-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJNL
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
            ENDDO

            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 1
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'

            DO NR=1,NRMAX
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
            ENDDO
            CALL FPNEWE
         ENDIF
C     +++++ end of toroidal electric field loop +++++
C
  250    CONTINUE
C
C         CALL FPGRAC('F -2',F,4)
C         CALL FPGRAC('F1-2',F1,4)

C     +++++ update velocity distribution function +++++

         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            FNS(NTH,NP,NR,NS)=FNS1(NTH,NP,NR,NS)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

C     +++++ calculate and save global data +++++

         TIMEFP=TIMEFP+DELT

         ISAVE=0
         IF (MOD(NT,NTSTP2).EQ.0) THEN
            CALL FPSGLB
            CALL FPWRT2
         ENDIF
         IF (MOD(NT,NTSTP1).EQ.0) THEN
            CALL FPSPRF
            CALL FPWRT1
         ENDIF

c      call FPWRT3

         IF(IERR.NE.0) RETURN
      ENDDO
C     +++++ end of time loop +++++
C
      RETURN
      END
C
C ************************************
C     PREDICTION OF ELECTRIC FIELD
C ************************************
C

      SUBROUTINE FPNEWE
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=2,NRMAX
         BP(NR)=BP(NR)+(E1(NR)-E1(NR-1))*DELT/(RA*DELR)
      ENDDO
C
      DO NR=1,NRMAX
         RJ2(NR)=(RG(NR+1)*BP(NR+1)-RG(NR)*BP(NR))
     &           /(RMU0*RM(NR)*DELR*RA)
      ENDDO
C
      DO NR=1,NRMAX
         E2(NR)=RJ2(NR)*E1(NR)/RJ1(NR)
      ENDDO
C
      RETURN
      END
