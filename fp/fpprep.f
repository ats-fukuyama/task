C     $Id$
C
C *****************************
C     PREPARATION OF FPLOOP
C *****************************

      SUBROUTINE FPPREP(IERR)

      USE plprof
      INCLUDE 'fpcomm.inc'

      EXTERNAL FPFN0U, FPFN0T, FPFN1A, FPFN2A
c      DIMENSION RCOEF(NRMAX)
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

c--------- normalize bounce averaged distribution function ---------

      IF(MODELA.eq.1)THEN
         DO NS=1,NSMAX
         DO NR=1,NRMAX
         RSUM1=0.D0
         RSUM2=0.D0
         RCOEF(NR)=0.D0
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM1 = RSUM1+VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
     &                 *RLAMDA(NTH,NR)
                  RSUM2 = RSUM2+VOL(NTH,NP)*FNS(NTH,NP,NR,NS)
               END DO
            END DO
            RCOEF(NR)=RSUM2/RSUM1
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FNS(NTH,NP,NR,NS) = FNS(NTH,NP,NR,NS)*RCOEF(NR)
               END DO
            END DO

         END DO
         END DO
      END IF
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
         CALL pl_prof(RHON,PLF)
C
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            RNFP(NR,NSA)=PLF(NS)%RN
            RTFP(NR,NSA)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
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
            ENDDO
         ENDDO
C
      ELSE
C
         DO NSA=1,NSAMAX
            THETA0(NSA)=RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)
            DO NR=1,NRMAX
               THETA(NR,NSA)=THETA0(NSA)*RTFP(NR,NSA)/RTFP0(NSA)
            ENDDO
         ENDDO
      ENDIF

c$$$      DO NSA=1,NSAMAX
c$$$         NS=NS_NSA(NSA)
c$$$      DO NR=1,NRMAX
c$$$         NTH=1
c$$$         DO NP=2,NPMAX-1
c$$$            IF(FNS(NTH,NP+1,NR,NSA).LE.1.D-70) THEN
c$$$               DFDP=0.D0
c$$$            ELSE
c$$$               DFDP=(FNS(NTH,NP+1,NR,NS)-FNS(NTH,NP-1,NR,NS))
c$$$     &                    /(2.D0*DELP*FNS(NTH,NP,NR,NS))
c$$$            ENDIF
c$$$            IF(FNS(NTH,NP,NR,NSA).LE.1.D-70) THEN
c$$$               DFDP1=0.D0
c$$$            ELSE
c$$$               DFDP1=(LOG(FNS(NTH,NP+1,NR,NS))-LOG(FNS(NTH,NP-1,NR,NS)))
c$$$     &                    /(2.D0*DELP)
c$$$            ENDIF
c$$$            DFDP2=-PM(NP)*RTFP0(NSA)/RTFP(NR,NSA)
c$$$            write(6,'(I5,1P5E12.4)')
c$$$     &           NP,PM(NP),DFDP,DFDP1,DFDP2,FNS(NTH,NP,NR,NS)
c$$$         ENDDO
c$$$      ENDDO
c$$$      ENDDO
      DO NSA=1,NSAMAX
         CALL FPCOEF(NSA)
      END DO
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
      DO NSA=NSFPMI,NSFPMA
      IF(MODELW(NSA).EQ.1.OR.MODELW(NSA).EQ.2) THEN
         CALL FPLDWR(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF
      IF(MODELW(NSA).EQ.4) THEN
         CALL FPWMREAD(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF
      ENDDO

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
      DO NR=1,NRMAX
         RHOL=RM(NR)
         RHOL1=RG(NR)
         RHOL2=RG(NR+1)
         VOLR(NR)=2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1))
     &           *2.D0*PI*RR
      ENDDO
      TVOLR=0.D0
      DO NR=1,NRMAX
         TVOLR=TVOLR+VOLR(NR)
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
            WRITE(6,*) 'NR,NTHC=',NR,NTH
            IF(NTH.EQ.NTHMAX/2) NTH=NTHMAX/2-1
C
            ITL(NR)=NTH
            ITU(NR)=NTHMAX-NTH+1
            EPSL=COSM(ITL(NR))**2/(2.D0-COSM(ITL(NR))**2)
c            EPSL=EPSR(NR)
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
                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0U,'FPFN0U')
C                  CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2U,'FPFN2U')
               ELSEIF(NTH.GT.ITL(NR)) THEN
                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0T,'FPFN0T')
C                  CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2T,'FPFN2T')
               ELSE
                  RINT0=0.D0
               ENDIF
               CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A,'FPFN2A')
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
ccc
c            DO NTH=1,NTHMAX
c               DELH=2.D0*ETAM(NTH,NR)/NAVMAX
c               sum11=0.D0  
c               DO NG=1,NAVMAX
c                  ETAL=DELH*(NG-0.5D0)
c                  X=EPSR(NR)*COS(ETAL)*RR
c                  PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
c                  PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
c                  sum11=sum11
c     &                 +1.D0/PCOS*ABS(COSM(NTH))
c               END DO
c               RLAMDA2(NTH,NR)=SUM11*DELH/PI
c               WRITE(*,*) NTH,RLAMDA(NTH,NR),RLAMDA2(nth,nr)
c ETAM(NTH,NR)*2.D0,ETAG(NTH,NR)*2.D0
c 
c            END DO
c            RLAMDA(ITL(NR),NR)=0.5D0*(RLAMDA(ITL(NR)-1,NR)
c     &                               +RLAMDA(ITL(NR)+1,NR))
c            RLAMDA(ITU(NR),NR)=0.5D0*(RLAMDA(ITU(NR)-1,NR)
c     &                               +RLAMDA(ITU(NR)+1,NR))
c            RLAMDA2(ITL(NR),NR)=0.5D0*(RLAMDA2(ITL(NR)-1,NR)
c     &                               +RLAMDA2(ITL(NR)+1,NR))
c            RLAMDA2(ITU(NR),NR)=0.5D0*(RLAMDA2(ITU(NR)-1,NR)
c     &                               +RLAMDA2(ITU(NR)+1,NR))

ccc            

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
      USE plprof
      INCLUDE 'fpcomm.inc'
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
C
      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.EQ.0) THEN
         RL=RM(1)-DELR
         RHON=RL
         CALL PL_PROF(RHON,PLF)
         RNFDL=RN(NS)/RNFD0L
         RTFDL=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=RL
         CALL PL_PROF(RHON,PLF)
         RNFDL=PLF(NS)%RN/RNFD0L
         RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
      ELSE
         RL=RM(NR)
         RHON=RL
         CALL PL_PROF(RHON,PLF)
         RNFDL=PLF(NS)%RN/RNFD0L
         RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
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
c$$$         IF(Z.LE.150.D0) THEN
            DKBSL=BESEKN(2,Z)
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL)
     &        *RTFD0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
c$$$         ELSE
c$$$            DAPPROX = SQRT( pi/(2.D0*Z) )*
c$$$     &           ( 1.D0+15.D0/(8.D0*Z) +15.D0*7.D0/(8.D0*Z)**2/2.D0 )
c$$$            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DAPPROX)
c$$$     &        *RTFD0L
c$$$            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
c$$$
c$$$         ENDIF
         IF(EX.LT.-100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(EX)
         ENDIF
      ENDIF

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
