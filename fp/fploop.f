C     $Id$
C
C *****************************
C     PREPARATION OF FPLOOP
C *****************************
C
      SUBROUTINE FPPREP(IERR)
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
            write(6,*) 'XX FPPREP:EQLOAD:IERR=',IERR
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
C     ----- set parameters for target species -----
C
      AEFP=PZ(NSFP)*AEE
      AMFP=PA(NSFP)*AMP
      RNFP0=PN(NSFP)
      RNFPS=PNS(NSFP)
      RTFP0=(PTPR(NSFP)+2.D0*PTPP(NSFP))/3.D0
      RTFPS=PTS(NSFP)
C
      PTFP0=SQRT(RTFP0*1.D3*AEE*AMFP)
      VTFP0=SQRT(RTFP0*1.D3*AEE/AMFP)
C
C     ----- set profile data -----
C
      DO NR=1,NRMAX
C
         RHON=RM(NR)
         CALL PLPROF(RHON)
C
         RNFP(NR)=RN(NSFP)
         RTFP(NR)=(RTPR(NSFP)+2.D0*RTPP(NSFP))/3.D0
         NTG3 = INT((NTG2+1)/(NSFPMA-NSFPMI+1)) + 1
         if(NSFPMA.eq.NSFPMI)NTG3 = INT((NTG2+1)/(NSFPMA-NSFPMI+1))
         IF(NTG3.eq.1)THEN
            DO NS=1,NSMAX
c               PTT2(1,NS)=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
            END DO
         END IF

         IF(NTG2.ge.NSFPMA+1)THEN
            IF(MODELC.eq.1.or.MODELC.eq.3)THEN
               RTFP(NR)=PTT2(NTG3-1,NSFP)
     &              *(RTPR(NSFP)+2.D0*RTPP(NSFP))/(3.D0*PTT2(1,NSFP))
            END IF
         END IF
c         write(*,*) "temperature",PTT(1), RTFP(1)

         PTFP(NR)=SQRT(RTFP(NR)*1.D3*AEE*AMFP)
         VTFP(NR)=SQRT(RTFP(NR)*1.D3*AEE/AMFP)
         RNE=RN(1)
         RTE=(RTPR(1)+2.D0*RTPP(1))/3.D0
         DO NS=1,NSMAX
            AEFD=PZ(NS)*AEE
            AMFD=PA(NS)*AMP
            RNFD(NR,NS)=RN(NS)
            RTFD(NR,NS)=(RTPR(NS)+2.D0*RTPP(NS))/3.D0
            IF(NTG2.ge.NSFPMA+1)THEN
               IF(MODELC.eq.1.or.MODELC.eq.3)THEN
                  IF(NS.ge.NSFPMI.and.NS.le.NSFPMA)
     &                 RTFD(NR,NS)=PTT2(NTG3-1,NS)
     &              *(RTPR(NSFP)+2.D0*RTPP(NSFP))/(3.D0*PTT2(1,NSFP))
               END IF
            END IF
               PTFD(NR,NS)=SQRT(RTFD(NR,NS)*1.D3*AEE*AMFD)
               VTFD(NR,NS)=SQRT(RTFD(NR,NS)*1.D3*AEE/AMFD)
            IF(NSFP.EQ.1.AND.NS.EQ.1) THEN
               RLNRL=14.9D0-0.5D0*LOG(RNE)+LOG(RTE)
            ELSEIF(NSFP.EQ.1.OR.NS.EQ.1) THEN
               RLNRL=15.2D0-0.5D0*LOG(RNE)+LOG(RTE)
            ELSE
               RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTFD(NR,NS))
            ENDIF
            FACT=AEFP**2*AEFD**2*RLNRL/(4.D0*PI*EPS0**2)
            RNUD(NR,NS)=FACT*RNFP0*1.D20
     &                 /(SQRT(2.D0)*VTFD(NR,NS)*PTFP0**2)
            RNUF(NR,NS)=FACT*RNFP0*1.D20
     &                 /(2*AMFD*VTFD(NR,NS)**2*PTFP0)
         ENDDO
      ENDDO

C
C     ----- set poloidal magneticl field -----
C
      DO NR=1,NRMAX+1
         RHON=RG(NR)
         CALL FPSETB(RHON,0.5D0*PI,BT,BP(NR))
         EPSR(NR)=RSRHON(RHON)/RR
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
C     ----- set relativistic parameters -----
C
      IF (MODELR.EQ.0) THEN
C
         THETA0=0.D0
         DO NR=1,NRMAX
            THETA(NR)=0.D0
            DKBSR(NR)=0.D0
         ENDDO
C
      ELSE
C
         THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
         DO NR=1,NRMAX
            THETA(NR)=THETA0*RTFP(NR)/RTFP0
            Z=1.D0/THETA(NR)
            IF(Z.LE.100.D0) THEN
               DKBSR(NR)=BESKN(2,Z)
            ELSE
               DKBSR(NR)=SQRT(PI/(2.D0*Z))*EXP(-Z)
            ENDIF
         ENDDO
      ENDIF
C
C     ----- set boundary distribution functions -----
C
      DO NP=1,NPMAX
         FL=FPMXWL(PM(NP),0,NSFP)
         DO NTH=1,NTHMAX
            FS1(NTH,NP)=FL
         ENDDO
      ENDDO
C
      DO NP=1,NPMAX
         FL=FPMXWL(PM(NP),NRMAX+1,NSFP)
         DO NTH=1,NTHMAX
            FS2(NTH,NP)=FL
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
C
      IERR=0
      RETURN
      END
C
C ****************************************
C     INITIALIZE VELOCITY DISTRIBUTION    
C ****************************************
C
      SUBROUTINE FPFINI
C
      INCLUDE 'fpcomm.inc'
C
c      open(7,file='testdelt.dat')
      DO NR=1,NRMAX
      DO NP=1,NPMAX
         FL=FPMXWL(PM(NP),NR,NSFP)
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=FL
            DO NS=1,NSMAX
               FLNS=FPMXWL(PM(NP),NR,NS)
               FNS(NTH,NP,NR,NS)=FLNS
               FNS2(NTH,NP,NR,NS)=FLNS
            END DO
         ENDDO
c         WRITE(6,'(I5,1P4E12.4)') NP,FNS(1,NP,NR,1)
c     &        ,FNS(NTHMAX,NP,NR,2)
c     &        ,FNS(1,NP,NR,2),PM(NP)
      ENDDO
      ENDDO

c-----distribution of beam ion 
      IF(NSBM.ne.0)THEN
      sum1=0.D0
      sum2=0.D0
      sum3=0.D0
      Do NP=1,NPMAX
         sum1=sum1+FNS(1,NP,1,1)*PM(NP)**2*DELP
         sum2=sum2+FNS(1,NP,1,3)*PM(NP)**2*DELP
      END DO
      Do NP=1,NPMAX
         PMP=PM(NP)-1.D0
         FLNS=0.D0
         if(NP.ne.1)PMM=PM(NP-1)-1.D0
         if(PMP.gt.0.D0.and.PMM.lt.0.D0)then
            FLNS=sum2/PM(NP)**2
         END if
         Do NTH=1,NTHMAX
            FNS(NTH,NP,1,NSBM)=0
            if(NSFP.eq.NSBM) F(NTH,NP,1)=0
            if(NTH.eq.1)then
               FNS(NTH,NP,1,NSBM)=FLNS*2.D0/SINM(NTH)
               if(NSFP.eq.NSBM) F(NTH,NP,1)=FLNS*2.D0/SINM(NTH)
            end if
         END DO
c         sum3 = sum3 +FNS(1,NP,1,3)*PM(NP)**2
c         WRITE(6,'(I5,1P4E12.4)') NP,FNS(1,NP,1,1)
c     &        ,FNS(1,NP,1,2)
c     &        ,FNS(1,NP,1,2),PM(NP)
      END DO
c      write(*,*)"NSBM",NSBM
c      write(*,*)sum1,sum2,sum3
      END IF
      RETURN
      END
C
C ****************************************
C     MAXWELLIAN VELOCITY DISTRIBUTION
C ****************************************
C
C      FUNCTION FPMXWL(PML,NR)
      FUNCTION FPMXWL(PML,NR,NS)
C
      INCLUDE 'fpcomm.inc'
C
      AMFD=PA(NS)*AMP
      AEFD=PZ(NS)*AEE
      RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)
      RNFD0=PN(NS)

      IF(NR.EQ.0) THEN
         RL=RM(1)-DELR
         RHON=RL
         CALL PLPROF(RHON)
         RNFDL=RN(NS)/RNFD0
         RTFDL=RTPR(NS)/RTFD0
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=RL
         CALL PLPROF(RHON)
         RNFDL=RN(NS)/RNFD0
         RTFDL=RTPR(NS)/RTFD0
      ELSE
         RNFDL=RNFD(NR,NS)
         RTFDL=RTFD(NR,NS)
      ENDIF

      IF (MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0)
         IF(EX.GT.100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(-EX)
         ENDIF

      ELSE
         THETA0=RTFD0*1.D3*AEE/(AMFD*VC*VC)
         THETAL=THETA0*RTFDL/RTFD0
         Z=1.D0/THETAL
         IF(Z.LE.150.D0) THEN
            DKBSL=BESKN(2,Z)
            FACT=RNFDL*SQRT(THETA0)/(4.D0*PI*RTFDL*DKBSL)
     &        *RTFD0*EXP(-1.D0/THETAL)
C            FACT=RNFDL*SQRT(THETA0**3)/(4.D0*PI*THETAL*DKBSL)
C     &           *EXP(-1.D0/THETAL)
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0))/THETAL
         ELSE
            FACT=RNFDL/(4.D0*PI)*SQRT(2.D0/PI)*SQRT(THETA0/THETAL)**3
            EX=-THETA0*PML**2/THETAL*0.5D0
            DKBSL=SQRT(PI/(2.D0*Z))*EXP(-Z)
         ENDIF
         IF(EX.LT.-100.D0) THEN
            FPMXWL=0.D0
         ELSE
            FPMXWL=FACT*EXP(EX)
         ENDIF
C         if(PML.le.2.D0)
C     &        WRITE(6,'(I5,1P4E12.4)') NS,PML,Z,FPMXWL,DKBSL
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
      SUBROUTINE FPSAVI
C
      INCLUDE 'fpcomm.inc'
C
      ISAVE=0
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         DCPP(NTH,NP,NR)=0.D0
         DCPT(NTH,NP,NR)=0.D0
         FCPP(NTH,NP,NR)=0.D0
         DWPP(NTH,NP,NR)=0.D0
         DWPT(NTH,NP,NR)=0.D0
         FEPP(NTH,NP,NR)=0.D0
         DWLHPP(NTH,NP,NR)=0.D0
         DWLHPT(NTH,NP,NR)=0.D0
         DWFWPP(NTH,NP,NR)=0.D0
         DWFWPT(NTH,NP,NR)=0.D0
         DWECPP(NTH,NP,NR)=0.D0
         DWECPT(NTH,NP,NR)=0.D0
      ENDDO
      ENDDO
      ENDDO
      CALL FPSPRF
      CALL FPSGLB
      CALL FPWRIT
C

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
      DIMENSION RJN(NRM),RJ3(NRM),E3(NRM),DELE(NRM)
      
C
      IF(MODELE.NE.0) CALL FPNEWE
C
      DO NT=1,NTMAX
C
      IF(MODELC.ne.-3)THEN
C----NSFP LOOP----------------------------------
      DO NSFP=NSFPMI,NSFPMA
         CALL FPPREP(IERR)
         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=FNS(NTH,NP,NR,NSFP)
         END DO
         END DO
         END DO
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
C
         IF (MOD(NT-1,NTSTPC).EQ.0) CALL FPCOEF
C
         CALL FPEXEC(IERR)
         IF(IERR.NE.0) GOTO 250
C
         IF(MODELE.NE.0) THEN
            DO NR=2,NRMAX
               RSUM=0.D0
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM=RSUM+VOL(NTH,NP)*F1(NTH,NP,NR)*PM(NP)
               ENDDO
               ENDDO
               RJN(NR)=AEFP*RNFP0*1.D20*PTFP0*DELP*RSUM/(AMFP*RM(NR)*RA)
            ENDDO
            RJN(1)=(4.D0*RJN(2)-RJN(3))/3.D0
C
            DELEM=0.D0
            DO NR=1,NRMAX
               IF(ABS(RJN(NR)-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJN(NR)-RJ2(NR))*(E2(NR)-E3(NR))
     &                    /(RJN(NR)-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJN(NR)
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
            ENDDO
C
            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 1
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'
C
            DO NR=1,NRMAX
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
            ENDDO
            CALL FPNEWE
         ENDIF
C
  250    CONTINUE
C         CALL FPGRAC('F -2',F,4)
C         CALL FPGRAC('F1-2',F1,4)
         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=F1(NTH,NP,NR)
            FNS2(NTH,NP,NR,NSFP)=F1(NTH,NP,NR) 
            if(NSFP.eq.NSFPMA)THEN
               DO NS2=NSFPMI,NSFPMA
                  FNS(NTH,NP,NR,NS2)=FNS2(NTH,NP,NR,NS2)
               END DO
            END IF
         ENDDO
         ENDDO
         ENDDO

C
         if(NT.eq.1.and.NSFP.eq.1)TIMEFP=TIMEFP+DELT
C
         ISAVE=0
         IF (MOD(NT,NTSTP1).EQ.0) THEN
            CALL FPSPRF
         ENDIF
         IF (MOD(NT,NTSTP2).EQ.0) THEN
            CALL FPSGLB
            CALL FPWRIT
         ENDIF

      END DO
      NSFP=1
      TIMEFP=TIMEFP+DELT
C-----END OF NSFP LOOP--------------------------
      ELSE
         CALL FPPREP(IERR)

         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=FNS(NTH,NP,NR,NSFP)
         END DO
         END DO
         END DO

         L=0
C
         IF(MODELE.NE.0) THEN
            DO NR=1,NRMAX
               E3(NR)=0.D0
               RJ3(NR)=0.D0
            ENDDO
         ENDIF
C
 2       L=L+1
C
         IF (MOD(NT-1,NTSTPC).EQ.0) CALL FPCOEF
C
         CALL FPEXEC(IERR)
         IF(IERR.NE.0) GOTO 251
C
         IF(MODELE.NE.0) THEN
            DO NR=2,NRMAX
               RSUM=0.D0
               DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM=RSUM+VOL(NTH,NP)*F1(NTH,NP,NR)*PM(NP)
               ENDDO
               ENDDO
               RJN(NR)=AEFP*RNFP0*1.D20*PTFP0*DELP*RSUM/(AMFP*RM(NR)*RA)
            ENDDO
            RJN(1)=(4.D0*RJN(2)-RJN(3))/3.D0
C
            DELEM=0.D0
            DO NR=1,NRMAX
               IF(ABS(RJN(NR)-RJ3(NR)).GT.1.D-20) THEN
                  DELE(NR)=(RJN(NR)-RJ2(NR))*(E2(NR)-E3(NR))
     &                    /(RJN(NR)-RJ3(NR))
                  E3(NR)=E2(NR)
                  RJ3(NR)=RJN(NR)
                  E2(NR)=E2(NR)-DELE(NR)
                  DELEM=MAX(ABS(DELE(NR))/MAX(ABS(E1(NR)),1.D-6),DELEM)
               ENDIF
            ENDDO
C
            IF (L.LT.LMAXE.AND.DELEM.GT.EPSE) GO TO 2
            IF (L.GE.LMAXE) WRITE(6,*) 'L IS LARGER THAN LMAXE'
C
            DO NR=1,NRMAX
               E1(NR)=E2(NR)
               RJ1(NR)=RJN(NR)
            ENDDO
            CALL FPNEWE
         ENDIF
C
 251     CONTINUE
C         CALL FPGRAC('F -2',F,4)
C         CALL FPGRAC('F1-2',F1,4)
         DO NR=1,NRMAX
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX
            F(NTH,NP,NR)=F1(NTH,NP,NR)
            FNS2(NTH,NP,NR,NSFP)=F1(NTH,NP,NR) 
c     tentative FNS
c            DO NS=1,NSMAX
c            F1(NTH,NP,NR)=FNS(NTH,NP,NR,NSFP)
               FNS(NTH,NP,NR,NSFP)=F1(NTH,NP,NR) 
c            END DO
         ENDDO
         ENDDO
         ENDDO

         TIMEFP=TIMEFP+DELT
         ISAVE=0
C
C
         IF (MOD(NT,NTSTP1).EQ.0) THEN
            CALL FPSPRF
         ENDIF
         IF (MOD(NT,NTSTP2).EQ.0) THEN
            CALL FPSGLB
            CALL FPWRIT
         ENDIF

      END IF
C------------------------------
         IF(IERR.NE.0) GOTO 1100
      ENDDO
 1100 CONTINUE
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
