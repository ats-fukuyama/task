C     $Id$

C     ****** RADIAL MESH AND METRIC TENSOR ******

      subroutine wmfem_setg(ierr)

      INCLUDE 'wmcomm.inc'

C     ****** Cylindrical coordinates ******

      IF(MODELG.EQ.0) THEN
         CALL wmmetric_cyl(ierr)

C     ****** Toroidal coordinates (circular tokamak) ******

      ELSEIF(MODELG.EQ.1) THEN
         CALL wmmetric_tor(IERR)

C     ****** Toroidal coordinates (simple equilibrium) ******
      ELSEIF(MODELG.EQ.2) THEN
         CALL wmmetric_tor(IERR)

C     ****** equilibrium file (TASK/EQ) ******

      ELSEIF(MODELG.EQ.3) THEN
         CALL wmmetric_eq(IERR)

C     ****** equilibrium file (VMEC) ******

      ELSEIF(MODELG.EQ.4) THEN
C         CALL wmmetric_vmec(IERR)

C     ****** equilibrium file (EQDSK) ******

      ELSEIF(MODELG.EQ.5) THEN
         CALL wmmetric_eq(IERR)

C     ****** equilibrium file (BOOZER) ******

      ELSEIF(MODELG.EQ.6) THEN
C         CALL wmmetric_booz(IERR)

C     ****** bpsd equilibrium ******

      ELSEIF(MODELG.EQ.9) THEN
C         CALL wmmetric_bpsd(IERR)

      ENDIF
C
      IF(IERR.NE.0) RETURN
C
      IF(MODELN.EQ.7) CALL WMDPRF(IERR)
      IF(MODELN.EQ.8) CALL WMXPRF(IERR)
      IF(MODELN.EQ.9) THEN
         IF(KNAMTR.NE.KNAMTR_SAVE) THEN
C            CALL WMTRLOAD(KNAMTR,IERR)
            KNAMTR_SAVE=KNAMTR
         ENDIF
      ENDIF
      IF(IERR.NE.0) RETURN
C
      IF(NPHMAX.EQ.1) THEN
         NDSIZ  = 1
         NDMIN  = 0
         NDMAX  = 0
         KDSIZ  = 1
         KDMIN  = 0
         KDMAX  = 0
         NDSIZX = 1
      ELSE
         NDSIZ  = NPHMAX
         NDMIN  =-NPHMAX/2+1
         NDMAX  = NPHMAX/2
         KDSIZ  = NPHMAX
         KDMIN  =-NPHMAX/2+1
         KDMAX  = NPHMAX/2
         NDSIZX = 3*NPHMAX/2
      ENDIF
C
      IF(NTHMAX.EQ.1) THEN
         MDSIZ  = 1
         MDMIN  = 0
         MDMAX  = 0
         LDSIZ  = 1
         LDMIN  = 0
         LDMAX  = 0
         MDSIZX = 1
      ELSE
         MDSIZ  = NTHMAX
         MDMIN  =-NTHMAX/2+1
         MDMAX  = NTHMAX/2
         LDSIZ  = NTHMAX
         LDMIN  =-NTHMAX/2+1
         LDMAX  = NTHMAX/2
         MDSIZX = 3*NTHMAX/2
      ENDIF
C
      IERR=0
      RETURN
      END
C
C     ****** RADIAL MESH (CYLINDRICAL COORDINATES) ******
C
      subroutine wmmetric_cyl(ierr)
C
      INCLUDE 'wmcomm.inc'

      ierr=0

      drho=(rb/ra)/nrmax
      dth=2.d0*pi/nthmax
      dph=2.d0*pi/nphmax
      dthg=2.d0*pi/nthgm

!     --- radial mesh, q profile ---
      do nr=1,nrmax+1
         rhol=drho*(nr-1)
         xrho(nr)=rhol
         xr(nr)=ra*rhol
         if(rhol.lt.1.d0) then
            qps(nr)=q0+(qa-q0)*rhol**2
         else
            qps(nr)=qa*rhol**2
         endif
      enddo

!     --- 2D grid ---
      do nr=1,nrmax+1
         do nth=1,nthmax
            rps(nth,nr)=ra*xrho(nr)*cos(dth*(nth-1))
            zps(nth,nr)=ra*xrho(nr)*sin(dth*(nth-1))
         enddo
         do nth=1,nthgm
            rpsg(nth,nr)=ra*xrho(nr)*cos(dthg*(nth-1))
            zpsg(nth,nr)=ra*xrho(nr)*sin(dthg*(nth-1))
         enddo
      enddo

!     --- 3D grid, metric, Bsup ---
      do nr=1,nrmax+1
         rhol=xrho(nr)
         do nph=1,nphmax
            do nth=1,nthmax
               rpst(nth,nph,nr)=rps(nth,nr)
               zpst(nth,nph,nr)=zps(nth,nr)

               rg11(nth,nph,nr)= ra**2
               rg12(nth,nph,nr)= 0.d0
               rg13(nth,nph,nr)= 0.d0
               rg22(nth,nph,nr)= (ra*rhol)**2
               rg23(nth,nph,nr)= 0.d0
               rg33(nth,nph,nr)= rr**2
               rj  (nth,nph,nr)= rr*ra**2*rhol

               bfld(2,nth,nph,nr)=0.d0
               bfld(3,nth,nph,nr)=bb/rr
               BPST(NTH,NPH,NR)=bb
            enddo
         enddo
      enddo

!     --- plasma boundary ---
         nsumax=31
         dthu=2.d0*pi/(nsumax-1)
         do nsu=1,nsumax
            rsu(nsu,1)=ra*cos(dthu*(nsu-1))
            zsu(nsu,1)=ra*sin(dthu*(nsu-1))
         enddo

!     --- wall boundary ---
         nswmax=31
         dthw=2.d0*pi/(nswmax-1)
         do nsw=1,nswmax
            rsw(nsu,1)=ra*cos(dthw*(nsw-1))
            zsw(nsu,1)=ra*sin(dthw*(nsw-1))
         enddo

!     --- graphic boundary ---
         rgmin=-rb*1.01d0
         rgmax= rb*1.01d0
         zgmin=-rb*1.01d0
         zgmax= rb*1.01d0

!     --- compatibility ---
         PSIPA=RA*RA*BB/(Q0+QA)
         PSIPB=SQRT(RB**2/RA**2+(RB**2/RA**2-1.D0)*Q0/QA)*PSIPA
C
         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6
C
         DO NR=1,NRMAX+1
            RHOL=XRHO(NR)
            IF(RHOL.LE.1.D0.AND.PN(1).NE.0.D0) THEN
               FEDGE=PNS(1)/PN(1)
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1)**PROFN2+FEDGE
               PT=(PTPR(1)+2*PTPP(1))/3.D0
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1)**PROFT2+FEDGE
               PPS(NR)=P0*FACTN*FACTT
            ELSE
               PPS(NR)=0.D0
            ENDIF
            RBPS(NR)=BB*RR
            VPS(NR)=2*PI*RR*PI*XR(NR)**2
            RLEN(NR)=2*PI*XR(NR)
         ENDDO
      RETURN
      END
C
C     ****** RADIAL MESH (TOROIDAL COORDINATES) ******
C
      SUBROUTINE wmmetric_tor(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      REAL(8) gm(3,3),gj
C
         IERR=0
C
         NSUMAX=31
         NSWMAX=31
         NPHMAX=1
C
         PSIPA=RA*RA*BB/(Q0+QA)
         DRHO=(RB/RA)/NRMAX
C
         DO NR=1,NRMAX+1
            RHOL=DRHO*(NR-1)
            XRHO(NR)=RHOL
            XR(NR)=RA*RHOL
C
C            WRITE(6,*) 'NR,RHO,XR=',NR,XRHO(NR),XR(NR)
C
            call wmfem_qprofile(xrho(nr),qinv)
            QPS(NR)  = 1.D0/qinv
         ENDDO
C
         DTH=2.D0*PI/NTHMAX
         DPH=2.D0*PI/NPHMAX

         DTHG=2.D0*PI/NTHGM
         DO NR=1,NRMAX+1
            DO NTH=1,NTHGM
               RCOS=COS(DTHG*(NTH-1))
               RSIN=SIN(DTHG*(NTH-1))
               RPSG(NTH,NR)  = RR + XR(NR)*RCOS
               ZPSG(NTH,NR)  =      XR(NR)*RSIN
            ENDDO
         ENDDO
C
         DTHU=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            RCOS=COS(DTHU*(NSU-1))
            RSIN=SIN(DTHU*(NSU-1))
            IF(MODELG.EQ.1) THEN
               RSU(NSU,1)=RA*RCOS
               RSW(NSU,1)=RB*RCOS
            ELSE
               RSU(NSU,1)=RR+RA*RCOS
               RSW(NSU,1)=RR+RB*RCOS
            ENDIF
            ZSU(NSU,1)=RA*RSIN
            ZSW(NSU,1)=RB*RSIN
         ENDDO
C
         RGMIN=RR-RB*1.01D0
         RGMAX=RR+RB*1.01D0
         ZGMIN=-RB*1.01D0
         ZGMAX= RB*1.01D0

C
         DO NR=1,NRMAX+1
         DO NTH=1,NTHMAX
            rhol=xrho(nr)
            th=dth*(nth-1)
            call wmfem_metrics(rhol,th,0.d0,gm,gj)
            call wmfem_magnetic(rhol,th,0.d0,babs,bsupth,bsupph)
         DO NPH=1,NPHMAX
            RPST(NTH,NPH,NR)=RR+RA*XRHO(NR)*COS(TH)
            ZPST(NTH,NPH,NR)=   RA*XRHO(NR)*SIN(TN)
C
            RG11(NTH,NPH,NR)= gm(1,1)
            RG12(NTH,NPH,NR)= gm(1,2)
            RG13(NTH,NPH,NR)= gm(1,3)
            RG22(NTH,NPH,NR)= gm(2,2)
            RG23(NTH,NPH,NR)= gm(2,3)
            RG33(NTH,NPH,NR)= gm(3,3)
            RJ  (NTH,NPH,NR)= gj
C
            BFLD(2,NTH,NPH,NR)=bsupth
            BFLD(3,NTH,NPH,NR)=bsupph
            BPST(NTH,NPH,NR)=BABS
         ENDDO
         ENDDO
         ENDDO
      RETURN
      END
C
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE wmmetric_eq(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      CHARACTER*(80) LINE
      SAVE NRMAXSV,NTHMAXSV,NSUMAXSV
      DATA NRMAXSV,NTHMAXSV,NSUMAXSV/0,0,0/
C
      IERR=0
C
      IF(NRMAXSV.EQ.NRMAX.AND.
     &   NTHMAXSV.EQ.NTHMAX.AND.
     &   NSUMAXSV.EQ.NSUMAX.AND.
     &   KNAMEQ_SAVE.EQ.KNAMEQ) RETURN
C
      NSUMAX=41
      NSWMAX=NSUMAX
      NPHMAX=1
C
      IF(NTHMAX.LT.4) THEN
         IF(MYRANK.EQ.0) 
     &        WRITE(6,*) 'XX WMXRZF: NTHMAX MUST BE GREATER THAN 4'
         IERR=1
         RETURN
      ENDIF
C
      IF(MYRANK.EQ.0) THEN
         CALL EQLOAD(MODELG,KNAMEQ,IERR)
      ENDIF
      IF(IERR.EQ.1) RETURN
C
      NRMAXSV=NRMAX
      NTHMAXSV=NTHMAX
      NSUMAXSV=NSUMAX
      KNAMEQ_SAVE=KNAMEQ
C
c$$$      IF(MYRANK.EQ.0) THEN
         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHMAX
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
C
         CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
C
C         WRITE(6,'(1P6E12.4)') BB,RR,RIP,RA,RKAP,RB
C
         CALL EQGETP(RHOT,PSIP,NRMAX+1)
         CALL EQGETR(RPS,DRPSI,DRCHI,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETZ(ZPS,DZPSI,DZCHI,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETBB(BPR,BPZ,BPT,BTP,NTHM,NTHMAX,NRMAX+1)
         CALL EQGETQ(PPS,QPS,RBPS,VPS,RLEN,NRMAX+1)
         CALL EQGETU(RSU,ZSU,RSW,ZSW,NSUMAX)
         CALL EQGETF(RGMIN,RGMAX,ZGMIN,ZGMAX)
         CALL EQGETA(RAXIS,ZAXIS,PSIPA,PSITA,Q0,QA)
C
C         WRITE(6,'(A,1P2E12.4)') 'PSIPA,PSITA=',PSIPA,PSITA
C         WRITE(6,'(1P5E12.4)') (RPS(NTH,NRMAX+1),NTH=1,NTHMAX)
C
         write(LINE,'(A,I5)') 'nrmax=',NRMAX+1
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHMAX
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
         CALL EQGETG(RPSG,ZPSG,NTHGM,NTHGM,NRMAX+1)
C
C         WRITE(6,'(1P5E12.4)') (RPSG(NTH,NRMAX+1),NTH=1,NTHGM)
c$$$      ENDIF
c$$$      CALL MPBCDA(BB)
c$$$      CALL MPBCDA(RR)
c$$$      CALL MPBCDA(RIP)
c$$$      CALL MPBCDA(RA)
c$$$      CALL MPBCDA(RKAP)
c$$$      CALL MPBCDA(RDLT)
c$$$      CALL MPBCDA(RB)
c$$$      CALL MPBCDN(RHOT,NRMAX+1)
c$$$      CALL MPBCDN(PSIP,NRMAX+1)
c$$$      CALL MPBCDN(RPS,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(DRPSI,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(DRCHI,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(ZPS,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(DZPSI,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(DZCHI,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(BPT,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(BTP,NTHM*(NRMAX+1))
c$$$      CALL MPBCDN(PPS,NRMAX+1)
c$$$      CALL MPBCDN(QPS,NRMAX+1)
c$$$      CALL MPBCDN(RBPS,NRMAX+1)
c$$$      CALL MPBCDN(VPS,NRMAX+1)
c$$$      CALL MPBCDN(RLEN,NRMAX+1)
c$$$      CALL MPBCDN(RSU,NSUMAX)
c$$$      CALL MPBCDN(ZSU,NSUMAX)
c$$$      CALL MPBCDN(RSW,NSUMAX)
c$$$      CALL MPBCDN(ZSW,NSUMAX)
c$$$      CALL MPBCDA(RGMIN)
c$$$      CALL MPBCDA(RGMAX)
c$$$      CALL MPBCDA(ZGMIN)
c$$$      CALL MPBCDA(ZGMAX)
c$$$      CALL MPBCDA(RAXIS)
c$$$      CALL MPBCDA(ZAXIS)
c$$$      CALL MPBCDA(PSIPA)
c$$$      CALL MPBCDA(PSITA)
c$$$      CALL MPBCDA(Q0)
c$$$      CALL MPBCDA(QA)
C
         RMAX=RSU(1,1)
         DO NSU=2,NSUMAX
            RMAX=MAX(RMAX,RSU(NSU,1))
         ENDDO
C
         DO NR=1,NRMAX+1
            ARG=RHOT(NR)
            IF(ARG.LE.0.D0) THEN
               XRHO(NR)=0.D0
            ELSE
               XRHO(NR)=ARG
            ENDIF
            XR(NR)=RA*XRHO(NR)
         ENDDO
C
C         WRITE(6,'(I5,1P5E12.4)') (NR,XRHO(NR),RHOT(NR),PSIP(NR),
C     &                     PSITA*RHOT(NR)**2,QPS(NR),NR=1,NRMAX+1)
C
         DO NR=1,NRMAX+1
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            RPST(NTH,NPH,NR)=RPS(NTH,NR)
            ZPST(NTH,NPH,NR)=ZPS(NTH,NR)
         ENDDO
         ENDDO
         ENDDO
C
         DO NR=2,NRMAX+1
         DO NTH=1,NTHMAX
         DO NPH=1,NPHMAX
            RG11(NTH,NPH,NR)= (DRPSI(NTH,NR)**2+DZPSI(NTH,NR)**2)
C     &                       *(2.D0*PI)**2
            RG12(NTH,NPH,NR)= (DRPSI(NTH,NR)*DRCHI(NTH,NR)
     &                        +DZPSI(NTH,NR)*DZCHI(NTH,NR))/XRHO(NR)
C     &                       * 2.D0*PI
            RG13(NTH,NPH,NR)=0.D0
            RG22(NTH,NPH,NR)= (DRCHI(NTH,NR)**2+DZCHI(NTH,NR)**2)
     &                        /(XRHO(NR)*XRHO(NR))
            RG23(NTH,NPH,NR)=0.D0
            RG33(NTH,NPH,NR)= RPS(NTH,NR)**2
            RJ  (NTH,NPH,NR)= RPS(NTH,NR)
     &                      *( DRPSI(NTH,NR)*DZCHI(NTH,NR)
     &                        -DRCHI(NTH,NR)*DZPSI(NTH,NR))/XRHO(NR)
C     &                      / 2.D0*PI
C
C            BFLD(2,NTH,NPH,NR)=1.D0/(2.D0*PI*RJ(NTH,NPH,NR))
C            BFLD(3,NTH,NPH,NR)=RBPS(NR)/RPS(NTH,NR)**2
C     &                        /(2.D0*PI)
C
            BPTL=(BPR(NTH,NR)*DZPSI(NTH,NR)
     &           -BPZ(NTH,NR)*DRPSI(NTH,NR))/XRHO(NR)
     &           /SQRT(RG11(NTH,NPH,NR))
     &           /SQRT(RG22(NTH,NPH,NR))
C     &                       * 2.D0*PI
C
            BFLD(2,NTH,NPH,NR)=BPTL
            BFLD(3,NTH,NPH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
C
C            IF(NTH.EQ.1) WRITE(6,'(2I3,1P6E12.4)') 
C     &           NR,NTH,1.D0/(2.D0*PI*RJ(NTH,NPH,NR)),
C     &           BPTL,
C     &           RBPS(NR)/RPS(NTH,NR)**2/(2.D0*PI),
C     &           BTP(NTH,NR)/RPS(NTH,NR),
C     &           BPZ(NTH,NR),XRHO(NR)
C
C            IF((NR.EQ.2).OR.(NR.EQ.3)) THEN
C            WRITE(6,*) 'NR,NTH,NPH=',NR,NTH,NPH
C            WRITE(6,'(1P3E21.4)') 
C     &           RG11(NTH,NPH,NR),RG12(NTH,NPH,NR),RG13(NTH,NPH,NR)
C            WRITE(6,'(1P3E21.4)') 
C     &           RG22(NTH,NPH,NR),RG23(NTH,NPH,NR),RG33(NTH,NPH,NR)
C            WRITE(6,'(1P3E21.4)') 
C     &           RJ(NTH,NPH,NR),BFLD(2,NTH,NPH,NR),BFLD(3,NTH,NPH,NR)
C            ENDIF
         ENDDO
         ENDDO
         ENDDO
C
         NR=1
         DO NPH=1,NPHMAX
         DO NTH=1,NTHMAX
            RG11(NTH,NPH,NR)= RG11(NTH,NPH,2)
            RG12(NTH,NPH,NR)= RG12(NTH,NPH,2)
            RG13(NTH,NPH,NR)= 0.D0
            RG22(NTH,NPH,NR)= RG22(NTH,NPH,2)
            RG23(NTH,NPH,NR)= 0.D0
            RG33(NTH,NPH,NR)= RPST(NTH,NPH,NR)**2
            RJ  (NTH,NPH,NR)= RJ(NTH,NPH,2)
            BPTL=0.D0
C
C            BFLD(2,NTH,NPH,NR)=1.D0/(2.D0*PI*RJ(NTH,NPH,NR))
C            BFLD(3,NTH,NPH,NR)=RBPS(NR)/RPS(NTH,NR)**2
C     &                        /(2.D0*PI)
C
            BFLD(2,NTH,NPH,NR)=BFLD(2,NTH,NPH,2)
            BFLD(3,NTH,NPH,NR)=BTP(NTH,NR)/RPS(NTH,NR)
C
C            WRITE(6,'(2I3,1P4E12.4)') NTH,NR,1.D0/RJ(NTH,NPH,NR),
C     &           RBPS(NR)/RPS(NTH,NR)**2,
C     &           BPTL,BPTL*RJ(NTH,NPH,NR)
C     &           BPT,BTP(NTH,NR)/RPS(NTH,NR)
C
         ENDDO
         ENDDO
      RETURN
      END
