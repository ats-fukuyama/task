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

      call wmfem_setr(ierr)

      dth=2.d0*pi/nthmax
      dph=2.d0*pi/nphmax
      dthg=2.d0*pi/nthgm

!     --- q profile ---
      do nr=1,nrmax
         if(xrho(nr).lt.1.d0) then
            qps(nr)=q0+(qa-q0)*xrho(nr)**2
         else
            qps(nr)=qa*xrho(nr)**2
         endif
      enddo

!     --- 2D grid ---
      do nr=1,nrmax
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
      do nr=1,nrmax
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
         dthu=2.d0*pi/nsumax
         do nsu=1,nsumax
            rsu(nsu,1)=ra*cos(dthu*(nsu-1))
            zsu(nsu,1)=ra*sin(dthu*(nsu-1))
         enddo

!     --- wall boundary ---
         dthw=2.d0*pi/nswmax
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

         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6

         DO NR=1,NRMAX
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
      REAL(8) gm(3,3),gj

      ierr=0

      call wmfem_setr(ierr)

      dth=2.d0*pi/nthmax
      dph=2.d0*pi/nphmax
      dthg=2.d0*pi/nthgm

      DO NR=1,NRMAX
         call wmfem_qprofile(xrho(nr),qinv)
         QPS(NR)  = 1.D0/qinv
      ENDDO

!     --- 2D grid ---
      DO NR=1,NRMAX
         DO NTH=1,NTHMAX
            RPS(NTH,NR)  = RR + XR(NR)*COS(DTH*(NTH-1))
            ZPS(NTH,NR)  =      XR(NR)*SIN(DTH*(NTH-1))
         ENDDO
         DO NTH=1,NTHGM
            RPSG(NTH,NR)  = RR + XR(NR)*COS(DTHG*(NTH-1))
            ZPSG(NTH,NR)  =      XR(NR)*SIN(DTHG*(NTH-1))
         ENDDO
      ENDDO

!     --- 3D grid, metric, Bsup ---
      DO NR=1,NRMAX
         DO NTH=1,NTHMAX
            rhol=xrho(nr)
            th=dth*(nth-1)
            call wmfem_metrics(rhol,th,0.d0,gm,gj)
            call wmfem_magnetic(rhol,th,0.d0,babs,bsupth,bsupph)
            DO NPH=1,NPHMAX
               RPST(NTH,NPH,NR)=RPS(NTH,NR)
               ZPST(NTH,NPH,NR)=ZPS(NTH,NR)

               RG11(NTH,NPH,NR)= gm(1,1)
               RG12(NTH,NPH,NR)= gm(1,2)
               RG13(NTH,NPH,NR)= gm(1,3)
               RG22(NTH,NPH,NR)= gm(2,2)
               RG23(NTH,NPH,NR)= gm(2,3)
               RG33(NTH,NPH,NR)= gm(3,3)
               RJ  (NTH,NPH,NR)= gj

               BFLD(2,NTH,NPH,NR)=bsupth
               BFLD(3,NTH,NPH,NR)=bsupph
               BPST(NTH,NPH,NR)=BABS
            ENDDO
         ENDDO
      ENDDO

!     --- plasma boundary ---
      DTHU=2.D0*PI/NSUMAX
      DO NSU=1,NSUMAX
         RSU(NSU,1)=RR+RA*COS(DTHU*(NSU-1))
         ZSU(NSU,1)=   RA*SIN(DTHU*(NSU-1))
      ENDDO

!     --- wall boundary ---
      DTHW=2.D0*PI/NSWMAX
      DO NSW=1,NSWMAX
         RSW(NSW,1)=RR+RA*COS(DTHW*(NSW-1))
         ZSW(NSW,1)=   RA*SIN(DTHW*(NSW-1))
      ENDDO

!     --- graphic boundary ---
      RGMIN=RR-RB*1.01D0
      RGMAX=RR+RB*1.01D0
      ZGMIN=-RB*1.01D0
      ZGMAX= RB*1.01D0

!     --- compatibility ---
         PSIPA=RA*RA*BB/(Q0+QA)
         PSIPB=SQRT(RB**2/RA**2+(RB**2/RA**2-1.D0)*Q0/QA)*PSIPA

         P0=0.D0
         DO NS=1,NSMAX
            P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
         ENDDO
         P0=P0*1.D20*AEE*1.D3/1.D6

         DO NR=1,NRMAX
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
C     ****** RADIAL MESH AND METRIC TENSOR ******
C
      SUBROUTINE wmmetric_eq(IERR)
C
      INCLUDE 'wmcomm.inc'
      real(8),dimension(3,3):: gm
      real(8):: gj
      CHARACTER*(80) LINE
C
      IERR=0
C
      IF(NTHMAX.LT.4) THEN
         IF(MYRANK.EQ.0) 
     &        WRITE(6,*) 'XX WMXRZF: NTHMAX MUST BE GREATER THAN 4'
         IERR=1
         RETURN
      ENDIF
C
      IF(KNAMEQ_SAVE.NE.KNAMEQ) THEN
C
         CALL EQLOAD(MODELG,KNAMEQ,IERR)
         IF(IERR.NE.0) RETURN
C
         KNAMEQ_SAVE=KNAMEQ
C
         write(LINE,'(A,I5)') 'nrmax=',NRMAX
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nthmax=',NTHMAX
         call eqparm(2,line,ierr)
         write(LINE,'(A,I5)') 'nsumax=',NSUMAX
         call eqparm(2,line,ierr)
         CALL EQCALQ(IERR)
         
         CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
         CALL EQGETQ(PPS,QPS,RBPS,VPS,RLEN,NRMAX)
         CALL EQGETA(RAXIS,ZAXIS,PSIPA,PSITA,Q0,QA)
         CALL EQGETU(RSU,ZSU,RSW,ZSW,NSUMAX)
         NSWMAX=NSUMAX
         CALL EQGETF(RGMIN,RGMAX,ZGMIN,ZGMAX)
      ENDIF

!     --- radial mesh after RA,RD,RB determined ---

      call wmfem_setr(ierr)

!     --- 3D grid, metric, Bsup ---
      dth=2.d0*pi/nthmax
      DO NR=1,NRMAX
         DO NTH=1,NTHMAX
            rhol=xrho(nr)
            th=dth*(nth-1)
            call wmeq_get_posrz(rhol,th,rrl,zzl,
     &                          drrrho,dzzrho,drrchi,dzzchi)
            call wmfem_metrics(rhol,th,0.d0,gm,gj)
            call wmfem_magnetic(rhol,th,0.d0,babs,bsupth,bsupph)
            DO NPH=1,NPHMAX
               RPST(NTH,NPH,NR)= rrl
               ZPST(NTH,NPH,NR)= zzl

               RG11(NTH,NPH,NR)= gm(1,1)
               RG12(NTH,NPH,NR)= gm(1,2)
               RG13(NTH,NPH,NR)= gm(1,3)
               RG22(NTH,NPH,NR)= gm(2,2)
               RG23(NTH,NPH,NR)= gm(2,3)
               RG33(NTH,NPH,NR)= gm(3,3)
               RJ  (NTH,NPH,NR)= gj

               BFLD(2,NTH,NPH,NR)=bsupth
               BFLD(3,NTH,NPH,NR)=bsupph
               BPST(NTH,NPH,NR)=  babs
!               write(6,'(A,2I5,1P5E12.4)') '%%%',NR,NTH,
!     &              babs,bsupth,bsupph,bsupth*sqrt(gm(2,2)),
!     &              bsupph*sqrt(gm(3,3))
            ENDDO
         ENDDO
      ENDDO

      dth=2.d0*pi/nthgm
      DO NR=1,NRMAX
         DO NTH=1,NTHGM
            rhol=xrho(nr)
            th=dth*(nth-1)
            call wmeq_get_posrz(rhol,th,rrl,zzl,
     &                          drrrho,dzzrho,drrchi,dzzchi)
            RPSG(NTH,NR)= rrl
            ZPSG(NTH,NR)= zzl
         end do
      end do

      RETURN
      END

!     ****** setup RADIAL MESH ******

      subroutine wmfem_setr(ierr)

      INCLUDE 'wmcomm.inc'
      integer:: nrmaxa,nrmaxd,nrmaxb
      real(8):: drho

      IF(rd.LE.ra.OR.rd.GE.rb) THEN
         rd=0.5D0*(ra+rb)
         WRITE(6,'(A,1P3E12.4)') 
     &        '!! wmfem_setr: rd = (ra+rb)/2: ra,rd,rb=',ra,rd,rb
      ENDIF

      drho=rb/(nrmax-1)
      nrmaxb=max(nint((rb-rd)/drho),2)
      nrmaxd=max(nint((rd-ra)/drho),2)
      nrmaxa=nrmax-nrmaxb-nrmaxd
      drho=1.d0/(nrmaxa-1)
      do nr=1,nrmaxa-1
         xrho(nr)=drho*(nr-1)
      enddo
      xrho(nrmaxa)=1.d0
      drho=(rd-ra)/(ra*nrmaxd)
      nrmaxd=nrmaxa+nrmaxd
      do nr=nrmaxa+1,nrmaxd-1
         xrho(nr)=1.d0+drho*(nr-nrmaxa)
      enddo
      xrho(nrmaxd)=rd/ra
      drho=(rb-rd)/(ra*nrmaxb)
      do nr=nrmaxd+1,nrmax-1
         xrho(nr)=rd/ra+drho*(nr-nrmaxd)
      enddo
      xrho(nrmax)=rb/ra

      do nr=1,nrmax
         xr(nr)=ra*xrho(nr)
      enddo
      return
      end
