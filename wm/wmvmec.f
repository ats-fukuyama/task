C     $Id$
C
C    ****** RADIAL MESH AND METRIC TENSOR ******
C
C
      SUBROUTINE WMXRZV(IERR)
C
      INCLUDE 'wmcomm.inc'
C     
      IERR=0
C
C
      IF(KNAMEQ.NE.KNAMEQ_SAVE) THEN
         CALL WMHGRD(IERR)
         IF(IERR.NE.0) RETURN
         KNAMEQ_SAVE=KNAMEQ
      ENDIF
      CALL WMHCRZ
      CALL WMHMTR
      CALL WMHCBB
C
      RETURN
      END
C
C    ****** READ WOUT FILE ******
C
      SUBROUTINE WMHGRD(IERR)
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      IERR=0
C
      NFL=8
      CALL FROPEN(NFL,KNAMEQ,1,0,'EQ',IERR)
C
C      open(8,FILE=FNAME,FORM='FORMATTED',STATUS='OLD')
c
c     write (nfort8,711) voli,gamma,1.0/nfp,Rmajor_p,bz_vmec,
c    >       ncurr,mnmax,ns,ntheta,nzeta,itfsq,niter/nstep+1,
c    >       mpol,ntor,ntor+1,1
      read  (     8,711) voli, rgam,     d2,      rc,     d3,
     >          n1,mnmax,nsrmax,ntheta_in,nzeta_in,n2,n3,
     >       mpol_in,ntor_in,n4,n5
 711  format(5e20.13,11i6)
c
            ntheta1 = 2*(ntheta/2)
            ntheta2 = 1 + ntheta1/2
            nznt = nzeta*ntheta2
            mpol1 = mpol_in - 1
            ntor0 = ntor_in + 1 ! or ntor_in not definite
            mn0   = 1           ! not definite
c
       print *,   ' nsrmax    = ',nsrmax
       print *,   '  mnmax    = ', mnmax
       print *,   '   mpol_in = ',  mpol_in
       print *,   '   ntor_in = ',  ntor_in
       print *,   ' ntheta    = ',ntheta
       print *,   '  nzeta    = ', nzeta
       print *,   '   nznt    = ',  nznt
c
      if( nsrmax .GT. nsd ) then
        print 600,' ns_in     = ',   nsrmax,'  > nsd    = ',nsd
        stop
      endif
      if( mnmax .GT. nmnm ) then
        print 600,' mnmax     = ',    mnmax,'  > nmnm   = ',nmnm
        stop
      endif
      if( mpol_in .NE. mpol ) then
        print 600,' mpol_in   = ',  mpol_in,' /= mpol   = ',mpol
        stop
      endif
      if( ntor_in .NE. nmax ) then
        print 600,' ntor_in   = ',  ntor_in,' /= nmax   = ',nmax
        stop
      endif
      if( ntheta_in .NE. ntheta ) then
        print 600,' ntheta_in = ',ntheta_in,' /= ntheta = ',ntheta
        stop
      endif
      if( nzeta_in .NE. nzeta ) then
        print 600,' nzeta_in  = ', nzeta_in,' /= nzeta  = ',nzeta
        stop
      endif
 600  format(2x,a13,i6,a13,i6)
c
            mn = 0
      do js = 1, nsrmax
         do m = 0, mpol1
            nmin0 = -ntor_in
            if (m .eq. 0) nmin0 = 0
            do n = nmin0, ntor_in
               mn = mn + 1
c                write (nfort8,722) xm(mn),xn(mn),
c    f			rmnc(mn),zmns(mn),lmns(mn),
c    h			bmn,gmn,
c    h			bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
c	 ----------------------------
                 read  (     8,722) xm(mn),xn(mn),
     f			rmnc(mn),zmns(mn),rlmns(mn),
     h			bmod(mn),rgmod(mn),
     h			     d1,      d2,      d3, bsu(mn), bsv(mn)
c	 ----------------------------
            enddo
         enddo
c          write (nfort8,723)
c    h          ( gsqrt(js,k),
c    h            bsubu(js,k),  bsubv(js,k),  bsubs(js,k),
c    h            bsupu(js,k),  bsupv(js,k),
c    h            bmod    (k),  k=1,nznt                  )
           read (      8,723) (d1,d2,d3,d4,d5,d6,d7,k=1,nznt)
      enddo
 722  format(12e20.13)
 723  format(7e20.13)
c
c          write (nfort8,732) (iotas(js),pres(js),vp(js),phips(js),
c    1            buco(js),bvco(js),phi(js),chi(js),jcuru(js),
c    2            jcurv(js),specw(js),js=2,ns)
c	 --------------------------------------------------------
           read  (     8,732) (riotas(js),pres(js),vp(js),phips(js),
     1            bpco(js),baco(js),phi(js),rchi(js),rjtheta(js),
     2            rjzeta(js),specw(js),js=2,nsrmax)
 732  format(11e20.13)
      print 995, (js,riotas(js),pres(js),vp(js),phips(js),
     >               bpco(js),baco(js),phi(js),rchi(js),rjtheta(js),
     >               rjzeta(js),specw(js), js=2,nsrmax)
 995  format(/1x,' is' ,1x,'    iota_h',1x,'    pres_h'
     >                 ,1x,'      vp_h',1x,'    phip_h'
     >                 ,1x,'    buco_h',1x,'    bvco_h'
     >                 ,1x,'     phi_f',1x,'     chi_f'
     >                 ,1x,'   jcuru_f',1x,'   jcurv_f'
     >                 ,1x,'     specw'/(1x,i3,11(1x,1pd10.3)) )

c       write (nfort8,734) (am(i),i=0,10)
c       write (nfort8,735) (fsqt(i),wdot(i),i=1,100)
        read  (     8,734) (amvm(i),i=0,10)
        read  (     8,735) (fsqt(i),wdot(i),i=1,100)
 734  format(11e20.13)
 735  format(2e20.13)
c
      close(8)
C
C     ----- Setup normalized radius for file data -----
C
         phi(1) = 0.0d0
      PSIA =-PHI(NSRMAX)
      DO NSR=1,NSRMAX
         S=-PHI(NSR)/PSIA
         IF(S.LT.0.D0) THEN
            XS(NSR)=0.D0
         ELSE
            XS(NSR)=SQRT(S)
         ENDIF
         IF(NSR.EQ.1) THEN
            XSH(NSR)=0.D0
         ELSE
            XSH(NSR)=0.5D0*(XS(NSR)+XS(NSR-1))
         ENDIF
C         WRITE(6,'(I5,2X,1P4E12.5)') 
C     &           NSR,PHI(NSR),XS(NSR),PRES(NSR),riotas(nsr)
      ENDDO
C
      RHOB=RB/RA
      DRHO=RHOB/NRMAX
      DO NR=1,NRMAX+1
         XRHO(NR)=DRHO*(NR-1)
         XR(NR)  =RB*XRHO(NR)
      ENDDO
C
      RETURN
      END
C
C    ****** CALCULATE R AND Z ******
C
      SUBROUTINE WMHCRZ
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C      ***** SPLINE PSIPS *****
C
      CALL SPL1D(XS,PHI,FX6,U6,NSRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMHCRZ: SPL1D: PHI: IEER=',IERR
C
      DO NR=1,NRMAX+1
         IF(XRHO(NR).GT.1.D0) THEN
            PSIPS(NR)=0.D0
         ELSE
            CALL SPL1DF(XRHO(NR),PSIPS(NR),XS,U6,NSRMAX,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX WMHCRZ: SPL1DF: PSIPS: IEER=',IERR
            ENDIF         
         ENDIF         
      ENDDO
C     
C      ***** SPLINE RMNC(S),ZMNS(S) DRMNC(S) DZMNS(S) *****
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            YRBS(NSR)=RMNC(MN+(NSR-1)*MNMAX)
            YZBS(NSR)=ZMNS(MN+(NSR-1)*MNMAX)
         ENDDO
C         IF(MYRANK.EQ.0) THEN
C            WRITE(6,'( I5,1P6E12.4)') MN,XM(MN),XN(MN),
C     &                                YRBS(1),YRBS(2),YRBS(3),YRBS(4)
C            WRITE(6,'(29X,1P4E12.4)') YZBS(1),YZBS(2),YZBS(3),YZBS(4)
C         ENDIF
         CALL SPL1D(XS,YRBS,FX1,U1(1,1,MN),NSRMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX WMHCRZ: SPL1D: YRBS: IEER=',IERR
         CALL SPL1D(XS,YZBS,FX2,U2(1,1,MN),NSRMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX WMHCRZ: SPL1D: YZBS: IEER=',IERR
C
         DO NR=1,NRMAX+1
            IF(XRHO(NR).GT.1.D0) THEN
               DRMNC(MN,NR)=(YRBS(NSRMAX)-YRBS(NSRMAX-1))
     &                     /(XS(NSRMAX)-XS(NSRMAX-1))
               SRMNC(MN,NR)=YRBS(NSRMAX)
     &                     +DRMNC(MN,NR)*(XRHO(NR)-XS(NSRMAX))
            ELSE
               CALL SPL1DD(XRHO(NR),SRMNC(MN,NR),DRMNC(MN,NR),
     &                     XS,U1(1,1,MN),NSRMAX,IERR)
               IF(IERR.NE.0) THEN
                  WRITE(6,*) 'XX WMHCRZ: SPL1DD: SRMNC: IEER=',IERR
C                 WRITE(6,'(3I5,1P2E12.4)') MN,NR,IERR,XRHO(NR),XS(NSRMAX)
               ENDIF
            ENDIF
C
            IF(XRHO(NR).GT.1.D0) THEN
               DZMNS(MN,NR)=(YZBS(NSRMAX)-YZBS(NSRMAX-1))
     &                     /(XS(NSRMAX)-XS(NSRMAX-1))
               SZMNS(MN,NR)=YZBS(NSRMAX)
     &                     +DZMNS(MN,NR)*(XRHO(NR)-XS(NSRMAX))
            ELSE
               CALL SPL1DD(XRHO(NR),SZMNS(MN,NR),DZMNS(MN,NR),
     &                     XS,U2(1,1,MN),NSRMAX,IERR)
               IF(IERR.NE.0) THEN
                  WRITE(6,*) 'XX WMHCRZ: SPL1DD: SSMNC: IEER=',IERR
C                  WRITE(6,'(3I5,1P2E12.4)') MN,NR,IERR,XRHO(NR),XS(NSRMAX)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C     ***** CALCULATE R,Z *****
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            RMNCC(MN,NSR)=RMNC(MN+(NSR-1)*MNMAX)
            ZMNSS(MN,NSR)=ZMNS(MN+(NSR-1)*MNMAX)
         ENDDO
      ENDDO
      RETURN
      END
C
C    ****** CALCULTE METRIC ******
C
      SUBROUTINE WMHMTR
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C      ***** CULCULATE METRIC TENSOR AND JACOBIAN*****
C 
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/(NHC*NPHMAX)
C
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C
         RPSS=0.D0
         ZPSS=0.D0
         DRS=0.D0
         DZS=0.D0
         DRTH=0.D0
         DZTH=0.D0
         DRPH=0.D0
         DZPH=0.D0
         DPHPH=0.D0
C
         TH=DTH*(NTH-1)
         PH=DPH*(NPH-1)
         DO MN=1,MNMAX
            RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
            RCOS=COS(XM(MN)*TH-XN(MN)*PH)
            RPSS =RPSS +       SRMNC(MN,NR)*RCOS
            ZPSS =ZPSS +       SZMNS(MN,NR)*RSIN
            DRS  =DRS  +       DRMNC(MN,NR)*RCOS
            DZS  =DZS  +       DZMNS(MN,NR)*RSIN
            DRTH =DRTH -XM(MN)*SRMNC(MN,NR)*RSIN
            DZTH =DZTH +XM(MN)*SZMNS(MN,NR)*RCOS
            DRPH =DRPH +XN(MN)*SRMNC(MN,NR)*RSIN
            DZPH =DZPH -XN(MN)*SZMNS(MN,NR)*RCOS
            DPHPH=DPHPH+       SRMNC(MN,NR)*RCOS
         ENDDO
C
C  *****  DRS,DZS,DPHS,DRTH,DZTH,DPHTH,DRPH,DZPH,DPHPH GRAPH *****
C            
         RPST(  NTH,NPH,NR)=RPSS
         ZPST(  NTH,NPH,NR)=ZPSS
C
         IF(NR.NE.1) THEN
            DRS=DRS/(2.D0*PSIA*XRHO(NR))
            DZS=DZS/(2.D0*PSIA*XRHO(NR))
            RG11(NTH,NPH,NR)=(DRS *DRS +DZS *DZS )*XRHO(NR)**2
            RG12(NTH,NPH,NR)=(DRS *DRTH+DZS *DZTH)
            RG13(NTH,NPH,NR)=(DRS *DRPH+DZS *DZPH)*XRHO(NR)
            RG22(NTH,NPH,NR)=(DRTH*DRTH+DZTH*DZTH)/XRHO(NR)**2
            RG23(NTH,NPH,NR)=(DRTH*DRPH+DZTH*DZPH)/XRHO(NR)
            RG33(NTH,NPH,NR)= DRPH*DRPH+DZPH*DZPH+DPHPH*DPHPH
            RJ  (NTH,NPH,NR)=DPHPH*(DRS*DZTH-DZS*DRTH)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         RG11(NTH,NPH,1)=RG11(NTH,NPH,2)
         RG12(NTH,NPH,1)=RG12(NTH,NPH,2)
         RG13(NTH,NPH,1)=RG13(NTH,NPH,2)
         RG22(NTH,NPH,1)=RG22(NTH,NPH,2)
         RG23(NTH,NPH,1)=RG23(NTH,NPH,2)
         RG33(NTH,NPH,1)=RG33(NTH,NPH,2)
         RJ(NTH,NPH,1)  =RJ(NTH,NPH,2)
      ENDDO
      ENDDO
C
      RETURN
      END
C
C    ***** SPLINE POLOIDAL AND TOROIDAL MAGNETIC FIELD *****
C
      SUBROUTINE WMHCBB
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
      DIMENSION BSUS(NSRM),BSVS(NSRM)
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            BSUS(NSR)=BSU(MN+(NSR-1)*MNMAX)
            BSVS(NSR)=BSV(MN+(NSR-1)*MNMAX)
         ENDDO
C         IF(MYRANK.EQ.0) THEN
C            WRITE(6,'( I5,1PE10.2,1P5E12.4)') MN,XM(MN),
C     &                      BSUS(1),BSUS(2),BSUS(3),BSUS(4),BSUS(5)
C            WRITE(6,'( 5X,1PE10.2,1P5E12.4)') XN(MN),
C     &                      BSVS(1),BSVS(2),BSVS(3),BSVS(4),BSVS(5)
C         ENDIF
         IF(XM(MN).EQ.0) THEN
            BSUS(1)=(4*BSUS(2)-BSUS(3))/3.D0
            BSVS(1)=(4*BSVS(2)-BSVS(3))/3.D0
C            BSUS(1)=(2*BSUS(2)-BSUS(3))/1.D0
C            BSVS(1)=(2*BSVS(2)-BSVS(3))/1.D0
            FX3(1)=0.D0
            FX4(1)=0.D0
            CALL SPL1D(XS,BSUS,FX3,U3(1,1,MN),NSRMAX,1,IERR)
            CALL SPL1D(XS,BSVS,FX4,U4(1,1,MN),NSRMAX,1,IERR)
         ELSEIF(XM(MN).EQ.1) THEN
C            BSUS(1)=(SQRT(2.D0)*BSUS(2)-BSUS(3))/(SQRT(2.D0)-1.D0)
C            BSVS(1)=(SQRT(2.D0)*BSVS(2)-BSVS(3))/(SQRT(2.D0)-1.D0)
            BSUS(1)=(2*BSUS(2)-BSUS(3))/1.D0
            BSVS(1)=(2*BSVS(2)-BSVS(3))/1.D0
C            BSUS(1)= 0.D0
C            BSVS(1)= 0.D0
            CALL SPL1D(XS,BSUS,FX3,U3(1,1,MN),NSRMAX,0,IERR)
            CALL SPL1D(XS,BSVS,FX4,U4(1,1,MN),NSRMAX,0,IERR)
         ELSE
            BSUS(1)= 0.D0
            BSVS(1)= 0.D0
C            FX1(1)=0.D0
C            FX2(1)=0.D0
            CALL SPL1D(XS,BSUS,FX3,U3(1,1,MN),NSRMAX,0,IERR)
            CALL SPL1D(XS,BSVS,FX4,U4(1,1,MN),NSRMAX,0,IERR)
         ENDIF
C
         BSTHSV(MN)=BSUS(NSRMAX)
         BSTHSD(MN)=(BSUS(NSRMAX)-BSUS(NSRMAX-1))
     &             /(XS(NSRMAX)  -XS(NSRMAX-1))
         BSPHSV(MN)=BSVS(NSRMAX)
         BSPHSD(MN)=(BSVS(NSRMAX)-BSVS(NSRMAX-1))
     &             /(XS(NSRMAX)  -XS(NSRMAX-1))
         DO NR=1,NRMAX+1
            IF(XRHO(NR)-XS(NSRMAX).GT.0.D0) THEN
               BSTHL=BSTHSV(MN)+BSTHSD(MN)*(XRHO(NR)-XS(NSRMAX))
               BSPHL=BSPHSV(MN)+BSPHSD(MN)*(XRHO(NR)-XS(NSRMAX))
            ELSE
               CALL SPL1DF(XRHO(NR),BSTHL,XS,U1(1,1,MN),NSRMAX,IERR)
               CALL SPL1DF(XRHO(NR),BSPHL,XS,U2(1,1,MN),NSRMAX,IERR)
            ENDIF
            BSTH(MN,NR)=BSTHL
            BSPH(MN,NR)=BSPHL
         ENDDO
      ENDDO
C
C      ***** CULCULATE MAGNETIC FIELD *****
C 
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/(NHC*NPHMAX)
C      DO NR=2,NRMAX+1
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
         TH=DTH*(NTH-1)
         PH=DPH*(NPH-1)
C
         SBTH=0.D0
         SBPH=0.D0
         DO MN=1,MNMAX
            RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
            RCOS=COS(XM(MN)*TH-XN(MN)*PH)
            SBTH=SBTH+BSTH(MN,NR)*RCOS
            SBPH=SBPH+BSPH(MN,NR)*RCOS
         ENDDO
         BFLD(2,NTH,NPH,NR)=SBTH
         BFLD(3,NTH,NPH,NR)=SBPH
      ENDDO
      ENDDO
      ENDDO
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C         BFLD(2,NTH,NPH,1)=(4*BFLD(2,NTH,NPH,2)-BFLD(2,NTH,NPH,3))/3
C         BFLD(3,NTH,NPH,1)=(4*BFLD(3,NTH,NPH,2)-BFLD(3,NTH,NPH,3))/3
      ENDDO
      ENDDO
C     
C    ***************************************
C
      P0=0.D0
      DO NS=1,NSMAX
         P0=P0+PN(NS)*(PTPR(NS)+2*PTPP(NS))/3.D0
      ENDDO
      P0=P0*1.D20*AEE*1.D3/1.D6
C
      DO NR=1,NRMAX+1
         RHOL=XRHO(NR)
         IF(RHOL.LE.1.D0) THEN
            IF(PN(1).LE.0.D0) THEN
               FACTN=0.D0
            ELSE
               FEDGE=PNS(1)/PN(1)
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1)**PROFN2+FEDGE
            ENDIF
            PT=(PTPR(1)+2*PTPP(1))/3.D0
            IF(PT.LE.0.D0) THEN
               FACTT=0.D0
            ELSE
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1)**PROFT2+FEDGE
            ENDIF
            PPS(NR)=P0*FACTN*FACTT
         ELSE
            PPS(NR)=0.D0
         ENDIF
      ENDDO
C
      NSUMAX=31
      DTHU=2.D0*PI/(NSUMAX-1)
      DO NSU=1,NSUMAX
         DO NPH=1,NPHMAX
            RSU(NSU,NPH)=0.D0
            ZSU(NSU,NPH)=0.D0
            TH=DTHU*(NSU-1)
            PH=DPH*(NPH-1)
            DO MN=1,MNMAX
               RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
               RCOS=COS(XM(MN)*TH-XN(MN)*PH)
               RSU(NSU,NPH)=RSU(NSU,NPH)+RMNCC(MN,NSRMAX)*RCOS
               ZSU(NSU,NPH)=ZSU(NSU,NPH)+ZMNSS(MN,NSRMAX)*RSIN
            ENDDO
         ENDDO
      ENDDO
C
      NSWMAX=31
      DTHW=2.D0*PI/(NSWMAX-1)
      DO NSW=1,NSWMAX
         DO NPH=1,NPHMAX
            RSW(NSW,NPH)=0.D0
            ZSW(NSW,NPH)=0.D0
            TH=DTHW*(NSW-1)
            PH=DPH*(NPH-1)
            DO MN=1,MNMAX
               RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
               RCOS=COS(XM(MN)*TH-XN(MN)*PH)
               RSW(NSW,NPH)=RSW(NSW,NPH)+SRMNC(MN,NRMAX+1)*RCOS
               ZSW(NSW,NPH)=ZSW(NSW,NPH)+SZMNS(MN,NRMAX+1)*RSIN
            ENDDO
         ENDDO
      ENDDO
C
C     ***** SPLINE IOTAS AND CULCULATE QPS *****
C
      RIOTAS(1)=2.D0*RIOTAS(2)-RIOTAS(3)
C
      CALL SPL1D(XS,RIOTAS,FX5,U5,NSRMAX,0,IERR)
C
      DO NR=1,NRMAX+1
         CALL SPL1DF(XRHO(NR),RIOTASL,XS,U5,NSRMAX,IERR)
         IF(IERR.EQ.2) THEN
            RIOTASL=RIOTAS(NSRMAX)+(RIOTAS(NSRMAX)-RIOTAS(NSRMAX-1))
     &                      /(XS(NSRMAX)-XS(NSRMAX-1))
     &                      *(XRHO(NR)-XS(NSRMAX))
         ENDIF
         QPS(NR)=2.D0*PI/RIOTASL
      ENDDO                     
C
C     ***** COMPUTE R,Z MAGNETIC AXES *****
C
      RGMIN=RSW(1,1)
      RGMAX=RSW(1,1)
      ZGMIN=ZSW(1,1)
      ZGMAX=ZSW(1,1)
      DO NPH=1,NPHMAX
         DO NSW=1,NSWMAX
            RGMIN=MIN(RGMIN,RSW(NSW,NPH))
            RGMAX=MAX(RGMAX,RSW(NSW,NPH))
            ZGMIN=MIN(ZGMIN,ZSW(NSW,NPH))
            ZGMAX=MAX(ZGMAX,ZSW(NSW,NPH))
         ENDDO
      ENDDO
      RR=0.5D0*(RGMIN+RGMAX)
C
      RETURN
      END
C
C    ***** LOCAL POLOIDAL AND TOROIDAL MAGNETIC FIELD *****
C
      SUBROUTINE WMHCBL(XRHOL,THL,PHL,BFLD2L,BFLD3L,QPSL)
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      SBTH=0.D0
      SBPH=0.D0
      DO MN=1,MNMAX
         IF(XRHOL-XS(NSRMAX).GT.0.D0) THEN
            BSTHL=BSTHSV(MN)+BSTHSD(MN)*(XRHOL-XS(NSRMAX))
            BSPHL=BSPHSV(MN)+BSPHSD(MN)*(XRHOL-XS(NSRMAX))
         ELSE
            CALL SPL1DF(XRHOL,BSTHL,XS,U1(1,1,MN),NSRMAX,IERR)
            CALL SPL1DF(XRHOL,BSPHL,XS,U2(1,1,MN),NSRMAX,IERR)
         ENDIF
C
C      ***** CULCULATE MAGNETIC FIELD *****
C 
         RSIN=SIN(XM(MN)*THL-XN(MN)*PHL)
         RCOS=COS(XM(MN)*THL-XN(MN)*PHL)
         SBTH=SBTH+BSTHL*RCOS
         SBPH=SBPH+BSPHL*RCOS
      ENDDO
      BFLD2L=SBTH
      BFLD3L=SBPH
C
      IF(XRHOL-XS(NSRMAX).GT.0.D0) THEN
         RIOTASL=RIOTAS(NSRMAX)+(RIOTAS(NSRMAX)-RIOTAS(NSRMAX-1))
     &                         /(XS(NSRMAX)-XS(NSRMAX-1))
     &                         *(XRHOL-XS(NSRMAX))
      ELSE
         CALL SPL1DF(XRHOL,RIOTASL,XS,U5,NSRMAX,IERR)
      ENDIF
      QPSL=2.D0*PI/RIOTASL
C
      RETURN
      END
C
C   ***** DRAW R,Z GRAPH *****
C
      SUBROUTINE WMGS1D(XG,YG,NXGMAX,KTITL,NGD)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION XG(NXGMAX),YG(NXGMAX)
      DIMENSION GX(NRM),GY(NRM)
      CHARACTER KTITL*6
C
      DIMENSION GP(4,4)
      DATA GP/ 3.0, 10.8,  9.5, 16.5,
     &         3.0, 10.8,  1.0,  8.0,
     &        13.8, 21.6,  9.5, 16.5,
     &        13.8, 21.6,  1.0,  8.0/
C
      DO NXG=1,NXGMAX
         GX(NXG)=GUCLIP(XG(NXG))
         GY(NXG)=GUCLIP(YG(NXG))
      ENDDO
      GXMIN=0.0
      GXMAX=GUCLIP(RB/RA)
C
      CALL WMGSUB(NXGMAX,GX,GXMIN,GXMAX,GY,1,GP(1,NGD),KTITL)
C
      RETURN
      END
