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
     f          rmnc(mn),zmns(mn),rlmns(mn),
     h          bmod(mn),rgmod(mn),
     h          d1,      d2,      d3, bsu(mn), bsv(mn)
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
c     PSIPA =-PHI(NSRMAX)
      PSIPA = PHI(NSRMAX)
      DO NSR=1,NSRMAX
c        S=-PHI(NSR)/PSIPA
         S= PHI(NSR)/PSIPA
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
          WRITE(6,'(I5,2X,1P4E12.5)') 
     &           NSR,PHI(NSR),XS(NSR),PRES(NSR),riotas(nsr)
      ENDDO
C
      RHOB=RB/RA
      DRHO=RHOB/NRMAX
      write(6,*) 'NR,XRHO,XR'
      DO NR=1,NRMAX+1
         XRHO(NR)=DRHO*(NR-1)
         XR(NR)  =RB*XRHO(NR)
      write(6,*) nr,xrho(nr),xr(nr)
      ENDDO
      write(6,*) 'RA,RB,DRHO=',RA,RB,DRHO
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
      integer np,i
      parameter (np=3)
      real*8 nra,psipax,srmnca,drmnca,szmnsa,dzmnsa
      dimension nra(np),psipax(np),srmnca(np),drmnca(np),szmnsa(np)
      dimension dzmnsa(np)
C
C      ***** SPLINE PSIP *****
C
      CALL SPL1D(XS,PHI,FX6,U6,NSRMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMHCRZ: SPL1D: PHI: IEER=',IERR
C
c      DPSIP=(PHI(NSRMAX)-PHI(NSRMAX-5))/(XS(NSRMAX)-XS(NSRMAX-5))
      DO NR=1,NRMAX+1
         IF(XRHO(NR).GT.1.D0) THEN
c            PSIP(NR)=PHI(NSRMAX)+DPSIP*(XRHO(NR)-XS(NSRMAX))
c
c
            do i=1,np
               nra(i)=nr-np-1+i
               psipax(i)=psip(nr-np-1+i)
            enddo
c     
            call polint(nra,psipax,np,nr,psip(nr),dy) 
c
         ELSE
            CALL SPL1DF(XRHO(NR),PSIP(NR),XS,U6,NSRMAX,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX WMHCRZ: SPL1DF: PSIP: IEER=',IERR
               WRITE(6,'(A,I5,1P3E12.4)') 
     &                    '   NR,XRHO(NR),XS(1),XS(NRMAX)=',
     &                        NR,XRHO(NR),XS(1),XS(NRMAX)
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
C         IF(NRANK.EQ.0) THEN
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
c               DRMNC(MN,NR)=(YRBS(NSRMAX)-YRBS(NSRMAX-1))
c     &                     /(XS(NSRMAX)-XS(NSRMAX-1))
c               SRMNC(MN,NR)=YRBS(NSRMAX)
c     &                     +DRMNC(MN,NR)*(XRHO(NR)-XS(NSRMAX))
c
c
            do i=1,np
               nra(i)=nr-np-1+i
               drmnca(i)=drmnc(mn,nr-np-1+i)
            enddo
c     
            call polint(nra,drmnca,np,nr,drmnc(mn,nr),dy) 
c
            do i=1,np
               nra(i)=nr-np-1+i
               srmnca(i)=srmnc(mn,nr-np-1+i)
            enddo
c     
            call polint(nra,srmnca,np,nr,srmnc(mn,nr),dy) 
c
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
c               DZMNS(MN,NR)=(YZBS(NSRMAX)-YZBS(NSRMAX-1))
c     &                     /(XS(NSRMAX)-XS(NSRMAX-1))
c               SZMNS(MN,NR)=YZBS(NSRMAX)
c     &                     +DZMNS(MN,NR)*(XRHO(NR)-XS(NSRMAX))
c
c
            do i=1,np
               nra(i)=nr-np-1+i
               dzmnsa(i)=dzmns(mn,nr-np-1+i)
            enddo
c     
            call polint(nra,dzmnsa,np,nr,dzmns(mn,nr),dy) 
c
            do i=1,np
               nra(i)=nr-np-1+i
               szmnsa(i)=szmns(mn,nr-np-1+i)
            enddo
c     
            call polint(nra,szmnsa,np,nr,szmns(mn,nr),dy) 
c
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
      DPH=2.D0*PI/(NHC*NHHMAX)
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
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
         PH=DPH*(NHH-1)
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
         RPST(  NTH,NHH,NR)=RPSS
         ZPST(  NTH,NHH,NR)=ZPSS
C
         IF(NR.NE.1) THEN
            DRS=DRS/(2.D0*PSIPA*XRHO(NR))
            DZS=DZS/(2.D0*PSIPA*XRHO(NR))
            RG11(NTH,NHH,NR)=(DRS *DRS +DZS *DZS )*XRHO(NR)**2
            RG12(NTH,NHH,NR)=(DRS *DRTH+DZS *DZTH)
            RG13(NTH,NHH,NR)=(DRS *DRPH+DZS *DZPH)*XRHO(NR)
            RG22(NTH,NHH,NR)=(DRTH*DRTH+DZTH*DZTH)/XRHO(NR)**2
            RG23(NTH,NHH,NR)=(DRTH*DRPH+DZTH*DZPH)/XRHO(NR)
            RG33(NTH,NHH,NR)= DRPH*DRPH+DZPH*DZPH+DPHPH*DPHPH
            RJ  (NTH,NHH,NR)=DPHPH*(DRS*DZTH-DZS*DRTH)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         RG11(NTH,NHH,1)=RG11(NTH,NHH,2)
         RG12(NTH,NHH,1)=RG12(NTH,NHH,2)
         RG13(NTH,NHH,1)=RG13(NTH,NHH,2)
         RG22(NTH,NHH,1)=RG22(NTH,NHH,2)
         RG23(NTH,NHH,1)=RG23(NTH,NHH,2)
         RG33(NTH,NHH,1)=RG33(NTH,NHH,2)
         RJ(NTH,NHH,1)  =RJ(NTH,NHH,2)
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
      integer np,i
      real*8 nra,bstha,bspha,qpsa
      parameter (np=3)
      dimension nra(np),bstha(np),bspha(np),qpsa(np)
      DIMENSION BSUS(NSRM),BSVS(NSRM)
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX
            BSUS(NSR)=BSU(MN+(NSR-1)*MNMAX)
            BSVS(NSR)=BSV(MN+(NSR-1)*MNMAX)
         ENDDO
C         IF(NRANK.EQ.0) THEN
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
         BSTHSV(MN)=BSUS(NSRMAX-1)
         BSTHSD(MN)=(BSUS(NSRMAX-1)-BSUS(NSRMAX-5))
     &             /(XS(NSRMAX-1)  -XS(NSRMAX-5))
         BSPHSV(MN)=BSVS(NSRMAX)
         BSPHSD(MN)=(BSVS(NSRMAX)-BSVS(NSRMAX-5))
     &             /(XS(NSRMAX)  -XS(NSRMAX-5))
         DO NR=1,NRMAX+1
C
C               BSTHL=BSTHSV(MN)+BSTHSD(MN)*(XRHO(NR)-XS(NSRMAX))
C               BSPHL=BSPHSV(MN)+BSPHSD(MN)*(XRHO(NR)-XS(NSRMAX))
C
            CALL SPL1DF(XRHO(NR),BSTHL,XS,U3(1,1,MN),NSRMAX,IERR)
            CALL SPL1DF(XRHO(NR),BSPHL,XS,U4(1,1,MN),NSRMAX,IERR)
            BSTH(MN,NR)=BSTHL
            BSPH(MN,NR)=BSPHL
         ENDDO
      ENDDO
c
      do mn=1,mnmax
         do nr=1,nrmax+1
            if((xrho(nr)-xs(nsrmax).gt.0.d0).or.(nr.eq.84)) then
c
                  do i=1,np
                     nra(i)=nr-np-1+i
                     bstha(i)=bsth(mn,nr-np-1+i)
                  enddo
c     
                  call polint(nra,bstha,np,nr,bsth(mn,nr),dy) 
c     
               do i=1,np
                  nra(i)=nr-np-1+i
                  bspha(i)=bsph(mn,nr-np-1+i)
               enddo
c     
               call polint(nra,bspha,np,nr,bsph(mn,nr),dy) 
c          
            else
            endif
         enddo
      enddo
c
c$$$      do nr=1,nrmax+1
c$$$            print*,nr,bsth(1,nr),bsph(1,nr)
c$$$      enddo
C
C      ***** CULCULATE MAGNETIC FIELD *****
C 
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/(NHC*NHHMAX)
C      DO NR=2,NRMAX+1
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
         TH=DTH*(NTH-1)
         PH=DPH*(NHH-1)
C
         SBTH=0.D0
         SBPH=0.D0
         DO MN=1,MNMAX
            RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
            RCOS=COS(XM(MN)*TH-XN(MN)*PH)
            SBTH=SBTH+BSTH(MN,NR)*RCOS
            SBPH=SBPH+BSPH(MN,NR)*RCOS
         ENDDO
         BFLD(2,NTH,NHH,NR)=SBTH
         BFLD(3,NTH,NHH,NR)=SBPH
      ENDDO
      ENDDO
      ENDDO
C
C      DO NHH=1,NHHMAX
C      DO NTH=1,NTHMAX
C         BFLD(2,NTH,NHH,1)=(4*BFLD(2,NTH,NHH,2)-BFLD(2,NTH,NHH,3))/3
C         BFLD(3,NTH,NHH,1)=(4*BFLD(3,NTH,NHH,2)-BFLD(3,NTH,NHH,3))/3
C      ENDDO
C      ENDDO
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
         DO NHH=1,NHHMAX
            RSU(NSU,NHH)=0.D0
            ZSU(NSU,NHH)=0.D0
            TH=DTHU*(NSU-1)
            PH=DPH*(NHH-1)
            DO MN=1,MNMAX
               RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
               RCOS=COS(XM(MN)*TH-XN(MN)*PH)
               RSU(NSU,NHH)=RSU(NSU,NHH)+RMNCC(MN,NSRMAX)*RCOS
               ZSU(NSU,NHH)=ZSU(NSU,NHH)+ZMNSS(MN,NSRMAX)*RSIN
            ENDDO
         ENDDO
      ENDDO
C
      NSWMAX=31
      DTHW=2.D0*PI/(NSWMAX-1)
      DO NSW=1,NSWMAX
         DO NHH=1,NHHMAX
            RSW(NSW,NHH)=0.D0
            ZSW(NSW,NHH)=0.D0
            TH=DTHW*(NSW-1)
            PH=DPH*(NHH-1)
            DO MN=1,MNMAX
               RSIN=SIN(XM(MN)*TH-XN(MN)*PH)
               RCOS=COS(XM(MN)*TH-XN(MN)*PH)
               RSW(NSW,NHH)=RSW(NSW,NHH)+SRMNC(MN,NRMAX+1)*RCOS
               ZSW(NSW,NHH)=ZSW(NSW,NHH)+SZMNS(MN,NRMAX+1)*RSIN
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
c         IF(IERR.EQ.2) THEN
c            RIOTASL=RIOTAS(NSRMAX)+(RIOTAS(NSRMAX)-RIOTAS(NSRMAX-1))
c     &                      /(XS(NSRMAX)-XS(NSRMAX-1))
c     &                      *(XRHO(NR)-XS(NSRMAX))
c         ENDIF
C         QPS(NR)=2.D0*PI/RIOTASL
            IF(XRHO(NR)-XS(NSRMAX).LT.0.D0) THEN!
               QPS(NR)=1.D0/RIOTASL
            else
c
                  do i=1,np
                     nra(i)=nr-np-1+i
                     qpsa(i)=qps(nr-np-1+i)
                  enddo
c     
                  call polint(nra,qpsa,np,nr,qps(nr),dy) 
c     
               endif
      ENDDO                     
C
C     ***** COMPUTE R,Z MAGNETIC AXES *****
C
      RGMIN=RSW(1,1)
      RGMAX=RSW(1,1)
      ZGMIN=ZSW(1,1)
      ZGMAX=ZSW(1,1)
      DO NHH=1,NHHMAX
         DO NSW=1,NSWMAX
            RGMIN=MIN(RGMIN,RSW(NSW,NHH))
            RGMAX=MAX(RGMAX,RSW(NSW,NHH))
            ZGMIN=MIN(ZGMIN,ZSW(NSW,NHH))
            ZGMAX=MAX(ZGMAX,ZSW(NSW,NHH))
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
            CALL SPL1DF(XRHOL,BSTHL,XS,U3(1,1,MN),NSRMAX,IERR)
            CALL SPL1DF(XRHOL,BSPHL,XS,U4(1,1,MN),NSRMAX,IERR)
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
C      QPSL=2.D0*PI/RIOTASL
      QPSL=1.D0/RIOTASL
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
