C     $Id$
C
C    ***** INTERFACE FOR BOOZER COORDINATE EQUILIBRIUM *****
C
      SUBROUTINE WMXRZB(IERR)
C
      INCLUDE 'wmcomm.inc'
C     
      IERR=0
C
      CALL WMHBRD
      CALL WMHBRZ
      CALL WMHBMT
      CALL WMHBBB
C
      RETURN
      END
C
C    ****** READ ASCII BOOZER FILE ******
C
      SUBROUTINE WMHBRD
C
C       Definition in Boozer coordinates on half-mesh
C
C	ns 	   number of Boozer radial grids
C	nmboz 	   number of Boozer harmonics
C	xmboz 	   poloidal mode number
C	xnboz 	   toroidal mode number
C
C	B spectrum
C	    B(i,t,z) =	   sum bbozh(m,i) cos(-t*mboz(m) + z*nboz(m))
C
C	Transformation from Boozer to cylindrical coordinates
C	  ( psi, theta, zeta ) -> ( R, Phi, Z )
C
C	    R(i,t,z) =	   sum rbozh(m,i) cos(-t*mboz(m) + z*nboz(m))
C	    P(i,t,z) = z + sum pbozh(m,i) sin(-t*mboz(m) + z*nboz(m))
C	    Z(i,t,z) =	   sum zbozh(m,i) sin(-t*mboz(m) + z*nboz(m))
C
C	Surface quantities ( based on vmec normalization )
C	    psi   : toroidal flux within a flux surface/(2*pi)
C	    shalf : normalized psi as psi/psi_edge*pi, D(shalf)=[0,pi]
C	    xiota : rotational transform
C	    wjs	  : toroidal current within a flux surface
C	    wis	  : poloidal current without a flux surface
C
C
      USE libfio
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      NFL=11
      CALL FROPEN(NFL,KNAMEQ,1,0,'EQ',IERR)
      IF(IERR.NE.0) STOP
C
 9997 format(3(i7))
 9996 format(3(1pd25.16))
 9995 format(2(i7))
      rewind NFL
      read(NFL,9997) MNMAX, NSRMAX, nfp
7001  format(3i4)
C
      IERR=0
      IF(MNMAX.GT.NMNM) THEN
         WRITE(6,'(A,2I8)') 
     &        'XX WMHBRD: MNMAX.GT.NMNM: MNMAX,NMNM=',MNMAX,NMNM
         IERR=1
      ENDIF
      IF(NSRMAX.GT.NSRM) THEN
         WRITE(6,'(A,2I8)') 
     &        'XX WMHBRD: NSRMAX.GT.NSRM: NSRMAX,NSRM=',NSRMAX,NSRM
         IERR=2
      ENDIF
      IF(IERR.NE.0) STOP
C

      read(NFL,9996) (shalf(i),i=2,nsrmax+1)
      read(NFL,9996) (RIOTAS(i),i=2,nsrmax+1)
      read(NFL,9996) (wjs(i),i=2,nsrmax+1)
      read(NFL,9996) (wis(i),i=2,nsrmax+1)
      
C
      DO MN=1,MNMAX
         read(NFL,9995) M,N
         read(NFL,9996,ERR=8001) (BBOZH(MN,i),i=2,nsrmax+1)
         XM(MN)=DBLE(M)
         XN(MN)=DBLE(N)
      ENDDO
C
      DO MN=1,MNMAX
         read(NFL,9995) M,N
         read(NFL,9996,ERR=8001) (rbozh(MN,i),i=2,nsrmax+1)
      ENDDO
C
      DO MN=1,MNMAX
         read(NFL,9995) M,N
         read(NFL,9996,ERR=8001) (zbozh(MN,i),i=2,nsrmax+1)
      ENDDO
C
      DO MN=1,MNMAX
         read(NFL,9995) M,N
         read(NFL,9996,ERR=8001) (pbozh(MN,i),i=2,nsrmax+1)
      ENDDO
C
      CLOSE(NFL)
C
      nsrchk = nsrmax/2
      rchichk = 2.d0*PI/90.0d0
      zsum  =  0.0d0
      do mn = 1, mnmax
C         print *, 'mn,xm(mn),zbozh=',mn, xm(mn),zbozh(mn,nsrchk)
         zsum  =  zsum + zbozh(mn,nsrchk)*sin(-xm(mn)*rchichk)
      enddo
c
C         print *,' nsrchk,chichk=', nsrchk, rchichk
C         print *,' zsum=', zsum
c
        if( zsum .gt. 0.0d0 ) then
c.....    counterclockwise
          drthta =  1.0d0
          print *,' --------------------------------------'
          print *,'     LHS ( theta is counterclockwise ) '
          print *,' --------------------------------------'
        else
c.....    clockwise
          drthta = -1.0d0
          print *,' --------------------------------------'
          print *,'     RHS ( theta is clockwise ) '
          print *,' --------------------------------------'
        endif
c

c+++++ updated on  9/29 94 by N^2 : for left -> right start
c
c      - theta -> theta in VMEC
c
      if(drthta .gt. 0.0d0 ) then
         do i = 2, nsrmax+1
            riotas(i) = - riotas(i)
            wis  (i) = - wis  (i)
         enddo
c
         do m = 1, mnmax
            xm(m) = - xm(m)
         enddo
c
          print *,' --------------------------------------'
          print *,'     LHS -> RHS ( theta is clockwise ) '
          print *,' --------------------------------------'
          drthta = -1.0d0
      endif
C
      WRITE(6,'(A,4I10)') 'NSRMAX,MNMAX,NSRM,NMNM=',
     &                     NSRMAX,MNMAX,NSRM,NMNM
c
c+++++ updated on  9/29 94 by N^2 : for left -> right end
c
C
C     ----- Setup normalized poloidal flux for file data -----
C
C        XS(1)=0                     XSH(1)=0
C        XS(2)=DELPSI                XSH(2)=0.5*DELPSI
C
C        XS(NSRMAX+1)=1.D0           XSH(NSRMAX+1)=1.D0-0.5*DELPSI
C        XS(NSRMAX+2)=(RB/RA)**2     XSH(NSRMAX+2)=1.D0+0.5*DELPSI
C                                    XSH(NSRMAX+3)=(RB/RA)**2
C
      DELPSI=1.D0/NSRMAX
      DO NSR=1,NSRMAX+1
         XS(NSR)=DELPSI*(NSR-1)
      ENDDO
      XS(NSRMAX+2)=(RB/RA)**2
      XSH(1)=0.D0
      DO NSR=2,NSRMAX+2
         XSH(NSR)=DELPSI*(NSR-1.5D0)
      ENDDO
      XSH(NSRMAX+3)=(RB/RA)**2
C
C     ----- Extraporate to axis and wall radius -----
C
      RR=0.D0
      DO MN=1,MNMAX
         IF(XM(MN).EQ.0.D0) THEN
            BBOZH(MN,1)=(3*BBOZH(MN,2)-BBOZH(MN,3))/2
            RBOZH(MN,1)=(3*RBOZH(MN,2)-RBOZH(MN,3))/2
            ZBOZH(MN,1)=(3*ZBOZH(MN,2)-ZBOZH(MN,3))/2
            PBOZH(MN,1)=(3*PBOZH(MN,2)-PBOZH(MN,3))/2
            IF(XN(MN).EQ.0.D0) RR=RBOZH(MN,1)
         ELSE
            BBOZH(MN,1)=0.D0
            RBOZH(MN,1)=0.D0
            ZBOZH(MN,1)=0.D0
            PBOZH(MN,1)=0.D0
         ENDIF
      ENDDO
      SHALF(1)=0.D0
      RIOTAS(1)=(3*RIOTAS(2)-RIOTAS(3))/2
      WIS(1)=(3*WIS(2)-WIS(3))/2
      WJS(1)=0.D0
C
      FACTOR=(XSH(NSRMAX+2)-XSH(NSRMAX+1))
     &      /(XSH(NSRMAX+1)-XSH(NSRMAX))
      DO MN=1,MNMAX
         BBOZH(MN,NSRMAX+2)=BBOZH(MN,NSRMAX+1)
     &             +FACTOR*(BBOZH(MN,NSRMAX+1)-BBOZH(MN,NSRMAX))
         RBOZH(MN,NSRMAX+2)=RBOZH(MN,NSRMAX+1)
     &             +FACTOR*(RBOZH(MN,NSRMAX+1)-RBOZH(MN,NSRMAX))
         ZBOZH(MN,NSRMAX+2)=ZBOZH(MN,NSRMAX+1)
     &             +FACTOR*(ZBOZH(MN,NSRMAX+1)-ZBOZH(MN,NSRMAX))
         PBOZH(MN,NSRMAX+2)=PBOZH(MN,NSRMAX+1)
     &             +FACTOR*(PBOZH(MN,NSRMAX+1)-PBOZH(MN,NSRMAX))
      ENDDO
      SHALF(NSRMAX+2)=SHALF(NSRMAX+1)
     &       +FACTOR*(SHALF(NSRMAX+1)-SHALF(NSRMAX))
      RIOTAS(NSRMAX+2)=RIOTAS(NSRMAX+1)
     &       +FACTOR*(RIOTAS(NSRMAX+1)-RIOTAS(NSRMAX))
C
      FACTOR=(1.D0-XSH(NSRMAX+1))
     &      /(XSH(NSRMAX+1)-XSH(NSRMAX))
      WJSSURF=WJS(NSRMAX+1)
     &       +FACTOR*(WJS(NSRMAX+1)-WJS(NSRMAX))
      WISSURF=WIS(NSRMAX+1)
     &       +FACTOR*(WIS(NSRMAX+1)-WIS(NSRMAX))
      WJS(NSRMAX+2)=WJSSURF
      WIS(NSRMAX+2)=WISSURF
C
      FACTOR=(XSH(NSRMAX+3)-XSH(NSRMAX+1))
     &      /(XSH(NSRMAX+1)-XSH(NSRMAX))
      DO MN=1,MNMAX
         BBOZH(MN,NSRMAX+3)=BBOZH(MN,NSRMAX+1)
     &             +FACTOR*(BBOZH(MN,NSRMAX+1)-BBOZH(MN,NSRMAX))
         RBOZH(MN,NSRMAX+3)=RBOZH(MN,NSRMAX+1)
     &             +FACTOR*(RBOZH(MN,NSRMAX+1)-RBOZH(MN,NSRMAX))
         ZBOZH(MN,NSRMAX+3)=ZBOZH(MN,NSRMAX+1)
     &             +FACTOR*(ZBOZH(MN,NSRMAX+1)-ZBOZH(MN,NSRMAX))
         PBOZH(MN,NSRMAX+3)=PBOZH(MN,NSRMAX+1)
     &             +FACTOR*(PBOZH(MN,NSRMAX+1)-PBOZH(MN,NSRMAX))
      ENDDO
      WJS(NSRMAX+3)=WJSSURF
      WIS(NSRMAX+3)=WISSURF
c
      BB0=wissurf/RR
      factor=DABS(RR*BB/wissurf)
      WRITE(6,'(A,1P4E12.4)')
     &     'RR,BB0,BB,factor=',RR,BB0,BB,factor
c
      do nsr=1,nsrmax+3
         do mn=1,mnmax
            bbozh(mn,nsr)=bbozh(mn,nsr)*factor
         enddo
      enddo
         
C
C      DO NSR=1,NSRMAX+3
C         WRITE(6,'(I5,2X,1P5E12.5)') NSR,SHALF(NSR),XS(NSR),XSH(NSR),
C     &                               WIS(NSR),WJS(NSR)
C      ENDDO
C
      RETURN
C
 8001 DO NSR=1,NSRMAX+1
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,BBOZH(MN,NSR)
      ENDDO
      STOP
 8002 DO NSR=1,NSRMAX+1
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,RBOZH(MN,NSR)
      ENDDO
      STOP
 8003 DO NSR=1,NSRMAX+1
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,ZBOZH(MN,NSR)
      ENDDO
      STOP
 8004 DO NSR=1,NSRMAX+1
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,PBOZH(MN,NSR)
      ENDDO
      STOP
      END
C
C    ****** CALCULATE R AND Z ******
C
      SUBROUTINE WMHBRZ
C
      USE libspl1d
      USE libpol
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      INTEGER I,NP
      REAL*8 SBMNCA,SRMNCA,SZMNSA,SPMNSA,DBMNCA,DRMNCA,DZMNSA,DPMNSA
      PARAMETER (NP=3)
      DIMENSION BSPL(NSRMP),RSPL(NSRMP),ZSPL(NSRMP),PSPL(NSRMP)
      DIMENSION NRA(NP),SBMNCA(NP),SRMNCA(NP),SZMNSA(NP),SPMNSA(NP)
      DIMENSION DBMNCA(NP),DRMNCA(NP),DZMNSA(NP),DPMNSA(NP)
C
C     ***** DEFINE XRHO, XR AND XSHRHO *****
C
C     XRHO   : NORMALIZED RADIUS         (poloidal flux)
C     XR     : AVERAGE RADIUS            (poloidal flux)
C     XSHRHO : NORMALIZED RADIUS OF DATA (poloidal flux)
C
      RHOB=RB/RA
      DRHO=RHOB/NRMAX
      DO NR=1,NRMAX+1
         XRHO(NR)=DRHO*(NR-1)
         XR(NR)  =RB*XRHO(NR)
      ENDDO
C
      DO NSR=1,NSRMAX+3
         XSHRHO(NSR)=SQRT(XSH(NSR))
      ENDDO
C
C      ***** PSIP (normalized) *****
C
      DO NR=1,NRMAX+1
         PSIP(NR)=XRHO(NR)**2
      ENDDO
C     
C      ***** SPLINE BBOZH(S),RBOZH(S),ZBOZH(S),PBOZH(S) *****
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX+3
            BSPL(NSR)=BBOZH(MN,NSR)
            RSPL(NSR)=RBOZH(MN,NSR)
            ZSPL(NSR)=ZBOZH(MN,NSR)
            PSPL(NSR)=PBOZH(MN,NSR)
         ENDDO
         IF(XM(MN).EQ.0.D0) THEN
            FX1(1)=0.D0
            CALL SPL1D(XSHRHO,BSPL,FX1,U1(1,1,MN),NSRMAX+3,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: BBOZH'
            FX2(1)=0.D0
            CALL SPL1D(XSHRHO,RSPL,FX2,U2(1,1,MN),NSRMAX+3,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: RBOZH'
            FX3(1)=0.D0
            CALL SPL1D(XSHRHO,ZSPL,FX3,U3(1,1,MN),NSRMAX+3,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: ZBOZH'
            FX3(1)=0.D0
            CALL SPL1D(XSHRHO,PSPL,FX4,U4(1,1,MN),NSRMAX+3,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: PBOZH'
         ELSE
            CALL SPL1D(XSHRHO,BSPL,FX1,U1(1,1,MN),NSRMAX+3,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: BBOZH'
            CALL SPL1D(XSHRHO,RSPL,FX2,U2(1,1,MN),NSRMAX+3,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: RBOZH'
            CALL SPL1D(XSHRHO,ZSPL,FX3,U3(1,1,MN),NSRMAX+3,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: ZBOZH'
            CALL SPL1D(XSHRHO,PSPL,FX4,U4(1,1,MN),NSRMAX+3,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: PBOZH'
         ENDIF
C
         DO NR=1,NSRMAX
            CALL SPL1DD(XRHO(NR),SBMNC(MN,NR),DBMNC(MN,NR),
     &                  XSHRHO,U1(1,1,MN),NSRMAX+3,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: BMNC: NR=',NR
            CALL SPL1DD(XRHO(NR),SRMNC(MN,NR),DRMNC(MN,NR),
     &                  XSHRHO,U2(1,1,MN),NSRMAX+3,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: RMNC: NR=',NR
            CALL SPL1DD(XRHO(NR),SZMNS(MN,NR),DZMNS(MN,NR),
     &                  XSHRHO,U3(1,1,MN),NSRMAX+3,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: ZMNC: NR=',NR
            CALL SPL1DD(XRHO(NR),SPMNS(MN,NR),DPMNS(MN,NR),
     &                  XSHRHO,U4(1,1,MN),NSRMAX+3,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: PMNC: NR=',NR
         ENDDO
         DO NR=NSRMAX+1,NRMAX+1
            DO I=1,NP
               NRA(I)=NR-2*NP-1+2*I
            ENDDO     
            DO I=1,NP
               SBMNCA(I)=SBMNC(MN,NR-2*NP-1+2*I)
            ENDDO     
            CALL POLINTN(NRA,SBMNCA,NP,NR,SBMNC(MN,NR)) 
            DO I=1,NP
               DBMNCA(I)=DBMNC(MN,NR-2*NP-1+2*I)
            ENDDO     
            CALL POLINTN(NRA,DBMNCA,NP,NR,DBMNC(MN,NR)) 
C     
            DO I=1,NP
               SRMNCA(I)=SRMNC(MN,NR-2*NP-1+2*I)
            ENDDO
            CALL POLINTN(NRA,SRMNCA,NP,NR,SRMNC(MN,NR)) 
            DO I=1,NP
               DRMNCA(I)=DRMNC(MN,NR-2*NP-1+2*I)
            ENDDO     
            CALL POLINTN(NRA,DRMNCA,NP,NR,DRMNC(MN,NR)) 
C     
            DO I=1,NP
               SZMNSA(I)=SZMNS(MN,NR-2*NP-1+2*I)
            ENDDO
            CALL POLINTN(NRA,SZMNSA,NP,NR,SZMNS(MN,NR)) 
            DO I=1,NP
               DZMNSA(I)=DZMNS(MN,NR-2*NP-1+2*I)
            ENDDO     
            CALL POLINTN(NRA,DZMNSA,NP,NR,DZMNS(MN,NR)) 
C     
            DO I=1,NP
               SPMNSA(I)=SPMNS(MN,NR-2*NP-1+2*I)
            ENDDO
            CALL POLINTN(NRA,SPMNSA,NP,NR,SPMNS(MN,NR)) 
            DO I=1,NP
               DPMNSA(I)=DPMNS(MN,NR-2*NP-1+2*I)
            ENDDO     
            CALL POLINTN(NRA,DPMNSA,NP,NR,DPMNS(MN,NR)) 
         ENDDO
      ENDDO
C
C      DO NR=1,NRMAX+1
C         WRITE(6,'(A,1P5E12.4)')
C     &        'XRHO,SBMNC=',XRHO(NR),SBMNC(1,NR),SBMNC(2,NR),
C     &        SBMNC(3,NR),SBMNC(4,NR)
C      ENDDO
C
      RETURN
      END
C
C    ****** CALCULTE METRIC ******
C
      SUBROUTINE WMHBMT
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C      ***** CULCULATE METRIC TENSOR AND JACOBIAN*****
C 
      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/(NHC*NHHMAX)
      PSIPA=1.D0
C
      DO NR=1,NRMAX+1
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
C
         TH=DTH*(NTH-1)
         PH=DPH*(NHH-1)
C
         BPSS=0.D0
         RPSS=0.D0
         ZPSS=0.D0
         PPSS=PH
         DRS=0.D0
         DZS=0.D0
         DPS=0.D0
         DRTH=0.D0
         DZTH=0.D0
         DPTH=0.D0
         DRPH=0.D0
         DZPH=0.D0
         DPPH=1.D0
C
         DO MN=1,MNMAX
            RSIN=SIN(-XM(MN)*TH+XN(MN)*PH)
            RCOS=COS(-XM(MN)*TH+XN(MN)*PH)
            BPSS =BPSS +       SBMNC(MN,NR)*RCOS
            RPSS =RPSS +       SRMNC(MN,NR)*RCOS
            ZPSS =ZPSS +       SZMNS(MN,NR)*RSIN
            PPSS =PPSS +       SPMNS(MN,NR)*RSIN
            DRS  =DRS  +       DRMNC(MN,NR)*RCOS
            DZS  =DZS  +       DZMNS(MN,NR)*RSIN
            DPS  =DPS  +       DPMNS(MN,NR)*RSIN
            DRTH =DRTH +XM(MN)*SRMNC(MN,NR)*RSIN
            DZTH =DZTH -XM(MN)*SZMNS(MN,NR)*RCOS
            DPTH =DPTH -XM(MN)*SPMNS(MN,NR)*RCOS
            DRPH =DRPH -XN(MN)*SRMNC(MN,NR)*RSIN
            DZPH =DZPH +XN(MN)*SZMNS(MN,NR)*RCOS
            DPPH =DPPH +XN(MN)*SPMNS(MN,NR)*RCOS
         ENDDO
C
C  *****  DRS,DZS,DPHS,DRTH,DZTH,DPHTH,DRPH,DZPH,DPHPH GRAPH *****
C            
         BPST(  NTH,NHH,NR)=BPSS
         RPST(  NTH,NHH,NR)=RPSS
         ZPST(  NTH,NHH,NR)=ZPSS
         PPST(  NTH,NHH,NR)=PPSS
C
         IF(NR.NE.1) THEN
            DRS=DRS/(2.D0*PSIPA*XRHO(NR))
            DZS=DZS/(2.D0*PSIPA*XRHO(NR))
            DPS=DPS/(2.D0*PSIPA*XRHO(NR))
            RG11(NTH,NHH,NR)=DRS *DRS +DZS *DZS +DPS *DPS *RPSS**2
            RG12(NTH,NHH,NR)=DRS *DRTH+DZS *DZTH+DPS *DPTH*RPSS**2
            RG13(NTH,NHH,NR)=DRS *DRPH+DZS *DZPH+DPS *DPPH*RPSS**2
            RG22(NTH,NHH,NR)=DRTH*DRTH+DZTH*DZTH+DPTH*DPTH*RPSS**2
            RG23(NTH,NHH,NR)=DRTH*DRPH+DZTH*DZPH+DPTH*DPPH*RPSS**2
            RG33(NTH,NHH,NR)=DRPH*DRPH+DZPH*DZPH+DPPH*DPPH*RPSS**2
C            RJ  (NTH,NHH,NR)=(DRS *DZTH-DZS *DRTH)*DPPH*RPSS
C     &                      +(DZS *DPTH-DPS *DZTH)*DRPH*RPSS
C     &                      +(DPS *DRTH-DRS *DPTH)*DZPH*RPSS
            RJ  (NTH,NHH,NR)=(DRS *DPTH-DPS *DRTH)*DZPH*RPSS
     &                      +(DPS *DZTH-DZS *DPTH)*DRPH*RPSS
     &                      +(DZS *DRTH-DRS *DZTH)*DPPH*RPSS
C
            RG11(NTH,NHH,NR)=RG11(NTH,NHH,NR)*XRHO(NR)**2
            RG13(NTH,NHH,NR)=RG13(NTH,NHH,NR)*XRHO(NR)
            RG22(NTH,NHH,NR)=RG22(NTH,NHH,NR)/XRHO(NR)**2
            RG23(NTH,NHH,NR)=RG23(NTH,NHH,NR)/XRHO(NR)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      do nr=1,nrmax+1
         print*,nr,rj(1,1,nr),drmnc(1,nr)
      enddo
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
      SUBROUTINE WMHBBB
C
      USE libspl1d
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
      DIMENSION SRMNCL(NMNM),SZMNSL(NMNM)
C
C     ***** INTERPORATE RIOTAS *****
C
      FX5(1)=0.D0
      CALL SPL1D(XSHRHO,RIOTAS,FX5,U5,NSRMAX+3,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMHBBB: SPL1D: RIOTAS'
C
      DO NR=1,NRMAX+1
         CALL SPL1DF(XRHO(NR),SIOTA(NR),
     &               XSHRHO,U5,NSRMAX+3,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX WMHBBB: SPL1DD: IOTA: NR=',NR
         QPS(NR)=1.D0/SIOTA(NR)
      ENDDO
C
      DO NHH=1,NHHMAX
      DO NTH=1,NTHMAX
      DO NR=1,NRMAX+1
         FACT=  RG22(NTH,NHH,NR)*XRHO(NR)**2*SIOTA(NR)**2
     &       +2*RG23(NTH,NHH,NR)*XRHO(NR)   *SIOTA(NR)
     &       +  RG33(NTH,NHH,NR)
         BFLD(3,NTH,NHH,NR)=BPST(NTH,NHH,NR)/SQRT(FACT)
         BFLD(2,NTH,NHH,NR)=SIOTA(NR)*BFLD(3,NTH,NHH,NR)
      ENDDO
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
               FACTN=(1.D0-FEDGE)*(1.D0-RHOL**PROFN1(1))**PROFN2(1)+FEDGE
            ENDIF
            PT=(PTPR(1)+2*PTPP(1))/3.D0
            IF(PT.LE.0.D0) THEN
               FACTT=0.D0
            ELSE
               FEDGE=PTS(1)/PT
               FACTT=(1.D0-FEDGE)*(1.D0-RHOL**PROFT1(1))**PROFT2(1)+FEDGE
            ENDIF
            PPS(NR)=P0*FACTN*FACTT
         ELSE
            PPS(NR)=0.D0
         ENDIF
      ENDDO
C
      NSUMAX=31
      DTHU=2.D0*PI/(NSUMAX-1)
      DPH=2.D0*PI/(NHC*NHHMAX)
C
      DO MN=1,MNMAX
         CALL SPL1DF(1.D0,SRMNCL(MN),XSHRHO,U2(1,1,MN),NSRMAX+3,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX WMHBBB: SPL1DF: SRMNCL'
         CALL SPL1DF(1.D0,SZMNSL(MN),XSHRHO,U3(1,1,MN),NSRMAX+3,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX WMHBBB: SPL1DF: SZMNCL'
      ENDDO
C
      DO NSU=1,NSUMAX
         DO NHH=1,NHHMAX
            RSU(NSU,NHH)=0.D0
            ZSU(NSU,NHH)=0.D0
            TH=DTHU*(NSU-1)
            PH=DPH*(NHH-1)
            DO MN=1,MNMAX
               RSIN=SIN(-XM(MN)*TH+XN(MN)*PH)
               RCOS=COS(-XM(MN)*TH+XN(MN)*PH)
               RSU(NSU,NHH)=RSU(NSU,NHH)+SRMNCL(MN)*RCOS
               ZSU(NSU,NHH)=ZSU(NSU,NHH)+SZMNSL(MN)*RSIN
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
               RSIN=SIN(-XM(MN)*TH+XN(MN)*PH)
               RCOS=COS(-XM(MN)*TH+XN(MN)*PH)
               RSW(NSW,NHH)=RSW(NSW,NHH)+SRMNC(MN,NRMAX+1)*RCOS
               ZSW(NSW,NHH)=ZSW(NSW,NHH)+SZMNS(MN,NRMAX+1)*RSIN
            ENDDO
         ENDDO
      ENDDO
C
C     ***** COMPUTE FIGURE BOUNDARY and MAJOR RADIUS *****
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
C
      RETURN
      END
C
C    ****** DRAW BOOZER DATA ******
C
      SUBROUTINE WMGBOOZ
C
      USE libgrf
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      DIMENSION GX(NSRMP),GY(NSRMP,NMNM)
      DIMENSION GXR(NRM),GYR(NRM,NMNM)
C
C     ----- graphics with psit axis -----
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(BBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/BBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(RBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/RBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(ZBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/ZBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(PBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/PBOZH/',0)
      CALL PAGEE
C
C     ----- graphics with rho axis -----
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(BBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/BBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(RBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/RBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(ZBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/ZBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+3
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(PBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRMP,NSRMAX+3,MNMAX,'/PBOZH/',0)
      CALL PAGEE
C
C     ----- graphics with xrho axis -----
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(SBMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/SBMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(DBMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/DBMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(SRMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/SRMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(DRMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/DRMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(SZMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/SZMNS/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(DZMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/DZMNS/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(SPMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/SPMNS/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GXR(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GYR(NR,MN)=GUCLIP(DPMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GXR,GYR,NRM,NRMAX+1,MNMAX,'/DPMNS/',0)
      CALL PAGEE
C
      RETURN
      END
