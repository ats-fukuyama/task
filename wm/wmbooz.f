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
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      DIMENSION shalf(NSRM),wjs(NSRM),wis(NSRM)
C
C 1002    WRITE(6,*) 'READ FILE NAME ?'
C         READ(5,'(A80)',ERR=1002,END=9000) FNAME
C          FNAME='wout'
C
      NFL=11
C      KNAMEQ='m08n06.X1.1320_b_0.00.asc'
      CALL FROPEN(NFL,KNAMEQ,1,0,'EQ',IERR)
C
      rewind NFL
      read(NFL,*) MNMAX, NSRMAX, nfp
      NSRMAX=NSRMAX+1
C
      read(NFL,*) (shalf(NSR),NSR=2,NSRMAX)
      read(NFL,*) (RIOTAS(NSR),NSR=2,NSRMAX)
C
      read(NFL,*) (wjs(NSR),NSR=2,NSRMAX)
      read(NFL,*) (wis(NSR),NSR=2,NSRMAX)
C
      DO MN=1,MNMAX
         read(NFL,*) M,N
         read(NFL,*,ERR=8001) (BBOZH(MN,NSR),NSR=2,NSRMAX)
         XM(MN)=DBLE(M)
         XN(MN)=DBLE(N)
      ENDDO
C
      DO MN=1,MNMAX
         read(NFL,*) M,N
         read(NFL,*,ERR=8002) (RBOZH(MN,NSR),NSR=2,NSRMAX)
      ENDDO
C
      DO MN=1,MNMAX
         read(NFL,*) M,N
         read(NFL,*,ERR=8003) (ZBOZH(MN,NSR),NSR=2,NSRMAX)
      ENDDO
C
      DO MN=1,MNMAX
         read(NFL,*) M,N
         read(NFL,*,ERR=8004) (PBOZH(MN,NSR),NSR=2,NSRMAX)
      ENDDO
C
      CLOSE(NFL)
C
C     ----- Setup normalized radius for file data -----
C
      DELPSI=1.D0/(NSRMAX-1)
      DO NSR=1,NSRMAX
         XS(NSR)=DELPSI*(NSR-1)
      ENDDO
      XSH(1)=0.D0
      DO NSR=2,NSRMAX
         XSH(NSR)=DELPSI*(NSR-1.5D0)
      ENDDO
      XSH(NSRMAX+1)=(RB/RA)**2
C
C     ----- Extraporate to axis and wall radius -----
C
      DO MN=1,MNMAX
         IF(XM(MN).EQ.0.D0) THEN
            BBOZH(MN,1)=(3*BBOZH(MN,2)-BBOZH(MN,3))/2
            RBOZH(MN,1)=(3*RBOZH(MN,2)-RBOZH(MN,3))/2
            ZBOZH(MN,1)=(3*ZBOZH(MN,2)-ZBOZH(MN,3))/2
            PBOZH(MN,1)=(3*PBOZH(MN,2)-PBOZH(MN,3))/2
         ELSE
            BBOZH(MN,1)=0.D0
            RBOZH(MN,1)=0.D0
            ZBOZH(MN,1)=0.D0
            PBOZH(MN,1)=0.D0
         ENDIF
         DEL1=XSH(NSRMAX+1)-XSH(NSRMAX)
         DEL2=XSH(NSRMAX)-XSH(NSRMAX-1)
         BBOZH(MN,NSRMAX+1)=((DEL1+DEL2)*BBOZH(MN,NSRMAX)
     &                       -DEL1*BBOZH(MN,NSRMAX-1))/DEL2
         RBOZH(MN,NSRMAX+1)=((DEL1+DEL2)*RBOZH(MN,NSRMAX)
     &                       -DEL1*RBOZH(MN,NSRMAX-1))/DEL2
         ZBOZH(MN,NSRMAX+1)=((DEL1+DEL2)*ZBOZH(MN,NSRMAX)
     &                       -DEL1*ZBOZH(MN,NSRMAX-1))/DEL2
         PBOZH(MN,NSRMAX+1)=((DEL1+DEL2)*PBOZH(MN,NSRMAX)
     &                       -DEL1*PBOZH(MN,NSRMAX-1))/DEL2
      ENDDO
C
C      DO NSR=1,NSRMAX+1
C         WRITE(6,'(I5,2X,1P3E12.5)') NSR,SHALF(NSR),XS(NSR),XSH(NSR)
C      ENDDO
C
      RETURN
C
 8001 DO NSR=1,NSRMAX
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,BBOZH(MN,NSR)
      ENDDO
      STOP
 8002 DO NSR=1,NSRMAX
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,RBOZH(MN,NSR)
      ENDDO
      STOP
 8003 DO NSR=1,NSRMAX
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,ZBOZH(MN,NSR)
      ENDDO
      STOP
 8004 DO NSR=1,NSRMAX
         WRITE(6,*) 'MN,M,N,NSR,B=',MN,M,N,NSR,PBOZH(MN,NSR)
      ENDDO
      STOP
      END
C
C    ****** CALCULATE R AND Z ******
C
      SUBROUTINE WMHBRZ
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      DIMENSION BSPL(NSRM),RSPL(NSRM),ZSPL(NSRM),PSPL(NSRM)
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
      DO NSR=1,NSRMAX+1
         XSHRHO(NSR)=SQRT(XSH(NSR))
      ENDDO
C
C      ***** PSIPS (normalized) *****
C
      DO NR=1,NRMAX+1
         PSIPS(NR)=XRHO(NR)**2
      ENDDO
C     
C      ***** SPLINE BBOZH(S),RBOZH(S),ZBOZH(S),PBOZH(S) *****
C
      DO MN=1,MNMAX
         DO NSR=1,NSRMAX+1
            BSPL(NSR)=BBOZH(MN,NSR)
            RSPL(NSR)=RBOZH(MN,NSR)
            ZSPL(NSR)=ZBOZH(MN,NSR)
            PSPL(NSR)=PBOZH(MN,NSR)
         ENDDO
         IF(XM(MN).EQ.0.D0) THEN
            FX1(1)=0.D0
            CALL SPL1D(XSHRHO,BSPL,FX1,U1,NSRMAX+1,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: BBOZH'
            FX2(1)=0.D0
            CALL SPL1D(XSHRHO,RSPL,FX2,U2,NSRMAX+1,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: RBOZH'
            FX3(1)=0.D0
            CALL SPL1D(XSHRHO,ZSPL,FX3,U3,NSRMAX+1,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: ZBOZH'
            FX3(1)=0.D0
            CALL SPL1D(XSHRHO,PSPL,FX4,U4,NSRMAX+1,1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: PBOZH'
         ELSE
            CALL SPL1D(XSHRHO,BSPL,FX1,U1,NSRMAX+1,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: BBOZH'
            CALL SPL1D(XSHRHO,RSPL,FX2,U2,NSRMAX+1,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: RBOZH'
            CALL SPL1D(XSHRHO,ZSPL,FX3,U3,NSRMAX+1,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: ZBOZH'
            CALL SPL1D(XSHRHO,PSPL,FX4,U4,NSRMAX+1,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1D: PBOZH'
         ENDIF
C
         DO NR=1,NRMAX+1
            CALL SPL1DD(XRHO(NR),SBMNC(MN,NR),DBMNC(MN,NR),
     &                  XSHRHO,U1,NSRMAX+1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: BMNC: NR=',NR
            CALL SPL1DD(XRHO(NR),SRMNC(MN,NR),DRMNC(MN,NR),
     &                  XSHRHO,U2,NSRMAX+1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: RMNC: NR=',NR
            CALL SPL1DD(XRHO(NR),SZMNS(MN,NR),DZMNS(MN,NR),
     &                  XSHRHO,U3,NSRMAX+1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: ZMNC: NR=',NR
            CALL SPL1DD(XRHO(NR),SPMNS(MN,NR),DPMNS(MN,NR),
     &                  XSHRHO,U4,NSRMAX+1,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX WMHBRZ: SPL1DD: PMNC: NR=',NR
         ENDDO
      ENDDO
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
      DPH=2.D0*PI/(NHC*NPHMAX)
      PSIA=1.D0
C
      DO NR=1,NRMAX+1
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
C
         TH=DTH*(NTH-1)
         PH=DPH*(NPH-1)
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
         BPST(  NTH,NPH,NR)=BPSS
         RPST(  NTH,NPH,NR)=RPSS
         ZPST(  NTH,NPH,NR)=ZPSS
         PPST(  NTH,NPH,NR)=PPSS
C
         IF(NR.NE.1) THEN
            DRS=DRS/(2.D0*PSIA*XRHO(NR))
            DZS=DZS/(2.D0*PSIA*XRHO(NR))
            DPS=DPS/(2.D0*PSIA*XRHO(NR))
            RG11(NTH,NPH,NR)=DRS *DRS +DZS *DZS +DPS *DPS *RPSS**2
            RG12(NTH,NPH,NR)=DRS *DRTH+DZS *DZTH+DPS *DPTH*RPSS**2
            RG13(NTH,NPH,NR)=DRS *DRPH+DZS *DZPH+DPS *DPPH*RPSS**2
            RG22(NTH,NPH,NR)=DRTH*DRTH+DZTH*DZTH+DPTH*DPTH*RPSS**2
            RG23(NTH,NPH,NR)=DRTH*DRPH+DZTH*DZPH+DPTH*DPPH*RPSS**2
            RG33(NTH,NPH,NR)=DRPH*DRPH+DZPH*DZPH+DPPH*DPPH*RPSS**2
            RJ  (NTH,NPH,NR)=(DRS *DZTH-DZS *DRTH)*DPPH*RPSS
     &                      +(DZS *DPTH-DPS *DZTH)*DRPH*RPSS
     &                      +(DPS *DRTH-DRS *DPTH)*DZPH*RPSS
C
            RG11(NTH,NPH,NR)=RG11(NTH,NPH,NR)*XRHO(NR)**2
            RG13(NTH,NPH,NR)=RG13(NTH,NPH,NR)*XRHO(NR)
            RG22(NTH,NPH,NR)=RG22(NTH,NPH,NR)/XRHO(NR)**2
            RG23(NTH,NPH,NR)=RG23(NTH,NPH,NR)/XRHO(NR)
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
      SUBROUTINE WMHBBB
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
C     ***** EXTRAPORATE RIOTAS *****
C
      RIOTAS(1)=(3*RIOTAS(2)-RIOTAS(3))/2
      DEL1=XSH(NSRMAX+1)-XSH(NSRMAX)
      DEL2=XSH(NSRMAX)-XSH(NSRMAX-1)
      RIOTAS(NSRMAX+1)=((DEL1+DEL2)*RIOTAS(NSRMAX)
     &                       -DEL1 *RIOTAS(NSRMAX-1))/DEL2
C
C     ***** INTERPORATE RIOTAS *****
C
      FX1(1)=0.D0
      CALL SPL1D(XSHRHO,RIOTAS,FX1,U1,NSRMAX+1,1,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX WMHBBB: SPL1D: RIOTAS'
C
      DO NR=1,NRMAX+1
         CALL SPL1DF(XRHO(NR),SIOTA(NR),
     &               XSHRHO,U1,NSRMAX+1,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX WMHBBB: SPL1DD: IOTA: NR=',NR
         QPS(NR)=1.D0/SIOTA(NR)
      ENDDO
C
      DO NPH=1,NPHMAX
      DO NTH=1,NTHMAX
      DO NR=1,NRMAX+1
         FACT=  RG22(NTH,NPH,NR)*XRHO(NR)**2*SIOTA(NR)**2
     &       +2*RG23(NTH,NPH,NR)*XRHO(NR)   *SIOTA(NR)
     &       +  RG33(NTH,NPH,NR)
         BFLD(3,NTH,NPH,NR)=BPST(NTH,NPH,NR)/SQRT(FACT)
         BFLD(2,NTH,NPH,NR)=SIOTA(NR)*BFLD(3,NTH,NPH,NR)
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
      DPH=2.D0*PI/(NHC*NPHMAX)
      DO NSU=1,NSUMAX
         DO NPH=1,NPHMAX
            RSU(NSU,NPH)=0.D0
            ZSU(NSU,NPH)=0.D0
            TH=DTHU*(NSU-1)
            PH=DPH*(NPH-1)
            DO MN=1,MNMAX
               RSIN=SIN(-XM(MN)*TH+XN(MN)*PH)
               RCOS=COS(-XM(MN)*TH+XN(MN)*PH)
               RSU(NSU,NPH)=RSU(NSU,NPH)+RBOZH(MN,NSRMAX)*RCOS
               ZSU(NSU,NPH)=ZSU(NSU,NPH)+ZBOZH(MN,NSRMAX)*RSIN
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
               RSIN=SIN(-XM(MN)*TH+XN(MN)*PH)
               RCOS=COS(-XM(MN)*TH+XN(MN)*PH)
               RSW(NSW,NPH)=RSW(NSW,NPH)+SRMNC(MN,NRMAX)*RCOS
               ZSW(NSW,NPH)=ZSW(NSW,NPH)+SZMNS(MN,NRMAX)*RSIN
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
      DO NPH=1,NPHMAX
         DO NSW=1,NSWMAX
            RGMIN=MIN(RGMIN,RSW(NSW,NPH))
            RGMAX=MAX(RGMAX,RSW(NSW,NPH))
            ZGMIN=MIN(ZGMIN,ZSW(NSW,NPH))
            ZGMAX=MAX(ZGMAX,ZSW(NSW,NPH))
         ENDDO
      ENDDO
      RR=0.5D0*(RGMIN+RGMAX)
      BB=BPST(1,1,1)
C
      RETURN
      END
C
C    ****** DRAW BOOZER DATA ******
C
      SUBROUTINE WMGBOOZ
C
      INCLUDE 'wmcomm.inc'
      INCLUDE 'vmcomm.inc'
C
      DIMENSION GX(NSRM),GY(NSRM,NMNM)
C
C     ----- graphics with psit axis -----
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(BBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/BBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(RBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/RBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(ZBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/ZBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(XSH(NSR))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(PBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/PBOZH/',0)
      CALL PAGEE
C
C     ----- graphics with rho axis -----
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(BBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/BBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(RBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/RBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(ZBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/ZBOZH/',0)
      CALL PAGEE
C
      DO NSR=1,NSRMAX+1
         GX(NSR)=GUCLIP(SQRT(XSH(NSR)))
         DO MN=1,MNMAX
            GY(NSR,MN)=GUCLIP(PBOZH(MN,NSR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NSRMAX+1,MNMAX,'/PBOZH/',0)
      CALL PAGEE
C
C     ----- graphics with xrho axis -----
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(SBMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/SBMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(DBMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/DBMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(SRMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/SRMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(DRMNC(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/DRMNC/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(SZMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/SZMNS/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(DZMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/DZMNS/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(SPMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/SPMNS/',0)
      CALL PAGEE
C
      DO NR=1,NRMAX+1
         GX(NR)=GUCLIP(XRHO(NR))
         DO MN=1,MNMAX
            GY(NR,MN)=GUCLIP(DPMNS(MN,NR))
         ENDDO
      ENDDO
      CALL PAGES
      CALL GRF1D(0,GX,GY,NSRM,NRMAX+1,MNMAX,'/DPMNS/',0)
      CALL PAGEE
C
      RETURN
      END
