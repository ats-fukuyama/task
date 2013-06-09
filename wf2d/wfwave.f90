!     $Id$

!     ********** WF WAVE SOLVER ( A VECTOR, FIRST ORDER ) **********

SUBROUTINE WFWAVE

  use libmtx
  use libmtx
  use wfcomm
  implicit none
  integer :: IERR
  real(4) :: GTMAIN,GTSOLV,GCPUT0,GCPUT1,GCPUT2,GCPUT3

  GTMAIN=0.0
  GTSOLV=0.0
  
  CALL GUTIME(GCPUT0)
  
  if (nrank.eq.0) WRITE(6,*) '--- WFWPRE started ---'
  CALL WFWPRE(IERR)
  IF(IERR.NE.0) GOTO 9000
  
  if (nrank.eq.0) WRITE(6,*) '--- CVCALC started ---'
  CALL CVCALC
  
  CALL GUTIME(GCPUT1)
  
  if(nrank.eq.0) WRITE(6,*) '--- CVSOLV started ---'
  CALL CVSOLV(IERR)
  
  CALL GUTIME(GCPUT2)
  IF(IERR.NE.0) GOTO 9000
  
  CALL CALFLD
  CALL PWRABS
  CALL PWRRAD
  CALL TERMEP
  CALL WFCALB
  CALL GUTIME(GCPUT3)
  GTSOLV=GTSOLV+GCPUT2-GCPUT1
  GTMAIN=GTMAIN+GCPUT3-GCPUT2+GCPUT1-GCPUT0

  if (nrank.eq.0) write (*,*) "GCPUT0,GCPUT1,GCPUT2,GCPUT3,GTSOLV",&
                               GCPUT0,GCPUT1,GCPUT2,GCPUT3,GTSOLV
  if (nrank.eq.0) CALL LPEFLD
  if (nrank.eq.0) WRITE(6,100) GTMAIN,GTSOLV

100 FORMAT(' ','****** CPU TIME : MAIN = ',F10.3,' SEC',5X,&
         &                     ': SOLV = ',F10.3,' SEC ******')
  
9000 continue

  RETURN
END SUBROUTINE WFWAVE

!     ********** WF WAVE PREPARATION **********

SUBROUTINE WFWPRE(IERR)

  use wfcomm
  implicit none
  integer :: IERR

  IERR=0

  CALL LPELMT

  if (nrank.eq.0) WRITE(6,*) '----- SETBDY start ---'
  CALL SETBDY(IERR)
  IF(IERR.NE.0) RETURN
  
  if (nrank.eq.0) WRITE(6,*) '----- SETSID start ---'
  CALL SETSID(IERR)
  IF(IERR.NE.0) RETURN
  
  if (nrank.eq.0) WRITE(6,*) '----- MODANT start ---'
  CALL MODANT(IERR)
  IF(IERR.NE.0) RETURN
  
  if (nrank.eq.0) WRITE(6,*) '----- SETKAN start ---'
  CALL SETKAN(IERR)
  IF(IERR.NE.0) RETURN

  if (nrank.eq.0) WRITE(6,*) '----- DEFBND start ---'
  CALL DEFBND(IERR)
  IF(IERR.NE.0) RETURN
  
  call wffld_allocate
  call wfpwr_allocate

  if (nrank.eq.0) then
     WRITE(6,'(6X,5A10)')  &
          &     '     NNMAX','     NEMAX','    NSDMAX',&
          &     '    NSFMAX','      MLEN'
     WRITE(6,'(A4,2X,6I10)') &
          &     'USED',NNMAX,NEMAX,NSDMAX,NSFMAX,MLEN
  end if
  RETURN
END SUBROUTINE WFWPRE

!     ******* SET KANOD *******

SUBROUTINE SETKAN(IERR)

  use wfcomm
  implicit none
  
  integer :: IERR,NE,NB,NBP,IN,IN1,NN,NSF,ISD,IN2,IN3
  integer :: KAN1,KAN2,KAN3,K1,K2,K3

  IERR=0
  IF(MODELN.EQ.0.OR.MODELN.EQ.4) THEN
     NBMAX=0
     NMMAX=0
     NKMAX=1
     NMKA(1)=0
     DO NE=1,NEMAX
        KAELM(NE)=1
     ENDDO
  ELSEIF(MODELN.EQ.1.OR.MODELN.EQ.5) THEN
     NBMAX=0
  ELSEIF(MODELN.EQ.2.OR.MODELN.EQ.6) THEN
     NMMAX=1
     EPSDM(1)=1.D0
     AMUDM(1)=1.D0
     SIGDM(1)=0.D0
     NKMAX=1
     NMKA(1)=1
     DO NE=1,NEMAX
        KAELM(NE)=1
     ENDDO
  ELSEIF(MODELN.GE.100) THEN
     CALL SETWGB(IERR)
     IF(IERR.NE.0) RETURN
     NMMAX=1
     EPSDM(1)=1.D0
     AMUDM(1)=1.D0
     SIGDM(1)=0.D0
     NKMAX=1
     NMKA(1)=1
     DO NE=1,NEMAX
        KAELM(NE)=1
     ENDDO
  ENDIF
  
  DO NB=1,NBMAX
     if (nrank.eq.0) WRITE(6,'(3I8,4X,A)') NB,KABDY(NB),NBPMAX(NB),KDBDY(NB)
     IF(KABDY(NB).LT.8) THEN
        DO NBP=1,NBPMAX(NB)
           NE=NENBP(NBP,NB)
           DO IN=1,4
              NN=NDELM(IN,NE)
              IF(NN.NE.NDNBP(NBP,NB)) THEN
                 IF(KANOD(NN).EQ.0) THEN
                    if (nrank.eq.0) then
                       WRITE(6,'(A)')     'XX INCONSISTENT KANOD=0 ON BOUNDARY:'
                       WRITE(6,'(A,6I8)') '   NB,NBP,NE,IN,NN,KABDY=',&
                                              NB,NBP,NE,IN,NN,KABDY(NB)
                       WRITE(6,'(A,1P4E12.4)')&
                        '   X,Y,Z,R=',XND(NN),YND(NN),ZND(NN),SQRT(XND(NN)**2+YND(NN)**2)
                    end if
                    IERR=1
                 ELSE
                    IF(KABDY(NB).EQ.4) THEN
                       if (nrank.eq.0) then
                          WRITE(6,'(A,5I8)') 'NB,NBP,NE,NN,KA=',&
                                              NB,NBP,NE,NN,KANOD(NN)
                       end if
                       IF(MODELN.EQ.4  .OR.MODELN.EQ.5.OR.&
                          MODELN.EQ.6  .OR.MODELN.EQ.7.OR.&
                          MODELN.GE.100                   ) THEN
                          KANOD(NN)=-NB
                       ELSE
                          KANOD(NN)=1
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ELSE
        DO NBP=1,NBPMAX(NB)
           NN=NDNBP(NBP,NB)
           IF(KANOD(NN).EQ.0) THEN
              if (nrank.eq.0) then
                 WRITE(6,'(A)') 'XX INCONSISTENT KANOD=0 ON BOUNDARY:'
                 WRITE(6,'(A,6I8)') '   NB,NBP,NE,IN,NN,KABDY=',&
                      &                 NB,NBP,NE,IN,NN,KABDY(NB)
                 WRITE(6,'(A,1P4E12.4)') '   X,Y,Z,R=',&
                      &       XND(NN),YND(NN),ZND(NN),SQRT(XND(NN)**2+YND(NN)**2)
              end if
              IERR=1
           ELSE
              KANOD(NN)=-NB
           ENDIF
        ENDDO
     ENDIF
     
     DO NSF=1,NSFMAX
        KAN1=KANOD(NDSRF(1,NSF))
        KAN2=KANOD(NDSRF(2,NSF))
        KAN3=KANOD(NDSRF(3,NSF))
        IF(KAN1.EQ.-NB.AND.&
          &KAN2.EQ.-NB.AND.&
          &KAN3.EQ.-NB) THEN
           DO ISD=1,3
              KASID(ABS(NSDSRF(ISD,NSF)))=-NB
           ENDDO
        ENDIF
     ENDDO
     
     DO NE=1,NEMAX
        DO IN=1,4
           IN1=MOD(IN,4)+1
           IN2=MOD(IN+1,4)+1
           IN3=MOD(IN+2,4)+1
           K1=KANOD(NDELM(IN1,NE))
           K2=KANOD(NDELM(IN2,NE))
           K3=KANOD(NDELM(IN3,NE))
           IF(K1.EQ.-NB.AND.&
              K2.EQ.-NB.AND.&
              K3.EQ.-NB) THEN
              KNELM(IN,NE)=-NB
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE SETKAN

!     ****** DIELECTRIC TENSOR ******

SUBROUTINE DTENSR(NN,CK)

  use wfcomm
  implicit none

  integer    :: NS,NN
  real(8)    :: WW,WP(NSM),WC(NSM),BABS,AL(3),RN(NSM),RTPR(NSM)
  real(8)    :: RTPP(NSM),RZCL(NSM),FX,FY,FZ
  complex(8) :: CWP,CWC,CDT0,CDX0,CDP0,CDT,CDP,CDX
  complex(8) :: CXX,CXY,CXZ,CYX,CYY,CYZ,CZX,CZY,CZZ,CK(3,3,NSM)

  WW=2.D0*PI*RF*1.D6
  
  DO NS=1,NSMAX
     WP(NS)=PZ(NS)*PZ(NS)*AEE*AEE*1.D20/(PA(NS)*AMP*EPS0*WW*WW)
     WC(NS)=PZ(NS)*AEE/(PA(NS)*AMP*WW)
  ENDDO
  
  CALL WFSMAG(NN,BABS,AL)
  CALL WFSDEN(NN,RN,RTPR,RTPP,RZCL)
  
  DO NS=1,NSMAX
     CWP=WP(NS)*RN(NS)/(1.D0+CII*RZCL(NS))
     CWC=WC(NS)*BABS/(1.D0+CII*RZCL(NS))
     CDT0= CWP/(1.D0-CWC**2)
     CDX0= CII*CWP*CWC/(1.D0-CWC**2)
     CDP0= CWP
     
     CDT=CDT0
     CDP=CDP0-CDT
     CDX=CDX0
     
     FX=AL(1)
     FY=AL(2)
     FZ=AL(3)
     CXX= CDT   +CDP*FX*FX
     CXY= CDX*FZ+CDP*FX*FY
     CXZ=-CDX*FY+CDP*FX*FZ
     CYX=-CDX*FZ+CDP*FY*FX
     CYY= CDT   +CDP*FY*FY
     CYZ= CDX*FX+CDP*FY*FZ
     CZX= CDX*FY+CDP*FZ*FX
     CZY=-CDX*FX+CDP*FZ*FY
     CZZ= CDT   +CDP*FZ*FZ
     
     CK(1,1,NS)=-CXX
     CK(1,2,NS)=-CXY
     CK(1,3,NS)=-CXZ
     CK(2,1,NS)=-CYX
     CK(2,2,NS)=-CYY
     CK(2,3,NS)=-CYZ
     CK(3,1,NS)=-CZX
     CK(3,2,NS)=-CZY
     CK(3,3,NS)=-CZZ
  end do
  
  return
END SUBROUTINE DTENSR

!     ****** CURRENT COEFFICIENT VECTOR CALCULATION ******

SUBROUTINE CVCALC

  use wfcomm
  implicit none
  integer    :: NE,N,NA,IJ,IN,ISD
  real(8)    :: W(4),F(4),DF(3,4),RWE(3,6),DWE(3,4,6),RW,PHASE
  real(8)    :: V,X1,Y1,Z1,X2,Y2,Z2,XM,YM,ZM
  complex(8) :: CJ(3),CVJ

  RW=2.D0*PI*RF*1.D6

  DO NE=1,NEMAX
     DO N=1,7
        CVTOT(N,NE)=0.D0
     ENDDO
  ENDDO
  
  DO NA=1,NAMAX
     PHASE =APH(NA)*PI/180.D0
     CVJ=CII*RW*RMU0*AJ(NA)*EXP(CII*PHASE)
     DO IJ=2,JNUM(NA)
        NE=JELMT(IJ,NA)
        CALL WFABCDX(NE,F,DF,RWE,DWE,V)
        X1=XJ(IJ-1,NA)
        Y1=YJ(IJ-1,NA)
        Z1=ZJ(IJ-1,NA)
        X2=XJ(IJ,NA)
        Y2=YJ(IJ,NA)
        Z2=ZJ(IJ,NA)
        CJ(1)=CVJ*(X2-X1)
        CJ(2)=CVJ*(Y2-Y1)
        CJ(3)=CVJ*(Z2-Z1)
        XM=0.5D0*(X1+X2)
        YM=0.5D0*(Y1+Y2)
        ZM=0.5D0*(Z1+Z2)
        DO IN=1,4
           W(IN)=F(IN)+DF(1,IN)*XM+DF(2,IN)*YM+DF(3,IN)*ZM
        ENDDO
        DO ISD=1,6
           DO IN=1,4
              CVTOT(ISD,NE)=CVTOT(ISD,NE) &
                           +(DWE(1,IN,ISD)*W(IN)*CJ(1)&
                            +DWE(2,IN,ISD)*W(IN)*CJ(2)&
                            +DWE(3,IN,ISD)*W(IN)*CJ(3))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE CVCALC

!     ****** LOCAL ELEMENT MATRIX WITH WG******

SUBROUTINE CMCALC(NE)

  use wfcomm
  implicit none
  integer    :: NME,NE,NBE,ISDMAX,IMDMAX,IM,IN,I,J,K,NK,NS,ID,N,M
  integer    :: K1,K2,K3,L,KN,NB,KA,K0S,K0,ND,K1S,K2S,NSD
  real(8)    :: RWE(3,6),DWE(3,4,6),F(4),DF(3,4),RW,WC,WC2,V,XMU,POS
  real(8)    :: EPSD,SIGD,AMUD,SUM,S
  complex(8) :: CKL(3,3,4),CK(3,3,NSM)
  complex(8) :: CEWG(0:NMDM),CBWG(3,0:NMDM,3),CEPS,CP

  RW=2.D0*PI*RF*1.D6
  WC=RW/VC
  WC2=WC**2
  
  NME=NMKA(KAELM(NE))
  NBE=NBELM(NE)
  
  IF(NBE.EQ.0) THEN
     ISDMAX=6
     IMDMAX=1
  ELSE
     ISDMAX=7
     IMDMAX=NMBDY(NBE)
  ENDIF
  
  CALL WFABCDX(NE,F,DF,RWE,DWE,V)
  
  DO IM=1,ISDMAX
     DO IN=1,ISDMAX
        DO J=1,IMDMAX
           DO I=1,IMDMAX
              CM(I,J,IN,IM)=(0.D0,0.D0)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  
  DO IN=1,ISDMAX
     DO I=1,IMDMAX
        CV(I,IN)=0.D0
     ENDDO
     CV(1,IN)=CVTOT(IN,NE)
  ENDDO

  IF(NME.EQ.0) THEN
     DO K=1,4
        NK=NDELM(K,NE)
        CALL DTENSR(NK,CK)
        DO I=1,3
           DO J=1,3
              IF(I.EQ.J) THEN
                 CKL(I,J,K)=1.D0
              ELSE
                 CKL(I,J,K)=0.D0
              ENDIF
              DO NS=1,NSMAX
                 CKL(I,J,K)=CKL(I,J,K)+CK(I,J,NS)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     XMU=1.D0
  ELSE
     DO K=1,4
        IF(MODELA.EQ.0) THEN
           ID=0
        ELSE
           NK=NDELM(K,NE)
           IF(MODELA.EQ.1) THEN
              POS=XND(NK)
           ELSEIF(MODELA.EQ.2) THEN
              POS=YND(NK)
           ELSEIF(MODELA.EQ.3) THEN
              POS=ZND(NK)
           ELSEIF(MODELA.EQ.4) THEN
              POS=SQRT(XND(NK)**2+YND(NK)**2)
           ELSE
              if (nrank.eq.0) WRITE(6,*) 'XX CMCALC: UNDEFINED MODELA=',MODELA
           ENDIF
           IF((POSRES-POS)*(POS-POSABS).GE.0.D0) THEN
              ID=1
           ELSE
              ID=0
           ENDIF
        ENDIF
        
        IF(ID.EQ.0) THEN
           EPSD=EPSDM(NME)
           SIGD=SIGDM(NME)
           AMUD=AMUDM(NME)
        ELSE
           CEPS=EPSDM(NME)+2.D0*EPSABS*(POS-POSABS)**2 &
                &              /(ABS(POSRES-POSABS)*(ABS(POSRES-POS)-CII*DLTABS))
           EPSD=REAL(CEPS)
           SIGD=REAL(-CII*CEPS)*RW*EPS0
           AMUD=AMUDM(NME)
        ENDIF
        
        DO I=1,3
           DO J=1,3
              IF(I.EQ.J) THEN
                 CKL(I,J,K)=EPSD+CII*SIGD/(RW*EPS0)
              ELSE
                 CKL(I,J,K)=0.D0
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     XMU=1.D0/AMUD
  ENDIF

  DO N=1,6
     DO M=1,6
        CP=RWE(1,N)*RWE(1,M)+RWE(2,N)*RWE(2,M)+RWE(3,N)*RWE(3,M)
        
        CM(1,1,N,M)=CM(1,1,N,M)+XMU*V*CP

        DO I=1,3
           DO J=1,3
              DO K2=1,4
                 SUM=0.D0
                 DO K1=1,4
                    DO K3=1,4
                       SUM=SUM+DWE(I,K1,N)*DWE(J,K3,M)*AIG3(K1,K2,K3)
                    ENDDO
                 ENDDO
                 CM(1,1,N,M)=CM(1,1,N,M)-WC2*V*SUM*CKL(I,J,K2)
              ENDDO
           ENDDO
        ENDDO
        
     ENDDO
  ENDDO
 
!     ***** BOUNDARY SURFACE *****

  DO L=1,4
     KN=KNELM(L,NE)
     IF(KN.LT.0) THEN
        NB=-KN
        KA=KABDY(NB)
        CALL WF2ABC(L,NE,S)
        
        IF(KA.GE.8) THEN
           DO K0S=1,3
              K0=MOD(L+K0S-1,4)+1
              ND=NDELM(K0,NE)
              CALL WFBFWG(ND,CBWG(1,0,K0S))
           ENDDO
         
           DO N=1,6
              DO K1S=1,3
                 K1=MOD(L+K1S-1,4)+1
                 DO K2S=1,3
                    K2=MOD(L+K2S-1,4)+1
                    DO I=1,3
                       CV(1,N)=&
                          CV(1,N)&
                             +CII*WC*XMU*S*DWE(I,K1,N)*CBWG(I,0,K2S)*AIF2(K1S,K2S)
                       DO J=1,IMDMAX
                          CM(1,J,N,7)=&
                             CM(1,J,N,7)&
                               -CII*WC*XMU*S*DWE(I,K1,N)*CBWG(I,J,K2S)*AIF2(K1S,K2S)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           
        ELSEIF(KA.EQ.4) THEN
           DO N=1,6
              DO M=1,6
                 DO K1S=1,3
                    K1=MOD(L+K1S-1,4)+1
                    DO K2S=1,3
                       K2=MOD(L+K2S-1,4)+1
                       DO I=1,3
                          CM(1,1,N,M)=&
                             CM(1,1,N,M)&
                                -CII*WC*XMU*S*DWE(I,K1,N)*DWE(I,K2,M)*AIF2(K1S,K2S)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDIF
  ENDDO
 
!     ----- boundary tangential electric field is given -----
  
  DO M=1,6
     NSD=ABS(NSDELM(M,NE))
     KA=KASID(NSD)
     IF(KA.LT.0) THEN
        NB=-KA
        IF(KABDY(NB).GE.8) THEN
           CALL WFEFWG(NSD,CEWG)
           DO N=1,6
              DO J=1,IMDMAX
                 CM(1,J,N,7)=CM(1,J,N,7) + CM(1,1,N,M)*CEWG(J)
              ENDDO
              CV(1,N)=CV(1,N) - CM(1,1,N,M)*CEWG(0)
           ENDDO
        ENDIF
     ENDIF
  ENDDO

!     ----- boundary weighting function are coupled -----
  
  DO N=1,6
     NSD=ABS(NSDELM(N,NE))
     KA=KASID(NSD)
     IF(KA.LT.0) THEN
        NB=-KA
        IF(KABDY(NB).GE.8) THEN
           CALL WFEFWG(NSD,CEWG)
           DO M=1,6
              DO I=1,IMDMAX
                 CM(I,1,7,M)=CM(I,1,7,M) + DCMPLX(CEWG(I))*CM(1,1,N,M)
              ENDDO
           ENDDO
           M=7
           DO I=1,IMDMAX
              DO J=1,IMDMAX
                 CM(I,J,7,M)=CM(I,J,7,M) + DCMPLX(CEWG(I))*CM(1,J,N,M)
              ENDDO
           ENDDO
           DO I=1,IMDMAX
              CV(I,7)=CV(I,7) + DCMPLX(CEWG(I))*CV(1,N)
           ENDDO
        ENDIF
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE CMCALC

!     ******* ELECTRIC FIELD CALCULATION *******

SUBROUTINE CALFLD

  use wfcomm
  implicit none

  integer    :: NB,NN,L,NSD,KA,I,NE,IN,ISD
  real(8)    :: F(4),DF(3,4),RWE(3,6),DWE(3,4,6),REFL,V
  complex(8) :: CE(3),CEWG(0:NMDM)

  DO NB=1,NBMAX
     IF(KABDY(NB).GE.8) THEN
        NN=NDBDY(NB)
        DO L=1,NMBDY(NB)
           CRFL(L,NB)=CSV(IMLEN(NN)+L-1)
           REFL=ABS(CRFL(L,NB))**2
           if (nrank.eq.0) WRITE(6,'(A,2I5,1P3E12.4)') &
                               'L,NB,CRFL,REFL=',L,NB,CRFL(L,NB),REFL
        ENDDO
     ENDIF
  ENDDO
  
  DO NSD=1,NSDMAX
     KA=KASID(NSD)
     IF(KA.EQ.0) THEN
        CESD(NSD)=CSV(IMLEN(NSD))
     ELSEIF(KA.EQ.1) THEN
        CESD(NSD)=0.D0
     ELSEIF(KA.EQ.4) THEN
        CESD(NSD)=CSV(IMLEN(NSD))
     ELSEIF(KA.LT.0) THEN
        NB=-KA
        CALL WFEFWG(NSD,CEWG)
        CESD(NSD)=CEWG(0)
        DO L=1,NMBDY(NB)
           CESD(NSD)=CESD(NSD)+CEWG(L)*CRFL(L,NB)
        ENDDO
     ELSE
        if (nrank.eq.0) WRITE(6,*) 'XX INVALID KASID: KASID(',NSD,')= ',KASID(NSD)
     ENDIF
  ENDDO

  DO I=1,3
     DO NN=1,NNMAX
        CEF(I,NN)=0.D0
     ENDDO
  ENDDO
  
  DO NE=1,NEMAX
     CALL WFABCDX(NE,F,DF,RWE,DWE,V)
     DO IN=1,4
        CE(1)=(0.D0,0.D0)
        CE(2)=(0.D0,0.D0)
        CE(3)=(0.D0,0.D0)
        DO ISD=1,6
           NSD=ABS(NSDELM(ISD,NE))
           CE(1)=CE(1)+DWE(1,IN,ISD)*CESD(NSD)
           CE(2)=CE(2)+DWE(2,IN,ISD)*CESD(NSD)
           CE(3)=CE(3)+DWE(3,IN,ISD)*CESD(NSD)
        ENDDO
        NN=NDELM(IN,NE)
        CEF(1,NN)=CEF(1,NN)+CE(1)*0.25D0*VELM(NE)/VNOD(NN)
        CEF(2,NN)=CEF(2,NN)+CE(2)*0.25D0*VELM(NE)/VNOD(NN)
        CEF(3,NN)=CEF(3,NN)+CE(3)*0.25D0*VELM(NE)/VNOD(NN)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE CALFLD

!     ****** POWER ABSORPTION ******

SUBROUTINE PWRABS

  use wfcomm
  implicit none

  integer    :: NS,NK,NN,NEDO,NE,NME,NKE,NSMAXL,K,I,J,ID,N,NSD,M,MSD,K1,K2,K3
  real(8)    :: PABSL(NSM),F(4),DF(3,4),RWE(3,6),DWE(3,4,6),XE(4),YE(4),ZE(4)
  real(8)    :: RW,V,POS,EPSD,SIGD,SUM
  complex(8) :: CK(3,3,NSM),CKL(3,3,4,NSM),CNST,CEPS
  
  RW=2.D0*PI*RF*1.D6
  CNST=0.5D0*CII*RW*EPS0
  
  PABST=0.D0
  DO NS=1,NSMAX
     PABSS(NS)=0.D0
  ENDDO
  DO NK=1,NKMAX
     PABSK(NK)=0.D0
  ENDDO
  
  DO NN=1,NNMAX
     PABSTN(NN)=0.D0
     DO NS=1,NSMAX
        PABSSN(NN,NS)=0.D0
     ENDDO
     DO NK=1,NKMAX
        PABSKN(NN,NK)=0.D0
     ENDDO
  ENDDO
  
  DO NEDO=1,NEMAX
     NE=NEDO
     NME=NMKA(KAELM(NE))
     NKE=KAELM(NE)
     
     IF(NME.EQ.0) THEN
        NSMAXL=NSMAX
     ELSE
        NSMAXL=1
     ENDIF
     DO NS=1,NSMAXL
        PABSL(NS)=0.D0
     ENDDO
     
     CALL WFNODE(NE,XE,YE,ZE)
     CALL WFABCDX(NE,F,DF,RWE,DWE,V)
     
     IF(NME.EQ.0) THEN
        DO K=1,4
           NK=NDELM(K,NE)
           CALL DTENSR(NK,CK)
           
           DO I=1,3
              DO J=1,3
                 DO NS=1,NSMAX
                    CKL(I,J,K,NS)=CK(I,J,NS)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ELSE
        DO K=1,4
           IF(MODELA.EQ.0) THEN
              ID=0
           ELSE
              NK=NDELM(K,NE)
              IF(MODELA.EQ.1) THEN
                 POS=XND(NK)
              ELSEIF(MODELA.EQ.2) THEN
                 POS=YND(NK)
              ELSEIF(MODELA.EQ.3) THEN
                 POS=ZND(NK)
              ELSEIF(MODELA.EQ.4) THEN
                 POS=SQRT(XND(NK)**2+YND(NK)**2)
              ELSE
                 if (nrank.eq.0) WRITE(6,*) 'XX PWRABS: UNDEFINED MODELA=',MODELA
              ENDIF
              IF((POSRES-POS)*(POS-POSABS).GE.0.D0) THEN
                 ID=1
              ELSE
                 ID=0
              ENDIF
           ENDIF
           
           IF(ID.EQ.0) THEN
              EPSD=EPSDM(NME)
              SIGD=SIGDM(NME)
           ELSE
              CEPS=EPSDM(NME)+2.D0*EPSABS*(POS-POSABS)**2 &
                   &                 /(ABS(POSRES-POSABS)*(ABS(POSRES-POS)-CII*DLTABS))
              EPSD=REAL(CEPS)
              SIGD=REAL(-CII*CEPS)*RW*EPS0
           ENDIF
           DO I=1,3
              DO J=1,3
                 IF(I.EQ.J) THEN
                    CKL(I,J,K,1)=EPSD+CII*SIGD/(RW*EPS0)
                 ELSE
                    CKL(I,J,K,1)=0.D0
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     
     DO N=1,6
        NSD=ABS(NSDELM(N,NE))
        DO M=1,6
           MSD=ABS(NSDELM(M,NE))
           DO K2=1,4
              DO I=1,3
                 DO J=1,3
                    
                    SUM=0.D0
                    DO K1=1,4
                       DO K3=1,4
                          SUM=SUM+DWE(I,K1,N)*DWE(J,K3,M)*AIG3(K1,K2,K3)
                       ENDDO
                    ENDDO
                    DO NS=1,NSMAXL
                       PABSL(NS)=PABSL(NS) &
                                 -DBLE(CNST*CONJG(CESD(NSD))*SUM &
                                  *CKL(I,J,K2,NS)*CESD(MSD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     
     IF(NME.EQ.0) THEN
        SUM=0.D0
        DO NS=1,NSMAX
           SUM      =SUM      +PABSL(NS)*V
           PABSS(NS)=PABSS(NS)+PABSL(NS)*V
           DO N=1,4
              NN=NDELM(N,NE)
              PABSTN(NN)    =PABSTN(NN)    +PABSL(NS)*V/4.D0
              PABSKN(NN,NKE)=PABSKN(NN,NKE)+PABSL(NS)*V/4.D0
              PABSSN(NN,NS) =PABSSN(NN,NS) +PABSL(NS)*V/4.D0
           ENDDO
        ENDDO
        PABST     =PABST     +SUM
        PABSK(NKE)=PABSK(NKE)+SUM
     ELSE
        PABST     =PABST     +PABSL(1)*V
        PABSK(NKE)=PABSK(NKE)+PABSL(1)*V
        DO N=1,4
           NN=NDELM(N,NE)
           PABSTN(NN)    =PABSTN(NN)    +PABSL(1)*V/4.D0
           PABSKN(NN,NKE)=PABSKN(NN,NKE)+PABSL(1)*V/4.D0
        ENDDO
     ENDIF
  ENDDO
  
  RETURN
END SUBROUTINE PWRABS

!     ****** RADIATED POWER ******

SUBROUTINE PWRRAD

  use wfcomm
  implicit none

  integer    :: NA,IJ,NE,IN,ISD,NSD
  real(8)    :: W(4),F(4),DF(3,4),RWE(3,6),DWE(3,4,6),PHASE,V
  real(8)    :: X1,Y1,Z1,X2,Y2,Z2,XM,YM,ZM
  complex(8) :: CJ(3),CVJ
  
  DO NA=1,NAMAX
     PHASE =APH(NA)*PI/180.D0
     CVJ=AJ(NA)*EXP(CII*PHASE)
     CIMP(NA)=0.D0
     DO IJ=2,JNUM(NA)
        NE=JELMT(IJ,NA)
        CALL WFABCDX(NE,F,DF,RWE,DWE,V)
        X1=XJ(IJ-1,NA)
        Y1=YJ(IJ-1,NA)
        Z1=ZJ(IJ-1,NA)
        X2=XJ(IJ,NA)
        Y2=YJ(IJ,NA)
        Z2=ZJ(IJ,NA)
        CJ(1)=CVJ*(X2-X1)
        CJ(2)=CVJ*(Y2-Y1)
        CJ(3)=CVJ*(Z2-Z1)
        XM=0.5D0*(X1+X2)
        YM=0.5D0*(Y1+Y2)
        ZM=0.5D0*(Z1+Z2)
        DO IN=1,4
           W(IN)=F(IN)+DF(1,IN)*XM+DF(2,IN)*YM+DF(3,IN)*ZM
        ENDDO
        DO ISD=1,6
           NSD=ABS(NSDELM(ISD,NE))
           DO IN=1,4
              CIMP(NA)=CIMP(NA)&
                       -0.5D0*CONJG(CESD(NSD)) &
                           *(DWE(1,IN,ISD)*W(IN)*CJ(1) &
                            +DWE(2,IN,ISD)*W(IN)*CJ(2) &
                            +DWE(3,IN,ISD)*W(IN)*CJ(3))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  CTIMP=0.D0
  DO NA=1,NAMAX
     CTIMP=CTIMP+CIMP(NA)
  ENDDO
  
  RETURN
END SUBROUTINE PWRRAD

!     ****** COMPLETE EFIELD AND POWER ******

SUBROUTINE TERMEP
  
  use wfcomm
  implicit none

  integer :: NN,NS,NK,I,J,MAX
  real(8) :: EABS(3),FACT,FACTR,ETOTL
  
  IF(PIN.GE.0.D0.OR.PABST.EQ.0.D0) THEN
     FACT=1.D0
  ELSE
     FACT=-PIN/ABS(PABST)
  ENDIF
  FACTR=SQRT(FACT)
  
  PABST=FACT*PABST
  DO NN=1,NNMAX
     PABSTN(NN)=FACT*PABSTN(NN)/VNOD(NN)
  ENDDO
  DO NS=1,NSMAX
     PABSS(NS)=FACT*PABSS(NS)
     DO NN=1,NNMAX
        PABSSN(NN,NS)=FACT*PABSSN(NN,NS)/VNOD(NN)
     ENDDO
  ENDDO
  DO NK=1,NKMAX
     PABSK(NK)=FACT*PABSK(NK)
     DO NN=1,NNMAX
        PABSKN(NN,NK)=FACT*PABSKN(NN,NK)/VNOD(NN)
     ENDDO
  ENDDO
  
  PNMAX=PABSTN(1)
  DO NN=2,NNMAX
     PNMAX=MAX(PNMAX,PABSTN(NN))
  ENDDO
  
  ETMAX=0.D0
  DO I=1,3
     EMAX(I)=0.D0
  ENDDO
  
  DO NN=1,NNMAX
     DO J=1,3
        CEF(J,NN)=FACTR*CEF(J,NN)
        EABS(J)=ABS(CEF(J,NN))
        EMAX(J)=MAX(EMAX(J),EABS(J))
     ENDDO

     ETOTL =SQRT(EABS(1)*EABS(1)+EABS(2)*EABS(2)+EABS(3)*EABS(3))
     ETMAX =MAX(ETMAX,ETOTL)
  ENDDO
  
  RETURN
END SUBROUTINE TERMEP
      
!     ****** CALCULATE MAGNETIC FIELD AND POLARIZED FIELD ******

SUBROUTINE WFCALB

  use wfcomm
  implicit none

  integer    :: I,NN,NE,ISD,NSD,IN
  real(8)    :: F(4),DF(3,4),RWE(3,6),DWE(3,4,6),AL(3)
  real(8)    :: RW,V,BABS,SUM,FACTOR
  complex(8) :: CB(3),COEFB,CB1,CB2,CE1,CE2

  RW=2.D0*PI*RF*1.D6
  COEFB=-CII*VC/RW
  
  DO I=1,3
     DO NN=1,NNMAX
        CBF(I,NN)=0.D0
     ENDDO
  ENDDO
  
  DO NE=1,NEMAX
     CALL WFABCDX(NE,F,DF,RWE,DWE,V)
     CB(1)=(0.D0,0.D0)
     CB(2)=(0.D0,0.D0)
     CB(3)=(0.D0,0.D0)
     DO ISD=1,6
        NSD=ABS(NSDELM(ISD,NE))
        CB(1)=CB(1)+COEFB*RWE(1,ISD)*CESD(NSD)
        CB(2)=CB(2)+COEFB*RWE(2,ISD)*CESD(NSD)
        CB(3)=CB(3)+COEFB*RWE(3,ISD)*CESD(NSD)
     ENDDO
     DO IN=1,4
        NN=NDELM(IN,NE)
        CBF(1,NN)=CBF(1,NN)+CB(1)*0.25D0*VELM(NE)/VNOD(NN)
        CBF(2,NN)=CBF(2,NN)+CB(2)*0.25D0*VELM(NE)/VNOD(NN)
        CBF(3,NN)=CBF(3,NN)+CB(3)*0.25D0*VELM(NE)/VNOD(NN)
     ENDDO
  ENDDO
  
  DO NN=1,NNMAX
     CALL WFSMAG(NN,BABS,AL)
     SUM=SQRT(AL(1)*AL(1)+AL(3)*AL(3))
     IF(SUM.EQ.0.D0) THEN
        CE1=CEF(3,NN)
        CE2=CEF(1,NN)
        CB1=CBF(3,NN)
        CB2=CBF(1,NN)
     ELSE
        CE1=(AL(3)*CEF(1,NN)-AL(1)*CEF(3,NN))/SUM
        CB1=(AL(3)*CBF(1,NN)-AL(1)*CBF(3,NN))/SUM
        CE2=SUM*CEF(2,NN) &
             &         -AL(2)*(AL(1)*CEF(1,NN)+AL(3)*CEF(3,NN))/SUM
        CB2=SUM*CBF(2,NN) &
             &         -AL(2)*(AL(1)*CBF(1,NN)+AL(3)*CBF(3,NN))/SUM
     ENDIF
     CEP(1,NN)=CE1-CII*CE2
     CEP(2,NN)=CE1+CII*CE2
     CEP(3,NN)=AL(1)*CEF(1,NN)+AL(2)*CEF(2,NN)+AL(3)*CEF(3,NN)
     CBP(1,NN)=CB1-CII*CB2
     CBP(2,NN)=CB1+CII*CB2
     CBP(3,NN)=AL(1)*CBF(1,NN)+AL(2)*CBF(2,NN)+AL(3)*CBF(3,NN)
  ENDDO
  
  FACTOR=0.5D0/(VC*RMU0)
  DO NN=1,NNMAX
     PFV(NN,1)=FACTOR*DBLE((DCONJG(CEF(2,NN))*CBF(3,NN) &
          &                -DCONJG(CEF(3,NN))*CBF(2,NN)))
     PFV(NN,2)=FACTOR*DBLE((DCONJG(CEF(3,NN))*CBF(1,NN) &
          &                -DCONJG(CEF(1,NN))*CBF(3,NN)))
     PFV(NN,3)=FACTOR*DBLE((DCONJG(CEF(1,NN))*CBF(2,NN) &
          &                -DCONJG(CEF(2,NN))*CBF(1,NN)))
!     WRITE(6,'(A,I8,1P3E12.4)') 'NN,PFV=', &
!          &        NN,PFV(NN,1),PFV(NN,2),PFV(NN,3)
  ENDDO
  
  RETURN
END SUBROUTINE WFCALB

!     ******* OUTPUT FIELD DATA *******

SUBROUTINE LPEFLD

  use wfcomm
  implicit none

  integer :: I,J,NS,NA,NK

  IF(NPRINT.LT.1) RETURN
     
  WRITE(6,110) (EMAX(I),I=1,3),ETMAX,PNMAX
110 FORMAT(' ','EXMAX  =',1PE12.4 &
         &,3X ,'EYMAX  =',1PE12.4 &
         &,3X ,'EZMAX  =',1PE12.4/&
         & ' ','EMAX   =',1PE12.4 &
         &,3X ,'PNMAX  =',1PE12.4)
  
  WRITE(6,120) DBLE(CTIMP),PABST
120 FORMAT(' ','RADIATED POWER =',1PE12.4/&
         & ' ','ABSORBED POWER =',1PE12.4)
  
  IF(ABS(PABST).GT.1.D-32) THEN
     DO NS=1,NSMAX
        WRITE(6,125) NS,PABSS(NS)
125     FORMAT(' ','      PABS(NS=',I2,') =',1PE12.4)
     ENDDO
     DO NK=1,NKMAX
        WRITE(6,126) NK,PABSK(NK)
126     FORMAT(' ','      PABS(NK=',I2,') =',1PE12.4)
     ENDDO
  ENDIF
  
  IF(NAMAX.GT.0) THEN
     WRITE(6,130)
130  FORMAT(' ',' I JNUM', ' AJ(I)','  APH(I)','  AWD(I)', &
          &    ' APOS(I)','  XJ(I)','   YJ(I)',&
          & 8X,'RADIATED POWER')
  ENDIF
  DO NA=1,NAMAX
     WRITE(6,140) NA,JNUM(NA),AJ(NA),APH(NA),AWD(NA),APOS(NA),&
          &                            XJ(1,NA),YJ(1,NA),CIMP(NA)
140  FORMAT(' ',I3,I3,0PF7.4,F7.2,4F7.4,2X,'(',1P2E12.4,')')
  ENDDO
  
  IF(NPRINT.LT.2) RETURN
  
  WRITE(6,150) (I,KANOD(I),PABSTN(I),(CEF(J,I),J=1,3),I=1,NNMAX)
150 FORMAT('ABSORBED POWER DENSITY AND ELECTRIC FIELD'/&
         &       'NOD ','KB',3X,'POWER',12X,'EX',19X,'EY',19X,'EZ'/&
         &      (I4,I2,1PE10.2,3(1X,1P2E10.2)))
  
  RETURN
END SUBROUTINE LPEFLD

!     ******* OUTPUT ELEMENT DATA *******

SUBROUTINE LPELMT

  use wfcomm
  implicit none

  integer :: I,J,NA

  IF(NPRINT.LT.3) RETURN
     
  WRITE(6,110) NNMAX
110 FORMAT(/' ','NODE DATA     : #### NNMAX =',I5,' ####'/&
         &       ' ',2('  NNMAX',' IMLEN',' KANOD',' INLEN',&
         &       9X,'X',14X,'Y',9X))
  WRITE(6,115) (I,IMLEN(I),KANOD(I),INLEN(I),XND(I),YND(I),&
       &              I=1,NNMAX)
115 FORMAT((' ',2(4I6,2X,1P2E15.7,2X)))
  
  WRITE(6,120) NEMAX,(I,(NDELM(J,I),J=1,4),I=1,NEMAX)
120 FORMAT(/' ','ELEMENT DATA  : #### NEMAX =',I5,' ####'/&
         &      (' ',5(I6,'(',4I5,')',2X)))
  
  WRITE(6,125) NEMAX,(I,(ISDELM(J,I),J=1,7),I=1,NEMAX)
125 FORMAT(/' ','NOP     DATA  : #### NEMAX =',I5,' ####'/&
         &      (' ',(I8,'(',7I8,')',2X)))
  
  WRITE(6,130) NBMAX
130 FORMAT(/' ','BOUNDARY DATA : #### NBMAX =',I5,' ####'/&
         &       ' ',3('  NO.',' NDBDY',2X))
  WRITE(6,135) (I,NDBDY(I),I=1,NBMAX)
135 FORMAT((' ',3(2I5,2X)))
  
  DO NA=1,NAMAX
     WRITE(6,140) NA,JNUM0(NA)
140  FORMAT(/' ','ORIGINAL ANTENNA DATA : NA =',I5,' JNUM0 =',I5/&
          &          ' ',3('  NO.',13X,' XJ0',11X,' YJ0',6X))
     WRITE(6,150) (I,XJ0(I,NA),YJ0(I,NA),I=1,JNUM0(NA))
150  FORMAT((' ',3(I5,8X,1P2E15.7)))
     
     WRITE(6,154) NA,JNUM(NA)
154  FORMAT(/' ','MODIFIED ANTENNA DATA : NA =',I5,' JNUM  =',I5/&
          &          ' ',3('  NO.',' JELM',8X,' JX ',11X,' JY ',6X))
     WRITE(6,156) (I,JELMT(I,NA),XJ(I,NA),YJ(I,NA),I=1,JNUM(NA))
156  FORMAT((' ',3(2I5,3X,1P2E15.7)))
  ENDDO
  
  RETURN
END SUBROUTINE LPELMT
