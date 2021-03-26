! wmdisp.f90

MODULE wmdisp

  PRIVATE
  PUBLIC wm_tnsr

CONTAINS

!     ****** CALCULATE DIELECTRIC TENSOR ******

  SUBROUTINE wm_tnsr(NR,NS,MD,ND)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS,MD,ND
    INTEGER:: NTH,NHH
    
    DO NHH=1,NHHMAX_F
    DO NTH=1,NTHMAX_F
       CTNSR(1,1,NTH,NHH)=0.D0
       CTNSR(1,2,NTH,NHH)=0.D0
       CTNSR(1,3,NTH,NHH)=0.D0
       CTNSR(2,1,NTH,NHH)=0.D0
       CTNSR(2,2,NTH,NHH)=0.D0
       CTNSR(2,3,NTH,NHH)=0.D0
       CTNSR(3,1,NTH,NHH)=0.D0
       CTNSR(3,2,NTH,NHH)=0.D0
       CTNSR(3,3,NTH,NHH)=0.D0
    ENDDO
    ENDDO


    IF((MOD(MODELA,2).EQ.1).AND.(NS.EQ.3)) THEN
       CALL WMTNAX(NR,MD,ND)
    ELSEIF((MOD(MODELA/2,2).EQ.1).AND.(NS.EQ.1)) THEN
       CALL WMTNEX(NR,MD,ND)
    ELSEIF(MODELP(NS).EQ.16) THEN
       CALL WMTNDK(NR,NS,MD,ND)
    ELSE
       IF(MODELP(NS).GE.11) THEN
          CALL WMTNSX(NR,NS,MD,ND)
       ELSE
          CALL WMDPIN(NR,NS,MD,ND)
       ENDIF
    ENDIF

!      IF(NR.EQ.1) THEN
!      WRITE(6,*) 'WMDISP: NR,NS=',NR,NS
!      WRITE(6,'(1P6E12.4)') 
!     &     CTNSR(1,1,0,0,1,1),CTNSR(1,2,0,0,1,1),CTNSR(1,3,0,0,1,1),
!     &     CTNSR(2,1,0,0,1,1),CTNSR(2,2,0,0,1,1),CTNSR(2,3,0,0,1,1),
!     &     CTNSR(3,1,0,0,1,1),CTNSR(3,2,0,0,1,1),CTNSR(3,3,0,0,1,1)
!      ENDIF
!      IF(NR.EQ.2) STOP

    RETURN
  END SUBROUTINE wm_tnsr

!     ****** CALCULATE DIELECTRIC TENSOR: MODELP=11..16 ******

  SUBROUTINE WMTNSX(NR,NS,MD,ND)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    USE wmprof
    USE libdsp,ONLY: DSPFN
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS,MD,ND
    REAL(rkind):: RN(NSMAX),RTPR(NSMAX),RTPP(NSMAX),RU(NSMAX)
    INTEGER:: NHH,NTH,NN,MM,NSS,NC
    REAL(rkind):: WW,BABS,BSUPTH,BSUPPH,RKTH,RKPH,RKPR,RKPP
    REAL(rkind):: RNPR,DTT,DTX,AM,AE,WP,WC,RNPP2,RKPP2,RKX2,RKT2,RKR2
    REAL(rkind):: UXX2,UYY2,RWC,RKTPR,RKTPP,RT,BTH,WTPR,RKPR_EFF
    COMPLEX(rkind):: CW,CWN,CWU,CWP,CWC,CPERP,CPARA,CCROS,CPERM
    COMPLEX(rkind):: CGZ0,CGZ,CADD,CPERP1,CPERP2,CPERM2,CPARA1,CPARA2
    COMPLEX(rkind):: CCROS1,CCROS2,CZ,CDZ

    CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
    WW=DBLE(CW)

      CALL WMCDEN(NR,RN,RTPR,RTPP,RU)

      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F

         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)

         IF(MDLWMK.EQ.1) THEN
            BTH=ABS(BSUPTH*RA*XRHO(NR))
            WTPR=RTPR(NS)*1.D3*AEE/(AMP*PA(NS))
!            RKPR_EFF=SQRT(WW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
            RKPR_EFF=SQRT(WW*BTHOBN(NR)/SQRT(8.D0*WTPR)/RR)
         ELSE
            RKPR_EFF=0.D0
         END IF

            NN=NPH0+NHC*ND
            MM=NTH0+MD

            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH

            IF(ABS(RKPR).LE.RKPR_EFF) RKPR=DSIGN(RKPR_EFF,RKPR)
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            RNPR=VC*RKPR/WW

            IF(MODELP(NS).EQ.14) THEN  ! evaluate cold k_perp
               DTT=1.D0
               DTX=0.D0
               DO NSS=1,NSMAX
                  AM=PA(NSS)*AMP
                  AE=PZ(NSS)*AEE
                  WP=AE*AE*RN(NSS)*1.D20/(AM*EPS0)
                  WC=AE*BABS/AM
                  DTT=DTT-WP/(WW*WW-WC*WC)
                  DTX=DTX-WP*WC/((WW*WW-WC*WC)*WW)
               ENDDO
               RNPP2=((DTT-RNPR**2)**2-DTX**2)/(DTT-RNPR**2)
               RKPP2=RNPP2*WW*WW/(VC*VC)
               RKX2=     MM*MM*RG22(NTH,NHH,NR)*XRHO(NR)**2 &
                   +2.D0*MM*NN*RG23(NTH,NHH,NR)*XRHO(NR) &
                   +     NN*NN*RG33(NTH,NHH,NR)
               IF(RKPP2.GT.0.D0) THEN
                  RKT2=RKX2-RKPR**2
                  IF(RKT2.GT.0.D0) THEN
                     IF(RKT2.LE.RKPP2) THEN
                        RKR2=RKPP2-RKT2
                     ELSE
                        RKR2=0.D0
                     ENDIF
                  ELSE
                     RKT2=0.D0
                     RKR2=RKPP2
                  ENDIF
                  UXX2= RKX2/RKPP2
                  UYY2= RKR2/RKPP2
               ELSE
                  RKPP2=0.D0
                  RKPP=0.D0
                  RKT2=0.D0
                  RKR2=0.D0
                  UXX2=0.D0
                  UYY2=0.D0
               ENDIF
            ELSE
               RKPP2=0.D0
               RKPP=0.D0
               RKT2=0.D0
               RKR2=0.D0
               UXX2=0.D0
               UYY2=0.D0
            ENDIF

            AM=PA(NS)*AMP
            AE=PZ(NS)*AEE

            IF(MODELP(NS).EQ.11) THEN
               CWN=CW-RU(NS)* RKPR+CI*PZCL(NS)*WW
               CWU=CW-RU(NS)* RKPR
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0)
               CWC=AE*BABS/AM
               IF(NS.EQ.1) THEN
                  CPERP=(0.D0,0.D0)
                  CPARA=CI*CWP/(CW*CW*PZCL(NS))
               ELSE
                  CPERP= CWP*CWU*CWN/(CW*CW*CWC**2)
                  CPARA=CI*CWP/(CW*CW*PZCL(NS))
               END IF
               CCROS=(0.D0,0.D0)
               CPERM=(0.D0,0.D0)
            ELSE IF(MODELP(NS).EQ.12) THEN
               CWN=CW-RU(NS)* RKPR+CI*PZCL(NS)*WW
               CWU=CW-RU(NS)* RKPR
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0)
               CWC=AE*BABS/AM
               CPERP=-   CWP*CWU*CWN/(CW*CW*(CWN**2-CWC**2))
               CCROS= CI*CWP*CWC*CWU/(CW*CW*(CWN**2-CWC**2))
               CPARA=-   CWP/(CW*CW) &
                         *(CW*CW/(CWU*CWN) &
                          +2*RKPP2*RU(NS)*RU(NS)/(CWN*CWN-CWC*CWC) &
                           *CWN/CWU)
               CPERM=(0.D0,0.D0)
            ELSE IF(MODELP(NS).EQ.13) THEN
               RT=(RTPR(NS)+2*RTPP(NS))/3.D0
               RKTPR=ABS(RKPR)*SQRT(2.D0*RT*AEE*1.D3/AM)
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0*CW*CW)
               CWC=AE*BABS/AM
               CGZ0=(CW       -RKPR*RU(NS))/RKTPR
               CPERP=0.D0
               CCROS=0.D0
               CPARA=0.D0
               DO NC=-1,1
                  CGZ=(CW-NC*CWC-RKPR*RU(NS))/RKTPR
                  CALL DSPFN(CGZ,CZ,CDZ)
                  IF(NC.EQ.-1) THEN
                     CPERP=CPERP+   CWP*CGZ0*CZ/2
                     CCROS=CCROS-CI*CWP*CGZ0*CZ/2
                  ELSEIF(NC.EQ.0) THEN
!                     CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS))
                     cadd=1+RKPR*RU(NS)*CW/(CW-RKPR*RU(NS))**2 !
!                     CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD*CADD
                     CPARA=CPARA-CWP*CDZ*CGZ*CGZ0*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP=CPERP+   CWP*CGZ0*CZ/2
                     CCROS=CCROS+CI*CWP*CGZ0*CZ/2
                  ENDIF
               ENDDO
               CPERM=(0.D0,0.D0)
            ELSE IF(MODELP(NS).EQ.14) THEN
               CWP=AE*AE*RN(NS)*1.D20/(AM*EPS0*CW*CW)
               RWC=AE*BABS/AM
               RKTPR=ABS(RKPR)*SQRT(2.D0*RTPR(NS)*AEE*1.D3/AM)
               RKTPP=RTPP(NS)*AEE*1.D3/(AM*RWC*RWC)
               CWC=DCMPLX(RWC,0.D0)
               CGZ0=(CW       -RKPR*RU(NS))/RKTPR
               CPERP1= 0.D0
               CPERP2= 0.D0
               CPERM2= 0.D0
               CPARA1= 0.D0
               CPARA2= 0.D0
               CCROS1= 0.D0
               CCROS2= 0.D0
               DO NC=-2,2
                  CGZ =(CW-NC*CWC-RKPR*RU(NS))/RKTPR
                  CADD=1+RKPR*RU(NS)/(CW-RKPR*RU(NS)-NC*CWC)
                  CALL DSPFN(CGZ,CZ,CDZ)
                  IF(NC.EQ.-2) THEN 
                     CPERP2=CPERP2+CWP*RKTPP*CGZ0*CZ
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*CZ
                  ELSEIF(NC.EQ.-1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
                     CPARA2=CPARA2 &
                           -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
                     CCROS1=CCROS1-CI*CWP*CGZ0*CZ/2
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*2*CZ
                  ELSEIF(NC.EQ.0) THEN
                     CPERM2=CPERM2+CWP*RKTPP*CGZ0*4*CZ
                     CPARA1=CPARA1-CWP*CDZ*CGZ*CGZ0*CADD*CADD
                     CPARA2=CPARA2 &
                           +CWP*RKTPP*CGZ0*2*CGZ*CDZ*CADD*CADD
                  ELSEIF(NC.EQ.1) THEN
                     CPERP1=CPERP1+CWP*CGZ0*CZ/2
                     CPERP2=CPERP2-CWP*RKTPP*CGZ0*CZ
                     CPERM2=CPERM2-CWP*RKTPP*CGZ0*2*CZ
                     CPARA2=CPARA2 &
                           -CWP*RKTPP*CGZ0*CGZ*CDZ*CADD*CADD
                     CCROS1=CCROS1+CI*CWP*CGZ0*CZ/2
                     CCROS2=CCROS2-CI*CWP*RKTPP*CGZ0*2*CZ
                  ELSEIF(NC.EQ.2) THEN
                     CPERP2=CPERM2+CWP*RKTPP*CGZ0*CZ
                     CCROS2=CCROS2+CI*CWP*RKTPP*CGZ0*CZ
                  ENDIF
               ENDDO

!               CPERP1= CWP*CGZ(0)*(CZ(1)+CZ(-1))/2
!               CPERP2=-CWP*RKTPP*CGZ(0)
!     &                    *(CZ(1)+CZ(-1)-CZ(2)-CZ(-2))    
!               CPERM2=-CWP*RKTPP*CGZ(0)
!     &                    *(2*CZ(1)+2*CZ(-1)-4*CZ(0))
!               CPARA1=-CWP*CDZ(0)*CGZ(0)*CGZ(0)*CADD(0)*CADD(0)
!               CPARA2= CWP*RKTPP*CGZ(0)
!     &                 *(2*CGZ( 0)*CDZ( 0)*CADD( 0)*CADD( 0)
!     &                    -CGZ( 1)*CDZ( 1)*CADD( 1)*CADD( 1)
!     &                    -CGZ(-1)*CDZ(-1)*CADD(-1)*CADD(-1))
!               CCROS1= CI*CWP*CGZ(0)*(CZ(1)-CZ(-1))/2
!               CCROS2=-CI*CWP*RKTPP*CGZ(0)
!     &                   *(2*CZ(1)-2*CZ(-1)-CZ(2)+CZ(-2))

               CPERP=CPERP1+0.5D0*CPERP2*RKPP2
               CPERM=       0.5D0*CPERM2*RKPP2
               CPARA=CPARA1+0.5D0*CPARA2*RKPP2
               CCROS=CCROS1+0.5D0*CCROS2*RKPP2
            ELSE
               CPERP=(0.D0,0.D0)
               CPARA=(0.D0,0.D0)
               CCROS=(0.D0,0.D0)
               CPERM=(0.D0,0.D0)
            ENDIF

!      IF(NR.EQ.1) THEN
!         WRITE(6,*) 'CPERP,CPERM,CPARA:',CPERP,CPERM,CPARA
!         WRITE(6,*) 'UXX2,UYY2:',UXX2,UYY2
!      ENDIF

!            CTNSR(1,1,NTH,NHH)= CPERP + UXX2*CPERM
!            CTNSR(1,2,NTH,NHH)=-CCROS
!            CTNSR(1,3,NTH,NHH)= 0.D0
!            CTNSR(2,1,NTH,NHH)= CCROS
!            CTNSR(2,2,NTH,NHH)= CPERP + UYY2*CPERM
!            CTNSR(2,3,NTH,NHH)= 0.D0
!            CTNSR(3,1,NTH,NHH)= 0.D0
!            CTNSR(3,2,NTH,NHH)= 0.D0
!            CTNSR(3,3,NTH,NHH)= CPARA

            CTNSR(1,1,NTH,NHH) &
           =CTNSR(1,1,NTH,NHH) + CPERP + UXX2*CPERM
            CTNSR(1,2,NTH,NHH) &
           =CTNSR(1,2,NTH,NHH) - CCROS
            CTNSR(2,1,NTH,NHH) &
           =CTNSR(2,1,NTH,NHH) + CCROS
            CTNSR(2,2,NTH,NHH) &
           =CTNSR(2,2,NTH,NHH) + CPERP + UYY2*CPERM
            CTNSR(3,3,NTH,NHH) &
           =CTNSR(3,3,NTH,NHH) + CPARA

!            IF(NR.EQ.30.AND.NTH.EQ.1.AND.MDX.EQ.1) THEN
!               WRITE(6,*) 'RN,RTPR,RTPP=',RN(1)/1.D20,
!     &                    RTPR(1)/(AEE*1.D3),RTPP(1)/(AEE*1.D3)
!               WRITE(6,*) 'BABS,QR=',BABS,QR
!               WRITE(6,*) 'RKPR,RKPP2=',RKPR,RKPP2
!               WRITE(6,*) 'CGZ(0)=',CGZ(0)
!               WRITE(6,*) 'CDZ(0)=',CDZ(0)
!               WRITE(6,*) 'CTNSR(3,3)=',CTNSR(3,3,NTH,MDX)
!            ENDIF
      ENDDO
      ENDDO

      RETURN
  END SUBROUTINE WMTNSX

!     ****** IMPORT FROM TASK/DP ******

  SUBROUTINE WMDPIN(NR,NS,MD,ND)

!           NR : NODE NUMBER (RADIAL POSITION)
!           NS : PARTICLE SPECIES 

    USE wmcomm
    USE plprofw
    USE plprof
    USE wmprof
    USE dptnsr0
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NS,MD,ND
    TYPE(pl_plfw_type),DIMENSION(NSMAX):: plfw
    TYPE(pl_mag_type):: mag
    TYPE(pl_grd_type),DIMENSION(NSMAX):: grd
    INTEGER:: NHH,NTH,MM,NN,NS1
    COMPLEX(rkind):: CLDISP(6),CDTNS(3,3),CW,CKPR,CKPP
    REAL(rkind):: WW,RHON,RKPR_EFF,RKPR,RNPR,RKPP2
    REAL(rkind):: BABS,BSUPTH,BSUPPH,BTH,WTPR,RKTH,RKPH

    DO ns1=1,nsmax
       grd(ns1)%grdn=0.D0
       grd(ns1)%grdtpr=0.D0
       grd(ns1)%grdtpp=0.D0
       grd(ns1)%grdu=0.D0
    END DO

      CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
      WW=DBLE(CW)

      RHON=XRHO(NR)
      CALL pl_profw(rhon,plfw)

      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         mag%babs=babs
         mag%bnx=0.D0
         mag%bny=0.D0
         mag%bnz=1.D0
         mag%rhon=rhon
         IF(MDLWMK.EQ.1) THEN
            BTH=ABS(BSUPTH*RA*XRHO(NR))
            WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
            RKPR_EFF=SQRT(WW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
         ELSE
            RKPR_EFF=0.D0
         END IF

            NN=NPH0+NHC*ND
            MM=NTH0+MD

            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH

            IF(ABS(RKPR).LE.RKPR_EFF) RKPR=SIGN(RKPR_EFF,RKPR)
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5
            RNPR=VC*RKPR/WW
            
            IF(NR.EQ.1) THEN
               CKPP=(0.D0,0.D0)
            ELSE
!               WRITE(6,'(A,5I4,3ES12.4)') &
!                    'wmdisp:',NR,NS,MM,NN,NTH,XR(NR),RPS(NTH,NR),RKPR
               RKPP2=(MM/XR(NR))**2+(NN/RPST(NTH,NHH,NR))**2-RKPR**2
               IF(RKPP2.GT.0.D0) THEN
                  CKPP=SQRT(DCMPLX(RKPP2,0.D0))
               ELSEIF(RKPP2.LT.0.D0) THEN
                  CKPP=CI*SQRT(DCMPLX(-RKPP2,0.D0))
               ELSE
                  CKPP=(0.D0,0.D0)
               END IF
            END IF
               
            CKPR=RKPR
            CALL dp_tnsr0(CW,CKPR,CKPP,ns,mag,plfw,grd,CLDISP)
            CDTNS(1,1)= CLDISP(1)
            CDTNS(1,2)= CLDISP(5)
            CDTNS(1,3)= CLDISP(4)
            CDTNS(2,1)=-CLDISP(5)
            CDTNS(2,2)= CLDISP(1)+CLDISP(3)
            CDTNS(2,3)= CLDISP(6)
            CDTNS(3,1)= CLDISP(4)
            CDTNS(3,2)=-CLDISP(6)
            CDTNS(3,3)= CLDISP(1)+CLDISP(2)
!            WRITE(6,'(A,3I6)') 'modelp:',modelp(1),modelp(2)
!            WRITE(6,'(A,6E12.4)') 'CW:,',CW,CKPR,CKPP
!            WRITE(6,'(A,5E12.4)') &
!                 'mag:',mag%babs,mag%bnx,mag%bny,mag%bnz,mag%rhon
!            WRITE(6,'(A,3E12.4)') &
!                 'plfw:ns=1:',plfw(1)%rn,plfw(1)%rtpr,plfw(1)%rtpp
!            WRITE(6,'(A,3E12.4)') &
!                 'plfw:ns=2:',plfw(2)%rn,plfw(2)%rtpr,plfw(2)%rtpp
!            WRITE(6,'(A,(6ES12.4))') 'CDTNS:',CDTNS(1:3,1:3)
!            IF(NR.EQ.3) STOP
            CTNSR(1,1,NTH,NHH) &
           =CTNSR(1,1,NTH,NHH) + CDTNS(1,1)
            CTNSR(1,2,NTH,NHH) &
           =CTNSR(1,2,NTH,NHH) + CDTNS(1,2)
            CTNSR(2,1,NTH,NHH) &
           =CTNSR(2,1,NTH,NHH) + CDTNS(2,1)
            CTNSR(2,2,NTH,NHH) &
           =CTNSR(2,2,NTH,NHH) + CDTNS(1,1)
            CTNSR(3,3,NTH,NHH) &
           =CTNSR(3,3,NTH,NHH) + CDTNS(3,3)

      ENDDO
      ENDDO

      RETURN
  END SUBROUTINE WMDPIN

!     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******

  SUBROUTINE WMTNAX(NR,MD,ND)

    USE wmcomm
    USE wmprof
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,MD,ND
    INTEGER:: NHH,NTH,NN,MM,NSEP
    REAL(rkind):: WW,XL,RNA,DRN,RTA,AM,AE,VTA,WP2
    REAL(rkind):: BABS,BSUPTH,BSUPPH,BTH,WTPR,RKPR,RKPR_EFF
    REAL(rkind):: ANGTH,WC,RHOA,RHOR,RKTH,RKPH
    COMPLEX(rkind):: CW,COEF,CWAST,CFN,CX,CEX,CPM,CQM,CRM

      CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
      WW=2.D0*PI*RF*1.D6

      NSEP=3
      CALL WMCPOS(NR,XL)
      IF(XL.LT.RA) THEN
         RNA=PNA*EXP(-(XL/PNAL)**2)*1.D20
         DRN=-2.D0*XL/(PNAL**2)
      ELSE
         RETURN
      ENDIF
      RTA=PTA*AEE*1.D3
      AM=PA(NSEP)*AMP
      AE=PZ(NSEP)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)

      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F

         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)

         IF(MDLWMK.EQ.1) THEN
            BTH=ABS(BSUPTH*RA*XRHO(NR))
            WTPR=RTA*1.D3*AEE/(AMP*PA(NSEP))
            RKPR_EFF=SQRT(WW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
         ELSE
            RKPR_EFF=0.D0
         END IF

         ANGTH=(NTH-1)*2.D0*PI/NTHMAX_F
!         ANGPH=(NHH-1)*2.D0*PI/NHHMAX_F
         WC=AE*BABS/AM
         RHOA=VTA/WC
         RHOR=RHOA/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTA*DRN/(BABS*CW*AE*XL)
         ENDIF

            NN=NPH0+NHC*ND
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.RKPR_EFF) RKPR=RKPR_EFF
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CX=CW/(ABS(RKPR)*VTA)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

            CTNSR(1,1,NTH,NHH) &
           =CTNSR(1,1,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,NTH,NHH) &
           =CTNSR(1,2,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,NTH,NHH) &
           =CTNSR(1,3,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,NTH,NHH) &
           =CTNSR(2,1,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,NTH,NHH) &
           =CTNSR(2,2,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,NTH,NHH) &
           =CTNSR(2,3,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,NTH,NHH) &
           =CTNSR(3,1,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,NTH,NHH) &
           =CTNSR(3,2,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,NTH,NHH) &
           =CTNSR(3,3,NTH,NHH)+CRM
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE WMTNAX

!     ****** CALCULATE ALPHA PARTICLE MAGNETIC DRIFT DIELECTRIC TENSOR ******

  SUBROUTINE WMTNEX(NR,MD,ND)

    USE wmcomm
    USE wmprof
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,MD,ND
    REAL(rkind):: RN(NSMAX),RTPR(NSMAX),RTPP(NSMAX),RU(NSMAX)
    INTEGER:: NSEL,NHH,NTH,NN,MM
    REAL(rkind):: WW,XL,XLP,RNX,RTX,RNXP,XLM,RNXM,RTE,AM,AE,VTE,WP2
    REAL(rkind):: BABS,BSUPTH,BSUPPH,BTH,WTPR,RKPR_EFF
    REAL(rkind):: ANGTH,WC,RHOE,RHOR,RKTH,RKPH,RKPR,RNE,DRN
    COMPLEX(rkind):: CW,COEF,CWAST,CFN,CX,CEX,CPM,CQM,CRM
    
      CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
      WW=2.D0*PI*RF*1.D6

      CALL WMCPOS(NR,XL)
      CALL WMCDEN(NR,RN,RTPR,RTPP,RU)
      NSEL=1
      RNX=RN(NSEL)
      RTX=RTPP(NSEL)
      IF(NR.LT.NRMAX) THEN
         CALL WMCPOS(NR+1,XLP)
         CALL WMCDEN(NR+1,RN,RTPR,RTPP,RU)
         RNXP=RN(NSEL)
      ELSE
         XLP=XL
         RNXP=RNX
      ENDIF
      IF(NR.GT.1) THEN
         CALL WMCPOS(NR-1,XLM)
         CALL WMCDEN(NR-1,RN,RTPR,RTPP,RU)
         RNXM=RN(NSEL)
      ELSE
         XLM=XL
         RNXM=RNX
      ENDIF
      IF(XL.LT.RA) THEN
         RNE= RNX
         DRN=(RNXP-RNXM)/(XLP-XLM)
      ELSE
         RETURN
      ENDIF
      RTE=RTX*AEE*1.D3
      AM=AME
      AE=-AEE
      VTE=SQRT(2.D0*RTE/AM)
      WP2=AE*AE*RNE*1.D20/(AM*EPS0)

      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)

         IF(MDLWMK.EQ.1) THEN
            BTH=ABS(BSUPTH*RA*XRHO(NR))
            WTPR=RTPR(NSEL)*1.D3*AEE/(AMP*PA(NSEL))
            RKPR_EFF=SQRT(WW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
         ELSE
            RKPR_EFF=0.D0
         END IF

         ANGTH=(NTH-1)*2.D0*PI/NTHMAX_F
!         ANGPH=(NHH-1)*2.D0*PI/NHHMAX_F
         WC=AE*BABS/AM
         RHOE=VTE/WC
         RHOR=RHOE/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTE*DRN/(BABS*CW*AE*XL)
         ENDIF

            NN=NPH0+NHC*ND
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.RKPR_EFF) RKPR=RKPR_EFF
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CX=CW/(ABS(RKPR)*VTE)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

            CTNSR(1,1,NTH,NHH) &
           =CTNSR(1,1,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,NTH,NHH) &
           =CTNSR(1,2,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,NTH,NHH) &
           =CTNSR(1,3,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,NTH,NHH) &
           =CTNSR(2,1,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,NTH,NHH) &
           =CTNSR(2,2,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,NTH,NHH) &
           =CTNSR(2,3,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,NTH,NHH) &
           =CTNSR(3,1,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,NTH,NHH) &
           =CTNSR(3,2,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,NTH,NHH) &
           =CTNSR(3,3,NTH,NHH)+CRM
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE WMTNEX

!     ****** CALCULATE MAGNETIC DRIFT DIELECTRIC TENSOR ******

  SUBROUTINE WMTNDK(NR,NS,MD,ND)

    USE wmcomm
    USE plprof
    USE plprofw
    USE wmprof
    IMPLICIT NONE
    TYPE(pl_plfw_type),DIMENSION(NSMAX):: plfw
    INTEGER,INTENT(IN):: NR,NS,MD,ND
    INTEGER:: NHH,NTH,NN,MM
    REAL(rkind):: WW,RHON,RNAP,RTAP,RNA,RTA,DRN
    REAL(rkind):: XL,AM,AE,VTA,WP2,BSUPTH,BSUPPH
    REAL(rkind):: BABS,BTH,WTPR,RKPR_EFF,ANGTH,WC,RHOA,RHOR
    REAL(rkind):: RKTH,RKPH,RKPR
    COMPLEX(rkind):: CW,COEF,CWAST,CFN,CX,CEX,CPM,CQM,CRM

      CW=2.D0*PI*DCMPLX(RF,RFI)*1.D6
      WW=2.D0*PI*RF*1.D6

      RHON=XRHO(NR)
      CALL pl_profw(RHON+0.01D0,plfw)
      RNAP=plfw(NS)%RN*1.D20
      RTAP=(plfw(NS)%RTPR+2.D0*plfw(NS)%RTPP)/3.D0*1.D3*AEE
      CALL pl_profw(RHON,plfw)
      RNA=plfw(NS)%RN*1.D20
      RTA=(plfw(NS)%RTPR+2.D0*plfw(NS)%RTPP)/3.D0*1.D3*AEE
      DRN=(RNAP-RNA)/(0.01D0*RA*RNA)
!      WRITE(6,'(I5,1P3E12.4)') NR,RNA,RTA,DRN

      XL=RHON*RA
      AM=PA(NS)*AMP
      AE=PZ(NS)*AEE
      VTA=SQRT(2.D0*RTA/AM)
      WP2=AE*AE*RNA/(AM*EPS0)

      DO NHH=1,NHHMAX_F
      DO NTH=1,NTHMAX_F
         CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
         IF(MDLWMK.EQ.1) THEN
            BTH=ABS(BSUPTH*RA*XRHO(NR))
            WTPR=plfw(NS)%RTPR*1.D3*AEE/(AMP*PA(NS))
            RKPR_EFF=SQRT(WW*BTH/SQRT(8.D0*WTPR)/RR/BABS)
         ELSE
            RKPR_EFF=0.D0
         END IF

         ANGTH=(NTH-1)*2.D0*PI/NTHMAX_F
!         ANGPH=(NHH-1)*2.D0*PI/NHHMAX_F
         WC=AE*BABS/AM
         RHOA=VTA/WC
         RHOR=RHOA/RR
         COEF=CI*WP2*SQRT(PI)/(CW*CW)
         IF(NR.EQ.1) THEN
            CWAST=0.D0
         ELSE
            CWAST=RTA*DRN/(BABS*CW*AE*XL)
         ENDIF

            NN=NPH0+NHC*ND
            MM=NTH0+MD
            CFN=1.D0-CWAST*MM
            RKTH=MM*BSUPTH/BABS
            RKPH=NN*BSUPPH/BABS
            RKPR=RKTH+RKPH
            IF(ABS(RKPR).LT.RKPR_EFF) RKPR=RKPR_EFF
            IF(ABS(RKPR).LT.1.D-5) RKPR=1.D-5

            CX=CW/(ABS(RKPR)*VTA)
            IF(ABS(CX).GT.5.D0) THEN
               CEX=(0.D0,0.D0)
            ELSE
               CEX=CX*EXP(-CX*CX)
            ENDIF
            CPM=CFN*COEF*RHOR*RHOR*CEX*CX*(1.D0+2.D0*CX*CX+CX*CX*CX*CX)
            CQM=CFN*COEF*RHOR     *CEX*CX*CX*(1.D0+2.D0*CX*CX)
            CRM=CFN*COEF          *CEX*CX*CX*CX*2.D0

            CTNSR(1,1,NTH,NHH) &
           =CTNSR(1,1,NTH,NHH)+CPM*COS(ANGTH)**2
            CTNSR(1,2,NTH,NHH) &
           =CTNSR(1,2,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(1,3,NTH,NHH) &
           =CTNSR(1,3,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(2,1,NTH,NHH) &
           =CTNSR(2,1,NTH,NHH)+CPM*COS(ANGTH)*SIN(ANGTH)
            CTNSR(2,2,NTH,NHH) &
           =CTNSR(2,2,NTH,NHH)+CPM*SIN(ANGTH)**2
            CTNSR(2,3,NTH,NHH) &
           =CTNSR(2,3,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,1,NTH,NHH) &
           =CTNSR(3,1,NTH,NHH)+CQM*COS(ANGTH)
            CTNSR(3,2,NTH,NHH) &
           =CTNSR(3,2,NTH,NHH)+CQM*SIN(ANGTH)
            CTNSR(3,3,NTH,NHH) &
           =CTNSR(3,3,NTH,NHH)+CRM
      ENDDO
      ENDDO
      RETURN
  END SUBROUTINE WMTNDK

END MODULE wmdisp
