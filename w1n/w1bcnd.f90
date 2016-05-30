C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 1) ******
C
C     ****** SET BOUNDARY CONDITIONS AT R=RA ******
C
      SUBROUTINE W1BCND
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1BND1/ CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
C
      DATA CI/(0.D0,1.D0)/
C
      RW=2.D6*PI*RF
      RKV=RW/VC
      RKPR=RKZ
      RNPR=VC*RKPR/RW
      RCE=VC*EPS0
      CKKV2=RNPR*RNPR-1.D0
      CKKV =SQRT(CKKV2)
      CFCTD=EXP(-2.D0*CKKV*RKV*(RD-RA))
      CFCTB=EXP(-2.D0*CKKV*RKV*(RB-RA))
      IF(ABS(WALLR).GT.1.D-12) THEN
         CKKW2=RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR))
         CKKW =SQRT(CKKW2)
         CKKU2=RNPR*RNPR-CKKW2
         CRE=CFCTB*(CKKV-CKKW)/(CKKV+CKKW)
         CRM=CFCTB*(CKKW-CKKV*CKKU2)/(CKKW+CKKV*CKKU2)
      ELSE 
         CRE=-CFCTB
         CRM=-CFCTB
      ENDIF
      CRJ=-CFCTD
      CFCJ=(0.D0,-0.5D0)*EXP(CKKV*RKV*(RD-RA))/(RCE*CKKV)
      CJY1=CFCJ*CFJY1
      CJY2=CFCJ*CFJY2
      CJZ1=CFCJ*CFJZ1
      CJZ2=CFCJ*CFJZ2
C
      DO 110 I=1,3
         DO 110 J=1,5
            CGIN(I,J)=(0.D0,0.D0)
            CGOT(I,J)=(0.D0,0.D0)
  110 CONTINUE
C
      IF(MOD(NMODEL,2).EQ.1) THEN
         CGIN(1,1)=     1.D0+CRE
         CGIN(1,2)=-CI*(1.D0-CRE)*CKKV
         CGIN(2,3)=     1.D0+CRM
         CGIN(2,4)=-CI*(1.D0-CRM)/CKKV
         CGIN(3,1)=    (1.D0+CRJ)     *CJY1
         CGIN(3,2)=-CI*(1.D0-CRJ)*CKKV*CJY1
         CGIN(3,3)=    (1.D0+CRJ)     *CJZ1
         CGIN(3,4)=-CI*(1.D0-CRJ)/CKKV*CJZ1
         CGOT(1,1)=     1.D0+CRE
         CGOT(1,2)= CI*(1.D0-CRE)*CKKV
         CGOT(2,3)=     1.D0+CRM
         CGOT(2,4)= CI*(1.D0-CRM)/CKKV
         CGOT(3,1)=    (1.D0+CRJ)     *CJY2
         CGOT(3,2)= CI*(1.D0-CRJ)*CKKV*CJY2
         CGOT(3,3)=    (1.D0+CRJ)     *CJZ2
         CGOT(3,4)= CI*(1.D0-CRJ)*CKKV*CJZ2
      ELSE
         CGIN(1,1)=  1.D0+CRE
         CGIN(1,2)=-(1.D0-CRE)*CKKV
         CGIN(2,3)=  1.D0+CRM
         CGIN(2,4)=-(1.D0-CRM)/CKKV
         CGIN(3,1)= (1.D0+CRJ)     *CJY1
         CGIN(3,2)=-(1.D0-CRJ)*CKKV*CJY1
         CGIN(3,3)= (1.D0+CRJ)     *CJZ1
         CGIN(3,4)=-(1.D0-CRJ)/CKKV*CJZ1
         CGOT(1,1)=  1.D0+CRE
         CGOT(1,2)=-(1.D0-CRE)*CKKV
         CGOT(2,3)=  1.D0+CRM
         CGOT(2,4)=-(1.D0-CRM)/CKKV
         CGOT(3,1)= (1.D0+CRJ)     *CJY2
         CGOT(3,2)=-(1.D0-CRJ)*CKKV*CJY2
         CGOT(3,3)= (1.D0+CRJ)     *CJZ2
         CGOT(3,4)=-(1.D0-CRJ)*CKKV*CJZ2
      ENDIF
      RETURN
      END
C
C     ******* ELECTROMAGNETIC FIELD IN VACUUM *******
C
      SUBROUTINE W1EVAC(NZ,NSYM)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1BND1/ CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
      COMMON /W1BND2/ CA(NXM)
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1PWR1/ PABSX(NXPM,ISM),FLUXX(NXTM),PABSXZ(ISM),PABSTT,
     &                RANT1(IAM),RANT2(IAM),XANT1(IAM),XANT2(IAM),
     &                PANT1,PANT2,PANT,PWALL1,PWALL2,PWALL,
     &                PIBW1,PIBW2,PIBW,PERR
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
C
      RW=2.D6*PI*RF
      RKV=RW/VC
      RCE=VC*EPS0
C
      IF(NSYM.NE.0.AND.NZ.NE.1.AND.NZ.NE.(NZP/2+1)) THEN
         FACT = 2.0*RZ
      ELSEIF(NSYM.EQ.-1.AND.NZ.EQ.(NZP/2+1)) THEN
         FACT = 0.0
      ELSE
         FACT =     RZ
      ENDIF
      RKPR=RKZ
      RNPR=VC*RKPR/RW
      CKKV2=RNPR*RNPR-1.D0
      CKKV =SQRT(CKKV2)
      CFCTD=EXP(-2.D0*CKKV*RKV*(RD-RA))
      CFCTB=EXP(-2.D0*CKKV*RKV*(RB-RA))
      IF(ABS(WALLR).GT.1.D-12) THEN
         CKKW2=RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR))
         CKKW =SQRT(CKKW2)
         CKKU2=RNPR*RNPR-CKKW2
         CRE=CFCTB*(CKKV-CKKW)/(CKKV+CKKW)
         CRM=CFCTB*(CKKW-CKKV*CKKU2)/(CKKW+CKKV*CKKU2)
      ELSE
         CRE=-CFCTB
         CRM=-CFCTB
      ENDIF
      CRJ=-CFCTD
      CFCJ=(0.D0,-0.5D0)*EXP(CKKV*RKV*(RD-RA))/(RCE*CKKV)
      CJY1=CFCJ*CFJY1
      CJY2=CFCJ*CFJY2
      CJZ1=CFCJ*CFJZ1
      CJZ2=CFCJ*CFJZ2
C
      IF(NMODEL.EQ.1.OR.NMODEL.EQ.3) THEN
         NSF=4*NXP+4
      ELSEIF(NMODEL.EQ.5) THEN
         NSF=6*NXP+4
      ELSE
         NSF=3*NXP+4
      ENDIF
      NXVH=NXV/2
C
      DO 7400 I=1,NXV
C
         NX1=NXP+NXV+I
         NX2=NXP+NXV+1-I
         X=XAM(NX1)
         CEAM=EXP( CKKV*RKV*(X+RA))
         CEAP=EXP(-CKKV*RKV*(X+RA))
C
         IF(I.LE.NXVH) THEN
            DX=RKV*(RB-RD)/DBLE(NXVH-1)
            CAY1=(0.D0,0.D0)
            CAY2=(0.D0,0.D0)
            CAZ1=(0.D0,0.D0)
            CAZ2=(0.D0,0.D0)
         ELSE
            DX=RKV*(RD-RA)/DBLE(NXVH-1)
            CAY1=CJY1
            CAY2=CJY2
            CAZ1=CJZ1
            CAZ2=CJZ2
         ENDIF
C
         CEX1  =CA(2    )*( 0.D0,-1.D0)*RNPR/CKKV*(CEAM-CRM*CEAP)
     &         +CAZ1     *( 0.D0,-1.D0)*RNPR/CKKV*(CEAM-CRJ*CEAP)
         CEY1  =CA(1    )                        *(CEAM+CRE*CEAP)
     &         +CAY1                             *(CEAM+CRJ*CEAP)
         CEZ1  =CA(2    )                        *(CEAM+CRM*CEAP)
     &         +CAZ1                             *(CEAM+CRJ*CEAP)
         CEX1DX=CA(2    )*(-1.D0, 0.D0)*RNPR     *(CEAM+CRM*CEAP)
     &         +CAZ1     *(-1.D0, 0.D0)*RNPR     *(CEAM+CRJ*CEAP)
         CEY1DX=CA(1    )*( 0.D0,-1.D0)*CKKV     *(CEAM-CRE*CEAP)
     &         +CAY1     *( 0.D0,-1.D0)*CKKV     *(CEAM-CRJ*CEAP)
         CEZ1DX=CA(2    )*( 0.D0,-1.D0)*CKKV     *(CEAM-CRM*CEAP)
     &         +CAZ1     *( 0.D0,-1.D0)*CKKV     *(CEAM-CRJ*CEAP)
         CEX2  =CA(NSF  )*( 0.D0, 1.D0)*RNPR/CKKV*(CEAM-CRM*CEAP)
     &         +CAZ2     *( 0.D0, 1.D0)*RNPR/CKKV*(CEAM-CRJ*CEAP)
         CEY2  =CA(NSF-1)                        *(CEAM+CRE*CEAP)
     &         +CAY2                             *(CEAM+CRJ*CEAP)
         CEZ2  =CA(NSF  )                        *(CEAM+CRM*CEAP)
     &         +CAZ2                             *(CEAM+CRJ*CEAP)
         CEX2DX=CA(NSF  )*(-1.D0, 0.D0)*RNPR     *(CEAM+CRM*CEAP)
     &         +CAZ2     *(-1.D0, 0.D0)*RNPR     *(CEAM+CRJ*CEAP)
         CEY2DX=CA(NSF-1)*( 0.D0, 1.D0)*CKKV     *(CEAM-CRE*CEAP)
     &         +CAY2     *( 0.D0, 1.D0)*CKKV     *(CEAM-CRJ*CEAP)
         CEZ2DX=CA(NSF  )*( 0.D0, 1.D0)*CKKV     *(CEAM-CRM*CEAP)
     &         +CAZ2     *( 0.D0, 1.D0)*CKKV     *(CEAM-CRJ*CEAP)
C
         CE2DA(NZ,NX1,1)=CEX1
         CE2DA(NZ,NX1,2)=CEY1
         CE2DA(NZ,NX1,3)=CEZ1
         CJ2DX          = 0.5D0*RNPR*CEZ1
         CJ2DY          =                -CEY1DX
         CJ2DZ          = 0.5D0*RNPR*CEX1-CEZ1DX
         FLUXX(NX1)=FLUXX(NX1) - FACT*( DCONJG(CEX1)*CJ2DX
     &                       + DCONJG(CEY1)*CJ2DY
     &                       + DCONJG(CEZ1)*CJ2DZ ) * RCE
C
         CE2DA(NZ,NX2,1)=CEX2
         CE2DA(NZ,NX2,2)=CEY2
         CE2DA(NZ,NX2,3)=CEZ2
         CJ2DX          = 0.5D0*RNPR*CEZ2
         CJ2DY          =                -CEY2DX
         CJ2DZ          = 0.5D0*RNPR*CEX2-CEZ2DX
         FLUXX(NX2)=FLUXX(NX2) - FACT*( DCONJG(CEX2)*CJ2DX
     &                       + DCONJG(CEY2)*CJ2DY
     &                       + DCONJG(CEZ2)*CJ2DZ ) * RCE
 7400 CONTINUE
C
      RETURN
      END
