!   wrsub.f90

MODULE wrsub

  PRIVATE
  PUBLIC wr_cold_rkr0,wrmodconv,wrnwtn,wrmodnwtn,dispfn,dispxr,dispxi, &
         wrcale,wrcale_i,wrcalk

CONTAINS

!  --- calculate radial wave number with cold plasma model ---
!      Input: RF,RPI,PHII,ZPI

  SUBROUTINE WR_COLD_RKR0(RKR0_1,RKR0_2)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(OUT):: RKR0_1,RKR0_2
    REAL(rkind):: OMG,PSIN,OMG_C,WP2
    REAL(rkind):: STIX_X,STIX_Y,STIX_S,STIX_P,STIX_D
    REAL(rkind):: W_A,W_B,W_C,RN_1,RN_2


    OMG=2.D6*PI*RF
    CALL PL_MAG_OLD(RPI*COS(PHII),RPI*SIN(PHII),ZPI,PSIN)
    CALL PL_PROF_OLD(PSIN)
    OMG_C = BABS*AEE/(AME)
    WP2=RN(1)*1.D20*AEE*AEE/(EPS0*AMP*PA(1))

!         SWNSNK1 = 1.D0 - WP2/(OMG*OMG-OMG_C*OMG_C)
!         SWNSNK2 = DCMPLX(0.D0, OMG_C*WP2/OMG/(OMG*OMG-OMG_C*OMG_C))
!         SWNSNK3 = 1.D0 - WP2/OMG/OMG
!
!         W_A = SWNSNK1
!         W_B =-((SWNSNK1 - RNPHII**2)*(SWNSNK1+SWNSNK3)+SWNSNK2*SWNSNK2)
!         W_C = SWNSNK3*((SWNSNK1 - RNPHII**2)**2 + SWNSNK2*SWNSNK2)

    STIX_X = WP2/OMG/OMG
    STIX_Y = OMG_C/OMG
    STIX_S = (1.D0-(STIX_X+STIX_Y**2))/(1.D0-STIX_Y**2)
    STIX_P = 1.D0 - STIX_X
    STIX_D = - STIX_X*STIX_Y/(1.D0-STIX_Y**2)

    W_A = STIX_S
    W_B = - STIX_S**2 + STIX_D**2 - STIX_S*STIX_P
    W_B = W_B + (STIX_S + STIX_P)*RNPHII**2
    W_C = STIX_P * ( (STIX_S - RNPHII**2)**2 - STIX_D**2 )


    RN_1 = (-W_B+SQRT(W_B**2-4.D0*W_A*W_C))/2.D0/W_A
    RN_2 = (-W_B-SQRT(W_B**2-4.D0*W_A*W_C))/2.D0/W_A

    IF(RN_1.GT.0.D0) THEN
       RKR0_1 = - SQRT(RN_1)*OMG/VC
    ELSE
       RKR0_1 = 0.D0
    END IF
    IF(RN_2.GT.0.D0) THEN
       RKR0_2 = - SQRT(RN_2)*OMG/VC
    ELSE
       RKR0_2 = 0.D0
    END IF

    RETURN
  END SUBROUTINE WR_COLD_RKR0

!  ----- mode conversion -----

  SUBROUTINE WRMODCONV(IOX, Y, F, OXEFF) 

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IOX
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: F(NEQ)
    REAL(rkind),INTENT(OUT):: OXEFF
    INTEGER:: I
    REAL(rkind):: OMG,OX_K0,OX_LN,OX_Y,OX_NZOPT,OX_NZ,OX_NY,RHON
    REAL(rkind):: BNX0,BNY0,BNZ0,RL0,Rr_IDEI,RKPARA0,S_O_X
    REAL(rkind):: DELTAB,Y10,Y20,Y30,OX_KC,Y4_OX,Y5_OX,Y6_OX
    REAL(rkind):: Y1_OX,Y2_OX,Y3_OX,DELTA

       OMG=2.D6*PI*RF

       CALL RAMBDA_N_OX(Y(1), Y(2), Y(3), OX_LN)
       CALL REFINDEX(Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)
       OX_K0 = OMG / VC
       WRITE(6,*)'N_OPT=',OX_NZOPT,'NZ=',OX_NZ,'NY=',OX_NY
       WRITE(6,*)'K0=',OX_K0,'N=', &
            SQRT((Y(4)**2+Y(5)**2+Y(6)**2))*(VC/OMG)

       OXEFF = ( 2.0*(1.0+OX_Y)*((OX_NZ-OX_NZOPT)**2) + OX_NY**2 )
       OXEFF = EXP(-PI*OX_K0*OX_LN*SQRT(0.5*OX_Y)*OXEFF)
       WRITE(6,*) 'OXEFF=',OXEFF 

       CALL PL_MAG_OLD(Y(1),Y(2),Y(3),RHON)
       CALL PL_PROF_OLD(RHON)

       BNX0 = BNX
       BNY0 = BNY
       BNZ0 = BNZ
       RL0  =SQRT(Y(1)**2+Y(2)**2)
       Rr_IDEI = RL0-RR
       RKPARA0=Y(4)*BNX0+Y(5)*BNY0+Y(6)*BNZ0
       S_O_X = 5.0D-5

       DELTAB =1.0D0
       DO I=1,1000000
          IF(I.EQ.1) THEN 
             DELTAB=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), OMG )
             Y10 = (Rr_IDEI/RL0)*Y(1)
             Y20 = (Rr_IDEI/RL0)*Y(2)
             Y30 = Y(3)
             Y10 = Y10 / SQRT(Y10**2+Y20**2+Y30**2)
             Y20 = Y20 / SQRT(Y10**2+Y20**2+Y30**2)
             Y30 = Y30 / SQRT(Y10**2+Y20**2+Y30**2)
             OX_KC = (Y(4)*Y10 + Y(5)*Y20 + Y(6)*Y30)
             Y4_OX = Y(4) - OX_KC * Y10 
             Y5_OX = Y(5) - OX_KC * Y20 
             Y6_OX = Y(6) - OX_KC * Y30
          END IF
          Y1_OX = Y(1) - IOX * S_O_X * Y10*Y(1)
          Y2_OX = Y(2) - IOX * S_O_X * Y20*Y(2)
          Y3_OX = Y(3) - IOX * S_O_X * Y30*Y(3)

          DELTA=DISPXR( Y1_OX, Y2_OX, Y3_OX, Y4_OX, Y5_OX, Y6_OX, OMG )
		
          IF ( DELTA*DELTAB.LT.0D0 ) THEN
             Y(1) =Y1_OX
             Y(2) =Y2_OX
             Y(3) =Y3_OX
             Y(4) =Y4_OX
             Y(5) =Y5_OX
             Y(6) =Y6_OX
             EXIT
          END IF
          IOX=I
       END DO
					    
       DO I =1,NEQ
          F(I) = Y(I)
       END DO

       RETURN
  END SUBROUTINE WRMODCONV

!  ----- calculate density scalelength -----

  SUBROUTINE  RAMBDA_N_OX(X,Y,Z, OX_LN)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y,Z
    REAL(rkind),INTENT(OUT):: OX_LN
    REAL(rkind):: RHON,OX_NE,RL0,Rr_IDEI
    REAL(rkind):: OX_X0,OX_Y0,OX_Z0,D_OX_X0,D_OX_Y0,D_OX_Z0,OX_NE_P,OX_NE_M

      CALL PL_MAG_OLD(X,Y,Z, RHON)
      CALL PL_PROF_OLD(RHON)
      OX_NE = RN(1) 
	  
      RL0  =SQRT(X**2+Y**2)
      Rr_IDEI = RL0-RR

      OX_X0 = (Rr_IDEI/RL0)*X
      OX_Y0 = (Rr_IDEI/RL0)*Y
      OX_Z0 = Z
      D_OX_X0 = - (DELDER) * OX_X0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
      D_OX_Y0 = - (DELDER) * OX_Y0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
      D_OX_Z0 = - (DELDER) * OX_Z0 / SQRT(OX_X0**2+OX_Y0**2+OX_Z0**2)
	  
!	  WRITE(6,*)'X=',X,'Y=',Y,'Z=',Z
!	  WRITE(6,*)'DX=',D_OX_X0,'DY=',D_OX_Y0,'DZ=',D_OX_Z0
	  
      CALL PL_MAG_OLD(X+D_OX_X0, Y+D_OX_Y0, Z+D_OX_Z0, RHON)
      CALL PL_PROF_OLD(RHON)
      OX_NE_P = RN(1)
      CALL PL_MAG_OLD(X-D_OX_X0, Y-D_OX_Y0, Z-D_OX_Z0, RHON)
      CALL PL_PROF_OLD(RHON)
      OX_NE_M = RN(1) 
	  
!	  WRITE(6,*) 'OX_NE_P(E18)=',OX_NE_P,'OX_NE_M(E18)=', OX_NE_M
	  
      OX_LN = OX_NE / ( (OX_NE_P-OX_NE_M)/(2.0D0*DELDER) )
!	  WRITE(6,*) 'NE(E18)=',OX_NE,'OX_LN=',OX_LN

      RETURN
  END SUBROUTINE RAMBDA_N_OX

!  ----- calculate refractive index after mode conversion -----

  SUBROUTINE  REFINDEX(Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: OX_Y,OX_NZOPT,OX_NZ,OX_NY
    REAL(rkind):: OMG,RHON,OMG_C_OX,RL,Rr_IDEI
    REAL(rkind):: Y10,Y20,Y30,OX_NX,OX_N2

       OMG=2.D6*PI*RF 
		
       CALL PL_MAG_OLD(Y(1), Y(2), Y(3), RHON)
       OMG_C_OX = BABS*AEE/(AME)
       OX_Y = OMG_C_OX / OMG
       OX_NZOPT = SQRT(OX_Y/(1.D0+OX_Y))
       OX_NZ = (Y(4)*BNX + Y(5)*BNY + Y(6)*BNZ)*VC/OMG
	
       RL  =SQRT(Y(1)**2+Y(2)**2)
       Rr_IDEI = RL-RR
       Y10 = (Rr_IDEI/RL)*Y(1)
       Y20 = (Rr_IDEI/RL)*Y(2)
       Y30 = Y(3)
       Y10 = Y10 / SQRT(Y10**2+Y20**2+Y30**2)
       Y20 = Y20 / SQRT(Y10**2+Y20**2+Y30**2)
       Y30 = Y30 / SQRT(Y10**2+Y20**2+Y30**2)
       OX_NX = (Y(4)*Y10 + Y(5)*Y20 + Y(6)*Y30)*VC/OMG
		
       OX_N2 = (Y(4)**2+Y(5)**2+Y(6)**2)*(VC/OMG)*(VC/OMG)
       OX_NY = OX_N2 - OX_NZ**2 - OX_NX**2
       IF (OX_NY.LT.0.D0) OX_NY=0.0D0
       OX_NY = SQRT(OX_NY)
			
       RETURN
  END SUBROUTINE REFINDEX

!  ----- calculate radial wave number satifying D=0 -----

  SUBROUTINE WRNWTN(IERR)

    USE wrcomm
    IMPLICIT NONE
!    REAL(rkind),INTENT(IN):: RKR0
!    REAL(rkind),INTENT(IN):: RKZI,RKPHII
!    REAL(rkind),INTENT(OUT):: RKRI
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ICOUNT
    REAL(rkind):: OMG,S,T,RKR,XP,YP,ZP,RKXP,RKXP1,RKXP2,RKYP,RKYP1,RKYP2,RKZP

      IERR=0
      OMG=2.D6*PI*RF
      RKRI=RKR0

    IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
       XP= RPI
       YP= PHII
       ZP= ZPI
       RKZP= RKZI
    ELSE
       XP= RPI*COS(PHII)
       YP= RPI*SIN(PHII)
       ZP= ZPI
       RKZP= RKZI
    ENDIF

      ICOUNT=0
 10   CONTINUE
      ICOUNT=ICOUNT+1

    IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
       RKXP= RKRI
       RKXP1= RKRI+DELKR
       RKXP2= RKRI-DELKR
       RKYP= RKPHII
       RKYP1= RKYP
       RKYP2= RKYP
    ELSE
       RKXP=  RKRI       *COS(PHII)-RKPHII*SIN(PHII)
       RKXP1=(RKRI+DELKR)*COS(PHII)-RKPHII*SIN(PHII)
       RKXP2=(RKRI-DELKR)*COS(PHII)-RKPHII*SIN(PHII)
       RKYP=  RKRI       *SIN(PHII)+RKPHII*COS(PHII)
       RKYP1=(RKRI+DELKR)*SIN(PHII)+RKPHII*COS(PHII)
       RKYP2=(RKRI-DELKR)*SIN(PHII)+RKPHII*COS(PHII)
    ENDIF

!      WRITE(6,'(1P6E12.4)') XP,YP,ZP,RKXP, RKYP, RKZP
!      WRITE(6,'(1PE12.4)' ) DISPXR(XP,YP,ZP,RKXP, RKYP, RKZP,OMG)
!      WRITE(6,'(1P6E12.4)') XP,YP,ZP,RKXP1, RKYP1,RKZP
!      WRITE(6,'(1PE12.4)' ) DISPXR(XP,YP,ZP,RKXP1,RKYP1,RKZP,OMG)
!      WRITE(6,'(1P6E12.4)') XP,YP,ZP,RKXP2, RKYP2,RKZP
!      WRITE(6,'(1PE12.4)' ) DISPXR(XP,YP,ZP,RKXP2,RKYP2,RKZP,OMG)

      S= DISPXR(XP,YP,ZP,RKXP, RKYP, RKZP,OMG)
      T=(DISPXR(XP,YP,ZP,RKXP1,RKYP1,RKZP,OMG) &
        -DISPXR(XP,YP,ZP,RKXP2,RKYP2,RKZP,OMG))/(2*DELKR)

      RKR=RKRI-S/T
      IF(MDLWRW.EQ.1) &
           WRITE(6,'(A,1P3E12.4)') 'RKR,RKRI,-S/T=',RKR,RKRI,-S/T

      IF(ABS((RKR-RKRI)/RKRI).LE.EPSNW) GOTO 9000
      IF(ICOUNT.GT.LMAXNW) GOTO 8000
      RKRI=RKR
      GOTO 10

 8000 WRITE(6,*) ' WRNWTN: DOES NOT CONVERGE'
      IERR=1000
 9000 RETURN   
  END SUBROUTINE WRNWTN

!  ----- calculate wave number satifying D=0 -----

  SUBROUTINE WRMODNWTN(Y_zh, YK) 

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y_zh(NEQ)
    REAL(rkind),INTENT(OUT):: YK(3)
    REAL(rkind):: Y(NEQ)
    INTEGER:: I,IMODNWTN
    REAL(rkind):: OMG,DELTA,XP,YP,ZP,VV,TT
    REAL(rkind):: RKXP,RKYP,RKZP,RRKXP,RRKYP,RRKZP,DKXP,DKYP,DKZP
    REAL(rkind):: DS2,FAC_NWTN,DDELTA

    Y(1:NEQ)=Y_zh(1:NEQ)
    OMG=2.D6*PI*RF
    DELTA=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), OMG )
    XP=Y(1)
    YP=Y(2)
    ZP=Y(3)
    VV=DELDER
    TT=DELDER

    DO I=1,LMAXNW
       RKXP=Y(4)
       RKYP=Y(5)
       RKZP=Y(6)
       RRKXP=MAX(ABS(RKXP)*VV,TT)
       RRKYP=MAX(ABS(RKYP)*VV,TT)
       RRKZP=MAX(ABS(RKZP)*VV,TT)

       DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,OMG) &
            -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,OMG))/(2.D0*RRKXP)
       DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,OMG) &
            -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,OMG))/(2.D0*RRKYP)
       DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,OMG) &
            -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,OMG))/(2.D0*RRKZP)

       DS2 = DKXP**2+DKYP**2+DKZP**2
       FAC_NWTN = 10.D0
       DO IMODNWTN=1,10
          FAC_NWTN = FAC_NWTN/10.D0
          Y(4) = RKXP - (DELTA*DKXP)/DS2*FAC_NWTN
          Y(5) = RKYP - (DELTA*DKYP)/DS2*FAC_NWTN
          Y(6) = RKZP - (DELTA*DKZP)/DS2*FAC_NWTN
          DDELTA = DISPXR(XP,YP,ZP,Y(4),Y(5),Y(6),OMG)
          IF (ABS(DDELTA) .LT. ABS(DELTA)) EXIT
       END DO
       IF (ABS(DDELTA).LT.1.0D-6) EXIT
       DELTA = DDELTA
    END DO

    DO I =1,3
       YK(I) = Y(I+3)
    END DO

    RETURN
  END SUBROUTINE WRMODNWTN

!  ----- calculate dispersion relation D -----

  FUNCTION DISPFN(RKR,RKPHI,RKZ,RP,ZP,PHI,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp,ONLY: cfdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RKR,RKPHI,RKZ,RP,ZP,PHI,OMG
    REAL(rkind):: DISPFN
    INTEGER:: MODELPS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

    CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
    IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
       CKX=DCMPLX(RKR,0.D0)
       CKY=DCMPLX(RKPHI,0.D0)
       CKZ=DCMPLX(RKZ,0.D0)
       X=RP
       Y=2.D0*PI*RR*SIN(PHI)
       Z=ZP
    ELSEIF(MODELG.EQ.11) THEN
       CKX=DCMPLX(RKR,0.D0)
       CKY=DCMPLX(RKPHI,0.D0)
       CKZ=DCMPLX(RKZ,0.D0)
       X=RP
       Y=PHI
       Z=ZP
    ELSE
       CKX=DCMPLX(RKR*COS(PHI)-RKPHI*SIN(PHI),0.D0)
       CKY=DCMPLX(RKR*SIN(PHI)+RKPHI*COS(PHI),0.D0)
       CKZ=DCMPLX(RKZ,0.D0)
       X=RP*COS(PHI)
       Y=RP*SIN(PHI)
       Z=ZP
    ENDIF
        
    DO NS=1,NSMAX
       MODELPS(NS)=MODELP(NS)
       IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
       IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
       IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
    ENDDO

    CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

    CWW=2.D0*PI*1.D6*CRF
    DO NS=1,NSMAX
       CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
       CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
    ENDDO

    DO NS=1,NSMAX
       MODELP(NS)=MODELPS(NS)
    ENDDO

    DISPFN=DBLE(CF)
    RETURN
  END FUNCTION DISPFN

!  ----- calculate real part of dispersion relation Re D -----

  FUNCTION DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp,ONLY: cfdispr
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    REAL(rkind):: DISPXR
    INTEGER:: MODELPS(NSMAX),MODELVS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
!            
      DO NS=1,NSMAX
         MODELPS(NS)=MODELP(NS)
         MODELVS(NS)=MODELV(NS)
         IF(MODELP(NS).GE.100) MODELV(NS)=0
         IF(MODELP(NS).GE.100.AND.MODELP(NS).LT.300) MODELP(NS)=0
         IF(MODELP(NS).GE.300.AND.MODELP(NS).LT.500) MODELP(NS)=4
         IF(MODELP(NS).GE.500.AND.MODELP(NS).LT.600) MODELP(NS)=6
      ENDDO

!      write(6,*) 'dispxr:',MODELP(1),MODELV(1)

      CF=CFDISPR(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         IF(NSDP(NS).EQ.1) THEN
            CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
         ENDIF
      ENDDO

      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
         MODELV(NS)=MODELVS(NS)
      ENDDO

      DISPXR=DBLE(CF)
      RETURN
  END FUNCTION DISPXR

!  ----- calculate imaginary part of dispersion relation Im D -----

  FUNCTION DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,OMG)

    USE wrcomm
    USE pllocal
    USE dpdisp,ONLY: cfdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    REAL(rkind):: DISPXI
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
            
      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         IF(NSDP(NS).EQ.1) THEN
            CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS))
            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
         ENDIF
      ENDDO

      DISPXI=DIMAG(CF)
      RETURN
  END FUNCTION DISPXI

!  ----- calculate eigen electric field and wave number -----

  SUBROUTINE WRCALE(YN,NSTPMAX_L,NRAY)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(IN):: NSTPMAX_L,NRAY
    COMPLEX(rkind):: CDET(3,3),CDETP(3,3),CDETM(3,3),CDETD(3,3)
    INTEGER:: NSTP,J,I
    REAL(rkind):: OMG,X1,Y1,Z1,VV,TT,ROMG,UE2,UE,RHON,EA
    COMPLEX(rkind):: CRF,CKX1,CKY1,CKZ1,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ,CRFP,CRFM,CUE2

      OMG=2.D6*PI*RF

      DO NSTP=0,NSTPMAX_L
         CRF =DCMPLX(OMG/(2.D6*PI),0.D0)
         X1  =YN(1,NSTP)
         Y1  =YN(2,NSTP)
         Z1  =YN(3,NSTP)
         CKX1=DCMPLX(YN(4,NSTP),0.D0)
         CKY1=DCMPLX(YN(5,NSTP),0.D0)
         CKZ1=DCMPLX(YN(6,NSTP),0.D0)
         CALL DP_DISP(CRF,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDET)
         CE1=CDET(2,2)*CDET(1,3)-CDET(1,2)*CDET(2,3)
         CE2=CDET(1,1)*CDET(2,3)-CDET(2,1)*CDET(1,3)
         CE3=CDET(2,2)*CDET(1,1)-CDET(1,2)*CDET(2,1)
         CE4=CDET(1,3)*CDET(2,1)-CDET(2,3)*CDET(1,1)
         IF(CE2.EQ.0.OR.CE4.EQ.0)THEN
            CEXY=0
            CEZY=0
         ELSE
            CEXY=CE1/CE2
            CEZY=CE3/CE4
         ENDIF

         EA=SQRT(ABS(CEXY)**2+1.D0+ABS(CEZY)**2)

         CUEX=CEXY/EA
         CUEY=1.D0/EA
         CUEZ=CEZY/EA

         VV=DELDER
         TT=DELDER

         ROMG=MAX(ABS(OMG)*VV,TT)
         CRFP =DCMPLX((OMG+ROMG)/(2.D6*PI),0.D0)
         CALL DP_DISP(CRFP,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETP)
         CRFM =DCMPLX((OMG-ROMG)/(2.D6*PI),0.D0)
         CALL DP_DISP(CRFM,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETM)
         
         DO J=1,3
         DO I=1,3
            CDETD(I,J)=OMG*(CDETP(I,J)-CDETM(I,J))/(2.D0*ROMG) &
                      +CDET(I,J)
         ENDDO
         ENDDO
         CUE2=DCONJG(CUEX)*(CDETD(1,1)*CUEX &
                           +CDETD(1,2)*CUEY &
                           +CDETD(1,3)*CUEZ) &
             +DCONJG(CUEY)*(CDETD(2,1)*CUEX &
                           +CDETD(2,2)*CUEY &
                           +CDETD(2,3)*CUEZ) &
             +DCONJG(CUEZ)*(CDETD(3,1)*CUEX &
                           +CDETD(3,2)*CUEY &
                           +CDETD(3,3)*CUEZ)
         UE2=DBLE(CUE2)
         UE=SQRT(ABS(YN(7,NSTP)/UE2))

         CEXS(NSTP,NRAY)=UE*CUEX
         CEYS(NSTP,NRAY)=UE*CUEY
         CEZS(NSTP,NRAY)=UE*CUEZ
         RKXS(NSTP,NRAY)=DBLE(CKX1)
         RKYS(NSTP,NRAY)=DBLE(CKY1)
         RKZS(NSTP,NRAY)=DBLE(CKZ1)
         RXS(NSTP,NRAY)=X1
         RYS(NSTP,NRAY)=Y1
         RZS(NSTP,NRAY)=Z1

         CALL pl_mag_old(X1,Y1,Z1,RHON)
         BNXS(NSTP,NRAY)=BNX
         BNYS(NSTP,NRAY)=BNY
         BNZS(NSTP,NRAY)=BNZ
         BABSS(NSTP,NRAY)=BABS
      ENDDO

      RETURN
  END SUBROUTINE WRCALE

!  ----- calculate wave number along ray -----

  SUBROUTINE WRCALK(NSTP,NRAY,RKPARA,RKPERP)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NSTP,NRAY
    REAL(rkind),INTENT(OUT):: RKPARA,RKPERP
    REAL(rkind):: X,Y,Z,RHON,RKX,RKY,RKZ

      X=RAYS(1,NSTP,NRAY)
      Y=RAYS(2,NSTP,NRAY)
      Z=RAYS(3,NSTP,NRAY)

      CALL PL_MAG_OLD(X,Y,Z,RHON)

      RKX=RAYS(4,NSTP,NRAY)
      RKY=RAYS(5,NSTP,NRAY)
      RKZ=RAYS(6,NSTP,NRAY)

      RKPARA=RKX*BNX+RKY*BNY+RKZ*BNZ
      RKPERP=SQRT(RKX**2+RKY**2+RKZ**2-RKPARA**2)
      RETURN
  END SUBROUTINE WRCALK

!  ----- calculate eigen electric field for initial position -----

  SUBROUTINE WRCALE_I(EPARA,EPERP)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    IMPLICIT NONE
    REAL(rkind),INTENT(OUT):: EPARA,EPERP
    COMPLEX(rkind):: CDET(3,3)
    REAL(rkind):: OMG,X1,Y1,Z1,RHON,EA
    COMPLEX(rkind):: CRF,CKX1,CKY1,CKZ1,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ

    OMG=2.D6*PI*RF
    CRF=DCMPLX(OMG/(2.D6*PI),0.D0)

    IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
       X1= RPI
       Y1= PHII
       Z1= ZPI
       CKX1= RKR0
       CKY1= RKPHII
       CKZ1= RKZI
    ELSE
       X1= RPI*COS(PHII)
       Y1= RPI*SIN(PHII)
       Z1= ZPI
       CKX1= RKR0*COS(PHII)-RKPHII*SIN(PHII)
       CKY1= RKR0*SIN(PHII)+RKPHII*COS(PHII)
       CKZ1= RKZI
    ENDIF
    CALL DP_DISP(CRF,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDET)
    CE1=CDET(2,2)*CDET(1,3)-CDET(1,2)*CDET(2,3)
    CE2=CDET(1,1)*CDET(2,3)-CDET(2,1)*CDET(1,3)
    CE3=CDET(2,2)*CDET(1,1)-CDET(1,2)*CDET(2,1)
    CE4=CDET(1,3)*CDET(2,1)-CDET(2,3)*CDET(1,1)
    IF(CE2.EQ.0.OR.CE4.EQ.0)THEN
       CEXY=0
       CEZY=0
    ELSE
       CEXY=CE1/CE2
       CEZY=CE3/CE4
    ENDIF

    EA=SQRT(ABS(CEXY)**2+1.D0+ABS(CEZY)**2)

    CUEX=CEXY/EA
    CUEY=1.D0/EA
    CUEZ=CEZY/EA

    CALL pl_mag_old(X1,Y1,Z1,RHON)
    EPARA=ABS(CUEX*BNX+CUEY*BNY+CUEZ*BNZ)
    EPERP=SQRT(ABS(CUEX)**2*(1.D0-BNX**2) &
              +ABS(CUEY)**2*(1.D0-BNY**2) &
              +ABS(CUEZ)**2*(1.D0-BNZ**2))
    RETURN
  END SUBROUTINE WRCALE_I

END MODULE wrsub
