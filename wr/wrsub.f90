!   wrsub.f90

MODULE wrsub

  PRIVATE
  PUBLIC wr_cold_rkperp
  PUBLIC cold_rkr0
  PUBLIC wrmodconv
  PUBLIC wrnwtn
  PUBLIC wrmodnwtn
  PUBLIC dispfn
  PUBLIC dispxr
  PUBLIC dispxi
  PUBLIC wrcale
  PUBLIC wrcalep
  PUBLIC wrcalk
  PUBLIC wrcale_i
  PUBLIC wrcale_xyz

CONTAINS

!  --- calculate wave number with cold plasma model ---
!      Input: RF,RPI,PHII,ZPI

  SUBROUTINE WR_COLD_RKPERP(R,Z,PHI,RKPARA,RKPERP_1,RKPERP_2)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: pl_mag_type,pl_mag
    USE plprofw,ONLY: pl_prfw_type,pl_profw
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: R,Z,PHI,RKPARA
    REAL(rkind),INTENT(OUT):: RKPERP_1,RKPERP_2
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    REAL(rkind):: OMG,OMG_C,WP2,RNPARA
    REAL(rkind):: STIX_X,STIX_Y,STIX_S,STIX_P,STIX_D
    REAL(rkind):: W_A,W_B,W_C,RN_1,RN_2

    OMG=2.D6*PI*RF
    CALL PL_MAG(R*COS(PHI),R*SIN(PHI),Z,mag)
    CALL PL_PROFW(mag%rhon,plfw)

    ! --- consider only electron contribution for ow density ---
    
    OMG_C = mag%BABS*AEE/(AME)
    WP2=plfw(1)%rn*1.D20*AEE*AEE/(EPS0*AMP*PA(1))
    RNPARA=RKPARA*VC/OMG

    STIX_X = WP2/OMG/OMG
    STIX_Y = OMG_C/OMG
    STIX_S = (1.D0-(STIX_X+STIX_Y**2))/(1.D0-STIX_Y**2)
    STIX_P = 1.D0 - STIX_X
    STIX_D = - STIX_X*STIX_Y/(1.D0-STIX_Y**2)

    W_A = STIX_S
    W_B = - STIX_S**2 + STIX_D**2 - STIX_S*STIX_P
    W_B = W_B + (STIX_S + STIX_P)*RNPARA**2
    W_C = STIX_P * ( (STIX_S - RNPARA**2)**2 - STIX_D**2 )


    RN_1 = (-W_B+SQRT(W_B**2-4.D0*W_A*W_C))/2.D0/W_A
    RN_2 = (-W_B-SQRT(W_B**2-4.D0*W_A*W_C))/2.D0/W_A

    WRITE(6,'(A,2ES12.4)') 'WR_COLD_RKPERP:',RN_1,RN_2

    IF(RN_1.GT.0.D0) THEN
       RKPERP_1 = - SQRT(RN_1)*OMG/VC
    ELSE
       RKPERP_1 = 0.D0
    END IF
    IF(RN_2.GT.0.D0) THEN
       RKPERP_2 = - SQRT(RN_2)*OMG/VC
    ELSE
       RKPERP_2 = 0.D0
    END IF

    RETURN
  END SUBROUTINE WR_COLD_RKPERP

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

  SUBROUTINE WRNWTN(RKR_initial,RKR_final,IERR)

    USE wrcomm
    USE dppola
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RKR_initial
    REAL(rkind),INTENT(OUT):: RKR_final
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ICOUNT
    REAL(rkind):: OMG,S,T,RKR,RKR_PRE
    REAL(rkind):: XP,YP,ZP,RKXP,RKXP1,RKXP2,RKYP,RKYP1,RKYP2,RKZP
    COMPLEX(rkind):: cdet(3,3),cepola(3),CRF,CKX,CKY,CKZ
    REAL(rkind):: err

    IERR=0
    OMG=2.D6*PI*RF
    RKR_PRE=RKR_initial

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
       RKXP1= RKR_PRE+DELKR
       RKXP2= RKR_PRE-DELKR
       RKYP= RKPHII
       RKYP1= RKYP
       RKYP2= RKYP
    ELSE
       RKXP=  RKR_PRE       *COS(PHII)-RKPHII*SIN(PHII)
       RKXP1=(RKR_PRE+DELKR)*COS(PHII)-RKPHII*SIN(PHII)
       RKXP2=(RKR_PRE-DELKR)*COS(PHII)-RKPHII*SIN(PHII)
       RKYP=  RKR_PRE       *SIN(PHII)+RKPHII*COS(PHII)
       RKYP1=(RKR_PRE+DELKR)*SIN(PHII)+RKPHII*COS(PHII)
       RKYP2=(RKR_PRE-DELKR)*SIN(PHII)+RKPHII*COS(PHII)
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

      RKR=RKR_PRE-S/T
      IF(MDLWRW.EQ.1) &
           WRITE(6,'(A,1P3E12.4)') 'RKR,RKR_PRE,-S/T=',RKR,RKR_PRE,-S/T

      IF(ABS((RKR-RKR_PRE)/RKR_PRE).LE.EPSNW) GOTO 9000
      IF(ICOUNT.GT.LMAXNW) GOTO 8000
      RKR_PRE=RKR
      GOTO 10

8000  WRITE(6,*) ' WRNWTN: DOES NOT CONVERGE'
      IERR=1000
9000  CONTINUE
      RKR_final=RKR
      CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      CALL dp_pola(crf,ckx,cky,ckz,xp,yp,zp,cdet,cepola,err)
      WRITE(6,'(A,1P2E12.4)') 'CRF=',crf
      WRITE(6,'(A,1P6E12.4)') 'K,X=',rkxp,rkyp,rkzp,xp,yp,zp
      WRITE(6,'(A,1P6E12.4)') 'cep=',cepola(1),cepola(2),cepola(3)
      WRITE(6,'(A,1PE12.4)')  'err=',err
      RETURN   
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
    USE plprof
    USE dpdisp,ONLY: cfdispr
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    TYPE(pl_mag_type):: MAG
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
      CALL pl_mag(X,Y,Z,MAG)
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
            CWC=MAG%BABS*PZ(NS)*AEE/(AMP*PA(NS))
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
    USE plprof
    USE dpdisp,ONLY: cfdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,OMG
    TYPE(pl_mag_type):: MAG
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
      CALL pl_mag(X,Y,Z,MAG)
            
      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

      CWW=2.D0*PI*1.D6*CRF
      DO NS=1,NSMAX
         IF(NSDP(NS).EQ.1) THEN
            CWC=MAG%BABS*PZ(NS)*AEE/(AMP*PA(NS))
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

!  ----- calculate polarization along ray -----

  SUBROUTINE wrcalep(nstp,nray,cepola,cenorm,err)

    USE wrcomm
    USE plprof,ONLY: pl_mag_type,pl_mag
    USE dppola,ONLY: dp_pola
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nstp,nray
    COMPLEX(rkind),INTENT(OUT):: cepola(3),cenorm(3)
    REAL(rkind),INTENT(OUT):: err
    TYPE(pl_mag_type):: mag
    REAL(rkind):: x,y,z,rkx,rky,rkz,rk
    REAL(rkind):: bnx,bny,bnz,knx,kny,knz,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z
    REAL(rkind):: en1,en2,en3
    COMPLEX(rkind):: crf,ckx,cky,ckz,cdet(3,3),cphase

    x=rays(1,nstp,nray)
    y=rays(2,nstp,nray)
    z=rays(3,nstp,nray)
    rkx=rays(4,nstp,nray)
    rky=rays(5,nstp,nray)
    rkz=rays(6,nstp,nray)
    CALL pl_mag(x,y,z,mag)

    crf=rfin(nray)
    ckx=rkx
    cky=rky
    ckz=rkz
    CALL dp_pola(crf,ckx,cky,ckz,x,y,z,cdet,cepola,err)

    bnx=mag%bnx
    bny=mag%bny
    bnz=mag%bnz
    rk=sqrt(rkx**2+rky**2+rkz**2)
    knx=rkx/rk
    kny=rky/rk
    knz=rkz/rk

    r3x=knx              ! r3 = bn
    r3y=kny
    r3z=knz
    r2x=bny*knz-bnz*kny  ! r2 = bn x kn
    r2y=bnz*knx-bnx*knz
    r2z=bnx*kny-bny*knx
    r1x=r2y*r3z-r2z*r3y  ! r1 = r2 x r3
    r1y=r2z*r3x-r2x*r3z
    r1z=r2x*r3y-r2y*r3x

    cenorm(1)=cepola(1)*r1x+cepola(2)*r1y+cepola(3)*r1z  ! O : plane on B and k
    cenorm(2)=cepola(1)*r2x+cepola(2)*r2y+cepola(3)*r2z  ! X : perp to B and k
    cenorm(3)=cepola(1)*r3x+cepola(2)*r3y+cepola(3)*r3z  ! P : parallel to B
    en1=ABS(cenorm(1))
    en2=ABS(cenorm(2))
    en3=ABS(cenorm(3))
    IF(en1.GE.en2) THEN
       IF(en1.GE.en3) THEN
          cphase=CONJG(cenorm(1)/en1)
       ELSE
          cphase=CONJG(cenorm(3)/en3)
       END IF
    ELSE
       IF(en2.GE.en3) THEN
          cphase=CONJG(cenorm(2)/en2)
       ELSE
          cphase=CONJG(cenorm(3)/en3)
       END IF
    END IF
    cenorm(1)=cenorm(1)*cphase
    cenorm(2)=cenorm(2)*cphase
    cenorm(3)=cenorm(3)*cphase
    RETURN
  END SUBROUTINE wrcalep

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

!  ----- calculate eigen electric field for initial position -----

  SUBROUTINE WRCALE_XYZ(X,Y,Z,RKX,RKY,RKZ,EPARA,EPERP)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    USE plcomm_type,ONLY: pl_mag_type
    USE plprof,ONLY: pl_mag
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,Y,Z,RKX,RKY,RKZ
    REAL(rkind),INTENT(OUT):: EPARA,EPERP
    TYPE(pl_mag_type):: mag
    COMPLEX(rkind):: CDET(3,3)
    REAL(rkind):: OMG,EA
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ

    OMG=2.D6*PI*RF
    CRF=DCMPLX(OMG/(2.D6*PI),0.D0)
    CKX=DCMPLX(RKX,0.D0)
    CKY=DCMPLX(RKY,0.D0)
    CKZ=DCMPLX(RKZ,0.D0)
    
    CALL DP_DISP(CRF,CKX,CKY,CKZ,X,Y,Z,CDET)
    
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

    CALL pl_mag(X,Y,Z,mag)
    EPARA=ABS(CUEX*mag%bnx+CUEY*mag%bny+CUEZ*mag%bnz)
    EPERP=SQRT(ABS(CUEX)**2*(1.D0-mag%bnx**2) &
              +ABS(CUEY)**2*(1.D0-mag%bny**2) &
              +ABS(CUEZ)**2*(1.D0-mag%bnz**2))
    RETURN
  END SUBROUTINE WRCALE_XYZ

END MODULE wrsub
