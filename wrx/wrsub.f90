
!   wrsub.f90

MODULE wrsub

  PRIVATE
  PUBLIC wr_cold_rkperp
  PUBLIC wr_newton
  PUBLIC cold_rkr0
  PUBLIC wrmodconv
  PUBLIC wrmodnwtn
  PUBLIC dispxr
  PUBLIC dispxi
  PUBLIC wrcale
  PUBLIC wrcalep
  PUBLIC wrcalk
  PUBLIC wrcale_i
  PUBLIC wrcale_xyz

CONTAINS

  ! *** calculate perpendicular wave numbers for parallel wave number ***
  !             (rkperp_1)**2 >= (rkperp_2)**2
  !              rkperp = 0 for evanescent mode
  
  SUBROUTINE WR_cold_rkperp(omega,R,Z,phi,rkpara,rkperp_1,rkperp_2)

    USE plcomm,ONLY: nsmax
    USE bpsd_constants
    USE pllocal
    USE plprof,ONLY: pl_mag_type,pl_mag
    USE plprofw,ONLY: pl_prfw_type,pl_profw
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: omega,R,Z,phi,rkpara
    REAL(rkind),INTENT(OUT):: rkperp_1,rkperp_2
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    REAL(rkind):: omega_c,omega_p2,rnpara
    REAL(rkind):: STIX_X,STIX_Y,STIX_S,STIX_P,STIX_D
    REAL(rkind):: W_A,W_B,W_C,RN_1,RN_2,RN_temp

    CALL pl_mag(R*cos(phi),R*SIN(phi),Z,mag)
    CALL PL_PROFW(mag%rhon,plfw)

    ! --- consider only electron contribution for ow density ---
    
    omega_c = mag%BABS*AEE/AME
    omega_p2=plfw(1)%rn*1.D20*AEE*AEE/(EPS0*AME)
    rnpara=rkpara*VC/omega

    STIX_X = omega_p2/omega**2
    STIX_Y = omega_c/omega
    STIX_S = (1.D0-(STIX_X+STIX_Y**2))/(1.D0-STIX_Y**2)
    STIX_P = 1.D0 - STIX_X
    STIX_D = - STIX_X*STIX_Y/(1.D0-STIX_Y**2)

    W_A = STIX_S
    W_B = - STIX_S**2 + STIX_D**2 - STIX_S*STIX_P
    W_B = W_B + (STIX_S + STIX_P)*RNPARA**2
    W_C = STIX_P * ( (STIX_S - RNPARA**2)**2 - STIX_D**2 )


    RN_1 = (-W_B+SQRT(W_B**2-4.D0*W_A*W_C))/2.D0/W_A
    RN_2 = (-W_B-SQRT(W_B**2-4.D0*W_A*W_C))/2.D0/W_A

!    WRITE(6,'(A,2ES12.4)') 'WR_COLD_RKPERP:',RN_1,RN_2

    IF(RN_2.GT.RN_1) THEN
       RN_temp=RN_1
       RN_1=RN_2
       RN_2=RN_temp
    END IF

    IF(RN_1.GT.0.D0) THEN
       RKPERP_1 = SQRT(RN_1)*omega/VC
    ELSE
       RKPERP_1 = 0.D0
    END IF
    IF(RN_2.GT.0.D0) THEN
       RKPERP_2 = SQRT(RN_2)*omega/VC
    ELSE
       RKPERP_2 = 0.D0
    END IF

    RETURN
  END SUBROUTINE WR_COLD_RKPERP

!  ----- calculate radial wave number satifying D=0 -----

  SUBROUTINE wr_newton(omega,R,PHI,Z,angp,angt,modew,rk_initial,rk_final,IERR)

    USE wrcomm
    USE dppola
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: omega,R,PHI,Z,angp,angt,rk_initial
    INTEGER,INTENT(IN):: modew
    REAL(rkind),INTENT(OUT):: rk_final
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ICOUNT
    REAL(rkind):: S,T,rk,rk_new,delrk
    REAL(rkind):: XP,YP,ZP,rk0,rk1,rk2,rkr0,rkr1,rkr2,rkz0,rkz1,rkz2
    REAL(rkind):: rkx0,rkx1,rkx2,rky0,rky1,rky2,rkph0,rkph1,rkph2
    COMPLEX(rkind):: cdet(3,3),cepola(3),CRF,CKX,CKY,CKZ
    REAL(rkind):: err,deg

    deg=PI/180.D0
    IERR=0

    XP=R*COS(PHI)
    YP=R*SIN(PHI)
    ZP=Z
    rk=rk_initial
    delrk=DELKR
    
    ICOUNT=0
10  CONTINUE
    ICOUNT=ICOUNT+1

    rk0=rk
    rk1=rk+delrk
    rk2=rk-delrk

    SELECT CASE(mdlwri)
    CASE(0,2,100,102)
       rkr0= -rk0*COS(angp*deg)*COS(angt*deg)
       rkph0= rk0*COS(angp*deg)*SIN(angt*deg)
       rkz0=  rk0*SIN(angp*deg)
       rkr1= -rk1*COS(angp*deg)*COS(angt*deg)
       rkph1= rk1*COS(angp*deg)*SIN(angt*deg)
       rkz1=  rk1*SIN(angp*deg)
       rkr2= -rk2*COS(angp*deg)*COS(angt*deg)
       rkph2= rk2*COS(angp*deg)*SIN(angt*deg)
       rkz2=  rk2*SIN(angp*deg)
    CASE(1,3,101,103)
       rkr0= -rk0*COS(angp*deg)*COS(angt*deg)
       rkph0= rk0              *SIN(angt*deg)
       rkz0=  rk0*SIN(angp*deg)*COS(angt*deg)
       rkr1= -rk1*COS(angp*deg)*COS(angt*deg)
       rkph1= rk1              *SIN(angt*deg)
       rkz1=  rk1*SIN(angp*deg)*COS(angt*deg)
       rkr2= -rk2*COS(angp*deg)*COS(angt*deg)
       rkph2= rk2              *SIN(angt*deg)
       rkz2=  rk2*SIN(angp*deg)*COS(angt*deg)
    END SELECT
    SELECT CASE(modew)
    CASE(2,3)
       rkr0 =-rkr0
       rkph0=-rkph0
       rkz0 =-rkz0
       rkr1 =-rkr1
       rkph1=-rkph1
       rkz1 =-rkz1
       rkr2 =-rkr2
       rkph2=-rkph2
       rkz2 =-rkz2
    END SELECT

    rkx0=rkr0*COS(PHI)-rkph0*SIN(PHI)
    rky0=rkr0*SIN(PHI)+rkph0*COS(PHI)
    rkx1=rkr1*COS(PHI)-rkph1*SIN(PHI)
    rky1=rkr1*SIN(PHI)+rkph1*COS(PHI)
    rkx2=rkr2*COS(PHI)-rkph2*SIN(PHI)
    rky2=rkr2*SIN(PHI)+rkph2*COS(PHI)

    S= DISPXR(XP,YP,ZP,rkx0,rky0,rkz0,omega)
    T=(DISPXR(XP,YP,ZP,rkx1,rky1,rkz1,omega) &
      -DISPXR(XP,YP,ZP,rkx2,rky2,rkz2,omega))/(2*delrk)

    rk_new=rk-S/T
    IF(idebug_wr(6).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(6):'
       WRITE(6,'(A,3ES12.4)') '      rk,rk_new,-S/T=',rk,rk_new,-S/T
    END IF

    IF(ABS((rk_new-rk)/rk).LE.EPSNW) GOTO 9000
    IF(ICOUNT.GT.LMAXNW) GOTO 8000
    rk=rk_new
    GOTO 10

8000 WRITE(6,*) 'XX wr_newton: DOES NOT CONVERGE'
    IERR=1000
9000 CONTINUE
    rk_final=rk_new
    crf=DCMPLX(omega/(2.D6*PI),0.D0)
    ckx=rkx0
    cky=rky0
    ckz=rkz0
    CALL dp_pola(crf,ckx,cky,ckz,xp,yp,zp,cdet,cepola,err)
    IF(idebug_wr(7).NE.0) THEN
       WRITE(6,'(A,2I4)') '*** idebug_wr(7):'
       WRITE(6,'(A,3ES12.4)') '   CRF=',crf,rk0
       WRITE(6,'(A,6ES12.4)') '   K,X=',rkx0,rky0,rkz0,xp,yp,zp
       WRITE(6,'(A,6ES12.4)') '   cep=',cepola(1),cepola(2),cepola(3)
       WRITE(6,'(A,ES12.4)')  '   err=',err
    END IF
    RETURN   
  END SUBROUTINE wr_newton

!  ----- mode conversion -----

  SUBROUTINE WRMODCONV(IOX, Y, F, OXEFF,omega) 

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IOX
    REAL(rkind),INTENT(IN):: omega
    REAL(rkind),INTENT(INOUT):: Y(NEQ)
    REAL(rkind),INTENT(OUT):: F(NEQ)
    REAL(rkind),INTENT(OUT):: OXEFF
    INTEGER:: I
    REAL(rkind):: OX_K0,OX_LN,OX_Y,OX_NZOPT,OX_NZ,OX_NY,RHON
    REAL(rkind):: BNX0,BNY0,BNZ0,RL0,Rr_IDEI,RKPARA0,S_O_X
    REAL(rkind):: DELTAB,Y10,Y20,Y30,OX_KC,Y4_OX,Y5_OX,Y6_OX
    REAL(rkind):: Y1_OX,Y2_OX,Y3_OX,DELTA

       CALL RAMBDA_N_OX(Y(1), Y(2), Y(3), OX_LN)
       CALL REFINDEX(omega, Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)
       OX_K0 = omega / VC
       WRITE(6,*)'N_OPT=',OX_NZOPT,'NZ=',OX_NZ,'NY=',OX_NY
       WRITE(6,*)'K0=',OX_K0,'N=', &
            SQRT((Y(4)**2+Y(5)**2+Y(6)**2))*(VC/omega)

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
             DELTAB=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), omega )
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

          DELTA=DISPXR( Y1_OX, Y2_OX, Y3_OX, Y4_OX, Y5_OX, Y6_OX, omega )
		
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

  SUBROUTINE  REFINDEX(omega,Y,OX_Y,OX_NZOPT,OX_NZ, OX_NY)

    USE wrcomm
    USE pllocal
    USE plprof,ONLY: PL_MAG_OLD
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: omega,Y(NEQ)
    REAL(rkind),INTENT(OUT):: OX_Y,OX_NZOPT,OX_NZ,OX_NY
    REAL(rkind):: RHON,omega_C_OX,RL,Rr_IDEI
    REAL(rkind):: Y10,Y20,Y30,OX_NX,OX_N2

       CALL PL_MAG_OLD(Y(1), Y(2), Y(3), RHON)
       omega_C_OX = BABS*AEE/(AME)
       OX_Y = omega_C_OX / omega
       OX_NZOPT = SQRT(OX_Y/(1.D0+OX_Y))
       OX_NZ = (Y(4)*BNX + Y(5)*BNY + Y(6)*BNZ)*VC/omega
	
       RL  =SQRT(Y(1)**2+Y(2)**2)
       Rr_IDEI = RL-RR
       Y10 = (Rr_IDEI/RL)*Y(1)
       Y20 = (Rr_IDEI/RL)*Y(2)
       Y30 = Y(3)
       Y10 = Y10 / SQRT(Y10**2+Y20**2+Y30**2)
       Y20 = Y20 / SQRT(Y10**2+Y20**2+Y30**2)
       Y30 = Y30 / SQRT(Y10**2+Y20**2+Y30**2)
       OX_NX = (Y(4)*Y10 + Y(5)*Y20 + Y(6)*Y30)*VC/omega
		
       OX_N2 = (Y(4)**2+Y(5)**2+Y(6)**2)*(VC/omega)*(VC/omega)
       OX_NY = OX_N2 - OX_NZ**2 - OX_NX**2
       IF (OX_NY.LT.0.D0) OX_NY=0.0D0
       OX_NY = SQRT(OX_NY)
			
       RETURN
  END SUBROUTINE REFINDEX

!  ----- calculate wave number satifying D=0 -----

  SUBROUTINE WRMODNWTN(Y_zh,omega,YK)

    USE wrcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Y_zh(NEQ),omega
    REAL(rkind),INTENT(OUT):: YK(3)
    REAL(rkind):: Y(0:NEQ)
    INTEGER:: I,IMODNWTN
    REAL(rkind):: DELTA,XP,YP,ZP,VV,TT
    REAL(rkind):: RKXP,RKYP,RKZP,RRKXP,RRKYP,RRKZP,DKXP,DKYP,DKZP
    REAL(rkind):: DS2,FAC_NWTN,DDELTA

    Y(1:NEQ)=Y_zh(1:NEQ)
    DELTA=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), omega )
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

       DKXP=(DISPXR(XP,YP,ZP,RKXP+RRKXP,RKYP,RKZP,omega) &
            -DISPXR(XP,YP,ZP,RKXP-RRKXP,RKYP,RKZP,omega))/(2.D0*RRKXP)
       DKYP=(DISPXR(XP,YP,ZP,RKXP,RKYP+RRKYP,RKZP,omega) &
            -DISPXR(XP,YP,ZP,RKXP,RKYP-RRKYP,RKZP,omega))/(2.D0*RRKYP)
       DKZP=(DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP+RRKZP,omega) &
            -DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP-RRKZP,omega))/(2.D0*RRKZP)

       DS2 = DKXP**2+DKYP**2+DKZP**2
       FAC_NWTN = 10.D0
       DO IMODNWTN=1,10
          FAC_NWTN = FAC_NWTN/10.D0
          Y(4) = RKXP - (DELTA*DKXP)/DS2*FAC_NWTN
          Y(5) = RKYP - (DELTA*DKYP)/DS2*FAC_NWTN
          Y(6) = RKZP - (DELTA*DKZP)/DS2*FAC_NWTN
          DDELTA = DISPXR(XP,YP,ZP,Y(4),Y(5),Y(6),omega)
          IF (ABS(DDELTA) .LT. ABS(DELTA)) EXIT
       END DO
       IF(idebug_wr(13).NE.0) THEN
          WRITE(6,'(A)') '*** idebug_wr(13):'
          WRITE(6,'(A,4ES12.4)') &
               '      Y(4,5,6),DDELTA=',Y(4),Y(5),Y(6),DDELTA
       END IF
       IF (ABS(DDELTA).LT.EPSD0) EXIT
       DELTA = DDELTA
    END DO

    DO I =1,3
       YK(I) = Y(I+3)
    END DO

    RETURN
  END SUBROUTINE WRMODNWTN

!  ----- calculate real part of dispersion relation Re D -----

  FUNCTION DISPXR(XP,YP,ZP,RKXP,RKYP,RKZP,omega)

    USE wrcomm
    USE plprof
    USE dpdisp,ONLY: cfdispr
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,omega
    TYPE(pl_mag_type):: MAG
    REAL(rkind):: DISPXR
    INTEGER:: MODELPS(NSMAX),MODELVS(NSMAX)
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

      CRF=DCMPLX(omega/(2.D6*PI),0.D0)
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

!      CWW=2.D0*PI*1.D6*CRF
!      DO NS=1,NSMAX
!         IF(NSDP(NS).EQ.1) THEN
!            CWC=MAG%BABS*PZ(NS)*AEE/(AMP*PA(NS))
!            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
!         ENDIF
!      ENDDO

      DO NS=1,NSMAX
         MODELP(NS)=MODELPS(NS)
         MODELV(NS)=MODELVS(NS)
      ENDDO

      DISPXR=DBLE(CF)
      RETURN
  END FUNCTION DISPXR

!  ----- calculate imaginary part of dispersion relation Im D -----

  FUNCTION DISPXI(XP,YP,ZP,RKXP,RKYP,RKZP,omega)

    USE wrcomm
    USE plprof
    USE dpdisp,ONLY: cfdisp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: XP,YP,ZP,RKXP,RKYP,RKZP,omega
    TYPE(pl_mag_type):: MAG
    REAL(rkind):: DISPXI
    INTEGER:: NS
    REAL(rkind):: X,Y,Z
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CF,CWW,CWC

      CRF=DCMPLX(omega/(2.D6*PI),0.D0)
      CKX=RKXP
      CKY=RKYP
      CKZ=RKZP
      X=XP
      Y=YP
      Z=ZP
      CALL pl_mag(X,Y,Z,MAG)
            
      CF=CFDISP(CRF,CKX,CKY,CKZ,X,Y,Z)

!      CWW=2.D0*PI*1.D6*CRF
!      DO NS=1,NSMAX
!         IF(NSDP(NS).EQ.1) THEN
!            CWC=MAG%BABS*PZ(NS)*AEE/(AMP*PA(NS))
!            CF=CF*(CWW-ABS(CWC))/(CWW+ABS(CWC))
!         ENDIF
!      ENDDO

      DISPXI=DIMAG(CF)
      RETURN
  END FUNCTION DISPXI

!  ----- calculate eigen electric field and wave number -----

  SUBROUTINE WRCALE(RF,YN,NSTPMAX_L,NRAY)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RF,YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(IN):: NSTPMAX_L,NRAY
    COMPLEX(rkind):: CDET(3,3),CDETP(3,3),CDETM(3,3),CDETD(3,3)
    INTEGER:: NSTP,J,I
    REAL(rkind):: X1,Y1,Z1,VV,TT,domega,UE2,UE,RHON,EA,omega
    COMPLEX(rkind):: CRF,CKX1,CKY1,CKZ1,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ,CRFP,CRFM,CUE2

    omega=2.D6*PI*RF
      DO NSTP=0,NSTPMAX_L
         CRF =DCMPLX(RF,0.D0)
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

         domega=MAX(ABS(omega)*VV,TT)
         CRFP =DCMPLX((omega+domega)/(2.D6*PI),0.D0)
         CALL DP_DISP(CRFP,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETP)
         CRFM =DCMPLX((omega-domega)/(2.D6*PI),0.D0)
         CALL DP_DISP(CRFM,CKX1,CKY1,CKZ1,X1,Y1,Z1,CDETM)
         
         DO J=1,3
         DO I=1,3
            CDETD(I,J)=omega*(CDETP(I,J)-CDETM(I,J))/(2.D0*domega) &
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

  SUBROUTINE WRCALE_I(RF,R,PHI,Z,RKR,RKPH,RKZ,EPARA,EPERP)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RF,R,PHI,Z,RKR,RKPH,RKZ
    REAL(rkind),INTENT(OUT):: EPARA,EPERP
    COMPLEX(rkind):: CDET(3,3)
    REAL(rkind):: omega,X1,Y1,Z1,RHON,EA
    COMPLEX(rkind):: CRF,CKX1,CKY1,CKZ1,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ

    omega=2.D6*PI*RF
    CRF=DCMPLX(RF,0.D0)

    IF(MODELG.EQ.0.OR.MODELG.EQ.1.OR.MODELG.EQ.11) THEN
       X1= R
       Y1= PHI
       Z1= Z
       CKX1= RKR
       CKY1= RKPH
       CKZ1= RKZ
    ELSE
       X1= R*COS(PHI)
       Y1= R*SIN(PHI)
       Z1= Z
       CKX1= RKR*COS(PHI)-RKPH*SIN(PHI)
       CKY1= RKR*SIN(PHI)+RKPH*COS(PHI)
       CKZ1= RKZ
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

  SUBROUTINE WRCALE_XYZ(RF,X,Y,Z,RKX,RKY,RKZ,EPARA,EPERP)

    USE wrcomm
    USE plprof,ONLY: PL_MAG_OLD
    USE pllocal
    USE dpdisp,ONLY: dp_disp
    USE plcomm_type,ONLY: pl_mag_type
    USE plprof,ONLY: pl_mag
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RF,X,Y,Z,RKX,RKY,RKZ
    REAL(rkind),INTENT(OUT):: EPARA,EPERP
    TYPE(pl_mag_type):: mag
    COMPLEX(rkind):: CDET(3,3)
    REAL(rkind):: omega,EA
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CE1,CE2,CE3,CE4,CEXY,CEZY
    COMPLEX(rkind):: CUEX,CUEY,CUEZ

    omega=2.D6*PI*RF
    CRF=DCMPLX(RF,0.D0)
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
