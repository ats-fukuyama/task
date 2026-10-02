! wroxb.f90

!     developed by R. Miyata (Kyushu U) 2023

MODULE wroxb

  PRIVATE
  PUBLIC WRRKFT_RKF_OXB

CONTAINS

! --- Auto-step-size Runge-Kutta-F method (OXB) ---

  SUBROUTINE WRRKFT_RKF_OXB(nstp_start,YN,nstp_end)

    USE wrcomm
    USE pllocal
    USE plprof
    USE wrfdrv
    USE wrsub,ONLY: wr_write_line
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: nstp_start
    REAL(rkind),INTENT(OUT):: YN(0:NEQ,0:NSTPMAX)
    INTEGER,INTENT(OUT):: nstp_end
    REAL(rkind):: Y(NEQ),Y_pre(NEQ),YK(NEQ)
    INTEGER:: NSTPLIM,NSTPMIN,NSTP,INIT,NDE,IER,I
    REAL(rkind):: RELERR,ABSERR,X0,XE,WORK0,PW,YM(NEQ),RHON,R
    REAL(rkind):: ESTERR(NEQ),WORK1(NEQ),WORK2(NEQ),WORK3(NEQ),WORK4(NEQ,11)
    REAL(rkind):: RF,omega,rkv,IOX,RL,RKRL,DDELS,OXEFF,DELTA,RK_b,X_pre
    RELERR = EPSRAY
    ABSERR = EPSRAY
    INIT = 1
    RF=wr_nray_status%RF
    omega=2.D6*PI*RF
    rkv = omega/VC
    NSTPLIM=MAX(INT(SMAX/DELS),NSTPMAX)
    NSTPMIN=INT(0.1D0*SMAX/DELS)

    NSTP=nstp_start
    X0 = YN(0,NSTP)
    XE = X0+DELS
    DO I=1,7
       Y(I)=YN(I,NSTP)
    ENDDO

    IOX = 0

    IF(idebug_wr(11).NE.0) THEN
       WRITE(6,'(A,I8)') '*** idebug_wr(12): wrrkft_rkf: nstp=',nstp
       WRITE(6,'(A,4ES12.4)') '      x0,y1,y2,y3 =',X0,Y(1),Y(2),Y(3)
       WRITE(6,'(A,4ES12.4)') '      y4,y5,y6,y7 =',Y(4),Y(5),Y(6),Y(7)
    END IF

       CALL wr_write_line(NSTP,X0,Y,YN(8,NSTP))

    DO NSTP = nstp_start+1,NSTPLIM
       PW=Y(7)
       CALL RKF(7,wr_fdrv,X0,XE,Y,INIT,RELERR,ABSERR,YM, &
                ESTERR,NDE,IER,WORK0,WORK1,WORK2,WORK3,WORK4)

       DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),omega)

       IF (ABS(DELTA).GT.1.0D-6) THEN
            CALL WRMODNWTN(YM,omega,YK)
            DO I=1,3
                YM(I+3) = YK(I)
            ENDDO
       ENDIF

!        RK_b = SQRT(YM(4)**2+YM(5)**2+YM(6)**2)
!
 !       IF (ABS(RK_b).gt.ABS(rkv)*2.AND.IOX.eq.1) THEN
  !          DELS=0.5*DELS
    !        IOX = 2
   !     ENDIF

       IF (IER .NE. 0) THEN
          WRITE(6,'(A,2I6)') 'XX wrrkft_rkf: NSTP,IER=',NSTP,IER
          RETURN
       ENDIF

       CALL VECR_2_PROJ(YM(1),YM(2),YM(3),Y(4),YM(5),YM(6),RKRL)
       RL  =SQRT(YM(1)**2+YM(2)**2+YM(3)**2)
       RKRL=RKRL/RL
        IF ( RKRL.GE.0.D0.AND.IOX.EQ.0 ) THEN

            DO I=1,7
                Y(I)=Y_pre(I)
            ENDDO
            DDELS = 0.1D0*DELS
            X0 = X0 - DELS
            XE = X0 + DDELS

            DO
                PW=Y(7)
                CALL RKF(7,wr_fdrv,X0,XE,Y,INIT,RELERR,ABSERR,YM, &
                         ESTERR,NDE,IER,WORK0,WORK1,WORK2,WORK3,WORK4)
                DELTA=DISPXR(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),omega)
                
                !IF (ABS(DELTA).GT.1.0D-6) THEN
                 !   CALL WRMODNWTN(YM,omega,YK)
                  !  DO I=1,3
                   !     YM(I+3) = YK(I)
                    !END DO
                !END IF

                CALL VECR_2_PROJ(YM(1),YM(2),YM(3),YM(4),YM(5),YM(6),RKRL)
                RL  = SQRT( YM(1)**2 + YM(2)**2 + YM(3)**2 )
                RKRL = RKRL / RL

                IF (RKRL.GT.0.D0) THEN
                    DO I=1,7
                        Y(I) = Y_pre(I)
                    ENDDO
                    X0 = X0 - DDELS
                    DDELS = 0.1D0*DDELS
                    XE = X0 + DDELS

                ELSEIF (RKRL.LT.0.D0 .AND. RKRL.GT.-rkv*1.0D-2) THEN
                    DO I=1,7
                        Y(I) = YM(I)
                    ENDDO
                    EXIT

                ELSE
                    DO I=1,7
                        Y(I) = YM(I)
                    ENDDO

                    DO I=1,7
                       Y_pre(I)=YM(I)
                    END DO

                    XE =X0+DDELS
                ENDIF

            END DO

            YN(0,NSTP)=XE
            DO I=1,7
                YN(I,NSTP) = Y(I)
            ENDDO
            YN(8,NSTP)=PW-YM(7)

            CALL wr_write_line(NSTP,XE,Y,YN(8,NSTP))
            CALL WRMODCONV(IOX,Y,YM,OXEFF,omega)

        ENDIF

        DO I=1,7
            Y_pre(I)=YM(I)
        END DO


       YN(0,NSTP)=XE
       DO I=1,7
          YN(I,NSTP)=YM(I)
       ENDDO
       YN(8,NSTP)=PW-YM(7)

       CALL wr_write_line(NSTP,XE,YM,YN(8,NSTP))

       DO I=1,7
          Y(I)=YM(I)
       ENDDO
       X0=XE
       XE=X0+DELS

       R=SQRT(Y(1)**2+Y(2)**2)
       IF(Y(7).LT.UUMIN.OR. &
         (NSTP.GT.NSTPMIN.AND. &
         (R.GT.Rmax_wr.OR. &
          R.LT.Rmin_wr.OR. &
          Y(3).GT.Zmax_wr.OR. &
          Y(3).LT.Zmin_wr))) THEN
          nstp_end = NSTP
          GOTO 8000
       ENDIF
       IF(MODELG.LE.10) THEN
          CALL pl_mag_old(Y(1),Y(2),Y(3),RHON)
          IF(RHON.GT.RB/RA) THEN
             nstp_end = NSTP
             GOTO 8000
          END IF
       ENDIF
    END DO
    nstp_end=NSTPLIM
     
8000 CONTINUE
    IF(YN(7,nstp_end).LT.0.D0) YN(7,nstp_end)=0.D0
    CALL wr_write_line(NSTP,X0,YM,YN(8,NSTP))

    RETURN
  END SUBROUTINE WRRKFT_RKF_OXB

  

!  ----- mode conversion -----
    SUBROUTINE WRMODCONV(IOX,Y,YM,OXEFF,omega)

      USE wrcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
      IMPLICIT NONE
      REAL(rkind),INTENT(OUT):: IOX
      REAL(rkind),INTENT(IN):: omega
      REAL(rkind),INTENT(INOUT):: Y(NEQ)
      REAL(rkind),INTENT(OUT):: YM(NEQ)
      REAL(rkind),INTENT(OUT):: OXEFF
      INTEGER:: I
      REAL(rkind):: Y10,Y20,Y30,R_a
      REAL(rkind):: OX_K0,OX_NE,OX_LN,OX_DNE_DX,OX_DNE_DY,OX_DNE_DZ,OX_DNE_DD,OX_Y,OX_NZOPT,OX_NZ,OX_NY,OX_NX
      REAL(rkind):: RL0,Rr_IDEI,S_O_X,DELTA,DELTAB
        
         RL0 = SQRT(Y(1)**2+Y(2)**2)
         Rr_IDEI = RL0-RR
         Y10 = (Rr_IDEI/RL0)*Y(1)
         Y20 = (Rr_IDEI/RL0)*Y(2)
         Y30 = Y(3)
         R_a = SQRT(Y10**2+Y20**2+Y30**2)
         CALL LAMBDA_N_OX(Y(1), Y(2), Y(3),OX_NE, OX_LN , OX_DNE_DX , OX_DNE_DY , OX_DNE_DZ , OX_DNE_DD)
         CALL REFINDEX(omega, Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY, OX_NX)
         OX_K0 = omega / VC
         WRITE(6,*)'N_OPT=',OX_NZOPT,'NZ=',OX_NZ,'NY=',OX_NY,'NX=',OX_NX, 'ne=',OX_NE
         WRITE(6,*)'K0=',OX_K0,'N=', &
              SQRT((Y(4)**2+Y(5)**2+Y(6)**2))*(VC/omega),'OX_LN=',OX_LN, 'wc/w=',OX_Y

         OXEFF = ( 2.0*(1.0+OX_Y)*((OX_NZ-OX_NZOPT)**2) + OX_NY**2 )
         OXEFF = EXP(-PI*OX_K0*OX_LN*SQRT(0.5*OX_Y)*OXEFF)
         WRITE(6,*) 'OXEFF=',OXEFF,'Ra=',R_a
         
         S_O_X = 5.0D-5

         DELTAB=ABS(DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), omega ))

         DO I=1,100000
            CALL LAMBDA_N_OX(Y(1), Y(2), Y(3),OX_NE, OX_LN,OX_DNE_DX,OX_DNE_DY,OX_DNE_DZ,OX_DNE_DD)
            Y(1) = Y(1) + I * S_O_X * OX_DNE_DX/OX_DNE_DD
            Y(2) = Y(2) + I * S_O_X * OX_DNE_DY/OX_DNE_DD
            Y(3) = Y(3) + I * S_O_X * OX_DNE_DZ/OX_DNE_DD

            DELTA=DISPXR( Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), omega )

            IF ( DELTA*DELTAB.LT.0D0 ) THEN
                EXIT
            END IF
        END DO

        DO I =1,NEQ
            YM(I) = Y(I)
        END DO
        IOX =1
        RETURN
    END SUBROUTINE WRMODCONV

  SUBROUTINE WRMODCONV_OLD(IOX, Y, F, OXEFF,omega)

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
       CALL REFINDEX_OLD(omega, Y, OX_Y, OX_NZOPT, OX_NZ, OX_NY)
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
  END SUBROUTINE WRMODCONV_OLD

!  ----- calculate density scalelength -----
    SUBROUTINE  LAMBDA_N_OX(X,Y,Z, OX_NE, OX_LN,OX_DNE_DX,OX_DNE_DY,OX_DNE_DZ,OX_DNE_DD)

      USE wrcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD,PL_PROF_OLD
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,Y,Z
      REAL(rkind),INTENT(OUT):: OX_NE,OX_LN,OX_DNE_DX,OX_DNE_DY,OX_DNE_DZ,OX_DNE_DD
      REAL(rkind):: RHON
      REAL(rkind):: OX_NE_XP,OX_NE_XM,OX_NE_YP,OX_NE_YM,OX_NE_ZP,OX_NE_ZM
      REAL(rkind):: OX_DNE_DR

        CALL PL_MAG_OLD(X,Y,Z, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE = RN(1)

        CALL PL_MAG_OLD(X+DELDER, Y, Z, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE_XP = RN(1)
        CALL PL_MAG_OLD(X-DELDER, Y, Z, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE_XM = RN(1)
        OX_DNE_DX = (OX_NE_XP - OX_NE_XM)/(2.D0*DELDER)

        CALL PL_MAG_OLD(X, Y+DELDER, Z, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE_YP = RN(1)
        CALL PL_MAG_OLD(X, Y-DELDER, Z, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE_YM = RN(1)
        OX_DNE_DY = (OX_NE_YP - OX_NE_YM)/(2.D0*DELDER)

        CALL PL_MAG_OLD(X, Y, Z+DELDER, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE_ZP = RN(1)
        CALL PL_MAG_OLD(X, Y, Z-DELDER, RHON)
        CALL PL_PROF_OLD(RHON)
        OX_NE_ZM = RN(1)
        OX_DNE_DZ = (OX_NE_ZP - OX_NE_ZM)/(2.D0*DELDER)

        OX_DNE_DD = SQRT(OX_DNE_DX**2 + OX_DNE_DY**2 + OX_DNE_DZ**2)

        CALL VECR_2_PROJ(X,Y,Z, OX_DNE_DX , OX_DNE_DY , OX_DNE_DZ , OX_DNE_DR )

        OX_LN = OX_NE / OX_DNE_DD
    !      WRITE(6,*) 'NE(E18)=',OX_NE,'OX_LN=',OX_LN

      RETURN
    END SUBROUTINE LAMBDA_N_OX

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

    SUBROUTINE  REFINDEX(omega,Y,OX_Y,OX_NZOPT,OX_NZ, OX_NY,OX_NX)
      USE wrcomm
      USE pllocal
      USE plprof,ONLY: PL_MAG_OLD
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: omega,Y(NEQ)
      REAL(rkind),INTENT(OUT):: OX_Y,OX_NZOPT,OX_NZ,OX_NY,OX_NX
      REAL(rkind):: RHON,omega_C_OX,NX,NY,NZ,OX_DNE_DX,OX_DNE_DY,OX_DNE_DZ,OX_DNE_DD
      REAL(rkind):: OX_N2,OX_LN,OX_NE

        CALL PL_MAG_OLD(Y(1), Y(2), Y(3), RHON)
        omega_C_OX = BABS*AEE/AME
        OX_Y = omega_C_OX / omega

        OX_NZ = (Y(4)*BNX + Y(5)*BNY + Y(6)*BNZ)*VC/omega

        IF ( OX_NZ.LT.0.D0 ) THEN
        OX_NZOPT = -SQRT(OX_Y/(1.D0+OX_Y))
        ELSE
        OX_NZOPT = SQRT(OX_Y/(1.D0+OX_Y))
        END IF

        NX = Y(4)*VC/omega
        NY = Y(5)*VC/omega
        NZ = Y(6)*VC/omega

        !CALL VECR_2_PROJ(Y(1),Y(2),Y(3), NX , NY , NZ , OX_NX )
        CALL LAMBDA_N_OX(Y(1),Y(2),Y(3),OX_NE, OX_LN,OX_DNE_DX,OX_DNE_DY,OX_DNE_DZ,OX_DNE_DD)
        
        OX_NX = (NX * OX_DNE_DX + NY * OX_DNE_DY + NZ * OX_DNE_DZ)/OX_DNE_DD

        OX_N2 = (Y(4)**2+Y(5)**2+Y(6)**2)*(VC/omega)*(VC/omega)
        OX_NY = OX_N2 - OX_NZ**2 - OX_NX**2

        IF ( OX_NY.LT.0.D0 ) THEN
        OX_NY=0.0D0
        ELSE
        OX_NY = SQRT(OX_NY)
        END IF
        RETURN
    END SUBROUTINE REFINDEX

  SUBROUTINE  REFINDEX_OLD(omega,Y,OX_Y,OX_NZOPT,OX_NZ, OX_NY)

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
  END SUBROUTINE REFINDEX_OLD

    SUBROUTINE VECR_2_PROJ(OX_X,OX_Y,OX_Z,A,B,C,OX_DR)

      USE wrcomm
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: OX_X,OX_Y,OX_Z,A,B,C
      REAL(rkind),INTENT(OUT):: OX_DR
      REAL(rkind):: RL,Rr_IDEI
      REAL(rkind):: OX_X0,OX_Y0,OX_Z0
      REAL(rkind):: PHI, THETA

        RL = SQRT(OX_X**2+OX_Y**2)
        Rr_IDEI = RL - RR
        OX_X0 = Rr_IDEI * (OX_X/RL)
        OX_Y0 = Rr_IDEI * (OX_Y/RL)
        OX_Z0 = OX_Z
        PHI = ATAN2(OX_Y0, OX_X0)
        THETA = ATAN2(SQRT(OX_X0**2 + OX_Y0**2),OX_Z0)

        OX_DR = SIN(THETA)*COS(PHI)*A
        OX_DR = OX_DR + SIN(THETA)*SIN(PHI)*B
        OX_DR = OX_DR + COS(THETA)*C

        RETURN
    END SUBROUTINE VECR_2_PROJ

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
          !WRITE(6,*) DDELTA,Y(4),Y(5),Y(6)
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

  SUBROUTINE RKF(NEQ,FUNC,X0,XE,YH0,INIT,RELERR,ABSERR, &
                 YHN,ESTERR,NDE,IER,H,YL0,YLN,ERR, &
                 WORK)

!*********************************************************************
!     SUBROUTINE RKF NUMERICALLY INTEGRATES A SYSTEM OF NEQ          *
!     FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM        *
!             DY(I)/DX = F(X, Y(1),..., Y(NEQ)),                     *
!     BY THE RUNGE-KUTTA-FEHLBERG (4,5) FORMULA.                     *
!                                                                    *
!     PARAMETERS                                                     *
!  === INPUT ===                                                     *
!     (1) NEQ: NUMBER OF EQUATIONS TO BE INTEGRATED                  *
!     (2) FUNC: SUBROUTINE FUNC(X,Y,F) TO EVALUATE DERIVATIVES       *
!                F(I)=DY(I)/DX                                       *
!     (3) X0: INITIAL VALUE OF INDEPENDENT VARIABLE                  *
!     (4) XE: OUTPUT POINT AT WHICH THE SOLUTION IS DESIRED          *
!     (5) YH0(I) (I=1,..,NEQ): INITIAL VALUE AT X0                   *
!     (6) INIT: INDICATOR TO INITIALIZE THE CODE                     *
!          INIT=1..THE CODE WILL BE INITIALIZED(FIRST CALL).         *
!          INIT=2..THE CODE WILL NOT BE INITIALIZED(SUBSEQUENT CALL).*
!     (7) RELERR: RELATIVE LOCAL ERROR TOLERANCE                     *
!     (8) ABSERR: ABSOLUTE LOCAL ERROR TOLERANCE                     *
!  === OUTPUT ===                                                    *
!     (9) YHN(I) (I=1,..,NEQ): APPROXIMATE SOLUTION AT XE            *
!    (10) ESTERR(I) (I=1,..,NEQ): ESTIMATE OF THE GLOBAL ERROR IN    *
!          THE APPROXIMATE SOLUTION YHN(I)                           *
!    (11) NDE: NUMBER OF DERIVATIVE EVALUATIONS                      *
!    (12) IER: INDICATOR FOR STATUS OF INTEGRATION                   *
!          IER=0..INTEGRATION REACHED XE. INDICATES SUCCESSFUL       *
!             RETURN AND IS THE NORMAL MODE FOR CONTINUING THE       *
!             INTEGRATION.                                           *
!          IER=10000..INTEGRATION WAS NOT COMPLETED BECAUSE TOO MANY *
!             DERIVATIVES EVALUATIONS WERE NEEDED.  IF THE USER WANTS*
!             TO CONTINUE THE INTEGRATION, HE JUST CALLS RKF AGAIN.  *
!          IER=20000..INTEGRATION WAS NOT COMPLETED BECAUSE ERROR    *
!             TOLERANCE WAS INAPPROPRIATE(=IT WAS ZERO).  MUST USE   *
!             NON-ZERO ABSERR TO CONTINUE THE INTEGRATION.           *
!          IER=30000..INTEGRATION WAS NOT COMPLETED BECAUSE REQUESTED*
!             ACCURACY COULD NOT BE ACHIEVED USING SMALLEST STEPSIZE.*
!             MUST INCREASE THE ERROR TOLERANCE TO CONTINUE THE      *
!             INTEGRATION.                                           *
!  === OTHERS ===                                                    *
!    (13) H: VARIABLE TO HOLD INFORMATION INTERNAL TO RKF WHICH IS   *
!           NECESSARY FOR SUBSEQUENT CALLS.                          *
!    (14) YL0(), YLN(): ARRAY (SIZE=NEQ) TO HOLD INFORMATION INTERNAL*
!           TO RKF WHICH IS NECESSARY FOR SUBSEQUENT CALLS.          *
!    (15) ERR(): ARRAY (SIZE=NEQ) TO BE USED INSIDE RKF              *
!    (16) WORK(): TWO-DIMENTIONAL ARRAY (SIZE=(NEQ,11)) TO BE        *
!                 USED INSIDE RKF                                    *
!    COPYRIGHT: M. SUGIHARA, NOVEMBER 15, 1989    V. 1               *
!*********************************************************************
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQ
    REAL(dp),INTENT(INOUT):: X0,YH0(NEQ),XE
    REAL(dp),INTENT(IN):: RELERR,ABSERR
    REAL(dp),INTENT(OUT):: YHN(NEQ),ESTERR(NEQ)
    INTEGER,INTENT(INOUT):: INIT
    INTEGER,INTENT(OUT):: NDE,IER
    REAL(dp),INTENT(INOUT):: H,YL0(NEQ),YLN(NEQ),ERR(NEQ),WORK(NEQ,11)
    EXTERNAL FUNC
!    INTEGER,PARAMETER:: ITEMAX=100000
    INTEGER,PARAMETER:: ITEMAX=1000
    REAl(dp),PARAMETER:: EPSMIN=1.0D-15
    INTEGER:: I,ITER
    REAL(dp):: TOLDY,TOL,HMIN,ERRET,ET,X

    IF (INIT .EQ. 1) THEN

!   -------------- INITIALIZATION (FIRST CALL)--------------------------

       DO I = 1,NEQ
          YL0(I) = YH0(I)
       END DO
       
!     -------- SET INITIAL STEP SIZE ---------

       CALL FUNC(X0,YL0,WORK(1:NEQ,1))
       TOLDY = 10.0D+74
       DO I = 1,NEQ
          IF (WORK(I,1) .NE. 0.0D0) THEN
             TOL = RELERR * DABS(YL0(I)) + ABSERR
             TOLDY = MIN(TOLDY, ABS(TOL / WORK(I,1)))
          END IF
       END DO
       HMIN = EPSMIN * (XE - X0)
       IF (TOLDY .EQ. 10.0D+74) THEN
          H = HMIN
       ELSE
          H = MIN((XE - X0) / 100.0D0, TOLDY ** 0.2D0)
       END IF
!    ----------------------------------------
       INIT = 2
       NDE = 1
    ELSE
       NDE = 0
    END IF
!  ------------------------------------------------------------------
    IER = 0
    X = X0
    ITER = 0

!********************* MAIN ITERATION *********************************

    !DO ITER = 1,ITEMAX
    DO
       ITER = ITER + 1
        !WRITE(6,*) ITER
       CALL FLSTEP(NEQ,FUNC,X,H,YL0,YLN,ERR,WORK)
       NDE = NDE + 6
!   ----------- COMPUTE ERROR TOLERANCES ---------------------
       ERRET = 0.0D0
       DO I = 1,NEQ
          ET = RELERR * 0.5D0 * (DABS(YL0(I)) + DABS(YLN(I))) + ABSERR
          IF (ET .EQ. 0.0D0) THEN
             GO TO 20000
          ELSE
             ERRET = MAX(ERRET, ABS(ERR(I)) / ET)
          END IF
       END DO
       IF (ERRET .GT. 1.0D0) THEN
          
!   ----------- UNSUCCESSFUL STEP ------------------------------

          IF (ERRET .GE. 59049.0D0) THEN
             H = 0.1D0 * H
          ELSE
             H = 0.9D0 * H / ERRET ** 0.2D0
          END IF
          IF (H .LE. EPSMIN) GO TO 30000
       ELSE IF ((X + H) .LT. XE) THEN
          !WRITE(6,*) ITER
          
!   ----------- SUCCESSFUL STEP (X+H < THE END POINT) -----------

          CALL FHSTEP(NEQ,FUNC,X,H,YH0,YHN,ERR,WORK)
          NDE = NDE + 12
          X = X + H
          X0 = X
          DO I = 1,NEQ
             YL0(I) = YLN(I)
             YH0(I) = YHN(I)
          END DO
          
          IF (ITER.GE.ITEMAX) THEN
             XE = X0
            RETURN
            ENDIF
          
!       ------- CHOOSE NEXT STEP --------

          IF (ERRET .LE. 1.889568D-4) THEN
             H = 5.0D0 * H
          ELSE
             H = 0.9D0 * H / ERRET ** 0.2D0
          END IF
!       ---------------------------------
       ELSE
          
!   ------------ SUCCESSFUL STEP (X+H > =  THE END POINT) ----------

          H = XE - X
          CALL FLSTEP(NEQ,FUNC,X,H,YL0,YLN,ERR,WORK)
          CALL FHSTEP(NEQ,FUNC,X,H,YH0,YHN,ERR,WORK)
          NDE = NDE + 18
          
!       -------- ESTIMATE GLOBAL ERROR -------

          DO I = 1,NEQ
             ESTERR(I) = (YLN(I) - YHN(I)) / 31.0D0
          END DO
!       --------------------------------------
          X0 = XE
          DO I = 1,NEQ
             YL0(I) = YLN(I)
             YH0(I) = YHN(I)
          END DO
          RETURN
       ENDIF
    END DO
    
!**********************************************************************

    IER = 10000
    WRITE(6,10001) X
10001 FORMAT(' ','(SUBR.-RKF) TROUBLE(TOO MANY ITERATIONS))', &
           ' AT X = ',1PE15.7)
    RETURN
20000 CONTINUE
    IER = 20000
    WRITE(6,20001) X
20001 FORMAT(' ','(SUBR.-RKF) TROUBLE(INAPPROPRIATE ERROR TOLERANCE)', &
           ' AT X = ',1PE15.7)
    RETURN
30000 CONTINUE
    IER = 30000
    WRITE(6,30001) X
30001 FORMAT(' ','(SUBR.-RKF) TROUBLE(TOO SMALL STEP SIZE)', &
           ' AT X = ',1PE15.7)
    RETURN
  END SUBROUTINE RKF

  SUBROUTINE FLSTEP(NEQ,FUNC,X,H,Y0,YN,ERR,WORK)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQ
    REAL(dp),INTENT(IN):: X,H,Y0(NEQ)
    REAL(dp),INTENT(OUT):: YN(NEQ)
    REAL(dp),INTENT(INOUT):: ERR(NEQ),WORK(NEQ,11)
    EXTERNAL FUNC
    INTEGER:: LINDEX

    LINDEX = 1
    CALL FEHL(LINDEX,NEQ,FUNC,X,H,Y0,YN,ERR,WORK(1,1),WORK(1,2), &
              WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7), &
              WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
    RETURN
  END SUBROUTINE FLSTEP

! --- internal one-step routine ---
    
  SUBROUTINE FHSTEP(NEQ,FUNC,X,H,Y0,YN,ERR,WORK)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQ
    REAL(dp),INTENT(IN):: X,H
    REAL(dp),INTENT(OUT):: YN(:)
    REAL(dp),INTENT(INOUT):: Y0(NEQ),ERR(NEQ),WORK(NEQ,11)
    EXTERNAL FUNC
    INTEGER:: LINDEX,I
    REAL(dp):: X1,H1

    X1 = X
    H1 = 0.5D0 * H
    LINDEX = 0
    CALL FEHL(LINDEX,NEQ,FUNC,X1,H1,Y0,YN,ERR,WORK(1,1),WORK(1,2), &
              WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7), &
              WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
    X1 = X1 + H1
    DO I = 1,NEQ
       Y0(I) = YN(I)
    END DO
    CALL FEHL(LINDEX,NEQ,FUNC,X1,H1,Y0,YN,ERR,WORK(1,1),WORK(1,2), &
              WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7), &
              WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
    RETURN
  END SUBROUTINE FHSTEP

! ---       
      
  SUBROUTINE FEHL(LINDEX,NEQ,FUNC,X,H,Y0,YN,ERR,AK1,AK2,AK3,AK4, &
                  AK5,AK6,W2,W3,W4,W5,W6)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: LINDEX,NEQ
    REAL(dp),INTENT(IN):: X,H,Y0(NEQ)
    REAL(dp),INTENT(OUT):: YN(NEQ)
    REAL(dp),INTENT(OUT):: ERR(NEQ)
    REAL(dp),INTENT(OUT):: AK1(NEQ),AK2(NEQ),AK3(NEQ)
    REAL(dp),INTENT(OUT):: AK4(NEQ),AK5(NEQ),AK6(NEQ)
    REAL(dp),INTENT(OUT):: W2(NEQ),W3(NEQ),W4(NEQ),W5(NEQ),W6(NEQ)
    INTEGER:: I
!    EXTERNAL FUNC
    INTERFACE
       SUBROUTINE FUNC(X,Y,F)
         USE task_kinds,ONLY: dp
         IMPLICIT NONE
         REAL(dp),INTENT(IN):: X,Y(*)
         REAL(dp),INTENT(OUT):: F(*)
       END SUBROUTINE FUNC
    END INTERFACE

    REAL(dp),PARAMETER:: ONE = 1.D0, TWO = 2.D0, THR = 3.D0, TWL = 12.D0
    REAL(dp),PARAMETER:: &
         AL2 = ONE / 4.D0, AL3 = THR / 8.D0, AL4 = TWL / 13.D0, &
         AL5 = ONE,        AL6 = ONE / 2.D0
    REAL(dp),PARAMETER:: &
         B21 = ONE / 4.D0, B31 = THR / 32.D0, B32 = 9.0D0 / 32.D0, &
         B41 = 1932.0D0 / 2197.D0, B42 = -7200.0D0 / 2197.D0,  &
         B43 = 7296.0D0 / 2197.D0, B51 = 439.0D0 / 216.D0, B52 = -8.D0, &
         B53 = 3680.0D0 / 513.D0,  B54 = -845.0D0 / 4104.D0, &
         B61 = -8.0D0 / 27.D0, B62 = 2.D0, B63 = -3544.0D0 / 2565.D0, &
         B64 = 1859.0D0 / 4104.D0, B65 = -11.0D0 / 40.D0
    REAL(dp),PARAMETER:: &
         GA1 = 16.0D0 / 135.D0, GA3 = 6656.0D0 / 12825.D0, &
         GA4 = 28561.0D0 / 56430.D0, &
         GA5 = -9.0D0 / 50.D0, GA6 = TWO / 55.D0
    REAL(dp),PARAMETER:: &
         DA1 = ONE / 360.D0, DA3 = -128.0D0 / 4275.D0, &
         DA4 = -2197.0D0 / 75240.D0, DA5 = ONE / 50.D0, DA6 = TWO / 55.D0

    CALL FUNC(X,Y0,AK1)
    DO I = 1,NEQ
       W2(I) = Y0(I) + H * B21 * AK1(I)
    END DO
    CALL FUNC(X + AL2 * H,W2,AK2)
    DO I = 1,NEQ
       W3(I) = Y0(I) + H * (B31 * AK1(I)+B32 * AK2(I))
    END DO
    CALL FUNC(X+AL3 * H,W3,AK3)
    DO I = 1,NEQ
       W4(I) = Y0(I) + H * (B41 * AK1(I) + B42 * AK2(I) + B43 * AK3(I))
    END DO
    CALL FUNC(X + AL4 * H,W4,AK4)
    DO I = 1,NEQ
       W5(I) = Y0(I) + H * (B51 * AK1(I) + B52 * AK2(I) &
                          + B53 * AK3(I) + B54 * AK4(I))
    END DO
    CALL FUNC(X + AL5 * H,W5,AK5)
    DO I = 1,NEQ
       W6(I) = Y0(I) + H * (B61 * AK1(I) + B62 * AK2(I) + B63 * AK3(I) &
                          + B64 * AK4(I) + B65 * AK5(I))
    END DO
    CALL FUNC(X + AL6 * H,W6,AK6)
    DO I = 1,NEQ
       YN(I) = Y0(I) + H * (GA1 * AK1(I) + GA3 * AK3(I) + GA4 * AK4(I) &
                          + GA5 * AK5(I) + GA6 * AK6(I))
    END DO
    IF (LINDEX .EQ. 1) THEN
       DO I = 1,NEQ 
          ERR(I) = H * (DA1 * AK1(I) + DA3 * AK3(I) + DA4 * AK4(I) &
                      + DA5 * AK5(I) + DA6 * AK6(I))
       END DO
    ENDIF
    RETURN
  END SUBROUTINE FEHL

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
  END MODULE wroxb
