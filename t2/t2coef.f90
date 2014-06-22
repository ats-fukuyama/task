!--------------------------------------------------------------------
!
! MODULE T2COEF
! 
!
!      CALCULATION OF COEFFICIENTS 
!      OF ADVECTION-DIFFUSION EQUATIONS 
!      FOR MULTI-FLUID TRANSPORT MODEL 
!      (ONLY ANOMALOUS TRANSPORT MODEL IS BASED ON TWO-FLUID MODEL)
!
!                       2014-06-22 H. Seto
!
!-------------------------------------------------------------------
MODULE T2COEF
  
  USE T2CNST, ONLY:ikind,rkind
  
  IMPLICIT NONE
  PRIVATE

  INTEGER(ikind),SAVE::i_m
    
  PUBLIC T2COEF_EXECUTE
  
CONTAINS
  
  SUBROUTINE T2COEF_EXECUTE
    
    USE T2COMM,ONLY:&
         NMMAX,CoordinateSwitch,&
         UsePotentialDescription,UseCoefficientCheck
    USE T2COUT,ONLY: T2COUT_EXECUTE
    USE T2CALV,ONLY: T2CALV_EXECUTE
    
    SELECT CASE (CoordinateSwitch)
    CASE (1)

       IF(.NOT.UsePotentialDescription)THEN 
          DO i_m = 1, NMMAX
             CALL T2CALV_EXECUTE(   i_m)
             CALL T2COEF_KV_EB(     i_m)
             CALL T2COEF_MS_COEF_EB(i_m)
             CALL T2COEF_AV_COEF_EB(i_m)
             CALL T2COEF_AT_COEF_EB(i_m)
             CALL T2COEF_DT_COEF_EB(i_m)
             CALL T2COEF_GV_COEF_EB(i_m)
             CALL T2COEF_GT_COEF_EB(i_m)
             CALL T2COEF_ES_COEF_EB(i_m)
             CALL T2COEF_EV_COEF_EB(i_m)
             CALL T2COEF_ET_COEF_EB(i_m)
             CALL T2COEF_SS_COEF_EB(i_m)
             
             !IF(i_m.EQ.1)THEN
             !   CALL T2COEF_CHECK(i_m)
             !   STOP
             !ENDIF
             IF(UseCoefficientCheck) CALL T2COUT_EXECUTE
             
          ENDDO
       ELSE
          
          DO i_m = 1, NMMAX
             
             CALL T2CALV_EXECUTE(     i_m)         
             CALL T2COEF_MS_COEF_PhiA(i_m)
             CALL T2COEF_AV_COEF_PhiA(i_m)
             CALL T2COEF_AT_COEF_PhiA(i_m)
             CALL T2COEF_DT_COEF_PhiA(i_m)
             CALL T2COEF_GV_COEF_PhiA(i_m)
             CALL T2COEF_GT_COEF_PhiA(i_m)
             CALL T2COEF_ES_COEF_PhiA(i_m)
             CALL T2COEF_EV_COEF_PhiA(i_m)
             CALL T2COEF_ET_COEF_PhiA(i_m)
             CALL T2COEF_SS_COEF_PhiA(i_m)
          ENDDO
       ENDIF
    CASE (2)
       DO i_m = 1, NMMAX
          CALL T2COEF_KV_TEST(     i_m)
          CALL T2COEF_MS_COEF_TEST(i_m)
          CALL T2COEF_AV_COEF_TEST(i_m)
          CALL T2COEF_AT_COEF_TEST(i_m)
          CALL T2COEF_DT_COEF_TEST(i_m)
          CALL T2COEF_GV_COEF_TEST(i_m)
          CALL T2COEF_GT_COEF_TEST(i_m)
          CALL T2COEF_ES_COEF_TEST(i_m)
          CALL T2COEF_EV_COEF_TEST(i_m)
          CALL T2COEF_ET_COEF_TEST(i_m)
          CALL T2COEF_SS_COEF_TEST(i_m)
       ENDDO
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_EXECUTE


  SUBROUTINE T2COEF_KV_EB(i_m)
    
    USE T2COMM,ONLY: NSMAX,R_rz,Bb,Ub,Wb,KnownVar
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_s,kOffsetX
    
    KnownVar(1,i_m) = Bb
    KnownVar(2,i_m) = R_rz
    
    DO i_s = 1,NSMAX
       kOffsetX = 2*(i_s-1)
       KnownVar(kOffsetX + 3,i_m) = Ub(i_s)
       KnownVar(kOffsetX + 4,i_m) = Wb(i_s)
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2COEF_KV_EB
  
  !---------------------------------------------------------
  ! 
  !       CALCULATION OF MASS SCALAR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-22 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_MS_COEF_EB(i_m)
    
    USE T2CNST, ONLY:VcSqRe
    USE T2COMM, ONLY:&
         & BpNF,  BtNF,  EtNF,  EpNF,        & 
         & NnNF,         FbNF,  FtNF,        &
         & PpNF,         QbNF,  QtNF,        &        
         & EqBpNF,EqBtNF,EqEtNF,EqEpNF,      &
         & EqNnNF,EqFrNF,EqFbNF,EqFtNF,      &
         & EqPpNF,EqQrNF,EqQbNF,EqQtNF,      &        
         !
         & NVMAX,NSMAX,NFMAX,GRt,G33Ct,Bb,R_mc,Mm,Dt,MassScaCoef,EqSet
    
    INTEGER(ikind),INTENT(IN)::i_m
    REAL(   rkind)::mmA
    INTEGER(ikind)::i_v,j_v,i_s,vOffsetA,vOffsetB
    
    ! INITIALIZATION
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       MassScaCoef(i_v,j_v,i_m) = 0.D0
    ENDDO
    ENDDO

    !
    ! variables as field (from i_v= 1 to i_v = NFMAX)
    !
    SELECT CASE (EqSet)
    CASE (1)
       i_v = 1                   ! Equation for psi'
       j_v = 1
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt              * BpNF/EqBpNF

       i_v = 2                   ! Equation for I
       j_v = 2 
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt              * BtNF/EqBtNF

       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 3
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt*VcSqRe       * EtNF/EqEtNF

       i_v = 4                   ! Equation for E_{\chi}
       j_v = 4
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt*VcSqRe*R_mc  * EpNF/EqEpNF
       
       i_v = 5                   ! Equation for E_{\rho}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          mmA = Mm(i_s)
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB  
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt              * NnNF/EqNnNF
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA*Bb       * FbNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA          * FtNF/EqFtNF
       
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*1.5D0        * PpNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  8 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA*Bb       * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  9 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA          * QtNF/EqQtNF
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    CASE (2)
       i_v = 1                   ! Equation for psi'
       j_v = 3
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt*VcSqRe*G33Ct * EtNF/EqBpNF
       
       i_v = 2                   ! Equation for I
       j_v = 4
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt*VcSqRe       * EpNF/EqBtNF
       
       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 1
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt              * BpNF/EqEtNF
       
       i_v = 4                   ! Equation for I
       j_v = 2 
       MassScaCoef(i_v,j_v,i_m) = GRt/Dt              * BtNF/EqEpNF
    
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
       
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
       
          mmA = Mm(i_s)
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB  
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt              * NnNF/EqNnNF
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA          * FtNF/EqFrNF
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA*Bb       * FbNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*1.5D0        * PpNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  9 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA          * QtNF/EqQrNF
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  8 + vOffsetB
          MassScaCoef(i_v,j_v,i_m) = GRt/Dt*mmA*Bb       * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
       
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_MS_COEF_EB
  
  !---------------------------------------------------------
  ! 
  !       CALCULATION OF ADVECTION VECTOR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-22 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_AV_COEF_EB(i_m)
    
    USE T2CNST,ONLY:VcSqRe
    USE T2COMM,ONLY:&
         & NSMAX,NFMAX,NVMAX,NDMAX,        &
         & BpNF,BtNF,EtNF,EpNF,ErNF, &
         & NnNF,     FbNF,FtNF,      &
         & PpNF,     QbNF,QtNF,      &
         & EqBpNF, EqBtNF, EqEtNF, EqEpNF, EqErNF, &
         & EqNnNF, EqFrNF, EqFbNF, EqFtNF,         &
         & EqPpNF, EqQrNF, EqQbNF, EqQtNF,         &
         !
         & Mm,Tt,Nn,UrCt,UpCt,Ub,UuSq,WrCt,WpCt,Pp,QrCt,QpCt,Wb,&
         & GRt,G11Ct,G12Ct,G22xCt,G33Ct,UgrCt,UgpCt,R_mc,Bb,BpCt,&
         & NVMAX,NDMAX,AdveVecCoef,EqSet
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_d,i_v,vOffsetA,&
         &         j_v,vOffsetB
    
    REAL(   rkind)::&
         & mmA,ttA,&
         & nnA,urCtA,upCtA,      ubA,uuSqA,&
         & ppA,qrCtA,qpCtA,      wbA,      &
         &     wrCtA,wpCtA,krCtA,kpCtA
    
    
    ! initialization    
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_d = 1, NDMAX
          AdveVecCoef(i_d,i_v,j_v,i_m) = 0.D0
       ENDDO
    ENDDO
    ENDDO
    
    SELECT CASE (EqSet)
    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       j_v = 1
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*UgrCt             * BpNF/EqBpNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*UgpCt             * BpNF/EqBpNF
       
       i_v = 2                   ! Equation for I
       j_v = 2
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*UgrCt             * BtNF/EqBtNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*UgpCt             * BtNF/EqBtNF

       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 1
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*G11Ct             * BpNF/EqEtNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*G12Ct             * BpNF/EqEtNF
       j_v = 3
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*ugrCt*VcSqRe      * EtNF/EqEtNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*ugpCt*VcSqRe      * EtNF/EqEtNF

       i_v = 4                   ! Equation for E_{\chi}
       j_v = 4
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*ugrCt*VcSqRe*R_mc * EpNF/EqEpNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*ugpCt*VcSqRe*R_mc * EpNF/EqEpNF

       i_v = 5                   ! Equation for E_{\rho}
       j_v = 4
       AdveVecCoef(1,i_v,j_v,i_m) =  GRt*G12Ct*R_mc        * EpNF/EqErNF
       AdveVecCoef(2,i_v,j_v,i_m) =  GRt*G22xCt            * EpNF/EqErNF
       j_v = 5
       AdveVecCoef(1,i_v,j_v,i_m) =  GRt*G11Ct             * ErNF/EqErNF
       AdveVecCoef(2,i_v,j_v,i_m) =  GRt*G12Ct             * ErNF/EqErNF
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          mmA = Mm(i_s)
          ttA = Tt(i_s)
          nnA = Nn(i_s)
          
          urCtA = UrCt(i_s)
          upCtA = UpCt(i_s)
          wrCtA = WrCt(i_s)
          wpCtA = WpCt(i_s)
          ubA   = Ub(  i_s)
          ppA   = Pp(  i_s)
          uuSqA = UuSq(i_s)
          wbA   = Wb(  i_s)
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)= GRt*(urCtA-UgrCt)     * NnNF/EqNnNF
          AdveVecCoef(2,i_v,j_v,i_m)= GRt*(upCtA-UgpCt)     * NnNF/EqNnNF
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)&
               &           = -GRt*mmA*          UgrCt*Bb    * FbNF/EqFbNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &           =  GRt*mmA*(ubA*BpCt-UgpCt*Bb)   * FbNF/EqFbNF
          j_v =  6 + vOffsetB
          AdveVecCoef(2,i_v,j_v,i_m)=  GRt*BpCt             * PpNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)=  GRt*mmA*(urCtA-UgrCt)* FtNF/EqFtNF
          AdveVecCoef(2,i_v,j_v,i_m)=  GRt*mmA*(upCtA-UgpCt)* FtNF/EqFtNF
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB      
          krCtA =  0.5D0*mmA*nnA*uuSqA*urCtA/ppA
          kpCtA =  0.5D0*mmA*nnA*uuSqA*upCtA/ppA
          AdveVecCoef(1,i_v,j_v,i_m)&
               &        = GRt*(wrCtA-krCtA-1.5D0*UgrCt)     * PpNF/EqPpNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        = GRt*(wpCtA-kpCtA-1.5D0*UgpCt)     * PpNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3 + vOffsetB
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        =  GRt*mmA*(wbA-1.5D0*ubA)*ttA*BpCt * FbNF/EqQbNF
          j_v =  6 + vOffsetB
          AdveVecCoef(2,i_v,j_v,i_m)=  2.5D0*GRt*ttA*BpCt   * PpNF/EqQbNF
          j_v =  8 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)&
               &        = -GRt*mmA          *UgrCt*Bb       * QbNF/EqQbNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        =  GRt*mmA*(ubA*BpCt-UgpCt*Bb)      * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  4 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)&
               &        =  GRt*mmA*ttA*(wrCtA-1.5D0*urCtA)  * FtNF/EqQtNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        =  GRt*mmA*ttA*(wpCtA-1.5D0*upCtA)  * FtNF/EqQtNF
          j_v =  9 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)=  GRt*mmA*(urCtA-UgrCt)* QtNF/EqQtNF
          AdveVecCoef(2,i_v,j_v,i_m)=  GRt*mmA*(upCtA-UgpCt)* QtNF/EqQtNF
       ENDDO

    CASE (2)
       i_v = 1                  ! Equation for psi'
       j_v = 1
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*G33Ct*G11Ct        * BpNF/EqBpNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*G33Ct*G12Ct        * BpNF/EqBpNF
       j_v = 3
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*G33Ct*ugrCt*VcSqRe * EtNF/EqBpNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*G33Ct*ugpCt*VcSqRe * EtNF/EqBpNF
       
       i_v = 2                   ! Equation for I
       j_v = 4
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*ugrCt*VcSqRe       * EpNF/EqBtNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*ugpCt*VcSqRe       * EpNF/EqBtNF
       
       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 1
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*UgrCt              * BpNF/EqEtNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*UgpCt              * BpNF/EqEtNF

       i_v = 4                   ! Equation for E_{\chi}
       j_v = 2
       AdveVecCoef(1,i_v,j_v,i_m) = -GRt*UgrCt              * BtNF/EqEpNF
       AdveVecCoef(2,i_v,j_v,i_m) = -GRt*UgpCt              * BtNF/EqEpNF

       i_v = 5                   ! Equation for E_{\rho}
       j_v = 4
       AdveVecCoef(1,i_v,j_v,i_m) =  GRt*G12Ct              * EpNF/EqErNF
       j_v = 5
       AdveVecCoef(1,i_v,j_v,i_m) =  GRt*G11Ct              * ErNF/EqErNF
       j_v = 6
       AdveVecCoef(2,i_v,j_v,i_m) =  GRt                    * EpNF/EqErNF
       
       i_v = 6                   ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          mmA = Mm(i_s)
          ttA = Tt(i_s)
          nnA = Nn(i_s)
          
          urCtA = UrCt(i_s)
          upCtA = UpCt(i_s)
          wrCtA = WrCt(i_s)
          wpCtA = WpCt(i_s)
          ubA   = Ub(  i_s)
          ppA   = Pp(  i_s)
          uuSqA = UuSq(i_s)
          wbA   = Wb(  i_s)
          krCtA =  0.5D0*mmA*nnA*uuSqA*urCtA/ppA
          kpCtA =  0.5D0*mmA*nnA*uuSqA*upCtA/ppA
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          j_v =  1 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)= GRt*(urCtA-UgrCt)     * NnNF/EqNnNF
          AdveVecCoef(2,i_v,j_v,i_m)= GRt*(upCtA-UgpCt)     * NnNF/EqNnNF
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)=  GRt*mmA*(urCtA-UgrCt)* FtNF/EqFrNF
          AdveVecCoef(2,i_v,j_v,i_m)=  GRt*mmA*(upCtA-UgpCt)* FtNF/EqFrNF
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  3 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)&
               &           = -GRt*mmA*          UgrCt*Bb    * FbNF/EqFbNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &           =  GRt*mmA*(ubA*BpCt-UgpCt*Bb)   * FbNF/EqFbNF
          j_v =  6 + vOffsetB
          AdveVecCoef(2,i_v,j_v,i_m)=  GRt*BpCt             * PpNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB      
          AdveVecCoef(1,i_v,j_v,i_m)&
               &        = GRt*(wrCtA-krCtA-1.5D0*UgrCt)     * PpNF/EqPpNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        = GRt*(wpCtA-kpCtA-1.5D0*UgpCt)     * PpNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  4 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)&
               &        =  GRt*mmA*ttA*(wrCtA-1.5D0*urCtA)  * FtNF/EqQrNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        =  GRt*mmA*ttA*(wpCtA-1.5D0*upCtA)  * FtNF/EqQrNF
          j_v =  9 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)=  GRt*mmA*(urCtA-UgrCt)* QtNF/EqQrNF
          AdveVecCoef(2,i_v,j_v,i_m)=  GRt*mmA*(upCtA-UgpCt)* QtNF/EqQrNF
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  3 + vOffsetB
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        =  GRt*mmA*(wbA-1.5D0*ubA)*ttA*BpCt * FbNF/EqQbNF
          j_v =  6 + vOffsetB
          AdveVecCoef(2,i_v,j_v,i_m)=  2.5D0*GRt*ttA*BpCt   * PpNF/EqQbNF
          j_v =  8 + vOffsetB
          AdveVecCoef(1,i_v,j_v,i_m)&
               &        = -GRt*mmA          *UgrCt*Bb       * QbNF/EqQbNF
          AdveVecCoef(2,i_v,j_v,i_m)&
               &        =  GRt*mmA*(ubA*BpCt-UgpCt*Bb)      * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
          
       ENDDO
       
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_AV_COEF_EB

  !---------------------------------------------------------
  ! 
  ! CALCULATION OF ADVECTION TENSOR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-22 H.Seto
  !
  !---------------------------------------------------------  
  SUBROUTINE T2COEF_AT_COEF_EB(i_m)
    
    USE T2CNST
    
    USE T2COMM, ONLY:&
         & NSMAX,NFMAX,NVMAX,NDMAX,NKMAX,&
         & FtNF,FpNF,QtNF,QpNF,&
         & EqFrNF,EqFbNF,EqFtNF,EqPpNF,EqQrNF,EqQbNF,EqQtNF,&
         & BNCXb1,BNCXb2,BNCXt1,BNCXt2,BNCPp1,BNCPp2,&
         & CNCV01,CNCV02,CNCV03,CNCV04,CNCV05,CNCV06,CNCV07,CNCV08,&
         & AdveTenCoef,EqSet    
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_d,i_k,i_v,vOffsetA,&
         &     j_d,    j_v,vOffsetB
    
    REAL(   rkind)::&
         & cncv01A,cncv02A,cncv03A,cncv04A,&
         & cncv05A,cncv06A,cncv07A,cncv08A
    
    ! INITIALIZATION
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_k = 1, NKMAX
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             AdveTenCoef(i_d,j_d,i_k,i_v,j_v,i_m) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    SELECT CASE (EqSet)
    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          cncv01A = CNCV01(i_s)
          cncv02A = CNCV02(i_s)
          cncv03A = CNCV03(i_s)
          cncv04A = CNCV04(i_s)
          cncv05A = CNCV05(i_s)
          cncv06A = CNCV06(i_s)
          cncv07A = CNCV07(i_s)
          cncv08A = CNCV08(i_s)
       
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}

          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv01A*FtNF/EqFbNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb2*cncv01A*FpNF/EqFbNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv02A*QtNF/EqFbNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb2*cncv02A*QpNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv01A*FtNF/EqFtNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv01A*FpNF/EqFtNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv02A*QtNF/EqFtNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv02A*QpNF/EqFtNF
       
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  4 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv05A*FtNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv06A*FtNF/EqPpNF
          j_v =  5 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv05A*FpNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv06A*FpNF/EqPpNF
          j_v =  9 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv07A*QtNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv08A*QtNF/EqPpNF
          j_v = 10 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv07A*QpNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv08A*QpNF/EqPpNF
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv03A*FtNF/EqQbNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb2*cncv03A*FpNF/EqQbNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv04A*QtNF/EqQbNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb1*cncv04A*QpNF/EqQbNF

          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv03A*FtNF/EqQtNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv03A*FpNF/EqQtNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv04A*QtNF/EqQtNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv04A*QpNF/EqQtNF
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a\zeta}
          
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          cncv01A = CNCV01(i_s)
          cncv02A = CNCV02(i_s)
          cncv03A = CNCV03(i_s)
          cncv04A = CNCV04(i_s)
          cncv05A = CNCV05(i_s)
          cncv06A = CNCV06(i_s)
          cncv07A = CNCV07(i_s)
          cncv08A = CNCV08(i_s)
       
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv01A*FtNF/EqFrNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv01A*FpNF/EqFrNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv02A*QtNF/EqFrNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv02A*QpNF/EqFrNF

          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv01A*FtNF/EqFbNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb2*cncv01A*FpNF/EqFbNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv02A*QtNF/EqFbNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb2*cncv02A*QpNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  4 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv05A*FtNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv06A*FtNF/EqPpNF
          j_v =  5 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv05A*FpNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv06A*FpNF/EqPpNF
          j_v =  9 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv07A*QtNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCPp1*cncv08A*QtNF/EqPpNF
          j_v = 10 + vOffsetB; i_k =  1
          AdveTenCoef(1,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv07A*QpNF/EqPpNF
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCPp2*cncv08A*QpNF/EqPpNF
       
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{rho}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv03A*FtNF/EqQrNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv03A*FpNF/EqQrNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXt1*cncv04A*QtNF/EqQrNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXt2*cncv04A*QpNF/EqQrNF
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  4 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv03A*FtNF/EqQbNF
          j_v =  5 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb2*cncv03A*FpNF/EqQbNF
          j_v =  9 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) =  BNCXb1*cncv04A*QtNF/EqQbNF
          j_v = 10 + vOffsetB; i_k = 1
          AdveTenCoef(2,2,i_k,i_v,j_v,i_m) = -BNCXb1*cncv04A*QpNF/EqQbNF

          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA    ! Equation for Q_{a\zeta}
       ENDDO
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_AT_COEF_EB
  
  !---------------------------------------------------------
  ! 
  ! CALCULATION OF DIFFUSION TENSOR COEFFICIENTS
  !
  !                 LAST UPDATE    2014-06-22 H.Seto
  !
  !---------------------------------------------------------  
  SUBROUTINE T2COEF_DT_COEF_EB(i_m)
    
    USE T2COMM,ONLY:&
         & NVMAX,NFMAX,NSMAX,NDMAX,                                 &
         & NnNF,FbNF,PpNF,QbNF,                                     &
         & EqFrNF,EqFbNF,EqFtNF,EqPpNF,EqQrNF,EqQbNF,EqQtNF,        &
         & Ub,Wb,                                                   &
         & BNCXb3,BNCXt3,BNCPp3,                                    &
         & CNCV01,CNCV02,CNCV03,CNCV04,CNCV05,CNCV06,CNCV07,CNCV08, &
         & DiffTenCoef,EqSet 
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_d,i_v,vOffsetA,&
         &     j_d,j_v,vOffsetB
    
    REAL(   rkind)::&
         ubA,wbA,&
         cncv01A,cncv02A,cncv03A,cncv04A,&
         cncv05A,cncv06A,cncv07A,cncv08A
    
    ! INITIALIZATION
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          DiffTenCoef(i_d,j_d,i_v,j_v,i_m) = 0.D0 
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    SELECT CASE (EqSet)
    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          ubA     = Ub(    i_s)
          wbA     = Wb(    i_s)
          cncv01A = CNCV01(i_s)
          cncv02A = CNCV02(i_s)
          cncv03A = CNCV03(i_s)
          cncv04A = CNCV04(i_s)
          cncv05A = CNCV05(i_s)
          cncv06A = CNCV06(i_s)
          cncv07A = CNCV07(i_s)
          cncv08A = CNCV08(i_s)
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv01A*ubA * NnNF/EqFbNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv01A     * FbNF/EqFbNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv02A*wbA * PpNF/EqFbNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv02A     * QbNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv01A*ubA * NnNF/EqFtNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv01A     * FbNF/EqFtNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv02A*wbA * PpNF/EqFtNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv02A     * QbNF/EqFtNF
       
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA    ! Equation for P_{a}
          j_v =  1 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)= -BNCPp3*cncv05A*ubA * NnNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCPp3*cncv06A*ubA * NnNF/EqPpNF
          j_v =  3 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)=  BNCPp3*cncv05A     * FbNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCPp3*cncv06A     * FbNF/EqPpNF
          j_v =  6 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)= -BNCPp3*cncv07A*wbA * PpNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCPp3*cncv08A*wbA * PpNF/EqPpNF
          j_v =  8 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)=  BNCPp3*cncv07A     * QbNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCPp3*cncv08A     * QbNF/EqPpNF

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv03A*ubA * NnNF/EqQbNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv03A     * FbNF/EqQbNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv04A*wbA * PpNF/EqQbNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv04A     * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv03A*ubA * NnNF/EqQtNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv03A     * FbNF/EqQtNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv04A*wbA * PpNF/EqQtNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv04A     * QbNF/EqQtNF
          
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          ubA     = Ub(    i_s)
          wbA     = Wb(    i_s)
          cncv01A = CNCV01(i_s)
          cncv02A = CNCV02(i_s)
          cncv03A = CNCV03(i_s)
          cncv04A = CNCV04(i_s)
          cncv05A = CNCV05(i_s)
          cncv06A = CNCV06(i_s)
          cncv07A = CNCV07(i_s)
          cncv08A = CNCV08(i_s)
          
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv01A*ubA * NnNF/EqFrNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv01A     * FbNF/EqFrNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv02A*wbA * PpNF/EqFrNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv02A     * QbNF/EqFrNF

          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv01A*ubA * NnNF/EqFbNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv01A     * FbNF/EqFbNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv02A*wbA * PpNF/EqFbNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv02A     * QbNF/EqFbNF
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}

          i_v =  6 + vOffsetA    ! Equation for P_{a}
          j_v =  1 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)= -BNCPp3*cncv05A*ubA * NnNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCPp3*cncv06A*ubA * NnNF/EqPpNF
          j_v =  3 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)=  BNCPp3*cncv05A     * FbNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCPp3*cncv06A     * FbNF/EqPpNF
          j_v =  6 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)= -BNCPp3*cncv07A*wbA * PpNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCPp3*cncv08A*wbA * PpNF/EqPpNF
          j_v =  8 + vOffsetB
          DiffTenCoef(1,2,i_v,j_v,i_m)=  BNCPp3*cncv07A     * QbNF/EqPpNF
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCPp3*cncv08A     * QbNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q^{\rho}_{a}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv03A*ubA * NnNF/EqQrNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv03A     * FbNF/EqQrNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXt3*cncv04A*wbA * PpNF/EqQrNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXt3*cncv04A     * QbNF/EqQrNF

          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv03A*ubA * NnNF/EqQbNF
          j_v =  3 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv03A     * FbNF/EqQbNF
          j_v =  6 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)= -BNCXb3*cncv04A*wbA * PpNF/EqQbNF
          j_v =  8 + vOffsetB
          DiffTenCoef(2,2,i_v,j_v,i_m)=  BNCXb3*cncv04A     * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_DT_COEF_EB

  !---------------------------------------------------------
  ! 
  !       CALCULATION OF GRADIENT VECTOR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-22 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_GV_COEF_EB(i_m)

    USE T2COMM, ONLY:&
         & UseAnomalousTransportFT, UseAnomalousTransportGT, & 
         & NSMAX,NFMAX,NVMAX,NDMAX,                          &
         & BtNF,EtNF,EpNF,ErNF,NnNF,PpNF,                    &
         & EqBpNF,EqBtNF,EqEtNF,EqEpNF,                      &
         & EqFrNF,EqFtNF,EqPpNF,EqQrNF,EqQtNF,               &
         & GRt,G11Ct,G12Ct,G22Co,G33Co,R_mc,BpCt,            &
         & Pp,Tt,UrCt,UpCt,FtAnom1,FtAnom2,GtAnom1,GtAnom2,  &
         & GradVecCoef,EqSet
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_d,i_v,vOffsetA,&
         &         j_v,vOffsetB
    
    REAL(   rkind)::&
         ttA,urCtA,upCtA,&
         ftAnom1A,ftAnom2A,gtAnom1A,gtAnom2A,gCoefrCt,gCoefpCt

    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_d = 1, NDMAX
          GradVecCoef(i_d,i_v,j_v,i_m) = 0.D0
       ENDDO
    ENDDO
    ENDDO

    SELECT CASE (EqSet)
    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v =  1                 ! Equation for psi'
       j_v =  3
       GradVecCoef(1,i_v,j_v,i_m) = -GRt                 * EtNF/EqBpNF
       
       i_v =  2                  ! Equation for I
       j_v =  4
       GradVecCoef(1,i_v,j_v,i_m) =  G33Co*R_mc          * EpNF/EqBtNF
       j_v =  5
       GradVecCoef(2,i_v,j_v,i_m) = -G33Co               * ErNF/EqBtNF
       
       i_v =  3                  ! Equation for E_{\zeta}
       
       i_v =  4                  ! Equation for E_{\chi}
       j_v =  2
       GradVecCoef(1,i_v,j_v,i_m) =  G22Co               * BtNF/EqEpNF
       
       i_v =  5                  ! Equation for E_{\rho}    
    
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       gCoefrCt = GRt*G11Ct*BpCt
       gCoefpCt = GRt*G12Ct*BpCt
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          ttA      = Tt(     i_s)
          urCtA    = UrCt(   i_s)
          upCtA    = UpCt(   i_s)
          ftAnom1A = FtAnom1(i_s)
          ftAnom2A = FtAnom2(i_s)
          GtAnom1A = GtAnom1(i_s)
          GtAnom2A = GtAnom2(i_s)
       
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          j_v =  6 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m) =  gCoefrCt         * PpNF/EqFrNF
          GradVecCoef(2,i_v,j_v,i_m) =  gCoefpCt         * PpNF/EqFrNF
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          IF(UseAnomalousTransportFT)THEN
             j_v =  1 + NFMAX
             GradVecCoef(1,i_v,j_v,i_m) = GRt*ftAnom1A   * NnNF/EqFtNF
             j_v =  6 + NFMAX
             GradVecCoef(1,i_v,j_v,i_m) = GRt*ftAnom2A   * PpNF/EqFtNF
          ENDIF
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m)= -GRt*urCtA         * PpNF/EqPpNF
          GradVecCoef(2,i_v,j_v,i_m)= -GRt*upCtA         * PpNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          j_v =  1 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m)&
               &            = -2.5D0*(ttA**2)*gCoefrCt   * NnNF/EqQrNF
          GradVecCoef(2,i_v,j_v,i_m)&
               &            = -2.5D0*(ttA**2)*gCoefpCt   * NnNF/EqQrNF
          j_v =  6 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m)=  5.0D0*ttA*gCoefrCt* PpNF/EqQrNF
          GradVecCoef(2,i_v,j_v,i_m)=  5.0D0*ttA*gCoefpCt* PpNF/EqQrNF
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          IF(UseAnomalousTransportGT)THEN
             j_v =  1 + vOffsetB 
             GradVecCoef(1,i_v,j_v,i_m) = GRt*gtAnom1A   * NnNF/EqQtNF
             j_v =  6 + vOffsetB
             GradVecCoef(1,i_v,j_v,i_m) = GRt*gtAnom2A   * PpNF/EqQtNF
          ENDIF
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
          
          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
          
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !       
       i_v =  1                  ! Equation for psi
       
       i_v =  2                  ! Equation for I
       j_v =  2
       GradVecCoef(1,i_v,j_v,i_m) =  G22Co               * BtNF/EqBtNF
       
       i_v =  3                  ! Equation for psi'
       j_v =  3
       GradVecCoef(1,i_v,j_v,i_m) = -GRt                 * EtNF/EqEtNF
       
       i_v =  4                  ! Equation for \bar{E}_{\chi}
       j_v =  4
       GradVecCoef(1,i_v,j_v,i_m) =  G33Co               * EpNF/EqEpNF
       j_v =  5
       GradVecCoef(2,i_v,j_v,i_m) = -G33Co               * ErNF/EqEpNF
       
       i_v =  5                  ! Equation for E_{\rho}
       i_v =  6                  ! Equation for \bar{E}_{\chi}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !    
       gCoefrCt = GRt*G11Ct*BpCt
       gCoefpCt = GRt*G12Ct*BpCt
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA
          
          ttA      = Tt(     i_s)
          urCtA    = UrCt(   i_s)
          upCtA    = UpCt(   i_s)
          ftAnom1A = FtAnom1(i_s)
          ftAnom2A = FtAnom2(i_s)
          GtAnom1A = GtAnom1(i_s)
          GtAnom2A = GtAnom2(i_s)
       
          i_v =  1 + vOffsetA    ! Equation for n_{a}
          
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          IF(UseAnomalousTransportFT)THEN
             j_v =  1 + NFMAX
             GradVecCoef(1,i_v,j_v,i_m) = GRt*ftAnom1A   * NnNF/EqFrNF
             j_v =  6 + NFMAX
             GradVecCoef(1,i_v,j_v,i_m) = GRt*ftAnom2A   * PpNF/EqFrNF
          ENDIF
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          j_v =  6 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m) =  gCoefrCt         * PpNF/EqFtNF
          GradVecCoef(2,i_v,j_v,i_m) =  gCoefpCt         * PpNF/EqFtNF
          
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  6 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m)= -GRt*urCtA         * PpNF/EqPpNF
          GradVecCoef(2,i_v,j_v,i_m)= -GRt*upCtA         * PpNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          ! >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
          IF(UseAnomalousTransportGT)THEN
             j_v =  1 + vOffsetB 
             GradVecCoef(1,i_v,j_v,i_m) = GRt*gtAnom1A   * NnNF/EqQrNF
             j_v =  6 + vOffsetB
             GradVecCoef(1,i_v,j_v,i_m) = GRt*gtAnom2A   * PpNF/EqQrNF
          ENDIF
          ! <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<<
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
          j_v =  1 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m)&
               &            = -2.5D0*(ttA**2)*gCoefrCt   * NnNF/EqQtNF
          GradVecCoef(2,i_v,j_v,i_m)&
               &            = -2.5D0*(ttA**2)*gCoefpCt   * NnNF/EqQtNF
          j_v =  6 + vOffsetB
          GradVecCoef(1,i_v,j_v,i_m)=  5.0D0*ttA*gCoefrCt* PpNF/EqQtNF
          GradVecCoef(2,i_v,j_v,i_m)=  5.0D0*ttA*gCoefpCt* PpNF/EqQtNF

          i_v = 10 + vOffsetA    ! Equation for Q^{\chi}_{a}
          
       END DO
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_GV_COEF_EB
  
  !---------------------------------------------------------
  ! 
  !       CALCULATION OF GRADIENT TENSOR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-22 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_GT_COEF_EB(i_m)
    
    USE T2COMM,ONLY:&
         & NSMAX,NVMAX,NFMAX,NKMAX,NDMAX,         &
         & NnNF,FbNF,PpNF,QbNF,                   &
         & EqFbNF,EqPpNF,EqQbNF,            &
         & Ub,Wb,&
         & BNCXb4,BNCPp4,BNCPp7,&
         & CNCV01,CNCV02,CNCV03,CNCV04,CNCV09,CNCV10,&
         & GradTenCoef,EqSet

    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_d,i_k,i_v,vOffsetA,kOffsetX,&
         &     j_d,    j_v,vOffsetB
    REAL(   rkind)::&
         ubA,wbA,cncv01A,cncv02A,cncv03A,cncv04A,cncv09A,cncv10A
         
    
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_k = 1, NKMAX
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             GradTenCoef(i_d,j_d,i_k,i_v,j_v,i_m) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    SELECT CASE (EqSet)

    CASE (1)
       !
       ! variables as field (from i_v= 1 to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
       
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA 
          kOffsetX =  2*(i_s-1)
          
          ubA     = Ub(  i_s)
          wbA     = Wb(  i_s)
          cncv01A = CNCV01(i_s)
          cncv02A = CNCV02(i_s)
          cncv03A = CNCV03(i_s)
          cncv04A = CNCV04(i_s)
          cncv09A = CNCV09(i_s)
          cncv10A = CNCV10(i_s)

          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv01A*ubA * NnNF/EqFbNF
          j_v =  3 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv01A     * FbNF/EqFbNF
          j_v =  6 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv02A*wbA * PpNF/EqFbNF
          j_v =  8 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv02A     * QbNF/EqFbNF
       
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  1 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp4*cncv09A*ubA * NnNF/EqPpNF
          j_v =  1 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp7*cncv01A*ubA * NnNF/EqPpNF
          j_v =  3 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp4*cncv09A     * FbNF/EqPpNF
          j_v =  3 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp7*cncv01A     * FbNF/EqPpNF
          j_v =  6 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp4*cncv10A*wbA * PpNF/EqPpNF
          j_v =  6 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp7*cncv02A*wbA * PpNF/EqPpNF
          j_v =  8 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp4*cncv10A     * QbNF/EqPpNF
          j_v =  8 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp7*cncv02A     * QbNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv03A*ubA * NnNF/EqQbNF
          j_v =  3 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv03A     * FbNF/EqQbNF
          j_v =  6 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv04A*wbA * PpNF/EqQbNF
          j_v =  8 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv04A     * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}      
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    CASE (2)
       !
       ! variables as field (from i_v= 1        to i_v = NFMAX)
       !
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
       i_v = 6                   ! Equation for E^{\chi}
       
       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1) + NFMAX
          vOffsetB = vOffsetA 
          kOffsetX =  2*(i_s-1)
          
          ubA     = Ub(  i_s)
          wbA     = Wb(  i_s)
          cncv01A = CNCV01(i_s)
          cncv02A = CNCV02(i_s)
          cncv03A = CNCV03(i_s)
          cncv04A = CNCV04(i_s)
          cncv09A = CNCV09(i_s)
          cncv10A = CNCV10(i_s)

          i_v =  1 + vOffsetA    ! Equation for n_{a}
          i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
          
          i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
          j_v =  1 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv01A*ubA * NnNF/EqFbNF
          j_v =  3 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv01A     * FbNF/EqFbNF
          j_v =  6 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv02A*wbA * PpNF/EqFbNF
          j_v =  8 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv02A     * QbNF/EqFbNF
       
          i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
          i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
          
          i_v =  6 + vOffsetA    ! Equation for p_{a}
          j_v =  1 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp4*cncv09A*ubA * NnNF/EqPpNF
          j_v =  1 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp7*cncv01A*ubA * NnNF/EqPpNF
          j_v =  3 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp4*cncv09A     * FbNF/EqPpNF
          j_v =  3 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp7*cncv01A     * FbNF/EqPpNF
          j_v =  6 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp4*cncv10A*wbA * PpNF/EqPpNF
          j_v =  6 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCPp7*cncv02A*wbA * PpNF/EqPpNF
          j_v =  8 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp4*cncv10A     * QbNF/EqPpNF
          j_v =  8 + vOffsetB; i_k =  3 + kOffsetX
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCPp7*cncv02A     * QbNF/EqPpNF
          
          i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
          
          i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
          j_v =  1 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv03A*ubA * NnNF/EqQbNF
          j_v =  3 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv03A     * FbNF/EqQbNF
          j_v =  6 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   = -BNCXb4*cncv04A*wbA * PpNF/EqQbNF
          j_v =  8 + vOffsetB; i_k =  1
          GradTenCoef(2,2,i_k,i_v,j_v,i_m)&
               &                   =  BNCXb4*cncv04A     * QbNF/EqQbNF
          
          i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}      
          i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
          
       ENDDO
    END SELECT

    RETURN
    
  END SUBROUTINE T2COEF_GT_COEF_EB

  !---------------------------------------------------------
  ! 
  !       CALCULATION OF EXCITATION SCALAR COEFFICIENTS
  !
  !                 LAST UPDATE    2014-06-22 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_ES_COEF_EB(i_m)
    
    USE T2CNST,ONLY:&
         Eps0,Rmu0
    USE T2COMM,ONLY:&
         & UseAnomalousTransportFT,UseAnomalousTransportGT,&
         & NSMAX,NVMAX,& 
         &           EtNF,EpNF,ErNF,&
         & NnNF,FrNF,FbNF,FtNF,FpNF,&
         & PpNF,QrNF,QbNF,QtNF,QpNF,&
         & EqBtNF,EqEtNF,EqEpNF,EqErNF,&
         &        EqFrNF,EqFbNF,EqFtNF,EqFpNF,&
         & EqPpNF,EqQrNF,EqQbNF,EqQtNF,EqQpNF,&
         & GRt,G12Ct,G11Ct,G12Co,G22Co,G33Co, &
         & R_mc,BpCt,BtCo,BtCt,Bb,BbSq,       &
         & Ee,Nn,Pp,Tt,Hex,CNCF01,CNCF02,CNCF03,CNCF04,&
         & FtAnom3,FtAnom4,GtAnom3,GtAnom4,&
         & ExciScaCoef,EqSet
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_v,vOffsetA,&
         & j_s,j_v,vOffsetB
    
    REAL(   rkind)::&
         & eeA,nnA,ppA,ttA,hexA,nnB,ppB,eeB,&
         & ftAnom3A,ftAnom4A,gtAnom3A,gtAnom4A,&
         & gCoef11,gCoef12,&
         & cncf01AB,cncf02AB,cncf03AB,cncf04AB
    
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       ExciScaCoef(i_v,j_v,i_m) = 0.D0
    ENDDO
    ENDDO
    
    !
    ! variables as field (from i_v= 1 to i_v = NFMAX)
    !
    SELECT CASE (EqSet)
    CASE (1)
       
       i_v = 1                   ! Equation for psi'
       
       i_v = 2                   ! Equation for I
       j_v = 4
       ExciScaCoef(i_v,j_v,i_m) =  G33Co                 * EpNF/EqBtNF
       
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}

       !
       ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
       !
       gCoef12 = GRt*G12Ct*BpCt*R_mc
       gCoef11 = GRt*G11Ct*BpCt
       DO j_s = 1, NSMAX
          
          vOffsetB = 10*(j_s-1) +NFMAX
          nnB = Nn(j_s)
          ppB = Pp(j_s)
          
          DO i_s = 1, NSMAX
             
             vOffsetA = 10*(i_s-1) + NFMAX
             cncf01AB = CNCF01(i_s,j_s) 
             cncf02AB = CNCF02(i_s,j_s) 
             cncf03AB = CNCF03(i_s,j_s) 
             cncf04AB = CNCF04(i_s,j_s)
             
             IF(i_s.EQ.j_s)THEN
                
                nnA      = Nn( i_s)
                ppA      = Pp( i_s)
                eeA      = Ee( i_s)
                hexA     = Hex(i_s)
                ftAnom3A = FtAnom3(i_s)
                ftAnom4A = FtAnom4(i_s)
                gtAnom3A = GtAnom3(i_s)
                gtAnom4A = GtAnom4(i_s)
                         
                i_v =  1 + vOffsetA    ! Equation for n_{a}
                
                i_v =  2 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
                j_v =  4
                ExciScaCoef(i_v,j_v,i_m)= -gCoef12*eeA*nnA*EpNF/EqFrNF
                j_v =  5
                ExciScaCoef(i_v,j_v,i_m)= -gCoef11*eeA*nnA*ErNF/EqFrNF
                j_v =  3 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -BtCo*Bb*eeA    *FbNF/EqFrNF
                j_v =  4 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  BbSq   *eeA    *FtNF/EqFrNF
                
                i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
                j_v =  3
                ExciScaCoef(i_v,j_v,i_m)&
                     &               = -GRt*BtCt*eeA*nnA * EtNF/EqFbNF
                j_v =  4
                ExciScaCoef(i_v,j_v,i_m)&
                     &          = -GRt*BpCt*eeA*nnA*R_mc * EpNF/EqFbNF
                j_v =  3 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf01AB*Bb   * FbNF/EqFbNF
                j_v =  8 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf02AB*Bb   * QbNF/EqFbNF
                
                i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
                j_v =  3
                ExciScaCoef(i_v,j_v,i_m)= -GRt*eeA*nnA   * EtNF/EqFtNF
                j_v =  2 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)&
                     &          = -GRt*GRt*BpCt*R_mc*eeA * FrNF/EqFtNF
                j_v =  4 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                     &                    -cncf01AB      * FtNF/EqFtNF
                j_v =  9 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf02AB      * QtNF/EqFtNF
                !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
                IF(UseAnomalousTransportFT)THEN
                   j_v =  3 + NFMAX
                   ExciScaCoef(i_v,j_v,i_m)&
                        &           =  GRt*ftAnom3A      * FbNF/EqFtNF
                   j_v =  4 + NFMAX
                   ExciScaCoef(i_v,j_v,i_m)&
                        &           =  ExciScaCoef(i_v,j_v,i_m) &
                        &             +GRt*ftAnom4A      * FtNF/EqFtNF
                ENDIF
                !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
             
                i_v =  5 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
                j_v =  3 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -GRt*Bb        * FbNF/EqFpNF
                j_v =  4 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  GRt*BtCt      * FtNF/EqFpNF
                j_v =  5 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  GRt*G22Co*BpCt* FpNF/EqFpNF
             
                i_v =  6 + vOffsetA    ! Equation for p_{a}
                j_v =  6 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  GRt*hexA      * PpNF/EqPpNF
                         
                i_v =  7 + vOffsetA    ! Equation for Q_{a}^{\rho}
                j_v =  4
                ExciScaCoef(i_v,j_v,i_m)&
                     &          = -2.5D0*gCoef12*eeA*ppA * EpNF/EqQrNF
                j_v =  5
                ExciScaCoef(i_v,j_v,i_m)&
                     &          = -2.5D0*gCoef11*eeA*ppA * ErNF/EqQrNF
                j_v =  8 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)&
                     &             = -BtCo*Bb*eeA        * QbNF/EqQrNF
                j_v =  9 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)&
                     &             =  BbSq*eeA           * QtNF/EqQrNF
                
                i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
                j_v =  3
                ExciScaCoef(i_v,j_v,i_m)&
                  &       = -2.5D0*GRt*BtCt     *eeA*ppA * EtNF/EqQbNF
                j_v =  4
                ExciScaCoef(i_v,j_v,i_m)&
                     &    = -2.5D0*GRt*BpCt*R_mc*eeA*ppA * EpNF/EqQbNF
                j_v =  3 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf03AB*Bb   * FbNF/EqQbNF
                j_v =  8 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf04AB*Bb   * QbNF/EqQbNF
                
                i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
                j_v =  3
                ExciScaCoef(i_v,j_v,i_m)&
                     &          = -2.5D0*GRt*eeA*ppA     * EtNF/EqQtNF
                j_v =  7 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)&
                     &          = -GRt*GRt*BpCt*R_mc*eeA * QrNF/EqQtNF
                !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
                IF(UseAnomalousTransportGT)THEN
                   j_v =  8 + vOffsetB
                   ExciScaCoef(i_v,j_v,i_m)&
                        &             =  GRt*gtAnom3A    * QbNF/EqQtNF
                   j_v =  9 + vOffsetB
                   ExciScaCoef(i_v,j_v,i_m)&
                        &             =  ExciScaCoef(i_v,j_v,i_m) &
                        &               +GRt*gtAnom4A    * QtNF/EqQtNF
                ENDIF
                !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
                j_v =  4 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf03AB      * FtNF/EqQtNF
                
                j_v =  9 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                     &                    -cncf04AB      * QtNF/EqQtNF
                
                i_v = 10 + vOffsetA    ! Equation for Q_{a}^{\chi}
                j_v =  8 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -GRt*Bb        * QbNF/EqQpNF
                j_v =  9 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  GRt*BtCt      * QtNF/EqQpNF
                j_v = 10 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  GRt*G22Co*BpCt* QpNF/EqQpNF
             ELSE
                i_v =  3 + vOffsetA    ! Equation for Gamma_{a\para}
                j_v =  3 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf01AB*Bb   * FbNF/EqFbNF
                j_v =  8 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf02AB*Bb   * QbNF/EqFbNF
                
                i_v =  4 + vOffsetA    ! Equation for Gamma_{a\zeta}
                j_v =  4 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                     &                    -cncf01AB      * FtNF/EqFtNF
                j_v =  9 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf02AB      * QtNF/EqFtNF
             
                i_v =  8 + vOffsetA    ! Equation for Q_{a\para}
                j_v =  3 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf03AB*Bb   * FbNF/EqQbNF
                j_v =  8 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf04AB*Bb   * QbNF/EqQbNF
                
                i_v =  9 + vOffsetA    ! Equation for Q_{a\zeta}
                j_v =  4 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)= -cncf03AB      * FtNF/EqQtNF
                j_v =  9 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                     &                    -cncf04AB      * QtNF/EqQtNF
             ENDIF
          ENDDO
       ENDDO
    CASE (2)
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for \bar{E}_{\chi}
       i_v = 5                   ! Equation for E_{\rho}

       i_v = 6                   ! Equation for \bar{E}_{\chi}
       ! KOKOMADE
    END SELECT
    
    !
    ! variables as fluid (from i_v = NFMAX+1 to i_v = NVMAX)
    !
    gCoef12 = GRt*G12Ct*BpCt*R_mc
    gCoef11 = GRt*G11Ct*BpCt
    DO j_s = 1, NSMAX
       
       vOffsetB = 10*(j_s-1)
       nnB = Nn(j_s)
       ppB = Pp(j_s)
       
       DO i_s = 1, NSMAX
          
          vOffsetA = 10*(i_s-1)
          cncf01AB = CNCF01(i_s,j_s) 
          cncf02AB = CNCF02(i_s,j_s) 
          cncf03AB = CNCF03(i_s,j_s) 
          cncf04AB = CNCF04(i_s,j_s)
          
          IF(i_s.EQ.j_s)THEN
             
             nnA      = Nn( i_s)
             ppA      = Pp( i_s)
             eeA      = Ee( i_s)
             hexA     = Hex(i_s)
             ftAnom3A = FtAnom3(i_s)
             ftAnom4A = FtAnom4(i_s)
             gtAnom3A = GtAnom3(i_s)
             gtAnom4A = GtAnom4(i_s)
             
             
             i_v =  6 + vOffsetA    ! Equation for n_{a}
             
             i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
             j_v =  4
             ExciScaCoef(i_v,j_v,i_m) = -gCoef12*eeA*nnA * EpNF/EqFrNF
             j_v =  5
             ExciScaCoef(i_v,j_v,i_m) = -gCoef11*eeA*nnA * ErNF/EqFrNF
             j_v =  8 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m) = -BtCo*Bb*eeA     * FbNF/EqFrNF
             j_v =  9 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m) =  BbSq   *eeA     * FtNF/EqFrNF
             
             i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
             j_v =  3
             ExciScaCoef(i_v,j_v,i_m)= -GRt*BtCt*eeA*nnA * EtNF/EqFbNF
             j_v =  4
             ExciScaCoef(i_v,j_v,i_m)&
                  &             = -GRt*BpCt*eeA*nnA*R_mc * EpNF/EqFbNF
             j_v =  8 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf01AB*Bb      * FbNF/EqFbNF
             j_v = 13 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf02AB*Bb  * QbNF/EqFbNF
             
             i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
             !j_v =  3
             ExciScaCoef(i_v,j_v,i_m)= -GRt*eeA*nnA      * EtNF/EqFtNF
             j_v =  7 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)&
                  &             = -GRt*GRt*BpCt*R_mc*eeA * FrNF/EqFtNF
             j_v =  9 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                  &                    -cncf01AB         * FtNF/EqFtNF
             j_v = 14 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf02AB         * QtNF/EqFtNF
             !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
             IF(UseAnomalousTransportFT)THEN
                j_v =  8
                ExciScaCoef(i_v,j_v,i_m)&
                     &             =  GRt*ftAnom3A       * FbNF/EqFtNF
                j_v =  9
                ExciScaCoef(i_v,j_v,i_m)&
                     &             =  ExciScaCoef(i_v,j_v,i_m) &
                     &               +GRt*ftAnom4A       * FtNF/EqFtNF
             ENDIF
             !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
             
             i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
             j_v =  8 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m) = -GRt*Bb          * FbNF/EqFpNF
             j_v =  9 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m) =  GRt*BtCt        * FtNF/EqFpNF
             j_v = 10 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m) =  GRt*G22Co*BpCt  * FpNF/EqFpNF
             
             i_v = 11 + vOffsetA    ! Equation for p_{a}
             j_v = 11 + vOffsetA
             ExciScaCoef(i_v,j_v,i_m) =  GRt*hexA        * PpNF/EqPpNF
                         
             i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}
             j_v =  4
             ExciScaCoef(i_v,j_v,i_m)&
                  &             = -2.5D0*gCoef12*eeA*ppA * EpNF/EqQrNF
             j_v =  5
             ExciScaCoef(i_v,j_v,i_m)&
                  &             = -2.5D0*gCoef11*eeA*ppA * ErNF/EqQrNF
             j_v = 13 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)&
                  &             = -BtCo*Bb*eeA           * QbNF/EqQrNF
             j_v = 14 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)&
                  &             =  BbSq*eeA              * QtNF/EqQrNF
             
             i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
             j_v =  3
             ExciScaCoef(i_v,j_v,i_m)&
                  &       = -2.5D0*GRt*BtCt     *eeA*ppA * EtNF/EqQbNF
             j_v =  4
             ExciScaCoef(i_v,j_v,i_m)&
                  &       = -2.5D0*GRt*BpCt*R_mc*eeA*ppA * EpNF/EqQbNF
             j_v =  8 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf03AB*Bb      * FbNF/EqQbNF
             j_v = 13 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf04AB*Bb      * QbNF/EqQbNF
             
             i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
             j_v =  3
             ExciScaCoef(i_v,j_v,i_m)= -2.5D0*GRt*eeA*ppA* EtNF/EqQtNF
             j_v = 12 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)&
                  &              = -GRt*GRt*BpCt*R_mc*eeA* QrNF/EqQtNF
             !  >>>> ANOMALOUS TRANSPORT * two-fluid model >>>>
             IF(UseAnomalousTransportGT)THEN
                j_v = 13 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)&
                     &             =  GRt*gtAnom3A       * QbNF/EqQtNF
                j_v = 14 + vOffsetB
                ExciScaCoef(i_v,j_v,i_m)&
                     &             =  ExciScaCoef(i_v,j_v,i_m) &
                     &               +GRt*gtAnom4A       * QtNF/EqQtNF
             ENDIF
             !  <<<< ANOMALOUS TRANSPORT * two-fluid model <<<<
             j_v =  9 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf03AB         * FtNF/EqQtNF
             
             j_v = 14 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                  &                    -cncf04AB         * QtNF/EqQtNF
       
             i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
             j_v = 13 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -GRt*Bb           * QbNF/EqQpNF
             j_v = 14 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)=  GRt*BtCt         * QtNF/EqQpNF
             j_v = 15 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)=  GRt*G22Co*BpCt   * QpNF/EqQpNF
          ELSE
             i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
             j_v =  8 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf01AB*Bb      * FbNF/EqFbNF
             j_v = 13 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf02AB*Bb      * QbNF/EqFbNF
             
             i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
             j_v =  9 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                  &                    -cncf01AB         * FtNF/EqFtNF
             j_v = 14 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf02AB         * QtNF/EqFtNF

             i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
             j_v =  8 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf03AB*Bb      * FbNF/EqQbNF
             j_v = 13 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf04AB*Bb      * QbNF/EqQbNF
             
             i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
             j_v =  9 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)= -cncf03AB         * FtNF/EqQtNF
             
             j_v = 14 + vOffsetB
             ExciScaCoef(i_v,j_v,i_m)=  ExciScaCoef(i_v,j_v,i_m) &
                  &                    -cncf04AB         * QtNF/EqQtNF
          ENDIF
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2COEF_ES_COEF_EB

  !---------------------------------------------------------
  ! 
  !       CALCULATION OF ECITATION VECTOR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-03 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_EV_COEF_EB(i_m)
    
    USE T2COMM,ONLY:&
         & NSMAX,NVMAX,NDMAX,NKMAX,                                &
         & BpNF,EtNF,EpNF,FbNF,QbNF,                               &
         & EqBpNF,EqEtNF,EqFbNF,EqQbNF,EqQtNF,                     & 
         & GRt,BpCt,Bb,G11Ct,G12Ct,R_mc,R_rz,Mm,Tt,Ub,Wb,          &
         & BNCQb1,BNCQb2,BNCQb3,BNCQb4,BNCQt1,BNCQt2,BNCQt3,BNCQt4,&
         & CNCV11,CNCV12,CNCV13,                                   &
         & ExciVecCoef,EqSet
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_v,i_k,i_d,vOffsetA,kOffsetX,&
         &     j_v,        vOffsetB   
    
    REAL(   rkind)::&
         mmA,ttA,ubA,wbA,b01,u01A,cncv11A,cncv12A,cncv13A
    
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_k = 1, NKMAX 
          DO i_d = 1, NDMAX
             ExciVecCoef(i_d,i_k,i_v,j_v,i_m) = 0.D0
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    SELECT CASE (EqSet)
    CASE (1)
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I

       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 1; i_k = 2
       ExciVecCoef(1,i_k,i_v,j_v,i_m)=  2.D0*GRt*G11Ct/R_rz * BpNF/EqEtNF
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  2.D0*GRt*G12Ct/R_rz * BpNF/EqEtNF
       
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
    CASE (2)
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}
       i_v = 5                   ! Equation for E_{\rho}
    END SELECT
        
    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !
    b01 = GRt*BpCt/Bb
    
    DO i_s = 1, NSMAX
       
       vOffsetA = 10*(i_s-1)
       vOffsetB = vOffsetA
       kOffsetX =  2*(i_s-1)
       
       mmA     = Mm(    i_s)
       ttA     = Tt(    i_s)
       ubA     = Ub(    i_s)
       wbA     = Wb(    i_s)
       cncv11A = CNCV11(i_s)
       cncv12A = CNCV12(i_s)
       cncv13A = CNCV13(i_s)
       
       u01A = wbA -1.5D0*ubA      
       
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}
       
       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  8 + vOffsetB; i_k =  1
       ExciVecCoef(2,i_k,i_v,j_v,i_m)= -b01*mmA*ubA      * FbNF/EqFbNF

       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       i_v = 12 + vOffsetA    ! Equation for Q_{a}^{\rho}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  3; i_k =  1
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQb1*cncv11A   * EtNF/EqQbNF
       j_v =  3; i_k =  3+kOffsetX
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQb2*cncv12A   * EtNF/EqQbNF
       j_v =  3; i_k =  4+kOffsetX
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQb2*cncv13A   * EtNF/EqQbNF
       j_v =  4; i_k =  1
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQb3*cncv11A   * EpNF/EqQbNF
       j_v =  4; i_k =  3+kOffsetX
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQb4*cncv12A   * EpNF/EqQbNF
       j_v =  4; i_k =  4+kOffsetX
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQb4*cncv13A   * EtNF/EqQbNF
       
       j_v =  8 + vOffsetB; i_k =  1
       ExciVecCoef(2,i_k,i_v,j_v,i_m)= -b01*mmA*u01A*ttA * FbNF/EqQbNF
       j_v = 13 + vOffsetB; i_k =  1
       ExciVecCoef(2,i_k,i_v,j_v,i_m)= -b01*mmA*ubA      * QbNF/EqQbNF
       
       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       j_v =  3; i_k =  1
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQt1*cncv11A   * EtNF/EqQtNF
       j_v =  3; i_k =  3+kOffsetX;
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQt2*cncv12A   * EtNF/EqQtNF
       j_v =  3; i_k =  4+kOffsetX;
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQt2*cncv13A   * EtNF/EqQtNF
       j_v =  4; i_k =  1;
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQt3*cncv11A   * EpNF/EqQtNF
       j_v =  4; i_k =  3+kOffsetX;
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQt4*cncv12A   * EpNF/EqQtNF
       j_v =  4; i_k =  4+kOffsetX;
       ExciVecCoef(2,i_k,i_v,j_v,i_m)=  BNCQt4*cncv13A   * EpNF/EqQtNF

       i_v = 15 + vOffsetA    ! Equation for Q_{a}^{\chi}
   
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2COEF_EV_COEF_EB
  
  !---------------------------------------------------------
  ! 
  !       CALCULATION OF ECITATION VECTOR COEFFICIENTS
  !
  !                 LAST UPDATE 2014-06-11 H.Seto
  !
  !---------------------------------------------------------
  SUBROUTINE T2COEF_ET_COEF_EB(i_m)
    

    USE T2COMM,ONLY:&
         & NSMAX,NDMAX,NVMAX,NKMAX,                   &
         & FtNF,FpNF,QtNF,QpNF,                       &
         & EqFbNF,EqPpNF,EqQbNF,                      &
         & BNCXb5,BNCXb6,BNCPp5,BNCPp6,BNCPp8,BNCPp9, &
         & CNCV01,CNCV02,CNCV03,CNCV04,CNCV09,CNCV10, &
         & ExciTenCoef
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         & i_s,i_d,i_k,i_v,vOffsetA,kOffsetX,&
         &     j_d,j_k,j_v,vOffsetB,kOffsetY
    
    REAL(   rkind)::&
         & cncv01A,cncv02A,cncv03A,cncv04A,cncv09A,cncv10A
        
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO j_k = 1, NKMAX
       DO i_k = 1, NKMAX
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             ExciTenCoef(i_d,j_d,i_k,j_k,i_v,j_v,i_m) = 0.D0
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

    !
    ! variables as field (from i_v= 1 to i_v = 5)
    !
    i_v = 1                   ! Equation for psi'
    i_v = 2                   ! Equation for I
    i_v = 3                   ! Equation for E_{\zeta}
    i_v = 4                   ! Equation for E_{\chi}
    i_v = 5                   ! Equation for E_{\rho}
    
    !
    ! variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
    !
    DO i_s = 1, NSMAX
       
       vOffsetA = 10*(i_s-1)
       vOffsetB = vOffsetA
       kOffsetX =  2*(i_s-1)
       kOffsetY = kOffsetX
       
       cncv01A = CNCV01(i_s)
       cncv02A = CNCV02(i_s)
       cncv03A = CNCV03(i_s)
       cncv04A = CNCV04(i_s)
       cncv09A = CNCV09(i_s)
       cncv10A = CNCV10(i_s)
 
       i_v =  6 + vOffsetA    ! Equation for n_{a}
       i_v =  7 + vOffsetA    ! Equation for Gamma_{a}^{\rho}

       i_v =  8 + vOffsetA    ! Equation for Gamma_{a\para}
       j_v =  9 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                          =  BNCXb5*cncv01A * FtNF/EqFbNF
       j_v = 10 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                          = -BNCXb6*cncv01A * FpNF/EqFbNF
       j_v = 14 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                          =  BNCXb5*cncv02A * QtNF/EqFbNF
       j_v = 15 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                          = -BNCXb6*cncv02A * QpNF/EqFbNF
       
       i_v =  9 + vOffsetA    ! Equation for Gamma_{a\zeta}
       i_v = 10 + vOffsetA    ! Equation for Gamma_{a}^{\chi}
       
       i_v = 11 + vOffsetA    ! Equation for p_{a}
       j_v =  9 + vOffsetB; i_k = 1           ; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         = -BNCPp5*cncv09A * FtNF/EqPpNF
       j_v = 10 + vOffsetB; i_k = 1           ; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         =  BNCPp6*cncv09A * FpNF/EqPpNF
       j_v = 14 + vOffsetB; i_k = 1           ; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         = -BNCPp5*cncv10A * QtNF/EqPpNF
       j_v = 15 + vOffsetB; i_k = 1           ; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         =  BNCPp6*cncv10A * QpNF/EqPpNF
       
       j_v =  9 + vOffsetB; i_k = 3 + kOffsetX; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         = -BNCPp8*cncv01A * FtNF/EqPpNF
       j_v = 10 + vOffsetB; i_k = 3 + kOffsetX; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         =  BNCPp9*cncv01A * FpNF/EqPpNF
       j_v = 14 + vOffsetB; i_k = 3 + kOffsetX; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         = -BNCPp8*cncv02A * QtNF/EqPpNF
       j_v = 15 + vOffsetB; i_k = 3 + kOffsetX; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         =  BNCPp9*cncv02A * QpNF/EqPpNF
       
       i_v = 12 + vOffsetA    ! Equation for \bar{Q}^{\rho}_{a}

       i_v = 13 + vOffsetA    ! Equation for Q_{a\para}
       j_v =  9 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         =  BNCXb5*cncv03A * FtNF/EqQbNF
       j_v = 10 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         = -BNCXb6*cncv03A * FpNF/EqQbNF
       j_v = 14 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         =  BNCXb5*cncv04A * QtNF/EqQbNF
       j_v = 15 + vOffsetB; i_k = 1; j_k = 1
       ExciTenCoef(2,2,i_k,j_k,i_v,j_v,i_m)&
            &                         = -BNCXb6*cncv04A * QpNF/EqQbNF

       i_v = 14 + vOffsetA    ! Equation for Q_{a\zeta}
       i_v = 15 + vOffsetA    ! Equation for Q^{\chi}_{a}

    ENDDO
    
    RETURN
    
  END SUBROUTINE T2COEF_ET_COEF_EB

  SUBROUTINE T2COEF_SS_COEF_EB(i_m)
    
    USE T2CNST, ONLY: Rmu0,EPS0
    USE T2COMM, ONLY: NVMAX,NSMAX,GRt,Ee,Nn,FpCo,FtCo,FtCt,SourScaCoef,&
         &            EqBpNF,EqBtNF,EqEtNF,EqEpNF,EqErNF,EqSet
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v,j_s
    REAL(   rkind)::eeB,ftCoB,ftCtB,fpCoB,nnB

    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       SourScaCoef(i_v,j_v,i_m) = 0.D0
    ENDDO
    ENDDO

    SELECT CASE (EqSet)
    CASE (1)
       i_v = 1                   ! Equation for psi'
       i_v = 2                   ! Equation for I
       
       i_v = 3                   ! Equation for E_{\zeta}
       j_v = 3
       DO j_s = 1, NSMAX
          eeB   = Ee(  j_s)
          ftCoB = FtCo(j_s)
          SourScaCoef(i_v,j_v,i_m) =  SourScaCoef(i_v,j_v,i_m)&
               &                     -GRt*Rmu0*eeB*ftCoB / EqEtNF
       ENDDO
       
       i_v = 4                   ! Equation for E_{\chi}
       j_v = 4
       DO j_s = 1, NSMAX
          eeB   = Ee(  j_s)
          fpCoB = FpCo(j_s)
          SourScaCoef(i_v,j_v,i_m) =  SourScaCoef(i_v,j_v,i_m)&
               &                     -GRt*Rmu0*eeB*fpCoB / EqEpNF
       ENDDO
       
       i_v = 5                   ! Equation for E_{\rho}
       j_v = 5
       DO j_s = 1, NSMAX
          eeB = Ee(j_s)
          nnB = Nn(j_s)
          SourScaCoef(i_v,j_v,i_m) =  SourScaCoef(i_v,j_v,i_m)&
               &                     +GRt*eeB*nnB/Eps0   / EqErNF
       ENDDO
    CASE (2)
       i_v = 1                   ! Equation for psi'
       j_v = 1
       DO j_s = 1, NSMAX
          eeB   = Ee(  j_s)
          ftCtB = FtCt(j_s)
          SourScaCoef(i_v,j_v,i_m) =  SourScaCoef(i_v,j_v,i_m)&
               &                     -GRt*Rmu0*eeB*ftCtB / EqBpNF
       ENDDO

       i_v = 2                   ! Equation for I
       j_v = 2
       DO j_s = 1, NSMAX
          eeB   = Ee(  j_s)
          fpCoB = FpCo(j_s)
          SourScaCoef(i_v,j_v,i_m) =  SourScaCoef(i_v,j_v,i_m)&
               &                     -GRt*Rmu0*eeB*fpCoB / EqBtNF
       ENDDO
    
       i_v = 3                   ! Equation for E_{\zeta}
       i_v = 4                   ! Equation for E_{\chi}

       i_v = 5                   ! Equation for E_{\rho}
       j_v = 5
       DO j_s = 1, NSMAX
          eeB = Ee(j_s)
          nnB = Nn(j_s)
          SourScaCoef(i_v,j_v,i_m) =  SourScaCoef(i_v,j_v,i_m)&
               &                     +GRt*eeB*nnB/Eps0   / EqErNF
       ENDDO
    END SELECT
    
    RETURN
    
  END SUBROUTINE T2COEF_SS_COEF_EB

  SUBROUTINE T2COEF_MS_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_MS_COEF_PhiA

  SUBROUTINE T2COEF_AV_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_AV_COEF_PhiA
  
  SUBROUTINE T2COEF_AT_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_AT_COEF_PhiA
  
  SUBROUTINE T2COEF_DT_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_DT_COEF_PhiA
  
  SUBROUTINE T2COEF_GV_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_GV_COEF_PhiA
  SUBROUTINE T2COEF_GT_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_GT_COEF_PhiA

  SUBROUTINE  T2COEF_ES_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_ES_COEF_PhiA

  SUBROUTINE T2COEF_EV_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_EV_COEF_PhiA

  SUBROUTINE T2COEF_ET_COEF_PhiA(i_m)
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_ET_COEF_PhiA

  SUBROUTINE  T2COEF_SS_COEF_PhiA(i_m)  
    INTEGER(ikind),INTENT(IN)::i_m
    RETURN
  END SUBROUTINE T2COEF_SS_COEF_PhiA

  SUBROUTINE T2COEF_CHECK(i_m)
    
    USE T2COMM, ONLY:&
         NVMAX,NKMAX,NDMAX,KnownVar,&
         MassScaCoef,AdveVecCoef,AdveTenCoef,DiffTenCoef,GradVecCoef,&
         GradTenCoef,ExciScaCoef,ExciVecCoef,ExciTenCoef,SourScaCoef
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v,i_k,j_k,i_d,j_d
    OPEN(31,FILE="TEST_T2COEF.txt")
    
100 FORMAT(A15,1X,'ik=',I2,1X,'jk=',I2,1X,'iv=',I2,1X,'jv=',I2,1X,'id=',I2,1X,'jd=',I2,1X,'val=',D15.8)

    DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       !DO i_v = 1, NVMAX
       !DO j_v = 1, NVMAX
          !DO i_d = 1, NDNAX
          !DO j_d = 1, NDMAX
             WRITE(31,100)"KnownVar",i_k,0,0,0,0,0,KnownVar(i_k,i_m)
          !ENDDO
          !ENDDO
       !ENDDO
       !ENDDO
    !ENDDO
    ENDDO

    !DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          !DO i_d = 1, NDNAX
          !DO j_d = 1, NDMAX
             WRITE(31,100)"MassScaCoef",0,0,i_v,j_v,0,0,MassScaCoef(i_v,j_v,i_m)
          !ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    !DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          DO i_d = 1,NDMAX
          !DO j_d = 1,NDMAX
             WRITE(31,100)"AdveVecCoef",0,0,i_v,j_v,i_d,0,AdveVecCoef(i_d,i_v,j_v,i_m)
          ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    !DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          !DO i_d = 1, NDNAX
          !DO j_d = 1, NDMAX
             WRITE(31,100)"MassScaCoef",0,0,i_v,j_v,0,0,MassScaCoef(i_v,j_v,i_m)
          !ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    !DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          DO i_d = 1,NDMAX
          !DO j_d = 1,NDMAX
             WRITE(31,100)"AdveVecCoef",0,0,i_v,j_v,i_d,0,AdveVecCoef(i_d,i_v,j_v,i_m)
          ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          DO i_d = 1,NDMAX
          DO j_d = 1,NDMAX
             WRITE(31,100)"AdveTenCoef",i_k,0,i_v,j_v,i_d,j_d,AdveTenCoef(i_d,j_d,i_k,i_v,j_v,i_m)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    !ENDDO
    
    !DO i_k = 1, NDMAX 
    !DO j_k = 1, NDMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          DO i_d = 1, NDMAX
          DO j_d = 1, NDMAX
             WRITE(31,100)"DiffTenCoef",0,0,i_v,j_v,i_d,j_d,DiffTenCoef(i_d,j_d,i_v,j_v,i_m)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    !DO i_k = 1, NKMAX
    !DO j_k = 1, NKMAX
       DO i_v = 1, NVMAX
       DO j_v = 1, NVMAX
          DO i_d = 1, NDMAX
          !DO j_d = 1, NDMAX
             WRITE(31,100)"GradVecCoef",0,0,i_v,j_v,i_d,0,GradVecCoef(i_d,i_v,j_v,i_m)
          ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO
    
    DO i_k = 1,NKMAX
    !DO j_k = 1,NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          DO i_d = 1,NDMAX
          DO j_d = 1,NDMAX
             WRITE(31,100)"GradTenCoef",i_k,0,i_v,j_v,i_d,j_d,GradTenCoef(i_d,j_d,i_k,i_v,j_v,i_m)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    !ENDDO
    ENDDO

    !DO i_k = 1,NKMAX
    !DO j_k = 1,NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          !DO i_d = 1,NDMAX
          !DO j_d = 1,NDMAX
             WRITE(31,100)"ExciScaCoef",0,0,i_v,j_v,0,0,ExciScaCoef(i_v,j_v,i_m)
          !ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    DO i_k = 1,NKMAX
    !DO j_k = 1,NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          DO i_d = 1,NDMAX
          !DO j_d = 1,NDMAX
             WRITE(31,100)"ExciVecCoef",i_k,0,i_v,j_v,i_d,0,ExciVecCoef(i_d,i_k,i_v,j_v,i_m)
          ENDDO
          ENDDO
       !ENDDO
       ENDDO
    !ENDDO
    ENDDO

    DO i_k = 1,NKMAX
    DO j_k = 1,NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          DO i_d = 1,NDMAX
          DO j_d = 1,NDMAX
             WRITE(31,100)"ExciTenCoef",i_k,j_k,i_v,j_v,i_d,j_d,ExciTenCoef(i_d,j_d,i_k,j_k,i_v,j_v,i_m)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

    !DO i_k = 1,NKMAX
    !DO j_k = 1,NKMAX
       DO i_v = 1,NVMAX
       DO j_v = 1,NVMAX
          !DO i_d = 1,NDMAX
          !DO j_d = 1,NDMAX
             WRITE(31,100)"SourScaCoef",0,0,i_v,j_v,0,0,SourScaCoef(i_v,j_v,i_m)
          !ENDDO
          !ENDDO
       ENDDO
       ENDDO
    !ENDDO
    !ENDDO

    CLOSE(31)
    
    RETURN
    
  END SUBROUTINE T2COEF_CHECK

  SUBROUTINE T2COEF_KV_TEST(     i_m)
    
    USE T2COMM,ONLY: GlobalCrd,KnownVar,TestCase 
    
    INTEGER(ikind),INTENT(IN)::i_m

    SELECT CASE (TestCase)
    CASE (1)
       KnownVar(1,i_m) = GlobalCrd(1,i_m)
    CASE (2)
       KnownVar(1,i_m) = GlobalCrd(1,i_m)
    CASE (3)
       KnownVar(1,i_m) = GlobalCrd(1,i_m)
    END SELECT

    RETURN
  END SUBROUTINE T2COEF_KV_TEST
  
  SUBROUTINE  T2COEF_MS_COEF_TEST(i_m)
    USE T2COMM,ONLY: NVMAX,TestMS,MassScaCoef,dt
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v

    ! INITIALIZATION
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       MassScaCoef(i_v,j_v,i_m) = 0.D0
    ENDDO
    ENDDO

    IF(TestMS) MassScaCoef(1,1,i_m) = 1.D0/dt

    RETURN

  END SUBROUTINE T2COEF_MS_COEF_TEST
 
  SUBROUTINE T2COEF_AV_COEF_TEST(i_m)
    
    USE T2COMM,ONLY: NVMAX,NDMAX,TestAV,AdveVecCoef,TestCase
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v,i_d
    ! initialization    
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_d = 1, NDMAX
          AdveVecCoef(i_d,i_v,j_v,i_m) = 0.D0
       ENDDO
    ENDDO
    ENDDO
    
    IF(TestAV)THEN
       SELECT CASE (TestCase)
       CASE (1)
          AdveVecCoef(1,1,1,i_m) = 1.D0
       CASE (2)
          AdveVecCoef(1,1,1,i_m) = 1.D0
       CASE (3)
          AdveVecCoef(1,1,1,i_m) = 1.D0
       END SELECT
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2COEF_AV_COEF_TEST

  SUBROUTINE  T2COEF_AT_COEF_TEST(i_m)
    USE T2COMM,ONLY: NVMAX,NKMAX,NDMAX,TestAT,TestCase,AdveTenCoef 
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v,i_d,j_d,i_k

    ! INITIALIZATION
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_k = 1, NKMAX
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             AdveTenCoef(i_d,j_d,i_k,i_v,j_v,i_m) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    IF(TestAT)THEN
       SELECT CASE (TestCase)
       CASE (1)
          AdveTenCoef(1,1,1,1,1,i_m) = -1.D0
       CASE (2)
          AdveTenCoef(2,1,1,1,1,i_m) = -1.D0
       END SELECT
    ENDIF

    RETURN

  END SUBROUTINE T2COEF_AT_COEF_TEST

  SUBROUTINE T2COEF_DT_COEF_TEST(i_m)
    USE T2COMM,ONLY: NVMAX,NDMAX,TestDT,TestCase,DiffTenCoef
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v,i_d,j_d

    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO j_d = 1, NDMAX
       DO i_d = 1, NDMAX
          DiffTenCoef(i_d,j_d,i_v,j_v,i_m) = 0.D0 
       ENDDO
       ENDDO
    ENDDO
    ENDDO

    IF(TestDT) DiffTenCoef(1,1,1,1,i_m) = 1.D0
    
    RETURN

  END SUBROUTINE T2COEF_DT_COEF_TEST
  
  SUBROUTINE  T2COEF_GV_COEF_TEST(i_m)

    USE T2COMM,ONLY: NVMAX,NDMAX,TestGV,TestCase,GradVecCoef
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v,i_d

    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_d = 1, NDMAX
          GradVecCoef(i_d,i_v,j_v,i_m) = 0.D0
       ENDDO
    ENDDO
    ENDDO
    
    IF(TestGV)THEN
       SELECT CASE (TestCase)
       CASE (1)
          GradVecCoef(1,1,1,i_m) = 1.D0
       CASE (2)
          GradVecCoef(2,1,1,i_m) = 1.D0
       END SELECT
    ENDIF

    RETURN
    
  END SUBROUTINE T2COEF_GV_COEF_TEST
  
  SUBROUTINE  T2COEF_GT_COEF_TEST(i_m)
    
    USE T2COMM,ONLY: NVMAX,NDMAX,NKMAX,TestGT,TestCase,GradTenCoef
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_d,i_k,i_v,j_d,    j_v
    
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_k = 1, NKMAX
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             GradTenCoef(i_d,j_d,i_k,i_v,j_v,i_m) = 0.D0
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    IF(TestGT)THEN
       SELECT CASE (TestCase)
       CASE (1)
          GradTenCoef(1,1,1,1,1,i_m) = 1.D0
       CASE (2)
          GradTenCoef(2,2,1,1,1,i_m) = 1.D0
       END SELECT
    END IF

    RETURN

  END SUBROUTINE T2COEF_GT_COEF_TEST
  
  SUBROUTINE T2COEF_ES_COEF_TEST(i_m)
    
    USE T2COMM, ONLY: NVMAX,TestES,ExciScaCoef
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v
    
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       ExciScaCoef(i_v,j_v,i_m) = 0.D0
    ENDDO
    ENDDO
    
    IF(TestES) ExciScaCoef(1,1,i_m) = LOG(2.D0)
    
    RETURN
  
  END SUBROUTINE T2COEF_ES_COEF_TEST
  
  SUBROUTINE  T2COEF_EV_COEF_TEST(i_m)

    USE T2COMM, ONLY: NVMAX,NKMAX,NDMAX,TestEV,ExciVecCoef

    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,i_k,i_d,j_v
  
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO i_k = 1, NKMAX 
          DO i_d = 1, NDMAX
             ExciVecCoef(i_d,i_k,i_v,j_v,i_m) = 0.D0
          ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    IF(TestEV) ExciVecCoef(1,1,1,1,i_m) = LOG(2.D0)
    
    RETURN
  
  END SUBROUTINE T2COEF_EV_COEF_TEST
  
  SUBROUTINE T2COEF_ET_COEF_TEST(i_m)

    USE T2COMM,ONLY:NVMAX,NKMAX,NDMAX,TestET,ExciTenCoef
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_d,i_k,i_v,j_d,j_k,j_v

    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       DO j_k = 1, NKMAX
       DO i_k = 1, NKMAX
          DO j_d = 1, NDMAX
          DO i_d = 1, NDMAX
             ExciTenCoef(i_d,j_d,i_k,j_k,i_v,j_v,i_m) = 0.D0
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    IF(TestET) ExciTenCoef(1,1,1,1,1,1,i_m) = LOG(2.D0)
    
    RETURN
    
  END SUBROUTINE T2COEF_ET_COEF_TEST
  
  SUBROUTINE T2COEF_SS_COEF_TEST(i_m)
    
    USE T2COMM, ONLY: NVMAX,TestSS,SourScaCoef
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::i_v,j_v
    
    ! initialization
    DO j_v = 1, NVMAX
    DO i_v = 1, NVMAX
       SourScaCoef(i_v,j_v,i_m) = 0.D0
    ENDDO
    ENDDO
    
    IF(TestSS) SourScaCoef(1,1,i_m) = 1.D0

    RETURN
    
  END SUBROUTINE T2COEF_SS_COEF_TEST
END MODULE T2COEF
