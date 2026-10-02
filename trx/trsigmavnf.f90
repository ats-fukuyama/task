! trsigmavnf.f90

MODULE trsigmavnf

  PRIVATE
  PUBLIC sigmav_DTmm

  CONTAINS

!     ***********************************************************

!           REACTION CROSS SECTION (MAXELLIAN) DT

!     ***********************************************************

    FUNCTION sigmav_DTmm(TD,TT)
      USE bpsd_kinds
      USE libsigma
      IMPLICIT NONE
      REAL(rkind),INTENT(IN)::  TD,TT
      REAL(rkind):: sigmav_DTmm
      REAL(rkind):: TI

         TI = (3.D0*ABS(TD)+2.D0*ABS(TT))/5.D0

         sigmav_DTmm=sigmavm_dt(TI)*(1E-6)
         
         RETURN
       END FUNCTION sigmav_DTmm

!     ***********************************************************

!           REACTION CROSS SECTION (MAXELLIAN) DD

!     ***********************************************************

         FUNCTION sigmav_DDmm(TD)

            ! not completed
    
          USE trcomm,ONLY: rkind
          IMPLICIT NONE
          REAL(rkind),INTENT(IN):: TD
          REAL(rkind):: sigmav_DDmm
          REAL(rkind) TI,H,ARG
    
          TI=TD
         !  H  = TI/37.D0 + 5.45D0/(3.D0+TI*(1.D0+(TI/37.5D0)**2.8D0))
         !  ARG= -20.D0/TI**(1.D0/3.D0)
          !selection of fusion reaction cross section approximation
         !  IF(MDLSS=0) THEN
            ! IF(ARG.GE.-100.D0)  THEN
         !        SIGMAM_DD = 3.7D-18*TI**(-2.D0/3.D0)*EXP(ARG)/H
         !    ELSE
         !        SIGMAM_DD = 0.D0
         !    ENDIF
          !  ELSE
          
            IF(TI.GE.4.D0) THEN 
               sigmav_DDmm = ((2.1069E-20) &
                     + (-3.7748E-20)*TI &
                     + (1.1242E-20)*TI**2 &
                     + (-1.9831E-22)*TI**3 &
                     + (1.3421E-24)*TI**4)*1E-6
            ELSE
               sigmav_DDmm = (4.428*1E-20)*1E-6
            ENDIF    
         !  ENDIF  
    
          RETURN
        END FUNCTION sigmav_DDmm

!     ***********************************************************

!      REACTION RATE : TAIL

!     ***********************************************************

      FUNCTION sigmav_DTbm(EB,EC,TI,PTNT)

!      APPROXIMATE FORMULA OF FUSION REACTION RATE
!         FOR SLOWING DOWN ION DISTRIBUTION
!      REF. TAKIZUKA AND YAMAGIWA, JAERI-M 87-066

!      EB   : BEAM ENERGY (KEV)
!      EC   : CRITICAL ENERGY (KEV)
!      TI   : TRITIUM TEMPERATURE (KEV)
!      PTNT : PB * TAUS / (ND * EB)

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: EB,EC,TI,PTNT
      REAL(rkind):: sigmav_DTbm
      REAL(rkind) XB,XC,AG1,AG2,AG3,AL1,AL2,AL3,X1,X2,X3,X4,SA

      XB=SQRT(EB/127.D0)
      XC=SQRT(EC/127.D0)

      AG1= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG2= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG3= 0.33D0
      AL1= 1.D0/(0.40D0+0.032D0*TI)
      AL2=-1.D0/(0.91D0+0.016D0*SQRT(ABS(TI)))
      AL3=-0.11D0
      X1=0.97D0-AG1/AL1
      X2=0.97D0
      X3=0.97D0+(AG2-AG3)/(AL3-AL2)
      X4=0.97D0+3.D0

      IF(XB.LT.X1) THEN
         SA=0.D0
      ELSE
         SA=-sigmav_DTbm_sub(X1,AG1,AL1,XC)
         IF(XB.LT.X2) THEN
            SA=SA+sigmav_DTbm_sub(XB,AG1,AL1,XC)
         ELSE
            SA=SA+sigmav_DTbm_sub(X2,AG1,AL1,XC) &
                 -sigmav_DTbm_sub(X2,AG2,AL2,XC)
            IF(XB.LT.X3) THEN
               SA=SA+sigmav_DTbm_sub(XB,AG2,AL2,XC)
            ELSE
               SA=SA+sigmav_DTbm_sub(X3,AG2,AL2,XC) &
                    -sigmav_DTbm_sub(X3,AG3,AL3,XC)
               IF(XB.LT.X4) THEN
                  SA=SA+sigmav_DTbm_sub(XB,AG3,AL3,XC)
               ELSE
                  SA=SA+sigmav_DTbm_sub(X4,AG3,AL3,XC)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      sigmav_DTbm=PTNT*1.67D-21*SA
      RETURN
    END FUNCTION sigmav_DTbm

    FUNCTION sigmav_DTbm_sub(XX,RGG,RGL,XC)

      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: XX,RGG,RGL,XC
      REAL(rkind):: sigmav_DTbm_sub
      REAL(rkind) :: X

      X=XX/XC
      sigmav_DTbm_sub=((RGG-0.97D0*RGL)/3.D0+XC*RGL/6.D0)*LOG(X*X*X+1.D0) &
           +XC*RGL*(X-LOG(X+1.D0)/2.D0 &
           -ATAN((2.D0*X-1.D0)/SQRT(3.D0))/SQRT(3.D0))
      RETURN
    END FUNCTION sigmav_DTbm_sub
  END MODULE trsigmavnf
