!     $Id: fpcoef.f90,v 1.38 2013/02/08 07:36:24 nuga Exp $
!
! ************************************************************
!
!                      CALCULATION OF D AND F
!
! ************************************************************
!
      MODULE fpfunc

      USE fpcomm
      USE libbes,ONLY: besekn

      INTEGER,parameter:: MODEL_T_IMP=2
      contains
!-------------------------------------------------------------
! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_prf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.eq.0)THEN
         RL=0.D0
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF

      IF(NR.eq.0.or.NR.eq.NRMAX+1)THEN
         IF(MODEL_EX_READ_Tn.eq.0)THEN 
            CALL PL_PROF(RHON,PLF)
            RNFDL=PLF(NS)%RN/RNFD0L
            RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
         ELSE
            IF(NR.eq.0)THEN
               RNFDL=RN_TEMP(1,NS)/RNFD0L
               RTFDL=RT_TEMP(1,NS)
            ELSEIF(NR.eq.NRMAX+1)THEN
               RNFDL=RNE_EXP_EDGE/RNFD0L
               RTFDL=RTE_EXP_EDGE
            END IF
         END IF
      ELSE
         RNFDL=RN_TEMP(NR,NS)/RNFD0L
         RTFDL=RT_TEMP(NR,NS)
      END IF
      
      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         FPMXWL=FACT*EXP(-EX)
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
         DKBSL=BESEKN(2,Z)
         FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
              *RTFD0L
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         FPMXWL=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL
!-------------------------------------------------------------
      FUNCTION FPMXWL_S(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_prf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_S

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(MODEL_EX_READ_Tn.eq.0)THEN
         CALL PL_PROF(RHON,PLF)
         RNFDL=PNS(NS)/RNFD0L*1.D-1
         RTFDL=PTS(NS)*1.D-2
      ELSEIF(MODEL_EX_READ_Tn.ne.0)THEN
         RNFDL=RNE_EXP_EDGE/RNFD0L*1.D-1
         RTFDL=RTE_EXP_EDGE*1.D-1
      END IF


      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL_S=0.D0
         ELSE
            FPMXWL_S=FACT*EXP(-EX)
         ENDIF
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
         DKBSL=BESEKN(2,Z)
         FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
              *RTFD0L
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         FPMXWL_S=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_S
!-------------------------------------------------------------

      SUBROUTINE FPMXWL_EDGE(NP,NSA,FL)

      implicit none
!      integer,intent(in):: NP, NSA
      integer:: NP, NSA
      integer:: NS
      real(kind8):: FL1, FL2
!      real(kind8),intent(out):: FL
      real(kind8):: FL

      NS=NS_NSA(NSA)

!      FL1=FPMXWL_S(PM(NP,NS),NRMAX,NS) ! at RM(NRMAX)
!      FL2=FPMXWL_S(PM(NP,NS),NRMAX+1,NS) ! at RG(NRMAX+1)

!     F at R=1.0+DELR/2
!      FL=FL2*1.D-1

      FL=FPMXWL_S(PM(NP,NS),NRMAX,NS) 

      RETURN
      END SUBROUTINE FPMXWL_EDGE
!-------------------------------------------------------------
      FUNCTION FPMXWL_LT(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_prf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_LT

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.eq.0)THEN
         RL=0.D0
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF

      CALL PL_PROF(RHON,PLF)
      RNFDL=PLF(NS)%RN/RNFD0L
!      RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
      RTFDL=PTS(NS)*1.D-1

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL_LT=0.D0
         ELSE
            FPMXWL_LT=FACT*EXP(-EX)
         ENDIF
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
            DKBSL=BESEKN(2,Z)
            FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
             *RTFD0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
            FPMXWL_LT=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_LT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-------------------------------------------------------------
      FUNCTION FPMXWL_IMP(PML,NR,NSA,N_ion)

      implicit none
      integer,intent(in):: NR,NSA,N_ion
      integer:: NS
      real(kind8),intent(in):: PML
      real(kind8):: rnfp0l,rtfp0l,ptfp0l,RTFD0L
      real(kind8):: amfpl,rnfpl,rtfpl,fact,ex,theta0l,thetal,z,dkbsl
      real(kind8):: FPMXWL_IMP, target_Z, target_ne, target_ni

      NS=NS_NSA(NSA)
      AMFPL=PA(NS)*AMP
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      RNFP0L=PN(NS)
      RTFP0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFP0L=SQRT(RTFD0L*1.D3*AEE*AMFPL)

      target_z=ABS(PZ(NS))
      target_ni=(SPITOT-RNFP0(1))*RNFP(NR,NSA)/RNFP0(NSA)/target_z

      RNFPL=target_ni
      
      IF(MODEL_T_IMP.eq.0)THEN
         RTFPL=RT_quench(NR)
      ELSEIF(MODEL_T_IMP.eq.1)THEN
         RTFPL=( 0.1D0*RT_quench(NR) + 0.9D0*RT_quench_f(NR) )
      ELSEIF(MODEL_T_IMP.eq.2)THEN
         RTFPL=( 0.01D0*RT_quench(NR) + 0.99D0*RT_quench_f(NR) )
      END IF

      IF(MODELR.EQ.0) THEN
         FACT=RNFPL/SQRT(2.D0*PI*RTFPL/RTFP0L)**3
         EX=PML**2/(2.D0*RTFPL/RTFP0L)
         IF(EX.GT.100.D0) THEN
            FPMXWL_IMP=0.D0
         ELSE
            FPMXWL_IMP=FACT*EXP(-EX)
         ENDIF
      ELSE
         THETA0L=THETA0(NS)
         THETAL=THETA0L*RTFPL/RTFP0L
         Z=1.D0/THETAL
            DKBSL=BESEKN(2,Z)
            FACT=RNFPL*SQRT(THETA0L)/(4.D0*PI*RTFPL*DKBSL) &
             *RTFP0L
            EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
            FPMXWL_IMP=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_IMP
!-------------------------------------------------------------
! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL_EXP(PML,NR,NS)

      USE plprof
      implicit none
      integer :: NR, NS, NS_beam
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      real(kind8):: FPMXWL_EXP, sum_beam_dens

      AMFDL=PA(NS)*AMP
      AEFDL=PZ(NS)*AEE
      RNFD0L=PN(NS)
      RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      PTFD0L=SQRT(RTFD0L*1.D3*AEE*AMFDL)

      IF(NR.eq.0)THEN
         RL=0.D0
         RHON=ABS(RL)
      ELSEIF(NR.EQ.NRMAX+1) THEN
         RL=RM(NRMAX)+DELR
         RHON=MIN(RL,1.D0)
      ELSE
         RL=RM(NR)
         RHON=RL
      ENDIF

      IF(MODEL_DELTA_F(NS).eq.1)THEN
         RNFDL=RN_BULK(NR,NS)/RNFD0L
      ELSE
         RNFDL=(RN_TEMP(NR,NS))/RNFD0L
      END IF

      RTFDL=RT_TEMP(NR,NS)

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         FPMXWL_EXP=FACT*EXP(-EX)
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
         DKBSL=BESEKN(2,Z)
         FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
              *RTFD0L
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         FPMXWL_EXP=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL_EXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE fpfunc
