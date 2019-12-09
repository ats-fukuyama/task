MODULE fpsub

  PRIVATE
!  PUBLIC fpbave_dpp
!  PUBLIC fpbave_dth
  PUBLIC FPMXWL
  PUBLIC FPMXWL_S
  PUBLIC FPMXWL_EDGE
  PUBLIC FPMXWL_LT
  PUBLIC update_fnsb_maxwell
  PUBLIC FPNEWTON
  PUBLIC PROF_OF_NF_REACTION_RATE
  
CONTAINS

  SUBROUTINE fpbave_dpp(DPP,NR,NSA,ID)
    USE fpcomm,ONLY: NTHMAX,NPSTART,NPENDWG,NRSTART,NRENDWM,NSAMAX, &
                     RLAMDA,ITL,ITU
    IMPLICIT NONE
    REAL(8),INTENT(INOUT):: &
         DPP(NTHMAX  ,NPSTART :NPENDWG,NRSTART:NRENDWM,NSAMAX)
    INTEGER,INTENT(IN):: NR,NSA,ID  ! ID=0 for DPP,FPP, ID=1 for DPT
    INTEGER:: NP,NTH

    DO NP=NPSTART,NPENDWG
       IF(ID.EQ.0) THEN
          DO NTH=ITL(NR)+1,NTHMAX/2
             DPP(NTH,NP,NR,NSA)    =0.5D0*(DPP(NTH,NP,NR,NSA) &
                                          +DPP(NTHMAX-NTH+1,NP,NR,NSA))
             DPP(NTHMAX-NTH+1,NP,NR,NSA)  =DPP(NTH,NP,NR,NSA)
          END DO
       ELSE
          DO NTH=ITL(NR)+1,NTHMAX/2
             DPP(NTH,NP,NR,NSA)    =0.5D0*(DPP(NTH,NP,NR,NSA) &
                                          -DPP(NTHMAX-NTH+1,NP,NR,NSA))
             DPP(NTHMAX-NTH+1,NP,NR,NSA) =-DPP(NTH,NP,NR,NSA)
          END DO
       END IF
       IF(ID.EQ.0) THEN
          DPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
                      *( DPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                        +DPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                        +DPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                        +DPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR)) 
          DPP(ITU(NR),NP,NR,NSA)=DPP(ITL(NR),NP,NR,NSA)
       ELSE
          DPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0 &
                      *( DPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR) &
                        +DPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR) &
                        -DPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR) &
                        -DPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR)) 
          DPP(ITU(NR),NP,NR,NSA)=-DPP(ITL(NR),NP,NR,NSA)
       END IF
    END DO
    RETURN
  END SUBROUTINE fpbave_dpp

  SUBROUTINE fpbave_dth(DTH,NR,NSA,ID)
    USE fpcomm,ONLY: NTHMAX,NPSTARTW,NPENDWM,NRSTART,NRENDWM,NSAMAX, &
                     RLAMDA,ITL,ITU
    IMPLICIT NONE
    REAL(8),INTENT(INOUT):: &
         DTH(NTHMAX+1,NPSTARTW:NPENDWM,NRSTART:NRENDWM,NSAMAX)
    INTEGER,INTENT(IN):: NR,NSA,ID  ! ID=0 for DTT, ID=1 for DTP,FTH
    INTEGER:: NP,NTH

    DO NP=NPSTARTW,NPENDWM
       DO NTH=ITL(NR)+1,NTHMAX/2
          IF(ID.EQ.0) THEN
             DTH(NTH,NP,NR,NSA)    =0.5D0*(DTH(NTH,         NP,NR,NSA) &
                                          +DTH(NTHMAX-NTH+2,NP,NR,NSA))
             DTH(NTHMAX-NTH+2,NP,NR,NSA)  =DTH(NTH,NP,NR,NSA)
          ELSE
             DTH(NTH,NP,NR,NSA)    =0.5D0*(DTH(NTH,         NP,NR,NSA) &
                                          -DTH(NTHMAX-NTH+2,NP,NR,NSA))
             DTH(NTHMAX-NTH+2,NP,NR,NSA) =-DTH(NTH,NP,NR,NSA)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE fpbave_dth

!-------------------------------------------------------------

! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL(PML,NR,NS)

      USE fpcomm
      USE plprof
      USE libbes
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
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

      CALL PL_PROF(RHON,PLF)
      RNFDL=PLF(NS)%RN/RNFD0L
      RTFDL=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0

      IF(MODELR.EQ.0) THEN
         FACT=RNFDL/SQRT(2.D0*PI*RTFDL/RTFD0L)**3
         EX=PML**2/(2.D0*RTFDL/RTFD0L)
         FPMXWL=FACT*EXP(-EX)
      ELSE
         THETA0L=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
         THETAL=THETA0L*RTFDL/RTFD0L
         Z=1.D0/THETAL
         DKBSL=BESEKNX(2,Z)
         FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
              *RTFD0L
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         FPMXWL=FACT*EXP(EX)
      END IF

      RETURN
      END FUNCTION FPMXWL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Subroutine FPNEWTON(NR,NSA,RNSL_,RWSL_,rtemp)

      USE fpcomm,ONLY: NS_NSA,AMFP,THETA0,PTFP0,VC,AEE
      USE libbes
      IMPLICIT NONE
      INTEGER,INTENT(IN)::NR,NSA
      double precision,intent(in):: RNSL_, RWSL_ 
      double precision,intent(out)::rtemp
      integer:: ncount, NS
      real(8):: xeave
      real(8):: xtemp, thetal, EAVE

      NS=NS_NSA(NSA)
!-----Average kinetic energy
!      EAVE=RWS(NR,NSA)*AMFP(NSA)*THETA0(NS) &
!           /(RNS(NR,NSA)*1.D20*PTFP0(NSA)**2*1.D-6)
      EAVE=RWSL_*AMFP(NSA)*THETA0(NS) & 
           /(RNSL_*1.D20*PTFP0(NSA)**2*1.D-6)
!-----initial value of THETAL
      THETAL=2.d0*EAVE/3.d0
      xtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)

      CALL XNEWTON(EAVE,THETAL,ncount)

      rtemp=AMFP(NSA)*VC**2*THETAL/(AEE*1.D3)
      xeave=AMFP(NSA)*VC**2*EAVE/(AEE*1.D3)

      RETURN

      CONTAINS

      SUBROUTINE xnewton(eave,thetal,ncount)
      IMPLICIT NONE
      REAL(8),intent(in):: eave
      REAL(8),intent(inout):: thetal
      INTEGer,intent(out):: ncount
      REAL(8),parameter:: eps=1.d-10
      REAL(8):: delthetal,thetalnew,epsthetal

!--------iteration start
      ncount=0
      DO while(ncount.le.100)
         ncount=ncount+1
         delthetal=-(rfunc(thetal)-eave)/dfunc(thetal)
         thetalnew=thetal+delthetal
         epsthetal=ABS(delthetal/thetal)

         thetal=thetalnew
         IF(epsthetal.le.eps) EXIT
      END DO
      RETURN
      END SUBROUTINE xnewton

      FUNCTION rfunc(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,rfunc
      REAL(8):: z,dkbsl1,dkbsl2
      z=1.D0/thetal
      dkbsl1=BESEKNX(1,Z)
      dkbsl2=BESEKNX(2,Z)
      rfunc= (dkbsl1 /dkbsl2 -1.D0+3.D0/Z)
      RETURN
      END FUNCTION rfunc

      FUNCTION rfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,rfuncp
      REAL(8):: z,dkbsl1,dkbsl2
      z=1.D0/thetal
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      rfuncp= dkbsl1 /dkbsl2 -1.D0+3.D0/Z
      RETURN
      END FUNCTION rfuncp
      
      FUNCTION dfunc(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,dfunc
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3
      z=1.D0/thetal
      dkbsl0=BESEKNX(0,z)
      dkbsl1=BESEKNX(1,z)
      dkbsl2=BESEKNX(2,z)
      dkbsl3=BESEKNX(3,z)
      dfunc =( (dkbsl0 +dkbsl2 )/dkbsl2                           &
                -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2 &
            +3.d0
      RETURN
      END FUNCTION dfunc

      FUNCTION dfuncp(thetal)
      IMPLICIT NONE
      REAL(8):: thetal,dfuncp
      REAL(8):: z,dkbsl0,dkbsl1,dkbsl2,dkbsl3
      z=1.D0/thetal
      dkbsl0=1.D0 -  1.D0/8.D0/z +   9.D0/128.D0/z**2
      dkbsl1=1.D0 +  3.D0/8.D0/z -  15.D0/128.D0/z**2
      dkbsl2=1.D0 + 15.D0/8.D0/z + 105.D0/128.D0/z**2
      dkbsl3=1.D0 + 35.D0/8.D0/z + 945.D0/128.D0/z**2
      dfuncp =( (dkbsl0 +dkbsl2 )/dkbsl2                          &
               -(dkbsl1 +dkbsl3 )*dkbsl1 /dkbsl2 **2)*0.5d0*z**2  & 
            +3.d0
      RETURN
      END FUNCTION dfuncp
      
      end Subroutine FPNEWTON

      SUBROUTINE PROF_OF_NF_REACTION_RATE(ID)

      USE fpcomm
      use libmpi
      use fpmpi
      IMPLICIT NONE
      integer,intent(in):: ID
      integer:: NR, ndata, NS
      double precision,dimension(NRSTART:NREND):: send1, send2
      double precision,dimension(NRMAX):: recv1, recv2, sigmav_mx
      double precision:: tot1, tot2, sum

      CALL mtx_set_communicator(comm_nr)
      ndata=NREND-NRSTART+1
      DO NR=NRSTART, NREND
         send1(NR)=RATE_NF(NR,ID)
         send2(NR)=RATE_NF_BB(NR,ID)
      END DO
      CALL mtx_allgather_real8(send1,ndata,recv1)
      CALL mtx_allgather_real8(send2,ndata,recv2)

      CALL mtx_reset_communicator

      tot1=0.D0
      tot2=0.D0
      sum=0.D0
      NS=NS_NSA(NSA_F1)
      DO NR=1, NRMAX
         tot1 = tot1 + recv1(NR)*VOLR(NR)
         tot2 = tot2 + recv2(NR)*VOLR(NR)
         sigmav_mx(NR) = 2.33D-14*RT_TEMP(NR,NS)**(-2.D0/3.D0)*EXP(-18.76D0*RT_TEMP(NR,NS)**(-1.D0/3.D0))*1.D-6 &
              *RN_TEMP(NR,NS)**2*1.D20*0.5D0*0.5D0 ! double count, only ID=1,2,5
         sum = sum + sigmav_mx(NR)
      END DO

      IF(NRANK.eq.0)THEN
         WRITE(25,'(99E14.6)') TIMEFP, (recv1(NR), NR=1, NRMAX), (recv2(NR), NR=1, NRMAX), (sigmav_mx(NR), NR=1, NRMAX)
         WRITE(26,'(10E14.6)') TIMEFP, tot1, tot2, sum
      END IF

      END SUBROUTINE PROF_OF_NF_REACTION_RATE

    END MODULE fpsub
