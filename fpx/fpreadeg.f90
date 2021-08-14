! fpreadeg

MODULE fpreadeg

      USE fpcomm

      PUBLIC READ_EXP_DATA
      PUBLIC MAKE_EXP_PROF
      PUBLIC FPMXWL_EXP
      
      PRIVATE
!      REAL(rkind),dimension(:,:),ALLOCATABLE:: read_tms_double, read_cx_double
!      integer,dimension(:,:),ALLOCATABLE:: read_tms_int
!      REAL(rkind),dimension(5):: cte_fit
!      REAL(rkind),dimension(6):: cne_fit
!      integer:: nend_tms, nend_cx
!      integer,parameter:: NRMAX=32
!      REAL(rkind),dimension(32):: RM
!      REAL(rkind),dimension(32):: RNE_EXP, RTE_EXP

      contains
!-----------------------------------
      SUBROUTINE READ_EXP_DATA ! EGDATA OF TMS and CX are READ

      IMPLICIT NONE

      CALL READ_EXP_TMS
      CALL READ_EXP_CX 

      END SUBROUTINE READ_EXP_DATA
!-----------------------------------
      SUBROUTINE MAKE_EXP_PROF(time) ! RNE_EXP and RTE_EXP are updated

      IMPLICIT NONE
      REAL(rkind),intent(in):: time
      REAL(rkind):: weight
      integer:: ntime1, NS, NR, NSA

      IF(MODEL_EX_READ_Tn.ne.0)THEN
         CALL time_interpolation_tms(time, ntime1, weight)
         CALL MAKE_PROF_FROM_TMS(ntime1,weight)
      END IF

      IF(MODEL_EX_READ_Tn.eq.2)THEN
         CALL time_interpolation_cx(time, ntime1, weight)
         CALL MAKE_PROF_FROM_CX(ntime1,weight)
      END IF

      DO NS=1, NSMAX
         IF(NS.eq.1)THEN
            DO NR=1, NRMAX
               RT_READ(NR,NS)=RTE_EXP(NR)
               RN_READ(NR,NS)=RNE_EXP(NR)
               
               RT_TEMP(NR,NS)=RTE_EXP(NR)
               RN_TEMP(NR,NS)=RNE_EXP(NR)
            END DO
         ELSE
            IF(MODEL_EX_READ_Tn.eq.1)THEN
               DO NR=1, NRMAX
                  RT_READ(NR,NS)=RTE_EXP(NR) ! assume to be same to electron
                  RN_READ(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
               
                  RT_TEMP(NR,NS)=RTE_EXP(NR)
                  RN_TEMP(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
               END DO
            ELSEIF(MODEL_EX_READ_Tn.eq.2)THEN
               DO NR=1, NRMAX
                  RT_READ(NR,NS)=RTI_EXP(NR)
                  RN_READ(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
               
                  RT_TEMP(NR,NS)=RTI_EXP(NR)
                  RN_TEMP(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
               END DO
            END IF
         END IF
      END DO
      DO NSA=1, NSAMAX
         NS=NS_NSA(NSA)
         IF(NS.eq.1)THEN
            DO NR=1,NRMAX
               RT_BULK(NR,NSA)=RTE_EXP(NR)
            END DO
         ELSE
            IF(MODEL_EX_READ_Tn.eq.1)THEN
               DO NR=1, NRMAX
                  RT_BULK(NR,NSA)=RTE_EXP(NR)
               END DO
            ELSEIF(MODEL_EX_READ_Tn.eq.2)THEN
               DO NR=1, NRMAX
                  RT_BULK(NR,NSA)=RTI_EXP(NR)
               END DO
            END IF
         END IF
      END DO
      END SUBROUTINE MAKE_EXP_PROF
!-----------------------------------
      SUBROUTINE READ_EXP_TMS

      IMPLICIT NONE
      integer:: i,j,k,m,nend1
!      REAL(rkind),dimension(27,1000):: read_temporary_double
      REAL(rkind),allocatable:: read_temporary_double(:,:)
      integer,dimension(2,1000):: read_temporary_int

      ALLOCATE(read_temporary_double(27,1000))

!      open(22,file='./dat/tswpe@119489.dat',status='old')
      open(22,file=EG_NAME_TMS,status='old')

      DO i=1,29
         READ(22,*) 
      END DO

      nend1=0
      DO i=1,1000
         READ(22,*,end=100) &
              (read_temporary_double(j,i),j=1,3), (read_temporary_int(k,i),k=1,2), (read_temporary_double(m,i),m=4,27)
         nend1=nend1+1
      END DO

100   close(22)
      nend_tms=nend1

      allocate(read_tms_double(27,nend1))
      allocate(read_tms_int(2,nend1))

      DO i=1,nend1
         DO j=1,27
            read_tms_double(j,i)=read_temporary_double(j,i)
         END DO
         DO j=1,2
            read_tms_int(j,i)=read_temporary_int(j,i)
         END DO
      END DO

!       WRITE(*,*) "END LINE", nend1, nend_tms

!       DO i=1, 2
!          WRITE(*,'(1P3E14.6,I8,I3,1P24E12.4)') &
!               (read_tms_double(j,i),j=1,3), (read_tms_int(k,i),k=1,2), (read_tms_double(m,i),m=4,27) 
!       END DO

      END SUBROUTINE READ_EXP_TMS
!--------------------------------
      SUBROUTINE READ_EXP_CX

      IMPLICIT NONE
      integer:: i,j,nend2, k
!      REAL(rkind),dimension(34,1000):: read_cx_temp, read_cx_temp2
      REAL(rkind),allocatable:: read_cx_temp(:,:), read_cx_temp2(:,:)

      ALLOCATE(read_cx_temp(34,1000), read_cx_temp2(34,1000))

      !      open(22,file='./dat/cxswpi7@119489.dat',status='old')
      open(22,file=EG_NAME_CX,status='old')

      DO i=1,21
         READ(22,*) 
      END DO
      
      nend2=0
      DO i=1,1000
         READ(22,*,end=200) (read_cx_temp(j,i),j=1,34)
!         READ(22,*,end=200) (read_cx_double(j,i),j=1,34)
         nend2=nend2+1
      END DO
      
200   close(22)
      nend_cx=nend2


      k=0
      DO i=1, nend_cx
         IF(read_cx_temp(29,i).ne.0.D0)THEN
            k=k+1
            DO j=1,34
               read_cx_temp2(j,k)=read_cx_temp(j,i)
            END DO
         END IF
      END DO

      nend_cx=k
      allocate(read_cx_double(34,nend_cx))

      DO i=1, nend_cx
         DO j=1,34
            read_cx_double(j,i)=read_cx_temp2(j,i)
         END DO
      END DO
      
!       WRITE(*,*) "END LINE", nend2

!       DO i=nend2-1,nend2
!          WRITE(*,'(1P34E14.6)') (read_cx_double(j,i),j=1,34)
!       END DO

      END SUBROUTINE READ_EXP_CX
!--------------------------------
      SUBROUTINE MAKE_PROF_FROM_TMS(ntime1,weight)

      IMPLICIT NONE
      integer,intent(in):: ntime1
      REAL(rkind),intent(in):: weight
      integer:: i,j,k
      integer:: NR
      REAL(rkind):: rte_ex, rne_ex, rho
      REAL(rkind),dimension(NRMAX,2):: prof_ne_temp, prof_te_temp
      
      DO k=1, 2
         i=k-1+ntime1
         
         DO j=1,5
            cte_fit(j)=read_tms_double(14+j,i)
         END DO
         DO j=1,5
            cne_fit(j)=read_tms_double(18+j,i)
         END DO

!          WRITE(*,*) " cte_fit(1)", i, cte_fit(1)
!          WRITE(*,*) " cne_fit(1)", i, cne_fit(1)
         
         DO NR=1, NRMAX
            rho = RM(NR)
            rte_ex = cte_fit(1)+cte_fit(2)*rho**2+cte_fit(3)*rho**4+cte_fit(4)*rho**6
            rne_ex = cne_fit(1)+cne_fit(2)*rho**2+cne_fit(3)*rho**4+cne_fit(4)*rho**6+cne_fit(5)*rho**8
            IF(rte_ex.gt.0.D0)THEN
               prof_te_temp(NR,k)=rte_ex
            ELSE
               prof_te_temp(NR,k)=1.D-1
            END IF
            IF(rne_ex.gt.0.D0)THEN
               prof_ne_temp(NR,k)=rne_ex
            ELSE
               prof_ne_temp(NR,k)=1.D-3               
            END IF
         END DO
      END DO
      
! time interpolation
      DO NR=1, NRMAX
         RTE_EXP(NR)=(1.D0-weight)*prof_te_temp(NR,1) + weight*prof_te_temp(NR,2)
         RNE_EXP(NR)=( (1.D0-weight)*prof_ne_temp(NR,1) + weight*prof_ne_temp(NR,2) )*1.D-1 ! 10^20
      END DO
      RTE_EXP_EDGE= read_tms_double(24,i) 
      RNE_EXP_EDGE= read_tms_double(25,i)*1.D-1
       
      open(22,file='prof_tms.dat')
      IF(NRANK.eq.0)THEN
         WRITE(22,'(A,3E14.6)') "# TIME, TE_EDGE, NE_EDGE", TIMEFP, RTE_EXP_EDGE, RNE_EXP_EDGE
         DO NR=1, NRMAX
            WRITE(22,'(I5,7E14.6)') NR, RM(NR), prof_te_temp(NR,1), RTE_EXP(NR), prof_te_temp(NR,2), &
                 prof_ne_temp(NR,1), RNE_EXP(NR), prof_ne_temp(NR,2)
         END DO
         WRITE(22,*) " "
         WRITE(22,*) " "
      END IF
      close(22)
      
      open(22,file='evol_tms.dat')
      IF(NRANK.eq.0)THEN
         DO i=1, nend_tms 
            IF(read_tms_double(15,i).ne.0.and.read_tms_double(19,i).ne.0)THEN ! time, Te, n_e
               WRITE(22,'(3E14.6)') read_tms_double(1,i), read_tms_double(15,i), read_tms_double(19,i)
            END IF
         END DO
      END IF
      close(22)

      END SUBROUTINE MAKE_PROF_FROM_TMS
!--------------------------------
      SUBROUTINE MAKE_PROF_FROM_CX(ntime1,weight)

      IMPLICIT NONE
      integer,intent(in):: ntime1
      REAL(rkind),intent(in):: weight
      integer:: i,j,k
      integer:: NR
      REAL(rkind):: rti_ex, rho
      REAL(rkind),dimension(NRMAX,2):: prof_ti_temp
      
      DO k=1, 2
         i=k-1+ntime1
         
         DO j=1,5
            cti_fit(j)=read_cx_double(28+j,i)
         END DO

!          WRITE(*,*) " cte_fit(1)", i, cte_fit(1)
!          WRITE(*,*) " cne_fit(1)", i, cne_fit(1)
         
         DO NR=1, NRMAX
            rho = RM(NR)
            rti_ex = cti_fit(1)+cti_fit(2)*rho**2+cti_fit(3)*rho**4+cti_fit(4)*rho**6
            
            prof_ti_temp(NR,k)=rti_ex
         END DO
      END DO
      
! time interpolation
      DO NR=1, NRMAX
         RTI_EXP(NR)=(1.D0-weight)*prof_ti_temp(NR,1) + weight*prof_ti_temp(NR,2)
      END DO
      RTI_EXP_EDGE= read_cx_double(34,i) 
       
      open(22,file='prof_cx.dat')
      IF(NRANK.eq.0)THEN
         WRITE(22,'(A,3E14.6)') "# TIME, TI_EDGE ", TIMEFP, RTI_EXP_EDGE
         DO NR=1, NRMAX
            WRITE(22,'(I5,7E14.6)') NR, RM(NR), prof_ti_temp(NR,1), RTI_EXP(NR), prof_ti_temp(NR,2)
         END DO
         WRITE(22,*) " "
         WRITE(22,*) " "
      END IF
      close(22)
      
      open(22,file='evol_cx.dat')
      IF(NRANK.eq.0)THEN
         DO i=1, nend_cx 
            IF(read_cx_double(28,i).ne.0)THEN ! time, Ti
               WRITE(22,'(2E14.6)') read_cx_double(1,i), read_cx_double(29,i)
            END IF
         END DO
      END IF
      close(22)

      END SUBROUTINE MAKE_PROF_FROM_CX
!------------------------------------
      SUBROUTINE time_interpolation_cx(time, ntime1, weight)

      IMPLICIT NONE
      REAL(rkind), intent(in):: time
      integer,intent(out):: ntime1
      integer:: ntime2
      REAL(rkind),intent(out):: weight ! y=(1-alpha)*y1 + alpha*y2
      REAL(rkind):: time_exp
      integer:: i

      time_exp= time + time_exp_offset
!      WRITE(*,*) "TEST_NEND_CX ", nend_cx
      IF(time_exp.lt.read_cx_double(1,1))THEN
         ntime1=1
         ntime2=1
         weight=1.D0
      ELSEIF(read_cx_double(1,nend_cx).lt.time_exp)THEN
         ntime1=nend_cx
         ntime2=nend_cx
         weight=1.D0
      ELSE
         DO i=1, nend_cx
            IF(read_cx_double(1,i).le.time_exp.and.time_exp.lt.read_cx_double(1,i+1))THEN
               ntime1=i
               ntime2=i+1
            END IF
!      IF(NRANK.eq.0) WRITE(*,*) "TEST", i, nend_cx, read_cx_double(1,1), read_cx_double(1,2)
         END DO
         weight=(time_exp-read_cx_double(1,ntime1))/(read_cx_double(1,ntime2)-read_cx_double(1,ntime1))
      END IF
 

!      IF(NRANK.eq.0)THEN
!         WRITE(*,'(A,2I5)') "ntime1, 2= ", ntime1, ntime2
!         WRITE(*,'(A,1P3E14.6)') "timefp + offset ", read_tms_double(1,ntime1), time+time_exp_offset, read_tms_double(1,ntime2)
!         WRITE(*,*) "weight ", weight
!      END IF

      END SUBROUTINE time_interpolation_cx
!------------------------------------
      SUBROUTINE time_interpolation_tms(time, ntime1, weight)

      IMPLICIT NONE
      REAL(rkind), intent(in):: time
      integer,intent(out):: ntime1
      integer:: ntime2
      REAL(rkind),intent(out):: weight ! y=(1-alpha)*y1 + alpha*y2
      REAL(rkind):: time_exp
      integer:: i

      time_exp= time + time_exp_offset
      DO i=1, nend_tms
         IF(read_tms_double(1,i).le.time_exp.and.time_exp.lt.read_tms_double(1,i+1))THEN
            ntime1=i
            ntime2=i+1
         END IF
      END DO
      
      weight=(time_exp-read_tms_double(1,ntime1))/(read_tms_double(1,ntime2)-read_tms_double(1,ntime1))

!      IF(NRANK.eq.0)THEN
!         WRITE(*,'(A,2I5)') "ntime1, 2= ", ntime1, ntime2
!         WRITE(*,'(A,1P3E14.6)') "timefp + offset ", read_tms_double(1,ntime1), time+time_exp_offset, read_tms_double(1,ntime2)
!         WRITE(*,*) "weight ", weight
!      END IF

      END SUBROUTINE time_interpolation_tms
!------------------------------------
! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL_EXP(PML,NR,NS)

      USE plprof
      USE libbes,ONLY: beseknx 
      implicit none
      integer :: NR, NS
      REAL(rkind) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      REAL(rkind) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      REAL(rkind):: FPMXWL_EXP

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

      IF(MODEL_DELTA_F(NS).eq.1)THEN
         RNFDL=(RN_TEMP(NR,NS)-RNS_DELF(NR,NS))/RNFD0L
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
         DKBSL=BESEKNX(2,Z)
         FACT=RNFDL*SQRT(THETA0L)/(4.D0*PI*RTFDL*DKBSL) &
              *RTFD0L
         EX=(1.D0-SQRT(1.D0+PML**2*THETA0L))/THETAL
         FPMXWL_EXP=FACT*EXP(EX)
      END IF

      RETURN
    END FUNCTION FPMXWL_EXP

  END MODULE fpreadeg
    
