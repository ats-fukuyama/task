!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE EG_READ

      USE fpcomm
      
!      double precision,dimension(:,:),pointer:: read_tms_double, read_cx_double
!      integer,dimension(:,:),pointer:: read_tms_int
!      double precision,dimension(5):: cte_fit
!      double precision,dimension(6):: cne_fit
!      integer:: nend_tms, nend_cx
!      integer,parameter:: NRMAX=32
!      double precision,dimension(32):: RM
!      double precision,dimension(32):: RNE_EXP, RTE_EXP

      contains
!-----------------------------------
      SUBROUTINE READ_EXP_DATA ! RNE_EXP and RTE_EXP are updated

      IMPLICIT NONE
      double precision:: weight
      integer:: ntime1, ntime2, NS, NR

!      CALL TEMPORARY_RHO_GRID 

      CALL READ_EXP_TMS
      CALL READ_EXP_CX 

      CALL time_interpolation(timefp, ntime1, ntime2, weight)
      CALL MAKE_PROF_FROM_TMS(ntime1,weight)

      DO NS=1, NSMAX
         DO NR=1, NRMAX
            RT_READ(NR,NS)=RTE_EXP(NR)
            RN_READ(NR,NS)=RNE_EXP(NR)

            RT_TEMP(NR,NS)=RTE_EXP(NR)
            RN_TEMP(NR,NS)=RNE_EXP(NR)
         END DO
      END DO

      END SUBROUTINE READ_EXP_DATA
!-----------------------------------
      SUBROUTINE READ_EXP_TMS

      IMPLICIT NONE
      integer:: i,j,k,m,nend1
      double precision,dimension(27,1000):: read_temporary_double
      integer,dimension(2,1000):: read_temporary_int

      open(22,file='./dat/tswpe@119489.dat',status='old')

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
      integer:: i,j,nend2

      allocate(read_cx_double(34,1000))

      open(22,file='./dat/cxswpi7@119489.dat',status='old')

      DO i=1,21
         READ(22,*) 
      END DO
      
      nend2=0
      DO i=1,1000
         READ(22,*,end=200) (read_cx_double(j,i),j=1,34)
         nend2=nend2+1
      END DO
      nend_cx=nend2
      
200   close(22)
      
!       WRITE(*,*) "END LINE", nend2

!       DO i=nend2-1,nend2
!          WRITE(*,'(1P34E14.6)') (read_cx_double(j,i),j=1,34)
!       END DO

      END SUBROUTINE READ_EXP_CX
!--------------------------------
      SUBROUTINE MAKE_PROF_FROM_TMS(ntime1,weight)

      IMPLICIT NONE
      integer,intent(in):: ntime1
      double precision,intent(in):: weight
      integer:: i,j,k
      integer:: NR
      double precision:: rte_ex, rne_ex, rho
      double precision,dimension(NRMAX,2):: prof_ne_temp, prof_te_temp
      
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
            
            prof_te_temp(NR,k)=rte_ex
            prof_ne_temp(NR,k)=rne_ex
         END DO
      END DO
      
! time interpolation
      DO NR=1, NRMAX
         RTE_EXP(NR)=(1.D0-weight)*prof_te_temp(NR,1) + weight*prof_te_temp(NR,2)
         RNE_EXP(NR)=( (1.D0-weight)*prof_ne_temp(NR,1) + weight*prof_ne_temp(NR,2) )*1.D-1 ! 10^20
      END DO
      RTE_EXP_EDGE= read_tms_double(24,i) 
      RNE_EXP_EDGE= read_tms_double(25,i)*1.D-1
       
      open(22,file='prof.dat')
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
      
      open(22,file='evol.dat')
      IF(NRANK.eq.0)THEN
         DO i=1, nend_tms 
            IF(read_tms_double(15,i).ne.0.and.read_tms_double(19,i).ne.0)THEN
               WRITE(22,'(3E14.6)') read_tms_double(1,i), read_tms_double(15,i), read_tms_double(19,i)
            END IF
         END DO
      END IF
      close(22)

      END SUBROUTINE MAKE_PROF_FROM_TMS
!------------------------------------
      SUBROUTINE time_interpolation(time, ntime1, ntime2, weight)

      IMPLICIT NONE
      double precision, intent(in):: time
      integer,intent(out):: ntime1, ntime2
      double precision,intent(out):: weight ! y=(1-alpha)*y1 + alpha*y2
      double precision:: time_exp
      integer:: i

      time_exp= time + time_exp_offset
      DO i=1, nend_tms
         IF(read_tms_double(1,i).le.time_exp.and.time_exp.lt.read_tms_double(1,i+1))THEN
            ntime1=i
            ntime2=i+1
         END IF
      END DO
      
      weight=(time_exp-read_tms_double(1,ntime1))/(read_tms_double(1,ntime2)-read_tms_double(1,ntime1))

      IF(NRANK.eq.0)THEN
         WRITE(*,'(A,2I5)') "ntime1, 2= ", ntime1, ntime2
         WRITE(*,'(A,1P3E14.6)') "timefp + offset ", read_tms_double(1,ntime1), time+time_exp_offset, read_tms_double(1,ntime2)
         WRITE(*,*) "weight ", weight
      END IF

      END SUBROUTINE time_interpolation
!------------------------------------
! ****************************************
!     MAXWELLIAN VELOCITY DISTRIBUTION
! ****************************************

      FUNCTION FPMXWL_EXP(PML,NR,NS)

      USE plprof
      USE libbes,ONLY: besekn 
      implicit none
      integer :: NR, NS
      real(kind8) :: PML,amfdl,aefdl,rnfd0l,rtfd0l,ptfd0l,rl,rhon
      real(kind8) :: rnfdl,rtfdl,fact,ex,theta0l,thetal,z,dkbsl
      TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
      real(kind8):: FPMXWL_EXP

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

      RNFDL=RN_TEMP(NR,NS)/RNFD0L
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
!------------------------------------
!      SUBROUTINE TEMPORARY_RHO_GRID
!
!      IMPLICIT NONE
!      double precision:: delr
!      INTEGER:: NR
!
!      delr=1.D0/NRMAX
!
!      DO NR=1, NRMAX
!         RM(NR)=(NR-1)*DELR + 0.5D0*DELR
!      END DO
!
!      END SUBROUTINE TEMPORARY_RHO_GRID
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE EG_READ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Program read_eg_file

       ! USE READ_EXP
       ! IMPLICIT NONE
       ! double precision:: timefp, weight
       ! integer:: ntime1, ntime2

       ! timefp=4.5D0
       ! CALL TEMPORARY_RHO_GRID 

       ! CALL READ_EXP_TMS
       ! CALL READ_EXP_CX 

       ! CALL time_interpolation(timefp, ntime1, ntime2, weight)
       ! CALL MAKE_PROF_FROM_TMS(ntime1,weight)
 
       ! END Program read_eg_file
