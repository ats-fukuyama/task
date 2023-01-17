!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE fpreadeg

      USE fpcomm
      USE libmpi

      PUBLIC READ_EXP_DATA
      PUBLIC MAKE_EXP_PROF
      PUBLIC MAKE_Tene_PROF_POLYNOMIAL
      PUBLIC MAKE_QLM_QLG_FROM_kspdiag

      PRIVATE
      integer:: TMS_DimSize, TMS_ValNo, CX_DimSize, CX_ValNo
      CHARACTER(len=100),dimension(:),allocatable:: TMS_ValName, CX_ValName
      double precision,dimension(:,:),allocatable:: TMS_DATA_ALL, CX_DATA_ALL
      double precision,dimension(:),allocatable:: TMS_TIME, CX_TIME
      double precision,dimension(:,:),allocatable:: TMS_te_poly_coef, CX_ti_poly_coef
      double precision,dimension(:,:),allocatable:: TMS_ne_poly_coef
      double precision,dimension(:),allocatable:: TMS_tedge, TMS_nedge, CX_ti_edge

      contains
!-----------------------------------
      SUBROUTINE READ_EXP_DATA ! EGDATA OF TMS and CX are READ

      IMPLICIT NONE

!      CALL READ_EXP_TMS
      CALL MAKE_TMS_DATA
      IF(MODEL_EX_READ_Tn.eq.2)THEN
!      CALL READ_EXP_CX 
         CALL MAKE_CXSWPI7_DATA
      END IF

      END SUBROUTINE READ_EXP_DATA
!-----------------------------------
      SUBROUTINE MAKE_EXP_PROF(time) ! RNE_EXP and RTE_EXP are updated

      IMPLICIT NONE
      double precision,intent(in):: time
      double precision:: weight, zeff, CN
      integer:: ntime1, NS, NR, NSA
      double precision,dimension(NSMAX-1,NSMAX-1):: Aji
      double precision,dimension(NSMAX-1):: bk,xsolved      

      IF(MODEL_EX_READ_Tn.ne.0)THEN
!         CALL time_interpolation_tms(time, ntime1, weight)
!         CALL MAKE_PROF_FROM_TMS(ntime1,weight)
         CALL time_interpolation_v2(time,TMS_TIME,TMS_DimSize, ntime1,weight)
         CALL MAKE_PROF_FROM_TMS_v2(ntime1,weight)
      END IF

      IF(MODEL_EX_READ_Tn.eq.2)THEN
         CALL time_interpolation_v2(time,CX_TIME,CX_DimSize, ntime1,weight)
!         CALL MAKE_PROF_FROM_CX(ntime1,weight)
         CALL MAKE_PROF_FROM_CX_v2(ntime1,weight)
      END IF

      DO NR=1, NRMAX
         RT_READ(NR,1)=RTE_EXP(NR)
         RN_READ(NR,1)=RNE_EXP(NR)
               
         RT_BULK(NR,1)=RTE_EXP(NR)
         RN_TEMP(NR,1)=RNE_EXP(NR)
         RN_BULK(NR,1)=RNE_EXP(NR)
      END DO
      IF(MODEL_DELTA_F_NI_RATIO.eq.1)THEN
         DO NS=2, NSMAX
            IF(MODEL_EX_READ_Tn.eq.1)THEN
               DO NR=1, NRMAX
                  RT_READ(NR,NS)=RTE_EXP(NR)
                  RN_READ(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
                  
                  RT_BULK(NR,NS)=RTE_EXP(NR)
                  RN_TEMP(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
               END DO
            ELSEIF(MODEL_EX_READ_Tn.eq.2)THEN
               DO NR=1, NRMAX
                  RT_READ(NR,NS)=RTI_EXP(NR)
                  RN_READ(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
                  
                  RT_BULK(NR,NS)=RTI_EXP(NR)
                  RN_TEMP(NR,NS)=RNE_EXP(NR)*NI_RATIO(NS)
               END DO
            END IF
         END DO
      ELSEIF(MODEL_DELTA_F_NI_RATIO.eq.2)THEN
         DO NR=1,NRMAX
            CALL make_matrix_bulk_dens(NR,Aji,bk)
            CALL solve_bulk_density_matrix(Aji,bk,xsolved)
            DO NS=2,NSMAX
               RN_BULK(NR,NS)=xsolved(NS-1)
               RN_TEMP(NR,NS)=RN_BULK(NR,NS)+RNS_DELF(NR,NS)
            END DO
         END DO
         IF(MODEL_EX_READ_Tn.eq.1)THEN
            DO NS=2,NSMAX
               DO NR=1,NRMAX
                  RT_READ(NR,NS)=RTE_EXP(NR)*Ti_Te_ratio(NS)
                  RT_BULK(NR,NS)=RTE_EXP(NR)*Ti_Te_ratio(NS)
!                  RT_READ(NR,NS)=RTE_EXP(NR)*0.5D0 ! temporal
               END DO
            END DO
         ELSEIF(MODEL_EX_READ_Tn.eq.2)THEN
            DO NS=2,NSMAX
               DO NR=1,NRMAX
                  RT_READ(NR,NS)=RTI_EXP(NR)
                  RT_BULK(NR,NS)=RTI_EXP(NR)
               END DO
            END DO
         END IF
      END IF

      DO NS=1, NSMAX
         IF(NS.eq.1)THEN
            DO NR=1,NRMAX
               RT_BULK(NR,NS)=RTE_EXP(NR)
            END DO
         ELSE
            IF(MODEL_EX_READ_Tn.eq.1)THEN
               DO NR=1, NRMAX
                  RT_BULK(NR,NS)=RTE_EXP(NR)*Ti_Te_ratio(NS)
               END DO
            ELSEIF(MODEL_EX_READ_Tn.eq.2)THEN
               DO NR=1, NRMAX
                  RT_BULK(NR,NS)=RTI_EXP(NR)
               END DO
            END IF
         END IF
      END DO

      !temporal
!      WRITE(*,'(4E14.6)') TIMEFP, RTE_EXP(1), RT_BULK(1,1), Ti_Te_ratio(NS_NSA(1))
!      WRITE(35,'(4E14.6)') TIMEFP, RN_TEMP(1,1), RN_TEMP(1,2), RN_TEMP(1,3)
      zeff=0.D0
      CN=0.D0
      DO NS=2, NSMAX
         zeff=zeff+PZ(NS)**2*RN_TEMP(1,NS)/RN_TEMP(1,1)
         CN = CN + PZ(NS)*RN_TEMP(1,NS)/RN_TEMP(1,1)
      END DO

      WRITE(37,'(99E14.6)') TIMEFP, (RN_BULK(1,NS), NS=1, NSMAX), &
           (RN_TEMP(1,NS), NS=1, NSMAX), zeff, CN

      END SUBROUTINE MAKE_EXP_PROF
!-----------------------------------
      SUBROUTINE READ_EXP_TMS

      IMPLICIT NONE
      integer:: i,j,k,m,nend1
      double precision,dimension(27,1000):: read_temporary_double
      integer,dimension(2,1000):: read_temporary_int

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
      double precision,dimension(34,1000):: read_cx_temp, read_cx_temp2

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
      double precision,intent(in):: weight
      integer:: i,j,k
      integer:: NR
      double precision:: rti_ex, rho
      double precision,dimension(NRMAX,2):: prof_ti_temp
      
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
      double precision, intent(in):: time
      integer,intent(out):: ntime1
      integer:: ntime2
      double precision,intent(out):: weight ! y=(1-alpha)*y1 + alpha*y2
      double precision:: time_exp
      integer:: i

      time_exp= time + time_exp_offset
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
      double precision, intent(in):: time
      integer,intent(out):: ntime1
      integer:: ntime2
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

!      IF(NRANK.eq.0)THEN
!         WRITE(*,'(A,2I5)') "ntime1, 2= ", ntime1, ntime2
!         WRITE(*,'(A,1P3E14.6)') "timefp + offset ", read_tms_double(1,ntime1), time+time_exp_offset, read_tms_double(1,ntime2)
!         WRITE(*,*) "weight ", weight
!      END IF

      END SUBROUTINE time_interpolation_tms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------
      SUBROUTINE make_matrix_bulk_dens(NR,Aji,bk) ! Ax=b, x means bulk density 

      IMPLICIT NONE
      integer,intent(in):: NR
      double precision,dimension(NSMAX-1,NSMAX-1),intent(out):: Aji
      double precision,dimension(NSMAX-1),intent(out):: bk
      double precision:: RSUMZ, RSUMZ2
      integer:: i_row,j_column,NS,judge,judge_zeff,imax, i, j
      integer,dimension(NSMAX-1):: ID_NS

      DO NS=2,NSMAX
         IF(PA(NS).eq.1)THEN ! Proton
            ID_NS(NS-1)=1
         ELSEIF(PA(NS).eq.2.and.PZ(NS).eq.1)THEN ! D
            ID_NS(NS-1)=2
         ELSEIF(PA(NS).eq.3.and.PZ(NS).eq.1)THEN ! T
            ID_NS(NS-1)=3
         ELSEIF(PA(NS).eq.4.and.PZ(NS).eq.2)THEN ! He4
            ID_NS(NS-1)=4
         ELSEIF(PA(NS).eq.12.and.PZ(NS).eq.6)THEN ! C
            ID_NS(NS-1)=6
         ELSE
            ID_NS(NS-1)=0
         END IF
      END DO

      Aji(:,:)=0.D0
      bk(:)=0.D0

      judge=0
      judge_zeff=0
      DO NS=2,NSMAX
         IF(PZ(NS).ge.judge)THEN
            judge=PZ(NS)
         END IF
      END DO
!      IF(judge.ge.2) judge_zeff=1
      IF(judge.ge.3) judge_zeff=1 

      i_row=0
! 1  eq. of charge neutrality
      i_row=i_row+1
      RSUMZ=0.D0
      RSUMZ2=0.D0
      DO j_column=1,NSMAX-1
         Aji(i_row,j_column) = PZ(j_column+1)
      END DO

      IF(MODEL_DELTA_F_CN.eq.0)THEN
         RSUMZ=0.D0
      ELSEIF(MODEL_DELTA_F_CN.eq.1)THEN
         DO NS=2,NSMAX
            RSUMZ =RSUMZ +RNS_DELF(NR,NS)*PZ(NS)
            RSUMZ2=RSUMZ2+RNS_DELF(NR,NS)*PZ(NS)**2
         END DO
      END IF
      bk(i_row) = RN_TEMP(NR,1)-RSUMZ

! 2  eq. of Zeff
      IF(judge_zeff.eq.1)THEN
         i_row=i_row+1
         DO j_column=1,NSMAX-1
            Aji(i_row,j_column) = PZ(j_column+1)**2
         END DO
         bk(i_row) = given_zeff*RN_TEMP(NR,1)-RSUMZ2
      END IF

! 3  eq. of DH_RATIO
      IF(DH_RATIO.ne.0.and.DH_RATIO.ne.1.D0)THEN
         i_row=i_row+1
         DO j_column=1,NSMAX-1
            IF(ID_NS(j_column).eq.1)THEN
               Aji(i_row,j_column) = -DH_RATIO
            ELSEIF(ID_NS(j_column).eq.2)THEN
               Aji(i_row,j_column) = 1.D0 - DH_RATIO
            ELSE
               Aji(i_row,j_column) = 0.D0
            END IF
         END DO
         bk(i_row)=0.D0
      END IF

! 4  eq. of HHe_RATIO ! Here, H means P and D. 
      IF(HHe_RATIO.ne.0.and.HHe_RATIO.ne.1.D0)THEN
         i_row=i_row+1
         DO j_column=1,NSMAX-1
            IF(ID_NS(j_column).eq.1)THEN
               Aji(i_row,j_column) = 1.D0 - HHe_RATIO
            ELSEIF(ID_NS(j_column).eq.2)THEN
               Aji(i_row,j_column) = 1.D0 - HHe_RATIO
            ELSEIF(ID_NS(j_column).eq.4)THEN
               Aji(i_row,j_column) = - HHe_RATIO
            ELSE
               Aji(i_row,j_column) = 0.D0
            END IF
         END DO
         bk(i_row)=0.D0
      END IF

! 5  eq. of DT_RATIO
      IF(DT_RATIO.ne.0.and.DT_RATIO.ne.1.D0)THEN
         i_row=i_row+1
         DO j_column=1,NSMAX-1
            IF(ID_NS(j_column).eq.1)THEN
               Aji(i_row,j_column) = 1.D0 - HHe_RATIO
            ELSEIF(ID_NS(j_column).eq.4)THEN
               Aji(i_row,j_column) = - HHe_RATIO
            ELSE
               Aji(i_row,j_column) = 0.D0
            END IF
         END DO
         bk(i_row)=0.D0
      END IF
      imax=i_row

!      IF(NRANK.eq.0)THEN
!         DO i=1,imax
!            WRITE(*,'(A,I2,10E14.6)') "A",i, (Aji(i,j),j=1,NSMAX-1)
!         END DO
!         WRITE(*,*) " "
!         DO i=1,imax
!            WRITE(*,'(A,I2,10E14.6)') "B",i, bk(i)
!         END DO
!      END IF

      END SUBROUTINE make_matrix_bulk_dens
!------------------------------------
      SUBROUTINE solve_bulk_density_matrix(Aji,bk,xsolved)

      IMPLICIT NONE
      double precision,dimension(NSMAX-1,NSMAX-1),intent(in)::Aji
      double precision,dimension(NSMAX-1),intent(in)::bk
      double precision,dimension(NSMAX-1),intent(out)::xsolved
      integer::NS,IERR
      integer,dimension(NSMAX-1):: IPIV

      CALL DGESV(NSMAX-1,1,Aji,NSMAX-1,IPIV,bk,NSMAX-1,IERR)
      IF(IERR.ne.0)THEN
         WRITE(*,'(A,I4)') "DGSEV ERR, IERR=", IERR
         WRITE(*,'(A,I4)') "DGSEV ERR, NRANK=", NRANK
      END IF

      DO NS=1,NSMAX-1
         xsolved(NS)=bk(NS)
      END DO

      END SUBROUTINE solve_bulk_density_matrix
!------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE EXTRACT_NUMBER_FROM_EG(STR_IN,INT_OUT)

      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)::STR_IN
      INTEGER,INTENT(OUT):: INT_OUT
      CHARACTER(LEN=1000):: trimmed
      INTEGER:: length, i_loc

      trimmed = trim(STR_IN)
      length = len(trimmed)
      i_loc = index(trimmed,"=")
      trimmed = trim(trimmed(i_loc+1:length))

      READ(trimmed,*) INT_OUT

      END SUBROUTINE EXTRACT_NUMBER_FROM_EG
!------------------------------------
      SUBROUTINE EXTRACT_NAME_FROM_EG(STR_IN,SIZE_IN,STR_OUT)
      
      IMPLICIT NONE
      INTEGER,INTENT(IN):: SIZE_IN
      CHARACTER(*),INTENT(IN)::STR_IN
      CHARACTER(len=100),dimension(SIZE_IN),INTENT(OUT)::STR_OUT
      CHARACTER(len=1000):: trimmed1, trimmed2
      INTEGER:: i_loc1, i_loc2, length, NoV

      trimmed1 = trim(STR_IN)
      trimmed1 = trim(trimmed1)

      DO NoV = 1, SIZE_IN
         length = len(trimmed1)
         i_loc1 = index(trimmed1,"'")
         trimmed2 = trimmed1(i_loc1+1:length)
         i_loc2 = index(trimmed2,"'")
         STR_OUT(NoV) = trimmed2(1:i_loc2-1)
!         IF(NRANK.eq.0) WRITE(*,'(A,4I5,A)') "TEST EXTRACT ", Nov, length, i_loc1, i_loc2, STR_OUT(NoV)
         length = len(trimmed2)
         trimmed1=trimmed2(i_loc2+1:length)
      END DO

      END SUBROUTINE EXTRACT_NAME_FROM_EG
!------------------------------------
      SUBROUTINE EXTRACT_DATA_FROM_EG(SIZE_IN, STRING_IN, OUT_VAL)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: SIZE_IN
      character(LEN=1000),intent(in):: STRING_IN
      character(LEN=100),dimension(SIZE_IN):: c_val
      double precision,dimension(SIZE_IN),intent(out):: OUT_VAL
      CHARACTER(len=1000):: trimmed1, trimmed2
      INTEGER:: NA, length, i_loc

      trimmed1 = trim(STRING_IN)
      DO NA=1, SIZE_IN-1
         length=len(trimmed1)
         i_loc = index(trimmed1,',')
         c_val(NA) = trimmed1(1:i_loc-1)
         trimmed2 = trimmed1(i_loc+1:length)
         trimmed1 = trim(trimmed2)
      END DO
      c_val(SIZE_IN) = trimmed1

      DO NA=1, SIZE_IN
         READ(c_val(NA),*) OUT_VAL(NA)
      END DO

      END SUBROUTINE EXTRACT_DATA_FROM_EG
!------------------------------------
      SUBROUTINE READ_EGFILE_1D(FILENAME,FILEPATH, DimSize,ValNo,ValName,DATA_ALL)

      IMPLICIT NONE
      character(len=1000),intent(in):: FILENAME, FILEPATH
      integer,intent(out):: DimSize, ValNo
      CHARACTER(len=100),dimension(:),allocatable,intent(out):: ValName
      double precision,dimension(:,:),allocatable,intent(out):: DATA_ALL
      character(len=1000):: path, name
      character(len=1000):: string
      integer:: ierr, nt, i_column
      double precision,dimension(:),allocatable:: TIME_EG
      double precision,dimension(:),allocatable:: data_row
      
      name=trim(FILENAME)//'@'//trim(SHOT_NUMBER)//'.dat'
      path=trim(FILEPATH)//trim(name)

!      IF(NRANK.eq.0) WRITE(*,*) "EG_PATH=", path

      OPEN(22,file=path,status='old',IOSTAT=ierr)
      IF(ierr.ne.0)THEN
         WRITE(*,'(A,I5,A)') "READ EG ERR ", nrank, path
         call mtx_abort(ierr)         
      END IF

      string='initialize'
      DO WHILE(index(string,'[Data]').eq.0.and.index(string,'[data]').eq.0.and.index(string,'[DATA]').eq.0)
         READ(22,'(A)') string
         IF(index(string,'DimSize').ne.0)THEN
            CALL EXTRACT_NUMBER_FROM_EG(string,DimSize) ! number of time grid
            allocate(TIME_EG(DimSize))
         END IF
         IF(index(string,'ValNo').ne.0)THEN
            CALL EXTRACT_NUMBER_FROM_EG(string,ValNo) ! number of data array-1
            allocate(ValName(ValNo))
         END IF
         IF(index(string,'ValName').ne.0)THEN
            CALL EXTRACT_NAME_FROM_EG(string,ValNo,ValName)
         END IF
      END DO

      allocate(data_row(ValNo+1))
      allocate(DATA_ALL(ValNo+1,DimSize))
      DO nt = 1, DimSize
         READ(22,'(A)') string
         CALL EXTRACT_DATA_FROM_EG(ValNo+1, string, data_row)
         DO i_column=1, ValNo+1
            DATA_ALL(i_column,nt)=data_row(i_column)
         END DO
      END DO

      CLOSE(22)
      END SUBROUTINE READ_EGFILE_1D
!------------------------------------
      SUBROUTINE EXTRACT_ne_te_POSITION_FROM_TSWPE(ValName,ValNo,INT_OUT)

      IMPLICIT NONE
      integer,intent(in):: ValNo
      character(len=100),dimension(ValNo),intent(in):: ValName
      integer,dimension(11),intent(out):: INT_OUT
! 1-4: cte0-cte4, 5-9: cne0-cne5, 10-11: te_edge, ne_edge
      integer:: i

      DO i=2, ValNo+1
!         IF(NRANK.eq.0) WRITE(*,'(I5,A,A)') i,", TEST ValName= ", ValName(i-1)
         IF(index(ValName(i-1),'cte0').ne.0) INT_OUT(1)=i
         IF(index(ValName(i-1),'cte2').ne.0) INT_OUT(2)=i
         IF(index(ValName(i-1),'cte4').ne.0) INT_OUT(3)=i
         IF(index(ValName(i-1),'cte6').ne.0) INT_OUT(4)=i

         IF(index(ValName(i-1),'cne0').ne.0) INT_OUT(5)=i
         IF(index(ValName(i-1),'cne2').ne.0) INT_OUT(6)=i
         IF(index(ValName(i-1),'cne4').ne.0) INT_OUT(7)=i
         IF(index(ValName(i-1),'cne6').ne.0) INT_OUT(8)=i
         IF(index(ValName(i-1),'cne8').ne.0) INT_OUT(9)=i

         IF(index(ValName(i-1),'tedge').ne.0) INT_OUT(10)=i
         IF(index(ValName(i-1),'nedge').ne.0) INT_OUT(11)=i
      END DO

      END SUBROUTINE EXTRACT_ne_te_POSITION_FROM_TSWPE
!------------------------------------
      SUBROUTINE EXTRACT_ti_POSITION_FROM_CXSWPI7(ValName,ValNo,INT_OUT)

      IMPLICIT NONE
      integer,intent(in):: ValNo
      character(len=100),dimension(ValNo),intent(in):: ValName
      integer,dimension(5),intent(out):: INT_OUT
! 1-4: cti0-cti4, 5: ti_edge
      integer:: i

      DO i=2, ValNo+1
         IF(index(ValName(i-1),'cti0').ne.0) INT_OUT(1)=i
         IF(index(ValName(i-1),'cti2').ne.0) INT_OUT(2)=i
         IF(index(ValName(i-1),'cti4').ne.0) INT_OUT(3)=i
         IF(index(ValName(i-1),'cti6').ne.0) INT_OUT(4)=i

         IF(index(ValName(i-1),'ti_edge').ne.0) INT_OUT(5)=i
      END DO

      END SUBROUTINE EXTRACT_ti_POSITION_FROM_CXSWPI7
!------------------------------------
      SUBROUTINE MAKE_TMS_DATA

      IMPLICIT NONE
      character(len=1000):: FILENAME, FILEPATH
      integer,dimension(11):: I_column
      integer:: nt

      FILENAME=trim(EG_NAME_TMS)
      FILEPATH=trim(EG_PATH)

      CALL READ_EGFILE_1D(FILENAME,FILEPATH, TMS_DimSize,TMS_ValNo,TMS_ValName,TMS_DATA_ALL)
      CALL EXTRACT_ne_te_POSITION_FROM_TSWPE(TMS_ValName,TMS_ValNo,I_column)

      allocate(TMS_TIME(TMS_DimSize))
      allocate(TMS_te_poly_coef(4,TMS_DimSize))
      allocate(TMS_ne_poly_coef(5,TMS_DimSize))
      allocate(TMS_tedge(TMS_DimSize))
      allocate(TMS_nedge(TMS_DimSize))

      DO nt=1, TMS_DimSize
         TMS_TIME(nt) = TMS_DATA_ALL(1,nt)
         TMS_te_poly_coef(1,nt) = TMS_DATA_ALL(I_column(1),nt)
         TMS_te_poly_coef(2,nt) = TMS_DATA_ALL(I_column(2),nt)
         TMS_te_poly_coef(3,nt) = TMS_DATA_ALL(I_column(3),nt)
         TMS_te_poly_coef(4,nt) = TMS_DATA_ALL(I_column(4),nt)

         TMS_ne_poly_coef(1,nt) = TMS_DATA_ALL(I_column(5),nt)
         TMS_ne_poly_coef(2,nt) = TMS_DATA_ALL(I_column(6),nt)
         TMS_ne_poly_coef(3,nt) = TMS_DATA_ALL(I_column(7),nt)
         TMS_ne_poly_coef(4,nt) = TMS_DATA_ALL(I_column(8),nt)
         TMS_ne_poly_coef(5,nt) = TMS_DATA_ALL(I_column(9),nt)

         TMS_tedge(nt) = TMS_DATA_ALL(I_column(10),nt)
         TMS_nedge(nt) = TMS_DATA_ALL(I_column(11),nt)
      END DO
!      IF(NRANK.eq.0)THEN
!         WRITE(*,'(A,11I5)') "TEST, I_column", (I_column(nt),nt=1,11)
!         DO nt=1, TMS_DimSize
!            WRITE(*,'(A,I5,5E14.6)') "TEST_CORE_EDGE, ", nt, TMS_DATA_ALL(1,nt), &
!                 TMS_te_poly_coef(1,nt), TMS_ne_poly_coef(1,nt), TMS_tedge(nt), TMS_nedge(nt)
!         END DO
!      END IF

      END SUBROUTINE MAKE_TMS_DATA
!------------------------------------
      SUBROUTINE MAKE_CXSWPI7_DATA

      IMPLICIT NONE
      character(len=1000):: FILENAME, FILEPATH
      integer,dimension(5):: I_column
      integer:: nt, I_sum, I_sum2

      FILENAME=trim(EG_NAME_CX)
      FILEPATH=trim(EG_PATH)

      CALL READ_EGFILE_1D(FILENAME,FILEPATH, CX_DimSize,CX_ValNo,CX_ValName,CX_DATA_ALL)
      CALL EXTRACT_ti_POSITION_FROM_CXSWPI7(CX_ValName,CX_ValNo,I_column)

      I_sum = 0
      DO nt=1, CX_DimSize
         IF(CX_DATA_ALL(I_column(1),nt).gt.0.D0)THEN
            I_sum = I_sum + 1 
         END IF
      END DO

      allocate(CX_TIME(I_sum))
      allocate(CX_ti_poly_coef(4,I_sum))
      allocate(CX_ti_edge(I_sum)) 

      CX_TIME(:)=0.D0
      CX_ti_poly_coef(:,:)=0.D0
      CX_ti_edge(:)=0.D0

      I_sum2=0
      DO nt=1, CX_DimSize
         IF(CX_DATA_ALL(I_column(1),nt).gt.0.D0)THEN
            I_sum2=I_sum2+1
            CX_TIME(I_sum2) = CX_DATA_ALL(1,nt)
            CX_ti_poly_coef(1,I_sum2) = CX_DATA_ALL(I_column(1),nt)
            CX_ti_poly_coef(2,I_sum2) = CX_DATA_ALL(I_column(2),nt)
            CX_ti_poly_coef(3,I_sum2) = CX_DATA_ALL(I_column(3),nt)
            CX_ti_poly_coef(4,I_sum2) = CX_DATA_ALL(I_column(4),nt)
            
            CX_ti_edge(I_sum2) = CX_DATA_ALL(I_column(5),nt)         
         END IF
      END DO
      CX_DimSize = I_sum

      END SUBROUTINE MAKE_CXSWPI7_DATA
!------------------------------------
      SUBROUTINE time_interpolation_v2(time, TIME_MEASURE, SIZE_IN, ntime1, weight)

      IMPLICIT NONE
      double precision, intent(in):: time ! timefp
      integer,intent(in):: SIZE_IN
      double precision,dimension(SIZE_IN),intent(in):: TIME_Measure
      integer,intent(out):: ntime1
      integer:: ntime2
      double precision,intent(out):: weight ! y=(1-alpha)*y1 + alpha*y2
      double precision:: time_exp
      integer:: nt

      time_exp= time + time_exp_offset
      DO nt=1, SIZE_IN-1
         IF(TIME_MEASURE(nt).le.time_exp.and.time_exp.lt.TIME_MEASURE(nt+1))THEN
            ntime1=nt
            ntime2=nt+1
         END IF
      END DO
      
      weight=(time_exp-TIME_MEASURE(ntime1))/(TIME_MEASURE(ntime2)-TIME_MEASURE(ntime1))

      END SUBROUTINE time_interpolation_v2
!------------------------------------
      SUBROUTINE MAKE_PROF_FROM_TMS_v2(ntime1,weight)

      IMPLICIT NONE
      integer,intent(in):: ntime1
      double precision,intent(in):: weight
      integer:: i,k
      integer:: NR
      double precision:: rte_ex, rne_ex, rho
      double precision,dimension(NRMAX,2):: prof_ne_temp, prof_te_temp

      
      DO k=1, 2
         i=k-1+ntime1
         DO NR=1, NRMAX
            rho = RM(NR)
            rte_ex = TMS_te_poly_coef(1,i)       +TMS_te_poly_coef(2,i)*rho**2 &
                    +TMS_te_poly_coef(3,i)*rho**4+TMS_te_poly_coef(4,i)*rho**6
            rne_ex = TMS_ne_poly_coef(1,i)       +TMS_ne_poly_coef(2,i)*rho**2 &
                    +TMS_ne_poly_coef(3,i)*rho**4+TMS_ne_poly_coef(4,i)*rho**6 &
                    +TMS_ne_poly_coef(5,i)*rho**8
            IF(rte_ex.gt.0.D0)THEN
               prof_te_temp(NR,k)=rte_ex
            ELSE
               prof_te_temp(NR,k)=1.D-1
            END IF
            IF(rne_ex.gt.0.D0)THEN
               prof_ne_temp(NR,k)=rne_ex
            ELSE
               prof_ne_temp(NR,k)=1.D-1
            END IF
         END DO
      END DO
      
! time interpolation
      DO NR=1, NRMAX
         RTE_EXP(NR)=(1.D0-weight)*prof_te_temp(NR,1) + weight*prof_te_temp(NR,2)
         RNE_EXP(NR)=( (1.D0-weight)*prof_ne_temp(NR,1) + weight*prof_ne_temp(NR,2) )*1.D-1 ! 10^20
         IF(RTE_EXP(NR).le.0) WRITE(*,'(A,I5,2E14.6)') "Negative RTE_EXP ", NR, RTE_EXP(NR), weight
         IF(RNE_EXP(NR).le.0) WRITE(*,'(A,I5,2E14.6)') "Negative RNE_EXP ", NR, RNE_EXP(NR), weight
      END DO
      RTE_EXP_EDGE= (1.D0-weight)*TMS_tedge(ntime1) + weight*TMS_tedge(ntime1+1)
      RNE_EXP_EDGE=((1.D0-weight)*TMS_nedge(ntime1) + weight*TMS_nedge(ntime1+1) )*1.D-1      
      IF(RTE_EXP_EDGE.le.0.D0) RTE_EXP_EDGE=1.D-1
      IF(RNE_EXP_EDGE.le.0.D0) RNE_EXP_EDGE=1.D-2
       
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
      
      END SUBROUTINE MAKE_PROF_FROM_TMS_v2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------
      SUBROUTINE MAKE_PROF_FROM_CX_v2(ntime1,weight)

      IMPLICIT NONE
      integer,intent(in):: ntime1
      double precision,intent(in):: weight
      integer:: i,k
      integer:: NR, ntime2
      double precision:: rti_ex1, rho
      double precision,dimension(NRMAX,2):: prof_ti_temp

      DO k=1, 2
         i=k-1+ntime1
         DO NR=1, NRMAX
            rho = RM(NR)
            rti_ex1 = CX_ti_poly_coef(1,i)       +CX_ti_poly_coef(2,i)*rho**2 &
                     +CX_ti_poly_coef(3,i)*rho**4+CX_ti_poly_coef(4,i)*rho**6
         IF(rti_ex1.gt.0.D0)THEN
            prof_ti_temp(NR,k)=rti_ex1
         ELSE
            prof_ti_temp(NR,k)=1.D-1
         END IF
      END DO
      END DO
      
! time interpolation
      DO NR=1, NRMAX
         RTI_EXP(NR)=(1.D0-weight)*prof_ti_temp(NR,1) + weight*prof_ti_temp(NR,2)
         IF(RTI_EXP(NR).le.0) WRITE(*,'(A,2I5,2E14.6)') "Negative RTI_EXP ", ntime1, NR, RTI_EXP(NR), weight
      END DO
      RTI_EXP_EDGE= (1.D0-weight)*CX_ti_edge(ntime1) + weight*CX_ti_edge(ntime1+1)
      IF(RTI_EXP_EDGE.le.0.D0) RTE_EXP_EDGE=1.D-1
       
      open(22,file='prof_cx.dat')
      IF(NRANK.eq.0)THEN
         WRITE(22,'(A,3E14.6)') "# TIME, TI_EDGE, ", TIMEFP, RTI_EXP_EDGE
         DO NR=1, NRMAX
            WRITE(22,'(I5,4E14.6)') NR, RM(NR), prof_ti_temp(NR,1), RTI_EXP(NR), prof_ti_temp(NR,2)
         END DO
         WRITE(22,*) " "
         WRITE(22,*) " "
      END IF
      close(22)
      
      END SUBROUTINE MAKE_PROF_FROM_CX_v2
!------------------------------------
      SUBROUTINE MAKE_Tene_PROF_POLYNOMIAL

      IMPLICIT NONE
      INTEGER:: NR, NS, NSA, i
      
      DO NR=1, NRMAX
         RT_INIT(NR,1)=te_poly(1) + te_poly(2)*RM(NR)**2 + &
              te_poly(3)*RM(NR)**4 + te_poly(4)*RM(NR)**6
         

         RN_INIT(NR,1)=(ne_poly(1) + ne_poly(2)*RM(NR)**2 + &
              ne_poly(3)*RM(NR)**4 + &
              ne_poly(4)*RM(NR)**6 + ne_poly(5)*RM(NR)**8)*1.D-1
         
         RN_BULK(NR,1)=RN_TEMP(NR,1)
      END DO
      
      DO NS=2, NSMAX
         DO NR=1, NRMAX
            RT_INIT(NR,NS)=RT_INIT(NR,1)
            RN_INIT(NR,NS)=RN_TEMP(NR,1)*NI_RATIO(NS)
               
            RN_BULK(NR,NS)=RN_INIT(NR,NS)
         END DO
      END DO

      END SUBROUTINE MAKE_Tene_PROF_POLYNOMIAL
!------------------------------------
      SUBROUTINE READ_kspdiag(FILEPATH, FILETIME, nrho_max, RM_ksp, iota_bar_ksp)

      IMPLICIT NONE
      character(len=1000),intent(in):: FILEPATH, FILETIME
      integer,intent(out):: nrho_max
      double precision,dimension(100),intent(out):: RM_ksp, iota_bar_ksp

      character(len=1000):: path, name
      character(len=1000):: string, trimmed
      integer:: ierr, i_pos_s, i_pos_e, nrho

      name='kspdiag_data_'//trim(SHOT_NUMBER)//'t'//trim(FILETIME)//'.flx'
      path=trim(FILEPATH)//trim(name)
      RM_ksp(:)=0.D0
      iota_bar_ksp(:)=0.D0

!      IF(NRANK.eq.0) WRITE(*,*) "VMEC_PATH=", path

      OPEN(22,file=path,status='old',IOSTAT=ierr)
      IF(ierr.ne.0)THEN
         WRITE(*,'(A,I5,A)') "READ EG ERR ", nrank, path
         call mtx_abort(ierr)         
      END IF

      string='initialize'
      DO WHILE(index(string,'nrho').eq.0)
         READ(22,'(A)') string
      END DO
      i_pos_s = index(string,'=') + 1
      i_pos_e = index(string,'modnum') - 2

      trimmed = trim(string(i_pos_s:i_pos_e))
      READ(trimmed,*) nrho_max

!      IF(NRANK.eq.0) WRITE(*,*) "TEST ", nrho_max

      string='initialize'
      DO nrho=1, nrho_max
         DO WHILE(index(string,'rho=').eq.0)
            READ(22,'(A)') string
         END DO
!         IF(NRANK.eq.0) WRITE(*,*) "TEST string", string

         i_pos_s = index(string,'rho=') + 4
         i_pos_e = index(string,'1/q') - 2
         trimmed = trim(string(i_pos_s:i_pos_e))
!         IF(NRANK.eq.0) WRITE(*,*) "TEST RM_ksp", trimmed
         READ(trimmed,*) RM_ksp(nrho)

         i_pos_s = index(string,'1/q=') + 4
         i_pos_e = index(string,'Vp') - 2
         trimmed = trim(string(i_pos_s:i_pos_e))
!         IF(NRANK.eq.0) WRITE(*,*) "TEST iotabar_ksp", trimmed
         READ(trimmed,*) iota_bar_ksp(nrho)
         string='initialize'
      END DO

      CLOSE(22)

      END SUBROUTINE READ_kspdiag
!------------------------------------
      SUBROUTINE MAKE_QLM_QLG_FROM_kspdiag(FILEPATH, FILETIME)

      IMPLICIT NONE
      character(len=1000),intent(in):: FILEPATH, FILETIME
      integer:: nrho_max, NR, nrho, nrho_m, nrho_p
      double precision,dimension(100):: RM_ksp, iota_bar_ksp
      double precision:: a1, a2, iota_bar

      CALL READ_kspdiag(FILEPATH, FILETIME, nrho_max, RM_ksp, iota_bar_ksp)

      DO nrho=1, nrho_max
         IF(NRANK.eq.0) WRITE(*,'(I4,3E14.6)') nrho, RM_ksp(nrho), iota_bar_ksp(nrho), &
              1.D0/SQRT(1.D0 + (RM_ksp(nrho)*RA*iota_bar_ksp(nrho)/RR)**2)
      END DO

      DO NR=1, NRMAX
         DO nrho=1, nrho_max-1
            IF(RM_ksp(nrho).lt.RM(NR).and.RM(NR).le.RM_ksp(nrho+1))THEN
               nrho_m = nrho
               nrho_p = nrho + 1
            END IF
         END DO
         a1= (RM(NR)-RM_ksp(nrho_m))/(RM_ksp(nrho_p)-RM_ksp(nrho_m))
         a2= (RM_ksp(nrho_p)-RM(NR))/(RM_ksp(nrho_p)-RM_ksp(nrho_m))
         iota_bar = a2*iota_bar_ksp(nrho_m) + a1*iota_bar_ksp(nrho_p)
         QLM(NR) = 1.D0 / iota_bar
         QLM_INIT(NR) = QLM(NR)
      END DO

      DO NR=1, NRMAX+1
         DO nrho=1, nrho_max-1
            IF(RM_ksp(nrho).lt.RG(NR).and.RG(NR).le.RM_ksp(nrho+1))THEN
               nrho_m = nrho
               nrho_p = nrho + 1
            END IF
         END DO
         a1= (RG(NR)-RM_ksp(nrho_m))/(RM_ksp(nrho_p)-RM_ksp(nrho_m))
         a2= (RM_ksp(nrho_p)-RG(NR))/(RM_ksp(nrho_p)-RM_ksp(nrho_m))
         iota_bar = a2*iota_bar_ksp(nrho_m) + a1*iota_bar_ksp(nrho_p)
         QLG(NR) = 1.D0 / iota_bar
         QLG_INIT(NR) = QLG(NR)
      END DO

      END SUBROUTINE MAKE_QLM_QLG_FROM_kspdiag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END MODULE fpreadeg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
