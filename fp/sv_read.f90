!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE FP_READ_FIT

      USE fpcomm      

       integer,dimension(:,:),pointer:: I_FIT_temp
       double precision,dimension(:),pointer:: D_FIT_temp

       contains
!-----------------------------------
       SUBROUTINE READ_FIT3D

       IMPLICIT NONE
       integer:: i,j,k,l_max 
       integer:: iempty, n_p_theta
       double precision:: mass_ratio, vmax_mks, rempty

!       open(23,file='./dat/sv_fp_119489.txt',status='old')
       open(23,file=SV_FILE_NAME,status='old')

       READ(23,*) 
       READ(23,*) iempty, iempty, iempty, ntmax_fit
       allocate(time_grid_fit(ntmax_fit))
       READ(23,*) 
       READ(23,*) vmax_mks, mass_ratio
       READ(23,*) 
       READ(23,*) ! v grid
       READ(23,*) 
       READ(23,*) ! theta grid
       READ(23,*) 
       READ(23,*) ! rho grid
       READ(23,*) 
       READ(23,*) (time_grid_fit(i),i=1,ntmax_fit) 

       allocate(I_FIT_temp(4,ntmax_fit*5000))
       allocate(D_FIT_temp(ntmax_fit*5000))
       allocate(number_of_lines_fit(2,ntmax_fit))
       I_FIT_temp(:,:)=0
       D_FIT_temp(:)=0.D0
       number_of_lines_fit(:,:)=0.D0

       k=0
       DO i=1, ntmax_fit
          READ(23,*) 
          READ(23,*) rempty, n_p_theta
          READ(23,*) 
          IF(n_p_theta.ne.0)THEN
             number_of_lines_fit(1,i)=k+1
             DO j=1, n_p_theta
                k=k+1
                READ(23,*) I_FIT_temp(1,k), rempty, I_FIT_temp(2,k), rempty, I_FIT_temp(3,k), rempty, D_FIT_temp(k)
                I_FIT_temp(4,K)=i
             END DO
             number_of_lines_fit(2,i)=k
          END IF
!      I_FIT(i,k)= NP, NTH, NR, NT, k=line number
       END DO

       close(23)

!       WRITE(*,'(5I8,E14.6)') k, (I_FIT_temp(i,k),i=1,4), D_FIT_temp(k)
       l_max=k

       allocate(I_FIT(4,k))
       allocate(D_FIT(k))
       
       DO j=1,l_max
          DO i=1,4
             I_FIT(i,j)=I_FIT_temp(i,j)
          END DO
          D_FIT(j)=D_FIT_temp(j)
       END DO

       deallocate(I_FIT_temp,D_FIT_temp)

!       WRITE(*,'(5I8,E14.6)') l_max, (I_FIT(i,l_max),i=1,4), D_FIT(l_max)

       END SUBROUTINE READ_FIT3D
!------------------------------------
       SUBROUTINE time_interpolation_fit(time, ntime1, ntime2, weight)

       IMPLICIT NONE
       double precision, intent(in):: time
       integer,intent(out):: ntime1, ntime2
       double precision,intent(out):: weight ! y=(1-alpha)*y1 + alpha*y2
       integer:: i

       DO i=1, ntmax_fit
          IF(time_grid_fit(i).le.time.and.time.lt.time_grid_fit(i+1))THEN
             ntime1=i
             ntime2=i+1
          END IF
       END DO
       weight=(time-time_grid_fit(ntime1))/(time_grid_fit(ntime2)-time_grid_fit(ntime1))

!       IF(time.lt.time_grid_fit(1).or.time_grid_fit(ntmax_fit).lt.time)THEN
!          ntime1=1
!          ntime2=1
!          weight=0.D0
!       END IF

!       WRITE(*,'(A,2I5)') "ntime1, 2= ", ntime1, ntime2
!       WRITE(*,'(A,1P3E14.6)') "time ", time_grid_fit(ntime1), time, time_grid_fit(ntime2)
!       WRITE(*,*) "weight ", weight

       END SUBROUTINE time_interpolation_fit
!------------------------------------
       SUBROUTINE MAKE_SPPB_FIT(weight, ntime1, ntime2, NSA)

       IMPLICIT NONE
       double precision,intent(in):: weight
       integer,intent(in):: ntime1, ntime2, NSA
       integer:: i
       integer:: NTH, NP, NR
       double precision,dimension(NTHMAX,NPMAX,NRMAX):: SPPB_FIT_TEMP1, SPPB_FIT_TEMP2
       double precision:: FACT
!       double precision,dimension(NRMAX):: power1, power2, power3

!       FACT = VTFP0(NSA)**3/(RNFP0(NSA)*1.D20)*1.D6
       FACT = AMFP(NSA)/(RNFP0(NSA)*1.D20*PTFP0(NSA)**2)*1.D6
!       FACT = 1.D0
       SPPB_FIT_TEMP1(:,:,:)=0.D0
       SPPB_FIT_TEMP2(:,:,:)=0.D0

!       power1(:)=0.D0
!       power2(:)=0.D0
!       power3(:)=0.D0
!       WRITE(*,*) "time1 ", NRSTART, NPSTART, number_of_lines_fit(1,ntime1), number_of_lines_fit(2,ntime1)
!       WRITE(*,*) "time2 ", NRSTART, NPSTART, number_of_lines_fit(1,ntime2), number_of_lines_fit(2,ntime2)

       DO i=number_of_lines_fit(1,ntime1), number_of_lines_fit(2,ntime1)
          NP=I_FIT(1,i)
          NTH=I_FIT(2,i)
          NR=I_FIT(3,i)
          SPPB_FIT_TEMP1(NTH,NP,NR)=D_FIT(i)*FACT/(VOLP(NTH,NP,NSA)*PM(NP,NSA)**2*0.5D0)
!          power1(NR)=power1(NR)+D_FIT(i)
       END DO

       DO i=number_of_lines_fit(1,ntime2), number_of_lines_fit(2,ntime2)
          NP=I_FIT(1,i)
          NTH=I_FIT(2,i)
          NR=I_FIT(3,i)
          SPPB_FIT_TEMP2(NTH,NP,NR)=D_FIT(i)*FACT/(VOLP(NTH,NP,NSA)*PM(NP,NSA)**2*0.5D0)
!          power2(NR)=power2(NR)+D_FIT(i)
       END DO

       IF(timefp+time_exp_offset.lt.time_grid_fit(1).or.time_grid_fit(ntmax_fit).lt.timefp+time_exp_offset)THEN
          DO NR=NRSTART, NREND
             DO NP=NPSTART, NPEND
                DO NTH=1, NTHMAX
                   SPPB(NTH,NP,NR,NSA)=0.D0
                END DO
             END DO
          END DO
       ELSE
          DO NR=NRSTART, NREND
             DO NP=NPSTART, NPEND
                DO NTH=1, NTHMAX
                   SPPB(NTH,NP,NR,NSA)=(1.D0-weight)*SPPB_FIT_TEMP1(NTH,NP,NR) &
                        + weight*SPPB_FIT_TEMP2(NTH,NP,NR)
                END DO
             END DO
          END DO
       END IF
!       DO NR=1, NRMAX
!          power3(NR)=(1.D0-weight)*power1(NR) + weight*power2(NR)
!          IF(NRANK.eq.0) WRITE(*,*) NR, power3(NR)
!       END DO

       END SUBROUTINE MAKE_SPPB_FIT
!------------------------------------
       SUBROUTINE NBI_SOURCE_FIT3D

       IMPLICIT NONE
       
       double precision:: weight, time
       integer:: ntime1, ntime2, NSA

       NSA=NSSPB(1)
!       timefp=4.5D0
!       CALL READ_FIT3D
       time = timefp + time_exp_offset
       CALL time_interpolation_fit(time, ntime1, ntime2, weight)
       CALL MAKE_SPPB_FIT(weight, ntime1, ntime2, NSA)


       END SUBROUTINE NBI_SOURCE_FIT3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     END MODULE FP_READ_FIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Program read_fit_file
!
!       USE READ_FIT
!       IMPLICIT NONE
!       double precision:: timefp, weight
!       integer:: ntime1, ntime2
!
!       timefp=4.5D0
!       CALL TEMPORARY_RHO_GRID 
!
!       CALL READ_FIT3D!
!
!       CALL time_interpolation_fit(timefp, ntime1, ntime2, weight)
!       CALL MAKE_SPPB_FIT(weight, ntime1, ntime2)
! 
!       END Program read_fit_file
