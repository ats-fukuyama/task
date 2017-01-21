!     $Id: fpfile.f90,v 1.2 2010/03/19 09:27:53 nuga Exp $
!
!     ***********************************************************
!
!           SAVE VELOCITY DISTRIBUTION DATA
!
!     ***********************************************************
!
      MODULE fpfile

      USE fpcomm

      contains
!----------------------------------------------------------------
      SUBROUTINE FP_SAVE

      CALL FWOPEN(21,KNAMFP,0,MODEFW,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)
      WRITE(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         NSBA=NSB_NSA(NSA)
         WRITE(21) NS_NSA(NSA)
         WRITE(21) DELP(NSBA)
         WRITE(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         WRITE(21) (((FNS(NTH,NP,NR,NSBA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

  900 RETURN
      END SUBROUTINE FP_SAVE
!
!     ***********************************************************
!
!           LOAD VELOCITY DISTRIBUTION DATA
!
!     ***********************************************************
!
      SUBROUTINE FP_LOAD

      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)
      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) DELP(NSBA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  900 RETURN
      END SUBROUTINE FP_LOAD
!------------------------------------------      
      SUBROUTINE FP_SAVE2

      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NSB, NS, IERR, NSBA

      CALL FWOPEN(21,KNAMFP,0,MODEFW,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF

! required for FP
! the values that update in each step
      REWIND(21)
      WRITE(21) TIMEFP, NT_init

      DO NSA=1, NSAMAX
         DO NR=1, NRMAX
            DO NP=1, NPMAX
               WRITE(21) (FNS(NTH,NP,NR,NSA), NTH=1, NTHMAX) ! -> FNSP, FNSB
            END DO
         END DO
      END DO

      WRITE(21) ( (RN_IMPL(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      WRITE(21) ( (RT_IMPL(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      WRITE(21) (E1(NR),NR=1,NRMAX) ! -> EP

      IF(MODELD.ne.0)THEN
         WRITE(21) ((((WEIGHR_G(NTH,NP,NR,NSA), NTH=1,NTHMAX),NP=1,NPMAX),NR=1, NRMAX+1),NSA=1,NSAMAX)
      END IF

      IF(MODEL_DISRUPT.ne.0)THEN
         WRITE(21) (conduct_sp(NR), NR=1, NRMAX)
         WRITE(21) (RN_disrupt(NR), NR=1, NRMAX)
         WRITE(21) (RN_runaway(NR), NR=1, NRMAX)
         WRITE(21) (RN_drei(NR), NR=1, NRMAX)
         WRITE(21) (RJ_ohm(NR), NR=1, NRMAX)
         WRITE(21) (RJ_runaway(NR), NR=1, NRMAX)
         WRITE(21) (RJ_bs(NR), NR=1, NRMAX)
         WRITE(21) (RT_quench(NR),NR=1,NRMAX)
         WRITE(21) (previous_rate_G(NR), NR=1, NRMAX)
         WRITE(21) (previous_rate_p_G(NR), NR=1, NRMAX)
         WRITE(21) (RFP(NR), NR=1, NRMAX)
         WRITE(21) (RFP_ava(NR), NR=1, NRMAX)


         IF(MODEL_IMPURITY.ne.0)THEN
            WRITE(21) ((RN_MGI_G(NR,NSB), NR=1, NRMAX),NSB=1,NSBMAX)
            WRITE(21) (RN0_MGI(NSB), NSB=1, NSBMAX)
         END IF
      END IF

! required for other components
      WRITE(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         NSBA=NSB_NSA(NSA)
         WRITE(21) NS_NSA(NSA)
         WRITE(21) DELP(NSBA),PMAX(NSA)
         WRITE(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

  900 RETURN
      END SUBROUTINE FP_SAVE2
!------------------------------------------      
      SUBROUTINE FP_LOAD2

      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NSB, IERR, NSBA, NS

      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF

! required for FP
! the values that update in each step
      REWIND(21)
      READ(21) TIMEFP, NT_init

      DO NSA=1, NSAMAX
         DO NR=1, NRMAX
            DO NP=1, NPMAX
               READ(21) (FNS(NTH,NP,NR,NSA), NTH=1, NTHMAX) ! -> FNSP, FNSB
            END DO
         END DO
      END DO
      
      READ(21) ( (RN_IMPL(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      READ(21) ( (RT_IMPL(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      READ(21) (E1(NR),NR=1,NRMAX) ! -> EP

      IF(MODELD.ne.0)THEN
         READ(21) ((((WEIGHR_G(NTH,NP,NR,NSA), NTH=1,NTHMAX),NP=1,NPMAX),NR=1, NRMAX+1),NSA=1,NSAMAX)
      END IF

      IF(MODEL_DISRUPT.ne.0)THEN
         READ(21) (conduct_sp(NR), NR=1, NRMAX)
         READ(21) (RN_disrupt(NR), NR=1, NRMAX)
         READ(21) (RN_runaway(NR), NR=1, NRMAX)
         READ(21) (RN_drei(NR), NR=1, NRMAX)
         READ(21) (RJ_ohm(NR), NR=1, NRMAX)
         READ(21) (RJ_runaway(NR), NR=1, NRMAX)
         READ(21) (RJ_bs(NR), NR=1, NRMAX)
         READ(21) (RT_quench(NR),NR=1,NRMAX)
         READ(21) (previous_rate_G(NR), NR=1, NRMAX)
         READ(21) (previous_rate_p_G(NR), NR=1, NRMAX)
         READ(21) (RFP(NR), NR=1, NRMAX)
         READ(21) (RFP_ava(NR), NR=1, NRMAX)

         IF(MODEL_IMPURITY.ne.0)THEN
            READ(21) ((RN_MGI_G(NR,NSB), NR=1, NRMAX),NSB=1,NSBMAX)
            READ(21) (RN0_MGI(NSB), NSB=1, NSBMAX)
         END IF
      END IF


! required for other compponents
      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         NSBA=NSB_NSA(NSA)
         READ(21) NS_NSA(NSA)
         READ(21) DELP(NSBA),PMAX(NSA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  900 RETURN
      END SUBROUTINE FP_LOAD2
!------------------------------------------      
      SUBROUTINE FP_PRE_LOAD

      USE fpprep
      USE libmpi
      USE fpmpi
      IMPLICIT NONE
      integer:: NSW, N, NSA, ierr, i
      real:: gut1, gut2
      integer,dimension(1:6):: idata
      integer,dimension(6*nsize):: idata2

      CALL GUTIME(gut1) 
      CALL fp_comm_setup
      CALL fp_allocate
      call fp_allocate_ntg1
      call fp_allocate_ntg2

      CALL mtx_set_communicator(comm_nr)
      allocate(MTXLEN(nsize),MTXPOS(nsize))
      CALL mtx_set_communicator(comm_nsanr)
      allocate(SAVLEN(nsize)) 
      allocate(SAVPOS(nsize,NSAEND-NSASTART+1)) 
      CALL mtx_reset_communicator
      allocate(Rank_Partition_Data(1:6,0:nsize-1))

      CALL mtx_reset_communicator

      idata(1)=NPSTARTW
      idata(2)=NPENDWM
      idata(3)=NRSTARTW
      idata(4)=NRENDWM
      idata(5)=NSASTART
      idata(6)=NSAEND

      CALL mtx_gather_integer(idata,6,idata2)

      IF(NRANK.eq.0)THEN
         DO N=0, nsize-1
            DO I=1,6
               Rank_Partition_Data(I,N)=idata2(I+N*6)
            END DO
         END DO
      END IF

      CALL mtx_set_communicator(comm_nr)

      CALL mtx_allgather1_integer(nrend-nrstart+1,mtxlen)
      CALL mtx_allgather1_integer(nrstart-1,mtxpos)

      CALL mtx_set_communicator(comm_nsanr)
      CALL mtx_allgather1_integer(nrend-nrstart+1,savlen)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL mtx_allgather1_integer((nsa-1)*NRMAX+NRSTART,savpos(1:nsize,N))
      END DO
      CALL mtx_reset_communicator

      CALL fp_set_nsa_nsb
!     ----- create meches -----
!      WRITE(6,*) "START MESH"
      CALL fp_mesh(ierr)
!      WRITE(6,*) "END MESH"
!     ----- Initialize diffusion coef. -----
      call FPCINI
!     ----- set parameters for target species -----
      CALL fp_set_normalize_param
      IF(MODEL_DISRUPT.ne.0)THEN
         CALL set_initial_disrupt_param
         CALL set_post_disrupt_Clog_f
         call mtx_broadcast_real8(POST_tau_ta0_f,nsamax)
      END IF

      END SUBROUTINE FP_PRE_LOAD
!------------------------------------------      
      SUBROUTINE FP_POST_LOAD

      USE fpprep
      USE libmpi
      USE fpmpi
      IMPLICIT NONE
      integer:: NSW, N, NSA, ierr, i, NR, NSB

      CALL mtx_reset_communicator
      CALL scatter_fns_to_fns0

      CALL mtx_set_communicator(comm_nsa)
      CALL update_fnsb
      CALL mtx_reset_communicator

      CALL Bcast_loaded_data
      IF(MODEL_DISRUPT.ne.0)THEN
         CALL set_post_disrupt_Clog
      END IF

      DO NR=NRSTART,NREND
         EP(NR)=E1(NR)
         EM(NR)=0.D0
      END DO

      CALL FNSP_INIT_EDGE
      CALL fp_continue(ierr)
      CALL fp_set_initial_value_from_f


      END SUBROUTINE FP_POST_LOAD
!------------------------------------------      
      SUBROUTINE Bcast_loaded_data

      USE libmpi
      IMPLICIT NONE
      INTEGER,DIMENSION(99):: idata
      real(8),DIMENSION(99):: rdata
      integer:: NR, NSB, NSA, NTH, NP

      idata(1)=NT_init
      CALL mtx_broadcast_integer(idata,1)
      NT_init = idata(1)

      DO NR=1,NRMAX
         rdata(NR)=E1(NR)
      END DO
      rdata(NRMAX+1)=TIMEFP
      CALL mtx_broadcast_real8(rdata,NRMAX+1)
      DO NR=1,NRMAX
         E1(NR)=rdata(NR)
      END DO
      TIMEFP=rdata(NRMAX+1)

      CALL mtx_broadcast_real8(RN_IMPL,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RT_IMPL,NRMAX*NSAMAX)

      IF(MODELD.ne.0)THEN
         CALL mtx_broadcast_real8(WEIGHR_G,NTHMAX*NPMAX*(NRMAX+1)*NSAMAX)
         DO NSA=NSASTART, NSAEND
            DO NR=NRSTART, NREND+1
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     WEIGHR(NTH,NP,NR,NSA)=WEIGHR_G(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
         END DO
      END IF

      IF(MODEL_DISRUPT.ne.0)THEN
         DO NR=1,NRMAX
            rdata(NR)=conduct_sp(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            conduct_sp(NR)=rdata(NR)
         END DO      
!
         DO NR=1,NRMAX
            rdata(NR)=RN_disrupt(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RN_disrupt(NR)=rdata(NR)
         END DO      
!
         DO NR=1,NRMAX
            rdata(NR)=RN_runaway(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RN_runaway(NR)=rdata(NR)
         END DO      
!
         DO NR=1,NRMAX
            rdata(NR)=RN_drei(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RN_drei(NR)=rdata(NR)
         END DO      
!
         DO NR=1,NRMAX
            rdata(NR)=RJ_ohm(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RJ_ohm(NR)=rdata(NR)
         END DO      
!
         DO NR=1,NRMAX
            rdata(NR)=RJ_runaway(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RJ_runaway(NR)=rdata(NR)
         END DO      
!
         DO NR=1,NRMAX
            rdata(NR)=RJ_bs(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RJ_bs(NR)=rdata(NR)
         END DO
!
         DO NR=1,NRMAX
            rdata(NR)=RT_quench(NR)
         END DO
         CALL mtx_broadcast_real8(rdata,NRMAX)
         DO NR=1,NRMAX
            RT_quench(NR)=rdata(NR)
         END DO

         CALL mtx_broadcast_real8(previous_rate_G,NRMAX)
         CALL mtx_broadcast_real8(previous_rate_p_G,NRMAX)
         DO NR=NRSTART, NREND
            previous_rate(NR)=previous_rate_G(NR)
            previous_rate_p(NR)=previous_rate_p_G(NR)
         END DO
         IF(MODEL_IMPURITY.ne.0)THEN
            CALL mtx_broadcast_real8(RN0_MGI,NSBMAX)
            CALL mtx_broadcast_real8(RN_MGI_G,NRMAX*NSBMAX)
            DO NSB=1,NSBMAX
               DO NR=NRSTART, NREND
                  RN_MGI(NR,NSB)=RN_MGI_G(NR,NSB)
               END DO
            END DO
         END IF
      END IF

      END SUBROUTINE Bcast_loaded_data
!------------------------------------------      
! GATHER DATA to RANK 0
      SUBROUTINE FP_PRE_SAVE

      USE libmpi
      IMPLICIT NONE
      double precision,dimension(NRMAX)::temp
      double precision,dimension(NRSTART:NREND)::temp_l
      double precision,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND):: dsend
      double precision,dimension(nthmax,npmax,nrmax):: drecv
      double precision,dimension(NTHMAX,NPMAX):: dsend_max, drecv_max
      integer,dimension(NTHMAX,NPMAX)::vloc_max

      double precision,dimension(nthmax,npmax,nrmax+1,nsastart:nsaend):: temp_l2
      double precision,dimension(nthmax,npstart:npend):: dsend2
      double precision,dimension(nthmax,npmax):: drecv2
      integer,dimension(NRMAX,NSBMAX)::vloc
      integer,dimension(NTHMAX*(NPEND-NPSTART+1))::vloc2
      integer:: nsend
      INTEGER:: NR, NSB, NSA, NTH, NP, dest, source, tag, nn, Isum, NRE
!!!!!!!!!!!!

      IF(MODELD.ne.0)THEN
         CALL mtx_set_communicator(comm_nrnp)
         nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)

         DO NSA=NSASTART, NSAEND
            DO NR=NRSTART, NREND
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     dsend(nth,np,nr)=WEIGHR(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
            CALL mtx_allgather_real8(dsend,nsend,drecv)
            DO NR=1, NRMAX
               DO NP=1, NPMAX
                  DO NTH=1, NTHMAX
                     temp_l2(nth,np,nr,nsa)=drecv(NTH,NP,NR)
                  END DO
               END DO
            END DO
            IF(NREND.eq.NRMAX)THEN
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     dsend_max(nth,np)=WEIGHR(NTH,NP,NRMAX+1,NSA)
                  END DO
               END DO
            ELSE
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     dsend_max(nth,np)=0.D0
                  END DO
               END DO
            END IF
            CALL mtx_allreduce_real8(dsend_max,NTHMAX*NPMAX,3,drecv_max,vloc_max)
            DO NP=1, NPMAX
               DO NTH=1, NTHMAX
                  temp_l2(nth,np,nrmax+1,nsa)=drecv_max(NTH,NP)
               END DO
            END DO
         END DO
         CALL mtx_set_communicator(comm_nsa)
         nsend=NTHMAX*NPMAX*(NRMAX+1)*(NSAEND-NSASTART+1) 
         CALL mtx_gather_real8(temp_l2,nsend,WEIGHR_G)
         CALL mtx_reset_communicator 
      END IF ! END MODELD

!!!!!!!!!!!!!
      IF(MODEL_DISRUPT.ne.0)THEN
         previous_rate_g(:)=0.D0
         previous_rate_p_g(:)=0.D0
! NS, NR
         IF(MODEL_IMPURITY.ne.0)THEN
            RN_MGI_G(:,:)=0.D0
            CALL mtx_set_communicator(comm_nr)
            DO NSB=1,NSBMAX
               DO NR=NRSTART, NREND
                  temp_l(NR)=RN_MGI(NR,NSB)
               END DO
               CALL mtx_gather_real8(temp_l,NREND-NRSTART+1,temp)
               DO NR=1, NRMAX
                  RN_MGI_G(NR,NSB)=temp(NR)
               END DO
            END DO

            CALL mtx_reset_communicator
         END IF
! NR
         CALL mtx_set_communicator(comm_nr)
         CALL mtx_gather_real8(previous_rate,NREND-NRSTART+1,previous_rate_g)
         CALL mtx_gather_real8(previous_rate_p,NREND-NRSTART+1,previous_rate_p_g)
      END IF

      CALL mtx_reset_communicator

      END SUBROUTINE FP_PRE_SAVE
!------------------------------------------      
      SUBROUTINE OPEN_EVOLVE_DATA_OUTPUT

      USE libmpi      
      IMPLICIT NONE

      IF(MODEL_DISRUPT.ne.0.and.NRANK.eq.0)THEN
         OPEN(9,file="f1_1.dat") 
         open(10,file='time_evol.dat') 
         open(11,file='efield_e1.dat') 
         open(12,file='dndt.dat') 
         open(13,file='radial.dat') 
         open(14,file='nth-re.dat')
         open(15,file='re_pitch.dat')
         open(18,file='efield_ref.dat')
         open(23,file='collision_time.dat')
      END IF

      open(19,file='p-T_bulk.dat')
!      open(20,file='DCPP.dat')
!      open(16,file='err_message_for_RT_BULK.dat')

      END SUBROUTINE OPEN_EVOLVE_DATA_OUTPUT
!------------------------------------------      
      SUBROUTINE CLOSE_EVOLVE_DATA_OUTPUT

      USE libmpi      
      IMPLICIT NONE

      IF(MODEL_DISRUPT.ne.0.and.NRANK.eq.0)THEN
         close(9)
         close(10)
         close(11)
         close(12)
         close(13)
         close(14)
         close(15)
         close(18)  
         close(23)  
      END IF
      close(19)
!      close(20)  
!      close(16)  

      END SUBROUTINE CLOSE_EVOLVE_DATA_OUTPUT
!------------------------------------------      

      END MODULE fpfile
