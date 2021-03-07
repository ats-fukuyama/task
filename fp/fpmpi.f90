      MODULE FPMPI

      USE libmpi
      USE fpcomm

      CONTAINS

      SUBROUTINE fp_gatherv_real8_sav(vsend,nscnt,vreturn,n,nsa)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: nscnt
      REAL(rkind),DIMENSION(nrstart:nrend,NSAMAX),INTENT(IN):: vsend
      double precision,dimension(nrstart:nrend):: vtemp
!      REAL(rkind),DIMENSION(nrmax*nsamax):: vrecv
      REAL(rkind),DIMENSION(nrmax*N_partition_s):: vrecv
      REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT):: vreturn
      INTEGER:: n, nsa, nn, ns, nr, nsw, nse
      integer,dimension(nsize):: idisp

      DO nn=1,nsize
         idisp(nn)=savpos(nn,n)-1
      END DO
      DO NR=NRSTART,NREND
         vtemp(NR)=vsend(nr,nsa)
      END DO
!      call mtx_gather_real8(vtemp,nscnt,vrecv)
      call mtx_gather_real8(vtemp,nscnt,vrecv)

      IF(NRANK.eq.0)THEN
         nsw = NSAEND-NSASTART+1
         DO NS=1,n_partition_s
            NSE=(NS-1)*NSW + N
            DO NR=1,NRMAX
               vreturn(NR,NSE)=vrecv( (NS-1)*NRMAX+NR )
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE fp_gatherv_real8_sav
!-----
      SUBROUTINE update_fnsb

      IMPLICIT NONE
      integer:: nsend, nth, np, nr, nsa, nsb
      double precision,dimension(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend):: dsend
      double precision,dimension(nthmax,npstart:npend,nrstart:nrend,nsamax):: drecv

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  dsend(nth,np,nr,nsa)=FNSP(nth,np,nr,nsa)
               END DO
            END DO
         END DO
      END DO

      nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)*(NSAEND-NSASTART+1)
!      CALL mtx_allgather_real8(dsend,nsend,FNSB(1:NTHMAX,NPSTART,NRSTART,1)) 
      CALL mtx_allgather_real8(dsend,nsend,drecv(1,NPSTART,NRSTART,1)) 


      DO NSA=1, NSAMAX
         NSB=NSB_NSA(NSA)
         IF(NSB.ne.0)THEN
            DO NR=NRSTART,NREND
               DO NP=NPSTART,NPEND
                  DO NTH=1,NTHMAX
                     FNSB(NTH,NP,NR,NSB)=drecv(nth,np,nr,nsa)
                  END DO
               END DO
            END DO
         END IF
      END DO

      END SUBROUTINE update_fnsb
!-----
      SUBROUTINE update_fns

      USE libmtx
      IMPLICIT NONE
      integer:: nsend, nth, np, nr, nsa, nsb, nsw, nswi
!      double precision,dimension(nthmax,npstart:npend,nrstart:nrend):: dsend
!      double precision,dimension(nthmax,npmax,nrmax,n_partition_s):: drecv
      double precision,dimension(nthmax,npstart:npend,nrstart:nrend):: dsend
      double precision,dimension(nthmax,npmax,nrmax):: drecv
      double precision,dimension(nthmax,npmax,nrmax,nsastart:nsaend):: dsend2

      CALL mtx_set_communicator(comm_nrnp)
      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  dsend(nth,np,nr)=fnsp(nth,np,nr,nsa)
               END DO
            END DO
         END DO
         nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)
         CALL mtx_allgather_real8(dsend,nsend,drecv) 
         DO NR=1,NRMAX
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  dsend2(nth,np,nr,nsa)=drecv(nth,np,nr)
               END DO
            END DO
         END DO
      END DO
      CALL mtx_set_communicator(comm_nsa)
      nsend=NTHMAX*NPMAX*NRMAX*(NSAEND-NSASTART+1)
      CALL mtx_gather_real8(dsend2,nsend,FNS) 
      CALL mtx_reset_communicator 

      END SUBROUTINE update_fns
!-----

      SUBROUTINE source_allreduce(array)

      IMPLICIT NONE
      DOUBLE PRECISION,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX), &
           INTENT(INOUT):: array
      DOUBLE PRECISION,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX):: &
           sendbuf, recvbuf
      INTEGER,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX):: &
           lloc
      INTEGER:: ierr, ncount

         ncount = NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)*NSAMAX
         sendbuf(:,:,:,:)=0.D0
         sendbuf(:,:,:,:)=array(:,:,:,:)

         CALL mtx_allreduce_real8(sendbuf, ncount, 3, recvbuf,lloc)

         array(:,:,:,:) = recvbuf(:,:,:,:)

      RETURN

      END SUBROUTINE source_allreduce
!-----
      SUBROUTINE p_theta_integration(vdata)

      IMPLICIT NONE
      double precision,intent(inout):: vdata
      double precision:: vrecv
      integer::vloc

      CALL mtx_allreduce1_real8(vdata,3,vrecv,vloc) 
      vdata = vrecv

      END SUBROUTINE p_theta_integration
!-----
      SUBROUTINE fpl_comm(vdata,vrecv)

      IMPLICIT NONE
      double precision,dimension(NPSTART:NPEND),intent(in)::vdata
      double precision,dimension(NPMAX),intent(out)::vrecv
      integer:: ndata

      ndata = (NPEND-NPSTART+1)
      CALL mtx_allgather_real8(vdata,ndata,vrecv)

      END SUBROUTINE fpl_comm
!-----
      SUBROUTINE shadow_comm_np(NR,NSA)

      IMPLICIT NONE
      double precision,dimension(nthmax)::sendbuf
      double precision,dimension(nthmax)::recvbuf
      integer,intent(in):: NR,NSA
      integer:: sendcount, recvcount, dest, source, nth

      DO NTH=1,NTHMAX
         sendbuf(nth)=FNS0(NTH,NPEND,NR,NSA)
         recvbuf(nth)=0.D0
      END DO
      
      sendcount=nthmax
      recvcount=sendcount
      dest=nrank+1
      source=nrank-1
      CALL mtx_sendrecv_real8(sendbuf,sendcount,dest, &
                              recvbuf,recvcount,source) 

      IF(NPSTART.ne.NPSTARTW)THEN
         DO NTH=1,NTHMAX
            FNS0(NTH,NPSTARTW,NR,NSA)=recvbuf(nth)
         END DO
      END IF
!============
      DO NTH=1,NTHMAX
         sendbuf(nth)=FNS0(NTH,NPSTART,NR,NSA)
         recvbuf(nth)=0.D0
      END DO
      
      dest=nrank-1
      source=nrank+1

      CALL mtx_sendrecv_real8(sendbuf,sendcount,dest, &
                              recvbuf,recvcount,source)

      IF(NPEND.ne.NPENDWM)THEN
         DO NTH=1,NTHMAX
            FNS0(NTH,NPENDWM,NR,NSA)=recvbuf(nth)
         END DO
      END IF

      END SUBROUTINE shadow_comm_np
!-----
      SUBROUTINE shadow_comm_nr(NSA)

      IMPLICIT NONE
      double precision,dimension(nthmax*(npendwm-npstartw+1))::sendbuf
      double precision,dimension(nthmax*(npendwm-npstartw+1))::recvbuf
      integer,intent(in):: NSA
      integer:: sendcount, recvcount, dest, source, nth, np, NM

      DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX
            NM=NTH+NTHMAX*(NP-NPSTARTW)
            sendbuf(NM)=FNS0(NTH,NP,NREND,NSA)
            recvbuf(NM)=0.D0
         END DO
      END DO

      sendcount=nthmax*(npendwm-npstartw+1)
      recvcount=sendcount
      dest=nrank+1
      source=nrank-1

      CALL mtx_sendrecv_real8(sendbuf,sendcount,dest, &
                              recvbuf,recvcount,source) 

      IF(NRSTART.ne.NRSTARTW)THEN
         DO NP=NPSTARTW, NPENDWM
            DO NTH=1,NTHMAX
               NM=NTH+NTHMAX*(NP-NPSTARTW)
               FNS0(NTH,NP,NRSTARTW,NSA)=recvbuf(NM)
            END DO
         END DO
      END IF
!===
      DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX
            NM=NTH+NTHMAX*(NP-NPSTARTW)
            sendbuf(NM)=FNS0(NTH,NP,NRSTART,NSA)
            recvbuf(NM)=0.D0
         END DO
      END DO

      sendcount=nthmax*(npendwm-npstartw+1)
      recvcount=sendcount
      dest=nrank-1
      source=nrank+1

      CALL mtx_sendrecv_real8(sendbuf,sendcount,dest, &
                              recvbuf,recvcount,source) 

      IF(NREND.ne.NRENDWM)THEN
         DO NP=NPSTARTW, NPENDWM
            DO NTH=1,NTHMAX
               NM=NTH+NTHMAX*(NP-NPSTARTW)
               FNS0(NTH,NP,NRENDWM,NSA)=recvbuf(NM)
            END DO
         END DO
      END IF

      END SUBROUTINE shadow_comm_nr
!-----------------------------------------------------------
      SUBROUTINE scatter_fns_to_fns0

      IMPLICIT NONE
      integer:: NTH,NP,NR,NSA, I, J
      integer:: NPS, NPE, NRS, NRE, NSAS, NSAE
      integer:: sendcount, recvcount
      integer:: dest, source, tag
      double precision,dimension(:,:,:,:),allocatable:: sendbuf, recvbuf

      CALL mtx_reset_communicator 

      IF(NRANK.eq.0)THEN
         DO I=0,nsize-1
            NPS=Rank_Partition_Data(1,I)
            NPE=Rank_Partition_Data(2,I)
            NRS=Rank_Partition_Data(3,I)
            NRE=Rank_Partition_Data(4,I)
            NSAS=Rank_Partition_Data(5,I)
            NSAE=Rank_Partition_Data(6,I)
            allocate(sendbuf(NTHMAX, 1-NPS+NPE, 1-NRS+NRE, 1-NSAS+NSAE))
            
            DO NSA=NSAS,NSAE
               DO NR=NRS,NRE
                  DO NP=NPS,NPE
                     DO NTH=1, NTHMAX
                        sendbuf(NTH, NP-NPS+1, NR-NRS+1, NSA-NSAS+1)=FNS(NTH,NP,NR,NSA)
                     END DO
                  END DO
               END DO
            END DO
            
            IF(I.eq.0)THEN
               DO NSA=NSASTART,NSAEND
                  DO NR=NRSTARTW,NRENDWM
                     DO NP=NPSTARTW,NPENDWM
                        DO NTH=1,NTHMAX
                           FNS0(NTH,NP,NR,NSA)=sendbuf(NTH,NP-NPSTARTW+1,NR-NRSTARTW+1,NSA-NSASTART+1)
                        END DO
                     END DO
                  END DO
               END DO
            ELSE
               dest = I
               tag=dest
               sendcount=nthmax*(NPE-NPS+1)*(NRE-NRS+1)*(NSAE-NSAS+1)
               CALL mtx_send_real8(sendbuf,sendcount,dest,tag)
            END IF
            deallocate(sendbuf)
         END DO
      ELSE
         source = 0
         recvcount=nthmax*(NPENDWM-NPSTARTW+1)*(NRENDWM-NRSTARTW+1)*(NSAEND-NSASTART+1)
         allocate(recvbuf(NTHMAX, NPSTARTW:NPENDWM, NRSTARTW:NRENDWM, NSASTART:NSAEND))

         tag=NRANK
         CALL mtx_recv_real8(recvbuf,recvcount,source,tag)
         DO NSA=NSASTART,NSAEND
            DO NR=NRSTARTW,NRENDWM
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX
                     FNS0(NTH,NP,NR,NSA)=recvbuf(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
         END DO
         deallocate(recvbuf)
      END IF

      FNSP(:,:,:,:)=FNS0(:,:,:,:)

      END SUBROUTINE scatter_fns_to_fns0
!-----------------------------------------------------------
      END MODULE FPMPI
