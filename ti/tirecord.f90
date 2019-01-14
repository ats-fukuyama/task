! tirecord.f90

MODULE tirecord

  USE ticomm,ONLY: rkind

  PRIVATE
  PUBLIC ti_snap,ti_record_ngt,ti_record_ngr

  INTEGER:: ngt_max,ngr_max,ngt_allocate_max,ngr_allocate_max
  REAL(rkind),ALLOCATABLE:: gt(:),gvt(:,:),gvta(:,:,:)
  REAL(rkind),ALLOCATABLE:: grt(:),gvrt(:,:,:),gvrta(:,:,:,:)

  PUBLIC ngt_max,ngr_max,ngt_allocate_max,ngr_allocate_max
  PUBLIC gt,gvt,gvta,grt,gvrt,gvrta

CONTAINS

  SUBROUTINE ti_snap
    USE ticomm
    IMPLICIT NONE

    IF(nrank.EQ.0) THEN
       WRITE(6,'(A,I5,1PE12.4,1PE12.4,2I5)') &
            '#NT,T,RD_LOOP,IC_LOOP,IC_MAT=', &
              NT,T,residual_loop_max,icount_loop_max,icount_mat_max
    END IF
    residual_loop_max=0.D0
    icount_loop_max=0
    icount_mat_max=0
    RETURN
  END SUBROUTINE ti_snap

  SUBROUTINE ti_integrate
    USE ticomm
    USE libmpi
    
    IMPLICIT NONE
    INTEGER:: nsa,nr,neq,nv,i
    REAL(rkind):: temp(nsmax)
    INTEGER:: loc(nsmax)

    DO nsa=1,nsa_max
       rnatot(nsa)=0.D0
       rnuatot(nsa)=0.D0
       rntatot(nsa)=0.D0
    END DO

    DO nr=nr_start,nr_end
       DO neq=1,neqmax
          i=(nr-1)*neqmax+neq
          IF(i.GE.istart.AND.i.LE.iend) THEN
             nsa=nsa_neq(neq)
             nv=nv_neq(neq)
             SELECT CASE(nv)
             CASE(1)
                rnatot(nsa) =rnatot(nsa) +rna(nsa,nr)*dvrho(nr)*dr
             CASE(2)
                rnuatot(nsa)=rnuatot(nsa)+rna(nsa,nr)*rua(nsa,nr)*dvrho(nr)*dr
             CASE(3)
                rntatot(nsa)=rntatot(nsa)+rna(nsa,nr)*rta(nsa,nr)*dvrho(nr)*dr
             END SELECT
          END IF
       END DO
    END DO

    CALL mtx_reduce_real8(rnatot,3, nsmax,temp,loc)
    IF(nrank.EQ.0) rnatot(1:nsmax)=temp(1:nsmax)
    CALL mtx_reduce_real8(rnuatot,3,nsmax,temp,loc)
    IF(nrank.EQ.0) rnuatot(1:nsmax)=temp(1:nsmax)
    CALL mtx_reduce_real8(rntatot,3,nsmax,temp,loc)
    IF(nrank.EQ.0) rntatot(1:nsmax)=temp(1:nsmax)

    IF(nrank.EQ.0) THEN
       DO nsa=1,nsa_max
          rnaave(nsa)=rnatot(nsa)/voltot
          IF(rnatot(nsa).NE.0.D0) THEN
             ruaave(nsa)=rnuatot(nsa)/rnatot(nsa)
             rtaave(nsa)=rntatot(nsa)/rnatot(nsa)
          ELSE
             ruaave(nsa)=0.D0
             rtaave(nsa)=0.D0
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE ti_integrate

  SUBROUTINE ti_allocate_ngt
    USE ticomm
    IMPLICIT NONE
    REAL(rkind),ALLOCATABLE:: gttemp(:),gvtatemp(:,:,:)
    
    IF(ngt_max.GT.0) THEN
       ALLOCATE(gttemp(ngt_max))
       ALLOCATE(gvtatemp(ngt_max,nsmax,9))
       gttemp(1:ngt_max)=gt(1:ngt_max)
       gvtatemp(1:ngt_max,1:nsmax,1:9)=gvta(1:ngt_max,1:nsmax,1:9)
    END IF

    IF(ALLOCATED(gt)) DEALLOCATE(gt)
    IF(ALLOCATED(gvta)) DEALLOCATE(gvta)

    ngt_allocate_max=ngt_allocate_max+ngt_allocate_step
    ALLOCATE(gt(ngt_allocate_max))
    ALLOCATE(gvta(ngt_allocate_max,nsmax,9))

    IF(ngt_max.GT.0) THEN
       gt(1:ngt_max)=gttemp(1:ngt_max)
       gvta(1:ngt_max,1:nsmax,1:9)=gvtatemp(1:ngt_max,1:nsmax,1:9)
       DEALLOCATE(gttemp,gvtatemp)
    END IF
    RETURN
  END SUBROUTINE ti_allocate_ngt

  SUBROUTINE ti_record_ngt
    USE ticomm
    IMPLICIT NONE
    INTEGER:: ngt,nsa

    CALL ti_integrate

    IF(ngt_max.GE.ngt_allocate_max) CALL ti_allocate_ngt
    ngt_max=ngt_max+1
    ngt=ngt_max

    IF(nrank.EQ.0) THEN
       gt(ngt)=t
       DO nsa=1,nsa_max
          gvta(ngt,nsa,1)=rna(nsa,1)
          gvta(ngt,nsa,2)=rua(nsa,1)
          gvta(ngt,nsa,3)=rta(nsa,1)
          gvta(ngt,nsa,4)=rna(nsa,nrmax)
          gvta(ngt,nsa,5)=rua(nsa,nrmax)
          gvta(ngt,nsa,6)=rta(nsa,nrmax)
          gvta(ngt,nsa,7)=rnaave(nsa)
          gvta(ngt,nsa,8)=ruaave(nsa)
          gvta(ngt,nsa,9)=rtaave(nsa)
       END DO
    END IF

    RETURN
  END SUBROUTINE ti_record_ngt

  SUBROUTINE ti_allocate_ngr
    USE ticomm
    IMPLICIT NONE
    REAL(rkind),ALLOCATABLE:: grttemp(:),gvrtatemp(:,:,:,:)
    
    IF(ngr_max.GT.0) THEN
       ALLOCATE(grttemp(ngr_max),gvrtatemp(nrmax,ngr_max,nsmax,3))
       grttemp(1:ngr_max)=grt(1:ngr_max)
       gvrtatemp(1:nrmax,1:ngr_max,1:nsa_max,1:3) &
            =gvrta(1:nrmax,1:ngr_max,1:nsa_max,1:3)
    END IF

    IF(ALLOCATED(grt)) DEALLOCATE(grt)
    IF(ALLOCATED(gvrta)) DEALLOCATE(gvrta)

    ngr_allocate_max=ngr_allocate_max+ngr_allocate_step
    ALLOCATE(grt(ngr_allocate_max),gvrta(nrmax,ngr_allocate_max,nsa_max,3))

    IF(ngr_max.GT.0) THEN
       grt(1:ngr_max)=grttemp(1:ngr_max)
       gvrta(1:nrmax,1:ngr_max,1:nsa_max,1:3) &
            =gvrtatemp(1:nrmax,1:ngr_max,1:nsa_max,1:3)
       DEALLOCATE(grttemp,gvrtatemp)
    END IF
    RETURN
  END SUBROUTINE ti_allocate_ngr

  SUBROUTINE ti_record_ngr
    USE ticomm
    IMPLICIT NONE
    INTEGER:: ngr,nsa,nr

    IF(ngr_max.GE.ngr_allocate_max) CALL ti_allocate_ngr
    ngr_max=ngr_max+1
    ngr=ngr_max

    IF(nrank.EQ.0) THEN
       grt(ngr)=t
       DO nsa=1,nsa_max
          DO nr=1,nrmax
             gvrta(nr,ngr,nsa,1)=rna(nsa,nr)
             gvrta(nr,ngr,nsa,2)=rua(nsa,nr)
             gvrta(nr,ngr,nsa,3)=rta(nsa,nr)
          END DO
       END DO
    END IF

    RETURN
  END SUBROUTINE ti_record_ngr
END MODULE tirecord
