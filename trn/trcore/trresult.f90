MODULE trresult

  USE trcomm, ONLY: rkind,ikind

  PUBLIC tr_calc_global,tr_status,tr_save_ngt
  PRIVATE

CONTAINS

! ***** calculate global values *****

  SUBROUTINE tr_calc_global

    RETURN
  END SUBROUTINE tr_calc_global

! ***** simple status report *****

  SUBROUTINE tr_status

    USE trcomm, ONLY : nsamax,t,qp,rt,kidnsa,nitmax
    IMPLICIT NONE
    REAL(rkind):: wp,taue
    INTEGER(ikind):: nsa

    wp=0.d0
    taue=0.d0

    WRITE(6,601) t,wp,taue,qp(0)
601 FORMAT(' # T: ',F7.3,'(s)     WP:',F7.2,'(MJ)  ', &
           '  TAUE:',F7.3,'(s)   Q0:',F7.3)
    WRITE(6,602) (kidnsa(nsa),rt(nsa,0),nsa=1,nsamax)
602 FORMAT(4(' T',A1,':',F7.3,'(keV)  ':))

    WRITE(6,'(A13,I4)') 'Iterations: ',nitmax
    nitmax = 0
    
    RETURN
  END SUBROUTINE tr_status

! ***** save data for time history *****

  SUBROUTINE tr_save_ngt

    USE trcomm, ONLY : &
         nrmax,nsamax,ngtmax,neqmax,               &
         ngt,gvt,gvts,gvti,gvrt,gvrts,gvrtj,gparts,     &
         t,rn,ru,rt,qp,jtot,joh,htr,rip
    USE trcoeftb, ONLY: Pereverzev_check
    IMPLICIT NONE
    INTEGER(ikind):: nsa,nr
    REAL(rkind),DIMENSION(neqmax,0:nrmax):: add_prv

    IF(ngt >= ngtmax) RETURN

    ngt=ngt+1
    gvt(ngt, 0) = t

    ! ----- values on the axis or the edge -----
    gvt(ngt, 1) = qp(0)
    gvt(ngt, 2) = qp(nrmax)

    DO nsa=1,nsamax
       gvts(ngt,nsa, 1) = rn(nsa,0)
       gvts(ngt,nsa, 2) = ru(nsa,0)
       gvts(ngt,nsa, 3) = rt(nsa,0)
    END DO

    ! plasma curent
    gvti(ngt,1) = rip
!    gvti(ngt,2) = 
!    gvti(ngt,3) = 

    ! ----- radial profile -----
    DO nr=0,nrmax
       gvrt(nr,ngt, 1) = qp(nr)
    END DO

    DO nsa=1,nsamax
          gvrts(0:nrmax,ngt,nsa, 1) = rn(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 2) = ru(nsa,0:nrmax)
          gvrts(0:nrmax,ngt,nsa, 3) = rt(nsa,0:nrmax)
    END DO

    gvrtj(0:nrmax,ngt,1) = jtot(0:nrmax) + htr(1,0:nrmax)
    gvrtj(0:nrmax,ngt,2) = joh(0:nrmax)
    gvrtj(0:nrmax,ngt,3) = htr(1,0:nrmax)


    ! for Pereverzev method
    ! numerically addtional term in nodal equation (relative value)
    IF(t == 0)THEN
       gparts(0:nrmax,ngt,1:nsamax,1) = 0.d0
       gparts(0:nrmax,ngt,1:nsamax,2) = 0.d0
       gparts(0:nrmax,ngt,1:nsamax,3) = 0.d0
    ELSE
       CALL Pereverzev_check(add_prv)
       DO nsa=1,nsamax
          DO nr=0,nrmax
             gparts(nr,ngt,nsa,1) = add_prv(1+3*nsa-2,nr)
             gparts(nr,ngt,nsa,2) = add_prv(1+3*nsa-1,nr)
             gparts(nr,ngt,nsa,3) = add_prv(1+3*nsa  ,nr)
          END DO
       END DO
    END IF

    RETURN
  END SUBROUTINE tr_save_ngt

END MODULE trresult
