MODULE trresult

  PUBLIC tr_calc_global,tr_status,tr_save_ngt
  PRIVATE

CONTAINS

! ***** calculate global values *****

  SUBROUTINE tr_calc_global

    RETURN
  END SUBROUTINE tr_calc_global

! ***** simple status report *****

  SUBROUTINE tr_status

    USE trcomm, ONLY : rkind,ikind,nsamax,t,qp,rt,kidnsa
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
    RETURN
  END SUBROUTINE tr_status

! ***** save data for time history *****

  SUBROUTINE tr_save_ngt

    USE trcomm, ONLY : &
         rkind,ikind,nrmax,nsamax,ngtmax, &
         ngt,gvt,gvts,gvrt,gvrts,t,rg,rn,ru,rt,qp
    IMPLICIT NONE
    INTEGER(ikind):: nsa,nr

    IF(ngt >= ngtmax) RETURN

    ngt=ngt+1

    gvt(ngt, 0) = t
    gvt(ngt, 1) = qp(0)
    gvt(ngt, 2) = qp(nrmax)

    DO nsa=1,nsamax
       gvts(ngt,nsa, 1) = rn(nsa,0)
       gvts(ngt,nsa, 2) = ru(nsa,0)
       gvts(ngt,nsa, 3) = rt(nsa,0)
    END DO

    DO nr=0,nrmax
       gvrt(nr,ngt, 1) = qp(nr)
    END DO
    DO nsa=1,nsamax
       DO nr=0,nrmax
          gvrts(nr,ngt,nsa, 1) = rn(nsa,nr)
          gvrts(nr,ngt,nsa, 2) = ru(nsa,nr)
          gvrts(nr,ngt,nsa, 3) = rt(nsa,nr)
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_save_ngt

END MODULE trresult
