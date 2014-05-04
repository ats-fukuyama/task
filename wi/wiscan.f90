!    $Id$

MODULE wiscan
  PRIVATE
  PUBLIC wi_scan
 
CONTAINS

  SUBROUTINE wi_scan(ierr)

    USE wicomm,ONLY: rkind,ikind,ntaumax,taumin,taumax,alfa,beta,any, &
         xmax,pn0,nxmax,nwmax 
    USE wiexec,ONLY: wi_exec
    USE libgrf,ONLY: grd1d
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: ntau
    REAL(rkind):: tau,dtau,rk0l,ratea,any_save
    REAL(rkind),DIMENSION(ntaumax):: taua,rateaa

    any_save=any
    rk0l=beta/alfa
    dtau=(taumax-taumin)/(ntaumax-1)

    DO ntau=1,ntaumax
       tau=taumin+dtau*(ntau-1)
       any=tau/rk0l**(1.D0/3.D0)
       IF(any < 1.0) THEN
          CALL wi_exec(0,ratea,ierr)
       ELSE
          ratea=0.D0
       END IF
       WRITE(6,'(A,I5,1P3E12.4)') 'ntau,tau,any,ratea=',ntau,tau,any,ratea
       taua(ntau)=tau
       rateaa(ntau)=ratea
    END DO
    any=any_save

    CALL PAGES
    CALL GRD1D(0,taua,rateaa,ntaumax,ntaumax,1,TITLE='@abs vs tau@')
    CALL PAGEE

    RETURN
  END SUBROUTINE wi_scan
END MODULE wiscan
