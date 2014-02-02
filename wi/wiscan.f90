!    $Id$

MODULE wiscan
  PRIVATE
  PUBLIC wi_scan
 
CONTAINS

  SUBROUTINE wi_scan(ierr)

    USE wicomm,ONLY: rkind,ikind,modewi,ntaumax,taumin,taumax,alfa,beta,aky, &
         xmax,pn0,nxmax,nwmax 
    USE wiunmag,ONLY: wi_unmag
    USE libgrf,ONLY: grd1d
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: ntau
    REAL(rkind):: tau,dtau,rk0l,ratea,aky_save
    REAL(rkind),DIMENSION(ntaumax):: taua,rateaa

    aky_save=aky
    rk0l=beta/alfa
    dtau=(taumax-taumin)/(ntaumax-1)
    write(6,'(A,1P4E12.4)') '## alfa,beta,rk0l,rk0l**(1/3)=', &
                             alfa,beta,rk0l,rk0l**(1.D0/3.D0)
    write(6,'(A,1P2E12.4,2I12)') '## xmax,pn0,nxmax,nwmax=', &
                             xmax,pn0,nxmax,nwmax
    DO ntau=1,ntaumax
       tau=taumin+dtau*(ntau-1)
       aky=tau/rk0l**(1.D0/3.D0)
       IF(aky < 1.0) THEN
          CALL wi_unmag(ratea,0,ierr)
       ELSE
          ratea=0.D0
       END IF
       WRITE(6,'(A,I5,1P3E12.4)') 'ntau,tau,aky,ratea=',ntau,tau,aky,ratea
       taua(ntau)=tau
       rateaa(ntau)=ratea
    END DO
    aky=aky_save

    CALL PAGES
    CALL GRD1D(0,taua,rateaa,ntaumax,ntaumax,1,TITLE='@abs vs tau@')
    CALL PAGEE

    RETURN
  END SUBROUTINE wi_scan
END MODULE wiscan
