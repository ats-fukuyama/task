MODULE wiscan
  PRIVATE
  PUBLIC wi_scan
 
CONTAINS

  SUBROUTINE wi_scan(ierr)

    USE wicomm,ONLY: ikind
    USE wiparm,ONLY: wi_parm
    INTEGER(ikind),INTENT(OUT):: ierr
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind):: mode

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## WI SCAN MENU: T/TAU  A/ALFA  X/EXIT'

    CALL TASK_KLIN(line,kid,mode,wi_parm)
    IF(mode /= 1) GOTO 1

    IF(kid.EQ.'T') THEN
       CALL wi_scan_tau(ierr)
    ELSEIF(kid.EQ.'A') THEN
       CALL wi_scan_alfa(ierr)
    ELSEIF(kid.EQ.'X') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX WI SCAN MENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wi_scan


  SUBROUTINE wi_scan_tau(ierr)

    USE wicomm,ONLY: rkind,ikind,ntaumax,taumin,taumax,alfa,beta,any, &
         xmax,pn0,nxmax,nwmax,kfscan
    USE wiexec,ONLY: wi_exec
    USE libgrf,ONLY: grd1d
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: ntau
    INTEGER(ikind),PARAMETER:: nfl=21
    REAL(rkind):: tau,dtau,rk0l,ratea,any_save
    REAL(rkind),DIMENSION(ntaumax):: taua,rateaa

    any_save=any
    rk0l=beta/alfa
    dtau=(taumax-taumin)/(ntaumax-1)
    
    IF(TRIM(kfscan)//'X'.NE.'X') CALL FWOPEN(nfl,kfscan,1,1,'SCAN',ierr)

    DO ntau=1,ntaumax
       tau=taumin+dtau*(ntau-1)
       any=tau/rk0l**(1.D0/3.D0)
       IF(any < 1.0) THEN
          CALL wi_exec(0,ratea,ierr)
       ELSE
          ratea=0.D0
       END IF
       WRITE(6,'(A,I5,1P3E12.4)') 'ntau,tau,any,ratea=',ntau,tau,any,ratea
       IF(TRIM(kfscan)//'X'.NE.'X') &
            WRITE(nfl,'(I5,1P3E12.4)') ntau,tau,any,ratea
       taua(ntau)=tau
       rateaa(ntau)=ratea
    END DO
    any=any_save

    IF(TRIM(kfscan)//'X'.NE.'X') CLOSE(nfl)

    CALL PAGES
    CALL GRD1D(0,taua,rateaa,ntaumax,ntaumax,1,TITLE='@abs vs tau@')
    CALL PAGEE

    RETURN
  END SUBROUTINE wi_scan_tau

  SUBROUTINE wi_scan_alfa(ierr)

    USE wicomm,ONLY: rkind,ikind,nalfamax,alfamin,alfamax,alfa,beta,any, &
         xmax,pn0,nxmax,nwmax,kfscan
    USE wiexec,ONLY: wi_exec
    USE libgrf,ONLY: grd1d
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: nalfa
    INTEGER(ikind),PARAMETER:: nfl=21
    REAL(rkind):: dalfa,rk0l,ratea,alfa_save
    REAL(rkind),DIMENSION(nalfamax):: rk0la,rateaa

    alfa_save=alfa
    dalfa=(log(alfamax)-log(alfamin))/(nalfamax-1)
    
    
    IF(TRIM(kfscan)//'X'.NE.'X') CALL FWOPEN(nfl,kfscan,1,1,'SCAN',ierr)

    DO nalfa=1,nalfamax
       alfa=exp(log(alfamin)+dalfa*(nalfa-1))
       rk0l=1.D0/alfa
       alfa=alfa*beta
       IF(any < 1.0) THEN
          CALL wi_exec(0,ratea,ierr)
       ELSE
          ratea=0.D0
       END IF
       WRITE(6,'(A,I5,1P3E12.4)') 'nalfa,alfa,rk0l,ratea=', &
                                   nalfa,alfa,rk0l,ratea
       IF(TRIM(kfscan)//'X'.NE.'X') &
            WRITE(nfl,'(I5,1P3E12.4)') nalfa,alfa,rk0l,rateaa
       rk0la(nalfa)=LOG10(rk0l)
       rateaa(nalfa)=ratea
    END DO
    alfa=alfa_save

    IF(TRIM(kfscan)//'X'.NE.'X') CLOSE(nfl)

    CALL PAGES
    CALL GRD1D(0,rk0la,rateaa,nalfamax,nalfamax,1,TITLE='@abs vs k0L@',&
               MODE_LS=1)
    CALL PAGEE

    RETURN
  END SUBROUTINE wi_scan_alfa
END MODULE wiscan
