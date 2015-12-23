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
         xmax,xmin,dx0,pn0,nxmax,nwmax,kfscan,pi,xwint,xgrid
    USE wiexec,ONLY: wi_exec
    USE wiprep,ONLY: wi_prep
    USE wigout,ONLY: wi_gra1
    USE libgrf,ONLY: grd1d
    
    IMPLICIT NONE
    INTEGER(ikind),INTENT(OUT):: ierr
    INTEGER(ikind):: nalfa
    INTEGER(ikind),PARAMETER:: nfl=21
    REAL(rkind):: dalfa,rk0l,ratea,anb
    REAL(rkind):: alfa_save,xmax_save,xmin_save,dx0_save,xwint_save
    REAL(rkind),DIMENSION(nalfamax):: rk0la,rateaa

    alfa_save=alfa
    IF(nalfamax.EQ.1) THEN
       dalfa=0.D0
    ELSE
       dalfa=(log(alfamax)-log(alfamin))/(nalfamax-1)
    END IF
    xmax_save=xmax
    xmin_save=xmin
    xwint_save=xwint
    dx0_save=dx0
    
    IF(TRIM(kfscan)//'X'.NE.'X') CALL FWOPEN(nfl,kfscan,1,1,'SCAN',ierr)

    WRITE(6,'(A)') 'nalfa,alfa,rk0l,ratea,xmin,xmax,dx0='
    DO nalfa=1,nalfamax
       alfa=exp(log(alfamin)+dalfa*(nalfa-1))
       rk0l=1.D0/alfa
       IF(ALFA.LT.0.01D0) THEN
          dx0=0.3D0
          xmax=200.D0/beta
          xmin=-50.0D0/BETA
       ELSEIF(ALFA.LT.1.D0) THEN
          dx0=dx0_save
          xmax=10.D0/(alfa*beta)
          xmin=-10.0D0/BETA
       ELSEIF(ALFA*BETA.LT.1.D0) THEN
          dx0=dx0_save
          xmax=10.D0/(beta)
          xmin=-5.0D0/(beta)
       ELSE
!          dx0=0.5D0*dx0_save/log(alfa)
          dx0=0.1D0*dx0_save/log(alfa)
!          dx0=0.03D0*dx0_save/log(alfa)
!          IF(dx0.GT.0.3D0*dx0_save/alfa) dx0=0.3D0*dx0_save/alfa 
          xmax=500.D0*dx0
          xmin=-500.D0*dx0
!          IF(xmin.gt.-10.d0) xmin=-10.d0
!          IF(xmax.lt. 10.d0) xmax= 10.d0
!          xmax=10.D0
!          xmin=-10.D0
!          xwint=100*dx0
       END IF
       alfa=alfa*beta
       WRITE(6,'(I5,1P6E12.4)') nalfa,alfa,rk0l,0.D0,xmin,xmax,dx0
       CALL wi_prep
       ANB=DEXP(-ALFA*xgrid(nxmax))
       IF(any**2 < 1.0-ANB) THEN
          CALL wi_exec(0,ratea,ierr)
       ELSE
          ratea=0.D0
       END IF

       WRITE(6,'(I5,1P6E12.4)') nalfa,alfa,rk0l,ratea,xmin,xmax,dx0
       IF(TRIM(kfscan)//'X'.NE.'X') &
            WRITE(nfl,'(I5,1P3E12.4)') nalfa,alfa,rk0l,ratea
       rk0la(nalfa)=LOG10(rk0l)
       rateaa(nalfa)=ratea
       if(NALFAMAX.eq.1) CALL wi_gra1
    END DO
    alfa=alfa_save
    xmax=xmax_save
    xmin=xmin_save
    xwint=xwint_save
    dx0=dx0_save

    IF(TRIM(kfscan)//'X'.NE.'X') CLOSE(nfl)

    IF(NALFAMAX.NE.1) THEN
       CALL PAGES
       CALL GRD1D(0,rk0la,rateaa,nalfamax,nalfamax,1,TITLE='@abs vs k0L@',&
                  MODE_LS=1,FMIN=0.D0)
       CALL PAGEE
    END IF

    RETURN
  END SUBROUTINE wi_scan_alfa
END MODULE wiscan
