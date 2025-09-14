! libpdkf_

MODULE libpdkf
  USE task_kinds,ONLY: rkind
  USE task_constants,ONLY: PI

  INTEGER:: N1
  REAL(rkind):: G1,G2,G3,G4,G5
  INTEGER:: np_tau
  REAL(rkind):: xi,eta,rnu

  PRIVATE
  PUBLIC pdkf_hift
  PUBLIC pdkf_ft

CONTAINS

  !   *****  REAL PART  *****

  FUNCTION func_pdkf_r(tau)
    
    !   real part fucntion of tau

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau
    REAL(rkind):: func_pdkf_r
    REAL(rkind):: arg

    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_r=tau**np_tau*DEXP(arg)*DCOS(rnu*tau)
    ELSE
       func_pdkf_r=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_r

  !   *****  IMAGINARY PART  *****

  FUNCTION func_pdkf_i(tau)
    
    !   imaginary part fucntion of tau

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau
    REAL(rkind):: func_pdkf_i
    REAL(rkind):: arg

    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_i=tau**np_tau*DEXP(arg)*DSIN(rnu*tau)
    ELSE
       func_pdkf_i=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_i

  !   *****  REAL PART  *****

  FUNCTION func_pdkf_rab(tau,Btau,tauA)
    
    ! real part fucntion of tau, tau_max-tau, tau-tau_min

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau,Btau,tauA
    REAL(rkind):: func_pdkf_rab
    REAL(rkind):: arg,dummy

    dummy=Btau
    dummy=tauA
    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_rab=tau**np_tau*DEXP(arg)*DCOS(rnu*tau)
    ELSE
       func_pdkf_rab=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_rab

  !     *****  IMAGINARY PART  *****

  FUNCTION func_pdkf_iab(tau,Btau,tauA)
    
    ! imaginary part fucntion of tau, tau_max-tau, tau-tau_min

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: tau,Btau,tauA
    REAL(rkind):: func_pdkf_iab
    REAL(rkind):: arg,dummy

    dummy=Btau
    dummy=tauA
    arg=0.5D0*(xi/tau)**2+0.5D0*(eta*tau)**2
    IF(arg.LT.176.D0) THEN
       func_pdkf_iab=tau**np_tau*DEXP(arg)*DSIN(rnu*tau)
    ELSE
       func_pdkf_iab=0.D0
    END IF
    RETURN
  END FUNCTION func_pdkf_iab

  !   *** plasma dispersion kernel function (zero-infinity DE integral) ***

  FUNCTION pdkf_hift(xi_,eta_,rnu_,np_tau_)
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi_,eta_,rnu_
    INTEGER,INTENT(IN):: np_tau_
    COMPLEX(rkind):: pdkf_hift
    REAL(rkind):: H0,EPS,SR,SI,ESR,ESI
    INTEGER:: ILST,M

    xi=xi_
    eta=eta_
    rnu=rnu_
    np_tau=np_tau_
    
    H0=0.5D0
    EPS=1.D-8
    ILST=0

    CALL DEHIFT(SR,ESR,H0,EPS,ILST,func_pdkf_r,'func_pdkf_r')
    CALL DEHIFT(SI,ESI,H0,EPS,ILST,func_pdkf_i,'func_pdkf_i')
    pdkf_hift=DCMPLX(SR,SI)
    RETURN
  END FUNCTION pdkf_hift

  !   *** plasma dispersion kernel function (finite DE integral) ***

  FUNCTION pdkf_ft(xi_,eta_,rnu_,np_tau_)
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi_,eta_,rnu_
    INTEGER,INTENT(IN):: np_tau_
    COMPLEX(rkind):: pdkf_ft
    REAL(rkind):: H0,EPS,tau_min,tau_max,SR,SI,ESR,ESI
    INTEGER:: ILST,M

    xi=xi_
    eta=eta_
    rnu=rnu_
    np_tau=np_tau_
    
    H0=0.5D0
    EPS=1.D-8
    ILST=0

    tau_min=xi/18.D0       ! xi/tau<18
    IF(eta.EQ.0.D0) THEN
       tau_max=0.D0
    ELSE
       tau_max=18.D0/eta      ! eta*tau<18
    END IF

    IF(tau_max.EQ.0.D0) THEN
       CALL DEHIFT(SR,ESR,H0,EPS,ILST,func_pdkf_r,'func_pdkf_r')
       CALL DEHIFT(SI,ESI,H0,EPS,ILST,func_pdkf_i,'func_pdkf_i')
    ELSE
       CALL DEFTAB(tau_min,tau_max,SR,ESR,H0,EPS,ILST,func_pdkf_rab, &
            'func_pdkf_rab')
       CALL DEFTAB(tau_min,tau_max,SI,ESI,H0,EPS,ILST,func_pdkf_iab, &
            'func_pdkf_iab')
    END IF
    pdkf_ft=DCMPLX(SR,SI)
    RETURN
  END FUNCTION pdkf_ft
END MODULE libpdkf
