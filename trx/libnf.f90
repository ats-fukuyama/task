! libnf.f90

MODULE libnf_local
  USE bpsd_kinds
  INTEGER:: id_nf_local,pm_local
  REAL(rkind):: temperature_local
END MODULE libnf_local

MODULE libnf
  USE bpsd_kinds
  USE bpsd_constants

  ! Fusion model
  !   model_pnf=0 : no fusion reaction
  !   model_pnf=1 : D + T -> He4 + n              nnfmax=1  nsmax=4 DT
  !   model_pnf=2 : D + D -> T + p                nnfmax=4  nsmax=6 DD1 DD2
  !                 D + D -> He3 + n                                DD3
  !                 D + T -> He4 + n                                DT
  !   model_pnf=3 : D + D -> T + p                nnfmax=6  nsmax=6 DD1 DD2
  !                 D + D -> He3 + n                                DD3
  !                 D + T -> He4 + n                                DT
  !                 D + He3 -> He4 + p                              DHe31 DHe32
  !   model_pnf=4 : D + D -> T + p                nnfmax=13 nsmax=7 DD1 DD2
  !                 D + D -> He3 + n                                DD3
  !                 D + T -> He4 + n                                DT
  !                 D + He3 -> He4 + p                              DH31 DHe32
  !                 T + T -> He4 + 2n                               TT
  !                 T + He3 -> He4 + p + n                          THe31 THe32
  !                 T + He3 -> He4 + D                              THe33 THe34
  !                 T + He3 -> He5 + p                              THe35 THe36
  
  !   model_pnf=12: D + D -> T + p   | T + He4/2  nnfmax=4  nsmax=4 DD1 DD2
  !                 D + D -> He3 + n ! He4 + n                      DD3
  !                 D + T -> He4 + n                                DT
  !   model_pnf=14: D + D -> T + p   ! T +He4/2   nnfmax=13 nsmax=4 DD1 DD2
  !                 D + D -> He3 + n ! He4 + n                      DD3
  !                 D + T -> He4 + n                                DT
  !                 D + He3 -> He4 + p                              DHe31 DHe32
  !                 T + T -> He4 + 2n                               TT
  !                 T + He3 -> He4 + p + n                          THe31 THe32
  !                 T + He3 -> He4 + D                              THe33 THe34
  !                 T + He3 -> He5 + p ! He4 + p                    THe35 THe36
  ! Fusion reaction id
  
  INTEGER,PARAMETER,PUBLIC:: id_nf_DT=    1 ! D + T   -> <He4> +  n
  INTEGER,PARAMETER,PUBLIC:: id_nf_DD1=   2 ! D + D   -> <T>   +  p
  INTEGER,PARAMETER,PUBLIC:: id_nf_DD2=   3 ! D + D   ->  T    + <p> 
  INTEGER,PARAMETER,PUBLIC:: id_nf_DD3=   4 ! D + D   -> <He3> +  n
  INTEGER,PARAMETER,PUBLIC:: id_nf_DHe31= 5 ! D + He3 -> <He4> +  p
  INTEGER,PARAMETER,PUBLIC:: id_nf_DHe32= 6 ! D + He3 ->  He4  + <p>
  INTEGER,PARAMETER,PUBLIC:: id_nf_TT=    7 ! T + T   -> <He4> + 2n
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe31= 8 ! T + He3 -> <He4> +  p  + n
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe32= 9 ! T + He3 ->  He4  + <p> + n
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe33=10 ! T + He3 -> <He4> +  D
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe34=11 ! T + He3 ->  He4  + <D>
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe35=12 ! T + He3 -> <He5> +  p
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe36=13 ! T + He3 ->  He5  + <p>

  INTEGER,DIMENSION(:),ALLOCATABLE:: id_nf_nnf
  
  Integer,DIMENSION(13),PUBLIC:: &
       ns1_idnf,ns2_idnf,nsp_idnf
  REAL(rkind),DIMENSION(13),PUBLIC::  &
       wgt_idnf,eng_idnf,enn_idnf

       ! Duane coef (NEL Formulary 2019)

  REAL(rkind),DIMENSION(5,6):: &
       Duane=reshape((/46.097D0, 372.D0, 4.36D-4,  1.220D0, 0.D0, &
                       47.88D0,  482.D0, 3.08D-4,  1.177D0, 0.D0, &
                       45.95D0,  5.02D4, 1.368D-2, 1.076D0, 409.D0, &
                       89.27D0,  2.59D4, 3.98D-3,  1.297D0, 647.D0, &
                       38.39D0,  448.D0, 1.02D-3,  2.09D0,  0.D0, &
                       123.1D0, 1.125D4, 0.D0,     0.D0,    0.D0/), &
                       (/5,6/))

  ! temperature range for reaction rate sigmav
    
  REAL(rkind),DIMENSION(10):: &
       tempa=(/ 1.D0, 2.D0, 5.D0, 10.D0, 20.D0, &
       50.D0, 100.D0, 200.D0, 500.D0, 1000.D0 /)
  
  ! reaction rate sigmav for DD
  REAL(rkind),DIMENSION(10):: &
       svnf_dd=(/ 1.5D-21, 5.4D-21, 1.8D-19, 1.2D-18, 5.2D-18, &
                  2.1D-17, 4.5D-17, 8.8D-17, 1.8D-16, 2.2D-16 /)
  ! reaction rate sigmav for DT
  REAL(rkind),DIMENSION(10):: &
       svnf_dt=(/ 5.5D-21, 2.6D-19, 1.3D-17, 1.1D-16, 4.2D-16, &
                  8.7D-16, 8.5D-16, 6.3D-16, 3.7D-16, 2.7D-16 /)
  ! reaction rate sigmav for DHe3
  REAL(rkind),DIMENSION(10):: &
       svnf_dhe3=(/ 1.0D-26, 1.4D-23, 6.7D-21, 2.3D-19, 3.8D-18, &
                    5.4D-17, 1.6D-16, 2.4D-16, 2.3D-16, 1.8D-16 /)
  ! reaction rate sigmav for TT
  REAL(rkind),DIMENSION(10):: &
       svnf_tt=(/ 3.3D-22, 7.1D-21, 1.4D-19, 7.2D-19, 2.5D-18, &
                  8.7D-18, 1.9D-17, 4.2D-17, 8.4D-17, 8.0D-17 /)
  ! reaction rate sigmav for THe3
  REAL(rkind),DIMENSION(10):: &
       svnf_the3=(/ 1.0D-28, 1.0D-25, 2.1D-22, 1.2D-20, 2.6D-19, &
                    5.3D-18, 2.7D-17, 9.2D-17, 2.9D-16, 5.2D-16 /)

  ! Mass of incident particle

  REAL(rkind),DIMENSION(4,10):: &
       usvnf_dd,usvnf_dt,usvnf_dhe3,usvnf_tt,usvnf_the3
  REAL(rkind),DIMENSION(10):: &
       dsvnf,tempa_log

  ! *** library subroutines ***

  PUBLIC set_usigmav_nf    ! set_usigmav_nf
  PUBLIC sigma_nf          ! sigma_nf(id_nf,energy)
  PUBLIC sigmav_nf         ! sigmav_nf(id_nf,temperature)

CONTAINS
  
  ! --- set spline coefficients for reaction rate sigmav

  SUBROUTINE set_usigmav_nf
    USE trcomm
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),DIMENSION(10):: dsvnf
    INTEGER:: ntemp,id,ierr

    SELECT CASE(model_pnf)
    CASE(0)
       nnfmax=0
       RETURN
    CASE(1)
       nnfmax=1
    CASE(2,12)
       nnfmax=4
    CASE(3)
       nnfmax=6
    CASE(4,14)
       nnfmax=13
    CASE DEFAULT
       WRITE(6,*) 'XX Error libnf: undefined model_pnf: model_pnf=',model_pnf
       STOP
    END SELECT

    IF(ALLOCATED(id_nf_nnf)) DEALLOCATE(id_nf_nnf)
    ALLOCATE(id_nf_nnf(nnfmax))

    SELECT CASE(model_pnf)
    CASE(1)
       id_nf_nnf(1)=id_nf_dt
    CASE(2,12)
       id_nf_nnf(1)=id_nf_dt
       id_nf_nnf(2)=id_nf_dd1
       id_nf_nnf(3)=id_nf_dd2
       id_nf_nnf(4)=id_nf_dd3
    CASE(3)
       id_nf_nnf(1)=id_nf_dt
       id_nf_nnf(2)=id_nf_dd1
       id_nf_nnf(3)=id_nf_dd2
       id_nf_nnf(4)=id_nf_dd3
       id_nf_nnf(5)=id_nf_dhe31
       id_nf_nnf(6)=id_nf_dhe32
    CASE(4,14)
       id_nf_nnf(1)=id_nf_dt
       id_nf_nnf(2)=id_nf_dd1
       id_nf_nnf(3)=id_nf_dd2
       id_nf_nnf(4)=id_nf_dd3
       id_nf_nnf(5)=id_nf_dhe31
       id_nf_nnf(6)=id_nf_dhe32
       id_nf_nnf(7)=id_nf_tt
       id_nf_nnf(8)=id_nf_the31
       id_nf_nnf(9)=id_nf_the32
       id_nf_nnf(10)=id_nf_the33
       id_nf_nnf(11)=id_nf_the34
       id_nf_nnf(12)=id_nf_the35
       id_nf_nnf(13)=id_nf_the36
    CASE DEFAULT
       WRITE(6,*) 'XX Error libnf: undefined model_pnf: model_pnf=',model_pnf
       STOP
    END SELECT

    SELECT CASE(model_pnf)
    CASE(1,12,13)
       IF(NS_D*NS_T*NS_He4.EQ.0) THEN
          IF(NS_D.EQ.0) WRITE(6,*)   'XX Error: libnf: NS_D=0'
          IF(NS_T.EQ.0) WRITE(6,*)   'XX Error: libnf: NS_T=0'
          IF(NS_He4.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He4=0'
          STOP
       END IF
    CASE(2,3)
       IF(NS_D*NS_T*NS_He4*NS_H*NS_He3.EQ.0) THEN
          IF(NS_D.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_D=0'
          IF(NS_T.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_T=0'
          IF(NS_H.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_H=0'
          IF(NS_He4.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He4=0'
          IF(NS_He3.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He3=0'
          STOP
       END IF
    CASE(4)
       IF(NS_D*NS_T*NS_He4*NS_H*NS_He3*NS_He5.EQ.0) THEN
          IF(NS_D.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_D=0'
          IF(NS_T.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_T=0'
          IF(NS_H.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_H=0'
          IF(NS_He4.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He4=0'
          IF(NS_He3.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He3=0'
          IF(NS_He5.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He5=0'
          STOP
       END IF
    CASE DEFAULT
       WRITE(6,*) 'XX Error libnb: undefined model_pnf: model_pnf=',model_pnf
       STOP
    END SELECT

    
    
    ns1_idnf(id_nf_dt)=NS_D
    ns2_idnf(id_nf_dt)=NS_T
    wgt_idnf(id_nf_dt)=1.0D0
    nsp_idnf(id_nf_dt)=NS_He4
    eng_idnf(id_nf_dt)=3.5D3*RKEV
    enn_idnf(id_nf_dt)=14.1D3*RKEV

    ns1_idnf(id_nf_dd1)=NS_D
    ns2_idnf(id_nf_dd1)=NS_D
    wgt_idnf(id_nf_dd1)=0.5D0
    nsp_idnf(id_nf_dd1)=NS_T
    eng_idnf(id_nf_dd1)=1.01D3*RKEV
    enn_idnf(id_nf_dd1)=0.D0
    
    ns1_idnf(id_nf_dd2)=NS_D
    ns2_idnf(id_nf_dd2)=NS_D
    wgt_idnf(id_nf_dd2)=0.5D0
    IF(model_pnf.EQ.12) THEN
       nsp_idnf(id_nf_dd2)=NS_He4
    ELSE
       nsp_idnf(id_nf_dd2)=NS_H
    END IF
    eng_idnf(id_nf_dd2)=3.02D3*RKEV
    enn_idnf(id_nf_dd3)=0.D0

    ns1_idnf(id_nf_dd3)=NS_D
    ns2_idnf(id_nf_dd3)=NS_D
    wgt_idnf(id_nf_dd3)=0.5D0
    IF(model_pnf.EQ.12) THEN
       nsp_idnf(id_nf_dd3)=NS_He3
    ELSE
       nsp_idnf(id_nf_dd3)=NS_He4
    END IF
    eng_idnf(id_nf_dd3)=0.82D3*RKEV
    enn_idnf(id_nf_dd3)=2.45D3*RKEV

    IF(model_pnf.GE.3) THEN
       ns1_idnf(id_nf_dhe31)=NS_D
       ns2_idnf(id_nf_dhe31)=NS_He3
       wgt_idnf(id_nf_dhe31)=1.D0
       nsp_idnf(id_nf_dhe31)=NS_He4
       eng_idnf(id_nf_dhe31)=3.6D3*RKEV
       enn_idnf(id_nf_dhe31)=0.D0

       ns1_idnf(id_nf_dhe32)=NS_D
       ns2_idnf(id_nf_dhe32)=NS_He3
       wgt_idnf(id_nf_dhe32)=1.D0
       nsp_idnf(id_nf_dhe32)=NS_H
       eng_idnf(id_nf_dhe32)=14.7D3*RKEV
       enn_idnf(id_nf_dhe32)=0.D0

       ns1_idnf(id_nf_tt)=NS_T
       ns2_idnf(id_nf_tt)=NS_T
       wgt_idnf(id_nf_tt)=1.D0
       nsp_idnf(id_nf_tt)=NS_He4
       eng_idnf(id_nf_tt)=1.25D3*RKEV  ! 11.3MeV*0.25/2.25
       enn_idnf(id_nf_tt)=10.05D3*RKEV ! 11.3Mev*2.00/2.25

       ns1_idnf(id_nf_the31)=NS_T
       ns2_idnf(id_nf_the31)=NS_He3
       wgt_idnf(id_nf_the31)=0.51D0
       nsp_idnf(id_nf_the31)=NS_He4
       eng_idnf(id_nf_the31)=1.34D3*RKEV ! 12.1MeV*0.25/2.25
       enn_idnf(id_nf_the31)=5.38D3*RKEV ! 12.1Mev*1.00/2.25

       ns1_idnf(id_nf_the32)=NS_T
       ns2_idnf(id_nf_the32)=NS_He3
       wgt_idnf(id_nf_the32)=0.51D0
       IF(model_pnf.EQ.14) THEN
          nsp_idnf(id_nf_the32)=NS_He4
       ELSE
          nsp_idnf(id_nf_the32)=NS_H
       END IF
       eng_idnf(id_nf_the32)=5.38D3*RKEV ! 12.1MeV*1.0/2.25
       enn_idnf(id_nf_the32)=0.D0

       ns1_idnf(id_nf_the33)=NS_T
       ns2_idnf(id_nf_the33)=NS_He3
       wgt_idnf(id_nf_the33)=0.43D0
       nsp_idnf(id_nf_the33)=NS_He4
       eng_idnf(id_nf_the33)=4.8D3*RKEV
       enn_idnf(id_nf_the33)=0.D0

       ns1_idnf(id_nf_the34)=NS_T
       ns2_idnf(id_nf_the34)=NS_He3
       wgt_idnf(id_nf_the34)=0.43D0
       nsp_idnf(id_nf_the34)=NS_D
       eng_idnf(id_nf_the34)=9.58D3*RKEV
       enn_idnf(id_nf_the34)=0.D0

       ns1_idnf(id_nf_the35)=NS_T
       ns2_idnf(id_nf_the35)=NS_He3
       wgt_idnf(id_nf_the35)=0.06D0
       nsp_idnf(id_nf_the35)=NS_He5
       eng_idnf(id_nf_the35)=1.89D3*RKEV
       enn_idnf(id_nf_the35)=0.D0

       ns1_idnf(id_nf_the36)=NS_T
       ns2_idnf(id_nf_the36)=NS_He3
       wgt_idnf(id_nf_the36)=0.06D0
       nsp_idnf(id_nf_the36)=NS_H
       eng_idnf(id_nf_the36)=9.46D3*RKEV
       enn_idnf(id_nf_the36)=0.D0
    END IF

    DO ntemp=1,10
       tempa_log(ntemp)=LOG10(tempa(ntemp))
    END DO
    
    CALL SPL1D(tempa_log,svnf_dt,  dsvnf,usvnf_dt,  10,0,ierr)
    id=1
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_dd,  dsvnf,usvnf_dd,  10,0,ierr)
       id=2
    ENDIF
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_dhe3,dsvnf,usvnf_dhe3,10,0,ierr)
       id=3
    END IF
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_tt,  dsvnf,usvnf_tt,  10,0,ierr)
       id=4
    END IF
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_the3,dsvnf,usvnf_the3,10,0,ierr)
       id=5
    END IF
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX SPL1D error in set_usvnf: id=',id
       STOP
    END IF
    RETURN
  END SUBROUTINE set_usigmav_nf

  ! --- cross section of nuclear fusion reaction ---
  ! ---     in barn (10^{-28}m^{-2})
  ! ---     as a function of energy in keV

  FUNCTION sigma_nf(id_nf,energy)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: id_nf
    REAL(rkind),INTENT(IN):: energy
    REAL(rkind):: sigma_nf
    
    IF(id_nf.LT.1.OR.id_nf.GT.6) THEN
       WRITE(6,'(A,I4)') &
            'XX sigma_nf: input error: undefined id_nf: ',id_nf
       STOP
    ENDIF
    
    IF(energy.LE.0.D0) THEN
       WRITE(6,*) 'XX sigma_duane: input error: non-positive energy: ',energy
       STOP
    ENDIF
    
    sigma_nf=(Duane(5,id_nf) &
         +Duane(2,id_nf) &
         /((Duane(4,id_nf)-Duane(3,id_nf)*energy)**2+1.D0)) &
         /(energy*(EXP(Duane(1,id_nf)/SQRT(energy))-1.D0))
    RETURN
  END FUNCTION sigma_nf
  
  ! --- reaction rate of nuclear fusion: sigmav  ---
  ! ---     as a function of temperature in keV

  FUNCTION sigmav_nf(id_nf,temperature)

    USE libspl1d
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id_nf
    REAL(rkind),INTENT(IN):: temperature
    REAL(rkind):: sigmav_nf,temperature_log
    INTEGER:: ierr
    
    IF(id_nf.LT.1.OR.id_nf.GT.13) THEN
       WRITE(6,'(A,I4)') 'XX sigmav_nf: input error: undefined id_nf: ',id_nf
       STOP
    END IF

    IF(temperature.LT.1.D0) THEN
       sigmav_nf=0.D0
       RETURN
    END IF

    IF(temperature.GT.1000.D0) THEN
       WRITE(6,'(A,ES12.4)') &
            'XX sigmav_nf: input error: Too high temperature: temperature=', &
            temperature
       STOP
    END IF

    temperature_log=LOG10(temperature)
    SELECT CASE(id_nf)
    CASE(id_nf_dt)
       CALL SPL1DF(temperature_log,sigmav_nf,tempa_log,usvnf_dt,10,ierr)
    CASE(id_nf_dd1,id_nf_dd2,id_nf_dd3)
       CALL SPL1DF(temperature_log,sigmav_nf,tempa_log,usvnf_dd,10,ierr)
    CASE(id_nf_dhe31,id_nf_dhe32)
       CALL SPL1DF(temperature_log,sigmav_nf,tempa_log,usvnf_dhe3,10,ierr)
    CASE(id_nf_tt)
       CALL SPL1DF(temperature_log,sigmav_nf,tempa_log,usvnf_tt,10,ierr)
    CASE(id_nf_the31,id_nf_the32,id_nf_the33, &
         id_nf_the34,id_nf_the35,id_nf_the36)
       CALL SPL1DF(temperature_log,sigmav_nf,tempa_log,usvnf_the3,10,ierr)
    END SELECT
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)')     'XX SPL1DF error in sigmav_nf: id_nf=',id_nf
       WRITE(6,'(A,ES12.4)') '       temperature=',temperature
       STOP
    END IF
  END FUNCTION sigmav_nf

  ! --- sigmav_nf by integeral over energy ---
  !         <sigmav>=\int_0^\infty 4\pi v^2 dv sigma v f(v)
  !                  f(v)=(m/2\pi T)^{3/2} exp(-mv^2/2T)   normalized to 1
  !                  E=mv^2/2
  !                  v=SQRT{2E/m}
  !                  dv=(1/2) SQRT{2/mE} dE
  !                  4\pi v^2 dv=4\pi (2E/m) (1/2) SQRT{2/mE} dE
  !                             =4\pi SQRT{E^2/m^2 2/mE} dE
  !                             =4\pi SQRT{2E/m^3} dE
  !                  f(E)=(m/2\pi T)^{3/2} exp(-E/T)
  !                <f(E)>=\int_0^\infty 4\pi v^2 dv f(v) 
  !                      =\int_0^\infty 4\pi SQRT{2E/m^3} dE
  !                               (m/2\pi T)^{3/2} exp(-E/T)
  !                      =\int_0^\infty SQRT{16 \pi^2 2E/m^3 m^3/8 \pi^3 T^3}
  !                               dE exp(-E/T)
  !                      =\int_0^\infty SQRT{4E/\pi T} exp(-E/T) dE/T
  !                      =\int_0^\infty SQRT{4X/\pi} exp(-X) dX
  !                      = (2/SQRT{\pi}) \int_0^\infty SQRT{X} exp(-X) dX
  !                      = 1
  !          <sigmav>=\int_0^\infty 4\pi SQRT{2E/m^3} dE sigma SQRT{2E/m} f(E)
  !                  =\int_0^\infty 4\pi 2E/m^2 dE sigma
  !                                   (m/2\pi T)^{3/2} exp(-E/T)
  !                  =\int_0^\infty SQRT{8 E^2/m\pi T} sigma exp(-E/T) dE/T
  !                  =\int_0^\infty SQRT{8T/m\pi} (E/T) sigma exp(-E/T) dE/T
  !                  =\int_0^\infty SQRT{8T/m\pi} X sigma exp(-X) dX

  ! --- sigmav for energy --- X=energy/temperature

  FUNCTION sigmav_nf_local(X)
    USE trcomm,ONLY: RKEV
    USE libnf_local
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: sigmav_nf_local
    REAL(rkind):: energy,velocity

    energy=temperature_local*X
    velocity=SQRT(2.D0*energy*RKEV/pm_local)
    sigmav_nf_local=sigma_nf(id_nf_local,energy)*velocity
    RETURN
  END FUNCTION sigmav_nf_local


  FUNCTION sigmav_nf_int(id_nf,temperature)

    USE plcomm
    USE libnf_local
    USE libde
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id_nf
    REAL(rkind),INTENT(IN):: temperature
    REAL(rkind):: sigmav_nf_int
    REAL(rkind):: error_int,H0,EPS
    INTEGER:: ILST

    IF(id_nf.LT.1.OR.id_nf.GT.6) THEN
       WRITE(6,'(A,I4)') 'XX sigmav_nf: input error: undefined id_nf: ',id_nf
       STOP
    END IF

    IF(temperature.LT.1.D0) THEN
       sigmav_nf_int=0.D0
       RETURN
    END IF

    IF(temperature.GT.1000.D0) THEN
       WRITE(6,'(A,ES12.4)') &
            'XX sigmav_nf: input error: Too high temperature: temperature=', &
            temperature
       STOP
    END IF

    id_nf_local=id_nf
    temperature_local=temperature
    pm_local=PA(nsp_idnf(id_nf))*AMP

    H0=1.D-4
    EPS=1.D-6
    
    CALL DEHIFE(sigmav_nf_int,error_int,H0,eps,ilst,sigmav_nf_local, &
         'sigmav_nf_int')

    RETURN
  END FUNCTION sigmav_nf_int
  
END MODULE libnf
      
