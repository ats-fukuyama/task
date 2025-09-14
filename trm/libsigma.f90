! libsigma.f90

MODULE libsigma
  USE bpsd_kinds
  USE bpsd_constants

  REAL(rkind),DIMENSION(10):: &
       tempa=(/ 1.D0, 2.D0, 5.D0, 10.D0, 20.D0, &
               50.D0, 100.D0, 200.D0, 500.D0, 1000.D0 /)
   
  ! for DT maxwellian
  REAL(rkind),DIMENSION(10):: &
       sigmavma_dt=(/ 5.5D-21, 2.6D-19, 1.3D-17, 1.1D-16, 4.2D-16, &
             8.7D-16, 8.5D-16, 6.3D-16, 3.7D-16, 2.7D-16 /)
  REAL(rkind),DIMENSION(10):: sigmavmal_dt

  ! for Dbulk T beam
  REAL(rkind),DIMENSION(10):: sigmaBdta, sigmaBdtal 
  REAL(rkind),DIMENSION(10):: dsigmaBdtal, usigmaBdtal
  

  REAL(rkind):: sigmaBDT 

  REAL(rkind),DIMENSION(10):: tempal
  REAL(rkind),DIMENSION(4,10):: usigmavmal_dt
  REAL(rkind),DIMENSION(4,10)::usigmaB_dt

  PRIVATE

  PUBLIC set_usigmavmal_dt
  PUBLIC sigmavm_dt
  PUBLIC set_usigmaBdt
  PUBLIC sigmaDbuTbm

CONTAINS

  ! set spline coefficients usvla for reaction rate for discrete energies

  SUBROUTINE set_usigmavmal_dt
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),DIMENSION(10):: dsigmavmal_dt
    INTEGER:: ierr,n

    DO n=1,10
       tempal(n)=LOG10(tempa(n))
       sigmavmal_dt(n)=LOG10(sigmavma_dt(n))
    END DO

    CALL SPL1D(tempal,sigmavmal_dt,dsigmavmal_dt,usigmavmal_dt,10,0,ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX spl1d sigmavmal_dt: ierr=',ierr
       STOP
    END IF

    RETURN
  END SUBROUTINE set_usigmavmal_dt

  ! interpolate sigma_v_maxwellian

  FUNCTION sigmavm_dt(temp)

    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: temp
    REAL(rkind):: templ
    REAL(rkind):: sigmavm_dt
    REAL(rkind):: sigmavml_dt
    INTEGER:: ierr

    templ=LOG10(temp)
    IF(templ.LE.tempal(1)) THEN
       sigmavm_dt=0.D0
    ELSE IF(templ.GE.tempal(10)) THEN
       sigmavm_dt=sigmavma_dt(10)
    ELSE
       CALL SPL1DF(templ,sigmavml_dt,tempal,usigmavmal_dt,10,ierr)
!       WRITE(6,'(5ES12.4)') temp,templ,sigmavml_dt,tempal(1),tempal(10)
       IF(ierr.NE.0) THEN
          WRITE(6,'(A,I4)') 'XX sigmavm_dt: SPL1DF error: ierr=',ierr
          STOP
       END IF
       sigmavm_dt=10.D0**sigmavml_dt
    END IF
    RETURN
  END FUNCTION sigmavm_dt
  

! !   calculation of spline for reactivity of DbulkTNBI
!   FUNCTION sigmaDbuTbm(temp,BENG)
!    Use trcomm !特定のものを使うのであればonlyを用いる
!    REAL(rkind),INTENT(IN):: temp, BENG
!    REAL(rkind):: ETNBI, A1, A2, A3, A4, A5, v0, v
!    INTEGER :: i, n, k
!    REAL(rkind),DIMENSION(10):: sigmaBdta,vcT
!    ! A1 = 45.95
!    ! A2 = 5.02E4
!    ! A3 = 1.368E-2
!    ! A4 = 1.076
!    ! A5 = 409

!    ! BENG = 1500
!    ETNBI = BENG !Beam energy in kev 
!    v0 = sqrt((2*ETNBI*1000*AEE)/AMT)!derive form ETNBI

!    DO n=1,10
!       vcT(n) = ([(3/4)*sqrt(pi)*(AME/AMT)]**(1/3))* sqrt(2*tempa(n)*AEE/AME) !function of tempa(n)
!    END DO

!    DO n=1,10
!       ! calculating integral of sigma*distribution function for certain electron temperatture
!       k = 10000
!       h = v0 / real(k)
!       sum = 0.0 !at v=0
!       do i = 1, k-1
!          v = i * h
!          ! sum = sum + sigBdt(v)*f(v)
!          sum = sum + sigBdt(v)*((v**3)/(v**3 + vcT(n)**3))
!       end do
!       sum = sum + 0.5 * sigBdt(v0)*((v0**3)/(v0**3 + vcT(n)**3))
!       Intf = sum * h
      
!       sigmaBdta(n) = Intf
!       sigmaBdtal(n) =LOG10(sigmaBdta(n))
!    END DO

!    contains
! !   function f(v)
! !     double precision f, v
! !     f = (v**3)/ (v**3 + vc(n)**3)   ! 式を定義
! !   end function f

!   function sigBdt(v)
!    ! real(8), intent(in) :: v
!    ! double precision SigBdt, A1, A2, A3, A4, A5, E
!    REAL(8), INTENT(IN) :: v
!    REAL(8) :: sigBdt
!    REAL(8) :: A1, A2, A3, A4, A5, E
!    !Using NRL fitting formula, substituing E with 1/2mv^2
!    ! sigBdt = (A5 + [(A4-A3*ETNBI)**2+1]**(-1)) / (ETNBI * [exp(A1/sqrt(ETNBI))-1])
!    A1 = 45.95
!    A2 = 5.02E4
!    A3 = 1.368E-2
!    A4 = 1.076
!    A5 = 409
!    E = 0.5d0*AMT*(v**2)/AEE/1000
!    ! sigBdt = (A5 + [(A4-A3*(1/2*AMT*v**2*(1/1.60218E16)))**2+1]**(-1)) / ((1/2*AMT*v**2*(1/1.60218E16)) * [exp(A1/sqrt(1/2*AMT*v**2*(1/1.60218E16)))-1])    ! 式を定義
!    sigBdt = (A5 + [(A4-A3*E)**2+1]**(-1)) / (E * [exp(A1/sqrt(E))-1])    ! 式を定義
! end function sigBdt

! ! set Spline coefficient
!    ! CALL SPL1D(tempal,sigmaBdtal,dsigmaBdtal,usigmaBdtal,10,0,ierr)

!   END FUNCTION sigmaDbuTbm

! Spline calculation
!   CALL SPL1DF(templ,sigmaB_dtl,tempal,usigmaBdtal,10,ierr)

  ! calculation of spline coefficient for reactivity of DbulkTNBI
  SUBROUTINE set_usigmaBdt
   Use trcomm !特定のものを使うのであればonlyを用いる
   USE libspl1d
   IMPLICIT NONE
   ! REAL(rkind),INTENT(IN):: temp, BENG
   REAL(rkind):: ETNBI, A1, A2, A3, A4, A5, v0, v
   ! REAL(rkind),DIMENSION(10):: sigmaBdta
   REAL(rkind),DIMENSION(10):: vcT
   REAL(rkind):: sum, h, Intf
   INTEGER :: i, n, k, ierr
   ! REAL(rkind),DIMENSION(10):: sigmaBdtal, dsigmaBdtal, usigmaBdtal, tempal
   REAL(rkind),DIMENSION(10):: dsigmaBdtal, usigmaBdtal, tempal

   ! 定数の定義
   ETNBI = 1000 !!!!!!!!!!!!!!!!!!!!!!!!!BENG !Beam energy in keV, 後ほどBENGに置き換える。BENGをtrcommで定義するようにする。!!!!!!!!!!!!!!!!!!!!!!!!
   v0 = sqrt((2*ETNBI*1000*AEE)/AMT) ! ETNBIから導出

   ! vcTの計算
   DO n = 1, 10
      vcT(n) = (3.0/4.0)*sqrt(pi)*(AME/AMT)**(1.0/3.0) * sqrt(2.0*tempa(n)*AEE/AME)
   END DO

   ! sigmaBdtaの計算
   DO n = 1, 10
      k = 10000
      h = v0 / real(k)
      sum = 0.0_rkind

      DO i = 1, k-1
         v = i * h
         sum = sum + sigBdt(v) * (v**3) / (v**3 + vcT(n)**3)
      END DO
      sum = sum + 0.5 * sigBdt(v0) * (v0**3) / (v0**3 + vcT(n)**3)

      Intf = sum * h
      sigmaBdta(n) = Intf
      sigmaBdtal(n) = LOG10(sigmaBdta(n))
   END DO

   ! SigmaBDT*distribution functionのSpline補間
   CALL SPL1D(tempal, sigmaBdtal, dsigmaBdtal, usigmaBdtal, 10, 0, ierr)
   contains
   ! sigBdt関数
   FUNCTION sigBdt(v)
      REAL(8), INTENT(IN) :: v
      REAL(8) :: sigBdt
      REAL(8) :: A1, A2, A3, A4, A5, E

      A1 = 45.95
      A2 = 5.02E4
      A3 = 1.368E-2
      A4 = 1.076
      A5 = 409
      E = 0.5d0 * AMT * v**2 / AEE / 1000.0d0
      sigBdt = (A5 + (A4 - A3*E)**2 + 1)**(-1) / (E * (EXP(A1/SQRT(E)) - 1))
   END FUNCTION sigBdt
END SUBROUTINE set_usigmaBdt

FUNCTION sigmaDbuTbm(temp)
   USE libspl1d
   IMPLICIT NONE
   REAL(rkind),INTENT(IN):: temp
   REAL(rkind):: templ
   REAL(rkind):: sigmaDbuTbm
   REAL(rkind):: sigmaDbuTbml
   INTEGER:: ierr
   templ=LOG10(temp)
   IF(templ.LE.tempal(1)) THEN
      sigmaDbuTbm=sigmaBdta(1)
   ELSE IF(templ.GE.tempal(10)) THEN
      sigmaDbuTbm=sigmaBdta(10)
   ELSE
      CALL SPL1DF(templ,sigmaDbuTbml,tempal,usigmaBdtal,10,ierr)
!       WRITE(6,'(5ES12.4)') temp,templ,sigmavml_dt,tempal(1),tempal(10)
      IF(ierr.NE.0) THEN
         WRITE(6,'(A,I4)') 'XX sigmavm_dt: SPL1DF error: ierr=',ierr
         STOP
      END IF
      sigmaDbuTbm = 10.D0**sigmaDbuTbml
   END IF
   RETURN
END FUNCTION sigmaDbuTbm


END MODULE libsigma
      
