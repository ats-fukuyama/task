! sak3.f90

MODULE sak3

  PRIVATE
  PUBLIC sak_3

CONTAINS

  SUBROUTINE sak_3

    USE sakcomm
    USE saksub
    USE libgrf
    IMPLICIT NONE
    REAL(dp):: wr,wi    ! real and imag parts of cw=omega/omegape
    REAL(dp):: rk       ! k vte/omegape = k lambda
    REAL(dp):: sg       ! sigma = l q_theta/k
    COMPLEX(dp):: cf    ! 
    INTEGER:: nxmax,nymax,nx,ny
    REAL(dp):: delta_nw,eps_nw,wwr,wwi,wra,wia,rd
    INTEGER:: lmax_nw,list_nw,mode_nw,ierr

    wr=1.D0
    wi=0.D0
    rk=0.1D0
    sg=0.D0

!        delta_nw: Step size in evaluating derivatives in Newton method
!        eps_nw  : Convergence criterion in Newton method
!        lmax_nw : Maximum iteration count in Newton method
!        list_nw : Listing in Newton method
!        mode_nw : Type of Newton method
    delta_nw = 1.D-6
    eps_nw = 1.D-8
    lmax_nw= 40
    list_nw= 1
    mode_nw= 0

1   CONTINUE
    WRITE(6,'(A/6ES12.4,2I4)') &
         '## INPUT: wr,wi,rk,sg,delta_nw,eps_nw,lmax_nw,list_nw?', &
         wr,wi,rk,sg,delta_nw,eps_nw,lmax_nw,list_nw
    READ(5,*,ERR=1,END=9000) wr,wi,rk,sg,delta_nw,eps_nw,lmax_nw,list_nw

    CALL set_rksg(rk,sg)
    CALL newtn0(subeps,wr,wi,wwr,wwi,rd,delta_nw,eps_nw,lmax_nw,list_nw,ierr)
    wra=SQRT((1.D0+sg**2)*(1.D0+3.D0*rk**2))
    wia=-SQRT(0.125D0*Pi)*SQRT(1.D0+sg**2)*(1.D0+3.D0*rk**2)**2/rk**3 &
         *EXP(-0.5D0/rk**2*(1.D0+3.D0*rk**2))
    WRITE(6,'(A,2F6.3,5ES12.4)') &
         'rk,sg:',rk,sg,wra,wia,wwr,wwi,rd
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE sak_3

END MODULE sak3
