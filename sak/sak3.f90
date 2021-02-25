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
    REAL(dp):: delta_nw,eps_nw,wr1,wi1,wi2,wim1,wim2,rd
    INTEGER:: lmax_nw,list_nw,mode_nw,ierr

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
    WRITE(6,'(A/4ES12.4,2I4)') &
         '## INPUT: rk,sg,delta_nw,eps_nw,lmax_nw,list_nw?', &
         rk,sg,delta_nw,eps_nw,lmax_nw,list_nw
    READ(5,*,ERR=1,END=9000) rk,sg,delta_nw,eps_nw,lmax_nw,list_nw

    CALL cwaprx(rk,sg,wr1,wi1,wi2,wim1,wim2)
    CALL set_rksg(rk,sg)
    CALL newtn0(subeps,wr1,wi2,wr,wi,rd,delta_nw,eps_nw,lmax_nw,list_nw,ierr)
    WRITE(6,'(A,2F6.3,5ES12.4)') &
         'rk,sg:',rk,sg,wr1,wi2,wr,wi,rd
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE sak_3

END MODULE sak3
