module tx_interface

  !****************!
  !   txmenu.f90   !
  !****************!

  interface
     SUBROUTINE TXKLIN(LINE,KID,MODE)
       integer(4), intent(out) :: MODE
       character(len=80), intent(out) :: LINE
       character(len=1), intent(out) :: KID
     end SUBROUTINE TXKLIN
  end interface

  !****************!
  !   txntv.f90    !
  !****************!

  interface
     subroutine perturb_mag(B_lambda)
       real(8), dimension(*), intent(out) :: B_lambda
     end subroutine perturb_mag
  end interface

  interface
     subroutine Wnm_spline(fmnq, wnm, umnq, nmnqm)
       integer(4), intent(in) :: nmnqm
       real(8), dimension(1:nmnqm), intent(in) :: fmnq, wnm
       real(8), dimension(1:4,1:nmnqm), intent(in) :: umnq
     end subroutine Wnm_spline
  end interface

  !****************!
  !   txmmm.f90    !
  !****************!

  interface
     subroutine txmmm95(dNedr,dNidr,dTedr,dTidr,dQdr,gamma)
       real(8), dimension(*), intent(in) :: dNedr,dNidr,dTedr,dTidr,dQdr
       real(8), dimension(*), intent(out), optional :: gamma
     end subroutine txmmm95
  end interface

  !****************!
  !   txlib.f90    !
  !****************!

  interface
     pure REAL(8) FUNCTION EXPV(X)
       REAL(8), INTENT(IN) :: X
     end FUNCTION EXPV
  end interface

  interface APTOS
     SUBROUTINE APITOS(STR, NSTR, I)
       character(len=*), INTENT(INOUT) :: STR
       INTEGER(4),       INTENT(INOUT) :: NSTR
       INTEGER(4),       INTENT(IN)    :: I
     end SUBROUTINE APITOS

     SUBROUTINE APSTOS(STR, NSTR, INSTR, NINSTR)
       character(len=*), INTENT(INOUT) :: STR
       INTEGER(4),       INTENT(INOUT) :: NSTR
       character(len=*), INTENT(IN)    :: INSTR
       INTEGER(4),       INTENT(IN)    :: NINSTR
     end SUBROUTINE APSTOS

     SUBROUTINE APDTOS(STR, NSTR, D, FORM)
       character(len=*), INTENT(INOUT) :: STR
       INTEGER(4),       INTENT(INOUT) :: NSTR
       REAL(8),          INTENT(IN)    :: D
       character(len=*), INTENT(IN)    :: FORM
     end SUBROUTINE APDTOS

     SUBROUTINE APRTOS(STR, NSTR, GR, FORM)
       character(len=*), INTENT(INOUT) :: STR
       INTEGER(4),       INTENT(INOUT) :: NSTR
       REAL(4),          INTENT(IN)    :: GR
       character(len=*), INTENT(IN)    :: FORM
     end SUBROUTINE APRTOS
  end interface

  interface
     SUBROUTINE TOUPPER(KTEXT)
       character(len=*), INTENT(INOUT) ::  KTEXT
     end SUBROUTINE TOUPPER
  end interface

  interface
     SUBROUTINE KSPLIT_TX(KKLINE,KID,KKLINE1,KKLINE2)
       CHARACTER(LEN=*),  INTENT(IN)  :: KKLINE
       CHARACTER(LEN=1),  INTENT(IN)  :: KID
       CHARACTER(LEN=*), INTENT(OUT) :: KKLINE1, KKLINE2
     end SUBROUTINE KSPLIT_TX
  end interface

  interface DERIVS
     SUBROUTINE DERIVS1(R,F,NRMAX,G)
       integer(4), intent(in) :: NRMAX
       real(8), dimension(0:NRMAX), intent(in)  :: R, F
       real(8), dimension(0:NRMAX), intent(out) :: G
     end SUBROUTINE DERIVS1

     SUBROUTINE DERIVS2(R,F,LQ,NQMAX,NRMAX,G)
       integer(4), intent(in) :: LQ, NRMAX, NQMAX
       real(8), dimension(0:NRMAX), intent(in)  :: R
       real(8), dimension(1:NQMAX,0:NRMAX), intent(in)  :: F
       real(8), dimension(0:NRMAX), intent(out) :: G
     end SUBROUTINE DERIVS2
  end interface

  interface
     pure REAL(8) FUNCTION DERIVF(NR,R,F,LQ,NQMAX,NRMAX)
       integer(4), intent(in) :: NR, LQ, NRMAX, NQMAX
       real(8), dimension(0:NRMAX), intent(in)  :: R
       real(8), dimension(1:NQMAX,0:NRMAX), intent(in)  :: F
     end FUNCTION DERIVF
  end interface

  interface
     function dfdx(x,f,nmax,mode)
       integer(4), intent(in) :: nmax, mode
       real(8), dimension(0:nmax), intent(in) :: x, f
       real(8), dimension(0:nmax) :: dfdx
     end function dfdx
  end interface

  interface
     REAL(8) FUNCTION INTG_F(X)
       real(8), dimension(*), intent(in) :: X
     end FUNCTION INTG_F
  end interface

  interface
     REAL(8) FUNCTION INTG_P(X,NR,ID)
       integer(4), intent(in) :: NR, ID
       real(8), dimension(*), intent(in) :: X
     end FUNCTION INTG_P
  end interface

  interface
     SUBROUTINE VALINT_SUB(X,NRLMAX,VAL,NR_START)
       integer(4), intent(in) :: NRLMAX
       integer(4), intent(in), optional :: NR_START
       real(8), dimension(*), intent(in) :: X
       real(8), intent(out) :: VAL
     end SUBROUTINE VALINT_SUB
  end interface

  interface
     SUBROUTINE INTDERIV3(X,R,intX,FVAL,NRMAX,ID)
       integer(4), intent(in) :: NRMAX, ID
       real(8), intent(in), dimension(0:NRMAX) :: X, R
       real(8), intent(in) :: FVAL
       real(8), intent(out), dimension(0:NRMAX) :: intX
     end SUBROUTINE INTDERIV3
  end interface

  interface
     pure REAL(8) FUNCTION LORENTZ(R,C1,C2,W1,W2,RC1,RC2,AMP)
       real(8), intent(in) :: r, c1, c2, w1, w2, rc1, rc2
       real(8), intent(in), optional :: AMP
     end FUNCTION LORENTZ
  end interface

  interface
     pure REAL(8) FUNCTION LORENTZ_PART(R,W1,W2,RC1,RC2,ID)
       real(8), intent(in) :: r, w1, w2, rc1, rc2
       integer(4), intent(in) :: ID
     end FUNCTION LORENTZ_PART
  end interface

  interface
     SUBROUTINE BISECTION(f,cl1,cl2,w1,w2,rc1,rc2,amp,s,valmax,val,valmin)
       real(8), external :: f
       real(8), intent(in) :: cl1, cl2, w1, w2, rc1, rc2, amp, s, valmax
       real(8), intent(in), optional :: valmin
       real(8), intent(out) :: val
     end SUBROUTINE BISECTION
  end interface

  interface
     subroutine inexpolate(nmax_in,r_in,dat_in,nmax_std,r_std,iedge,dat_out,ideriv,nrbound,idx)
       integer(4), intent(in) :: nmax_in, nmax_std, iedge
       integer(4), intent(in),  optional :: ideriv, idx
       integer(4), intent(out), optional :: nrbound
       real(8), dimension(1:nmax_in), intent(in) :: r_in, dat_in
       real(8), dimension(0:nmax_std), intent(in) :: r_std
       real(8), dimension(0:nmax_std), intent(out) :: dat_out
     end subroutine inexpolate
  end interface

  interface
     pure real(8) function fgaussian(x,mu,sigma,norm)
       real(8), intent(in) :: x, mu, sigma
       integer(4), intent(in), optional :: norm
     end function fgaussian
  end interface

  interface
     pure real(8) function moving_average(i,f,imax,iend)
       integer(4), intent(in) :: i, imax
       integer(4), intent(in), optional :: iend
       real(8), dimension(0:imax), intent(in) :: f
     end function moving_average
  end interface

  !*****************!
  !   txmisc.f90    !
  !*****************!

  interface
     pure REAL(8) FUNCTION TRCOFS(S,ALFA,RKCV)
       real(8), intent(in) :: S, ALFA, RKCV
     end FUNCTION TRCOFS
  end interface

  interface
     pure REAL(8) FUNCTION CORR(X)
       real(8), intent(in) :: X
     end FUNCTION CORR
  end interface

  !**********************!
  !   coulomb_log.f90    !
  !**********************!

  interface
     function coulog( zeff, ne, te, ti, A1, Z1, A2, Z2, tb ) result( lambda )
       real(8), intent(in) :: zeff, ne, te, ti, A1, Z1, A2, Z2
       real(8), intent(in), optional :: tb
       real(8) :: lambda
     end function coulog
  end interface

  interface
     function coulog_gen( ne, te, CDi, A1, Z1, t1, A2, Z2, t2 ) result( lambda )
       real(8), intent(in) :: ne, te, CDi, A1, Z1, t1, A2, Z2, t2
       real(8) :: lambda
     end function coulog_gen
  end interface

  interface
     real(8) function coulog_NRL(imodel, Ne, Te, Ni, Ti, PA, PZ) result(f)
       integer(4), intent(in) :: imodel
       real(8), intent(in) :: Ne, Te
       real(8), intent(in), optional :: Ni, Ti, PA, PZ
     end function coulog_NRL
  end interface

  !****************!
  !   txfile.f90   !
  !****************!

  interface
     SUBROUTINE TXLOAD(IST)
       integer(4), intent(out) :: IST
     end SUBROUTINE TXLOAD
  end interface

  interface
     SUBROUTINE TXGLOD(IST)
       integer(4), intent(out) :: IST
     end SUBROUTINE TXGLOD
  end interface

  interface
     REAL(8) FUNCTION rLINEAVE(Rho)
       REAL(8), INTENT(IN) :: Rho
     end FUNCTION rLINEAVE
  end interface

  interface
     integer(4) function detect_datatype(kchar)
       character(len=*), intent(in) :: kchar
     end function detect_datatype
  end interface

  interface
     subroutine initprof_input(nr, idx, out)
       integer(4), optional :: nr, idx
       real(8), optional :: out
     end subroutine initprof_input
  end interface

  !***************************!
  !   txg3d.f90, txg2d.f90    !
  !***************************!

  interface
     subroutine TXGRUR(GX,GTX,GYL,NRMAX,NGT,NGTM)!,STR,KV,INQ)
       integer(4),                  intent(in) :: NRMAX, NGT, NGTM!, INQ
       real(4), dimension(0:NRMAX), intent(in) :: GX
       real(4), dimension(0:NGT),   intent(in) :: GTX
       real(4), dimension(0:NRMAX,0:NGTM), intent(in) :: GYL
!       CHARACTER(LEN=80),INTENT(IN):: STR, KV
     end subroutine TXGRUR
  end interface

  interface
     subroutine TXGRURA(GX,GTX,GYL,NRMAX,NGT,NGTM)!,STR,KV,INQ)
       integer(4),                  intent(in) :: NRMAX, NGT, NGTM!, INQ
       real(4), dimension(0:NRMAX), intent(in) :: GX
       real(4), dimension(0:NGT),   intent(in) :: GTX
       real(4), dimension(0:NRMAX,0:NGTM), intent(in) :: GYL
!       CHARACTER(LEN=80),INTENT(IN):: STR, KV
     end subroutine TXGRURA
  end interface

  interface
     subroutine TXGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KV,MODE)
       REAL(4),    INTENT(IN) :: GX1, GX2, GY1, GY2
       INTEGER(4), INTENT(IN) :: NXM, NXMAX, NYMAX, MODE
       REAL(4), DIMENSION(NXMAX),     INTENT(IN) :: GX
       REAL(4), DIMENSION(NYMAX),     INTENT(IN) :: GY
       REAL(4), DIMENSION(NXM,NYMAX), INTENT(IN) :: GZ
       CHARACTER(LEN=80) :: STR, KV
     end subroutine TXGR3D
  end interface

end module tx_interface
