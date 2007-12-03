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

  interface DERIVS
     SUBROUTINE DERIVS1D(R,F,NRMAX,G)
       real(8), dimension(0:NRMAX), intent(in)  :: R, F
       real(8), dimension(0:NRMAX), intent(out) :: G
       integer(4), intent(in) :: NRMAX
     end SUBROUTINE DERIVS1D

     SUBROUTINE DERIVS2D(R,F,LQ,NQMAX,NRMAX,G)
       real(8), dimension(0:NRMAX), intent(in)  :: R
       real(8), dimension(1:NQMAX,0:NRMAX), intent(in)  :: F
       real(8), dimension(0:NRMAX), intent(out) :: G
     end SUBROUTINE DERIVS2D
  end interface

  interface
     pure REAL(8) FUNCTION DERIVF(NR,R,F,LQ,NQMAX,NRMAX)
       real(8), dimension(0:NRMAX), intent(in)  :: R
       real(8), dimension(1:NQMAX,0:NRMAX), intent(in)  :: F
       integer(4), intent(in) :: NR, LQ, NRMAX, NQMAX
     end FUNCTION DERIVF
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
     pure REAL(8) FUNCTION TRCOFS(S,ALFA,RKCV)
       real(8), intent(in) :: S, ALFA, RKCV
     end FUNCTION TRCOFS
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

  !****************!
  !   txfile.f90   !
  !****************!

  interface
     REAL(8) FUNCTION rLINEAVE(Rho)
       REAL(8), INTENT(IN) :: Rho
     end FUNCTION rLINEAVE
  end interface

  !****************!
  !   txg3d.f90    !
  !****************!

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
     subroutine TXGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KV,MODE)
       REAL(4),    INTENT(IN) :: GX1, GX2, GY1, GY2
       INTEGER(4), INTENT(IN) :: NXM, NXMAX, NYMAX, MODE
       REAL(4), DIMENSION(NXMAX),     INTENT(IN) :: GX
       REAL(4), DIMENSION(NYMAX),     INTENT(IN) :: GY
       REAL(4), DIMENSION(NXM,NYMAX), INTENT(IN) :: GZ
     end subroutine TXGR3D
  end interface

end module tx_interface
