! libpdkf2

MODULE libpdkf2

  USE task_kinds,ONLY: rkind
  USE task_constants,ONLY: PI

  REAL(rkind):: xi_func,alpha_func,beta_func,N_func
  INTEGER:: ntau_func
  INTEGER:: N0_eul=5
  INTEGER:: LMAX_eul=200
  REAL(rkind):: EPS_de=1.D-5
  REAL(rkind):: H0_de=0.5D0
  INTEGER:: ILST_de=0

  PRIVATE
  PUBLIC pdkf2

CONTAINS

  !     --- PDKF2: plasma dispersion kernel function ---
  !     \int_0^\infty \rd\tau \tau^{ntau-1}
  !          \exp^{-0.5D0*xi^2/\tau^2-0.5D0*\beta\tau^2-\imi\alpha\tau}

  FUNCTION pdkf2(xi,alpha,beta,ntau)
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi,alpha,beta
    INTEGER,INTENT(IN):: ntau
    COMPLEX(rkind):: pdkf2

    CALL pdkf_eul(xi,alpha,beta,ntau,pdkf2)

    RETURN
  END FUNCTION pdkf2

!     *****  EULER TRANSFOMATION of pdkf  *****

  SUBROUTINE pdkf_eul(xi,alpha,beta,ntau,pdkf)

    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi,alpha,beta
    INTEGER,INTENT(IN):: ntau
    COMPLEX(rkind),INTENT(OUT):: pdkf
    REAL(rkind):: SR,SI,SR1,SI1,SR2,SI2,SKR,SKI,PARITY
    INTEGER:: N,L,K
    REAL(rkind),ALLOCATABLE:: A(:),B(:)

    N0_eul=5
    ALLOCATE(A(LMAX_eul),B(LMAX_eul))

    SR1=0.D0
    SI1=0.D0
    IF(N0_eul.NE.0) THEN 
       DO N=N0_eul-1,0,-1
          CALL pdkf_elm(xi,alpha,beta,ntau,N,SR,SI)
          SR1=SR1+SR
          SI1=SI1+SI
       END DO
    ENDIF

    CALL pdkf_elm(xi,alpha,beta,ntau,N0_eul,SR,SI)
    A(1)=SR
    SR2=0.5D0*SR
    B(1)=SI
    SI2=0.5D0*SI
    
    PARITY=-1.D0
    L=0

30  CONTINUE
    L=L+1
    IF(L.GE.LMAX_eul) GO TO 9000
    CALL pdkf_elm(xi,alpha,beta,ntau,N0_eul+L,SR,SI)
    A(L+1)=SR*PARITY
    B(L+1)=SI*PARITY
    DO K=L,1,-1    
       A(K)=A(K+1)-A(K)
       B(K)=B(K+1)-B(K)
    END DO
    SKR=A(1)*PARITY*0.5D0**(L+1)
    SR2=SR2+SKR
    SKI=B(1)*PARITY*0.5D0**(L+1)
    SI2=SI2+SKI
    PARITY=-PARITY
    IF(DABS(SKR).GT.EPS_de.OR.DABS(SKI).GT.EPS_de) GOTO 30

    SR=(SR1+SR2)/SQRT(2.D0*PI)
    SI=(SI1+SI2)/SQRT(2.D0*PI)
    pdkf=DCMPLX(SR,SI)

    DEALLOCATE(A,B)
    RETURN

9000 CONTINUE
    WRITE(6,*) ' ## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
    RETURN
  END SUBROUTINE pdkf_eul

  !     *****  INTEGRAL over half period  *****

  SUBROUTINE pdkf_elm(xi,alpha,beta,ntau,N,SR,SI)

    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: xi,alpha,beta
    INTEGER,INTENT(IN):: ntau,N
    REAL(rkind),INTENT(OUT):: SR,SI
    REAL(rkind):: ER,EI
    
    xi_func=xi
    alpha_func=alpha
    beta_func=beta
    ntau_func=ntau
    N_func=N

    CALL DEFT(SR,ER,H0_de,EPS_de,ILST_de,FUNR_pdkf,'FUNR_pdkf')
    CALL DEFT(SI,EI,H0_de,EPS_de,ILST_de,FUNI_pdkf,'FUNI_pdkf')

    RETURN
  END SUBROUTINE pdkf_elm

!     *****  REAL PART  *****

    FUNCTION FUNR_pdkf(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNR_pdkf
      REAL(rkind):: Y1,T,T2,YY

      Y1=XM ! dummy
      IF(N_func.EQ.0) THEN 
         T=0.25D0*PI*XP/alpha_func   ! XP=0~2  T=0~0.5Pi
         T2=T*T
         YY=-0.5D0*(xi_func*xi_func/T2+beta_func*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(alpha_func*T).LT.1.D4) THEN
            FUNR_pdkf=T**ntau_func*0.25D0*PI*DEXP(YY)*DCOS(alpha_func*T)
         ELSE
            FUNR_pdkf=0.D0
         ENDIF 
      ELSE
         T=(0.5D0*PI*X+N_func*PI)/alpha_func ! X=-1~1,T=(N-0.5)*Pi~(N+0.5)*Pi
         T2=T*T
         YY=-0.5D0*(xi_func*xi_func/T2+beta_func*T2)
         IF(ABS(YY).LT.352.D0.AND.ABS(alpha_func).LT.1.D4) THEN
            FUNR_pdkf=T**ntau_func*0.5D0*PI*DEXP(YY)*DCOS(alpha_func*T)
         ELSE
            FUNR_pdkf=0.D0
         ENDIF
      ENDIF
      RETURN
    END FUNCTION FUNR_pdkf

!     *****  IMAG PART  *****

    FUNCTION FUNI_pdkf(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,XM,XP
      REAL(rkind):: FUNI_pdkf
      REAL(rkind):: Y1,Y2,T,T2,YY,AN2

      Y1=X ! dummy
      Y2=XM ! dummy
      T=(0.5D0*PI*XP+N_func*PI)/alpha_func ! XP=0-2, T=N*PI-(N+1)*PI
      T2=T*T
      YY=-0.5D0*(xi_func*xi_func/T2+beta_func*T2)
      IF(ABS(YY).LT.352.D0.AND.ABS(alpha_func*T).LT.1.D4) THEN
         FUNI_pdkf=T**ntau_func*0.5D0*PI*DEXP(YY)*DSIN(alpha_func*T)
      ELSE
         FUNI_pdkf=0.D0
      ENDIF
      RETURN
    END FUNCTION funi_pdkf

END MODULE libpdkf2
