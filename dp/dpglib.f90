! dpglib.f90

MODULE dpglib

  USE dpcomm,ONLY: rkind
  INTEGER:: nmax_a,nmax_b,nmax_c
  REAL(rkind),ALLOCATABLE:: f_a(:),urgb_ar(:,:),urgb_ag(:,:),urgb_ab(:,:)
  REAL(rkind),ALLOCATABLE:: f_b(:),urgb_br(:,:),urgb_bg(:,:),urgb_bb(:,:)
  REAL(rkind),ALLOCATABLE:: f_c(:),urgb_cr(:,:),urgb_cg(:,:),urgb_cb(:,:)

  PRIVATE
  PUBLIC set_rgbd,rgb_bard,rgbf_a,rgbf_b,rgbf_c,dpclip

CONTAINS

  SUBROUTINE set_rgbd(rgb)
    USE dpcomm,ONLY: rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rgb(3)
    EXTERNAL setrgb,POLY,LINES

    CALL setrgb(REAL(RGB(1)),REAL(RGB(2)),REAL(RGB(3)))
  END SUBROUTINE set_rgbd


!     ****** RGB COLOR PATERN ******

  SUBROUTINE rgb_bard(X1,X2,Y1,Y2,RGB,nmax,ind)

    USE dpcomm,ONLY: rkind
    IMPLICIT NONE
    REAL,INTENT(IN):: X1,X2,Y1,Y2
    INTEGER,INTENT(IN):: nmax,ind
    REAL(rkind),INTENT(IN):: RGB(3,nmax)
    REAL(rkind):: DX,DY,DXL,DYL
    REAL:: X(5),Y(5)
    INTEGER:: n
    EXTERNAL POLY,SETRGB,LINES

    IF(IND.EQ.0) THEN
       DX=(X2-X1)/DBLE(nmax)
       DY=0.0
       DXL=0.0
       DYL=Y2-Y1
    ELSE
       DX=0.0
       DY=(Y2-Y1)/DBLE(nmax)
       DXL=X2-X1
       DYL=0.0
    ENDIF
    DO n=1,nmax
       X(1)=dpclip(X1+DX*(n-1))
       Y(1)=dpclip(Y1+DY*(n-1))
       X(2)=dpclip(X1+DX*n+DXL)
       Y(2)=Y(1)
       X(3)=X(2)
       Y(3)=dpclip(Y1+DY*n+DYL)
       X(4)=X(1)
       Y(4)=Y(3)
       X(5)=X(1)
       Y(5)=Y(1)
       CALL set_rgbd(rgb(1:3,n))
       CALL POLY(X,Y,5)
    END DO

    X(1)=X1
    Y(1)=Y1
    X(2)=X2
    Y(2)=Y1
    X(3)=X2
    Y(3)=Y2
    X(4)=X1
    Y(4)=Y2
    X(5)=X1
    Y(5)=Y1
    CALL SETRGB(0.0,0.0,0.0)
    CALL LINES(X,Y,5)

    RETURN
  END SUBROUTINE rgb_bard

! --- rgb color pattern - positive ---

  SUBROUTINE rgbf_a(f,rgb)
    USE dpcomm,ONLY: rkind
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: f
    REAL(rkind),INTENT(OUT):: rgb(3)
    REAL(rkind),ALLOCATABLE:: rgb_a(:,:),dummy(:)
    INTEGER,SAVE:: init=0
    INTEGER:: ierr1,ierr2,ierr3

    IF(init.EQ.0) THEN
       nmax_a=11
       ALLOCATE(f_a(nmax_a),rgb_a(nmax_a,3),dummy(nmax_a))
       ALLOCATE(urgb_ar(4,nmax_a),urgb_ag(4,nmax_a),urgb_ab(4,nmax_a))
       f_a( 1)=0.0D0; rgb_a( 1,1:3)=(/0.4D0,0.9D0,0.4D0/)
       f_a( 2)=0.1D0; rgb_a( 2,1:3)=(/0.6D0,1.0D0,0.2D0/)
       f_a( 3)=0.2D0; rgb_a( 3,1:3)=(/0.8D0,1.0D0,0.2D0/)
       f_a( 4)=0.3D0; rgb_a( 4,1:3)=(/0.9D0,1.0D0,0.2D0/)
       f_a( 5)=0.4D0; rgb_a( 5,1:3)=(/1.0D0,1.0D0,0.2D0/)
       f_a( 6)=0.5D0; rgb_a( 6,1:3)=(/1.0D0,0.8D0,0.2D0/)
       f_a( 7)=0.6D0; rgb_a( 7,1:3)=(/1.0D0,0.6D0,0.2D0/)
       f_a( 8)=0.7D0; rgb_a( 8,1:3)=(/1.0D0,0.4D0,0.2D0/)
       f_a( 9)=0.8D0; rgb_a( 9,1:3)=(/1.0D0,0.2D0,0.4D0/)
       f_a(10)=0.9D0; rgb_a(10,1:3)=(/1.0D0,0.6D0,0.7D0/)
       f_a(11)=1.0D0; rgb_a(11,1:3)=(/1.0D0,1.0D0,1.0D0/)
       CALL SPL1D(f_a,rgb_a(1:nmax_a,1),dummy,urgb_ar,nmax_a,0,ierr1)
       CALL SPL1D(f_a,rgb_a(1:nmax_a,2),dummy,urgb_ag,nmax_a,0,ierr2)
       CALL SPL1D(f_a,rgb_a(1:nmax_a,3),dummy,urgb_ab,nmax_a,0,ierr3)
       IF(ierr1+ierr2+ierr3.NE.0) THEN
          WRITE(6,'(A,3I5)') &
               'XX rgb_a: SPL1D error: ierr1/2/3=',ierr1,ierr2,ierr3
          RETURN
       END IF
       DEALLOCATE(rgb_a,dummy)
       INIT=1
    END IF

    CALL SPL1DF(f,rgb(1),f_a,urgb_ar,nmax_a,ierr1)
    CALL SPL1DF(f,rgb(2),f_a,urgb_ag,nmax_a,ierr2)
    CALL SPL1DF(f,rgb(3),f_a,urgb_ab,nmax_a,ierr3)
    IF(ierr1+ierr2+ierr3.NE.0) THEN
       WRITE(6,'(A,3I5)') &
            'XX rgb_a: SPL1DF error: ierr1/2/3=',ierr1,ierr2,ierr3
       RETURN
    END IF
  END SUBROUTINE rgbf_a


! --- rgb color pattern - negative ---

  SUBROUTINE rgbf_b(f,rgb)
    USE dpcomm,ONLY: rkind
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: f
    REAL(rkind),INTENT(OUT):: rgb(3)
    REAL(rkind),ALLOCATABLE:: rgb_b(:,:),dummy(:)
    INTEGER,SAVE:: init=0
    INTEGER:: ierr1,ierr2,ierr3

    IF(init.EQ.0) THEN
       nmax_b=11
       ALLOCATE(f_b(nmax_b),rgb_b(nmax_b,3),dummy(nmax_b))
       ALLOCATE(urgb_br(4,nmax_b),urgb_bg(4,nmax_b),urgb_bb(4,nmax_b))
       f_b( 1)=0.0D0; rgb_b( 1,1:3)=(/1.0D0,1.0D0,1.0D0/)
       f_b( 2)=0.1D0; rgb_b( 2,1:3)=(/0.7D0,0.6D0,1.0D0/)
       f_b( 3)=0.2D0; rgb_b( 3,1:3)=(/0.4D0,0.2D0,1.0D0/)
       f_b( 4)=0.3D0; rgb_b( 4,1:3)=(/0.2D0,0.4D0,1.0D0/)
       f_b( 5)=0.4D0; rgb_b( 5,1:3)=(/0.2D0,0.6D0,1.0D0/)
       f_b( 6)=0.5D0; rgb_b( 6,1:3)=(/0.2D0,0.8D0,1.0D0/)
       f_b( 7)=0.6D0; rgb_b( 7,1:3)=(/0.2D0,1.0D0,1.0D0/)
       f_b( 8)=0.7D0; rgb_b( 8,1:3)=(/0.2D0,1.0D0,0.9D0/)
       f_b( 9)=0.8D0; rgb_b( 9,1:3)=(/0.2D0,1.0D0,0.8D0/)
       f_b(10)=0.9D0; rgb_b(10,1:3)=(/0.2D0,1.0D0,0.6D0/)
       f_b(11)=1.0D0; rgb_b(11,1:3)=(/0.4D0,0.9D0,0.4D0/)
       CALL SPL1D(f_b,rgb_b(1:nmax_b,1),dummy,urgb_br,nmax_b,0,ierr1)
       CALL SPL1D(f_b,rgb_b(1:nmax_b,2),dummy,urgb_bg,nmax_b,0,ierr2)
       CALL SPL1D(f_b,rgb_b(1:nmax_b,3),dummy,urgb_bb,nmax_b,0,ierr3)
       IF(ierr1+ierr2+ierr3.NE.0) THEN
          WRITE(6,'(A,3I5)') &
               'XX rgb_b: SPL1D error: ierr1/2/3=',ierr1,ierr2,ierr3
          RETURN
       END IF
       DEALLOCATE(rgb_b,dummy)
       INIT=1
    END IF

    CALL SPL1DF(f,rgb(1),f_b,urgb_br,nmax_b,ierr1)
    CALL SPL1DF(f,rgb(2),f_b,urgb_bg,nmax_b,ierr2)
    CALL SPL1DF(f,rgb(3),f_b,urgb_bb,nmax_b,ierr3)
    IF(ierr1+ierr2+ierr3.NE.0) THEN
       WRITE(6,'(A,3I5)') &
            'XX rgb_b: SPL1DF error: ierr1/2/3=',ierr1,ierr2,ierr3
       RETURN
    END IF
  END SUBROUTINE rgbf_b


! --- rgb color pattern - positive and negative ---

  SUBROUTINE rgbf_c(f,rgb)
    USE dpcomm,ONLY: rkind
    USE libspl1d
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: f
    REAL(rkind),INTENT(OUT):: rgb(3)
    REAL(rkind),ALLOCATABLE:: rgb_c(:,:),dummy(:)
    INTEGER,SAVE:: init=0
    INTEGER:: ierr1,ierr2,ierr3

    IF(init.EQ.0) THEN
       nmax_c=21
       ALLOCATE(f_c(nmax_c),rgb_c(nmax_c,3),dummy(nmax_c))
       ALLOCATE(urgb_cr(4,nmax_c),urgb_cg(4,nmax_c),urgb_cb(4,nmax_c))
       f_c( 1)=0.00D0; rgb_c( 1,1:3)=(/1.0D0,1.0D0,1.0D0/)
       f_c( 2)=0.05D0; rgb_c( 2,1:3)=(/0.7D0,0.6D0,1.0D0/)
       f_c( 3)=0.10D0; rgb_c( 3,1:3)=(/0.4D0,0.2D0,1.0D0/)
       f_c( 4)=0.15D0; rgb_c( 4,1:3)=(/0.2D0,0.4D0,1.0D0/)
       f_c( 5)=0.20D0; rgb_c( 5,1:3)=(/0.2D0,0.6D0,1.0D0/)
       f_c( 6)=0.25D0; rgb_c( 6,1:3)=(/0.2D0,0.8D0,1.0D0/)
       f_c( 7)=0.30D0; rgb_c( 7,1:3)=(/0.2D0,1.0D0,1.0D0/)
       f_c( 8)=0.35D0; rgb_c( 8,1:3)=(/0.2D0,1.0D0,0.9D0/)
       f_c( 9)=0.40D0; rgb_c( 9,1:3)=(/0.2D0,1.0D0,0.8D0/)
       f_c(10)=0.45D0; rgb_c(10,1:3)=(/0.2D0,1.0D0,0.6D0/)
       f_c(11)=0.50D0; rgb_c(11,1:3)=(/0.4D0,0.9D0,0.4D0/)
       f_c(12)=0.55D0; rgb_c(12,1:3)=(/0.6D0,1.0D0,0.2D0/)
       f_c(13)=0.60D0; rgb_c(13,1:3)=(/0.8D0,1.0D0,0.2D0/)
       f_c(14)=0.65D0; rgb_c(14,1:3)=(/0.9D0,1.0D0,0.2D0/)
       f_c(15)=0.70D0; rgb_c(15,1:3)=(/1.0D0,1.0D0,0.2D0/)
       f_c(16)=0.75D0; rgb_c(16,1:3)=(/1.0D0,0.8D0,0.2D0/)
       f_c(17)=0.8D00; rgb_c(17,1:3)=(/1.0D0,0.6D0,0.2D0/)
       f_c(18)=0.85D0; rgb_c(18,1:3)=(/1.0D0,0.4D0,0.2D0/)
       f_c(19)=0.90D0; rgb_c(19,1:3)=(/1.0D0,0.2D0,0.4D0/)
       f_c(20)=0.95D0; rgb_c(20,1:3)=(/1.0D0,0.6D0,0.7D0/)
       f_c(21)=1.00D0; rgb_c(21,1:3)=(/1.0D0,1.0D0,1.0D0/)
       CALL SPL1D(f_c,rgb_c(1:nmax_c,1),dummy,urgb_cr,nmax_c,0,ierr1)
       CALL SPL1D(f_c,rgb_c(1:nmax_c,2),dummy,urgb_cg,nmax_c,0,ierr2)
       CALL SPL1D(f_c,rgb_c(1:nmax_c,3),dummy,urgb_cb,nmax_c,0,ierr3)
       IF(ierr1+ierr2+ierr3.NE.0) THEN
          WRITE(6,'(A,3I5)') &
               'XX rgb_c: SPL1D error: ierr1/2/3=',ierr1,ierr2,ierr3
          RETURN
       END IF
       DEALLOCATE(rgb_c,dummy)
       INIT=1
    END IF

    CALL SPL1DF(f,rgb(1),f_c,urgb_cr,nmax_c,ierr1)
    CALL SPL1DF(f,rgb(2),f_c,urgb_cg,nmax_c,ierr2)
    CALL SPL1DF(f,rgb(3),f_c,urgb_cb,nmax_c,ierr3)
    IF(ierr1+ierr2+ierr3.NE.0) THEN
       WRITE(6,'(A,3I5)') &
            'XX rgb_c: SPL1DF error: ierr1/2/3=',ierr1,ierr2,ierr3
       RETURN
    END IF
  END SUBROUTINE rgbf_c


! ****** AVOID REAL*4 UNDERFLOW ******

  FUNCTION dpclip(D)

    USE dpcomm,ONLY: rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: D
    REAL:: dpclip

    IF(ABS(D).LT.1.D-30) then
       dpclip=0.0
    ELSEIF(D.GT. 1.D30) THEN
       dpclip= 1.E30
    ELSEIF(D.lT.-1.D30) THEN
       dpclip=-1.E30
    ELSE
       dpclip=SNGL(D)
    ENDIF
    RETURN
  END FUNCTION DPCLIP

END MODULE dpglib
