      FUNCTION U_ERF(x)
!***********************************************************************
!U_ERF evaluates the error function
!  W.A.Houlberg 1/99
!Input:
!  x-argument of error function
!***********************************************************************
      USE task_kinds,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind)    U_ERF
!Declaration of input variables
      REAL(rkind)    x
!Declaration of local variables
      INTEGER        i
      REAL(rkind)    a,                       b,
     #               c,                       d,
     #               del,                     gln,
     #               h,                       x2
      gln=5.723649D-1
      x2=x**2
      U_ERF=0.D0
      IF(x2.lt.1.5D0) THEN
        a=0.5D0
        b=2.0D0
        del=b
        DO i=1,100
          a=a+1.D0
          del=del*x2/a
          b=b+del
          IF(ABS(del).lt.ABS(b)*3.0D-7) THEN
            U_ERF=b*EXP(-x2+0.5D0*LOG(x2)-gln)
            IF(x.lt.0.D0) U_ERF=-U_ERF
            GOTO 1000
          ENDIF
        ENDDO
      ELSE
        b=x2+0.5D0
        c=1.D0/1.0D-30
        d=1.D0/b
        h=d
        DO i=1,100
          a=-i*(i-0.5D0)
          b=b+2.D0
          d=a*d+b
          IF(ABS(d).lt.1.0D-30) d=1.0D-30
          c=b+a/c
          IF(ABS(c).lt.1.0D-30) c=1.0D-30
          d=1.D0/d
          del=d*c
          h=h*del
          IF(ABS(del-1.D0).lt.3.D-7) THEN
            IF(-x2+0.5D0*LOG(x2)-gln.LT.-100.D0) THEN
               U_ERF=0.D0
            ELSE
               U_ERF=1.D0-EXP(-x2+0.5D0*LOG(x2)-gln)*h
               IF(x.lt.0.D0) U_ERF=-U_ERF
            END IF
            GOTO 1000
          ENDIF
        ENDDO
      ENDIF
 1000 RETURN
      END
