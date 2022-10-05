! libpol.f90

MODULE libpol

  USE TASK_KINDS,ONLY: dp
  PRIVATE
  PUBLIC polintx    ! real version
  PUBLIC polintn    ! integer version
  PUBLIC polintn1
  PUBLIC polintn2
  PUBLIC polintn_old

CONTAINS

  !  *** polynomial interpolation ***

  !       Neville's algorithm
  !         f(x_i)=y_i for i=0..n
  !       f_i_j: polynomial of order j-i satisfies for i..j
  !         f_i_i(x)=y_i  for 0<=i<=n
  !         f_i_j(x)=((x-x_j)p_i_j-1(x)-(x-x_i)p_i+1_j(x))/(x_i-x_j)
  !                       for 0<=i<j<=n
  !       f_0_0
  !       f_1_1  f_0_1
  !       f_2_2  f_1_2  f_0_2
  !       f_3_3  f_2_3  f_1_3  f_0_3
  !       f_4_4  f_3_4  f_2_4  f_1_4  f_0_4

  SUBROUTINE polintx(xa,fa,n,x,f)
    IMPLICIT NONE
    INTEGER,INTENT(IN):: n
    REAL(dp),INTENT(IN):: xa(n),fa(n),x
    REAL(dp),INTENT(OUT):: f
    REAL(dp):: den,dif,dift,ho,hp,w,c(n),d(n),dy
    INTEGER:: i,m,ns

    ns=1
    dif=ABS(x-xa(1))
    DO i=1,n
       dift=ABS(x-xa(i))
       IF (dift.LT.dif) THEN
          ns=i
          dif=dift
       END IF
       c(i)=fa(i)
       d(i)=fa(i)
    END DO

    f=fa(ns)
    ns=ns-1
    DO m=1,n-1
       DO i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          IF(den.eq.0.D0) THEN
             WRITE(6,*) 'failure in polintx'
             STOP
          ENDIF
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       END DO
       IF (2*ns.LT.n-m) THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       f=f+dy
    END DO
    RETURN
  END SUBROUTINE polintx

  ! **********************************************
  
  SUBROUTINE polintn(nra,psa,n,nr,ps)
    IMPLICIT NONE
    INTEGER,INTENT(IN):: n,nr
    INTEGER,INTENT(IN):: nra(n)
    REAL(dp),INTENT(IN):: psa(n)
    REAL(dp),INTENT(OUT):: ps
    REAL(dp):: den,dif,dift,w,c(n),d(n),dy
    INTEGER:: i,m,ns,nho,nhp,ndif,ndift,nden

    ns=1
    ndif=ABS(nr-nra(1))
    DO i=1,n
       ndift=ABS(nr-nra(i))
       IF (ndift.LT.ndif) THEN
          ns=i
          ndif=ndift
       END IF
       c(i)=psa(i)
       d(i)=psa(i)
    END DO

    ps=psa(ns)
    ns=ns-1
    DO m=1,n-1
       DO i=1,n-m
          nho=nra(i)-nr
          nhp=nra(i+m)-nr
          w=c(i+1)-d(i)
          nden=nho-nhp
          IF(nden.eq.0.D0) THEN
             WRITE(6,*) 'failure in polint'
             STOP
          ENDIF
          den=w/DBLE(nden)
          d(i)=nhp*den
          c(i)=nho*den
       END DO
       IF (2*ns.LT.n-m) THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       END IF
       ps=ps+dy
    END DO
    RETURN
  END SUBROUTINE polintn

  SUBROUTINE polintn1(nr,npmax,nrm,data)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr,npmax,nrm
    REAL(dp),DIMENSION(nrm),INTENT(INOUT):: data
    INTEGER,DIMENSION(npmax):: nra
    REAL(dp),DIMENSION(npmax):: datapa
    INTEGER:: np

    DO np=1,npmax
       nra(np)=nr-npmax-1+np
       datapa(np)=data(nr-npmax-1+np)
    ENDDO

    CALL polintn(nra,datapa,npmax,nr,data(nr)) 
    RETURN
  END SUBROUTINE polintn1

! **********************************************

  SUBROUTINE polintn2(nr,nthmax,npmax,nthm,nrm,data)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr,nthmax,npmax,nthm,nrm
    REAL(dp),DIMENSION(nthm,nrm),INTENT(INOUT):: data
    INTEGER,DIMENSION(npmax):: nra
    REAL(dp),DIMENSION(npmax):: datapa
    INTEGER:: nth,np

    DO nth=1,nthmax
       DO np=1,npmax
          nra(np)=nr-npmax-1+np
          datapa(np)=data(nth,nr-npmax-1+np)
       END DO

       CALL polintn(nra,datapa,npmax,nr,data(nth,nr)) 
    END DO
    RETURN
  END SUBROUTINE polintn2

! **********************************************

  SUBROUTINE polintn1_old(nr,npmax,nrm,data)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr,npmax,nrm
    REAL(dp),DIMENSION(nrm),INTENT(INOUT):: data
    REAL(dp),DIMENSION(npmax):: nra,datapa
    INTEGER:: np
    REAL(dp):: dy

    DO np=1,npmax
       nra(np)=nr-npmax-1+np
       datapa(np)=data(nr-npmax-1+np)
    ENDDO
      
    CALL polintn_old(nra,datapa,npmax,nr,data(nr),dy) 
    RETURN
  END SUBROUTINE polintn1_old

  SUBROUTINE polintn2_old(nr,nthmax,npmax,nthm,nrm,data)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr,nthmax,npmax,nthm,nrm
    REAL(dp),DIMENSION(nthm,nrm),INTENT(INOUT):: data
    REAL(dp),DIMENSION(npmax):: nra,datapa
    INTEGER:: nth,np
    REAL(dp):: dy

    DO nth=1,nthmax
       DO np=1,npmax
          nra(np)=nr-npmax-1+np
          datapa(np)=data(nth,nr-npmax-1+np)
       ENDDO

       CALL polintn_old(nra,datapa,npmax,nr,data(nth,nr),dy) 
    ENDDO
    RETURN
  END SUBROUTINE polintn2_old

  SUBROUTINE polintn_old(nra,psa,n,nr,ps,dy)
    IMPLICIT NONE
    INTEGER:: n,nr
    REAL(dp):: dy,nra(n),ps,psa(n)
    INTEGER:: i,m,ns
    REAL(dp):: den,dif,dift,ho,hp,w,c(n),d(n)

    ns=1
    dif=ABS(nr-nra(1))
    DO i=1,n
       dift=ABS(nr-nra(i))
       IF (dift.LT.dif) THEN
          ns=i
          dif=dift
       ENDIF
       c(i)=psa(i)
       d(i)=psa(i)
    END DO
    ps=psa(ns)
    ns=ns-1
    DO m=1,n-1
       DO i=1,n-m
          ho=nra(i)-nr
          hp=nra(i+m)-nr
          w=c(i+1)-d(i)
          den=ho-hp
          IF(den.EQ.0.) THEN
             WRITE(6,*) 'XX polint: failure'
             STOP
          ENDIF
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       END DO
       IF (2*ns.LT.n-m) THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       ENDIF
       ps=ps+dy
    END DO
    RETURN
  END SUBROUTINE polintn_old
END MODULE libpol
