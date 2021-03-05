  ! testpol.f90

PROGRAM testpol

  USE task_kinds,ONLY: dp
  USE libpol
  USE libgrf
  IMPLICIT NONE

  INTEGER,ALLOCATABLE:: na(:)
  REAL(dp),ALLOCATABLE:: fa(:),xa(:)
  INTEGER:: n,nb,i,mode
  REAL(dp):: xb,fb1,fb2,fb3,dy
  CHARACTER(LEN=20):: K
  EXTERNAL GSOPEN,PAGES,PAGEE,GSCLOS

  CALL GSOPEN
  mode=1
  n=3

1 CONTINUE
  WRITE(6,'(A,2I4)') '## INPUT mode,n:',mode,n
  READ(5,*,ERR=1,END=9000) mode,n
  IF(mode.EQ.0) GOTO 9000
  IF(n.LT.3) n=3
  
  IF(ALLOCATED(na)) DEALLOCATE(na,fa,xa)
  ALLOCATE(na(n),fa(n+1),xa(n+1))
  DO i=1,n
     na(i)=i-1
     xa(i)=DBLE(na(i))
     fa(i)=DBLE(i-1)
  END DO
  nb=n
  xb=DBLE(n)

2 CONTINUE
  SELECT CASE(mode)
  CASE(1)
     WRITE(K,'(A,I2,A)') '(A,',n,'I4)'
     WRITE(6,K) '## INPUT na:',(na(i),i=1,n)
     READ(5,*,ERR=2,END=1) (na(i),i=1,n)
     IF(na(1).EQ.na(2)) GO TO 1
     DO i=1,n
        xa(i)=DBLE(na(i))
     END DO
  CASE(2)
     WRITE(K,'(A,I2,A)') '(A,',n,'ES12.4)'
     WRITE(6,K) '## INPUT xa:',(xa(i),i=1,n)
     READ(5,*,ERR=2,END=1) (xa(i),i=1,n)
     IF(xa(1).EQ.xa(2)) GO TO 1
  END SELECT
3 CONTINUE
  WRITE(K,'(A,I2,A)') '(A,',n,'ES12.4)'
  WRITE(6,K) '## INPUT fa:',(fa(i),i=1,n)
  READ(5,*,ERR=3,END=2) (fa(i),i=1,n)
  IF(fa(1).EQ.0.D0.AND.fa(2).EQ.0.D0.AND.fa(3).EQ.0.D0) GO TO 2
4 CONTINUE
  SELECT CASE(mode)
  CASE(1)
     WRITE(6,'(A,I2)') '## INPUT nb:',nb
     READ(5,*,ERR=4,END=2) nb
     IF(nb.EQ.0) GO TO 2
     xb=DBLE(nb)
  CASE(2)
     WRITE(6,'(A,ES12.4)') '## INPUT xb:',xb
     READ(5,*,ERR=4,END=2) xb
     IF(xb.EQ.0.D0) GO TO 2
  END SELECT

  SELECT CASE(mode)
  CASE(1)
     CALL polintn(na,fa,n,nb,fb1)
     CALL polintx(xa,fa,n,xb,fb2)
     CALL polintn_old(xa,fa,n,nb,fb3,dy)
  CASE(2)
     CALL polintx(xa,fa,n,xb,fb2)
     fb1=fb2
     fb3=fb2
  END SELECT
  
  DO i=1,n
     WRITE(6,'(A,I4,2ES12.4)') 'n,x,f=',i,xa(i),fa(i)
  END DO
  WRITE(6,'(A,4X,4ES12.4)') 'n,x,f=',xb,fb1,fb2,fb3
  
  CALL PAGES

  xa(n+1)=xb
  fa(n+1)=fb1
  CALL grd1d(1,xa,fa,n+1,n+1,1,'@polintn@')
  fa(n+1)=fb2
  CALL grd1d(2,xa,fa,n+1,n+1,1,'@polintx@')
  fa(n+1)=fb3
  CALL grd1d(3,xa,fa,n+1,n+1,1,'@polintn_old@')

  CALL PAGEE

  GO TO 2

9000 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM testpol
