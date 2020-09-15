! wqgout.f90

MODULE wqgout

  PRIVATE
  PUBLIC wq_gout

CONTAINS

  SUBROUTINE wq_gout
    USE wqcomm
    USE wqparm,ONLY: wq_parm,wq_broadcast
    USE wqview,ONLY: wq_view
    USE libmpi
    IMPLICIT NONE
    CHARACTER(LEN=80):: line
    CHARACTER(LEN=1):: kid
    INTEGER:: mode

1   CONTINUE
    
    IF(nrank.EQ.0) THEN
       WRITE(6,'(A)') &
            '## wqgout menu: A,B  parm  X/end'
       CALL TASK_KLIN(line,kid,mode,wq_parm)
    END IF
    CALL ToUpper(kid)
    CALL mtx_broadcast1_character(kid)
    CALL mtx_broadcast1_integer(mode)

    IF(mode.EQ.2) CALL wq_broadcast
    IF(mode.NE.1) GO TO 1

    SELECT CASE(kid)
    CASE('A')
       CALL wq_gsub(W,dx,RR,omega,EX,EY,EZ,ne,OCE,OPE,OUH,AP,OR,OL,0)
    CASE('B')
       CALL wq_gsub(W,dx,RR,omega,EX,EY,EZ,ne,OCE,OPE,OUH,AP,OR,OL,1)
    CASE('X')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE wq_gout

  SUBROUTINE wq_gsub(W,dx,RR,omega,EX,EY,EZ,ne,OCE,OPE,OUH,AP,OR,OL,flag)

  use bpsd
  USE libgrf
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: W
  REAL(8),INTENT(IN)    :: dx,RR,omega,ne(W,W),OCE(W,W)
  REAL(8),INTENT(IN)    :: OPE(W,W),OUH(W,W),AP(W,W),OR(W,W),OL(W,W)
  COMPLEX(8),INTENT(IN) :: EX(W,W),EY(W,W),EZ(W,W)
  INTEGER(4)            :: i,j,k,NGP,NXM,NXMAX,NYMAX,MODE,IPRD,IG2D,flag
  CHARACTER(LEN=80)     :: STR
  REAL(8)               :: FX(W),FY(W),FZ(W,W)


  WRITE(6,'(A,I6)') 'wa_gsub: flag=',flag
  DO i=1,W
     FX(i)=RR+dx*DBLE(i-W/2)
  ENDDO
  DO j=1,W
     FY(j)=dx*DBLE(j-1-W/2)
  ENDDO

  SELECT CASE(flag)
  CASE(0)

  call PAGES
  
  NXM=W
  NXMAX=W
  NYMAX=W
  MODE=0
  IPRD=0
  IG2D=1

  NGP=5
  STR='/AP/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=AP(i,j)
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=6
  STR='/OUH/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=OUH(i,j)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=7
  STR='/R/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=OR(i,j)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=8
  STR='/L/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=OL(i,j)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=9
  STR='/OCE/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=OCE(i,j)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=10
  STR='/OPE/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=OPE(i,j)**2/omega**2-1.d0
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  NGP=11
  STR='/ne/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=ne(i,j)
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D,ASPECT=0.D0)

  call PAGEE

  CASE(1)

  call PAGES

  NXM=W
  NXMAX=W
  NYMAX=W
  MODE=0
  IPRD=0
  IG2D=1

  NGP=5
  STR='/Real(EX)/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=REAL(EX(i,j))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D)

  NGP=8
  STR='/Image(EX)/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=IMAG(EX(i,j))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D)
  
  NGP=6
  STR='/Real(EY)/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=Real(EY(i,j))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D)

  NGP=9
  STR='/Image(EY)/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=IMAG(EY(i,j))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D)
  
  NGP=7
  STR='/Real(EZ)/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=Real(EZ(i,j))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D)

  NGP=10
  STR='/Image(EZ)/'
  DO j=1,W
     DO i=1,W
        FZ(i,j)=IMAG(EZ(i,j))
     END DO
  END DO
  CALL GRD2D(NGP,FX,FY,FZ,NXM,NXMAX,NYMAX,STR,MODE,IPRD,IG2D)
  
  CALL PAGEE

   END SELECT
  
  RETURN
  END SUBROUTINE wq_gsub
END MODULE wqgout
