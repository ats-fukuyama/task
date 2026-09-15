MODULE w1gout

CONTAINS

  SUBROUTINE w1_gout
    USE w1comm
    USE w1grf1,ONLY: w1gr2dw,w1gr2dr,w1grud,w1gr1dj
    USE w1grf2,ONLY: w1gr1b,w1gr1d,w1gr1f,w1gr1h,w1gruf
    IMPLICIT NONE
    INTEGER:: NID

1   WRITE(6,*) '## Input choice of plot: 1-6 for 1D, 1-4 for 2D, 0 for end'
    READ(5,*,END=9000,ERR=1) NID
    IF(NID.EQ.0) GO TO 9000
    IF(NZMAX.EQ.1) THEN
       SELECT CASE(NID)
       CASE(1)
          CALL W1GR1D
       CASE(2)          
          CALL W1GR1B
       CASE(3)          
          CALL W1GR1F
       CASE(4)          
          CALL W1GR1H
       CASE(5)          
          CALL W1GRUD
       CASE(6)          
          CALL W1GRUF
       END SELECT
    ELSE
       SELECT CASE(NID)
       CASE(1)
          CALL W1GR2DW(1)
       CASE(2)
          CALL W1GR2DW(2)
       CASE(3)
          CALL W1GR2DW(3)
       CASE(4)
          CALL W1GR2DR
       END SELECT
    ENDIF
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE w1_gout
END MODULE w1gout
