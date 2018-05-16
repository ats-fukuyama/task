MODULE w1gout

CONTAINS

  SUBROUTINE w1_gout
    USE w1comm
    USE w1grf1,ONLY: w1gr2d,w1grud
    USE w1grf2,ONLY: w1gr1b,w1gr1d,w1gr1f,w1gr1h,w1gruf
    IMPLICIT NONE

    IF(NGRAPH.GT.0) THEN
       IF(NZPMAX.EQ.1) THEN
          CALL W1GR1D
          CALL W1GR1B
          CALL W1GR1F
          CALL W1GR1H
          CALL W1GRUD
          CALL W1GRUF
       ELSE
          CALL W1GR2D
          CALL W1GR1D
       ENDIF
    ENDIF
  END SUBROUTINE w1_gout
END MODULE w1gout
