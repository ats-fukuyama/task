MODULE T2PREC

  USE T2CNST,ONLY:ikind,rkind

  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE T2PREC_EXECUTE

    USE T2COMM,ONLY:NXMAX,NVMAX,NFMAX,Xvec,XvecIn,XvecOut,d2xout
    USE T2EXEC,ONLY:T2EXEC_EXECUTE
    USE T2GOUT,ONLY:T2_GOUT
    USE T2CONV,ONLY:T2CONV_EXECUTE
    USE T2COEF,ONLY:T2COEF_EXECUTE

    INTEGER(ikind)::cnt,ix,iv
    REAL(   rkind)::res

    cnt = 0
    
    DO 

       cnt = cnt + 1

       DO ix = 1, NXMAX
          DO iv = 1, NVMAX
             XvecIn(iv,ix) = Xvec(iv,ix)
          ENDDO
       ENDDO

       CALL T2COEF_EXECUTE
       !CALL T2EXEC_EXECUTE(NFMAX)
       !CALL T2CONV_EXECUTE(NFMAX,res)
       CALL T2EXEC_EXECUTE(6)
       CALL T2CONV_EXECUTE(6,res)

       DO ix = 1, NXMAX
          DO iv = 1, NVMAX
             Xvec(iv,ix) = XvecOut(iv,ix)
          ENDDO
       ENDDO
       
       IF(res.LE.1.D-7) EXIT

       IF(cnt.GE.100)THEN
          WRITE(6,*)'T2PREC did not converge'
          d2xout = Xvec
          CALL T2_GOUT
          STOP
       END IF
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2PREC_EXECUTE

END MODULE T2PREC
