MODULE w1fflr

  USE w1comm,ONLY: rkind,PI
  IMPLICIT NONE
  REAL(rkind):: G1
  INTEGER:: NF,NN,NM

CONTAINS

!     ****** SLAVE FUNTION FOR DE INTEGRATION ******

  FUNCTION W1FNFG(X,XM,XP)
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,XM,XP
    REAL(rkind):: W1FNFG
    REAL(rkind):: SP=2.506628275D0
    REAL(rkind):: DUMMY,T,S,XX

    DUMMY=X
    DUMMY=XM
    T=0.5D0*PI*XP
    S=SIN(T)
    XX=G1/(2.82843D0*S)
    IF(ABS(XX).LT.10.D0) THEN
       IF(NF.EQ.0) THEN
          W1FNFG = 2.D0*S*S*S*(-XX*ERFC(XX)+EXP(-XX*XX)/SP)
       ELSEIF(NF.EQ.1) THEN
          W1FNFG = 0.5D0*S**NM*EXP(-XX*XX)*COS(2*NN*T)/SP
       ELSEIF(NF.EQ.2) THEN
          W1FNFG = 0.5D0*COS(T)*S**NM*EXP(-XX*XX)*SIN(2*NN*T)/SP
       ELSEIF(NF.EQ.3) THEN
          W1FNFG = 0.5D0*S**(NM+1)*ERFC(XX)*COS(2*NN*T)
       ELSEIF(NF.EQ.4) THEN
          W1FNFG = 0.5D0*COS(T)*S**(NM+1)*ERFC(XX)*SIN(2*NN*T)
       ENDIF
    ELSE
       W1FNFG = 0.D0
    ENDIF
  END FUNCTION W1FNFG
  
!     ****** FUNCTION FOR INTEGRO-DIFFERENTIAL ANALYSIS ******

  FUNCTION W1FNMN(X,NF1,NN1,NM1)
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    INTEGER,INTENT(IN):: NF1,NN1,NM1
    REAL(rkind):: W1FNMN
    REAL(rkind),PARAMETER:: SP=2.506628275D0
    INTEGER:: ILST
    REAL(rkind):: H0,EPS,CS,ES

    NF=NF1
    NN=NN1
    NM=NM1
    G1=X

    IF(ABS(X).LT.1.D-6) THEN
       IF(NF.EQ.0) THEN
          W1FNMN=16.D0/(3.D0*PI*SP)
       ELSEIF(NF.EQ.1) THEN
          IF(NM.EQ.1) W1FNMN=2.D0/(PI*SP*(1-4*NN*NN))
          IF(NM.EQ.3) W1FNMN=3.D0/(2.D0*PI*SP) &
                            *(1.D0/(1-4*NN*NN) &
                             -1.D0/(9-4*NN*NN))
          IF(NM.EQ.5) W1FNMN=5.D0/(8.D0*PI*SP) &
                            *(1.D0/(25-4*NN*NN) &
                             -3.D0/( 9-4*NN*NN) &
                             +2.D0/( 1-4*NN*NN))
       ELSEIF(NF.EQ.2) THEN
          IF(NM.EQ.0) W1FNMN=4.D0*NN/(PI*SP*(4.D0*NN*NN-1.D0))
          IF(NM.EQ.2) W1FNMN=1.D0/(2.D0*PI*SP) &
                            *( 4.D0*NN   /(4*NN*NN-1) &
                             -(2.D0*NN+1)/((2*NN+1)**2-4) &
                             -(2.D0*NN-1)/((2*NN-1)**2-4))
          IF(NM.EQ.4) W1FNMN=1.D0/(8.D0*PI*SP) &
                            *((2.D0*NN+1)/((2*NN+1)**2-16) &
                             +(2.D0*NN-1)/((2*NN-1)**2-16) &
                             -4*(2.D0*NN+1)/((2*NN+1)**2-4) &
                             -4*(2.D0*NN-1)/((2*NN-1)**2-4) &
                             +3.D0/(2*NN+1) &
                             +3.D0/(2*NN-1))
       ELSEIF(NF.EQ.3) THEN
          W1FNMN=0.D0
          IF(NM.EQ.1) THEN
             IF(NN.EQ.0) W1FNMN= 0.5D0
             IF(NN.EQ.1) W1FNMN=-0.25D0
          ELSEIF(NM.EQ.3) THEN
             IF(NN.EQ.0) W1FNMN= 0.375D0
             IF(NN.EQ.1) W1FNMN=-0.25D0
             IF(NN.EQ.2) W1FNMN= 0.0625D0
          ELSEIF(NM.EQ.5) THEN
             IF(NN.EQ.0) W1FNMN= 5.D0/16.D0
             IF(NN.EQ.1) W1FNMN=-15.D0/64.D0
             IF(NN.EQ.2) W1FNMN= 3.D0/32.D0
             IF(NN.EQ.3) W1FNMN=-1.D0/64.D0
          ENDIF
       ELSEIF(NF.EQ.4) THEN
          W1FNMN=0.D0
          IF(NM.EQ.0) THEN
             IF(NN.EQ.1) W1FNMN= 0.25D0
          ELSEIF(NM.EQ.2) THEN
             IF(NN.EQ.1) W1FNMN= 0.125D0
             IF(NN.EQ.2) W1FNMN=-0.0625D0
          ELSEIF(NM.EQ.4) THEN
             IF(NN.EQ.1) W1FNMN= 5.D0/64.D0
             IF(NN.EQ.2) W1FNMN=-1.D0/16.D0
             IF(NN.EQ.3) W1FNMN= 1.D0/64.D0
          ENDIF
       ENDIF
    ELSE
       H0=0.5
       EPS=1.D-6
       ILST=0
       CALL DEFT(CS,ES,H0,EPS,ILST,W1FNFG,'w1fnfg')
       W1FNMN=CS
    ENDIF
    RETURN
  END FUNCTION W1FNMN
END MODULE w1fflr

