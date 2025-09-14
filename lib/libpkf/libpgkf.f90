MODULE libpgkf

  USE task_kinds,ONLY: rkind
  USE task_constants,ONLY: PI
  IMPLICIT NONE

  PRIVATE
  PUBLIC pgkf
  
  REAL(rkind):: G1
  INTEGER:: NF,NN,NM

CONTAINS

!     ****** SLAVE FUNTION FOR DE INTEGRATION ******

  FUNCTION pgkf_slave(X,XM,XP)
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,XM,XP
    REAL(rkind):: pgkf_slave
    REAL(rkind):: SP=2.506628275D0
    REAL(rkind):: DUMMY,T,S,XX

    DUMMY=X
    DUMMY=XM
    T=0.5D0*PI*XP
    S=SIN(T)
    XX=G1/(2.82843D0*S)
    IF(ABS(XX).LT.10.D0) THEN
       IF(NF.EQ.0) THEN
          pgkf_slave = 2.D0*S*S*S*(-XX*ERFC(XX)+EXP(-XX*XX)/SP)
       ELSEIF(NF.EQ.1) THEN
          pgkf_slave = 0.5D0*S**NM*EXP(-XX*XX)*COS(2*NN*T)/SP
       ELSEIF(NF.EQ.2) THEN
          pgkf_slave = 0.5D0*COS(T)*S**NM*EXP(-XX*XX)*SIN(2*NN*T)/SP
       ELSEIF(NF.EQ.3) THEN
          pgkf_slave = 0.5D0*S**(NM+1)*ERFC(XX)*COS(2*NN*T)
       ELSEIF(NF.EQ.4) THEN
          pgkf_slave = 0.5D0*COS(T)*S**(NM+1)*ERFC(XX)*SIN(2*NN*T)
       ELSE
          WRITE(6,'(A,I6)') 'XX pgkf_slave: unknown NF: NF=',NF
          STOP
       ENDIF
    ELSE
       pgkf_slave = 0.D0
    ENDIF
  END FUNCTION pgkf_slave
  
  !     ****** plasma gyro kernel function ******
  !         --- NF=[0:4] ---
  !         --- NM=1,3,5 for NF=1
  !         --- NM=0,2,4 for NF=2
  !         --- NM=1,3,5 for NF=3
  !         --- NM=0,2,4 for NF=4
  

  FUNCTION pgkf(X,NF1,NN1,NM1)
    USE libde
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    INTEGER,INTENT(IN):: NF1,NN1,NM1
    REAL(rkind):: pgkf
    REAL(rkind),PARAMETER:: SP=2.506628275D0
    INTEGER:: ILST
    REAL(rkind):: H0,EPS,CS,ES

    NF=NF1
    NN=NN1
    NM=NM1
    G1=X

    IF(ABS(X).LT.1.D-6) THEN
       pgkf=0.D0
       IF(NF.EQ.0) THEN
          pgkf=16.D0/(3.D0*PI*SP)
       ELSEIF(NF.EQ.1) THEN
          IF(NM.EQ.1) THEN
             pgkf=2.D0/(PI*SP*(1-4*NN*NN))
          ELSE IF(NM.EQ.3) THEN
             pgkf=3.D0/(2.D0*PI*SP) &
                  *(1.D0/(1-4*NN*NN) &
                  -1.D0/(9-4*NN*NN))
          ELSE IF(NM.EQ.5) THEN
             pgkf=5.D0/(8.D0*PI*SP) &
                  *(1.D0/(25-4*NN*NN) &
                  -3.D0/( 9-4*NN*NN) &
                  +2.D0/( 1-4*NN*NN))
          ENDIF
       ELSEIF(NF.EQ.2) THEN
          IF(NM.EQ.0) THEN
             pgkf=4.D0*NN/(PI*SP*(4.D0*NN*NN-1.D0))
          ELSE IF(NM.EQ.2) THEN
             pgkf=1.D0/(2.D0*PI*SP) &
                            *( 4.D0*NN   /(4*NN*NN-1) &
                             -(2.D0*NN+1)/((2*NN+1)**2-4) &
                             -(2.D0*NN-1)/((2*NN-1)**2-4))
          ELSE IF(NM.EQ.4) THEN
             pgkf=1.D0/(8.D0*PI*SP) &
                            *((2.D0*NN+1)/((2*NN+1)**2-16) &
                             +(2.D0*NN-1)/((2*NN-1)**2-16) &
                             -4*(2.D0*NN+1)/((2*NN+1)**2-4) &
                             -4*(2.D0*NN-1)/((2*NN-1)**2-4) &
                             +3.D0/(2*NN+1) &
                             +3.D0/(2*NN-1))
          ENDIF
       ELSEIF(NF.EQ.3) THEN
          pgkf=0.D0
          IF(NM.EQ.1) THEN
             IF(NN.EQ.0) pgkf= 0.5D0
             IF(NN.EQ.1) pgkf=-0.25D0
          ELSEIF(NM.EQ.3) THEN
             IF(NN.EQ.0) pgkf= 0.375D0
             IF(NN.EQ.1) pgkf=-0.25D0
             IF(NN.EQ.2) pgkf= 0.0625D0
          ELSEIF(NM.EQ.5) THEN
             IF(NN.EQ.0) pgkf= 5.D0/16.D0
             IF(NN.EQ.1) pgkf=-15.D0/64.D0
             IF(NN.EQ.2) pgkf= 3.D0/32.D0
             IF(NN.EQ.3) pgkf=-1.D0/64.D0
          ENDIF
       ELSEIF(NF.EQ.4) THEN
          pgkf=0.D0
          IF(NM.EQ.0) THEN
             IF(NN.EQ.1) pgkf= 0.25D0
          ELSEIF(NM.EQ.2) THEN
             IF(NN.EQ.1) pgkf= 0.125D0
             IF(NN.EQ.2) pgkf=-0.0625D0
          ELSEIF(NM.EQ.4) THEN
             IF(NN.EQ.1) pgkf= 5.D0/64.D0
             IF(NN.EQ.2) pgkf=-1.D0/16.D0
             IF(NN.EQ.3) pgkf= 1.D0/64.D0
          ENDIF
       ELSE
          WRITE(6,'(A,I6)') 'XX pgkf: out of range : NF=',NF
          STOP
       ENDIF
    ELSE
       H0=0.5D0
       EPS=1.D-6
       ILST=0
       CALL DEFT(CS,ES,H0,EPS,ILST,pgkf_slave,'pgkf')
       pgkf=CS
    ENDIF
    WRITE(6,'(A,ES12.4,3I6,ES12.4)') 'pgkf:',X,NF1,NN1,NM1,pgkf
    RETURN
  END FUNCTION pgkf
END MODULE libpgkf

