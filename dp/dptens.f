C     $Id$
C
C     ****** CALCULATE DIELECTRIC TENSOR ******
C
      SUBROUTINE DPTENS(CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION CLDISP(6),CLDISP1(6)
C
      IF(RN(NS).LE.0.D0) THEN
         CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP)
      ELSEIF(RHON_LOC.LT.RHON_MIN.OR.RHON_LOC.GT.RHON_MAX) THEN
         ID1=MOD(MODELP(NS),100)
         CALL DPTENS_AN(ID1,CW,CKPR,CKPP,NS,CLDISP)
      ELSE
         ID1=MOD(MODELP(NS),100)
         ID2=MODELP(NS)/100
         IDV=MODELV(NS)
         SELECT CASE(IDV)
         CASE(0)
            CALL DPTENS_AN(ID1,CW,CKPR,CKPP,NS,CLDISP)
         CASE(1)
            CALL DPFMFL(NS,0)
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DPHOTFI(CW,CKPR,CKPP,NS,CLDISP)
            ELSE
               CALL DPHOTF(CW,CKPR,CKPP,NS,CLDISP)
            ENDIF
         CASE(2)
            CALL DPFPFL(NS)
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DPHOTFI(CW,CKPR,CKPP,NS,CLDISP)
            ELSE
               CALL DPHOTF(CW,CKPR,CKPP,NS,CLDISP)
            ENDIF
         CASE(3)
            CALL DPFMFL(NS,1)
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DPHOTRI(CW,CKPR,CKPP,NS,CLDISP)
            ELSE
               CALL DPHOTR(CW,CKPR,CKPP,NS,CLDISP)
            ENDIF
         CASE(4)
            CALL DPFPFL(NS)
            IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
               CALL DPHOTRI(CW,CKPR,CKPP,NS,CLDISP)
            ELSE
               CALL DPHOTR(CW,CKPR,CKPP,NS,CLDISP)
            ENDIF
         END SELECT
C
         IF(ID2.EQ.2) THEN
            CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP1)
         ELSEIF(ID2.EQ.3) THEN
            CALL DPTNKP(CW,CKPR,CKPP,NS,CLDISP1)
         ENDIF
C
         IF(ID2.EQ.2.OR.ID2.EQ.3) THEN
            CLDISP(1)=DBLE(CLDISP1(1))+CI*DIMAG(CLDISP(1))
            CLDISP(2)=DBLE(CLDISP1(2))+CI*DIMAG(CLDISP(2))
            CLDISP(3)=DBLE(CLDISP1(3))+CI*DIMAG(CLDISP(3))
            CLDISP(4)=DBLE(CLDISP1(4))+CI*DIMAG(CLDISP(4))
            CLDISP(5)=CI*DIMAG(CLDISP1(5))+DBLE(CLDISP(5))
            CLDISP(6)=CI*DIMAG(CLDISP1(6))+DBLE(CLDISP(6))
         ENDIF
      ENDIF
C         WRITE(6,'(A,3I5,1PE12.4)') 
C     &        'NS,ID1,ID2,RN =',NS,ID1,ID2,RN(NS)
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CW,R,P=',CW,CKPR,CKRR
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CLDISP=',CLDISP(1),CLDISP(2),CLDISP(3)
C         WRITE(6,'(A,1P6E12.4)') 
C     &        'CLDISP=',CLDISP(4),CLDISP(5),CLDISP(6)
      RETURN
      END
C
C     ****** CALCULATE ANALYTIC DIELECTRIC TENSOR ******
C
      SUBROUTINE DPTENS_AN(ID1,CW,CKPR,CKPP,NS,CLDISP)
C
      INCLUDE '../dp/dpcomm.inc'
      INCLUDE '../pl/plcom2.inc'
C
      DIMENSION CLDISP(6)
C
      SELECT CASE(ID1)
      CASE(0) ! collisionless cold model
         CALL DPTNCL(CW,CKPR,CKPP,NS,CLDISP)
      CASE(1) ! collisional cold model
         CALL DPTNCC(CW,CKPR,CKPP,NS,CLDISP)
      CASE(2) ! idial MHD model
         CALL DPTNIM(CW,CKPR,CKPP,NS,CLDISP)
      CASE(3) ! resistive MHD model
         CALL DPTNRM(CW,CKPR,CKPP,NS,CLDISP)
      CASE(4) ! kinetic without FLR model
         CALL DPTNHP(CW,CKPR,CKPP,NS,CLDISP)
      CASE(5) ! kinetic with FLR model
         CALL DPTNKP(CW,CKPR,CKPP,NS,CLDISP)
      CASE(6) ! weakly relativistic model
         CALL DPTNKR(CW,CKPR,CKPP,NS,CLDISP)
      CASE(11:15) ! old WM models model
         CALL DPTNFK2(CW,CKPR,CKPP,NS,CLDISP)
      CASE(16) ! old WM drift kinetic model
         CALL DPTNDK0(CW,CKPR,CKPP,NS,CLDISP)
      CASE DEFAULT
         WRITE(6,*) 'XX WRONG MODELP IN DPTENS: ID1=',ID1
      END SELECT
      RETURN
      END
