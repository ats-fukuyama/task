C     $Id$
C
C     ***** TASK/TR/GOUT MENU *****
C
      SUBROUTINE TRGOUT
C
      INCLUDE 'trcomm.inc'
C
      CHARACTER KIG*5,K1*1,K2*1,K3*1,K4*1,K5*1,KK*3
      SAVE INIT,INQG
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         INQG=0
         INIT=1
      ENDIF
C
  1   CONTINUE
         WRITE(6,*) '# SELECT : R1-R9, T1-T9, G1-G7, P1-P5, Z1, Y1,',
     &                        ' A1-A2, E1-E9, D1-D59, M1-M5'
         WRITE(6,*) '           N1-N6, S/SAVE  L/LOAD  H/HELP  ',
     &              'C/CLEAR  I/INQ  X/EXIT'
         READ(5,'(A5)',ERR=1,END=9000) KIG
         K1=KIG(1:1)
         K2=KIG(2:2)
         K3=KIG(3:3)
         K4=KIG(4:4)
         K5=KIG(5:5)
         CALL GUCPTL(K1)
         CALL GUCPTL(K2)
         CALL GUCPTL(K3)
         CALL GUCPTL(K4)
         CALL GUCPTL(K5)
         KK=K3//K4//K5
C
         IF(K1.EQ.'C') THEN
            NGR=0
            NGT=0
         ELSEIF(K1.EQ.'I') THEN
            IF(INQG.EQ.0) THEN
               INQG=4      
               WRITE(6,*) '## GRAPHIC SCALE INQUIRE MODE : ON'
            ELSE
               INQG=0      
               WRITE(6,*) '## GRAPHIC SCALE INQUIRE MODE : OFF'
            ENDIF
         ELSEIF(K1.EQ.'S') THEN
            CALL TRGRSV
         ELSEIF(K1.EQ.'L') THEN
            CALL TRGRLD
         ELSEIF(K1.EQ.'R') THEN
            CALL TRGRR0(K2,INQG)
         ELSEIF(K1.EQ.'Y') THEN
            CALL TRGRY0(K2,INQG)
         ELSEIF(K1.EQ.'T') THEN
            CALL TRGRT0(K2,INQG)
         ELSEIF(K1.EQ.'Z') THEN
            CALL TRGRX0(K2,INQG)
         ELSEIF(K1.EQ.'G') THEN
            CALL TRGRG0(K2,INQG)
         ELSEIF(K1.EQ.'P') THEN
            CALL TRGRP0(K2,INQG)
         ELSEIF(K1.EQ.'A') THEN
            CALL TRGRA0(K2,INQG)
         ELSEIF(K1.EQ.'E') THEN
            CALL TRGRE0(K2,INQG)
         ELSEIF(K1.EQ.'D') THEN
            CALL TRGRD0(K2,KK,INQG)
         ELSEIF(K1.EQ.'H') THEN
            CALL TRHELP('G')
         ELSEIF(K1.EQ.'M') THEN
            CALL TRCOMP(K2,INQG)
         ELSEIF(K1.EQ.'N') THEN
            CALL TRGRN0(K2,INQG)
         ELSEIF(K1.EQ.'X') THEN
            GOTO 9000
         ELSE
            WRITE(6,*) 'UNSUPPORTED GRAPH ID'
         ENDIF
      GOTO 1
C
 9000 RETURN
      END
