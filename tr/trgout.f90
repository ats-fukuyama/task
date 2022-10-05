! trgout.f

MODULE trgout

  PRIVATE
  PUBLIC tr_gout

CONTAINS

!     ***** TASK/TR/GOUT MENU *****

      SUBROUTINE tr_gout

      USE TRCOMM, ONLY : NGR, NGT
      USE trfile
      USE libchar
      IMPLICIT NONE
      INTEGER, SAVE :: INIT =0, INQG
      CHARACTER(LEN=5) :: KIG
      CHARACTER(LEN=3) :: KK
      CHARACTER(LEN=1) :: K1, K2, K3, K4, K5
      integer       :: ist

      IF(INIT.EQ.0) THEN
         INQG=0
         INIT=1
      ENDIF

      DO
         WRITE(6,*) '# SELECT : R1-R9, T1,2,5-8, G1-G7, Z1, Y1,',' A1-A2, E1-E6, EQ, D1-D67, M1-M5'
         WRITE(6,*) '           N1-N2, S/SAVE  L/LOAD  H/HELP  ','C/CLEAR  I/INQ  X/EXIT'
         READ(5,'(A5)',iostat=ist) KIG
         if(ist > 0) then
            cycle
         elseif(ist < 0) then
            exit
         end if
         K1=KIG(1:1)
         K2=KIG(2:2)
         K3=KIG(3:3)
         K4=KIG(4:4)
         K5=KIG(5:5)
         CALL toupper(K1)
         CALL toupper(K2)
         CALL toupper(K3)
         CALL toupper(K4)
         CALL toupper(K5)
         KK=K3//K4//K5

         select case(K1)
         case('C')
            NGR=0
            NGT=0
         case('I')
            IF(INQG.EQ.0) THEN
               INQG=4
               WRITE(6,*) '## GRAPHIC SCALE INQUIRE MODE : ON'
            ELSE
               INQG=0
               WRITE(6,*) '## GRAPHIC SCALE INQUIRE MODE : OFF'
            ENDIF
         case('S')
            CALL tr_grsv
         case('L')
            CALL tr_grld
         case('R')
            CALL TRGRR0(K2,INQG)
         case('Y')
            CALL TRGRY0(K2,INQG)
         case('T')
            CALL TRGRT0(K2,INQG)
         case('Z')
            CALL TRGRX0(K2,INQG)
         case('G')
            CALL TRGRG0(K2,INQG)
         case('A')
            CALL TRGRA0(K2,INQG)
         case('E')
            CALL TRGRE0(K2,INQG)
         case('D')
            CALL TRGRD0(K2,KK,INQG)
         case('H')
            CALL TRHELP('G')
         case('M')
            CALL TRCOMP(K2,INQG)
         case('N')
            CALL TRGRN0(K2,INQG)
         case('X')
            GOTO 9000
         case default
            WRITE(6,*) 'UNSUPPORTED GRAPH ID'
         end select
      END DO

 9000 RETURN
      END SUBROUTINE tr_gout
    END MODULE trgout
