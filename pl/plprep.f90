! plprep.f90

MODULE plprep

  PRIVATE
  PUBLIC pl_prep

CONTAINS

  ! --- Initialize pl ---'
  
  SUBROUTINE pl_prep
    IMPLICIT NONE

  ! --- Initialize ID_NS and KID_NS ---'
    CALL pl_prep_ns

  END SUBROUTINE pl_prep
  
  ! --- Initialize ID_NS and KID_NS ---'
  
  SUBROUTINE pl_prep_ns

    USE plcomm
    IMPLICIT NONE
    INTEGER:: ns,npm,id_stop

    ID_STOP=0
    DO NS=1,NSMAX
       IF(NPA(NS).EQ.0) THEN  ! electron
          KID_NS(NS)='e   '
          ID_NS(NS)=-1
          PA(NS)=AME/AMP
          PZ(NS)=-1.D0
       ELSE                   ! ion
          ID_NS(NS)=1
          NPM=NINT(PA(NS))
          SELECT CASE(NPA(NS))
          CASE(1)
             SELECT CASE(NPM)
             CASE(1)
                KID_NS(NS)='H   '
             CASE(2)
                KID_NS(NS)='D   '
             CASE(3)
                KID_NS(NS)='T   '
             CASE DEFAULT
                ID_NS(NS)=0
             END SELECT
          CASE(2)
             SELECT CASE(NPM)
             CASE(3)
                KID_NS(NS)='He3 '
             CASE(4)
                KID_NS(NS)='He4 '
             CASE DEFAULT
                ID_NS(NS)=0
             END SELECT
          CASE(6)
             SELECT CASE(NPM)
             CASE(12)
                KID_NS(NS)='C   '
             CASE DEFAULT
                ID_NS(NS)=0
             END SELECT
          CASE(12)
             SELECT CASE(NPM)
             CASE(12)
                KID_NS(NS)='C   '
             CASE DEFAULT
                ID_NS(NS)=0
             END SELECT
          CASE(26)
             SELECT CASE(NPM)
             CASE(50)
                KID_NS(NS)='Fe  '
             CASE DEFAULT
                ID_NS(NS)=0
             END SELECT
          END SELECT
          IF(ID_NS(NS).EQ.0) THEN
             WRITE(6,*) 'XX tr_prep_ns: undefined NPA and PM FOR NS !!'
             WRITE(6,*) 'XX      NS,NPA,PM=',NS,NPA(NS),PA(NS)
             ID_STOP=1
             STOP
          END IF
          IF(PZ(NS).EQ.0) ID_NS(NS)=0 ! neutral
       END IF
    END DO
    IF(ID_STOP.EQ.1) STOP
  END SUBROUTINE pl_prep_ns
END MODULE plprep

