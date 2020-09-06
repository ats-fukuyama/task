! libcharx.f90

MODULE libcharx

  PRIVATE
  PUBLIC ksplitn

CONTAINS

!     ***** SEPARATE STRING AT CHAR *****

  SUBROUTINE ksplitn(line,kidn,lword,nword_max,kworda)

    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)  :: line      ! target string
    CHARACTER(LEN=*),INTENT(IN)  :: kidn      ! list of separating chars
    INTEGER,INTENT(IN)  :: lword              ! length of a word
    INTEGER,INTENT(OUT):: nword_max          ! number of words
    CHARACTER(LEN=lword),INTENT(OUT),ALLOCATABLE:: &
         kworda(:) ! array of separated words
    CHARACTER(LEN=1),ALLOCATABLE:: kida(:)    ! array of separating chars
    INTEGER:: nkid ! number of separating chars
    INTEGER,allocatable:: isepa(:) ! position of separated char
    INTEGER:: i,n,nword,id
    CHARACTER(LEN=1):: kn

    nkid=LEN_TRIM(kidn)
    IF(ALLOCATED(kida)) DEALLOCATE(kida)
    ALLOCATE(kida(nkid))
    DO i=1,nkid
       kida(i)=kidn(i:i)
    END DO

    ! --- calculate number of separating chars ---
    
    nword=0
    DO n=1,LEN_TRIM(line)
       kn=line(n:n)
       id=0
       DO i=1,nkid
          IF(kn.EQ.kidn(i:i)) THEN
             id=1
             EXIT
          END IF
       END DO
       IF(id.EQ.1) nword=nword+1
    END DO
    nword_max=nword+1  ! number of words

    ! --- allocate kworda and isepa ---

    IF(ALLOCATED(kworda)) DEALLOCATE(kworda,isepa)
    ALLOCATE(kworda(nword_max),isepa(nword_max))

    ! --- make a list of separator position ---

    nword=0
    DO n=1,LEN_TRIM(line)
       kn=line(n:n)
       id=0
       DO i=1,nkid
          IF(kn.EQ.kida(i)) THEN
             id=1
             EXIT
          END IF
       END DO
       IF(id.EQ.1) THEN
          nword=nword+1
          isepa(nword)=n
       END IF
    END DO
    nword=nword+1
    isepa(nword)=LEN_TRIM(line)+1

    ! --- make a list of separated words ---

    n=1
    DO nword=1,nword_max
       kworda(nword)=line(n:isepa(nword)-1)
       n=isepa(nword)+1
    END DO
    RETURN
  END SUBROUTINE ksplitn
END MODULE libcharx
