!
!   Data library based on 
!       Atomic data and Nuclear data Tables, 20 (1977) 397-439
!       by D. E. Post, R. V. Jensen, C. B. Tarter, W. H. Grasberger, 
!          and W. A. Lokke
!
MODULE ADPOST

  USE bpsd
  USE libfio

  PRIVATE
  PUBLIC read_adpost,broadcast_adpost,func_adpost,nd_npa_adpost

  INTEGER,PARAMETER:: IZDIMD = 140   ! maximum number of NPAs

  INTEGER:: NDMAX=0                      ! Maximum number of atomic data
  INTEGER:: LMAXM                        ! Maximum of number of data blocks
  CHARACTER(LEN=8),ALLOCATABLE:: KZ0A(:) ! Element name of data ND
  INTEGER,ALLOCATABLE:: NPAA(:)          ! Atomic Number of data ND
  INTEGER,ALLOCATABLE:: IM0A(:)          ! Mass Number of data ND
  INTEGER,ALLOCATABLE:: LMAXA(:)         ! Number of data blocks of data ND
  REAL(rkind),ALLOCATABLE:: TMINA(:,:)   ! Minimum Te of a block L [keV]
  REAL(rkind),ALLOCATABLE:: TMAXA(:,:)   ! Minimum Te of a block L [keV]
  REAL(rkind),ALLOCATABLE:: TMINB(:,:)   ! Minimum Te of a block L [keV]
  REAL(rkind),ALLOCATABLE:: TMAXB(:,:)   ! Maximum Te of a block L [keV]
  REAL(rkind),ALLOCATABLE:: TMINC(:,:)   ! Minimum Te of a block L [keV]
  REAL(rkind),ALLOCATABLE:: TMAXC(:,:)   ! Maximum Te of a block L [keV]
  REAL(rkind),ALLOCATABLE:: COEFA(:,:,:) ! Cooling rate coef of ..
  REAL(rkind),ALLOCATABLE:: COEFB(:,:,:) ! Zav coef ..
  REAL(rkind),ALLOCATABLE:: COEFC(:,:,:) ! Z2av coef of (log_10 Te [keV])^N
  INTEGER:: ND_TABLE(IZDIMD)             ! Table of ND for NPA
  INTEGER,ALLOCATABLE:: NPA_TABLE(:)     ! Table of available NPA 

CONTAINS

  SUBROUTINE read_adpost(adpost_filename,IERR)

    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: adpost_filename
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN,NPA,ND,NL,N

!   --- deallocate data array, if allocatd ---

    IF(ALLOCATED(KZ0A)) DEALLOCATE(KZ0A)
    IF(ALLOCATED(NPAA)) DEALLOCATE(NPAA)
    IF(ALLOCATED(IM0A)) DEALLOCATE(IM0A)
    IF(ALLOCATED(LMAXA)) DEALLOCATE(LMAXA)
    IF(ALLOCATED(TMINA)) DEALLOCATE(TMINA)
    IF(ALLOCATED(TMINB)) DEALLOCATE(TMINB)
    IF(ALLOCATED(TMINC)) DEALLOCATE(TMINC)
    IF(ALLOCATED(TMAXA)) DEALLOCATE(TMAXA)
    IF(ALLOCATED(TMAXB)) DEALLOCATE(TMAXB)
    IF(ALLOCATED(TMAXC)) DEALLOCATE(TMAXC)
    IF(ALLOCATED(COEFA)) DEALLOCATE(COEFA)
    IF(ALLOCATED(COEFB)) DEALLOCATE(COEFB)
    IF(ALLOCATED(COEFC)) DEALLOCATE(COEFC)
    IF(ALLOCATED(NPA_TABLE)) DEALLOCATE(NPA_TABLE)

!   --- initialize conversion table between NPA and ND ---

    DO NPA=1,IZDIMD
       ND_TABLE(NPA)=0
    END DO
    LUN=12

    CALL FROPEN(LUN,adpost_filename,1,0,'ADPOST',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX READ_ADPOST: FROPEN: IERR =',IERR
       IERR=1
       RETURN
    END IF

!   --- start reading data file ---

    READ(LUN,*,ERR=9001,END=9011) NDMAX
    ALLOCATE(KZ0A(NDMAX),NPAA(NDMAX),IM0A(NDMAX),LMAXA(NDMAX))
    READ(LUN,*,ERR=9002,END=9012) LMAXM
    ALLOCATE(TMINA(LMAXM,NDMAX),TMAXA(LMAXM,NDMAX))
    ALLOCATE(TMINB(LMAXM,NDMAX),TMAXB(LMAXM,NDMAX))
    ALLOCATE(TMINC(LMAXM,NDMAX),TMAXC(LMAXM,NDMAX))
    ALLOCATE(COEFA(0:5,LMAXM,NDMAX))
    ALLOCATE(COEFB(0:5,LMAXM,NDMAX))
    ALLOCATE(COEFC(0:5,LMAXM,NDMAX))
    ALLOCATE(NPA_TABLE(NDMAX))

    DO ND=1,NDMAX
       READ(LUN,*,ERR=9003,END=9013) KZ0A(ND)
       READ(LUN,*,ERR=9004,END=9014) NPAA(ND)
       READ(LUN,*,ERR=9005,END=9015) IM0A(ND)
       READ(LUN,*,ERR=9006,END=9016) LMAXA(ND)
       DO NL=1,LMAXA(ND)
          READ(LUN,*,ERR=9007,END=9017) &
               TMINA(NL,ND),TMAXA(NL,ND),(COEFA(N,NL,ND),N=0,5)
       END DO
       DO NL=1,LMAXA(ND)
          READ(LUN,*,ERR=9008,END=9018) &
               TMINB(NL,ND),TMAXB(NL,ND),(COEFB(N,NL,ND),N=0,5)
       END DO
       DO NL=1,LMAXA(ND)
          READ(LUN,*,ERR=9009,END=9019) &
               TMINC(NL,ND),TMAXC(NL,ND),(COEFC(N,NL,ND),N=0,5)
       END DO
       ND_TABLE(NPAA(ND))=ND
    END DO

    ND=0
    DO NPA=1,IZDIMD
       IF(ND_TABLE(NPA).NE.0) THEN
          ND=ND+1
          NPA_TABLE(ND)=NPA
       END IF
    END DO
    WRITE(6,'(A)') '## List of available NPA'
    WRITE(6,'(20I4)') (NPA_TABLE(ND),ND=1,NDMAX)

    REWIND LUN
    CLOSE(LUN)
    IERR=0
    RETURN

9001 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading NDMAX'
    IERR=1
    RETURN
9011 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading NDMAX'
    IERR=11
    RETURN
9002 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading LMAXM'
    IERR=2
    RETURN
9012 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading LMAXM'
    IERR=12
    RETURN
9003 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading KZ0A'
    IERR=3
    RETURN
9013 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading KZ0A'
    IERR=13
    RETURN
9004 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading NPAA'
    IERR=4
    RETURN
9014 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading NPAA'
    IERR=14
    RETURN
9005 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading IM0A'
    IERR=5
    RETURN
9015 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading IM0A'
    IERR=15
    RETURN
9006 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading LMAXA'
    IERR=6
    RETURN
9016 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading LMAXA'
    IERR=16
    RETURN
9007 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading TMINA,TMAXA,COEFA'
    IERR=7
    RETURN
9017 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading TMINA,TMAXA,COEFA'
    IERR=17
    RETURN
9008 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading TMINB,TMAXB,COEFB'
    IERR=78
    RETURN
9018 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading TMINB,TMAXB,COEFB'
    IERR=18
    RETURN
9009 WRITE(6,'(A)') 'XX READ_ADPOST: READ ERROR in reading TMINC,TMAXC,COEFC'
    IERR=9
    RETURN
9019 WRITE(6,'(A)') 'XX READ_ADPOST: END OF FILE in reading TMINC,TMAXC,COEFC'
    IERR=19
    RETURN
  END SUBROUTINE read_adpost

! --- broadcast loaded adpost data ---

  SUBROUTINE broadcast_adpost(IERR)

    USE libmpi
    USE commpi
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NPA,ND,NL,N,i
    REAl(rkind),ALLOCATABLE:: work(:)

!   --- deallocate data array, if allocatd ---

    IF(nrank.NE.0) THEN
       IF(ALLOCATED(KZ0A)) DEALLOCATE(KZ0A)
       IF(ALLOCATED(NPAA)) DEALLOCATE(NPAA)
       IF(ALLOCATED(IM0A)) DEALLOCATE(IM0A)
       IF(ALLOCATED(LMAXA)) DEALLOCATE(LMAXA)
       IF(ALLOCATED(TMINA)) DEALLOCATE(TMINA)
       IF(ALLOCATED(TMINB)) DEALLOCATE(TMINB)
       IF(ALLOCATED(TMINC)) DEALLOCATE(TMINC)
       IF(ALLOCATED(TMAXA)) DEALLOCATE(TMAXA)
       IF(ALLOCATED(TMAXB)) DEALLOCATE(TMAXB)
       IF(ALLOCATED(TMAXC)) DEALLOCATE(TMAXC)
       IF(ALLOCATED(COEFA)) DEALLOCATE(COEFA)
       IF(ALLOCATED(COEFB)) DEALLOCATE(COEFB)
       IF(ALLOCATED(COEFC)) DEALLOCATE(COEFC)
       IF(ALLOCATED(NPA_TABLE)) DEALLOCATE(NPA_TABLE)
    END IF
!   --- initialize conversion table between NPA and ND ---

    DO NPA=1,IZDIMD
       ND_TABLE(NPA)=0
    END DO

!   --- start broadcasting data ---

    CALL mtx_broadcast1_integer(NDMAX)
    CALL mtx_broadcast1_integer(LMAXM)

    IF(nrank.NE.0) ALLOCATE(KZ0A(NDMAX),NPAA(NDMAX),IM0A(NDMAX),LMAXA(NDMAX))

    DO ND=1,NDMAX
       CALL mtx_broadcast_character(KZ0A(ND),8)
    END DO
    CALL mtx_broadcast_integer(NPAA,NDMAX)
    CALL mtx_broadcast_integer(IM0A,NDMAX)
    CALL mtx_broadcast_integer(LMAXA,NDMAX)

    IF(nrank.NE.0) THEN
       ALLOCATE(TMINA(LMAXM,NDMAX),TMAXA(LMAXM,NDMAX))
       ALLOCATE(TMINB(LMAXM,NDMAX),TMAXB(LMAXM,NDMAX))
       ALLOCATE(TMINC(LMAXM,NDMAX),TMAXC(LMAXM,NDMAX))
       ALLOCATE(COEFA(0:5,LMAXM,NDMAX))
       ALLOCATE(COEFB(0:5,LMAXM,NDMAX))
       ALLOCATE(COEFC(0:5,LMAXM,NDMAX))
       ALLOCATE(NPA_TABLE(NDMAX))
    END IF

    ALLOCATE(work(6*LMAXM))
    DO ND=1,NDMAX
       IF(nrank.EQ.0) work(1:LMAXA(ND))=TMINA(1:LMAXA(ND),ND)
       CALL mtx_broadcast_real8(work,LMAXA(ND))
       IF(nrank.NE.0) TMINA(1:LMAXA(ND),ND)=work(1:LMAXA(ND))
       IF(nrank.EQ.0) work(1:LMAXA(ND))=TMAXA(1:LMAXA(ND),ND)
       CALL mtx_broadcast_real8(work,LMAXA(ND))
       IF(nrank.NE.0) TMAXA(1:LMAXA(ND),ND)=work(1:LMAXA(ND))
       IF(nrank.EQ.0) work(1:LMAXA(ND))=TMINB(1:LMAXA(ND),ND)
       CALL mtx_broadcast_real8(work,LMAXA(ND))
       IF(nrank.NE.0) TMINB(1:LMAXA(ND),ND)=work(1:LMAXA(ND))
       IF(nrank.EQ.0) work(1:LMAXA(ND))=TMAXB(1:LMAXA(ND),ND)
       CALL mtx_broadcast_real8(work,LMAXA(ND))
       IF(nrank.NE.0) TMAXB(1:LMAXA(ND),ND)=work(1:LMAXA(ND))
       IF(nrank.EQ.0) work(1:LMAXA(ND))=TMINC(1:LMAXA(ND),ND)
       CALL mtx_broadcast_real8(work,LMAXA(ND))
       IF(nrank.NE.0) TMINC(1:LMAXA(ND),ND)=work(1:LMAXA(ND))
       IF(nrank.EQ.0) work(1:LMAXA(ND))=TMAXC(1:LMAXA(ND),ND)
       CALL mtx_broadcast_real8(work,LMAXA(ND))
       IF(nrank.NE.0) TMAXC(1:LMAXA(ND),ND)=work(1:LMAXA(ND))

       IF(nrank.EQ.0) THEN
          i=0
          DO NL=1,LMAXA(ND)
             DO N=0,5
                i=i+1
                work(i)=COEFA(N,NL,ND)
             END DO
          END DO
       END IF

       CALL mtx_broadcast_real8(work,6*LMAXA(ND))

       IF(nrank.NE.0) THEN
          i=0
          DO NL=1,LMAXA(ND)
             DO N=0,5
                i=i+1
                COEFA(N,NL,ND)=work(i)
             END DO
          END DO
       ELSE
          i=0
          DO NL=1,LMAXA(ND)
             DO N=0,5
                i=i+1
                work(i)=COEFB(N,NL,ND)
             END DO
          END DO
       END IF

       CALL mtx_broadcast_real8(work,6*LMAXA(ND))

       IF(nrank.ne.0) THEN
          i=0
          DO NL=1,LMAXA(ND)
             DO N=0,5
                i=i+1
                COEFB(N,NL,ND)=work(i)
             END DO
          END DO
       ELSE
          i=0
          DO NL=1,LMAXA(ND)
             DO N=0,5
                i=i+1
                work(i)=COEFC(N,NL,ND)
             END DO
          END DO
       END IF

       CALL mtx_broadcast_real8(work,6*LMAXA(ND))

       IF(nrank.ne.0) THEN
          i=0
          DO NL=1,LMAXA(ND)
             DO N=0,5
                i=i+1
                COEFC(N,NL,ND)=work(i)
             END DO
          END DO
       END IF
    END DO
    DEALLOCATE(work)

    DO ND=1,NDMAX
       ND_TABLE(NPAA(ND))=ND
    END DO

    ND=0
    DO NPA=1,IZDIMD
       IF(ND_TABLE(NPA).NE.0) THEN
          ND=ND+1
          NPA_TABLE(ND)=NPA
       END IF
    END DO

    IERR=0
    RETURN

  END SUBROUTINE broadcast_adpost

!--- fucntion to calculate Z_ave, Z^2_ave, and D_rad.

  FUNCTION func_adpost(NPA,ID,PT)

    IMPLICIT NONE
    REAL(rkind):: func_adpost
    INTEGER,INTENT(IN):: NPA    ! Atomic number
    INTEGER,INTENT(IN):: ID     ! Data type 1=Z_ave,2=Z^2_ave, 3=I_rad
    REAL(rkind),INTENT(IN):: PT ! Electron temperature [keV]
    REAL(rkinD):: PTL,SUMX
    INTEGER:: ND,N,L,LMAX,NL

    ND=ND_TABLE(NPA)

    IF(ND.EQ.0) THEN
       WRITE(6,'(A,I4)') 'XX func_adpost: Undefined NPA: NPA =',NPA
       func_adpost=0.D0
       RETURN
    END IF

    IF(PT.LE.0.D0) THEN
       WRITE(6,'(A,1PE12.4)') &
            'XX func_adpost: Negatic electron temperature: PT =',PT
       func_adpost=0.D0
       RETURN
    END IF

    PTL=PT
    LMAX=LMAXA(ND)
    SELECT CASE(ID)
    CASE(3) ! Radiated power
       IF(PTL.LT.TMINA(1,ND)) THEN
          PTL=TMINA(1,ND)
          L=1
       ELSE IF(PTL.GT.TMAXA(LMAX,ND)) THEN
          PTL=TMAXA(LMAX,ND)
          L=LMAX
       ELSE
          DO NL=1,LMAX
             L=NL
             IF(PTL.LE.TMAXA(NL,ND)) EXIT
          END DO
       END IF
       SUMX=0.D0
       PTL=LOG10(PTL)
       DO N=5,0,-1
          SUMX=COEFA(N,L,ND)+PTL*SUMX
       END DO
       func_adpost=SUMX
    CASE(1)  ! Z
       IF(PTL.LT.TMINB(1,ND)) THEN
          PTL=TMINB(1,ND)
          L=1
       ELSE IF(PTL.GT.TMAXB(LMAX,ND)) THEN
          PTL=TMAXB(LMAX,ND)
          L=LMAX
       ELSE
          DO NL=1,LMAX
             L=NL
             IF(PTL.LE.TMAXB(NL,ND)) EXIT
          END DO
       END IF
       SUMX=0.D0
       PTL=LOG10(PTL)
       DO N=5,0,-1
          SUMX=COEFB(N,L,ND)+PTL*SUMX
       END DO
       func_adpost=SUMX
    CASE(2) ! Z^2
       IF(PTL.LT.TMINC(1,ND)) THEN
          PTL=TMINC(1,ND)
          L=1
       ELSE IF(PTL.GT.TMAXC(LMAX,ND)) THEN
          PTL=TMAXC(LMAX,ND)
          L=LMAX
       ELSE
          DO NL=1,LMAX
             L=NL
             IF(PTL.LE.TMAXC(NL,ND)) EXIT
          END DO
       END IF
       SUMX=0.D0
       PTL=LOG10(PTL)
       DO N=5,0,-1
          SUMX=COEFC(N,L,ND)+PTL*SUMX
       END DO
       func_adpost=SUMX
    CASE DEFAULT
       func_adpost=0.D0
    END SELECT
    RETURN
  END FUNCTION func_adpost

  FUNCTION nd_npa_adpost(NPA)
    IMPLICIT NONE
    INTEGER:: nd_npa_adpost
    INTEGER,INTENT(IN):: NPA

    IF(NDMAX.EQ.0) THEN
       WRITE(6,'(A)') 'XX ADPOST data is empty'
       nd_npa_adpost=0
    ELSE
       nd_npa_adpost=ND_TABLE(NPA)
    END IF
    RETURN
  END FUNCTION ND_NPA_ADPOST
    

END MODULE ADPOST

