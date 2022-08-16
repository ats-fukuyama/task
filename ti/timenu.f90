! timenu.f90

! *** ti menu ***

MODULE timenu

  PRIVATE
  PUBLIC ti_menu

CONTAINS

  SUBROUTINE ti_menu

      USE bpsd
      USE ticomm
      USE tiinit
      USE tiparm
      USE tiprep
      USE tiexec
      USE tigout
      USE plload,ONLY: pl_load_trdata
      USE plprof,ONLY: pl_prf_type,pl_prof
      USE libkio
      USE libchar
      USE libmpi
      
      IMPLICIT NONE
      INTEGER(ikind) :: IERR, MODE, nid,nr,nsa,ns
      TYPE(pl_prf_type),ALLOCATABLE:: plf(:)
      REAL(rkind):: RHON
      INTEGER(ikind), SAVE :: INIT=0
      REAL:: cputime1,cputime2
      CHARACTER(LEN=1) :: KID
      CHARACTER(LEN=80):: LINE

!     ------ SELECTION OF TASK TYPE ------

      IERR=0

    1 CONTINUE
      IF(nrank.EQ.0) THEN
         IF(INIT.EQ.0) THEN
            WRITE(6,601)
601         FORMAT('## TI MENU: P,V/PARM  R/RUN  L/LOAD  ', &
                               'W/WRITE  H/HELP  Q/QUIT')
         ELSEIF(INIT.EQ.1) THEN
            WRITE(6,602)
602         FORMAT('## TI MENU: G/GRAPH  '/ &
                               'P,V/PARM  R/RUN  L/LOAD  ', &
                               'W/WRITE  S/SAVE  Q/QUIT')
         ELSE
            WRITE(6,603)
603         FORMAT('## TI MENU: C/CONT  G/GRAPH  ', &
                               'P,V/PARM  R/RUN  L/LOAD  ',&
                               'W/WRITE  S/SAVE  Q/QUIT')
         ENDIF

         CALL TASK_KLIN(LINE,KID,MODE,ti_parm)
      END IF
      CALL mtx_broadcast1_character(KID)
      CALL mtx_broadcast1_integer(MODE)

      IF(MODE.EQ.2) CALL ti_broadcast
      IF(MODE.NE.1) GOTO 1

      SELECT CASE(KID)
      CASE('P')                     ! parameter input
         IF(nrank.eq.0) CALL ti_parm(0,'TI',IERR)
         CALL ti_broadcast
      CASE('V')                     ! view input parameters
         IF(nrank.eq.0) CALL ti_view
      CASE('L')                     ! load bulk component profile from TR
         WRITE(6,*) '## input nid: 0-6'
         READ(5,*,ERR=1,END=1) nid
         CALL pl_load_trdata(nid,ierr)
         if(ierr.ne.0) GO TO 1
         IF(nid.ne.0) THEN
            ALLOCATE(plf(NSMAX))
            DO nr=1,nrmax
               RHON=RM(NR)/RA
               CALL pl_prof(RHON,plf)
               DO NSA=1,nsa_max
                  NS=NS_NSA(NSA)
                  IF(NS.EQ.1.OR.NS.EQ.2) THEN
                     RNA(NSA,NR)=PLF(NS)%RN
                     RTA(NSA,NR)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
                     RUA(NSA,NR)=PLF(NS)%RU
                  END IF
               END DO
            END DO
            DEALLOCATE(plf)
         END IF
      CASE('R','C')                 ! run or continue simulation
         IF(nrank.EQ.0) CALL CPU_TIME(cputime1)
         IF(KID.EQ.'R') THEN
            CALL ti_prep(ierr)
            if(ierr.ne.0) GO TO 1
            INIT=2
         END IF
         IF(INIT.EQ.2) THEN
            CALL ti_exec(ierr)
            if(ierr.ne.0) GO TO 1
            IF(nrank.EQ.0) THEN
               CALL CPU_TIME(cputime2)
               WRITE(6,'(A,F12.3)') '--cpu time =',cputime2-cputime1
            END IF
         END IF
      CASE('G')                     ! graphic output
         IF(INIT.GE.1.AND.nrank.EQ.0) CALL ti_gout
      CASE('W')                     ! show simulation results in text
         IF(INIT.EQ.2.and.nrank.EQ.0) THEN
102         WRITE(6,*) '# SELECT ',  ': PRINT TYPE (1..9)  N/NAMELIST X/EXIT'
            READ(5,'(A1)',ERR=102,END=1) KID
            CALL toupper(KID)
            IF(KID.EQ.'H') THEN
               !            CALL ti_help('W')
            ELSEIF(KID.EQ.'X') THEN
               GOTO 1
            ELSE
               !            CALL ti_print(KID)
            ENDIF
            GOTO 102
         END IF
      CASE('H')                     ! show helpful information
!         CALL ti_help('M')
      CASE('S')                     ! save results
!         CALL ti_save(ierr)
      CASE('Q')                     ! quit 
         GOTO 9000
      CASE('X','#','!')
         CONTINUE
      CASE default
         IF(nrank.EQ.0) WRITE(6,*) 'XX timenu: UNKNOWN KID'
      END SELECT

      GOTO 1

 9000 CONTINUE
      RETURN
    END SUBROUTINE ti_menu
  END MODULE timenu
