!--------------------------------------------------------------------- 
!
!    MODULE FOR GENERATING GLOBAL PARAMETERS AND INITIAL PROFILES
!
!
!                   LAST UPDATE 2014-05-28 H.Seto
!
!---------------------------------------------------------------------
MODULE T2PREP
  
  USE T2CNST,ONLY: ikind,rkind
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC T2_PREP
  
CONTAINS

  !-------------------------------------------------------------------
  !
  !
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2_PREP
    
    USE T2COMM,ONLY: time_t2,time_init
    USE T2DIV, ONLY: T2_DIV
    USE T2PROF,ONLY: T2_PROF
    
    REAL(4):: e0time_0,e0time_1
    
    !d0iar   = d0rmnr/d0rmjr
    
    ! set up normalization factors for T2CALV
    CALL T2PREP_SETUP_NORMALIZATION_FACTOR 
    
    ! set up global parameters for T2COMM
    CALL T2PREP_SETUP_GLOBAL_PARAMETER
    
    ! set up initial profiles
    CALL T2PREP_SETUP_INITIAL_PROFILES
    
    CALL T2_DIV
    
    CALL CPU_TIME(e0time_0)
    CALL T2_PROF
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- Initial profile calculated: cpu=', &
                         e0time_1-e0time_0,' [s]'

 !   IF(idfile.ge.1) CALL T2_WRIT_MAIN(d1guv,0,c10rname)
 !   IF(idfile.ge.2) CALL T2_WRIT_GPT(20,0,d1guv)
 !   IF(idfile.ge.3) CALL T2_WRIT_GP1(22,0,d1guv)

    time_t2=time_init

    RETURN

  END SUBROUTINE T2_PREP
  
  !--------------------------------------------------------------
  ! 
  !       SET NORMALIZATION FACTORS FOR EQAUTIONS 
  !       IN ORDER TO IMPROVE CONVERGENCE CHARACTERISTIC
  !
  !                           LAST UPDATE 2014-05-28 H.Seto          
  !
  !--------------------------------------------------------------
  SUBROUTINE T2PREP_SETUP_NORMALIZATION_FACTOR
    USE T2CNST,ONLY:Aee
    USE T2COMM,ONLY:UsePotentialDescription,&
         &          UseNormalization,&
         !
         &          BpNF, BtNF, EtNF, EpNF, ErNF,&
         &          NnNF, FrNF, FbNF, FtNF, FpNF,&
         &          PpNF, QrNF, QbNF, QtNF, QpNF
    
    IF(.NOT.UsePotentialDescription)THEN
       IF(UseNormalization)THEN
          BpNF = 1.D0                  ! [        T    ]
          BtNF = 1.D0                  ! [        T    ]
          EtNF = 1.D0                  ! [        V/m  ]
          EpNF = 1.D0                  ! [        V/m  ]
          ErNF = 1.D3                  ! [       kV/m  ] 
          NnNF = 1.D20                 ! [10^20    /m3 ]
          FrNF = NnNF                  ! [10^20    /m2s]
          FbNF = NnNF*1.D3             ! [10^23    /m2s]
          FtNF = NnNF*1.D3             ! [10^23    /m2s]
          FpNF = NnNF*1.D3             ! [10^23    /m2s]
          PpNF = NnNF*1.D3*Aee         ! [10^20 keV/m3 ]
          QrNF = FrNF*1.D3*Aee         ! [10^20 keV/m2s]
          QbNF = FbNF*1.D3*Aee         ! [10^23 keV/m2s]
          QtNF = FtNF*1.D3*Aee         ! [10^23 keV/m2s]
          QpNF = FpNF*1.D3*Aee         ! [10^23 keV/m2s]
       ELSE
          BpNF = 1.D0         ! 
          BtNF = 1.D0         ! 
          EtNF = 1.D0         !
          EpNF = 1.D0         !
          ErNF = 1.D0         !
          NnNF = 1.D0         !
          FrNF = 1.D0         !
          FbNF = 1.D0         !
          FtNF = 1.D0         !
          FpNF = 1.D0         !
          PpNF = 1.D0         !
          QrNF = 1.D0         !
          QbNF = 1.D0         !
          QtNF = 1.D0         !
          QpNF = 1.D0         !
       ENDIF
    ELSE
       WRITE(6,*)"T2PREP_SETUP_NORMALIZATION_FACTOR              "
       WRITE(6,*)"Potential Description Ver. is underconstruction"
       STOP
    ENDIF
    
    RETURN
  
  END SUBROUTINE T2PREP_SETUP_NORMALIZATION_FACTOR

  !--------------------------------------------------------------
  ! 
  !       SET GLOBAL PARAMETERS FOR T2COMM_ALLOCATE
  !
  !                           LAST UPDATE 2014-05-28 H.Seto          
  !
  !--------------------------------------------------------------
  SUBROUTINE T2PREP_SETUP_GLOBAL_PARAMETER
    
    USE T2COMM,ONLY:&
         NLMAX, NVMAX, NKMAX, NQMAX, NNMAX, NMMAX, NXMAX, &
         NBMAX, NRMAX, NEMAX, NHMAX, NAMAX, NNRMX, NERMX, &
         NECMX, NDMAX, NSMAX, NPMAX, NPMIN,        &
         !
         UsePotentialDescription,&
         !
         i1mlvl,i1pdn1,i1pdn2,i1rdn1,i1rdn2,i1mmax,i1bmax,i1emax,&
         !
         T2NGRA_ALLOCATE, T2COMM_ALLOCATE

    INTEGER(ikind)::&
         i0mlva,i0mlvb,i0mlvc,i0lidi,i0lidj,i0mlvl
    
    CALL T2NGRA_ALLOCATE
    
    DO i0lidi = 1, NLMAX
       i0mlvl=i1mlvl(i0lidi)-1
       i1pdn2(i0lidi) = NPMIN*(2**i0mlvl)   
    ENDDO
    
    !C
    !C CALCULATE NUMBER OF NODES IN EACH MESH LEVEL
    !C 
    
    NMMAX = 0
    NBMAX = 1
    NXMAX = 0
    NRMAX = 1
    
    DO i0lidi = 1, NLMAX
       i1rdn1(i0lidi) = i1rdn2(i0lidi) + 1
       i1pdn1(i0lidi) = i1pdn2(i0lidi) + 1
    ENDDO
    
    SELECT CASE (NNMAX)
       !C  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !C  9: QADRADIC RECTANGULAR ELEMENT (  9 POINTS)
       !C 16: CUBIC    RECTANGULAR ELEMENT ( 16 POINTS)
    CASE(4)
       DO i0lidi = 1, NLMAX
          i0mlva = i1mlvl(i0lidi-1)
          i0mlvb = i1mlvl(i0lidi  )
          i0mlvc = i1mlvl(i0lidi+1)
          i1mmax(i0lidi) = i1rdn1(i0lidi)*i1pdn1(i0lidi)
          IF(i0mlva.EQ.0)THEN
             i1bmax(i0lidi) = i1rdn2(i0lidi)*i1pdn2(i0lidi)+1
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             i1bmax(i0lidi) = i1rdn2(i0lidi)*i1pdn2(i0lidi)+i1pdn2(i0lidi-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             i1bmax(i0lidi) = i1rdn2(i0lidi)*i1pdn2(i0lidi)+i1pdn2(i0lidi)
          ENDIF
          NMMAX = NMMAX + i1rdn1(i0lidi)*i1pdn1(i0lidi)
          NBMAX = NBMAX + i1rdn2(i0lidi)*i1pdn2(i0lidi)
       ENDDO
       
       NXMAX = NBMAX
       
       DO i0lidi = 2, NLMAX
          i0mlva=i1mlvl(i0lidi-1)
          i0mlvb=i1mlvl(i0lidi  )
          IF(i0mlva.NE.i0mlvb) NXMAX = NXMAX + i1pdn2(i0lidi-1)
       ENDDO
       
       DO i0lidi = 1, NLMAX
          NRMAX = NRMAX + i1rdn2(i0lidi)
       ENDDO
       
       NERMX = NXMAX + 1
       NECMX = i1pdn2(1) + 4*(NBMAX - i1pdn1(NLMAX))&
            + 2*(NXMAX - NBMAX +i1pdn2(NLMAX))
    CASE (8)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE (12)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE DEFAULT
       WRITE(6,*)'T2_GRPH: IMPROPER IMPUT NNMAX0=',NNMAX
       STOP
    END SELECT
    
    !C
    !C CALCULATE NUMBER OF ELEMENTS IN EACH MESH LEVEL
    !C                  NEMAX, I1EMAX

    i1emax(0:NLMAX)= 0
    
    DO i0lidi = 1, NLMAX
       DO i0lidj = 1, i0lidi
          i1emax(i0lidi) = i1emax(i0lidi) + i1rdn2(i0lidj)*i1pdn2(i0lidj)
       ENDDO
    ENDDO
    
    NEMAX = i1emax(NLMAX)

    !C
    !C CALCULATE NUMBER OF HANGED NODE
    !C                  NHMAX

    NHMAX = 0

    DO i0lidi = 2, NLMAX
       i0mlva=i1mlvl(i0lidi-1)
       i0mlvb=i1mlvl(i0lidi  )
       IF(i0mlva.NE.i0mlvb) NHMAX = NHMAX + i1pdn2(i0lidi-1)
    ENDDO
    

    !C CALCULATE ARRAY SIZE OF MATRIX INDICES FOR NODE
    !C NNRMX,I0NCMX
    
    NNRMX = 1 + NBMAX
    NAMAX = 0
    SELECT CASE (NNMAX)
       !C  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !C  9: QADRADIC RECTANGULAR ELEMENT (  9 POINTS)
       !C 16: CUBIC    RECTANGULAR ELEMENT ( 16 POINTS)
    CASE( 4)
       
       DO i0lidi = 1, NLMAX
          i0mlva = i1mlvl(i0lidi-1)
          i0mlvb = i1mlvl(i0lidi  )
          i0mlvc = i1mlvl(i0lidi+1)
          !C POINTS ON LEFT-SIDE-EDGES 
          IF(    i0mlva.EQ.0)THEN
             !C FIRST AND SECOND EDGE POINTS ON LEFT-SIDE IN Lv-1 MESH
             NAMAX = NAMAX + 1+i1pdn2(i0lidi  ) +  7*i1pdn2(i0lidi  )
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             !C SECOND EDGE POINTS ON LEFT-SIDES IN Lv-2~LMAX MESH
             NAMAX = NAMAX + 8*i1pdn2(i0lidi-1) + 11*i1pdn2(i0lidi-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             NAMAX = NAMAX + 9*i1pdn2(i0lidi  )
          ENDIF
          !C CONTSRAINT FREE POINTS ON MIDDLE AREA
          IF(i1rdn2(i0lidi).GE.2)THEN
             NAMAX = NAMAX + 9*i1pdn2(i0lidi)*(i1rdn2(i0lidi)-2)
          ENDIF
          !C FIRST EDGE POINTS ON RIGHT-SIDE
          IF(    i0mlvc.EQ.0     )THEN  
             NAMAX = NAMAX +  6*i1pdn2(i0lidi)
          ELSEIF(i0mlvb.NE.i0mlvc)THEN
             NAMAX = NAMAX + 11*i1pdn2(i0lidi)
          ELSEIF(i0mlvb.EQ.i0mlvc)THEN
             NAMAX = NAMAX +  9*i1pdn2(i0lidi)
          ELSE
             WRITE(6,*)'ERROR IN SET_MATRIX_INDICES_FOR_N-CRS'
             STOP
          ENDIF
       ENDDO
    CASE (9)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE (16)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE DEFAULT
       WRITE(6,*)'T2PREP: IMPROPER IMPUT NNMAX=',NNMAX 
       STOP
    END SELECT
    
    IF(.NOT.UsePotentialDescription)THEN
       NVMAX  = 10*NSMAX + 5
       NKMAX  = 2*NSMAX + 2 
    ELSE
       WRITE(6,*)"Potential Description Ver. is underconstruction"
       STOP
    END IF
    
    WRITE(6,*)'NDMAX=',NDMAX,'NNMAX=',NNMAX,'NSMAX=',NSMAX
    WRITE(6,*)'NLMAX=',NLMAX,'NPMAX=',NPMAX,'NQMAX=',NQMAX
    WRITE(6,*)'NMMAX=',NMMAX,'NRMAX=',NRMAX,'NEMAX=',NEMAX
    WRITE(6,*)'NAMAX=',NAMAX,'NXMAX=',NXMAX,'NBMAX=',NBMAX
    WRITE(6,*)'NNRMX=',NNRMX,'NERMX=',NERMX,'NECMX=',NECMX
    WRITE(6,*)'NKMAX=',NKMAX,'NHMAX=',NHMAX

    CALL T2COMM_ALLOCATE
    
    RETURN
  END SUBROUTINE T2PREP_SETUP_GLOBAL_PARAMETER

END MODULE T2PREP
