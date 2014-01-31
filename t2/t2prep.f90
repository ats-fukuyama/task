!C 
!C
!C
MODULE T2PREP
  
  PRIVATE
  PUBLIC T2_PREP
  
CONTAINS

  !C
  !C
  !C

  SUBROUTINE T2_PREP
    
    USE T2CNST, ONLY: &
         i0ikind,i0rkind,i0lmaxm,i0spcsm,d0aee,d0ame,d0amp
    USE T2COMM
    USE T2DIV, ONLY: T2_DIV
    USE T2PROF, ONLY: T2_PROF
!    USE T2WRIT, ONLY: T2_WRIT_MAIN,T2_WRIT_GPT,T2_WRIT_GP1
    IMPLICIT NONE
    INTEGER(i0ikind)::i0mlva,i0mlvb,i0mlvc,i1,j1,i0mesh_level
    REAL(4):: e0time_0,e0time_1
        
    d0iar = d0rmnr/d0rmjr
    d0mfcst = 1.D0
    d0btcst = 1.D0/d0iar
    d0ercst = 1.D0
    d0epcst = 1.D0
    d0etcst = 1.D0
    d0nncst = 1.D20
    d0frcst = 1.D20 
    d0fbcst = 1.D20
    d0ftcst = 1.D20/d0iar
    d0ppcst = d0aee*1.D23
    d0qrcst = d0aee*1.D23
    d0qbcst = d0aee*1.D23
    d0qtcst = d0aee*1.D23/d0iar

    CALL T2NGRA_ALLOCATE

    DO i0lidi = 1, i0lmax
       i0mesh_level=i1mlvl(i0lidi)-1
       i1pdn2(i0lidi) = i0pdiv_number*(2**i0mesh_level)   
    ENDDO
    
    !C CALCULATE NUMBER OF NODES IN EACH MESH LEVEL
    !C 
        
    !i0nmax1=0
    !i0nmax2=1
    !i0nmax3=0
    !i0nmax4=1

    i0mmax = 0
    i0bmax = 1
    i0xmax = 0
    i0rmax = 1
    
    DO i0lidi = 1, i0lmax
       i1rdn1(i0lidi) = i1rdn2(i0lidi) + 1
       i1pdn1(i0lidi) = i1pdn2(i0lidi) + 1
    ENDDO
    
    SELECT CASE (i0nmax)
       !C  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !C  8: QADRADIC RECTANGULAR ELEMENT (  8 POINTS)
       !C 12: CUBIC    RECTANGULAR ELEMENT ( 12 POINTS)
    CASE(4)
       DO i0lidi = 1, i0lmax
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
          i0mmax = i0mmax + i1rdn1(i0lidi)*i1pdn1(i0lidi)
          i0bmax = i0bmax + i1rdn2(i0lidi)*i1pdn2(i0lidi)
       ENDDO
       
       i0nmax3 = i0nmax2
       
       DO i0lidi = 2, i0lmax
          i0mlva=i1mlvl(i0lidi-1)
          i0mlvb=i1mlvl(i0lidi  )
          IF(i0mlva.NE.i0mlvb) i0xmax = i0xmax + i1pdn2(i0lidi-1)
       ENDDO
       
       DO i0lidi = 1, i0lmax
          i0rmax = i0rmax + i1rdn2(i0lidi)
       ENDDO
       
       i0ermx = i0xmax + 1
       i0ecmx = i1pdn2(1) + 4*(i0bmax - i1pdn1(i0lmax))&
            + 2*(i0xmax - i0bmax +i1pdn2(i0lmax))
    CASE (8)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE (12)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE DEFAULT
       WRITE(6,*)'T2_GRPH: IMPROPER IMPUT I0NMAX0=',i0nmax
       STOP
    END SELECT
    
    !C
    !C CALCULATE NUMBER OF ELEMENTS IN EACH MESH LEVEL
    !C                  I0EMAX, I1EMAX

    i1emax(0:i0lmax)= 0
    
    DO i0lidi = 1, i0lmax
       DO i0lidj = 1, i0lidi
          i1emax(i0lidi) = i1emax(i0lidi) + i1rdn2(i0lidj)*i1pdn2(i0lidj)
       ENDDO
    ENDDO
    
    i0emax = i1emax(i0lmax)

    !C
    !C CALCULATE NUMBER OF HANGED NODE
    !C                  I0HMAX

    i0hmax = 0

    DO i0lidi = 2, i0lmax
       i0mlva=i1mlvl(i0lidi-1)
       i0mlvb=i1mlvl(i0lidi  )
       IF(i0mlva.NE.i0mlvb) i0hmax = i0hmax + i1pdn2(i0lidi-1)
    ENDDO
    

    !C CALCULATE ARRAY SIZE OF MATRIX INDICES FOR NODE
    !C I0NRMX,I0NCMX
    
    i0nrmx = 1 + i0bmax
    i0ncmx = 0
    SELECT CASE (i0nmax)
       !C  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !C  8: QADRADIC RECTANGULAR ELEMENT (  8 POINTS)
       !C 12: CUBIC    RECTANGULAR ELEMENT ( 12 POINTS)
    CASE( 4)
       
       DO i0lidi = 1, i0lmax
          i0mlva = i1mlvl(i0lid-1)
          i0mlvb = i1mlvl(i0lid  )
          i0mlvc = i1mlvl(i0lid+1)
          !C POINTS ON LEFT-SIDE-EDGES 
          IF(    i0mlva.EQ.0)THEN
             !C FIRST AND SECOND EDGE POINTS ON LEFT-SIDE IN Lv-1 MESH
             i0ncmx = i0ncmx + 1+i1pdn2(i0lidi  ) +  7*i1pdn2(i0lidi  )
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             !C SECOND EDGE POINTS ON LEFT-SIDES IN Lv-2~LMAX MESH
             i0ncmx = i0ncmx + 8*i1pdn2(i0lidi-1) + 11*i1pdn2(i0lidi-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             i0ncmx = i0ncmx + 9*i1pdn2(i0lidi  )
          ENDIF
          !C CONTSRAINT FREE POINTS ON MIDDLE AREA
          IF(i1rdn2(i0lidi).GE.2)THEN
             i0ncmx = i0ncmx + 9*i1pdn2(i0lidi)*(i1rdn2(i0lidi)-2)
          ENDIF
          !C FIRST EDGE POINTS ON RIGHT-SIDE
          IF(    i0mlvc.EQ.0     )THEN  
             i0ncmx = i0ncmx +  6*i1pdn2(i0lidi)
          ELSEIF(i0mlvb.NE.i0mlvc)THEN
             i0ncmx = i0ncmx + 11*i1pdn2(i0lidi)
          ELSEIF(i0mlvb.EQ.i0mlvc)THEN
             i0ncmx = i0ncmx +  9*i1pdn2(i0lidi)
          ELSE
             WRITE(6,*)'ERROR IN SET_MATRIX_INDICES_FOR_N-CRS'
             STOP
          ENDIF
       ENDDO
    CASE (8)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE (12)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE DEFAULT
       WRITE(6,*)'T2_GRPH: IMPROPER IMPUT I0NMAX0=',i0nmax0 
       STOP
    END SELECT
    
    !C
    !C
    !C
    
    i0vmax  = 8*i0spcs + 5
    i0wmax  = 2*i0spcs + 2 
    i0cmax = i0vmax*i0ncmx! DUMMY
    i0bmax = i0vmax *i0nmax2
    i0xmax = i0vmax *i0nmax3

    CALL T2COMM_ALLOCATE

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
END MODULE T2PREP
