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
    IMPLICIT NONE

    INTEGER(i0ikind)::i0mlva,i0mlvb,i0mlvc,i1,j1,i0mesh_level
        
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

    CALL T2NGRA_ALLOCATE_1

    DO i1=1,i0lmax
       i0mesh_level=i1mlvl(i1)-1
       i1pdn2(i1) = i0pdiv_number*(2**i0mesh_level)   
    ENDDO
    
    !C CALCULATE NUMBER OF NODES IN EACH MESH LEVEL
    !C          I0NMAX1, i0NMAX2, I0NMAX3
        
    i0nmax1=0
    i0nmax2=1
    i0nmax3=0
    
    DO i1 =1,i0lmax
       i1rdn1(i1)= i1rdn2(i1)+1
       i1pdn1(i1)= i1pdn2(i1)+1
    ENDDO
    
    SELECT CASE (i0nmax0)
       !C  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !C  8: QADRADIC RECTANGULAR ELEMENT (  8 POINTS)
       !C 12: CUBIC    RECTANGULAR ELEMENT ( 12 POINTS)
    CASE(4)
       DO i1=1,i0lmax
          i0mlva=i1mlvl(i1-1)
          i0mlvb=i1mlvl(i1)
          i0mlvc=i1mlvl(i1+1)
          i1nmax1(i1) = i1rdn1(i1)*i1pdn1(i1)
          IF(i0mlva.EQ.0)THEN
             i1nmax2(i1)= i1rdn2(i1)*i1pdn2(i1)+1
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             i1nmax2(i1)= i1rdn2(i1)*i1pdn2(i1)+i1pdn2(i1-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             i1nmax2(i1)= i1rdn2(i1)*i1pdn2(i1)+i1pdn2(i1)
          ENDIF
          i0nmax1=i0nmax1+i1rdn1(i1)*i1pdn1(i1)
          i0nmax2=i0nmax2+i1rdn2(i1)*i1pdn2(i1)
       ENDDO
       
       i0nmax3=i0nmax2
       
       DO i1=2,i0lmax
          i0mlva=i1mlvl(i1-1)
          i0mlvb=i1mlvl(i1)
          IF(i0mlva.NE.i0mlvb)i0nmax3=i0nmax3+i1pdn2(i1-1)
       ENDDO
       
       i0ermx = i0nmax3+1
       i0ecmx = i1pdn2(1) + 4*(i0nmax2 - i1pdn1(i0lmax))&
            + 2*(i0nmax3-i0nmax2 +i1pdn2(i0lmax))
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
    

    !C CALCULATE NUMBER OF ELEMENTS IN EACH MESH LEVEL
    !C                  I0EMAX, I1EMAX

    
    i1emax(0:i0lmax)= 0
    DO i1=1,i0lmax
       DO j1=1,i1
          i1emax(i1)=i1emax(i1)+i1rdn2(j1)*i1pdn2(j1)
       ENDDO
    ENDDO
    
    i0emax=i1emax(i0lmax)


    !C CALCULATE NUMBER OF HANGED NODE
    !C                  I0HMAX

    i0hmax=0

    DO i1=2,i0lmax
       i0mlva=i1mlvl(i1-1)
       i0mlvb=i1mlvl(i1)
       IF(i0mlva.NE.i0mlvb)i0hmax=i0hmax+i1pdn2(i1-1)
    ENDDO
    

    !C CALCULATE ARRAY SIZE OF MATRIX INDICES FOR NODE
    !C I0NRMX,I0NCMX
    
    i0nrmx = 1 + i0nmax2
    i0ncmx = 0
    SELECT CASE (i0nmax0)
       !C  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !C  8: QADRADIC RECTANGULAR ELEMENT (  8 POINTS)
       !C 12: CUBIC    RECTANGULAR ELEMENT ( 12 POINTS)
    CASE( 4)
       
       DO i1=1,i0lmax
          i0mlva=i1mlvl(i1-1)
          i0mlvb=i1mlvl(i1)
          i0mlvc=i1mlvl(i1+1)
          !C POINTS ON LEFT-SIDE-EDGES 
          IF(    i0mlva.EQ.0)THEN
             !C FIRST AND SECOND EDGE POINTS ON LEFT-SIDE IN Lv-1 MESH
             i0ncmx= i0ncmx + 1+i1pdn2(i1)   +  7*i1pdn2(i1)
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             !C SECOND EDGE POINTS ON LEFT-SIDES IN Lv-2~LMAX MESH
             i0ncmx= i0ncmx + 8*i1pdn2(i1-1) + 11*i1pdn2(i1-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             i0ncmx= i0ncmx + 9*i1pdn2(i1)
          ENDIF
          !C CONTSRAINT FREE POINTS ON MIDDLE AREA
          IF(i1rdn2(i1).GE.2)THEN
             i0ncmx= i0ncmx + 9*i1pdn2(i1)*(i1rdn2(i1)-2)
          ENDIF
          !C FIRST EDGE POINTS ON RIGHT-SIDE
          IF(    i0mlvc.EQ.0     )THEN  
             i0ncmx= i0ncmx +  6*i1pdn2(i1)
          ELSEIF(i0mlvb.NE.i0mlvc)THEN
             i0ncmx= i0ncmx + 11*i1pdn2(i1)
          ELSEIF(i0mlvb.EQ.i0mlvc)THEN
             i0ncmx= i0ncmx + 9*i1pdn2(i1)
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
    
    i0vmax  = 8*i0spcs + 5
    i0vgrmx = i0vmax + 1
    i0vgcmx = 9*(i0spcs**2) + 48*i0spcs + 11
    
    i0vmax  = 8*i0spcs + 5
    i0wmax  = 5*i0spcs + 1
    i0vgrmx = i0vmax   + 1
    i0msrmx = i0vgrmx
    i0avrmx = i0vgrmx
    i0atrmx = i0vgrmx
    i0dtrmx = i0vgrmx
    i0gvrmx = i0vgrmx
    i0gtrmx = i0vgrmx
    i0esrmx = i0vgrmx
    i0evrmx = i0vgrmx
    i0etrmx = i0vgrmx
    i0ssrmx = i0vgrmx
    
    i0vgcmx = 11 + 49*i0spcs + 8*(i0spcs**2)
    i0mscmx =  5 +  8*i0spcs    
    i0avcmx =  7 + 10*i0spcs
    i0atcmx =      20*i0spcs
    i0dtcmx =      20*i0spcs
    i0gvcmx =  4 +  4*i0spcs
    i0gtcmx =      16*i0spcs
    i0escmx =      21*i0spcs + 8*(i0spcs**2)
    i0evcmx =      23*i0spcs
    i0etcmx =      16*i0spcs
    i0sscmx =  5 +  8*i0spcs
     
    i0cmax = i0vgcmx*i0ncmx
    i0bmax = i0vmax *i0nmax2
    i0xmax = i0vmax *i0nmax3
    
    CALL T2COMM_ALLOCATE

    RETURN

  END SUBROUTINE T2_PREP
END MODULE T2PREP
