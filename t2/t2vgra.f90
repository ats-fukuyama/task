!C------------------------------------------------------------------
!C
!C         MODULE T2VGRA
!C       
!C         VARIABLE GRAPH GENERATOR FOR TASK/T2
!C
!C                       2014-02-22 H.SETO
!C
!C------------------------------------------------------------------
MODULE T2VGRA
  
  USE T2CNST,ONLY:&
       i0ikind,i0rkind
  
  IMPLICIT NONE
  
  PUBLIC T2_VGRA,T2_VGRA_OUTPUT
  
  PRIVATE 
  
CONTAINS 
  
  SUBROUTINE T2_VGRA
    
    USE T2COMM, ONLY: idfile
    
    CALL T2VGRA_VV
    CALL T2VGRA_MS
    CALL T2VGRA_AV
    CALL T2VGRA_AT
    CALL T2VGRA_DT
    CALL T2VGRA_GV
    CALL T2VGRA_GT
    CALL T2VGRA_ES
    CALL T2VGRA_EV
    CALL T2VGRA_ET
    CALL T2VGRA_SS

    IF(idfile.ge.5) CALL T2_VGRA_OUTPUT
    CALL T2_VGRA_OUTPUT
    
    RETURN
    
  END SUBROUTINE T2_VGRA
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR VALIABLE MATRIX ARRAY
  !C 
  !C          MODIFIED 2014-03-06
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_VV
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2vvvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj,&
         i0vofi,i0vofj
    
    !C
    !C INITIALIZATION
    !C
    
    i2vvvt(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_002
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2vvvt(i0vidi,i0vidj) = 1

    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1

    !C Er 
    i0vidj = 5
    i2vvvt(i0vidi,i0vidj) = 1

    !C
    !C EQ_003
    !C
    
    i0vidi = 3

    !C PSI'
    i0vidj = 1
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2vvvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax

       i0vofj = 10*i0sidj

       !C Ft
       i0vidj = i0vofj - 1
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQ_004
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1


    DO i0sidj = 1, i0smax
    
       i0vofj = 10*i0sidj

       !C Fr
       i0vidj = i0vofj - 3
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Fp
       i0vidj = i0vofj 
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQ_005
    !C

    i0vidi = 5
    
    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1

    !C Er
    i0vidj = 5
    i2vvvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax
       
       i0vofj = 10*i0sidj
       
       !C N
       i0vidj = i0vofj - 4
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi

       !C
       !C EQ_006
       !C
       
       i0vidi = i0vofi - 4
       
       !C N
       i0vidj = i0vofi - 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_007
       !C
       
       i0vidi= i0vofi - 3
       
       !C Ep
       i0vidj = 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fr
       i0vidj = i0vofi - 3
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = i0vofi - 2 
       i2vvvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C P 
       i0vidj = i0vofi + 1
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
      
             !C Ft
             i0vidj = i0vofj - 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj 
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qt
             i0vidj = i0vofj + 4
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1
          
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2vvvt(i0vidi,i0vidj) = 1

          !C Qb
          i0vidj = i0vofj + 3
          i2vvvt(i0vidi,i0vidj) = 1

       ENDDO
       
       !C
       !C EQ_009
       !C

       i0vidi = i0vofi - 1
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fr
             i0vidj = i0vofj - 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fb
             i0vidj = i0vofj - 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qb
             i0vidj = i0vofj + 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C 
          !C FOR ANORMALOUS TRANSPORT 
          !C 
          
          !C Ne
          i0vidj = 6
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Fbe
          i0vidj = 8
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Fte 
          i0vidj = 9
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Pe
          i0vidj = 11
          i2vvvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQ_010
       !C

       i0vidi = i0vofi

       !C Fb
       i0vidj = i0vofi - 2 
       i2vvvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fp
       i0vidj = i0vofi
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C N
       i0vidj = i0vofi - 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Ft
       i0vidj = i0vofi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fp
       i0vidj = i0vofi
       i2vvvt(i0vidi,i0vidj) = 1

       !C P
       i0vidj = i0vofi + 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = i0vofi + 3
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qp
       i0vidj = i0vofi + 5
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi + 2
       
       !C Ep
       i0vidj = 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C N
       i0vidj = i0vofi - 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C P
       i0vidj = i0vofi + 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qr
       i0vidj = i0vofi + 2
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = i0vofi + 3
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qt
       i0vidj = i0vofi + 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1             

             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Ft
             i0vidj = i0vofj - 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qt
             i0vidj = i0vofj + 4
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1

          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = i0vofj + 3
          i2vvvt(i0vidi,i0vidj) = 1

       ENDDO

       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4

       DO i0sidj = 1, i0smax
       
          i0vofj = 10*i0sidj

          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = i0vofj - 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fb
             i0vidj = i0vofj - 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fp
             i0vidj = i0vofj
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = i0vofj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qr
             i0vidj = i0vofj + 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qb
             i0vidj = i0vofj + 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qp
             i0vidj = i0vofj + 5
             i2vvvt(i0vidi,i0vidj) = 1

          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2vvvt(i0vidi,i0vidj) = 1              
          
          !C 
          !C FOR ANORMALOUS TRANSPORT 
          !C 
          
          !C Ne
          i0vidj = 6
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Pe
          i0vidj = 11
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qbe
          i0vidj = 13
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qte 
          i0vidj = 14
          i2vvvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQ_015
       !C

       i0vidi = i0vofi + 5
       
       !C Qb
       i0vidj = i0vofi + 3 
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qp
       i0vidj = i0vofi + 5
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_VV
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR MASS SCALAR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_MS
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2msvt
    
    INTEGER(i0ikind)::&
         & i0sidi,i0vidi,i0vofi,&
         &        i0vidj
    
    !C VARIABLE-VARIABLE GRAPH
    
    !C INITIALIZE
    
    i2msvt(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2msvt(i0vidi,i0vidj) = 1

    !C
    !C EQ_002
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2msvt(i0vidi,i0vidj) = 1

    !C
    !C EQ_003
    !C
    
    i0vidi = 3

    !C Et
    i0vidj = 3
    i2msvt(i0vidi,i0vidj) = 1
        
    !C
    !C EQ_004
    !C
    
    i0vidi = 4
    
    !C Ep
    i0vidj = 4
    i2msvt(i0vidi,i0vidj) = 1
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi

       !C
       !C EQ_006
       !C
       
       i0vidi = i0vofi - 4
       
       !C N
       i0vidj = i0vofi - 4
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C Fb
       i0vidj = i0vofi - 2
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_009
       !C
       
       i0vidi = i0vofi - 1
       
       !C Ft
       i0vidj = i0vofi - 1
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C Qb
       i0vidj = i0vofi + 3
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4

       !C Qt
       i0vidj = i0vofi + 4
       i2msvt(i0vidi,i0vidj) = 1              
       
       !C
       !C EQUATION FOR Qp
       !C
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_MS

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR ADVECTION VECTOR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AV
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2avvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0vidj,i0vofi
        
    !C
    !C INITIALIZE
    !C
    
    i2avvt(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2avvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_002
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2avvt(i0vidi,i0vidj) = 1

    !C
    !C EQ_003
    !C
    
    i0vidi = 3

    !C PSI'
    i0vidj = 1
    i2avvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2avvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_004
    !C
    
    i0vidi = 4
        
    !C Ep
    i0vidj = 4
    i2avvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_005
    !C

    i0vidi = 5
    
    !C Ep
    i0vidj = 4
    i2avvt(i0vidi,i0vidj) = 1

    !C Er
    i0vidj = 5
    i2avvt(i0vidi,i0vidj) = 1

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax

       i0vofi  = 10*i0sidi

       !C
       !C EQUATION FOR N
       !C

       i0vidi = i0vofi - 4
              
       !C N
       i0vidj = i0vofi - 4
       i2avvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = i0vofi - 2
       
       !C Fb
       i0vidj = i0vofi - 2
       i2avvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2avvt(i0vidi,i0vidj) = 1
                     
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = i0vofi - 1
       
       !C Ft
       i0vidj = i0vofi - 1
       i2avvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       i2avvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = i0vofi + 3

       !C Fb
       i0vidj = i0vofi - 2
       i2avvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2avvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2avvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = i0vofi + 4
       
       !C Ft
       i0vidj = i0vofi - 1
       i2avvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2avvt(i0vidi,i0vidj) = 1           
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_AV

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR ADVECTION TENSOR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AT
  
    USE T2COMM,ONLY:&
         i0smax,i0wmax,i0vmax,i2atvt,i3atwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2atvt(         1:i0vmax,1:i0vmax) = 0
    i3atwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi

       !C
       !C EQ_008
       !C
       i0vidi = i0vofi - 2
       
       !C Ft (B)
       i0vidj = i0vofi - 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Fp (B)
       i0vidj = i0vofi
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qt (B)
       i0vidj = i0vofi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qp (B)
       i0vidj = i0vofi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_009
       !C
       
       i0vidi = i0vofi - 1
      
       !C Ft (B)
       i0vidj = i0vofi - 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Fp (B)
       i0vidj = i0vofi 
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qt (B)
       i0vidj = i0vofi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qp (B)
       i0vidj = i0vofi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1 
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C Ft (B)
       i0vidj = i0vofi - 1
       i2atvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Fp (B)
       i0vidj = i0vofi
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qt (B)
       i0vidj = i0vofi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qp (B)
       i0vidj = i0vofi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C Ft (B)
       i0vidj = i0vofi - 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Fp (B)
       i0vidj = i0vofi 
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qt (B)
       i0vidj = i0vofi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qp (B)
       i0vidj = i0vofi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4
       
       !C Ft (B)
       i0vidj = i0vofi - 1
       i2atvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Fp (B)
       i0vidj = i0vofi 
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qt (B)
       i0vidj = i0vofi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qp (B)
       i0vidj = i0vofi + 5
       i2atvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_AT
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR DIFFUSION TENSOR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_DT
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2dtvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2dtvt(1:i0vmax,1:i0vmax) = 0
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
      
       !C N
       i0vidj = i0vofi - 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_009
       !C
       
       i0vidi = i0vofi - 1
       
       !C N
       i0vidj = i0vofi - 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C N
       i0vidj = i0vofi - 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C N
       i0vidj = i0vofi - 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4
       
       !C N
       i0vidj = i0vofi - 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2dtvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_DT

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR GRADIENT VECOTR ARRAY
  !C 
  !C                     2014-03-15 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GV
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2gvvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0vidj,i0vofi
    
    !C VARIABLE-VARIABLE GRAPH
    
    !C
    !C INITIALIZATION
    !C

    i2gvvt(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQ_001
    !C
    
    i0vidi = 1
    
    !C Et
    i0vidj = 3
    i2gvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_002
    !C
    
    i0vidi = 2
    
    !C Ep
    i0vidj = 4
    i2gvvt(i0vidi,i0vidj) = 1

    !C Er 
    i0vidj = 5
    i2gvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQ_004
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    i2gvvt(i0vidi,i0vidj) = 1
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_007
       !C
       
       i0vidi = i0vofi - 3
       
       !C P 
       i0vidj = i0vofi + 1
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_009
       !C
       
       i0vidi = i0vofi - 1
       
       !C Ne (ANORMALOUS)
       i0vidj = 6
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C Pe (ANORMALOUS)
       i0vidj = 11
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi + 2
       
       !C N
       i0vidj = i0vofi - 4
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = i0vofi + 1
       i2gvvt(i0vidi,i0vidj) = 1

       !C
       !C EQ_014
       !C
       
       i0vidi = i0vofi + 4
       
       !C Ne (ANORMALOUS)
       i0vidj = 6
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C Pe (ANORMALOUS)
       i0vidj = 11
       i2gvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GV

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR GRADIENT TENSOR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GT
    
    USE T2COMM,ONLY:&
         i0smax,i0wmax,i0vmax,i2gtvt,i3gtwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2gtvt(         1:i0vmax,1:i0vmax) = 0
    i3gtwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_008
       !C

       i0vidi = i0vofi - 2
       
       !C N : (B)
       i0vidj = i0vofi - 4
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb : (B)
       i0vidj = i0vofi - 2
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C P : (B)
       i0vidj = i0vofi + 1
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb: (B)
       i0vidj = i0vofi + 3
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C N: (B), (Ub)
       i0vidj = i0vofi - 4
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb: (B), (Ub)
       i0vidj = i0vofi - 2
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C P: (B), (Ub)
       i0vidj = i0vofi + 1
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Qb: (B), (Ub)
       i0vidj = i0vofi + 3
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C N : (B)
       i0vidj = i0vofi - 4
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb: (B)
       i0vidj = i0vofi - 2
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C P : (B)
       i0vidj = i0vofi + 1
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb: (B)
       i0vidj = i0vofi + 3
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GT

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION SCALAR ARRAY
  !C 
  !C                     2014-03-06
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2VGRA_ES
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2esvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj,&
         i0vofi,i0vofj
    
    !C
    !C INITIALIZATION
    !C
    
    i2esvt(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C EQ_002
    !C
    
    i0vidi = 2
 
    !C Ep
    
    i0vidj = 4
    i2esvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax
       
       i0vofj = 10*i0sidj
       
       !C
       !C EQ_003
       !C
       
       i0vidi = 3
       
       !C Ft
       i0vidj = i0vofj - 1
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_004
       !C
       
       i0vidi = 4
       
       !C Fr
       i0vidj = i0vofj - 3
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Fp
       i0vidj = i0vofj
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_005
       !C
       
       i0vidi = 5
       
       !C N
       i0vidj = i0vofj - 4
       i2esvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_007
       !C
       
       i0vidi= i0vofi - 3
       
       !C Ep
       i0vidj = 4
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Er
       i0vidj = 5
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = i0vofi - 2
       i2esvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2esvt(i0vidi,i0vidj) = 1
       
       
       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj          
          
          !C
          !C EQ_008
          !C
          
          i0vidi = i0vofi - 2
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2esvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = i0vofj + 3
          i2esvt(i0vidi,i0vidj) = 1
          
          !C
          !C EQ_009
          !C
          
          i0vidi = i0vofi - 1
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Fr
             i0vidj = i0vofj - 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Fbe (ANORMALOUS)
             i0vidj = 8
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Fte (ANORMALOUS)
             i0vidj = 9
             i2esvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2esvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQ_010
       !C

       i0vidi = i0vofi

       !C Fb
       i0vidj = i0vofi - 2
       i2esvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = i0vofi - 1
       i2esvt(i0vidi,i0vidj) = 1

       !C Fp
       i0vidj = i0vofi
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C P
       i0vidj = i0vofi + 1
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQ_012
       !C
       
       i0vidi = i0vofi + 2
       
       !C Ep
       i0vidj = 4
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Er
       i0vidj = 5
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = i0vofi + 3
       i2esvt(i0vidi,i0vidj) = 1

       !C Qt
       i0vidj = i0vofi + 4
       i2esvt(i0vidi,i0vidj) = 1

       DO i0sidj = 1, i0smax
          
          i0vofj = 10*i0sidj

          !C
          !C EQ_013
          !C
          
          i0vidi = i0vofi + 3
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2esvt(i0vidi,i0vidj) = 1             
             
          ENDIF
          
          !C Fb
          i0vidj = i0vofj - 2
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = i0vofj + 3
          i2esvt(i0vidi,i0vidj) = 1
          
          !C
          !C EQ_014
          !C
          
          i0vidi = i0vofi + 4
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Qr
             i0vidj = i0vofj + 2
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Qbe (ANORMALOUS)
             i0vidj = 13
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Qte (ANORMALOUS)
             i0vidj = 14
             i2esvt(i0vidi,i0vidj) = 1
          ENDIF
          
          !C Ft
          i0vidj = i0vofj - 1
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = i0vofj + 4
          i2esvt(i0vidi,i0vidj) = 1              
          
       ENDDO
       
       !C
       !C EQ_015
       !C
       
       i0vidi = i0vofi + 5
       
       !C Qb
       i0vidj = i0vofi + 3
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = i0vofi + 4
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qp
       i0vidj = i0vofi + 5
       i2esvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_ES

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION VECTOR ARRAY
  !C 
  !C                     2014-03-06
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EV
    
    USE T2COMM,ONLY:&
         i0smax,i0wmax,i0vmax,i2evvt,i3evwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2evvt(         1:i0vmax,1:i0vmax) = 0
    i3evwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
    !C
    !C EQ_003
    !C
    
    i0vidi = 3
    
    !C PSI': (R)
    i0vidj = 1
    i2evvt(i0vidi,i0vidj) = 1
    
    i0widi = 2
    i3evwt(i0widi,i0vidi,i0vidj) = 1
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_008
       !C

       i0vidi = i0vofi - 2
       
       !C Fb: (B)
       i0vidj = i0vofi - 2
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQ_013
       !C
       
       i0vidi = i0vofi + 3
       
       !C Et: (B) (Ub) (Qb/P)
       i0vidj = 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Ep: (B) (Ub) (Qb/P)
       i0vidj = 4
       i2evvt(i0vidi,i0vidj) = 1             
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb: (B)
       i0vidj = i0vofi - 2
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb: (B)
       i0vidj = i0vofi + 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = i0vofi + 4
       
       !C Et: B, Ub, Wb
       i0vidj = 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Ep: B, Ub, Wb
       i0vidj = 4
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_EV

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION TENSOR ARRAY
  !C 
  !C                     2014-02-20 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ET

    USE T2COMM,ONLY:&
         i0smax,i0wmax,i0vmax,i2etvt,i4etwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,&
                i0widj,i0vidj,i0vofi
    
    !C
    !C INITIALIZATION
    !C
    
    i2etvt(                  1:i0vmax,1:i0vmax) = 0
    i4etwt(1:i0wmax,1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
    DO i0sidi = 1, i0smax
       
       i0vofi = 10*i0sidi
       
       !C
       !C EQ_008
       !C
       
       i0vidi = i0vofi - 2
       
       !C Ft: (B,B)
       i0vidj = i0vofi - 1
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Fp: (B,B)
       i0vidj = i0vofi
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qt (B,B)
       i0vidj = i0vofi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qp (B,B)
       i0vidj = i0vofi + 5
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C
       !C EQ_011
       !C
       
       i0vidi = i0vofi + 1
       
       !C Ft: (B,B), (Ub,B)
       i0vidj = i0vofi - 1
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Fp: (B,B), (Ub,B)
       i0vidj = i0vofi
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qt: (B,B), (Ub,B)
       i0vidj = i0vofi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Qp: (B,B), (Ub,B)
       i0vidj = i0vofi + 5
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = i0vofi + 3
       
       !C Ft: (B,B)
       i0vidj = i0vofi - 1
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Fp: (B,B)
       i0vidj = i0vofi
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Qt: (B,B)
       i0vidj = i0vofi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qp: (B,B)
       i0vidj = i0vofi + 5
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN 

  END SUBROUTINE T2VGRA_ET
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR SOURCE SCALAR ARRAY
  !C 
  !C                     2014-03-06 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_SS
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2ssvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2ssvt(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2ssvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2ssvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vidi = 3

    !C Et
    i0vidj = 3
    i2ssvt(i0vidi,i0vidj) = 1
        
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    !C Ep
    i0vidj = 4
    i2ssvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Er
    !C

    i0vidi = 5
    
    !C Ep
    i0vidj = 5
    i2ssvt(i0vidi,i0vidj) = 1
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       !C
       !C EQUATION FOR N
       !C

       i0vidi = 10*i0sidi - 4
              
       !C N
       i0vidj = 10*i0sidi - 4
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi = 10*i0sidi - 3
       
       !C Qb
       i0vidj = 10*i0sidi - 3
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 10*i0sidi - 2
       
       !C Fb
       i0vidj = 10*i0sidi - 2
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = 10*i0sidi - 1
       
       !C Ft
       i0vidj = 10*i0sidi - 1
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Fp
       !C

       i0vidi = 10*i0sidi
       
       !C Fp
       i0vidj = 10*i0sidi 
       i2ssvt(i0vidi,i0vidj) = 1

       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 10*i0sidi + 1
       
       !C P
       i0vidj = 10*i0sidi + 1
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 10*i0sidi + 2
       
       !C Qr
       i0vidj = 10*i0sidi + 2
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 10*i0sidi + 3
       
       !C Qb
       i0vidj = 10*i0sidi + 3
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 10*i0sidi + 4

       !C Qt
       i0vidj = 10*i0sidi + 4
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qp
       !C
       
       i0vidi = 10*i0sidi + 5

       !C Qp
       i0vidj = 10*i0sidi + 5
       i2ssvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    RETURN

  END SUBROUTINE T2VGRA_SS

    


  SUBROUTINE T2_VGRA_OUTPUT

    USE T2COMM
    INTEGER(i0ikind):: i0sidi,i0sidj,i0widi,i0widj

    OPEN(10,FILE='TEST_VV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'VV=',i2vvvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_MS.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'MS=',i2msvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'AV=',i2avvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'AT=',i2atvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_DT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'DT=',i2dtvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'GV=',i2gvvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'GT=',i2gtvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ES.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'ES=',i2esvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_EV.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'EV=',i2evvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ET.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'ET=',i2etvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_SS.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'SS=',i2ssvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I3ATWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,&
               'ATWT=',i3atwt(i0widi,i0sidi,i0sidj)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I3GTWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,&
               'GTWT=',i3gtwt(i0widi,i0sidi,i0sidj)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I3EVWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,&
               'EVWT=',i3evwt(i0widi,i0sidi,i0sidj)
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_I4ETWT.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       DO i0widj = 1, i0wmax
       DO i0widi = 1, i0wmax
          WRITE(10,*)'vi=',i0sidi,'vj=',i0sidj,'wi=',i0widi,'wj=',i0widj,&
               'ETWT=',i4etwt(i0widi,i0widj,i0sidi,i0sidj)
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_SS.dat')
    DO i0sidi = 1, i0vmax
    DO i0sidj = 1, i0vmax
       WRITE(10,*)'i=',i0sidi,'j=',i0sidj,'SS=',i2ssvt(i0sidi,i0sidj)
    ENDDO
    ENDDO
    CLOSE(10)

    RETURN
  END SUBROUTINE T2_VGRA_OUTPUT
    
END MODULE T2VGRA
