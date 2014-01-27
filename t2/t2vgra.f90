!C------------------------------------------------------------------
!C
!C         MODULE T2VGRA
!C       
!C         VARIABLE GRAPH GENERATOR FOR TASK/T2
!C
!C         MODIFIED 2014-01-27
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
    
    IF(idfile.ge.5) CALL T2_VGRA_OUTPUT
    
    RETURN
    
  END SUBROUTINE T2_VGRA
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR VALIABLE MATRIX ARRAY
  !C 
  !C          MODIFIED 2014-01-27
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_VV
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2vvtbl
    
    INTEGER(i0ikind)::&
         i0sida,i0sidb,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2vvtbl(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    
    i0vida = 1
    
    !C PSI'   
    i0vidb = 1
    i2vvtbl(i0vida,i0vidb) = 1
    
    !C Et
    i0vidb = 3
    i2vvtbl(i0vida,i0vidb) = 1
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vida = 2
    
    !C I
    i0vidb = 2
    i2vvtbl(i0vida,i0vidb) = 1

    !C Ep
    i0vidb = 4
    i2vvtbl(i0vida,i0vidb) = 1

    !C Er 
    i0vidb = 5
    i2vvtbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vida = 3

    !C PSI'
    i0vidb = 1
    i2vvtbl(i0vida,i0vidb) = 1
    
    !C Et
    i0vidb = 3
    i2vvtbl(i0vida,i0vidb) = 1
    
    DO i0sidb = 1, i0spcs

       !C Ft
       i0vidb = 8*i0sidb + 1
       i2vvtbl(i0vida,i0vidb) = 1

    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vida = 4
    
    !C I
    i0vidb = 2
    i2vvtbl(i0vida,i0vidb) = 1
    
    !C Ep
    i0vidb = 4
    i2vvtbl(i0vida,i0vidb) = 1


    DO i0sidb = 1, i0spcs
    
       !C Fb
       i0vidb = 8*i0sidb
       i2vvtbl(i0vida,i0vidb) = 1

       !C Ft
       i0vidb = 8*i0sidb + 1
       i2vvtbl(i0vida,i0vidb) = 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C

    i0vida = 5
    
    !C Ep
    i0vidb = 4
    i2vvtbl(i0vida,i0vidb) = 1

    !C Er
    i0vidb = 5
    i2vvtbl(i0vida,i0vidb) = 1
    
    DO i0sidb = 1, i0spcs
       
       !C N
       i0vidb = 8*i0sidb - 2
       i2vvtbl(i0vida,i0vidb) = 1

    ENDDO

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C

       i0vida = 8*i0sida - 2
              
       !C N
       i0vidb = 8*i0sida - 2
       i2vvtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vida= 8*i0sida - 1
       
       !C Ep
       i0vidb = 4
       i2vvtbl(i0vida,i0vidb) = 1

       !C Er
       i0vidb = 5
       i2vvtbl(i0vida,i0vidb) = 1

       !C Fr
       i0vidb = 8*i0sida - 1
       i2vvtbl(i0vida,i0vidb) = 1

       !C Fb
       i0vidb = 8*i0sida 
       i2vvtbl(i0vida,i0vidb) = 1

       !C Ft
       i0vidb = 8*i0sida + 1
       i2vvtbl(i0vida,i0vidb) = 1

       !C P 
       i0vidb = 8*i0sida + 2
       i2vvtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN
             
             !C Et
             i0vidb = 3
             i2vvtbl(i0vida,i0vidb) = 1
             
             !C Ep
             i0vidb = 4
             i2vvtbl(i0vida,i0vidb) = 1
             
             !C N
             i0vidb = 8*i0sidb - 2
             i2vvtbl(i0vida,i0vidb) = 1

             
             !C Ft
             i0vidb = 8*i0sidb + 1
             i2vvtbl(i0vida,i0vidb) = 1
             
             !C P
             i0vidb = 8*i0sidb + 2
             i2vvtbl(i0vida,i0vidb) = 1

             !C Qt
             i0vidb = 8*i0sidb + 5
             i2vvtbl(i0vida,i0vidb) = 1
          
          ENDIF
          
          !C Fb
          i0vidb = 8*i0sidb
          i2vvtbl(i0vida,i0vidb) = 1

          !C Qb
          i0vidb = 8*i0sidb + 4
          i2vvtbl(i0vida,i0vidb) = 1

       ENDDO
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vida = 8*i0sida + 1
       
       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN
             
             !C Et
             i0vidb = 3
             i2vvtbl(i0vida,i0vidb) = 1

             
             !C N
             i0vidb = 8*i0sidb - 2
             i2vvtbl(i0vida,i0vidb) = 1

             !C Fr
             i0vidb = 8*i0sidb - 1
             i2vvtbl(i0vida,i0vidb) = 1

             !C Fb
             i0vidb = 8*i0sidb 
             i2vvtbl(i0vida,i0vidb) = 1

             !C P
             i0vidb = 8*i0sidb + 2
             i2vvtbl(i0vida,i0vidb) = 1

             !C Qb
             i0vidb = 8*i0sidb + 4
             i2vvtbl(i0vida,i0vidb) = 1
             
          ENDIF
          
          !C Ft
          i0vidb = 8*i0sidb + 1
          i2vvtbl(i0vida,i0vidb) = 1
          
          !C Qt
          i0vidb = 8*i0sidb + 5
          i2vvtbl(i0vida,i0vidb) = 1
          
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C N
       i0vidb = 8*i0sida - 2
       i2vvtbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2vvtbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2vvtbl(i0vida,i0vidb) = 1

       !C P
       i0vidb = 8*i0sida + 2
       i2vvtbl(i0vida,i0vidb) = 1

       !C Qb
       i0vidb = 8*i0sida + 4
       i2vvtbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2vvtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vida = 8*i0sida + 3
       
       !C Ep
       i0vidb = 4
       i2vvtbl(i0vida,i0vidb) = 1

       !C Er
       i0vidb = 5
       i2vvtbl(i0vida,i0vidb) = 1

       !C N
       i0vidb = 8*i0sida - 2
       i2vvtbl(i0vida,i0vidb) = 1

       !C P
       i0vidb = 8*i0sida + 2
       i2vvtbl(i0vida,i0vidb) = 1

       !C Qr
       i0vidb = 8*i0sida + 3
       i2vvtbl(i0vida,i0vidb) = 1

       !C Qb
       i0vidb = 8*i0sida + 4
       i2vvtbl(i0vida,i0vidb) = 1

       !C Qt
       i0vidb = 8*i0sida + 5
       i2vvtbl(i0vida,i0vidb) = 1

       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN

             !C Et
             i0vidb = 3
             i2vvtbl(i0vida,i0vidb) = 1

             !C Ep
             i0vidb = 4
             i2vvtbl(i0vida,i0vidb) = 1             

             !C N
             i0vidb = 8*i0sidb - 2
             i2vvtbl(i0vida,i0vidb) = 1

             !C Ft
             i0vidb = 8*i0sidb + 1
             i2vvtbl(i0vida,i0vidb) = 1

             !C P
             i0vidb = 8*i0sidb + 2
             i2vvtbl(i0vida,i0vidb) = 1

             !C Qt
             i0vidb = 8*i0sidb + 5
             i2vvtbl(i0vida,i0vidb) = 1

          ENDIF
          
          !C Fb
          i0vidb = 8*i0sidb 
          i2vvtbl(i0vida,i0vidb) = 1
          
          !C Qb
          i0vidb = 8*i0sidb + 4
          i2vvtbl(i0vida,i0vidb) = 1

       ENDDO

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vida = 8*i0sida + 5

       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN
             
             !C Et
             i0vidb = 3
             i2vvtbl(i0vida,i0vidb) = 1

             !C Ep
             i0vidb = 4
             i2vvtbl(i0vida,i0vidb) = 1
             
             !C N
             i0vidb = 8*i0sidb - 2
             i2vvtbl(i0vida,i0vidb) = 1
             
             !C Fb
             i0vidb = 8*i0sidb 
             i2vvtbl(i0vida,i0vidb) = 1

             !C P
             i0vidb = 8*i0sidb + 2
             i2vvtbl(i0vida,i0vidb) = 1

             !C Qr
             i0vidb = 8*i0sidb + 3
             i2vvtbl(i0vida,i0vidb) = 1

             !C Qb
             i0vidb = 8*i0sidb + 4
             i2vvtbl(i0vida,i0vidb) = 1


          ENDIF
          
          !C Ft
          i0vidb = 8*i0sidb + 1
          i2vvtbl(i0vida,i0vidb) = 1
          
          !C Qt
          i0vidb = 8*i0sidb + 5
          i2vvtbl(i0vida,i0vidb) = 1              
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_VV
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR MASS SCALAR ARRAY
  !C 
  !C                     2014-01-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_MS
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2mstbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C VARIABLE-VARIABLE GRAPH

    !C INITIALIZE

    i2mstbl(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    
    i0vida = 1
    
    !C PSI'   
    i0vidb = 1
    i2mstbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR I
    !C
    
    i0vida = 2
    
    !C I
    i0vidb = 2
    i2mstbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vida = 3

    !C Et
    i0vidb = 3
    i2mstbl(i0vida,i0vidb) = 1
        
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vida = 4
    
    !C Ep
    i0vidb = 4
    i2mstbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR Er
    !C

    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs
       
       !C
       !C EQUATION FOR N
       !C

       i0vida = 8*i0sida - 2
              
       !C N
       i0vidb = 8*i0sida - 2
       i2mstbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       
       !C Fb
       i0vidb = 8*i0sida
       i2mstbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vida = 8*i0sida + 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2mstbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C P
       i0vidb = 8*i0sida + 2
       i2mstbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2mstbl(i0vida,i0vidb) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vida = 8*i0sida + 5

       !C Qt
       i0vidb = 8*i0sida + 5
       i2mstbl(i0vida,i0vidb) = 1              
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_MS

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR ADVECTION VECTOR ARRAY
  !C 
  !C                     2014-01-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AV
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2avtbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    

    !C INITIALIZE
    
    i2avtbl(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    
    i0vida = 1
    
    !C PSI'   
    i0vidb = 1
    i2avtbl(i0vida,i0vidb) = 1
    
    !C
    !C EQUATION FOR Bt
    !C
    
    i0vida = 2
    
    !C I
    i0vidb = 2
    i2avtbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vida = 3

    !C PSI'
    i0vidb = 1
    i2avtbl(i0vida,i0vidb) = 1
    
    !C Et
    i0vidb = 3
    i2avtbl(i0vida,i0vidb) = 1
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vida = 4
        
    !C Ep
    i0vidb = 4
    i2avtbl(i0vida,i0vidb) = 1
    
    !C
    !C EQUATION FOR Er
    !C

    i0vida = 5
    
    !C Ep
    i0vidb = 4
    i2avtbl(i0vida,i0vidb) = 1

    !C Er
    i0vidb = 5
    i2avtbl(i0vida,i0vidb) = 1

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C

       i0vida = 8*i0sida - 2
              
       !C N
       i0vidb = 8*i0sida - 2
       i2avtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       !C Fb
       i0vidb = 8*i0sida
       i2avtbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2avtbl(i0vida,i0vidb) = 1
                     
       !C
       !C EQUATION FOR Ft
       !C

       i0vida = 8*i0sida + 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2avtbl(i0vida,i0vidb) = 1
   
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C P
       i0vidb = 8*i0sida + 2
       i2avtbl(i0vida,i0vidb) = 1

       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4

       !C Fb
       i0vidb = 8*i0sida 
       i2avtbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2avtbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2avtbl(i0vida,i0vidb) = 1
       
       !C EQUATION FOR Qt
       
       i0vida = 8*i0sida + 5
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2avtbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2avtbl(i0vida,i0vidb) = 1              
       
    ENDDO
    
    RETURN

  END SUBROUTINE T2VGRA_AV

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR ADVECTION TENSOR ARRAY
  !C 
  !C                     2014-01-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_AT
  
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2attbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2attbl(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    
    !C
    !C EQUATION FOR I
    !C
    
    !C
    !C EQUATION FOR Et
    !C
    
    !C
    !C EQUATION FOR Ep
    !C

    !C
    !C EQUATION FOR Er
    !C

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C

       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       !C Fb
       i0vidb = 8*i0sida
       i2attbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2attbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vida = 8*i0sida + 1
       
       !C Fb
       i0vidb = 8*i0sida
       i2attbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2attbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C Fb
       i0vidb = 8*i0sida 
       i2attbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2attbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       !C Fb
       i0vidb = 8*i0sida 
       i2attbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2attbl(i0vida,i0vidb) = 1

       !C Qb
       i0vidb = 8*i0sida + 5
       i2attbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vida = 8*i0sida + 5
       
       !C Fb
       i0vidb = 8*i0sida 
       i2attbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2attbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2attbl(i0vida,i0vidb) = 1              
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_AT
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR DIFFUSION TENSOR ARRAY
  !C 
  !C                     2014-01-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_DT
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2dttbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2dttbl(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C

    !C
    !C EQUATION FOR I
    !C

    !C
    !C EQUATION FOR Et
    !C
        
    !C
    !C EQUATION FOR Ep
    !C
    
    !C
    !C EQUATION FOR Er
    !C
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C

       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
      
       !C N
       i0vidb = 8*i0sida - 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida
       i2dttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2dttbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vida = 8*i0sida + 1
       
       !C N
       i0vidb = 8*i0sida - 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2dttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2dttbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C N
       i0vidb = 8*i0sida - 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2dttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2dttbl(i0vida,i0vidb) = 1

       !C Qb
       i0vidb = 8*i0sida + 4
       i2dttbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       !C N
       i0vidb = 8*i0sida - 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2dttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2dttbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vida = 8*i0sida + 5
       
       !C N
       i0vidb = 8*i0sida - 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2dttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2dttbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2dttbl(i0vida,i0vidb) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_DT

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR GRADIENT VECOTR ARRAY
  !C 
  !C                     2014-12-05 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GV
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2gvtbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C VARIABLE-VARIABLE GRAPH
    
    !C
    !C INITIALIZATION
    !C

    i2gvtbl(1:i0vmax,1:i0vmax) = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    
    i0vida = 1
    
    !C Et
    i0vidb = 3
    i2gvtbl(i0vida,i0vidb) = 1
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vida = 2
    
    !C Ep
    i0vidb = 4
    i2gvtbl(i0vida,i0vidb) = 1

    !C Er 
    i0vidb = 5
    i2gvtbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vida = 4
    
    !C I
    i0vidb = 2
    i2gvtbl(i0vida,i0vidb) = 1
    
    !C
    !C EQUATION FOR Er
    !C

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C

       !C
       !C EQUATION FOR Fr
       !C
       
       i0vida = 8*i0sida - 1
       
       !C P 
       i0vidb = 8*i0sida + 2
       i2gvtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       !C
       !C EQUATION FOR Ft
       !C
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C P
       i0vidb = 8*i0sida + 2
       i2gvtbl(i0vida,i0vidb) = 1

       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vida = 8*i0sida + 3
       
       !C N
       i0vidb = 8*i0sida - 2
       i2gvtbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2gvtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       !C
       !C EQUATION FOR Qt
       !C
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GV

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR GRADIENT TENSOR ARRAY
  !C 
  !C                     2014-01-27
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_GT
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2gttbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2gttbl(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    
    !C
    !C EQUATION FOR I
    !C

    !C
    !C EQUATION FOR Et
    !C
    
    !C
    !C EQUATION FOR Ep
    !C
    
    !C
    !C EQUATION FOR Er
    !C

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       !C N
       i0vidb = 8*i0sida - 2
       i2gttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida
       i2gttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2gttbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2gttbl(i0vida,i0vidb) = 1
          
       !C
       !C EQUATION FOR Ft
       !C

       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C N
       i0vidb = 8*i0sida - 2
       i2gttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2gttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2gttbl(i0vida,i0vidb) = 1

       !C Qb
       i0vidb = 8*i0sida + 4
       i2gttbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       !C N
       i0vidb = 8*i0sida - 2
       i2gttbl(i0vida,i0vidb) = 1
       
       !C Fb
       i0vidb = 8*i0sida 
       i2gttbl(i0vida,i0vidb) = 1
       
       !C P
       i0vidb = 8*i0sida + 2
       i2gttbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2gttbl(i0vida,i0vidb) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_GT

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION SCALAR ARRAY
  !C 
  !C                     2014-01-27
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2VGRA_ES
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2estbl
    
    INTEGER(i0ikind)::&
         i0sida,i0sidb,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2estbl(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    
    !C
    !C EQUATION FOR I
    !C
    
    !C
    !C EQUATION FOR Et
    !C
    
    i0vida = 3
    
    DO i0sidb = 1, i0spcs
       
       !C Ft
       i0vidb = 8*i0sidb + 1
       i2estbl(i0vida,i0vidb) = 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vida = 4

    DO i0sidb = 1, i0spcs
    
       !C Fb
       i0vidb = 8*i0sidb
       i2estbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sidb + 1
       i2estbl(i0vida,i0vidb) = 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C

    i0vida = 5
    
    DO i0sidb = 1, i0spcs
       
       !C N
       i0vidb = 8*i0sidb - 2
       i2estbl(i0vida,i0vidb) = 1

    ENDDO
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vida= 8*i0sida - 1
       
       !C Ep
       i0vidb = 4
       i2estbl(i0vida,i0vidb) = 1

       !C Er
       i0vidb = 5
       i2estbl(i0vida,i0vidb) = 1

       !C Fb
       i0vidb = 8*i0sida 
       i2estbl(i0vida,i0vidb) = 1

       !C Ft
       i0vidb = 8*i0sida + 1
       i2estbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Fb
       !C
       
       i0vida = 8*i0sida
       
       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN
             
             !C Et
             i0vidb = 3
             i2estbl(i0vida,i0vidb) = 1
             
             !C Ep
             i0vidb = 4
             i2estbl(i0vida,i0vidb) = 1
             
          ENDIF
          
          !C Fb
          i0vidb = 8*i0sidb
          i2estbl(i0vida,i0vidb) = 1
          
          !C Qb
          i0vidb = 8*i0sidb + 4
          i2estbl(i0vida,i0vidb) = 1
          
       ENDDO
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vida = 8*i0sida + 1
       
       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN
             
             !C Et
             i0vidb = 3
             i2estbl(i0vida,i0vidb) = 1

             !C Fr
             i0vidb = 8*i0sidb - 1
             i2estbl(i0vida,i0vidb) = 1

          ENDIF
          
          !C Ft
          i0vidb = 8*i0sidb + 1
          i2estbl(i0vida,i0vidb) = 1
          
          !C Qt
          i0vidb = 8*i0sidb + 5
          i2estbl(i0vida,i0vidb) = 1
          
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C P
       i0vidb = 8*i0sida + 2
       i2estbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vida = 8*i0sida + 3
       
       !C Ep
       i0vidb = 4
       i2estbl(i0vida,i0vidb) = 1

       !C Er
       i0vidb = 5
       i2estbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2estbl(i0vida,i0vidb) = 1

       !C Qt
       i0vidb = 8*i0sida + 5
       i2estbl(i0vida,i0vidb) = 1

       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN

             !C Et
             i0vidb = 3
             i2estbl(i0vida,i0vidb) = 1

             !C Ep
             i0vidb = 4
             i2estbl(i0vida,i0vidb) = 1             

          ENDIF
          
          !C Fb
          i0vidb = 8*i0sidb 
          i2estbl(i0vida,i0vidb) = 1
          
          !C Qb
          i0vidb = 8*i0sidb + 4
          i2estbl(i0vida,i0vidb) = 1

       ENDDO

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vida = 8*i0sida + 5

       DO i0sidb = 1, i0spcs
          
          IF(i0sida.EQ.i0sidb)THEN
             
             !C Et
             i0vidb = 3
             i2estbl(i0vida,i0vidb) = 1
             
             !C Qr
             i0vidb = 8*i0sidb + 3
             i2estbl(i0vida,i0vidb) = 1
             
          ENDIF
          
          !C Ft
          i0vidb = 8*i0sidb + 1
          i2estbl(i0vida,i0vidb) = 1
          
          !C Qt
          i0vidb = 8*i0sidb + 5
          i2estbl(i0vida,i0vidb) = 1              
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_ES

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION VECTOR ARRAY
  !C 
  !C                     2014-01-27
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_EV
        
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2evtbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2evtbl(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    
    !C
    !C EQUATION FOR I
    !C

    !C
    !C EQUATION FOR Et
    !C
    
    i0vida = 3

    !C PSI'
    i0vidb = 1
    i2evtbl(i0vida,i0vidb) = 1

    !C
    !C EQUATION FOR Ep
    !C
    
    !C
    !C EQUATION FOR Er
    !C

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C

       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       !C Fb
       i0vidb = 8*i0sida
       i2evtbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Ft
       !C
       
       !C
       !C EQUATION FOR P
       !C
       
       !C
       !C EQUATION FOR Qr
       !C

       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       
       !C Et
       i0vidb = 3
       i2evtbl(i0vida,i0vidb) = 1
       
       !C Ep
       i0vidb = 4
       i2evtbl(i0vida,i0vidb) = 1             
       
       !C Fb
       i0vidb = 8*i0sida
       i2evtbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2evtbl(i0vida,i0vidb) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vida = 8*i0sida + 5
       
       !C Et
       i0vidb = 3
       i2evtbl(i0vida,i0vidb) = 1
       
       !C Ep
       i0vidb = 4
       i2evtbl(i0vida,i0vidb) = 1
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2VGRA_EV

  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR EXCITATION TENSOR ARRAY
  !C 
  !C                     2014-01-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_ET

    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i2ettbl
    
    INTEGER(i0ikind)::&
         i0sida,i0vida,i0vidb
    
    !C
    !C INITIALIZATION
    !C
    
    i2ettbl(1:i0vmax,1:i0vmax) = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    
    !C
    !C EQUATION FOR I
    !C
    
    !C
    !C EQUATION FOR Et
    !C
    
    !C
    !C EQUATION FOR Ep
    !C
    
    !C
    !C EQUATION FOR Er
    !C

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sida = 1, i0spcs

       !C
       !C EQUATION FOR N
       !C

       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vida = 8*i0sida
       
       !C Fb
       i0vidb = 8*i0sida
       i2ettbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2ettbl(i0vida,i0vidb) = 1
       
       !C Qb
       i0vidb = 8*i0sida + 4
       i2ettbl(i0vida,i0vidb) = 1

       !C Qt
       i0vidb = 8*i0sida + 5
       i2ettbl(i0vida,i0vidb) = 1
              
       !C
       !C EQUATION FOR Ft
       !C
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vida = 8*i0sida + 2
       
       !C Fb
       i0vidb = 8*i0sida 
       i2ettbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2ettbl(i0vida,i0vidb) = 1

       !C Qb
       i0vidb = 8*i0sida + 4
       i2ettbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2ettbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vida = 8*i0sida + 4
       
       !C Fb
       i0vidb = 8*i0sida 
       i2ettbl(i0vida,i0vidb) = 1
       
       !C Ft
       i0vidb = 8*i0sida + 1
       i2ettbl(i0vida,i0vidb) = 1
          
       !C Qb
       i0vidb = 8*i0sida + 4
       i2ettbl(i0vida,i0vidb) = 1
       
       !C Qt
       i0vidb = 8*i0sida + 5
       i2ettbl(i0vida,i0vidb) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
    ENDDO
    
    RETURN 

  END SUBROUTINE T2VGRA_ET
  

  SUBROUTINE T2_VGRA_OUTPUT

    USE T2COMM
    IMPLICIT NONE
    INTEGER(i0ikind):: i0sida,i0sidb

    OPEN(10,FILE='TEST_VV.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'VV=',i2vvtbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_MS.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'MS=',i2mstbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AV.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'AV=',i2avtbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AT.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'AT=',i2attbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_DT.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'DT=',i2dttbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GV.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'GV=',i2gvtbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GT.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'GT=',i2gttbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ES.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'ES=',i2estbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_EV.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'EV=',i2evtbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ET.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'ET=',i2ettbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_SS.dat')
    DO i0sida = 1, i0vmax
    DO i0sidb = 1, i0vmax
       WRITE(10,*)'i=',i0sida,'j=',i0sidb,'SS=',i2sstbl(i0sida,i0sidb)
    ENDDO
    ENDDO
    CLOSE(10)

    RETURN
  END SUBROUTINE T2_VGRA_OUTPUT
    
END MODULE T2VGRA
