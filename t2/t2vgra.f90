!C------------------------------------------------------------------
!C
!C         MODULE T2VGRA
!C       
!C         VARIABLE GRAPH GENERATOR FOR TASK/T2
!C
!C                       2014-01-28 H.SETO
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
  !C          MODIFIED 2014-01-27
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_VV
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2vvvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj
    
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
    !C EQUATION FOR PSI
    !C
    
    i0vidi = 1
    
    !C PSI'   
    i0vidj = 1
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR I
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
    !C EQUATION FOR Et
    !C
    
    i0vidi = 3

    !C PSI'
    i0vidj = 1
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2vvvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax

       !C Ft
       i0vidj = 8*i0sidj + 1
       i2vvvt(i0vidi,i0vidj) = 1

    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    i2vvvt(i0vidi,i0vidj) = 1
    
    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1


    DO i0sidj = 1, i0smax
    
       !C Fb
       i0vidj = 8*i0sidj
       i2vvvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = 8*i0sidj + 1
       i2vvvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C

    i0vidi = 5
    
    !C Ep
    i0vidj = 4
    i2vvvt(i0vidi,i0vidj) = 1

    !C Er
    i0vidj = 5
    i2vvvt(i0vidi,i0vidj) = 1
    
    DO i0sidj = 1, i0smax
       
       !C N
       i0vidj = 8*i0sidj - 2
       i2vvvt(i0vidi,i0vidj) = 1

    ENDDO

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C

       i0vidi = 8*i0sidi - 2
              
       !C N
       i0vidj = 8*i0sidi - 2
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi= 8*i0sidi - 1
       
       !C Ep
       i0vidj = 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fr
       i0vidj = 8*i0sidi - 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = 8*i0sidi 
       i2vvvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = 8*i0sidi + 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C P 
       i0vidj = 8*i0sidi + 2
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = 8*i0sidj - 2
             i2vvvt(i0vidi,i0vidj) = 1

             
             !C Ft
             i0vidj = 8*i0sidj + 1
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C P
             i0vidj = 8*i0sidj + 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qt
             i0vidj = 8*i0sidj + 5
             i2vvvt(i0vidi,i0vidj) = 1
          
          ENDIF
          
          !C Fb
          i0vidj = 8*i0sidj
          i2vvvt(i0vidi,i0vidj) = 1

          !C Qb
          i0vidj = 8*i0sidj + 4
          i2vvvt(i0vidi,i0vidj) = 1

       ENDDO
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = 8*i0sidi + 1
       
       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             
             !C N
             i0vidj = 8*i0sidj - 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fr
             i0vidj = 8*i0sidj - 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C Fb
             i0vidj = 8*i0sidj 
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = 8*i0sidj + 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qb
             i0vidj = 8*i0sidj + 4
             i2vvvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Ft
          i0vidj = 8*i0sidj + 1
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = 8*i0sidj + 5
          i2vvvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = 8*i0sidi 
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Ft
       i0vidj = 8*i0sidi + 1
       i2vvvt(i0vidi,i0vidj) = 1

       !C P
       i0vidj = 8*i0sidi + 2
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = 8*i0sidi + 4
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = 8*i0sidi + 5
       i2vvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 8*i0sidi + 3
       
       !C Ep
       i0vidj = 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C N
       i0vidj = 8*i0sidi - 2
       i2vvvt(i0vidi,i0vidj) = 1

       !C P
       i0vidj = 8*i0sidi + 2
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qr
       i0vidj = 8*i0sidi + 3
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = 8*i0sidi + 4
       i2vvvt(i0vidi,i0vidj) = 1

       !C Qt
       i0vidj = 8*i0sidi + 5
       i2vvvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN

             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1             

             !C N
             i0vidj = 8*i0sidj - 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ft
             i0vidj = 8*i0sidj + 1
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = 8*i0sidj + 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qt
             i0vidj = 8*i0sidj + 5
             i2vvvt(i0vidi,i0vidj) = 1

          ENDIF
          
          !C Fb
          i0vidj = 8*i0sidj 
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = 8*i0sidj + 4
          i2vvvt(i0vidi,i0vidj) = 1

       ENDDO

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5

       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C N
             i0vidj = 8*i0sidj - 2
             i2vvvt(i0vidi,i0vidj) = 1
             
             !C Fb
             i0vidj = 8*i0sidj 
             i2vvvt(i0vidi,i0vidj) = 1

             !C P
             i0vidj = 8*i0sidj + 2
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qr
             i0vidj = 8*i0sidj + 3
             i2vvvt(i0vidi,i0vidj) = 1

             !C Qb
             i0vidj = 8*i0sidj + 4
             i2vvvt(i0vidi,i0vidj) = 1


          ENDIF
          
          !C Ft
          i0vidj = 8*i0sidj + 1
          i2vvvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = 8*i0sidj + 5
          i2vvvt(i0vidi,i0vidj) = 1              
          
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
         i0smax,i0vmax,i2msvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
                i0vidj
    
    !C VARIABLE-VARIABLE GRAPH

    !C INITIALIZE

    i2msvt(1:i0vmax,1:i0vmax) = 0

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
    i2msvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR I
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2msvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vidi = 3

    !C Et
    i0vidj = 3
    i2msvt(i0vidi,i0vidj) = 1
        
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    !C Ep
    i0vidj = 4
    i2msvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Er
    !C

    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax
       
       !C
       !C EQUATION FOR N
       !C

       i0vidi = 8*i0sidi - 2
              
       !C N
       i0vidj = 8*i0sidi - 2
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       
       !C Fb
       i0vidj = 8*i0sidi
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = 8*i0sidi + 1
       
       !C Ft
       i0vidj = 8*i0sidi + 1
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2msvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2msvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5

       !C Qt
       i0vidj = 8*i0sidi + 5
       i2msvt(i0vidi,i0vidj) = 1              
       
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
         i0smax,i0vmax,i2avvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
                i0vidj
    
    

    !C INITIALIZE
    
    i2avvt(1:i0vmax,1:i0vmax) = 0

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
    i2avvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR Bt
    !C
    
    i0vidi = 2
    
    !C I
    i0vidj = 2
    i2avvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    i0vidi = 3

    !C PSI'
    i0vidj = 1
    i2avvt(i0vidi,i0vidj) = 1
    
    !C Et
    i0vidj = 3
    i2avvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
        
    !C Ep
    i0vidj = 4
    i2avvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR Er
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

       !C
       !C EQUATION FOR N
       !C

       i0vidi = 8*i0sidi - 2
              
       !C N
       i0vidj = 8*i0sidi - 2
       i2avvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb
       i0vidj = 8*i0sidi
       i2avvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2avvt(i0vidi,i0vidj) = 1
                     
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = 8*i0sidi + 1
       
       !C Ft
       i0vidj = 8*i0sidi + 1
       i2avvt(i0vidi,i0vidj) = 1
   
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2avvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4

       !C Fb
       i0vidj = 8*i0sidi 
       i2avvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2avvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2avvt(i0vidi,i0vidj) = 1
       
       !C EQUATION FOR Qt
       
       i0vidi = 8*i0sidi + 5
       
       !C Ft
       i0vidj = 8*i0sidi + 1
       i2avvt(i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = 8*i0sidi + 5
       i2avvt(i0vidi,i0vidj) = 1              
       
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
         i0smax,i0wmax,i0vmax,i2atvt,i3atwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,&
                       i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2atvt(         1:i0vmax,1:i0vmax) = 0
    i3atwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
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
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C

       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb (B)
       i0vidj = 8*i0sidi
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Ft (B)
       i0vidj = 8*i0sidi + 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qb (B)
       i0vidj = 8*i0sidi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qt (B)
       i0vidj = 8*i0sidi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vidi = 8*i0sidi + 1
       
       !C Fb (B)
       i0vidj = 8*i0sidi
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Ft (B)
       i0vidj = 8*i0sidi + 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qb (B)
       i0vidj = 8*i0sidi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qt (B)
       i0vidj = 8*i0sidi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C Fb (B)
       i0vidj = 8*i0sidi 
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Ft (B)
       i0vidj = 8*i0sidi + 1
       i2atvt(i0vidi,i0vidj) = 1 
       
       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb (B)
       i0vidj = 8*i0sidi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qt (B)
       i0vidj = 8*i0sidi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C Fb (B)
       i0vidj = 8*i0sidi 
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Ft (B)
       i0vidj = 8*i0sidi + 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb (B)
       i0vidj = 8*i0sidi + 4
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Qb (B)
       i0vidj = 8*i0sidi + 5
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5
       
       !C Fb
       i0vidj = 8*i0sidi 
       i2atvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = 8*i0sidi + 1
       i2atvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2atvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3atwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qt
       i0vidj = 8*i0sidi + 5
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
  !C                     2014-01-27 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2VGRA_DT
    
    USE T2COMM,ONLY:&
         i0smax,i0vmax,i2dtvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
                i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2dtvt(1:i0vmax,1:i0vmax) = 0

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
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C

       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
      
       !C N
       i0vidj = 8*i0sidi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = 8*i0sidi
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vidi = 8*i0sidi + 1
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = 8*i0sidi 
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = 8*i0sidi 
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2dtvt(i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = 8*i0sidi + 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = 8*i0sidi 
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vidi = 8*i0sidi + 5
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Fb
       i0vidj = 8*i0sidi 
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2dtvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2dtvt(i0vidi,i0vidj) = 1
       
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
         i0smax,i0vmax,i2gvvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
                i0vidj
    
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
    !C EQUATION FOR PSI
    !C
    
    i0vidi = 1
    
    !C Et
    i0vidj = 3
    i2gvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vidi = 2
    
    !C Ep
    i0vidj = 4
    i2gvvt(i0vidi,i0vidj) = 1

    !C Er 
    i0vidj = 5
    i2gvvt(i0vidi,i0vidj) = 1

    !C
    !C EQUATION FOR Et
    !C
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4
    
    !C I
    i0vidj = 2
    i2gvvt(i0vidi,i0vidj) = 1
    
    !C
    !C EQUATION FOR Er
    !C

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C

       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi = 8*i0sidi - 1
       
       !C P 
       i0vidj = 8*i0sidi + 2
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       !C
       !C EQUATION FOR Ft
       !C
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2gvvt(i0vidi,i0vidj) = 1

       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 8*i0sidi + 3
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2gvvt(i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2gvvt(i0vidi,i0vidj) = 1
       
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
         i0smax,i0wmax,i0vmax,i2gtvt,i3gtwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,&
                       i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2gtvt(         1:i0vmax,1:i0vmax) = 0
    i3gtwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0

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
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C N (B)
       i0vidj = 8*i0sidi - 2
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = 8*i0sidi
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C P
       i0vidj = 8*i0sidi + 2
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
          
       !C
       !C EQUATION FOR Ft
       !C
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C N (B), (Ub)
       i0vidj = 8*i0sidi - 2
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = 8*i0sidi 
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C P
       i0vidj = 8*i0sidi + 2
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Qb
       i0vidj = 8*i0sidi + 4
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C N
       i0vidj = 8*i0sidi - 2
       i2gtvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = 8*i0sidi 
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2gtvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3gtwt(i0widi,i0vidi,i0vidj) = 1

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
         i0smax,i0vmax,i2esvt
    
    INTEGER(i0ikind)::&
         i0sidi,i0vidi,&
         i0sidj,i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2esvt(1:i0vmax,1:i0vmax) = 0
    
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
    
    i0vidi = 3
    
    DO i0sidj = 1, i0smax
       
       !C Ft
       i0vidj = 8*i0sidj + 1
       i2esvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vidi = 4

    DO i0sidj = 1, i0smax
    
       !C Fb
       i0vidj = 8*i0sidj
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Ft
       i0vidj = 8*i0sidj + 1
       i2esvt(i0vidi,i0vidj) = 1
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C

    i0vidi = 5
    
    DO i0sidj = 1, i0smax
       
       !C N
       i0vidj = 8*i0sidj - 2
       i2esvt(i0vidi,i0vidj) = 1

    ENDDO
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi= 8*i0sidi - 1
       
       !C Ep
       i0vidj = 4
       i2esvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2esvt(i0vidi,i0vidj) = 1

       !C Fb
       i0vidj = 8*i0sidi 
       i2esvt(i0vidi,i0vidj) = 1

       !C Ft
       i0vidj = 8*i0sidi + 1
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C
       
       i0vidi = 8*i0sidi
       
       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Ep
             i0vidj = 4
             i2esvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Fb
          i0vidj = 8*i0sidj
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = 8*i0sidj + 4
          i2esvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQUATION FOR Ft
       !C
       
       i0vidi = 8*i0sidi + 1
       
       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1

             !C Fr
             i0vidj = 8*i0sidj - 1
             i2esvt(i0vidi,i0vidj) = 1

          ENDIF
          
          !C Ft
          i0vidj = 8*i0sidj + 1
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = 8*i0sidj + 5
          i2esvt(i0vidi,i0vidj) = 1
          
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2esvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 8*i0sidi + 3
       
       !C Ep
       i0vidj = 4
       i2esvt(i0vidi,i0vidj) = 1

       !C Er
       i0vidj = 5
       i2esvt(i0vidi,i0vidj) = 1
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2esvt(i0vidi,i0vidj) = 1

       !C Qt
       i0vidj = 8*i0sidi + 5
       i2esvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN

             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1

             !C Ep
             i0vidj = 4
             i2esvt(i0vidi,i0vidj) = 1             

          ENDIF
          
          !C Fb
          i0vidj = 8*i0sidj 
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qb
          i0vidj = 8*i0sidj + 4
          i2esvt(i0vidi,i0vidj) = 1

       ENDDO

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5

       DO i0sidj = 1, i0smax
          
          IF(i0sidi.EQ.i0sidj)THEN
             
             !C Et
             i0vidj = 3
             i2esvt(i0vidi,i0vidj) = 1
             
             !C Qr
             i0vidj = 8*i0sidj + 3
             i2esvt(i0vidi,i0vidj) = 1
             
          ENDIF
          
          !C Ft
          i0vidj = 8*i0sidj + 1
          i2esvt(i0vidi,i0vidj) = 1
          
          !C Qt
          i0vidj = 8*i0sidj + 5
          i2esvt(i0vidi,i0vidj) = 1              
          
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
         i0smax,i0wmax,i0vmax,i2evvt,i3evwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,&
                       i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2evvt(         1:i0vmax,1:i0vmax) = 0
    i3evwt(1:i0wmax,1:i0vmax,1:i0vmax) = 0
    
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
    
    i0vidi = 3

    !C PSI' (R)
    i0vidj = 1
    i2evvt(i0vidi,i0vidj) = 1

    i0widi = 2
    i3evwt(i0widi,i0vidi,i0vidj) = 1

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
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C
       
       !C
       !C EQUATION FOR Fr
       !C

       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb (B)
       i0vidj = 8*i0sidi
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

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
       
       i0vidi = 8*i0sidi + 4
       
       
       !C Et (B) (Ub) (Qb/P)
       i0vidj = 3
       i2evvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Ep (B) (Ub) (Qb/P)
       i0vidj = 4
       i2evvt(i0vidi,i0vidj) = 1             
       
       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 2
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C Fb (B)
       i0vidj = 8*i0sidi
       i2evvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       
       !C Qb (B)
       i0vidj = 8*i0sidi + 4
       i2evvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i3evwt(i0widi,i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5
       
       !C Et
       i0vidj = 3
       i2evvt(i0vidi,i0vidj) = 1
       
       !C Ep
       i0vidj = 4
       i2evvt(i0vidi,i0vidj) = 1
       
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
         i0smax,i0wmax,i0vmax,i2etvt,i4etwt
    
    INTEGER(i0ikind)::&
         i0sidi,i0widi,i0vidi,&
                i0widj,i0vidj
    
    !C
    !C INITIALIZATION
    !C
    
    i2etvt(                  1:i0vmax,1:i0vmax) = 0
    i4etwt(1:i0wmax,1:i0wmax,1:i0vmax,1:i0vmax) = 0

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
    
    DO i0sidi = 1, i0smax

       !C
       !C EQUATION FOR N
       !C

       !C
       !C EQUATION FOR Fr
       !C
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb (B,B)
       i0vidj = 8*i0sidi
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Ft (B,B)
       i0vidj = 8*i0sidi + 1
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qb (B,B)
       i0vidj = 8*i0sidi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qt (B,B)
       i0vidj = 8*i0sidi + 5
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Ft
       !C
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C Fb (B,B), (Ub,B)
       i0vidj = 8*i0sidi 
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Ft (B,B), (Ub,B)
       i0vidj = 8*i0sidi + 1
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qb (B,B), (Ub,B)
       i0vidj = 8*i0sidi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Qt (B,B), (Ub,B)
       i0vidj = 8*i0sidi + 5
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       i0widi = 2*i0sidi + 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C Fb (B,B)
       i0vidj = 8*i0sidi 
       i2etvt(i0vidi,i0vidj) = 1
       
       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Ft (B,B)
       i0vidj = 8*i0sidi + 1
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C Qb (B,B)
       i0vidj = 8*i0sidi + 4
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1

       !C Qt (B,B)
       i0vidj = 8*i0sidi + 5
       i2etvt(i0vidi,i0vidj) = 1

       i0widi = 1
       i0widj = 1
       i4etwt(i0widi,i0widj,i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qt
       !C
       
    ENDDO
    
    RETURN 

  END SUBROUTINE T2VGRA_ET
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR SOURCE SCALAR ARRAY
  !C 
  !C                     2014-01-28 H.SETO
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

       i0vidi = 8*i0sidi - 2
              
       !C N
       i0vidj = 8*i0sidi - 2
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fr
       !C
       
       i0vidi = 8*i0sidi - 1
       
       !C Qb
       i0vidj = 8*i0sidi - 1
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Fb
       !C

       i0vidi = 8*i0sidi
       
       !C Fb
       i0vidj = 8*i0sidi
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Ft
       !C

       i0vidi = 8*i0sidi + 1
       
       !C Ft
       i0vidj = 8*i0sidi + 1
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR P
       !C
       
       i0vidi = 8*i0sidi + 2
       
       !C P
       i0vidj = 8*i0sidi + 2
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qr
       !C
       
       i0vidi = 8*i0sidi + 3
       
       !C Qr
       i0vidj = 8*i0sidi + 3
       i2ssvt(i0vidi,i0vidj) = 1
       
       !C
       !C EQUATION FOR Qb
       !C
       
       i0vidi = 8*i0sidi + 4
       
       !C Qb
       i0vidj = 8*i0sidi + 4
       i2ssvt(i0vidi,i0vidj) = 1

       !C
       !C EQUATION FOR Qt
       !C
       
       i0vidi = 8*i0sidi + 5

       !C Qt
       i0vidj = 8*i0sidi + 5
       i2ssvt(i0vidi,i0vidj) = 1              
       
    ENDDO
    
    RETURN

  END SUBROUTINE T2VGRA_SS

    


  SUBROUTINE T2_VGRA_OUTPUT

    USE T2COMM
    INTEGER(i0ikind):: i0sidi,i0sidj

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

    RETURN
  END SUBROUTINE T2_VGRA_OUTPUT
    
END MODULE T2VGRA
