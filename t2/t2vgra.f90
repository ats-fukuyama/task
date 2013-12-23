!C------------------------------------------------------------------
!C
!C         MODULE T2VGRA
!C       
!C         VARIABLE GRAPH GENERATOR FOR TASK/T2
!C
!C         MODIFIED 2013-12-05
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
    
    CALL T2VGRA_VM
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
  
  !C
  !C
  !C SUBROUTINE FOR VALIABLE MATRIX ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C
  SUBROUTINE T2VGRA_VM
    
    USE T2COMM,ONLY:&
         i0spcs,i0vmax,i0vgrmx,i0vgcmx,&
         i1vgidr,i1vgidc,i2vtbl
    
    INTEGER(i0ikind)::&
         i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vofs
    
    !C VARIABLE-VARIABLE GRAPH

    !C INITIALIZE
    i0vofs                       = 0
    i2vtbl(1:i0vmax,1:i0vmax)    = 0
    i1vgidr(1:i0vgrmx)           = 0
    i1vgidc(1:i0vgcmx)           = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR Psi
    !C
    i0vr                       = 1
    i1vgidr(i0vr)              = i0vofs+1
    
    !C PSI'
    i0vc                       = 1
    i0vg                       = i0vofs+1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    !C Et
    i0vc                       = 3
    i0vg                       = i0vofs+2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    i0vofs                     = i0vofs + 2
    
    !C
    !C EQUATION FOR Bt
    !C
    i0vr                       = 2
    i1vgidr(i0vr)              = i0vofs + 1
    
    !C I
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    !C Er 
    i0vc                       = 5
    i0vg                       = i0vofs + 3
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    i0vofs                     = i0vofs + 3

    !C
    !C EQUATION FOR Et
    !C
    i0vr                       = 3
    i1vgidr(i0vr)              = i0vofs + 1
    
    !C PSI'
    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    !C Et
    i0vc                       = 3
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    DO i1 = 1, i0spcs
       !C Ft
       i0vc                    = 8*i1 + 1
       i0vg                    = i0vofs + 2 + i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO
    
    i0vofs = i0vofs + 2 + i0spcs

    !C
    !C EQUATION FOR Ep
    !C
    i0vr                       = 4
    i1vgidr(i0vr)              = i0vofs + 1
    
    !C I
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    DO i1 = 1, i0spcs
       !C Fb
       i0vc                    = 8*i1
       i0vg                    = i0vofs + 1 + 2*i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Ft
       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + 2 + 2*i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO
    
    i0vofs = i0vofs + 2 + 2*i0spcs
    
    !C
    !C EQUATION FOR Er
    !C
    i0vr                       = 5
    i1vgidr(i0vr)              = i0vofs + 1
    
    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    !C Er
    i0vc                       = 5
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    DO i1 = 1, i0spcs
       !C N
       i0vc                    = 8*i1 - 2
       i0vg                    = i0vofs + 2 + i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO

    i0vofs = i0vofs + 2 + i0spcs

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    DO i1 = 1, i0spcs
       !C
       !C EQUATION FOR N
       !C
       i0vr                    = 8*i1   - 2
       i1vgidr(i0vr)           = i0vofs + 1
       
       !C N
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
       
       i0vofs                  = i0vofs + 1
       
       !C
       !C EQUATION FOR Fr
       !C
       i0vr                    = 8*i1   - 1
       i1vgidr(i0vr)           = i0vofs + 1
       
       !C Ep
       i0vc                    = 4
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Er
       i0vc                    = 5
       i0vg                    = i0vofs + 2
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Fr
       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Fb
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 4
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Ft
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 5
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C P
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 6
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vofs                  = i0vofs + 6
       
       !C
       !C EQUATION FOR Fb
       !C
       i0vr                    = 8*i1
       i1vgidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs
          
          i0vcx   = 8*j1
          
          IF(i1.EQ.j1)THEN
             
             !C Et
             i0vc              = 3
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Ep
             i0vc              = 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C N
             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Fb
             i0vc              = i0vcx
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Ft
             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C P
             i0vc              = i0vcx  + 2
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qb
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qt
             i0vc              = i0vcx  + 5
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vofs = i0vofs + 8
         
          ELSE
             
             !C Fb
             i0vc              = i0vcx
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qb
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vofs = i0vofs + 2
             
          END IF
       ENDDO
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr                    = 8*i1   + 1
       i1vgidr(i0vr)           = i0vofs + 1
       
       DO j1=1,i0spcs
          
          i0vcx                = 8*j1 + 1
          
          IF(i1.EQ.j1)THEN
             !C Et
             i0vc              = 3
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C N
             i0vc              = i0vcx  - 3
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Fr
             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Fb
             i0vc              = i0vcx  - 1
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Ft
             i0vc              = i0vcx  
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C P
             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Qb
             i0vc              = i0vcx  + 3
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Qt
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 8
             
          ELSE
             
             !C Ft
             i0vc              = i0vcx  
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qt
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vofs = i0vofs + 2
             
          END IF
       ENDDO
       
       !C
       !C EQUATION FOR P
       !C
       i0vr                    = 8*i1   + 2
       i1vgidr(i0vr)           = i0vofs + 1
       
       !C N
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Fb
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Ft
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C P
       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Qb
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 5
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Qt
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 6
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vofs                  = i0vofs + 6

       !C
       !C EQUATION FOR Qr
       !C
       i0vr                    = 8*i1   + 3
       i1vgidr(i0vr)           = i0vofs + 1
       
       !C Ep
       i0vc                    = 4
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
       
       !C Er
       i0vc                    = 5
       i0vg                    = i0vofs + 2
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
       
       !C N
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 3
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
       
       !C P
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 4
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Qr
       i0vc                    = i0vr  
       i0vg                    = i0vofs + 5
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Qb
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 6
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C Qt
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 7
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
       
       i0vofs                  = i0vofs + 7

       !C
       !C EQUATION FOR Qb
       !C
       i0vr                    = 8*i1   + 4
       i1vgidr(i0vr)           = i0vofs + 1
       
       DO j1=1,i0spcs
          
          i0vcx   = 8*j1+4
          
          IF(i1.EQ.j1)THEN
             !C Et
             i0vc              = 3
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Ep
             i0vc              = 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C N
             i0vc              = i0vcx  - 6
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Fb
             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Ft
             i0vc              = i0vcx  - 3
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C P
             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qb
             i0vc              = i0vcx
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qt
             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 8

          ELSE
             !C Fb
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Qb
             i0vc              = i0vcx
             i0vg              = i0vofs +  2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 2

          END IF
       ENDDO
       
       !C EQUATION FOR Qt
       
       i0vr                    = 8*i1   + 5
       i1vgidr(i0vr)           = i0vofs + 1
       
       DO j1 = 1, i0spcs
          
          i0vcx   = 8*j1+5
          
          IF(i1.EQ.j1)THEN
             !C Et
             i0vc              = 3
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Ep
             i0vc              = 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C N
             i0vc              = i0vcx  - 7
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Fb
             i0vc              = i0vcx  - 5
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Ft
             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C P
             i0vc              = i0vcx  - 3
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             !C Qr
             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Qb
             i0vc              = i0vcx  - 1
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Qt
             i0vc              = i0vcx
             i0vg              = i0vofs + 9
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 9

          ELSE
             !C Ft
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             !C Qt
             i0vc              = i0vcx
             i0vg              = i0vofs +  2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vofs = i0vofs + 2
             
          END IF
       ENDDO
    ENDDO

    i0vr          =  6 + 8*i0spcs
    i1vgidr(i0vr) =  i0vofs + 1
    
    RETURN
    
  END SUBROUTINE T2VGRA_VM
  
  !C
  !C
  !C SUBROUTINE FOR MASS SCALAR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C
  SUBROUTINE T2VGRA_MS
    
    USE T2COMM,ONLY:&
         i0msrmx,i0mscmx,i0spcs,i1msidr,i1msidc,i0dbg

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1msidr(1:i0msrmx)         = 0
    i1msidc(1:i0mscmx)         = 0
    i0vofs                     = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI'
    !C
    i0vr                       = 1
    i1msidr(i0vr)              = i0vofs + 1
    
    !C PSI' 
    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc
    
    i0vofs                     = i0vofs + 1
    
    !C
    !C EQUATION FOR I
    !C
    i0vr                       = 2
    i1msidr(i0vr)              = i0vofs + 1

    !C I
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C
    !C EQUATION FOR Et
    !C
    i0vr                       = 3
    i1msidr(i0vr)              = i0vofs + 1

    !C Et
    i0vc                       = 3
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc
    
    i0vofs = i0vofs + 1

    !C
    !C EQUATION FOR Ep
    !C
    i0vr                       = 4
    i1msidr(i0vr)              = i0vofs + 1
    
    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc

    i0vofs = i0vofs + 1

    !C
    !C EQUATION FOR Er
    !C
    i0vr                        = 5
    i1msidr(i0vr)               = i0vofs + 1

    !C Er
    i0vc                        = 5
    i0vg                        = i0vofs + 1
    i1msidc(i0vg)               = i0vc

    i0vofs                      = i0vofs + 1

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID
    !C
    !C
    DO i1 = 1,i0spcs
       !C
       !C EQUATION FOR N
       !C
       i0vr                    = 8*i1   - 2
       i1msidr(i0vr)           = i0vofs + 1

       !C N
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc
       
       i0vofs = i0vofs + 1

       !C
       !C EQUATION FOR Fr
       !C
       i0vr                    = 8*i1   - 1
       i1msidr(i0vr)           = i0vofs + 1
       
       !C Fr
       i0vc                       = i0vr
       i0vg                       = i0vofs + 1
       i1msidc(i0vg)              = i0vc
       
       i0vofs = i0vofs + 1
       
       !C
       !C EQUATION FOR Fb
       !C
       i0vr                    = 8*i1
       i1msidr(i0vr)           = i0vofs + 1

       !C Fb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR Ft
       !C 
       i0vr                    = 8*i1   + 1
       i1msidr(i0vr)           = i0vofs + 1
       
       !C Ft
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1
       
       !C
       !C EQUATION FOR P
       !C
       i0vr                    = 8*i1   + 2
       i1msidr(i0vr)           = i0vofs + 1

       !C P
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR Qr
       !C
       i0vr                    = 8*i1   + 3
       i1msidr(i0vr)           = i0vofs + 1

       !C Qr
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       i0vr                    = 8*i1   + 4
       i1msidr(i0vr)           = i0vofs + 1

       !C Qb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1msidr(i0vr)           = i0vofs + 1
       
       !C Qt
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1
    ENDDO
    
    i0vr                       = i0msrmx
    i1msidr(i0vr)              = i0vofs + 1
    
    
    RETURN

  END SUBROUTINE T2VGRA_MS

  !C
  !C
  !C SUBROUTINE FOR ADVECTION VECTOR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C
  SUBROUTINE T2VGRA_AV

    USE T2COMM,ONLY:&
         i0spcs, i0avrmx, i0avcmx, i1avidr, i1avidc

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1avidr(1:i0avrmx)         = 0
    i1avidc(1:i0avcmx)         = 0
    i0vofs                     = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    i0vr                       = 1
    i1avidr(i0vr)              = i0vofs + 1
    
    !C PSI
    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C
    !C EQUATION FOR I
    !C
    i0vr                       = 2
    i1avidr(i0vr)              = i0vofs + 1

    !C I
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1
    
    !C
    !C EQUATION FOR Et
    !C
    i0vr                       = 3
    i1avidr(i0vr)              = i0vofs + 1

    !C PSI
    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    !C Et
    i0vc                       = 3
    i0vg                       = i0vofs + 2
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 2

    !C EQUATION FOR Ep
    i0vr                       = 4
    i1avidr(i0vr)              = i0vofs + 1

    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C EQUATION FOR Er
    i0vr                       = 5
    i1avidr(i0vr)              = i0vofs + 1

    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    !C Er
    i0vc                       = 5
    i0vg                       = i0vofs + 2
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 2

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID
    !C
    !C
    DO i1 = 1,i0spcs
       !C
       !C EQUATION FOR N
       !C
       i0vr                    = 8*i1   - 2
       i1avidr(i0vr)           = i0vofs + 1

       !C N
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1
       
       !C
       !C EQUATION FOR Fr
       !C
       i0vr                    = 8*i1   - 1
       i1avidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR Fb
       !C
       i0vr                    = 8*i1
       i1avidr(i0vr)           = i0vofs + 1
       
       !C Fb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 2
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 2
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr                    = 8*i1   + 1
       i1avidr(i0vr)           = i0vofs + 1

       !C Ft
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR P
       !C
       i0vr                    = 8*i1   + 2
       i1avidr(i0vr)           = i0vofs + 1

       !C P
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR Qr
       !C
       i0vr                    = 8*i1   + 3
       i1avidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       i0vr                    = 8*i1   + 4
       i1avidr(i0vr)           = i0vofs + 1
       
       !C Fb
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1avidc(i0vg)           = i0vc

       !C Qb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1avidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 3
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1avidr(i0vr)           = i0vofs + 1
       
       !C Ft
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       !C Qt
       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 2
       
    ENDDO
    

    i0vr                       = i0avrmx
    i1avidr(i0vr)              = i0vofs + 1

    RETURN
    
  END SUBROUTINE T2VGRA_AV
  
  !C
  !C
  !C SUBROUTINE FOR ADVECTION TENSOR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C
  
  SUBROUTINE T2VGRA_AT
    
    USE T2COMM,ONLY:&
         i0spcs,i0atrmx,i0atcmx,i1atidr,i2atidc    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    
    !C INITIALIZE
    i1atidr(1:i0atrmx)     = 0
    i2atidc(1:i0atcmx,1:2) = 0
    i0vofs                 = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    i0vr                       = 1
    i1atidr(i0vr)              = i0vofs + 1

    !C
    !C EQUATION FOR I
    !C
    i0vr                       = 2
    i1atidr(i0vr)              = i0vofs + 1

    !C
    !C EQUATION FOR Et
    !C
    i0vr                       = 3
    i1atidr(i0vr)              = i0vofs + 1
    
    !C
    !C EQUATION FOR Ep
    !C
    i0vr                       = 4
    i1atidr(i0vr)              = i0vofs + 1

    !C
    !C EQUATION FOR Er
    !C
    i0vr                        = 5
    i1atidr(i0vr)               = i0vofs + 1
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    DO i1 = 1,i0spcs
       !C
       !C EQUATION FOR N
       !C
       i0vr                    = 8*i1   - 2
       i1atidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR Fr
       !C
       i0vr                    = 8*i1   - 1
       i1atidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR Fb
       !C
       i0vr                    = 8*i1
       i1atidr(i0vr)           = i0vofs + 1
       
       !C Fb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Ft
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 2
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Qb
       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 3
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       !C Qt
       i0vc                    = i0vr   + 5
       i0vg                    = i0vofs + 4
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr                    = 8*i1   + 1
       i1atidr(i0vr)           = i0vofs + 1
       
       !C Fb
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 1
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Ft
       i0vc                    = i0vr   
       i0vg                    = i0vofs + 2
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       !C Qb
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 3
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       !C Qt
       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 4
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR P
       !C
       i0vr                    = 8*i1   + 2
       i1atidr(i0vr)           = i0vofs + 1

       !C Fb
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       !C Ft
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Qb
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Qt
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 4
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       i0vofs                  = i0vofs + 4

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1atidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1atidr(i0vr)           = i0vofs + 1
       
       !C Fb
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Ft
       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 2
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       !C Qb
       i0vc                    = i0vr   
       i0vg                    = i0vofs + 3
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Qt
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 4
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1atidr(i0vr)           = i0vofs + 1
       
       !C Fb
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 1
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Ft
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1
       
       !C Qb
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       !C Qt
       i0vc                    = i0vr   
       i0vg                    = i0vofs + 4
       i2atidc(i0vg,1)         = i0vc
       i2atidc(i0vg,2)         = 1

       i0vofs                  = i0vofs + 4
       
    ENDDO
    
    i0vr                       = i0atrmx
    i1atidr(i0vr)              = i0vofs + 1
    
    RETURN
    
  END SUBROUTINE T2VGRA_AT

  !C
  !C
  !C SUBROUTINE FOR DIFFUSION TENSOR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C  
  SUBROUTINE T2VGRA_DT
    
    USE T2COMM,ONLY:&
         i0spcs,i0dtrmx,i0dtcmx,i1dtidr,i1dtidc
    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1dtidr(1:i0dtrmx)         = 0
    i1dtidc(1:i0dtcmx)         = 0
    i0vofs                     = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    i0vr                       = 1
    i1dtidr(i0vr)              = i0vofs+1

    !C
    !C EQUATION FOR I
    !C
    i0vr                       = 2
    i1dtidr(i0vr)              = i0vofs+1

    !C
    !C EQUATION FOR Et
    !C
    i0vr                       = 3
    i1dtidr(i0vr)              = i0vofs+1

    !C
    !C EQUATION FOR Ep
    !C
    i0vr                       = 4
    i1dtidr(i0vr)              = i0vofs+1

    !C
    !C EQUATION FOR Er
    !C
    i0vr                       = 5
    i1dtidr(i0vr)              = i0vofs+1

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    DO i1 = 1,i0spcs
       
       !C
       !C EQUATION FOR N
       !C
       i0vr                    = 8*i1   - 2
       i1dtidr(i0vr)           = i0vofs + 1

       !C
       !C EQUATION FOR Fr
       !C
       i0vr                    = 8*i1   - 1
       i1dtidr(i0vr)           = i0vofs + 1

       !C
       !C EQUATION FOR Fb
       !C
       i0vr                    = 8*i1
       i1dtidr(i0vr)           = i0vofs + 1
       
       !C N
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       !C Fb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       !C Qb
       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr                    = 8*i1   + 1
       i1dtidr(i0vr)           = i0vofs + 1

       !C N
       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       !C Fb
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       !C Qb
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR P
       !C
       i0vr                    = 8*i1   + 2
       i1dtidr(i0vr)           = i0vofs + 1

       !C N
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       !C Fb
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       !C Qb
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR Qr
       !C
       i0vr                    = 8*i1   + 3
       i1dtidr(i0vr)           = i0vofs + 1

       !C
       !C EQUATION FOR Qb
       !C
       i0vr                    = 8*i1   + 4
       i1dtidr(i0vr)           = i0vofs + 1
       
       !C N
       i0vc                    = i0vr   - 6
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       !C Fb
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       !C Qb
       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1dtidr(i0vr)           = i0vofs + 1

       !C N
       i0vc                    = i0vr   - 7
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       !C Fb
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       !C P
       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       !C Qb
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
    ENDDO
    
    i0vr                       = i0dtrmx
    i1dtidr(i0vr)              = i0vofs + 1
    
    RETURN
    
  END SUBROUTINE T2VGRA_DT

  !C
  !C
  !C SUBROUTINE FOR GRADIENT VECOTR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C  
  SUBROUTINE T2VGRA_GV

    USE T2COMM,ONLY:&
         i0spcs,i0gvrmx,i0gvcmx,i1gvidr,i1gvidc
  
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1gvidr(1:i0gvrmx) = 0
    i1gvidc(1:i0gvcmx) = 0
    i0vofs             = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI'
    !C
    i0vr                       = 1
    i1gvidr(i0vr)              = i0vofs + 1

    !C Et
    i0vc                       = 3
    i0vg                       = i0vofs + 1
    i1gvidc(i0vg)              = i0vc
    
    i0vofs                     = i0vofs + 1
    
    !C
    !C EQUATION FOR I
    !C
    i0vr                       = 2
    i1gvidr(i0vr)              = i0vofs + 1

    !C Ep
    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1gvidc(i0vg)              = i0vc
    
    !C Er
    i0vc                       = 5
    i0vg                       = i0vofs + 2
    i1gvidc(i0vg)              = i0vc
    
    i0vofs                     = i0vofs + 2

    !C
    !C EQUATION FOR Et
    !C
    i0vr                       = 3
    i1gvidr(i0vr)              = i0vofs + 1
    
    !C
    !C EQUATION FOR Ep
    !C
    i0vr                       = 4
    i1gvidr(i0vr)              = i0vofs + 1

    !C I
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1gvidc(i0vg)              = i0vc
    
    i0vofs                     = i0vofs + 1
    
    !C
    !C EQUATION FOR Er
    !C
    i0vr                       = 5
    i1gvidr(i0vr)              = i0vofs + 1

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    DO i1 = 1,i0spcs
       !C
       !C EQUATION FOR N
       !C 
       i0vr                    = 8*i1   - 2
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR Fr
       !C 
       i0vr                    = 8*i1   - 1
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C P
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 1
       i1gvidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR Fb
       !C
       i0vr                    = 8*i1
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr                    = 8*i1   + 1
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C
       !C EQUATION FOR P
       !C
       i0vr                    = 8*i1   + 2
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C P
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1gvidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1

       !C
       !C EQUATION FOR Qr
       !C
       i0vr                    = 8*i1   + 3
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C N
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 1
       i1gvidc(i0vg)           = i0vc
       
       !C P
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1gvidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 2

       !C
       !C EQUATION FOR Qb
       !C
       i0vr                    = 8*i1   + 4
       i1gvidr(i0vr)           = i0vofs + 1

       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1gvidr(i0vr)           = i0vofs + 1
       
    ENDDO
    
    i0vr                       = i0gvrmx
    i1gvidr(i0vr)              = i0vofs + 1

    RETURN

  END SUBROUTINE T2VGRA_GV

  !C
  !C
  !C SUBROUTINE FOR GRADIENT TENSOR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C  
  SUBROUTINE T2VGRA_GT
    
    USE T2COMM,ONLY:&
         i0spcs,i0gtrmx,i0gtcmx,i1gtidr,i2gtidc
    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C VARIABLE-VARIABLE GRAPH
    
    !C INITIALIZE
    i1gtidr(1:i0gtrmx)     = 0
    i2gtidc(1:i0gtcmx,1:2) = 0
    i0vofs                 = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C

    !C
    !C EQUATION FOR PSI
    !C
    i0vr          = 1
    i1gtidr(i0vr) = i0vofs+1
    
    !C
    !C EQUATION FOR I
    !C
    i0vr          = 2
    i1gtidr(i0vr) = i0vofs+1

    !C
    !C EQUATION FOR Et
    !C
    i0vr          = 3
    i1gtidr(i0vr) = i0vofs+1
  
    !C
    !C EQUATION FOR Ep
    !C
    i0vr          = 4
    i1gtidr(i0vr) = i0vofs+1

    !C
    !C EQUATION FOR Er
    !C
    i0vr          = 5
    i1gtidr(i0vr) = i0vofs+1

    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    DO i1 = 1,i0spcs
       !C
       !C EQUATION FOR N
       !C
       i0vr          = 8*i1   - 2
       i1gtidr(i0vr) = i0vofs + 1

       !C
       !C EQUATION FOR Fr
       !C
       i0vr          = 8*i1   - 1
       i1gtidr(i0vr) = i0vofs + 1

       !C
       !C EQUATION FOR Fb
       !C
       i0vr          = 8*i1
       i1gtidr(i0vr) = i0vofs + 1
       
       !C N
       i0vc            = i0vr   - 2
       i0vg            = i0vofs + 1
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1
       
       !C Fb
       i0vc            = i0vr
       i0vg            = i0vofs + 2
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       !C P
       i0vc            = i0vr   + 2
       i0vg            = i0vofs + 3
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1
       
       !C Qb
       i0vc            = i0vr   + 4
       i0vg            = i0vofs + 4
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1
       
       i0vofs = i0vofs + 4
  
       !C
       !C EQUATION FOR Ft
       !C
       i0vr            = 8*i1   + 1
       i1gtidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR P
       !C
       i0vr            = 8*i1   + 2
       i1gtidr(i0vr)   = i0vofs + 1

       !C N (B)
       i0vc            = i0vr   - 4
       i0vg            = i0vofs + 1
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       !C N (Ub)
       i0vc            = i0vr   - 4
       i0vg            = i0vofs + 2
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 5*i1   - 2

       !C Fb (B)
       i0vc            = i0vr   - 2
       i0vg            = i0vofs + 3
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       !C Fb (Ub)
       i0vc            = i0vr   - 2
       i0vg            = i0vofs + 4
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 5*i1   - 2

       !C P (B)
       i0vc            = i0vr
       i0vg            = i0vofs + 5
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       !C P (Ub)
       i0vc            = i0vr
       i0vg            = i0vofs + 6
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 5*i1   - 2

       !C Qb (B)
       i0vc            = i0vr   + 2
       i0vg            = i0vofs + 7
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1
       
       !C Qb (Ub)
       i0vc            = i0vr   + 2
       i0vg            = i0vofs + 8
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 5*i1   - 2

       i0vofs = i0vofs + 8
       
       !C
       !C EQUATION FOR Qr
       !C
       i0vr            = 8*i1   + 3
       i1gtidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       i0vr            = 8*i1   + 4
       i1gtidr(i0vr)   = i0vofs + 1
       
       !C N (B)
       i0vc            = i0vr   - 6
       i0vg            = i0vofs + 1
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       !C Fb (B)
       i0vc            = i0vr   - 4
       i0vg            = i0vofs + 2
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       !C P (B)
       i0vc            = i0vr   - 2
       i0vg            = i0vofs + 3
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1
       
       !C Qb (B)
       i0vc            = i0vr
       i0vg            = i0vofs + 4
       i2gtidc(i0vg,1) = i0vc
       i2gtidc(i0vg,2) = 1

       i0vofs = i0vofs + 4
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1gtidr(i0vr)           = i0vofs + 1
       
    ENDDO
    
    i0vr                       = i0gtrmx
    i1gtidr(i0vr)              = i0vofs + 1

    RETURN
    
  END SUBROUTINE T2VGRA_GT

  !C
  !C
  !C SUBROUTINE FOR EXCITATION SCALAR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C  
  SUBROUTINE T2VGRA_ES
    
    USE T2COMM, ONLY:&
         i0spcs,i0esrmx,i0escmx,i1esidr,i1esidc
    
    INTEGER(i0ikind)::&
         i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vofs
    
    !C INITIALIZE
    i1esidr(1:i0esrmx)         = 0
    i1esidc(1:i0escmx)         = 0
    i0vofs                     = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    i0vr          = 1
    i1esidr(i0vr) = i0vofs + 1
    
    !C
    !C EQUATION FOR I
    !C
    i0vr          = 2
    i1esidr(i0vr) = i0vofs + 1

    !C
    !C EQUATION FOR Et
    !C
    i0vr          = 3
    i1esidr(i0vr) = i0vofs + 1
    
    DO i1 = 1, i0spcs
       i0vc          = 8*i1+1
       i0vg          = i0vofs + 1
       i1esidc(i0vg) = i0vc
       
       i0vofs = i0vofs + 1
    ENDDO
    
    !C
    !C EQUATION FOR Ep
    !C
    i0vr          = 4
    i1esidr(i0vr) = i0vofs + 1

    DO i1 = 1, i0spcs
       
       i0vc          = 8*i1
       i0vg          = i0vofs + 1
       i1esidc(i0vg) = i0vc
       
       i0vc          = 8*i1+1
       i0vg          = i0vofs + 2
       i1esidc(i0vg) = i0vc
       
       i0vofs = i0vofs + 2
       
    ENDDO
    
    !C
    !C EQUATION FOR Er
    !C
    i0vr          = 5
    i1esidr(i0vr) = i0vofs + 1
    
    DO i1 = 1, i0spcs
  
       i0vc          = 8*i1-2
       i0vg          = i0vofs + 1
       i1esidc(i0vg) = i0vc
       
       i0vofs = i0vofs + 1
       
    ENDDO



    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID 
    !C
    !C
    DO i1 = 1, i0spcs
       
       !C
       !C EQUATION FOR N
       !C
       i0vr          = 8*i1   - 2
       i1esidr(i0vr) = i0vofs + 1

       !C
       !C EQUATION FOR Fr
       !C
       i0vr          = 8*i1   - 1
       i1esidr(i0vr) = i0vofs + 1

       !C Ep
       i0vc          = 4
       i0vg          = i0vofs + 1
       i1esidc(i0vg) = i0vc
       
       !C Er
       i0vc          = 5
       i0vg          = i0vofs + 2
       i1esidc(i0vg) = i0vc

       !C Fb
       i0vc          = i0vr   + 1
       i0vg          = i0vofs + 3
       i1esidc(i0vg) = i0vc

       !C Ft
       i0vc          = i0vr   + 2
       i0vg          = i0vofs + 4
       i1esidc(i0vg) = i0vc

       i0vofs = i0vofs + 4
       
       !C
       !C EQUATION FOR Fb
       !C
       i0vr          = 8*i1
       i1esidr(i0vr) = i0vofs + 1
       
       DO j1 = 1, i0spcs
          
          i0vcx   = 8*j1
          
          IF(i1.EQ.j1)THEN
             !C Et
             i0vc          = 3
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc

             !C Ep
             i0vc          = 4
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc

             !C Fb
             i0vc          = i0vcx
             i0vg          = i0vofs + 3
             i1esidc(i0vg) = i0vc

             !C Qb
             i0vc          = i0vcx  + 4
             i0vg          = i0vofs + 4
             i1esidc(i0vg) = i0vc

             i0vofs = i0vofs + 4
             
          ELSE
             !C Fb
             i0vc          = i0vcx
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc

             !C Qb
             i0vc          = i0vcx  + 4
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc
             
             i0vofs            = i0vofs + 2
             
          END IF
       ENDDO

       !C
       !C EQUATION FOR Ft
       !C
       
       i0vr          = 8*i1   + 1
       i1esidr(i0vr) = i0vofs + 1
       
       DO j1=1,i0spcs
       
          i0vcx = 8*j1 + 1
          
          IF(i1.EQ.j1)THEN
             !C Et
             i0vc          = 3
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc
             
             !C Fr
             i0vc          = i0vcx  - 2
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc

             !C Ft
             i0vc          = i0vcx
             i0vg          = i0vofs + 3
             i1esidc(i0vg) = i0vc

             !C Qt
             i0vc          = i0vcx  + 4
             i0vg          = i0vofs + 4
             i1esidc(i0vg) = i0vc

             i0vofs = i0vofs + 4
             
          ELSE
             !C Ft
             i0vc          = i0vcx
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc
             
             !C Qt
             i0vc          = i0vcx  + 4
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc

             i0vofs = i0vofs + 2
             
          END IF
       ENDDO

       !C
       !C EQUATION FOR P
       !C
       i0vr          = 8*i1   + 2
       i1esidr(i0vr) = i0vofs + 1

       !C P
       i0vc          = i0vr
       i0vg          = i0vofs + 1
       i1esidc(i0vg) = i0vc
       
       i0vofs = i0vofs + 1
       
       !C
       !C EQUATION FOR Qr
       !C
       i0vr          = 8*i1   + 3
       i1esidr(i0vr) = i0vofs + 1

       !C Ep
       i0vc          = 4
       i0vg          = i0vofs + 1
       i1esidc(i0vg) = i0vc
       
       !C Er
       i0vc          = 5
       i0vg          = i0vofs + 2
       i1esidc(i0vg) = i0vc

       i0vc          = i0vr   + 1
       i0vg          = i0vofs + 3
       i1esidc(i0vg) = i0vc


       i0vc          = i0vr   + 2
       i0vg          = i0vofs + 4
       i1esidc(i0vg) = i0vc

       i0vofs = i0vofs + 4

       !C
       !C EQUATION FOR Qb
       !C
       i0vr          = 8*i1   + 4
       i1esidr(i0vr) = i0vofs + 1
       
       DO j1=1,i0spcs

          i0vcx = 8*j1 + 4

          IF(i1.EQ.j1)THEN
             !C Et
             i0vc          = 3
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc

             !C Ep
             i0vc          = 4
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc
             
             !C Fb
             i0vc          = i0vcx  - 4
             i0vg          = i0vofs + 3
             i1esidc(i0vg) = i0vc
             
             !C Qb
             i0vc          = i0vcx
             i0vg          = i0vofs + 4
             i1esidc(i0vg) = i0vc
             
             i0vofs = i0vofs + 4
          
          ELSE
             !C Fb
             i0vc          = i0vcx  - 4
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc

             !C Qb
             i0vc          = i0vcx
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc

             i0vofs = i0vofs + 2
             
          END IF
       ENDDO
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vr                    = 8*i1   + 5
       i1esidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs

          i0vcx = 8*j1 + 5
          
          IF(i1.EQ.j1)THEN
             !C Et
             i0vc          = 3
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc
             
             !C Ft
             i0vc          = i0vcx  - 4
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc
             
             !C Qr
             i0vc          = i0vcx  - 2
             i0vg          = i0vofs + 3
             i1esidc(i0vg) = i0vc

             !C Qt
             i0vc          = i0vcx
             i0vg          = i0vofs + 4
             i1esidc(i0vg) = i0vc
             
             i0vofs = i0vofs + 4
             
          ELSE
             !C Ft
             i0vc          = i0vcx  - 4
             i0vg          = i0vofs + 1
             i1esidc(i0vg) = i0vc

             !C Qt
             i0vc          = i0vcx
             i0vg          = i0vofs + 2
             i1esidc(i0vg) = i0vc
             
             i0vofs = i0vofs + 2
             
          END IF
       ENDDO
    ENDDO
    
    i0vr          = i0esrmx
    i1esidr(i0vr) = i0vofs + 1
    
    RETURN
    
  END SUBROUTINE T2VGRA_ES

  !C
  !C
  !C SUBROUTINE FOR EXCITATION VECTOR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C
  SUBROUTINE T2VGRA_EV
    
    USE T2COMM, ONLY:&
         i0spcs,i0evrmx,i0evcmx,i1evidr,i2evidc
    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1evidr(1:i0evrmx)     = 0
    i2evidc(1:i0evcmx,1:2) = 0
    i0vofs                 = 0
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    
    i0vr            = 1
    i1evidr(i0vr)   = i0vofs + 1
    
    !C
    !C EQUATION FOR I
    !C
    
    i0vr            = 2
    i1evidr(i0vr)   = i0vofs + 1
    
    !C
    !C EQUATION FOR Et
    !C
    
    i0vr            = 3
    i1evidr(i0vr)   = i0vofs + 1
    
    !C PSI' ( ln R)
    i0vc            = 1
    i0vg            = i0vofs + 1
    i2evidc(i0vg,1) = i0vc
    i2evidc(i0vg,2) = 2
    
    i0vofs = i0vofs + 1
    
    !C
    !C EQUATION FOR Ep
    !C
    
    i0vr          = 4
    i1evidr(i0vr) = i0vofs + 1
    
    !C
    !C EQUATION FOR Er
    !C
    
    i0vr          = 5
    i1evidr(i0vr) = i0vofs + 1
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID
    !C
    !C
    DO i1 = 1,i0spcs
      
       !C
       !C EQUATION FOR N
       !C
       i0vr            = 8*i1   - 2
       i1evidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR Fr
       !C
       i0vr            = 8*i1   - 1
       i1evidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR Fb
       !C
       i0vr            = 8*i1
       i1evidr(i0vr)   = i0vofs + 1
       
       !C Fb (B)
       i0vc            = i0vr
       i0vg            = i0vofs + 1
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1
       
       i0vofs = i0vofs + 1
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr            = 8*i1   + 1
       i1evidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR P
       !C
       i0vr            = 8*i1   + 2
       i1evidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR Qr
       !C
       i0vr            = 8*i1   + 3
       i1evidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       i0vr            = 8*i1   + 4
       i1evidr(i0vr)   = i0vofs + 1
       
       !C Et (B)
       i0vc            = 3
       i0vg            = i0vofs + 1
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1
       
       !C Et (N)
       i0vc            = 3
       i0vg            = i0vofs + 2
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   - 1
       
       !C Et (Fb)
       i0vc            = 3
       i0vg            = i0vofs + 3
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1
       
       !C Et (P)
       i0vc            = 3
       i0vg            = i0vofs + 4
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 1
       
       !C Et (Qb)
       i0vc            = 3
       i0vg            = i0vofs + 5
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 2
       
       !C Ep (B)
       i0vc            = 4
       i0vg            = i0vofs + 6
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1
       
       !C Ep (N)
       i0vc            = 4
       i0vg            = i0vofs + 7
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   - 1

       !C Ep (Fb)
       i0vc            = 4
       i0vg            = i0vofs + 8
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1  
       
       !C Ep (P)
       i0vc            = 4
       i0vg            = i0vofs + 9
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 1
       
       !C Ep (Qb)
       i0vc            = 4
       i0vg            = i0vofs + 10
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 2
       
       !C Fb (B)
       i0vc            = i0vr   - 4
       i0vg            = i0vofs + 11
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1

       !C Qb (B)
       i0vc            = i0vr   
       i0vg            = i0vofs + 12
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1
       
       i0vofs = i0vofs + 12
       
       !C
       !C EQUATION FOR Qt
       !C
       
       i0vr            = 8*i1   + 5
       i1evidr(i0vr)   = i0vofs + 1
       
       !C Et (B)
       i0vc            = 3
       i0vg            = i0vofs + 1
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1
       
       !C Et (N)
       i0vc            = 3
       i0vg            = i0vofs + 2
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   - 1
       
       !C Et (Fb)
       i0vc            = 3
       i0vg            = i0vofs + 3
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1
       
       !C Et (P)
       i0vc            = 3
       i0vg            = i0vofs + 4
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 1
       
       !C Et (Qb)
       i0vc            = 3
       i0vg            = i0vofs + 5
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 2
       
       !C Ep (B)
       i0vc            = 4
       i0vg            = i0vofs + 6
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 1
       
       !C Ep (N)
       i0vc            = 4
       i0vg            = i0vofs + 7
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   - 1
       
       !C Ep (Fb)
       i0vc            = 4
       i0vg            = i0vofs + 8
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1
       
       !C Ep (P)
       i0vc            = 4
       i0vg            = i0vofs + 9
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 1
       
       !C Ep (Qb)
       i0vc            = 4
       i0vg            = i0vofs + 10
       i2evidc(i0vg,1) = i0vc
       i2evidc(i0vg,2) = 5*i1   + 2
       
       i0vofs = i0vofs + 10
       
    ENDDO
    
    i0vr          = i0evrmx
    i1evidr(i0vr) = i0vofs + 1
    
    RETURN
    
  END SUBROUTINE T2VGRA_EV

  !C
  !C
  !C SUBROUTINE FOR EXCITATION TENSOR ARRAY
  !C 
  !C          MODIFIED 2013-12-05
  !C
  !C  
  SUBROUTINE T2VGRA_ET
    
    USE T2COMM, ONLY:&
         i0spcs,i0etrmx,i0etcmx,i1etidr,i2etidc
    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C VARIABLE-VARIABLE GRAPH

    !C INITIALIZE
    i1etidr(1:i0etrmx)     = 0
    i2etidc(1:i0etcmx,1:3) = 0
    i0vofs                 = 0

    !C
    !C
    !C VARIABLES OF PLASMA AS FIELD
    !C
    !C
    
    !C
    !C EQUATION FOR PSI
    !C
    i0vr            = 1
    i1etidr(i0vr)   = i0vofs+1

    !C
    !C EQUATION FOR I
    !C
    i0vr            = 2
    i1etidr(i0vr)   = i0vofs+1

    !C
    !C EQUATION FOR Et
    !C
    i0vr            = 3
    i1etidr(i0vr)   = i0vofs+1

    !C
    !C EQUATION FOR Ep
    !C
    i0vr            = 4
    i1etidr(i0vr)   = i0vofs+1

    !C
    !C EQUATION FOR Er
    !C
    i0vr            = 5
    i1etidr(i0vr)   = i0vofs+1
    
    !C
    !C
    !C VARIABLES OF PLASMA AS FLUID
    !C
    !C
    DO i1 = 1,i0spcs
       
       !C
       !C EQUATION FOR N
       !C
       i0vr            = 8*i1   - 2
       i1etidr(i0vr)   = i0vofs + 1

       !C
       !C EQUATION FOR Fr
       !C
       i0vr            = 8*i1   - 1
       i1etidr(i0vr)   = i0vofs + 1

       !C
       !C EQUATION FOR Fb
       !C
       i0vr            = 8*i1
       i1etidr(i0vr)   = i0vofs + 1
       
       !C Fb (B,B)
       i0vc            = i0vr
       i0vg            = i0vofs + 1
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       !C Ft (B,B)
       i0vc            = i0vr   + 1
       i0vg            = i0vofs + 2
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1
       
       !C Qb (B,B)
       i0vc            = i0vr   + 4
       i0vg            = i0vofs + 3
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1
       
       !C Qt (B,B)
       i0vc            = i0vr   + 5
       i0vg            = i0vofs + 4
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       i0vofs = i0vofs + 4
       
       !C
       !C EQUATION FOR Ft
       !C
       i0vr            = 8*i1   + 1
       i1etidr(i0vr)   = i0vofs + 1

       !C
       !C EQUATION FOR P
       !C
       i0vr            = 8*i1   + 2
       i1etidr(i0vr)   = i0vofs + 1
       
       !C Fb (B,B)
       i0vc            = i0vr   - 2
       i0vg            = i0vofs + 1
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       !C Fb (Ub,B)
       i0vc            = i0vr   - 2
       i0vg            = i0vofs + 2
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 5*i1   - 2
       i2etidc(i0vg,3) = 1

       !C Ft (B,B)
       i0vc            = i0vr   - 1
       i0vg            = i0vofs + 3
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       !C Ft (Ub,B)
       i0vc            = i0vr   - 1
       i0vg            = i0vofs + 4
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 5*i1   - 2
       i2etidc(i0vg,3) = 1

       !C Qb (B,B)
       i0vc            = i0vr   + 2
       i0vg            = i0vofs + 5
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1
       
       !C Qb (Ub,B)
       i0vc            = i0vr   + 2
       i0vg            = i0vofs + 6
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 5*i1   - 2
       i2etidc(i0vg,3) = 1
       
       !C Qt (B,B)
       i0vc            = i0vr   + 3
       i0vg            = i0vofs + 7
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       !C Qt (Ub,B)
       i0vc            = i0vr   + 3
       i0vg            = i0vofs + 8
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 5*i1   - 2
       i2etidc(i0vg,3) = 1

       i0vofs = i0vofs + 8
       
       !C
       !C EQUATION FOR Qr
       !C
       i0vr            = 8*i1   + 3
       i1etidr(i0vr)   = i0vofs + 1
       
       !C
       !C EQUATION FOR Qb
       !C
       i0vr            = 8*i1   + 4
       i1etidr(i0vr)   = i0vofs + 1
       
       !C Fb (B,B)
       i0vc            = i0vr   - 4
       i0vg            = i0vofs + 1
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       !C Ft (B,B)
       i0vc            = i0vr   - 3
       i0vg            = i0vofs + 2
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1
       
       !C Qb (B,B)
       i0vc            = i0vr
       i0vg            = i0vofs + 3
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       !C Qt (B,B)
       i0vc            = i0vr   + 1
       i0vg            = i0vofs + 4
       i2etidc(i0vg,1) = i0vc
       i2etidc(i0vg,2) = 1
       i2etidc(i0vg,3) = 1

       i0vofs = i0vofs + 4
       
       !C
       !C EQUATION FOR Qt
       !C
       i0vr             = 8*i1   + 5
       i1etidr(i0vr)    = i0vofs + 1
       
    ENDDO
    
    i0vr          = i0etrmx
    i1etidr(i0vr) = i0vofs + 1
    
    RETURN
    
  END SUBROUTINE T2VGRA_ET
  

  SUBROUTINE T2_VGRA_OUTPUT

    USE T2COMM
    IMPLICIT NONE
    INTEGER(i0ikind):: i1,j1

    OPEN(10,FILE='TEST_VMR.dat')
    DO i1=1,i0vgrmx
       WRITE(10,*)'i1=',i1,'I1VGIDR=',i1vgidr(i1)
    ENDDO
    
    OPEN(10,FILE='TEST_VMC.dat')
    DO i1=1,i0vgcmx
       WRITE(10,*)'i1=',i1,'I1VGIDC=',i1vgidc(i1)
    ENDDO
    
    OPEN(10,FILE='TEST_VMX.dat')
    DO i1=1,i0vmax
    DO j1=1,i0vmax
       WRITE(10,*)'i1=',i1,'j1=',j1,'I2VTBL=',i2vtbl(i1,j1)
    ENDDO
    ENDDO

    OPEN(10,FILE='TEST_MSR.dat')
    DO i1=1,i0msrmx
       WRITE(10,*)'i1=',i1,'I1MSIDR=',i1msidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_MSC.dat')
    DO i1=1,i0mscmx
       WRITE(10,*)'i1=',i1,'I1MSIDC=',i1msidc(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_AVR.dat')
    DO i1=1,i0avrmx
       WRITE(10,*)'i1=',i1,'I1AVIDR=',i1avidr(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_AVC.dat')
    DO i1=1,i0avcmx
       WRITE(10,*)'i1=',i1,'I1AVIDC=',i1avidc(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ATR.dat')
    DO i1=1,i0atrmx
       WRITE(10,*)'i1=',i1,'I1ATIDR=',i1atidr(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ATC.dat')
    DO i1=1,i0atcmx
       WRITE(10,*)'i1=',i1,&
            'I1ATIDC1=',i2atidc(i1,1),'I1ATIDC1=',i2atidc(i1,2)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_DTR.dat')
    DO i1=1,i0dtrmx
       WRITE(10,*)'i1=',i1,'I1DTIDR=',i1dtidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_DTC.dat')
    DO i1=1,i0dtcmx
       WRITE(10,*)'i1=',i1,'I1DTIDC=',i1dtidc(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GVR.dat')
    DO i1=1,i0gvrmx
       WRITE(10,*)'i1=',i1,'I1GVIDR=',i1gvidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_GVC.dat')
    DO i1=1,i0gvcmx
       WRITE(10,*)'i1=',i1,'I1GVIDC=',i1gvidc(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GTR.dat')
    DO i1=1,i0gtrmx
       WRITE(10,*)'i1=',i1,'I1GTIDR=',i1gtidr(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GTC.dat')
    DO i1=1,i0gtcmx
       WRITE(10,*)'i1=',i1,&
            'I2GTIDC1=',i2gtidc(i1,1),'I2GTIDC2=',i2gtidc(i1,2)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_ESR.dat')
    DO i1=1,i0esrmx
       WRITE(10,*)'i1=',i1,'I1ESIDR=',i1esidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_ESC.dat')
    DO i1=1,i0escmx
       WRITE(10,*)'i1=',i1,'I1ESIDC=',i1esidc(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_EVR.dat')
    DO i1=1,i0evrmx
       WRITE(10,*)'i1=',i1,'I1EVIDR=',i1evidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_EVC.dat')
    DO i1=1,i0evcmx
       WRITE(10,*)'i1=',i1,&
            'I2EVIDC1=',i2evidc(i1,1),'I2EVIDC2=',i2evidc(i1,2)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_ETR.dat')
    DO i1=1,i0etrmx
       WRITE(10,*)'i1=',i1,'I1ETIDR=',i1etidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_ETC.dat')
    DO i1=1,i0etcmx
       WRITE(10,*)'i1=',i1,'I2ETIDC1=',i2etidc(i1,1),&
            'I2ETIDC2=',i2etidc(i1,2),'I2ETIDC3=',i2etidc(i1,3)
    ENDDO
    CLOSE(10)
    
    RETURN
  END SUBROUTINE T2_VGRA_OUTPUT
    
END MODULE T2VGRA
