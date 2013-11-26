!C------------------------------------------------------------------
!C
!C         MODULE T2VGRA
!C       
!C         VARIABLE GRAPH GENERATOR FOR TASK/T2
!C
!C
!C------------------------------------------------------------------

MODULE T2VGRA
  USE T2CNST,ONLY:&
       i0ikind,i0rkind
  
  IMPLICIT NONE
  
  PUBLIC T2_VGRA
  
PRIVATE 
  
CONTAINS 
  
  SUBROUTINE T2_VGRA

    CALL T2_VGRA_VM
    CALL T2_VGRA_MS
    CALL T2_VGRA_AV
    CALL T2_VGRA_AT
    CALL T2_VGRA_DT
    CALL T2_VGRA_GV
    CALL T2_VGRA_GT
    CALL T2_VGRA_ES
    CALL T2_VGRA_EV
    CALL T2_VGRA_ET
    
    RETURN
    
  END SUBROUTINE T2_VGRA

  !C
  !C SUBROUTINE FOR VALIABLE MATRIX 
  !C

  SUBROUTINE T2_VGRA_VM
    
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

    !C EQUATION FOR Bp

    i0vr                       = 1
    i1vgidr(i0vr)              = i0vofs+1
    
    i0vc                       = 1
    i0vg                       = i0vofs+1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 5
    i0vg                       = i0vofs+2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vofs                     = i0vofs + 2

    !C EQUATION FOR Bt

    i0vr                       = 2
    i1vgidr(i0vr)              = i0vofs + 1
    
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 3
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 4
    i0vg                       = i0vofs + 3
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vofs                     = i0vofs + 3

    !C EQUATION FOR Er

    i0vr                       = 3
    i1vgidr(i0vr)              = i0vofs + 1

    i0vc                       = 3
    i0vg                       = i0vofs+1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 4
    i0vg                       = i0vofs+2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1-2
       i0vg                    = i0vofs+2+i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO

    i0vofs                     = i0vofs + 2 + i0spcs

    !C EQUATION FOR Ep

    i0vr                       = 4
    i1vgidr(i0vr)              = i0vofs + 1
    
    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 4
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1
       i0vg                    = i0vofs + 1 + 2*i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + 2 + 2*i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO

    i0vofs                     = i0vofs + 2 + 2*i0spcs

    !C EQUATION FOR Et
    
    i0vr                       = 5
    i1vgidr(i0vr)               = i0vofs + 1
    
    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    i0vc                       = 5
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + 2 + i1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO
    
    i0vofs = i0vofs + 2 + i0spcs
    
    DO i1 = 1, i0spcs
       
       !C EQUATION FOR N

       i0vr                    = 8*i1   - 2
       i1vgidr(i0vr)            = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Fr
       
       i0vr                    = 8*i1   - 1
       i1vgidr(i0vr)           = i0vofs + 1
       
       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 4
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 5
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 6
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vofs                  = i0vofs + 6
       
       !C EQUATION FOR Fb

       i0vr                    = 8*i1
       i1vgidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs
          
          i0vcx   = 8*j1
          
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = 5
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 2
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  + 5
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 8
         
          ELSE
             
             i0vc              = i0vcx
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 2
          
          END IF
       ENDDO
       
       !C EQUATION FOR Ft
       
       i0vr                    = 8*i1   + 1
       i1vgidr(i0vr)           = i0vofs + 1
       
       DO j1=1,i0spcs
          
          i0vcx                = 8*j1 + 1
          
          IF(i1.EQ.j1)THEN
             i0vc              = 5
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  - 3
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 1
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 3
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 8
          ELSE
             
             i0vc              = i0vcx  
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vofs = i0vofs + 2
             
          END IF
       ENDDO

       !C EQUATION FOR P
       
       i0vr                    = 8*i1   + 2
       i1vgidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 5
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 6
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vofs                  = i0vofs + 6

       !C EQUATION FOR Qr

       i0vr                    = 8*i1   + 3
       i1vgidr(i0vr)           = i0vofs + 1
       
       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 3
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 4
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr  
       i0vg                    = i0vofs + 5
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 6
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 7
       i1vgidc(i0vg)           = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
       
       i0vofs                  = i0vofs + 7

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1vgidr(i0vr)           = i0vofs + 1
       
       DO j1=1,i0spcs
          
          i0vcx   = 8*j1+4
          
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = 5
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 6
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 3
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 8

          ELSE
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

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
       
       DO j1=1,i0spcs
          i0vcx   = 8*j1+5
          IF(i1.EQ.j1)THEN
             
             i0vc              = 4
             i0vg              = i0vofs + 1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = 5
             i0vg              = i0vofs + 2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 7
             i0vg              = i0vofs + 3
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 5
             i0vg              = i0vofs + 4
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 5
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 3
             i0vg              = i0vofs + 6
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 7
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 1
             i0vg              = i0vofs + 8
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx
             i0vg              = i0vofs + 9
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 9

          ELSE
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx
             i0vg              = i0vofs +  2
             i1vgidc(i0vg)     = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 2
          END IF
       ENDDO
    ENDDO

    i0vr                       =  6 + 8*i0spcs
    i1vgidr(i0vr)              =  i0vofs + 1

    !CCCC
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
    
    RETURN
    
  END SUBROUTINE T2_VGRA_VM

  !C
  !C
  !C
  SUBROUTINE T2_VGRA_MS
    
    USE T2COMM,ONLY:&
         i0msrmx,i0mscmx,i0spcs,i1msidr,i1msidc,i0dbg

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1msidr(1:i0msrmx)         = 0
    i1msidc(1:i0mscmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1msidr(i0vr)              = i0vofs + 1
    
    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc
    
    i0vofs                     = i0vofs+1
    
    !C EQUATION FOR Bt
    i0vr                       = 2
    i1msidr(i0vr)              = i0vofs + 1

    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc

    i0vofs                     = i0vofs+1

    !C EQUATION FOR Er
    i0vr                       = 3
    i1msidr(i0vr)              = i0vofs + 1

    i0vc                       = 3
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc
    
    i0vofs = i0vofs+1

    !C EQUATION FOR Ep
    i0vr                       = 4
    i1msidr(i0vr)              = i0vofs + 1
    
    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1msidc(i0vg)              = i0vc

    i0vofs = i0vofs+1

    !C EQUATION FOR Et
    i0vr                        = 5
    i1msidr(i0vr)               = i0vofs + 1

    i0vc                        = 5
    i0vg                        = i0vofs + 1
    i1msidc(i0vg)               = i0vc

    i0vofs                      = i0vofs+1

    DO i1 = 1,i0spcs

       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1msidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc
       
       i0vofs = i0vofs +1

       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1msidr(i0vr)           = i0vofs + 1
       
       i0vc                       = i0vr
       i0vg                       = i0vofs + 1
       i1msidc(i0vg)              = i0vc
       
       i0vofs = i0vofs+1
       
       !C EQUATION FOR Fb
       
       i0vr                    = 8*i1
       i1msidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs +1

       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1msidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1msidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1msidr(i0vr)           = i0vofs + 1

       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs+1
       
       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1msidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1msidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1msidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1
    ENDDO
    
    i0vr                       = i0msrmx
    i1msidr(i0vr)              = i0vofs + 1
    
    !CCCC
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
    RETURN

  END SUBROUTINE T2_VGRA_MS
  
  !C
  !C
  !C
  SUBROUTINE T2_VGRA_AV

    USE T2COMM,ONLY:&
         i0spcs, i0avrmx, i0avcmx, i1avidr, i1avidc

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1avidr(1:i0avrmx)         = 0
    i1avidc(1:i0avcmx)         = 0
    i0vofs                     = 0
    
    !C EQUATION FOR Bp
    i0vr                       = 1
    i1avidr(i0vr)              = i0vofs + 1

    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C EQUATION FOR Bt
    i0vr                       = 2
    i1avidr(i0vr)              = i0vofs + 1

    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C EQUATION FOR Er
    i0vr                       = 3
    i1avidr(i0vr)              = i0vofs + 1

    i0vc                       = 3
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vc                       = 4
    i0vg                       = i0vofs + 2
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 2

    !C EQUATION FOR Ep
    i0vr                       = 4
    i1avidr(i0vr)              = i0vofs + 1

    i0vc                       = 4
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C EQUATION FOR Et
    i0vr                       = 5
    i1avidr(i0vr)              = i0vofs + 1

    i0vc                       = 1
    i0vg                       = i0vofs + 1
    i1avidc(i0vg)              = i0vc

    i0vc                       = 5
    i0vg                       = i0vofs + 2
    i1avidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 2

    DO i1 = 1,i0spcs

       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1avidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1
       
       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1avidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1avidr(i0vr)           = i0vofs + 1
       

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc


       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 2
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 2
       
       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1avidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1avidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1avidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1avidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1avidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1avidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 3
       
       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1avidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1avidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1avidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 2
       
    ENDDO
    

    i0vr                       = i0avrmx
    i1avidr(i0vr)              = i0vofs + 1
    !CCCC
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

    RETURN
    
  END SUBROUTINE T2_VGRA_AV

  !C
  !C
  !C
  SUBROUTINE T2_VGRA_AT
    
    USE T2COMM,ONLY:&
         i0spcs,i0atrmx,i0atcmx,i1atidr,i1atidc,i1atws

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    
    !C INITIALIZE
    i1atidr(1:i0atrmx)         = 0
    i1atidc(1:i0atcmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1atidr(i0vr)              = i0vofs + 1

    !C EQUATION FOR Bt
    i0vr                       = 2
    i1atidr(i0vr)              = i0vofs + 1


    !C EQUATION FOR Er
    i0vr                       = 3
    i1atidr(i0vr)              = i0vofs + 1
    
    !C EQUATION FOR Ep
    i0vr                       = 4
    i1atidr(i0vr)              = i0vofs + 1

    !C EQUATION FOR Et
    i0vr                        = 5
    i1atidr(i0vr)               = i0vofs + 1
    
    DO i1 = 1,i0spcs
       
       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1atidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1atidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1atidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 2
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 3
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
       
       i0vc                    = i0vr   + 5
       i0vg                    = i0vofs + 4
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1atidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 1
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   
       i0vg                    = i0vofs + 2
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
        
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 3
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
       
       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 4
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
       
       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1atidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
       
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 4
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vofs                  = i0vofs + 4

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1atidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1atidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 2
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   
       i0vg                    = i0vofs + 3
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
       
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 4
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1atidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 1
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1
       
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vc                    = i0vr   
       i0vg                    = i0vofs + 4
       i1atidc(i0vg)           = i0vc
       i1atws( i0vg)           = 1

       i0vofs                  = i0vofs + 4
       
    ENDDO
    
    i0vr                       = i0atrmx
    i1atidr(i0vr)              = i0vofs + 1

    !C FOR DEBUG
    OPEN(10,FILE='TEST_ATR.dat')
    DO i1=1,i0atrmx
       WRITE(10,*)'i1=',i1,'I1ATIDR=',i1atidr(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ATC.dat')
    DO i1=1,i0atcmx
       WRITE(10,*)'i1=',i1,'I1ATIDC=',i1atidc(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ATW.dat')
    DO i1=1,i0atcmx
       WRITE(10,*)'i1=',i1,'I1ATWS=',i1atws(i1)
    ENDDO
    CLOSE(10)

    RETURN

  END SUBROUTINE T2_VGRA_AT
  
  SUBROUTINE T2_VGRA_DT
    
    USE T2COMM,ONLY:&
         i0spcs,i0dtrmx,i0dtcmx,i1dtidr,i1dtidc
    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1dtidr(1:i0dtrmx)         = 0
    i1dtidc(1:i0dtcmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1dtidr(i0vr)              = i0vofs+1


    !C EQUATION FOR Bt
    i0vr                       = 2
    i1dtidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Er
    i0vr                       = 3
    i1dtidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Ep
    i0vr                       = 4
    i1dtidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Et
    i0vr                        = 5
    i1dtidr(i0vr)               = i0vofs+1

    DO i1 = 1,i0spcs

       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1dtidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1dtidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1dtidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR FT
       i0vr                    = 8*i1   + 1
       i1dtidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1dtidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1dtidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1dtidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 6
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1dtidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 7
       i0vg                    = i0vofs + 1
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 2
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 3
       i1dtidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 4
       i1dtidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
    ENDDO

    i0vr                       = i0dtrmx
    i1dtidr(i0vr)              = i0vofs + 1

    !CCCC
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

    RETURN

  END SUBROUTINE T2_VGRA_DT

  SUBROUTINE T2_VGRA_GV

    USE T2COMM,ONLY:&
         i0spcs,i0gvrmx,i0gvcmx,i1gvidr,i1gvidc
  
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1gvidr(1:i0gvrmx)         = 0
    i1gvidc(1:i0gvcmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1gvidr(i0vr)              = i0vofs + 1

    i0vc                       = 5
    i0vg                       = i0vofs + 1
    i1gvidc(i0vg)              = i0vc
    
    i0vofs                     = i0vofs + 1

    !C EQUATION FOR Bt
    i0vr                       = 2
    i1gvidr(i0vr)              = i0vofs + 1

    i0vc                       = 3
    i0vg                       = i0vofs + 1
    i1gvidc(i0vg)              = i0vc
    
    i0vc                       = 4
    i0vg                       = i0vofs + 2
    i1gvidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 2

    !C EQUATION FOR Er
    i0vr                       = 3
    i1gvidr(i0vr)              = i0vofs + 1
    
    !C EQUATION FOR Ep
    i0vr                       = 4
    i1gvidr(i0vr)              = i0vofs + 1

    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1gvidc(i0vg)              = i0vc

    i0vofs                     = i0vofs + 1

    !C EQUATION FOR Et
    i0vr                       = 5
    i1gvidr(i0vr)              = i0vofs + 1
    
    DO i1 = 1,i0spcs
       
       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1gvidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 1
       i1gvidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1gvidr(i0vr)           = i0vofs + 1
              
       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1gvidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1gvidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1gvidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1gvidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 1
       i1gvidc(i0vg)           = i0vc
       
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1gvidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 2

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1gvidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1gvidr(i0vr)           = i0vofs + 1
       
    ENDDO
    
    i0vr                       = i0gvrmx
    i1gvidr(i0vr)              = i0vofs + 1

    !CCCC
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

    RETURN

  END SUBROUTINE T2_VGRA_GV

  SUBROUTINE T2_VGRA_GT
    
    USE T2COMM,ONLY:&
         i0spcs,i0gtrmx,i0gtcmx,i1gtidr,i1gtidc,i1gtws

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C VARIABLE-VARIABLE GRAPH
    
    !C INITIALIZE
    i1gtidr(1:i0gtrmx)         = 0
    i1gtidc(1:i0gtcmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1gtidr(i0vr)              = i0vofs+1


    !C EQUATION FOR Bt
    i0vr                       = 2
    i1gtidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Er
    i0vr                       = 3
    i1gtidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Ep
    i0vr                       = 4
    i1gtidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Et
    i0vr                        = 5
    i1gtidr(i0vr)               = i0vofs+1

    DO i1 = 1,i0spcs

       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1gtidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1gtidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1gtidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 4
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1gtidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1gtidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 5*i1   - 3

       
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 3
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 4
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 5*i1   - 3

       i0vc                    = i0vr
       i0vg                    = i0vofs + 5
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 6
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 5*i1   - 3

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 7
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 8
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 5*i1   - 3

       i0vofs                  = i0vofs + 8
       
       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1gtidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1gtidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 6
       i0vg                    = i0vofs + 1
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 3
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1gtidc(i0vg)           = i0vc
       i1gtws( i0vg)           = 1

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1gtidr(i0vr)           = i0vofs + 1

    ENDDO
    
    i0vr                       = i0gtrmx
    i1gtidr(i0vr)              = i0vofs + 1

    !CCCC
    OPEN(10,FILE='TEST_GTR.dat')
    DO i1=1,i0gtrmx
       WRITE(10,*)'i1=',i1,'I1GTIDR=',i1gtidr(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_GTC.dat')
    DO i1=1,i0gtcmx
       WRITE(10,*)'i1=',i1,'I1GTIDC=',i1gtidc(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_GTW.dat')
    DO i1=1,i0gtcmx
       WRITE(10,*)'i1=',i1,'I1GTWS=',i1gtws(i1)
    ENDDO
    CLOSE(10)
    
    RETURN
  END SUBROUTINE T2_VGRA_GT
  
  SUBROUTINE T2_VGRA_ES
    
    USE T2COMM, ONLY:&
         i0spcs,i0esrmx,i0escmx,i1esidr,i1esidc
    
    INTEGER(i0ikind)::&
         i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vofs
    
    !C INITIALIZE
    i1esidr(1:i0esrmx)         = 0
    i1esidc(1:i0escmx)         = 0
    i0vofs                     = 0
    
    !C EQUATION FOR BP
    i0vr                       = 1
    i1esidr(i0vr)              = i0vofs + 1
    
    !C EQUATION FOR BT
    i0vr                       = 2
    i1esidr(i0vr)              = i0vofs + 1
    
    !C EQUATION FOR ER
    i0vr                       = 3
    i1esidr(i0vr)              = i0vofs + 1
    
    DO i1 = 1, i0spcs
  
       i0vc                    = 8*i1-2
       i0vg                    = i0vofs + 1
       i1esidc(i0vg)           = i0vc
       
       i0vofs                  = i0vofs + 1
       
    ENDDO

    !C EQUATION FOR EP
    i0vr                       = 4
    i1esidr(i0vr)              = i0vofs + 1

    DO i1 = 1, i0spcs

       i0vc                    = 8*i1
       i0vg                    = i0vofs + 1
       i1esidc(i0vg)           = i0vc

       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + 2
       i1esidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 2

    ENDDO

    !C EQUATION FOR ET
    i0vr                       = 5
    i1esidr(i0vr)              = i0vofs + 1

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + 1
       i1esidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1
    ENDDO

    DO i1 = 1,i0spcs

       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1esidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR FR
       i0vr                    = 8*i1   - 1
       i1esidr(i0vr)           = i0vofs + 1

       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1esidc(i0vg)           = i0vc


       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1esidc(i0vg)           = i0vc


       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1esidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1esidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR FB
       i0vr                    = 8*i1
       i1esidr(i0vr)           = i0vofs + 1
       
       DO j1 = 1, i0spcs
          
          i0vcx   = 8*j1
          
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc

             i0vc              = 5
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx
             i0vg              = i0vofs + 3
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 4
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 4

          ELSE
             i0vc              = i0vcx
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 2

          END IF
       ENDDO

       !C EQUATION FOR FT
       i0vr                    = 8*i1   + 1
       i1esidr(i0vr)           = i0vofs + 1
       
       DO j1=1,i0spcs
       
          i0vcx                = 8*j1 + 1
          
          IF(i1.EQ.j1)THEN
             i0vc              = 5
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc
             
             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx
             i0vg              = i0vofs + 3
             i1esidc(i0vg)     = i0vc


             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 4
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 4
          ELSE
             i0vc              = i0vcx
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 2
          END IF
       ENDDO

       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1esidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1esidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 1

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1esidr(i0vr)           = i0vofs + 1

       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1esidc(i0vg)           = i0vc


       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1esidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1esidc(i0vg)           = i0vc


       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1esidc(i0vg)           = i0vc

       i0vofs                  = i0vofs + 4

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1esidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs

          i0vcx = 8*j1 + 4

          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc

             i0vc              = 5
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 3
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx
             i0vg              = i0vofs + 4
             i1esidc(i0vg)     = i0vc
             
             i0vofs            = i0vofs + 4
          
          ELSE
             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 2
             
          END IF
       ENDDO

       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1esidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs

          i0vcx = 8*j1 + 5
          
          IF(i1.EQ.j1)THEN
             i0vc              = 5
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc
             
             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 3
             i1esidc(i0vg)     = i0vc

             i0vc              = i0vcx
             i0vg              = i0vofs + 4
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 4

          ELSE

             i0vc              = i0vcx  - 4
             i0vg              = i0vofs + 1
             i1esidc(i0vg)     = i0vc


             i0vc              = i0vcx
             i0vg              = i0vofs + 2
             i1esidc(i0vg)     = i0vc

             i0vofs            = i0vofs + 2

          END IF
       ENDDO
    ENDDO
    
    i0vr                       = i0esrmx
    i1esidr(i0vr)              = i0vofs + 1
    
    !CCCC
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

    RETURN
    
  END SUBROUTINE T2_VGRA_ES
  
  SUBROUTINE T2_VGRA_EV
    
    USE T2COMM, ONLY:&
         i0spcs,i0evrmx,i0evcmx,i1evidr,i1evidc,i1evws

    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C INITIALIZE
    i1evidr(1:i0evrmx)         = 0
    i1evidc(1:i0evcmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1evidr(i0vr)              = i0vofs + 1

    !C EQUATION FOR Bt
    i0vr                       = 2
    i1evidr(i0vr)              = i0vofs + 1


    !C EQUATION FOR Er
    i0vr                       = 3
    i1evidr(i0vr)              = i0vofs + 1
    
    !C EQUATION FOR Ep
    i0vr                       = 4
    i1evidr(i0vr)              = i0vofs + 1

    !C EQUATION FOR Et
    i0vr                        = 5
    i1evidr(i0vr)               = i0vofs + 1
    
    DO i1 = 1,i0spcs
       
       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1evidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1evidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1evidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1

       i0vofs                  = i0vofs + 1
       
       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1evidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1evidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1evidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1evidr(i0vr)           = i0vofs + 1
       
       i0vc                    = 4
       i0vg                    = i0vofs + 1
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1

       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 2

       i0vc                    = 4
       i0vg                    = i0vofs + 3
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 1
       
       i0vc                    = 4
       i0vg                    = i0vofs + 4
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   
       
       i0vc                    = 4
       i0vg                    = i0vofs + 5
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   + 1
       

       i0vc                    = 5
       i0vg                    = i0vofs + 6
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1

       i0vc                    = 5
       i0vg                    = i0vofs + 7
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 2

       i0vc                    = 5
       i0vg                    = i0vofs + 8
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 1

       i0vc                    = 5
       i0vg                    = i0vofs + 9
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   

       i0vc                    = 5
       i0vg                    = i0vofs + 10
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   + 1

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 11
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1

       i0vc                    = i0vr   
       i0vg                    = i0vofs + 12
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1       
       
       i0vofs                  = i0vofs + 12
       
       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1evidr(i0vr)           = i0vofs + 1
       
       i0vc                    = 4
       i0vg                    = i0vofs + 1
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1

       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 2

       i0vc                    = 4
       i0vg                    = i0vofs + 3
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 1

       i0vc                    = 4
       i0vg                    = i0vofs + 4
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1

       i0vc                    = 4
       i0vg                    = i0vofs + 5
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   + 1

       i0vc                    = 5
       i0vg                    = i0vofs + 6
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 1

       i0vc                    = 5
       i0vg                    = i0vofs + 7
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 2
 
       i0vc                    = 5
       i0vg                    = i0vofs + 8
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   - 1

       i0vc                    = 5
       i0vg                    = i0vofs + 9
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1

       i0vc                    = 5
       i0vg                    = i0vofs + 10
       i1evidc(i0vg)           = i0vc
       i1evws( i0vg)           = 5*i1   + 1

       i0vofs                  = i0vofs + 10
       
    ENDDO
    
    i0vr                       = i0evrmx
    i1evidr(i0vr)              = i0vofs + 1
    
    !C FOR DEBUG
    OPEN(10,FILE='TEST_EVR.dat')
    DO i1=1,i0evrmx
       WRITE(10,*)'i1=',i1,'I1EVIDR=',i1evidr(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_EVC.dat')
    DO i1=1,i0evcmx
       WRITE(10,*)'i1=',i1,'I1EVIDC=',i1evidc(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_EVW.dat')
    DO i1=1,i0evcmx
       WRITE(10,*)'i1=',i1,'I1EVWS=',i1evws(i1)
    ENDDO
    CLOSE(10)

    RETURN

  END SUBROUTINE T2_VGRA_EV

  SUBROUTINE T2_VGRA_ET
    
    USE T2COMM, ONLY:&
         i0spcs,i0etrmx,i0etcmx,i1etidr,i1etidc,i2etws
    
    INTEGER(i0ikind)::&
         i1,i0vr,i0vc,i0vg,i0vofs
    
    !C VARIABLE-VARIABLE GRAPH

    !C INITIALIZE
    i1etidr(1:i0etrmx)         = 0
    i1etidc(1:i0etcmx)         = 0
    i0vofs                     = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i1etidr(i0vr)              = i0vofs+1


    !C EQUATION FOR Bt
    i0vr                       = 2
    i1etidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Er
    i0vr                       = 3
    i1etidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Ep
    i0vr                       = 4
    i1etidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Et
    i0vr                        = 5
    i1etidr(i0vr)               = i0vofs+1

    DO i1 = 1,i0spcs

       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i1etidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i1etidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i1etidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 2
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 3
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 5
       i0vg                    = i0vofs + 4
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i1etidr(i0vr)           = i0vofs + 1
       
       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i1etidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 5*i1   - 3
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 4
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 5*i1   - 3
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 5
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 6
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 5*i1   - 3
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 7
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 8
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 5*i1   - 3
       i2etws(2,i0vg)          = 1

       i0vofs                  = i0vofs + 8
       
       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i1etidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i1etidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 2
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 4
       i1etidc( i0vg)          = i0vc
       i2etws(1,i0vg)          = 1
       i2etws(2,i0vg)          = 1

       i0vofs                  = i0vofs + 4
       
       !C EQUATION FOR Qt
       i0vr                    = 8*i1   + 5
       i1etidr(i0vr)           = i0vofs + 1
    ENDDO
    
    i0vr                       = i0etrmx
    i1etidr(i0vr)              = i0vofs + 1
    
    !CCCC
    OPEN(10,FILE='TEST_ETR.dat')
    DO i1=1,i0etrmx
       WRITE(10,*)'i1=',i1,'I1ETIDR=',i1etidr(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_ETC.dat')
    DO i1=1,i0etcmx
       WRITE(10,*)'i1=',i1,'I1ETIDC=',i1etidc(i1)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='TEST_ETW.dat')
    DO i1=1,i0etcmx
          WRITE(10,*)'i1=',i1,'ETWS1=',i2etws(1,i1),'ETWS1=',i2etws(2,i1)
    ENDDO
    CLOSE(10)

    RETURN
  END SUBROUTINE T2_VGRA_ET

END MODULE T2VGRA
