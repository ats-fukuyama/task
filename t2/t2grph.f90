
!C------------------------------------------------------------------
!C
!C         MODULE T2GPRH
!C       
!C         NODE AND VARIABLE GRAPH GENERATOR FOR TASK/T2
!C
!C        (ARGUMANT DECLARATION PART WILL BE FEEDBACKED TO T2COMM)
!C
!C------------------------------------------------------------------
MODULE T2GRPH
  USE T2CNST,ONLY:&
       i0ikind,i0rkind,i0lmax_tmp
  USE T2COMM

  IMPLICIT NONE

CONTAINS 
  
  !C------------------------------------------------------------------
  !C
  !C          T2GRPH_MAIN
  !C
  !C          MAIN SUBROUTINE OF T2GRPH
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2GRPH_MAIN
    
    CALL T2INIT_GRPH
    
    CALL T2GRPH_SETUP
    
    
    CALL T2GRPH_VGRAPH
    CALL T2GRPH_VMGRAPH
    CALL T2GRPH_VAGRAPH
    CALL T2GRPH_VDGRAPH
    CALL T2GRPH_VGGRAPH    
    CALL T2GRPH_VEGRAPH
    
    SELECT CASE (i0nmax0)
    CASE( 4)
       CALL T2GRPH_NGRAPH1
    CASE( 8)
       CALL T2GRPH_NGRAPH2
    CASE(12)
       CALL T2GRPH_NGRAPH3
    END SELECT

    CALL T2GRPH_MFCS
    
    CALL T2GRPH_OUTPUT
    
    CALL T2GRPH_TERM
    RETURN
  END SUBROUTINE T2GRPH_MAIN
  
  !C------------------------------------------------------------------
  !C
  !C         PARAMETER SETTING ROUTINE WITH NAMELIST FOR
  !C  
  !C        (THIS WILL BE FEEDBACKED TO T2INIT.F90)
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2INIT_GRPH
    INTEGER(i0ikind)::i1,&
         i1mesh_level(  0:i0lmax_tmp+1),&
         i1rdiv_number(-1:i0lmax_tmp  ),&
         i0pdiv_number,i0mesh_level
    REAL(   i0rkind)::&
         d1rec_tmp(  0:i0lmax_tmp  )
    NAMELIST /T2GRPH/ &
         i0spcs, i0nmax0, i0lmax, i1mesh_level, i0pdiv_number, i1rdiv_number, d1rec_tmp
     
    i1mesh_level( 0:i0lmax_tmp+1) = 0
    i1rdiv_number(-1:i0lmax_tmp  ) = 0

    OPEN(10,file='t2cprm.nl')
    READ(10,T2GRPH)
    CLOSE(10)
    print*,i0lmax
  
    IF(i0lmax.GT.i0lmax_tmp)THEN
       WRITE(6,*)'STACK OVERFLOW IN T2_INIT: I0LMAX'
       WRITE(6,*)'I0LMAX=',i0lmax 
       STOP
    ENDIF
    
    CALL T2GRPH_ALLOCATE_1
    
    i1mlvl( 0:i0lmax+1) = i1mesh_level(  0:i0lmax+1)
    i1rdn2(-1:i0lmax  ) = i1rdiv_number(-1:i0lmax  )
    i1pdn2(-1:i0lmax  ) = 0
    d1rec(  0:i0lmax  ) = d1rec_tmp(  0:i0lmax  )
    
    DO i1=1,i0lmax
       i0mesh_level=i1mlvl(i1)-1
       i1pdn2(i1) = i0pdiv_number*(2**i0mesh_level)   
    ENDDO

    RETURN

  END SUBROUTINE T2INIT_GRPH

  !C------------------------------------------------------------------
  !C
  !C 
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2GRPH_SETUP
    INTEGER(i0ikind)::i0mlva,i0mlvb,i0mlvc,i1,j1

    !C
    !C CALCULATE NUMBER OF NODES IN EACH MESH LEVEL
    !C
    !C          I0NMAX1, i0NMAX2, I0NMAX3
    !C
    
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
    !C CALCULATE NUMBER OF ELEMENTS IN EACH MESH LEVEL
    !C 
    !C                  I0EMAX, I1EMAX
    !C
    
    i1emax(0:i0lmax)= 0
    DO i1=1,i0lmax
       DO j1=1,i1
          i1emax(i1)=i1emax(i1)+i1rdn2(j1)*i1pdn2(j1)
       ENDDO
    ENDDO
    
    i0emax=i1emax(i0lmax)

    !C
    !C CALCULATE NUMBER OF HANGED NODE
    !C 
    !C                  I0HMAX
    !C
    
    i0hmax=0

    DO i1=2,i0lmax
       i0mlva=i1mlvl(i1-1)
       i0mlvb=i1mlvl(i1)
       IF(i0mlva.NE.i0mlvb)i0hmax=i0hmax+i1pdn2(i1-1)
    ENDDO
    
    !C
    !C CALCULATE ARRAY SIZE OF MATRIX INDICES FOR NODE
    !C I0NRMX: 
    !C I0NCMX: 
    !C
    
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
             i0ncmx= i0ncmx + 1+i1pdn2(i1) + 7*i1pdn2(i1)
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             !C SECOND EDGE POINTS ON LEFT-SIDES IN Lv-2~LMAX MESH
             i0ncmx= i0ncmx + 8*i1pdn2(i1-1) + 9*i1pdn2(i1-1)
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
    
    
    !C CALCULATE ARRAY SIZE OF MATRIX INDICES FOR VARIABLE SPECIES
    !C I0VRMX: 
    !C I0VCMX: 
    SELECT CASE (i0spcs)
       
    CASE(0)
       i0vmax  = 1
       i0vcmx  = 1
       i0vrmx  = 2
    CASE(2:)
       i0vmax  = 8*i0spcs + 5
       i0vrmx  = i0vmax + 1
       i0vcmx  = 9*(i0spcs**2) + 48*i0spcs + 11
       !C
       i0vmcmx = 4 +  6*i0spcs
       !C
       i0vacmx = 7 + 37*i0spcs
       !C
       i0vdcmx =     20*i0spcs
       !C
       i0vgcmx = 4 + 15*i0spcs
       !C
       i0vecmx =     38*i0spcs + 8*(i0spcs**2)
       !C
       i0vscmx = 5 +  8*i0spcs
       
    END SELECT
    
    i0cmax = i0vcmx*i0ncmx
    i0xmax = i0vmax*i0nmax2
    
    CALL T2GRPH_ALLOCATE_2
    
    RETURN
    
  END SUBROUTINE T2GRPH_SETUP
  
  !C------------------------------------------------------------------
  !C
  !C        VARIABLE GRAPH GENERATING ROUTINE 
  !C                 FOR  MULTI-FLUID EQUATIONS FEM SOLVER
  !C        
  !C
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2GRPH_VGRAPH
    
    INTEGER(i0ikind)::&
         i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vgx,i0vofs
    
    !C VARIABLE-VARIABLE GRAPH 
    
    !C INITIALIZE
    i2vtbl(1:i0vmax,1:i0vmax)  = 0
    i1vidr(1:i0vrmx)           = 0
    i1vidc(1:i0vcmx)           = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i0vofs                     = 0
    i1vidr(i0vr)               = i0vofs+1
    
    i0vc                       = 1
    i0vg                       = i0vofs+1
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 5
    i0vg                       = i0vofs+2 
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    !C EQUATION FOR Bt
    i0vr                       = 2
    i0vofs                     = 2
    i1vidr(i0vr)               = i0vofs+1
    
    i0vc                       = 2 
    i0vg                       = i0vofs+1
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 3
    i0vg                       = i0vofs+2 
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 4
    i0vg                       = i0vofs+3 
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg



    !C EQUATION FOR Er
    i0vr                       = 3 
    i0vofs                     = 5
    i1vidr(i0vr)               = i0vofs+1

    i0vc                       = 3
    i0vg                       = i0vofs+1 
    i1vidc(i0vg)               = i0vc 
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 4
    i0vg                       = i0vofs+2 
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1-2
       i0vg                    = i0vofs+2+i1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO

    !C EQUATION FOR Ep
    i0vr                       = 4
    i0vofs                     = 7+i0spcs
    i1vidr(i0vr)               = i0vofs+1

    i0vc                       = 2
    i0vg                       = i0vofs+1
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
   
    i0vc                       = 4
    i0vg                       = i0vofs+2
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg
    
    DO i1 = 1, i0spcs
       i0vc                    = 8*i1
       i0vg                    = i0vofs+2 + 2*i1 -1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = 8*i1+1
       i0vg                    = i0vofs+2 + 2*i1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO
    
    !C EQUATION FOR Et
    i0vr                       = 5
    i0vofs                     = 9+3*i0spcs
    i1vidr(i0vr)               = i0vofs+1

    i0vc                       = 1
    i0vg                       = i0vofs+1 
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    i0vc                       = 5
    i0vg                       = i0vofs+2
    i1vidc(i0vg)               = i0vc
    i2vtbl(i0vr,i0vc)          = i0vg

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1+1
       i0vg                    = i0vofs+2 + i1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
    ENDDO
    
    DO i1 = 1,i0spcs
       i0vgx                   = 4*i0spcs+11+(i1-1)*(54+8*i0spcs)
       !C EQUATION FOR N
       i0vr                    = 8*i1-2
       i0vofs                  = i0vgx
       i1vidr(i0vr)            = i0vofs+1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr + 1
       i0vg                    = i0vofs + 2
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr + 2
       i0vg                    = i0vofs + 3
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr + 3
       i0vg                    = i0vofs + 4
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       
       !C EQUATION FOR Fr
       i0vofs                  = i0vgx  + 4
       i0vr                    = 8*i1   - 1
       i1vidr(i0vr)            = i0vofs + 1

       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr+1
       i0vg                    = i0vofs + 3
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr+2
       i0vg                    = i0vofs + 4
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr+3
       i0vg                    = i0vofs + 5
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       !C EQUATION FOR Fb
       i0vofs                  = i0vgx  + 9
       i0vr                    = 8*i1
       i1vidr(i0vr)            = i0vofs + 1
       DO j1=1,i0spcs
          
          i0vcx   = 8*j1
          
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs +  1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = 5
             i0vg              = i0vofs +  2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  -  2
             i0vg              = i0vofs +  3
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  1
             i0vg              = i0vofs +  4
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx
             i0vg              = i0vofs +  5
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  +  1
             i0vg              = i0vofs +  6
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  +  2
             i0vg              = i0vofs +  7
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  +  3
             i0vg              = i0vofs +  8
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  +  4
             i0vg              = i0vofs +  9
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  +  5
             i0vg              = i0vofs + 10
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs+10
          ELSE
             i0vc              = i0vcx
             i0vg              = i0vofs + 1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vofs = i0vofs+2
          END IF
       ENDDO

       !C EQUATION FOR Ft
       i0vofs                  = i0vgx  + 2*i0spcs + 17
       i0vr                    = 8*i1   + 1
       i1vidr(i0vr)            = i0vofs + 1

       DO j1=1,i0spcs
          i0vcx                = 8*j1
          IF(i1.EQ.j1)THEN
             i0vc              = 5
             i0vg              = i0vofs + 1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  - 1
             i0vg              = i0vofs + 3
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  
             i0vg              = i0vofs + 4
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 5
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 2
             i0vg              = i0vofs + 6
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  + 3
             i0vg              = i0vofs + 7
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 8
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  + 5
             i0vg              = i0vofs + 9
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 9
          ELSE
             i0vc              = i0vcx  + 1
             i0vg              = i0vofs + 1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  + 5
             i0vg              = i0vofs + 2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 2
          END IF
       ENDDO

       !C EQUATION FOR P
       i0vofs                  = i0vgx  + 4*i0spcs + 24
       i0vr                    = 8*i1   + 2
       i1vidr(i0vr)            = i0vofs + 1
       
       i0vc                    = i0vr-4
       i0vg                    = i0vofs + 1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr-3
       i0vg                    = i0vofs + 2
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr-2
       i0vg                    = i0vofs + 3
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr-1
       i0vg                    = i0vofs + 4
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr
       i0vg                    = i0vofs + 5
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr+1
       i0vg                    = i0vofs + 6
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr+2
       i0vg                    = i0vofs + 7
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr+3
       i0vg                    = i0vofs + 8
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg


       !C EQUATION FOR Qr
       i0vofs                  = i0vgx  + 4*i0spcs + 32
       i0vr                    = 8*i1   + 3
       i1vidr(i0vr)            = i0vofs + 1

       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr - 5
       i0vg                    = i0vofs + 3
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr - 1
       i0vg                    = i0vofs + 4
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr +1
       i0vg                    = i0vofs + 5
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg

       i0vc                    = i0vr + 2 
       i0vg                    = i0vofs + 6
       i1vidc(i0vg)            = i0vc
       i2vtbl(i0vr,i0vc)       = i0vg
      
       !C EQUATION FOR Qb
       i0vofs                  = i0vgx  + 4*i0spcs + 38
       i0vr                    = 8*i1   + 4
       i1vidr(i0vr)            = i0vofs + 1

       DO j1=1,i0spcs
          i0vcx   = 8*j1+4
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs +  1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = 5
             i0vg              = i0vofs +  2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  6
             i0vg              = i0vofs +  3
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  5
             i0vg              = i0vofs +  4
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  5
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  3
             i0vg              = i0vofs +  6
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  2
             i0vg              = i0vofs +  7
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  1
             i0vg              = i0vofs +  8
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  
             i0vg              = i0vofs +  9
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  +  1
             i0vg              = i0vofs + 10
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs+10
          ELSE
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  
             i0vg              = i0vofs +  2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs+2
          END IF
       ENDDO

       !C EQUATION FOR Qt
       i0vofs                  = i0vgx  + 6*i0spcs + 46
       i0vr                    = 8*i1   + 5
       i1vidr(i0vr)            = i0vofs + 1
       
       DO j1=1,i0spcs
          i0vcx   = 8*j1+5
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs +  1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = 5
             i0vg              = i0vofs +  2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  7
             i0vg              = i0vofs +  3
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  6
             i0vg              = i0vofs +  4
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  5
             i0vg              = i0vofs +  5
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  6
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  3
             i0vg              = i0vofs +  7
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  2
             i0vg              = i0vofs +  8
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx  -  1
             i0vg              = i0vofs +  9
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vc              = i0vcx 
             i0vg              = i0vofs + 10
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 10
          ELSE
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg
             
             i0vc              = i0vcx  
             i0vg              = i0vofs +  2
             i1vidc(i0vg)      = i0vc
             i2vtbl(i0vr,i0vc) = i0vg

             i0vofs = i0vofs + 2
          END IF
       ENDDO
    ENDDO
    
    i0vr                       = 5+ 8*i0spcs + 1
    i0vofs                     = 8*(i0spcs**2)+58*i0spcs+11
    i1vidr(i0vr)               = i0vofs+1

    RETURN
    
  END SUBROUTINE T2GRPH_VGRAPH
  
  !C
  !C 
  !C FOR MASS SUBMATRIX
  !C
  !C
  SUBROUTINE T2GRPH_VMGRAPH
    
    INTEGER(i0ikind)::&
         i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vgx,i0vofs
    
    !C INITIALIZE
    i1vmidr(1:i0vrmx )         = 0
    i1vmidc(1:i0vmcmx)         = 0
    
    !C EQUATION FOR Bp
    i0vr                       = 1
    i0vofs                     = 0
    i1vmidr(i0vr)              = i0vofs+1
    
    i0vc                       = 1
    i0vg                       = i0vofs+1
    i1vmidc(i0vg)              = i0vc

    !C EQUATION FOR Bt
    i0vr                       = 2
    i0vofs                     = 1
    i1vmidr(i0vr)              = i0vofs+1
    
    i0vc                       = 2 
    i0vg                       = i0vofs+1
    i1vmidc(i0vg)              = i0vc

    !C EQUATION FOR Er
    i0vr                       = 3 
    i0vofs                     = 2
    i1vmidr(i0vr)              = i0vofs+1


    !C EQUATION FOR Ep
    i0vr                       = 4
    i0vofs                     = 2
    i1vmidr(i0vr)              = i0vofs+1

    i0vc                       = 4
    i0vg                       = i0vofs+1
    i1vmidc(i0vg)              = i0vc

   
    !C EQUATION FOR Et
    i0vr                        = 5
    i0vofs                      = 3
    i1vmidr(i0vr)               = i0vofs+1

    i0vc                        = 5
    i0vg                        = i0vofs+1 
    i1vmidc(i0vg)               = i0vc
    
    DO i1 = 1,i0spcs
       i0vgx                   = 6*i1   - 2
       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i0vofs                  = i0vgx
       i1vmidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vmidc(i0vg)           = i0vc

      
       !C EQUATION FOR Fr
       i0vofs                  = i0vgx  + 1
       i0vr                    = 8*i1   - 1
       i1vmidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fb
       i0vofs                  = i0vgx  + 1
       i0vr                    = 8*i1
       i1vmidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vmidc(i0vg)           = i0vc

       !C EQUATION FOR Ft
       i0vofs                  = i0vgx  + 2
       i0vr                    = 8*i1   + 1
       i1vmidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vmidc(i0vg)           = i0vc

       !C EQUATION FOR P
       i0vofs                  = i0vgx  + 3
       i0vr                    = 8*i1   + 2
       i1vmidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vmidc(i0vg)           = i0vc


       !C EQUATION FOR Qr
       i0vofs                  = i0vgx  + 4
       i0vr                    = 8*i1   + 3
       i1vmidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qb
       i0vofs                  = i0vgx  + 4
       i0vr                    = 8*i1   + 4
       i1vmidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vmidc(i0vg)           = i0vc

       !C EQUATION FOR Qt
       i0vofs                  = i0vgx  + 5
       i0vr                    = 8*i1   + 5
       i1vmidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vmidc(i0vg)           = i0vc
    ENDDO
    
    i0vr                       = 5+ 8*i0spcs + 1
    i0vofs                     = 6*i0spcs+4
    i1vmidr(i0vr)              = i0vofs+1

    RETURN
    
  END SUBROUTINE T2GRPH_VMGRAPH
  
  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2GRPH_VAGRAPH

  INTEGER::&
       i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vgx,i0vofs

    !C VARIABLE-VARIABLE GRAPH 
    
    !C INITIALIZE
    i1vaidr(1:i0vrmx)           = 0
    i1vaidc(1:i0vacmx)          = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i0vofs                     = 0
    i1vaidr(i0vr)              = i0vofs+1
    
    i0vc                       = 1
    i0vg                       = i0vofs+1
    i1vaidc(i0vg)              = i0vc

    !C EQUATION FOR Bt
    i0vr                       = 2
    i0vofs                     = 1
    i1vaidr(i0vr)              = i0vofs+1
    
    i0vc                       = 2 
    i0vg                       = i0vofs+1
    i1vaidc(i0vg)              = i0vc

    !C EQUATION FOR Er
    i0vr                       = 3 
    i0vofs                     = 2
    i1vaidr(i0vr)              = i0vofs+1

    i0vc                       = 3
    i0vg                       = i0vofs+1
    i1vaidc(i0vg)              = i0vc

    i0vc                       = 4
    i0vg                       = i0vofs+2
    i1vaidc(i0vg)              = i0vc

    !C EQUATION FOR Ep
    i0vr                       = 4
    i0vofs                     = 4
    i1vaidr(i0vr)              = i0vofs+1

    i0vc                       = 1
    i0vg                       = i0vofs+2
    i1vaidc(i0vg)              = i0vc
    
    i0vc                       = 4
    i0vg                       = i0vofs+1
    i1vaidc(i0vg)              = i0vc
   
    !C EQUATION FOR Et
    i0vr                        = 5
    i0vofs                      = 5
    i1vaidr(i0vr)               = i0vofs+1

    i0vc                        = 1
    i0vg                        = i0vofs+1 
    i1vaidc(i0vg)               = i0vc

    i0vc                        = 5
    i0vg                        = i0vofs+2 
    i1vaidc(i0vg)               = i0vc
    
    DO i1 = 1,i0spcs
       i0vgx                   = 37*i1  - 30
       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i0vofs                  = i0vgx
       i1vaidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr
       i0vg                    = i0vofs + 1
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 2
       i1vaidc(i0vg)           = i0vc
       
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i1vaidc(i0vg)           = i0vc
       
       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 4
       i1vaidc(i0vg)           = i0vc
      
       !C EQUATION FOR Fr
       i0vofs                  = i0vgx  + 4
       i0vr                    = 8*i1   - 1
       i1vaidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fb
       i0vofs                  = i0vgx  + 4
       i0vr                    = 8*i1
       i1vaidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 1
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1vaidc(i0vg)           = i0vc
       
       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 5
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 6
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 5
       i0vg                    = i0vofs + 7
       i1vaidc(i0vg)           = i0vc

       !C EQUATION FOR Ft
       i0vofs                  = i0vgx  + 11
       i0vr                    = 8*i1   + 1
       i1vaidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 5
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 6
       i1vaidc(i0vg)           = i0vc

       !C EQUATION FOR P
       i0vofs                  = i0vgx  + 17
       i0vr                    = 8*i1   + 2
       i1vaidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 1
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 5
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 6
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 7
       i1vaidc(i0vg)           = i0vc


       !C EQUATION FOR Qr
       i0vofs                  = i0vgx  + 24
       i0vr                    = 8*i1   + 3
       i1vaidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Qb
       i0vofs                  = i0vgx  + 24
       i0vr                    = 8*i1   + 4
       i1vaidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 1
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 3
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 4
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 5
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 6
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 7
       i1vaidc(i0vg)           = i0vc


       !C EQUATION FOR Qt
       i0vofs                  = i0vgx  + 31
       i0vr                    = 8*i1   + 5
       i1vaidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 6
       i0vg                    = i0vofs + 1
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 2
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 3
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 4
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 5
       i1vaidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 6
       i1vaidc(i0vg)           = i0vc
    ENDDO
    
    i0vr                       =  8*i0spcs + 6
    i0vofs                     = 37*i0spcs + 7
    i1vaidr(i0vr)              = i0vofs    + 1

    RETURN

  END SUBROUTINE T2GRPH_VAGRAPH

  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2GRPH_VDGRAPH
      INTEGER::&
       i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vgx,i0vofs

    !C VARIABLE-VARIABLE GRAPH 
    
    !C INITIALIZE
    i1vdidr(1:i0vrmx)           = 0
    i1vdidc(1:i0vdcmx)          = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i0vofs                     = 0
    i1vdidr(i0vr)              = i0vofs+1
    

    !C EQUATION FOR Bt
    i0vr                       = 2
    i0vofs                     = 0
    i1vdidr(i0vr)              = i0vofs+1
    
    !C EQUATION FOR Er
    i0vr                       = 3 
    i0vofs                     = 0
    i1vdidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Ep
    i0vr                       = 4
    i0vofs                     = 0
    i1vdidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Et
    i0vr                        = 5
    i0vofs                      = 0
    i1vdidr(i0vr)               = i0vofs+1
    
    DO i1 = 1,i0spcs
       
       i0vgx                   = 20*i1  - 20
       
       !C EQUATION FOR N
       i0vr                    = 8*i1   - 2
       i0vofs                  = i0vgx
       i1vdidr(i0vr)           = i0vofs + 1
             
       !C EQUATION FOR Fr
       i0vr                    = 8*i1   - 1
       i0vofs                  = i0vgx 
       i1vdidr(i0vr)           = i0vofs + 1

       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i0vofs                  = i0vgx 
       i1vdidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 1
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 2
       i1vdidc(i0vg)           = i0vc
       
       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 3
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 4
       i0vg                    = i0vofs + 4
       i1vdidc(i0vg)           = i0vc


       !C EQUATION FOR Ft
       i0vr                    = 8*i1   + 1
       i0vofs                  = i0vgx  + 4
       i1vdidr(i0vr)           = i0vofs + 1

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 1
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 2
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 4
       i1vdidc(i0vg)           = i0vc

       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i0vofs                  = i0vgx  + 8
       i1vdidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 1
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 3
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 4
       i1vdidc(i0vg)           = i0vc

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   +  3
       i0vofs                  = i0vgx  + 12
       i1vdidr(i0vr)           = i0vofs +  1

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   +  4
       i0vofs                  = i0vgx  + 12
       i1vdidr(i0vr)           = i0vofs +  1
       
       i0vc                    = i0vr   - 6
       i0vg                    = i0vofs + 1
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 4
       i0vg                    = i0vofs + 2
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 3
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr 
       i0vg                    = i0vofs + 4
       i1vdidc(i0vg)           = i0vc

       !C EQUATION FOR Qt
       i0vr                    = 8*i1   +  5
       i0vofs                  = i0vgx  + 16
       i1vdidr(i0vr)           = i0vofs +  1
       
       i0vc                    = i0vr   - 7
       i0vg                    = i0vofs + 1
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 5
       i0vg                    = i0vofs + 2
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 3
       i1vdidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 4
       i1vdidc(i0vg)           = i0vc
    ENDDO
    
    i0vr                       =  8*i0spcs + 6
    i0vofs                     = 20*i0spcs 
    i1vdidr(i0vr)              = i0vofs    + 1
    
    RETURN
  END SUBROUTINE T2GRPH_VDGRAPH

  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2GRPH_VGGRAPH
      INTEGER::&
       i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vgx,i0vofs

    !C VARIABLE-VARIABLE GRAPH 
    
    !C INITIALIZE
    i1vgidr(1:i0vrmx)           = 0
    i1vgidc(1:i0vgcmx)          = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i0vofs                     = 0
    i1vgidr(i0vr)              = i0vofs+1
    
    i0vc                       = 5
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    
    !C EQUATION FOR Bt
    i0vr                       = 2
    i0vofs                     = 1
    i1vgidr(i0vr)              = i0vofs+1
    
    i0vc                       = 3
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc

    i0vc                       = 4
    i0vg                       = i0vofs + 2
    i1vgidc(i0vg)              = i0vc
    !C EQUATION FOR Er

    i0vr                       = 3 
    i0vofs                     = 3
    i1vgidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Ep
    i0vr                       = 4
    i0vofs                     = 3
    i1vgidr(i0vr)              = i0vofs+1

    i0vc                       = 2
    i0vg                       = i0vofs + 1
    i1vgidc(i0vg)              = i0vc
    !C EQUATION FOR Et
    i0vr                        = 5
    i0vofs                      = 4
    i1vgidr(i0vr)               = i0vofs+1
    
    DO i1 = 1,i0spcs
       
       i0vgx                   = 15*i1  - 11
       
       !C EQUATION FOR N
       i0vr                    = 8*i1   -  2
       i0vofs                  = i0vgx
       i1vgidr(i0vr)           = i0vofs +  1
             
       !C EQUATION FOR Fr
       i0vr                    = 8*i1   -  1
       i0vofs                  = i0vgx 
       i1vgidr(i0vr)           = i0vofs +  1

       i0vc                    = i0vr   +  3
       i0vg                    = i0vofs +  1
       i1vgidc(i0vg)           = i0vc

       !C EQUATION FOR Fb
       i0vr                    = 8*i1
       i0vofs                  = i0vgx  +  1
       i1vgidr(i0vr)           = i0vofs +  1
       
       i0vc                    = i0vr   -  2
       i0vg                    = i0vofs +  1
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs +  2
       i1vgidc(i0vg)           = i0vc
       
       i0vc                    = i0vr   +  2
       i0vg                    = i0vofs +  3
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr   +  4 
       i0vg                    = i0vofs +  4
       i1vgidc(i0vg)           = i0vc


       !C EQUATION FOR Ft
       i0vr                    = 8*i1   +  1
       i0vofs                  = i0vgx  +  5
       i1vgidr(i0vr)           = i0vofs +  1

       !C EQUATION FOR P
       i0vr                    = 8*i1   +  2
       i0vofs                  = i0vgx  +  5
       i1vgidr(i0vr)           = i0vofs +  1
       
       i0vc                    = i0vr   -  4
       i0vg                    = i0vofs +  1
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr   -  2
       i0vg                    = i0vofs +  2
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs +  3
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr   +  2
       i0vg                    = i0vofs +  4
       i1vgidc(i0vg)           = i0vc

       !C EQUATION FOR Qr
       i0vr                    = 8*i1   +  3
       i0vofs                  = i0vgx  +  9
       i1vgidr(i0vr)           = i0vofs +  1

       i0vc                    = i0vr   -  5
       i0vg                    = i0vofs +  1
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr   -  1
       i0vg                    = i0vofs +  2
       i1vgidc(i0vg)           = i0vc

       !C EQUATION FOR Qb
       i0vr                    = 8*i1   +  4
       i0vofs                  = i0vgx  + 11
       i1vgidr(i0vr)           = i0vofs +  1
       
       i0vc                    = i0vr   -  6
       i0vg                    = i0vofs +  1
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr   -  4
       i0vg                    = i0vofs +  2
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr   -  2
       i0vg                    = i0vofs +  3
       i1vgidc(i0vg)           = i0vc

       i0vc                    = i0vr 
       i0vg                    = i0vofs +  4
       i1vgidc(i0vg)           = i0vc

       !C EQUATION FOR Qt
       i0vr                    = 8*i1   +  5
       i0vofs                  = i0vgx  + 15
       i1vgidr(i0vr)           = i0vofs +  1
    ENDDO
    
    i0vr                       =  8*i0spcs + 6
    i0vofs                     = 15*i0spcs + 4
    i1vgidr(i0vr)              = i0vofs    + 1

    RETURN

  END SUBROUTINE T2GRPH_VGGRAPH

  !C
  !C
  !C
  !C
  !C

  SUBROUTINE T2GRPH_VEGRAPH
      INTEGER::&
       i1,j1,i0vr,i0vc,i0vg,i0vcx,i0vgx,i0vofs

    !C INITIALIZE
    i1veidr(1:i0vrmx)           = 0
    i1veidc(1:i0vecmx)           = 0

    !C EQUATION FOR Bp
    i0vr                       = 1
    i0vofs                     = 0
    i1veidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Bt
    i0vr                       = 2
    i0vofs                     = 0
    i1veidr(i0vr)              = i0vofs+1

    !C EQUATION FOR Er
    i0vr                       = 3 
    i0vofs                     = 0
    i1veidr(i0vr)              = i0vofs+1

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1-2
       i0vg                    = i0vofs+i1
       i1veidc(i0vg)           = i0vc
    ENDDO

    !C EQUATION FOR Ep
    i0vr                       = 4
    i0vofs                     = i0spcs
    i1veidr(i0vr)              = i0vofs+1
    DO i1 = 1, i0spcs
       i0vc                    = 8*i1
       i0vg                    = i0vofs + 2*i1 -1
       i1veidc(i0vg)           = i0vc
       
       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + 2*i1
       i1veidc(i0vg)           = i0vc
    ENDDO
    
    !C EQUATION FOR Et
    i0vr                       = 5
    i0vofs                     = 3*i0spcs
    i1veidr(i0vr)              = i0vofs+1

    DO i1 = 1, i0spcs
       i0vc                    = 8*i1+1
       i0vg                    = i0vofs + i1
       i1veidc(i0vg)           = i0vc
    ENDDO
    
    DO i1 = 1,i0spcs
       i0vgx                   = 4*i0spcs + (i1-1)*(32+8*i0spcs)
       !C EQUATION FOR N
       i0vr                    = 8*i1-2
       i0vofs                  = i0vgx
       i1veidr(i0vr)           = i0vofs+1
       
       !C EQUATION FOR Fr
       i0vofs                  = i0vgx 
       i0vr                    = 8*i1   - 1
       i1veidr(i0vr)           = i0vofs + 1

       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1veidc(i0vg)           = i0vc


       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1veidc(i0vg)           = i0vc


       i0vc                    = i0vr+1
       i0vg                    = i0vofs + 3
       i1veidc(i0vg)           = i0vc

       i0vc                    = i0vr+2
       i0vg                    = i0vofs + 4
       i1veidc(i0vg)           = i0vc

       !C EQUATION FOR Fb
       i0vofs                  = i0vgx  + 4
       i0vr                    = 8*i1
       i1veidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs
          
          i0vcx   = 8*j1
          
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs +  1
             i1veidc(i0vg)     = i0vc
             
             i0vc              = 5
             i0vg              = i0vofs +  2
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  -  1
             i0vg              = i0vofs +  3
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx
             i0vg              = i0vofs +  4
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  +  1
             i0vg              = i0vofs +  5
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  +  3
             i0vg              = i0vofs +  6
             i1veidc(i0vg)      = i0vc

             i0vc              = i0vcx  +  4
             i0vg              = i0vofs +  7
             i1veidc(i0vg)      = i0vc

             i0vc              = i0vcx  +  5
             i0vg              = i0vofs +  8
             i1veidc(i0vg)     = i0vc

             i0vofs = i0vofs+8
          ELSE
             i0vc              = i0vcx
             i0vg              = i0vofs + 1
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1veidc(i0vg)     = i0vc
             
             i0vofs = i0vofs+2
          END IF
       ENDDO

       !C EQUATION FOR Ft
       i0vofs                  = i0vgx  + 2*i0spcs + 10
       i0vr                    = 8*i1   + 1
       i1veidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs
          i0vcx                = 8*j1 + 1
          IF(i1.EQ.j1)THEN
             i0vc              = 5
             i0vg              = i0vofs + 1
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  - 2
             i0vg              = i0vofs + 2
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  
             i0vg              = i0vofs + 3
             i1veidc(i0vg)     = i0vc
             
             
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 4
             i1veidc(i0vg)     = i0vc

             i0vofs = i0vofs + 4
          ELSE
             i0vc              = i0vcx 
             i0vg              = i0vofs + 1
             i1veidc(i0vg)     = i0vc
             
             i0vc              = i0vcx  + 4
             i0vg              = i0vofs + 2
             i1veidc(i0vg)     = i0vc

             i0vofs = i0vofs + 2
          END IF
       ENDDO

       !C EQUATION FOR P
       i0vr                    = 8*i1   + 2
       i0vofs                  = i0vgx  + 4*i0spcs + 12
       i1veidr(i0vr)           = i0vofs + 1
       
       i0vc                    = i0vr   - 3
       i0vg                    = i0vofs + 1
       i1veidc(i0vg)           = i0vc


       i0vc                    = i0vr   - 2
       i0vg                    = i0vofs + 2
       i1veidc(i0vg)           = i0vc

       i0vc                    = i0vr   - 1
       i0vg                    = i0vofs + 3
       i1veidc(i0vg)           = i0vc

       i0vc                    = i0vr
       i0vg                    = i0vofs + 4
       i1veidc(i0vg)           = i0vc


       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 5
       i1veidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 2
       i0vg                    = i0vofs + 6
       i1veidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 3
       i0vg                    = i0vofs + 7
       i1veidc(i0vg)           = i0vc



       !C EQUATION FOR Qr
       i0vr                    = 8*i1   + 3
       i0vofs                  = i0vgx  + 4*i0spcs + 19
       i1veidr(i0vr)           = i0vofs + 1

       i0vc                    = 3
       i0vg                    = i0vofs + 1
       i1veidc(i0vg)           = i0vc


       i0vc                    = 4
       i0vg                    = i0vofs + 2
       i1veidc(i0vg)           = i0vc

       i0vc                    = i0vr   + 1
       i0vg                    = i0vofs + 3
       i1veidc(i0vg)           = i0vc


       i0vc                    = i0vr   + 2 
       i0vg                    = i0vofs + 4
       i1veidc(i0vg)           = i0vc

      
       !C EQUATION FOR Qb
       i0vr                    = 8*i1   + 4
       i0vofs                  = i0vgx  + 4*i0spcs + 23
       i1veidr(i0vr)           = i0vofs + 1

       DO j1=1,i0spcs

          i0vcx   = 8*j1+4

          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs +  1
             i1veidc(i0vg)     = i0vc
          
             i0vc              = 5
             i0vg              = i0vofs +  2
             i1veidc(i0vg)     = i0vc
          
             i0vc              = i0vcx  -  5
             i0vg              = i0vofs +  3
             i1veidc(i0vg)     = i0vc
          

             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  4
             i1veidc(i0vg)     = i0vc
          

             i0vc              = i0vcx  -  3
             i0vg              = i0vofs +  5
             i1veidc(i0vg)     = i0vc
             
             i0vc              = i0vcx  -  1
             i0vg              = i0vofs +  6
             i1veidc(i0vg)     = i0vc
          

             i0vc              = i0vcx  
             i0vg              = i0vofs +  7
             i1veidc(i0vg)     = i0vc
             
             i0vc              = i0vcx  +  1
             i0vg              = i0vofs +  8
             i1veidc(i0vg)     = i0vc
             
             i0vofs = i0vofs + 8
          ELSE
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1veidc(i0vg)     = i0vc
             
             i0vc              = i0vcx  
             i0vg              = i0vofs +  2
             i1veidc(i0vg)     = i0vc
             
             i0vofs = i0vofs + 2
          END IF
       ENDDO

       !C EQUATION FOR Qt
       i0vofs                  = i0vgx  + 6*i0spcs + 29
       i0vr                    = 8*i1   + 5
       i1veidr(i0vr)           = i0vofs + 1
       
       DO j1=1,i0spcs
          i0vcx   = 8*j1+5
          IF(i1.EQ.j1)THEN
             i0vc              = 4
             i0vg              = i0vofs +  1
             i1veidc(i0vg)     = i0vc
             
             i0vc              = 5
             i0vg              = i0vofs +  2
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  3
             i1veidc(i0vg)     = i0vc

             i0vc              = i0vcx  -  2
             i0vg              = i0vofs +  4
             i1veidc(i0vg)      = i0vc

             i0vc              = i0vcx  
             i0vg              = i0vofs +  5
             i1veidc(i0vg)     = i0vc

             i0vofs = i0vofs + 5
          ELSE
             i0vc              = i0vcx  -  4
             i0vg              = i0vofs +  1
             i1veidc(i0vg)     = i0vc

             
             i0vc              = i0vcx  
             i0vg              = i0vofs +  2
             i1veidc(i0vg)     = i0vc

             i0vofs = i0vofs + 2
          END IF
       ENDDO
    ENDDO
    
    i0vr                       = 5+ 8*i0spcs + 1
    i0vofs                     = 8*(i0spcs**2)+36*i0spcs
    i1veidr(i0vr)              = i0vofs + 1

    RETURN

  END SUBROUTINE T2GRPH_VEGRAPH

  !C------------------------------------------------------------------
  !C
  !C            NODE GRAPH GENERATING ROUTINE 
  !C                 FOR  MULTI-FLUID EQUATIONS FEM SOLVER
  !C             
  !C             COMPRESSED ROW STRAGE FORMAT
  !C
  !C      FOR LINEAR RECTANGULAR ELEMENT   4------3 
  !C                                       |      |
  !C                                       |      |
  !C                                       1------2
  !C
  !C------------------------------------------------------------------  
  SUBROUTINE T2GRPH_NGRAPH1
    INTEGER(i0ikind)::&
         i0ecnt,i0hcnt,i0rcnt,i0ccnt,i0ncnt,&
         i0ll,i0lr,i0ul,i0ur,&
         i0ppc1,i0ppl2,i0ppc2,i0ppr2,&
         i0rdc1,i0rdc2,&
         i0jl,i0jc,i0jr,i0il,&
         i0stc1,i0nd,i0nc,i0nu,&
         i0stl2,i0stc2,i0str2,&
         i0mlva,i0mlvb,i0mlvc,&
         i1tmp06(6),i1tmp07(7),i1tmp08(8),i1tmp09(9),i1tmp11(11),&
         i1subtot(0:i0lmax),i0offset
    INTEGER(i0ikind)::&
         i1,j1,i2,j2,i3
    !C--------------------------------------------
    
    !C CONSTRUCT NON-DEGENERATED NODE
    !C                    - DEGENERATED NODE GRAPH
    
    !C SUBTOTALS UP TO Nth-DOMAIN: I1SUBTOT[N]    
    DO i1=0,i0lmax
       i1subtot(i1)=1
    ENDDO
    DO i1=1,i0lmax
       DO j1=1,i1
          i1subtot(i1)=i1subtot(i1)+i1rdn2(j1)*i1pdn2(j1)
       ENDDO
    ENDDO
    
    !C SET NODE NUMBER 
    i0ncnt=0
    i0hcnt=0
    DO i1=1,i0lmax
       i0mlva=i1mlvl(i1-1)
       i0mlvb=i1mlvl(i1)
       i0ppc2 = i1pdn2(i1)
       i0ppl2 = INT(i1pdn2(i1)/2)
       i0ppc1 = i1pdn1(i1)
       i0rdc1 = i1rdn1(i1)
       
       !C SET OFFSET
       i0offset = 0
       IF(i1mlvl(i1).GE.3)THEN
          DO i2=0,i1mlvl(i1)-3
             i0offset = i0offset+i1pdn2(1)*(2**i2)
          ENDDO
       END IF
       
       !C SET SUBTOTALS
       i0stl2  = i1subtot(i1-1) - i1pdn2(i1-1)
       i0stc2  = i1subtot(i1-1)
       i0stm2  = i1subtot(i0lmax)
       IF(i0mlva.EQ.0)THEN
          !C DOMAIN WITH AXIS BOUNDARY
          DO i2=1,i0rdc1
             IF(i2.NE.1)THEN
                !C NODES EXCEPT AXIS BOUNADRY 
                i0il=i2-2
                DO j2=1,i0ppc1
                   i0ncnt=i0ncnt+1
                   i2crt(i0ncnt,1)= i0ncnt
                   i2crt(i0ncnt,2)&
                        = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                ENDDO
             ELSEIF(    i2.EQ.1)THEN
                !C NODES ON AXIS BOUNADRY 
                DO j2=1,i0ppc1
                   i0ncnt=i0ncnt+1
                   i2crt(i0ncnt,1)= i0ncnt
                   i2crt(i0ncnt,2)= i0stc2
                ENDDO
             ENDIF
          ENDDO
       ELSEIF(i0mlva.NE.i0mlvb)THEN
          !C DOMAIN WITH INTERFACE BOUNDARY
          DO i2=1,i0rdc1
             IF(    i2.NE.1)THEN
                !C NODES EXCEPT INTERFACE BOUNADRY 
                i0il=i2-2
                DO j2=1,i0ppc1
                   i0ncnt=i0ncnt+1
                   i2crt(i0ncnt,1)= i0ncnt
                   i2crt(i0ncnt,2)&
                        = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                ENDDO
             ELSEIF(i2.EQ.1)THEN
                !C NODES ON INTERFACE BOUNDARY
                DO j2=1,i0ppc1
                   i0ncnt=i0ncnt+1
                   IF(MOD(j2,2).EQ.1)THEN
                      i2crt(i0ncnt,1)= i0ncnt
                      i2crt(i0ncnt,2)&
                           = i0stl2 + MOD(INT((j2+1)/2)-1,i0ppl2)+1
                   ELSEIF(MOD(j2,2).EQ.0)THEN
                      i0hcnt=i0hcnt+1
                      i2crt(i0ncnt,1)= i0ncnt
                      i2crt(i0ncnt,2)&
                           = i0stm2 + i0offset + INT(j2/2)
                      i2hbc2(i0hcnt,1)&
                           = i0stl2 + MOD(INT((j2-1)/2),i0ppl2)+1
                      i2hbc2(i0hcnt,2)&
                           = i0stl2 + MOD(INT((j2+1)/2),i0ppl2)+1
                   ENDIF
                ENDDO
             END IF
          ENDDO
       ELSEIF(i0mlva.EQ.i0mlvb)THEN
          !C DOMAIN WITHOUT INTERFACE BOUNDARY
          DO i2=1,i0rdc1
             i0il=i2-2
             DO j2=1,i0ppc1
                i0ncnt=i0ncnt+1
                i2crt(i0ncnt,1)= i0ncnt
                i2crt(i0ncnt,2)&
                     = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
             ENDDO
          ENDDO
       ENDIF
    ENDDO


    !C CONSTRUCT ELEMENT - NODE GRAPH
    !C LAST UPDATE  2013/02/13
    i0ecnt = 0
    DO i1=1,i0lmax
       i0ppc1 = i1pdn1(i1)
       i0ppc2 = i1pdn2(i1)
       i0rdc2 = i1rdn2(i1)
       i0stc1=0
       DO i2=0,i1-1
          i0stc1 = i0stc1+i1rdn1(i2)*i1pdn1(i2)
       ENDDO
       DO i2=1,i0rdc2
       DO j2=1,i0ppc2
          i0ecnt = i0ecnt+1
          i0ll =i0stc1+ i0ppc1*(i2-1)+j2
          i0lr =i0stc1+ i0ppc1*i2    +j2
          i0ur =i0stc1+ i0ppc1*i2    +j2+1
          i0ul =i0stc1+ i0ppc1*(i2-1)+j2+1
          
          i3enr(i0ecnt,1,1) = i2crt(i0ll,1)
          i3enr(i0ecnt,1,2) = i2crt(i0lr,1)
          i3enr(i0ecnt,1,3) = i2crt(i0ur,1)
          i3enr(i0ecnt,1,4) = i2crt(i0ul,1)

          i3enr(i0ecnt,2,1) = i2crt(i0ll,2)
          i3enr(i0ecnt,2,2) = i2crt(i0lr,2)
          i3enr(i0ecnt,2,3) = i2crt(i0ur,2)
          i3enr(i0ecnt,2,4) = i2crt(i0ul,2)
       ENDDO
       ENDDO
    ENDDO
    
    !C CONSTRUCT NODE - NODE RELATION TABLE
    !C CHEKED 2012/11/14
    i0ccnt = 0
    i0rcnt = 0
    DO i1=1,i0lmax
       i0mlva=i1mlvl(i1-1)
       i0mlvb=i1mlvl(i1)
       i0mlvc=i1mlvl(i1+1)
       DO i2=1,i1rdn1(i1)
          IF(    (i0mlva.EQ.0).AND.(i2.EQ.1))THEN
             !C PATTERN 1
             i0rcnt=i0rcnt+1
             i1nidr(i0rcnt)=1
             i0ppc2=i1pdn2(i1)
             DO i0jc=1,i0ppc2+1
                i0ccnt=i0ccnt+1
                i1nidc(i0ccnt)=i0jc
             ENDDO
             i0rcnt          = i0rcnt+1
             i1nidr(i0rcnt) = i0ccnt+1    
          ELSEIF((i0mlva.EQ.0).AND.(i2.EQ.2))THEN
             !C PATTERN 2
             i0stc2 = 1
             i0ppc2 = i1pdn2(i1)
             DO i0jc = 1,i0ppc2
                i0nd = i0stc2 + i0ppc2 - MOD(i0ppc2+1-i0jc,i0ppc2)
                i0nc = i0stc2 + i0jc 
                i0nu = i0stc2 + MOD(i0jc,i0ppc2) + 1
                i1tmp07(1) = 1
                i1tmp07(2) = i0nd
                i1tmp07(3) = i0nc
                i1tmp07(4) = i0nu 
                i1tmp07(5) = i0nd + i0ppc2
                i1tmp07(6) = i0nc + i0ppc2
                i1tmp07(7) = i0nu + i0ppc2
                CALL SORT_ARRAY( 7,i1tmp07)
                DO i3=1,7
                   i0ccnt         = i0ccnt+1
                   i1nidc(i0ccnt)= i1tmp07(i3)
                ENDDO
                i0rcnt        = i0rcnt+1
                i1nidr(i0rcnt)= i0ccnt+1
             ENDDO
          ELSEIF((i0mlvb.NE.i0mlvc).AND.&
                 (i0mlvc.NE.0).AND.&
                 (i2.EQ.i1rdn1(i1)))THEN
             !C PATTERN 4
             i0ppc2 = i1pdn2(i1)
             i0ppr2 = i1pdn2(i1+1)
             i0stc2 = i1subtot(i1-1)+i0ppc2*(i2-2)
             i0str2 = i1subtot(i1)
             DO i0jc=1,i0ppc2
                i0jr = 2*i0jc-1
                i0nd = i0stc2 + i0ppc2 - MOD(i0ppc2+1-i0jc,i0ppc2) 
                i0nc = i0stc2 + i0jc
                i0nu = i0stc2 + MOD(i0jc,i0ppc2) + 1
                i1tmp11( 1) = i0nd - i0ppc2
                i1tmp11( 2) = i0nc - i0ppc2
                i1tmp11( 3) = i0nu - i0ppc2
                i1tmp11( 4) = i0nd
                i1tmp11( 5) = i0nc
                i1tmp11( 6) = i0nu
                i1tmp11( 7) = i0str2 + i0ppr2 - MOD(i0ppr2+2-i0jr,i0ppr2)  
                i1tmp11( 8) = i0str2 + i0ppr2 - MOD(i0ppr2+1-i0jr,i0ppr2)
                i1tmp11( 9) = i0str2 + i0jr
                i1tmp11(10) = i0str2 + MOD(i0jr  ,i0ppr2) + 1
                i1tmp11(11) = i0str2 + MOD(i0jr+1,i0ppr2) + 1
                CALL SORT_ARRAY(11,i1tmp11)
                DO i3=1,11
                   i0ccnt         = i0ccnt+1
                   i1nidc(i0ccnt) = i1tmp11(i3)
                ENDDO
                i0rcnt         = i0rcnt+1
                i1nidr(i0rcnt)= i0ccnt+1
             ENDDO
          ELSEIF((i0mlva.NE.i0mlvb).AND.&
                 (i0mlva.NE.0).AND.&
                 (i2.EQ.2))THEN
             i0ppc2 = i1pdn2(i1)
             i0ppl2 = i1pdn2(i1-1)
             i0stc2 = i1subtot(i1-1)
             i0stl2 = i0stc2 -i0ppl2
             DO i0jc=1,i0ppc2
                IF(    MOD(i0jc,2).EQ.0)THEN
                   !C PATTERN 5
                   i0jl = INT(i0jc/2)
                   i0nd = i0stc2 + i0ppc2 - MOD(i0ppc2+1-i0jc,i0ppc2)
                   i0nc = i0stc2 + i0jc
                   i0nu = i0stc2 + MOD(i0jc,i0ppc2) + 1
                   i1tmp08(1) = i0stl2 + i0jl
                   i1tmp08(2) = i0stl2 + MOD(i0jl,i0ppl2) + 1
                   i1tmp08(3) = i0nd
                   i1tmp08(4) = i0nc
                   i1tmp08(5) = i0nu
                   i1tmp08(6) = i0nd + i0ppc2
                   i1tmp08(7) = i0nc + i0ppc2
                   i1tmp08(8) = i0nu + i0ppc2
                   CALL SORT_ARRAY( 8,i1tmp08)
                   DO i3=1,8
                      i0ccnt          = i0ccnt+1
                      i1nidc(i0ccnt) = i1tmp08(i3)
                   ENDDO
                   i0rcnt          = i0rcnt+1
                   i1nidr(i0rcnt) = i0ccnt+1
                ELSEIF(MOD(i0jc,2).EQ.1)THEN
                   !C PATTERN 6
                   i0jl = INT((i0jc+1)/2)
                   i0nd = i0stc2 + i0ppc2 - MOD(i0ppc2+1-i0jc,i0ppc2)
                   i0nc = i0stc2 + i0jc
                   i0nu = i0stc2 + MOD(i0jc,i0ppc2) + 1
                   i1tmp09(1) = i0stl2&
                        + i0ppl2 - MOD(i0ppl2+1-i0jl,i0ppl2)
                   i1tmp09(2) = i0stl2 + i0jl
                   i1tmp09(3) = i0stl2 + MOD(i0jl,i0ppl2) + 1
                   i1tmp09(4) = i0nd
                   i1tmp09(5) = i0nc
                   i1tmp09(6) = i0nu
                   i1tmp09(7) = i0nd + i0ppc2
                   i1tmp09(8) = i0nc + i0ppc2
                   i1tmp09(9) = i0nu + i0ppc2
                   CALL SORT_ARRAY( 9,i1tmp09)
                   DO i3=1,9
                      i0ccnt          = i0ccnt+1
                      i1nidc(i0ccnt) = i1tmp09(i3)
                   ENDDO
                   i0rcnt         = i0rcnt+1
                   i1nidr(i0rcnt) = i0ccnt+1
                ENDIF
             ENDDO
          ELSEIF((i0mlvb.NE.i0mlvc).AND.&
                 (i0mlvc.EQ.0).AND.&
                 (i2.EQ.i1rdn1(i1)))THEN
             ! PATTERN 7
             i0ppc2 = i1pdn2(i1)
             i0stc2 = i1subtot(i1-1)+i0ppc2*(i2-2)
             DO i0jc=1,i0ppc2
                i0nd = i0stc2 + i0ppc2 - MOD(i0ppc2+1-i0jc,i0ppc2)
                i0nc = i0stc2 + i0jc
                i0nu = i0stc2 + MOD(i0jc,i0ppc2) + 1
                i1tmp06(1) = i0nd - i0ppc2
                i1tmp06(2) = i0nc - i0ppc2
                i1tmp06(3) = i0nu - i0ppc2
                i1tmp06(4) = i0nd
                i1tmp06(5) = i0nc
                i1tmp06(6) = i0nu
                CALL SORT_ARRAY( 6,i1tmp06)
                DO i3=1,6
                   i0ccnt          = i0ccnt+1
                   i1nidc(i0ccnt) = i1tmp06(i3)
                ENDDO
                i0rcnt         = i0rcnt+1
                i1nidr(i0rcnt)= i0ccnt+1
             ENDDO
          ELSEIF((i0mlva.NE.0).AND.(i2.EQ.1))THEN
             
          ELSE
             !C PATTERN 3
             i0ppc2 = i1pdn2(i1)
             i0stc2 = i1subtot(i1-1)+i0ppc2*(i2-2)
             DO i0jc=1,i0ppc2
                i0nd = i0stc2 + i0ppc2 - MOD(i0ppc2+1-i0jc,i0ppc2)
                i0nc = i0stc2 + i0jc
                i0nu = i0stc2 + MOD(i0jc,i0ppc2) + 1
                i1tmp09(1) = i0nd - i0ppc2
                i1tmp09(2) = i0nc - i0ppc2
                i1tmp09(3) = i0nu - i0ppc2
                i1tmp09(4) = i0nd
                i1tmp09(5) = i0nc
                i1tmp09(6) = i0nu
                i1tmp09(7) = i0nd + i0ppc2
                i1tmp09(8) = i0nc + i0ppc2
                i1tmp09(9) = i0nu + i0ppc2
                CALL SORT_ARRAY( 9,i1tmp09)
                DO i3=1,9
                   i0ccnt          = i0ccnt+1
                   i1nidc(i0ccnt) = i1tmp09(i3)
                ENDDO
                i0rcnt         = i0rcnt+1
                i1nidr(i0rcnt)= i0ccnt+1
             ENDDO
          ENDIF
       ENDDO
    ENDDO
        
    RETURN
 
  END SUBROUTINE T2GRPH_NGRAPH1

  !C------------------------------------------------------------------
  !C
  !C            NODE GRAPH GENERATING ROUTINE 
  !C                 FOR  MULTI-FLUID EQUATIONS FEM SOLVER
  !C             
  !C             COMPRESSED ROW STRAGE FORMAT
  !C
  !C      FOR QUADRADICR RECTANGULAR ELEMENT  
  !C                                            7---6---5 
  !C                                            |       |
  !C                                            8       4
  !C                                            |       |
  !C                                            1---2---3
  !C      (WILL BE IMPLEMENTATED)  
  !C------------------------------------------------------------------  
  SUBROUTINE T2GRPH_NGRAPH2
    RETURN
  END SUBROUTINE T2GRPH_NGRAPH2

  !C------------------------------------------------------------------
  !C
  !C            NODE GRAPH GENERATING ROUTINE 
  !C                 FOR  MULTI-FLUID EQUATIONS FEM SOLVER
  !C             
  !C             COMPRESSED ROW STRAGE FORMAT
  !C
  !C      FOR CUBIC RECTANGULAR ELEMENT  
  !C                                            10--09--08--07
  !C                                            |            |
  !C                                            11          06 
  !C                                            |            |
  !C                                            12          05
  !C                                            |            |
  !C                                            01--02--03--04
  !C      (WILL BE IMPLEMENTATED)  
  !C------------------------------------------------------------------  
  SUBROUTINE T2GRPH_NGRAPH3
    RETURN
  END SUBROUTINE T2GRPH_NGRAPH3
  
  
  SUBROUTINE T2GRPH_MFCS
    INTEGER(i0ikind)::&
         i1,i2,j2,&
         i0ncnt1,i0rdn1,i0pdn1,i0rdn2,i0pdn2
    REAL(   i0rkind)::&
         d0rsiz,d0psiz
    REAL(i0rkind),DIMENSION(:),ALLOCATABLE::&
         d1mcr1,d1mcp1
    !C------------------------------------------------------     
    
    !C SET MAGNETIC COORDINATES
    i0ncnt1 = 0
    DO i1=1,i0lmax
       
       i0rdn1=i1rdn1(i1)
       i0pdn1=i1pdn1(i1)
       i0rdn2=i1rdn2(i1)
       i0pdn2=i1pdn2(i1)
       
       ALLOCATE(d1mcr1(i0rdn1),d1mcp1(i0pdn1))

       DO i2=1,i0rdn1
          d1mcr1(i2) = 0.D0
       ENDDO
       DO i2=1,i0pdn1
          d1mcp1(i2) = 0.D0
       ENDDO
       
       d1rsiz(i1) = (d1rec(i1)-d1rec(i1-1))/DBLE(i0rdn2)
       d1psiz(i1) = 1.d0/DBLE(i0pdn2)
       d1msiz(i1) = SQRT(d1rsiz(i1)*d1psiz(i1))
       
       d0rsiz=d1rsiz(i1)
       d0psiz=d1psiz(i1)
       
       DO i2= 1,i0rdn1
          d1mcr1(i2) = d0rsiz*DBLE(i2-1)+d1rec(i1-1)
       ENDDO
       DO i2=1,i0pdn1
          d1mcp1(i2) = d0psiz*DBLE(i2-1)
       ENDDO
    
       !C CONSTRUCT GEOMETRY
       DO j2=1,i0rdn1 
       DO i2=1,i0pdn1
          i0ncnt1=i0ncnt1+1
          !C MAGNETIC FLUX COORDINATE (RHO, CHI)
          d2mfc1(i0ncnt1,1)=d1mcr1(j2)
          d2mfc1(i0ncnt1,2)=d1mcp1(i2)
       ENDDO
       ENDDO
       
       DEALLOCATE(d1mcr1,d1mcp1)
   
    ENDDO
    
    IF(i0ncnt1.NE.i0nmax1)THEN
       WRITE(6,*)'CHECK SUM ERROR'
       STOP
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2GRPH_MFCS

  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2GRPH_OUTPUT
    INTEGER(i0ikind)::i1
    
    OPEN(10,FILE='I2CRT_TEST.dat')
    DO i1=1,i0nmax1
       WRITE(10,'("NUM1=",I9,1X,"CRT1=",I9,1X,"CRT2=",I9)')&
            i1,i2crt(i1,1),i2crt(i1,2)
    ENDDO
     
    CLOSE(10)
    OPEN(10,FILE='I2HBC2_TEST.dat')
    DO i1=1,i0hmax
       WRITE(10,'("HNUM=",I9,1X,"LHN=",I9,1X,"UHN=",I9)')&
            i1,i2hbc2(i1,1),i2hbc2(i1,2)
    ENDDO

    CLOSE(10)
    OPEN(10,FILE='I3ENR1_TEST.dat')
    DO i1=1,i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i1,i3enr(i1,1,1),i3enr(i1,1,2),i3enr(i1,1,3),i3enr(i1,1,4)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR2_TEST.dat')
    DO i1=1,i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i1,i3enr(i1,2,1),i3enr(i1,2,2),i3enr(i1,2,3),i3enr(i1,2,4)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='CRS_TEST_INDR.dat')
    DO i1=1,i0nrmx
       IF(i1.EQ.1)THEN
          WRITE(10,'("I1=",I9,1X,"INDNR=",I9)')i1,i1nidr(i1)
       ELSE
          WRITE(10,'("I1=",I9,1X,"INDNR=",I9,1X,"J1=",I9,1X,"NOC=",I9)')&
               i1,i1nidr(i1),i1-1,i1nidr(i1)-i1nidr(i1-1)
       ENDIF
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='CRS_TEST_INDC.dat')
    DO i1=1,i0ncmx
       WRITE(10,'("I1=",I9,1X,"INDNC=",I9)')i1,i1nidc(i1)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_NMAX.dat')
    DO i1=1,i0lmax
       WRITE(10,'("MESH_NUMBER=",I3,1X,"MESH_LEVEL=",I3)')&
            i1,i1mlvl(i1)
       WRITE(10,'("NMAX1=",I9,1X,"NMAX2=",I9)')&
            i1nmax1(i1),i1nmax2(i1)
    ENDDO
    WRITE(10,'("TOTAL_NMAX1=",I9,1X,"TOTAL_NMAX2=",I9)')&
         i0nmax1,i0nmax2
    CLOSE(10)
    OPEN(10,FILE='TEST_EMAX.dat')
    DO i1=0,i0lmax
       WRITE(10,'("MESH_NUMBER=",I3,1X,"EMAX=",I9)')i1,i1emax(i1)
    ENDDO
    WRITE(10,'("TOTAL_EMAX1=",I9)')i0emax
    CLOSE(10)
    !C CALCURATE NUMBER OF HANGED-NODE
    OPEN(10,FILE='MATRIX_INFO1.dat')
    WRITE(10,'("I1NIDR_ARRAY_SIZE=",I9)')i0nrmx
    WRITE(10,'("I1NIDC_ARRAY_SIZE=",I9)')i0ncmx
    WRITE(10,'("D1GSM2_ARRAY_SIZE=",I9)')i0cmax
    WRITE(10,'("D1GUV2_ARRAY_SIZE=",I9)')i0nmax3
    CLOSE(10)
    
    OPEN(10,FILE='MFC1_CHECK.dat')
    DO i1=1,i0nmax1
       WRITE(10,'("i1=",I5,1X,"RHO=",D15.8,1X,"CHI=",D15.8)')&
            i1,d2mfc1(i1,1),d2mfc1(i1,2)
    ENDDO
    CLOSE(10)
    
    RETURN
  END SUBROUTINE T2GRPH_OUTPUT

  !C------------------------------------------------------------------
  !C
  !C        ROUTINE FOR TERMINATE UNNECESSARY WORKING ARRAYS
  !C
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2GRPH_TERM
    RETURN
  END SUBROUTINE T2GRPH_TERM

  !C------------------------------------------------------------------
  !C 
  !C         BUBLE SORT ALGORITHM
  !C 
  !C------------------------------------------------------------------
  SUBROUTINE SORT_ARRAY(i0asiz, i1array)
    INTEGER(i0ikind),INTENT(IN)::i0asiz
    INTEGER(i0ikind),DIMENSION(i0asiz),INTENT(INOUT)::i1array
    INTEGER(i0ikind)::&
         i0mtmp, i0mloc,i1
    
    DO i1 = 1, i0asiz-1
       i0mtmp = minval(i1array(i1+1:i0asiz))
       i0mloc = minloc(i1array(i1+1:i0asiz),1)+i1
       IF(i1array(i1).GT.i0mtmp)THEN
          i1array(i0mloc) = i1array(i1)
          i1array(i1)     = i0mtmp
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE SORT_ARRAY
END MODULE T2GRPH
