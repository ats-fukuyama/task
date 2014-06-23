!------------------------------------------------------------------
!
!         MODULE T2NGRA
!       
!         NODE GRAPH GENERATOR FOR TASK/T2
! 
!                   LAST UPDATE 2014-06-20 H.Seto
!
!    NodeRowCRS NodeColCRS NodeDiaCRS ElementNodeGraph  
!    HangedNodeTable
!
!    i1eidr, i1eidc i2crt i1mc1d,d1mc1d,i1mfc1 
!   
!------------------------------------------------------------------

MODULE T2NGRA
  
  USE T2CNST,ONLY:&
       ikind,rkind
  
  IMPLICIT NONE
  
  PUBLIC T2_NGRA,T2_NGRA_OUTPUT
  
  PRIVATE
  
CONTAINS 
  
  !------------------------------------------------------------------
  !
  !          T2_NGRA
  !
  !          MAIN SUBROUTINE OF T2NGRA
  !
  !------------------------------------------------------------------
  SUBROUTINE T2_NGRA
    
    USE T2COMM,ONLY:NNMAX,idfile
    
    CALL T2NGRA_COORDINATE
    
    SELECT CASE (NNMAX)
       
    CASE( 4) !C LINIEAR   ELEMENT
       CALL T2NGRA_NGRAPH1
    CASE( 9)
       CALL T2NGRA_NGRAPH2
    CASE(16) !C CUBIC     ELEMENT
       CALL T2NGRA_NGRAPH3
    END SELECT
    
    IF(idfile.ge.4) CALL T2_NGRA_OUTPUT
    CALL T2_NGRA_OUTPUT
    
    CALL T2NGRA_TERM
    
    RETURN
    
  END SUBROUTINE T2_NGRA
  
  !------------------------------------------------------------------
  !
  !            NODE GRAPH GENERATING ROUTINE 
  !                 FOR  MULTI-FLUID EQUATIONS FEM SOLVER
  !             
  !             COMPRESSED ROW STRAGE FORMAT
  !
  !      FOR LINEAR RECTANGULAR ELEMENT   4------3 
  !                                       |      |
  !                                       |      |
  !                                       1------2
  !
  !------------------------------------------------------------------  
  
  SUBROUTINE T2NGRA_NGRAPH1
    
    ! CONSTRUCT NON-DEGENERATED NODE 
    !             - DEGENERATED NODE GRAPH! CHECKED 
    CALL T2NGRA_CRT1    ! checked i2crt
        
    ! CONSTRUCT ELEMENT - NODE GRAPH      ! CHECKED 
    CALL T2NGRA_ENR1    ! checked i3enr
    
    ! CONSTRUCT NODE - ELEMENT  GRAPH     ! CHECKED
    CALL T2NGRA_NER1    ! i1eidr,i1eidc
    
    ! CONSTRUCT NODE - NODE GRAPH ! CHECKED 
    CALL T2NGRA_NNR1    ! NodeColCRS, NodeRowCRS, NodeDiaCRS
    
    CALL T2NGRA_FSA1    ! 
    RETURN
    
  END SUBROUTINE T2NGRA_NGRAPH1
  
  SUBROUTINE T2NGRA_FSA1
    
    USE T2COMM,ONLY:NMMAX,NFAMX,i2crt,NodeFSA
    INTEGER(ikind)::im,cnt,crtb,crtc,crtb_tmp
    
    cnt = 0
    crtb_tmp = 1
    DO im = 1, NMMAX
       crtb = i2crt(2,im)
       crtc = i2crt(3,im)
       IF((crtb.GT.crtb_tmp).AND.(crtb.NE.crtc))THEN
          cnt = cnt + 1
          NodeFSA(1,cnt) = crtb-1
          NodeFSA(2,cnt) = crtc-1
       ENDIF
       crtb_tmp = crtb
    ENDDO
    
    NFAMX = cnt
    
    RETURN

  END SUBROUTINE T2NGRA_FSA1
  !
  !
  !
  !                     2014-06-18 H. Seto
  !
  SUBROUTINE T2NGRA_COORDINATE
    
    USE T2CNST,ONLY:PI
    USE T2COMM,ONLY:&
         NMMAX,NLMAX,NRMAX,CoordinateSwitch,&
         i1rdn1,i1pdn1,i1rdn2,i1pdn2,i1mc1d,&
         d1mc1d,d1rec,&
         GlobalCrd,i1mfc1
    
    INTEGER(ikind)::&
         i_l
    INTEGER(ikind)::&
         i1,j1,&
         i0mcnt,i0bcnt,i0rcnt,i0rdn1,i0pdn1,i0rdn2,i0pdn2
    REAL(   rkind)::&
         d0rsiz,d0psiz
    REAL(rkind),DIMENSION(:),ALLOCATABLE::&
         d1mcr1,d1mcp1
    
    SELECT CASE (CoordinateSwitch)
    
    CASE (1) ! polar system
       
       i0mcnt  = 0
       i0rcnt  = 1
       i0bcnt  = 1
    
       i1mc1d(i0rcnt) = 1
       d1mc1d(i0rcnt) = 0.D0
       
       ! s = r^{2} constant width grid 
       DO i_l = 1, NLMAX
          i0rdn1 = i1rdn1(i_l)
          i0pdn1 = i1pdn1(i_l)
          i0rdn2 = i1rdn2(i_l)
          i0pdn2 = i1pdn2(i_l)
          
          ALLOCATE(d1mcr1(i0rdn1),d1mcp1(i0pdn1))
          
          DO i1 = 1, i0rdn1
             d1mcr1(i1) = 0.D0
          ENDDO
          
          DO i1 = 1, i0pdn1
             d1mcp1(i1) = 0.D0
          ENDDO
       
          d0rsiz = (d1rec(i_l)-d1rec(i_l-1))/DBLE(i0rdn2)
          d0psiz = 2.d0*PI/DBLE(i0pdn2)
          
          DO i1 = 1, i0rdn1
             d1mcr1(i1) = (d0rsiz*DBLE(i1-1)+d1rec(i_l-1))**2
          ENDDO
          
          DO i1 = 1, i0pdn1
             d1mcp1(i1) = d0psiz*DBLE(i1-1)
          ENDDO
          
          ! CONSTRUCT GEOMET
          DO j1 = 1, i0rdn1 
             IF(j1.NE.1)THEN
                i0rcnt         = i0rcnt + 1
                i0bcnt         = i0bcnt + i0pdn2
                i1mc1d(i0rcnt) = i0bcnt
                d1mc1d(i0rcnt) = d1mcr1(j1)
             ENDIF
             
             DO i1 = 1, i0pdn1
                i0mcnt = i0mcnt + 1
                !C MAGNETIC FLUX COORDINATE (RHO, CHI)
                GlobalCrd(1,i0mcnt) = d1mcr1(j1)
                GlobalCrd(2,i0mcnt) = d1mcp1(i1)
                i1mfc1(     i0mcnt) = i0bcnt
             ENDDO
          ENDDO
          
          DEALLOCATE(d1mcr1,d1mcp1)
          
       ENDDO
       
    CASE (2)! cartesian ! y pediodic
       
       i0mcnt  = 0
       i0rcnt  = 0
       i0bcnt  = 0
       
       DO i_l = 1, NLMAX
          i0rdn1 = i1rdn1(i_l)
          i0pdn1 = i1pdn1(i_l)
          i0rdn2 = i1rdn2(i_l)
          i0pdn2 = i1pdn2(i_l)
          
          ALLOCATE(d1mcr1(i0rdn1),d1mcp1(i0pdn1))
          
          DO i1 = 1, i0rdn1
             d1mcr1(i1) = 0.D0
          ENDDO
          
          DO i1 = 1, i0pdn1
             d1mcp1(i1) = 0.D0
          ENDDO
          
          d0rsiz = (d1rec(i_l)-d1rec(i_l-1))/DBLE(i0rdn2)
          d0psiz = 1.D0/DBLE(i0pdn2)
          
          DO i1 = 1, i0rdn1
             d1mcr1(i1) = d0rsiz*DBLE(i1-1)+d1rec(i_l-1)
          ENDDO
          
          DO i1 = 1, i0pdn1
             d1mcp1(i1) = d0psiz*DBLE(i1-1)
          ENDDO
          
          ! construct geometry
          DO j1 = 1, i0rdn1 
             IF((i_l.EQ.1).OR.(j1.NE.1))THEN
                i0rcnt         = i0rcnt + 1
                i0bcnt         = i0bcnt + i0pdn2
                i1mc1d(i0rcnt) = i0bcnt
                d1mc1d(i0rcnt) = d1mcr1(j1)
             ENDIF
             
             DO i1 = 1, i0pdn1
                i0mcnt = i0mcnt + 1
                ! cartesian coordinate
                GlobalCrd(1,i0mcnt) = d1mcr1(j1)
                GlobalCrd(2,i0mcnt) = d1mcp1(i1)
                i1mfc1(     i0mcnt) = i0bcnt
             ENDDO
          ENDDO
          
          DEALLOCATE(d1mcr1,d1mcp1)
          
       ENDDO
       
    END SELECT
    
    !DO i1 = 1,NRMAX
    !   print*,i1,i1mc1d(i1),d1mc1d(i1)
    !ENDDO

    IF(i0mcnt.NE.NMMAX)THEN
       WRITE(6,*)'CHECK SUM ERROR'
       STOP
    ENDIF
       
    RETURN
    
  END SUBROUTINE T2NGRA_COORDINATE
  
  
  !-------------------------------------------------------------------
  !  SUBROUTINE FOR NODE-MAPPING TABLE I2CRT
  !   
  !  I2CRT[1,1:NMMAX1]: NODE-NUMBER FOR COEF. CALC.
  !  I2CRT[2,1:NMMAX1]: NODE-NUMBER FOR 2D GRID
  !  I2CRT[3,1:NMMAX1]: NODE-NUMBER FOR 1D GRID
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2NGRA_CRT1

    USE T2COMM,ONLY:&
         NLMAX,NPMIN,CoordinateSwitch,i1rdn2,i1pdn2,i1mlvl,&
         i1pdn1,i1rdn1,i2crt,HangedNodeTable

    INTEGER(ikind)::&
         i0hcnt,i0mcnt,i0stm2,&
         i0ppc1,i0ppl2,i0ppc2,i0rdc1,&
         i0il,i0stl2,i0stc2,i0mlva,i0mlvb,&
         i1subtot(0:NLMAX),i0offset
    INTEGER(ikind)::&
         i_l,i0lidj,&
         i2,j2,i0stack1d,i0stack2d
    
    !C SUBTOTALS UP TO Nth-DOMAIN: I1SUBTOT[N]    
    SELECT CASE (CoordinateSwitch)
    CASE (1)
       DO i_l = 0, NLMAX
          i1subtot(i_l) = 1
       ENDDO
    CASE (2)
       DO i_l = 0, NLMAX
          i1subtot(i_l) = NPMIN
       ENDDO
    END SELECT

    DO i_l=1,NLMAX
       DO i0lidj = 1, i_l
          i1subtot(i_l)=i1subtot(i_l)+i1rdn2(i0lidj)*i1pdn2(i0lidj)
       ENDDO
    ENDDO
    
    ! SET NODE NUMBER 
    i0mcnt = 0
    i0hcnt = 0
    
    DO i_l= 1, NLMAX
       
       i0mlva=i1mlvl(i_l-1)
       i0mlvb=i1mlvl(i_l  )
       
       i0ppc2 = i1pdn2(i_l)
       i0ppl2 = INT(i0ppc2/2)
       i0ppc1 = i1pdn1(i_l)
       i0rdc1 = i1rdn1(i_l)
       ! 
       ! SET OFFSET
       !
       i0offset = 0
       
       IF(i1mlvl(i_l).GE.3)THEN
          DO i2 = 0, i1mlvl(i_l) - 3
             i0offset = i0offset+i1pdn2(1)*(2**i2)
          ENDDO
       END IF

       ! SET SUBTOTALS
       i0stl2  = i1subtot(i_l-1) - i1pdn2(i_l-1)
       i0stc2  = i1subtot(i_l-1)
       i0stm2  = i1subtot(NLMAX)
       
       IF(i0mlva.EQ.0)THEN
          ! DOMAIN WITH x = 0 AXIS BOUNDARY
          DO i2=1,i0rdc1
             IF(i2.EQ.1)THEN
                !C FOR AXIS IN POLAR COODRINATE
                SELECT CASE(CoordinateSwitch)
                CASE (1)
                   DO j2 = 1, i0ppc1
                      i0mcnt = i0mcnt + 1
                      i2crt(1,i0mcnt) = i0mcnt
                      i2crt(2,i0mcnt) = 1
                      i2crt(3,i0mcnt) = 1 
                   ENDDO
                CASE (2)
                   DO j2=1,i0ppc1
                      i0mcnt = i0mcnt + 1
                      i2crt(1,i0mcnt) = i0mcnt
                      i2crt(2,i0mcnt) = MOD(j2-1,i0ppc2)+1 
                      i2crt(3,i0mcnt) = i0ppc2
                   ENDDO
                END SELECT
             ELSE
                i0il = i2 - 2
                i0stack1d = i0stc2 + i0ppc2*i0il+i0ppc2               
                DO j2=1,i0ppc1
                   i0mcnt = i0mcnt + 1
                   i0stack2d =  i0stc2 + i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                   i2crt(1,i0mcnt) = i0mcnt
                   i2crt(2,i0mcnt) = i0stack2d 
                   i2crt(3,i0mcnt) = i0stack1d 
                ENDDO
             ENDIF
          ENDDO
       ELSEIF(i0mlva.NE.i0mlvb)THEN
          !C
          !C DOMAIN WITH INTERFACE BOUNDARY
          !C
          DO i2 = 1, i0rdc1
             IF(    i2.NE.1)THEN
                !C
                !C NODES EXCEPT INTERFACE BOUNADRY 
                !C
                i0il = i2 - 2
                i0stack1d = i0stc2 + i0ppc2*i0il + i0ppc2

                DO j2 = 1, i0ppc1
                   
                   i0mcnt    = i0mcnt+1
                   i0stack2d = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                   
                   i2crt(1,i0mcnt) = i0mcnt
                   i2crt(2,i0mcnt) = i0stack2d
                   i2crt(3,i0mcnt) = i0stack1d
             
                ENDDO
             
             ELSEIF(i2.EQ.1)THEN
                !C
                !C NODES ON INTERFACE BOUNDARY
                !C
                DO j2 = 1, i0ppc1
                   i0mcnt = i0mcnt + 1
                   i0stack1d = i0stl2 + i0ppl2 
                   IF(MOD(j2,2).EQ.1)THEN
                      i0stack2d = i0stl2 + MOD(INT((j2+1)/2)-1,i0ppl2)+1
                      i2crt(1,i0mcnt) = i0mcnt
                      i2crt(2,i0mcnt) = i0stack2d
                      i2crt(3,i0mcnt) = i0stack1d
                   ELSEIF(MOD(j2,2).EQ.0)THEN
                      i0hcnt    = i0hcnt+1
                      i0stack2d = i0stm2 + i0offset + INT(j2/2)

                      i2crt(1,i0mcnt)= i0mcnt
                      i2crt(2,i0mcnt)= i0stack2d
                      i2crt(3,i0mcnt)= i0stack1d

                      HangedNodeTable(1,i0hcnt)&
                           = i0stl2 + MOD(INT((j2-1)/2),i0ppl2)+1
                      HangedNodeTable(2,i0hcnt)&
                           = i0stl2 + MOD(INT((j2+1)/2),i0ppl2)+1
                   ENDIF
                ENDDO
             END IF
          ENDDO
       ELSEIF(i0mlva.EQ.i0mlvb)THEN
          !C
          !C DOMAIN WITHOUT INTERFACE BOUNDARY
          !C
          DO i2 = 1, i0rdc1
             
             i0il = i2 - 2
             i0stack1d = i0stc2+i0ppc2*i0il+i0ppc2
             DO j2 = 1, i0ppc1
                i0mcnt = i0mcnt + 1
                i0stack2d = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                i2crt(1,i0mcnt) = i0mcnt
                i2crt(2,i0mcnt) = i0stack2d
                i2crt(3,i0mcnt) = i0stack1d
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2NGRA_CRT1

  SUBROUTINE T2NGRA_ENR1
    
    USE T2COMM,ONLY:&
         NLMAX,NBMAX,&
         i1pdn1,i1pdn2,i1rdn1,i1rdn2,i2crt,&
         HangedNodeTable,ElementNodeGraph

    INTEGER(ikind)::&
         i2,    j2,    i_l,i0ecnt,&
         i0ll , i0lr , i0ul , i0ur,&
         i0ll1, i0lr1, i0ul1, i0ur1,&
         i0ll2, i0lr2, i0ul2, i0ur2,&
         i0ll3, i0lr3, i0ul3, i0ur3,&
         i0ll4, i0lr4, i0ul4, i0ur4,&
         i0ppc1,i0ppc2,i0rdc2,i0stc1
         
    i0ecnt = 0
    
    DO i_l = 1, NLMAX
       i0ppc1 = i1pdn1(i_l)
       i0ppc2 = i1pdn2(i_l)
       i0rdc2 = i1rdn2(i_l)
       i0stc1 = 0
!
!       DO i2=0,i1-1
!                     modified by AF: i1rdn1(0) and i1pdn1(0) are not defined
!
       DO i2 = 1, i_l - 1
          i0stc1 = i0stc1+i1rdn1(i2)*i1pdn1(i2)
       ENDDO
       
       DO i2 = 1,i0rdc2
       DO j2 = 1,i0ppc2
          i0ecnt = i0ecnt + 1
          
          i0ll = i0stc1 + i0ppc1*(i2-1) + j2
          i0lr = i0stc1 + i0ppc1*i2     + j2
          i0ur = i0stc1 + i0ppc1*i2     + j2 + 1
          i0ul = i0stc1 + i0ppc1*(i2-1) + j2 + 1
          
          i0ll1 = i2crt(1,i0ll)
          i0lr1 = i2crt(1,i0lr)
          i0ur1 = i2crt(1,i0ur)
          i0ul1 = i2crt(1,i0ul)

          i0ll2 = i2crt(2,i0ll)
          i0lr2 = i2crt(2,i0lr)
          i0ur2 = i2crt(2,i0ur)
          i0ul2 = i2crt(2,i0ul)

          i0ll3 = i2crt(2,i0ll)
          i0lr3 = i2crt(2,i0lr)
          IF(i0ll3.GT.NBMAX) i0ll3 = HangedNodeTable(1,i0ll3-NBMAX)
          IF(i0lr3.GT.NBMAX) i0lr3 = HangedNodeTable(1,i0lr3-NBMAX)
          i0ur3 = i0lr3
          i0ul3 = i0ll3
          
          i0ll4 = i2crt(3,i0ll)
          i0lr4 = i2crt(3,i0lr)
          i0ur4 = i2crt(3,i0ur)
          i0ul4 = i2crt(3,i0ul)
          
          ElementNodeGraph(1,1,i0ecnt) = i0ll1
          ElementNodeGraph(2,1,i0ecnt) = i0lr1
          ElementNodeGraph(3,1,i0ecnt) = i0ur1
          ElementNodeGraph(4,1,i0ecnt) = i0ul1
          
          ElementNodeGraph(1,2,i0ecnt) = i0ll2
          ElementNodeGraph(2,2,i0ecnt) = i0lr2
          ElementNodeGraph(3,2,i0ecnt) = i0ur2
          ElementNodeGraph(4,2,i0ecnt) = i0ul2

          ElementNodeGraph(1,3,i0ecnt) = i0ll3
          ElementNodeGraph(2,3,i0ecnt) = i0lr3
          ElementNodeGraph(3,3,i0ecnt) = i0ur3
          ElementNodeGraph(4,3,i0ecnt) = i0ul3
          
          ElementNodeGraph(1,4,i0ecnt) = i0ll4
          ElementNodeGraph(2,4,i0ecnt) = i0lr4
          ElementNodeGraph(3,4,i0ecnt) = i0ur4
          ElementNodeGraph(4,4,i0ecnt) = i0ul4

          !i1mlel(i0ecnt) = i_l
          
       ENDDO
       ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE T2NGRA_ENR1
  
  SUBROUTINE T2NGRA_NER1

    USE T2COMM,ONLY:&
         CoordinateSwitch,NLMAX,NECMX,i1eidr,i1eidc,i1pdn2,i1rdn2,i1mlvl

    INTEGER(ikind)::&
         i_l,i0ecnt,i0ncnt,&
         i0ll,i0lr,i0ul,i0ur,&
         i0mlva,               i0rda,i0pda,i0sbta,&
         i0mlvb,i0rb,i0pb,i0xb,i0rdb,i0pdb,i0sbtb,&
         i0mlvc,i0pc,     i0xc,i0rdc,i0pdc
    !C SUBTOTALS UP TO Nth-DOMAIN: I1SUBTOT[N]    
    
    !C SET NODE NUMBER 
    i0ecnt = 0
    i0ncnt = 1
    i1eidr(i0ncnt) = 1
    i0sbta = 0
    i0sbtb = 0
    DO i_l = 1, NLMAX
       
       i0mlva = i1mlvl(i_l-1)
       i0mlvb = i1mlvl(i_l  )
       i0mlvc = i1mlvl(i_l+1)

       i0pda = i1pdn2(i_l-1)       
       i0pdb = i1pdn2(i_l  )
       i0rda = i1rdn2(i_l-1)
       i0rdb = i1rdn2(i_l  )
       
       i0sbta = i0sbta + i0rda*i0pda
       
       !C PROCESS FOR THE POINT OF ORIGIN

       IF(i0mlva.EQ.0)THEN
          SELECT CASE(CoordinateSwitch)

          CASE (1)
             i0ncnt = i0ncnt+1
             DO i0pb=1,i0pdb
                i1eidc(i0ecnt+1)= i0pb
                i0ecnt = i0ecnt + 1
             ENDDO
             i1eidr(i0ncnt)= i0ecnt+1
          CASE (2)
             DO i0pb = 0, i0pdb-1
                
                i0ncnt = i0ncnt + 1
                i0xb = i0pb+i0pdb
                i0lr =  MOD(i0xb-1,i0pdb)+1 + i0sbta
                i0ur =  MOD(i0xb  ,i0pdb)+1 + i0sbta
             
                i1eidc(i0ecnt+1)= i0lr
                i1eidc(i0ecnt+2)= i0ur
                
                i0ecnt = i0ecnt + 2
                i1eidr(i0ncnt)= i0ecnt+1
             ENDDO
          END SELECT
       ENDIF
       
       DO i0rb = 0, i0rdb-2
          DO i0pb = 0, i0pdb-1
             
             i0ncnt = i0ncnt + 1

             i0xb = i0pb+i0pdb
             
             i0ll = i0pdb* i0rb     + MOD(i0xb-1,i0pdb)+1 + i0sbta
             i0lr = i0pdb*(i0rb +1) + MOD(i0xb-1,i0pdb)+1 + i0sbta
             i0ur = i0pdb*(i0rb +1) + MOD(i0xb  ,i0pdb)+1 + i0sbta
             i0ul = i0pdb* i0rb     + MOD(i0xb  ,i0pdb)+1 + i0sbta
             
             i1eidc(i0ecnt+1)= i0ll
             i1eidc(i0ecnt+2)= i0lr
             i1eidc(i0ecnt+3)= i0ur
             i1eidc(i0ecnt+4)= i0ul
             
             i0ecnt = i0ecnt + 4
             i1eidr(i0ncnt)= i0ecnt+1
             
          ENDDO
       ENDDO

       IF(i0mlvc.EQ.0)THEN
          
          DO i0pb = 0, i0pdb-1
             
             i0ncnt = i0ncnt + 1
             
             i0xb = i0pb+i0pdb
             
             i0ll = i0pdb*(i0rdb-1) + MOD(i0xb-1,i0pdb)+1 + i0sbta
             i0ul = i0pdb*(i0rdb-1) + MOD(i0xb  ,i0pdb)+1 + i0sbta
             
             i1eidc(i0ecnt+1)= i0ll
             i1eidc(i0ecnt+2)= i0ul
             
             i0ecnt = i0ecnt + 2
             i1eidr(i0ncnt)= i0ecnt+1
          ENDDO
          
       ELSEIF(i0mlvc.EQ.i0mlvb)THEN
          
          i0pdc = i1pdn2(i_l+1)
          i0rdc = i1rdn2(i_l+1)
          
          DO i0pb = 0,i0pdb-1
             
             i0ncnt = i0ncnt + 1
             
             i0xb = i0pb+i0pdb
             
             i0ll = i0pdb*(i0rdb-1) + MOD(i0xb-1,i0pdb)+1 + i0sbta
             i0lr = i0pdb* i0rdb    + MOD(i0xb-1,i0pdb)+1 + i0sbta
             i0ur = i0pdb* i0rdb    + MOD(i0xb  ,i0pdb)+1 + i0sbta
             i0ul = i0pdb*(i0rdb-1) + MOD(i0xb  ,i0pdb)+1 + i0sbta
             
             i1eidc(i0ecnt+1)= i0ll
             i1eidc(i0ecnt+2)= i0lr
             i1eidc(i0ecnt+3)= i0ur
             i1eidc(i0ecnt+4)= i0ul
             
             i0ecnt = i0ecnt + 4
             i1eidr(i0ncnt)= i0ecnt+1
             
          ENDDO
          
       ELSEIF(i0mlvc.EQ.(i0mlvb+1))THEN
       
          i0pdc = i1pdn2(i_l+1)
          i0rdc = i1rdn2(i_l+1)
          
          DO i0pb = 0, i0pdb-1
             
             i0ncnt = i0ncnt + 1
             i0pc = 2*i0pb

             i0xb = i0pb + i0pdb
             i0xc = i0pc + i0pdc

             i0ll = i0pdb*(i0rdb-1) + MOD(i0xb-1,i0pdb)+1 + i0sbta
             i0lr = i0pdb* i0rdb    + MOD(i0xc-1,i0pdc)+1 + i0sbta
             i0ur = i0pdb* i0rdb    + MOD(i0xc  ,i0pdc)+1 + i0sbta
             i0ul = i0pdb*(i0rdb-1) + MOD(i0xb  ,i0pdb)+1 + i0sbta
             
             i1eidc(i0ecnt+1)= i0ll
             i1eidc(i0ecnt+2)= i0lr
             i1eidc(i0ecnt+3)= i0ur
             i1eidc(i0ecnt+4)= i0ul
             
             i0ecnt = i0ecnt + 4
             i1eidr(i0ncnt)= i0ecnt+1
             
          ENDDO
          
       ELSE
          WRITE(6,*)'ERROR IN T2NGRA_NER1: WRONG MESH DATA'
          STOP
       END IF
    ENDDO

    !C VIRTUAL POINTS

    i0sbtb = 0
    
    DO i_l = 1, NLMAX
       
       i0mlvb = i1mlvl(i_l  )
       i0mlvc = i1mlvl(i_l+1)
       
       i0rdb  = i1rdn2(i_l  )   
       i0pdb  = i1pdn2(i_l  )
       i0sbtb = i0sbtb + i0rdb*i0pdb
       
       IF(i0mlvc.EQ.(i0mlvb+1))THEN
          
          i0pdc = i1pdn2(i_l+1)            
          i0rdc = i1rdn2(i_l+1)
          
          DO i0pb = 1,i0pdb
             
             i0ncnt = i0ncnt + 1
           
             i0pc = 2*i0pb-1
             i0xc = i0pc + i0pdc
             
             i0lr =  MOD(i0xc-1,i0pdc)+1 + i0sbtb
             i0ur =  MOD(i0xc  ,i0pdc)+1 + i0sbtb
             
             i1eidc(i0ecnt+1)= i0lr
             i1eidc(i0ecnt+2)= i0ur

             i0ecnt = i0ecnt + 2
             i1eidr(i0ncnt)= i0ecnt+1
             
          ENDDO
       ENDIF
    ENDDO
    IF(i0ecnt.NE.NECMX)THEN
       print*,i0ecnt,NECMX
       WRITE(6,*)'ERROR IN NER'
       STOP
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2NGRA_NER1

  SUBROUTINE T2NGRA_NNR1
    
    USE T2COMM, ONLY:&
         NNMAX,NBMAX,NLMAX,NodeRowCRS,NodeColCRS,NodeDiaCRS,&
         i1eidr,i1eidc,ElementNodeGraph
    INTEGER(ikind)::&
         i0nidi,i0nidj,&
         i0nsiz,i0esiz,i0nidr,i1,i0ncnt,&
         i0na,i0nb,i0eg,i0ec,i0nidc,&
         i0ngx,i0ex,i0ecx,i0ng
    INTEGER(ikind),DIMENSION(:),ALLOCATABLE::i1nstk

    NodeRowCRS(1) = 1
    i0nidc = 0

    !DO i0nidr = 2,NBMAX+1
    DO i0nidr = 1,NBMAX
       !i0esiz = i1eidr(i0nidr)-i1eidr(i0nidr-1)
       i0esiz = i1eidr(i0nidr+1)-i1eidr(i0nidr)
       i0nsiz = MAX(i0esiz,1000)
       ALLOCATE(i1nstk(1:i0nsiz))
       i1nstk(1:i0nsiz) = 0
       i0ncnt = 0       
       !DO i0eg = i1eidr(i0nidr-1),i1eidr(i0nidr)-1
       DO i0eg = i1eidr(i0nidr),i1eidr(i0nidr+1)-1
          
          i0ec = i1eidc(i0eg)
          
          DO i0nidi = 1,NNMAX
             
             i0ng = ElementNodeGraph(i0nidi,2,i0ec)
             
             IF(    i0ng.LE.NBMAX)THEN
                i0ncnt = i0ncnt + 1 
                i1nstk(i0ncnt) = i0ng
             ELSEIF(i0ng.GT.NBMAX)THEN
                DO i0ex  = i1eidr(i0ng),i1eidr(i0ng+1)-1
                   i0ecx = i1eidc(i0ex)
                   DO i0nidj = 1, NNMAX
                      i0ngx = ElementNodeGraph(i0nidj,2,i0ecx)
                      IF(i0ngx.LE.NBMAX)THEN
                         i0ncnt = i0ncnt + 1 
                         i1nstk(i0ncnt) = i0ngx
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       
       CALL SORT_ARRAY(i0nsiz,i1nstk)
       
       DO i1 = 1, i0nsiz
          
          IF(i1.EQ.1)THEN
             i0na = 0 
             i0nb = i1nstk(1)
          ELSE
             i0na = i1nstk(i1-1)
             i0nb = i1nstk(i1  )
          ENDIF
          
          IF(i0na.NE.i0nb)THEN
             i0nidc = i0nidc + 1
             NodeColCRS(i0nidc) = i0nb
             IF(i0nidr.EQ.i0nb) NodeDiaCRS(i0nidr)=i0nidc
          ENDIF
       ENDDO
       NodeRowCRS(i0nidr+1) = i0nidc+1
       
       DEALLOCATE(i1nstk)
       
    ENDDO
    
    RETURN

  END SUBROUTINE T2NGRA_NNR1
  !C------------------------------------------------------------------
  !C
  !C        ROUTINE FOR TERMINATE UNNECESSARY WORKING ARRAYS
  !C
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2NGRA_TERM
    RETURN
  END SUBROUTINE T2NGRA_TERM

  SUBROUTINE T2NGRA_NGRAPH2
    RETURN
  END SUBROUTINE T2NGRA_NGRAPH2

  SUBROUTINE T2NGRA_NGRAPH3
    RETURN
  END SUBROUTINE T2NGRA_NGRAPH3
  !C------------------------------------------------------------------
  !C 
  !C         BUBLE SORT ALGORITHM
  !C 
  !C------------------------------------------------------------------
  SUBROUTINE SORT_ARRAY(i0asiz, i1array)
  
    INTEGER(ikind),INTENT(IN)::i0asiz
    INTEGER(ikind),DIMENSION(i0asiz),INTENT(INOUT)::i1array
    INTEGER(ikind)::&
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

  SUBROUTINE T2_NGRA_OUTPUT

    USE T2COMM,ONLY:&
         NAMAX,NBMAX,NXMAX,NHMAX,NEMAX,NMMAX,NLMAX,NNMAX,&
         NNRMX,NERMX,NECMX,&
         HangedNodeTable, ElementNodeGraph,&
         i2crt, NodeRowCRS,NodeColCRS,NodeDiaCRS,i1eidr,i1eidc,&
         i1mlvl,i1emax,GlobalCrd,i1mfc1
       
         INTEGER(ikind)::i1
    INTEGER(ikind)::&
         i0midi,i0hidi,i0eidi,i_l,i0aidi,i0nidi

    OPEN(10,FILE='T2NGRA_I2CRT.dat')
    DO i0midi = 1, NMMAX
       WRITE(10,'("NUM1=",I7,1X,I7,1X,I7,1X,I7)')&
            i0midi,i2crt(1,i0midi),i2crt(2,i0midi),&
            i2crt(3,i0midi)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2NGRA_HANGEDNODETABLE.dat')
    DO i0hidi = 1, NHMAX
       WRITE(10,'("HNUM=",I9,1X,"LHN=",I9,1X,"UHN=",I9)')&
            i0hidi,HangedNodeTable(1,i0hidi),HangedNodeTable(2,i0hidi)
    ENDDO    
    CLOSE(10)

    OPEN(10,FILE='T2NGRA_ELEMENTNODEGRAPH1_TEST.dat')
    DO i0eidi = 1, NEMAX
       WRITE(10,'("ELM_NUMBER=",I9,1X,"ELEMENTNODEGRAPH=",4I9)')&
            i0eidi,(ElementNodeGraph(i0nidi,1,i0eidi),i0nidi = 1, NNMAX)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2NGRA_ELEMENTNODEGRAPH2_TEST.dat')
    DO i0eidi = 1, NEMAX
       WRITE(10,'("ELM_NUMBER=",I9,1X,"ELEMENTNODEGRAPH=",4I9)')&
            i0eidi,(ElementNodeGraph(i0nidi,2,i0eidi),i0nidi = 1, NNMAX)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2NGRA_ELEMENTNODEGRAPH3_TEST.dat')
    DO i0eidi = 1, NEMAX
       WRITE(10,'("ELM_NUMBER=",I9,1X,"ELEMENTNODEGRAPH=",4I9)')&
            i0eidi,(ElementNodeGraph(i0nidi,3,i0eidi),i0nidi = 1, NNMAX)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2NGRA_ELEMENTNODEGRAPH4_TEST.dat')
    DO i0eidi = 1, NEMAX
       WRITE(10,'("ELM_NUMBER=",I9,1X,"ELEMENTNODEGRAPH=",4I9)')&
            i0eidi,(ElementNodeGraph(i0nidi,4,i0eidi),i0nidi = 1, NNMAX)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2NGRA_NodeRowCRS.dat')
    DO i1=1,NNRMX
       IF(i1.EQ.1)THEN
          WRITE(10,'("I1=",I9,1X,"INDNR=",I9)')i1,NodeRowCRS(i1)
       ELSE
          WRITE(10,'("I1=",I9,1X,"INDNR=",I9,1X,"J1=",I9,1X,"NOC=",I9)')&
               i1,NodeRowCRS(i1),i1-1,NodeRowCRS(i1)-NodeRowCRS(i1-1)
       ENDIF
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2NGRA_NodeColCRS.dat')
    DO i0aidi = 1, NAMAX
       WRITE(10,'("I1=",I9,1X,"COL=",I9)')i0aidi,NodeColCRS(i0aidi)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='T2NGRA_NodeDiaCRS.dat')
    DO i0aidi = 1, NBMAX
       WRITE(10,'("I1=",I9,1X,"DIA=",I9)')i0aidi,NodeDiaCRS(i0aidi)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='T2NGRA_I1MLVL.dat')
    DO i_l = 1, NLMAX
       WRITE(10,'("MESH_NUMBER=",I3,1X,"MESH_LEVEL=",I3)')&
            i_l, i1mlvl(i_l)
    ENDDO
    WRITE(10,'("TOTAL_MMAX=",I9,1X,"TOTAL_BMAX=",I9)')&
         NMMAX,NBMAX
    CLOSE(10)
    
    OPEN(10,FILE='T2NGRA_I1EMAX.dat')
    DO i_l = 0, NLMAX
       WRITE(10,'("MESH_NUMBER=",I3,1X,"EMAX=",I9)')&
            i_l,i1emax(i_l)
    ENDDO
    WRITE(10,'("TOTAL_EMAX=",I9)')NEMAX
    CLOSE(10)
    
    !C CALCURATE NUMBER OF HANGED-NODE
    OPEN(10,FILE='T2NGRA_MATRIX.dat')
    WRITE(10,'("NODEROWCRS_ARRAY_SIZE=",I9)')NNRMX
    WRITE(10,'("NODECOLCRS_ARRAY_SIZE=",I9)')NAMAX
    WRITE(10,'("D2XVEC_ARRAY_SIZE=",I9)')NXMAX
    WRITE(10,'("D2BVEC_ARRAY_SIZE=",I9)')NBMAX
    CLOSE(10)
    
    OPEN(10,FILE='T2NGRA_MFC1.dat')
    DO i0midi = 1, NMMAX
       WRITE(10,'("i1=",I5,1X,"RHO=",D15.8,1X,"CHI=",D15.8,1X,"SN=",I5)')&
            i0midi,GlobalCrd(1,i0midi),GlobalCrd(2,i0midi),i1mfc1(i0midi)
    ENDDO
    CLOSE(10)
    
    OPEN(30,FILE='T2NGRA_I1EIDR.dat')
    DO i1 = 1,NERMX
       write(30,*)i1,i1eidr(i1)
    ENDDO
    CLOSE(30)

    OPEN(30,FILE='T2NGRA_I1EIDC.dat')
    DO i1 = 1,NECMX
       write(30,*)i1,i1eidc(i1)
    ENDDO
    CLOSE(30)

    
    RETURN
  END SUBROUTINE T2_NGRA_OUTPUT
END MODULE T2NGRA
