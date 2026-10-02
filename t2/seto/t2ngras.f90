!C------------------------------------------------------------------
!C
!C         MODULE T2NGRA
!C       
!C         NODE GRAPH GENERATOR FOR TASK/T2
!C
!C
!C------------------------------------------------------------------

MODULE T2NGRA
  
  USE T2CNST,ONLY:&
       i0ikind,i0rkind
  
  IMPLICIT NONE
  
  PUBLIC T2_NGRA,T2_NGRA_OUTPUT
  
  PRIVATE
  
CONTAINS 
  
  !C------------------------------------------------------------------
  !C
  !C          T2NGRA_MAIN
  !C
  !C          MAIN SUBROUTINE OF T2NGRA
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2_NGRA
    
    USE T2COMM,ONLY:i0nmax,idfile
    
    CALL T2NGRA_MSC
    
    SELECT CASE (i0nmax)
       
    CASE( 4) !C LINIEAR   ELEMENT
       CALL T2NGRA_NGRAPH1
    !CASE( 8) !C QUADRATIC ELEMENT 
    CASE( 9)
       CALL T2NGRA_NGRAPH2
    !CASE(12) !C CUBIC     ELEMENT
    CASE(16) !C CUBIC     ELEMENT
       CALL T2NGRA_NGRAPH3
    END SELECT
    
    IF(idfile.ge.4) CALL T2_NGRA_OUTPUT
    CALL T2_NGRA_OUTPUT
    
    CALL T2NGRA_TERM
    
    RETURN
    
  END SUBROUTINE T2_NGRA
  
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
  
  SUBROUTINE T2NGRA_NGRAPH1
    
    !C--------------------------------------------
    
    !C
    !C CONSTRUCT NON-DEGENERATED NODE
    !C                    - DEGENERATED NODE GRAPH
    
    CALL T2NGRA_CRT1    
    
    !C
    !C CONSTRUCT ELEMENT - NODE GRAPH
    !C
    
    CALL T2NGRA_ENR1
    
    !C
    !C CONSTRUCT NODE - ELEMENT  GRAPH
    !C
    
    CALL T2NGRA_NER1
    
    !C
    !C CONSTRUCT NODE - NODE GRAPH
    !C
    
    CALL T2NGRA_NNR1
    
    RETURN
 
  END SUBROUTINE T2NGRA_NGRAPH1
  
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
  SUBROUTINE T2NGRA_NGRAPH2
    RETURN
  END SUBROUTINE T2NGRA_NGRAPH2
  
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
  SUBROUTINE T2NGRA_NGRAPH3
    RETURN
  END SUBROUTINE T2NGRA_NGRAPH3
  
  !C
  !C
  !C
  !C                     2014-01-30 H.SETO
  !C
  SUBROUTINE T2NGRA_MSC
    
    USE T2CNST,ONLY:d0pi
    USE T2COMM,ONLY:&
         i0mmax,i0lmax,&
         i1rdn1,i1pdn1,i1rdn2,i1pdn2,i1mc1d,&
         d1mc1d,d1rsiz,d1psiz,d1msiz,d1rec,&
         d2mfc1,i1mfc1
    
    INTEGER(i0ikind)::&
         i0lidi
    INTEGER(i0ikind)::&
         i1,j1,&
         i0mcnt,i0bcnt,i0rcnt,i0rdn1,i0pdn1,i0rdn2,i0pdn2
    REAL(   i0rkind)::&
         d0rsiz,d0psiz
    REAL(i0rkind),DIMENSION(:),ALLOCATABLE::&
         d1mcr1,d1mcp1
    
    i0mcnt  = 0
    i0rcnt  = 1
    i0bcnt  = 1
    
    i1mc1d(i0rcnt) = 1
    d1mc1d(i0rcnt) = 0.D0
    !C
    !C s = r^{2} constant width grid 
    !C 
    DO i0lidi = 1, i0lmax
       
       i0rdn1 = i1rdn1(i0lidi)
       i0pdn1 = i1pdn1(i0lidi)
       i0rdn2 = i1rdn2(i0lidi)
       i0pdn2 = i1pdn2(i0lidi)
       
       ALLOCATE(d1mcr1(i0rdn1),d1mcp1(i0pdn1))
       
       DO i1 = 1, i0rdn1
          d1mcr1(i1) = 0.D0
       ENDDO
       
       DO i1 = 1, i0pdn1
          d1mcp1(i1) = 0.D0
       ENDDO
       
       
       d0rsiz = (d1rec(i0lidi)**2-d1rec(i0lidi-1)**2)/DBLE(i0rdn2)
       d0psiz = 2.d0*d0pi/DBLE(i0pdn2)
       
       d1rsiz(i0lidi) = d0rsiz
       d1psiz(i0lidi) = d0psiz
       d1msiz(i0lidi) = SQRT(d0rsiz*d0psiz)
       
       
       DO i1 = 1, i0rdn1
          d1mcr1(i1) = d0rsiz*DBLE(i1-1)+d1rec(i0lidi-1)**2
          !d1mcr1(i1) = d0rsiz*DBLE(i1-1)+d1rec(i0lidi-1)
       ENDDO
       
       DO i1 = 1, i0pdn1
          d1mcp1(i1) = d0psiz*DBLE(i1-1)
       ENDDO
    
       !C CONSTRUCT GEOMET
       DO j1 = 1, i0rdn1 
          
          IF(j1.NE.1)THEN
             i0rcnt         = i0rcnt + 1
             i0bcnt         = i0bcnt + i0pdn2
             i1mc1d(i0rcnt) = i0bcnt
             d1mc1d(i0rcnt) = d1mcr1(j1)
          ENDIF
          
          DO i1 = 1, i0pdn1
             
             i0mcnt = i0mcnt + 1
             
             !C
             !C MAGNETIC FLUX COORDINATE (RHO, CHI)
             !C
             d2mfc1(1,i0mcnt) = d1mcr1(j1)
             d2mfc1(2,i0mcnt) = d1mcp1(i1)
             i1mfc1(  i0mcnt) = i0bcnt
          ENDDO
       ENDDO
       
       DEALLOCATE(d1mcr1,d1mcp1)
       
    ENDDO
    
    IF(i0mcnt.NE.i0mmax)THEN
       WRITE(6,*)'CHECK SUM ERROR'
       STOP
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2NGRA_MSC
  
  !C
  !C
  !C
  !C
  !C
  SUBROUTINE T2_NGRA_OUTPUT

    USE T2COMM,ONLY:&
         i0amax,i0bmax,i0xmax,i0hmax,i0emax,i0mmax,i0lmax,i0nmax,&
         i0nrmx,i0ermx,i0ecmx,&
         i2hbc, i3enr, i2crt, i1nidr,i1nidc,i1eidr,i1eidc,&
         i1mlvl,i1mmax,i1bmax,i1emax,d2mfc1,i1mfc1

    INTEGER(i0ikind)::i1
    INTEGER(i0ikind)::&
         i0midi,i0hidi,i0eidi,i0lidi,i0aidi,i0nidi

    OPEN(10,FILE='I2CRT_TEST.dat')
    DO i0midi = 1, i0mmax
       WRITE(10,'("NUM1=",I7,1X,I7,1X,I7,1X,I7)')&
            i0midi,i2crt(1,i0midi),i2crt(2,i0midi),&
            i2crt(3,i0midi)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='I2HBC_TEST.dat')
    DO i0hidi = 1, i0hmax
       WRITE(10,'("HNUM=",I9,1X,"LHN=",I9,1X,"UHN=",I9)')&
            i0hidi,i2hbc(1,i0hidi),i2hbc(2,i0hidi)
    ENDDO    
    CLOSE(10)

    OPEN(10,FILE='I3ENR1_TEST.dat')
    DO i0eidi = 1, i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i0eidi,(i3enr(i0nidi,1,i0eidi),i0nidi = 1, i0nmax)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR2_TEST.dat')
    DO i0eidi = 1, i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i0eidi,(i3enr(i0nidi,2,i0eidi),i0nidi = 1, i0nmax)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR3_TEST.dat')
    DO i0eidi = 1, i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i0eidi,(i3enr(i0nidi,3,i0eidi),i0nidi = 1, i0nmax)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR4_TEST.dat')
    DO i0eidi = 1, i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i0eidi,(i3enr(i0nidi,4,i0eidi),i0nidi = 1, i0nmax)
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
    DO i0aidi = 1, i0amax
       WRITE(10,'("I1=",I9,1X,"INDNC=",I9)')i0aidi,i1nidc(i0aidi)
    ENDDO
    CLOSE(10)
    
    OPEN(10,FILE='TEST_NMAX.dat')
    DO i0lidi = 1, i0lmax
       WRITE(10,'("MESH_NUMBER=",I3,1X,"MESH_LEVEL=",I3)')&
            i0lidi, i1mlvl(i0lidi)
       WRITE(10,'("MMAX=",I9,1X,"BMAX=",I9)')&
            i1mmax(i0lidi),i1bmax(i0lidi)
    ENDDO
    WRITE(10,'("TOTAL_MMAX=",I9,1X,"TOTAL_BMAX=",I9)')&
         i0mmax,i0bmax
    CLOSE(10)
    
    OPEN(10,FILE='TEST_EMAX.dat')
    DO i0lidi = 0, i0lmax
       WRITE(10,'("MESH_NUMBER=",I3,1X,"EMAX=",I9)')&
            i0lidi,i1emax(i0lidi)
    ENDDO
    WRITE(10,'("TOTAL_EMAX=",I9)')i0emax
    CLOSE(10)
    
    !C CALCURATE NUMBER OF HANGED-NODE
    OPEN(10,FILE='MATRIX_INFO1.dat')
    WRITE(10,'("I1NIDR_ARRAY_SIZE=",I9)')i0nrmx
    WRITE(10,'("I1NIDC_ARRAY_SIZE=",I9)')i0amax
    WRITE(10,'("D2XVEC_ARRAY_SIZE=",I9)')i0xmax
    WRITE(10,'("D2BVEC_ARRAY_SIZE=",I9)')i0bmax
    CLOSE(10)
    
    OPEN(10,FILE='MFC1_CHECK.dat')
    DO i0midi = 1, i0mmax
       WRITE(10,'("i1=",I5,1X,"RHO=",D15.8,1X,"CHI=",D15.8,1X,"SN=",I5)')&
            i0midi,d2mfc1(1,i0midi),d2mfc1(2,i0midi),i1mfc1(i0midi)
    ENDDO
    CLOSE(10)
    
    OPEN(30,FILE='I1EIDR.dat')
    DO i1 = 1,i0ermx
       write(30,*)i1,i1eidr(i1)
    ENDDO
    CLOSE(30)

    OPEN(30,FILE='I1EIDC.dat')
    DO i1 = 1,i0ecmx
       write(30,*)i1,i1eidc(i1)
    ENDDO
    CLOSE(30)

    RETURN
  END SUBROUTINE T2_NGRA_OUTPUT

  !C-------------------------------------------------------------------
  !C SUBROUTINE FOR NODE-MAPPING TABLE I2CRT
  !C   
  !C I2CRT[1,1:I0MMAX1]: NODE-NUMBER FOR COEF. CALC.
  !C I2CRT[2,1:I0MMAX1]: NODE-NUMBER FOR 2D GRID
  !C I2CRT[3,1:I0MMAX1]: NODE-NUMBER FOR 1D GRID
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2NGRA_CRT1

    USE T2COMM,ONLY:&
         i0lmax,i1rdn2,i1pdn2,i1mlvl,i1pdn1,i1rdn1,i2crt,i2hbc

    INTEGER(i0ikind)::&
         i0hcnt,i0mcnt,i0stm2,&
         i0ppc1,i0ppl2,i0ppc2,i0rdc1,&
         i0il,i0stl2,i0stc2,i0mlva,i0mlvb,&
         i1subtot(0:i0lmax),i0offset
    INTEGER(i0ikind)::&
         i0lidi,i0lidj,&
         i2,j2,i0stack1d,i0stack2d
    
    !C SUBTOTALS UP TO Nth-DOMAIN: I1SUBTOT[N]    
    DO i0lidi = 0, i0lmax
       i1subtot(i0lidi) = 1
    ENDDO
    
    DO i0lidi=1,i0lmax
       DO i0lidj = 1, i0lidi
          i1subtot(i0lidi)=i1subtot(i0lidi)+i1rdn2(i0lidj)*i1pdn2(i0lidj)
       ENDDO
    ENDDO
    !C
    !C SET NODE NUMBER 
    !C
    i0mcnt = 0
    i0hcnt = 0
    
    DO i0lidi= 1, i0lmax
       
       i0mlva=i1mlvl(i0lidi-1)
       i0mlvb=i1mlvl(i0lidi  )
       
       i0ppc2 = i1pdn2(i0lidi)
       i0ppl2 = INT(i0ppc2/2)
       i0ppc1 = i1pdn1(i0lidi)
       i0rdc1 = i1rdn1(i0lidi)
       !C 
       !C SET OFFSET
       !C
       i0offset = 0
       
       IF(i1mlvl(i0lidi).GE.3)THEN
          DO i2 = 0, i1mlvl(i0lidi) - 3
             i0offset = i0offset+i1pdn2(1)*(2**i2)
          ENDDO
       END IF
       !C
       !C SET SUBTOTALS
       !C
       i0stl2  = i1subtot(i0lidi-1) - i1pdn2(i0lidi-1)
       i0stc2  = i1subtot(i0lidi-1)
       i0stm2  = i1subtot(i0lmax)
       
       IF(i0mlva.EQ.0)THEN
          !C
          !C DOMAIN WITH AXIS BOUNDARY
          !C
          DO i2=1,i0rdc1
             IF(i2.NE.1)THEN
                !C
                !C NODES EXCEPT AXIS BOUNADRY 
                !C
                i0il = i2 - 2
                i0stack1d = i0stc2 + i0ppc2*i0il+i0ppc2               
                DO j2=1,i0ppc1
                   i0mcnt = i0mcnt + 1
                   i0stack2d =  i0stc2 + i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                   i2crt(1,i0mcnt) = i0mcnt
                   i2crt(2,i0mcnt) = i0stack2d 
                   i2crt(3,i0mcnt) = i0stack1d 
                ENDDO
             ELSEIF(    i2.EQ.1)THEN
                !C
                !C NODES ON AXIS BOUNADRY 
                !C
                DO j2 = 1, i0ppc1
                   
                   i0mcnt = i0mcnt + 1

                   i2crt(1,i0mcnt) = i0mcnt
                   i2crt(2,i0mcnt) = 1
                   i2crt(3,i0mcnt) = 1 
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

                      i2hbc(1,i0hcnt)&
                           = i0stl2 + MOD(INT((j2-1)/2),i0ppl2)+1
                      i2hbc(2,i0hcnt)&
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
         i0lmax,i0bmax,&
         i1pdn1,i1pdn2,i1rdn1,i1rdn2,i2crt,i2hbc,i3enr,i1mlel

    INTEGER(i0ikind)::&
         i2,    j2,    i0lidi,i0ecnt,&
         i0ll , i0lr , i0ul , i0ur,&
         i0ll1, i0lr1, i0ul1, i0ur1,&
         i0ll2, i0lr2, i0ul2, i0ur2,&
         i0ll3, i0lr3, i0ul3, i0ur3,&
         i0ll4, i0lr4, i0ul4, i0ur4,&
         i0ppc1,i0ppc2,i0rdc2,i0stc1
         
    i0ecnt = 0
    
    DO i0lidi = 1, i0lmax
       i0ppc1 = i1pdn1(i0lidi)
       i0ppc2 = i1pdn2(i0lidi)
       i0rdc2 = i1rdn2(i0lidi)
       i0stc1 = 0
!
!       DO i2=0,i1-1
!                     modified by AF: i1rdn1(0) and i1pdn1(0) are not defined
!
       DO i2 = 1, i0lidi - 1
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
          IF(i0ll3.GT.i0bmax) i0ll3 = i2hbc(1,i0ll3-i0bmax)
          IF(i0lr3.GT.i0bmax) i0lr3 = i2hbc(1,i0lr3-i0bmax)
          i0ur3 = i0lr3
          i0ul3 = i0ll3
          
          i0ll4 = i2crt(3,i0ll)
          i0lr4 = i2crt(3,i0lr)
          i0ur4 = i2crt(3,i0ur)
          i0ul4 = i2crt(3,i0ul)
          
          i3enr(1,1,i0ecnt) = i0ll1
          i3enr(2,1,i0ecnt) = i0lr1
          i3enr(3,1,i0ecnt) = i0ur1
          i3enr(4,1,i0ecnt) = i0ul1
          
          i3enr(1,2,i0ecnt) = i0ll2
          i3enr(2,2,i0ecnt) = i0lr2
          i3enr(3,2,i0ecnt) = i0ur2
          i3enr(4,2,i0ecnt) = i0ul2

          i3enr(1,3,i0ecnt) = i0ll3
          i3enr(2,3,i0ecnt) = i0lr3
          i3enr(3,3,i0ecnt) = i0ur3
          i3enr(4,3,i0ecnt) = i0ul3
          
          i3enr(1,4,i0ecnt) = i0ll4
          i3enr(2,4,i0ecnt) = i0lr4
          i3enr(3,4,i0ecnt) = i0ur4
          i3enr(4,4,i0ecnt) = i0ul4

          i1mlel(i0ecnt) = i0lidi
          
       ENDDO
       ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE T2NGRA_ENR1
  
  SUBROUTINE T2NGRA_NER1

    USE T2COMM,ONLY:&
         i0lmax,i0ecmx,i1eidr,i1eidc,i1pdn2,i1rdn2,i1mlvl

    INTEGER(i0ikind)::&
         i0lidi,i0ecnt,i0ncnt,&
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
    DO i0lidi = 1, i0lmax
       
       i0mlva = i1mlvl(i0lidi-1)
       i0mlvb = i1mlvl(i0lidi  )
       i0mlvc = i1mlvl(i0lidi+1)

       i0pda = i1pdn2(i0lidi-1)       
       i0pdb = i1pdn2(i0lidi  )
       i0rda = i1rdn2(i0lidi-1)
       i0rdb = i1rdn2(i0lidi  )
       
       i0sbta = i0sbta + i0rda*i0pda
       
       !C PROCESS FOR THE POINT OF ORIGIN

       IF(i0mlva.EQ.0)THEN
          
          i0ncnt = i0ncnt+1
          
          DO i0pb=1,i0pdb
             i1eidc(i0ecnt+1)= i0pb
             i0ecnt = i0ecnt + 1
          ENDDO
          
          i1eidr(i0ncnt)= i0ecnt+1
          
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
          
          i0pdc = i1pdn2(i0lidi+1)
          i0rdc = i1rdn2(i0lidi+1)
          
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
       
          i0pdc = i1pdn2(i0lidi+1)
          i0rdc = i1rdn2(i0lidi+1)
          
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
    
    DO i0lidi = 1, i0lmax
       
       i0mlvb = i1mlvl(i0lidi  )
       i0mlvc = i1mlvl(i0lidi+1)
       
       i0rdb  = i1rdn2(i0lidi  )   
       i0pdb  = i1pdn2(i0lidi  )
       i0sbtb = i0sbtb + i0rdb*i0pdb
       
       IF(i0mlvc.EQ.(i0mlvb+1))THEN
          
          i0pdc = i1pdn2(i0lidi+1)            
          i0rdc = i1rdn2(i0lidi+1)
          
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
    IF(i0ecnt.NE.i0ecmx)THEN
       print*,i0ecnt,i0ecmx
       WRITE(6,*)'ERROR IN NER'
       STOP
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2NGRA_NER1

  SUBROUTINE T2NGRA_NNR1
    
    USE T2COMM, ONLY:&
         i0nmax,i0bmax,i0lmax,i1nidr,i1nidc,i1eidr,i1eidc,i3enr
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,&
         i0nsiz,i0esiz,i0nidr,i1,i0ncnt,&
         i0na,i0nb,i0eg,i0ec,i0nidc,&
         i0ngx,i0ex,i0ecx,i0ng
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1nstk

    i1nidr(1) = 1
    i0nidc = 0

    DO i0nidr = 2,i0bmax+1
       
       i0esiz = i1eidr(i0nidr)-i1eidr(i0nidr-1)
       !i0nsiz = MAX(i0esiz,i0nmax*8)
       i0nsiz = MAX(i0esiz,1000)
       ALLOCATE(i1nstk(1:i0nsiz))
       i1nstk(1:i0nsiz) = 0
       i0ncnt = 0       
       
       DO i0eg = i1eidr(i0nidr-1),i1eidr(i0nidr)-1
          
          i0ec = i1eidc(i0eg)
          
          DO i0nidi = 1,i0nmax
             
             i0ng = i3enr(i0nidi,2,i0ec)
             
             IF(    i0ng.LE.i0bmax)THEN
                i0ncnt = i0ncnt + 1 
                i1nstk(i0ncnt) = i0ng
             ELSEIF(i0ng.GT.i0bmax)THEN
                DO i0ex  = i1eidr(i0ng),i1eidr(i0ng+1)-1
                   i0ecx = i1eidc(i0ex)
                   DO i0nidj = 1, i0nmax
                      i0ngx = i3enr(i0nidj,2,i0ecx)
                      IF(i0ngx.LE.i0bmax)THEN
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
             i1nidc(i0nidc) = i0nb
          ENDIF
       ENDDO
       
       i1nidr(i0nidr) = i0nidc+1
       
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
END MODULE T2NGRA
