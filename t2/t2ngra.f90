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
    
    USE T2COMM,ONLY:i0nmax0,idfile
    
    CALL T2NGRA_MSC
    
    SELECT CASE (i0nmax0)
       
    CASE( 4) !C LINIEAR   ELEMENT
       CALL T2NGRA_NGRAPH1
    CASE( 8) !C QUADRATIC ELEMENT 
       CALL T2NGRA_NGRAPH2
    CASE(12) !C CUBIC     ELEMENT
       CALL T2NGRA_NGRAPH3
    END SELECT

    IF(idfile.ge.4) CALL T2_NGRA_OUTPUT
    
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

  SUBROUTINE T2NGRA_MSC
    
    USE T2CNST,ONLY:d0pi
    USE T2COMM

    INTEGER(i0ikind)::&
         i1,i2,j2,&
         i0ncnt1,i0ncnt4,i0mcnt,i0rdn1,i0pdn1,i0rdn2,i0pdn2
    REAL(   i0rkind)::&
         d0rsiz,d0psiz
    REAL(i0rkind),DIMENSION(:),ALLOCATABLE::&
         d1mcr1,d1mcp1
    !C------------------------------------------------------     
    !C SET MAGNETIC COORDINATES
    
    i0ncnt1  = 0
    i0ncnt4  = 1
    i0mcnt   = 1

    i1mfc4(i0ncnt4)= 1
    d1mfc4(i0ncnt4)= 0.D0
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
       d1psiz(i1) = 2.d0*d0pi/DBLE(i0pdn2)
       d1msiz(i1) = SQRT(d1rsiz(i1)*d1psiz(i1))
       
       d0rsiz=d1rsiz(i1)
       d0psiz=d1psiz(i1)
       
       DO i2= 1,i0rdn1
          d1mcr1(i2) = d0rsiz*DBLE(i2-1)+d1rec(i1-1)
       ENDDO
       
       DO i2=1,i0pdn1
          d1mcp1(i2) = d0psiz*DBLE(i2-1)
       ENDDO
    
       !C CONSTRUCT GEOMET
       DO j2=1,i0rdn1 
          
          IF(j2.NE.1)THEN
             i0ncnt4         = i0ncnt4 + 1
             i0mcnt          = i0mcnt  + i0pdn2
             i1mfc4(i0ncnt4) = i0mcnt
             d1mfc4(i0ncnt4) = d1mcr1(j2)
          ENDIF
          
          DO i2=1,i0pdn1
             
             i0ncnt1=i0ncnt1+1
             
             !C
             !C MAGNETIC FLUX COORDINATE (RHO, CHI)
             !C
             
             d2mfc1(i0ncnt1,1) = d1mcr1(j2)
             d2mfc1(i0ncnt1,2) = d1mcp1(i2)
             i1mfc1(i0ncnt1  ) = i0mcnt
       
          ENDDO
       ENDDO
       
       DEALLOCATE(d1mcr1,d1mcp1)
       
    ENDDO
    
    IF(i0ncnt1.NE.i0nmax1)THEN
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

    USE T2COMM

    INTEGER(i0ikind)::i1,i2
    
    OPEN(10,FILE='I2CRT_TEST.dat')
    DO i1=1,i0nmax1
       WRITE(10,'("NUM1=",I9,1X,"CRT1=",I9,1X,"CRT2=",I9)')&
            i1,i2crt(i1,1),i2crt(i1,2)
    ENDDO

     
    CLOSE(10)
    OPEN(10,FILE='I2HBC_TEST.dat')
    DO i1=1,i0hmax
       WRITE(10,'("HNUM=",I9,1X,"LHN=",I9,1X,"UHN=",I9)')&
            i1,i2hbc(i1,1),i2hbc(i1,2)
    ENDDO

    CLOSE(10)
    OPEN(10,FILE='I3ENR1_TEST.dat')
    DO i1=1,i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i1,(i3enr(1,i2,i1),i2=1,i0nmax0)
    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR2_TEST.dat')
    DO i1=1,i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i1,(i3enr(2,i2,i1),i2=1,i0nmax0)

    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR3_TEST.dat')
    DO i1=1,i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i1,(i3enr(3,i2,i1),i2=1,i0nmax0)

    ENDDO
    CLOSE(10)

    OPEN(10,FILE='I3ENR4_TEST.dat')
    DO i1=1,i0emax
       WRITE(10,'("ELM_NUMBER=",I9,1X,"I3ENR=",4I9)')&
            i1,(i3enr(4,i2,i1),i2=1,i0nmax0)

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
       WRITE(10,'("i1=",I5,1X,"RHO=",D15.8,1X,"CHI=",D15.8,1X,"SN=",I5)')&
            i1,d2mfc1(i1,1),d2mfc1(i1,2),i1mfc1(i1)
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
  !C I2CRT[1:I0NMAX1,1]: NODE-NUMBER FOR COEF. CALC.
  !C I2CRT[1:I0NMAX1,2]: NODE-NUMBER FOR 2D GRID
  !C I2CRT[1:I0NMAX1,3]: NODE-NUMBER FOR 1D GRID
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2NGRA_CRT1

    USE T2COMM

    INTEGER(i0ikind)::&
         i0ecnt,i0hcnt,i0rcnt,i0ccnt,i0ncnt,&
         i0ll,i0lr,i0ul,i0ur,&
         i0ppc1,i0ppl2,i0ppc2,i0ppr2,&
         i0rdc1,i0rdc2,&
         i0jl,i0jc,i0jr,i0il,&
         i0stc1,i0nd,i0nc,i0nu,&
         i0stl2,i0stc2,i0str2,&
         i0mlva,i0mlvb,i0mlvc,&
         i1subtot(0:i0lmax),i0offset
    INTEGER(i0ikind)::&
         i1,j1,i2,j2,i3,i0stack1d,i0stack2d
    
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
                i0stack1d = i0stc2 + i0ppc2*i0il+i0ppc2

                DO j2=1,i0ppc1
                   
                   i0ncnt=i0ncnt+1
                   i0stack2d =  i0stc2 + i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                   
                   i2crt(i0ncnt,1) = i0ncnt
                   i2crt(i0ncnt,2) = i0stack2d 
                   i2crt(i0ncnt,3) = i0stack1d 

                ENDDO
             ELSEIF(    i2.EQ.1)THEN
                !C NODES ON AXIS BOUNADRY 

                DO j2=1,i0ppc1

                   i0ncnt=i0ncnt+1

                   i2crt(i0ncnt,1) = i0ncnt
                   i2crt(i0ncnt,2) = 1
                   i2crt(i0ncnt,3) = 1 
                ENDDO
             ENDIF
          ENDDO
       ELSEIF(i0mlva.NE.i0mlvb)THEN
          !C DOMAIN WITH INTERFACE BOUNDARY
          DO i2=1,i0rdc1
             IF(    i2.NE.1)THEN
                !C NODES EXCEPT INTERFACE BOUNADRY 

                i0il=i2-2
                i0stack1d = i0stc2 + i0ppc2*i0il + i0ppc2

                DO j2=1,i0ppc1
                   
                   i0ncnt    = i0ncnt+1
                   i0stack2d = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                   
                   i2crt(i0ncnt,1) = i0ncnt
                   i2crt(i0ncnt,2) = i0stack2d
                   i2crt(i0ncnt,3) = i0stack1d
             
                ENDDO
             
             ELSEIF(i2.EQ.1)THEN
                !C NODES ON INTERFACE BOUNDARY
                DO j2=1,i0ppc1
                   
                   i0ncnt=i0ncnt+1
                   i0stack1d = i0stl2 + i0ppl2 
                   IF(MOD(j2,2).EQ.1)THEN
                      
                      i0stack2d = i0stl2 + MOD(INT((j2+1)/2)-1,i0ppl2)+1

                      i2crt(i0ncnt,1) = i0ncnt
                      i2crt(i0ncnt,2) = i0stack2d
                      i2crt(i0ncnt,3) = i0stack1d
                      
                   ELSEIF(MOD(j2,2).EQ.0)THEN
                      
                      i0hcnt    = i0hcnt+1
                      i0stack2d = i0stm2 + i0offset + INT(j2/2)
                      
                      i2crt(i0ncnt,1)= i0ncnt
                      i2crt(i0ncnt,2)= i0stack2d
                      i2crt(i0ncnt,3)= i0stack1d
                      
                      i2hbc(i0hcnt,1)&
                           = i0stl2 + MOD(INT((j2-1)/2),i0ppl2)+1
                      i2hbc(i0hcnt,2)&
                           = i0stl2 + MOD(INT((j2+1)/2),i0ppl2)+1
                   ENDIF
                ENDDO
             END IF
          ENDDO
       ELSEIF(i0mlva.EQ.i0mlvb)THEN
          
          !C DOMAIN WITHOUT INTERFACE BOUNDARY
          DO i2=1,i0rdc1
             
             i0il=i2-2
             !i0stack1d = i0stc2+i0ppc2*i0il+1
             i0stack1d = i0stc2+i0ppc2*i0il+i0ppc2
             DO j2=1,i0ppc1
                
                i0ncnt=i0ncnt+1
                i0stack2d = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1
                
                !i2crt(i0ncnt,1)= i0ncnt
                !i2crt(i0ncnt,2)&
                !     = i0stc2+i0ppc2*i0il+MOD(j2-1,i0ppc2)+1

                i2crt(i0ncnt,1) = i0ncnt
                i2crt(i0ncnt,2) = i0stack2d
                i2crt(i0ncnt,3) = i0stack1d
                
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    OPEN(30,FILE='CRT.dat')
    DO i1 =1,i0nmax1
       WRITE(30,*)'i=',i1,'CRT1=',i2crt(i1,1),&
            'CRT2=',i2crt(i1,2),'CRT3=',i2crt(i1,3)
    ENDDO
    CLOSE(30)
    RETURN
    
  END SUBROUTINE T2NGRA_CRT1

  SUBROUTINE T2NGRA_ENR1
    
    USE T2COMM

    INTEGER(i0ikind)::&
         i0ecnt,i0hcnt,i0rcnt,i0ccnt,i0ncnt,&
         i0ll ,i0lr ,i0ul ,i0ur,&
         i0ll1,i0lr1,i0ul1,i0ur1,&
         i0ll2,i0lr2,i0ul2,i0ur2,&
         i0ll3,i0lr3,i0ul3,i0ur3,&
         i0ll4,i0lr4,i0ul4,i0ur4,&
         i0ppc1,i0ppl2,i0ppc2,i0ppr2,&
         i0rdc1,i0rdc2,&
         i0jl,i0jc,i0jr,i0il,&
         i0stc1,i0nd,i0nc,i0nu,&
         i0stl2,i0stc2,i0str2,&
         i0mlva,i0mlvb,i0mlvc,&
         i1subtot(0:i0lmax),i0offset
    INTEGER(i0ikind)::&
         i1,j1,i2,j2,i3
    
    i0ecnt = 0
    
    DO i1=1,i0lmax
       i0ppc1 = i1pdn1(i1)
       i0ppc2 = i1pdn2(i1)
       i0rdc2 = i1rdn2(i1)
       i0stc1=0
!
!       DO i2=0,i1-1
!                     modified by AF: i1rdn1(0) and i1pdn1(0) are not defined
!
       DO i2=1,i1-1
          i0stc1 = i0stc1+i1rdn1(i2)*i1pdn1(i2)
       ENDDO
       DO i2=1,i0rdc2
       DO j2=1,i0ppc2
          i0ecnt = i0ecnt+1
          
          i0ll = i0stc1 + i0ppc1*(i2-1) + j2
          i0lr = i0stc1 + i0ppc1*i2     + j2
          i0ur = i0stc1 + i0ppc1*i2     + j2 + 1
          i0ul = i0stc1 + i0ppc1*(i2-1) + j2 + 1
          
          i0ll1 = i2crt(i0ll,1)
          i0lr1 = i2crt(i0lr,1)
          i0ur1 = i2crt(i0ur,1)
          i0ul1 = i2crt(i0ul,1)

          i0ll2 = i2crt(i0ll,2)
          i0lr2 = i2crt(i0lr,2)
          i0ur2 = i2crt(i0ur,2)
          i0ul2 = i2crt(i0ul,2)

          i0ll3 = i2crt(i0ll,2)
          i0lr3 = i2crt(i0lr,2)
          IF(i0ll3.GT.i0nmax2) i0ll3 = i2hbc(i0ll3-i0nmax2,1)
          IF(i0lr3.GT.i0nmax2) i0lr3 = i2hbc(i0lr3-i0nmax2,1)
          i0ur3 = i0lr3
          i0ul3 = i0ll3
          
          !i0ll4 = i1mfc1(i0ll)
          !i0lr4 = i1mfc1(i0lr)
          !i0ur4 = i1mfc1(i0ur)
          !i0ul4 = i1mfc1(i0ul)

          i0ll4 = i2crt(i0ll,3)
          i0lr4 = i2crt(i0lr,3)
          i0ur4 = i2crt(i0ur,3)
          i0ul4 = i2crt(i0ul,3)
          
          i3enr(1,1,i0ecnt) = i0ll1
          i3enr(1,2,i0ecnt) = i0lr1
          i3enr(1,3,i0ecnt) = i0ur1
          i3enr(1,4,i0ecnt) = i0ul1
          
          i3enr(2,1,i0ecnt) = i0ll2
          i3enr(2,2,i0ecnt) = i0lr2
          i3enr(2,3,i0ecnt) = i0ur2
          i3enr(2,4,i0ecnt) = i0ul2

          i3enr(3,1,i0ecnt) = i0ll3
          i3enr(3,2,i0ecnt) = i0lr3
          i3enr(3,3,i0ecnt) = i0ur3
          i3enr(3,4,i0ecnt) = i0ul3
          
          i3enr(4,1,i0ecnt) = i0ll4
          i3enr(4,2,i0ecnt) = i0lr4
          i3enr(4,3,i0ecnt) = i0ur4
          i3enr(4,4,i0ecnt) = i0ul4

          i1mlel(i0ecnt) = i1
          
       ENDDO
       ENDDO
    ENDDO

    RETURN

  END SUBROUTINE T2NGRA_ENR1

  SUBROUTINE T2NGRA_NER1

    USE T2COMM

    INTEGER(i0ikind)::&
         i1,i0ecnt,i0ncnt,&
         i0ll,i0lr,i0ul,i0ur,&
         i0mlva,i0ra,i0pa,     i0rda,i0pda,i0sbta,&
         i0mlvb,i0rb,i0pb,i0xb,i0rdb,i0pdb,i0sbtb,&
         i0mlvc,i0pc,i0rc,i0xc,i0rdc,i0pdc
    !C SUBTOTALS UP TO Nth-DOMAIN: I1SUBTOT[N]    
    
    !C SET NODE NUMBER 
    i0ecnt = 0
    i0ncnt = 1
    i1eidr(i0ncnt) = 1
    i0sbta = 0
    i0sbtb = 0
    DO i1=1,i0lmax
       
       i0mlva = i1mlvl(i1-1)
       i0mlvb = i1mlvl(i1)
       i0mlvc = i1mlvl(i1+1)

       i0pda = i1pdn2(i1-1)       
       i0pdb = i1pdn2(i1  )

       
       i0rda = i1rdn2(i1-1)
       i0rdb = i1rdn2(i1  )
       
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
       
       DO i0rb=0,i0rdb-2
          DO i0pb =0,i0pdb-1
             
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
          
          DO i0pb =0,i0pdb-1
             
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
          
          i0pdc = i1pdn2(i1+1)
          i0rdc = i1rdn2(i1+1)
          
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
       
          i0pdc = i1pdn2(i1+1)
          i0rdc = i1rdn2(i1+1)
          
          DO i0pb =0,i0pdb-1
             
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
    
    DO i1 = 1,i0lmax
       
       i0mlvb = i1mlvl(i1)
       i0mlvc = i1mlvl(i1+1)
       
       i0rdb  = i1rdn2(i1  )   
       i0pdb  = i1pdn2(i1  )
       i0sbtb = i0sbtb + i0rdb*i0pdb
       
       IF(i0mlvc.EQ.(i0mlvb+1))THEN
          
          i0pdc = i1pdn2(i1+1)            
          i0rdc = i1rdn2(i1+1)
          
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
    
    USE T2COMM
    INTEGER(i0ikind)::&
         i0nsiz,i0esiz,i0nidr,i1,i0ecnt,i0ncnt,&
         i0na,i0nb,i0eg,i0ec,i0nidc,&
         i0ngx,i0ex,i0ecx,i2,i0ng
    INTEGER(i0ikind),DIMENSION(:),ALLOCATABLE::i1nstk

    i1nidr(1) = 1
    i0nidc = 0

    DO i0nidr = 2,i0nmax2+1
       
       i0esiz = i1eidr(i0nidr)-i1eidr(i0nidr-1)
       i0nsiz = MAX(i0esiz,i0nmax0*8)
       
       ALLOCATE(i1nstk(1:i0nsiz))
       i1nstk(1:i0nsiz) = 0
       i0ncnt = 0       
       
       DO i0eg = i1eidr(i0nidr-1),i1eidr(i0nidr)-1
          
          i0ec = i1eidc(i0eg)
          
          DO i1 = 1,i0nmax0
             
             i0ng = i3enr(2,i1,i0ec)
             
             IF(    i0ng.LE.i0nmax2)THEN
                i0ncnt = i0ncnt + 1 
                i1nstk(i0ncnt) = i0ng
             ELSEIF(i0ng.GT.i0nmax2)THEN
                DO i0ex  = i1eidr(i0ng),i1eidr(i0ng+1)-1
                   i0ecx = i1eidc(i0ex)
                   DO i2 = 1,i0nmax0
                      i0ngx = i3enr(2,i2,i0ecx)
                      IF(i0ngx.LE.i0nmax2)THEN
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
