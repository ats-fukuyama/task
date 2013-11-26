!C-------------------------------------------------------------------- 
!C INTERPORATION AND INTEGRATION ROUTINE
!C ALGORITHM: GAUSSIAN QUADRATURE
!C
!C
!C LAST UPDATE 2013-11-11
!C -------------------------------------------------------------------
MODULE T2INTG
  
  USE T2CNST, ONLY: i0rkind,i0ikind

  IMPLICIT NONE
  
  PUBLIC T2_INTG
  
  PRIVATE  
  
CONTAINS
  
  SUBROUTINE T2_INTG
    
    USE T2COMM, ONLY:&
         i0nmax0,i0amax0,i0dmax0,&
         d3imsn0,d4iavn0,d6iatn0,d5idtn0,d4igvn0,d6igtn0,&
         d3iesn0,d5ievn0,d7ietn0,d2issn0,&
         d1wfct0,d1absc0,d2wfct0,d4ifnc0
    
    INTEGER(i0ikind)::&
         i2,j2,&
         i3,j3,k3,l3,m3,&
         i4,j4
    !C------------------------------------------------------
    
    CALL T2INTG_IFUNC
    
    !C MS: D3IMSN0 
    
    d3imsn0(1:i0nmax0,1:i0nmax0,1:i0nmax0)=0.D0
    
    DO k3=1,i0nmax0
    DO j3=1,i0nmax0
    DO i3=1,i0nmax0
       DO j4= 1,i0amax0
       DO i4= 1,i0amax0
          d3imsn0(       i3,j3,k3)&
               = d3imsn0(i3,j3,k3        )&
               + d4ifnc0(i3,      0,i4,j4)&
               * d4ifnc0(   j3,   0,i4,j4)&
               * d4ifnc0(      k3,0,i4,j4)&
               * d2wfct0(           i4,j4)
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C AV: D4IAVN0 

    d4iavn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0) = 0.D0
    
    DO i2=1,i0dmax0
       DO k3=1,i0nmax0
       DO j3=1,i0nmax0
       DO i3=1,i0nmax0
          DO j4=1,i0amax0
          DO i4=1,i0amax0
             d4iavn0(       i3,j3,k3,i2        )&
                  = d4iavn0(i3,j3,k3,i2        )&
                  + d4ifnc0(i3,         0,i4,j4)&
                  * d4ifnc0(   j3,   i2,  i4,j4)&
                  * d4ifnc0(      k3,   0,i4,j4)&
                  * d2wfct0(              i4,j4)&
                  + d4ifnc0(i3,         0,i4,j4)&
                  * d4ifnc0(   j3,      0,i4,j4)&
                  * d4ifnc0(      k3,i2,  i4,j4)&
                  * d2wfct0(              i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO

    !C AT: D6IATN0 
    
    d6iatn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0,1:i0dmax0) = 0.D0
    
    DO j2 = 1, i0dmax0
    DO i2 = 1, i0dmax0
       DO l3 = 1, i0nmax0
       DO k3 = 1, i0nmax0
       DO j3 = 1, i0nmax0
       DO i3 = 1, i0nmax0
          DO j4 = 1, i0amax0
          DO i4 = 1, i0amax0
             d6iatn0(       i3,j3,k3,l3,i2,j2        )&
                  = d6iatn0(i3,j3,k3,l3,i2,j2        )&
                  + d4ifnc0(i3,         i2,     i4,j4)&
                  * d4ifnc0(   j3,            0,i4,j4)&
                  * d4ifnc0(      k3,      j2,  i4,j4)&
                  * d4ifnc0(         l3,      0,i4,j4)&
                  * d2wfct0(                    i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C DT: D5IDTN0 
    
    d5idtn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0,1:i0dmax0) = 0.D0
    
    DO i2=1,i0dmax0
    DO j2=1,i0dmax0
       DO k3=1,i0nmax0
       DO j3=1,i0nmax0
       DO i3=1,i0nmax0
          DO j4=1,i0amax0
          DO i4=1,i0amax0
             d5idtn0(       i3,j3,k3,i2,j2        )&
                  = d5idtn0(i3,j3,k3,i2,j2        )&
                  + d4ifnc0(i3,      i2,     i4,j4)&
                  * d4ifnc0(   j3,      j2,  i4,j4)&
                  * d4ifnc0(      k3,      0,i4,j4)&
                  * d2wfct0(                 i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    !C GV: D4IGVN0 
    
    d4igvn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0) = 0.D0

    DO i2=1,i0dmax0
       DO k3=1,i0nmax0
       DO j3=1,i0nmax0
       DO i3=1,i0nmax0
          DO j4=1,i0amax0
          DO i4=1,i0amax0
             d4igvn0(       i3,j3,k3,i2        )&
                  = d4igvn0(i3,j3,k3,i2        )&
                  + d4ifnc0(i3,         0,i4,j4)&
                  * d4ifnc0(   j3,   i2,  i4,j4)&                       
                  * d4ifnc0(      k3,   0,i4,j4)&
                  * d2wfct0(              i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    
    !C GT: D6IGTN0 
    
    d6igtn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0,1:i0dmax0) = 0.D0

    DO j2=1,i0dmax0
    DO i2=1,i0dmax0
       DO l3=1,i0nmax0
       DO k3=1,i0nmax0
       DO j3=1,i0nmax0
       DO i3=1,i0nmax0
          DO j4=1,i0amax0
          DO i4=1,i0amax0
             d6igtn0(       i3,j3,k3,l3,i2,j2        )&
                  = d6igtn0(i3,j3,k3,l3,i2,j2        )&
                  + d4ifnc0(i3,               0,i4,j4)&
                  * d4ifnc0(   j3,         j2,  i4,j4)&                       
                  * d4ifnc0(      k3,   i2,     i4,j4)&
                  * d4ifnc0(         l3,      0,i4,j4)&
                  * d2wfct0(                    i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C ES: D3IESN0 
    
    d3iesn0(1:i0nmax0,1:i0nmax0,1:i0nmax0) = 0.D0
    
    DO k3=1,i0nmax0
    DO j3=1,i0nmax0
    DO i3=1,i0nmax0
       DO j4=1,i0amax0
       DO i4=1,i0amax0
          d3iesn0(       i3,j3,k3         )&
               = d3iesn0(i3,j3,k3         )&
               + d4ifnc0(i3,       0,i4,j4)&
               * d4ifnc0(    j3,   0,i4,j4)&
               * d4ifnc0(       k3,0,i4,j4)&
               * d2wfct0(            i4,j4)
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C EV: D5IESN0 
    
    d5ievn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0) = 0.D0
    
    DO i2=1,i0dmax0
       DO l3=1,i0nmax0
       DO k3=1,i0nmax0
       DO j3=1,i0nmax0
       DO i3=1,i0nmax0
          DO j4=1,i0amax0
          DO i4=1,i0amax0
             d5ievn0(       i3,j3,k3,l3,i2        )&
                  = d5ievn0(i3,j3,k3,l3,i2        )&
                  + d4ifnc0(i3,            0,i4,j4)&
                  * d4ifnc0(   j3,         0,i4,j4)&
                  * d4ifnc0(      k3,   i2,  i4,j4)&
                  * d4ifnc0(         l3,   0,i4,j4)&
                  * d2wfct0(                 i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO

    !C ET: D7IETN0 
    
    d7ietn0(1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,1:i0nmax0,&
            1:i0dmax0,1:i0dmax0) = 0.D0
    
    DO j2 = 1, i0dmax0
    DO i2 = 1, i0dmax0
       DO m3 = 1, i0nmax0
       DO l3 = 1, i0nmax0
       DO k3 = 1, i0nmax0
       DO j3 = 1, i0nmax0
       DO i3 = 1, i0nmax0
          DO j4 =1, i0amax0
          DO i4 =1, i0amax0
             d7ietn0(       i3,j3,k3,l3,m3,i2,j2        )&
                  = d7ietn0(i3,j3,k3,l3,m3,i2,j2        )&
                  + d4ifnc0(i3,                  0,i4,j4)&
                  * d4ifnc0(   j3,               0,i4,j4)&
                  * d4ifnc0(      k3,      i2,     i4,j4)&
                  * d4ifnc0(         l3,      j2,  i4,j4)&
                  * d4ifnc0(            m3,      0,i4,j4)&
                  * d2wfct0(                       i4,j4)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C SS: D2ISSN0
    
    d2issn0(1:i0nmax0,1:i0nmax0) = 0.D0
    
    DO j3=1,i0nmax0
    DO i3=1,i0nmax0
       DO j4=1,i0amax0
       DO i4=1,i0amax0
          d2issn0(       i3,j3        )&
               = d2issn0(i3,j3        )&
               + d4ifnc0(i3,   0,i4,j4)&
               * d4ifnc0(   j3,0,i4,j4)&
               * d2wfct0(        i4,j4)
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2_INTG
  
  !C------------------------------------------------------------------
  !C SUBROUTINE SET_GAUSSIAN_QUADRATURE_ARRAYS
  !C CHECKED 2012/06/19
  !C------------------------------------------------------------------
  SUBROUTINE T2INTG_IFUNC
    
    USE T2CNST, ONLY: d1absc32,d1wfct32
    USE T2COMM, ONLY:&
         i0nmax0,i0amax0,i0dmax0,&
         d1wfct0,d1absc0,d2wfct0,d4ifnc0

    INTEGER(i0ikind):: i1,j1,i0iodd,i0ieve,i0jodd,i0jeve
    REAL(   i0rkind):: d0x,d0y,d0wfct
    !C
    !C SET NUMBER OF ABSCISSAS FOR GAUSS INTEGARATION
    !C
    SELECT CASE (i0amax0)

    CASE(32)!C 32*32 POINTS 2D GAUSS INTEGRATION
       
       DO j1=1,16
          
          i0jodd=2*j1-1
          i0jeve=2*j1
          
          DO i1=1,16
             
             i0iodd=2*i1-1
             i0ieve=2*i1
             
             d0wfct=d1wfct32(i1)*d1wfct32(j1)
             d2wfct0(i0iodd,i0jodd)=d0wfct
             d2wfct0(i0ieve,i0jodd)=d0wfct
             d2wfct0(i0iodd,i0jeve)=d0wfct
             d2wfct0(i0ieve,i0jeve)=d0wfct
                          
          ENDDO
          
          d0x=d1absc32(j1)
          d1absc0(i0jodd)= - d0x
          d1absc0(i0jeve)=   d0x
          
       ENDDO
       
    CASE DEFAULT
       WRITE(6,*)'-------------------------------------------------'
       WRITE(6,*)'SUBROUTINE SET_INTEGRATED_INTERPOLATION FUNCTIONS'
       WRITE(6,*)'ERROR: ILLEGAL I0AMAX0'
       WRITE(6,*)'-------------------------------------------------'
       STOP
    END SELECT
    
    
    !C
    !C SET INTERPOLATION FUNCTION 
    !C
    SELECT CASE (i0nmax0)
       
    CASE(4)
       !C FOR 4-NODE LINEAR LAGANGIAN RECTANGULAR ELEMENT 
       !C  
       !C                 4------3
       !C                 |      |
       !C                 |      |
       !C                 1------2
       !C
       
       !C
       !C PHI (x_i,y_i)
       !C
       DO i1 = 1, i0amax0
       DO j1 = 1, i0amax0
          d0x=d1absc0(i1)
          d0y=d1absc0(j1)
          d4ifnc0(1,0,i1,j1) = (1.D0-d0x)*(1.D0-d0y)/4.D0
          d4ifnc0(2,0,i1,j1) = (1.D0+d0x)*(1.D0-d0y)/4.D0
          d4ifnc0(3,0,i1,j1) = (1.D0+d0x)*(1.D0+d0y)/4.D0
          d4ifnc0(4,0,i1,j1) = (1.D0-d0x)*(1.D0+d0y)/4.D0
       ENDDO
       ENDDO
       
       !C
       !C dPHI/dx (x_i,y_i)
       !C
       DO i1 = 1, i0amax0
       DO j1 = 1, i0amax0
          d0y=d1absc0(j1)
          d4ifnc0(1,1,i1,j1) = -(1.D0-d0y)/4.D0 
          d4ifnc0(2,1,i1,j1) =  (1.D0-d0y)/4.D0 
          d4ifnc0(3,1,i1,j1) =  (1.D0+d0y)/4.D0
          d4ifnc0(4,1,i1,j1) = -(1.D0+d0y)/4.D0
       ENDDO
       ENDDO

       !C
       !C dPHI/dy (x_i,y_i)
       !C
       DO i1 = 1, i0amax0
       DO j1 = 1, i0amax0
          d0x=d1absc0(i1)
          d4ifnc0(1,2,i1,j1) = -(1.D0-d0x)/4.D0 
          d4ifnc0(2,2,i1,j1) = -(1.D0+d0x)/4.D0 
          d4ifnc0(3,2,i1,j1) =  (1.D0+d0x)/4.D0
          d4ifnc0(4,2,i1,j1) =  (1.D0-d0x)/4.D0
       ENDDO
       ENDDO

    CASE(8)
       !C FOR 8-NODE QUADRADIC SERENDIPITY RECTANGULAR ELEMENT
       !C
       !C    7---6---5
       !C    |       |
       !C    8       4
       !C    |       |
       !C    1---2---3
       !C
       
       !C
       !C PHI(x_i,y_i) 
       !C
       DO i1 = 1, i0amax0
       DO j1 = 1, i0amax0
          d0x=d1absc0(i1)
          d0y=d1absc0(j1)
          d4ifnc0(1,0,i1,j1)&
               = -(1.D0-d0x)*(1.D0-d0y)*(1.D0+d0x+d0y)/4.D0
          d4ifnc0(2,0,i1,j1)&
               =  (1.D0-d0x**2)*(1.D0-d0y   )/2.D0
          d4ifnc0(3,0,i1,j1)&
               = -(1.D0+d0x)*(1.D0-d0y)*(1.D0-d0x+d0y)/4.D0
          d4ifnc0(4,0,i1,j1)&
               =  (1.D0+d0x   )*(1.D0-d0y**2)/2.D0
          d4ifnc0(5,0,i1,j1)&
               = -(1.D0+d0x)*(1.D0+d0y)*(1.D0-d0x-d0y)/4.D0
          d4ifnc0(6,0,i1,j1)&
               =  (1.D0-d0x**2)*(1.D0+d0y   )/2.D0
          d4ifnc0(7,0,i1,j1)&
               = -(1.D0-d0x)*(1.D0+d0y)*(1.D0+d0x-d0y)/4.D0
          d4ifnc0(8,0,i1,j1)&
               =  (1.D0-d0x   )*(1.D0-d0y**2)/2.D0 
       ENDDO
       ENDDO
 
       !C
       !C dPHI/dx (x_i,y_i)
       !C
       DO i1 = 1, i0amax0
       DO j1 = 1, i0amax0
          d0x=d1absc0(i1)
          d0y=d1absc0(j1)
          d4ifnc0(1,1,i1,j1) =  (1.D0-d0y)*(2.D0*d0x+d0y)/4.D0
          d4ifnc0(2,1,i1,j1) = -d0x*(1.D0-d0y)
          d4ifnc0(3,1,i1,j1) =  (1.D0-d0y)*(2.D0*d0x-d0y)/4.D0 
          d4ifnc0(4,1,i1,j1) =  (1.D0-d0y**2)/2.D0 
          d4ifnc0(5,1,i1,j1) =  (1.D0+d0y)*(2.D0*d0x+d0y)/4.D0
          d4ifnc0(6,1,i1,j1) = -d0x*(1.D0+d0y)
          d4ifnc0(7,1,i1,j1) =  (1.D0+d0y)*(2.D0*d0x-d0y)/4.D0
          d4ifnc0(8,1,i1,j1) = -(1.D0-d0y**2)/2.D0  
       ENDDO
       ENDDO
       
       !C
       !C dPHI/dy (x_i,y_i)
       !C
       DO i1 = 1, i0amax0
       DO j1 = 1, i0amax0
          d0x=d1absc0(i1)
          d0y=d1absc0(j1)
          d4ifnc0(1,2,i1,j1) =  (1.D0-d0x)*(2.D0*d0y+d0x)/4.D0
          d4ifnc0(2,2,i1,j1) = -(1.D0-d0x**2)/2.D0
          d4ifnc0(3,2,i1,j1) =  (1.D0+d0x)*(2.D0*d0y-d0x)/4.D0
          d4ifnc0(4,2,i1,j1) = -d0y*(1.D0+d0x)
          d4ifnc0(5,2,i1,j1) =  (1.D0+d0x)*(2.D0*d0y+d0x)/4.D0 
          d4ifnc0(6,2,i1,j1) =  (1.D0-d0x**2)/2.D0
          d4ifnc0(7,2,i1,j1) =  (1.D0-d0x)*(2.D0*d0y-d0x)/4.D0
          d4ifnc0(8,2,i1,j1) = -d0y*(1.D0-d0x) 
       ENDDO
       ENDDO
    !CASE(12)
       !C FOR 12-NODE CUBIC SERENDIPITY RECTANGULAR ELEMENT
       !C
       !C    10--09--08--07
       !C     |           |
       !C    11          06
       !C     |           |
       !C    12          05
       !C     |           |
       !C    01--02--03--04
       !C
       
    CASE DEFAULT
       WRITE(6,*)'-------------------------------------------------'
       WRITE(6,*)'SUBROUTINE SET_INTEGRATED_INTERPOLATION FUNCTIONS'
       WRITE(6,*)'ERROR: ILLEGAL I0NMAX0'
       WRITE(6,*)'-------------------------------------------------'
       STOP
    END SELECT
    RETURN
  END SUBROUTINE T2INTG_IFUNC
END MODULE T2INTG
