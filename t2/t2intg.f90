!C-------------------------------------------------------------------- 
!C INTERPORATION AND INTEGRATION ROUTINE
!C ALGORITHM: GAUSSIAN QUADRATURE
!C
!C
!C LAST UPDATE 2014-02-08
!C
!C -------------------------------------------------------------------
MODULE T2INTG
  
  USE T2CNST, ONLY: i0rkind,i0ikind

  IMPLICIT NONE
  
  PUBLIC T2_INTG
  
  PRIVATE  
  
CONTAINS

  !C-------------------------------------------------------------------
  !C
  !C  INTEGRAL ARRAYS
  !C
  !C                     2014-01-29 H.SETO     
  !C 
  !C-------------------------------------------------------------------
  SUBROUTINE T2_INTG
    
    USE T2COMM, ONLY:&
         i0nmax,i0qmax,i0dmax,&
         d3imsn,d4iavn,d6iatn,d5idtn,d4igvn,d6igtn,&
         d3iesn,d5ievn,d7ietn,d2issn,&
         !d5imss,d6iavs,d8iats,d7idts,d6igvs,d8igts,&
         !d5iess,d7ievs,d9iets,d4isss,&
         d2wfct,d4ifnc
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,i0nidn,&
         i0didi,i0didj,i0didk,&
         i0qidi,i0qidj
    REAL(   i0rkind)::&
         d0ifnci,d0ifncj,d0ifnck,d0ifncl,d0ifncm,d0ifncn,&
         d0wfct, d0temp
    !C------------------------------------------------------
    
    CALL T2INTG_IFUNC
    
    !C
    !C
    !C STANDARD INTEGRAL ARRAYS FOR PG-FEM
    !C
    !C

    !C
    !C MS: D3IMSN
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax    
    DO i0nidk = 1, i0nmax
       d0temp = 0.D0
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax
          d0ifnci = d4ifnc(i0qidi,i0qidj,0,i0nidi) 
          d0ifncj = d4ifnc(i0qidi,i0qidj,0,i0nidj) 
          d0ifnck = d4ifnc(i0qidi,i0qidj,0,i0nidk) 
          d0wfct  = d2wfct(i0qidi,i0qidj         )
          d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
       ENDDO
       ENDDO
       d3imsn(i0nidk,i0nidi,i0nidj) = d0temp
    ENDDO
    ENDDO
    ENDDO

    !C
    !C AV: D4IAVN 
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax    
    DO i0nidk = 1, i0nmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,0,     i0nidi) 
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didi,i0nidj) 
             d0ifnck = d4ifnc(i0qidi,i0qidj,0,     i0nidk) 
             d0wfct  = d2wfct(i0qidi,i0qidj              ) 
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
             d0ifnci = d4ifnc(i0qidi,i0qidj,0,     i0nidi) 
             d0ifncj = d4ifnc(i0qidi,i0qidj,0,     i0nidj) 
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didi,i0nidk) 
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
          ENDDO
          ENDDO
          d4iavn(i0didi,i0nidk,i0nidi,i0nidj) = d0temp
       ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C
    !C AT: D6IATN 
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0,     i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didj,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0,     i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d6iatn(i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C DT: D5IDTN 
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi) 
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didj,i0nidj) 
             d0ifnck = d4ifnc(i0qidi,i0qidj,0,     i0nidk) 
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
          ENDDO
          ENDDO
          d5idtn(i0didi,i0didj,i0nidk,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C GV: D4IGVN
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,0,     i0nidi) 
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didi,i0nidj) 
             d0ifnck = d4ifnc(i0qidi,i0qidj,0,     i0nidk) 
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
          ENDDO
          ENDDO
          d4igvn(i0didi,i0nidk,i0nidi,i0nidj) = d0temp
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C GT: D6IGTN
    !C 
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax    
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,0,     i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didj,i0nidj) 
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didi,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0,     i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d6igtn(i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C
    !C ES: D3IESN
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidk = 1, i0nmax
       d0temp = 0.D0
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax
          d0ifnci = d4ifnc(i0qidi,i0qidj,0,i0nidi)
          d0ifncj = d4ifnc(i0qidi,i0qidj,0,i0nidj)
          d0ifnck = d4ifnc(i0qidi,i0qidj,0,i0nidk)
          d0wfct  = d2wfct(i0qidi,i0qidj         )
          d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
       ENDDO
       ENDDO
       d3iesn(i0nidk,i0nidi,i0nidj) = d0temp
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C EV: D5IEVN
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,0,     i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0,     i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didi,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0,     i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d5ievn(i0didi,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C
    !C ET: D7IETN
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidm = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,0,     i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0,     i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didi,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,i0didj,i0nidl)
             d0ifncm = d4ifnc(i0qidi,i0qidj,0,     i0nidm)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0ifncm*d0wfct
          ENDDO
          ENDDO
          d7ietn(i0didi,i0didj,i0nidk,i0nidl,i0nidm,i0nidi,i0nidj)&
               = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C
    !C SS: D2ISSN
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d0temp = 0.D0
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax
          d0ifnci = d4ifnc(i0qidi,i0qidj,0,i0nidi)
          d0ifncj = d4ifnc(i0qidi,i0qidj,0,i0nidj)
          d0wfct  = d2wfct(i0qidi,i0qidj         )
          d0temp  = d0temp + d0ifnci*d0ifncj*d0wfct
       ENDDO
       ENDDO
       d2issn(i0nidi,i0nidj) = d0temp
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2_INTG
  
  !C------------------------------------------------------------------
  !C
  !C GAUSSIAN QUADRATURE ARRAYS
  !C
  !C                     2014-01-29 H.SETO
  !C
  !C------------------------------------------------------------------
  SUBROUTINE T2INTG_IFUNC
    
    USE T2CNST, ONLY: d1absc32,d1wfct32
    USE T2COMM, ONLY:&
         i0nmax,i0qmax,i0dmax,&
         d1wfct,d1absc,d2wfct,d4ifnc
    
    INTEGER(i0ikind)::&
         i0qidi,i0qidi_odd,i0qidi_eve,&
         i0qidj,i0qidj_odd,i0qidj_eve
    
    REAL(   i0rkind):: d0x,d0y,d0wfct
    
    !C
    !C SET NUMBER OF ABSCISSAS FOR GAUSS INTEGARATION
    !C
    
    SELECT CASE (i0qmax)
       
    CASE(32) !C 32*32 POINTS 2D GAUSS INTEGRATION
       
       DO i0qidj = 1, 16
          
          i0qidj_odd = 2*i0qidj - 1
          i0qidj_eve = 2*i0qidj
          
          DO i0qidi = 1, 16
             
             i0qidi_odd = 2*i0qidi - 1
             i0qidi_eve = 2*i0qidi
             
             d0wfct = d1wfct32(i0qidi)*d1wfct32(i0qidj)

             d2wfct(i0qidi_odd,i0qidj_odd) = d0wfct
             d2wfct(i0qidi_eve,i0qidj_odd) = d0wfct
             d2wfct(i0qidi_odd,i0qidj_eve) = d0wfct
             d2wfct(i0qidi_eve,i0qidj_eve) = d0wfct
             
          ENDDO
          
          d0x = d1absc32(i0qidj)
          
          d1absc(i0qidj_odd) = - d0x
          d1absc(i0qidj_eve) =   d0x
          
       ENDDO
       
    CASE DEFAULT
       WRITE(6,*)'-------------------------------------------------'
       WRITE(6,*)'SUBROUTINE SET_INTEGRATED_INTERPOLATION FUNCTIONS'
       WRITE(6,*)'ERROR: ILLEGAL I0QMAX0'
       WRITE(6,*)'-------------------------------------------------'
       STOP
    END SELECT
    
    
    !C
    !C SET INTERPOLATION FUNCTION 
    !C
    SELECT CASE (i0nmax)
       
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
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax

          d0x = d1absc(i0qidi)
          d0y = d1absc(i0qidj)
          
          d4ifnc(i0qidi,i0qidj,0,1) = (1.D0-d0x)*(1.D0-d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,2) = (1.D0+d0x)*(1.D0-d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,3) = (1.D0+d0x)*(1.D0+d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,4) = (1.D0-d0x)*(1.D0+d0y)/4.D0
          
       ENDDO
       ENDDO
       
       !C
       !C dPHI/dx (x_i,y_i)
       !
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax
          
          d0y = d1absc(i0qidj)
          
          d4ifnc(i0qidi,i0qidj,1,1) = -(1.D0-d0y)/4.D0 
          d4ifnc(i0qidi,i0qidj,1,2) =  (1.D0-d0y)/4.D0 
          d4ifnc(i0qidi,i0qidj,1,3) =  (1.D0+d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,1,4) = -(1.D0+d0y)/4.D0

       ENDDO
       ENDDO

       !C
       !C dPHI/dy (x_i,y_i)
       !C

       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax

          d0x = d1absc(i0qidi)

          d4ifnc(i0qidi,i0qidj,2,1) = -(1.D0-d0x)/4.D0 
          d4ifnc(i0qidi,i0qidj,2,2) = -(1.D0+d0x)/4.D0 
          d4ifnc(i0qidi,i0qidj,2,3) =  (1.D0+d0x)/4.D0
          d4ifnc(i0qidi,i0qidj,2,4) =  (1.D0-d0x)/4.D0

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
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax
       
          d0x = d1absc(i0qidi)
          d0y = d1absc(i0qidj)

          d4ifnc(i0qidi,i0qidj,0,1) &
               = -(1.D0-d0x)*(1.D0-d0y)*(1.D0+d0x+d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,2) &
               =  (1.D0-d0x**2)*(1.D0-d0y   )/2.D0
          d4ifnc(i0qidi,i0qidj,0,3) &
               = -(1.D0+d0x)*(1.D0-d0y)*(1.D0-d0x+d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,4) &
               =  (1.D0+d0x   )*(1.D0-d0y**2)/2.D0
          d4ifnc(i0qidi,i0qidj,0,5) &
               = -(1.D0+d0x)*(1.D0+d0y)*(1.D0-d0x-d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,6) &
               =  (1.D0-d0x**2)*(1.D0+d0y   )/2.D0
          d4ifnc(i0qidi,i0qidj,0,7) &
               = -(1.D0-d0x)*(1.D0+d0y)*(1.D0+d0x-d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,0,8) &
               =  (1.D0-d0x   )*(1.D0-d0y**2)/2.D0 
       ENDDO
       ENDDO
 
       !C
       !C dPHI/dx (x_i,y_i)
       !C
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax

          d0x = d1absc(i0qidi)
          d0y = d1absc(i0qidj)
          
          d4ifnc(i0qidi,i0qidj,1,1) =  (1.D0-d0y)*(2.D0*d0x+d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,1,2) = -d0x*(1.D0-d0y)
          d4ifnc(i0qidi,i0qidj,1,3) =  (1.D0-d0y)*(2.D0*d0x-d0y)/4.D0 
          d4ifnc(i0qidi,i0qidj,1,4) =  (1.D0-d0y**2)/2.D0 
          d4ifnc(i0qidi,i0qidj,1,5) =  (1.D0+d0y)*(2.D0*d0x+d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,1,6) = -d0x*(1.D0+d0y)
          d4ifnc(i0qidi,i0qidj,1,7) =  (1.D0+d0y)*(2.D0*d0x-d0y)/4.D0
          d4ifnc(i0qidi,i0qidj,1,8) = -(1.D0-d0y**2)/2.D0

       ENDDO
       ENDDO
       
       !C
       !C dPHI/dy (x_i,y_i)
       !C
       
       DO i0qidj = 1, i0qmax
       DO i0qidi = 1, i0qmax
          

          d0x = d1absc(i0qidi)
          d0y = d1absc(i0qidj)

          d4ifnc(i0qidi,i0qidj,2,1) =  (1.D0-d0x)*(2.D0*d0y+d0x)/4.D0
          d4ifnc(i0qidi,i0qidj,2,2) = -(1.D0-d0x**2)/2.D0
          d4ifnc(i0qidi,i0qidj,2,3) =  (1.D0+d0x)*(2.D0*d0y-d0x)/4.D0
          d4ifnc(i0qidi,i0qidj,2,4) = -d0y*(1.D0+d0x)
          d4ifnc(i0qidi,i0qidj,2,5) =  (1.D0+d0x)*(2.D0*d0y+d0x)/4.D0 
          d4ifnc(i0qidi,i0qidj,2,6) =  (1.D0-d0x**2)/2.D0
          d4ifnc(i0qidi,i0qidj,2,7) =  (1.D0-d0x)*(2.D0*d0y-d0x)/4.D0
          d4ifnc(i0qidi,i0qidj,2,8) = -d0y*(1.D0-d0x)

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
