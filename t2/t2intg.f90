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

  !C-------------------------------------------------------------------
  !C
  !C  INTEGRAL ARRAYS
  !C
  !C                     2014-01-29 H.SETO     
  !C 
  !C-------------------------------------------------------------------
  SUBROUTINE T2_INTG
    
    USE T2COMM, ONLY:&
         i0nmax,i0amax,i0dmax,&
         d3imsn,d4iavn,d6iatn,d5idtn,d4igvn,d6igtn,&
         d3iesn,d5ievn,d7ietn,d2issn,&
         d2wfct,d4ifnc
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,&
         i0didi,i0didj,&
         i0aidi,i0aidj
    !C------------------------------------------------------
    
    CALL T2INTG_IFUNC
    
    !C
    !C MS: D3IMSN
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax    
    DO i0nidk = 1, i0nmax
       d3imsn(i0nidk,i0nidi,i0nidj) = 0.D0
       DO i0aidj = 1,i0amax
       DO i0aidi = 1,i0amax
          d3imsn(                       i0nidk,i0nidi,i0nidj) &
               = d3imsn(                i0nidk,i0nidi,i0nidj) &
               + d4ifnc(i0aidi,i0aidj,0,       i0nidi       ) &
               * d4ifnc(i0aidi,i0aidj,0,              i0nidj) &
               * d4ifnc(i0aidi,i0aidj,0,i0nidk              ) &
               * d2wfct(i0aidi,i0aidj                       )
       ENDDO
       ENDDO
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
          d4iavn(i0didi,i0nidk,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d4iavn(                       i0didi,i0nidk,i0nidi,i0nidj) &
                  = d4iavn(                i0didi,i0nidk,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,0,              i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,  i0didi,              i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,0,       i0nidk              ) &
                  * d2wfct(i0aidi,i0aidj                              ) &
                  + d4ifnc(i0aidi,i0aidj,0,              i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,0,                     i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,  i0didi,i0nidk              ) &
                  * d2wfct(i0aidi,i0aidj                              )
          ENDDO
          ENDDO
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
          d6iatn(i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d6iatn(                       i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) &
                  = d6iatn(                i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,  i0didi,                     i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,0,                                   i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,         i0didj,i0nidk                     ) &
                  * d4ifnc(i0aidi,i0aidj,0,                     i0nidl              ) &
                  * d2wfct(i0aidi,i0aidj                                            ) 
          ENDDO
          ENDDO
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
          d5idtn(i0didi,i0didj,i0nidk,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d5idtn(                       i0didi,i0didj,i0nidk,i0nidi,i0nidj) &
                  = d5idtn(                i0didi,i0didj,i0nidk,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,  i0didi,              i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,         i0didj,              i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,0,              i0nidk              ) &
                  * d2wfct(i0aidi,i0aidj                                     )
          ENDDO
          ENDDO
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
          d4igvn(i0didi,i0nidk,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d4igvn(                       i0didi,i0nidk,i0nidi,i0nidj) &
                  = d4igvn(                i0didi,i0nidk,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,0,              i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,  i0didi,              i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,0,       i0nidk              ) &
                  * d2wfct(i0aidi,i0aidj                              )
          ENDDO
          ENDDO
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
          d6igtn(i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d6igtn(                       i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) &
                  = d6igtn(                i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,0,                            i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,         i0didj,                     i0nidj) & 
                  * d4ifnc(i0aidi,i0aidj,  i0didi,       i0nidk                     ) &
                  * d4ifnc(i0aidi,i0aidj,0,                     i0nidl              ) &
                  * d2wfct(i0aidi,i0aidj                                            )
          ENDDO
          ENDDO
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
       d3iesn(i0nidk,i0nidi,i0nidj) = 0.D0
       DO i0aidj =1, i0amax
       DO i0aidi =1, i0amax
          d3iesn(                       i0nidi,i0nidj,i0nidk) &
               = d3iesn(                i0nidi,i0nidj,i0nidk) &
               + d4ifnc(i0aidi,i0aidj,0,i0nidi              ) &
               * d4ifnc(i0aidi,i0aidj,0,       i0nidj       ) &
               * d4ifnc(i0aidi,i0aidj,0,              i0nidk) &
               * d2wfct(i0aidi,i0aidj                       )
       ENDDO
       ENDDO
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
          d5ievn(i0didi,i0nidk,i0nidl,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d5ievn(                       i0didi,i0nidk,i0nidl,i0nidi,i0nidj) &
                  = d5ievn(                i0didi,i0nidk,i0nidl,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,0,                     i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,0,                            i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,  i0didi,i0nidk                     ) &
                  * d4ifnc(i0aidi,i0aidj,0,              i0nidl              ) &
                  * d2wfct(i0aidi,i0aidj                                     )
          ENDDO
          ENDDO
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
          d7ietn(i0didi,i0didj,i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) = 0.D0
          DO i0aidj = 1, i0amax
          DO i0aidi = 1, i0amax
             d7ietn(                       i0didi,i0didj,i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) &
                  = d7ietn(                i0didi,i0didj,i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) &
                  + d4ifnc(i0aidi,i0aidj,0,                                   i0nidi       ) &
                  * d4ifnc(i0aidi,i0aidj,0,                                          i0nidj) &
                  * d4ifnc(i0aidi,i0aidj,  i0didi,       i0nidk                            ) &
                  * d4ifnc(i0aidi,i0aidj,         i0didj,       i0nidl                     ) &
                  * d4ifnc(i0aidi,i0aidj,0,                            i0nidm              ) &
                  * d2wfct(i0aidi,i0aidj                                                   )
          ENDDO
          ENDDO
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
       d2issn(i0nidi,i0nidj) = 0.D0
       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax
          d2issn(                       i0nidi,i0nidj) &
               = d2issn(                i0nidi,i0nidj) &
               + d4ifnc(i0aidi,i0aidj,0,i0nidi       ) &
               * d4ifnc(i0aidi,i0aidj,0,       i0nidj) &
               * d2wfct(i0aidi,i0aidj                )
       ENDDO
       ENDDO
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
         i0nmax,i0amax,i0dmax,&
         d1wfct,d1absc,d2wfct,d4ifnc
    
    INTEGER(i0ikind)::&
         i0aidi,i0aidi_odd,i0aidi_eve,&
         i0aidj,i0aidj_odd,i0aidj_eve
    
    REAL(   i0rkind):: d0x,d0y,d0wfct
    
    !C
    !C SET NUMBER OF ABSCISSAS FOR GAUSS INTEGARATION
    !C
    
    SELECT CASE (i0amax)
       
    CASE(32) !C 32*32 POINTS 2D GAUSS INTEGRATION
       
       DO i0aidj = 1, 16
          
          i0aidj_odd = 2*i0aidj - 1
          i0aidj_eve = 2*i0aidj
          
          DO i0aidi = 1, 16
             
             i0aidi_odd = 2*i0aidi - 1
             i0aidi_eve = 2*i0aidi
             
             d0wfct = d1wfct32(i0aidi)*d1wfct32(i0aidj)

             d2wfct(i0aidi_odd,i0aidj_odd) = d0wfct
             d2wfct(i0aidi_eve,i0aidj_odd) = d0wfct
             d2wfct(i0aidi_odd,i0aidj_eve) = d0wfct
             d2wfct(i0aidi_eve,i0aidj_eve) = d0wfct
             
          ENDDO
          
          d0x = d1absc32(i0aidj)
          
          d1absc(i0aidj_odd) = - d0x
          d1absc(i0aidj_eve) =   d0x
          
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
       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax

          d0x = d1absc(i0aidi)
          d0y = d1absc(i0aidj)
          
          d4ifnc(i0aidi,i0aidj,0,1) = (1.D0-d0x)*(1.D0-d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,2) = (1.D0+d0x)*(1.D0-d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,3) = (1.D0+d0x)*(1.D0+d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,4) = (1.D0-d0x)*(1.D0+d0y)/4.D0
          
       ENDDO
       ENDDO
       
       !C
       !C dPHI/dx (x_i,y_i)
       !
       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax
          
          d0y = d1absc(i0aidj)
          
          d4ifnc(i0aidi,i0aidj,1,1) = -(1.D0-d0y)/4.D0 
          d4ifnc(i0aidi,i0aidj,1,2) =  (1.D0-d0y)/4.D0 
          d4ifnc(i0aidi,i0aidj,1,3) =  (1.D0+d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,1,4) = -(1.D0+d0y)/4.D0

       ENDDO
       ENDDO

       !C
       !C dPHI/dy (x_i,y_i)
       !C

       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax

          d0x = d1absc(i0aidi)

          d4ifnc(i0aidi,i0aidj,2,1) = -(1.D0-d0x)/4.D0 
          d4ifnc(i0aidi,i0aidj,2,2) = -(1.D0+d0x)/4.D0 
          d4ifnc(i0aidi,i0aidj,2,3) =  (1.D0+d0x)/4.D0
          d4ifnc(i0aidi,i0aidj,2,4) =  (1.D0-d0x)/4.D0

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
       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax
       
          d0x = d1absc(i0aidi)
          d0y = d1absc(i0aidj)

          d4ifnc(i0aidi,i0aidj,0,1) &
               = -(1.D0-d0x)*(1.D0-d0y)*(1.D0+d0x+d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,2) &
               =  (1.D0-d0x**2)*(1.D0-d0y   )/2.D0
          d4ifnc(i0aidi,i0aidj,0,3) &
               = -(1.D0+d0x)*(1.D0-d0y)*(1.D0-d0x+d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,4) &
               =  (1.D0+d0x   )*(1.D0-d0y**2)/2.D0
          d4ifnc(i0aidi,i0aidj,0,5) &
               = -(1.D0+d0x)*(1.D0+d0y)*(1.D0-d0x-d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,6) &
               =  (1.D0-d0x**2)*(1.D0+d0y   )/2.D0
          d4ifnc(i0aidi,i0aidj,0,7) &
               = -(1.D0-d0x)*(1.D0+d0y)*(1.D0+d0x-d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,0,8) &
               =  (1.D0-d0x   )*(1.D0-d0y**2)/2.D0 
       ENDDO
       ENDDO
 
       !C
       !C dPHI/dx (x_i,y_i)
       !C
       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax

          d0x = d1absc(i0aidi)
          d0y = d1absc(i0aidj)
          
          d4ifnc(i0aidi,i0aidj,1,1) =  (1.D0-d0y)*(2.D0*d0x+d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,1,2) = -d0x*(1.D0-d0y)
          d4ifnc(i0aidi,i0aidj,1,3) =  (1.D0-d0y)*(2.D0*d0x-d0y)/4.D0 
          d4ifnc(i0aidi,i0aidj,1,4) =  (1.D0-d0y**2)/2.D0 
          d4ifnc(i0aidi,i0aidj,1,5) =  (1.D0+d0y)*(2.D0*d0x+d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,1,6) = -d0x*(1.D0+d0y)
          d4ifnc(i0aidi,i0aidj,1,7) =  (1.D0+d0y)*(2.D0*d0x-d0y)/4.D0
          d4ifnc(i0aidi,i0aidj,1,8) = -(1.D0-d0y**2)/2.D0

       ENDDO
       ENDDO
       
       !C
       !C dPHI/dy (x_i,y_i)
       !C
       
       DO i0aidj = 1, i0amax
       DO i0aidi = 1, i0amax
          

          d0x = d1absc(i0aidi)
          d0y = d1absc(i0aidj)

          d4ifnc(i0aidi,i0aidj,2,1) =  (1.D0-d0x)*(2.D0*d0y+d0x)/4.D0
          d4ifnc(i0aidi,i0aidj,2,2) = -(1.D0-d0x**2)/2.D0
          d4ifnc(i0aidi,i0aidj,2,3) =  (1.D0+d0x)*(2.D0*d0y-d0x)/4.D0
          d4ifnc(i0aidi,i0aidj,2,4) = -d0y*(1.D0+d0x)
          d4ifnc(i0aidi,i0aidj,2,5) =  (1.D0+d0x)*(2.D0*d0y+d0x)/4.D0 
          d4ifnc(i0aidi,i0aidj,2,6) =  (1.D0-d0x**2)/2.D0
          d4ifnc(i0aidi,i0aidj,2,7) =  (1.D0-d0x)*(2.D0*d0y-d0x)/4.D0
          d4ifnc(i0aidi,i0aidj,2,8) = -d0y*(1.D0-d0x)

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
