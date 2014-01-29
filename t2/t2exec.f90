!C
!C
!C
!C
!C
!C
!C
MODULE T2EXEC
  
  USE T2CNST, ONLY:&
       i0ikind,i0rkind

  IMPLICIT NONE
  
  INTEGER(i0ikind)::&
       i0vidi,i0vidj,i0widi,i0widj,i0eidi
  
  PUBLIC T2_EXEC
  
  PRIVATE  
  
CONTAINS
  
  !C
  !C MAIN ROUTINE OF T2EXEC
  !C FEM SOLVER FOR SIMULTANEOUS ADVECTION-DIFFUSION EQUATIONS
  !C 
  
  SUBROUTINE T2_EXEC
    
    USE T2COMM,ONLY:&
         i0emax,i0wmax,i0vmax,&
         i2msvt,i2avvt,i2atvt,i2dtvt,i2gvvt,i2gtvt,&
         i2esvt,i2evvt,i2etvt,i2ssvt,&
         i3atwt,i3gtwt,i3evwt,i4etwt
    
    INTEGER(i0ikind)::&
         i0msvt,i0avvt,i0atvt,i0dtvt,i0gvvt,i0gtvt,&
         i0esvt,i0evvt,i0etvt,i0ssvt,&
         i0atwt,i0gtwt,i0evwt,i0etwt
    
    REAL(4)::e0time_0,e0time_1
    
    CALL CPU_TIME(e0time_0)
    
    DO i0eidi = 1, i0emax
       
       CALL T2EXEC_LV 
       
       DO i0vidj = 1, i0vmax
       DO i0vidi = 1, i0vmax
          
          i0msvt = i2msvt(i0vidi,i0vidj)
          i0avvt = i2avvt(i0vidi,i0vidj)
          i0atvt = i2atvt(i0vidi,i0vidj)
          i0dtvt = i2dtvt(i0vidi,i0vidj)
          i0gvvt = i2gvvt(i0vidi,i0vidj)
          i0gtvt = i2gtvt(i0vidi,i0vidj)
          i0esvt = i2esvt(i0vidi,i0vidj)
          i0evvt = i2evvt(i0vidi,i0vidj)
          i0etvt = i2etvt(i0vidi,i0vidj)
          i0ssvt = i2ssvt(i0vidi,i0vidj)
          
          IF(i0msvt.EQ.1)       CALL T2EXEC_MS
          
          IF(i0avvt.EQ.1)       CALL T2EXEC_AV
          
          IF(i0atvt.EQ.1)THEN
             DO i0widi = 1, i0wmax
                i0atwt = i3atwt(i0widi,i0vidi,i0vidj)
                IF(i0atwt.EQ.1) CALL T2EXEC_AT
             ENDDO
          ENDIF
          
          IF(i0dtvt.EQ.1)       CALL T2EXEC_DT
          
          IF(i0gvvt.EQ.1)       CALL T2EXEC_GV
          
          IF(i0gtvt.EQ.1)THEN
             DO i0widi = 1, i0wmax
                i0gtwt = i3gtwt(i0widi,i0vidi,i0vidj)
                IF(i0gtwt.EQ.1) CALL T2EXEC_GT
             ENDDO
          ENDIF
          
          IF(i0esvt.EQ.1)       CALL T2EXEC_ES
          
          IF(i0evvt.EQ.1)THEN
             DO i0widi = 1, i0wmax
                i0evwt = i3evwt(i0widi,i0vidi,i0vidj)
                IF(i0evwt.EQ.1) CALL T2EXEC_EV
             ENDDO
          ENDIF
          
          IF(i0etvt.EQ.1)THEN          
             DO i0widj = 1, i0wmax
             DO i0widi = 1, i0wmax
                i0etwt = i4etwt(i0widi,i0widj,i0vidi,i0vidj)
                IF(i0etwt.EQ.1) CALL T2EXEC_ET
             ENDDO
             ENDDO
          ENDIF
          
          IF(i0ssvt.EQ.1)       CALL T2EXEC_SS
          
       ENDDO
       ENDDO
       
       CALL T2EXEC_STORE

    ENDDO
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC: completed:          cpu=', &
         e0time_1-e0time_0,' [s]'
    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_BC
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_BC completed:        cpu=', &
         e0time_1-e0time_0,' [s]'
    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_SOLVE
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_SOLVE completed:     cpu=', &
         e0time_1-e0time_0,' [s]'
    
    RETURN
    
  END SUBROUTINE T2_EXEC
  
  !C------------------------------------------------------------------ 
  !C  SUBROUTINE T2EXEC_SET_LOCAL_VARIABLES
  !C  　
  !C    * SET LOCAL NODE-ELEMENT GRAPH
  !C    * CALCULATE JACOBIAN OF PARAMETRIC SPACE
  !C    * STORE VARIABLES ST N-th TIMESTEP   
  !C    * SET WORKING ARRAY FOR DIFFERENTIAL
  !C
  !C                     2014-01-27 H.SETO
  !C
  !C------------------------------------------------------------------
  
  SUBROUTINE T2EXEC_LV
    
    USE T2COMM, ONLY:&
         i0nmax,i0dmax,i0vmax,i0smax,i0lid,&
         d1rsiz,d1psiz,d1guv,d2ws,d2wrks,d2kwns,&
         i2enr0,i3enr,i1mlel,&
         d2jmpm,d0jdmp,i0nmax2,i0nmax3
    
    INTEGER::&
         i0ridi,i0nidi,i0sidi,i0node,i0val
    
    !C
    !C SET LOCAL NODE-ELEMENT GRAPH
    !C
    !C   I2ENR(1:i0nmax,1) : FOR COEFFICIENT CALCULATION
    !C   I2ENR(1:i0nmax,2) : FOR 2D-2D 
    !C   I2ENR(1:i0nmax,3) : FOR 1D-2D,1D-1D 
    !C   I2ENR(1:i0nmax,4) : FOR 2D-1D
    !C
    
    DO i0ridi = 1, 4
       DO i0nidi = 1, i0nmax
          i2enr0(i0nidi,i0ridi)=i3enr(i0nidi,i0ridi,i0eidi)
       ENDDO
    ENDDO
    
    !C
    !C CALCULATE JACOBIAN OF PARAMETRIC SPACE
    !C 
    !C D0JDMP: JACOBIAN OF PARAMETRIC SPACE
    !C D2JMPM: INVERSE JACOBI MATRIX OF PARAMETRIC SPACE
    !C
    
    i0lid = i1mlel(i0eidi) 
    d0jdmp= d1rsiz(i0lid)*d1psiz(i0lid)/4.D0
    
    IF(d0jdmp.LE.0.D0)THEN
       WRITE(6,'("ERROR:: D1JDMP IS SINGULAR")')
       STOP
    ENDIF
    
    d2jmpm(1,1)= 2.D0/d1rsiz(i0lid)
    d2jmpm(1,2)= 0.D0
    d2jmpm(2,1)= 0.D0
    d2jmpm(2,2)= 2.D0/d1psiz(i0lid)
    
    !C
    !C STORE VARIABLES AT N-th TIMESTEP 
    !C xyz

    DO i0vidi = 1, i0vmax
       
       IF(    i0vidi.GT.3)THEN
          DO i0nidi = 1, i0nmax
             i0node = i2enr0(i0nidi,2)
             i0val  = i0vmax*(i0node - 1) + i0vidi
             d2kwns(i0nidi,i0vidi) = d1guv(i0val)
          ENDDO
       ELSEIF(i0vidi.LE.3)THEN
          DO i0nidi = 1, i0nmax
             i0node = i2enr0(i0nidi,4)
             i0val  = i0vmax*(i0node - 1) + i0vidi
             d2kwns(i0nidi,i0vidi) = d1guv(i0val)
          ENDDO
       ENDIF
    ENDDO
    
    !C SET WORKING ARRAY FOR DIFFERENTIAL

    !C 
    !C D2WRKS(1   ,:) : B    AT L-TH PICARD ITERATION
    !C D2WRKS(2   ,:) : R    AT L-TH PICARD ITERATION
    !C D2WRKS(2N+1,:) : Ub   AT L-TH PICARD ITERATION 
    !C D2WRKS(2N+2,:) : Qb/P AT L-TH PICARD ITERATION 
    !C
    
    DO i0nidi = 1, i0nmax
       
       i0node = i2enr0(i0nidi,1)
       
       d2wrks(1,i0nidi) = d2ws(1,i0node)
       d2wrks(2,i0nidi) = d2ws(2,i0node)
       
       DO i0sidi = 1, i0smax
          
          i0widi = 2*i0sidi
          
          d2wrks(i0widi+1,i0nidi) = d2ws(i0widi+1,i0node)
          d2wrks(i0widi+2,i0nidi) = d2ws(i0widi+2,i0node)
          
       ENDDO
    ENDDO
        
    RETURN
    
  END SUBROUTINE T2EXEC_LV
  
  SUBROUTINE T2EXEC_BC
    RETURN
  END SUBROUTINE T2EXEC_BC
  
  !C
  !C MODIFIED 2013-12-02
  !C
  !C
  SUBROUTINE T2EXEC_SOLVE
    
    USE T2COMM
    USE LIBMPI
    USE COMMPI
    USE LIBMTX
    
    RETURN
    
  END SUBROUTINE T2EXEC_SOLVE
  
  !C
  !C SUBROUTINES FOR CALCULATING  MATRCES 
  !C
  
  SUBROUTINE T2EXEC_MS
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d0jdmp,dt,&
         d3ms,d3imsn,d4smat,d2svec,d2kwns
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,&
         i0node
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d1svec(1:i0nmax),&
         d1mass(1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
    
    d2smat(1:i0nmax,1:i0nmax) = 0.D0
    
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       d1mass(i0nidi)&
            = d3ms(i0vidi,i0vidj,i0node)*d0jdmp/dt
    ENDDO
    
    !C
    !C MAIN LOOP
    !C
    
    DO i0nidj = 1, i0nmax       
    DO i0nidi = 1, i0nmax         
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidk = 1, i0nmax         
          d2smat(              i0nidi,i0nidj) &
               = d2smat(       i0nidi,i0nidj) &
               + d3imsn(i0nidk,i0nidi,i0nidj) &
               * d1mass(i0nidk              )
       ENDDO       
    ENDDO
    ENDDO
    
    DO i0nidi = 1, i0nmax
       d1svec(i0nidi) = 0.D0
       DO i0nidj = 1, i0nmax  
          d1svec(       i0nidi              ) &
               = d1svec(i0nidi              ) &
               + d2smat(i0nidi,i0nidj       ) &
               * d2kwns(       i0nidj,i0vidj)
       ENDDO
    ENDDO
    
    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO

    !C
    !C STORE SUBVECTOR
    !C

    DO i0nidi = 1, i0nmax
       d2svec(       i0nidi,i0vidi) &
            = d2svec(i0nidi,i0vidi) & 
            + d1svec(i0nidi       )
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_MS
  
  
  SUBROUTINE T2EXEC_AV
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d0jdmp,d2jmpm,&
         d4av,d4iavn,d4smat
    
    INTEGER(i0ikind)::&
         i0didi,i0didj,&
         i0nidi,i0nidj,i0nidk,&
         i0node
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d2velo(1:i0dmax,1:i0nmax),&
         d2temp(1:i0dmax,1:i0nmax)
    
    !C
    !C INTITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       DO i0didi = 1, i0dmax
          d2velo(     i0didi,i0nidi                     ) &
               = d4av(i0didi,       i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
    ENDDO
    
    !C
    !C MAIN LOOP
    !C
    
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
          d2temp(i0didj,i0nidk) = 0.D0
          DO i0didi = 1, i0dmax
             d2temp(              i0didj,i0nidk) &
                  = d2temp(       i0didj,i0nidk) &
                  + d2velo(i0didi,       i0nidk) &
                  * d2jmpm(i0didi,i0didj       ) 
          ENDDO
       ENDDO
    ENDDO
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidk = 1, i0nmax
          DO i0didj = 1, i0dmax
             d2smat(                     i0nidi,i0nidj) &
                  = d2smat(              i0nidi,i0nidj) &
                  + d4iavn(i0didj,i0nidk,i0nidi,i0nidj) &
                  * d2temp(i0didj,i0nidk              )
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AV
  
  
  
  SUBROUTINE T2EXEC_AT
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,&
         d6at,d6iatn,d4smat,d2wrks
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,&
         i0didi,i0didj,i0didk,i0didl,&
         i0node
    
    REAL(   i0rkind)::&
         d3velo(1:i0dmax,1:i0dmax,1:i0nmax),&
         d2smat(1:i0nmax,1:i0nmax),&
         d1atwi(1:i0nmax),&
         d3temp(1:i0nmax,1:i0dmax,1:i0dmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       d1atwi(i0nidi) = d2wrks(i0nidi,i0widi)
    ENDDO
    
    DO i0nidi = 1, i0nmax       
       i0node = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3velo(     i0didi,i0didj,i0nidi)&
               = d6at(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
       ENDDO
    ENDDO
         
    !C
    !C MAIN LOOP
    !C
    
    DO i0nidl = 1, i0nmax
       DO i0didl = 1, i0dmax
       DO i0didk = 1, i0dmax
          d3temp(i0didk,i0didl,i0nidl) = 0.D0          
          DO i0didi = 1, i0dmax
          DO i0didj = 1, i0dmax
             d3temp(                     i0didk,i0didl,i0nidl) &
                  = d3temp(              i0didk,i0didl,i0nidl) &
                  + d3velo(i0didi,i0didj,              i0nidl) &
                  * d2jmpm(i0didi,       i0didk              ) &
                  * d2jmpm(       i0didj,       i0didl       )
          END DO
          END DO
       END DO
       END DO
    END DO
       
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidl = 1, i0nmax
       DO i0nidk = 1, i0nmax
          DO i0didl = 1, i0dmax
          DO i0didk = 1, i0dmax
             d2smat(                                   i0nidi,i0nidj)&
                  = d2smat(                            i0nidi,i0nidj)&
                  + d6iatn(i0didk,i0didl,i0nidk,i0nidl,i0nidi,i0nidj)&
                  * d1atwi(              i0nidk                     )&
                  * d3temp(i0didk,i0didl,       i0nidl              )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AT

  
  
  SUBROUTINE T2EXEC_DT
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,&
         d5dt,d5idtn
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,&
         i0didi,i0didj,i0didk,i0didl,&
         i0node

    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d3diff(1:i0nmax,1:i0dmax,1:i0dmax), &
         d3temp(1:i0nmax,1:i0dmax,1:i0dmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3diff(     i0didi,i0didj,i0nidi)&
               = d5dt(i0didi,i0didj,i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
       ENDDO
    ENDDO

    !C
    !C MAIN LOOP
    !C

    DO i0nidk = 1, i0nmax
       DO i0didl = 1, i0dmax
       DO i0didk = 1, i0dmax
          d3temp(i0didk,i0didl,i0nidk) = 0.D0
          DO i0didj = 1, i0dmax
          DO i0didi = 1, i0dmax
             d3temp(                     i0didk,i0didl,i0nidk) &
                  = d3temp(              i0didk,i0didl,i0nidk) &
                  + d3diff(i0didi,i0didj,              i0nidk) &
                  * d2jmpm(i0didi,       i0didk              ) &
                  * d2jmpm(       i0didj,       i0didl       )
          END DO
          END DO     
       END DO
       END DO
    END DO

    DO i0nidi = 1, i0nmax
    DO i0nidj = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidk = 1, i0nmax
          DO i0didl = 1, i0dmax
          DO i0didk = 1, i0dmax
             d2smat(                            i0nidi,i0nidj) &
                  = d2smat(                     i0nidi,i0nidj) &
                  + d5idtn(i0didk,i0didl,i0nidk,i0nidi,i0nidj) &
                  * d3temp(i0didk,i0didl,i0nidk              )
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_DT

  SUBROUTINE T2EXEC_GV

    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,&
         d4gv,d4igvn
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0node,&
         i0didi,i0didj
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax), &
         d2grad(1:i0dmax,1:i0nmax), &
         d2temp(1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
          
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       DO i0didi = 1, i0dmax
          d2grad(     i0didi,i0nidi) &
               = d4gv(i0didi,       i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
    ENDDO
    
    !C
    !C MAIN LOOP
    !C
    
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
          d2temp(i0didj,i0nidk) = 0.D0
          DO i0didi = 1, i0dmax
             d2temp(              i0didj,i0nidk) &
                  = d2temp(       i0didj,i0nidk) &
                  + d2grad(i0didi,       i0nidk) &
                  * d2jmpm(i0didi,i0didj       ) 
          ENDDO
       ENDDO
    ENDDO

    DO i0nidi =1, i0nmax
    DO i0nidj =1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidk =1, i0nmax
          DO i0didj = 1, i0dmax
             d2smat(                     i0nidi,i0nidj) &
                  = d2smat(              i0nidi,i0nidj) &
                  + d4igvn(i0didj,i0nidk,i0nidi,i0nidj) &
                  * d2temp(i0didj,i0nidk              )
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GV
  
  SUBROUTINE T2EXEC_GT
    
    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,&
         d6gt,d6igtn,d4smat,d2wrks
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,&
         i0didi,i0didj,i0didk,i0didl,&
         i0node
    
    REAL(   i0rkind)::&
         d3grad(1:i0dmax,1:i0dmax,1:i0nmax),&
         d2smat(1:i0nmax,1:i0nmax),&
         d1gtwi(1:i0nmax),&
         d3temp(1:i0dmax,1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
   
    d2smat(1:i0nmax,1:i0nmax) = 0.D0
    
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3grad(     i0didi,i0didj,i0nidi                     ) &
               = d6gt(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
       ENDDO
    ENDDO
    
    DO i0nidi = 1, i0nmax
       d1gtwi(i0nidi) = d2wrks(i0nidi,i0widi)
    ENDDO
    
    !C
    !C MAIN LOOP
    !C

    DO i0nidl = 1, i0nmax
       DO i0didl = 1, i0dmax
       DO i0didk = 1, i0dmax
          d3temp(i0didk,i0didl,i0nidl) = 0.D0
          DO i0didj = 1, i0dmax
          DO i0didi = 1, i0dmax
             d3temp(                     i0didk,i0didl,i0nidl) &
                  = d3temp(              i0didk,i0didl,i0nidl) &
                  + d3grad(i0didi,i0didj,              i0nidl) &
                  * d2jmpm(i0didi,       i0didk              ) &
                  * d2jmpm(       i0didj,       i0didl       )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didl = 1, i0dmax
       DO i0didk = 1, i0dmax
          d2smat(                                   i0nidi,i0nidj) &
               = d2smat(                            i0nidi,i0nidj) &
               + d6igtn(i0didk,i0didl,i0nidk,i0nidl,i0nidi,i0nidj) &
               * d1gtwi(              i0nidk                     ) &   
               * d3temp(i0didk,i0didl,       i0nidl              )         
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO 

    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO       
       
    RETURN
    
  END SUBROUTINE T2EXEC_GT
  
  SUBROUTINE T2EXEC_ES

    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,&
         d3es,d3iesn
         
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,&
         i0node

    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d1exct(1:i0nmax)

    !C
    !C INITIALIZATION
    !C
        
    DO i0nidi = 1,i0nmax
       i0node = i2enr0(i0nidi,1)
       d1exct(i0nidi)&
            = d3es(i0vidi,i0vidj,i0node) &
            * d0jdmp
    ENDDO

    !C
    !C MAIN LOOP
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidk= 1, i0nmax
          d2smat(              i0nidi,i0nidj) & 
               = d2smat(       i0nidi,i0nidj) &
               + d3iesn(i0nidk,i0nidi,i0nidj) &
               * d1exct(i0nidk              )
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO

    RETURN    
        
  END SUBROUTINE T2EXEC_ES
  
  SUBROUTINE T2EXEC_EV
    
    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,&
         d5ev,d5ievn,d2wrks,d4smat
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,&
         i0didi,i0didj,&
         i0node

    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d2exct(1:i0nmax,1:i0dmax),&
         d1evwi(1:i0nmax),&
         d2temp(1:i0nmax,1:i0dmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       DO i0didi = 1, i0dmax
          d2exct(     i0didi,i0nidi)&
               = d5ev(i0didi,i0widi,i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
    ENDDO
    
    DO i0nidi = 1, i0nmax
       d1evwi(i0nidi) = d2wrks(i0nidi,i0widi)
    ENDDO

    !C
    !C MAIN LOOP
    !C 
    
    DO i0nidl = 1, i0nmax   
       DO i0didj = 1, i0dmax
          d2temp(i0didj,i0nidl) = 0.D0     
          DO i0didi = 1, i0dmax
             d2temp(              i0didj,i0nidl) &
                  = d2temp(i0didj,       i0nidl) &
                  + d2exct(i0didi,       i0nidl) &
                  * d2jmpm(i0didi,i0didj       )
          END DO
       END DO
    END DO
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidl = 1, i0nmax
       DO i0nidk = 1, i0nmax
          DO i0didj = 1, i0dmax
             d2smat(                            i0nidi,i0nidj) &
                  = d2smat(                     i0nidi,i0nidj) &
                  + d5ievn(i0didj,i0nidk,i0nidl,i0nidi,i0nidj) &
                  * d1evwi(       i0nidk                     ) &
                  * d2temp(i0didj,       i0nidl              )
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
   
    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_EV
  
  SUBROUTINE T2EXEC_ET
  
    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,&
         d7et,d7ietn,d2wrks
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,&
         i0didi,i0didj,i0didk,i0didl,&
         i0node
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d3exct(1:i0dmax,1:i0dmax,1:i0nmax),&
         d1etwi(1:i0nmax),&
         d1etwj(1:i0nmax),&
         d3temp(1:i0dmax,1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
  
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3exct(     i0didi,i0didj,i0nidi)&
               = d7et(i0didi,i0didj,i0widi,i0widj,i0vidi,i0vidj,i0node) &
               * d0jdmp
       ENDDO
       ENDDO
    ENDDO
    
    DO i0nidi = 1, i0nmax
       d1etwi(i0nidi) = d2wrks(i0nidi,i0widi)
       d1etwj(i0nidi) = d2wrks(i0nidi,i0widi)
    ENDDO
    
    !C
    !C MAIN LOOP
    !C 
    
    DO i0nidm = 1, i0nmax   
       DO i0didl = 1, i0dmax
       DO i0didk = 1, i0dmax
          d3temp(i0didk,i0didl,i0nidm) = 0.D0
          DO i0didj = 1, i0dmax
          DO i0didi = 1, i0dmax
             d3temp(                     i0didk,i0didl,i0nidm) &
                  = d3temp(              i0didk,i0didl,i0nidm) &
                  + d3exct(i0didi,i0didj,              i0nidm) &
                  * d2jmpm(i0didi,       i0didk              ) &
                  * d2jmpm(       i0didj,       i0didl       )
          END DO
          END DO
       END DO
       END DO
    END DO

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidm = 1, i0nmax
       DO i0nidl = 1, i0nmax
       DO i0nidk = 1, i0nmax
          DO i0didl = 1, i0dmax
          DO i0didk = 1, i0dmax
             d2smat(                                          i0nidi,i0nidj) &
                  = d2smat(                                   i0nidi,i0nidj) &
                  + d7ietn(i0didk,i0didl,i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) &
                  * d1etwi(              i0nidk                            ) &
                  * d1etwj(                     i0nidl                     ) &
                  * d3temp(i0didk,i0didl,              i0nidm              )
          ENDDO
          ENDDO
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C STORE SUBMATRIX
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d4smat(       i0nidi,i0nidj,i0vidi,i0vidj) &
            = d4smat(i0nidi,i0nidj,i0vidi,i0vidj) &
            + d2smat(i0nidi,i0nidj              )
    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_ET
  
  SUBROUTINE T2EXEC_SS

    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,&
         d3ss,d2issn,d2svec
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,&
         i0node
 
    REAL(   i0rkind)::&
         d1sour(1:i0nmax),&
         d1svec(1:i0nmax)
    
    !C
    !C INITIALIZATION    
    !C
    
    DO i0nidi = 1, i0nmax
       i0node = i2enr0(i0nidi,1)
       d1sour(i0nidi)&
            = d3ss(i0vidi,i0vidj,i0node)&
            * d0jdmp
    ENDDO
    
    DO i0nidi = 1, i0nmax
       d1svec(i0nidi) = 0.D0
       DO i0nidj = 1, i0nmax
          d1svec(              i0nidi) &
               = d1svec(       i0nidi) &
               + d2issn(i0nidj,i0nidi) &
               * d1sour(i0nidj       )
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_SS
  
  !C
  !C SUBROUTINE FOR STORE SUBMATRIX 
  !C FOR BI-LINEAR RECTANGULAR ELEMENT
  !C
  
  SUBROUTINE T2EXEC_STORE
    
    USE T2COMM,ONLY:&
         i0nmax,i0nmax2,i0vmax,&
         i1nidr,i1nidc,i2hbc,i2enr0,d4smat,d3gmat
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0ng,&
         i0nrl,i0nrc,i0nru,&
         i0nc,i0ncl,i0ncc,i0ncu
    
    !C
    !C
    !C 1Dx1D
    !C
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax         
       i0nrc  = i2enr0(i0nidi,4)
       i0ncc  = i2enr0(i0nidj,4)
       DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
          i0nc = i1nidc(i0ng)
          IF(i0nc.EQ.i0ncc)THEN
             DO i0vidj = 1,3
             DO i0vidi = 1,3
                d3gmat(                     i0vidi,i0vidj,i0ng) &
                     = d3gmat(              i0vidi,i0vidj,i0ng) &
                     + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     )
             ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C
    !C 1Dx2D
    !C
    !C
   
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       
       i0nrc  = i2enr0(i0nidi,3)
       i0ncc  = i2enr0(i0nidj,2)
       
       IF(    i0ncc.LE.i0nmax2)THEN
          
          DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
             i0nc = i1nidc(i0ng)
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 1, 3
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF(i0ncc.GT.i0nmax2)THEN
          
          i0ncc  = i0ncc-i0nmax2
          i0ncl  = i2hbc(i0ncc,1)
          i0ncu  = i2hbc(i0ncc,2)
          
          DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
             i0nc = i1nidc(i0ng )
             IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 1, 3
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ENDIF
       
    ENDDO
    ENDDO
    
    !C
    !C
    !C 2Dx1D
    !C
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax

       i0nrc  = i2enr0(i0nidi,2)
       i0ncc  = i2enr0(i0nidj,4)

       IF(    i0nrc.LE.i0nmax2)THEN

          DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
             i0nc = i1nidc(i0ng)
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     )
                ENDDO
                ENDDO
             ENDIF
          ENDDO

       ELSEIF(i0nrc.GT.i0nmax2)THEN
       
          i0nrc  = i0nrc-i0nmax2
          i0nrl  = i2hbc(i0nrc,1)
          i0nru  = i2hbc(i0nrc,2)
          
          DO i0ng = i1nidr(i0nrl), i1nidr(i0nrl+1)-1
             i0nc = i1nidc(i0ng) 
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i0ng = i1nidr(i0nru), i1nidr(i0nru+1)-1
             i0nc = i1nidc(i0ng) 
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO

       ENDIF
    ENDDO
    ENDDO
      
    !C
    !C 
    !C 2Dx2D
    !C
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       
       i0nrc  = i2enr0(i0nidi,2)
       i0ncc  = i2enr0(i0nidj,2)
       
       IF((i0nrc.LE.i0nmax2).AND.(i0ncc.LE.i0nmax2))THEN
          
          DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
             i0nc = i1nidc(i0ng )
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
       ELSEIF((i0nrc.LE.i0nmax2).AND.(i0ncc.GT.i0nmax2))THEN
          
          i0ncc  = i0ncc-i0nmax2
          i0ncl  = i2hbc(i0ncc,1)
          i0ncu  = i2hbc(i0ncc,2)
          
          DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
             i0nc = i1nidc(i0ng )
             IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i0nrc.GT.i0nmax2).AND.(i0ncc.LE.i0nmax2))THEN
          
          i0nrc  = i0nrc-i0nmax2
          i0nrl  = i2hbc(i0nrc,1)
          i0nru  = i2hbc(i0nrc,2)
          
          DO i0ng = i1nidr(i0nrl), i1nidr(i0nrl+1)-1
             i0nc = i1nidc(i0ng ) 
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
          DO i0ng = i1nidr(i0nru), i1nidr(i0nru+1)-1
             i0nc = i1nidc(i0ng ) 
             IF(i0nc.EQ.i0ncc)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i0nrc.GT.i0nmax2).AND.(i0ncc.GT.i0nmax2))THEN
          
          i0nrc  = i0nrc-i0nmax2
          i0nrl  = i2hbc(i0nrc,1)
          i0nru  = i2hbc(i0nrc,2)
          
          i0ncc  = i0ncc-i0nmax2
          i0ncl  = i2hbc(i0ncc,1)
          i0ncu  = i2hbc(i0ncc,2)
          
          DO i0ng = i1nidr(i0nrl), i1nidr(i0nrl+1)-1
             i0nc = i1nidc(i0ng )
             IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i0ng = i1nidr(i0nru), i1nidr(i0nru+1)-1
             i0nc = i1nidc(i0ng )
             IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3gmat(                     i0vidi,i0vidj,i0ng) &
                        = d3gmat(              i0vidi,i0vidj,i0ng) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj     ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ENDIF

    ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_STORE

END MODULE T2EXEC
