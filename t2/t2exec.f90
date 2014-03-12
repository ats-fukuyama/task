!C--------------------------------------------------------------------
!C
!C T2EXEC
!C
!C                       2024-02-22 H.SETO
!C
!C
!C--------------------------------------------------------------------
MODULE T2EXEC
  
  USE T2CNST, ONLY:&
       i0ikind,i0rkind
  
  IMPLICIT NONE
  
  INTEGER(i0ikind)::&
       i0vidi,i0vidj,i0widi,i0widj,i0eidi
  
  PUBLIC T2_EXEC
  
  PRIVATE  
  
CONTAINS
  
  !C-------------------------------------------------------------------
  !C
  !C MAIN ROUTINE OF T2EXEC
  !C FEM SOLVER FOR SIMULTANEOUS ADVECTION-DIFFUSION EQUATIONS
  !C 
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------
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
          
          !IF(i0ssvt.EQ.1)       CALL T2EXEC_SS
          
       ENDDO
       ENDDO

       CALL T2EXEC_STORE
       
    ENDDO
    
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC: completed:          cpu=', &
         e0time_1-e0time_0,' [s]'

    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_BCOND
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_BCOND completed:     cpu=', &
         e0time_1-e0time_0,' [s]'

    CALL CPU_TIME(e0time_0)
    CALL T2EXEC_SOLVE
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- T2EXEC_SOLVE completed:     cpu=', &
         e0time_1-e0time_0,' [s]'

    RETURN
    
  END SUBROUTINE T2_EXEC
  
  !C
  !C MODIFIED 2014-03-09
  !C
  !C
  SUBROUTINE T2EXEC_SOLVE
    
    USE T2COMM, ONLY:&
         i0vmax,i0bmax,i0xmax,i0amax,i0lmax,i0bvmax,i0avmax,&
         i0dbg,i2vvvt,i2hbc,idebug,&
         i1nidr,i1nidc,i1pdn2,i1rdn2,&
         d3amat,d2bvec,d2xvec,d2xvec_befor,d2xvec_after
    
    USE LIBMPI
    USE COMMPI
    USE LIBMTX
    

    INTEGER(i0ikind)::istart,iend,its
    INTEGER(i0ikind)::itype, m1, m2
    REAL(   i0rkind)::tolerance,d0val
    REAL(   i0rkind),POINTER,SAVE::x(:)
    INTEGER(i0ikind)::&
         i0nr,i0nc,i0pdn2,i0rdn2
    INTEGER(i0ikind)::&
         i0lidi,i0ridi,i0pidi,i0aidi,i0bidi,&
         i0xidi,i0xidc,i0xidd,i0xidu,i0xrd,i0xru,&
         i0vvvt,i0offset,&
         i0br,i0xr,i0ar,i0ac,&
         i0arc,i0arc1,i0arc2,i0arc3,&
         i0acc,i0acc1,i0acc2,i0acc3,i0acl,i0acl1,i0acl2,i0acl3

100 FORMAT(A5,I3,A5,I3,A5,I3,A5,D15.6,A5,D15.6)

    !C
    !C
    !C MTX SETUP
    !C
    !C

    itype = 0
    m1    = 4
    
    IF(nsize.EQ.1)THEN
       m2 = 5
    ELSE
       m2 = 0
    ENDIF
    
    tolerance=1.D-7
    
    i0bvmax = i0bmax*i0vmax
    i0avmax = i0amax*i0vmax*i0vmax
    
    ALLOCATE(x(i0bvmax))
    
    CALL MTX_SETUP(i0bvmax,istart,iend,idebug=0)
    !CALL MTX_SETUP(i0bvmax,istart,iend,nzmax=i0avmax,idebug=0)
    
    !C 
    !C STORE GLOBAL STIFFNESS MATRIX  
    !C 
    DO i0nr = 1, i0bmax   
       DO i0aidi = i1nidr(i0nr), i1nidr(i0nr+1)-1
          i0nc = i1nidc(i0aidi) 
          DO i0vidj = 1, i0vmax
          DO i0vidi = 1, i0vmax
             i0vvvt = i2vvvt(i0vidi,i0vidj)
             IF(i0vvvt.EQ.1) THEN
                d0val = d3amat(i0vidi,i0vidj,i0aidi)
                i0ar  = i0vmax*(i0nr-1) + i0vidi
                i0ac  = i0vmax*(i0nc-1) + i0vidj
                CALL MTX_SET_MATRIX(i0ar,i0ac,d0val)
                IF(IDEBUG.EQ.1) THEN
                   WRITE(18,'(2I5,I10,2I3,2I7,1PE12.4)') &
                        i0nr,i0nc,i0aidi,i0vidi,i0vidj,i0ar,i0ac,d0val
                END IF

             END IF
             
          ENDDO
          ENDDO
       ENDDO
    ENDDO

    !C
    !C ADDITIONAL COMPONENTS FOR FLUX SURFACE AVERAGING
    !C
    
    i0offset  = 1

    DO i0lidi = 1, i0lmax
       
       i0pdn2 = i1pdn2(i0lidi)
       i0rdn2 = i1rdn2(i0lidi)
       
       DO i0ridi = 1, i0rdn2
          DO i0pidi = 1, i0pdn2
                
             i0arc  = i0vmax*i0offset
             i0arc1 = i0arc + 1
             i0arc2 = i0arc + 2
             i0arc3 = i0arc + 3
                
             i0acc  = i0arc
             i0acc1 = i0arc1
             i0acc2 = i0arc2
             i0acc3 = i0arc3
                
             i0acl  = i0vmax*(i0offset-1)
             i0acl1 = i0acl + 1
             i0acl2 = i0acl + 2
             i0acl3 = i0acl + 3
                
             IF((i0pidi.GE.1).AND.(i0pidi.LT.i0pdn2))THEN
                CALL MTX_SET_MATRIX(i0arc1,i0acc1,-1.D0)
                CALL MTX_SET_MATRIX(i0arc2,i0acc2,-1.D0)
                CALL MTX_SET_MATRIX(i0arc3,i0acc3,-1.D0)
             ENDIF
             
             IF((i0pidi.GT.1).AND.(i0pidi.LE.i0pdn2))THEN!
                CALL MTX_SET_MATRIX(i0arc1,i0acl1, 1.D0)
                CALL MTX_SET_MATRIX(i0arc2,i0acl2, 1.D0)
                CALL MTX_SET_MATRIX(i0arc3,i0acl3, 1.D0)
            ENDIF
             
             i0offset = i0offset + 1
             
          ENDDO
       ENDDO
    ENDDO
        
    !C
    !C SET GLOBAL RIGHT HAND SIDE VECTOR
    !C
    
    DO i0bidi = 1, i0bmax
       DO i0vidi = 1, i0vmax
          i0br  = i0vmax*(i0bidi-1) + i0vidi
          d0val = d2bvec(i0vidi,i0bidi)
          CALL MTX_SET_SOURCE(i0br,d0val)
       ENDDO
    ENDDO
    
    !C
    !C SET GLOBAL RIGHT HAND SIDE VECTOR
    !C
        
    DO i0xidi = 1, i0bmax
       DO i0vidi = 1, i0vmax
          i0xr  = i0vmax*(i0xidi-1) + i0vidi
          d0val = d2xvec_befor(i0vidi,i0xidi)
          CALL MTX_SET_VECTOR(i0xr,d0val)
       ENDDO
    ENDDO
    
    CALL MTX_SOLVE(itype,tolerance,its,&
         methodKSP=m1,methodPC=m2,max_steps = 999)
    
    CALL MTX_GATHER_VECTOR(x)
    
    DO i0xidi = 1, i0xmax
       IF(    i0xidi.LE.i0bmax)THEN
          DO i0vidi = 1, i0vmax
             i0xr = i0vmax*(i0xidi - 1) + i0vidi
             d2xvec_after(i0vidi,i0xidi) = x(i0xr)
          ENDDO
       ELSEIF(i0xidi.GT.i0bmax)THEN
          i0xidc = i0xidi - i0bmax
          i0xidd = i2hbc(1,i0xidc)
          i0xidu = i2hbc(2,i0xidc)
          DO i0vidi = 1, i0vmax 
             i0xrd = i0vmax*(i0xidd-1) + i0vidi
             i0xru = i0vmax*(i0xidu-1) + i0vidi
             d2xvec_after(i0vidi,i0xidi) = 0.5D0*(x(i0xrd)+x(i0xru))
          ENDDO
       ENDIF
    ENDDO
    
    CALL MTX_CLEANUP
    
    DEALLOCATE(x)
    
    RETURN
    
  END SUBROUTINE T2EXEC_SOLVE


  !C------------------------------------------------------------------ 
  !C  SUBROUTINE T2EXEC_SET_LOCAL_VARIABLES
  !C  　
  !C    * SET LOCAL NODE-ELEMENT GRAPH
  !C    * CALCULATE JACOBIAN OF PARAMETRIC SPACE
  !C    * STORE VARIABLES ST N-th TIMESTEP   
  !C    * SET WORKING ARRAY FOR DIFFERENTIAL
  !C    * SET STABILIZATION FACTORS FOR SUPG
  !C
  !C                     2014-02-07 H.SETO
  !C
  !C------------------------------------------------------------------  
  SUBROUTINE T2EXEC_LV
    
    USE T2COMM, ONLY:&
         i0nmax,i0dmax,i0vmax,i0smax,&
         d1rsiz,d1psiz,d2ws,d2xvec,d2wrks,d2kwns,d2mtrc,&
         d4smat,d2svec,d4av,&
         i2enr0,i3enr,i1mlel,d0rmnr,&
         d2jmpm,d0jdmp,dt,d2jm1,d2mfc1!,i2avvt,d3eafv,i0supg
    
    INTEGER(i0ikind)::&
         i0midi,i0didi,i0didj,i0ridi,i0nidi,i0nidj,i0nidk,&
         i0sidi,i0lidi,i0bidi,i0xidi!,i0avvt
    
    REAL(   i0rkind)::&
         d0supg,d0rsiz,d0psiz,d0sqrtg,d0sqrtgi,&
         d0mtrc,d0eafv,d0area,&
         d1temp(1:i0dmax)
    
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
    !C checked 2014-02-04 H.SETO
    !C 
    i0nidi = i2enr0(1,1)
    i0nidj = i2enr0(2,1)
    i0nidk = i2enr0(4,1)

    d0rsiz = d2mfc1(1,i0nidj)-d2mfc1(1,i0nidi)
    d0psiz = d2mfc1(2,i0nidk)-d2mfc1(2,i0nidi)
    d0rsiz = ABS(d0rsiz)
    d0psiz = ABS(d0psiz)
    d0jdmp = d0rsiz*d0psiz*0.25D0

    IF(d0jdmp.LE.0.D0)THEN
       WRITE(6,'("ERROR:: D1JDMP IS SINGULAR")')
       STOP
    ENDIF

    d2jmpm(1,1) = 2.D0/d0rsiz
    d2jmpm(1,2) = 0.D0
    d2jmpm(2,1) = 0.D0
    d2jmpm(2,2) = 2.D0/d0psiz
    
    !C
    !C STORE DEPENDENT VARIABLES AT N-th TIMESTEP 
    !C 
    
    DO i0nidi = 1, i0nmax
       i0bidi = i2enr0(i0nidi,4)
       DO i0vidi = 1, 3
          d2kwns(i0nidi,i0vidi) = d2xvec(i0vidi,i0bidi)
       ENDDO
       
       i0xidi = i2enr0(i0nidi,2)
       DO i0vidi = 4, i0vmax
          d2kwns(i0nidi,i0vidi) = d2xvec(i0vidi,i0xidi)
       ENDDO

    ENDDO
    
    !C
    !C STORE WORKING ARRAY FOR DIFFERENTIAL
    !C

    !C 
    !C D2WRKS(1   ,:) : B    AT L-TH PICARD ITERATION
    !C D2WRKS(2   ,:) : R    AT L-TH PICARD ITERATION
    !C D2WRKS(2N+1,:) : Ub   AT L-TH PICARD ITERATION 
    !C D2WRKS(2N+2,:) : Qb/P AT L-TH PICARD ITERATION 
    !C
    
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       d2wrks(i0nidi,1) = d2ws(1,i0midi)
       d2wrks(i0nidi,2) = d2ws(2,i0midi)
       DO i0sidi = 1, i0smax
          i0widi = 2*i0sidi        
          d2wrks(i0nidi,i0widi+1) = d2ws(i0widi+1,i0midi)
          d2wrks(i0nidi,i0widi+2) = d2ws(i0widi+2,i0midi)
       ENDDO
    ENDDO
    
    !C
    !C INITIALIZATION OF SUBMATRIX AND SUBVETOR
    !C
    
    d4smat(1:i0nmax,1:i0nmax,1:i0vmax,1:i0vmax) = 0.D0
    d2svec(1:i0nmax,1:i0vmax) = 0.D0
    
    RETURN

  END SUBROUTINE T2EXEC_LV
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF MASS SUBMATRCES : M
  !C
  !C                     2014-02-08 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2EXEC_MS
    
    USE T2COMM,ONLY:&
         i0supg,i0nmax,i0dmax,i0vmax,i2enr0,&
         d0jdmp,d2jmpm,dt,d4smat,d2svec,d2kwns,&
         d3ms,d3imsn!,d3eafv,d5imss,i2avvt
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,&
         i0midi
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d1svec(1:i0nmax),&
         d1mass(1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       d1mass(                   i0nidi) &
            = d3ms(i0vidi,i0vidj,i0midi) &
            * d0jdmp/dt
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
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF ADVECTION SUBMATRCES: V1
  !C
  !C                     2014-02-12 H.SETO
  !C
  !C-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_AV
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,d0jdmp,d2jmpm,&
         d4av,d4iavn,d4smat
    
    INTEGER(i0ikind)::&
         i0didi,i0didj,&
         i0nidi,i0nidj,i0nidk,&
         i0midi
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d2velo(1:i0dmax,1:i0nmax),&
         d2temp(1:i0dmax,1:i0nmax)
    
    !C
    !C INTITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       DO i0didi = 1, i0dmax
          d2velo(     i0didi,              i0nidi) &
               = d4av(i0didi,i0vidi,i0vidj,i0midi) &
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
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF ADVECTION SUBMATRCES: V2
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_AT
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,d2wrks,d6at,d6iatn
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,&
         i0didi,i0didj,i0didk,i0didl,&
         i0midi
    
    REAL(   i0rkind)::&
         d3velo(1:i0dmax,1:i0dmax,1:i0nmax),&
         d2smat(1:i0nmax,1:i0nmax),&
         d1atwi(1:i0nmax),&
         d3temp(1:i0dmax,1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
      d1atwi(i0nidi) = d2wrks(i0nidi,i0widi)
    ENDDO
    
    DO i0nidi = 1, i0nmax       
       i0midi = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3velo(     i0didi,i0didj,                     i0nidi) &
               = d6at(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0midi) &
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
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
       
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidl = 1, i0nmax
       DO i0nidk = 1, i0nmax
          DO i0didl = 1, i0dmax
          DO i0didk = 1, i0dmax
             d2smat(                                   i0nidi,i0nidj) &
                  = d2smat(                            i0nidi,i0nidj) &
                  + d6iatn(i0didk,i0didl,i0nidk,i0nidl,i0nidi,i0nidj) &
                  * d3temp(i0didk,i0didl,       i0nidl              ) &
                  * d1atwi(              i0nidk                     )
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

  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF DIFFUSION SUBMATRCES: V1
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_DT
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,d5dt,d5idtn
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,&
         i0didi,i0didj,i0didk,i0didl,&
         i0midi

    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d3diff(1:i0dmax,1:i0dmax,1:i0nmax), &
         d3temp(1:i0dmax,1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3diff(     i0didi,i0didj,              i0nidi)&
               = d5dt(i0didi,i0didj,i0vidi,i0vidj,i0midi) &
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

  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF GRADIENT SUBMATRCES: A1
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GV
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,d4igvn,d4gv
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,&
         i0didi,i0didj,&
         i0midi
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax), &
         d2grad(1:i0dmax,1:i0nmax), &
         d2temp(1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
          
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       DO i0didi = 1, i0dmax
          d2grad(     i0didi,              i0nidi) &
               = d4gv(i0didi,i0vidi,i0vidj,i0midi) &
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

    DO i0nidi = 1, i0nmax
    DO i0nidj = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidk = 1, i0nmax
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
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF GRADIENT SUBMATRCES: A2
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------  
  SUBROUTINE T2EXEC_GT
    
    USE T2COMM,ONLY:&
         i0nmax,i0dmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,d2wrks,d6igtn,d6gt
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,&
         i0didi,i0didj,i0didk,i0didl,&
         i0midi
    
    REAL(   i0rkind)::&
         d3grad(1:i0dmax,1:i0dmax,1:i0nmax),&
         d2smat(1:i0nmax,1:i0nmax),&
         d1gtwi(1:i0nmax),&
         d3temp(1:i0dmax,1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
       
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3grad(     i0didi,i0didj,                     i0nidi) &
               = d6gt(i0didi,i0didj,i0widi,i0vidi,i0vidj,i0midi) &
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
       d2smat(i0nidi,i0nidj) = 0.D0
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
  
  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF EXCITATION SUBMATRCES: C1
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------    
  SUBROUTINE T2EXEC_ES

    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,d3iesn,d3es
         
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,&
         i0didi,&
         i0midi

    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d1exct(1:i0nmax)

    !C
    !C INITIALIZATION
    !C
        
    DO i0nidi = 1,i0nmax
       i0midi = i2enr0(i0nidi,1)
       d1exct(                   i0nidi) &
            = d3es(i0vidi,i0vidj,i0midi) &
            * d0jdmp
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

  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF EXCITATION SUBMATRCES: C2
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_EV
    
    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d2wrks,d4smat,d5ievn,d5ev
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,&
         i0didi,i0didj,&
         i0midi!,i0avvt
    
    REAL(   i0rkind)::&
         d2smat(1:i0nmax,1:i0nmax),&
         d2exct(1:i0dmax,1:i0nmax),&
         d1evwi(1:i0nmax),&
         d2temp(1:i0dmax,1:i0nmax)
    
    !C
    !C INITIALIZATION
    !C
    
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       DO i0didi = 1, i0dmax
          d2exct(     i0didi,                     i0nidi) &
               = d5ev(i0didi,i0widi,i0vidi,i0vidj,i0midi) &
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
                  = d2temp(       i0didj,i0nidl) &
                  + d2exct(i0didi,       i0nidl) &
                  * d2jmpm(i0didi,i0didj       )
          ENDDO
       ENDDO
    ENDDO
    
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

  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF EXCITATION SUBMATRCES: C3
  !C
  !C                     2014-01-30 H.SETO
  !C
  !C-------------------------------------------------------------------   
  SUBROUTINE T2EXEC_ET
    
    USE T2COMM,ONLY:&
         i0dmax,i0nmax,i0vmax,i2enr0,&
         d2jmpm,d0jdmp,d4smat,d2wrks,d7ietn,d7et
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,i0nidl,i0nidm,i0nidn,&
         i0didi,i0didj,i0didk,i0didl,&
         i0midi
    
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
       i0midi = i2enr0(i0nidi,1)
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d3exct(     i0didi,i0didj,i0nidi)&
               = d7et(i0didi,i0didj,i0widi,i0widj,&
               &      i0vidi,i0vidj,i0midi)&
               * d0jdmp
       ENDDO
       ENDDO
    ENDDO
    
    DO i0nidi = 1, i0nmax
       d1etwi(i0nidi) = d2wrks(i0nidi,i0widi)
       d1etwj(i0nidi) = d2wrks(i0nidi,i0widj)
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
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       d2smat(i0nidi,i0nidj) = 0.D0
       DO i0nidm = 1, i0nmax
       DO i0nidl = 1, i0nmax
       DO i0nidk = 1, i0nmax
          DO i0didl = 1, i0dmax
          DO i0didk = 1, i0dmax
             d2smat(                            i0nidi,i0nidj) &
                  = d2smat(                     i0nidi,i0nidj) &
                  + d7ietn(i0didk,i0didl,&
                  &        i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) &
                  * d1etwi(i0nidk                            ) &
                  * d1etwj(       i0nidl                     ) &
                  * d3temp(i0didk,i0didl,i0nidm              )
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
         d2jmpm,d0jdmp,d2svec,d2issn,d3ss

    INTEGER(i0ikind)::&
         i0nidi,i0nidj,i0nidk,&
         i0didi,&
         i0midi
 
    REAL(   i0rkind)::&
         d1sour(1:i0nmax),&
         d1svec(1:i0nmax)

    !C
    !C INITIALIZATION    
    !C
    
    DO i0nidi = 1, i0nmax
       i0midi = i2enr0(i0nidi,1)
       d1sour(i0nidi)&
            = d3ss(i0vidi,i0vidj,i0midi) * d0jdmp
    ENDDO

    !C
    !C MAIN LOOP
    !C

    DO i0nidi = 1, i0nmax
       d1svec(i0nidi) = 0.D0
       DO i0nidj = 1, i0nmax
          d1svec(              i0nidi) &
               = d1svec(       i0nidi) &
               + d2issn(i0nidj,i0nidi) &
               * d1sour(i0nidj       )
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
    
  END SUBROUTINE T2EXEC_SS
  
  !C-------------------------------------------------------------------
  !C
  !C SUBROUTINE FOR STORE SUBMATRIX 
  !C FOR BI-LINEAR RECTANGULAR ELEMENT
  !C
  !C
  !C
  !C-------------------------------------------------------------------
  
  SUBROUTINE T2EXEC_STORE
    
    USE T2COMM,ONLY:&
         i0nmax,i0vmax,i0bmax,&
         i1nidr,i1nidc,i2hbc,i2enr0,&
         d4smat,d3amat,d2svec,d2bvec
    
    INTEGER(i0ikind)::&
         i0nidi,i0nidj,&
         i0aidi,&
         i0bidi,i0bidj,i0bidk,&
         i0hidi,i0hidj,i0hidk,i0hidl,&
         i0xidi,i0xidj


    !C
    !C
    !C MATRIX
    !C
    !C
    
    !C
    !C 1Dx1D
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax         
       i0bidi  = i2enr0(i0nidi,4)
       i0bidj  = i2enr0(i0nidj,4)
       DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
          i0bidk = i1nidc(i0aidi)
          IF(i0bidk.EQ.i0bidj)THEN
             DO i0vidj = 1,3
             DO i0vidi = 1,3
                d3amat(                     i0vidi,i0vidj,i0aidi) &
                     = d3amat(              i0vidi,i0vidj,i0aidi) &
                     + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
             ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    
    !C
    !C 1Dx2D
    !C
   
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
              
       i0bidi = i2enr0(i0nidi,3)
       i0xidj = i2enr0(i0nidj,2)
       
       IF(    i0xidj.LE.i0bmax)THEN
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 1, 3
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF(i0xidj.GT.i0bmax)THEN
          i0xidj = i0xidj - i0bmax
          i0hidi = i2hbc(1,i0xidj)
          i0hidj = i2hbc(2,i0xidj)
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidi).OR.(i0bidk.EQ.i0hidj))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 1, 3
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.5D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    
    !C
    !C 2Dx1D
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       i0xidi = i2enr0(i0nidi,2)
       i0bidj = i2enr0(i0nidj,4)
       IF(    i0xidi.LE.i0bmax)THEN
          DO i0aidi = i1nidr(i0xidi), i1nidr(i0xidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF(i0bidk.EQ.i0bidj)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       ELSEIF(i0xidi.GT.i0bmax)THEN
       
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          
          DO i0aidi = i1nidr(i0hidi), i1nidr(i0hidi+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0bidj)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i0aidi = i1nidr(i0hidj), i1nidr(i0hidj+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0bidj)THEN
                DO i0vidj = 1, 3
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
       
       ENDIF
    ENDDO
    ENDDO
      
    !C 
    !C 2Dx2D
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
       
       i0xidi = i2enr0(i0nidi,2)
       i0xidj = i2enr0(i0nidj,2)
       
       IF((i0xidi.LE.i0bmax).AND.(i0xidj.LE.i0bmax))THEN
          
          DO i0aidi = i1nidr(i0xidi), i1nidr(i0xidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       )
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
       ELSEIF((i0xidi.LE.i0bmax).AND.(i0xidj.GT.i0bmax))THEN
          
          i0xidj = i0xidj - i0bmax
          i0hidi = i2hbc(1,i0xidj)
          i0hidj = i2hbc(2,i0xidj)
          
          DO i0aidi = i1nidr(i0xidi), i1nidr(i0xidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidi).OR.(i0bidk.EQ.i0hidj))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i0xidi.GT.i0bmax).AND.(i0xidj.LE.i0bmax))THEN
          
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          
          DO i0aidi = i1nidr(i0hidi), i1nidr(i0hidi+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
             
          DO i0aidi = i1nidr(i0hidj), i1nidr(i0hidj+1)-1
             i0bidk = i1nidc(i0aidi) 
             IF(i0bidk.EQ.i0xidj)THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.50D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ELSEIF((i0xidi.GT.i0bmax).AND.(i0xidj.GT.i0bmax))THEN
          
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          
          i0xidj = i0xidj - i0bmax
          i0hidk = i2hbc(1,i0xidj)
          i0hidl = i2hbc(2,i0xidj)
          
          DO i0aidi = i1nidr(i0hidi), i1nidr(i0hidi+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidk).OR.(i0bidk.EQ.i0hidl))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
          DO i0aidi = i1nidr(i0hidj), i1nidr(i0hidj+1)-1
             i0bidk = i1nidc(i0aidi)
             IF((i0bidk.EQ.i0hidk).OR.(i0bidk.EQ.i0hidl))THEN
                DO i0vidj = 4, i0vmax
                DO i0vidi = 4, i0vmax
                   d3amat(                     i0vidi,i0vidj,i0aidi) &
                        = d3amat(              i0vidi,i0vidj,i0aidi) &
                        + d4smat(i0nidi,i0nidj,i0vidi,i0vidj       ) &
                        * 0.25D0
                ENDDO
                ENDDO
             ENDIF
          ENDDO
          
       ENDIF
       
    ENDDO
    ENDDO
    
    !C
    !C
    !C RIGHT HANDSIDE VECTOR
    !C
    !C
    
    !C
    !C 1D
    !C
    DO i0nidi = 1, i0nmax
       i0bidi = i2enr0(i0nidi,4)
       DO i0vidi = 1, 3
          d2bvec(              i0vidi,i0bidi) &
               = d2bvec(       i0vidi,i0bidi) &
               + d2svec(i0nidi,i0vidi       )
       ENDDO
    ENDDO
    
    !C
    !C 2D
    !C
    
    DO i0nidi = 1, i0nmax
       i0xidi = i2enr0(i0nidi,2)
       IF(    i0xidi.LE.i0bmax)THEN
          DO i0vidi = 4, i0vmax
             d2bvec(              i0vidi,i0xidi) &
                  = d2bvec(       i0vidi,i0xidi) &
                  + d2svec(i0nidi,i0vidi       )
          ENDDO
       ELSEIF(i0xidi.GT.i0bmax)THEN
          i0xidi = i0xidi - i0bmax
          i0hidi = i2hbc(1,i0xidi)
          i0hidj = i2hbc(2,i0xidi)
          DO i0vidi = 4, i0vmax
             d2bvec(              i0vidi,i0hidi) &
                  = d2bvec(       i0vidi,i0hidi) &
                  + d2svec(i0nidi,i0vidi       ) &
                  * 0.50D0
             d2bvec(              i0vidi,i0hidj) &
                  = d2bvec(       i0vidi,i0hidj) &
                  + d2svec(i0nidi,i0vidi       ) &
                  * 0.50D0
          ENDDO
       ENDIF
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_STORE
 
  !C-------------------------------------------------------------------
  !C 
  !C BOUNDARY CONDITIONS
  !C
  !C                     2014-02-22 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2EXEC_BCOND
    
    USE T2COMM, ONLY:&
         i0solv,i0lmax,i0bmax,i0vmax,&
         i1pdn2,i1nidr,i1nidc,&
         d3amat,d2bvec,d2xvec
    
    INTEGER::&
         i0aidi,&
         i0bidi,i0bidj,&
         i0vidi,i0vidj,&
         i0bsta,i0bend,i0pdn2
    
    i0pdn2 = i1pdn2(i0lmax)
    i0bsta = i0bmax - i0pdn2 + 1
    i0bend = i0bmax

    !C
    !C
    !C SET FIXED VARIALES 
    !C 
    !C 
    
    SELECT CASE(i0solv)
       
       !C i0solv = 1
       !C SOLVE ONLY ELECTRON DENSITY AND MOMENTUM
       !C
       
    CASE(1)
       
       !C
       !C
       !C SET FIXED VALUES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,11:)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,11:)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
       ENDDO
       
       !C
       !C i0solv = 2
       !C SOLVE ONLY ELECTRON AND
       !C            ION DENSITIES AND MOMENTUMS
       !C
       
    CASE(2)
       
       !C
       !C
       !C SET FIXED VALUES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,11:15,21:25)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,11:15,21:25)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
          
       ENDDO

       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
       
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:16,20:25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:16,20:)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
          
       ENDDO
       
       !C
       !C i0slov = 3
       !C SOLVE ELECTRON DENSITY, MOMENTUMS,
       !C            AND PRESSURE
       !C
       
    CASE(3)

       !C
       !C
       !C SET FIXED VALUES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:11,16:21)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:11,16:21)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:11,15:21,25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                END SELECT
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:11,15:21,25)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
       ENDDO
       
       !C
       !C i0solv = 4: 
       !C SOLVE ELECTRON AND ION: MOMENTUMS AND HEATFLUXES
       !C
       
    CASE(4)
       
       !C
       !C
       !C SET FIXED VARIABLES
       !C
       !C
       
       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,11,16,21)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,11,16,21)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C SET DIRICHLET CONDITION 
       !C
    
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:11,15:16,20:21,25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
        
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:11,15:16,20:21,25)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             END SELECT
          ENDDO
       ENDDO

       !C
       !C i0solv = 5: 
       !C SOLVE ELECTRON AND ION DENSITIES, 
       !C                MOMENTUMS AND HEATFLUXES
       !C

    CASE(5)

       !C
       !C
       !C SET FIXED VARIABLES
       !C
       !C

       DO i0bidi= 1, i0bmax
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidj = 1, i0vmax
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:5,11,21)
                   IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                      d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                   ELSE
                      d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                   ENDIF
                CASE DEFAULT
                   CYCLE
                ENDSELECT
             ENDDO
             ENDDO
          ENDDO
          
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:5,11,21)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION 
       !C
       !C
       
       DO i0bidi = i0bsta, i0bend
          
          !C
          !C STIFFNESS MATRIX
          !C
          
          DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
             i0bidj = i1nidc(i0aidi)
             DO i0vidi = 1, i0vmax
                SELECT CASE(i0vidi)
                CASE(1:6,10:11,15:16,20:21,25)
                   CYCLE
                CASE DEFAULT
                   DO i0vidj = 1, i0vmax
                      IF((i0bidi.EQ.i0bidj).AND.(i0vidi.EQ.i0vidj))THEN
                         d3amat(i0vidi,i0vidj,i0aidi) = 1.D0
                      ELSE
                         d3amat(i0vidi,i0vidj,i0aidi) = 0.D0
                      ENDIF
                   ENDDO
                ENDSELECT
             ENDDO
          ENDDO
        
          !C
          !C RHS VECTOR 
          !C
          
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,10:11,15:16,20:21,25)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
       ENDDO
    CASE DEFAULT
       WRITE(6,*)'INCORRECT I0SOLV'
    ENDSELECT
    
    RETURN
    
  END SUBROUTINE T2EXEC_BCOND
END MODULE T2EXEC
