!C--------------------------------------------------------------------
!C
!C T2EXEC
!C
!C                       2024-03-26 H.SETO
!C
!C
!C--------------------------------------------------------------------
MODULE T2EXEC
  
  USE T2CNST, ONLY:&
       i0ikind,i0rkind,ikind,rkind
  
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
  !C                     2014-05-20 H.SETO
  !C
  !C-------------------------------------------------------------------
  SUBROUTINE T2_EXEC
    
    USE T2COMM,ONLY:&
         nnmax,nmmax,nemax,nkmax,nvmax,&
         i2msvt,i2avvt,i2atvt,i2dtvt,i2gvvt,i2gtvt,&
         i2esvt,i2evvt,i2etvt,i2ssvt,&
         i3atwt,i3gtwt,i3evwt,i4etwt
    
    INTEGER(ikind)::&
         i_v,j_v,i_k,j_k,i_e
    
    INTEGER(ikind)::&
         i0msvt,i0avvt,i0atvt,i0dtvt,i0gvvt,i0gtvt,&
         i0esvt,i0evvt,i0etvt,i0ssvt,&
         i0atwt,i0gtwt,i0evwt,i0etwt

!    LOGICAL::&
!         hasMassSca, hasAdveVec, hasAdveTen, hasDiffTen,&
!         hasGradVec, hasGradTen, hasExctSca, hasExctVec,&
!         hasExctTen, hasSourSca
    
    REAL(4)::e0time_0,e0time_1
    

    CALL CPU_TIME(e0time_0)

    DO i_e = 1, nemax
             
       CALL T2EXEC_INITIALIZE(i_e)
       
       DO j_v = 1, nvmax
       DO i_v = 1, nvmax
          
          i0msvt = i2msvt(i_v,j_v)
          i0avvt = i2avvt(i_v,j_v)
          i0atvt = i2atvt(i_v,j_v)
          i0dtvt = i2dtvt(i_v,j_v)
          i0gvvt = i2gvvt(i_v,j_v)
          i0gtvt = i2gtvt(i_v,j_v)
          i0esvt = i2esvt(i_v,j_v)
          i0evvt = i2evvt(i_v,j_v)
          i0etvt = i2etvt(i_v,j_v)
          i0ssvt = i2ssvt(i_v,j_v)


          
          IF(i0msvt.EQ.1)       CALL T2EXEC_MS(i_v,j_v)
          
!          IF(i0avvt.EQ.1)       CALL T2EXEC_AV(i_v,j_v)
          
          IF(i0atvt.EQ.1)THEN
             DO i_k = 1, nkmax
                i0atwt = i3atwt(i_k,i_v,j_v)
!                IF(i0atwt.EQ.1) CALL T2EXEC_AT(i_v,j_v,i_k)
             ENDDO
          ENDIF
          
!          IF(i0dtvt.EQ.1)       CALL T2EXEC_DT(i_v,j_v)
          
!          IF(i0gvvt.EQ.1)       CALL T2EXEC_GV(i_v,j_v)
          
          IF(i0gtvt.EQ.1)THEN
             DO i_k = 1, nkmax
                i0gtwt = i3gtwt(i_k,i_v,j_v)
!                IF(i0gtwt.EQ.1) CALL T2EXEC_GT(i_v,j_v,i_k)
             ENDDO
          ENDIF
          
!          IF(i0esvt.EQ.1)       CALL T2EXEC_ES(i_v,j_v)
          
          IF(i0evvt.EQ.1)THEN
             DO i_k = 1, nkmax
                i0evwt = i3evwt(i_k,i_v,j_v)
!                IF(i0evwt.EQ.1) CALL T2EXEC_EV(i_v,j_v,i_k)
             ENDDO
          ENDIF
          
          IF(i0etvt.EQ.1)THEN          
             DO j_k = 1, nkmax
             DO i_k = 1, nkmax
                i0etwt = i4etwt(i_k,j_k,i_v,j_v)
!                IF(i0etwt.EQ.1) CALL T2EXEC_ET(i_v,j_v,i_k,j_k)
             ENDDO
             ENDDO
          ENDIF
          
          !IF(i0ssvt.EQ.1)       CALL T2EXEC_SS
          
       ENDDO
       ENDDO

!       CALL T2EXEC_STORE(i_v,j_v)
       
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

    itype = 0
    m1    = 4
    
    IF(nsize.EQ.1)THEN
       m2 = 5
    ELSE
       m2 = 0
    ENDIF
    
    tolerance=1.D-7
    idebug = 0

    i0bvmax = i0bmax*i0vmax
    i0avmax = i0amax*i0vmax*i0vmax
    
       
    ALLOCATE(x(i0bvmax))
    
    CALL MTX_SETUP(i0bvmax,istart,iend,nzmax=i0avmax,idebug=0)
    
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
    !C ADDITIONAL COMPONENTS
    !C FOR FLUX SURFACE AVERAGING
    !C

    GOTO 2000

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
    
2000 CONTINUE

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
    !C SET INITIAL VALUE IN X
    !C
    
    DO i0xidi = 1, i0bmax
       DO i0vidi = 1, i0vmax
          i0xr  = i0vmax*(i0xidi-1) + i0vidi
          d0val = d2xvec_befor(i0vidi,i0xidi)
          CALL MTX_SET_VECTOR(i0xr,d0val)
       ENDDO
    ENDDO
    
    CALL MTX_SOLVE(itype,tolerance,its,methodKSP=m1,methodPC=m2)
    
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
  !C  SUBROUTINE T2EXEC_LV
  !C  　
  !C    * SET LOCAL NODE-ELEMENT GRAPH
  !C    * CALCULATE JACOBIAN OF LOCAL COORDINATES
  !C    * SET KNOWN VARIABLES FOR DIFFERENTIAL
  !C    * SET STABILIZATION FACTORS FOR SUPG (killed)
  !C
  !C                     2014-05-20 H.SETO
  !C
  !C------------------------------------------------------------------  
  SUBROUTINE T2EXEC_INITIALIZE(i_e)
    
    USE T2COMM, ONLY:&
         nnmax,i0dmax,i0vmax,i0smax,&
         d2ws,d2xvec,d2wrks,d2kwns,d2mtrc,&
         d4smat,d2svec,&
         i2enr0,i3enr,i1mlel,d0rmnr,&
         d2jmpm,d0jdmp,dt,d2jm1,d2mfc1,&
         !C
         nsmax,&       ! number of particle species
         nvmax,&       ! number of dependent variables
         grphNEElem,&  ! node-element graph 
         jacDetLocCrd,&! jacobian of local coordinates
         jacInvLocCrd,&! inverse jacobi matrix of local coordinates
         valKnownElem,&! element known variables at nodal points 
         bvecElem,&    ! element right hand side vector  (Ax=b)
         amatElem      ! element right hand side vector  (Ax=b)
    
    INTEGER(ikind),INTENT(IN)::i_e ! element number
    
    INTEGER(ikind)::i_r,i_m,i_k,i_s,i_n,j_n,k_n
    REAL(   rkind)::gridSizeRad,gridSizePol
    
    !C
    !C SET LOCAL NODE-ELEMENT GRAPH
    !C
    !C   I2ENR(1:i0nmax,1) : FOR COEFFICIENT CALCULATION
    !C   I2ENR(1:i0nmax,2) : FOR 2D-2D 
    !C   I2ENR(1:i0nmax,3) : FOR 1D-2D,1D-1D 
    !C   I2ENR(1:i0nmax,4) : FOR 2D-1D
    !C
    
    DO i_r = 1, 4
       DO i_n = 1, nnmax
          !grphNEElem(i_n,i_r)=grphNE(i_n,i_r,i_e)
          grphNEElem(i_n,i_r)=i3enr(i_n,i_r,i_e)
       ENDDO
    ENDDO
    
    
    !C
    !C CALCULATE JACOBIAN OF LOCAL COORDINATE 
    !C 
    !C jacDetLocCrd: JACOBIAN OF LOCAL COORDINATES
    !C jacInvLocCrd: INVERSE JACOBI MATRIX OF LOCAL COORDINATES
    !C
    !C modified 2014-05-20 H.SETO
    !C 
    
    i_n = grphNEElem(1,1)
    j_n = grphNEElem(2,1)
    k_n = grphNEElem(4,1)

    gridSizeRad = d2mfc1(1,j_n)-d2mfc1(1,i_n)
    gridSizePol = d2mfc1(2,k_n)-d2mfc1(2,i_n)

    gridSizeRad  = ABS(gridSizeRad)
    gridSizePol  = ABS(gridSizePol)
    
    jacDetLocCrd = gridSizeRad*gridSizePol*0.25D0

    IF(jacDetLocCrd.LE.0.D0)THEN
       WRITE(6,'("ERROR:: Jacobian of Local Coords is singular")')
       STOP
    ENDIF

    jacInvLocCrd(1,1) = 2.D0/gridSizeRad
    jacInvLocCrd(1,2) = 0.D0
    jacInvLocCrd(2,1) = 0.D0
    jacInvLocCrd(2,2) = 2.D0/gridSizePol
    
    !C
    !C STORE WORKING ARRAY FOR DIFFERENTIAL
    !C
    
    !C 
    !C valKnownElem(1   ,:) : B    AT L-TH PICARD ITERATION
    !C valKnownElem(2   ,:) : R    AT L-TH PICARD ITERATION
    !C valKnownElem(2N+1,:) : Ub   AT L-TH PICARD ITERATION 
    !C valKnownElem(2N+2,:) : Qb/P AT L-TH PICARD ITERATION 
    !C
    
    DO i_n = 1, nnmax
       i_m = grphNEElem(i_n,1)
       !valKnownElem(i_n,1) = valKnown(1,i_m)
       !valKnownElem(i_n,2) = valKnown(2,i_m)
       valKnownElem(i_n,1) = d2ws(1,i_m)
       valKnownElem(i_n,2) = d2ws(2,i_m)

       DO i_s = 1, nsmax
          i_k = 2*i_s
          !valKnownElem(i_n,i_k+1) = valKnown(i_k+1,i_m)
          !valKnownElem(i_n,i_k+2) = valKnown(i_k+2,i_m)
          valKnownElem(i_n,i_k+1) = d2ws(i_k+1,i_m)
          valKnownElem(i_n,i_k+2) = d2ws(i_k+2,i_m)
       ENDDO
    ENDDO
    
    amatElem(1:nnmax,1:nnmax,1:nvmax,1:nvmax) = 0.D0 
    bvecElem(1:nnmax,        1:nvmax        ) = 0.D0

    RETURN
    
  END SUBROUTINE T2EXEC_INITIALIZE

  !C------------------------------------------------------------------
  !C
  !C CALCULATION OF SUBMATRCES FROM MASS SCALAR COEFFICIENTS
  !C
  !C                     2014-05-20 H.SETO
  !C
  !C-------------------------------------------------------------------

  SUBROUTINE T2EXEC_MS(i_v,j_v)
    
    USE T2COMM,ONLY:&
         nnmax,ndmax,nvmax,dt,xvec,&
         massScaIntgPG,&
         massScaCoe,&
         grphNEElem,  amatElem,bvecElem,&
         jacInvLocCrd,jacDetLocCrd
    

    INTEGER(ikind), INTENT(IN)::i_v,j_v 
    INTEGER(ikind)::i_n,j_n,k_n,i_m,i_b,i_x
   
    REAL(   i0rkind)::&
         massScaMatElem(1:nnmax,1:nnmax),&
         massScaVecElem(1:nnmax        ),&
         massScaCoeElem(1:nnmax        ),&
         valsPriElem(   1:nnmax        )
    
    !C
    !C INITIALIZATION
    !C
    
    DO i_n = 1, nnmax
       i_m = grphNEElem(i_n,1)
       massScaCoeElem(       i_n        ) &
            = massScaCoe(    i_v,j_v,i_m) * jacDetLocCrd/dt
    ENDDO
    
    SELECT CASE(j_v)
       
    CASE (1:3)
       DO i_n = 1, nnmax
          i_b = grphNEElem(i_n,4)
          valsPriElem(i_n) = xvec(j_v,i_b)
       ENDDO
    CASE DEFAULT
       DO i_n = 1, nnmax
          i_x = grphNEElem(i_n,2)
          valsPriElem(i_n) = xvec(j_v,i_x)          
       ENDDO
    END SELECT
    
    !C
    !C MAIN LOOP 
    !C
    
    DO i_n = 1, nnmax
    DO j_n = 1, nnmax
       massScaMatElem(i_n,j_n) = 0.D0
       DO k_n = 1, nnmax
          massScaMatElem(           i_n,j_n) &
               = massScaMatElem(    i_n,j_n) &
               + massScaIntgPG( k_n,i_n,j_n) &
               * massScaCoeElem(k_n        ) 
       ENDDO
    ENDDO
    ENDDO

    DO i_n = 1, nnmax
       massScaVecElem(i_n) = 0.D0
       DO j_n = 1, nnmax  
          massScaVecElem(       i_n    ) &
               = massScaVecElem(i_n    ) &
               + massScaMatElem(i_n,j_n) &
               * valsPriElem(       j_n)
       ENDDO
    ENDDO
    
    !C
    !C STORE SUBMATRIX
    !C
    
    DO j_n = 1, nnmax
    DO i_n = 1, nnmax
       amatElem(             i_n,j_n,i_v,j_v) &
            = amatElem(      i_n,j_n,i_v,j_v) &
            + massScaMatElem(i_n,j_n        )
    ENDDO
    ENDDO
    
    !C
    !C STORE SUBVECTOR
    !C
    
    DO i_n = 1, nnmax
       bvecElem(             i_n,i_v) &
            = bvecElem(      i_n,i_v) & 
            + massScaVecElem(i_n    )
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
  !C                     2014-03-27 H.SETO
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
       !C SOLVE ONLY ELECTRON DENSITY 
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
                CASE(1:5,7:)
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
             CASE(1:5,7:)
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
                CASE(1:5,7:)
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
             CASE(1:5,7:)
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
                CASE(1:6,11,16,21)
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
             CASE(1:6,11,16,21)
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
                END SELECT
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
       !C i0solv = 4: 
       !C SOLVE ELECTRON 
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
                CASE(1:5,16:)
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
             CASE(1:5,16:)
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
                CASE(1:6,10:11,15:)
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
             CASE(1:6,10:11,15:)
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
                CASE(1:5)
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
             CASE(1:5)
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             CASE DEFAULT
                CYCLE
             ENDSELECT
          ENDDO
          
       ENDDO
       

       !C
       !C
       !C SET DIRICHLET CONDITION (MAGNETIC AXIS)
       !C
       !C
       
       i0bidi = 1
       
       !C
       !C STIFFNESS MATRIX
       !C
       
       DO i0aidi = i1nidr(i0bidi), i1nidr(i0bidi+1)-1
          i0bidj = i1nidc(i0aidi)
          DO i0vidi = 1, i0vmax
             SELECT CASE(i0vidi)
             CASE(1:6,8:11,13:16,18:21,23:25)
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
          CASE(1:6,8:11,13:16,18:21,23:25)
             CYCLE
          CASE DEFAULT
             d2bvec(i0vidi,i0bidi) = 0.D0
          ENDSELECT
       ENDDO
       
       !C
       !C
       !C SET DIRICHLET CONDITION (FIRST WALL)
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
                !CASE(1:6,10:11,15:16,20:21,25)
                CASE(1:5,7:8,10,12:13,15,17:18,20,22:23,25)
                !CASE(1:6,8,10:11,13,15:16,18,20:21,23,25)
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
             !CASE(1:6,10:11,15:16,20:21,25)
             !CASE(1:7,10:12,15:17,20:22,25)
             CASE(1:5,7:8,10,12:13,15,17:18,20,22:23,25)
             !CASE(1:6,8,10:11,13,15:16,18,20:21,23,25)
                CYCLE
             CASE DEFAULT
                d2bvec(i0vidi,i0bidi) = d2xvec(i0vidi,i0bidi)
             ENDSELECT
          ENDDO
       ENDDO

    CASE(10)

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
                CASE(6:)
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
             CASE(6:)
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
                CASE(6:)
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
             CASE(6:)
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
