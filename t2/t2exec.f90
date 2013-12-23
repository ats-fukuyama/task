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
  
  PUBLIC T2_EXEC
  
  PRIVATE  

CONTAINS
  
  !C
  !C MAIN ROUTINE OF T2EXEC
  !C FEM SOLVER FOR SIMULTANEOUS ADVECTION-DIFFUSION EQUATIONS
  !C 
  
  SUBROUTINE T2_EXEC
    
    USE T2COMM,ONLY:i0emax,i0eid
    
    INTEGER(i0ikind):: i1,j1,&
         i0vr,i0vc,i0vg,i0nr,i0nc,i0ng,i0tr,i0tc,i0tg
    
    REAL(4)::e0time1,e0time2
    
    CALL CPU_TIME(e0time1)
    
    DO i0eid=1,i0emax
       
       CALL T2EXEC_LV
       CALL T2EXEC_MS
       CALL T2EXEC_AV
       CALL T2EXEC_AT
       CALL T2EXEC_DT
       CALL T2EXEC_GV
       CALL T2EXEC_GT
       CALL T2EXEC_ES 
       CALL T2EXEC_EV    
       CALL T2EXEC_ET
       CALL T2EXEC_SS
       
    ENDDO
    
    CALL CPU_TIME(e0time2)
    WRITE(6,*)'TIME IN STIFF',e0time2-e0time1,'[s]'     
    CALL T2EXEC_DIRICHLET
    print*,'SOLVE_MATRIX'
    CALL T2EXEC_SOLVE
    
    RETURN
    
  END SUBROUTINE T2_EXEC
  
  !C------------------------------------------------------------------ 
  !C  SUBROUTINE T2EXEC_SET_LOCAL_VARIABLES
  !C  　
  !C    * SET LOCAL NODE-ELEMENT GRAPH
  !C    * CALCULATE JACOBIAN OF PARAMETRIC SPACE
  !C    * STORE VARIABLES ST N-th TIMESTEP   
  !C    * SET WORKING ARRAY FOR DIFFERENTIAL
  !C------------------------------------------------------------------
  
  SUBROUTINE T2EXEC_LV
    
    USE T2COMM, ONLY:&
         i0nmax0,i0dmax0,i0vmax,i0spcs,i0eid,i0lid,&
         d1rsiz,d1psiz,d2knv0,d1guv,d2ws,d2ws0,&
         i2enr0,i3enr,i1mlel,&
         d2jmpm,d0jdmp,i0nmax2,i0nmax3
    
    INTEGER::&
         i1,i3,i4,i0ln,i0lv,i0wid
    
    !C
    !C SET LOCAL NODE-ELEMENT GRAPH
    !C
    !C   I2ENR(1,:) : NON-DEGENERATED 
    !C   I2ENR(2,:) :     DEGENERATED
    !C   I2ENR(3,:) :     DEGENERATED
    !C   I2ENR(4,:) :     DEGENERATED
    !C
    
    DO i3=1,i0nmax0
       DO i1=1,4
          i2enr0(i1,i3)=i3enr(i0eid,i1,i3)
       ENDDO
    ENDDO
    
    !C
    !C CALCULATE JACOBIAN OF PARAMETRIC SPACE
    !C 
    !C D0JDMP: JACOBIAN OF PARAMETRIC SPACE
    !C D2JMPM: INVERSE JACOBI MATRIX OF PARAMETRIC SPACE
    !C
    
    i0lid = i1mlel(i0eid) 
    d0jdmp= d1rsiz(i0lid)*d1psiz(i0lid)/4.D0
    
    IF(d0jdmp.LE.0.D0)THEN
       WRITE(6,'("ERROR:: D1JDMP IS SINGULAR")')
       STOP
    ENDIF
    
    d2jmpm(1,1)= 2.D0/d1rsiz(i0lid)
    d2jmpm(1,2)= 0.D0
    d2jmpm(2,1)= 0.D0
    d2jmpm(2,2)= 2.D0/d1psiz(i0lid)
    
    !C STORE VARIABLES AT N-th TIMESTEP 

    DO i4=1,i0vmax    
       
       IF(    i4.GT.3)THEN
          DO i3=1,i0nmax0
             i0ln = i2enr0(2,i3)
             i0lv = i0vmax*(i0ln-1) + i4
             d2knv0(i3,i4) = d1guv(i0lv)
          ENDDO
       ELSEIF(i4.LE.3)THEN
          DO i3=1,i0nmax0
             i0ln = i2enr0(4,i3)
             i0lv = i0vmax*(i0ln-1) + i4
             d2knv0(i3,i4) = d1guv(i0lv)
          ENDDO
       ENDIF
    ENDDO
    
    !C SET WORKING ARRAY FOR DIFFERENTIAL

    !C 
    !C D2WS0(1   ,:) : B   AT L-TH PICARD ITERATION
    !C D2WS0(2   ,:) : lnR AT L-TH PICARD ITERATION
    !C D2WS0(5N-2,:) : Ub  AT L-TH PICARD ITERATION 
    !C D2WS0(5N-1,:) : N   AT L-TH PICARD ITERATION 
    !C D2WS0(5N  ,:) : Fb  AT L-TH PICARD ITERATION 
    !C D2WS0(5N+1,:) : P   AT L-TH PICARD ITERATION 
    !C D2WS0(5N+2,:) : Qb  AT L-TH PICARD ITERATION 
    !C
    
    DO i3 = 1, i0nmax0
       
       i0ln = i2enr0(1,i3)
       
       d2ws0(1,i3) = d2ws(1,i0ln)
       d2ws0(2,i3) = d2ws(2,i0ln)
       
       DO i4 = 1, i0spcs
          
          i0wid = 5*i0spcs - 3
          d2ws0(i0wid+1,i3) = d2ws(i0wid+1,i0ln)
          d2ws0(i0wid+2,i3) = d2ws(i0wid+2,i0ln)
          d2ws0(i0wid+3,i3) = d2ws(i0wid+3,i0ln)
          d2ws0(i0wid+4,i3) = d2ws(i0wid+4,i0ln)
          d2ws0(i0wid+5,i3) = d2ws(i0wid+5,i0ln)
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_LV
  
  SUBROUTINE T2EXEC_DIRICHLET
    
    USE T2COMM, ONLY:&
         i0lmax,i0nmax2,i0vmax,i0vgcmx,&
         i1pdn2,i1nidr,i1nidc,i1vgidr,i1vgidc,&
         d1gsm,d1grv,d1guv,i0dbg
    INTEGER::&
         i0nr,i0nc,i0ng,i0vr,i0vc,i0vg,i0tr,i0tc,i0tg,&
         i0dsta,i0dend,i0pdn2
    !C------------------------------------------------------
    
    i0pdn2 = i1pdn2(i0lmax)
    i0dsta = i0nmax2-i0pdn2
    i0dend = i0nmax2
    
    !C
    !C SET FIXED VARIALES 
    !C 

    DO i0nr= 1, i0dsta-1  
       
       !C STIFFNESS MATRIX
       
       DO i0ng = i1nidr(i0nr), i1nidr(i0nr+1)-1
          i0nc = i1nidc(i0ng)
          DO i0vr = 1, i0vmax
             !IF((i0vr.GE.6).AND.(MOD(i0vr-6,8).GE.i0dbg))THEN
             IF(i0vr.NE.3)THEN
                DO i0vg = i1vgidr(i0vr), i1vgidr(i0vr+1)-1
                   i0vc = i1vgidc(i0vg)
                   i0tr = i0vmax*( i0nr-1)+i0vr
                   i0tc = i0vmax*( i0nc-1)+i0vc
                   i0tg = i0vgcmx*(i0ng-1)+i0vg
                   IF((i0nr.EQ.i0nc).AND.(i0vr.EQ.i0vc))THEN
                      d1gsm(i0tg) = 1.D0
                   ELSE
                      d1gsm(i0tg) = 0.D0
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       
       !C RHS VECTOR 
       
       DO i0vr = 1, i0vmax
          !IF((i0vr.GE.6).AND.(MOD(i0vr-6,8).GE.i0dbg))THEN
          IF(i0vr.NE.3)THEN
             i0tr        = i0vmax*( i0nr-1)+i0vr
             d1grv(i0tr) = d1guv(i0tr)
          ENDIF
       ENDDO
    ENDDO
    
    !C SET DIRICHLET CONDITION 
    
    DO i0nr = i0dsta, i0dend
       
       !C STIFFNESS MATRIX
       
       DO i0ng = i1nidr(i0nr), i1nidr(i0nr+1)-1
          i0nc = i1nidc(i0ng)
          DO i0vr = 1, i0vmax
             !IF((i0vr.EQ.4).OR.(i0vr.EQ.5)) CYCLE
             !IF(  ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)))CYCLE
             !IF(  ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)).OR.&
             !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.3)))CYCLE
             !IF(  ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)).OR.&
             !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.3)).OR.&
             !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.6)).OR.&
             !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.7))) CYCLE
             DO i0vg = i1vgidr(i0vr),i1vgidr(i0vr+1)-1
                i0vc = i1vgidc(i0vg)
                i0tr = i0vmax*( i0nr-1)+i0vr
                i0tc = i0vmax*( i0nc-1)+i0vc
                i0tg = i0vgcmx*(i0ng-1)+i0vg
                IF((i0nr.EQ.i0nc).AND.(i0vr.EQ.i0vc))THEN
                   d1gsm(i0tg)=1.D0
                ELSE
                   d1gsm(i0tg)=0.D0
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
       !C RHS VECTOR 
       
       DO i0vr = 1, i0vmax
          
          !IF((i0vr.EQ.4).OR.(i0vr.EQ.5)) CYCLE
          !IF(  ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)))CYCLE
          !IF(  ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)).OR.&
          !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.3)))CYCLE
          !IF(  (i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)) CYCLE
          !IF(  ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.2)).OR.&
          !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.3)).OR.&
          !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.6)).OR.&
          !     ((i0vr.GE.6).AND.(MOD(i0vr-6,8).EQ.7))) CYCLE
          i0tr        = i0vmax*(i0nr-1)+i0vr
          d1grv(i0tr) = d1guv(i0tr)
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_DIRICHLET
  
  !C
  !C MODIFIED 2013-12-02
  !C
  !C
  SUBROUTINE T2EXEC_SOLVE
    
    USE T2COMM, ONLY:&
         i0vmax,i0vgcmx,i0cmax,i0bmax,i0xmax,&
         i0nmax2,i0nmax3,i0dbg,i0lmax,&
         i1nidr,i1nidc,i1vgidr,i1vgidc,i1pdn2,i1rdn2,&
         d1gsm, d1grv,d1guv_befor,d1guv_after,i2hbc,i0cmax
    
    USE LIBMPI
    USE COMMPI
    USE LIBMTX
    

    INTEGER(i0ikind)::istart,iend,its
    INTEGER(i0ikind)::itype, m1, m2,i2val(i0vmax,i0vmax)
    REAL(   i0rkind)::tolerance,d0val
    REAL(   i0rkind),POINTER,SAVE::x(:)
    REAL(4)::cputime1, cputime2 
    INTEGER(i0ikind)::&
         i1,i2,j2,i4,i5,i0cnt,&
         i0vr,i0vc,i0vg,i0nr,i0nc,i0ng,i0tr,i0tc,i0tg,&
         i0hnu,i0hnd,i0hind,i0xid,i0xidd,i0xidu,&
         i0pdn2,i0rdn2,i0ofs,&
         i0trc,i0trc1,i0trc2,i0trc3,&
         i0tcc,i0tcc1,i0tcc2,i0tcc3,&
         i0tcl,i0tcl1,i0tcl2,i0tcl3

100 FORMAT(A5,I3,A5,I3,A5,I3,A5,D15.6,A5,D15.6)
    
    !C MTX SETUP

    itype = 0
    m1    = 4
    
    IF(nsize.EQ.1)THEN
       m2 = 5
    ELSE
       m2 = 0
    ENDIF

    tolerance=1.D-7

    IF(nrank.EQ.0) CALL CPU_TIME(cputime1)
    
    ALLOCATE(x(i0bmax))
    
    CALL MTX_SETUP(i0bmax,istart,iend,nzmax=i0cmax)
    
    i0cnt=0
    
    DO i1=1,i0cmax
       d1gsm(i1)=0.d0
    ENDDO
       
    DO i1=1,i0bmax
       d1grv(i1)=0.d0
    ENDDO
       
       
    !C
    !C 
    !C
    DO i0nr=1,i0nmax2
       DO i0ng = i1nidr(i0nr), i1nidr(i0nr+1)-1
          i0nc = i1nidc(i0ng)
          DO i0vr = 1, i0vmax
             DO i0vg  = i1vgidr(i0vr), i1vgidr(i0vr+1)-1
                i0vc  = i1vgidc(i0vg)
                i0tr  = i0vmax*( i0nr - 1) + i0vr
                i0tc  = i0vmax*( i0nc - 1) + i0vc
                i0tg  = i0vgcmx*(i0ng - 1) + i0vg
                d0val = d1gsm(i0tg)
                
                IF(d0val.NE.0.D0) CALL MTX_SET_MATRIX(i0tr,i0tc,d0val)
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    !C
    !C
    !C
    
    i0ofs  = 1
    
    DO i1=1,i0lmax
       
       i0pdn2 = i1pdn2(i1)
       i0rdn2 = i1rdn2(i1)
       
       DO i2=1,i0rdn2
       DO j2=1,i0pdn2
          
          i0trc  = i0vmax*i0ofs
          i0trc1 = i0trc + 1
          i0trc2 = i0trc + 2
          i0trc3 = i0trc + 3
          
          i0tcc  = i0trc
          i0tcc1 = i0tcc + 1
          i0tcc2 = i0tcc + 2
          i0tcc3 = i0tcc + 3
          
          i0tcl  = i0trc - i0vmax
          i0tcl1 = i0tcl + 1
          i0tcl2 = i0tcl + 2
          i0tcl3 = i0tcl + 3
          
          IF((i1.EQ.i0lmax).AND.(i2.EQ.i0rdn2)) EXIT
          
          IF((j2.GE.1).AND.(j2.LT.i0pdn2))THEN
             CALL MTX_SET_MATRIX(i0trc1,i0tcc1,-1.D0)
             CALL MTX_SET_MATRIX(i0trc2,i0tcc2,-1.D0)
             CALL MTX_SET_MATRIX(i0trc3,i0tcc3,-1.D0)
          ENDIF
          IF((j2.GT.1).AND.(j2.LE.i0pdn2))THEN
             CALL MTX_SET_MATRIX(i0trc1,i0tcl1, 1.D0)
             CALL MTX_SET_MATRIX(i0trc2,i0tcl2, 1.D0)
             CALL MTX_SET_MATRIX(i0trc3,i0tcl3, 1.D0)
          ENDIF
          
          i0ofs = i0ofs + 1
          
       ENDDO
       ENDDO
    ENDDO
    
    !IF(i0ofs.NE.i0nmax2)THEN
    !   WRITE(6,*)'ERROR IN T2EXEC_SOLVE'
    !   STOP
    !ENDIF
    
    !C
    !C
    !C
    DO i4=1,i0bmax
       d0val = d1grv(i4)
       CALL MTX_SET_SOURCE(i4,d0val)
    ENDDO
    
    DO i4=1,i0bmax
       d0val = d1guv_befor(i4)
       CALL MTX_SET_VECTOR(i4,d0val)
    ENDDO
    
    CALL MTX_SOLVE(itype,tolerance,its,&
         methodKSP=m1,methodPC=m2,max_steps = 999)
    
    CALL MTX_GATHER_VECTOR(x)
    
    DO i4=1,i0nmax3
       
       IF(i4.LE.i0nmax2)THEN
          DO i5=1,i0vmax
             i0xid              = i0vmax*(i4-1)+i5
             d1guv_after(i0xid) = x(i0xid)
          ENDDO
       ELSEIF((i4.GT.i0nmax2).AND.(i4.LE.i0nmax3))THEN
          i0hind          = i4-i0nmax2
          i0hnd           = i2hbc(i0hind,1)
          i0hnu           = i2hbc(i0hind,2)
          
          DO i5=1,i0vmax 
             i0xid              = i0vmax*(i4-1)+i5
             i0xidd             = i0vmax*(i0hnd-1)+i5
             i0xidu             = i0vmax*(i0hnu-1)+i5
             d1guv_after(i0xid) = 0.5D0*(x(i0xidd)+x(i0xidu))
          ENDDO
       ENDIF
       
    ENDDO
    
    IF(nrank.eq.0)then
       CALL CPU_TIME(cputime2)
       WRITE(6,'(A,F12.3)')'--cputime',cputime2 - cputime1
    ENDIF
    
    CALL MTX_CLEANUP
    
    DEALLOCATE(x)
    
!    IF(i0tstp.eq.i0tmax) CALL MTX_FINALIZE
    
    RETURN
    
  END SUBROUTINE T2EXEC_SOLVE
  
  !C
  !C SUBROUTINES FOR CALCULATING  MATRCES 
  !C
  SUBROUTINE T2EXEC_MS
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d0jdmp,d2smat0,d1svec0,d2knv0,&
         i0msid,i1msidr,i1msidc,d2ms,d3imsn0, &
         dt
    
    INTEGER(i0ikind):: i3,j3,k3
    REAL(   i0rkind):: d1mass0(1:i0nmax0)
    
    !C INITIALIZE
    DO i0vida = 1, i0vmax
       
       !IF((i0vida.GE.6).AND.(MOD(i0vida-6,8).GT.i0dbg)) CYCLE
       
       DO i0msid = i1msidr(i0vida),i1msidr(i0vida+1)-1
          i0vidb = i1msidc(i0msid)
          
          d2smat0(1:i0nmax0,1:i0nmax0) = 0.D0
          d1svec0(1:i0nmax0          ) = 0.D0
          d1mass0(1:i0nmax0          ) = 0.D0
          
          DO i3=1,i0nmax0
             d1mass0(i3)=d2ms(i0msid,i2enr0(1,i3))
          ENDDO
          
          DO k3=1,i0nmax0
          DO j3=1,i0nmax0
          DO i3=1,i0nmax0
             d2smat0(       i3,j3   )&
                  = d2smat0(i3,j3   )&
                  + d3imsn0(i3,j3,k3)&
                  * d1mass0(      k3)&
                  * d0jdmp
          ENDDO
          ENDDO
          ENDDO 
          
          DO j3=1,i0nmax0
          DO i3=1,i0nmax0
             d2smat0(i3,j3) = d2smat0(i3,j3)/dt
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRM
    
          DO j3=1,i0nmax0
          DO i3=1,i0nmax0
             d1svec0(i3)=d1svec0(i3)+d2smat0(i3,j3)*d2knv0(j3,i0vidb)
          ENDDO
          ENDDO

          CALL T2EXEC_STRV
          
       ENDDO

    ENDDO
    RETURN
    
  END SUBROUTINE T2EXEC_MS
  
  SUBROUTINE T2EXEC_AV
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d0jdmp,d2jmpm,d2smat0,&
         i0avid,i1avidr,i1avidc,d3av,d4iavn0         
    INTEGER(i0ikind)::&
         i2,j2,i3,j3,k3
    
    REAL(   i0rkind)::&
         d2velo0(1:i0nmax0,1:i0dmax0)
    
    DO i0vida = 1, i0vmax   
       
       DO i0avid = i1avidr(i0vida),i1avidr(i0vida+1)-1
          i0vidb = i1avidc(i0avid)
          
          !C INTITIALIZE
                    
          d2velo0(1:i0nmax0,1:i0dmax0)=0.D0
          d2smat0(1:i0nmax0,1:i0nmax0)=0.D0
          
          DO i2=1,i0dmax0
             DO i3=1,i0nmax0
                d2velo0(i3,i2)=d3av(i2,i0avid,i2enr0(1,i3))
             ENDDO
          ENDDO

          DO j2=1,i0dmax0
          DO i2=1,i0dmax0
             DO k3=1,i0nmax0
             DO j3=1,i0nmax0
             DO i3=1,i0nmax0
                d2smat0(       i3,j3         )&
                     = d2smat0(i3,j3         )&
                     + d4iavn0(i3,j3,k3,   j2)&
                     * d2velo0(      k3,i2   )&
                     * d2jmpm(          i2,j2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRM
          
       ENDDO
    ENDDO
    RETURN

  END SUBROUTINE T2EXEC_AV
  
  SUBROUTINE T2EXEC_AT
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0atid,i1atidr,i2atidc,d4at,d6iatn0,d2ws0
    
    INTEGER(i0ikind)::&
         i3,j3,k3,l3,i2,j2,k2,l2,&
         i0nn,i0atsa
    REAL(   i0rkind)::&
         d3velo0(1:i0nmax0,1:i0dmax0,1:i0dmax0),&
         d1atsa0(1:i0nmax0)
    
    DO i0vida = 1, i0vmax   
       
       DO i0atid = i1atidr(i0vida),i1atidr(i0vida+1)-1
          
          i0vidb = i2atidc(i0atid,1)
          i0atsa = i2atidc(i0atid,2)
          
          !C INITIALIZATION
          
          d2smat0(1:i0nmax0,1:i0nmax0          ) = 0.D0
          d3velo0(1:i0nmax0,1:i0dmax0,1:i0dmax0) = 0.D0
          d1atsa0(1:i0nmax0                    ) = 0.D0
   
          DO i3 = 1, i0nmax0
             i0nn = i2enr0(1,i3)
             DO j2 = 1, i0dmax0
             DO i2 = 1, i0dmax0
                d3velo0(i3,i2,j2)=d4at(i2,j2,i0atid,i0nn)
             ENDDO
             ENDDO
          ENDDO
                    
          DO i3 = 1,i0nmax0
             d1atsa0(i3) = d2ws0(i0atsa,i3)
          ENDDO
          
          !C MAIN LOOP
          
          DO l2 = 1, i0dmax0
          DO k2 = 1, i0dmax0
          DO j2 = 1, i0dmax0
          DO i2 = 1, i0dmax0
             DO l3 = 1, i0nmax0
             DO k3 = 1, i0nmax0
             DO j3 = 1, i0nmax0
             DO i3 = 1, i0nmax0
                d2smat0(       i3,j3                  )&
                     = d2smat0(i3,j3                  )&
                     + d6iatn0(i3,j3,k3,l3,      k2,l2)&
                     * d1atsa0(      k3               )&   
                     * d3velo0(         l3,i2,j2      )&
                     * d2jmpm(             i2,   k2   )&
                     * d2jmpm(                j2,   l2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRM
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_AT
  
  SUBROUTINE T2EXEC_DT
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0dtid,i1dtidr,i1dtidc,d4dt,d5idtn0
    
    INTEGER(i0ikind)::&
         i2,j2,k2,l2,i3,j3,k3

    REAL(   i0rkind)::&
         d3diff0(1:i0nmax0,1:i0dmax0,1:i0dmax0)

    DO i0vida = 1, i0vmax
       
       DO i0dtid = i1dtidr(i0vida),i1dtidr(i0vida+1)-1
          
          i0vidb = i1dtidc(i0dtid)

          !C INITIALIZATION
          
          d2smat0(1:i0nmax0,1:i0nmax0          ) = 0.D0
          d3diff0(1:i0nmax0,1:i0dmax0,1:i0dmax0) = 0.D0
    
          DO j2=1,i0dmax0
          DO i2=1,i0dmax0
             DO i3=1,i0nmax0
                d3diff0(i3,i2,j2)=d4dt(i2,i2,i0dtid,i2enr0(1,i3))
             ENDDO
          ENDDO
          ENDDO
    
          !C MAIN LOOP

          DO l2=1,i0dmax0
          DO k2=1,i0dmax0
          DO j2=1,i0dmax0
          DO i2=1,i0dmax0
             DO k3=1,i0nmax0
             DO j3=1,i0nmax0
             DO i3=1,i0nmax0
                d2smat0(       i3,j3               )&
                     = d2smat0(i3,j3               )&
                     + d5idtn0(i3,j3,k3,      k2,l2)&
                     * d3diff0(      k3,i2,j2      )&
                     * d2jmpm(          i2,   k2   )&
                     * d2jmpm(             j2,   l2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO
          ENDDO
          ENDDO
    
          CALL T2EXEC_STRM

       ENDDO
       
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_DT

  SUBROUTINE T2EXEC_GV
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0gvid,i1gvidr,i1gvidc,d3gv,d4igvn0
    
    INTEGER(i0ikind)::&
         i3,j3,k3,i2,j2
    
    REAL(   i0rkind)::&
         d2grad0(1:i0nmax0,1:i0dmax0)
    
    DO i0vida = 1,i0vmax
       
       DO i0gvid = i1gvidr(i0vida),i1gvidr(i0vida+1)-1
          
          i0vidb = i1gvidc(i0gvid)
          
          !C INITIALIZATION
          
          d2grad0(1:i0nmax0,1:i0dmax0) = 0.D0
          d2smat0(1:i0nmax0,1:i0nmax0) = 0.D0
          
          DO i2=1,i0dmax0
             DO i3=1,i0nmax0
                d2grad0(i3,i2)=d3gv(i2,i0gvid,i2enr0(1,i3))
             ENDDO
          ENDDO
    
          !C MAIN LOOP
    
          DO j2=1,i0dmax0
          DO i2=1,i0dmax0
             DO k3=1,i0nmax0
             DO j3=1,i0nmax0
             DO i3=1,i0nmax0
                d2smat0(       i3,j3         )&
                     = d2smat0(i3,j3         )&
                     + d4igvn0(i3,j3,k3,   j2)&
                     * d2grad0(      k3,i2   )&
                     * d2jmpm(          i2,j2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO

          CALL T2EXEC_STRM
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GV
  
  SUBROUTINE T2EXEC_GT
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i0eid,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0gtid,i1gtidr,i2gtidc,d4gt,d6igtn0,d2ws0
    
    INTEGER(i0ikind)::&
         i3,j3,k3,l3,i2,j2,k2,l2,&
         i0gtsa,i0nn
    
    REAL(   i0rkind)::&
         d3grad0(1:i0nmax0,1:i0dmax0,1:i0dmax0),&
         d1gtsa0(1:i0nmax0)
    
    DO i0vida = 1, i0vmax
       
       DO i0gtid = i1gtidr(i0vida),i1gtidr(i0vida+1)-1
          
          i0vidb = i2gtidc(i0gtid,1)
          i0gtsa = i2gtidc(i0gtid,2)

          !C INITIALIZATION
          
          d2smat0(1:i0nmax0,1:i0nmax0          ) = 0.D0
          d3grad0(1:i0nmax0,1:i0dmax0,1:i0dmax0) = 0.D0
          d1gtsa0(1:i0nmax0                    ) = 0.D0
          
          DO i3 = 1, i0nmax0
             i0nn = i2enr0(1,i3)
             DO j2 = 1, i0dmax0
             DO i2 = 1, i0dmax0
                d3grad0(i3,i2,j2)=d4gt(i2,j2,i0gtid,i0nn)
             ENDDO
             ENDDO
          ENDDO
          
          DO i3 = 1,i0nmax0
             d1gtsa0(i3) = d2ws0(i0gtsa,i3)
          ENDDO
          
          !C MAIN LOOP
    
          DO l2 = 1, i0dmax0
          DO k2 = 1, i0dmax0
          DO j2 = 1, i0dmax0
          DO i2 = 1, i0dmax0
             DO l3 = 1, i0nmax0
             DO k3 = 1, i0nmax0
             DO j3 = 1, i0nmax0
             DO i3 = 1, i0nmax0
                d2smat0(       i3,j3                  )&
                     = d2smat0(i3,j3                  )&
                     + d6igtn0(i3,j3,k3,l3,      k2,l2)&
                     * d1gtsa0(      k3               )&   
                     * d3grad0(         l3,i2,j2      )&
                     * d2jmpm(             i2,   k2   )&
                     * d2jmpm(                j2,   l2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
             ENDDO 
          ENDDO
          ENDDO
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRM
       
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_GT
  
  SUBROUTINE T2EXEC_ES

    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0esid,i1esidr,i1esidc,d2es,d3iesn0
         
    INTEGER(i0ikind)::&
         i3,j3,k3
    REAL(   i0rkind)::&
         d1exct0(1:i0nmax0)

    DO i0vida = 1, i0vmax
       
       DO i0esid = i1esidr(i0vida),i1esidr(i0vida+1)-1
       
          i0vidb = i1esidc(i0esid)    
          
          !C INITIALIZATION
          d2smat0(1:i0nmax0,1:i0nmax0) = 0.D0
          d1exct0(1:i0nmax0)           = 0.D0
          
          DO i3=1,i0nmax0
             d1exct0(i3)=d2es(i0esid,i2enr0(1,i3))
          ENDDO

          !C MAIN LOOP
          
          DO k3=1,i0nmax0
          DO j3=1,i0nmax0
          DO i3=1,i0nmax0
             d2smat0(       i3,j3   )&
                  = d2smat0(i3,j3   )&
                  + d3iesn0(i3,j3,k3)&
                  * d1exct0(      k3)&
                  * d0jdmp
          ENDDO
          ENDDO
          ENDDO
    
          CALL T2EXEC_STRM
          
       ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE T2EXEC_ES
  
  SUBROUTINE T2EXEC_EV
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0evid,i1evidr,i2evidc,d3ev,d5ievn0,d2ws0
    
    INTEGER(i0ikind)::&
         i3,j3,k3,l3,i2,j2,&
         i0evsa
    REAL(   i0rkind)::&
         d2exct0(1:i0nmax0,1:i0dmax0),&
         d1evsa0(1:i0nmax0)
    
    DO i0vida = 1, i0vmax
       DO i0evid = i1evidr(i0vida),i1evidr(i0vida+1)-1

          i0vidb = i2evidc(i0evid,1)
          i0evsa = i2evidc(i0evid,2)
          
          !C INITIALIZATION
          d2smat0(1:i0nmax0,1:i0nmax0) = 0.D0
          d2exct0(1:i0nmax0,1:i0dmax0) = 0.D0
          d1evsa0(1:i0nmax0          ) = 0.D0
          
          DO i2 = 1, i0dmax0
             DO i3 = 1, i0nmax0
                d2exct0(i3,i2)=d3ev(i2,i0evid,i2enr0(1,i3))
             ENDDO
          ENDDO
          
          DO i3 = 1,i0nmax0
             d1evsa0(i3) = d2ws0(i0evsa,i3)
          ENDDO
          
          !C MAIN LOOP
          
          DO j2 = 1, i0dmax0
          DO i2 = 1, i0dmax0
             DO l3 = 1, i0nmax0
             DO k3 = 1, i0nmax0
             DO j3 = 1, i0nmax0
             DO i3 = 1, i0nmax0
                d2smat0(       i3,j3            )&
                     = d2smat0(i3,j3            )&
                     + d5ievn0(i3,j3,k3,l3,   j2)&
                     * d1evsa0(      k3         )&   
                     * d2exct0(         l3,i2   )&
                     * d2jmpm(             i2,j2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRM
          
       ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE T2EXEC_EV

  SUBROUTINE T2EXEC_ET
    
    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d2smat0,&
         i0etid,i1etidr,i2etidc,d4et,d7ietn0,i2etws,d2ws0
    
    INTEGER(i0ikind)::&
         i3,j3,k3,l3,m3,i2,j2,k2,l2,&
         i0etsa,i0etsb
    
    REAL(   i0rkind)::&
         d3exct0(1:i0nmax0,1:i0dmax0,1:i0dmax0),&
         d1etsa0(1:i0nmax0),&
         d1etsb0(1:i0nmax0)
    
    DO i0vida = 1, i0vmax
       
       DO i0etid = i1etidr(i0vida),i1etidr(i0vida+1)-1
          
          i0vidb = i2etidc(i0etid,1)
          i0etsa = i2etidc(i0etid,2)
          i0etsb = i2etidc(i0etid,3)
          
          !C INITIALIZATION
          d2smat0(1:i0nmax0,1:i0nmax0          ) = 0.D0
          d3exct0(1:i0nmax0,1:i0dmax0,1:i0dmax0) = 0.D0
          d1etsa0(1:i0nmax0                    ) = 0.D0
          d1etsb0(1:i0nmax0                    ) = 0.D0
    
          DO j2 = 1, i0dmax0
          DO i2 = 1, i0dmax0
             DO i3 = 1, i0nmax0
                d3exct0(i3,i2,j2)=d4et(i2,j2,i0etid,i2enr0(1,i3))
             ENDDO
          ENDDO
          ENDDO

          DO i3 = 1,i0nmax0
             d1etsa0(i3) = d2ws0(i0etsa,i3)
             d1etsb0(i3) = d2ws0(i0etsb,i3)
          ENDDO

          !C MAIN LOOP
          
          DO l2 = 1, i0dmax0
          DO k2 = 1, i0dmax0
          DO j2 = 1, i0dmax0
          DO i2 = 1, i0dmax0
          DO m3 = 1, i0nmax0
             DO l3 = 1, i0nmax0
             DO k3 = 1, i0nmax0
             DO j3 = 1, i0nmax0
             DO i3 = 1, i0nmax0
                d2smat0(       i3,j3                     )&
                     = d2smat0(i3,j3                     )&
                     + d7ietn0(i3,j3,k3,l3,m3,      k2,l2)&
                     * d1etsa0(      k3                  )&
                     * d1etsb0(         l3               )&
                     * d3exct0(            m3,i2,j2      )&
                     * d2jmpm(                i2,   k2   )&
                     * d2jmpm(                   j2,   l2)&
                     * d0jdmp
             ENDDO
             ENDDO
             ENDDO
             ENDDO
             ENDDO
          ENDDO
          ENDDO
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRM
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_ET
  
  SUBROUTINE T2EXEC_SS

    USE T2COMM,ONLY:&
         i0nmax0,i0dmax0,i0vmax,i2enr0,i0vida,i0vidb,&
         d2jmpm,d0jdmp,d1svec0,&
         i0ssid,i1ssidr,i1ssidc,d2ss,d2issn0
    
    INTEGER(i0ikind)::&
         i3,j3
    
    REAL(   i0rkind)::&
         d1sour0(1:i0nmax0)
    
    DO i0vida = 1, i0vmax
       
       DO i0ssid = i1ssidr(i0vida),i1ssidr(i0vida+1)-1
          i0vidb = i1ssidc(i0ssid)
          
          !C INITIALIZATION    
          d1svec0(1:i0nmax0) = 0.D0
          d1sour0(1:i0nmax0) = 0.D0
          
          DO i3=1,i0nmax0
             d1sour0(i3)=d2ss(i0ssid,i2enr0(1,i3))
          ENDDO
          
          DO j3=1,i0nmax0
          DO i3=1,i0nmax0
             d1svec0(       i3   )&
                  = d1svec0(i3   )&
                  + d2issn0(i3,j3)&
                  * d1sour0(   j3)&
                  * d0jdmp
          ENDDO
          ENDDO
          
          CALL T2EXEC_STRV
          
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2EXEC_SS
  
  !C
  !C SUBROUTINE FOR STORE SUBMATRIX 
  !C FOR BI-LINEAR RECTANGULAR ELEMENT
  !C
  
  SUBROUTINE T2EXEC_STRM
    
    USE T2COMM,ONLY:&
         i0nmax0,i0nmax2,i0vmax,i0vida,i0vidb,i0vgcmx,&
         i1nidr,i1nidc,i2hbc,i2vtbl,i2enr0,d2smat0,d1gsm
    
    INTEGER(i0ikind)::&
         i3,j3,i0ng,i0nc,i0tg,i0vofs,&
         i0nrl,i0nrc,i0nru,i0ncl,i0ncc,i0ncu
    REAL(   i0rkind)::&
         d0smat
    !C------------------------------------------------------
    
    i0vofs = i2vtbl(i0vida,i0vidb)
    
    
    IF(    (i0vida.GT.3).AND.(i0vidb.GT.3))THEN
       !C
       !C 2Dx2D
       !C
       DO j3 = 1, i0nmax0
       DO i3 = 1, i0nmax0
          
          i0nrc  = i2enr0(  2,i3)
          i0ncc  = i2enr0(  2,j3)
          d0smat = d2smat0(i3,j3)
          
          IF(    (i0nrc.LE.i0nmax2).AND.(i0ncc.LE.i0nmax2))THEN
             
             DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
                i0nc = i1nidc(i0ng )
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ELSEIF((i0nrc.LE.i0nmax2).AND.(i0ncc.GT.i0nmax2))THEN
             
             i0ncc  = i0ncc-i0nmax2
             i0ncl  = i2hbc(i0ncc,1)
             i0ncu  = i2hbc(i0ncc,2)
             d0smat = d0smat/2.D0
 
             DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
                i0nc = i1nidc(i0ng )
                IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ELSEIF((i0nrc.GT.i0nmax2).AND.(i0ncc.LE.i0nmax2))THEN
             
             i0nrc  = i0nrc-i0nmax2
             i0nrl  = i2hbc(i0nrc,1)
             i0nru  = i2hbc(i0nrc,2)
             d0smat = d0smat/2.D0
             
             DO i0ng = i1nidr(i0nrl), i1nidr(i0nrl+1)-1
                i0nc = i1nidc(i0ng ) 
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
             DO i0ng = i1nidr(i0nru), i1nidr(i0nru+1)-1
                i0nc = i1nidc(i0ng ) 
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ELSEIF((i0nrc.GT.i0nmax2).AND.(i0ncc.GT.i0nmax2))THEN
             
             i0nrc  = i0nrc-i0nmax2
             i0nrl  = i2hbc(i0nrc,1)
             i0nru  = i2hbc(i0nrc,2)
             
             i0ncc  = i0ncc-i0nmax2
             i0ncl  = i2hbc(i0ncc,1)
             i0ncu  = i2hbc(i0ncc,2)
             
             d0smat = d0smat/4.D0
             
             DO i0ng = i1nidr(i0nrl), i1nidr(i0nrl+1)-1
                i0nc = i1nidc(i0ng )
                IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
             DO i0ng = i1nidr(i0nru), i1nidr(i0nru+1)-1
                i0nc = i1nidc(i0ng )
                IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ENDIF
       ENDDO
       ENDDO
       
    ELSEIF((i0vida.GT.3).AND.(i0vidb.LE.3))THEN
       !C
       !C 2Dx1D
       !C
       DO j3 = 1, i0nmax0
       DO i3 = 1, i0nmax0
          
          i0nrc  = i2enr0(  2,i3)
          i0ncc  = i2enr0(  4,j3)
          d0smat = d2smat0(i3,j3)
          
          IF(    i0nrc.LE.i0nmax2)THEN
             
             DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
                i0nc = i1nidc(i0ng)
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ELSEIF(i0nrc.GT.i0nmax2)THEN
             
             i0nrc  = i0nrc-i0nmax2
             i0nrl  = i2hbc(i0nrc,1)
             i0nru  = i2hbc(i0nrc,2)
             d0smat = d0smat/2.D0
             
             DO i0ng = i1nidr(i0nrl), i1nidr(i0nrl+1)-1
                i0nc = i1nidc(i0ng) 
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
             DO i0ng = i1nidr(i0nru), i1nidr(i0nru+1)-1
                i0nc = i1nidc(i0ng) 
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ENDIF
       ENDDO
       ENDDO
       
    ELSEIF((i0vida.LE.3).AND.(i0vidb.GT.3))THEN
       
       !C
       !C 1Dx2D
       !C
       
       DO j3 = 1, i0nmax0
       DO i3 = 1, i0nmax0
          
          i0nrc  = i2enr0(  3,i3)
          i0ncc  = i2enr0(  2,j3)
          d0smat = d2smat0(i3,j3)
          
          IF(    i0ncc.LE.i0nmax2)THEN
             
             DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
                i0nc = i1nidc(i0ng)
                IF(i0nc.EQ.i0ncc)THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ELSEIF(i0ncc.GT.i0nmax2)THEN
             
             i0ncc  = i0ncc-i0nmax2
             i0ncl  = i2hbc(i0ncc,1)
             i0ncu  = i2hbc(i0ncc,2)
             d0smat = d0smat/2.D0
             
             DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
                i0nc = i1nidc(i0ng )
                IF((i0nc.EQ.i0ncl).OR.(i0nc.EQ.i0ncu))THEN
                   i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                   d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
                ENDIF
             ENDDO
             
          ENDIF
          
       ENDDO
       ENDDO
       
    ELSEIF((i0vida.LE.3).AND.(i0vidb.LE.3))THEN
       !C
       !C 1Dx1D
       !C
       DO j3 = 1, i0nmax0
       DO i3 = 1, i0nmax0
          
          i0nrc  = i2enr0(  4,i3)
          i0ncc  = i2enr0(  4,j3)
          d0smat = d2smat0(i3,j3)
          
          DO i0ng = i1nidr(i0nrc), i1nidr(i0nrc+1)-1
             i0nc = i1nidc(i0ng)
             IF(i0nc.EQ.i0ncc)THEN
                i0tg        = i0vgcmx*(i0ng-1) + i0vofs
                d1gsm(i0tg) = d1gsm(i0tg)      + d0smat
             ENDIF
          ENDDO
          
       ENDDO
       ENDDO
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2EXEC_STRM
  
  SUBROUTINE T2EXEC_STRV
    
    USE T2COMM,ONLY:&
         i0nmax0,i0nmax2,i0vmax,i0vida,i2enr0,&
         i1nidr,i1nidc,i2hbc,i2vtbl,d1grv,d1svec0
    
    REAL(   i0rkind)::&
         d0lvec
    INTEGER(i0ikind)::&
         i3,i0nrl,i0nrc,i0nru,i0trl,i0trc,i0tru
    
    IF(    i0vida.GT.3)THEN
       DO i3=1,i0nmax0
          
          i0nrc=i2enr0(2,i3)
          d0lvec = d1svec0(i3)
          
          IF(    i0nrc.LE.i0nmax2)THEN
             
             i0trc        = i0vmax*(i0nrc-1) + i0vida
             d1grv(i0trc) = d1grv(i0trc) + d0lvec
             
          ELSEIF(i0nrc.GT.i0nmax2)THEN
             
             d0lvec = d0lvec/2.D0
             i0nrc  = i0nrc - i0nmax2
             i0nrl  = i2hbc(i0nrc,1)
             i0nru  = i2hbc(i0nrc,2)
             i0trl  = i0vmax*(i0nrl-1) + i0vida
             i0tru  = i0vmax*(i0nru-1) + i0vida
             
             d1grv(i0trl)=d1grv(i0trl) + d0lvec
             d1grv(i0tru)=d1grv(i0tru) + d0lvec
             
          ENDIF
       ENDDO
    ELSEIF(i0vida.LE.3)THEN
       
       DO i3=1,i0nmax0
          
          i0nrc  = i2enr0(4,i3)
          d0lvec = d1svec0(i3)
          
          i0trc        = i0vmax*(i0nrc-1) + i0vida
          d1grv(i0trc) = d1grv(i0trc)     + d0lvec
          
       ENDDO
       
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2EXEC_STRV
END MODULE T2EXEC
