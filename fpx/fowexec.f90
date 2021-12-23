!     $Id: fpexec.f90,v 1.29 2013/02/08 07:36:24 nuga Exp $

!     **************************
!        EXECUTE TIME ADVANCE
!     **************************

MODULE fowexec

  use fpcomm
  use fowcomm

  PRIVATE
  PUBLIC fow_exec

contains

  subroutine fow_exec(NSA,IERR,its)

    USE libmpi
    USE libmtx
    USE fpmpi
    IMPLICIT NONE
    integer:: NSA, NP, NTH, NR, NL, NN, NS, imtx
    integer:: IERR,its
    integer:: imtxstart1,imtxend1
    REAL(rkind),dimension(imtxend-imtxstart+1):: BM_L

    integer :: imtx_pnc, imtx_D, imtx_Xstg, nl_lost, imtx_Ostg
    REAL(rkind) :: begin_time, end_time

    call cpu_time(begin_time)

    NS=NS_NSA(NSA)

    do nr = nrstart, nrend
      do np = npstart, npend
        do nth = 1, nthmax
          f(nth,np,nr) = fnsp(nth,np,nr,nsa)
        end do
      end do
    end do

    CALL mtx_set_communicator(comm_nrnp) !3D

    !     ----- Set up matrix solver -----
    CALL mtx_setup(imtxsize,imtxstart1,imtxend1,imtxwidth)
    IF(imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend) THEN
        WRITE(6,*) 'XX fp_exec: '
        WRITE(6,*) '   imtxstart1.NE.imtxstart.OR.imtxend1.NE.imtxend'
        WRITE(6,*) '   imtxstart1,imtxstart = ',imtxstart1,imtxstart
        WRITE(6,*) '   imtxend1,imtxend     = ',imtxend1,imtxend
        STOP
    ENDIF

    !     ----- Set up weight array -----

    CALL fowweight(NSA,IERR)

    !     ----- Set up index array IMTXA -----

    CALL SET_FM_IMTXA(NSA,FNSM)

    DO imtx=imtxstart,imtxend
        NLMAX(imtx)=0
        BM(imtx)=0.D0
        DO NL=1,NLMAXM
          LL(imtx,NL)=0
          AL(imtx,NL)=0.D0
        ENDDO
    ENDDO


    !     ----- Calculate matrix coefficients in a row -----
    do nr = nrstart, nrend
       do np = npstart, npend
          do nth = 1, nthmax
             imtx = imtxa(nth,np,nr)
             call fowsetm(nth, np, nr, nsa, nlmax(imtx))
          end do
       end do
    end do

    !     ----- internal boundary condition-----
    do nr = nrstart, nrend
       do np = npstart, npend
          
          if ( theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
             if ( pz(ns) >= 0.d0 ) then
                imtx_D = imtxa(nth_pnc(nsa)-1,np,nr)
                imtx_Xstg = imtxa(nth_stg(nsa),np,nr_pnc_point(np,nr,nsa))
             else
                imtx_D = imtxa(nth_pnc(nsa),np,nr)
                imtx_Xstg = imtxa(nth_stg(nsa)-1,np,nr_pnc_point(np,nr,nsa))
             end if
             call IBC_pinch(np,nr,nsa,nlmax(imtx_D),nlmax(imtx_Xstg))  
          end if


          if ( pz(ns) >= 0.d0 ) then
             imtx_Xstg = imtxa(nth_stg(nsa),np,nr)
             if ( nr_rhom_pinch(np,nr,nsa) <= nrmax ) then
                imtx_pnc = imtxa(nth_pnc(nsa)-1,np,nr_rhom_pinch(np,nr,nsa))
                call IBC_X_stagnation(np,nr,nsa, &
                     nlmax(imtx_Xstg),nlmax(imtx_pnc))
             else
                nl_lost = -1
                call IBC_X_stagnation(np,nr,nsa,nlmax(imtx_Xstg),nl_lost)
             end if

          else
             imtx_Xstg = imtxa(nth_stg(nsa)-1,np,nr)
             if ( nr_rhom_pinch(np,nr,nsa) <= nrmax ) then
                imtx_pnc = imtxa(nth_pnc(nsa),np,nr_rhom_pinch(np,nr,nsa))
                call IBC_X_stagnation(np,nr,nsa, &
                     nlmax(imtx_Xstg),nlmax(imtx_pnc))
             else
                nl_lost = -1
                call IBC_X_stagnation(np,nr,nsa,nlmax(imtx_Xstg),nl_lost)
             end if

          end if

          if ( pz(ns) >= 0.d0 ) then
             imtx_Ostg = imtxa(nth_pnc(nsa)-1,np,nr)
          else
             imtx_Ostg = imtxa(nth_pnc(nsa),np,nr)
          end if
          call IBC_O_stagnation(np,nr,nsa,nlmax(imtx_Ostg))

       end do
    end do


    !     ----- Diagonal term -----
    
    do nr = nrstart, nrend
       do np = npstart, npend
          do nth = 1, nthmax
             imtx = imtxa(nth,np,nr)
             bm(imtx) = (1.d0+(1.d0-rimpl)*delt*dl(imtx))*fm(imtx) &
                      +delt*spp(nth,np,nr,nsa)
             if ( imtx.ge.imtxstart.and.imtx.le.imtxend ) then
                call mtx_set_matrix(imtx, imtx, 1.d0-rimpl*delt*dl(imtx))
                call mtx_set_vector(imtx, fm(imtx))
             end if
          end do
       end do
    end do

    !     ----- Off diagonal term -----

    DO imtx=imtxstart,imtxend ! LHS
       DO NL=1,NLMAX(imtx)
          IF(LL(imtx,NL).NE.0) THEN
             CALL mtx_set_matrix(imtx,LL(IMTX,NL),-RIMPL*DELT*AL(imtx,NL))
          ENDIF
       ENDDO
    ENDDO
    !     ----- Source vector: contribution from off-diagonal term -----

    DO imtx=imtxstart,imtxend ! RHS
       DO NL=1,NLMAX(imtx)
          NN=LL(imtx,NL)
          IF(NN.NE.0) THEN
             IF(NN.ge.imtxstart-NTHMAX.and.NN.le.imtxend+NTHMAX)THEN
                BM(imtx)=BM(imtx)+(1.D0-RIMPL)*DELT*AL(imtx,NL)*FM(NN)
             ELSEIF(NN.lt.imtxstart-NTHMAX)THEN
                BM(imtx)=BM(imtx)+(1.D0-RIMPL)*DELT*AL(imtx,NL)*FM_shadow_m(NN)
             ELSE
                BM(imtx)=BM(imtx)+(1.D0-RIMPL)*DELT*AL(imtx,NL)*FM_shadow_p(NN)
             END IF
          ENDIF
       ENDDO
       CALL mtx_set_source(imtx,BM(imtx))
    ENDDO


    !     ----- Solve matrix equation -----

    CALL mtx_solve(imtx,epsm,its,MODEL_KSP,MODEL_PC)
    ierr=0

    !     ----- Get solution vector -----

    CALL mtx_get_vector(BM_L)

    DO nr=nrstart, nrend
       DO np=npstart, npend
          DO nth=1,nthmax
             imtx=imtxa(nth,np,nr)
             fns0(nth,np,nr,nsa)=bm_l(imtx-imtxstart+1)
          END DO
       END DO
    END DO
    
    !     shadow requires to communicate
    CALL mtx_set_communicator(comm_np)
    DO nr=nrstart, nrend
        CALL shadow_comm_np(nr,nsa)
    END DO
    CALL mtx_set_communicator(comm_nr)
    CALL shadow_comm_nr(nsa)
    CALL mtx_set_communicator(comm_nrnp) !3D


    !     ----- Clean up matrix solver -----
    CALL mtx_cleanup

    CALL mtx_reset_communicator

    call cpu_time(end_time)
    write(6,'(A,I0,A,ES10.3,A)') 'fowexec time (nsa=',nsa,') : ', &
         end_time-begin_time,'[sec]'

    RETURN
  END subroutine fow_exec

  !     ---------------------------------

  SUBROUTINE SET_FM_IMTXA(NSA,func_in)

    IMPLICIT NONE
    integer:: NTH, NP, NR, NSA, imtx, NRS
    REAL(rkind),intent(IN):: &
         func_in(NTHMAX,NPSTARTW:NPENDWM,NRSTARTW:NRENDWM,NSASTART:NSAEND)

    IF(NRSTART.eq.1)THEN
        NRS=1
    ELSE
        NRS=NRSTART-1
    END IF

    DO NR=NRSTARTW,NRENDWM
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              imtx=NTH+NTHMAX*(NP-1)+NPMAX*NTHMAX*(NR-1)
              IMTXA(NTH,NP,NR)=imtx
          END DO
        END DO
    END DO

    DO NR=NRSTART,NREND
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              IMTX=IMTXA(NTH,NP,NR)
              FM(IMTX)=func_in(NTH,NP,NR,NSA)
          ENDDO
        ENDDO
    ENDDO
    NR=NRSTARTW
    IF(NR.ne.NRSTART)THEN
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              IMTX=IMTXA(NTH,NP,NR)
              FM_shadow_m(IMTX)=func_in(NTH,NP,NR,NSA)
          ENDDO
        ENDDO
    END IF
    NR=NRENDWM
    IF(NR.ne.NREND)THEN
        DO NP=NPSTARTW,NPENDWM
          DO NTH=1,NTHMAX
              IMTX=IMTXA(NTH,NP,NR)
              FM_shadow_p(IMTX)=func_in(NTH,NP,NR,NSA)
          ENDDO
        ENDDO
    END IF


  END SUBROUTINE SET_FM_IMTXA

  !
  !     ***************************
  !        CALCULATION OF WEIGHT
  !     ***************************
  !
  subroutine fowweight(NSA,IERR) ! proposed by Chang and Cooper [30] in Karney

    use fpcomm
    use fowcomm

    implicit none
    integer:: nsa, np, nth, nr, nthl, nthr, ns
    integer:: ierr
    REAL(rkind):: dfdth, fvel, dfdp
    
    integer :: nrl, nrr, npl, npr
    REAL(rkind) :: dfdrm, width_p, width_r, width_t, dth

    !     +++++ calculation of weigthing (including off-diagonal terms) +++++

    REAL(rkind)::epswt=1.d-70


    ierr=0
    ns = ns_nsa(nsa)

    do nr = nrstart, nrend
      do np = npstart, npendwg
        do nth = 1, nthmax
            dfdth = 0.d0
            dfdrm  = 0.d0
            if ( np /= 1 ) then
              nthl = min(nth+1,nthmax)
              nthr = max(nth-1,1)
              nrl  = min(nr+1,nrmax)
              nrr  = max(nr-1,1)

              if ( np == npmax+1 ) then
                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                       /(2.d0*pg(np,ns)*delthm(nth,np-1,nr,nsa) &
                       *f(nth,np-1,nr))

                  dfdrm  = (f(nth,np-1,nrl)-f(nth,np-1,nrr)) &
                       /(2.d0*delr*f(nth,np-1,nr))
                end if
              else

                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  if ( abs(f(nth,np,nr)) > epswt ) then
                  
                    dth =(delthm(nth,np,nr,nsa)+delthm(nth,np-1,nr,nsa))*0.5d0
                    if ( dth > 0.d0 ) then
                      dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                              /(4.d0*pg(np,ns)*dth*f(nth,np-1,nr)) &
                            + (f(nthl,np  ,nr)-f(nthr,np  ,nr)) &
                              /(4.d0*pg(np,ns)*dth*f(nth,np  ,nr))
                    end if

                    dfdrm  = (f(nth,np-1,nrl)-f(nthr,np-1,nrr)) &
                            /(4.d0*delr*f(nth,np-1,nr)) &
                          + (f(nthl,np  ,nrl)-f(nthr,np  ,nrr)) &
                            /(4.d0*delr*f(nth,np  ,nr))

                  else
                    dth=(delthm(nth,np,nr,nsa)+delthm(nth,np-1,nr,nsa))*0.5d0
                    if ( dth > 0.d0 ) then
                      dfdth = (f(nthl,np-1,nr)-f(nthr,np-1,nr)) &
                           /(2.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)* &
                           f(nth,np-1,nr)) 
                    end if
                    
                    dfdrm  = (f(nth,np-1,nrl)-f(nth,np-1,nrr)) &
                         /(2.d0*delr*f(nth,np-1,nr)) 
                  end if
                else
                  if ( abs(f(nth,np,nr)) > epswt ) then
                    dth=(delthm(nth,np,nr,nsa)+delthm(nth,np-1,nr,nsa))*0.5d0
                    if ( dth >= 0.d0 ) then
                      dfdth = (f(nthl,np,nr)-f(nthr,np,nr)) &
                            /(2.d0*pg(np,ns)*dth*f(nth,np  ,nr))
                    end if

                    dfdrm  = (f(nth,np  ,nrl)-f(nth,np  ,nrr)) &
                          /(2.d0*delr*f(nth,np  ,nr))
                  end if
                end if
              end if
            end if
            fvel = Fppfow(nth,np,nr,nsa) &
                 -Dptfow(nth,np,nr,nsa)*dfdth &
                 -Dprfow(nth,np,nr,nsa)*dfdrm
            weighp(nth,np,nr,nsa)=fowwegh(-delp(ns)*fvel,dppfow(nth,np,nr,nsa))

        end do
      end do
    end do

    do nr = nrstart, nrend
      do np = npstart, npend
        do nth = 1, nthmax+1
          dfdp = 0.d0          
          npl = min(npmax,np+1)
          npr = max(1, np-1)
          width_p = dble(npl-npr)

          dfdrm = 0.d0
          nrl = min(nrmax,nr+1)
          nrr = max(1,nr-1)
          width_r = dble(nrl-nrr)
          
          if ( nth == 1 ) then
            if ( abs(f(nth,np,nr)) > epswt ) then
               dfdp = (f(nth,npl,nr)-f(nth,npr,nr)) &
                    /(width_p*delp(ns)*f(nth,np,nr))
               dfdrm = (f(nth,np,nrl)-f(nth,np,nrr)) &
                    /(width_r*delr   *f(nth,np,nr))  
            end if
          else if ( nth == nthmax+1 ) then
            if ( abs(f(nthmax,np,nr)) > epswt ) then
               dfdp = (f(nthmax,npl,nr)-f(nthmax,npr,nr)) &
                    /(width_p*delp(ns)*f(nthmax,np,nr))
               dfdrm = (f(nthmax,np,nrl)-f(nthmax,np,nrr)) &
                    /(width_r*delr*f(nthmax,np,nr))  
            end if
          else
            if(abs(f(nth,np,nr))>epswt.and.abs(f(nth-1,np,nr))>epswt) then
              dfdp = ( &
                   (f(nth-1,npl,nr)-f(nth-1,npr,nr)) &
                   /(width_p*delp(ns)*f(nth-1,np,nr)) &
                  +(f(nth  ,npl,nr)-f(nth  ,npr,nr)) &
                   /(width_p*delp(ns)*f(nth  ,np,nr)) &
                   )/2.d0

              dfdrm = ( &
                   (f(nth-1,np,nrl)-f(nth-1,np,nrr)) &
                   /(width_r*delr*f(nth-1,np,nr))+&
                   (f(nth  ,np,nrl)-f(nth  ,np,nrr)) &
                   /(width_r*delr*f(nth  ,np,nr)) &
              )
           else if( abs(f(nth,np,nr)) > epswt .and. &
                    abs(f(nth-1,np,nr)) <= epswt ) then
              dfdp = (f(nth  ,npl,nr)-f(nth  ,npr,nr)) &
                   /(width_p*delp(ns)*f(nth  ,np,nr))
              dfdrm = (f(nth  ,np,nrl)-f(nth  ,np,nrr)) &
                   /(width_r*delr*f(nth  ,np,nr))
           else if ( abs(f(nth,np,nr)) <= epswt .and. &
                     abs(f(nth-1,np,nr)) > epswt ) then
              dfdp = (f(nth-1,npl,nr)-f(nth-1,npr,nr)) &
                   /(width_p*delp(ns)*f(nth-1,np,nr))
              dfdrm = (f(nth-1,np,nrl)-f(nth-1,np,nrr)) &
                   /(width_r*delr*f(nth-1,np,nr))
            end if
          end if
          fvel = Fthfow(nth,np,nr,nsa) &
               -Dtpfow(nth,np,nr,nsa)*dfdp &
               -Dtrfow(nth,np,nr,nsa)*dfdrm
          if ( nth <= nthmax ) then
             weight(nth,np,nr,nsa) &
                  =fowwegh(-delthm(nth,np,nr,nsa)*fvel,Dttfow(nth,np,nr,nsa))
          else
             weight(nth,np,nr,nsa) &
                  =fowwegh(-delthm(nth-1,np,nr,nsa)*fvel,Dttfow(nth,np,nr,nsa))
          end if
        end do
      end do
    end do

    do nr = nrstart, nrendwg
      do np = npstart, npend
        do nth = 1, nthmax
          dfdp = 0.d0          
          npl = min(npmax,np+1)
          npr = max(1, np-1)
          width_p = dble(npl-npr)

          dfdth = 0.d0
          nthl = min(nth+1,nthmax)
          nthr = max(nth-1,1)
          width_t = dble(nthl-nthr)

          if ( nr == 1 ) then
            if ( abs(f(nth,np,1)) > epswt ) then
               dfdp = (f(nth,npl,1)-f(nth,npr,1)) &
                    /(width_p*delp(ns)*f(nth,np,1))
              if ( delthm_rg(nth,np,nr,nsa) /= 0.d0 ) then
                 dfdth = (f(nthl,np,1)-f(nthr,np,1)) &
                      /(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,1))
              else
                 dfdth = 0.d0
              end if
            end if
          else if ( nr == nrmax+1 ) then
             if ( abs(f(nth,np,nrmax)) > epswt ) then
                dfdp = (f(nth,npl,nrmax)-f(nth,npr,nrmax)) &
                     /(width_p*delp(ns)*f(nth,np,nrmax))
                if ( delthm_rg(nth,np,nr,nsa) /= 0.d0 ) then
                   dfdth= (f(nthl,np,nrmax)-f(nthr,np,nrmax)) &
                        /(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nrmax))
                end if
             end if
          else
             if ( abs(f(nth,np,nr)) > epswt .and. &
                  abs(f(nth,np,nr-1)) > epswt ) then
              dfdp = ( &
                   (f(nth,npl,nr-1)-f(nth,npr,nr-1)) &
                   /(width_p*delp(ns)*f(nth,np,nr-1)) &
                  +(f(nth,npl,nr  )-f(nth,npr,nr  )) &
                   /(width_p*delp(ns)*f(nth,np,nr  )) &
                   )/2.d0

              if ( nr == 2 ) then
                 if ( delthm_rg(nth,np,nr,nsa) > 0.d0 ) then
                    dfdth= (f(nthl,np,nr)-f(nthr,np,nr)) &
                         /(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nr))
                 end if
              else
                 if ( delthm_rg(nth,np,nr-1,nsa) > 0.d0 .and. &
                      delthm_rg(nth,np,nr  ,nsa) > 0.d0 ) then
                    if ( delthm_rg(nth,np,nr-1,nsa) > 0.d0 .and. &
                         delthm_rg(nth,np,nr  ,nsa) > 0.d0 ) then
                       dfdth=( &
                            (f(nthl,np,nr-1)-f(nthr,np,nr-1)) &
                            /(width_t*delthm_rg(nth,np,nr-1,nsa)*f(nth,np,nr))&
                           +(f(nthl,np,nr  )-f(nthr,np,nr  )) &
                           /(width_t*delthm_rg(nth,np,nr  ,nsa)*f(nth,np,nr))&
                           )/2.d0                  
                    end if
                 end if
              end if

           else if ( abs(f(nth,np,nr)) > epswt .and. &
                     abs(f(nth,np,nr-1)) <= epswt ) then
              dfdp = (f(nth,npl,nr  )-f(nth,npr,nr  )) &
                   /(width_p*delp(ns)*f(nth,np,nr))
              if ( delthm_rg(nth,np,nr,nsa) > 0.d0 ) then
                 dfdth= (f(nthl,np,nr  )-f(nthr,np,nr  )) &
                      /(width_t*delthm_rg(nth,np,nr,nsa)*f(nth,np,nr))
              end if
              
           else if ( abs(f(nth,np,nr)) <= epswt .and. &
                abs(f(nth,np,nr-1)) > epswt ) then
              dfdp = (f(nth,npl,nr-1)-f(nth,npr,nr-1)) &
                   /(width_p*delp(ns)*f(nth,np,nr-1))
              if ( delthm_rg(nth,np,nr-1,nsa) > 0.d0  ) then
                 dfdth= (f(nthl,np,nr-1)-f(nthr,np,nr-1)) &
                      /(width_t*delthm_rg(nth,np,nr-1,nsa)*f(nth,np,nr-1))
              end if

            end if
          end if
          fvel = Frrfow(nth,np,nr,nsa) &
               -Drpfow(nth,np,nr,nsa)*dfdp &
               -Drtfow(nth,np,nr,nsa)*dfdth
          weighr(nth,np,nr,nsa) = fowwegh(-delr*fvel,Drrfow(nth,np,nr,nsa))
        end do
      end do
    end do

    return
  end subroutine fowweight

  ! ************************************************
  !     WEIGHTING FUNCTION FOR CONVECTION EFFECT
  ! ************************************************

  FUNCTION fowwegh(X,Y)

    IMPLICIT NONE
    REAL(rkind):: X, Y, Z
    REAL(rkind):: fowwegh

    IF(ABS(Y).LT.1.D-70) THEN
        IF(X.GT.0.D0) THEN
          fowwegh=0.D0
        ELSEIF(X.LT.0.D0) THEN
          fowwegh=1.D0
        ELSE
          fowwegh=0.5D0
        end if
    ELSE
        Z=X/Y
        IF(ABS(Z).LT.1.D-5)THEN
          fowwegh=0.5D0-Z/12.D0+Z**3/720.D0
        ELSE IF(Z.GE.100.D0)THEN
          fowwegh=1.D0/Z
        ELSE IF(Z.LE.-100.D0)THEN
          fowwegh=1.D0/Z+1.D0
        ELSE
          fowwegh=1.D0/Z-1.D0/(EXP(Z)-1.D0)
        END IF
    ENDIF
    RETURN
  END FUNCTION fowwegh

  ! ******************************************
  !     Calculation of matrix coefficients
  ! ******************************************

  subroutine fowsetm(NTH,NP,NR,NSA,NL)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NTH,NP,NR,NSA
    INTEGER,INTENT(OUT):: NL
    INTEGER:: IMTX

    integer:: ns
    REAL(rkind):: pl, sl, rl
    REAL(rkind) :: deltath
    REAL(rkind) :: DIVD(3,3), DIVF(3), Ffow(2,3)
    REAL(rkind) :: D_term, F_term
    integer :: si, sj, sign_to_index(-1:1), alpha, beta, loc(4)
    integer :: boundary_flag


    NS=NS_NSA(NSA)

    NL=0
    IMTX=IMTXA(NTH,NP,NR)

    Ffow(1,1) = Fppfow(nth,np,nr,nsa)*pg(np,nsa)**2
    Ffow(2,1) = Fppfow(nth,np+1,nr,nsa)*pg(np+1,nsa)**2

    Ffow(1,2) = Fthfow(nth,np,nr,nsa)*sin(thetamg(nth,np,nr,nsa))
    Ffow(2,2) = Fthfow(nth+1,np,nr,nsa)*sin(thetamg(nth+1,np,nr,nsa))

    Ffow(1,3) = Frrfow(nth,np,nr,nsa)*rg(nr)
    Ffow(2,3) = Frrfow(nth,np,nr+1,nsa)*rg(nr+1)


    ! discretized (div(d/dX))_Y
    pl = pm(np,nsa)
    sl = SIN( thetam(nth,np,nr,nsa) )
    rl = rm(nr)
    deltath= delthm(nth,np,nr,nsa)

    DIVD(1,1) = 1.d0/(pl**2 * delp(nsa)**2)
    DIVD(1,2) = 1.d0/(pl**2 * delp(nsa) * deltath * 2.d0)
    DIVD(1,3) = 1.d0/(pl**2 * delp(nsa) * delr * 2.d0)
    DIVD(2,1) = 1.d0/(pl * sl * delp(nsa) * deltath * 2.d0)
    DIVD(2,2) = 1.d0/(pl * sl * deltath**2 )
    DIVD(2,3) = 1.d0/(pl * sl * deltath * delr * 2.d0)
    DIVD(3,1) = 1.d0/(rl * delp(nsa)  * delr * 2.d0)
    DIVD(3,2) = 1.d0/(rl * deltath * delr * 2.d0)
    DIVD(3,3) = 1.d0/(rl * delr**2)
    DIVF(1)   = 1.d0/(pl**2 * delp(nsa))
    DIVF(2)   = 1.d0/(pl * sl * deltath)
    DIVF(3)   = 1.d0/(rl * delr)

    sign_to_index(-1) = 1
    sign_to_index(1)  = 2
    sign_to_index(0)  = 0
    loc = [nth,np,nr,nsa]

    ! term of f(nth, np, nr)
    DL(IMTX) = 0.d0
    do alpha = 1, 3
       do si = -1, 1, 2
          D_term = -1.d0*Dfow(alpha,alpha,sign_to_index(si)-1,0,loc) &
               *DIVD(alpha,alpha)
          F_term = -1.d0*si*Ffow(sign_to_index(si),alpha) &
               *w(si,alpha,0,sign_to_index(si)-1,0,loc)*DIVF(alpha)
          dl(imtx) = dl(imtx) + (D_term + F_term)/JI(nth,np,nr,nsa)
       end do
    end do
    if ( nth == nth_pnc(nsa) .and. &
         theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
       D_term = Dfow(2,2,0,0,loc)*DIVD(2,2)*IBCflux_ratio(np,nr,nsa)
       F_term = -1.d0*Ffow(1,2)*w(-1,2,0,0,0,loc)*DIVF(2) &
            *IBCflux_ratio(np,nr,nsa)
       dl(imtx) = dl(imtx) + (D_term + F_term)/JI(nth,np,nr,nsa)
    end if

    ! terms of f(i+si, j, k) 
    do alpha = 1, 3
       do si = -1, 1, 2
          boundary_flag = check_external_boundary(alpha,0,si,0,loc)
          if ( boundary_flag == 1  ) cycle

          D_term = Dfow(alpha,alpha,sign_to_index(si)-1,0,loc) &
               *DIVD(alpha,alpha)
          F_term = -1.d0*si*Ffow(sign_to_index(si),alpha) &
               *w(-1*si,alpha,0,sign_to_index(si)-1,0,loc)*DIVF(alpha)

          do beta = 1, 3
             if ( alpha == beta ) cycle
             do sj = -1, 1, 2
                D_term = D_term &
                     + si * sj*Dfow(beta,alpha,sign_to_index(sj)-1,0,loc) &
                     *w(sj,beta,alpha,sign_to_index(sj)-1,si,loc) &
                     *DIVD(beta,alpha)
             end do
          end do

          nl = nl+1
          ll(imtx,nl) = get_imtxa(alpha,0,si,0,loc)
          al(imtx,nl) = al(imtx,nl) + (D_term + F_term)/JI(nth,np,nr,nsa)

          if ( ABS(al(imtx,nl)) < 1.d-70 ) then
             ll(imtx,nl) = 0
             al(imtx,nl) = 0.d0
             nl = nl-1
          end if

       end do
    end do

    ! terms of f(i+si, j+sj, k) 
    do alpha = 1, 3
       do beta = 1, 3
          if ( alpha >= beta) cycle
          do si = -1, 1, 2
             do sj = -1, 1, 2

                boundary_flag = check_external_boundary(alpha,beta,si,sj,loc)
                if ( boundary_flag == 1 ) cycle
            
                D_term = si*sj*( &
                     Dfow(alpha,beta,sign_to_index(si)-1,0,loc) &
                     *DIVD(alpha,beta) &
                     *w(-1*si,alpha,beta,sign_to_index(si)-1,sj,loc) &
                    +Dfow(beta,alpha,sign_to_index(sj)-1,0,loc) &
                     *DIVD(beta,alpha) &
                     *w(-1*sj,beta,alpha,sign_to_index(sj)-1,si,loc) &
                     )

                nl = nl+1
                ll(imtx,nl) = get_imtxa(alpha,beta,si,sj,loc)
                al(imtx,nl) = al(imtx,nl) + D_term/JI(nth,np,nr,nsa)

                if ( ABS(al(imtx,nl)) < 1.d-70 ) then
                   ll(imtx,nl) = 0
                   al(imtx,nl) = 0.d0
                   nl = nl-1
                end if

             end do
          end do
       end do
    end do

    spp(nth,np,nr,nsa) = sppb(nth,np,nr,nsa)

    RETURN
  END subroutine fowsetm

  function Dfow(alpha,beta,si,sj,loc)
    implicit none
    REAL(rkind) :: Dfow
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    if ( alpha == 1 ) then
       if ( beta == 1 ) then
          Dfow = Dppfow(nth,np+si,nr,nsa)*pg(np+si,nsa)**2
       else if ( beta == 2 ) then
          Dfow = Dptfow(nth+sj,np+si,nr,nsa)*pg(np+si,nsa)
       else if ( beta == 3 ) then
          Dfow = Dprfow(nth,np+si,nr+sj,nsa)*pg(np+si,nsa)**2
       end if

    else if ( alpha == 2 ) then
       if ( beta == 1 ) then
          Dfow = Dtpfow(nth+si,np+sj,nr,nsa)*sin(thetamg(nth+si,np+sj,nr,nsa))
       else if ( beta == 2 ) then
          Dfow = Dttfow(nth+si,np,nr,nsa)*sin(thetamg(nth+si,np+sj,nr,nsa)) &
               /pm(np,nsa)
       else if ( beta == 3 ) then
          Dfow = Dtrfow(nth+si,np,nr+sj,nsa)*sin(thetamg(nth+si,np,nr+sj,nsa))
       end if

    else if ( alpha == 3 ) then
       if ( beta == 1 ) then
          Dfow = Drpfow(nth,np+sj,nr+si,nsa)*rg(nr+si)
       else if ( beta == 2 ) then
          Dfow = Drtfow(nth+sj,np,nr+si,nsa)*rg(nr+si)/pm(np,nsa)
       else if ( beta == 3 ) then
          Dfow = Drrfow(nth,np,nr+si,nsa)*rg(nr+si)
       end if

    end if

  end function Dfow

  function w(sign,alpha,beta,si,sj,loc)
    ! weight function for alpha, beta, gamma
    implicit none
    REAL(rkind) :: w
    integer,intent(in) :: sign,alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)
    
    if ( alpha == 1 ) then
       if ( beta == 2 ) then
          w = WEIGHP(nth+sj, np+si, nr, nsa)
       else if ( beta == 3 ) then
          w = WEIGHP(nth, np+si, nr+sj, nsa)
       else
          w = WEIGHP(nth, np+si, nr, nsa)
       end if

    else if ( alpha == 2 ) then
       if ( beta == 1 ) then
          w = WEIGHT(nth+si, np+sj, nr, nsa)
       else if ( beta == 3 ) then
          w = WEIGHT(nth+si, np, nr+sj, nsa)
      else 
         w = WEIGHT(nth+si, np, nr, nsa)
      end if

    else
       if ( beta == 1 ) then
          w = WEIGHR(nth, np+sj, nr+si, nsa)
       else if ( beta == 2 ) then
          w = WEIGHR(nth+sj, np, nr+si, nsa)
       else
          w = WEIGHR(nth, np, nr+si, nsa)
       end if
    end if

    if ( sign < 0 ) w = 1.d0-w

  end function w

  function check_external_boundary(alpha,beta,si,sj,loc) result(flag)
    ! flag == 0 -> inside external boundary
    !         1 -> outside external boundary
    implicit none
    integer :: flag
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa
    
    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    flag = 0

    ! check external boundary
    if ( alpha == 1 ) then
      if ( beta == 2 ) then
        if ( (np+si <= 0 .or. npmax+1 <= np+si)&
            .or. (nth+sj <= 0 .or. nthmax+1 <= nth+sj) )&
            flag = 1

      else if ( beta == 3 ) then
        if ( (np+si <= 0 .or. npmax+1 <= np+si)&
            .or. (nr+sj <= 0 .or. nrmax+1 <= nr+sj) )&
          flag = 1

      else ! beta == 0
        if ( np+si <= 0 .or. npmax+1 <= np+si ) flag = 1

      end if

    else if ( alpha == 2 ) then
      if ( beta == 1 ) then
        if ( (np+sj <= 0 .or. npmax+1 <= np+sj)&
            .or. (nth+si <= 0 .or. nthmax+1 <= nth+si) )&
          flag = 1

      else if ( beta == 3 ) then
        if ( (nth+si <= 0 .or. nthmax+1 <= nth+si)&
            .or. (nr+sj <= 0 .or. nrmax+1 <= nr+sj) )&
          flag = 1

      else ! beta == 0
        if ( nth+si <= 0 .or. nthmax+1 <= nth+si ) flag = 1

      end if

    else
      if ( beta == 1 ) then
        if ( (np+sj <= 0 .or. npmax+1 <= np+sj)&
            .or. (nr+si <= 0 .or. nrmax+1 <= nr+si) )&
          flag = 1

      else if ( beta == 2 ) then
          if ( (nth+sj <= 0 .or. nthmax+1 <= nth+sj)&
            .or. (nr+si <= 0 .or. nrmax+1 <= nr+si) )&
          flag = 1

      else ! beta == 0
        if ( nr+si <= 0 .or. nrmax+1 <= nr+si ) flag = 1

      end if
    end if

    ! check internal boundary
    

  end function check_external_boundary

  function get_imtxa(alpha,beta,si,sj,loc) result(n)
    implicit none
    integer :: n
    integer,intent(in) :: alpha,beta,si,sj,loc(4)
    integer :: nth, np, nr, nsa

    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    SELECT CASE(alpha)
    CASE(1)
       SELECT CASE (beta)
       CASE(2)
          n = imtxa(nth+sj,np+si,nr)
       CASE(3)
          n = imtxa(nth,np+si,nr+sj)
       CASE DEFAULT ! beta == 0
          n = imtxa(nth,np+si,nr)
       END SELECT
    
    CASE(2)
       SELECT CASE(beta)
       CASE(1)
          n = imtxa(nth+si,np+sj,nr)
       CASE(3)
        n = imtxa(nth+si,np,nr+sj)
       CASE DEFAULT ! beta == 0
          n = imtxa(nth+si,np,nr)
       END SELECT
       
    CASE(3)
       SELECT CASE(beta)
       CASE(1)
          n = imtxa(nth,np+sj,nr+si)
       CASE(2)
          n = imtxa(nth+sj,np,nr+si)
       CASE DEFAULT ! beta == 0
          n = imtxa(nth,np,nr+si)
       END SELECT

    END SELECT

  end function get_imtxa

  subroutine IBC_pinch(np,nr,nsa,nl_D,nl_Xstg)
    implicit none
    integer,intent(in) :: np,nr,nsa
    integer,intent(inout) :: nl_D, nl_Xstg
    REAL(rkind) :: S_pnc_to_trp, S_pnc_to_pas

    integer:: ns, nth, imtx_D, imtx_Xstg, nl
    REAL(rkind):: pl, sl, rl
    REAL(rkind) :: deltath, fact1, fact2
    REAL(rkind) :: Fthp, Dtx(3)
    REAL(rkind) :: dfnspdp, dfnspdth, dfnspdr, fnspg
    integer :: loc_pnc(4), loc_Xstg(4), loc_D(4)

    if ( theta_pnc(np,nr,nsa) == NO_PINCH_ORBIT ) then
      return
    end if

    ns = ns_nsa(nsa)

    nth = nth_pnc(nsa)

    pl = pm(np,nsa)
    sl = SIN( thetam(nth-1,np,nr,nsa) )
    rl = rm(nr)
    deltath= delthm(nth-1,np,nr,nsa)

    Dtx(1)  = Dtpfow(nth,np,nr,nsa)*sl/(pl*sl*delp(nsa)*deltath*2.d0)
    Dtx(2)  = Dttfow(nth,np,nr,nsa)*sl/(pl**2*sl*deltath**2)
    Dtx(3)  = Dtrfow(nth,np,nr,nsa)*sl/(pl*sl*deltath*delr*2.d0)
    Fthp    = Fthfow(nth,np,nr,nsa)*sl

    loc_pnc = [nth,np,nr,nsa]
    loc_D   = [nth,np,nr,nsa]
    loc_Xstg= [nth_stg(nsa),np,nr_pnc_point(np,nr,nsa),nsa]

    if ( pz(ns) >= 0.d0 ) then
       ! calcurate sign of S_pnc_to_trp
       if ( np == npmax ) then
          dfnspdp = ( f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) ) &
                     / ( delp(ns) )
       else if ( np == 1 ) then
          dfnspdp = ( f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                     / ( delp(ns) )
       else
          dfnspdp = ( f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                     / ( 2.d0*delp(ns) )
       end if

       if ( nr == nrmax ) then
          dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                     / ( delr )  
       else if ( nr == 1 ) then
          dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                     / ( delr )
       else
          dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                     / ( 2.d0*delr )
       end if

       dfnspdth= ( fnsp(nth,np,nr,nsa)-fnsp(nth-1,np,nr,nsa) ) &
               / ( 0.5d0*pl*(delthm(nth,np,nr,nsa)+delthm(nth-1,np,nr,nsa)) )


       fnspg   = f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0)

       S_pnc_to_trp = -1.d0*( &
            Dtpfow(nth,np,nr,nsa)*dfnspdp  &
          + Dttfow(nth,np,nr,nsa)*dfnspdth &
          + Dtrfow(nth,np,nr,nsa)*dfnspdr  &
        ) + Fthfow(nth,np,nr,nsa)*fnspg

       if ( S_pnc_to_trp >= 0.d0 ) then
          return
       end if


       fact1 = (1.d0-IBCflux_ratio(np,nr,nsa))/JI(nth-1,np,nr,nsa)
       fact2 = IBCflux_ratio(np,nr,nsa) &
            /JI(nth_stg(nsa),np,nr_pnc_point(np,nr,nsa),nsa)
       imtx_D    = imtxa(nth-1,np,nr)
       imtx_Xstg = imtxa(nth_stg(nsa), np, nr_pnc_point(np,nr,nsa))
       dl(imtx_D) = dl(imtx_D) + fact1*( Dtx(2)+Fthp*weight(nth,np,nr,nsa) )

       do nl = 1, nl_D
          if ( ll(imtx_D,nl) == imtxa_boundary(nth,np+1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np+1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(1)*weight(nth,np+1,nr,nsa)*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np-1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np-1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(1)*weight(nth,np-1,nr,nsa)*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  -( Dtx(2) - Fthp*(1.d0-weight(nth,np,nr,nsa)) )*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np,nr+1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact1
          
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np,nr+1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(3)*weight(nth,np,nr+1,nsa)*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np,nr-1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np,nr-1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(3)*weight(nth,np,nr-1,nsa)*fact1
  
          end if
       end do
  
       nl_Xstg = nl_Xstg+1 
       ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np,nr)
       al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
            + fact2*( Dtx(2)+Fthp*weight(nth,np,nr,nsa) )
  
       if ( imtxa_boundary(nth,np+1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np+1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact2  
       end if
  
       if ( imtxa_boundary(nth-1,np+1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np+1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(1)*weight(nth,np+1,nr,nsa)*fact2  
       end if
  
       if ( imtxa_boundary(nth,np-1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np-1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact2  
       end if
  
       if (  imtxa_boundary(nth-1,np-1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np-1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(1)*weight(nth,np-1,nr,nsa)*fact2  
       end if
  
       if ( imtxa_boundary(nth,np,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               -( Dtx(2) + Fthp*(1.d0-weight(nth,np,nr,nsa)) )*fact2  
       end if
  
       if ( imtxa_boundary(nth,np,nr+1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np,nr+1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact2  
       end if
  
       if ( imtxa_boundary(nth-1,np,nr+1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np,nr+1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(3)*weight(nth,np,nr+1,nsa)*fact2  
       end if
  
       if ( imtxa_boundary(nth,np,nr-1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np,nr-1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact2  
       end if
  
       if ( imtxa_boundary(nth-1,np,nr-1)  /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np,nr-1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(3)*weight(nth,np,nr-1,nsa)*fact2  
       end if
  
    else ! PZ(ns) < 0
       if ( np == npmax ) then
          dfnspdp = ( f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                     / ( delp(ns) )
       else if ( np == 1 ) then
          dfnspdp = ( f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                     / ( delp(ns) )
       else
          dfnspdp = ( f_grid(loc_pnc, 0.d0, 1.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, -0.5d0, 0.5d0) )&
                     / ( 2.d0*delp(ns) )
       end if

       if ( nr == nrmax ) then
          dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                     / ( delr )  
       else if ( nr == 1 ) then
          dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0) )&
                     / ( delr )
       else
          dfnspdr = ( f_grid(loc_pnc, 0.d0, 0.5d0, 1.5d0) &
                     -f_grid(loc_pnc, 0.d0, 0.5d0, -0.5d0) )&
                     / ( 2.d0*delr )
       end if

       dfnspdth= ( fnsp(nth,np,nr,nsa)-fnsp(nth-1,np,nr,nsa) )&
            / ( 0.5d0*pl*(delthm(nth,np,nr,nsa)+delthm(nth-1,np,nr,nsa)) )

       fnspg   = f_grid(loc_pnc, 0.d0, 0.5d0, 0.5d0)

       S_pnc_to_pas = -1.d0*( &
            Dtpfow(nth,np,nr,nsa)*dfnspdp  &
          + Dttfow(nth,np,nr,nsa)*dfnspdth &
          + Dtrfow(nth,np,nr,nsa)*dfnspdr  &
          ) + Fthfow(nth,np,nr,nsa)*fnspg

       if ( S_pnc_to_pas <= 0.d0 ) then
          return 
       end if

       fact1 = -1.d0*IBCflux_ratio(np,nr,nsa)/JI(nth,np,nr,nsa)
       fact2 = -1.d0*(1.d0-IBCflux_ratio(np,nr,nsa)) &
            /JI(nth_stg(nsa)-1,np,nr_pnc_point(np,nr,nsa),nsa)
       imtx_D    = imtxa(nth,np,nr)
       imtx_Xstg = imtxa(nth_stg(nsa)-1, np, nr_pnc_point(np,nr,nsa))

       dl(imtx_D) = dl(imtx_D) + fact1*( -1.d0*Dtx(2) &
            +Fthp*(1.d0-weight(nth,np,nr,nsa)) )

       do nl = 1, nl_D
          if ( ll(imtx_D,nl) == imtxa_boundary(nth,np+1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np+1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(1)*weight(nth,np+1,nr,nsa)*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np-1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np-1,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(1)*weight(nth,np-1,nr,nsa)*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np,nr) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  +( Dtx(2) + Fthp*weight(nth,np,nr,nsa) )*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np,nr+1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact1
          
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np,nr+1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  - Dtx(3)*weight(nth,np,nr+1,nsa)*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth,np,nr-1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact1
  
          else if ( ll(imtx_D,nl) == imtxa_boundary(nth-1,np,nr-1) ) then
             al(imtx_D,nl) = al(imtx_D,nl) &
                  + Dtx(3)*weight(nth,np,nr-1,nsa)*fact1
  
          end if
       end do
  
       nl_Xstg = nl_Xstg+1 
       ll(imtx_Xstg,nl_Xstg) = imtx_D
       al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
            + fact2*( -1.d0*Dtx(2)+Fthp*(1.d0-weight(nth,np,nr,nsa)) )
  
       if ( imtxa_boundary(nth,np+1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np+1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact2  
       end if
  
       if ( imtxa_boundary(nth-1,np+1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np+1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(1)*weight(nth,np+1,nr,nsa)*fact2  
       end if
  
       if ( imtxa_boundary(nth,np-1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np-1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact2  
       end if
  
       if (  imtxa_boundary(nth-1,np-1,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np-1,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(1)*weight(nth,np-1,nr,nsa)*fact2  
       end if
  
       if ( imtxa_boundary(nth,np,nr) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np,nr)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               -( Dtx(2)  - Fthp*weight(nth,np,nr,nsa) )*fact2  
       end if
  
       if ( imtxa_boundary(nth,np,nr+1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np,nr+1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact2  
       end if
  
       if ( imtxa_boundary(nth-1,np,nr+1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np,nr+1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               - Dtx(3)*weight(nth,np,nr+1,nsa)*fact2  
       end if
  
       if ( imtxa_boundary(nth,np,nr-1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth,np,nr-1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact2  
       end if
  
       if ( imtxa_boundary(nth-1,np,nr-1) /= -1 ) then
          nl_Xstg = nl_Xstg+1
          ll(imtx_Xstg,nl_Xstg) = imtxa(nth-1,np,nr-1)
          al(imtx_Xstg,nl_Xstg) = al(imtx_Xstg,nl_Xstg) &
               + Dtx(3)*weight(nth,np,nr-1,nsa)*fact2  
       end if

    end if

  end subroutine IBC_pinch

  subroutine IBC_X_stagnation(np,nr,nsa,nl_Xstg,nl_pnc)
    implicit none
    integer,intent(in) :: np,nr,nsa
    integer,intent(inout) :: nl_Xstg, nl_pnc

    integer:: ns, nth, imtx_Xstg, imtx_pnc, nl
    REAL(rkind):: pl, sl, rl
    REAL(rkind) :: deltath, S_x, fact1, fact2
    REAL(rkind) :: Fthp, Dtx(3)
    REAL(rkind) :: dfnspdp, dfnspdth, dfnspdr, fnspg
    integer :: loc_Xstg(4), loc_pnc(4)

    ns = ns_nsa(nsa)

    if ( pz(ns) >= 0.d0 ) then
      
      nth = nth_stg(nsa)
      loc_Xstg = [nth,np,nr,nsa]
      if ( nl_pnc /= -1 ) then
        loc_pnc = [nth_pnc(nsa),np,nr_rhom_pinch(np,nr,nsa),nsa]
      end if

      pl = pm(np,nsa)
      sl = SIN( thetam(nth,np,nr,nsa) )
      rl = rm(nr)
      deltath= delthm(nth,np,nr,nsa)
  
      Dtx(1)  = Dtpfow(nth,np,nr,nsa)*sl/(pl*sl*delp(nsa)*deltath*2.d0)
      Dtx(2)  = Dttfow(nth,np,nr,nsa)*sl/(pl**2*sl*deltath**2)
      Dtx(3)  = Dtrfow(nth,np,nr,nsa)*sl/(pl*sl*deltath*delr*2.d0)
      Fthp    = Fthfow(nth,np,nr,nsa)*sl  

      if ( np == npmax ) then
         dfnspdp = ( f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 0.d0, -0.5d0, 0.5d0) )&
                    / ( delp(ns) )
      else if ( np == 1 ) then
         dfnspdp = ( f_grid(loc_Xstg, 0.d0, 1.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) )&
                    / ( delp(ns) )
      else
         dfnspdp = ( f_grid(loc_Xstg, 0.d0, 1.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 0.d0, -0.5d0, 0.5d0) )&
                    / ( 2.d0*delp(ns) )
      end if

      if ( nr == nrmax ) then
         dfnspdr = ( f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 0.d0, 0.5d0, -0.5d0) )&
                    / ( delr )  
      else if ( nr == 1 ) then
         dfnspdr = ( f_grid(loc_Xstg, 0.d0, 0.5d0, 1.5d0) &
                    -f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0) )&
                    / ( delr )
      else
         dfnspdr = ( f_grid(loc_Xstg, 0.d0, 0.5d0, 1.5d0) &
                    -f_grid(loc_Xstg, 0.d0, 0.5d0, -0.5d0) )&
                    / ( 2.d0*delr )
      end if

      dfnspdth= ( fnsp(nth,np,nr,nsa)-fnsp(nth-1,np,nr,nsa) )&
              / ( 0.5d0*pl*(delthm(nth,np,nr,nsa)+delthm(nth-1,np,nr,nsa)) )

      fnspg   = f_grid(loc_Xstg, 0.d0, 0.5d0, 0.5d0)

      S_x = -1.d0*( &
          Dtpfow(nth,np,nr,nsa)*dfnspdp  &
        + Dttfow(nth,np,nr,nsa)*dfnspdth &
        + Dtrfow(nth,np,nr,nsa)*dfnspdr  &
      ) + Fthfow(nth,np,nr,nsa)*fnspg

      imtx_Xstg = imtxa(nth, np, nr)

      if ( S_x <= 0.d0 .and. nl_pnc == -1 ) then ! orbital loss
        return

      else if ( S_x > 0.d0 .and. nl_pnc == -1 ) then

        fact1 = -1.d0/JI(nth,np,nr,nsa)

        dl(imtx_Xstg) = dl(imtx_Xstg) &
             + fact1*( -1.d0*Dtx(2)+Fthp*(1.d0-weight(nth,np,nr,nsa)) )

        do nl = 1, nl_Xstg
          if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np+1,nr) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth-1,np+1,nr) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  - Dtx(1)*weight(nth,np+1,nr,nsa)*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np-1,nr) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth-1,np-1,nr) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  + Dtx(1)*weight(nth,np-1,nr,nsa)*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np,nr) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  +( Dtx(2) + Fthp*weight(nth,np,nr,nsa) )*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np,nr+1) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact1
            
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth-1,np,nr+1) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  - Dtx(3)*weight(nth,np,nr+1,nsa)*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np,nr-1) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact1
    
          else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth-1,np,nr-1) ) then
             al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                  + Dtx(3)*weight(nth,np,nr-1,nsa)*fact1
    
          end if
        end do
    
        return

      else

        imtx_pnc=imtxa(loc_pnc(1),loc_pnc(2),loc_pnc(3))

        fact2 = 1.d0/JI(loc_pnc(1),loc_pnc(2),loc_pnc(3),loc_pnc(4))
        nl_pnc = nl_pnc+1 
        ll(imtx_pnc,nl_pnc) = imtx_Xstg
        al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
             + fact2*( -1.d0*Dtx(2)+Fthp*(1.d0-weight(nth,np,nr,nsa)) )
    
        if ( imtxa_boundary(nth,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np+1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact2  
        end if
    
        if ( imtxa_boundary(nth-1,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth-1,np+1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(1)*weight(nth,np+1,nr,nsa)*fact2  
        end if
    
        if ( imtxa_boundary(nth,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np-1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact2  
        end if
    
        if (  imtxa_boundary(nth-1,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth-1,np-1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(1)*weight(nth,np-1,nr,nsa)*fact2  
        end if
    
        if ( imtxa_boundary(nth,np,nr) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               -( Dtx(2) - Fthp*weight(nth,np,nr,nsa) )*fact2  
        end if
    
        if ( imtxa_boundary(nth,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np,nr+1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact2  
        end if
    
        if ( imtxa_boundary(nth-1,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth-1,np,nr+1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(3)*weight(nth,np,nr+1,nsa)*fact2  
        end if
    
        if ( imtxa_boundary(nth,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np,nr-1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact2  
        end if
    
        if ( imtxa_boundary(nth-1,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1
          ll(imtx_pnc,nl_pnc) = imtxa(nth-1,np,nr-1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(3)*weight(nth,np,nr-1,nsa)*fact2  
        end if  
  
      end if

    else ! pz(ns) < 0.d0
       
      nth = nth_stg(nsa)-1
      loc_Xstg = [nth,np,nr,nsa]
      if ( nl_pnc /= -1 ) then
        loc_pnc = [nth_pnc(nsa)-1,np,nr_rhom_pinch(np,nr,nsa),nsa]
      end if

      imtx_Xstg = imtxa(nth, np, nr)

      pl = pm(np,nsa)
      sl = SIN( thetam(nth,np,nr,nsa) )
      rl = rm(nr)
      deltath= delthm(nth,np,nr,nsa)
  
      Dtx(1)  = Dtpfow(nth+1,np,nr,nsa)*sl/(pl*sl*delp(nsa)*deltath*2.d0)
      Dtx(2)  = Dttfow(nth+1,np,nr,nsa)*sl/(pl**2*sl*deltath**2)
      Dtx(3)  = Dtrfow(nth+1,np,nr,nsa)*sl/(pl*sl*deltath*delr*2.d0)
      Fthp    = Fthfow(nth+1,np,nr,nsa)*sl  

      if ( np == npmax ) then
         dfnspdp = ( f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 1.d0, -0.5d0, 0.5d0) )&
                    / ( delp(ns) )
      else if ( np == 1 ) then
         dfnspdp = ( f_grid(loc_Xstg, 1.d0, 1.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) )&
                    / ( delp(ns) )
      else
         dfnspdp = ( f_grid(loc_Xstg, 1.d0, 1.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 1.d0, -0.5d0, 0.5d0) )&
                    / ( 2.d0*delp(ns) )
      end if

      if ( nr == nrmax ) then
         dfnspdr = ( f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) &
                    -f_grid(loc_Xstg, 1.d0, 0.5d0, -0.5d0) )&
                    / ( delr )  
      else if ( nr == 1 ) then
         dfnspdr = ( f_grid(loc_Xstg, 1.d0, 0.5d0, 1.5d0) &
                    -f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0) )&
                    / ( delr )
      else
         dfnspdr = ( f_grid(loc_Xstg, 1.d0, 0.5d0, 1.5d0) &
                    -f_grid(loc_Xstg, 1.d0, 0.5d0, -0.5d0) )&
                    / ( 2.d0*delr )
      end if

      dfnspdth= ( fnsp(nth+1,np,nr,nsa)-fnsp(nth,np,nr,nsa) )&
              / ( 0.5d0*pl*(delthm(nth+1,np,nr,nsa)+delthm(nth,np,nr,nsa)) )

      fnspg   = f_grid(loc_Xstg, 1.d0, 0.5d0, 0.5d0)

      S_x = -1.d0*( &
          Dtpfow(nth+1,np,nr,nsa)*dfnspdp  &
        + Dttfow(nth+1,np,nr,nsa)*dfnspdth &
        + Dtrfow(nth+1,np,nr,nsa)*dfnspdr  &
      ) + Fthfow(nth+1,np,nr,nsa)*fnspg

      if ( S_x >= 0.d0 .and. nl_pnc == -1 ) then ! orbital loss
        return

      else if ( S_x < 0.d0 .and. nl_pnc == -1 ) then

        fact1 = 1.d0/JI(nth,np,nr,nsa)
        imtx_Xstg = imtxa(nth, np, nr)

        dl(imtx_Xstg) = dl(imtx_Xstg) &
             + fact1*( Dtx(2)+Fthp*weight(nth+1,np,nr,nsa) )
        do nl = 1, nl_Xstg
           if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth+1,np+1,nr) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   - Dtx(1)*(1.d0-weight(nth+1,np+1,nr,nsa))*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np+1,nr) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   - Dtx(1)*weight(nth+1,np+1,nr,nsa)*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth+1,np-1,nr) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   + Dtx(1)*(1.d0-weight(nth+1,np-1,nr,nsa))*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np-1,nr) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   + Dtx(1)*weight(nth+1,np-1,nr,nsa)*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth+1,np,nr) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   -( Dtx(2) - Fthp*(1.d0-weight(nth+1,np,nr,nsa)) )*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth+1,np,nr+1) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   - Dtx(3)*(1.d0-weight(nth+1,np,nr+1,nsa))*fact1
            
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np,nr+1) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   - Dtx(3)*weight(nth+1,np,nr+1,nsa)*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth+1,np,nr-1) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   + Dtx(3)*(1.d0-weight(nth+1,np,nr-1,nsa))*fact1
    
           else if ( ll(imtx_Xstg,nl) == imtxa_boundary(nth,np,nr-1) ) then
              al(imtx_Xstg,nl) = al(imtx_Xstg,nl) &
                   + Dtx(3)*weight(nth+1,np,nr-1,nsa)*fact1
    
           end if
        end do
      
        return

     else

        imtx_pnc = imtxa(loc_pnc(1),loc_pnc(2),loc_pnc(3))

        fact1 = 1.d0/JI(loc_pnc(1),loc_pnc(2),loc_pnc(3),loc_pnc(4))

        nl_pnc = nl_pnc+1 
        ll(imtx_pnc,nl_pnc) = imtx_Xstg
        al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
             + fact1*( Dtx(2)+Fthp*weight(nth+1,np,nr,nsa) )

        if ( imtxa_boundary(nth+1,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth+1,np+1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(1)*(1.d0-weight(nth+1,np+1,nr,nsa))*fact1
  
        else if ( imtxa_boundary(nth,np+1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np+1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(1)*weight(nth+1,np+1,nr,nsa)*fact1
  
        else if ( imtxa_boundary(nth+1,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth+1,np-1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(1)*(1.d0-weight(nth+1,np-1,nr,nsa))*fact1
  
        else if ( imtxa_boundary(nth,np-1,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np-1,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(1)*weight(nth+1,np-1,nr,nsa)*fact1
  
        else if ( imtxa_boundary(nth+1,np,nr) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth+1,np,nr)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               -( Dtx(2) - Fthp*(1.d0-weight(nth+1,np,nr,nsa)) )*fact1
  
        else if ( imtxa_boundary(nth+1,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth+1,np,nr+1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(3)*(1.d0-weight(nth+1,np,nr+1,nsa))*fact1
          
        else if ( imtxa_boundary(nth,np,nr+1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np,nr+1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               - Dtx(3)*weight(nth+1,np,nr+1,nsa)*fact1
  
        else if ( imtxa_boundary(nth+1,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth+1,np,nr-1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(3)*(1.d0-weight(nth+1,np,nr-1,nsa))*fact1
  
        else if ( imtxa_boundary(nth,np,nr-1) /= -1 ) then
          nl_pnc = nl_pnc+1 
          ll(imtx_pnc,nl_pnc) = imtxa(nth,np,nr-1)
          al(imtx_pnc,nl_pnc) = al(imtx_pnc,nl_pnc) &
               + Dtx(3)*weight(nth+1,np,nr-1,nsa)*fact1
  
        end if
      
      end if

    end if

  end subroutine IBC_X_stagnation

  subroutine IBC_O_stagnation(np,nr,nsa,nl_Ostg)
    implicit none
    integer,intent(in) :: np,nr,nsa
    integer,intent(inout) :: nl_Ostg
    integer :: ns, nl, imtx_Ostg,nth
    real(rkind) :: fact
    REAL(rkind) :: Fthp, Dtx(3),pl,sl,deltath,rl

    ns = ns_nsa(nsa)

    if ( pz(ns) >= 0.d0 ) then
       nth = nth_stg(nsa)-1
       pl = pm(np,nsa)
       sl = SIN( thetam(nth,np,nr,nsa) )
       rl = rm(nr)
       deltath= delthm(nth,np,nr,nsa)
       imtx_Ostg = imtxa(nth,np,nr)
  
       Dtx(1)  = Dtpfow(nth+1,np,nr,nsa)*sl/(pl*sl*delp(nsa)*deltath*2.d0)
       Dtx(2)  = Dttfow(nth+1,np,nr,nsa)*sl/(pl**2*sl*deltath**2)
       Dtx(3)  = Dtrfow(nth+1,np,nr,nsa)*sl/(pl*sl*deltath*delr*2.d0)
       Fthp    = Fthfow(nth+1,np,nr,nsa)*sl

       fact = 1.d0/JI(nth,np,nr,nsa)

       dl(imtx_Ostg) = dl(imtx_Ostg) &
            + (Dtx(2)+Fthp*weight(nth+1,np,nr,nsa))*fact
      
       do nl = 1, nl_Ostg
          if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth+1,np+1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(1)*(1-weight(nth+1,np+1,nr,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np+1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(1)*weight(nth+1,np+1,nr,nsa)*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth+1,np-1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(1)*(1.d0-weight(nth+1,np-1,nr,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np+1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(1)*weight(nth+1,np-1,nr,nsa)*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth+1,np,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - (Dtx(2)+Fthp*(1.d0-weight(nth+1,np,nr,nsa)))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth+1,np,nr+1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(3)*(1.d0-weight(nth+1,np,nr+1,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np,nr+1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(3)*weight(nth+1,np,nr+1,nsa)*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth+1,np,nr-1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(3)*(1.d0-weight(nth+1,np,nr-1,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np,nr-1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(3)*weight(nth+1,np,nr-1,nsa)*fact
          end if
       end do

    else ! pz(ns) < 0.d0
       nth = nth_stg(nsa)
       pl = pm(np,nsa)
       sl = SIN( thetam(nth,np,nr,nsa) )
       rl = rm(nr)
       deltath= delthm(nth,np,nr,nsa)
       imtx_Ostg = imtxa(nth,np,nr)
  
       Dtx(1)  = Dtpfow(nth,np,nr,nsa)*sl/(pl*sl*delp(nsa)*deltath*2.d0)
       Dtx(2)  = Dttfow(nth,np,nr,nsa)*sl/(pl**2*sl*deltath**2)
       Dtx(3)  = Dtrfow(nth,np,nr,nsa)*sl/(pl*sl*deltath*delr*2.d0)
       Fthp    = Fthfow(nth,np,nr,nsa)*sl

       fact = -1.d0/JI(nth,np,nr,nsa)

       dl(imtx_Ostg) = dl(imtx_Ostg) &
            + fact*( -1.d0*Dtx(2)+Fthp*(1.d0-weight(nth,np,nr,nsa)) )
       do nl = 1, nl_Ostg
          if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np+1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(1)*(1.d0-weight(nth,np+1,nr,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth-1,np+1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(1)*weight(nth,np+1,nr,nsa)*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np-1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(1)*(1.d0-weight(nth,np-1,nr,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth-1,np-1,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(1)*weight(nth,np-1,nr,nsa)*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np,nr) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  +( Dtx(2)+Fthp*weight(nth,np,nr,nsa) )*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np,nr+1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(3)*(1.d0-weight(nth,np,nr+1,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth-1,np,nr+1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  - Dtx(3)*weight(nth,np,nr+1,nsa)*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth,np,nr-1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(3)*(1.d0-weight(nth,np,nr-1,nsa))*fact
          else if ( ll(imtx_Ostg,nl) == imtxa_boundary(nth-1,np,nr-1) ) then
             al(imtx_Ostg,nl) = al(imtx_Ostg,nl) &
                  + Dtx(3)*weight(nth,np,nr-1,nsa)*fact
          end if
       end do

    end if

  end subroutine IBC_O_stagnation

  function f_grid(loc,dth,dp,dr)
    ! dx = -1/2, 1/2, 3/2    -> half integer grid point
    ! dx = 0, 1              -> integer grid point
    implicit none
    real(rkind) :: f_grid
    integer,intent(in) :: loc(4)
    real(rkind),intent(in) :: dth, dp, dr
    integer :: nth, np, nr, nsa
    integer :: ith, ip, ir
    real(rkind) :: f_plus, f_minus

    nth = loc(1)
    np  = loc(2)
    nr  = loc(3)
    nsa = loc(4)

    if ( INT(2.d0*dth) == 2*INT(dth) ) then
      ith = nth+ INT(dth)
      ip  = np + INT(dp-0.5d0)
      ir  = nr + INT(dr-0.5d0)

      f_plus  = (1.d0-WEIGHT(ith, ip, ir, nsa)) * fnsp(ith,ip,ir,nsa)
      f_minus = WEIGHT(ith, ip, ir, nsa) * fnsp(ith-1,ip,ir,nsa)

    else if ( INT(2.d0*dp) == 2*INT(dp) ) then
      ith = nth + INT(dth-0.5d0)
      ip  = np  + INT(dp)
      ir  = nr  + INT(dr -0.5d0)

      f_plus  = (1.d0-WEIGHP(ith, ip, ir, nsa)) * fnsp(ith,ip,ir,nsa)
      f_minus = WEIGHP(ith, ip, ir, nsa) * fnsp(ith,ip-1,ir,nsa)

    else if ( INT(2.d0*dr) == 2*INT(dr) ) then
      ith = nth + INT(dth-0.5d0)
      ip  = np  + INT(dp -0.5d0)
      ir  = nr  + INT(dr)

      f_plus  = (1.d0-WEIGHR(ith, ip, ir, nsa)) * fnsp(ith,ip,ir,nsa)
      f_minus = WEIGHR(ith, ip, ir, nsa) * fnsp(ith,ip,ir-1,nsa)

    end if

    f_grid = f_plus + f_minus

  end function f_grid

  function imtxa_boundary(nth,np,nr) result(n)
    implicit none
    integer :: n
    integer,intent(in) :: nth,np,nr
    logical :: nth_,np_,nr_

    nth_ = 1 <= nth .and. nth <= nthmax 
    np_  = 1 <= np  .and. np  <= npmax 
    nr_  = 1 <= nr  .and. nr  <= nrmax 

    if ( nth_ .and. np_ .and. nr_ ) then
      n = imtxa(nth,np,nr)
    else
      n = -1
    end if

  end function imtxa_boundary

END MODULE fowexec
