!***************************************************************
!
!   Command loop
!
!***************************************************************

subroutine TXMENU

  use tx_commons, only : allocate_txcomm, deallocate_txcomm, &
       &              rMU0, IERR, rMUb1, rMUb2, TMAX, DT, NTMAX, T_TX, PNBHT1, &
       &              PNBHT2, PNBHP, PNBHex, X, LQm4, NRMAX, TPRE, PNBCD, &
       &              iflag_file, gkilo, &
       &              DelRho, DelN, Rho, LQe1, LQi1, LQz1, achg, ICONT, IRPIN, array_init_NR
  use tx_main, only : TXEXEC
  use tx_graphic, only : TX_GRAPH_SAVE, TXSTGR, TXGOUT, allocate_txgraf, deallocate_txgraf, &
       &                 NGR, NGRM, GT, GY
  use tx_variables, only : TXCALV
  use tx_parameter_control, only : TXPARM_CHECK, TXPARM, TXVIEW
  use tx_interface, only : TXKLIN, TOUPPER, TXLOAD
  use tx_ripple, only : ripple_input, ripple_spl, deallocate_ripple
  use tx_ntv, only : dealloc_ntv
  use eqread_mod, only : alloc_equ
  use mod_cross_section, only : deallocate_spline_table_carbon_rate_coef_adas, &
                                deallocate_spline_table_beam_rate_coef
  implicit none
  integer(4) :: MODE, I, IST, ier, NR, IER_RP
  character(len=80) :: LINE
  character(len=1)  :: KID, KID2

  IERR  = 0
  ICONT = 0
  IRPIN = 0 ! Read ripple file if IRPIN /= 0
  if((PNBHP + PNBHT1 + PNBHT2 + PNBHex) /= 0) then
     rMUb1 = 1.D0
     rMUb2 = rMU0
     GY%gnrm(4) = gkilo
  end if

  !  *** MENU ***

  do
     ier = 0
     if(ICONT == 0) TMAX = DT*NTMAX

     write(6,'(3(A,ES12.4))') &
          &   ' ## TIME=',T_TX,'  DT=',DT,'  NEXT TIME =',TMAX
     write(6,*) '## INPUT: ', &
          &   'R:RUN  C:CONT  P,V:PARM  G:GRAPH  '// &
          &   'W:STAT  S:SAVE  L:LOAD  I:INIT '
     write(6,'(11X,A)') 'F,FR:FILE  N:PTRB  M:ITG  O:OUT  Q:QUIT'
     flush(6)

     call TXKLIN(LINE,KID,MODE)
     if(MODE /= 1) then
        call TXPARM_CHECK
        TMAX=T_TX+DT*NTMAX
        if(rMUb1 == rMU0 .and. &
        & (PNBHT1 /= 0.D0 .or. PNBHT2 /= 0.D0 .or. PNBHP /= 0.D0 .or. PNBHex /= 0.D0)) then
           rMUb1 = 1.D0
           rMUb2 = rMU0
           GY%gnrm(4) = gkilo
           if(allocated(X)) X(0:NRMAX,LQm4) = X(0:NRMAX,LQm4) / rMUb2
        end if
        cycle
     end if

     select case(KID)
     case('R')
        if (ICONT /= 0) then
           write(6,*) '# Would you like to restart? [y/N]'
           read(5,'(A1)',iostat=IST) KID
           if(IST /= 0) cycle
           call TOUPPER(KID)
           if(KID /= 'Y') cycle
        end if
        call allocate_txcomm(ier, icont)
        call allocate_txgraf(ier, icont)
        if(ier /= 0) cycle
        T_TX = 0.D0
        TPRE = 0.D0
        IERR = 0
        ICONT = 1
        call TXPROF
        call TX_GRAPH_SAVE
        if(IRPIN /= 0) call ripple_spl
        call TXEXEC
        TMAX=T_TX+DT*NTMAX
     case('C')
        if (ICONT == 0) then
           write(6,*) 'XX RUN or LOAD before CONTINUE !'
           cycle
        end if
        NGR=-1
        call TXSTGR(NGR,GT,GY,NGRM)
        if(IRPIN /= 0) call ripple_spl
        call TXEXEC
        TMAX=T_TX+DT*NTMAX
     case('P')
        call TXPARM(KID)
        if(KID == 'Q') exit
     case('V')
        call TXVIEW
     case('I')
        call TXINIT
     case('Q')
        exit
     case('W')
        if (ICONT == 0) then
           write(6,*) 'XX RUN or LOAD before CONTINUE !'
           cycle
        end if
        call TXSTAT
     case('S')
        call TXSAVE
     case('L')
        call TXLOAD(IER)
        if(IER /= 0) cycle
        IERR = 0
        NGR = -1
        ICONT = 1
        call TX_GRAPH_SAVE
     case('G')
        call TXGOUT
     case('N')
        if (ICONT == 0) then
           write(6,*) 'XX RUN or LOAD before CONTINUE !'
           cycle
        end if
        do NR = 0, NRMAX-1
           if(Rho(NR) <= DelRho  .and.  Rho(NR+1) >= DelRho) then
              I = NR
              exit
           end if
        end do
        X(LQe1,I) = X(LQe1,I) + DelN
        X(LQi1,I) = X(LQi1,I) + DelN * achg(2)
        X(LQz1,I) = X(LQz1,I) + DelN * achg(3)
        call TXCALV(X)
     case('B')
        KID2=LINE(2:2)
        call TOUPPER(KID2)
        select case(KID2)
        case('P')
           PNBCD =  1.D0
        case('0')
           PNBCD =  0.D0
        case('M')
           PNBCD = -1.D0
        case default
           write(6,*) 'XX Unknown beam command'
        end select
     case('F')
        KID2=LINE(2:2)
        call TOUPPER(KID2)
        select case(KID2)
        case('R')
           call ripple_input(IER_RP)
           if(IER_RP == 0) IRPIN = 1
        case default
           call ascii_input

           if(iflag_file == 1) then
              if(rMUb2 == 1.D0) then
                 rMUb1 = 1.D0
                 rMUb2 = rMU0
                 GY%gnrm(4) = gkilo

                 if(allocated(X)) X(LQm4,0:NRMAX) = X(LQm4,0:NRMAX) / rMUb2
                 if(T_TX /= 0.D0) call TXCALV(X)
              end if
           end if
        end select
     case('O')
        if (ICONT == 0) then
           write(6,*) 'XX RUN or LOAD before CONTINUE !'
           cycle
        end if
        call outfile
     case('M')
        if (ICONT == 0) then
           write(6,*) 'XX RUN or LOAD before CONTINUE !'
           cycle
        end if
        call ITG_growthrate
     case('#')
        continue
     case default
        write(6,*) 'XX Unknown command'
     end select
  end do

  ! *** deallocate dynamic arrays before termination ***

  call deallocate_ripple
  call deallocate_spline_table_carbon_rate_coef_adas
  call deallocate_spline_table_beam_rate_coef
  call dealloc_ntv
  call deallocate_txcomm
  call deallocate_txgraf
  call alloc_equ(-1)

end subroutine TXMENU

!***** INPUT KID or LINE *****
!           MODE=0: LINE INPUT 
!                1: KID INPUT
!                2: PARM INPUT
!                3: NEW PROMPT

subroutine TXKLIN(LINE,KID,MODE)

  use tx_parameter_control, only : TXPARL
  use tx_interface, only : TOUPPER
  implicit none

  integer(4), intent(out) :: MODE
  character(len=80), intent(out) :: LINE
  character(len=1), intent(out) :: KID
  integer(4) :: ID, I, IST

  !  ----- read line input input -----

  read(5,'(A80)',iostat=IST) LINE

  if(IST == 0) then
     !  ----- parameter input -----

     ID=0
     do I=1,80
        if(LINE(I:I) == '=') ID=1
     end do
     if(ID == 1) then
        call TXPARL(LINE)
        KID=' '
        MODE=2
        return
     end if

     !  ----- command input -----

     KID=LINE(1:1)
     call TOUPPER(KID)
     if(KID >= 'A' .and. KID <= 'Z') then
        MODE=1
        return
     end if

     !  ----- line input -----

     KID=' '
     MODE=0

  else if(IST < 0) then

     !  ----- input end -----

     KID='Q'
     MODE=1

  else if(IST > 0) then

     !  ----- input error -----

     write(6,*) 'XX INPUT ERROR !'
     MODE=3

  end if

end subroutine TXKLIN
