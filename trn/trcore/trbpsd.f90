MODULE trbpsd

  use bpsd
  TYPE(bpsd_device_type),  PRIVATE,SAVE :: device
  TYPE(bpsd_species_type), PRIVATE,SAVE :: species
  TYPE(bpsd_equ1D_type),   PRIVATE,SAVE :: equ1D
  TYPE(bpsd_metric1D_type),PRIVATE,SAVE :: metric1D
  TYPE(bpsd_plasmaf_type), PRIVATE,SAVE :: plasmaf
  LOGICAL, PRIVATE, SAVE :: tr_bpsd_init_flag = .TRUE.

  PUBLIC

CONTAINS
!=======================================================================
  SUBROUTINE tr_bpsd_init(ierr)

    USE trcomm
    ! local variables
    INTEGER(ikind),INTENT(OUT) :: ierr
    INTEGER(ikind) :: ns

      IF(tr_bpsd_init_flag) THEN
         species%nsmax = 0
         plasmaf%nsmax = 0
         plasmaf%nrmax = 0
         tr_bpsd_init_flag = .FALSE.
      ENDIF

!      write(6,*) 'top of tr_bpsd_init:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_init:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      ! --- deivece ---
      device%rr   = RR
      device%zz   = 0.d0
      device%ra   = RA
      device%rb   = RA+0.1D0
      device%bb   = BB
      device%ip   = RIP
      device%elip = RKAP
      device%trig = RDLT
      CALL bpsd_set_data(device,ierr)

      ! --- species ---
      IF(species%nsmax.ne.nsmax) THEN
         IF(ASSOCIATED(species%data)) THEN
            DEALLOCATE(species%data)
         ENDIF
         species%nsmax = nsmax
         ALLOCATE(species%data(species%nsmax))
      ENDIF

      DO ns=1,species%nsmax
         species%data(ns)%pa  = pa(ns)
         species%data(ns)%pz  = pz(ns)
         species%data(ns)%pz0 = pz(ns)
      ENDDO
      CALL bpsd_set_data(species,ierr)

!!$      ! --- equ1D ---
!!$      IF((equ1D%nrmax.ne.nrmax+1)) THEN
!!$         IF(ASSOCIATED(equ1D%rho)) THEN
!!$            DEALLOCATE(equ1D%rho)
!!$         ENDIF
!!$         IF(ASSOCIATED(equ1D%data)) THEN
!!$            DEALLOCATE(equ1D%data)
!!$         ENDIF
!!$         equ1D%nrmax = nrmax+1
!!$         ALLOCATE(equ1D%rho(1:equ1D%nrmax))
!!$         ALLOCATE(equ1D%data(1:equ1D%nrmax))
!!$      ENDIF
!!$
!!$      ! --- metric1D ---
!!$      IF((metric1D%nrmax.ne.nrmax+1)) THEN
!!$         IF(ASSOCIATED(metric1D%rho)) THEN
!!$            DEALLOCATE(metric1D%rho)
!!$         ENDIF
!!$         IF(ASSOCIATED(metric1D%data)) THEN
!!$            DEALLOCATE(metric1D%data)
!!$         ENDIF
!!$         metric1D%nrmax = nrmax+1
!!$         ALLOCATE(metric1D%rho(1:metric1D%nrmax))
!!$         ALLOCATE(metric1D%data(1:metric1D%nrmax))
!!$      ENDIF

      ! --- plasmaf ---
      IF((plasmaf%nsmax.ne.nsmax).or. &
     &   (plasmaf%nrmax.ne.nrmax)) THEN
         IF(ASSOCIATED(plasmaf%rho)) THEN
            DEALLOCATE(plasmaf%rho)
         ENDIF
         IF(ASSOCIATED(plasmaf%data)) THEN
            DEALLOCATE(plasmaf%data)
         ENDIF
         IF(ASSOCIATED(plasmaf%qinv)) THEN
            DEALLOCATE(plasmaf%qinv)
         ENDIF
         plasmaf%nsmax = nsmax
         plasmaf%nrmax = nrmax+1
         ALLOCATE(plasmaf%rho(plasmaf%nrmax))
         ALLOCATE(plasmaf%data(plasmaf%nrmax,plasmaf%nsmax))
         ALLOCATE(plasmaf%qinv(plasmaf%nrmax))
      ENDIF

      CALL bpsd_set_data(plasmaf,ierr)

      RETURN
    END SUBROUTINE tr_bpsd_init

!=======================================================================
    SUBROUTINE tr_bpsd_set(ierr)

      USE trcomm
      INTEGER(ikind),INTENT(out) :: ierr

      ! local variables
      INTEGER(ikind) :: ns,nr

!      write(6,*) 'top of tr_bpsd_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_set: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      ! --- device ---
      device%rr   = RR
      device%zz   = 0.d0
      device%ra   = RA
      device%rb   = RA+0.1D0
      device%bb   = BB
      device%ip   = RIP
      device%elip = RKAP
      device%trig = RDLT
      CALL bpsd_set_data(device,ierr)

      ! --- plasmaf ---
      plasmaf%time = t

      DO nr=1,plasmaf%nrmax
         plasmaf%rho(nr)  = rhog(nr-1)
         plasmaf%qinv(nr) = 1.d0 / qp(nr-1)

         DO ns=1, nsamax
            plasmaf%data(nr,ns)%pn   = rn(ns,nr-1)*1.d20
            plasmaf%data(nr,ns)%pt   = rt(ns,nr-1)*1.d3
            plasmaf%data(nr,ns)%ptpr = rt(ns,nr-1)*1.d3
            plasmaf%data(nr,ns)%ptpp = rt(ns,nr-1)*1.d3
            plasmaf%data(nr,ns)%pu   = ru(ns,nr-1)
         ENDDO
      ENDDO

!      write(6,*) 'END of tr_bpsd_set:qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'END of tr_bpsd_set:rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      CALL bpsd_set_data(plasmaf,ierr)
      RETURN
    END SUBROUTINE tr_bpsd_set

!=======================================================================
    SUBROUTINE tr_bpsd_get(ierr)

      USE trcomm
      INTEGER(ikind),INTENT(out) :: ierr

!      write(6,*) 'top of tr_bpsd_get: qp'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      write(6,*) 'top of tr_bpsd_get: rt'
!      write(6,'(1P5E12.4)') (rt(nr,1),nr=1,nrmax)
!      pause

      CALL bpsd_get_data(device,ierr)
      RR   = device%rr
      RA   = device%ra
      BB   = device%bb
      RIP  = device%ip
      RKAP = device%elip
      RDLT = device%trig

!      IF(modelg.eq.3.or.modelg.eq.5) THEN
!         RIPS = RIP
!         RIPE = RIP
!      ENDIF

      modelg = mdlgmt
      ! TASK/EQ or EQDSK output geometry, call TOPICS/EQU or TASK/EQ
      IF(modelg.eq.9) THEN

         CALL bpsd_get_data(equ1D,ierr)
         CALL bpsd_get_data(metric1D,ierr)

!         write(*,*) nrmax+1
!         write(*,*) equ1D%nrmax
         IF(nrmax+1 /= equ1D%nrmax .OR. nrmax+1 /= metric1D%nrmax) THEN
            CALL get_spl_on(ierr) ! spline interpolation
         ELSE
            CALL get_spl_off(ierr)!
         END IF

      ENDIF

      RETURN
    END SUBROUTINE tr_bpsd_get

! **************************************************************************

    SUBROUTINE get_spl_off(ierr)
      USE trcomm
      INTEGER(ikind),INTENT(OUT) :: ierr
      ! local variables
      INTEGER(ikind) :: nr
      REAL(rkind)    :: psita

      ierr = 0 ! dummy

      DO nr = 1, equ1D%nrmax
         psitrho(nr-1) = equ1D%data(nr)%psit
         psiprho(nr-1) = equ1D%data(nr)%psip
         ppprho(nr-1)  = equ1D%data(nr)%ppp
         piqrho(nr-1)  = equ1D%data(nr)%piq
         ttrho(nr-1)   = equ1D%data(nr)%pip * rmu0/(2.d0*pi)
         pirho(nr-1)   = equ1D%data(nr)%pit
      ENDDO
      psita = equ1D%data(equ1D%nrmax)%psit
      qp(0:nrmax) = 1.d0/piqrho(0:nrmax)

      DO nr = 1, metric1D%nrmax
         rhog(nr-1)    = equ1D%rho(nr)
         pvolrho(nr-1) = metric1D%data(nr)%pvol      ! Plasma volume
         psurrho(nr-1) = metric1D%data(nr)%psur      ! Plasma surface
         dvrho(nr-1)   = metric1D%data(nr)%dvpsit * (2.d0*psita*rhog(nr-1))
         ! dV/drho
         rdpvrho(nr-1) = 1.d0 / (metric1D%data(nr)%dvpsip * 2.d0*pi)
         ! dpsi/dV
         !             metric1D%data(nr)%aver2       ! <R^2>
         arrho(nr-1)   = metric1D%data(nr)%aver2i    ! <1/R^2>
         abb2rho(nr-1) = metric1D%data(nr)%aveb2     ! <B^2>
         aib2rho(nr-1) = metric1D%data(nr)%aveb2i    ! <1/B^2>
         !             metric1D%data(nr)%avegv       ! <|grad V|>
         !             metric1D%data(nr)%avegv2      ! <|grad V|^2>
         abvrho(nr-1)  = metric1D%data(nr)%avegvr2   ! <|grad V|^2/R^2>
         ar1rho(nr-1)  = metric1D%data(nr)%avegr     ! <|grad rho|>
         ar2rho(nr-1)  = metric1D%data(nr)%avegr2    ! <|grad rho|^2>
         abrho(nr-1)   = metric1D%data(nr)%avegrr2   ! <|grad rho|^2/R^2>
         !             metric1D%data(nr)%avegpp2     ! <|grad Psip|^2>
         rmjrho(nr-1)  = metric1D%data(nr)%rr        ! local R
         rmnrho(nr-1)  = metric1D%data(nr)%rs        ! local r
         rkprho(nr-1)  = metric1D%data(nr)%elip      ! elipticity
         !             metric1D%data(nr)%trig        ! triangularity
         abb1rho(nr-1) = metric1D%data(nr)%aveb      ! <B>
         
         epsrho(nr-1)  = rmnrho(nr-1) / rmjrho(nr-1) ! rs/rr
      ENDDO

      rhom(1:nrmax) = 0.5d0*(rhog(0:nrmax-1)*rhog(1:nrmax))

      RETURN
    END SUBROUTINE get_spl_off

! --------------------------------------------------------------------------

    SUBROUTINE get_spl_on(ierr)
      USE trcomm
      INTEGER(ikind),INTENT(OUT) :: ierr
      ! local variables
      INTEGER(ikind) :: nr
      REAL(rkind)    :: psita

      REAL(rkind),DIMENSION(metric1D%nrmax) :: rhoeg,deriv,temp
      REAL(rkind),DIMENSION(4,metric1D%nrmax) :: urhog

      REAL(rkind),DIMENSION(0:nrmax) :: rhoeg_tr

      ! --- get 'rhog' ---
      DO nr = 1, metric1D%nrmax
         rhoeg(nr) = dble(nr-1)/dble(metric1D%nrmax-1)
      END DO
      DO nr = 0, nrmax
         rhoeg_tr(nr) = dble(nr)/dble(nrmax)
      END DO

      deriv=0.d0
      CALL SPL1D(rhoeg,metric1D%rho,deriv,urhog,metric1D%nrmax,0,ierr)
      DO nr = 0, nrmax
         CALL SPL1DF(rhoeg_tr(nr),rhog(nr),metric1D%rho,urhog, &
                     metric1D%nrmax,ierr)
      END DO
      rhog(nrmax) = 1.d0
      ! ------------------

      temp(1:equ1D%nrmax) = equ1D%data(1:equ1D%nrmax)%psit
      CALL spl1d_tr_array(equ1D%rho,temp,rhog,psitrho,equ1D%nrmax,1,ierr)

      temp(1:equ1D%nrmax) = equ1D%data(1:equ1D%nrmax)%psip
      CALL spl1d_tr_array(equ1D%rho,temp,rhog,psiprho,equ1D%nrmax,1,ierr)

      temp(1:equ1D%nrmax) = equ1D%data(1:equ1D%nrmax)%ppp
      CALL spl1d_tr_array(equ1D%rho,temp,rhog,ppprho,equ1D%nrmax,1,ierr)

      temp(1:equ1D%nrmax) = equ1D%data(1:equ1D%nrmax)%piq
      CALL spl1d_tr_array(equ1D%rho,temp,rhog,piqrho,equ1D%nrmax,1,ierr)

      temp(1:equ1D%nrmax) = equ1D%data(1:equ1D%nrmax)%pip
      CALL spl1d_tr_array(equ1D%rho,temp,rhog,ttrho,equ1D%nrmax,1,ierr)
      ttrho(0:nrmax) = ttrho(0:nrmax) * rmu0/(2.d0*pi)

      temp(1:equ1D%nrmax) = equ1D%data(1:equ1D%nrmax)%pit
      CALL spl1d_tr_array(equ1D%rho,temp,rhog,pirho,equ1D%nrmax,1,ierr)

      psita       = psitrho(nrmax)
      qp(0:nrmax) = 1.d0/piqrho(0:nrmax)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%pvol
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,pvolrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%psur
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,psurrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%dvpsit
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,dvrho,metric1D%nrmax,1,ierr)
      dvrho(0:nrmax) = dvrho(0:nrmax) * (2.d0*psita*rhog(0:nrmax))

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%dvpsip
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,rdpvrho,metric1D%nrmax,1,ierr)
      rdpvrho(0:nrmax) = 1.d0/(rdpvrho(0:nrmax)*2.d0*pi)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%aver2i
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,arrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%aveb2
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,abb2rho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%aveb2i
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,aib2rho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%avegvr2
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,abvrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%avegr
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,ar1rho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%avegr2
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,ar2rho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%avegrr2
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,abrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%rr
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,rmjrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%rs
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,rmnrho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%elip
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,rkprho,metric1D%nrmax,1,ierr)

      temp(1:metric1D%nrmax) = metric1D%data(1:metric1D%nrmax)%aveb
      CALL spl1d_tr_array(metric1D%rho,temp,rhog,abb1rho,metric1D%nrmax,1,ierr)

      epsrho(0:nrmax)  = rmnrho(0:nrmax) / rmjrho(0:nrmax)

      rhom(1:nrmax) = 0.5d0*(rhog(0:nrmax-1)*rhog(1:nrmax))

      RETURN
    END SUBROUTINE get_spl_on


    SUBROUTINE spl1d_tr_array(X,F,rhog,FU,NXMAX,ID,IERR)
! -------------------------------------------------------------------------
!      INPUT : X(NXMAX)  : Coordinates
!              F(NXMAX)  : Value
!              NXMAX     : Number of variables
!              ID        : 0 : Second derivatives = 0 at X(1) and X(NXMAX)
!                          1 : Derivative FX(1) is 0.
!
!      OUTPUT: FU(0:nrmax) : Value on grid of TR 
!              IERR        : Error indicator
! -------------------------------------------------------------------------
      USE trcomm, ONLY: nrmax

      IMPLICIT NONE
      INTEGER(ikind),                INTENT(IN)   :: NXMAX, ID
      REAL(rkind), DIMENSION(NXMAX), INTENT(IN)   :: X, F
      REAL(rkind), DIMENSION(0:nrmax), INTENT(IN) :: rhog
      REAL(rkind), DIMENSION(0:nrmax), INTENT(OUT):: FU
      INTEGER(ikind),                  INTENT(OUT):: IERR
      
      REAL(rkind),DIMENSION(NXMAX)    :: FX ! EDGE DERIVATIVE FOR 1<= ID <=3      
      REAL(rkind), DIMENSION(4,NXMAX) :: U  ! SPLINE COEFICIENTS
      INTEGER(ikind) :: nr
     
      FX = 0.d0      
      CALL SPL1D(X,F,FX,U,NXMAX,ID,IERR)
      
      DO nr = 0, nrmax
         CALL SPL1DF(rhog(nr),FU(nr),X,U,NXMAX,IERR)
      END DO
      
      RETURN
    END SUBROUTINE spl1d_tr_array


         ! TASK/EQ or EQDSK output geometry
         ! Calibration in order to keep consistency 
         !                          between metrics and plasma current
!         IF(modelg.eq.3.or.modelg.eq.5) THEN
!            RIP  = ABVRHOG(NRMAX)*RDPVRHOG(NRMAX)/(2.D0*PI*RMU0)*1.D-6
!            RIPS = RIP
!            RIPE = RIP
!         ENDIF

END MODULE trbpsd
