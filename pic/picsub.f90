!  ***** TASK/PIC PREPARATION *****

MODULE picsub
  PRIVATE
  !INCLUDE 'fftw3.f'
  PUBLIC poisson_f,poisson_m,efield,bfield,wave,kine,pote,absorb_phi

CONTAINS

  !***********************************************************************
  SUBROUTINE poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
       rho,phi,rhof,phif,awk,afwk,cform,ipssn)
    !***********************************************************************
    IMPLICIT NONE
    !INCLUDE 'fftw3.f'
    REAL(8), DIMENSION(nxmax1,nymax1) :: rho,phi
    REAL(8), DIMENSION(nxmax,nymax) :: awk
    COMPLEx(8), DIMENSION(nxmaxh1,nymax) :: afwk
    COMPLEx(8), DIMENSION(nxmaxh1,nymax) :: rhof, phif
    REAL(8), DIMENSION(nxmaxh1,nymax) :: cform
    INTEGER(4) :: nxmax, nymax,nxmaxh1,nxmax1,nymax1,ipssn
    INTEGER(4) :: ifset

    IF(ipssn.EQ.0) THEN
       CALL poisson_sub(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)
       ifset = 0
       CALL fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,rho,rhof,awk,afwk,ifset)
    ELSE

       !.......... fourier transform rho
       ifset = -1
       CALL fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,rho,rhof,awk,afwk,ifset)

       !.......... calculate phi from rho in fourier space
       ipssn = 1
       CALL poisson_sub(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)

       !.......... inverse fourier transform phi
       ifset = 1
       CALL fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,phi,phif,awk,afwk,ifset)

    END IF
  END SUBROUTINE poisson_f

  !***********************************************************************
  SUBROUTINE poisson_sub(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)
    !***********************************************************************
    IMPLICIT NONE
    COMPLEx(8), DIMENSION(nxmaxh1,nymax) :: rhof, phif
    REAL(8), DIMENSION(nxmaxh1,nymax) :: cform
    INTEGER(4) :: nxmax, nymax, nxmaxh1
    INTEGER(4) :: nx, ny, nymaxh,ipssn
    REAL(8) :: alxi, alyi, pi, twopi, am, an, amn2, afsp, afsp2

    afsp  = 1.5 !+++ size of finite-size-particle
    afsp2 = afsp * afsp

    IF( ipssn .EQ. 0 ) THEN
       pi    = 3.14159265358979d0
       twopi = 2.d0 * pi
       alxi   = 1.d0 / DBLE(nxmax)
       alyi   = 1.d0 / DBLE(nymax)
       nymaxh    = nymax / 2
       !$omp parallel do private(am,an,amn2)
       DO ny = 1, nymaxh+1
          DO nx = 1, nxmaxh1
             am = twopi * DBLE(nx-1) * alxi
             an = twopi * DBLE(ny-1) * alyi
             amn2 = am*am + an*an

             IF( nx .EQ. 1 .AND. ny .EQ. 1 ) THEN
                cform(nx,ny) = 0.0
             ELSE
                cform(nx,ny) = 1.d0 / amn2 * ExP( - amn2 * afsp2 )
             ENDIF

             IF( ny .GE. 2 .AND. ny .LE. nymaxh ) THEN
                cform(nx,nymax-ny+2) = cform(nx,ny)
             ENDIF
          END DO
       END DO
       !$omp end parallel do
    ELSE

       !----- solve poisson equation
       !$omp parallel do
       DO ny = 1, nymax
          DO nx = 1, nxmaxh1
             phif(nx,ny) = cform(nx,ny) * rhof(nx,ny)
          END DO
       END DO
       !$omp end parallel do
    ENDIF

  END SUBROUTINE poisson_sub

  !***********************************************************************
  SUBROUTINE fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,a,af,awk,afwk,ifset)
    !***********************************************************************
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(8), DIMENSION(nxmax1,nymax1) :: a
    REAL(8), DIMENSION(nxmax,nymax) :: awk
    REAL(8) :: alx, aly
    COMPLEx(8), DIMENSION(nxmaxh1,nymax) :: af, afwk
    INTEGER(4) :: nxmax, nymax, nxmaxh1, nxmax1, nymax1
    INTEGER(4) :: ifset, nx, ny
    !....integer(8), save :: FFTW_ESTIMATE
    INTEGER(8), SAVE :: plan1, plan2

    alx = DBLE(nxmax)
    aly = DBLE(nymax)

    IF( ifset .EQ. 0 ) THEN
       !----- initialization of fourier transform
       CALL dfftw_plan_dft_r2c_2d(plan1,nxmax,nymax,awk,afwk,FFTW_ESTIMATE)
       CALL dfftw_plan_dft_c2r_2d(plan2,nxmax,nymax,afwk,awk,FFTW_ESTIMATE)

    ELSEIF( ifset .EQ. -1 ) THEN

       !----- fourier transform
       !$omp parallel do
       DO ny = 1, nymax
          DO nx = 1, nxmax
             awk(nx,ny) = a(nx,ny)
          END DO
       END DO
       !$omp end parallel do
       CALL dfftw_Execute(plan1)
       !$omp parallel do
       DO ny = 1, nymax
          DO nx = 1, nxmaxh1
             af(nx,ny) = afwk(nx,ny) / ( alx * aly )
          END DO
       END DO
       !$omp end parallel do
    ELSE

       !----- inverse fourier transform
       !$omp parallel do
       DO ny = 1, nymax
          DO nx = 1, nxmaxh1
             afwk(nx,ny) = af(nx,ny)
          END DO
       END DO
       !$omp end parallel do
       CALL dfftw_Execute(plan2)
       !$omp parallel do
       DO ny = 1, nymax
          DO nx = 1, nxmax
             a(nx,ny) = awk(nx,ny)
          END DO
       END DO
       !$omp end parallel do
       !$omp parallel do
       DO ny = 1, nymax
          a(nxmax1,ny) = a(1,ny)
       END DO
       !$omp end parallel do
       !$omp parallel do
       DO nx = 1, nxmax1
          a(nx,nymax1) = a(nx,1)
       END DO
       !$omp end parallel do
    ENDIF

  END SUBROUTINE fftpic

  !***********************************************************************
  SUBROUTINE poisson_m(nxmax1,nymax1,rho,phi,ipssn, &
       model_matrix0,model_matrix1,model_matrix2, &
       tolerance_matrix,model_boundary,dlen)
    !***********************************************************************
    USE libmpi
    USE commpi
    USE libmtx
    IMPLICIT NONE
    REAL(8), DIMENSION(nxmax1,nymax1) :: rho,phi
    REAL(8), DIMENSION(:),ALLOCATABLE :: x
    REAL(8):: tolerance_matrix,dlen,inv,xd,yd,xdmax,ydmax
    INTEGER :: nxmax1,nymax1,nxymax,ipssn
    INTEGER :: model_matrix0,model_matrix1,model_matrix2
    INTEGER :: nxmax,nymax,mode,imax,isize,jwidth,ileng
    INTEGER,SAVE:: status=0,istart,iend,irange
    INTEGER :: i,nx,ny,l,m,its,ilen
    INTEGER :: model_boundary

    nxmax=nxmax1-2
    nymax=nymax1-2
    imax=nxmax*nymax
    IF(nxmax.LE.nymax) THEN
       mode=0
       isize=nxmax
       ileng=nymax
       jwidth=4*nxmax-1
    ELSE
       mode=1
       isize=nymax
       ileng=nxmax
       jwidth=4*nymax-1
    END IF
    CALL mtx_setup(imax,istart,iend,nzmax=5*imax)
    irange=iend-istart+1
    status=1
    ALLOCATE(x(imax))
    DO i=istart,iend
       l=MOD(i-1,isize)+1
       m=(i-1)/isize+1
       IF(m.GT.1) CALL mtx_set_matrix(i,i-isize,1.d0)
       IF(l.GT.1) CALL mtx_set_matrix(i,i-1,1.d0)
       CALL mtx_set_matrix(i,i,-4.d0)
       IF(l.LT.isize) CALL mtx_set_matrix(i,i+1,1.d0)
       IF(m.LT.ileng) CALL mtx_set_matrix(i,i+isize,1.d0)
    ENDDO
    IF(mode.EQ.0) THEN
      !$omp parallel do private(i,nx,ny)
       DO i=istart,iend
          nx=MOD(i-1,isize)+1
          ny=(i-1)/isize+1
          CALL mtx_set_source(i,-rho(nx+1,ny+1))
       ENDDO
       !$omp end parallel do
    ELSE
      !$omp parallel do private(i,nx,ny)
       DO i=istart,iend
          ny=MOD(i-1,isize)+1
          nx=(i-1)/isize+1
          CALL mtx_set_source(i,-rho(nx+1,ny+1))
       ENDDO
       !$omp end parallel do
    END IF
    CALL mtx_solve(model_matrix0,tolerance_matrix,its, &
         methodKSP=model_matrix1,methodPC=model_matrix2)
    !      IF((its.ne.0).and.(nrank.eq.0)) write(6,*) 'mtx_solve: its=',its
    CALL mtx_gather_vector(x)
    IF(mode.EQ.0) THEN
      !$omp parallel do private(i,nx,ny)
       DO i=1,imax
          nx=MOD(i-1,isize)+1
          ny=(i-1)/isize+1
          phi(nx+1,ny+1)=x(i)
       ENDDO
       !$omp end parallel do
    ELSE
      !$omp parallel do private(i,nx,ny)
       DO i=1,imax
          ny=MOD(i-1,isize)+1
          nx=(i-1)/isize+1
          phi(nx+1,ny+1)=x(i)
       ENDDO
       !$omp end parallel do
    END IF
    DEALLOCATE(x)
    status=2
    CALL mtx_cleanup

  END SUBROUTINE poisson_m

  !***********************************************************************
  SUBROUTINE efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb,&
       Ex,Ey,Ez,Exb,Eyb,Ezb,Exbb,Eybb,Ezbb,Bxb,Byb,Bzb,Esx,Esy,Esz,Emx,Emy,Emz,&
       jx,jy,jz,vcfact,model_push,model_boundary)
  !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) ::Ex,Ey,Ez,Exb,Eyb,Ezb,Exbb,Eybb,Ezbb,&
    Bxb,Byb,Bzb,Esx,Esy,Esz,Emx,Emy,Emz
    REAL(8), DIMENSION(0:nxmax,0:nymax) ::phi,Ax,Ay,Az,Axb,Ayb,Azb,jx,jy,jz
    REAL(8):: dt,vcfact,dl,dm,dn
    INTEGER:: nxmax, nymax, nx, ny, nxm, nxp, nym, nyp
    INTEGER:: model_push, model_boundary,model_wg
    IF(model_boundary .EQ. 0) THEN
      !$omp parallel do private(nxm,nym,nxp,nyp)
       DO nx = 0, nxmax
          DO ny = 0, nymax

             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             IF( nx .EQ. 0  )    nxm = nxmax - 1
             IF( nx .EQ. nxmax ) nxp = 1
             IF( ny .EQ. 0  )    nym = nymax - 1
             IF( ny .EQ. nymax ) nyp = 1
            ! Esx(nx,ny) = phi(nx,ny) - phi(nxp,ny)
            ! Esy(nx,ny) = phi(nx,ny) - phi(nx,nyp)
            ! Esz(nx,ny) = 0.d0
            ! Emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
            ! Emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
            ! Emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt

             Esx(nx,ny)=Exb(nx,ny) - dt * jx(nx,ny)
             Esy(nx,ny)=Eyb(nx,ny) - dt * jy(nx,ny)
             Esz(nx,ny)=Ezb(nx,ny) - dt * jz(nx,ny)
             Emx(nx,ny)=dt*vcfact**2*(Bzb(nx,ny)-Bzb(nx,nym))
             Emy(nx,ny)=dt*vcfact**2*(Bzb(nxm,ny)-Bzb(nx,ny))
             Emz(nx,ny)=dt*vcfact**2*(Byb(nx,ny)-Byb(nxm,ny)-Bxb(nx,ny)+Bxb(nx,nym))
          END DO
       END DO
      !$omp end parallel do
    ELSE
      !$omp parallel do private(nxm,nym,nxp,nyp)
       DO nx = 1, nxmax
          DO ny = 1, nymax

             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             !IF( nx .EQ. -1 ) nxm = 0
             !IF( nx .EQ. nxmax ) nxp = nxmax
             !IF( ny .EQ. -1 ) nym = 0
             !IF( ny .EQ. nymax ) nyp = nymax
            !  Esx(nx,ny) = phi(nx,ny) - phi(nxp,ny)
            !  Esy(nx,ny) = phi(nx,ny) - phi(nx,nyp)
            !  Esz(nx,ny) = 0.d0
            !  Emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
            !  Emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
            !  Emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt
            Esx(nx,ny)=Exb(nx,ny) - dt * jx(nx,ny)
            Esy(nx,ny)=Eyb(nx,ny) - dt * jy(nx,ny)
            Esz(nx,ny)=Ezb(nx,ny) - dt * jz(nx,ny)
            Emx(nx,ny)=dt*vcfact**2*(Bzb(nx,ny)-Bzb(nx,nym))
            Emy(nx,ny)=dt*vcfact**2*(Bzb(nxm,ny)-Bzb(nx,ny))
            Emz(nx,ny)=dt*vcfact**2*(Byb(nx,ny)-Byb(nxm,ny)-Bxb(nx,ny)+Bxb(nx,nym))
            Ex(nx,ny) = Esx(nx,ny) + Emx(nx,ny)
            Ey(nx,ny) = Esy(nx,ny) + Emy(nx,ny)
            Ez(nx,ny) = Esz(nx,ny) + Emz(nx,ny)
        END DO
       END DO
      !$omp end parallel do
      !boundary condition for electro static field
       DO ny = 1, nymax-1
         nym = ny - 1
         IF(model_boundary .eq. 1) THEN !reflective
          Esx(0,ny) = Esx(1,ny)!Exb(0,ny) - dt * jx(0,ny)
          Esy(0,ny) = - Esy(1,ny)!Eyb(0,ny) - dt * jy(0,ny)
          Esz(0,ny) = - Esz(1,ny)
          Emx(0,ny) = Emx(1,ny)!dt*vcfact**2*(Bzb(0,ny)-Bzb(0,nym))
          Emy(0,ny) = - Emy(1,ny)!dt*vcfact**2*(-Bzb(0,ny))
          Emz(0,ny) = - Emz(1,ny)
          Esx(nxmax,ny)=Esx(nxmax-1,ny)!Exb(nxmax,ny) - dt * jx(nxmax,ny)
          Esy(nxmax,ny)=0.d0!Eyb(nxmax,ny) - dt * jy(nxmax,ny)
          Esz(nxmax,ny)=0.d0
          Emx(nxmax,ny)=Emx(nxmax-1,ny)!dt*vcfact**2*(Bzb(nxmax,ny)-Bzb(nxmax,nym))
          Emy(nxmax,ny)=0.d0!dt*vcfact**2*(Bzb(nxmax-1,ny)-Bzb(nxmax,ny))
          Emz(nxmax,ny)=0.d0
          Ex(0,ny) = Esx(0,ny) + Emx(0,ny)
          Ey(0,ny) = Esy(0,ny) + Emy(0,ny)
          Ez(0,ny) = Esz(0,ny) + Emz(0,ny)
          Ex(nxmax,ny) = Esx(nxmax,ny) + Emx(nxmax,ny)
          Ey(nxmax,ny) = Esy(nxmax,ny) + Emy(nxmax,ny)
          Ez(nxmax,ny) = Esz(nxmax,ny) + Emz(nxmax,ny)
        ELSE IF(model_boundary .eq. 2) THEN !absorbing
            Ex(0,ny)=Exb(0,ny)-dt*jx(0,ny)&
                    +dt*vcfact**2*(Bzb(0,ny)-Bzb(0,nym))
            Ey(0,ny) = 2.d0*Ey(1,ny)-Ey(2,ny)
            Ey(nxmax,ny) = 2.d0*Ey(nxmax-1,ny)-Ey(nxmax-2,ny)

          Ez(0,ny)=-Ezbb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                  *(Ez(1,ny)+Ezbb(0,ny)) &
                  +2.d0/(vcfact*dt+1.d0)*(Ezb(0,ny)+Ezb(1,ny))&
                  +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                  *(Ezb(0,ny+1)-2.d0*Ezb(0,ny)+Ezb(0,ny-1)+Ezb(1,ny+1)&
                  -2.d0*Ezb(1,ny)+Ezb(1,ny-1))

          Ez(nxmax,ny)=-Ezbb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                  *(Ez(nxmax-1,ny)+Ezbb(nxmax,ny)) &
                  +2.d0/(vcfact*dt+1.d0)*(Ezb(nxmax,ny)+Ezb(nxmax-1,ny))&
                  +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                  *(Ezb(nxmax,ny+1)-2.d0*Ezb(nxmax,ny)&
                  +Ezb(nxmax,ny-1)+Ezb(nxmax-1,ny+1)&
                  -2.d0*Ezb(nxmax-1,ny)+Ezb(nxmax-1,ny-1))
        ENDIF
       ENDDO
       DO nx = 1, nxmax-1
         nxm = nx - 1
         IF(model_boundary .eq. 1) then
         Esx(nx,0) = - Esx(nx,1)
         Esy(nx,0) = Esy(nx,1)!Eyb(nx,-1) - dt * jy(nx,-1)
         Esz(nx,0) = - Esz(nx,1)
         Emx(nx,0) = - Emx(nx,1)!dt*vcfact**2*Bzb(nx,-1)
         Emy(nx,0) = Emy(nx,1)!dt*vcfact**2*(Bzb(nxm,-1)-Bzb(nx,-1))
         Emz(nx,0) = - Emz(nx,1)
         Esx(nx,nymax) = 0.d0!Exb(nx,nymax) - dt * jx(nx,nymax)
         Esy(nx,nymax) = Esy(nx,nymax-1)!Eyb(nx,nymax) - dt * jy(nx,nymax)
         Esz(nx,nymax) = 0.d0
         Emx(nx,nymax) = 0.d0!dt*vcfact**2*(Bzb(nx,nymax)-Bzb(nx,nymax-1))
         Emy(nx,nymax) = Emy(nx,nymax-1)
         !dt*vcfact**2*(Bzb(nxm,nymax)-Bzb(nx,nymax))
         Emz(nx,nymax) = 0.d0
         Ex(nx,0) = Esx(nx,0) + Emx(nx,0)
         Ey(nx,0) = Esy(nx,0) + Emy(nx,0)
         Ez(nx,0) = Esz(nx,0) + Emz(nx,0)
         Ex(nx,nymax) = Esx(nx,nymax) + Emx(nx,nymax)
         Ey(nx,nymax) = Esy(nx,nymax) + Emy(nx,nymax)
         Ez(nx,nymax) = Esz(nx,nymax) + Emz(nx,nymax)
       ELSE IF(model_boundary .eq. 2) THEN
           Ex(nx,0) = 2.d0*Ex(nx,1)-Ex(nx,2)
           Ey(nx,0) = Eyb(nx,0)-dt*jy(nx,0)&
                     +dt*vcfact**2*(Bzb(nxm,0)-Bzb(nx,0))
           Ex(nx,nymax) = 2.d0*Ex(nx,nymax-1)-Ex(nx,nymax-2)

            Ez(nx,0)=-Ezbb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                    *(Ez(nx,1)+Ezbb(nx,0))&
                    +2.d0/(vcfact*dt+1.d0)*(Ezb(nx,0)+Ezb(nx,1))&
                    +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                    *(Ezb(nx+1,0)-2.d0*Ezb(nx,0)+Ezb(nx-1,0)+Ezb(nx+1,1)&
                     -2.d0*Ezb(nx,1)+Ezb(nx-1,1))

            Ez(nx,nymax)=-Ezbb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                    *(Ez(nx,nymax-1)+Ezbb(nx,nymax))&
                    +2.d0/(vcfact*dt+1.d0)*(Ezb(nx,nymax)+Ezb(nx,nymax-1))&
                    +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                    *(Ezb(nx+1,nymax)-2.d0*Ezb(nx,nymax)&
                    +Ezb(nx-1,nymax)+Ezb(nx+1,nymax-1)&
                    -2.d0*Ezb(nx,nymax-1)+Ezb(nx-1,nymax-1))

          ENDIF
       ENDDO
       IF(model_boundary .eq. 1) THEN
          Ex(:,0) = 0.d0
          Ex(:,nymax) = 0.d0
          Ey(0,:) = 0.d0
          Ey(nxmax,:) = 0.d0
          Ez(:,0) = 0.d0
          Ez(:,nymax) = 0.d0
          Ez(0,:) = 0.d0
          Ez(nxmax,:) = 0.d0
      ENDIF
           IF(model_boundary .eq. 2) THEN
            dl=(vcfact*dt-sqrt(2.d0))/(vcfact*dt+sqrt(2.d0))
            dm=2.d0*sqrt(2.d0)/(vcfact*dt+sqrt(2.d0))
            dn=4.d0*vcfact**2*dt**2/(vcfact*dt*sqrt(2.d0)+2.d0)

             Ez(0,0)=-Ezbb(1,1)+dl*(Ez(1,1)+Ezbb(0,0))+dm*(Ezb(1,1)+Ezb(0,0))&
                   +dn*(Ezb(0,1)-Ezb(1,1)-Ezb(0,0)+Ezb(1,0))

            Ez(nxmax,0)=-Ezbb(nxmax-1,1)+dl*(Ez(nxmax-1,1)+Ezbb(nxmax,0))+dm*(Ezb(nxmax-1,1)+Ezb(nxmax,0))&
                          +dn*(Ezb(nxmax,1)-Ezb(nxmax-1,1)-Ezb(nxmax,0)+Ezb(nxmax-1,0))

            Ez(0,nymax)=-Ezbb(1,nymax-1)+dl*(Ez(1,nymax-1)+Ezbb(0,nymax))+dm*(Ezb(1,nymax-1)+Ezb(0,nymax))&
                         +dn*(Ezb(0,nymax-1)-Ezb(1,nymax-1)-Ezb(0,nymax)+Ezb(1,nymax))

             Ez(nxmax,nymax)=-Ezbb(nxmax-1,nymax-1)&
                           +dl*(Ez(nxmax-1,nymax-1)+Ezbb(nxmax,nymax))&
                           +dm*(Ezb(nxmax-1,nymax-1)+Ezb(nxmax,nymax))&
                           +dn*(Ezb(nxmax,nymax-1)-Ezb(nxmax-1,nymax-1)&
                           -Ezb(nxmax,nymax)+Ezb(nxmax-1,nymax))
           ENDIF
ENDIF
      IF(model_boundary .eq. 0) THEN
       IF(MOD(model_push/2,2).EQ.0) THEN
          IF(MOD(model_push,2).EQ.0) THEN
             Ex(0:nxmax,0:nymax) = 0.D0
             Ey(0:nxmax,0:nymax) = 0.D0
             Ez(0:nxmax,0:nymax) = 0.D0
          ELSE
             Ex(0:nxmax,0:nymax) = Esx(0:nxmax,0:nymax)
             Ey(0:nxmax,0:nymax) = Esy(0:nxmax,0:nymax)
             Ez(0:nxmax,0:nymax) = Esz(0:nxmax,0:nymax)
          END IF
       ELSE
          IF(MOD(model_push,2).EQ.0) THEN
             Ex(0:nxmax,0:nymax) = Emx(0:nxmax,0:nymax)
             Ey(0:nxmax,0:nymax) = Emy(0:nxmax,0:nymax)
             Ez(0:nxmax,0:nymax) = Emz(0:nxmax,0:nymax)
          ELSE
            Ex(0:nxmax,0:nymax) = Esx(0:nxmax,0:nymax) + Emx(0:nxmax,0:nymax)
            Ey(0:nxmax,0:nymax) = Esy(0:nxmax,0:nymax) + Emy(0:nxmax,0:nymax)
            Ez(0:nxmax,0:nymax) = Esz(0:nxmax,0:nymax) + Emz(0:nxmax,0:nymax)
          END IF
       END IF
     ENDIF
     IF(model_boundary .eq. 2) THEN
      !  Ex(0,nymax) = Exb(0,nymax)-dt*jx(0,nymax)&
      !              +dt*vcfact**2*(Bzb(0,nymax)-Bzb(0,nymax-1))
      !  Ey(nxmax,0) = Eyb(nxmax,0)-dt*jy(nxmax,0)&
      !              +dt*vcfact**2*(Bzb(nxmax-1,0)-Bzb(nxmax,0))
       Ex(nxmax,:) = 0.d0
       Ey(:,nymax) = 0.d0
       Ex(0,0) = 2.d0*Ex(0,2) - Ex(0,1)
       Ey(0,0) = 2.d0*Ey(2,0) - Ey(1,0)
     ENDIF
  END SUBROUTINE efield

  !***********************************************************************
  SUBROUTINE bfield(nxmax,nymax,dt,Ax,Ay,Az,Axb,Ayb,Azb,&
       Ex,Ey,Ez,Bx,By,Bz,Bxb,Byb,Bzb,Bxbb,Bybb,Bzbb,Bxbg,Bybg,Bzbg,bb,vcfact, &
       model_push,model_boundary)
  !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nymax) :: Bxnab,Bynab,Bznab
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: Ex,Ey,Ez,Bx,By,Bz,Bxb,Byb,Bzb, &
                                           Bxbb,Bybb,Bzbb,Bxbg,Bybg,Bzbg,bb
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: Ax,Ay,Az,Axb,Ayb,Azb
    INTEGER :: nxmax, nymax, nx, ny, nxp, nyp, nxm, nym
    INTEGER :: model_push, model_boundary
    REAL(8) :: Bxx,Byy,Bzz,vcfact,dt,dl,dm,dn
    IF(model_boundary .EQ. 0) THEN
      !$omp parallel do private(nxm,nym,nxp,nyp)
       DO ny = 0, nymax
          DO nx = 0, nxmax
             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             IF( nx .EQ. 0  )    nxm = nxmax - 1
             IF( nx .EQ. nxmax ) nxp = 1
             IF( ny .EQ. 0  )    nym = nymax - 1
             IF( ny .EQ. nymax ) nyp = 1

              ! Bx(nx,ny) = 0.5d0 * (Az(nx,nyp) + Azb(nx,nyp) &
              !      - Az(nx,ny) - Azb(nx,ny))
              ! By(nx,ny) = - 0.5d0 * (Az(nxp,ny) + Azb(nxp,ny) &
              !      - Az(nx,ny) - Azb(nx,ny))
              ! Bz(nx,ny) = 0.5d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
              !      - Ay(nx,ny) - Ayb(nx,ny) &
              !      - (Ax(nx,nyp) + Axb(nx,nyp) &
              !      - Ax(nx,ny) - Axb(nx,ny)))

            Bx(nx,ny)=dt*(-Ez(nx,nyp)+Ez(nx,ny))+Bxb(nx,ny)
            By(nx,ny)=dt*(Ez(nxp,ny)-Ez(nx,ny))+Byb(nx,ny)
            Bz(nx,ny)=dt*(-Ey(nxp,ny)+Ey(nx,ny)+Ex(nx,nyp)-Ex(nx,ny))+Bzb(nx,ny)
            Bxbb(nx,ny)=Bxb(nx,ny)
            Bybb(nx,ny)=Byb(nx,ny)
            Bzbb(nx,ny)=Bzb(nx,ny)
            Bxx=Bx(nx,ny)
            Byy=By(nx,ny)
            Bzz=Bz(nx,ny)
            Bx(nx,ny)=0.5d0*(Bx(nx,ny)+Bxb(nx,ny))
            By(nx,ny)=0.5d0*(By(nx,ny)+Byb(nx,ny))
            Bz(nx,ny)=0.5d0*(Bz(nx,ny)+Bzb(nx,ny))
            Bxb(nx,ny)=Bxx
            Byb(nx,ny)=Byy
            Bzb(nx,ny)=Bzz
          END DO
       END DO
       !$omp end parallel do
    ELSE
      !$omp parallel do private(nxm,nym,nxp,nyp)
       DO ny = 0, nymax-1
          DO nx = 0, nxmax-1
             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             !IF( nx .EQ. 0  )    nxm = 0
             !IF( nx .EQ. nxmax ) nxp = nxmax
             !IF( ny .EQ. 0  )    nym = 0
             !IF( ny .EQ. nymax ) nyp = nymax

                ! Bx(nx,ny) =   0.5d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                !      - Az(nx,ny) - Azb(nx,ny))
                ! By(nx,ny) = - 0.5d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                !      - Az(nx,ny) - Azb(nx,ny))
                ! Bz(nx,ny) =   0.5d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                !      - Ay(nx,ny) - Ayb(nx,ny) &
                !      - (Ax(nx,nyp) + Axb(nx,nyp) &
                !      - Ax(nx,ny) - Axb(nx,ny)))
            Bx(nx,ny)=dt*(-Ez(nx,nyp)+Ez(nx,ny))+Bxb(nx,ny)
            By(nx,ny)=dt*(Ez(nxp,ny)-Ez(nx,ny))+Byb(nx,ny)
            Bz(nx,ny)=dt*(-Ey(nxp,ny)+Ey(nx,ny)+Ex(nx,nyp)-Ex(nx,ny))+Bzb(nx,ny)
          END DO
       END DO
       !$omp end parallel do
        DO ny = 1, nymax-1
          nyp = ny + 1
          IF(model_boundary .eq. 1) then  !reflective
              Bx(nxmax,ny)=0.d0!dt*(-Ez(nxmax,nyp)+Ez(nxmax,ny))+Bxb(nxmax,ny)
              By(nxmax,ny)=By(nxmax-1,ny)!By(nxmax,ny-1)+Bx(nxmax-1,ny)
              Bz(nxmax,ny)=2.d0*Bz(nxmax-1,ny)-Bz(nxmax-2,ny)!dt*(Ey(nxmax,ny)+Ex(nxmax,nyp)-Ex(nxmax,ny))+Bzb(nxmax,ny)
              !Bx(-1,ny)=-Bx(1,ny)!-dt*(Ez(0,nyp)-Ez(0,ny))+Bxb(0,ny)
              !By(-1,ny)=By(0,ny)!dt*(Ez(1,ny)-Ez(0,ny))+Byb(0,ny)
              !Bz(-1,ny)=Bz(0,ny)!dt*(-Ey(1,ny)+Ey(0,ny)+Ex(0,nyp)-Ex(0,ny))+Bzb(0,ny)
          ELSE IF(model_boundary .eq. 2) then !Mur's abosorbing boundary condition
            Bx(nxmax,ny) = dt*(-Ez(nxmax,nyp)+Ez(nxmax,ny))+Bxb(nxmax,ny)
            !By(nxmax,ny) = 2.d0*By(nxmax-1,ny)-By(nxmax-2,ny)

              Bz(0,ny)=-Bzbb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                      *(Bz(1,ny)+Bzbb(0,ny)) &
                      +2.d0/(vcfact*dt+1.d0)*(Bzb(0,ny)+Bzb(1,ny))&
                      +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                      *(Bzb(0,ny+1)-2.d0*Bzb(0,ny)+Bzb(0,ny-1)+Bzb(1,ny+1)&
                      -2.d0*Bzb(1,ny)+Bzb(1,ny-1))

              Bz(nxmax-1,ny)=-Bzbb(nxmax-2,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                          *(Bz(nxmax-2,ny)+Bzbb(nxmax-1,ny)) &
                          +2.d0/(vcfact*dt+1.d0)*(Bzb(nxmax-1,ny)+Bzb(nxmax-2,ny))&
                          +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                          *(Bzb(nxmax-1,ny+1)-2.d0*Bzb(nxmax-1,ny)&
                          +Bzb(nxmax-1,ny-1)+Bzb(nxmax-2,ny+1)&
                          -2.d0*Bzb(nxmax-2,ny)+Bzb(nxmax-2,ny-1))
              !Bz(nxmax,ny)=2.d0*Bz(nxmax-1,ny)-Bz(nxmax-2,ny)

           ENDIF
         ENDDO
         DO nx = 1, nxmax-1
           nxp = nx + 1
           IF(model_boundary .eq. 1) then
               Bx(nx,nymax)=Bx(nx,nymax-1)!Bx(nx-1,nymax)+By(nx,nymax-1)
               By(nx,nymax)=0.d0!dt*(Ez(nxp,nymax)-Ez(nx,nymax))+Byb(nx,nymax)
               Bz(nx,nymax)=2.d0*Bz(nx,nymax-1)-Bz(nx,nymax-2)!dt*(-Ex(nx,nymax)-Ey(nxp,nymax)+Ey(nx,nymax))+Bzb(nx,nymax)
               !Bx(nx,-1)=Bx(nx,0)!-dt*(Ez(nx,1)-Ez(nx,0))+Bxb(nx,0)
               !By(nx,-1)=-By(nx,1)!dt*(Ez(nxp,0)-Ez(nx,0))+Byb(nx,0)
               !Bz(nx,-1)=Bz(nx,0)
               !dt*(-Ey(nxp,-1)+Ey(nx,-1)+Ex(nx,0)-Ex(nx,-1))+Bzb(nx,-1)
           ELSE IF(model_boundary .eq. 2) then
              !Bx(nx,nymax) = 2.d0*Bx(nx,nymax-1)-Bx(nx,nymax-2)
              By(nx,nymax) = dt*(Ez(nxp,nymax)-Ez(nx,nymax))+Byb(nx,nymax)
              Bz(nx,0)=-Bzbb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                       *(Bz(nx,1)+Bzbb(nx,0))&
                       +2.d0/(vcfact*dt+1.d0)*(Bzb(nx,0)+Bzb(nx,1))&
                       +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                       *(Bzb(nx+1,0)-2.d0*Bzb(nx,0)+Bzb(nx-1,0)+Bzb(nx+1,1)&
                       -2.d0*Bzb(nx,1)+Bzb(nx-1,1))

              Bz(nx,nymax-1)=-Bzbb(nx,nymax-2)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                          *(Bz(nx,nymax-2)+Bzbb(nx,nymax-1))&
                          +2.d0/(vcfact*dt+1.d0)*(Bzb(nx,nymax-1)+Bzb(nx,nymax-2))&
                          +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                          *(Bzb(nx+1,nymax-1)-2.d0*Bzb(nx,nymax-1)&
                          +Bzb(nx-1,nymax-1)+Bzb(nx+1,nymax-2)&
                          -2.d0*Bzb(nx,nymax-2)+Bzb(nx-1,nymax-2))
              !Bz(nx,nymax)=2.d0*Bz(nx,nymax-1)-Bz(nx,nymax-2)

            ENDIF
        ENDDO
        IF(model_boundary .eq.2) THEN
          dl=(vcfact*dt-sqrt(2.d0))/(vcfact*dt+sqrt(2.d0))
          dm=2.d0*sqrt(2.d0)/(vcfact*dt+sqrt(2.d0))
          dn=4.d0*vcfact**2*dt**2/(vcfact*dt*sqrt(2.d0)+2.d0)
          Bz(0,0)=-Bzbb(1,1)+dl*(Bz(1,1)+Bzbb(0,0))+dm*(Bzb(1,1)+Bzb(0,0))&
                +dn*(Bzb(0,1)-Bzb(1,1)-Bzb(0,0)+Bzb(1,0))
          Bz(nxmax-1,0)=-Bzbb(nxmax-2,1)+dl*(Bz(nxmax-2,1)+Bzbb(nxmax-1,0))+dm*(Bzb(nxmax-2,1)+Bzb(nxmax-1,0))&
                       +dn*(Bzb(nxmax-1,1)-Bzb(nxmax-2,1)-Bzb(nxmax-1,0)+Bzb(nxmax-2,0))
          Bz(0,nymax-1)=-Bzbb(1,nymax-2)+dl*(Bz(1,nymax-2)+Bzbb(0,nymax-1))+dm*(Bzb(1,nymax-2)+Bzb(0,nymax-1))&
                      +dn*(Bzb(0,nymax-2)-Bzb(1,nymax-2)-Bzb(0,nymax-1)+Bzb(1,nymax-1))
          Bz(nxmax-1,nymax-1)=-Bzbb(nxmax-2,nymax-2)+dl*(Bz(nxmax-2,nymax-2)&
                        +Bzbb(nxmax-1,nymax-1))+dm*(Bzb(nxmax-2,nymax-2)&
                        +Bzb(nxmax-1,nymax-1))&
                        +dn*(Bzb(nxmax-1,nymax-2)-Bzb(nxmax-2,nymax-2)&
                        -Bzb(nxmax-1,nymax-1)+Bzb(nxmax-2,nymax-1))
        ENDIF
        IF(model_boundary .eq. 1) THEN
          Bx(0,:) = 0.d0
          Bx(nxmax,:) = 0.d0
          !Bx(:,nymax) = 0.d0
          By(:,0) = 0.d0
          By(:,nymax) = 0.d0
      ENDIF
      IF(model_boundary .EQ. 2) THEN
        Bx(nxmax,0) = dt*(-Ez(nxmax,1)+Ez(nxmax,0))+Bxb(nxmax,0)
        Bx(:,nymax) = 0.d0
        By(0,nymax) = dt*(Ez(1,nymax)-Ez(0,nymax))+Byb(0,nymax)
        By(nxmax,:) = 0.d0
        Bz(nxmax,:) = 0.d0
        Bz(:,nymax) = 0.d0
        !Bz(0,0) = 0.d0
        !Bz(0,nymax-1) = 0.d0
        !Bz(nxmax-1,0) = 0.d0
        !Bz(nxmax-1,nymax-1) = 0.d0
      ENDIF
        DO nx = 0,nxmax
          DO ny = 0,nymax
        Bxbb(nx,ny)=Bxb(nx,ny)
        Bybb(nx,ny)=Byb(nx,ny)
        Bzbb(nx,ny)=Bzb(nx,ny)
        Bxx=Bx(nx,ny)
        Byy=By(nx,ny)
        Bzz=Bz(nx,ny)
        Bx(nx,ny)=0.5d0*(Bx(nx,ny)+Bxb(nx,ny))
        By(nx,ny)=0.5d0*(By(nx,ny)+Byb(nx,ny))
        Bz(nx,ny)=0.5d0*(Bz(nx,ny)+Bzb(nx,ny))
        Bxb(nx,ny)=Bxx
        Byb(nx,ny)=Byy
        Bzb(nx,ny)=Bzz
          ENDDO
        ENDDO
      ENDIF
    IF(MOD(model_push/8,2).EQ.0) THEN
       IF(MOD(model_push/4,2).EQ.0) THEN
          Bx(0:nxmax,0:nymax) = 0.D0
          By(0:nxmax,0:nymax) = 0.D0
          Bz(0:nxmax,0:nymax) = 0.D0
       ELSE
          Bx(0:nxmax,0:nymax) = Bxbg(0:nxmax,0:nymax)
          By(0:nxmax,0:nymax) = Bybg(0:nxmax,0:nymax)
          Bz(0:nxmax,0:nymax) = Bzbg(0:nxmax,0:nymax)
       END IF
    ELSE
       IF(MOD(model_push/4,2).EQ.0) THEN
          CONTINUE
       ELSE
          Bx(0:nxmax,0:nymax) = Bx(0:nxmax,0:nymax) + Bxbg(0:nxmax,0:nymax)
          By(0:nxmax,0:nymax) = By(0:nxmax,0:nymax) + Bybg(0:nxmax,0:nymax)
          Bz(0:nxmax,0:nymax) = Bz(0:nxmax,0:nymax) + Bzbg(0:nxmax,0:nymax)
       END IF
    END IF

    bb(0:nxmax,0:nymax) = SQRT(Bx(0:nxmax,0:nymax)**2 &
         +By(0:nxmax,0:nymax)**2 &
         +Bz(0:nxmax,0:nymax)**2)

  END SUBROUTINE bfield

  !***********************************************************************
  SUBROUTINE wave(nxmax,nymax,dt,Ex,Ey,Ez,vcfact,xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi)
  !***********************************************************************
     IMPLICIT NONE
     REAL(8), DIMENSION(0:nxmax,0:nymax) :: ex,ey,ez
     INTEGER :: nxmax,nymax,ny
     REAL(8) :: dt,xmin_wg,xmax_wg,ymin_wg,ymax_wg,yc,ylen,dph,y,factor,&
     amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi,vcfact,amp_start
     yc=0.5d0*(ymin_wg+ymax_wg)
     ylen=(ymax_wg-ymin_wg)
     IF(ylen .NE. 0) dph=ph_wg/ylen
     amp_start = time * vcfact/10.d0
     IF(amp_start .GE. 1.d0) amp_start=1.0d0
     DO ny=0,nymax
        y=DBLE(ny)
        IF(y.GE.ymin_wg.AND.y.LE.ymax_wg) THEN
           factor=ExP(-12.D0*(y-yc)**2/(ylen)**2)
           Ey(1,ny)=amp_wg*amp_start &
                *factor*COS(rot_wg*pi/180.D0) &
                *SIN(omega*time-pi*dph*(y-ymin_wg)/180.D0)
           Ez(0,ny)=amp_wg*amp_start &
                *factor*SIN(rot_wg*pi/180.D0) &
                *SIN(omega*time-pi*dph*(y-ymin_wg)/180.D0)
        END IF
    END DO
  END SUBROUTINE wave
  !***********************************************************************
  SUBROUTINE kine(npmax,vx,vy,vz,akin,mass,vcfact)
  !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(npmax) :: vx, vy, vz
    REAL(8) :: akin, mass, vcfact, gamma
    INTEGER(4) :: npmax, np
    akin = 0.d0
    !$omp parallel do Private(gamma)&
    !$omp reduction(+:akin)
    DO np = 1, npmax
      gamma = 1.d0/sqrt(1.d0 - (vx(np)**2+vy(np)**2+vz(np)**2)/vcfact**2)
       akin = akin + vcfact ** 2 * (gamma - 1)
    END DO
    !$omp end parallel do
    IF(npmax.EQ.0) THEN
       akin=0.D0
    ELSE
       akin = akin * mass  / DBLE(npmax)
    END IF
  END SUBROUTINE kine

  !***********************************************************************
  SUBROUTINE pote(nxmax,nymax,Ex,Ey,Ez,Bx,By,Bz,Bxbg,Bybg,Bzbg,vcfact, &
       apote,apotm)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: Ex,Ey,Ez,Bx,By,Bz,Bxbg,Bybg,Bzbg
    REAL(8) :: apote,apotm,vcfact
    INTEGER(4) :: nxmax, nymax, nx, ny

    apote = 0.d0
    apotm = 0.d0
    !$omp parallel do reduction(+:apote,apotm)
    DO ny = 1, nymax-1
       DO nx = 1, nxmax-1
          apote = apote + (Ex(nx,ny)**2 + Ey(nx,ny)**2 + Ez(nx,ny)**2)
          apotm = apotm + ((Bx(nx,ny)-Bxbg(nx,ny))**2 &
               + (By(nx,ny)-Bybg(nx,ny))**2 &
               + (Bz(nx,ny)-Bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    !$omp parallel do reduction(+:apote,apotm)
    DO nx = 0, nxmax,nxmax
       DO ny = 1, nymax-1
          apote = apote + 0.5D0*(Ex(nx,ny)**2+Ey(nx,ny)**2+Ez(nx,ny)**2)
          apotm = apotm + 0.5D0*((Bx(nx,ny)-Bxbg(nx,ny))**2 &
               + (By(nx,ny)-Bybg(nx,ny))**2 &
               + (Bz(nx,ny)-Bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    !$omp parallel do reduction(+:apote,apotm)
    DO ny = 0, nymax,nymax
       DO nx = 1, nxmax-1
          apote = apote + 0.5D0*(Ex(nx,ny)**2+Ey(nx,ny)**2+Ez(nx,ny)**2)
          apotm = apotm + 0.5D0*((Bx(nx,ny)-Bxbg(nx,ny))**2 &
               + (By(nx,ny)-Bybg(nx,ny))**2 &
               + (Bz(nx,ny)-Bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    !$omp parallel do reduction(+:apote,apotm)
    DO ny = 0, nymax,nymax
       DO nx = 0, nxmax,nxmax
          apote = apote + 0.25D0*(Ex(nx,ny)**2+Ey(nx,ny)**2+Ez(nx,ny)**2)
          apotm = apotm + 0.25D0*((Bx(nx,ny)-Bxbg(nx,ny))**2 &
               + (By(nx,ny)-Bybg(nx,ny))**2 &
               + (Bz(nx,ny)-Bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    apote = 0.5D0 * apote / (DBLE(nxmax)*DBLE(nymax))
    apotm = 0.5D0 * vcfact**2 * apotm / (DBLE(nxmax)*DBLE(nymax))
 END SUBROUTINE pote

END MODULE picsub
