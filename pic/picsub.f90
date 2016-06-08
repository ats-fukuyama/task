!  ***** TASK/PIC PREPARATION *****

MODULE picsub
  PRIVATE
  !INCLUDE 'fftw3.f'
  PUBLIC poisson_f,poisson_m,efield,bfield,kine,pote,absorb_phi

CONTAINS

  !***********************************************************************
  SUBROUTINE poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
       rho,phi,rhof,phif,awk,afwk,cform,ipssn)
    !***********************************************************************
    IMPLICIT NONE
    !INCLUDE 'fftw3.f'
    REAL(8), DIMENSION(nxmax1,nymax1) :: rho,phi
    REAL(8), DIMENSION(nxmax,nymax) :: awk
    COMPLEX(8), DIMENSION(nxmaxh1,nymax) :: afwk
    COMPLEX(8), DIMENSION(nxmaxh1,nymax) :: rhof, phif
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
    COMPLEX(8), DIMENSION(nxmaxh1,nymax) :: rhof, phif
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
                cform(nx,ny) = 1.d0 / amn2 * EXP( - amn2 * afsp2 )
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
    COMPLEX(8), DIMENSION(nxmaxh1,nymax) :: af, afwk
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
       CALL dfftw_execute(plan1)
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
       CALL dfftw_execute(plan2)
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
  SUBROUTINE absorb_phi(nxmax,nymax,phi,phib,phibb,dt,vcfact)
  !***********************************************************************
  IMPLICIT NONE
  REAL(8),DIMENSION(0:nxmax,0:nymax) :: phi,phib,phibb
  REAL(8) :: dt,vcfact
  INTEGER :: nxmax,nymax,nx,ny
  !aborobing boundary condition for phi

  DO nx = 1, nxmax-1
    phi(nx,0)=-phibb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            *(phi(nx,1)+phibb(nx,0)) &
            +2.d0/(vcfact*dt+1.d0)*(phib(nx,0)+phib(nx,1))&
            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            *(phib(nx+1,0)-2.d0*phib(nx,0)+phib(nx-1,0)+phib(nx+1,1)&
             -2.d0*phib(nx,1)+phib(nx-1,1))
    phi(nx,nymax)=-phibb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            *(phi(nx,nymax-1)+phibb(nx,nymax)) &
            +2.d0/(vcfact*dt+1.d0)*(phib(nx,nymax)+phib(nx,nymax-1))&
            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            *(phib(nx+1,nymax)-2.d0*phib(nx,nymax)+phib(nx-1,nymax)&
            +phib(nx+1,nymax-1)-2.d0*phib(nx,nymax-1)+phib(nx-1,nymax-1))
  END DO
  DO ny = 1, nymax-1
    phi(nxmax,ny)=-phibb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            *(phi(nxmax-1,ny)+phibb(nxmax,ny)) &
            +2.d0/(vcfact*dt+1.d0)*(phib(nxmax,ny)+phib(nxmax-1,ny))&
            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            *(phib(nxmax,ny+1)-2.d0*phib(nxmax,ny)&
            +phib(nxmax,ny-1)+phib(nxmax-1,ny+1)&
            -2.d0*phib(nxmax-1,ny)+phib(nxmax-1,ny-1))
  ENDDO

  END SUBROUTINE absorb_phi

  !***********************************************************************
  SUBROUTINE efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb,&
       ex,ey,ez,exb,eyb,ezb,exbb,eybb,ezbb,bxb,byb,bzb,esx,esy,esz,emx,emy,emz,&
       jx,jy,jz,vcfact,model_wg,xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi,model_push,model_boundary)
  !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) ::phi,Ax,Ay,Az,Axb,Ayb,Azb,ex,ey,ez,exb,eyb,ezb,exbb,eybb,ezbb,&
    bxb,byb,bzb,esx,esy,esz,emx,emy,emz,jx,jy,jz
    REAL(8):: dt,vcfact,dl,dm,dn,amp_wg,rot_wg,eli_wg,&
    omega,time,pi,amp_start,y,yc,ylen,factor,dph,ymin_wg,ymax_wg,xmin_wg,xmax_wg,ph_wg
    INTEGER :: nxmax, nymax, nx, ny, nxm, nxp, nym, nyp
    INTEGER:: model_push, model_boundary,model_wg
    IF(model_boundary .EQ. 0) THEN
      !$omp parallel do private(nx,ny,nxm,nym,nxp,nyp)
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
            ! esx(nx,ny) = phi(nx,ny) - phi(nxp,ny)
            ! esy(nx,ny) = phi(nx,ny) - phi(nx,nyp)
            ! esz(nx,ny) = 0.d0
            ! emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
            ! emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
            ! emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt

             esx(nx,ny)=ex(nx,ny) - dt * jx(nx,ny)
             esy(nx,ny)=ey(nx,ny) - dt * jy(nx,ny)
             esz(nx,ny)=ez(nx,ny) - dt * jz(nx,ny)
             emx(nx,ny)=dt*vcfact**2*(bzb(nx,ny)-bzb(nx,nym))
             emy(nx,ny)=dt*vcfact**2*(bzb(nxm,ny)-bzb(nx,ny))
             emz(nx,ny)=dt*vcfact**2*(byb(nx,ny)-byb(nxm,ny)-bxb(nx,ny)+bxb(nx,nym))
          END DO
       END DO
      !$omp end parallel do
    ELSE IF (model_boundary .NE. 0) THEN
      !$omp parallel do private(nx,ny,nxm,nym,nxp,nyp)
       DO nx = 1, nxmax-1
          DO ny = 1, nymax-1

             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             !IF( nx .EQ. 0  )    nxm = 0
             !IF( nx .EQ. nxmax ) nxp = nxmax
             !IF( ny .EQ. 0  )    nym = 0
             !IF( ny .EQ. nymax ) nyp = nymax
            !  esx(nx,ny) = phi(nx,ny) - phi(nxp,ny)
            !  esy(nx,ny) = phi(nx,ny) - phi(nx,nyp)
            !  esz(nx,ny) = 0.d0
            !  emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
            !  emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
            !  emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt
            esx(nx,ny)=ex(nx,ny) - dt * jx(nx,ny)
            esy(nx,ny)=ey(nx,ny) - dt * jy(nx,ny)
            esz(nx,ny)=ez(nx,ny) - dt * jz(nx,ny)
            emx(nx,ny)=dt*vcfact**2*(bzb(nx,ny)-bzb(nx,nym))
            emy(nx,ny)=dt*vcfact**2*(bzb(nxm,ny)-bzb(nx,ny))
            emz(nx,ny)=dt*vcfact**2*(byb(nx,ny)-byb(nxm,ny)-bxb(nx,ny)+bxb(nx,nym))
            ex(nx,ny) = esx(nx,ny) + emx(nx,ny)
            ey(nx,ny) = esy(nx,ny) + emy(nx,ny)
            ez(nx,ny) = esz(nx,ny) + emz(nx,ny)
        END DO
       END DO
      !$omp end parallel do
      !boundary condition for electro static field
       DO ny = 1, nymax-1
         nym = ny - 1
         IF(model_boundary .eq. 1) THEN !reflective
          esx(0,ny) = ex(0,ny) - dt * jx(0,ny)
          esy(0,ny) = ey(0,ny) - dt * jy(0,ny)
          esz(0,ny) = ez(0,ny) - dt * jz(0,ny)
          emx(0,ny) = dt*vcfact**2*(bzb(0,ny)-bzb(0,nym))
          emy(0,ny) = 0.d0
          emz(0,ny) = dt*vcfact**2*(-bxb(0,ny)+bxb(0,nym))
          esx(nxmax,ny) = esx(nxmax-1,ny)
          esy(nxmax,ny) = esy(nxmax-1,ny)
          esz(nxmax,ny) = esz(nxmax-1,ny)
          emx(nxmax,ny) = emx(nxmax-1,ny)
          emy(nxmax,ny) = emy(nxmax-1,ny)
          emz(nxmax,ny) = emz(nxmax-1,ny)
          ex(0,ny) = esx(0,ny) + emx(0,ny)
          ey(0,ny) = esy(0,ny) + emy(0,ny)
          ez(0,ny) = esz(0,ny) + emz(0,ny)
          ex(nxmax,ny) = esx(nxmax,ny) + emx(nxmax,ny)
          ey(nxmax,ny) = esy(nxmax,ny) + emy(nxmax,ny)
          ez(nxmax,ny) = esz(nxmax,ny) + emz(nxmax,ny)
        ELSE IF(model_boundary .eq. 2) THEN
          !  Ex(nxmax,ny)=-Exbb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !          *(Ex(nxmax-1,ny)+Exbb(nxmax,ny)) &
          !          +2.d0/(vcfact*dt+1.d0)*(Exb(nxmax,ny)+Exb(nxmax-1,ny))&
          !          +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !          *(Exb(nxmax,ny+1)-2.d0*Exb(nxmax,ny)&
          !          +Exb(nxmax,ny-1)+Exb(nxmax-1,ny+1)&
          !          -2.d0*Exb(nxmax-1,ny)+Exb(nxmax-1,ny-1))
           !
          !  Ey(nxmax,ny)=-Eybb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !          *(Ey(nxmax-1,ny)+Eybb(nxmax,ny)) &
          !          +2.d0/(vcfact*dt+1.d0)*(Eyb(nxmax,ny)+Eyb(nxmax-1,ny))&
          !          +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !          *(Eyb(nxmax,ny+1)-2.d0*Eyb(nxmax,ny)&
          !          +Eyb(nxmax,ny-1)+Eyb(nxmax-1,ny+1)&
          !          -2.d0*Eyb(nxmax-1,ny)+Eyb(nxmax-1,ny-1))

          Ez(nxmax,ny)=-Ezbb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                  *(Ez(nxmax-1,ny)+Ezbb(nxmax,ny)) &
                  +2.d0/(vcfact*dt+1.d0)*(Ezb(nxmax,ny)+Ezb(nxmax-1,ny))&
                  +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                  *(Ezb(nxmax,ny+1)-2.d0*Ezb(nxmax,ny)&
                  +Ezb(nxmax,ny-1)+Ezb(nxmax-1,ny+1)&
                  -2.d0*Ezb(nxmax-1,ny)+Ezb(nxmax-1,ny-1))

        !  Ex(0,ny)=-Exbb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
        !           *(Ex(1,ny)+Exbb(0,ny)) &
        !           +2.d0/(vcfact*dt+1.d0)*(Exb(0,ny)+Exb(1,ny))&
        !           +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
        !           *(Exb(0,ny+1)-2.d0*Exb(0,ny)&
        !           +Exb(0,ny-1)+Exb(1,ny+1)&
        !           -2.d0*Exb(1,ny)+Exb(1,ny-1))
         !
        !  Ey(0,ny)=-Eybb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
        !           *(Ey(1,ny)+Eybb(0,ny)) &
        !           +2.d0/(vcfact*dt+1.d0)*(Eyb(0,ny)+Eyb(1,ny))&
        !           +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
        !           *(Eyb(0,ny+1)-2.d0*Eyb(0,ny)&
        !           +Eyb(0,ny-1)+Eyb(1,ny+1)&
        !           -2.d0*Eyb(1,ny)+Eyb(1,ny-1))

        Ez(0,ny)=-Ezbb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                 *(Ez(1,ny)+Ezbb(0,ny)) &
                 +2.d0/(vcfact*dt+1.d0)*(Ezb(0,ny)+Ezb(1,ny))&
                 +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                 *(Ezb(0,ny+1)-2.d0*Ezb(0,ny)&
                 +Ezb(0,ny-1)+Ezb(1,ny+1)&
                 -2.d0*Ezb(1,ny)+Ezb(1,ny-1))
        ENDIF
       ENDDO
       DO nx = 1, nxmax-1
         nxm = nx - 1
         IF(model_boundary .eq. 1) then
         esx(nx,0) = ex(nx,0) - dt * jx(nx,0)
         esy(nx,0) = ey(nx,0) - dt * jy(nx,0)
         esz(nx,0) = ez(nx,0) - dt * jz(nx,0)
         emx(nx,0) = 0.d0
         emy(nx,0) = dt*vcfact**2*(bzb(nxm,0)-bzb(nx,0))
         emz(nx,0) = dt*vcfact**2*(byb(nx,0)-byb(nxm,0))
         esx(nx,nymax) = esx(nx,nymax-1)
         esy(nx,nymax) = esy(nx,nymax-1)
         esz(nx,nymax) = esz(nx,nymax-1)
         emx(nx,nymax) = emx(nx,nymax-1)
         emy(nx,nymax) = emy(nx,nymax-1)
         emz(nx,nymax) = emz(nx,nymax-1)
         ex(nx,0) = esx(nx,0) + emx(nx,0)
         ey(nx,0) = esy(nx,0) + emy(nx,0)
         ez(nx,0) = esz(nx,0) + emz(nx,0)
         ex(nx,nymax) = esx(nx,nymax) + emx(nx,nymax)
         ey(nx,nymax) = esy(nx,nymax) + emy(nx,nymax)
         ez(nx,nymax) = esz(nx,nymax) + emz(nx,nymax)
         ELSE IF(model_boundary .eq. 2) then
            ! Ex(nx,0)=-Exbb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !         *(Ex(nx,1)+Exbb(nx,0)) &
            !         +2.d0/(vcfact*dt+1.d0)*(Exb(nx,0)+Exb(nx,1))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !         *(Exb(nx+1,0)-2.d0*Exb(nx,0)+Exb(nx-1,0)+Exb(nx+1,1)&
            !          -2.d0*Exb(nx,1)+Exb(nx-1,1))
            !
            ! Ey(nx,0)=-Eybb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !         *(Ey(nx,1)+Eybb(nx,0))&
            !         +2.d0/(vcfact*dt+1.d0)*(Eyb(nx,0)+Eyb(nx,1))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !         *(Eyb(nx+1,0)-2.d0*Eyb(nx,0)+Eyb(nx-1,0)+Eyb(nx+1,1)&
            !          -2.d0*Eyb(nx,1)+Eyb(nx-1,1))

           Ez(nx,0)=-Ezbb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
                   *(Ez(nx,1)+Ezbb(nx,0))&
                   +2.d0/(vcfact*dt+1.d0)*(Ezb(nx,0)+Ezb(nx,1))&
                   +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
                   *(Ezb(nx+1,0)-2.d0*Ezb(nx,0)+Ezb(nx-1,0)+Ezb(nx+1,1)&
                    -2.d0*Ezb(nx,1)+Ezb(nx-1,1))

            ! Ex(nx,nymax)=-Exbb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !          *(Ex(nx,nymax-1)+Exbb(nx,nymax)) &
            !         +2.d0/(vcfact*dt+1.d0)*(Exb(nx,nymax)+Exb(nx,nymax-1))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !         *(Exb(nx+1,nymax)-2.d0*Exb(nx,nymax) &
            !         +Exb(nx-1,nymax)+Exb(nx+1,nymax-1)&
            !         -2.d0*Exb(nx,nymax-1)+Exb(nx-1,nymax-1))
            !
            ! Ey(nx,nymax)=-Eybb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !         *(Ey(nx,nymax-1)+Eybb(nx,nymax))&
            !         +2.d0/(vcfact*dt+1.d0)*(Eyb(nx,nymax)+Eyb(nx,nymax-1))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0)) &
            !         *(Eyb(nx+1,nymax)-2.d0*Eyb(nx,nymax)&
            !         +Eyb(nx-1,nymax)+Eyb(nx+1,nymax-1)&
            !         -2.d0*Eyb(nx,nymax-1)+Eyb(nx-1,nymax-1))

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
          ex(:,0) = 0.d0
          ex(:,nymax) = 0.d0
          ey(0,:) = 0.d0
          ey(nxmax,:) = 0.d0
          ez(:,0) = 0.d0
          ez(:,nymax) = 0.d0
          ez(0,:) = 0.d0
          ez(nxmax,:) = 0.d0
      ENDIF
      ex(0,0)=0.d0
      ex(0,nymax)=0.d0
      ex(nxmax,0)=0.d0
      ex(nxmax,nymax)=0.d0
      ey(0,0)=0.d0
      ey(0,nymax)=0.d0
      ey(nxmax,0)=0.d0
      ey(nxmax,nymax)=0.d0
      ez(0,0)=0.d0
      ez(0,nymax)=0.d0
      ez(nxmax,0)=0.d0
      ez(nxmax,nymax)=0.d0
           IF(model_boundary .eq. 2) THEN
          !   dl=(vcfact*dt-sqrt(2.d0))/(vcfact*dt+sqrt(2.d0))
          !   dm=2.d0*sqrt(2.d0)/(vcfact*dt+sqrt(2.d0))
          !   dn=4.d0*vcfact**2*dt**2/(vcfact*dt*sqrt(2.d0)+2.d0)
          !    Ex(0,0)=-Exbb(1,1)+dl*(Ex(1,1)+Exbb(0,0))+dm*(Exb(1,1)+Exb(0,0))&
          !           +dn*(Exb(0,1)-Exb(1,1)-Exb(0,0)+Exb(1,0))
          !    Ey(0,0)=-Eybb(1,1)+dl*(Ey(1,1)+Eybb(0,0))+dm*(Eyb(1,1)+Eyb(0,0))&
          !           +dn*(Eyb(0,1)-Eyb(1,1)-Eyb(0,0)+Eyb(1,0))
          !    Ez(0,0)=-Ezbb(1,1)+dl*(Ez(1,1)+Ezbb(0,0))+dm*(Ezb(1,1)+Ezb(0,0))&
          !          +dn*(Ezb(0,1)-Ezb(1,1)-Ezb(0,0)+Ezb(1,0))
          !    Ex(nxmax,0)=-Exbb(nxmax-1,1)+dl*(Ex(nxmax-1,1)+Exbb(nxmax,0))+dm*(Exb(nxmax-1,1)+Exb(nxmax,0))&
          !                   +dn*(Exb(nxmax,1)-Exb(nxmax-1,1)-Exb(nxmax,0)+Exb(nxmax-1,0))
          !    Ey(nxmax,0)=-Eybb(nxmax-1,1)+dl*(Ey(nxmax-1,1)+Eybb(nxmax,0))+dm*(Eyb(nxmax-1,1)+Eyb(nxmax,0))&
          !                  +dn*(Eyb(nxmax,1)-Eyb(nxmax-1,1)-Eyb(nxmax,0)+Eyb(nxmax-1,0))
          !   Ez(nxmax,0)=-Ezbb(nxmax-1,1)+dl*(Ez(nxmax-1,1)+Ezbb(nxmax,0))+dm*(Ezb(nxmax-1,1)+Ezb(nxmax,0))&
          !                 +dn*(Ezb(nxmax,1)-Ezb(nxmax-1,1)-Ezb(nxmax,0)+Ezb(nxmax-1,0))
          !    Ex(0,nymax)=-Exbb(1,nymax-1)+dl*(Ex(1,nymax-1)+Exbb(0,nymax))+dm*(Exb(1,nymax-1)+Exb(0,nymax))&
          !                 +dn*(Exb(0,nymax-1)-Exb(1,nymax-1)-Exb(0,nymax)+Exb(1,nymax))
          !    Ey(0,nymax)=-Eybb(1,nymax-1)+dl*(Ey(1,nymax-1)+Eybb(0,nymax))+dm*(Eyb(1,nymax-1)+Eyb(0,nymax))&
          !                +dn*(Eyb(0,nymax-1)-Eyb(1,nymax-1)-Eyb(0,nymax)+Eyb(1,nymax))
          !   Ez(0,nymax)=-Ezbb(1,nymax-1)+dl*(Ez(1,nymax-1)+Ezbb(0,nymax))+dm*(Ezb(1,nymax-1)+Ezb(0,nymax))&
          !                +dn*(Ezb(0,nymax-1)-Ezb(1,nymax-1)-Ezb(0,nymax)+Ezb(1,nymax))
          !
          !    Ex(nxmax,nymax)=-Exbb(nxmax-1,nymax-1)&
          !                   +dl*(Ex(nxmax-1,nymax-1)+Exbb(nxmax,nymax))&
          !                   +dm*(Exb(nxmax-1,nymax-1)+Exb(nxmax,nymax))&
          !                   +dn*(Exb(nxmax,nymax-1)-Exb(nxmax-1,nymax-1)&
          !                   -Exb(nxmax,nymax)+Exb(nxmax-1,nymax))
          !
          !    Ey(nxmax,nymax)=-Eybb(nxmax-1,nymax-1)&
          !                    +dl*(Ey(nxmax-1,nymax-1)+Eybb(nxmax,nymax))&
          !                    +dm*(Eyb(nxmax-1,nymax-1)+Eyb(nxmax,nymax))&
          !                    +dn*(Eyb(nxmax,nymax-1)-Eyb(nxmax-1,nymax-1)&
          !                    -Eyb(nxmax,nymax)+Eyb(nxmax-1,nymax))
          !
          !    Ez(nxmax,nymax)=-Ezbb(nxmax-1,nymax-1)&
          !                  +dl*(Ez(nxmax-1,nymax-1)+Ezbb(nxmax,nymax))&
          !                  +dm*(Ezb(nxmax-1,nymax-1)+Ezb(nxmax,nymax))&
          !                  +dn*(Ezb(nxmax,nymax-1)-Ezb(nxmax-1,nymax-1)&
          !                  -Ezb(nxmax,nymax)+Ezb(nxmax-1,nymax))
           ENDIF
ENDIF
      IF(model_boundary .eq. 0) THEN
       IF(MOD(model_push/2,2).EQ.0) THEN
          IF(MOD(model_push,2).EQ.0) THEN
             ex(0:nxmax,0:nymax) = 0.D0
             ey(0:nxmax,0:nymax) = 0.D0
             ez(0:nxmax,0:nymax) = 0.D0
          ELSE
             ex(0:nxmax,0:nymax) = esx(0:nxmax,0:nymax)
             ey(0:nxmax,0:nymax) = esy(0:nxmax,0:nymax)
             ez(0:nxmax,0:nymax) = esz(0:nxmax,0:nymax)
          END IF
       ELSE
          IF(MOD(model_push,2).EQ.0) THEN
             ex(0:nxmax,0:nymax) = emx(0:nxmax,0:nymax)
             ey(0:nxmax,0:nymax) = emy(0:nxmax,0:nymax)
             ez(0:nxmax,0:nymax) = emz(0:nxmax,0:nymax)
          ELSE
         ex(0:nxmax,0:nymax) = esx(0:nxmax,0:nymax) + emx(0:nxmax,0:nymax)
           ey(0:nxmax,0:nymax) = esy(0:nxmax,0:nymax) + emy(0:nxmax,0:nymax)
           ez(0:nxmax,0:nymax) = esz(0:nxmax,0:nymax) + emz(0:nxmax,0:nymax)
          END IF
       END IF
     ENDIF
      IF(model_boundary.ne.0) THEN
      SELECT CASE(model_wg)
      CASE(0)
         yc=0.5d0*(ymin_wg+ymax_wg)
         ylen=(ymax_wg-ymin_wg)
         IF(ylen .NE. 0) dph=ph_wg/ylen
         amp_start = time * vcfact/100.d0
         IF(amp_start .GE. 1.d0) amp_start=1.0d0
         DO ny=1,nymax
            y=DBLE(ny)
            IF(y.GE.ymin_wg.AND.y.LE.ymax_wg) THEN
               factor=EXP(-12.D0*(y-yc)**2/(ylen)**2)
               Ey(1,ny)=amp_wg*amp_start &
                    *factor*COS(rot_wg*pi/180.D0) &
                    *SIN(omega*time-pi*dph*(y-ymin_wg)/180.D0)
               Ez(1,ny)=amp_wg*amp_start &
                    *factor*SIN(rot_wg*pi/180.D0) &
                    *SIN(omega*time-pi*dph*(y-ymin_wg)/180.D0)
            END IF
         END DO
      END SELECT
    ENDIF

  END SUBROUTINE efield

  !***********************************************************************
  SUBROUTINE bfield(nxmax,nymax,dt,Ax,Ay,Az,Axb,Ayb,Azb, &
       ex,ey,ez,bx,by,bz,bxb,byb,bzb,bxbb,bybb,bzbb,bxbg,bybg,bzbg,bb,vcfact, &
       model_push,model_boundary)
  !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nymax) :: bxnab,bynab,bznab
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: ex,ey,ez,bx,by,bz,bxb,byb,bzb, &
                                           bxbb,bybb,bzbb,bxbg,bybg,bzbg,bb
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: Ax,Ay,Az,Axb,Ayb,Azb
    INTEGER :: nxmax, nymax, nx, ny, nxp, nyp, nxm, nym
    INTEGER:: model_push, model_boundary
    REAL(8):: bxx,byy,bzz,vcfact,dt,dl,dm,dn
    IF(model_boundary .EQ. 0) THEN
      !$omp parallel do private(nx,ny,nxm,nym,nxp,nyp)
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

              ! bx(nx,ny) = 0.5d0 * (Az(nx,nyp) + Azb(nx,nyp) &
              !      - Az(nx,ny) - Azb(nx,ny))
              ! by(nx,ny) = - 0.5d0 * (Az(nxp,ny) + Azb(nxp,ny) &
              !      - Az(nx,ny) - Azb(nx,ny))
              ! bz(nx,ny) = 0.5d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
              !      - Ay(nx,ny) - Ayb(nx,ny) &
              !      - (Ax(nx,nyp) + Axb(nx,nyp) &
              !      - Ax(nx,ny) - Axb(nx,ny)))

            bx(nx,ny)=dt*(-ez(nx,nyp)+ez(nx,ny))+bxb(nx,ny)
            by(nx,ny)=dt*(ez(nxp,ny)-ez(nx,ny))+byb(nx,ny)
            bz(nx,ny)=dt*(-ey(nxp,ny)+ey(nx,ny)+ex(nx,nyp)-ex(nx,ny))+bzb(nx,ny)
            bxbb(nx,ny)=bxb(nx,ny)
            bybb(nx,ny)=byb(nx,ny)
            bzbb(nx,ny)=bzb(nx,ny)
            bxx=bx(nx,ny)
            byy=by(nx,ny)
            bzz=bz(nx,ny)
            bx(nx,ny)=0.5d0*(bx(nx,ny)+bxb(nx,ny))
            by(nx,ny)=0.5d0*(by(nx,ny)+byb(nx,ny))
            bz(nx,ny)=0.5d0*(bz(nx,ny)+bzb(nx,ny))
            bxb(nx,ny)=bxx
            byb(nx,ny)=byy
            bzb(nx,ny)=bzz
          END DO
       END DO
       !$omp end parallel do
    ELSEIF(model_boundary .NE. 0) THEN
      !$omp parallel do private(nx,ny,nxm,nym,nxp,nyp)
       DO ny = 1, nymax-1
          DO nx = 1, nxmax-1
             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             IF( nx .EQ. 0  )    nxm = 0
             IF( nx .EQ. nxmax ) nxp = nxmax
             IF( ny .EQ. 0  )    nym = 0
             IF( ny .EQ. nymax ) nyp = nymax

                ! bx(nx,ny) =   0.5d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                !      - Az(nx,ny) - Azb(nx,ny))
                ! by(nx,ny) = - 0.5d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                !      - Az(nx,ny) - Azb(nx,ny))
                ! bz(nx,ny) =   0.5d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                !      - Ay(nx,ny) - Ayb(nx,ny) &
                !      - (Ax(nx,nyp) + Axb(nx,nyp) &
                !      - Ax(nx,ny) - Axb(nx,ny)))
            bx(nx,ny)=dt*(-ez(nx,nyp)+ez(nx,ny))+bxb(nx,ny)
            by(nx,ny)=dt*(ez(nxp,ny)-ez(nx,ny))+byb(nx,ny)
            bz(nx,ny)=dt*(-ey(nxp,ny)+ey(nx,ny)+ex(nx,nyp)-ex(nx,ny))+bzb(nx,ny)
          END DO
       END DO
       !$omp end parallel do
        DO ny = 1, nymax-1
          nyp = ny + 1
          IF(model_boundary .eq. 1) then  !reflective
            !  bx(nxmax,ny)=dt*(-ez(nxmax,nyp)+ez(nxmax,ny))+bxb(nxmax,ny)
            !  by(nxmax,ny)=byb(nxmax,ny)
            !  bz(nxmax,ny)=dt*(ex(nxmax,nyp)-ex(nxmax,ny))+bzb(nxmax,ny)
            !  bx(0,ny)=bx(1,ny)
            !  by(0,ny)=by(1,ny)
            !  bz(0,ny)=bz(1,ny)
          ELSE IF(model_boundary .eq. 2) then !Mur's abosorbing boundary condition
            !  bx(nxmax,ny)=-bxbb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !          *(bx(nxmax-1,ny)+bxbb(nxmax,ny)) &
            !          +2.d0/(vcfact*dt+1.d0)*(bxb(nxmax,ny)+bxb(nxmax-1,ny))&
            !          +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !          *(bxb(nxmax,ny+1)-2.d0*bxb(nxmax,ny)&
            !          +bxb(nxmax,ny-1)+bxb(nxmax-1,ny+1)&
            !          -2.d0*bxb(nxmax-1,ny)+bxb(nxmax-1,ny-1))
            !
            !  by(nxmax,ny)=-bybb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !          *(by(nxmax-1,ny)+bybb(nxmax,ny)) &
            !          +2.d0/(vcfact*dt+1.d0)*(byb(nxmax,ny)+byb(nxmax-1,ny))&
            !          +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !          *(byb(nxmax,ny+1)-2.d0*byb(nxmax,ny)&
            !          +byb(nxmax,ny-1)+byb(nxmax-1,ny+1)&
            !          -2.d0*byb(nxmax-1,ny)+byb(nxmax-1,ny-1))
            !
            ! bz(nxmax,ny)=-bzbb(nxmax-1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !         *(bz(nxmax-1,ny)+bzbb(nxmax,ny)) &
            !         +2.d0/(vcfact*dt+1.d0)*(bzb(nxmax,ny)+bzb(nxmax-1,ny))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !         *(bzb(nxmax,ny+1)-2.d0*bzb(nxmax,ny)&
            !         +bzb(nxmax,ny-1)+bzb(nxmax-1,ny+1)&
            !         -2.d0*bzb(nxmax-1,ny)+bzb(nxmax-1,ny-1))
            !
            !   bx(0,ny)=-bxbb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !         *(bx(1,ny)+bxbb(0,ny)) &
            !         +2.d0/(vcfact*dt+1.d0)*(bxb(0,ny)+bxb(1,ny))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !         *(bxb(0,ny+1)-2.d0*bxb(0,ny)&
            !         +bxb(0,ny-1)+bxb(1,ny+1)&
            !         -2.d0*bxb(1,ny)+bxb(1,ny-1))
            !
            !   by(0,ny)=-bybb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !         *(by(1,ny)+bybb(0,ny)) &
            !         +2.d0/(vcfact*dt+1.d0)*(byb(0,ny)+byb(1,ny))&
            !         +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !         *(byb(0,ny+1)-2.d0*byb(0,ny)&
            !         +byb(0,ny-1)+byb(1,ny+1)&
            !         -2.d0*byb(1,ny)+byb(1,ny-1))
            !
            !  bz(0,ny)=-bzbb(1,ny)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
            !        *(bz(1,ny)+bzbb(0,ny)) &
            !        +2.d0/(vcfact*dt+1.d0)*(bzb(0,ny)+bzb(1,ny))&
            !        +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
            !        *(bzb(0,ny+1)-2.d0*bzb(0,ny)&
            !        +bzb(0,ny-1)+bzb(1,ny+1)&
            !        -2.d0*bzb(1,ny)+bzb(1,ny-1))
           ENDIF
         ENDDO
         DO nx = 1, nxmax-1
           nxp = nx + 1
           IF(model_boundary .eq. 1) then
              ! bx(nx,nymax)=bxb(nx,nymax)
              ! by(nx,nymax)=dt*(ez(nxp,nymax)-ez(nx,nymax))+byb(nx,nymax)
              ! bz(nx,nymax)=dt*(-ey(nxp,nymax)+ey(nx,nymax))+bzb(nx,nymax)
              ! bx(nx,0)=bx(nx,1)
              ! by(nx,0)=by(nx,1)
              ! bz(nx,0)=bz(nx,1)
           ELSEIF(model_boundary .eq. 2) then
          !    bx(nx,0)=-bxbb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !            *(bx(nx,1)+bxbb(nx,0)) &
          !            +2.d0/(vcfact*dt+1.d0)*(bxb(nx,0)+bxb(nx,1))&
          !            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !            *(bxb(nx+1,0)-2.d0*bxb(nx,0)+bxb(nx-1,0)+bxb(nx+1,1)&
          !             -2.d0*bxb(nx,1)+bxb(nx-1,1))
          !
          !    by(nx,0)=-bybb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !            *(by(nx,1)+bybb(nx,0))&
          !            +2.d0/(vcfact*dt+1.d0)*(byb(nx,0)+byb(nx,1))&
          !            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !            *(byb(nx+1,0)-2.d0*byb(nx,0)+byb(nx-1,0)+byb(nx+1,1)&
          !             -2.d0*byb(nx,1)+byb(nx-1,1))
          !
          !   bz(nx,0)=-bzbb(nx,1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !           *(bz(nx,1)+bzbb(nx,0))&
          !           +2.d0/(vcfact*dt+1.d0)*(bzb(nx,0)+bzb(nx,1))&
          !           +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !           *(bzb(nx+1,0)-2.d0*bzb(nx,0)+bzb(nx-1,0)+bzb(nx+1,1)&
          !            -2.d0*bzb(nx,1)+bzb(nx-1,1))
          !
          !    bx(nx,nymax)=-bxbb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !             *(bx(nx,nymax-1)+bxbb(nx,nymax)) &
          !            +2.d0/(vcfact*dt+1.d0)*(bxb(nx,nymax)+bxb(nx,nymax-1))&
          !            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !            *(bxb(nx+1,nymax)-2.d0*bxb(nx,nymax) &
          !            +bxb(nx-1,nymax)+bxb(nx+1,nymax-1)&
          !            -2.d0*bxb(nx,nymax-1)+bxb(nx-1,nymax-1))
          !
          !    by(nx,nymax)=-bybb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !            *(by(nx,nymax-1)+bybb(nx,nymax))&
          !            +2.d0/(vcfact*dt+1.d0)*(byb(nx,nymax)+byb(nx,nymax-1))&
          !            +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0)) &
          !            *(byb(nx+1,nymax)-2.d0*byb(nx,nymax)&
          !            +byb(nx-1,nymax)+byb(nx+1,nymax-1)&
          !            -2.d0*byb(nx,nymax-1)+byb(nx-1,nymax-1))
          !
          !   bz(nx,nymax)=-bzbb(nx,nymax-1)+(vcfact*dt-1.d0)/(vcfact*dt+1.d0)&
          !           *(bz(nx,nymax-1)+bzbb(nx,nymax))&
          !           +2.d0/(vcfact*dt+1.d0)*(bzb(nx,nymax)+bzb(nx,nymax-1))&
          !           +(vcfact*dt)**2/(2.d0*(vcfact*dt+1.d0))&
          !           *(bzb(nx+1,nymax)-2.d0*bzb(nx,nymax)&
          !           +bzb(nx-1,nymax)+bzb(nx+1,nymax-1)&
          !           -2.d0*bzb(nx,nymax-1)+bzb(nx-1,nymax-1))
            ENDIF
        ENDDO
        IF(model_boundary .eq. 2) THEN
          dl=(vcfact*dt-sqrt(2.d0))/(vcfact*dt+sqrt(2.d0))
          dm=2.d0*sqrt(2.d0)/(vcfact*dt+sqrt(2.d0))
          dn=4.d0*vcfact**2*dt**2/(vcfact*dt*sqrt(2.d0)+2.d0)
          ! Bx(0,0)=-Bxbb(1,1)+dl*(Bx(1,1)+Bxbb(0,0))+dm*(Bxb(1,1)+Bxb(0,0))&
          !        +dn*(Bxb(0,1)-Bxb(1,1)-Bxb(0,0)+Bxb(1,0))
          ! By(0,0)=-Bybb(1,1)+dl*(By(1,1)+Bybb(0,0))+dm*(Byb(1,1)+Byb(0,0))&
          !        +dn*(Byb(0,1)-Byb(1,1)-Byb(0,0)+Byb(1,0))
          !Bz(0,0)=-Bzbb(1,1)+dl*(Bz(1,1)+Bzbb(0,0))+dm*(Bzb(1,1)+Bzb(0,0))&
          !       +dn*(Bzb(0,1)-Bzb(1,1)-Bzb(0,0)+Bzb(1,0))
          ! Bx(nxmax,0)=-Bxbb(nxmax-1,1)+dl*(Bx(nxmax-1,1)+Bxbb(nxmax,0))+dm*(Bxb(nxmax-1,1)+Bxb(nxmax,0))&
          !               +dn*(Bxb(nxmax,1)-Bxb(nxmax-1,1)-Bxb(nxmax,0)+Bxb(nxmax-1,0))
          ! By(nxmax,0)=-Bybb(nxmax-1,1)+dl*(By(nxmax-1,1)+Bybb(nxmax,0))+dm*(Byb(nxmax-1,1)+Byb(nxmax,0))&
          !               +dn*(Byb(nxmax,1)-Byb(nxmax-1,1)-Byb(nxmax,0)+Byb(nxmax-1,0))
          !Bz(nxmax,0)=-Bzbb(nxmax-1,1)+dl*(Bz(nxmax-1,1)+Bzbb(nxmax,0))+dm*(Bzb(nxmax-1,1)+Bzb(nxmax,0))&
          !              +dn*(Bzb(nxmax,1)-Bzb(nxmax-1,1)-Bzb(nxmax,0)+Bzb(nxmax-1,0))
          ! Bx(0,nymax)=-Bxbb(1,nymax-1)+dl*(Bx(1,nymax-1)+Bxbb(0,nymax))+dm*(Bxb(1,nymax-1)+Bxb(0,nymax))&
          !              +dn*(Bxb(0,nymax-1)-Bxb(1,nymax-1)-Bxb(0,nymax)+Bxb(1,nymax))
          ! By(0,nymax)=-Bybb(1,nymax-1)+dl*(By(1,nymax-1)+Bybb(0,nymax))+dm*(Byb(1,nymax-1)+Byb(0,nymax))&
          !              +dn*(Byb(0,nymax-1)-Byb(1,nymax-1)-Byb(0,nymax)+Byb(1,nymax))
          !Bz(0,nymax)=-Bzbb(1,nymax-1)+dl*(Bz(1,nymax-1)+Bzbb(0,nymax))+dm*(Bzb(1,nymax-1)+Bzb(0,nymax))&
          !             +dn*(Bzb(0,nymax-1)-Bzb(1,nymax-1)-Bzb(0,nymax)+Bzb(1,nymax))
          ! Bx(nxmax,nymax)=-Bxbb(nxmax-1,nymax-1)&
          !                +dl*(Bx(nxmax-1,nymax-1)+Bxbb(nxmax,nymax))&
          !                +dm*(Bxb(nxmax-1,nymax-1)+Bxb(nxmax,nymax))&
          !                +dn*(Bxb(nxmax,nymax-1)-Bxb(nxmax-1,nymax-1)&
          !                -Bxb(nxmax,nymax)+Bxb(nxmax-1,nymax))
          ! By(nxmax,nymax)=-Bybb(nxmax-1,nymax-1)+dl*(By(nxmax-1,nymax-1)&
          !                +Bybb(nxmax,nymax))+dm*(Byb(nxmax-1,nymax-1)&
          !                +Byb(nxmax,nymax))&
          !                +dn*(Byb(nxmax,nymax-1)-Byb(nxmax-1,nymax-1)&
          !                -Byb(nxmax,nymax)+Byb(nxmax-1,nymax))
        !  Bz(nxmax,nymax)=-Bzbb(nxmax-1,nymax-1)+dl*(Bz(nxmax-1,nymax-1)&
        !                 +Bzbb(nxmax,nymax))+dm*(Bzb(nxmax-1,nymax-1)&
        !                 +Bzb(nxmax,nymax))&
        !                 +dn*(Bzb(nxmax,nymax-1)-Bzb(nxmax-1,nymax-1)&
        !                 -Bzb(nxmax,nymax)+Bzb(nxmax-1,nymax))
        ENDIF
          bx(0,:) = 0.d0
          bx(nxmax,:) = 0.d0
          bx(:,0) = 0.d0
          bx(:,nymax) = 0.d0
          by(0,:) = 0.d0
          by(nxmax,:) = 0.d0
          by(:,0) = 0.d0
          by(:,nymax) = 0.d0
          bz(0,:) = 0.d0
          bz(nxmax,:) = 0.d0
          bz(:,0) = 0.d0
          bz(:,nymax) = 0.d0
        ! bxb(0,:) = 0.d0
        ! bxb(nxmax,:) = 0.d0
        ! bxb(:,nymax) = 0.d0
        ! byb(:,0) = 0.d0
        ! byb(:,nymax) = 0.d0
        ! byb(nxmax,:) = 0.d0
        ! bzb(nxmax,:) = 0.d0
        ! bzb(:,nymax) = 0.d0

        DO nx = 0,nxmax
          DO ny = 0,nymax
        bxbb(nx,ny)=bxb(nx,ny)
        bybb(nx,ny)=byb(nx,ny)
        bzbb(nx,ny)=bzb(nx,ny)
        bxx=bx(nx,ny)
        byy=by(nx,ny)
        bzz=bz(nx,ny)
        bx(nx,ny)=0.5d0*(bx(nx,ny)+bxb(nx,ny))
        by(nx,ny)=0.5d0*(by(nx,ny)+byb(nx,ny))
        bz(nx,ny)=0.5d0*(bz(nx,ny)+bzb(nx,ny))
        bxb(nx,ny)=bxx
        byb(nx,ny)=byy
        bzb(nx,ny)=bzz
      ENDDO
    ENDDO
      ENDIF
    IF(MOD(model_push/8,2).EQ.0) THEN
       IF(MOD(model_push/4,2).EQ.0) THEN
          bx(0:nxmax,0:nymax) = 0.D0
          by(0:nxmax,0:nymax) = 0.D0
          bz(0:nxmax,0:nymax) = 0.D0
       ELSE
          bx(0:nxmax,0:nymax) = bxbg(0:nxmax,0:nymax)
          by(0:nxmax,0:nymax) = bybg(0:nxmax,0:nymax)
          bz(0:nxmax,0:nymax) = bzbg(0:nxmax,0:nymax)
       END IF
    ELSE
       IF(MOD(model_push/4,2).EQ.0) THEN
          CONTINUE
       ELSE
          bx(0:nxmax,0:nymax) = bx(0:nxmax,0:nymax) + bxbg(0:nxmax,0:nymax)
          by(0:nxmax,0:nymax) = by(0:nxmax,0:nymax) + bybg(0:nxmax,0:nymax)
          bz(0:nxmax,0:nymax) = bz(0:nxmax,0:nymax) + bzbg(0:nxmax,0:nymax)
       END IF
    END IF

    bb(0:nxmax,0:nymax) = SQRT(bx(0:nxmax,0:nymax)**2 &
         +by(0:nxmax,0:nymax)**2 &
         +bz(0:nxmax,0:nymax)**2)


  END SUBROUTINE bfield

  !***********************************************************************
  SUBROUTINE kine(npmax,vx,vy,vz,akin,mass,vcfact)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(npmax) :: vx, vy, vz
    REAL(8) :: akin, mass, vcfact, gamma
    INTEGER(4) :: npmax, np
    akin = 0.d0
    DO np = 1, npmax
      gamma = 1.d0/sqrt(1.d0 - (vx(np)**2+vy(np)**2+vz(np)**2)/vcfact**2)
       akin = akin + vcfact ** 2 * (gamma - 1)
    END DO

    IF(npmax.EQ.0) THEN
       akin=0.D0
    ELSE
       akin = akin * mass  / DBLE(npmax)
    END IF
  END SUBROUTINE kine

  !***********************************************************************
  SUBROUTINE pote(nxmax,nymax,ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg,vcfact, &
       apote,apotm)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg
    REAL(8) :: apote,apotm,vcfact
    INTEGER(4) :: nxmax, nymax, nx, ny

    apote = 0.d0
    apotm = 0.d0
    !$omp parallel do reduction(+:apote,apotm)
    DO ny = 1, nymax-1
       DO nx = 1, nxmax-1
          apote = apote + (ex(nx,ny)**2 + ey(nx,ny)**2 + ez(nx,ny)**2)
          apotm = apotm + ((bx(nx,ny)-bxbg(nx,ny))**2 &
               + (by(nx,ny)-bybg(nx,ny))**2 &
               + (bz(nx,ny)-bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    !$omp parallel do reduction(+:apote,apotm)
    DO nx = 0, nxmax,nxmax
       DO ny = 1, nymax-1
          apote = apote + 0.5D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
          apotm = apotm + 0.5D0*((bx(nx,ny)-bxbg(nx,ny))**2 &
               + (by(nx,ny)-bybg(nx,ny))**2 &
               + (bz(nx,ny)-bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    !$omp parallel do reduction(+:apote,apotm)
    DO ny = 0, nymax,nymax
       DO nx = 1, nxmax-1
          apote = apote + 0.5D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
          apotm = apotm + 0.5D0*((bx(nx,ny)-bxbg(nx,ny))**2 &
               + (by(nx,ny)-bybg(nx,ny))**2 &
               + (bz(nx,ny)-bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    !$omp parallel do reduction(+:apote,apotm)
    DO ny = 0, nymax,nymax
       DO nx = 0, nxmax,nxmax
          apote = apote + 0.25D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
          apotm = apotm + 0.25D0*((bx(nx,ny)-bxbg(nx,ny))**2 &
               + (by(nx,ny)-bybg(nx,ny))**2 &
               + (bz(nx,ny)-bzbg(nx,ny))**2)
       END DO
    END DO
    !$omp end parallel do
    apote = 0.5D0 * apote / (DBLE(nxmax)*DBLE(nymax))
    apotm = 0.5D0 * vcfact**2 * apotm / (DBLE(nxmax)*DBLE(nymax))
  END SUBROUTINE pote

END MODULE picsub
