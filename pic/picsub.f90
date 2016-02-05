!  ***** TASK/PIC PREPARATION *****

MODULE picsub
  PRIVATE
  PUBLIC poisson_f,poisson_m,efield,bfield,kine,pote

CONTAINS

  !***********************************************************************
  SUBROUTINE poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
       rho,phi,rhof,phif,awk,afwk,cform,ipssn)
    !***********************************************************************
    IMPLICIT NONE
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

    ELSE

       !----- solve poisson equation
       DO ny = 1, nymax
          DO nx = 1, nxmaxh1
             phif(nx,ny) = cform(nx,ny) * rhof(nx,ny)
          END DO
       END DO

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
       DO ny = 1, nymax
          DO nx = 1, nxmax
             awk(nx,ny) = a(nx,ny)
          END DO
       END DO

       CALL dfftw_execute(plan1)

       DO ny = 1, nymax
          DO nx = 1, nxmaxh1
             af(nx,ny) = afwk(nx,ny) / ( alx * aly )
          END DO
       END DO

    ELSE

       !----- inverse fourier transform

       DO ny = 1, nymax
          DO nx = 1, nxmaxh1
             afwk(nx,ny) = af(nx,ny)
          END DO
       END DO

       CALL dfftw_execute(plan2)

       DO ny = 1, nymax
          DO nx = 1, nxmax
             a(nx,ny) = awk(nx,ny)
          END DO
       END DO

       DO ny = 1, nymax
          a(nxmax1,ny) = a(1,ny)
       END DO

       DO nx = 1, nxmax1
          a(nx,nymax1) = a(nx,1)
       END DO

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

    CALL mtx_setup(imax,istart,iend)
    irange=iend-istart+1
    status=1
    ALLOCATE(x(irange))

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
       DO i=istart,iend
          nx=MOD(i-1,isize)+1
          ny=(i-1)/isize+1
          CALL mtx_set_source(i,-rho(nx+1,ny+1))
       ENDDO
    ELSE
       DO i=istart,iend
          ny=MOD(i-1,isize)+1
          nx=(i-1)/isize+1
          CALL mtx_set_source(i,-rho(nx+1,ny+1))
       ENDDO
    END IF
    CALL mtx_solve(model_matrix0,tolerance_matrix,its, &
         methodKSP=model_matrix1,methodPC=model_matrix2)
    !      IF((its.ne.0).and.(nrank.eq.0)) write(6,*) 'mtx_solve: its=',its
    CALL mtx_get_vector(x)
    IF(mode.EQ.0) THEN
       DO i=istart,iend
          nx=MOD(i-1,isize)+1
          ny=(i-1)/isize+1
          phi(nx+1,ny+1)=x(i-istart+1)
       ENDDO
    ELSE
       DO i=istart,iend
          ny=MOD(i-1,isize)+1
          nx=(i-1)/isize+1
          phi(nx+1,ny+1)=x(i-istart+1)
       ENDDO
    END IF
    DEALLOCATE(x)
    status=2
    CALL mtx_cleanup
    IF(model_boundary .EQ. 2) THEN !damping phi in absorbing boundary
       inv = 1.0d0 / dlen
       ilen = int(dlen)
       DO ny = 1,nymax
          DO nx = nxmax-ilen,nxmax
             xd=DBLE(nx)
             xdmax=DBLE(nxmax)
             phi(nx,ny) =phi(nx,ny)*(-1.0d0*inv**2*xd**2 &
                                     +2.0d0*inv**2*(xdmax-dlen)*xd&
                                     +1.0d0-1.0d0*inv**2*(xdmax-dlen)**2)
          ENDDO
       ENDDO
       DO nx = 1,nxmax-ilen
          DO ny = nymax-ilen,nymax
             yd=DBLE(ny)
             ydmax=DBLE(nymax)
             phi(nx,ny) =phi(nx,ny)*(-1.0d0*inv**2*yd**2 &
                                     +2.0d0*inv**2*(ydmax-dlen)*yd&
                                     +1.0d0-1.0d0*inv**2*(ydmax-dlen)**2)
          ENDDO
       ENDDO
       DO nx = 1, nxmax-ilen
          DO ny = 1, ilen
             yd=DBLE(ny)
             phi(nx,ny) = phi(nx,ny)*(-1.0d0*inv**2*yd**2+2.0d0*inv*yd)
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE poisson_m

  !***********************************************************************
  SUBROUTINE efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
       ex,ey,ez,esx,esy,esz,emx,emy,emz, &
       model_push,model_boundary)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) ::  &
         phi,Ax,Ay,Az,Axb,Ayb,Azb,ex,ey,ez,esx,esy,esz,emx,emy,emz
    REAL(8):: dt
    INTEGER :: nxmax, nymax, nx, ny, nxm, nxp, nym, nyp
    INTEGER:: model_push, model_boundary

    IF(model_boundary .EQ. 0) THEN
       DO nx = 0, nymax
          DO ny = 0, nxmax

             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             IF( nx .EQ. 0  )    nxm = nxmax - 1
             IF( nx .EQ. nxmax ) nxp = 1
             IF( ny .EQ. 0  )    nym = nymax - 1
             IF( ny .EQ. nymax ) nyp = 1
             esx(nx,ny) = 0.5d0 * ( phi(nxm,ny) - phi(nxp,ny))
             esy(nx,ny) = 0.5d0 * ( phi(nx,nym) - phi(nx,nyp))
             esz(nx,ny) = 0.d0
             emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
             emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
             emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt

          END DO
       END DO
    ELSE IF (model_boundary .NE. 0) THEN
       DO nx = 0, nymax
          DO ny = 0, nxmax

             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             IF( nx .EQ. 0  )    nxm = 0
             IF( nx .EQ. nxmax ) nxp = nxmax
             IF( ny .EQ. 0  )    nym = 0
             IF( ny .EQ. nymax ) nyp = nymax

              IF(nx .EQ. 0 .OR. nx .EQ. nxmax) THEN
                 esx(nx,ny) = phi(nxm,ny) - phi(nxp,ny)
              ELSE
                 esx(nx,ny) = 0.5d0 * ( phi(nxm,ny) - phi(nxp,ny))
              ENDIF

             IF(ny .EQ. 0 .OR. ny .EQ. nymax) THEN
                 esy(nx,ny) = phi(nx,nym) - phi(nx,nyp)
             ELSE
                 esy(nx,ny) = 0.5d0 * ( phi(nx,nym) - phi(nx,nyp))
             END IF

             esz(nx,ny) = 0.d0
             emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
             emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
             emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt

          END DO
       END DO
       esx(:,0) = 0.d0
       esx(:,nymax) = 0.d0
       esy(0,:) = 0.d0
       esy(nxmax,:) = 0.d0
       emx(:,0) = 0.d0
       emx(:,nymax) = 0.d0
       emy(0,:) = 0.d0
       emy(nxmax,:) = 0.d0
       emz(:,0) = 0.d0
       emz(:,nymax) = 0.d0
       emz(0,:) = 0.d0
       emz(nxmax,:) = 0.d0

    END IF

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
    !  if(model_boundary .ne. 0 .and. nxmax .ge. 10 .and. nymax .ge. 10) then
    !   do ny = 10,nymax-10
    !   do nx = nxmax-10,nxmax
    !     Ex(nx,ny) = Ex(nx,ny) * (-0.01d0 * dble(nx) ** 2 &
    !                           + 2.0d0 * 0.01d0 * (dble(nxmax - 10)) * dble(nx) &
    !                           + 1.0d0 - 0.01d0 * (dble(nxmax - 10))**2)
    !     Ey(nx,ny) = Ey(nx,ny) * (-0.01d0 * dble(nx) ** 2 &
    !                           + 2.0d0 * 0.01d0 * (dble(nxmax - 10)) * dble(nx) &
    !                           + 1.0d0 - 0.01d0 * (dble(nxmax - 10))**2)
    !     Ez(nx,ny) = Ez(nx,ny) * (-0.01d0 * dble(nx) ** 2 &
    !                           + 2.0d0 * 0.01d0 * (dble(nxmax - 10)) * dble(nx) &
    !                           + 1.0d0 - 0.01d0 * (dble(nxmax - 10))**2)
    !   enddo
    !   enddo
    !   do nx = 0,nxmax-10
    !   do ny = nymax-10,nymax
    !     Ex(nx,ny) = Ex(nx,ny) * (-0.01d0 * dble(ny) ** 2 &
    !                           + 2.0d0 * 0.01d0 * (dble(nymax - 10)) * dble(ny) &
    !                           + 1.0d0 - 0.01d0 * (dble(nymax - 10))**2)
    !     Ey(nx,ny) = Ey(nx,ny) * (-0.01d0 * dble(ny) ** 2 &
    !                           + 2.0d0 * 0.01d0 * (dble(nymax - 10)) * dble(ny) &
    !                           + 1.0d0 - 0.01d0 * (dble(nymax - 10))**2)
    !     Ez(nx,ny) = Ez(nx,ny) * (-0.01d0 * dble(ny) ** 2 &
    !                           + 2.0d0 * 0.01d0 * (dble(nymax - 10)) * dble(ny) &
    !                           + 1.0d0 - 0.01d0 * (dble(nymax - 10))**2)
    !   enddo
    !   enddo
    !   do nx = 0, nxmax-10
    !   do ny = 0, 10
    !     Ex(nx,ny) = Ex(nx,ny) * (-0.01d0 * dble(ny) ** 2 &
    !                           + 2.0d0 * 0.1d0 * dble(ny))
    !     Ey(nx,ny) = Ey(nx,ny) * (-0.01d0 * dble(ny) ** 2 &
    !                           + 2.0d0 * 0.1d0 * dble(ny))
    !     Ez(nx,ny) = Ez(nx,ny) * (-0.01d0 * dble(ny) ** 2 &
    !                           + 2.0d0 * 0.1d0 * dble(ny))
    !   enddo
    !   enddo
    ! endif
  END SUBROUTINE efield

  !***********************************************************************
  SUBROUTINE bfield(nxmax,nymax,Ax,Ay,Az,Axb,Ayb,Azb, &
       bx,by,bz,bxbg,bybg,bzbg,bb, &
       model_push,model_boundary)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nymax) :: bxnab,bznab
    REAL(8), DIMENSION(0:nxmax) :: bynab
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: bx,by,bz,bxbg,bybg,bzbg,bb
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: Ax,Ay,Az,Axb,Ayb,Azb
    INTEGER :: nxmax, nymax, nx, ny, nxp, nyp, nxm, nym
    INTEGER:: model_push, model_boundary

    IF(model_boundary .EQ. 0) THEN
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

             bx(nx,ny) = 0.25d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                  - Az(nx,nym) - Azb(nx,nym))
             by(nx,ny) = - 0.25d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                  - Az(nxm,ny) - Azb(nxm,ny))
             bz(nx,ny) = 0.25d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                  - Ay(nxm,ny) - Ayb(nxm,ny) &
                  - (Ax(nx,nyp) + Axb(nx,nyp) &
                  - Ax(nx,nym) - Axb(nx,nym)))
          END DO
       END DO
    ELSEIF(model_boundary .NE. 0) THEN
       DO ny = 0, nymax
          DO nx = 0, nxmax
             nxm = nx - 1
             nxp = nx + 1
             nym = ny - 1
             nyp = ny + 1

             IF( nx .EQ. 0  )    nxm = 0
             IF( nx .EQ. nxmax ) nxp = nxmax
             IF( ny .EQ. 0  )    nym = 0
             IF( ny .EQ. nymax ) nyp = nymax

             IF(nx .EQ. 0 .OR. nx .EQ. nxmax .OR. &
                  ny .EQ. 0 .OR. ny .EQ. nymax) THEN

                bx(nx,ny) =   0.5d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                     - Az(nx,nym) - Azb(nx,nym))
                by(nx,ny) = - 0.5d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                     - Az(nxm,ny) - Azb(nxm,ny))
                bz(nx,ny) =   0.5d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                     - Ay(nxm,ny) - Ayb(nxm,ny) &
                     -(Ax(nx,nyp) + Axb(nx,nyp) &
                     - Ax(nx,nym) - Axb(nx,nym)))
             ELSE

                bx(nx,ny) =   0.25d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                     - Az(nx,nym) - Azb(nx,nym))
                by(nx,ny) = - 0.25d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                     - Az(nxm,ny) - Azb(nxm,ny))
                bz(nx,ny) =   0.25d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                     - Ay(nxm,ny) - Ayb(nxm,ny) &
                     - (Ax(nx,nyp) + Axb(nx,nyp) &
                     - Ax(nx,nym) - Axb(nx,nym)))
             END IF
          END DO
       END DO
       bx(0,:) = 0.d0
       bx(nxmax,:) = 0.d0
       by(:,0) = 0.d0
       by(:,nymax) = 0.d0
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
  SUBROUTINE kine(npmax,vx,vy,vz,akin,mass)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(npmax) :: vx, vy, vz
    REAL(8) :: akin, mass
    INTEGER(4) :: npmax, np
    akin = 0.d0
    DO np = 1, npmax
       akin = akin + vx(np)**2 + vy(np)**2 + vz(np)**2
    END DO

    IF(npmax.EQ.0) THEN
       akin=0.D0
    ELSE
       akin = 0.5 * akin * mass /DBLE(npmax)
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

    DO ny = 1, nymax-1
       DO nx = 1, nxmax-1
          apote = apote + (ex(nx,ny)**2 + ey(nx,ny)**2 + ez(nx,ny)**2)
          apotm = apotm + ((bx(nx,ny)-bxbg(nx,nx))**2 &
               + (by(nx,ny)-bybg(nx,nx))**2 &
               + (bz(nx,ny)-bzbg(nx,nx))**2)
       END DO
    END DO

    DO nx = 0, nxmax,nxmax
       DO ny = 1, nymax-1
          apote = apote + 0.5D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
          apotm = apotm + 0.5D0*((bx(nx,ny)-bxbg(nx,nx))**2 &
               + (by(nx,ny)-bybg(nx,nx))**2 &
               + (bz(nx,ny)-bzbg(nx,nx))**2)
       END DO
    END DO

    DO ny = 0, nymax,nymax
       DO nx = 1, nxmax-1
          apote = apote + 0.5D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
          apotm = apotm + 0.5D0*((bx(nx,ny)-bxbg(nx,nx))**2 &
               + (by(nx,ny)-bybg(nx,nx))**2 &
               + (bz(nx,ny)-bzbg(nx,nx))**2)
       END DO
    END DO

    DO ny = 0, nymax,nymax
       DO nx = 0, nxmax,nxmax
          apote = apote + 0.25D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
          apotm = apotm + 0.25D0*((bx(nx,ny)-bxbg(nx,nx))**2 &
               + (by(nx,ny)-bybg(nx,nx))**2 &
               + (bz(nx,ny)-bzbg(nx,nx))**2)
       END DO
    END DO
    apote = 0.5D0 * apote / (DBLE(nxmax)*DBLE(nymax))
    apotm = 0.5D0 * vcfact**2 * apotm / (DBLE(nxmax)*DBLE(nymax))

  END SUBROUTINE pote

END MODULE picsub
