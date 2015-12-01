!  ***** TASK/PIC PREPARATION *****

MODULE picsub
  PRIVATE
  PUBLIC poisson_f,poisson_m,efield,bfield,kine,pote
 
CONTAINS

!***********************************************************************
    subroutine poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
                         rho,phi,rhof,phif,awk,afwk,cform,ipssn)
!***********************************************************************
      implicit none
      real(8), dimension(nxmax1,nymax1) :: rho,phi
      real(8), dimension(nxmax,nymax) :: awk
      complex(8), dimension(nxmaxh1,nymax) :: afwk
      complex(8), dimension(nxmaxh1,nymax) :: rhof, phif
      real(8), dimension(nxmaxh1,nymax) :: cform
      integer(4) :: nxmax, nymax,nxmaxh1,nxmax1,nymax1,ipssn
      integer(4) :: ifset

      IF(ipssn.EQ.0) THEN
         call poisson_sub(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)
         ifset = 0
         call fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,rho,rhof,awk,afwk,ifset)
      ELSE

       !.......... fourier transform rho
       ifset = -1
       call fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,rho,rhof,awk,afwk,ifset)

       !.......... calculate phi from rho in fourier space
       ipssn = 1
       call poisson_sub(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)

       !.......... inverse fourier transform phi
       ifset = 1
       call fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,phi,phif,awk,afwk,ifset)

    END IF
  end subroutine poisson_f

!***********************************************************************
    subroutine poisson_sub(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)
!***********************************************************************
      implicit none
      complex(8), dimension(nxmaxh1,nymax) :: rhof, phif
      real(8), dimension(nxmaxh1,nymax) :: cform
      integer(4) :: nxmax, nymax, nxmaxh1
      integer(4) :: nx, ny, nymaxh,ipssn
      real(8) :: alxi, alyi, pi, twopi, am, an, amn2, afsp, afsp2

         afsp  = 1.5 !+++ size of finite-size-particle
         afsp2 = afsp * afsp 

      if( ipssn .eq. 0 ) then
         pi    = 3.14159265358979d0
         twopi = 2.d0 * pi 
         alxi   = 1.d0 / dble(nxmax)
         alyi   = 1.d0 / dble(nymax)
         nymaxh    = nymax / 2

         do ny = 1, nymaxh+1
         do nx = 1, nxmaxh1
            am = twopi * dble(nx-1) * alxi
            an = twopi * dble(ny-1) * alyi
            amn2 = am*am + an*an
           
            if( nx .eq. 1 .and. ny .eq. 1 ) then
               cform(nx,ny) = 0.0
            else 
               cform(nx,ny) = 1.d0 / amn2 * exp( - amn2 * afsp2 ) 
            endif

            if( ny .ge. 2 .and. ny .le. nymaxh ) then
               cform(nx,nymax-ny+2) = cform(nx,ny)
            endif
         end do
         end do

      else

         !----- solve poisson equation
         do ny = 1, nymax
         do nx = 1, nxmaxh1
            phif(nx,ny) = cform(nx,ny) * rhof(nx,ny)
         end do
         end do 

      endif
      
    end subroutine poisson_sub

!***********************************************************************
    subroutine fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,a,af,awk,afwk,ifset)
!***********************************************************************
      implicit none
      include 'fftw3.f'
      real(8), dimension(nxmax1,nymax1) :: a
      real(8), dimension(nxmax,nymax) :: awk
      real(8) :: alx, aly 
      complex(8), dimension(nxmaxh1,nymax) :: af, afwk
      integer(4) :: nxmax, nymax, nxmaxh1, nxmax1, nymax1
      integer(4) :: ifset, nx, ny
      !....integer(8), save :: FFTW_ESTIMATE
      integer(8), save :: plan1, plan2

      alx = dble(nxmax)
      aly = dble(nymax)

      if( ifset .eq. 0 ) then
         !----- initialization of fourier transform
         call dfftw_plan_dft_r2c_2d(plan1,nxmax,nymax,awk,afwk,FFTW_ESTIMATE)
         call dfftw_plan_dft_c2r_2d(plan2,nxmax,nymax,afwk,awk,FFTW_ESTIMATE)

      elseif( ifset .eq. -1 ) then

         !----- fourier transform
         do ny = 1, nymax
         do nx = 1, nxmax
            awk(nx,ny) = a(nx,ny)
         end do
         end do

         call dfftw_execute(plan1)

         do ny = 1, nymax
         do nx = 1, nxmaxh1
         af(nx,ny) = afwk(nx,ny) / ( alx * aly )
         end do
         end do

      else

         !----- inverse fourier transform

         do ny = 1, nymax
         do nx = 1, nxmaxh1
         afwk(nx,ny) = af(nx,ny)
         end do
         end do

         call dfftw_execute(plan2)

         do ny = 1, nymax
         do nx = 1, nxmax
            a(nx,ny) = awk(nx,ny)
         end do
         end do

         do ny = 1, nymax
            a(nxmax1,ny) = a(1,ny)
         end do

         do nx = 1, nxmax1
            a(nx,nymax1) = a(nx,1)
         end do

      endif

    end subroutine fftpic

!***********************************************************************
    subroutine poisson_m(nxmax1,nymax1,rho,phi,ipssn, &
                         model_matrix0,model_matrix1,model_matrix2, &
                         tolerance_matrix)
!***********************************************************************
      USE libmpi
      USE commpi
      USE libmtx
      implicit none
      real(8), dimension(nxmax1,nymax1) :: rho,phi
      real(8), dimension(:),allocatable :: x
      real(8):: tolerance_matrix
      integer :: nxmax1,nymax1,nxymax,ipssn
      integer :: model_matrix0,model_matrix1,model_matrix2
      integer :: nxmax,nymax,mode,imax,isize,jwidth,ileng
      integer,save:: status=0,istart,iend,irange
      integer :: i,nx,ny,l,m,its

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
         l=mod(i-1,isize)+1
         m=(i-1)/isize+1
         if(m.gt.1) CALL mtx_set_matrix(i,i-isize,1.d0)
         if(l.gt.1) CALL mtx_set_matrix(i,i-1,1.d0)
         CALL mtx_set_matrix(i,i,-4.d0)
         if(l.lt.isize) CALL mtx_set_matrix(i,i+1,1.d0)
         if(m.lt.ileng) CALL mtx_set_matrix(i,i+isize,1.d0)
      ENDDO

      IF(mode.EQ.0) THEN
         DO i=istart,iend
            nx=mod(i-1,isize)+1
            ny=(i-1)/isize+1
            CALL mtx_set_source(i,-rho(nx+1,ny+1))
         ENDDO
      ELSE
         DO i=istart,iend
            ny=mod(i-1,isize)+1
            nx=(i-1)/isize+1
            CALL mtx_set_source(i,-rho(nx+1,ny+1))
         ENDDO
      END IF
      CALL mtx_solve(model_matrix0,tolerance_matrix,its, &
                     methodKSP=model_matrix1,methodPC=model_matrix2)
      IF((its.ne.0).and.(nrank.eq.0)) write(6,*) 'mtx_solve: its=',its
      CALL mtx_get_vector(x)
      IF(mode.EQ.0) THEN
         DO i=istart,iend
            nx=mod(i-1,isize)+1
            ny=(i-1)/isize+1
            phi(nx+1,ny+1)=x(i)
         ENDDO
      ELSE
         DO i=istart,iend
            ny=mod(i-1,isize)+1
            nx=(i-1)/isize+1
            phi(nx+1,ny+1)=x(i)
         ENDDO
      END IF
      DEALLOCATE(x)
      status=2
      CALL mtx_cleanup
    END subroutine poisson_m

!***********************************************************************
    subroutine efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
                      ex,ey,ez,esx,esy,esz,emx,emy,emz,model_boundary)
!***********************************************************************
      implicit none
      real(8), dimension(0:nxmax,0:nymax) ::  &
           phi,Ax,Ay,Az,Axb,Ayb,Azb,ex,ey,ez,esx,esy,esz,emx,emy,emz
      real(8):: dt
      integer :: nxmax, nymax, nx, ny, nxm, nxp, nym, nyp,model_boundary
      if(model_boundary .eq. 0) then
         do nx = 0, nymax
         do ny = 0, nxmax

            nxm = nx - 1
            nxp = nx + 1
            nym = ny - 1
            nyp = ny + 1

            if( nx .eq. 0  )    nxm = nxmax - 1
            if( nx .eq. nxmax ) nxp = 1
            if( ny .eq. 0  )    nym = nymax - 1
            if( ny .eq. nymax ) nyp = 1
            esx(nx,ny) = 0.5d0 * ( phi(nxm,ny) - phi(nxp,ny))
            esy(nx,ny) = 0.5d0 * ( phi(nx,nym) - phi(nx,nyp))
            esz(nx,ny) = 0.d0
            emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
            emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
            emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt
            ex(nx,ny) = esx(nx,ny) + emx(nx,ny)
            ey(nx,ny) = esy(nx,ny) + emy(nx,ny)
            ez(nx,ny) = esz(nx,ny) + emz(nx,ny)

       end do
       end do
     else if (model_boundary .ne. 0) then
         do nx = 0, nymax
         do ny = 0, nxmax

            nxm = nx - 1
            nxp = nx + 1
            nym = ny - 1
            nyp = ny + 1

            if( nx .eq. 0  )    nxm = 0 
            if( nx .eq. nxmax ) nxp = nxmax
            if( ny .eq. 0  )    nym = 0
            if( ny .eq. nymax ) nyp = nymax

            if(nx .eq. 0 .or. nx .eq. nxmax) then
               esx(nx,ny) = phi(nxm,ny) - phi(nxp,ny)
            else
               esx(nx,ny) = 0.5d0 * ( phi(nxm,ny) - phi(nxp,ny))
            endif
            
            if(ny .eq. 0 .or. ny .eq. nymax) then
               esy(nx,ny) = phi(nx,nym) - phi(nx,nyp)
            else
               esy(nx,ny) = 0.5d0 * ( phi(nx,nym) - phi(nx,nyp))
            end if

            esz(nx,ny) = 0.d0
            emx(nx,ny) = - ( Ax(nx,ny) - Axb(nx,ny) ) / dt
            emy(nx,ny) = - ( Ay(nx,ny) - Ayb(nx,ny) ) / dt
            emz(nx,ny) = - ( Az(nx,ny) - Azb(nx,ny) ) / dt
            ex(nx,ny) = esx(nx,ny) + emx(nx,ny)
            ey(nx,ny) = esy(nx,ny) + emy(nx,ny)
            ez(nx,ny) = esz(nx,ny) + emz(nx,ny)

       end do
       end do
    end if
     end subroutine efield

!***********************************************************************
    subroutine bfield(nxmax,nymax,Ax,Ay,Az,Axb,Ayb,Azb, &
                                  bx,by,bz,bxbg,bybg,bzbg,bb)
!***********************************************************************
      implicit none
      real(8), dimension(0:nymax) :: bxnab,bznab
      real(8), dimension(0:nxmax) :: bynab
      real(8), dimension(0:nxmax,0:nymax) :: bx,by,bz,bxbg,bybg,bzbg,bb
      real(8), dimension(0:nxmax,0:nymax) :: Ax,Ay,Az,Axb,Ayb,Azb
      integer :: nxmax, nymax, nx, ny, nxp, nyp, nxm, nym, model_boundary
      if(model_boundary .eq. 0) then
         do ny = 0, nymax
         do nx = 0, nxmax
            nxm = nx - 1
            nxp = nx + 1
            nym = ny - 1
            nyp = ny + 1

            if( nx .eq. 0  )    nxm = nxmax - 1
            if( nx .eq. nxmax ) nxp = 1
            if( ny .eq. 0  )    nym = nymax - 1
            if( ny .eq. nymax ) nyp = 1
         
            bx(nx,ny) = 0.25d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                      - Az(nx,nym) - Azb(nx,nym)) + bxbg(nx,ny)
            by(nx,ny) = 0.25d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                      - Az(nxm,ny) - Azb(nxm,ny)) + bybg(nx,ny)
            bz(nx,ny) = 0.25d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                      - Ay(nxm,ny) - Ayb(nxm,ny) &
                      - (Ax(nx,nyp) + Axb(nx,nyp) &
                      -Ax(nx,nym) - Axb(nx,nym)))+ bzbg(nx,ny)
            bb(nx,ny) = SQRT(bx(nx,ny)**2+by(nx,ny)**2+bz(nx,ny)**2)
         end do
         end do
         elseif(model_boundary .ne. 0) then
       do ny = 0, nymax
       do nx = 0, nxmax
            nxm = nx - 1
            nxp = nx + 1
            nym = ny - 1
            nyp = ny + 1

            if( nx .eq. 0  )    nxm = 0
            if( nx .eq. nxmax ) nxp = nxmax
            if( ny .eq. 0  )    nym = 0
            if( ny .eq. nymax ) nyp = nymax
         
            bx(nx,ny) = 0.25d0 * (Az(nx,nyp) + Azb(nx,nyp) &
                      - Az(nx,nym) - Azb(nx,nym)) + bxbg(nx,ny)
            by(nx,ny) = 0.25d0 * (Az(nxp,ny) + Azb(nxp,ny) &
                      - Az(nxm,ny) - Azb(nxm,ny)) + bybg(nx,ny)
            if(nx .eq. 0 .or. nx .eq. nxmax .or. ny .eq. 0 .or. ny .eq. nymax) then
            bz(nx,ny) = 0.5d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                      - Ay(nxm,ny) - Ayb(nxm,ny) &
                      - (Ax(nx,nyp) + Axb(nx,nyp) &
                      -Ax(nx,nym) - Axb(nx,nym)))+ bzbg(nx,ny)
            else
            bz(nx,ny) = 0.25d0 * (Ay(nxp,ny) + Ayb(nxp,ny) &
                      - Ay(nxm,ny) - Ayb(nxm,ny) &
                      - (Ax(nx,nyp) + Axb(nx,nyp) &
                      -Ax(nx,nym) - Axb(nx,nym)))+ bzbg(nx,ny)
            end if
            bb(nx,ny) = SQRT(bx(nx,ny)**2+by(nx,ny)**2+bz(nx,ny)**2)
        end do
        end do
     endif 
    end subroutine bfield

!***********************************************************************
    subroutine kine(npmax,vx,vy,vz,akin,mass)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: vx, vy, vz
      real(8) :: akin, mass
      integer(4) :: npmax, np
      akin = 0.d0
      do np = 1, npmax
         akin = akin + vx(np)**2 + vy(np)**2 + vz(np)**2
      end do

      akin = 0.5 * akin * mass /dble(npmax)
    end subroutine kine

!***********************************************************************
    subroutine pote(nxmax,nymax,ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg,vcfact, &
                    apote,apotm)
!***********************************************************************
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg
      real(8) :: apote,apotm,vcfact
      integer(4) :: nxmax, nymax, nx, ny

      apote = 0.d0
      apotm = 0.d0

         do ny = 1, nymax-1
         do nx = 1, nxmax-1
            apote = apote + 0.5d0*(ex(nx,ny)**2 + ey(nx,ny)**2 + ez(nx,ny)**2)
            apotm = apotm + 0.5d0*((bx(nx,ny)-bxbg(nx,nx))**2 &
                  + (by(nx,ny)-bybg(nx,nx))**2 &
                  + (bz(nx,ny)-bzbg(nx,nx))**2)
         end do
         end do

         do nx = 0, nxmax,nxmax
         do ny = 1, nymax-1
            apote = apote + 0.5D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
            apotm = apotm + 0.5D0*((bx(nx,ny)-bxbg(nx,nx))**2 &
                  + (by(nx,ny)-bybg(nx,nx))**2 &
                  + (bz(nx,ny)-bzbg(nx,nx))**2)
         end do
         end do

         do ny = 0, nymax,nymax
         do nx = 1, nxmax-1
            apote = apote + 0.5D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
            apotm = apotm + 0.5D0*((bx(nx,ny)-bxbg(nx,nx))**2 &
                  + (by(nx,ny)-bybg(nx,nx))**2 &
                  + (bz(nx,ny)-bzbg(nx,nx))**2)
         end do
         end do

         do ny = 0, nymax,nymax
         do nx = 0, nxmax,nxmax
            apote = apote + 0.25D0*(ex(nx,ny)**2+ey(nx,ny)**2+ez(nx,ny)**2)
            apotm = apotm + 0.25D0*((bx(nx,ny)-bxbg(nx,nx))**2 &
                  + (by(nx,ny)-bybg(nx,nx))**2 &
                  + (bz(nx,ny)-bzbg(nx,nx))**2)
         end do
         end do
      apote = 0.5D0 * apote / (dble(nxmax)*dble(nymax))
      apotm = 0.5D0 * vcfact**2 * apotm / (dble(nxmax)*dble(nymax))

    end subroutine pote

end MODULE picsub
