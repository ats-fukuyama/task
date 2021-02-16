! sak2.f90

MODULE sak2

  PRIVATE
  PUBLIC sak_2

CONTAINS

  SUBROUTINE sak_2

    USE sakcomm
    USE saksub
    USE libgrf
    IMPLICIT NONE
    COMPLEX(dp):: cw    ! omega/omegape
    REAL(dp):: wr,wi    ! real and imag parts of cw
    REAL(dp):: rk2      ! (k vte/omegape)**2 = (k lambda)**2
    REAL(dp):: sg2      ! sigma**2 = (l q_theta/k)**2
    COMPLEX(dp):: cf    ! 
    INTEGER:: nxmax,nymax,nx,ny
    REAL(dp):: wrmin,wrmax,wimin,wimax,delwr,delwi
    REAL(dp),ALLOCATABLE:: wra(:),wia(:),fra(:,:),fia(:,:),faa(:,:)
    EXTERNAL GSOPEN,GSCLOS,PAGES,PAGEE
    
    rk2=0.1D0
    sg2=0.D0
    wrmin=1.0D0
    wrmax=1.5D0
    wimin=-0.1D0
    wimax= 0.1D0
    nxmax=101
    nymax=101

1   CONTINUE
    WRITE(6,'(A/6ES12.4,2I4)') &
         '## INPUT: rk2,sg2,wrmin,wrmax,wimin,wimax,nxmax,nymax?', &
         rk2,sg2,wrmin,wrmax,wimin,wimax,nxmax,nymax
    READ(5,*,ERR=1,END=9000) rk2,sg2,wrmin,wrmax,wimin,wimax,nxmax,nymax

    ALLOCATE(wra(nxmax),wia(nymax))
    ALLOCATE(fra(nxmax,nymax),fia(nxmax,nymax),faa(nxmax,nymax))
    
    delwr=(wrmax-wrmin)/(nxmax-1)
    delwi=(wimax-wimin)/(nymax-1)
    DO nx=1,nxmax
       wra(nx)=wrmin+delwr*(nx-1)
    END DO
    DO ny=1,nymax
       wia(ny)=wimin+delwi*(ny-1)
    END DO
    DO ny=1,nymax
       DO nx=1,nxmax
          cw=CMPLX(wra(nx),wia(ny))
          cf=cfeps(cw,rk2,sg2)
          fra(nx,ny)=REAL(cf)
          fia(nx,ny)=AIMAG(cf)
          faa(nx,ny)=fra(nx,ny)**2+fia(nx,ny)**2
       END DO
    END DO

    CALL PAGES
    CALL grd2d(1,wra,wia,fra,nxmax,nxmax,nymax,'@f_real(wr,wi)@')
    CALL grd2d(2,wra,wia,fia,nxmax,nxmax,nymax,'@f_imag(wr,wi)@')
    CALL grd2d(3,wra,wia,faa,nxmax,nxmax,nymax,'@f_abs(wr,wi)@')
    CALL PAGEE
    DEALLOCATE(wra,wia,fra,fia,faa)
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE sak_2

END MODULE sak2
