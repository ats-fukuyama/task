!  ***** TASK/PIC COMMON *****
! parameter

MODULE piccomm_parm

  USE bpsd_kinds
  USE bpsd_constants

  INTEGER:: npx,npy,npz,nx,ny,nz,iend,nhmod
  REAL(rkind):: dt,me,mi,chrge,chrgi,te,ti,eps

END MODULE piccomm_parm

MODULE piccomm !
		
  USE piccomm_parm

  INTEGER:: np,nxh1,nx1,ny1,nz1,nxy
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: ex,ey,rho,phi !
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: awk
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: xe,ye,ze,vxe,vye,vze, &
                                         xi,yi,zi,vxi,vyi,vzi
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: cform
  COMPLEX(rkind),ALLOCATABLE,DIMENSION(:,:):: rhof,phif,afwk

  REAL(rkind),ALLOCATABLE,DIMENSION(:):: timet,akinet,akinit,aktott,apott,atott

  REAL(8) :: ctome, ctomi,       &
             vte, vti, cfact, cfacti,              &
             akine , akini , aktot , apot , atot ,         &
             akine0, akini0, aktot0, apot0, atot0,         &
             akine1, akine2, akini1, akini2, time,         &
             x1, x2, y1, y2, z1, z2 ,alx, aly, alz,                &
             wkword, wtime, wtime1, wtime2
  integer :: iloop, ifset, ipssn, iran, iene, ienemax
  integer :: ierr, myid, nodes

CONTAINS

  SUBROUTINE pic_allocate

    INTEGER,SAVE:: nx_save=0, ny_save=0, nz_save=0
    INTEGER,SAVE:: npx_save=0, npy_save=0, npz_save=0

    IF(nx  == nx_save  .AND. &
       ny  == ny_save  .AND. &
       nz  == nz_save  .AND. &  
       npx == npx_save .AND. &
       npy == npy_save .AND. &
       npz == npz_save )  RETURN
       
    IF(ALLOCATED(ex)) CALL pic_deallocate

    ALLOCATE(ex(0:nx,0:ny),ey(0:nx,0:ny))
    ALLOCATE(rho(0:nx,0:ny),phi(0:nx,0:ny))
    ALLOCATE(awk(nx,ny))
    ALLOCATE(xe(np),ye(np),vxe(np),vye(np))
    ALLOCATE(xi(np),yi(np),vxi(np),vyi(np))
    ALLOCATE(cform(nxh1,ny))
    ALLOCATE(rhof(nxh1,ny),phif(nxh1,ny),afwk(nxh1,ny))

    nx_save=nx
    ny_save=ny
    nz_save=nz
    npx_save=npx
    npy_save=npy
    npz_save=npz

    RETURN
  END SUBROUTINE pic_allocate

  SUBROUTINE pic_deallocate

    IF(ALLOCATED(ex)) THEN
       DEALLOCATE(ex,ey,rho,phi,awk)
       DEALLOCATE(xe,ye,vxe,vye,xi,yi,vxi,vyi)
       DEALLOCATE(cform,rhof,phif,afwk)
    END IF

    IF(ALLOCATED(timet)) THEN
       DEALLOCATE(timet,akinet,akinit,aktott,apott,atott)
    END IF


    RETURN
  END SUBROUTINE pic_deallocate
END MODULE piccomm
