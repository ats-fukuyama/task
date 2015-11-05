!  ***** TASK/PIC COMMON *****
! parameter

MODULE piccomm_parm

  USE bpsd_kinds
  USE bpsd_constants

  INTEGER:: npx,npy,nx,ny,iend,nhmod
  REAL(rkind):: dt,me,mi,chrge,chrgi,te,ti,bxbg,bybg,bzbg,eps,c,omega

END MODULE piccomm_parm

MODULE piccomm 
		
  USE piccomm_parm

  INTEGER:: np,nxh1,nx1,ny1,nxy
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: ex,ey,ez,esx,esy,esz,emx,emy,emz
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: bx,by,bz
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: rho,phi,phib,awk,jx,jy,jz,&
                                           Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: bb,AA
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: xe,ye,ze,vxe,vye,vze, &
                                         xi,yi,zi,vxi,vyi,vzi, &
                                         xeb,yeb,zeb,xib,yib,zib
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: vparae,vperpe,vparai,vperpi
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: expi,eypi,ezpi,bxpi,bypi,bzpi
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

    INTEGER,SAVE:: nx_save=0, ny_save=0
    INTEGER,SAVE:: npx_save=0, npy_save=0

    IF(nx  == nx_save  .AND. &
       ny  == ny_save  .AND. &
       npx == npx_save .AND. &
       npy == npy_save )  RETURN
       
    IF(ALLOCATED(ex)) CALL pic_deallocate

    ALLOCATE(ex(0:nx,0:ny),ey(0:nx,0:ny),ez(0:nx,0:ny))
    ALLOCATE(esx(0:nx,0:ny),esy(0:nx,0:ny),esz(0:nx,0:ny))
    ALLOCATE(emx(0:nx,0:ny),emy(0:nx,0:ny),emz(0:nx,0:ny))
    ALLOCATE(bx(0:nx,0:ny),by(0:nx,0:ny),bz(0:nx,0:ny))
    ALLOCATE(rho(0:nx,0:ny),phi(0:nx,0:ny),phib(0:nx,0:ny))
    ALLOCATE(jx(0:nx,0:ny),jy(0:nx,0:ny),jz(0:nx,0:ny))
    ALLOCATE(awk(nx,ny))
    ALLOCATE(xe(np),ye(np),ze(np),vxe(np),vye(np),vze(np))
    ALLOCATE(xi(np),yi(np),zi(np),vxi(np),vyi(np),vzi(np))
    ALLOCATE(xeb(np),yeb(np),zeb(np),xib(np),yib(np),zib(np))
    ALLOCATE(vparae(np),vperpe(np),vparai(np),vperpi(np))
    ALLOCATE(cform(nxh1,ny))
    ALLOCATE(rhof(nxh1,ny),phif(nxh1,ny),afwk(nxh1,ny))
    ALLOCATE(Ax(0:nx,0:ny),Ay(0:nx,0:ny),Az(0:nx,0:ny))
    ALLOCATE(Axb(0:nx,0:ny),Ayb(0:nx,0:ny),Azb(0:nx,0:ny))
    ALLOCATE(Axbb(0:nx,0:ny),Aybb(0:nx,0:ny),Azbb(0:nx,0:ny))
    ALLOCATE(bb(0:nx,0:ny),AA(0:nx,0:ny))

    nx_save=nx
    ny_save=ny
    npx_save=npx
    npy_save=npy

    RETURN
  END SUBROUTINE pic_allocate

  SUBROUTINE pic_deallocate

    IF(ALLOCATED(ex)) THEN
       DEALLOCATE(ex,ey,ez,rho,phi,phib,awk)
       DEALLOCATE(bx,by,bz)
       DEALLOCATE(xe,ye,ze,vxe,vye,vze)
       DEALLOCATE(xi,yi,zi,vxi,vyi,vzi)
       DEALLOCATE(jx,jy,jz)
       DEALLOCATE(xeb,yeb,zeb,xib,yib,zib)
       DEALLOCATE(Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
       DEALLOCATE(bb,AA)
       DEALLOCATE(vparae,vperpe,vparai,vperpi)
!       DEALLOCATE(expi,eypi,ezpi,bxpi,bypi,bzpi)
       DEALLOCATE(cform,rhof,phif,afwk)
    END IF

    IF(ALLOCATED(timet)) THEN
       DEALLOCATE(timet,akinet,akinit,aktott,apott,atott)
    END IF


    RETURN
  END SUBROUTINE pic_deallocate
END MODULE piccomm
