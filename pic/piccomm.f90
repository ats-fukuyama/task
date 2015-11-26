!  ***** TASK/PIC COMMON *****
! parameter

MODULE piccomm_parm

  USE bpsd_kinds
  USE bpsd_constants

  INTEGER:: npxmax,npymax,nxmax,nymax,ntmax,ntstep,ntgstep,ntpstep
  INTEGER:: model_boundary,model_antenna,model_wg
  INTEGER:: model_matrix0,model_matrix1,model_matrix2
  REAL(rkind):: dt,me,mi,chrge,chrgi,te,ti,&
       bxmin,bxmax,bymin,bymax,bzmin,bzmax,vcfact,omega,eps
  REAL(rkind):: jxant,jyant,jzant,phxant,phyant,phzant
  REAL(rkind):: xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg,ph_wg,rot_wg,eli_wg
  REAL(rkind):: tolerance_matrix

END MODULE piccomm_parm

MODULE piccomm 
		
  USE piccomm_parm
  USE commpi

  INTEGER:: npmax,nxmaxh1,nxmax1,nymax1,nxymax
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: ex,ey,ez,esx,esy,esz,emx,emy,emz
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: bx,by,bz,bxbg,bybg,bzbg
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

  REAL(rkind),ALLOCATABLE,DIMENSION(:):: timet,akinet,akinit,aktott
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: atott,apotet,apotmt,aptott

  REAL(8) :: ctome, ctomi, &
             vte, vti,     &
             akine , akini , aktot , apote,  apotm , aptot , atot , &
             akine0, akini0, aktot0, apote0, apotm0, aptot0, atot0, &
             akine1, akine2, akini1, akini2, time,  &
             x1, x2, y1, y2, z1, z2 ,alx, aly, alz, &
             wkword, wtime, wtime1, wtime2
  integer :: ntcount, ntgcount, ntpcount, ntgmax, ntpmax
  integer :: ifset, ipssn, iran
  integer :: ierr

CONTAINS

  SUBROUTINE pic_allocate

    INTEGER,SAVE:: nxmax_save=0, nymax_save=0
    INTEGER,SAVE:: npxmax_save=0, npymax_save=0

    IF(nxmax  == nxmax_save  .AND. &
       nymax  == nymax_save  .AND. &
       npxmax == npxmax_save .AND. &
       npymax == npymax_save )  RETURN
       
    IF(ALLOCATED(ex)) CALL pic_deallocate

    ALLOCATE(ex(0:nxmax,0:nymax),ey(0:nxmax,0:nymax),ez(0:nxmax,0:nymax))
    ALLOCATE(esx(0:nxmax,0:nymax),esy(0:nxmax,0:nymax),esz(0:nxmax,0:nymax))
    ALLOCATE(emx(0:nxmax,0:nymax),emy(0:nxmax,0:nymax),emz(0:nxmax,0:nymax))
    ALLOCATE(bx(0:nxmax,0:nymax),by(0:nxmax,0:nymax),bz(0:nxmax,0:nymax))
    ALLOCATE(bxbg(0:nxmax,0:nymax),bybg(0:nxmax,0:nymax),bzbg(0:nxmax,0:nymax))
    ALLOCATE(rho(0:nxmax,0:nymax),phi(0:nxmax,0:nymax),phib(0:nxmax,0:nymax))
    ALLOCATE(jx(0:nxmax,0:nymax),jy(0:nxmax,0:nymax),jz(0:nxmax,0:nymax))
    ALLOCATE(awk(nxmax,nymax))
    ALLOCATE(xe(npmax),ye(npmax),ze(npmax),vxe(npmax),vye(npmax),vze(npmax))
    ALLOCATE(xi(npmax),yi(npmax),zi(npmax),vxi(npmax),vyi(npmax),vzi(npmax))
    ALLOCATE(xeb(npmax),yeb(npmax),zeb(npmax),xib(npmax),yib(npmax),zib(npmax))
    ALLOCATE(vparae(npmax),vperpe(npmax),vparai(npmax),vperpi(npmax))
    ALLOCATE(cform(nxmaxh1,nymax))
    ALLOCATE(rhof(nxmaxh1,nymax),phif(nxmaxh1,nymax),afwk(nxmaxh1,nymax))
    ALLOCATE(Ax(0:nxmax,0:nymax),Ay(0:nxmax,0:nymax),Az(0:nxmax,0:nymax))
    ALLOCATE(Axb(0:nxmax,0:nymax),Ayb(0:nxmax,0:nymax),Azb(0:nxmax,0:nymax))
    ALLOCATE(Axbb(0:nxmax,0:nymax),Aybb(0:nxmax,0:nymax),Azbb(0:nxmax,0:nymax))
    ALLOCATE(bb(0:nxmax,0:nymax),AA(0:nxmax,0:nymax))

    nxmax_save=nxmax
    nymax_save=nymax
    npxmax_save=npxmax
    npymax_save=npymax

    RETURN
  END SUBROUTINE pic_allocate

  SUBROUTINE pic_deallocate

    IF(ALLOCATED(ex)) THEN
       DEALLOCATE(ex,ey,ez,esx,esy,esz,emx,emy,emz)
       DEALLOCATE(bx,by,bz,bxbg,bybg,bzbg)
       DEALLOCATE(rho,phi,phib,awk)
       DEALLOCATE(jx,jy,jz)
       DEALLOCATE(xe,ye,ze,vxe,vye,vze)
       DEALLOCATE(xi,yi,zi,vxi,vyi,vzi)
       DEALLOCATE(xeb,yeb,zeb,xib,yib,zib)
       DEALLOCATE(Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
       DEALLOCATE(bb,AA)
       DEALLOCATE(vparae,vperpe,vparai,vperpi)
       DEALLOCATE(cform,rhof,phif,afwk)
    END IF

    IF(ALLOCATED(timet)) THEN
       DEALLOCATE(timet,akinet,akinit,aktott)
       DEALLOCATE(atott,apotet,apotmt,aptott)
    END IF


    RETURN
  END SUBROUTINE pic_deallocate
END MODULE piccomm
