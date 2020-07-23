! femmeshprep.f90

MODULE femmeshprep

  USE wfcomm,ONLY: rkind
  LOGICAL:: fem_mesh_allocated=.FALSE.
  INTEGER(rkind):: nseg_max,nside_max
  INTEGER(rkind):: nelm_node_max
  REAL(rkind):: xnode_min,xnode_max,ynode_min,ynode_max
  REAL(rkind):: vol_tot
  INTEGER(rkind),ALLOCATABLE:: nside_max_nelm(:)
  REAL(rkind),ALLOCATABLE:: xcenter_nelm(:),ycenter_nelm(:)
  REAL(rkind),ALLOCATABLE:: vol_nelm(:)
  REAL(rkind),ALLOCATABLE:: xcenter_nseg(:),ycenter_nseg(:)
  INTEGER(rkind),ALLOCATABLE:: node_nseg(:,:),nelm_nseg(:,:)
  INTEGER(rkind),ALLOCATABLE:: idir_nside_nelm(:,:),idseg_nseg(:)
  INTEGER(rkind),ALLOCATABLE:: nseg_nside_nelm(:,:)
  INTEGER(rkind),ALLOCATABLE:: nelm1_nside_nelm(:,:),nside1_nside_nelm(:,:)
  INTEGER(rkind),ALLOCATABLE:: nelm_nangle_node(:,:),nelm_max_node(:)

! --- prepartion of mesh related variables ---

  PUBLIC fem_meshprep

CONTAINS

  SUBROUTINE fem_mesh_allocate
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    IMPLICIT NONE
    INTEGER,SAVE:: nseg_max_save=0,nside_max_save=0,nelm_max_save=0

    IF(fem_mesh_allocated) THEN
       IF((nseg_max .EQ.nseg_max_save ).AND. &
          (nside_max.EQ.nside_max_save).AND. &
          (nelm_max .EQ.nelm_max_save )) RETURN
       CALL fem_mesh_deallocate
    ENDIF

    ALLOCATE(xcenter_nelm(nelm_max),ycenter_nelm(nelm_max))
    ALLOCATE(vol_nelm(nelm_max))
    ALLOCATE(xcenter_nseg(nseg_max),ycenter_nseg(nseg_max))
    ALLOCATE(idseg_nseg(nseg_max))
    ALLOCATE(nelm_nseg(2,nseg_max),node_nseg(2,nseg_max))
    ALLOCATE(nseg_nside_nelm(nside_max,nelm_max))
    ALLOCATE(idir_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nelm1_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nside1_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nelm_max_node(nelm_max))

    IF(nrank.eq.0) WRITE(6,'(A)') '## fem_mesh_allocated'
    IF(nrank.eq.0) WRITE(6,'(A,3I8)') '   nelm_max,nseg_max,nside_max=', &
                                          nelm_max,nseg_max,nside_max

    nseg_max_save =nseg_max
    nside_max_save=nside_max
    nelm_max_save =nelm_max

    fem_mesh_allocated=.TRUE.
    RETURN
  END SUBROUTINE fem_mesh_allocate

  SUBROUTINE fem_mesh_deallocate
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    IMPLICIT NONE

    DEALLOCATE(xcenter_nelm,ycenter_nelm)
    DEALLOCATE(vol_nelm)
    DEALLOCATE(xcenter_nseg,ycenter_nseg)
    DEALLOCATE(idseg_nseg)
    DEALLOCATE(nelm_nseg,node_nseg)
    DEALLOCATE(nseg_nside_nelm)
    DEALLOCATE(idir_nside_nelm)
    DEALLOCATE(nelm1_nside_nelm)
    DEALLOCATE(nside1_nside_nelm)
    DEALLOCATE(nelm_max_node)

    fem_mesh_allocated=.FALSE.
    RETURN
  END SUBROUTINE fem_mesh_deallocate

!**********************************************************************
  SUBROUTINE fem_meshprep
!**********************************************************************

    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    IMPLICIT NONE
    INTEGER,ALLOCATABLE:: node_nsega(:,:),nsega_nside_nelm(:,:)
    INTEGER,ALLOCATABLE:: nsega_pair(:),nsega_nseg(:),nseg_nsega(:)
    INTEGER,ALLOCATABLE:: ncount_max_nxzone_nyzone(:,:)
    INTEGER,ALLOCATABLE:: nsega_ncount_nxzone_nyzone(:,:,:)
    INTEGER,ALLOCATABLE:: nelm_nsega(:),nside_nsega(:)
    REAL(rkind),ALLOCATABLE:: xcenter_nsega(:),ycenter_nsega(:)
    INTEGER:: nseg_all
    INTEGER:: nelm,nseg,nside,n1,n2,n3
    INTEGER:: nseg1,nsega,nelm1,nelm2,nside1,k
    INTEGER:: node,ncount,ncount_zone_max
    INTEGER:: nsega1,nx,ny
    REAL(rkind):: x,y,xc,yc,xc1,yc1,xlen_zone,ylen_zone
    REAL(rkind):: x1,y1,x2,y2,x3,y3,xg,yg,s,v,xgsum,ygsum,vsum
    REAL(rkind):: sfactor,vfactor,weight,xce,yce,xcs,ycs,dl,ds,xn,yn

    nxzone_max=100
    nyzone_max=100

! --- initialize nside_nelm ---

    IF(ALLOCATED(nside_max_nelm)) DEALLOCATE(nside_max_nelm)
    ALLOCATE(nside_max_nelm(nelm_max))
    DO nelm=1,nelm_max
       nside_max_nelm(nelm)=3
    END DO
    nside_max=3
    
! --- Evaluate min and max of xnode and ynode of active elements ---

    node=node_nside_nelm(1,1)
    xnode_min=xnode(node)
    xnode_max=xnode(node)
    ynode_min=ynode(node)
    ynode_max=ynode(node)
    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm(nelm)
          node=node_nside_nelm(nside,nelm)
          x=xnode(node)
          xnode_min=MIN(x,xnode_min)
          xnode_max=MAX(x,xnode_max)
          y=ynode(node)
          ynode_min=MIN(y,ynode_min)
          ynode_max=MAX(y,ynode_max)
       END DO
    END DO

! --- Define nseg_nelm, node_nseg, and nelm_nseg ---
!    --- First, define each side of elements has different segment number
!       --- therefore most of segments are counted twice

    nseg_all=nside_max*nelm_max

    ALLOCATE(node_nsega(2,nseg_all),nsega_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nelm_nsega(nseg_all),nside_nsega(nseg_all))
    nelm_nsega(1:nseg_all)=0.D0
    nside_nsega(1:nseg_all)=0.D0
    
    nsega=0
    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm(nelm)
          nsega=nsega+1
          nsega_nside_nelm(nside,nelm)=nsega
          node_nsega(1,nsega)=node_nside_nelm(nside,nelm)
          IF(nside.EQ.nside_max_nelm(nelm)) THEN
             node_nsega(2,nsega)=node_nside_nelm(1,nelm)
          ELSE
             node_nsega(2,nsega)=node_nside_nelm(nside+1,nelm)
          END IF
          nelm_nsega(nsega)=nelm
          nside_nsega(nsega)=nside
       END DO
    END DO
    nseg_all=nsega

!    --- zoning segment with the position of segment center ---
    
    ALLOCATE(xcenter_nsega(nseg_all),ycenter_nsega(nseg_all))
    DO nsega=1,nseg_all
       x1=xnode(node_nsega(1,nsega))
       y1=ynode(node_nsega(1,nsega))
       x2=xnode(node_nsega(2,nsega))
       y2=ynode(node_nsega(2,nsega))
       xcenter_nsega(nsega)=0.5D0*(x1+x2)
       ycenter_nsega(nsega)=0.5D0*(y1+y2)
    END DO

  ! --- setup zone data ---
  
    ALLOCATE(ncount_max_nxzone_nyzone(nxzone_max,nyzone_max))
    ncount_max_nxzone_nyzone(1:nxzone_max,1:nyzone_max)=0

    xlen_zone=(xnode_max-xnode_min)/nxzone_max
    ylen_zone=(ynode_max-ynode_min)/nyzone_max

    ! Count number of segments in a zone

    DO nsega=1,nseg_all
       xc=xcenter_nsega(nsega)
       yc=ycenter_nsega(nsega)
       nx=INT((xc-xnode_min)/xlen_zone)+1
       ny=INT((yc-ynode_min)/ylen_zone)+1
       IF(nx.GT.nxzone_max) nx=nxzone_max  ! adjust at maximum xmax
       IF(ny.GT.nyzone_max) ny=nyzone_max  ! adjust at maximum ymax
       ncount_max_nxzone_nyzone(nx,ny)=ncount_max_nxzone_nyzone(nx,ny)+1
    END DO

    ! Count maximum number of elements in a zone

    ncount_zone_max=0
    DO ny=1,nyzone_max
       DO nx=1,nxzone_max
          ncount_zone_max=MAX(ncount_zone_max,ncount_max_nxzone_nyzone(nx,ny))
       END DO
    ENDDO

    WRITE(29,'(A)') 'ncount_max_nxzone_nyzone(nx,ny)'
    WRITE(29,'(3I5)') &
         ((nx,ny,ncount_max_nxzone_nyzone(nx,ny), &
         ny=1,nyzone_max),nx=1,nxzone_max)

    ! Set nelm of elements in a zone

    ALLOCATE(nsega_ncount_nxzone_nyzone(ncount_zone_max,nxzone_max,nyzone_max))
    
    ncount_max_nxzone_nyzone(1:nxzone_max,1:nyzone_max)=0
    DO nsega=1,nseg_all
       xc=xcenter_nsega(nsega)
       yc=ycenter_nsega(nsega)
       nx=INT((xc-xnode_min)/xlen_zone)+1
       ny=INT((yc-ynode_min)/ylen_zone)+1
       IF(nx.GT.nxzone_max) nx=nxzone_max  ! adjust at maximum xmax
       IF(ny.GT.nyzone_max) ny=nyzone_max  ! adjust at maximum ymax

       ncount_max_nxzone_nyzone(nx,ny)=ncount_max_nxzone_nyzone(nx,ny)+1
       nsega_ncount_nxzone_nyzone(ncount_max_nxzone_nyzone(nx,ny),nx,ny) &
                  =nsega
    END DO
    
!    --- find two segments which have same pair of nodes

    ALLOCATE(nsega_nseg(nseg_all),nseg_nsega(nseg_all))
    ALLOCATE(nsega_pair(nseg_all))

    DO nsega=1,nseg_all
       nsega_pair(nsega)=0
    END DO

    DO nsega=1,nseg_all
       xc=xcenter_nsega(nsega)
       yc=ycenter_nsega(nsega)
       nx=INT((xc-xnode_min)/xlen_zone)+1
       ny=INT((yc-ynode_min)/ylen_zone)+1
       IF(nx.GT.nxzone_max) nx=nxzone_max  ! adjust at maximum xmax
       IF(ny.GT.nyzone_max) ny=nyzone_max  ! adjust at maximum ymax

       DO ncount=1,ncount_max_nxzone_nyzone(nx,ny)
          nsega1=nsega_ncount_nxzone_nyzone(ncount,nx,ny)          
          IF(nsega_pair(nsega1).EQ.0) THEN
             xc1=xcenter_nsega(nsega1)
             yc1=ycenter_nsega(nsega1)
             IF(ABS(xc-xc1).LE.1.D-80.AND.ABS(yc-yc1).LE.1.D-80) THEN
                nsega_pair(nsega)=nsega1
                nsega_pair(nsega1)=-nsega
                EXIT
             END IF
          ENDIF
       END DO
    END DO

!    --- asign new nseg ---

    nseg=0
    DO nsega=1,nseg_all
       IF(nsega_pair(nsega).GT.0) THEN
          nseg=nseg+1
          nseg_nsega(nsega)=nseg
          nsega_nseg(nseg)=nsega
          nseg_nsega(nsega_pair(nsega))=-nseg
       ELSE IF(nsega_pair(nsega).EQ.0) THEN
          nseg=nseg+1
          nseg_nsega(nsega)=nseg
          nsega_nseg(nseg)=nsega
       END IF
    END DO
    nseg_max=nseg

!   --- allocate fluid mesh variables ---

    CALL fem_mesh_allocate

!   --- Setup segment variables ---

    DO nseg=1,nseg_max
       nsega=nsega_nseg(nseg)
       xcenter_nseg(nseg)=xcenter_nsega(nsega)
       ycenter_nseg(nseg)=ycenter_nsega(nsega)
       node_nseg(1,nseg)=node_nsega(1,nsega)
       node_nseg(2,nseg)=node_nsega(2,nsega)
       nelm=nelm_nsega(nsega)
       nside=nside_nsega(nsega)
       nelm_nseg(1,nseg)=nelm
       nseg_nside_nelm(nside,nelm)=nseg
       nsega1=ABS(nsega_pair(nsega))
       IF(nsega1.NE.0) THEN
          nelm1=nelm_nsega(nsega1)
          nside1=nside_nsega(nsega1)
          nelm_nseg(2,nseg)=nelm1
          nelm1_nside_nelm(nside,nelm)=nelm1
          nseg_nside_nelm(nside1,nelm1)=nseg
          nelm1_nside_nelm(nside1,nelm1)=nelm
       ELSE
          nelm_nseg(2,nseg)=0
          nelm1_nside_nelm(nside,nelm)=0
       END IF
    END DO
    
    DEALLOCATE(nseg_nsega,nsega_nseg,nsega_pair)
    DEALLOCATE(nsega_ncount_nxzone_nyzone)
    DEALLOCATE(ncount_max_nxzone_nyzone)
    DEALLOCATE(xcenter_nsega,ycenter_nsega)
    DEALLOCATE(nelm_nsega,nside_nsega)
    DEALLOCATE(node_nsega,nsega_nside_nelm)

! --- modelg = 0: rectanglar coordinates: 2D (volume means area)
! --- modelg = 1: cylindrical and toroidal coordinates: 3D

    IF(modelg.EQ.0) THEN
       sfactor=1.D0
       vfactor=0.D0
    ELSE
       sfactor=0.D0
       vfactor=1.D0
    END IF

! --- define elemnt center position: xcenter_nelm,ycenter_belm
! --- define volume of element: vol_nelm

    vol_tot=0.D0
    DO nelm=1,nelm_max
       xgsum=0.D0
       ygsum=0.D0
       vsum=0.D0
       DO nside=2,nside_max_nelm(nelm)-1
          n1=node_nside_nelm(1,nelm)
          x1=xnode(n1)
          y1=ynode(n1)
          n2=node_nside_nelm(nside,nelm)
          x2=xnode(n2)
          y2=ynode(n2)
          n3=node_nside_nelm(nside+1,nelm)
          x3=xnode(n3)
          y3=ynode(n3)
          s=0.5D0*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
          v=2.D0*pi*(x1+x2+x3)*s/3.D0
          xg=sfactor*(x1+x2+x3)*s/3.D0 &
            +vfactor*2.D0*pi*(x1*x1+x1*x2+x2*x2+x2*x3+x3*x3+x3*x1)*s/6.D0
          yg=sfactor*(y1+y2+y3)*s/3.D0 &
            +vfactor*2.D0*pi*(2.D0*x1*y1     +x1*y2     +x1*y3 &
                                  +x2*y1+2.D0*x2*y2     +x2*y3 &
                                  +x3*y1     +x3*y2+2.D0*x3*y3)*s/12.D0
          xgsum=xgsum+xg
          ygsum=ygsum+yg
          vsum=vsum+sfactor*s+vfactor*v
       END DO
       xcenter_nelm(nelm)=xgsum/vsum
       ycenter_nelm(nelm)=ygsum/vsum
!       IF(nprint.EQ.4.AND.MOD(nelm,10).EQ.0) THEN
!          WRITE(6,'(A,I8,1P3E12.4)') 'nelm,vol_disk,vol_calc,diff=', &
!               nelm,vol_nelm(nelm),vsum,vol_nelm(nelm)-vsum
!       END IF
       vol_nelm(nelm)=vsum
       vol_tot=vol_tot+vsum
    END DO

    IF(nprint.EQ.3) THEN
       DO nelm=1,nelm_max
          WRITE(6,'(A,I6,2X,4I6,2X,4I6)') 'nelm,node,nseg=', &
               nelm, &
               (node_nside_nelm(nside,nelm),nside=1,4), &
               (nseg_nside_nelm(nside,nelm),nside=1,4)
       END DO
       DO nelm=1,nelm_max
          WRITE(6,'(A,I8,1P3E12.4)') 'nelm,xc,yc,v=', &
               nelm,xcenter_nelm(nelm),ycenter_nelm(nelm),vol_nelm(nelm)
       END DO
       WRITE(6,'(A,1PE12.4)') 'vol_tot=', vol_tot
    END IF

    CALL fem_setup_nelm_node
    RETURN
  END SUBROUTINE fem_meshprep

  !   --- setup nelm_node ---

  SUBROUTINE fem_setup_nelm_node
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    IMPLICIT NONE
    INTEGER:: nelm,nside,node

    nelm_max_node(1:node_max)=0
    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm(nelm)
          node=node_nside_nelm(nside,nelm)
          nelm_max_node(node)=nelm_max_node(node)+1
       END DO
    END DO

    nelm_node_max=0
    DO node=1,node_max
       nelm_node_max=MAX(nelm_node_max,nelm_max_node(node))
    END DO
    WRITE(6,*) '-- fem_setup_nelm_node: nelm_node_max=',nelm_node_max

    IF(ALLOCATED(nelm_nangle_node)) DEALLOCATE(nelm_nangle_node)
    ALLOCATE(nelm_nangle_node(nelm_node_max,node_max))

    nelm_max_node(1:node_max)=0
    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm(nelm)
          node=node_nside_nelm(nside,nelm)
          nelm_max_node(node)=nelm_max_node(node)+1
          nelm_nangle_node(nelm_max_node(node),node)=nelm
       END DO
    END DO
    RETURN
  END SUBROUTINE fem_setup_nelm_node

END MODULE femmeshprep
