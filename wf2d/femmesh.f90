! femmesh.f90

MODULE femmesh

    USE wfcomm, &
         nelm_max=>NEMAX, &
         node_max=>NNMAX, &
         node_nside_nelm=>NDELM, &
         xnode=>RNODE, &
         ynode=>ZNODE, &
         xnode_min=>BDRMIN, &
         xnode_max=>BDRMAX, &
         ynode_min=>BDZMIN, &
         ynode_max=>BDZMAX, &
         nseg_max=>NSDMAX, &
         node_nseg=>NDSID, &
         nelm_nseg=>NESID, &
         nside_nseg=>INSID, &
         nseg_nside_nelm=>NSDELM, &
         mode_nseg=>KASID, &
         mode_node=>KANOD, &
         nelm1_nside_nelm=>KNELM, &
         nside1_nside_nelm=>KSELM

  LOGICAL:: fem_mesh_allocated=.FALSE.
  INTEGER:: nside_max
  INTEGER,ALLOCATABLE:: nside_max_nelm(:)
  INTEGER:: nelm_node_max
  REAL(rkind):: vol_tot
  REAL(rkind),ALLOCATABLE:: xcenter_nelm(:),ycenter_nelm(:)
  REAL(rkind),ALLOCATABLE:: vol_nelm(:)
  REAL(rkind),ALLOCATABLE:: xcenter_nseg(:),ycenter_nseg(:)
  INTEGER,ALLOCATABLE:: idir_nside_nelm(:,:)
  INTEGER,ALLOCATABLE:: nelm_nangle_node(:,:),nelm_max_node(:)

! --- prepartion of mesh related variables ---

  PUBLIC fem_meshprep
  PUBLIC fem_set_nseg
  PUBLIC fem_set_nbdy

CONTAINS

  SUBROUTINE fem_mesh_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: nseg_max_save=0,nside_max_save=0,nelm_max_save=0

    IF(fem_mesh_allocated) THEN
       IF((nseg_max .EQ.nseg_max_save ).AND. &
          (nside_max.EQ.nside_max_save).AND. &
          (nelm_max .EQ.nelm_max_save ).AND. &
          (nseg_max .EQ.nseg_max_save )) RETURN
       CALL fem_mesh_deallocate
    ENDIF

    ALLOCATE(xcenter_nelm(nelm_max),ycenter_nelm(nelm_max))
    ALLOCATE(vol_nelm(nelm_max))
    ALLOCATE(xcenter_nseg(nseg_max),ycenter_nseg(nseg_max))
    ALLOCATE(idir_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nelm_max_node(nelm_max))

    IF(nrank.eq.0) WRITE(6,'(A)') '## fem_mesh_allocated'
    IF(nrank.eq.0) WRITE(6,'(A,3I8)') '   nelm_max,nseg_max,nside_max=', &
                                          nelm_max,nseg_max,nside_max

    nseg_max_save =nseg_max
    nside_max_save=nside_max
    nelm_max_save =nelm_max
    nseg_max_save =nseg_max

    fem_mesh_allocated=.TRUE.
    RETURN
  END SUBROUTINE fem_mesh_allocate

  SUBROUTINE fem_mesh_deallocate
    IMPLICIT NONE

    DEALLOCATE(xcenter_nelm,ycenter_nelm)
    DEALLOCATE(vol_nelm)
    DEALLOCATE(xcenter_nseg,ycenter_nseg)
    DEALLOCATE(idir_nside_nelm)
    DEALLOCATE(nelm_max_node)

    fem_mesh_allocated=.FALSE.
    RETURN
  END SUBROUTINE fem_mesh_deallocate

!**********************************************************************
  SUBROUTINE fem_mesh_prep
!**********************************************************************

    IMPLICIT NONE
    INTEGER,ALLOCATABLE:: node_nsega(:,:)
    INTEGER,ALLOCATABLE:: nsega_nside_nelm(:,:)
    INTEGER,ALLOCATABLE:: nsega_pair(:),nsega_nseg(:),nseg_nsega(:)
    INTEGER,ALLOCATABLE:: ncount_max_nxzone_nyzone(:,:)
    INTEGER,ALLOCATABLE:: nsega_ncount_nxzone_nyzone(:,:,:)
    INTEGER,ALLOCATABLE:: nelm_nsega(:),nside_nsega(:)
    REAL(rkind),ALLOCATABLE:: xcenter_nsega(:),ycenter_nsega(:)
    INTEGER:: nseg,nsega,nseg_all,nsega1
    INTEGER:: nelm,nside,n1,n2,n3
    INTEGER:: nelm1,nside1
    INTEGER:: node,ncount,ncount_zone_max
    INTEGER:: nx,ny
    REAL(rkind):: x,y,xc,yc,xc1,yc1,xlen_zone,ylen_zone
    REAL(rkind):: x1,y1,x2,y2,x3,y3,xg,yg,s,v,xgsum,ygsum,vsum
    REAL(rkind):: sfactor,vfactor

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

!    WRITE(29,'(A)') 'ncount_max_nxzone_nyzone(nx,ny)'
!    WRITE(29,'(3I5)') &
!         ((nx,ny,ncount_max_nxzone_nyzone(nx,ny), &
!         ny=1,nyzone_max),nx=1,nxzone_max)

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
       WRIte(66,*) nseg,nsega,nsega_pair(nsega)
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
       nside_nseg(1,nseg)=nside
       IF(nside.EQ.0) &
            WRITE(6,*) 'nseg,nsega,nelm,nside=',nseg,nsega,nelm,nside
       nseg_nside_nelm(nside,nelm)=nseg
       nsega1=ABS(nsega_pair(nsega))
       IF(nsega1.NE.0) THEN
          nelm1=nelm_nsega(nsega1)
          nside1=nside_nsega(nsega1)
          nelm_nseg(2,nseg)=nelm1
          nelm1_nside_nelm(nside,nelm)=nelm1
          nside1_nside_nelm(nside,nelm)=nside1
          nseg_nside_nelm(nside1,nelm1)=nseg
          nelm1_nside_nelm(nside1,nelm1)=nelm
          nside1_nside_nelm(nside1,nelm1)=nside
       ELSE
          nelm_nseg(2,nseg)=0
          nelm1_nside_nelm(nside,nelm)=0
          nside1_nside_nelm(nside,nelm)=0
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
               (node_nside_nelm(nside,nelm),nside=1,3), &
               (nseg_nside_nelm(nside,nelm),nside=1,3)
       END DO
       DO nelm=1,nelm_max
          WRITE(6,'(A,I8,1P3E12.4)') 'nelm,xc,yc,v=', &
               nelm,xcenter_nelm(nelm),ycenter_nelm(nelm),vol_nelm(nelm)
       END DO
       WRITE(6,'(A,1PE12.4)') 'vol_tot=', vol_tot
    END IF

    CALL fem_setup_nelm_node
    RETURN
  END SUBROUTINE fem_mesh_prep

  !   --- setup nelm_node ---

  SUBROUTINE fem_setup_nelm_node
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
    IF(nrank.EQ.0) &
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


! *** convert nsega to nseg ***

  SUBROUTINE fem_set_nseg

    IMPLICIT NONE
    INTEGER,ALLOCATABLE:: node_nsega(:,:),nsega_nelm(:,:)
    INTEGER,ALLOCATABLE:: nsega_pair(:),nsega_nseg(:),nseg_nsega(:)
    INTEGER:: nelm,nseg,nside,n1
    INTEGER:: nseg1,nsega,nelm1,nelm2,nside1,k
    INTEGER:: node
    REAL(rkind):: x,y
    INTEGER,PARAMETER:: nside_max_nelm=3
    INTEGER:: nseg_all

    ! --- Evaluate min and max of xnode and ynode of active elements ---

    node=node_nside_nelm(1,1)
    xnode_min=xnode(node)
    xnode_max=xnode(node)
    ynode_min=ynode(node)
    ynode_max=ynode(node)
    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm
          node=node_nside_nelm(nside,nelm)
          x=xnode(node)
          xnode_min=MIN(x,xnode_min)
          xnode_max=MAX(x,xnode_max)
          y=ynode(node)
          ynode_min=MIN(y,ynode_min)
          ynode_max=MAX(y,ynode_max)
       END DO
    END DO

    ! --- Define nsega_nelm, node_nsega, and nelm_nseg ---
    !    --- First, define nsega for a segment
    !       --- unique nsega is given for all sides of a element
    !       --- therefore most of segments are counted twice
    !    --- node_nsega(1) and node_nsega(2) are generated
    !       --- according to the order of node_nside_nelm

    nside_max=3
    nseg_all=nside_max_nelm*nelm_max
    
    ALLOCATE(node_nsega(2,nseg_all),nsega_nelm(nside_max,nelm_max))

    nsega=0
    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm
          nsega=nsega+1
          nsega_nelm(nside,nelm)=nsega
          node_nsega(1,nsega)=node_nside_nelm(nside,nelm)
          IF(nside.EQ.nside_max_nelm) THEN
             node_nsega(2,nsega)=node_nside_nelm(1,nelm)
          ELSE
             node_nsega(2,nsega)=node_nside_nelm(nside+1,nelm)
          END IF
       END DO
    END DO
    nseg_all=nsega

    !    --- find two segments which have the same pair of nodes
    !           nsega_pair(lower nsega) = higher_nseg
    !           nsega_pair(higher nsega)=-lower_nseg

    ALLOCATE(nsega_pair(nseg_all))
    ALLOCATE(nsega_nseg(nseg_all),nseg_nsega(nseg_all))

    DO nsega=1,nseg_all
       nsega_pair(nsega)=0
    END DO

    DO nsega=1,nseg_all
       n1=node_nsega(1,nsega)
       DO nseg1=nsega+1,nseg_all
          IF(nsega_pair(nseg1).EQ.0) THEN
             IF(node_nsega(2,nseg1).EQ.n1) THEN
                IF(node_nsega(1,nseg1).EQ.node_nsega(2,nsega)) THEN
                   nsega_pair(nsega)=nseg1
                   nsega_pair(nseg1)=-nsega
                   EXIT
                ENDIF
             END IF
          END IF
       END DO
    END DO

    !    --- asign new nseg: removing negative (higher) nsega  ---
    !             nseg_nsega: positive: lower nsega of a segment
    !             nseg_nsega: zero:     boundary segment
    !             nseg_nsega: negative: higher nsega of a seg (to be removed)

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

!    !   --- allocate fluid mesh variables using nseg_max---
!
    CALL wfsid_allocate

    !   --- convert node_nsega to node_nseg ---

    DO nsega=1,nseg_all
       nseg=nseg_nsega(nsega)
       IF(nseg.GT.0) THEN
          IF(nseg.LE.nsega) THEN
             node_nseg(1,nseg)=node_nsega(1,nsega)
             node_nseg(2,nseg)=node_nsega(2,nsega)
          ELSE IF(nseg.GE.nsega) THEN
             WRITE(6,*) 'XX sdx_meshprep error in rewriting node_nseg'
             WRITE(6,*) '   nsega,nseg=',nsega,nseg
             STOP
          END IF
       END IF
    END DO
       
    !    --- set nseg_nside_nelm, idir_nside_nelm and nelm_nseg---
    !        --- nseg_nsega=positive, idir_nside_nlem=0: positive direction
    !        --- nseg_nsega=negative, idir_nside_nlem=1: negative direction
    !        --- nseg_nside_nelm is now all positive or 0

    DO nseg=1,nseg_max
       nelm_nseg(1,nseg)=0
       nelm_nseg(2,nseg)=0
    END DO

    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm
          nsega=nsega_nelm(nside,nelm)
          nseg=nseg_nsega(nsega)
          IF(nseg.GT.0) THEN
             nseg_nside_nelm(nside,nelm)=nseg
!             idir_nside_nelm(nside,nelm)=0
             nelm_nseg(1,nseg)=nelm
          ELSE IF(nseg.LT.0) THEN
             nseg_nside_nelm(nside,nelm)=-nseg
!             idir_nside_nelm(nside,nelm)=1
             nelm_nseg(2,-nseg)=nelm
          END IF
       END DO
    END DO

    !    --- set nelm1_nside_nelm to identify adjacent element ---
    !       --- first, find nelm different from original nelm and set it nelm1

    DO nelm=1,nelm_max
       DO nside=1,nside_max_nelm
          nseg=nseg_nside_nelm(nside,nelm)
          nelm1=nelm_nseg(1,nseg)
          nelm2=nelm_nseg(2,nseg)
          IF(nelm1.NE.nelm.AND.nelm2.EQ.nelm) THEN
             nelm1_nside_nelm(nside,nelm)=nelm1
          ELSE IF(nelm2.NE.nelm.AND.nelm1.EQ.nelm) THEN
             nelm1_nside_nelm(nside,nelm)=nelm2
          ELSE
             WRITE(6,'(A)')     'XX sdx_meshprep: CONFILCT in nseg_nside_nelm'
             WRITE(6,'(A,5I8)') '   nelm,nside,nseg,nelm1,nelm2=', &
                  nelm,nside,nseg,nelm1,nelm2
             STOP
          END IF
          nelm1=nelm1_nside_nelm(nside,nelm)

          ! --- if nelm1=0, the segment is a boundary segment.
          !     if not, find nside1 for nelm1 corresponds to the nseg
          !     and set it to nside1_nside_nelm

          IF(nelm1.EQ.0) THEN ! boundary
             mode_nseg(nseg)=1
             nside1=0
          ELSE ! internal boundary
             mode_nseg(nseg)=0
             nside1=0
             DO k=1,nside_max_nelm
                IF(nseg.EQ.nseg_nside_nelm(k,nelm1)) nside1=k
             END DO
             IF(nside1.EQ.0) THEN
                WRITE(6,*) 'XX sdx_meshprep: nside1 not found'
                WRITE(6,*) '   nelm,nside,nseg,nelm1=',nelm,nside,nseg,nelm1
                STOP
             END IF
             nside1_nside_nelm(nside,nelm)=nside1
          END IF
       END DO
    END DO

    DEALLOCATE(node_nsega,nsega_nelm)
    DEALLOCATE(nseg_nsega,nsega_nseg,nsega_pair)

  END SUBROUTINE fem_set_nseg

!     ****** Set Boundary Attribute for Side and Node ******

SUBROUTINE fem_set_nbdy

  implicit none
  INTEGER:: nbdy,nseg1,nseg2,nseg3,nelm,node
  INTEGER:: nbdy_xmax1,nbdy_xmax2,nseg,node1,node2,node3
  INTEGER:: node_start,node_next,i,j
  REAL(rkind):: xmax1,xmax2,x,y1,y2
  INTEGER,ALLOCATABLE:: nbdy_order(:),nseg_nbdy_temp(:)

! ----- SET nseg_max & nelm1_nside_nelm -----
! NBSID :: Number of Boundary Side
! nelm1_nside_nelm(ISD,nelm) :: THE ELEMENT WHICH SHARE THE SIDE (ISD,nelm)
!                  IF nelm1_nside_nelm.EQ.0, SIDE IS BOUNDARY

     DO nelm=1,nelm_max
        WRITE(61,'(A,I6,A,3I6,2X,3I6)') 'nelm1_nside_nelm:',nelm,' : ', &
             nelm1_nside_nelm(1,nelm),nelm1_nside_nelm(2,nelm), &
             nelm1_nside_nelm(3,nelm), &
             node_nside_nelm(1,nelm),node_nside_nelm(2,nelm), &
             node_nside_nelm(3,nelm)
     END DO
     DO nelm=1,nelm_max
        WRITE(62,'(A,I6,A,3I6)') 'nseg_nside_nelm:',nelm,' : ', &
             nseg_nside_nelm(1,nelm),nseg_nside_nelm(2,nelm), &
             nseg_nside_nelm(3,nelm)
     END DO
     DO nseg=1,nseg_max
        WRITE(63,'(A,I6,A,3I6)') 'node_nseg: ',nseg,' : ', &
             node_nseg(1,nseg),node_nseg(2,nseg)
     END DO

     ! --- SIDE IS ON THE BOUNDARY, SIDE IS nelmW ---

     DO node=1,node_max
        mode_node(node)=0
     END DO

     DO nelm=1,nelm_max
        node1=node_nside_nelm(1,nelm)
        node2=node_nside_nelm(2,nelm)
        node3=node_nside_nelm(3,nelm)
        IF(nelm1_nside_nelm(1,nelm).EQ.0) THEN
           mode_node(node1)=1
           mode_node(node2)=1
        ENDIF
        IF(nelm1_nside_nelm(2,nelm).EQ.0) THEN
           mode_node(node2)=1
           mode_node(node3)=1
        ENDIF
        IF(nelm1_nside_nelm(3,nelm).EQ.0) THEN
           mode_node(node3)=1
           mode_node(node1)=1
        ENDIF
     ENDDO

  
!  DO nseg=1,nseg_max
!     write(*,*) "nelm_nseg,nside_nseg,nseg",nelm_nseg(nseg), &
!     nside_nseg(nseg),nseg
!  ENDDO

!  DO nelm=1,nelm_max
!     DO nside=1,3
!        write(*,*) "nelm,nside,nseg_nside_nelm", &
!                    nelm,nside,nseg_nside_nelm(nside,nelm)
!     end DO
!  end DO

!    DO nelm=1,nelm_max
!     DO nside=1,3
!        write(*,*) "nelm,nside,nelm1_nside_nelm", &
!                    nelm,nside,nelm1_nside_nelm(nside,nelm)
!     end DO
!  end DO

 ! ----- SET mode_nseg -----

  if(nrank.eq.0) WRITE(6,*) '------- SETBDY set mode_nseg start ---'

  DO nseg=1,nseg_max
     mode_nseg(nseg)=0
  ENDDO

  DO nelm=1,nelm_max
     nseg1=ABS(nseg_nside_nelm(1,nelm))
     nseg2=ABS(nseg_nside_nelm(2,nelm))
     nseg3=ABS(nseg_nside_nelm(3,nelm))
     IF(nelm1_nside_nelm(1,nelm).EQ.0) mode_nseg(nseg1)=1
     IF(nelm1_nside_nelm(2,nelm).EQ.0) mode_nseg(nseg2)=1
     IF(nelm1_nside_nelm(3,nelm).EQ.0) mode_nseg(nseg3)=1
  ENDDO

  ! === count boundary segment ===

  nbdy=0
  DO nelm=1,nelm_max
     nseg1=ABS(nseg_nside_nelm(1,nelm))
     nseg2=ABS(nseg_nside_nelm(2,nelm))
     nseg3=ABS(nseg_nside_nelm(3,nelm))
     IF(nelm1_nside_nelm(1,nelm).EQ.0) nbdy=nbdy+1
     IF(nelm1_nside_nelm(2,nelm).EQ.0) nbdy=nbdy+1
     IF(nelm1_nside_nelm(3,nelm).EQ.0) nbdy=nbdy+1
  END DO
  nbdy_max=nbdy
  WRITE(6,'(A,I8)') '## nbdy_max=',nbdy_max

  ! === allocate nbdy_nseg and nseg_nbdy ===

  IF(ALLOCATED(nbdy_nseg)) THEN
     IF(SIZE(nbdy_nseg).NE.nseg_max) THEN
        DEALLOCATE(nbdy_nseg)
        ALLOCATE(nbdy_nseg(nseg_max))
     END IF
  ELSE
     ALLOCATE(nbdy_nseg(nseg_max))
  END IF
  IF(ALLOCATED(nseg_nbdy)) THEN
     IF(SIZE(nseg_nbdy).NE.nbdy_max) THEN
        DEALLOCATE(nseg_nbdy)
        ALLOCATE(nseg_nbdy(nbdy_max))
     END IF
  ELSE
     ALLOCATE(nseg_nbdy(nbdy_max))
  END IF

  DO nseg=1,nseg_max
     nbdy_nseg(nseg)=0
  END DO
  
  nbdy=0
  DO nelm=1,nelm_max
     nseg1=ABS(nseg_nside_nelm(1,nelm))
     nseg2=ABS(nseg_nside_nelm(2,nelm))
     nseg3=ABS(nseg_nside_nelm(3,nelm))
     IF(nelm1_nside_nelm(1,nelm).EQ.0) THEN
        nbdy=nbdy+1
        nseg_nbdy(nbdy)=nseg1
        nbdy_nseg(nseg1)=nbdy
     END IF
     IF(nelm1_nside_nelm(2,nelm).EQ.0) THEN
        nbdy=nbdy+1
        nseg_nbdy(nbdy)=nseg2
        nbdy_nseg(nseg2)=nbdy
     END IF
     IF(nelm1_nside_nelm(3,nelm).EQ.0) THEN
        nbdy=nbdy+1
        nseg_nbdy(nbdy)=nseg3
        nbdy_nseg(nseg3)=nbdy
     END IF
  END DO
  IF(nbdy.NE.nbdy_max) THEN
     WRITE(6,*) 'XX ERROR in fem_set_nby. inconsistent nbdy,nbdy_max=', &
          nbdy,nbdy_max
     STOP
  END IF

  ! --- find two boundaries with maximum x for nbdy=1 ---

  nbdy_xmax1=1
  nbdy_xmax2=1
  xmax1=xnode_min
  xmax2=xnode_min
  
  DO nbdy=1,nbdy_max
     nseg=nseg_nbdy(nbdy)
     node1=node_nseg(1,nseg)
     node2=node_nseg(2,nseg)
     x=0.5D0*(xnode(node1)+xnode(node2))
     IF(x.GT.xmax1) THEN
        nbdy_xmax2=nbdy_xmax1
        nbdy_xmax1=nbdy
        xmax2=xmax1
        xmax1=x
     ELSE
        IF(x.GT.xmax2) THEN
           nbdy_xmax2=nbdy
           xmax2=x
        END IF
     END IF
     IF(idebuga(25).EQ.1) WRITE(6,'(A,I6,ES12.4,2I6,2ES12.4)') &
          '## nbdy:',nbdy,x,nbdy_xmax1,nbdy_xmax2,xmax1,xmax2
  END DO
  IF(idebuga(25).EQ.1) &
       WRITE(6,'(A,4I8)') '## nbdy1,nseg1,nbdy2,nseg2:', &
       nbdy_xmax1,nbdy_xmax2, &
       nseg_nbdy(nbdy_xmax1),nseg_nbdy(nbdy_xmax2)
  IF(idebuga(25).EQ.1) &
       WRITE(6,'(A,3ES12.4)') '## xmax1,xmax2,diff=',xmax1,xmax2,xmax1-xmax2

  ! --- choose larger y if xmax1 and xmax2 are almost same ---

  ALLOCATE(nbdy_order(nbdy_max))
  
  IF(ABS(xmax1-xmax2).LT.1.D-16) THEN
     nseg=nseg_nbdy(nbdy_xmax1)
     node1=node_nseg(1,nseg)
     node2=node_nseg(2,nseg)
     y1=0.5D0*(ynode(node1)+ynode(node2))
     nseg=nseg_nbdy(nbdy_xmax2)
     node1=node_nseg(1,nseg)
     node2=node_nseg(2,nseg)
     y2=0.5D0*(ynode(node1)+ynode(node2))
     IF(y2.GT.y1) THEN
        nbdy_order(1)=nbdy_xmax2
     ELSE
        nbdy_order(1)=nbdy_xmax1
     END IF
     IF(idebuga(25).EQ.1) &
          WRITE(6,'(A,2I6,ES12.4,2I6,ES12.4)') &
          '## Y: ',nbdy_xmax1,nseg_nbdy(nbdy_xmax1),y1, &
                   nbdy_xmax2,nseg_nbdy(nbdy_xmax2),y2
     IF(idebuga(25).EQ.1) &
          WRITE(6,'(A,I6)') '## nbdy_order(1):',nbdy_order(1)
  ELSE
     nbdy_order(1)=nbdy_xmax1
  END IF

  ! --- choose upper node to find next boundary segment ---

  nseg=nseg_nbdy(nbdy_order(1))
  node1=node_nseg(1,nseg)
  node2=node_nseg(2,nseg)
  y1=ynode(node1)
  y2=ynode(node2)
  IF(y2.GT.y1) THEN
     node_start=node1
     node_next=node2
     IF(idebuga(25).EQ.1) &
          WRITE(6,'(A,3I6,2ES12.4)') &
          '# nseg,node1,node2,y1,y2=',nseg,node1,node2,y1,y2
  ELSE
     node_start=node2
     node_next=node1
     IF(idebuga(25).EQ.1) &
          WRITE(6,'(A,3I6,2ES12.4)') &
          '# nseg,node2,node1,y2,y1=',nseg,node2,node1,y2,y1
  END IF

  ! --- find next boundary segment ---
  
  IF(idebuga(25).EQ.1) &
       WRITE(6,'(A,3I8,2ES12.4)') &
       '## next nbdy:',0,nbdy_order(1),node_next, &
       xnode(node_next),ynode(node_next)
  DO i=1,nbdy_max-1
     nbdy=nbdy_order(i)
     DO j=1,nbdy_max
        IF(j.NE.nbdy) THEN
           nseg=nseg_nbdy(j)
           node1=node_nseg(1,nseg)
           node2=node_nseg(2,nseg)
           IF(node1.EQ.node_next) THEN
              node_next=node2
              nbdy_order(i+1)=j
              EXIT
           ENDIF
           IF(node2.EQ.node_next) THEN
              node_next=node1
              nbdy_order(i+1)=j
              EXIT
           ENDIF
        END IF
     END DO
     IF(idebuga(25).EQ.1) &
          WRITE(6,'(A,3I8,2ES12.4)') &
          '## next nbdy:',i,j,node_next,xnode(node_next),ynode(node_next)
  END DO

  IF(idebuga(25).EQ.1) &
       WRITE(6,'(A,2I8)') '## node_next,node_start: ',node_next,node_start

  ALLOCATE(nseg_nbdy_temp(nbdy_max))

  DO nbdy=1,nbdy_max
     nseg_nbdy_temp(nbdy)=nseg_nbdy(nbdy_order(nbdy))
  END DO
  DO nbdy=1,nbdy_max
     nseg_nbdy(nbdy)=nseg_nbdy_temp(nbdy)
  END DO

  DO nbdy=1,nbdy_max
     nseg=nseg_nbdy(nbdy)
     node1=node_nseg(1,nseg)
     IF(idebuga(25).EQ.1) &
          WRITE(6,'(A,3I6,2ES12.4)') '## nbdy,nseg,node,x,y:', &
          nbdy,nseg,node1,xnode(node1),ynode(node1)
  END DO

  DEALLOCATE(nbdy_order,nseg_nbdy_temp)
           
  RETURN
END SUBROUTINE fem_set_nbdy

END MODULE femmesh
