! feminterpolate.f90

MODULE feminterpolate

! Defined in femcomm.f90
! INTEGER:: nxzone_max,nyzone_max   ! maximum division of zoning in x and y

  USE wfcomm,ONLY: rkind
  INTEGER:: ncount_zone_max         ! maximum number of elements in a zone
  REAL(rkind):: xlen_zone,ylen_zone ! length of rectangular zone in x or y
  INTEGER,ALLOCATABLE:: ncount_max_nxzone_nyzone(:,:)
                                 ! number of elements in a zone
  INTEGER,ALLOCATABLE:: nelm_ncount_nxzone_nyzone(:,:,:)
                                 ! nelm for ncount in a zone
  
  PRIVATE
  PUBLIC fem_setup_zone,fem_find_nelm_for_xy,fem_interpolate_xy

CONTAINS

  ! --- setup zone data ---
  
  SUBROUTINE fem_setup_zone

    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    USE femmeshprep
    IMPLICIT NONE
    INTEGER:: nelm,node,nside,nx_min,nx_max,ny_min,ny_max,nx,ny
    REAL(rkind):: xmin,xmax,ymin,ymax
    
    IF(ALLOCATED(ncount_max_nxzone_nyzone)) THEN
       DEALLOCATE(ncount_max_nxzone_nyzone)
       DEALLOCATE(nelm_ncount_nxzone_nyzone)
    END IF
    ALLOCATE(ncount_max_nxzone_nyzone(nxzone_max,nyzone_max))
    ncount_max_nxzone_nyzone(1:nxzone_max,1:nyzone_max)=0

    xlen_zone=(xnode_max-xnode_min)/nxzone_max
    ylen_zone=(ynode_max-ynode_min)/nyzone_max

    ! Count number of elements in a zone

    DO nelm=1,nelm_max
       CALL xyrange_nelm(nelm,xmin,xmax,ymin,ymax)
       nx_min=INT((xmin-xnode_min)/xlen_zone)+1
       nx_max=nxzone_max-INT((xnode_max-xmax)/xlen_zone)+1
       ny_min=INT((ymin-ynode_min)/ylen_zone)+1
       ny_max=nyzone_max-INT((ynode_max-ymax)/ylen_zone)+1
       IF(nx_max.GT.nxzone_max) nx_max=nxzone_max  ! adjust at maximum xmax
       IF(ny_max.GT.nyzone_max) ny_max=nyzone_max  ! adjust at maximum ymax
       DO ny=ny_min,ny_max
          DO nx=nx_min,nx_max
             ncount_max_nxzone_nyzone(nx,ny)=ncount_max_nxzone_nyzone(nx,ny)+1
          END DO
       END DO
    END DO

    ! Count maximum number of elements in a zone

    ncount_zone_max=0
    DO ny=1,nyzone_max
       DO nx=1,nxzone_max
          ncount_zone_max=MAX(ncount_zone_max,ncount_max_nxzone_nyzone(nx,ny))
       END DO
    ENDDO
!       WRITE(29,'(A)') 'ncount_max_nxzone_nyzone(nx,ny)'
!       WRITE(29,'(3I5)') &
!            ((nx,ny,ncount_max_nxzone_nyzone(nx,ny), &
!            ny=1,nyzone_max),nx=1,nxzone_max)

    ! Set nelm of elements in a zone

    ALLOCATE(nelm_ncount_nxzone_nyzone(ncount_zone_max,nxzone_max,nyzone_max))
    
    ncount_max_nxzone_nyzone(1:nxzone_max,1:nyzone_max)=0
    DO nelm=1,nelm_max
       node=node_nside_nelm(1,nelm)
       xmin=xnode(node)
       xmax=xnode(node)
       ymin=ynode(node)
       ymax=ynode(node)
       DO nside=2,nside_max_nelm(nelm)
          node=node_nside_nelm(nside,nelm)
          xmin=MIN(xmin,xnode(node))
          xmax=MAX(xmax,xnode(node))
          ymin=MIN(ymin,ynode(node))
          ymax=MAX(ymax,ynode(node))
       END DO
       nx_min=INT((xmin-xnode_min)/xlen_zone)+1
       nx_max=nxzone_max-INT((xnode_max-xmax)/xlen_zone)+1
       ny_min=INT((ymin-ynode_min)/ylen_zone)+1
       ny_max=nyzone_max-INT((ynode_max-ymax)/ylen_zone)+1
       IF(nx_max.GT.nxzone_max) nx_max=nxzone_max  ! adjust at maximum xmax
       IF(ny_max.GT.nyzone_max) ny_max=nyzone_max  ! adjust at maximum ymax
!       WRITE(29,'(A,5I5)') 'nelm,nx_min/max,ny_min/max=', &
!            nelm,nx_min,nx_max,ny_min,ny_max
       DO ny=ny_min,ny_max
          DO nx=nx_min,nx_max
             ncount_max_nxzone_nyzone(nx,ny)=ncount_max_nxzone_nyzone(nx,ny)+1
             nelm_ncount_nxzone_nyzone(ncount_max_nxzone_nyzone(nx,ny),nx,ny) &
                  =nelm
          END DO
       END DO
    END DO
!    DO ny=1,2
!    DO nx=1,2
!    DO ncount=1,ncount_max_nxzone_nyzone(nx,ny)
!       nelm=nelm_ncount_nxzone_nyzone(ncount,nx,ny)
!       WRITE(29,'(A,4I5,1P2E12.4)') 'nx,ny,ncount,nelm,xc,yc=',&
!            nx,ny,ncount,nelm,xcenter_nelm(nelm),ycenter_nelm(nelm)
!    END DO
!    END DO
!    END DO
!    DO ny=7,8
!    DO nx=7,8
!    DO ncount=1,ncount_max_nxzone_nyzone(nx,ny)
!       nelm=nelm_ncount_nxzone_nyzone(ncount,nx,ny)
!       WRITE(29,'(A,4I5,1P2E12.4)') 'nx,ny,ncount,nelm,xc,yc=',&
!            nx,ny,ncount,nelm,xcenter_nelm(nelm),ycenter_nelm(nelm)
!    END DO
!    END DO
!    END DO
    RETURN
       
  END SUBROUTINE fem_setup_zone

  SUBROUTINE xyrange_nelm(nelm,xmin,xmax,ymin,ymax)
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    USE femmeshprep
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nelm
    REAL(rkind),INTENT(OUT):: xmin,xmax,ymin,ymax
    INTEGER:: node,nside
    
    node=node_nside_nelm(1,nelm)
    xmin=xnode(node)
    xmax=xnode(node)
    ymin=ynode(node)
    ymax=ynode(node)
    DO nside=2,nside_max_nelm(nelm)
       node=node_nside_nelm(nside,nelm)
       xmin=MIN(xmin,xnode(node))
       xmax=MAX(xmax,xnode(node))
       ymin=MIN(ymin,ynode(node))
       ymax=MAX(ymax,ynode(node))
    END DO
    RETURN
  END SUBROUTINE xyrange_nelm
   
  ! --- find nelm for x,y:  if nelm=0, start from the nearest nelm
  !                         otherwize start from given nelm
  !                         if (x,y) in the nelm, save the nelm and exit
  !                         otherwize look for the nside nearest and move there
  !                         if boundary side exits, record the nelm_boundary
  !                         if no boundary side exits, set the nelm_boundary 0
  !                         if new elm is the same as the nelm_boundary,
  !                         set nelm=0 and exit
  
  SUBROUTINE fem_find_nelm_for_xy(x,y,nelm)
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    USE femmeshprep
    USE libmpi
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x,y
    INTEGER,INTENT(INOUT):: nelm
    INTEGER:: nside,nseg,node1,node2,nelm1,nside1,nseg1,nelm2
    INTEGER:: nxzone,nyzone,ncount
    INTEGER:: nelm_ncount
    REAL(rkind):: xc,yc,x1,y1,x2,y2,xc1,yc1

    IF(nelm.NE.0) THEN   ! If nelm is given, look for arround
    
       IF(xy_in_nelm(x,y,nelm)) THEN ! no need to move from the present elm
!          WRITE(6,'(A,I5)') '-- fem_find_nelm_for_xy: level 0: nelm=',nelm
          RETURN
       END IF
       
       ! level1
       
       xc=xcenter_nelm(nelm)
       yc=ycenter_nelm(nelm)
       DO nside=1,nside_max_nelm(nelm)
          nseg=nseg_nside_nelm(nside,nelm)
          IF(nseg.EQ.0) CYCLE   ! boundary segment
          node1=node_nside_nelm(nside,nelm)
          IF(nside.EQ.nside_max_nelm(nelm)) THEN
             node2=node_nside_nelm(1,nelm)
          ELSE
             node2=node_nside_nelm(nside+1,nelm)
          END IF
          x1=xnode(node1)-xc
          y1=ynode(node1)-yc
          x2=xnode(node2)-xc
          y2=ynode(node2)-yc
          IF(xy_in_range(x-xc,y-yc,x1,y1,x2,y2)) THEN
             nelm1=nelm1_nside_nelm(nside,nelm)
             IF(nelm1.EQ.0) CYCLE
             IF(xy_in_nelm(x,y,nelm1)) THEN   ! found in level 1
                nelm=nelm1
!                WRITE(6,'(A,I5)') &
!                     '-- fem_find_nelm_for_xy: level 1: nelm=',nelm
                RETURN
             END IF

             ! level 2
             xc1=xcenter_nelm(nelm1)
             yc1=ycenter_nelm(nelm1)
             DO nside1=1,nside_max_nelm(nelm1)
                nseg1=nseg_nside_nelm(nside1,nelm1)
                IF(nseg1.EQ.0) CYCLE      ! boundary segment
                IF(nseg1.EQ.nseg) CYCLE   ! original segment
                node1=node_nside_nelm(nside1,nelm)
                IF(nside1.EQ.nside_max_nelm(nelm1)) THEN
                   node2=node_nside_nelm(1,nelm1)
                ELSE
                   node2=node_nside_nelm(nside1+1,nelm1)
                END IF
                x1=xnode(node1)-xc1
                y1=ynode(node1)-yc1
                x2=xnode(node2)-xc1
                y2=ynode(node2)-yc1
                IF(xy_in_range(x-xc1,y-yc1,x1,y1,x2,y2)) THEN
                   nelm2=nelm1_nside_nelm(nside1,nelm1)
                   IF(nelm2.EQ.0) CYCLE
                   IF(xy_in_nelm(x,y,nelm2)) THEN    ! found in level 2
                      nelm=nelm2
!                      WRITE(6,'(A,I5)') &
!                           '-- fem_find_nelm_for_xy: level 2: nelm=',nelm
                      RETURN
                   END IF
                END IF
             END DO
          END IF
       END DO
    END IF  ! not found in level 1 and level 2 adjacent elements
    
    ! find nelm by check all elements in a zone

    nxzone=INT((x-xnode_min)/xlen_zone)+1
    nyzone=INT((y-ynode_min)/ylen_zone)+1
    IF(nxzone.GT.nxzone_max) nxzone=nxzone_max
    IF(nyzone.GT.nyzone_max) nyzone=nyzone_max
    DO ncount=1,ncount_max_nxzone_nyzone(nxzone,nyzone)
       nelm_ncount=nelm_ncount_nxzone_nyzone(ncount,nxzone,nyzone)
       IF(xy_in_nelm(x,y,nelm_ncount)) THEN
          nelm=nelm_ncount
!          WRITE(6,'(A,I5)') &
!               '-- fem_find_nelm_for_xy: level 3: nelm=',nelm
          RETURN
       END IF
    END DO
!    WRITE(6,'(A)') &
!      'XX fem_find_nelm_for_xy: out of range: x,y,xnode_min/max,ynode_min/max'
!    WRITE(6,'(A,1P6E12.4)') &
!         '      ',x,y,xnode_min,xnode_max,ynode_min,ynode_max
    nelm=0
  END SUBROUTINE fem_find_nelm_for_xy

  ! Check (x,y) in nelm or not: left-hand-side of all sides or on any side
  
  FUNCTION xy_in_nelm(x,y,nelm)
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    USE femmeshprep
    IMPLICIT NONE
    LOGICAL:: xy_in_nelm
    REAL(rkind),INTENT(IN):: x,y
    INTEGER,INTENT(IN):: nelm
    INTEGER:: node,nside
    REAL(rkind):: x1,y1,x2,y2,f
!    REAL(rkind):: x3,y3
!    INTEGER,SAVE:: i=0

!    WRITE(6,'(A,2I5,1P2E12.4)') '1:',nelm,1, &
!         xnode(node_nside_nelm(1,nelm)), &
!         ynode(node_nside_nelm(1,nelm))
!    WRITE(6,'(A,2I5,1P2E12.4)') '2:',nelm,1, &
!         xnode(node_nside_nelm(2,nelm)), &
!         ynode(node_nside_nelm(2,nelm))
!    WRITE(6,'(A,2I5,1P2E12.4)') '3:',nelm,1, &
!         xnode(node_nside_nelm(3,nelm)), &
!         ynode(node_nside_nelm(3,nelm))
!    node=node_nside_nelm(1,nelm)
!    x1=xnode(node)
!    y1=ynode(node)
!    node=node_nside_nelm(2,nelm)
!    x2=xnode(node)
!    y2=ynode(node)
!    node=node_nside_nelm(3,nelm)
!    x3=xnode(node)
!    y3=ynode(node)
!    f=(y2-y1)*(x3-x1)-(x2-x1)*(y3-y1)
!    WRITE(6,'(A,1PE12.4)') 'f:',f
!    i=i+1
!    IF(i.EQ.5) STOP
    
    xy_in_nelm=.FALSE.
    DO nside=1,nside_max_nelm(nelm)
       node=node_nside_nelm(nside,nelm)
       x1=xnode(node)
       y1=ynode(node)
       IF(nside.EQ.nside_max_nelm(nelm)) THEN
          node=node_nside_nelm(1,nelm)
       ELSE
          node=node_nside_nelm(nside+1,nelm)
       END IF
       x2=xnode(node)
       y2=ynode(node)
       f=(y2-y1)*(x-x1)-(x2-x1)*(y-y1)
       IF(f.GT.1.D-8) THEN
!          WRITE(6,'(A,1P6E12.4)') &
!               'xy:',x,y,x1,y1,x2,y2
!          WRITE(6,'(A,2I6,1P5E12.4)') &
!               '.F:',nelm,nside,f,(y2-y1),(x-x1),(x2-x1),(y-y1)
          RETURN
       END IF
!       IF(f.GT.0.D0) RETURN
!       IF(f.GT.-1.D-8) RETURN
    END DO
    xy_in_nelm=.TRUE.
!    WRITE(6,'(A,1P6E12.4)') &
!         'xy:T:',x,y,x1,y1,x2,y2
    RETURN
  END FUNCTION xy_in_nelm
    
  ! Check (x,y) in the angle rage between (x1,y1) and (x2,y2) from origin

  FUNCTION xy_in_range(x,y,x1,y1,x2,y2)
    USE wfcomm,ONLY: rkind
    IMPLICIT NONE
    LOGICAL:: xy_in_range
    REAL(rkind),INTENT(IN):: x,y,x1,y1,x2,y2
    REAL(rkind):: f1,f2

    f1=y1*x-x1*y
    f2=y2*x-x2*y
    IF(f1.LE.0.D0.AND.f2.GE.0.D0) THEN
       xy_in_range=.TRUE.
    ELSE
       xy_in_range=.FALSE.
    END IF
    RETURN
  END FUNCTION xy_in_range


  SUBROUTINE fem_interpolate_xy(x,y,f,f_nelm,id)
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x,y
    REAL(rkind),INTENT(OUT):: f
    REAL(rkind),INTENT(IN):: f_nelm(nelm_max)
    INTEGER,SAVE:: nelm_save=0
    INTEGER,INTENT(IN):: id  ! 0 for new search, 1: use nelm previous search
    INTEGER:: nelm

    IF(id.EQ.0) THEN
       nelm=0
    ELSE
       nelm=nelm_save
    END IF
    
    CALL fem_find_nelm_for_xy(x,y,nelm)
    IF(nelm.EQ.0) THEN ! (x,y) is out of computation region
       f=0.D0
       RETURN
    ELSE
       nelm_save=nelm
    END IF

!    IF(model_interpolation.GT.2.OR.model_interpolation.EQ.0) THEN
       f=f_nelm(nelm)
!    ELSE
!       SELECT CASE(model_interpolation)
!       CASE(1) ! linear interpolation (discontinuous)
!          CALL fem_grad_f(nelm,f_nelm,dfx,dfy)
!!          WRITE(29,'(A,I5,1P5E12.4)') 'nelm,f,dfx,dfy=', &
!!               nelm,f_nelm(nelm),dfx,x-xcenter_nelm(nelm), &
!!                                 dfy,y-ycenter_nelm(nelm)
!          f=f_nelm(nelm)+dfx*(x-xcenter_nelm(nelm)) &
!                        +dfy*(y-ycenter_nelm(nelm))
!       CASE(2) ! linear interpolation (continuouas)
!          CALL fem_linear_interporate(x,y,nelm,f_nelm,f)
!       END SELECT
!    END IF
    RETURN
  END SUBROUTINE fem_interpolate_xy

  SUBROUTINE fem_linear_interporate(x,y,nelm,f_nelm,f)
    USE wfcomm, &
         nelm_max=>nemax,node_max=>nnmax,node_nside_nelm=>ndelm, &
         xnode=>rnode,ynode=>znode
    USE femmeshprep
    USE libinv
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x,y,f_nelm(nelm_max)
    INTEGER,INTENT(IN):: nelm
    REAL(rkind),INTENT(OUT):: f
    INTEGER:: nside,nangl,node1,node2,nelm0,nelm1,nelm2,ierr
    REAL(rkind):: x0,y0,f0,x1,y1,f1,x2,y2,f2
    REAL(rkind):: weight,sum_of_weight,sum_of_f_weight
    REAL(rkind):: cmat(3,3),a,b,c

    ! confirm (x,y) in nelm

    nelm0=nelm
    CALL fem_find_nelm_for_xy(x,y,nelm0)
    x0=xcenter_nelm(nelm0)
    y0=ycenter_nelm(nelm0)
    
    ! find triangle including (x,y)

    DO nside=1,nside_max_nelm(nelm0)
       node1=node_nside_nelm(nside,nelm0)
       IF(nside.EQ.nside_max_nelm(nelm0)) THEN
          node2=node_nside_nelm(1,nelm0)
       ELSE
          node2=node_nside_nelm(nside+1,nelm0)
       END IF
       x1=xnode(node1)
       y1=ynode(node1)
       x2=xnode(node2)
       y2=ynode(node2)
       IF(xy_in_range(x-x0,y-y0,x1-x0,y1-y0,x2-x0,y2-y0)) GO TO 100
    END DO
    WRITE(6,'(A/I5,1P4E12.4)') &
         'XX fem_linear_interporate: Logical error: nelm,x0,y0,x,y=', &
         nelm,x0,y0,x,y
    STOP

100 CONTINUE

    ! evaluate values at node1 and node2

    sum_of_weight=0.D0
    sum_of_f_weight=0.D0
    DO nangl=1,nelm_max_node(node1)
       nelm1=nelm_nangle_node(nangl,node1)
       weight=1.D0/((x1-xcenter_nelm(nelm1))**2 &
                   +(y1-ycenter_nelm(nelm1))**2)
       sum_of_weight=sum_of_weight+weight
       sum_of_f_weight=sum_of_f_weight+f_nelm(nelm1)*weight
    END DO
    f1=sum_of_f_weight/sum_of_weight

    sum_of_weight=0.D0
    sum_of_f_weight=0.D0
    DO nangl=1,nelm_max_node(node2)
       nelm2=nelm_nangle_node(nangl,node2)
       weight=1.D0/((x2-xcenter_nelm(nelm2))**2 &
                   +(y2-ycenter_nelm(nelm2))**2)
       sum_of_weight=sum_of_weight+weight
       sum_of_f_weight=sum_of_f_weight+f_nelm(nelm2)*weight
    END DO
    f2=sum_of_f_weight/sum_of_weight

    f0=f_nelm(nelm0)
    cmat(1,1)=x1
    cmat(1,2)=y1
    cmat(1,3)=1.D0
    cmat(2,1)=x2
    cmat(2,2)=y2
    cmat(2,3)=1.D0
    cmat(3,1)=x0
    cmat(3,2)=y0
    cmat(3,3)=1.D0
    CALL INVMRD(cmat,3,3,ierr)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX fem_linear_interporate: INVMRD error: IERR=',ierr
       STOP
    END IF

    a=cmat(1,1)*f1+cmat(1,2)*f2+cmat(1,3)*f0
    b=cmat(2,1)*f1+cmat(2,2)*f2+cmat(2,3)*f0
    c=cmat(3,1)*f1+cmat(3,2)*f2+cmat(3,3)*f0
    f=a*x+b*y+c
    RETURN
  END SUBROUTINE fem_linear_interporate

END MODULE feminterpolate
