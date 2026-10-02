! femcomm.f90

MODULE femcomm

  USE bpsd_kinds
  USE commpi
  
  IMPLICIT NONE

  PUBLIC

  ! === Basic mesh parameters ===  (Definition of a mesh)

  ! --- mesh type parameter ---  0:XY 1:RZ 2:ZR
  
  INTEGER:: mode_mesh=0

  ! --- mesh size parameters ---
  
  INTEGER:: node_max=0  ! number of nodes
  INTEGER:: nelm_max=0  ! number of element
  INTEGER,PARAMETER:: nside_max=3
                        ! maximum number of nodes in an element
                        ! maximum number of segments in an element
  ! --- node ---
  
  REAL(rkind),ALLOCATABLE:: xnode(:) ! x coordinate of node position (node_max)
  REAL(rkind),ALLOCATABLE:: ynode(:) ! y coordinate of node position (node_max)
  
  ! --- relation between nodes and elements ---
  
  INTEGER,ALLOCATABLE:: nside_max_nelm(:)
                      ! number of nodes/segments in the element
  INTEGER,ALLOCATABLE:: node_nside_nelm(:,:)
                      ! node number (=segment number) in the element 
                      !     (nside_max,nelm_max)

  ! === Variables calculated from the mesh ===

  ! --- mesh size variables ---
  
  INTEGER:: nseg_max=0  ! number of segments
  INTEGER:: nbdy_max=0  ! number of boundary segments

  ! --- node variables ---

  REAL(rkind):: xnode_min,xnode_max  ! minimum and maximum of xnode
  REAL(rkind):: ynode_min,ynode_max  ! minimum and maximum of xnode
  
  INTEGER,ALLOCATABLE:: mode_node(:) ! node prperty: 0:internal 2:boundary
  INTEGER,ALLOCATABLE:: nbdy_node(:) ! boundary number of node, 0 for internal

  INTEGER,ALLOCATABLE:: nvar_node(:) ! variable number of node E

  INTEGER,ALLOCATABLE:: nelm_max_node(:)
                                     ! number of elements allocated to a node
  INTEGER:: nelm_node_max            ! maximum of nelm_max_node  
  INTEGER,ALLOCATABLE:: ndir_max_node(:)
                                     ! number of elements related to a node 
  INTEGER,ALLOCATABLE:: nelm_ndir_node(:,:)
                                     ! list of elements related to a node
!  REAL(rkind),ALLOCATABLE:: weight_ndir_node(:,:)
!                                     ! weight of elements related to a node
!  INTEGER:: nseg_node_max              ! maximum of nseg_max_node
!  INTEGER,ALLOCATABLE:: nseg_max_node(:)
!                                       ! number of segments related to a node
!  INTEGER,ALLOCATABLE:: nseg_ndir_node(:,:)
!                                       ! list of segments related to a node     
  ! --- nelm variables ---

  REAL(rkind),ALLOCATABLE:: xcenter_nelm(:)  ! x coordinate of element center
  REAL(rkind),ALLOCATABLE:: ycenter_nelm(:)  ! y coordinate of element center

  INTEGER,ALLOCATABLE:: nseg_nside_nelm(:,:) ! segment number of the side
  INTEGER,ALLOCATABLE:: nelm1_nside_nelm(:,:)  ! element number for the side
  INTEGER,ALLOCATABLE:: nside1_nside_nelm(:,:) ! side number for the side
  
  REAL(rkind),ALLOCATABLE:: area_nelm(:)  ! area of element
  REAL(rkind),ALLOCATABLE:: vol_nelm(:)   ! volume of element 
  REAL(rkind):: area_tot ! total area
  REAL(rkind):: vol_tot  ! total volume

  REAL(rkind),ALLOCATABLE:: xmin_nelm(:),xmax_nelm(:) ! x range of element
  REAL(rkind),ALLOCATABLE:: ymin_nelm(:),ymax_nelm(:) ! y range of element

  INTEGER,ALLOCATABLE:: nmed_nelm(:) ! medium number of element
                                     !   0: vacuum, 1: plasma
  
  ! === segment variables ===
  
  REAL(rkind),ALLOCATABLE:: xcenter_nseg(:) ! x coordinate of segment
  REAL(rkind),ALLOCATABLE:: ycenter_nseg(:) ! y coordinate of segment
  REAL(rkind),ALLOCATABLE:: len_nseg(:)  ! length of segment
  INTEGER,ALLOCATABLE:: node_nseg(:,:) ! node_number of segment
  INTEGER,ALLOCATABLE:: nelm_nseg(:,:) ! element number of segment
  INTEGER,ALLOCATABLE:: nside_nseg(:,:) ! side number of segment
  INTEGER,ALLOCATABLE:: mode_nseg(:)   ! 0: inside, 1: exteranal boundary
  INTEGER,ALLOCATABLE:: nvar_nseg(:)   ! variable number of segment
        
  ! === boundary variables ===

  INTEGER,ALLOCATABLE:: nbdy_nseg(:)  ! boundary number of segment, 0:internal
  INTEGER,ALLOCATABLE:: nseg_nbdy(:)  ! segment number of boundary
  INTEGER,ALLOCATABLE:: nelm_nbdy(:)  ! element number of boundary
  INTEGER,ALLOCATABLE:: nside_nbdy(:) ! side number of boundary
  INTEGER,ALLOCATABLE:: node_nbdy(:)  ! node number of boundary

  ! === index variables for element sorting ===

  REAL(rkind),ALLOCATABLE :: sindex_nelm(:)
  INTEGER,ALLOCATABLE :: iv_nelm(:)
  INTEGER,ALLOCATABLE :: iw_nelm(:)
  INTEGER,ALLOCATABLE :: id_nelm(:)
  REAL(rkind),ALLOCATABLE :: sindex_min_nelm(:)
  REAL(rkind),ALLOCATABLE :: sindex_max_nelm(:)
        
  ! === surface volume ratio ===

  real(rkind),dimension(3,3,3):: AIF3,AIE3
  real(rkind),dimension(3,3)  :: AIF2,AIE2
  real(rkind),dimension(3)    :: AIF1,AIE1
        
  ! === allocation status ===

  LOGICAL:: fem_base_allocated=.FALSE.
  LOGICAL:: fem_mesh_allocated=.FALSE.

CONTAINS

  ! --- base allocation ---
  
  SUBROUTINE fem_base_allocate
    
    IMPLICIT NONE
    INTEGER,SAVE :: node_max_save=0
    INTEGER,SAVE :: nelm_max_save=0
    INTEGER,SAVE :: nside_max_save=0

    IF(fem_base_allocated) THEN
       IF((node_max.EQ.node_max_save).AND. &
          (nelm_max.EQ.nelm_max_save).AND. &
          (nside_max.EQ.nside_max_save)) RETURN
       CALL fem_base_deallocate
    END IF

    ALLOCATE(xnode(node_max),ynode(node_max))
    ALLOCATE(nside_max_nelm(nelm_max))
    ALLOCATE(node_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nelm1_nside_nelm(nside_max,nelm_max))
    ALLOCATE(nside1_nside_nelm(nside_max,nelm_max))

    node_max_save=node_max
    nelm_max_save=nelm_max
    nside_max_save=nside_max

    fem_base_allocated=.TRUE.
    
  END SUBROUTINE fem_base_allocate

  SUBROUTINE fem_base_deallocate

    IMPLICIT NONE
    
    DEALLOCATE(xnode,ynode)
    DEALLOCATE(nside_max_nelm)
    DEALLOCATE(node_nside_nelm)
    DEALLOCATE(nelm1_nside_nelm)
    DEALLOCATE(nside1_nside_nelm)

    fem_base_allocated=.FALSE.

  END SUBROUTINE fem_base_deallocate

  ! --- mesh allocation ---
  
  SUBROUTINE fem_mesh_allocate
    
    IMPLICIT NONE
    INTEGER,SAVE :: node_max_save=0
    INTEGER,SAVE :: nelm_max_save=0
    INTEGER,SAVE :: nseg_max_save=0
    INTEGER,SAVE :: nside_max_save=0
    INTEGER,SAVE :: nbdy_max_save=0

    IF(fem_mesh_allocated) THEN
       IF((node_max .EQ.node_max_save ).AND. &
          (nelm_max .EQ.nelm_max_save ).AND. &
          (nseg_max .EQ.nseg_max_save ).AND. &
          (nside_max.EQ.nside_max_save).AND. &
          (nbdy_max .EQ.nbdy_max_save )) RETURN
       CALL fem_mesh_deallocate
    END IF

    ALLOCATE(mode_node(node_max),nbdy_node(node_max))
    ALLOCATE(nvar_node(node_max))
    
    ALLOCATE(xcenter_nelm(nelm_max),ycenter_nelm(nelm_max))
    ALLOCATE(nseg_nside_nelm(nside_max,nelm_max))
    ALLOCATE(area_nelm(nelm_max),vol_nelm(nelm_max))
    ALLOCATE(xmin_nelm(nelm_max),xmax_nelm(nelm_max))
    ALLOCATE(ymin_nelm(nelm_max),ymax_nelm(nelm_max))
    ALLOCATE(nmed_nelm(nelm_max))

    ALLOCATE(xcenter_nseg(nseg_max),ycenter_nseg(nseg_max),len_nseg(nseg_max))
    ALLOCATE(node_nseg(2,nseg_max))
    ALLOCATE(nelm_nseg(2,nseg_max))
    ALLOCATE(nside_nseg(2,nseg_max))
    ALLOCATE(mode_nseg(nseg_max))
    ALLOCATE(nvar_nseg(nseg_max))

    ALLOCATE(nbdy_nseg(nseg_max))
    ALLOCATE(nseg_nbdy(nbdy_max))
    ALLOCATE(nelm_nbdy(nbdy_max))
    ALLOCATE(nside_nbdy(nbdy_max))
    ALLOCATE(node_nbdy(nbdy_max))

    ALLOCATE(sindex_nelm(nelm_max))
    ALLOCATE(iv_nelm(nelm_max),iw_nelm(nelm_max),id_nelm(nelm_max))
    ALLOCATE(sindex_min_nelm(nelm_max),sindex_max_nelm(nelm_max))

    IF(nrank.eq.0) WRITE(6,'(A)') '## fem_mesh_allocated'
    IF(nrank.eq.0) WRITE(6,'(A,4I8)') &
         '   node_max,nelm_max,nseg_max,nbdy_max=', &
             node_max,nelm_max,nseg_max,nbdy_max

    node_max_save =node_max
    nelm_max_save =nelm_max
    nseg_max_save =nseg_max
    nside_max_save=nside_max
    nbdy_max_save =nbdy_max

    fem_mesh_allocated=.TRUE.

  END SUBROUTINE fem_mesh_allocate

  SUBROUTINE fem_mesh_deallocate

    IMPLICIT NONE

    DEALLOCATE(mode_node,nbdy_node,nvar_node)
    
    DEALLOCATE(xcenter_nelm,ycenter_nelm)
    DEALLOCATE(nseg_nside_nelm)
    DEALLOCATE(area_nelm,vol_nelm)
    DEALLOCATE(xmin_nelm,xmax_nelm,ymin_nelm,ymax_nelm)
    DEALLOCATE(nmed_nelm)
    
    DEALLOCATE(xcenter_nseg,ycenter_nseg,len_nseg)
    DEALLOCATE(node_nseg,nelm_nseg,nside_nseg)
    DEALLOCATE(mode_nseg,nvar_nseg)
    
    DEALLOCATE(nbdy_nseg)
    DEALLOCATE(nseg_nbdy)
    DEALLOCATE(nelm_nbdy)
    DEALLOCATE(nside_nbdy)
    DEALLOCATE(node_nbdy)

    DEALLOCATE(sindex_nelm)
    DEALLOCATE(iv_nelm,iw_nelm,id_nelm)
    DEALLOCATE(sindex_min_nelm,sindex_max_nelm)

    fem_mesh_allocated=.FALSE.

    RETURN
  END SUBROUTINE fem_mesh_deallocate
END MODULE femcomm
