MODULE trlib
  USE trcomm, ONLY: ikind, rkind

  PUBLIC

CONTAINS

  SUBROUTINE lin_itp(x0,y0,ax,ay,nxmax,nmax)
! -------------------------------------------------------------------------
!   *** Linear interpolation                                               
! -------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(4),               INTENT(IN) :: nmax,nxmax
    REAL(8),                  INTENT(IN) :: x0
    REAL(8),DIMENSION(1:nmax),INTENT(IN) :: ax,ay
    REAL(8),INTENT(OUT) :: y0

    INTEGER(4) :: nx, xid, nxl, nxr
    REAL(8)    :: grd

    IF(nxmax == 1)THEN
       y0 = ay(1)
       RETURN
    END IF

    IF(x0 < ax(1))THEN
       nxl = 0
       nxr = 1
    ELSE IF(x0 > ax(nxmax))THEN
       nxl = nxmax-1
       nxr = nxmax
    ELSE
       DO nx = 2, nxmax
          IF(x0 < ax(nx))THEN
             nxl = nx-1
             nxr = nx
             EXIT
          ELSE IF(x0 == ax(nx))THEN
             y0 = ay(nxmax)
             RETURN
          END IF
       END DO
    END IF

    grd = (ay(nxr) - ay(nxl)) / (ax(nxr) - ax(nxl))
    y0  = grd * (x0 - ax(nxl)) + ay(nxl)

    RETURN
  END SUBROUTINE lin_itp

! **************************************************************************

  SUBROUTINE mesh_rem2g(rhog,datam,datag,nrmax,idcntr,idedge)
! -------------------------------------------------------------------------
!   REconvert from HALF integer grid to integer grid.
! -------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind),INTENT(IN) :: nrmax, idcntr, idedge
    REAL(rkind),DIMENSION(0:nrmax),INTENT(IN)  :: rhog
    REAL(rkind),DIMENSION(1:nrmax),INTENT(IN)  :: datam
    REAL(rkind),DIMENSION(0:nrmax),INTENT(OUT) :: datag

    INTEGER(ikind) :: nr
    REAL(rkind)    :: cm, cp, rhom1,rhom2, FCTR ! TASK/lib
    REAL(rkind),DIMENSION(1:nrmax) :: drhog
    
    drhog(1:nrmax) = rhog(1:nrmax) - rhog(0:nrmax-1)
    ! This is the assumption of the subroutine.
    rhom1 = 0.5d0*(rhog(0) + rhog(1))
    rhom2 = 0.5d0*(rhog(1) + rhog(2))

    DO nr = 1, nrmax-1
       cm = drhog(nr+1) / (drhog(nr)+drhog(nr+1))
       cp = drhog(nr  ) / (drhog(nr)+drhog(nr+1))

       datag(nr) = cm*datam(nr) + cp*datam(nr+1)
    END DO

    SELECT CASE(idcntr)
    CASE(0)
       datag(0) = 0.d0
    CASE(1) ! derivative = 0 at the axis
       datag(1) = FCTR(rhom1,rhom2,datam(1),datam(2))
    END SELECT

    SELECT CASE(idedge)
    CASE(0)
       datag(nrmax) = 0.d0
    CASE(1) ! linear extrapolation
       datag(nrmax) = ((2.d0*drhog(nrmax)+drhog(nrmax-1))*datam(nrmax  )  &
                                            -drhog(nrmax)*datam(nrmax-1)) &
                     /(drhog(nrmax)+drhog(nrmax-1))
    END SELECT

    RETURN
  END SUBROUTINE mesh_rem2g


  REAL(rkind) FUNCTION FEDG(Rm,Rp,Fm,Fp)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: Rm, Rp, Fm, Fp

    FEDG = Fp + (Fp - Fm)/(Rp - Rm)*(1.d0 - Rp)

    RETURN
  END FUNCTION FEDG

! **************************************************************************
! **************************************************************************

  SUBROUTINE mesh_convert_mtog(datam,datag,nrmax)
!----------------------------------------------------------------------------
!          *** convert half mesh to origin + grid mesh ***
!
!  Suppose that rho derivative of data be zero at the axis and the boundary.
!  datag(0) is at rho = 0, and datag(nrmax) is at rho = 1.
!----------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind), INTENT(in)  :: nrmax
    REAL(rkind),    INTENT(in)  :: datam(1:nrmax)
    REAL(rkind),    INTENT(out) :: datag(0:nrmax)
      
    datag(0)       = (9.d0*datam(1)-datam(2))/8.d0
    datag(1:nrmax-1) = 0.5d0 * (datam(1:nrmax-1) + datam(2:nrmax))
    datag(nrmax)   = (4.d0*datam(nrmax)-datam(nrmax-1))/3.d0

    RETURN
  END SUBROUTINE mesh_convert_mtog

  SUBROUTINE mesh_convert_gtom(datag,datam,nrmax)
!----------------------------------------------------------------------------
!          *** convert origin + grid mesh to half mesh ***
!                ( just invert mesh_convert_gtom )
!
!  << CAUTION >>
!  This routine should be used only in case that one reconverts       
!   data converted by "mesh_convert_mtog" routine.                    
!  In any other case, one should use "data_interpolate_gtom" routine. 
!----------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind), INTENT(in)  :: nrmax
    REAL(rkind),    INTENT(in)  :: datag(0:nrmax)
    REAL(rkind),    INTENT(out) :: datam(1:nrmax)
    REAL(rkind) :: c11=9.d0/8.d0,c12=-1.d0/8.d0,c21=0.5d0,c22=0.5d0
    REAL(rkind) :: det,a11,a12,a21,a22

    det = c11*c22-c12*c21
    a11 = c22/det
    a12 =-c12/det
    a21 =-c21/det
    a22 = c11/det
    datam(1) = a11*datag(0) + a12*datag(1)
    datam(2) = a21*datag(0) + a22*datag(1)
    datam(3:nrmax) = 2.d0*datag(2:nrmax-1) - datam(2:nrmax-1)
    RETURN
  END SUBROUTINE mesh_convert_gtom

!****************************************************************************

  SUBROUTINE mesh_convert_mtog0(datam,datag,nrmax)
!----------------------------------------------------------------------------
!          *** convert half mesh to origin + grid mesh ***
!
! Suppose that data be zero at the axis, 
!  and rho derivative of data be zero at the boundary.
!----------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind), INTENT(in)  :: nrmax
    REAL(rkind),    INTENT(in)  :: datam(nrmax)
    REAL(rkind),    INTENT(out) :: datag(nrmax+1)

    datag(1)       = 0.d0
    datag(2:nrmax) = 0.5d0 * (datam(1:nrmax-1) + datam(2:nrmax))
    datag(nrmax+1) = (4.d0*datam(nrmax)-datam(nrmax-1))/3.d0

    RETURN
  END SUBROUTINE mesh_convert_mtog0

  SUBROUTINE mesh_convert_gtom0(datag,datam,nrmax)
!----------------------------------------------------------------------------
!           *** convert origin + grid mesh to half mesh ***
!                 ( just invert mesh_convert_mtog0 )
!
! << CAUTION >>
!  This routine should be used only in case that one reconverts        
!   data converted by "mesh_convert_mtog" routine.                     
!  In any other case, one should use "data_interpolate_gtom0" routine. 
!----------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind), INTENT(in)  :: nrmax
    REAL(rkind),    INTENT(in)  :: datag(nrmax+1)
    REAL(rkind),    INTENT(out) :: datam(nrmax)
    REAL(rkind) :: c11=9.d0/8.d0,c12=-1.d0/8.d0,c21=0.5d0,c22=0.5d0
    REAL(rkind) :: det,a11,a12,a21,a22

    det=c11*c22-c12*c21
    a11= c22/det
    a12=-c12/det
    a21=-c21/det
    a22= c11/det
    datam(1)=0.5d0*datag(2)
    datam(2:nrmax) = 2.d0 * datag(2:nrmax) - datam(1:nrmax-1)
    RETURN
  END SUBROUTINE mesh_convert_gtom0

!****************************************************************************

  SUBROUTINE data_interpolate_gtom_full(datag,datam,nrmax)
!----------------------------------------------------------------------------
!      *** interpolate data on half mesh from full data on grid ***
!----------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind), INTENT(in)  :: nrmax
    REAL(rkind),    INTENT(in)  :: datag(nrmax+1)
    REAL(rkind),    INTENT(out) :: datam(nrmax)

    datam(1:nrmax) = 0.5d0 * (datag(1:nrmax) + datag(2:nrmax+1))

    RETURN
  END SUBROUTINE data_interpolate_gtom_full

  SUBROUTINE data_interpolate_gtom(datag,datam,nrmax)
!----------------------------------------------------------------------------
! *** interpolate data on half mesh from data on grid except the axis ***
!----------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(ikind), INTENT(in)  :: nrmax
    REAL(rkind),    INTENT(in)  :: datag(0:nrmax)
    REAL(rkind),    INTENT(out) :: datam(nrmax)

    ! linear extrapolation
    datam(1) = 1.5d0 * datag(1) - 0.5d0 * datag(2)

    datam(2:nrmax) = 0.5d0 * (datag(1:nrmax-1) + datag(2:nrmax))

    RETURN
  END SUBROUTINE data_interpolate_gtom

END MODULE trlib
