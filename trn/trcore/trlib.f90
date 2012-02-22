MODULE trlib
  USE trcomm, ONLY: ikind, rkind

  PUBLIC

  CONTAINS

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

      subroutine mesh_convert_gtom0(datag,datam,nrmax)
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
