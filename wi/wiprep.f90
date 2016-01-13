! $Id$

MODULE wiprep

  PRIVATE
  PUBLIC wi_prep

CONTAINS

  SUBROUTINE wi_prep

    USE wicomm,ONLY: wi_allocate
    IMPLICIT NONE
    
    CALL init_fem    ! initiaize FEM coefficients
    CALL set_grid    ! setup grid and density profile
    CALL wi_allocate
    CALL set_profile ! setup density profile
    RETURN

  END SUBROUTINE wi_prep

!     *****  INITIALIZE D0,D1,D2,D3  *****

    SUBROUTINE init_fem

      USE wicomm,ONLY: ikind,D0,D1,D2,D3
      IMPLICIT NONE
      INTEGER(ikind):: L,M,N

      D0(0,0)=1.D0/3.D0
      D0(0,1)=1.D0/6.D0
      D0(1,0)=1.D0/6.D0
      D0(1,1)=1.D0/3.D0
      D1(0,0)=-0.5D0
      D1(0,1)=-0.5D0
      D1(1,0)=0.5D0
      D1(1,1)=0.5D0
      D2(0,0)=1.D0
      D2(0,1)=-1.D0
      D2(1,0)=-1.D0
      D2(1,1)=1.D0
      DO L=0,1
         DO M=0,1
            DO N=0,1
               D3(L,M,N)=1.D0/12.D0
            END DO
         END DO
      END DO
      D3(0,0,0)=1.D0/4.D0
      D3(1,1,1)=1.D0/4.D0
      RETURN
    END SUBROUTINE init_fem

!     *****  set_grid

    SUBROUTINE set_grid

      USE libgrf,ONLY: GRD1D
      USE wicomm,ONLY: &
           ikind,rkind,dxmin,dx0,xwmin,xwint,alfa,pn0,xmin,xmax,nxmax,nwmax, &
           xgrid,modelp,idebug
      IMPLICIT NONE

      REAL(rkind):: xres,x,factor,dx1,dx2,range
      INTEGER(ikind):: nx,nx1,nx2
      INTEGER(ikind),SAVE:: nxmax_save=0
      REAL(rkind),DIMENSION(:),ALLOCATABLE:: xid

      IF(dxmin .LE. 0.D0) THEN
         nxmax=(xmax-xmin)/dx0
      ELSE
         xres=LOG(pn0)    ! resonance position (omegape=omega)
         x=xmin
         nx=0
         DO WHILE (x.LT.xmax)
            factor=((x-xres)/xwmin)**2
            IF(FACTOR.GT.100) THEN
               dx1=dx0
            ELSE
               dx1=dx0*dxmin/(dx0*EXP(-FACTOR)+dxmin)
            ENDIF
            x=x+dx1
            factor=((x-xres)/xwmin)**2
            IF(FACTOR.GT.100) THEN
               dx2=dx0
            ELSE
               dx2=dx0*dxmin/(dx0*EXP(-FACTOR)+dxmin)
            ENDIF
            x=x+0.5D0*(dx2-dx1)
            nx=nx+1
         END DO
         nxmax=nx
      END IF

      IF(nxmax.NE.nxmax_save) THEN
         IF(ALLOCATED(xgrid)) DEALLOCATE(xgrid)
         ALLOCATE(xgrid(0:nxmax))
      END IF

      IF(dxmin .LE. 0.D0) THEN
         DO nx=0,nxmax
            xgrid(nx)=xmin+dx0*nx
         END DO
      ELSE
         xres=LOG(pn0)/alfa    ! resonance position (omegape=omega)
         x=xmin
         nx=0
         xgrid(nx)=x
         DO WHILE (x.LT.xmax)
            factor=((x-xres)/xwmin)**2
            IF(FACTOR.GT.100) THEN
               dx1=dx0
            ELSE
               dx1=dx0*dxmin/(dx0*EXP(-FACTOR)+dxmin)
            ENDIF
            x=x+dx1
            factor=((x-xres)/xwmin)**2
            IF(FACTOR.GT.100) THEN
               dx2=dx0
            ELSE
               dx2=dx0*dxmin/(dx0*EXP(-FACTOR)+dxmin)
            ENDIF
            x=x+0.5D0*(dx2-dx1)
            nx=nx+1
            xgrid(nx)=x
         END DO
      END IF

      SELECT CASE(MODELP)
      CASE(0,1)
         nwmax=1
      CASE(2)
         nx1=0
         nx2=0
         nwmax=0
         DO WHILE(nx2.LT.nxmax)
            range=xgrid(nx2)-xgrid(nx1)
            IF(nx2-nx1.GT.nwmax) nwmax=nx2-nx1
            IF(range.LT.xwint) THEN
               nx2=nx2+1
            ELSE
               nx1=nx1+1
            ENDIF
         END DO
      END SELECT

      RETURN
    END SUBROUTINE set_grid

!     *****  SET PROFILE  ***** 

    SUBROUTINE set_profile

      USE wicomm,ONLY: ikind,rkind,nxmax,xgrid,alfa,pn0,cwe,cwp
      USE wigcom
      IMPLICIT NONE

      INTEGER(ikind):: nx

      DO NX=0,NXMAX
         CWE(NX)=DEXP(-0.5D0*ALFA*xgrid(nx))
         CWP(NX)=PN0
      END DO
      RETURN
    END SUBROUTINE set_profile

  END MODULE wiprep
