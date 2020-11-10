! ptexec1.f90

MODULE ptexec1

  PRIVATE
  PUBLIC pt_exec1

CONTAINS

  SUBROUTINE pt_exec1

    USE ptcomm
    USE libgrf
    IMPLICIT NONE
    REAL(rkind),ALLOCATABLE:: xc(:),yc(:),xe(:),ye(:),xd(:),yd(:)
    REAL(rkind):: dth,th,xmin,xmax,ymin,ymax
    REAL(4):: GUCLIP
    INTEGER:: nth

    ALLOCATE(xc(nthmax),yc(nthmax))
    ALLOCATE(xe(nthmax),ye(nthmax))
    ALLOCATE(xd(nthmax),yd(nthmax))
    
    dth=2.D0*Pi/(nthmax-1)

    DO nth=1,nthmax
       th=dth*(nth-1)
       xc(nth)=rr+ra*COS(th)
       yc(nth)=   ra*SIN(th)
       xe(nth)=rr+ra*COS(th)
       ye(nth)=   rkap*ra*SIN(th)
       xd(nth)=rr+ra*COS(th+rdlt*SIN(th))
       yd(nth)=   rkap*ra*SIN(th)
    END DO

    xmin=rr-1.1D0*ra
    xmax=rr+1.1D0*ra
    ymin=-1.1D0*rkap*ra
    ymax= 1.1D0*rkap*ra

    CALL PAGES
    CALL GRD2D_FRAME_START(0,xmin,xmax,ymin,ymax,'@coordinates@', &
                           ASPECT=0.D0,NOFRAME=1,NOINFO=1)
    CALL setlin(1,0,0)
    CALL setrgb(0.0,0.7,0.0)
    CALL setlnw(0.05)
    CALL MOVE2D(GUCLIP(xc(1)),GUCLIP(yc(1)))
    DO nth=2,nthmax
       CALL DRAW2D(GUCLIP(xc(nth)),GUCLIP(yc(nth)))
    END DO
    CALL setlin(2,0,0)
    CALL setrgb(0.0,0.5,1.0)
    CALL setlnw(0.05)
    CALL MOVE2D(GUCLIP(xe(1)),GUCLIP(ye(1)))
    DO nth=2,nthmax
       CALL DRAW2D(GUCLIP(xe(nth)),GUCLIP(ye(nth)))
    END DO
    CALL setlin(0,1,0)
    CALL setrgb(1.0,0.0,0.0)
    CALL setlnw(0.05)
    CALL MOVE2D(GUCLIP(xd(1)),GUCLIP(yd(1)))
    DO nth=2,nthmax
       CALL DRAW2D(GUCLIP(xd(nth)),GUCLIP(yd(nth)))
    END DO
    CALL setlin(0,0,0)
    CALL setrgb(0.0,0.0,0.0)
    CALL setlnw(0.03)
    CALL MOVE2D(GUCLIP(rr-ra),GUCLIP(0.D0))
    CALL DRAW2D(GUCLIP(rr+ra),GUCLIP(0.D0))
    CALL MOVE2D(GUCLIP(rr),GUCLIP(-rkap*ra))
    CALL DRAW2D(GUCLIP(rr),GUCLIP( rkap*ra))
    th=0.25D0*Pi
    CALL MOVEPT2D(GUCLIP(rr),GUCLIP(0.D0),1)
    CALL DRAWPT2D(GUCLIP(rr+ra*COS(th)),GUCLIP(ra*SIN(th)))
    CALL MOVEPT2D(GUCLIP(rr+ra*COS(th)),GUCLIP(0.D0),1)
    CALL DRAWPT2D(GUCLIP(rr+ra*COS(th)),GUCLIP(rkap*ra*SIN(th)))
    CALL MOVEPT2D(GUCLIP(rr),GUCLIP(0.D0),2)
    CALL DRAWPT2D(GUCLIP(rr+ra*COS(th)),GUCLIP(rkap*ra*SIN(th)))
    CALL DRAWPT2D(GUCLIP(rr+ra*COS(th+rdlt*SIN(th))),GUCLIP(rkap*ra*SIN(th)))
    CALL setlnw(0.05)
    CALL MOVE2D(GUCLIP(rr),GUCLIP(0.D0))
    CALL DRAW2D(GUCLIP(rr+ra*COS(th+rdlt*SIN(th))),GUCLIP(rkap*ra*SIN(th)))
    CALL GRD2D_FRAME_END
    CALL PAGEE

    CALL PAGES
    CALL GRD2D_FRAME_START(0,xmin,xmax,ymin,ymax,'@fixed boundary@', &
                           ASPECT=0.D0,NOFRAME=1,NOINFO=1)
    CALL setlin(0,1,0)
    CALL setrgb(1.0,0.0,0.0)
    CALL setlnw(0.05)
    CALL MOVE2D(GUCLIP(xd(1)),GUCLIP(yd(1)))
    DO nth=2,nthmax
       CALL DRAW2D(GUCLIP(xd(nth)),GUCLIP(yd(nth)))
    END DO
    CALL setlin(0,0,0)
    CALL setrgb(0.0,0.0,0.0)
    CALL setlnw(0.03)
    CALL MOVE2D(GUCLIP(rr-ra),GUCLIP(0.D0))
    CALL DRAW2D(GUCLIP(rr+ra),GUCLIP(0.D0))
    th=0.25D0*Pi
    CALL setlnw(0.05)
    CALL MOVE2D(GUCLIP(rr),GUCLIP(0.D0))
    CALL DRAW2D(GUCLIP(rr+ra*COS(th+rdlt*SIN(th))),GUCLIP(rkap*ra*SIN(th)))
    CALL GRD2D_FRAME_END
    CALL PAGEE

    RETURN
  END SUBROUTINE pt_exec1
END MODULE ptexec1
       
