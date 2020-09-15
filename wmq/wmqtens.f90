subroutine wmqtens(CD,ne,OCE,OPE,OUH,OR,OL,W,dx,omega,B0,nu,RR,RA,q0,qa,n0)

  use bpsd
  implicit none

  INTEGER(4),INTENT(IN)  :: W
  COMPLEX(8),INTENT(OUT) :: CD(W,W,3,3)
  REAL(8)   ,INTENT(OUT) :: ne(W,W),OCE(W,W),OPE(W,W),OUH(W,W),OR(W,W),OL(W,W)
  REAL(8)   ,INTENT(IN)  :: dx,omega,B0,nu,RR,RA,q0,qa,n0

  COMPLEX(8) :: Comega,COEF
  integer(4) :: i,j
  REAL(8)    :: OCEX(W,W),OCEY(W,W),OCEZ(W,W)
  REAL(8)    :: BX(W,W),BY(W,W),BZ(W,W),r,rx,ry,q
  
  DO j=1,W
     DO i=1,W
       
        rx  = (DBLE(i)-DBLE(W)*0.5d0)*dx
        ry  = (DBLE(j)-DBLE(W)*0.5d0)*dx
        r   =  SQRT(rx**2+ry**2)

        IF(r.le.RA)THEN

           q=q0+(qa-q0)*(r/RA)**2

           BX(i,j)=-B0*ry/(RR*q)
           BY(i,j)= B0*rx/(RR*q)
           BZ(i,j)= B0*RR/(RR+rx)
           ne(i,j)= n0*(1.d0-(r/RA)**2)

        ELSE

           q=qa*RA/r

           BX(i,j)=-B0*ry/(RR*q)
           BY(i,j)= B0*rx/(RR*q)
           BZ(i,j)= B0*RR/(RR+rx)
           ne(i,j)= 0.0d0

        END IF

     END DO
  END DO

!  DO j=1,W
!    DO i=1,W
!        BX(i,j) = 0.d0
!        BY(i,j) = 0.d0
!        BZ(i,j) = B0
!        ne(i,j) = n0*DBLE(2*W-i-j)/DBLE(W-2)
!     END DO
!  END DO

  DO j=1,W
     DO i=1,W

        OPE  (i,j) = sqrt(ne(i,j)*AEE**2/(AME*EPS0))
        OCEX (i,j) = AEE*BX(i,j)/AME
        OCEY (i,j) = AEE*BY(i,j)/AME
        OCEZ (i,j) = AEE*BZ(i,j)/AME
        OCE  (i,j) = sqrt(OCEX(i,j)**2+OCEY(i,j)**2+OCEZ(i,j)**2)
        OUH  (i,j) = sqrt(OCE(i,j)**2+OPE(i,j)**2)
        OR   (i,j) = 0.5d0*(sqrt(OCE(i,j)**2+4.d0*OPE(i,j)**2)+ABS(OCE(i,j)))
        OL   (i,j) = 0.5d0*(sqrt(OCE(i,j)**2+4.d0*OPE(i,j)**2)-ABS(OCE(i,j)))

     END DO
  END DO

  Comega=-CI*omega+nu !collision version of omega

  DO j=1,W
     DO i=1,W
        
        COEF        = ne(i,j)*AME**2/(AME*Comega*(OCE(i,j)**2+Comega**2))

        CD(i,j,1,1) = COEF*(    Comega**2    +OCEX(i,j)**2       ) 
        CD(i,j,1,2) = COEF*( OCEZ(i,j)*Comega+OCEX(i,j)*OCEY(i,j))
        CD(i,j,1,3) = COEF*(-OCEY(i,j)*Comega+OCEZ(i,j)*OCEX(i,j))
        CD(i,j,2,1) = COEF*(-OCEZ(i,j)*Comega+OCEX(i,j)*OCEY(i,j))
        CD(i,j,2,2) = COEF*(    Comega**2    +OCEY(i,j)**2       )
        CD(i,j,2,3) = COEF*( OCEX(i,j)*Comega+OCEY(i,j)*OCEZ(i,j))
        CD(i,j,3,1) = COEF*( OCEY(i,j)*Comega+OCEZ(i,j)*OCEX(i,j))
        CD(i,j,3,2) = COEF*(-OCEX(i,j)*Comega+OCEY(i,j)*OCEZ(i,j))
        CD(i,j,3,3) = COEF*(    Comega**2    +OCEZ(i,j)**2       )

    END DO
  END DO

  RETURN
END SUBROUTINE wmqtens
