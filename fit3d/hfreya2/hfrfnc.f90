!#####################################################################
      SUBROUTINE fcnv(x,n,fvec)
!#####################################################################

  USE hfrmod, ONLY : rr, zz, pphi
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: n
  REAL(8),INTENT(INOUT) :: x(n)
  REAL(8),INTENT(OUT):: fvec(n)
  REAL(8):: r,z,p

      IF(x(1).lt.0.0) THEN
         x(1)=-x(1)
         x(2)=x(2)+3.1415926535897
      ENDIF
!
!     x1=x(1)
!     x2=x(2)
!     x3=x(3)

      CALL rzpfun(x,r,z,p)

      fvec(1)=r-rr
      fvec(2)=z-zz
      fvec(3)=p-pphi

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE fcnj(x,n,fjac)
!#####################################################################

  USE hfrmod, ONLY : abmnum, cp, cr, cz, dp, dr, dz, ep, er, ez, loop, mnumbr,    &
    modmax, nnumbr, pnmn, psino, rnmn, znmn
  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: n
  REAL(8),INTENT(INOUT):: x(n)
  REAL(8),INTENT(OUT)  :: fjac(n,n)
  INTEGER(4):: i, ii, j, k, ldfjac
  REAL(8):: cosf, dpharm, drharm, dxabm, dzharm, pharm, rharm, sinf, x1, x2, x3,   &
    xabm, xx, zharm

      ldfjac = n

      IF(x(1).lt.0.0) THEN
         x(1)=-x(1)
         x(2)=x(2)+3.1415926535897
      ENDIF

      x1=x(1)
      x2=x(2)
      x3=x(3)

      IF(x1.gt.1.0) GO TO 500
      IF(x1.eq.0.0) THEN
        ii=0
        xx=0.
        GO TO 50
      ENDIF

      DO i=1,loop
        IF(x1.gt.psino(i-1).and.x1.le.psino(i)) THEN
          ii=i-1
          xx=x1-psino(i-1)
          GO TO 50
        ENDIF
      END DO

 50   CONTINUE

      DO i=1,ldfjac
        DO j=1,n
          fjac(i,j)=0.
        END DO
      END DO

      DO k=1,modmax
        xabm=x1**abmnum(k)
        dxabm=x1**(abmnum(k)-1)*abmnum(k)
        cosf=dcos(mnumbr(k)*x2-nnumbr(k)*x3)
        sinf=dsin(mnumbr(k)*x2-nnumbr(k)*x3)
        rharm=rnmn(k,ii)+(cr(k,ii)+(dr(k,ii)+er(k,ii)*xx)*xx)*xx
        zharm=znmn(k,ii)+(cz(k,ii)+(dz(k,ii)+ez(k,ii)*xx)*xx)*xx
        pharm=pnmn(k,ii)+(cp(k,ii)+(dp(k,ii)+ep(k,ii)*xx)*xx)*xx
        drharm=cr(k,ii)+(2.d0*dr(k,ii)+3.d0*er(k,ii)*xx)*xx
        dzharm=cz(k,ii)+(2.d0*dz(k,ii)+3.d0*ez(k,ii)*xx)*xx
        dpharm=cp(k,ii)+(2.d0*dp(k,ii)+3.d0*ep(k,ii)*xx)*xx
        fjac(1,1)=fjac(1,1)+(drharm*xabm+rharm*dxabm)*cosf
        fjac(2,1)=fjac(2,1)+(dzharm*xabm+zharm*dxabm)*sinf
        fjac(3,1)=fjac(3,1)+(dpharm*xabm+pharm*dxabm)*sinf
        fjac(1,2)=fjac(1,2)-mnumbr(k)*rharm*xabm*sinf
        fjac(1,3)=fjac(1,3)+nnumbr(k)*rharm*xabm*sinf
        fjac(2,2)=fjac(2,2)+mnumbr(k)*zharm*xabm*cosf
        fjac(2,3)=fjac(2,3)-nnumbr(k)*zharm*xabm*cosf
        fjac(3,2)=fjac(3,2)+mnumbr(k)*pharm*xabm*cosf
        fjac(3,3)=fjac(3,3)-nnumbr(k)*pharm*xabm*cosf
      END DO
      fjac(3,3)=1.+fjac(3,3)
      RETURN

 500  CONTINUE
 600  FORMAT(1h //,' **** x **** out of range****',1p3e14.4)
 800  CONTINUE
      DO i=1,ldfjac
        DO j=1,n
          fjac(i,j)=0.
        END DO
      END DO
!
      DO k=1,modmax
        cosf=dcos(mnumbr(k)*x2-nnumbr(k)*x3)
        sinf=dsin(mnumbr(k)*x2-nnumbr(k)*x3)
        IF(k.eq.1) GO TO 910
        fjac(1,1)=fjac(1,1)+rnmn(k,loop)*cosf
        fjac(1,2)=fjac(1,2)-mnumbr(k)*rnmn(k,loop)*sinf*x1
        fjac(1,3)=fjac(1,3)+nnumbr(k)*rnmn(k,loop)*sinf*x1
 910    fjac(2,1)=fjac(2,1)+znmn(k,loop)*sinf
        fjac(2,2)=fjac(2,2)+mnumbr(k)*znmn(k,loop)*cosf*x1
        fjac(2,3)=fjac(2,3)-nnumbr(k)*znmn(k,loop)*cosf*x1
        fjac(3,2)=fjac(3,2)+mnumbr(k)*pnmn(k,loop)*cosf
        fjac(3,3)=fjac(3,3)-nnumbr(k)*pnmn(k,loop)*cosf
      END DO

      fjac(3,3)=fjac(3,3)+1.0
      RETURN
      END SUBROUTINE

!##################################################################
      REAL(8) FUNCTION ran3(idum)
!##################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(INOUT):: idum
  INTEGER(4),PARAMETER:: mseed=161803398, mbig=1000000000, mz=0
  REAL(8),PARAMETER:: fac=1.d-9
  INTEGER(4):: i, ii, inext, inextp, k, mj, mk, ma(55)
  INTEGER(4):: iff=0
!!!     parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
!      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)

      IF(idum.lt.0.or.iff.eq.0)THEN
        iff=1
        mj=mseed-iabs(idum)
        mj=mod(mj,mbig)
        ma(55)=mj
        mk=1
        DO i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          IF(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
        END DO
        DO k=1,4
          DO i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            IF(ma(i).lt.mz)ma(i)=ma(i)+mbig
          END DO
        END DO
        inext=0
        inextp=31
        idum=1
      ENDIF
      inext=inext+1
      IF(inext.eq.56)inext=1
      inextp=inextp+1
      IF(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      IF(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      RETURN
      END FUNCTION

!#####################################################################
      SUBROUTINE rzpfun(x,rr,zz,pp)

!     Find the position in the cylinder coordinates
!     form the position in the Boozer coordinates.

!--------------------------------------------------------------------c

  USE hfrmod, ONLY : abmnum, cp, cr, cz, dp, dr, dz, ep, er, ez, loop, mnumbr,    &
    modmax, nnumbr, pnmn, psino, rnmn, znmn
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: x(3)
  REAL(8),INTENT(OUT):: rr, zz, pp
  INTEGER(4):: i, ii, k
  REAL(8):: cosf, sinf, xabm, xx, xx2, xx3

      IF(x(1).lt.0.d0) THEN
!!         iflag=0                             200710
      ELSE IF(x(1).eq.0.0d0) THEN
        pp=x(3)
        rr=rnmn(1,0)
        zz=0.0d0
      ELSE IF(x(1).gt.1.0d0) THEN
        rr=rnmn(1,loop)
        zz=0.0d0
        pp=x(3)
        DO k=1,modmax
          cosf=dcos(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          sinf=dsin(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          zz=zz+znmn(k,loop)*sinf*x(1)
          pp=pp+pnmn(k,loop)*sinf
          rr=rr+rnmn(k,loop)*cosf*x(1)
        END DO
        rr=rr-rnmn(1,loop)*dcos(mnumbr(1)*x(2)-nnumbr(1)*x(3))*x(1)
      ELSE
        DO i=1,loop
          IF(x(1).gt.psino(i-1).and.x(1).le.psino(i)) THEN
            ii=i-1
            xx=x(1)-psino(ii)
            GO TO 50
          ENDIF
        END DO
 50     CONTINUE

        rr=0.0d0
        zz=0.0d0
        pp=x(3)

        xx2=xx*xx
        xx3=xx2*xx
        DO k=1,modmax
          xabm=x(1)**abmnum(k)
          cosf=dcos(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          sinf=dsin(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          rr=rr+(rnmn(k,ii)+cr(k,ii)*xx+dr(k,ii)*xx2+er(k,ii)*xx3)*xabm*cosf
          zz=zz+(znmn(k,ii)+cz(k,ii)*xx+dz(k,ii)*xx2+ez(k,ii)*xx3)*xabm*sinf
          pp=pp+(pnmn(k,ii)+cp(k,ii)*xx+dp(k,ii)*xx2+ep(k,ii)*xx3)*xabm*sinf
        END DO
      ENDIF

      RETURN
      END SUBROUTINE
