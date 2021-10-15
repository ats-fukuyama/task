!     $Id: fpbroadcast.f90,v 1.3 2013/01/14 16:48:26 fukuyama Exp $

MODULE fpbroadcast

CONTAINS

!     ***** BROADCAST EQ DATA *****

    SUBROUTINE fp_eq_broadcast

      USE libmpi
      USE libmtx
      INCLUDE '../eq/eqcomc.inc'
      INTEGER,DIMENSION(3):: idata
      REAL(RKIND),DIMENSION(11):: ddata
      REAL(RKIND),DIMENSION(:),ALLOCATABLE:: temp
      INTEGER:: nrg,nzg

      idata( 1)=NRGMAX
      idata( 2)=NZGMAX
      idata( 3)=NPSMAX
      CALL mtx_broadcast_integer(idata,3)
      NRGMAX=idata( 1)
      NZGMAX=idata( 2)
      NPSMAX=idata( 3)

      ddata( 1)=RR
      ddata( 2)=BB
      ddata( 3)=RIP
      ddata( 4)=RA
      ddata( 5)=RKAP
      ddata( 6)=RDLT
      ddata( 7)=RB
      ddata( 8)=RAXIS
      ddata( 9)=ZAXIS
      ddata(10)=PSIPA
      ddata(11)=PSI0
      CALL mtx_broadcast_real8(ddata,11)
      RR=   ddata( 1)
      BB=   ddata( 2)
      RIP=  ddata( 3)
      RA=   ddata( 4)
      RKAP= ddata( 5)
      RDLT= ddata( 6)
      RB=   ddata( 7)
      RAXIS=ddata( 8)
      ZAXIS=ddata( 9)
      PSIPA=ddata(10)
      PSI0= ddata(10)

      CALL mtx_broadcast_real8(RG,NRGMAX)
      CALL mtx_broadcast_real8(ZG,NRGMAX)
      ALLOCATE(temp(nrgmax*nzgmax))
      DO nzg=1,nzgmax
         DO nrg=1,nrgmax
            temp(nrgmax*(nzg-1)+nrg)=PSIRZ(nrg,nzg)
         ENDDO
      ENDDO
      CALL mtx_broadcast_real8(temp,NRGMAX*NZGMAX)
      DO nzg=1,nzgmax
         DO nrg=1,nrgmax
            PSIRZ(nrg,nzg)=temp(nrgmax*(nzg-1)+nrg)
         ENDDO
      ENDDO
      DO nzg=1,nzgmax
         DO nrg=1,nrgmax
            temp(nrgmax*(nzg-1)+nrg)=HJTRZ(nrg,nzg)
         ENDDO
      ENDDO
      CALL mtx_broadcast_real8(temp,NRGMAX*NZGMAX)
      DO nzg=1,nzgmax
         DO nrg=1,nrgmax
            HJTRZ(nrg,nzg)=temp(nrgmax*(nzg-1)+nrg)
         ENDDO
      ENDDO
      DEALLOCATE(temp)
      CALL mtx_broadcast_real8(PSIPS,NPSMAX)
      CALL mtx_broadcast_real8(PPPS,NPSMAX)
      CALL mtx_broadcast_real8(TTPS,NPSMAX)
      CALL mtx_broadcast_real8(TEPS,NPSMAX)
      CALL mtx_broadcast_real8(OMPS,NPSMAX)

      RETURN
    END SUBROUTINE fp_eq_broadcast

END MODULE fpbroadcast
