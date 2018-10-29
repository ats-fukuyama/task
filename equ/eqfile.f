C     $Id$
C
C     ***** SAVE TASK/EQ DATA *****
C
      SUBROUTINE EQUSAVE(KNAMEQ,MODEP,IERR)
C
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use eqv_mod
      use eqt_mod
      use trn_mod
      use r2d_mod
      use com_mod
      use imp_mod
      use eqsub_mod
      use tpxssl_mod
      USE libfio
      implicit none
      character(len=80),intent(in)::  knameq
      integer,intent(in):: modep
      integer,intent(out):: ierr
      integer i,n
      integer NSGMAX,NTGMAX,NUGMAX,NRMAX,NTHMAX,NRVMAX,NTVMAX
      integer NSG,NTG,NRV
C
      CALL FWOPEN(21,KNAMEQ,0,MODEP,'EQU',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(21)
      WRITE(21) rmaj, btv/rmaj,tcu*1.D-6 !RR,BB,RIP
      WRITE(21) nr,nz                    !NRGMAX,NZGMAX
      WRITE(21) (rg(i),i=1,nr)
      WRITE(21) (zg(i),i=1,nz)
      WRITE(21) (psi(i)*(2.D0*cnpi),i=1,nr*nz)       
C                               !PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      WRITE(21) nv                       !NPSMAX
      WRITE(21) ((siv(n)-siv(1))*(2.D0*cnpi),n=1,nv) !(PSIPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (muv(n)/cnmu*hdv(n)**gam,n=1,nv) !(PPPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (hdv(n)/aav(n)*(2.D0*cnpi),n=1,nv) !(TTPS(NPS),NPS=1,NPSMAX)
      WRITE(6,'(I8,1pE12.4)') (n,hdv(n)/aav(n),n=1,nv) !(TTPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (0.d0,n=1,nv) !(TEPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (0.d0,n=1,nv) !(OMPS(NPS),NPS=1,NPSMAX)
C
      NSGMAX=64
      NTGMAX=64
      NUGMAX=64
      NRMAX =100
      NTHMAX=64
      NRVMAX=50
      NTVMAX=200
C
      WRITE(21) NSGMAX,NTGMAX,NUGMAX,NRMAX,NTHMAX,nsu,NRVMAX,NTVMAX
      WRITE(21) ((0.d0,NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) ((0.d0,NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) ((0.d0,NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) raxis,zaxis,hiv(nv),
     &     (siv(nv)-siv(1))*(2.D0*cnpi),(siv(1)-siv(nv)*(2.D0*cnpi)) 
C                                       !RAXIS,ZAXIS,PSITA,PSIPA,PSI0
      WRITE(21) (0.d0,NRV=1,NRVMAX)
      WRITE(21) (0.d0,NRV=1,NRVMAX)
      WRITE(21) (0.d0,NRV=1,NRVMAX)
      WRITE(21) (0.d0,NRV=1,NRVMAX)
      WRITE(21) (0.d0,NRV=1,NRVMAX)
      WRITE(21) rpla,elip,trig,rpla*1.1,1.D0    !RA,RKAP,RDLT,RB,FRBIN
      WRITE(21) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 !PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      WRITE(21) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 !PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      WRITE(21) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 !PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
      WRITE(21) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 !PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      WRITE(21) 0.d0,0.d0,0.d0 !PROFR0,PROFR1,PROFR2
C      WRITE(21) PTS,PN0,HM
      WRITE(21) 0.d0
      CLOSE(21)
C
C      WRITE(6,*) 'HJTRZ=',HJTRZ(10,10)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
