C     $Id$
C
C     ***** SAVE TASK/EQ DATA *****
C
      SUBROUTINE EQSAVE
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL FWOPEN(21,KNAMEQ,0,MODEFW,'EQ',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(21)
      WRITE(21) RR,BB,RIP
      WRITE(21) NRGMAX,NZGMAX
      WRITE(21) (RG(NRG),NRG=1,NRGMAX)
      WRITE(21) (ZG(NZG),NZG=1,NZGMAX)
      WRITE(21) ((PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      WRITE(21) NPSMAX
      WRITE(21) (PSIPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (PPPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (TTPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (TEPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (OMPS(NPS),NPS=1,NPSMAX)
C
      WRITE(21) NSGMAX,NTGMAX,NUGMAX,NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX
      WRITE(21) ((PSI(NSG,NTG),NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) ((DELPSI(NSG,NTG),NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) ((HJT(NSG,NTG),NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) RAXIS,ZAXIS,PSITA,PSIPA,PSI0
      WRITE(21) (PSIPNV(NRV),NRV=1,NRVMAX)
      WRITE(21) (PSIPV(NRV),NRV=1,NRVMAX)
      WRITE(21) (PSITV(NRV),NRV=1,NRVMAX)
      WRITE(21) (QPV(NRV),NRV=1,NRVMAX)
      WRITE(21) (TTV(NRV),NRV=1,NRVMAX)
      WRITE(21) RA,RKAP,RDLT,RB,FRBIN
      WRITE(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      WRITE(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      WRITE(21) PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
      WRITE(21) PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      WRITE(21) PROFR0,PROFR1,PROFR2
C      WRITE(21) PTS,PN0,HM
      WRITE(21) ((HJTRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      CLOSE(21)
C
C      WRITE(6,*) 'HJTRZ=',HJTRZ(10,10)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
C
C     ***** LOAD EQUILIBRIUM DATA *****
C
      SUBROUTINE EQLOAD(MODELG1,KNAMEQ1,IERR)
C
      USE eq_bpsd_mod
      INCLUDE '../eq/eqcomc.inc'
      CHARACTER KNAMEQ1*80
      INTEGER ierr
C
      MODELG=MODELG1
      KNAMEQ=KNAMEQ1
      CALL EQREAD(IERR)
      CALL eq_bpsd_set(ierr)
      RETURN
      END
C
C     ***** LOAD EQUILIBRIUM DATA *****
C
      SUBROUTINE EQREAD(IERR)
C
      USE eqread_mod
      INCLUDE '../eq/eqcomc.inc'
C
      IF(MODELG.EQ.3.OR.MODELG.EQ.9) THEN
         CALL EQRTSK(IERR)
      ELSEIF(MODELG.EQ.5) THEN
         CALL EQDSKR(IERR)
         CALL EQCALQ(IERR)
      ELSEIF(MODELG.EQ.8) THEN
         CALL EQJAEAR(IERR)
      ELSEIF(MODELG.EQ.15) THEN
         CALL EQDSK
      ELSE
         WRITE(6,*) 'XX EQLOAD: UNKNOWN MODELG: MODELG=',MODELG
      ENDIF
C
      RETURN
      END
C
C     ***** LOAD TASK/EQ DATA *****
C
      SUBROUTINE EQRTSK(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION DERIV(NRVM)
C
      CALL FROPEN(21,KNAMEQ,0,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
C
      READ(21) RR,BB,RIP
      READ(21) NRGMAX,NZGMAX
      READ(21) (RG(NRG),NRG=1,NRGMAX)
      READ(21) (ZG(NZG),NZG=1,NZGMAX)
      READ(21) ((PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      READ(21) NPSMAX
      READ(21) (PSIPS(NPS),NPS=1,NPSMAX)
      READ(21) (PPPS(NPS),NPS=1,NPSMAX)
      READ(21) (TTPS(NPS),NPS=1,NPSMAX)
      READ(21) (TEPS(NPS),NPS=1,NPSMAX)
      READ(21) (OMPS(NPS),NPS=1,NPSMAX)
C
      READ(21) NSGMAX,NTGMAX,NUGMAX,NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX
      READ(21) ((PSI(NTG,NSG),NTG=1,NTGMAX),NSG=1,NSGMAX)
      READ(21) ((DELPSI(NTG,NSG),NTG=1,NTGMAX),NSG=1,NSGMAX)
      READ(21) ((HJT(NTG,NSG),NTG=1,NTGMAX),NSG=1,NSGMAX)
      READ(21) RAXIS,ZAXIS,PSITA,PSIPA,PSI0
      READ(21) (PSIPNV(NRV),NRV=1,NRVMAX)
      READ(21) (PSIPV(NRV),NRV=1,NRVMAX)
      READ(21) (PSITV(NRV),NRV=1,NRVMAX)
      READ(21) (QPV(NRV),NRV=1,NRVMAX)
      READ(21) (TTV(NRV),NRV=1,NRVMAX)
      READ(21) RA,RKAP,RDLT,RB,FRBIN
      READ(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      READ(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      READ(21) PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
      READ(21) PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      READ(21) PROFR0,PROFR1,PROFR2
C      READ(21) PTS,PN0,HM
      READ(21,ERR=1000) ((HJTRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      GOTO 1001
 1000 CONTINUE
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            HJTRZ(NRG,NZG)=0.D0
         ENDDO
         ENDDO
 1001 CONTINUE
C      
      CLOSE(21)
C
      IF(MODELG.EQ.9) THEN
         CALL EQMESH
         CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
         CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
         CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
         CALL EQDEFB
C         DO NSG=1,NSGMAX
C            WRITE(6,'(A,I5)') 'PSI NSG=',NSG
C            WRITE(6,'(1P5E12.4)') (PSI(NTG,NSG),NTG=1,NTGMAX)
C         ENDDO
C         DO NSG=1,NSGMAX
C            WRITE(6,'(A,I5)') 'HJT NSG=',NSG
C            WRITE(6,'(1P5E12.4)') (HJT(NTG,NSG),NTG=1,NTGMAX)
C         ENDDO
C         WRITE(6,'(A,I5,1P4E12.4)')
C     &        ('NV:',NV,PSIPNV(NV),PSITV(NV),QPV(NV),TTV(NV),
C     &         NV=1,NRVMAX)
      ENDIF
C     
C     WRITE(6,*) 'HJTRZ=',HJTRZ(10,10)
C
      RETURN
      END
C
C     ***** READ EQDSK FORMAT FILE *****
C
      SUBROUTINE EQDSKR(IERR)
C
      USE libgrf
      INCLUDE '../eq/eqcomq.inc'
C
      character case(6)*10
      dimension rlim(NSUM),zlim(NSUM),pressw(NPSM),pwprim(NPSM),
     &          dmion(NSUM),rhovn(NSUM),ajtor(NPSM)
C
      ierr=0

      neqdsk=31
      CALL FROPEN(neqdsk,KNAMEQ,1,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
c
      REWIND(neqdsk)
      read (neqdsk,2000) (case(i),i=1,6),idum,NRGMAX,NZGMAX
      write(6,'(A,2I5)') 'NRGMAX,NZGMAX=',NRGMAX,NZGMAX
      NPSMAX=NRGMAX
      read (neqdsk,2020) rdim,zdim,Rctr,rleft,zmid
      write(6,'(A/1P5E12.4)') 'rdim,zdim,Rctr,rleft,zmid=',
     &                         rdim,zdim,Rctr,rleft,zmid
      read (neqdsk,2020) RAXIS,ZAXIS,PSI0,PSIA,Bctr
      write(6,'(A/1P5E12.4)') 'RAXIS,ZAXIS,PSI0,PSIA,Bctr=',
     &                         RAXIS,ZAXIS,PSI0,PSIA,Bctr
C
      read (neqdsk,2020) RIP,simag,xdum,rmaxis,xdum
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      read (neqdsk,2020) (TTPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (PPPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (TTDTTPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (DPPPS(i),i=1,NPSMAX)
      read (neqdsk,2020) ((PSIRZ(i,j),i=1,NRGMAX),j=1,NZGMAX)
      read (neqdsk,2020) (QQPS(i),i=1,NPSMAX)

C Radial and vertical mesh

      DR=rdim/(NRGMAX-1)
      DZ=zdim/(NZGMAX-1)
      DO NRG=1,NRGMAX
         RG(NRG)=rleft+DR*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=zmid-0.5D0*zdim+DZ*(NZG-1)
      ENDDO

C for negative Ip and negative BB
      IF(RIP.LT.0.D0) RIP=-RIP
      IF(Bctr.LT.0.D0) Bctr=-Bctr
      IF(TTPS(1).LT.0.D0) THEN
         DO NPS=1,NPSMAX
            TTPS(NPS)=-TTPS(NPS)
         ENDDO
      ENDIF

C  Change normalization

      RIP=RIP/1.D6

      PSI0=2.D0*PI*PSI0
      PSIA=2.D0*PI*PSIA
C
      DPS=(PSIA-PSI0)/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)+PSI0
      ENDDO
      PSIPA=PSIA-PSI0
      PSI0=PSI0-PSIA

      DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZ(NRG,NZG)=2.D0*PI*PSIRZ(NRG,NZG)
         ENDDO
      ENDDO
      DO i=1,NPSMAX
         TTPS(i)   =2.D0*PI*TTPS(i)
         TTDTTPS(i)=2.D0*PI*TTDTTPS(i)
         DPPPS(i)=DPPPS(i)/(2.D0*PI)
         DTTPS(i)  =TTDTTPS(i)/TTPS(i)
      ENDDO
      PSIA=0.D0
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         HJTRZ(NRG,NZG)=0.D0
      ENDDO
      ENDDO

      read (neqdsk,2022,END=1800) NSUMAX,limitr
      write(6,'(A,2I5)') 'NSUMAX,limitr=',NSUMAX,limitr
      read (neqdsk,2020) (RSU(i),ZSU(i),i=1,NSUMAX)
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)

      RSUMAX = Rctr
      RSUMIN = Rctr
      ZSUMAX = 0.d0
      ZSUMIN = 0.d0
      R_ZSUMAX = 0.d0
      R_ZSUMIN = 0.d0
      do i = 1, nsumax
         if(RSU(i) .GT. RSUMAX) RSUMAX = RSU(i)
         if(RSU(i) .LT. RSUMIN) RSUMIN = RSU(i)
         if(ZSU(i) .GT. ZSUMAX) then
            ZSUMAX   = ZSU(i)
            R_ZSUMAX = RSU(i)
         end if
         if(ZSU(i) .LT. ZSUMIN) then
            ZSUMIN   = ZSU(i)
            R_ZSUMIN = RSU(i)
         end if
      enddo

C *** The following variable defined in Tokamaks 3rd, Sec. 14.14 ***
      RR   = 0.5d0 * (RSUMAX - RSUMIN) + RSUMIN
      RA   = RR - RSUMIN
      IF(RBRA.LE.1.D0) THEN
         RB = 1.1D0*RA
      ELSE
         RB = RBRA * RA
      END IF
      RKAP = (ZSUMAX - ZSUMIN) / (RSUMAX - RSUMIN)
      RDLT = 0.5d0 * ((RR-R_ZSUMIN) + (RR-R_ZSUMAX)) / RA
      BB   = Bctr
      RIPX = RIP

      GO TO 1900
C
C     without surface data
C
 1800 CONTINUE

      NSUMAX = 0
      limitr = 0
      RR   = Rctr
      RB   = RG(NRGMAX)-Rctr
      RKAP = (ZG(NRGMAX)-ZG(1))/(RG(NRGMAX)-RG(1))
      RDLT = 0.D0
      BB   = Bctr
      RIPX = RIP

 1900 CONTINUE

      RETURN
      return
c     
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
c 2024 format (i5,e16.9,i5)
       end
C
C     ***** READ JAEA EQ FORMAT FILE *****
C
      SUBROUTINE EQJAEAR(IERR)
C
      USE libgrf
      INCLUDE '../eq/eqcomq.inc'
C
      integer(4):: ir,iz
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: psi_temp
      REAL(8),DIMENSION(:,:),ALLOCATABLE::  PSIRG,PSIZG,PSIRZG
      REAL(8),DIMENSION(:,:),ALLOCATABLE::  PSIx1,psix2,PSIx3
      REAL(8):: DERIV(NPSM)
      REAL(8),DIMENSION(1)::  THIN=(/0.035/)
      REAL(8),DIMENSION(:),ALLOCATABLE:: rc_xp,zc_xp,psic_xp
      REAL(8),DIMENSION(NSUM):: XA
      INTEGER:: icount,icountmax
      INTEGER,PARAMETER:: icountm=10
      INTEGER:: ic_min1,ic_min2,ic_min3
      REAL(8):: psic_max,psic_min
      EXTERNAL RBB,PSIGZ0
C
      ALLOCATE(PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM))
      ALLOCATE(PSIx1(NRGM,NZGM),psix2(NRGM,NZGM),PSIx3(NRGM,NZGM))

      neqdsk=21
      CALL FROPEN(neqdsk,KNAMEQ,0,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
c
      REWIND(neqdsk)
      READ (neqdsk) NRGMAX,NZGMAX
      write(6,*) 'nrgmax,nzgmax=',nrgmax,nzgmax
      READ (neqdsk) rmin,rmax,zmin,zmax
      write(6,*) 'rmin,rmax=',rmin,rmax
      write(6,*) 'zmin,zmax=',zmin,zmax
      DR=(rmax-rmin)/(NRGMAX-1)
      DZ=(zmax-zmin)/(NZGMAX-1)
      DO NRG=1,NRGMAX
         RG(NRG)=rmin+DR*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=zmin+DZ*(NZG-1)
      ENDDO

      ALLOCATE(psi_temp(NRGMAX,NZGMAX))
      READ(neqdsk) psi_temp
      psirz(1:nrgmax,1:nzgmax)=psi_temp(1:nrgmax,1:nzgmax)
      DEALLOCATE(psi_temp)

      read (neqdsk) BtR
      write(6,*) 'read BtR=',BtR
      CLOSE(neqdsk)

      neqdsk=22
      CALL FROPEN(neqdsk,KNAMEQ2,1,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
      write(6,*) 'open 2'
      npsmax=41
      REWIND(neqdsk)
      READ(neqdsk,*)
      DO nps=1,npsmax
         READ(neqdsk,*) i,ttps(nps),ppps(nps)
      ENDDO
      CLOSE(neqdsk)


! ----- calculate dpsirz/dr)**2+(dpsirz/dz)**2 -----

      psix1(1:nrgmax,1:nzgmax)=0.0
      DO nz=2,nzgmax-1
         DO nr=2,nrgmax-1
            psix1(nr,nz)=(psirz(nr+1,nz)-psirz(nr-1,nz))**2 
     &                  /(rg(nr+1)-rg(nr-1))**2 
     &                  +(psirz(nr,nz+1)-psirz(nr,nz-1))**2 
     &                  /(zg(nz+1)-zg(nz-1))**2
         END DO
      END DO

! ----- find local minimum -----
      ALLOCATE(rc_xp(icountm),zc_xp(icountm),psic_xp(icountm))
      icount=0
      DO nz=3,nzgmax-2
         DO nr=3,nrgmax-2
            IF((psix1(nr,nz) < psix1(nr-1,nz)) .AND. 
     &         (psix1(nr,nz) < psix1(nr+1,nz)) .AND.
     &         (psix1(nr,nz) < psix1(nr,nz-1)) .AND.
     &         (psix1(nr,nz) < psix1(nr,nz+1))) THEN
               icount=icount+1
               rc_xp(icount)=rg(nr)
               zc_xp(icount)=zg(nz)
               psic_xp(icount)=psirz(nr,nz)
               IF(icount == icountm) GOTO 1000
            END IF
         END DO
      END DO
 1000 CONTINUE
      icountmax=icount

      SELECT CASE(icountmax)
      CASE(0)
         WRITE(6,*) 'XX NO MAG AXIS FOUND'
         STOP
      CASE(1)
         RAXIS=rc_xp(1)
         ZAXIS=zc_xp(1)
         PSIAXIS=psic_xp(1)
         NXPOINT=0
         WRITE(6,'(A,1P3E12.4)') 'Axis:',RAXIS,ZAXIS,PSIAXIS
      CASE(2:icountm)
         ic_min1=0
         ic_min2=0
         ic_min3=0
         psic_max=psic_xp(1)
         DO icount=2,icountmax
            psic_max=MAX(psic_xp(icount),psic_max)
         END DO
         psic_min=psic_max
         ic_min1=0
         DO icount=1,icountmax
            IF(psic_xp(icount) < psic_min) THEN
               ic_min1=icount
               psic_min=psic_xp(icount)
            END IF
         END DO
         psic_min=psic_max
         ic_min2=0
         DO icount=1,icountmax
            IF((icount /= ic_min1).AND.
     &         (psic_xp(icount) < psic_min)) THEN
               ic_min2=icount
               psic_min=psic_xp(icount)
            END IF
         END DO
         psic_min=psic_max
         ic_min3=0
         DO icount=1,icountmax
            IF((icount /= ic_min1).AND.
     &         (icount /= ic_min2).AND.
     &         (psic_xp(icount) < psic_min)) THEN
               ic_min3=icount
               psic_min=psic_xp(icount)
            END IF
         END DO

         IF(ic_min1 /= 0) THEN
            RAXIS=rc_xp(ic_min1)
            ZAXIS=zc_xp(ic_min1)
            PSIAXIS=psic_xp(ic_min1)
            NXPOINT=0
            WRITE(6,'(A,1P3E12.4)') 'Axis:',RAXIS,ZAXIS,PSIAXIS
         ENDIF
         IF(ic_min2 /= 0) THEN
            RXPNT1=rc_xp(ic_min2)
            ZXPNT1=zc_xp(ic_min2)
            PSIXPNT1=psic_xp(ic_min2)
            NXPOINT=1
            WRITE(6,'(A,1P3E12.4)') 'Xp1: ',RXPNT1,ZXPNT1,PSIXPNT1
         ENDIF
         IF(ic_min3 /= 0) THEN
            RXPNT2=rc_xp(ic_min3)
            ZXPNT2=zc_xp(ic_min3)
            PSIXPNT2=psic_xp(ic_min3)
            NXPOINT=2
            WRITE(6,'(A,1P3E12.4)') 'Xp2: ',RXPNT2,ZXPNT2,PSIXPNT2
         ENDIF
      END SELECT

! -----
      CALL setup_psig
      CALL find_axis
      WRITE(6,'(A,1P3E12.4)') 'RAXIS,ZAXIS,PSI_AXIS=',
     &                         RAXIS,ZAXIS,PSIG(RAXIS,ZAXIS)
      IF(NXPOINT >= 1) THEN
         CALL find_xpoint1
         WRITE(6,'(A,1P3E12.4)') 'RXPNT1,ZXPNT1,PSI_XPNT1=',
     &                            RXPNT1,ZXPNT1,PSIG(RXPNT1,ZXPNT1)
      END IF
      IF(NXPOINT >= 2) THEN
         CALL find_xpoint2
         WRITE(6,'(A,1P3E12.4)') 'RXPNT2,ZXPNT2,PSI_XPNT2=',
     &                            RXPNT2,ZXPNT2,PSIG(RXPNT2,ZXPNT2)
      END IF

      REDGE=FBRENT(PSIGZ0,RAXIS+0.1D0,2*RAXIS,1.D-8)
      WRITE(6,'(A,1P3E12.4)') 'REDGE,ZAXIS,PSI_EDGE=',
     &                         REDGE,ZAXIS,PSIG(REDGE,ZAXIS)

      NMAX=400
      H=16.D0*(REDGE-RAXIS)/NMAX
      CALL calc_separtrix(REDGE,ZAXIS,RXPNT1,ZXPNT1,H,NMAX,
     &                    XA,RSU,ZSU,NSUMAX,IERR) 

!      CALL PAGES
!      CALL GRD2D(0,rg,zg,psirz,nrgm,nrgmax,nzgmax,'@psirz@',0,0,1,
!     &           NLMAX=31,ASPECT=0.D0,
!     &           LINE_thickness=THIN)
!      CALL SETRGB(1.0,0.0,0.0)
!      CALL SETLNW(0.035)
!      CALL draw_cross(RAXIS,ZAXIS,0.2D0)
!      IF(nxpoint >= 1) CALL draw_cross(RXPNT1,ZXPNT1,0.2D0)
!      IF(nxpoint >= 2) CALL draw_cross(RXPNT2,ZXPNT2,0.2D0)
!      IF(nxpoint >= 1) THEN
!         CALL SETRGB(0.0,0.0,0.0)
!         DO NSU=1,NSUMAX
!            CALL draw_cross(RSU(NSU),ZSU(NSU),0.1D0)
!C            WRITE(6,'(A,I5,1P3E12.4)') 'SU:',NSU,RSU(NSU),ZSU(NSU),
!C     &                                 PSIG(RSU(NSU),ZSU(NSU))
!         END DO
!      END IF
!      CALL PAGEE

      RSUMAX = RAXIS
      RSUMIN = RAXIS
      ZSUMAX = ZAXIS
      ZSUMIN = ZAXIS
      R_ZSUMAX = 0.d0
      R_ZSUMIN = 0.d0
      DO i = 1, nsumax
         IF(RSU(i) > RSUMAX) RSUMAX = RSU(i)
         IF(RSU(i) < RSUMIN) RSUMIN = RSU(i)
         IF(ZSU(i) > ZSUMAX) THEN
            ZSUMAX   = ZSU(i)
            R_ZSUMAX = RSU(i)
         END IF
         IF(ZSU(i) < ZSUMIN) THEN
            ZSUMIN   = ZSU(i)
            R_ZSUMIN = RSU(i)
         END IF
      END DO

C *** The following variable defined in Tokamaks 3rd, Sec. 14.14 ***
      RR = RAXIS
      RA = REDGE - RAXIS
      !==  RB: wall minor radius  ======================
      !    Multiplication factor 1.1 is tentatively set.
      RB   = 1.1d0 * RA
!      RB   = 1.2d0 * RA
      !=================================================
      RKAP = (ZSUMAX - ZSUMIN) / (RSUMAX - RSUMIN)

!  ---- corrected on 2010/01/18 for negative triangularity ----
!      RDLT = 0.5d0 * (ABS(RR-R_ZSUMIN) + ABS(RR-R_ZSUMAX)) / RA

      RDLT = 0.5d0 * ((RR-R_ZSUMIN) + (RR-R_ZSUMAX)) / RA
      BB   = BtR / RR

      write(6,'(A,1PE12.4)') 'RR    =',RR
      write(6,'(A,1PE12.4)') 'RA    =',RA
      write(6,'(A,1PE12.4)') 'RB    =',RB
      write(6,'(A,1PE12.4)') 'RKAP  =',RKAP
      write(6,'(A,1PE12.4)') 'RDLT  =',RDLT
      write(6,'(A,1PE12.4)') 'BB    =',BB

      PSI0 = 2.D0*PI*PSIG(RAXIS,ZAXIS)
      PSIA = 0.D0
      PSIPA=-PSI0
      DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZ(NRG,NZG)=2.D0*PI*PSIRZ(NRG,NZG)
         ENDDO
      ENDDO
      call setup_psig

      DPS=PSIPA/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)
         TTPS(NPS)=2.D0*PI*TTPS(NPS)
      ENDDO

      TTDTTPS(1:npsmax)=0.D0
      DPPPS(1:npsmax)=0.D0

      CALL SPL1D(PSIPS,PPPS,  DERIV,UPPPS, NPSMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PPPS: IERR=',IERR
      CALL SPL1D(PSIPS,TTPS,  DERIV,UTTPS, NPSMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTPS: IERR=',IERR
C
C      DO NPS=1,NPSMAX
C         X=PPFUNC(PSIPS(NPS))
C         WRITE(6,'(A,I5,1P5E12.4)') 
C     &        'NPS:',NPS,PSIPS(NPS),PPPS(NPS),TTPS(NPS),
C     &        PPFUNC(PSIPS(NPS)),TTFUNC(PSIPS(NPS))
C      END DO
C
      CALL EQCALQP(IERR)
      IF(IERR.NE.0) RETURN
C
      IF(NSUMAX.GT.0) THEN
         CALL EQCALQV(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF
C
      CALL EQSETS_RHO(IERR)
      CALL EQSETS(IERR)

      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         HJTRZ(NRG,NZG)=0.D0
      ENDDO
      ENDDO
      DEALLOCATE(PSIRG,PSIZG,PSIRZG)
      DEALLOCATE(PSIx1,psix2,PSIx3)
      RETURN
      END

      SUBROUTINE RBB(X,RGB)
      REAL(4),INTENT(IN):: X
      REAL(4),INTENT(OUT):: RGB(3)

      RGB(1)=0.0
      RGB(2)=0.0
      RGB(3)=0.0

      IF(X > 0.5) THEN
         RGB(1)=SQRT(2*(X-0.5))
      ELSE IF(X < 0.5) THEN
         RGB(3)=SQRT(2*(0.5-X))
      ENDIF
      RETURN
      END
C
C     ***** SAVE METRICS *****
C
      SUBROUTINE EQMETRIC(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      character KNAMET*80
      data KNAMET /'eq_metric.dat'/ 
C
      nmetric=21
      CALL FWOPEN(nmetric,KNAMET,1,MODEFW,'EQ',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(nmetric)
      WRITE(nmetric,'(A,2X,A,6X,A,8X,A,3X,A)') '#','rho_tor','dV/drho',
     &     '<1/R^2>',"<|grad rho|^2/R^2>"
      DO NR = 1, NRPMAX
         CALL SPL1DF(FNPSIP(RHOT(NR)),DAT1,PSIP,UDVDRHO ,NRMAX,IERR)
         WRITE(nmetric,'(1X,F10.7,1P3E15.7)') RHOT(NR),DAT1,
     &        fnavir2(rhot(nr)),fnavgrr2(nr)
      ENDDO
      WRITE(nmetric,'(80X)')
      WRITE(nmetric,'(A,2X,A,3X,A,3X,A,5X,A,9X,A)') '#','rho_tor',
     &     "<|grad rho|>","<|grad rho|^2>","<B^2>","<1/B^2>"
      WRITE(nmetric,'(1X,0PF10.7,1P4E15.7)') (RHOT(NR),
     &     fnavgr(rhot(nr)),fnavgr2(rhot(nr)),
     &     fnavbb2(rhot(nr)),fnavib2(rhot(nr)),NR=1,NRPMAX)
      CLOSE(nmetric)
C
      WRITE(6,*) '# METRIC DATA WAS SUCCESSFULLY SAVED TO "',
     &           KNAMET(1:13),'".'
C
      RETURN
      END
C
C     ***** READ RIPPLE CONTOUR DATA FROM OFMC ****
C
      subroutine read_rppl(ierr)
C
      INCLUDE '../eq/eqcomq.inc'
C
      character kfile*20, kline*130
C
      kfile='ripple.profile'
      nrppl=21
      CALL FROPEN(nrppl,kfile,1,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
c
      rewind(nrppl)
c
c     *** R-coordinates ***
      idx = 0
      do
         if(idx == 0) then
            read(nrppl,'(A130)',iostat=ist) kline
            if(index(kline,"R-coordinate") /= 0) then ! detect the start position of the data chunk
               idx = 1
            end if
            cycle
         end if
c
         read(nrppl,'(1x,10e13.5)') (Rrp(i),i=1,NRrpM)
         exit
      end do
c
c     *** Z-coordinates ***
      idx = 0
      do
         if(idx == 0) then
            read(nrppl,'(A130)',iostat=ist) kline
            if(index(kline,"Z-coordinate") /= 0) then ! detect the start position of the data chunk
               idx = 1
            end if
            cycle
         end if
c
         read(nrppl,'(1x,10e13.5)') (Zrp(i),i=1,NZrpM)
         exit
      end do
c
c     *** Ripple contour ***
      idx = 0
      j   = 0
      do
         if(idx == 0) then
            read(nrppl,'(A130)',iostat=ist) kline
            if(index(kline,"at") /= 0) then ! detect the start position of the data chunk
               idx = 1
               j = j + 1
            end if
            cycle
         end if
c
         read(nrppl,'(1x,10e13.5)') (RpplRZ(i,j),i=1,NRrpM)
         idx = 0
         if(j == NZrpM) then
            exit
         else
            cycle
         end if
      end do
c
      close(nrppl)
c
      return
      end

      SUBROUTINE draw_cross(x,y,len)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: x,y,len
      INTERFACE
         FUNCTION GUCLIP(x)
            REAL(8),INTENT(IN):: x
            REAL(4):: GUCLIP
         END FUNCTION GUCLIP
      END INTERFACE

      CALL MOVE2D(GUCLIP(x-0.5D0*len),GUCLIP(y))
      CALL DRAW2D(GUCLIP(x+0.5D0*len),GUCLIP(y))
      CALL MOVE2D(GUCLIP(x),GUCLIP(y-0.5D0*len))
      CALL DRAW2D(GUCLIP(x),GUCLIP(y+0.5D0*len))
      RETURN
      END
