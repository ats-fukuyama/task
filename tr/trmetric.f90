! trmetric.f90

MODULE trmetric

  PRIVATE
  PUBLIC tr_set_metric,trgfrg

CONTAINS

!     ***********************************************************

!           SET GEOMETRICAL FACTOR

!     ***********************************************************

      SUBROUTINE tr_set_metric(ierr)

      USE trcomm, ONLY : modelg, nrmax, knameq, knameq2
      USE trbpsd, ONLY: tr_bpsd_init,tr_bpsd_put,tr_bpsd_get
      USE equnit_mod, ONLY: eq_parm,eq_prof,eq_calc,eq_load
      USE plvmec, ONLY: pl_vmec
      implicit none
      integer, intent(out):: ierr
      character(len=80):: line

      CALL trstgf
      CALL trgfrg

      if(modelg.eq.3.or.modelg.eq.5.or.modelg.eq.8) then
         write(line,'(A,I5)') 'nrmax=',nrmax+1
         call eq_parm(2,line,ierr)
         write(line,'(A,I5)') 'nthmax=',64
         call eq_parm(2,line,ierr)
         write(line,'(A,I5)') 'nsumax=',0
         call eq_parm(2,line,ierr)
         if(modelg.eq.8) then
            write(line,'(A,A,A,A)') 'knameq2=','"',TRIM(knameq2),'"'
            write(6,'(A,A)') 'line=',line
            call eq_parm(2,line,ierr)
         end if
         call eq_load(modelg,knameq,ierr) ! load eq data and calculate eq
         IF(ierr.NE.0) THEN
            WRITE(6,*) 'XX eq_load: ierr=',ierr
            RETURN
         ENDIF
         call tr_bpsd_get(ierr)  ! 
         if(ierr.ne.0) write(6,*) 'XX tr_bpsd_get: ierr=',ierr
!         call trgout
      elseif(modelg.eq.7) then
         call pl_vmec(knameq,ierr) ! load vmec data
         call tr_bpsd_get(ierr)  ! 
!         call trgout
      elseif(modelg.eq.9) then
         call eq_prof ! initial calculation of eq
         call eq_calc         ! recalculate eq
         call tr_bpsd_get(ierr)  ! 
!         call trgout
      endif

      return
    end subroutine tr_set_metric
      
!     ***********************************************************

!           SET GEOMETRIC FACTOR AT HALF MESH

!     ***********************************************************

      SUBROUTINE TRSTGF

      USE TRCOMM, ONLY : &
           ABRHO, ABRHOU, AR1RHO, AR1RHOU, AR2RHO, AR2RHOU, ARRHO, &
           ARRHOU, BB, BP, BPRHO, DVRHO, DVRHOU, EPSRHO, MDLUF, &
           MDPHIA, NRMAX, PHIA, PI, QP, QRHO, RA, RG, &
           RHOG, RHOM, RJCB, RKAP, RKPRHO, RKPRHOU, RM, RMJRHO, &
           RMJRHOU, RMNRHO, RMNRHOU, RR, TTRHO, TTRHOU, VOLAU, &
           ABVRHO, PVOLRHOG, PSURRHOG, ABB1RHO
      IMPLICIT NONE
      INTEGER(4) :: NR
      REAL(8)    :: RKAPS, RHO_A

      RKAPS=SQRT(RKAP)
      IF(MDLUF.NE.0) THEN
         DO NR=1,NRMAX
            TTRHO(NR)=TTRHOU(1,NR)
            DVRHO(NR)=DVRHOU(1,NR)
            ABRHO(NR)=ABRHOU(1,NR)
            ABVRHO(NR)=DVRHO(NR)**2*ABRHO(NR)
            ARRHO(NR)=ARRHOU(1,NR)
            AR1RHO(NR)=AR1RHOU(1,NR)
            AR2RHO(NR)=AR2RHOU(1,NR)
            RMJRHO(NR)=RMJRHOU(1,NR)
            RMNRHO(NR)=RMNRHOU(1,NR)
            RKPRHO(NR)=RKPRHOU(1,NR)
            IF(MDPHIA.EQ.0) THEN
!     define rho_a from phi_a data
               RHO_A=SQRT(PHIA/(PI*BB))
               RJCB(NR)=1.D0/RHO_A
               RHOM(NR)=RM(NR)*RHO_A
               RHOG(NR)=RG(NR)*RHO_A
            ELSE
               RHO_A=SQRT(VOLAU(1)/(2.D0*PI**2*RMJRHOU(1,NRMAX)))
               RJCB(NR)=1.D0/RHO_A
               RHOM(NR)=RM(NR)/RJCB(NR)
               RHOG(NR)=RG(NR)/RJCB(NR)
            ENDIF
            EPSRHO(NR)=RMNRHO(NR)/RMJRHO(NR)
         ENDDO
         CALL FLUX
      ELSE
         DO NR=1,NRMAX
            BPRHO(NR)=BP(NR)
            QRHO(NR)=QP(NR)

            TTRHO(NR)=BB*RR
            DVRHO(NR)=2.D0*PI*RKAP*RA*RA*2.D0*PI*RR*RM(NR)
            ABRHO(NR)=1.D0/(RKAPS*RA*RR)**2
            ABVRHO(NR)=DVRHO(NR)**2*ABRHO(NR)
            ARRHO(NR)=1.D0/RR**2
            AR1RHO(NR)=1.D0/(RKAPS*RA)
            AR2RHO(NR)=1.D0/(RKAPS*RA)**2
            RMJRHO(NR)=RR
            RMNRHO(NR)=RA*RG(NR)
            RKPRHO(NR)=RKAP
            RJCB(NR)=1.D0/(RKAPS*RA)
            RHOM(NR)=RM(NR)/RJCB(NR)
            RHOG(NR)=RG(NR)/RJCB(NR)
            EPSRHO(NR)=RMNRHO(NR)/RMJRHO(NR)
            ABB1RHO(NR)=BB*(1.D0+0.25D0*EPSRHO(NR)**2)
            PVOLRHOG(NR)=PI*RKAP*(RA*RG(NR))**2*2.D0*PI*RR
            PSURRHOG(NR)=PI*(RKAP+1.D0)*RA*RG(NR)*2.D0*PI*RR
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE TRSTGF


!     ***********************************************************

!           CALCULATING FLUX FROM SOURCE TERM FILES

!     ***********************************************************

      SUBROUTINE FLUX

      USE TRCOMM, ONLY : DR, DVRHO, MDLFLX, NRMAX, NSM, PZ, RGFLX, SNBU, SWLU
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8),DIMENSION(NRMAX)::  SALEL, SALIL


      IF(MDLFLX.EQ.0) THEN
         RGFLX(1:NRMAX,1:NSM)=0.D0
      ELSE
         DO NR=1,NRMAX
            SALEL(NR)=SNBU(1,NR,1)+SWLU(1,NR)/PZ(2)
            RGFLX(NR,1)=SUM(SALEL(1:NR)*DVRHO(1:NR))*DR
            SALIL(NR)=SNBU(1,NR,2)+SWLU(1,NR)
            RGFLX(NR,2)=SUM(SALIL(1:NR)*DVRHO(1:NR))*DR
            RGFLX(NR,3:NSM)=0.D0
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE FLUX

!     ***********************************************************

!           GEOMETRIC QUANTITIES AT GRID MESH

!     ***********************************************************

      SUBROUTINE TRGFRG

      USE TRCOMM, ONLY : &
           ABB2RHOG, ABRHO, ABRHOG, AIB2RHOG, AR1RHO, AR1RHOG, AR2RHO, &
           AR2RHOG, ARHBRHOG, ARRHO, ARRHOG, BB, DVRHO, DVRHOG, EPSRHO, &
           NRMAX, RG, RKPRHO, RKPRHOG, RM, TTRHO, TTRHOG, ABVRHOG
      IMPLICIT NONE
      INTEGER(4) :: NR
      REAL(8)    :: RGL

      DO NR=1,NRMAX-1
         AR1RHOG(NR)=0.5D0*(AR1RHO(NR)+AR1RHO(NR+1))
         AR2RHOG(NR)=0.5D0*(AR2RHO(NR)+AR2RHO(NR+1))
         RKPRHOG(NR)=0.5D0*(RKPRHO(NR)+RKPRHO(NR+1))
         TTRHOG (NR)=0.5D0*(TTRHO (NR)+TTRHO (NR+1))
         DVRHOG (NR)=0.5D0*(DVRHO (NR)+DVRHO (NR+1))
         ARRHOG (NR)=0.5D0*(ARRHO (NR)+ARRHO (NR+1))
         ABRHOG (NR)=0.5D0*(ABRHO (NR)+ABRHO (NR+1))

         ABVRHOG(NR)=DVRHOG(NR)**2*ABRHOG(NR)
         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
      ENDDO
      NR=NRMAX
         RGL=RG(NR)

         CALL AITKEN(RGL,AR1RHOG(NR),RM,AR1RHO,2,NRMAX)
         CALL AITKEN(RGL,AR2RHOG(NR),RM,AR2RHO,2,NRMAX)
         CALL AITKEN(RGL,RKPRHOG(NR),RM,RKPRHO,2,NRMAX)
         CALL AITKEN(RGL,TTRHOG (NR),RM,TTRHO ,2,NRMAX)
         CALL AITKEN(RGL,DVRHOG (NR),RM,DVRHO ,2,NRMAX)
         CALL AITKEN(RGL,ARRHOG (NR),RM,ARRHO ,2,NRMAX)
         CALL AITKEN(RGL,ABRHOG (NR),RM,ABRHO ,2,NRMAX)

         ABVRHOG(NR)=DVRHOG(NR)**2*ABRHOG(NR)
         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
         ARHBRHOG(NR)=AR2RHOG(NR)*AIB2RHOG(NR)

      RETURN
      END SUBROUTINE TRGFRG

    END MODULE trmetric
