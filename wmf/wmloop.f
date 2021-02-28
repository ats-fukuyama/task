C
      SUBROUTINE WMLOOP(IERR,IFLG)
C
      INCLUDE 'wmcomm.inc'
      integer(4), intent(inout) :: IERR
      integer(4), intent(in), optional :: IFLG
      COMPLEX(8),DIMENSION(:),ALLOCATABLE:: CRADTTS
      REAL(8),DIMENSION(:),ALLOCATABLE:: PABSTTS,PCURTS
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: PABSTS,PCURRS
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE:: PABSRS
      REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE:: PCURS
      REAL(8),DIMENSION(:,:,:,:,:),ALLOCATABLE:: PABSS,PABSKS
      COMPLEX(8),DIMENSION(:,:,:,:,:),ALLOCATABLE:: 
     &     CEFLDKS,CBFLDKS,CEFLDS,CBFLDS,CENS,CEPS
      INTEGER(4):: NPH0_SAVE,NPHSMAX_SAVE
      INTEGER(4),DIMENSION(NPHSM):: NPH0S_SAVE
      REAL(8),DIMENSION(NPHSM):: PFRACS_SAVE
!mh      REAL(8),DIMENSION(NPHSM):: PFRAC
      REAL(8),DIMENSION(NPHSM):: PFACT1S
C
      NPH0_SAVE=NPH0
      NPHSMAX_SAVE=NPHSMAX
      DO NPHS=1,NPHSM
         NPH0S_SAVE(NPHS)=NPH0S(NPHS)
         PFRACS_SAVE(NPHS)=PFRACS(NPHS)
      ENDDO

      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/NHHMAX

      SELECT CASE(MDLWM_NPHS)
      CASE(0)
         NPHSMAX=1
         NPH0S(1)=NPH0
         PFRACS(1)=1.D0
      CASE(2:3)
         DO NPHS=1,NPHSMAX
            NPH0S(NPHS)=-NPHSMAX/2+(NPHS-1)
            PFRACS(NPHS)=1.D0
         ENDDO
      END SELECT

      ALLOCATE(CRADTTS(NPHSMAX))
      ALLOCATE(PABSTTS(NPHSMAX))
      ALLOCATE(PABSTS(NSMAX,NPHSMAX))
      ALLOCATE(PABSRS(NRMAX,NSMAX,NPHSMAX))
      ALLOCATE(PABSS(NTHMAX,NHHMAX,NRMAX,NSMAX,NPHSMAX))
      ALLOCATE(PABSKS(NTHMAX,NHHMAX,NRMAX,NSMAX,NPHSMAX))
      ALLOCATE(PCURTS(NPHSMAX))
      ALLOCATE(PCURRS(NRMAX,NPHSMAX))
      ALLOCATE(PCURS(NTHMAX,NHHMAX,NRMAX,NPHSMAX))
      ALLOCATE(CEFLDKS(3,NTHMAX,NHHMAX,NRMAX,NPHSMAX))
      ALLOCATE(CBFLDKS(3,NTHMAX,NHHMAX,NRMAX,NPHSMAX))
      ALLOCATE(CEFLDS(3,NTHMAX,NHHMAX,NRMAX,NPHSMAX))
      ALLOCATE(CBFLDS(3,NTHMAX,NHHMAX,NRMAX,NPHSMAX))
      ALLOCATE(CENS(3,NTHMAX,NHHMAX,NRMAX,NPHSMAX))
      ALLOCATE(CEPS(3,NTHMAX,NHHMAX,NRMAX,NPHSMAX))

!     ===== Main calculation loop =====

      DO NPHS=1,NPHSMAX

!     ----- calculate for NPHS -----

         NPH0=NPH0S(NPHS)
         CALL WMEXEC(IERR)
         IF(IERR.NE.0) GOTO 8000

!     ----- integrate power and current -----

         CALL WMDVOL

         PABSTT=0.D0
         DO NS=1,NSMAX
            PABST(NS)=0.D0
            DO NR=1,NRMAX
               PABSR(NR,NS)=0.D0
               DO NHH=1,NHHMAX
                  DO NTH=1,NTHMAX
                     PABSR(NR,NS)=PABSR(NR,NS)
     &                    +PABS(NTH,NHH,NR,NS)*DTH*DPH
                  END DO
               END DO
               PABST(NS)=PABST(NS)+PABSR(NR,NS)
            END DO
            PABSTT=PABSTT+PABST(NS)
         END DO
C
         PCURT=0.D0
         DO NR=1,NRMAX
            PCURR(NR)=0.D0
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PCURR(NR)=PCURR(NR)+PCUR(NTH,NHH,NR)*DTH
               END DO
            END DO
            PCURT=PCURT+PCURR(NR)
         END DO

!     ----- calculate power and current density -----

         DO NS=1,NSMAX
            DO NR=1,NRMAX
               PABSR(NR,NS)=PABSR(NR,NS)*DVOLS(NR)
               DO NHH=1,NHHMAX
                  DO NTH=1,NTHMAX
                     PABS(NTH,NHH,NR,NS)=PABS(NTH,NHH,NR,NS)
     &                                  *DVOL(NTH,NHH,NR)
                  END DO
               END DO
               DO ND=NDMIN,NDMAX
                  NDX=ND-NDMIN+1
                  DO MD=MDMIN,MDMAX
                     MDX=MD-MDMIN+1
                     PABSK(MDX,NDX,NR,NS)=PABSK(MDX,NDX,NR,NS)
     &                                   *DVOLS(NR)
                  END DO
               END DO
            END DO
         END DO
C
         DO NR=1,NRMAX
            PCURR(NR)=PCURR(NR)*DVOLS(NR)
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PCUR(NTH,NHH,NR)=PCUR(NTH,NHH,NR)*DVOL(NTH,NHH,NR)
               ENDDO
            ENDDO
         ENDDO

         IF(NRANK.EQ.0) THEN
            CALL WMPOUT
            IF(MODELW.EQ.1) CALL WMDOUT(IERR)
         ENDIF

!     ----- save data for NPHS ----

         CRADTTS(NPHS)=CRADTT
         PABSTTS(NPHS)=PABSTT
         DO NS=1,NSMAX
            PABSTS(NS,NPHS)=PABST(NS)
            DO NR=1,NRMAX
               PABSRS(NR,NS,NPHS)=PABSR(NR,NS)
               DO NHH=1,NHHMAX
                  DO NTH=1,NTHMAX
                     PABSS(NTH,NHH,NR,NS,NPHS)=PABS(NTH,NHH,NR,NS)
                     PABSKS(NTH,NHH,NR,NS,NPHS)=PABSK(NTH,NHH,NR,NS)
                  END DO
               END DO
            END DO
         END DO

         PCURTS(NPHS)=PCURT
         DO NR=1,NRMAX
            PCURRS(NR,NPHS)=PCURR(NR)
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PCURS(NTH,NHH,NR,NPHS)=PCUR(NTH,NHH,NR)
               END DO
            END DO
         END DO

         DO NR=1,NRMAX
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  DO I=1,3
                     CEFLDKS(I,NTH,NHH,NR,NPHS)=CEFLDK(I,NTH,NHH,NR)
                     CBFLDKS(I,NTH,NHH,NR,NPHS)=CBFLDK(I,NTH,NHH,NR)
                     CEFLDS(I,NTH,NHH,NR,NPHS)=CEFLD(I,NTH,NHH,NR)
                     CBFLDS(I,NTH,NHH,NR,NPHS)=CBFLD(I,NTH,NHH,NR)
                     CENS(I,NTH,NHH,NR,NPHS)=CEN(I,NTH,NHH,NR)
                     CEPS(I,NTH,NHH,NR,NPHS)=CEP(I,NTH,NHH,NR)
                  END DO
               END DO
            END DO
         END DO
      ENDDO
!     ===== END: Main calculation loop =====

!     ===== Post processing  =====

!     ----- calculate antenna weighting factor-----

      SELECT CASE(MDLWM_NPHS)
      CASE(0)
         PFACT1S(1)=1.D0
      CASE(1)
         DO NPHS=1,NPHSMAX
            IF(PABSTTS(NPHS).EQ.0.D0) THEN
               PFACT1S(NPHS)=0.D0
            ELSE
               PFACT1S(NPHS)=PFRACS(NPHS)/PABSTTS(NPHS)
            ENDIF
         END DO
      CASE(2)
         DO NPHS=1,NPHSMAX
            PFACT1S(NPHS)=1.D0
         END DO
      CASE(3)
         DO NPHS=1,NPHSMAX
            IF(ABS(CRADTTS(NPHS)).EQ.0.D0) THEN
               PFACT1S(NPHS)=0.D0
            ELSE
               PFACT1S(NPHS)=1.D0/ABS(CRADTTS(NPHS))**2
            ENDIF
         ENDDO
      END SELECT

!     ----- calculate total absorbed power ----

      PABSTOT=0.D0
      DO NPHS=1,NPHSMAX
         PABSTOT=PABSTOT+PFACT1S(NPHS)*PABSTTS(NPHS)
      ENDDO

!     ----- calculated power normalizing factor -----

      IF(PRFIN.EQ.0.D0.OR.PABSTOT.EQ.0.D0) THEN
         PFACT2=1.D0
      ELSE
         PFACT2=PRFIN/PABSTOT
      ENDIF

!     ----- normalize field quantities -----

      DO NPHS=1,NPHSMAX
         PWFACT=PFACT1S(NPHS)*PFACT2
         CEFACT=SQRT(PWFACT)

         CRADTTS(NPHS)=CEFACT*CRADTTS(NPHS)
!         write(6,*) NPHS,PWFACT,PABSTTS(NPHS)
         PABSTTS(NPHS)=PWFACT*PABSTTS(NPHS)
         DO NS=1,NSMAX
            PABSTS(NS,NPHS)=PWFACT*PABSTS(NS,NPHS)
            DO NR=1,NRMAX
               PABSRS(NR,NS,NPHS)=PWFACT*PABSRS(NR,NS,NPHS)
               DO NHH=1,NHHMAX
                  DO NTH=1,NTHMAX
                     PABSS(NTH,NHH,NR,NS,NPHS)
     &                    =PWFACT*PABSS(NTH,NHH,NR,NS,NPHS)
                     PABSKS(NTH,NHH,NR,NS,NPHS)
     &                    =PWFACT*PABSKS(NTH,NHH,NR,NS,NPHS)
                  END DO
               END DO
            END DO
         END DO

         PCURTS(NPHS)=PWFACT*PCURTS(NPHS)
         DO NR=1,NRMAX
            PCURRS(NR,NPHS)=PWFACT*PCURRS(NR,NPHS)
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PCURS(NTH,NHH,NR,NPHS)=PWFACT*PCURS(NTH,NHH,NR,NPHS)
               END DO
            END DO
         END DO

         DO NR=1,NRMAX
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  DO I=1,3
                     CEFLDKS(I,NTH,NHH,NR,NPHS)
     &                    =CEFACT*CEFLDKS(I,NTH,NHH,NR,NPHS)
                     CBFLDKS(I,NTH,NHH,NR,NPHS)
     &                    =CEFACT*CBFLDKS(I,NTH,NHH,NR,NPHS)
                     CEFLDS(I,NTH,NHH,NR,NPHS)
     &                    =CEFACT*CEFLDS(I,NTH,NHH,NR,NPHS)
                     CBFLDS(I,NTH,NHH,NR,NPHS)
     &                    =CEFACT*CBFLDS(I,NTH,NHH,NR,NPHS)
                     CENS(I,NTH,NHH,NR,NPHS)
     &                    =CEFACT*CENS(I,NTH,NHH,NR,NPHS)
                     CEPS(I,NTH,NHH,NR,NPHS)
     &                    =CEFACT*CEPS(I,NTH,NHH,NR,NPHS)
                  END DO
               END DO
            END DO
         END DO
      ENDDO

!     ----- calculate total field -----

      CRADTT=0.D0
      PABSTT=0.D0
      DO NS=1,NSMAX
         PABST(NS)=0.D0
         DO NR=1,NRMAX
            PABSR(NR,NS)=0.D0
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PABS(NTH,NHH,NR,NS)=0.D0
               END DO
            END DO
         END DO
      END DO

      DO NPHS=1,NPHSMAX
         CRADTT=CRADTT+CRADTTS(NPHS)
         PABSTT=PABSTT+PABSTTS(NPHS)
!         write(6,*) NPHS,PABSTT,PABSTTS(NPHS)
         DO NS=1,NSMAX
            PABST(NS)=PABST(NS)+PABSTS(NS,NPHS)
!            write(6,*) NS,NPHS,PABST(NS),PABSTS(NS,NPHS)
            DO NR=1,NRMAX
               PABSR(NR,NS)=PABSR(NR,NS)+PABSRS(NR,NS,NPHS)
               DO NHH=1,NHHMAX
                  DO NTH=1,NTHMAX
                     PABS(NTH,NHH,NR,NS)=PABS(NTH,NHH,NR,NS)
     &                    +PABSS(NTH,NHH,NR,NS,NPHS)
                     PABSK(NTH,NHH,NR,NS)=PABSK(NTH,NHH,NR,NS)
     &                    +PABSKS(NTH,NHH,NR,NS,NPHS)
                  END DO
               END DO
            END DO
         END DO
      END DO

      PCURT=0.D0
      DO NR=1,NRMAX
         PCURR(NR)=0.D0
         DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               PCUR(NTH,NHH,NR)=0.D0
            END DO
         END DO
      END DO

      DO NPHS=1,NPHSMAX
         PCURT=PCURT+PCURTS(NPHS)
         DO NR=1,NRMAX
            PCURR(NR)=PCURR(NR)+PCURRS(NR,NPHS)
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  PCUR(NTH,NHH,NR)=PCUR(NTH,NHH,NR)
     &                 +PCURS(NTH,NHH,NR,NPHS)
               END DO
            END DO
         END DO
      END DO

      DO NR=1,NRMAX
         DO NHH=1,NHHMAX
            DO NTH=1,NTHMAX
               DO I=1,3
                  CEFLDK(I,NTH,NHH,NR)=0.D0
                  CBFLDK(I,NTH,NHH,NR)=0.D0
                  CEFLD(I,NTH,NHH,NR)=0.D0
                  CBFLD(I,NTH,NHH,NR)=0.D0
                  CEN(I,NTH,NHH,NR)=0.D0
                  CEP(I,NTH,NHH,NR)=0.D0
               END DO
            END DO
         END DO
      END DO


      DO NPHS=1,NPHSMAX
         DO NR=1,NRMAX
            DO NHH=1,NHHMAX
               DO NTH=1,NTHMAX
                  DO I=1,3
                     CEFLDK(I,NTH,NHH,NR)=CEFLDK(I,NTH,NHH,NR)
     &                                   +CEFLDKS(I,NTH,NHH,NR,NPHS)
                     CBFLDK(I,NTH,NHH,NR)=CBFLDK(I,NTH,NHH,NR)
     &                                   +CBFLDKS(I,NTH,NHH,NR,NPHS)
                     CEFLD(I,NTH,NHH,NR)=CEFLD(I,NTH,NHH,NR)
     &                                  +CEFLDS(I,NTH,NHH,NR,NPHS)
                     CBFLD(I,NTH,NHH,NR)=CBFLD(I,NTH,NHH,NR)
     &                                  +CBFLDS(I,NTH,NHH,NR,NPHS)
                     CEN(I,NTH,NHH,NR)=CEN(I,NTH,NHH,NR)
     &                                +CENS(I,NTH,NHH,NR,NPHS)
                     CEP(I,NTH,NHH,NR)=CEP(I,NTH,NHH,NR)
     &                                +CEPS(I,NTH,NHH,NR,NPHS)
                  END DO
               END DO
            END DO
         END DO
      ENDDO

!     ----- Output P_abs(r,s) and J_CD(r) for TOPICS/ACCOME -----

      IF(present(IFLG)) CALL WM_TOPICS_OUT(PABSTS,IERR)

!     -----------------------------------------------------------

      IF(NRANK.EQ.0) THEN
         IF(PABSTT.NE.0.D0) THEN
            WRITE(6,*)
            WRITE(6,'(3A)') '  NPH','    PFRACS ', '   PABSTS(NS) [W]'
            DO NPHS=1,NPHSMAX
               WRITE(6,'(I5,1P6E12.4)') 
     &              NPH0S(NPHS),
     &              SUM(PABSTS(1:MIN(NSMAX,5),NPHS))/PABSTT,
     &              (PABSTS(NS,NPHS),NS=1,MIN(NSMAX,5))
            ENDDO
            WRITE(6,'(A)') 
     &           '--------------------------------------------'//
     &           '---------------------------------'
            WRITE(6,'(17X,1P5E12.4)') (PABST(NS),NS=1,MIN(NSMAX,5))
         ENDIF
         NPH0=0
         CALL WMPOUT
         IF(MODELW.EQ.1) CALL WMDOUT(IERR)
      ENDIF

 8000 CONTINUE

      DEALLOCATE(CEFLDKS)
      DEALLOCATE(CBFLDKS,CEFLDS,CBFLDS,CENS,CEPS)
      DEALLOCATE(PCURTS,PCURRS,PCURS)
      DEALLOCATE(PABSTTS,PABSTS,PABSRS,PABSS,PABSKS)
      DEALLOCATE(CRADTTS)

      NPH0=NPH0_SAVE
      NPHSMAX=NPHSMAX_SAVE
      DO NPHS=1,NPHSM
         NPH0S(NPHS)=NPH0S_SAVE(NPHS)
         PFRACS(NPHS)=PFRACS_SAVE(NPHS)
      ENDDO

      RETURN
      END

!     ----- calculate element volume and shell volume -----

      SUBROUTINE WMDVOL

      INCLUDE 'wmcomm.inc'
      real(8),dimension(3,3)::  gm
      real(8):: gj1,gj2

      DTH=2.D0*PI/NTHMAX
      DPH=2.D0*PI/NHHMAX

      DO NR=1,NRMAX-1
         DVOLS(NR)=0.D0
         DRHO=XRHO(NR+1)-XRHO(NR)
         DO NHH=1,NHHMAX
         DO NTH=1,NTHMAX
            th=dth*(nth-1)
            ph=dph*(nhh-1)
            call wmfem_metrics(xrho(nr),th,ph,gm,gj1)
            call wmfem_metrics(xrho(nr+1),th,ph,gm,gj2)
!            write(6,'(A,1P6E12.4)') 
!     &           'gj:',xrho(nr),xrho(nr+1),th,ph,gj1,gj2
            DSSS=DRHO*0.5d0*(gj1+gj2)
            DVOL(NTH,NHH,NR)=1.D0/DSSS
            DVOLS(NR)=DVOLS(NR)+DSSS*DTH*DPH
         ENDDO
         ENDDO
         DVOLS(NR)=1.D0/DVOLS(NR)
      ENDDO
      DVOLS(NRMAX)=0.D0
      DO NTH=1,NTHMAX
         DO NHH=1,NHHMAX
            DVOL(NTH,NHH,NRMAX)=0.d0
         ENDDO
      ENDDO

      RETURN
      END

