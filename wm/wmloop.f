C***********************************************************************
C
C     Multi-mode calculation
C
C***********************************************************************
c
      SUBROUTINE WM_LOOP(IERR)
c
      include '../wm/wmcomm.inc'
c
      IF(NPHMAX.EQ.1) THEN
         CALL WMEXEC(IERR)
         RETURN
      ENDIF

      NPH0_SV  = NPH0
      NPHMAX_SV = NPHMAX
      NHHMAX_SV = NHHMAX

      IF(NHHMAX.GT.1) THEN
         NHHMAX=NPHMAX/NHC
         NPHMAX=NHC
      END IF
C
C     Axisymmetric case : NHHMAX=1
C
      DO NR=1,NRMAX+1
         DO NPH=1,NPHMAX_SV
            DO NTH=1,NTHMAX
               CEFLD3D(1,NTH,NPH,NR)=(0.D0,0.D0)
               CEFLD3D(2,NTH,NPH,NR)=(0.D0,0.D0)
               CEFLD3D(3,NTH,NPH,NR)=(0.D0,0.D0)
            END DO
         END DO
      END DO
      DO NS=1,NSMAX
         DO NR=1,NRMAX+1
            DO NPH=1,NPHMAX_SV
               DO NTH=1,NTHMAX
                  CPABS3D(NTH,NPH,NR,NS)=0.D0
               END DO
            END DO
         END DO
      END DO

      DO NPH = 1,NPHMAX
         NPH0 = NPH-NPHMAX/2-1
         WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0
C
         CALL WMEXEC(IERR)
         CALL mtx_barrier
         IF(IERR.NE.0) EXIT

         DO NR=1,NRMAX+1
            DO NTH=1,NTHMAX
               CEFLD3D(1,NTH,NPH,NR)=CEFLD(1,NTH,1,NR)
               CEFLD3D(2,NTH,NPH,NR)=CEFLD(2,NTH,1,NR)
               CEFLD3D(3,NTH,NPH,NR)=CEFLD(3,NTH,1,NR)
            END DO
         END DO

         DO NS=1,NSMAX
            DO NR=1,NRMAX+1
               DO NTH=1,NTHMAX
                  CPABS3D(NTH,NPH,NR,NS)=CPABS(NTH,1,NR,NS)
               END DO
            END DO
         END DO

         DO NS=1,NSMAX
            PABST3D(NPH,NS)=PABST(NS)
         END DO
         PABSTT3D(NPH)=PABSTT

C         write(6,'(A,1P2E12.4)') 
C     &        'CEFLD3D(1,5,NPH,5)=',CEFLD3D(1,5,NPH,5)
C
      ENDDO
C
      NPH0  = NPH0_SV
C
      return
      end
