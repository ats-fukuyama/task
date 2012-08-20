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
C
C     Axisymmetric case : NHHMAX=1
C
      DO NR=1,NRMAX
         DO NPH=1,NPHMAX
            DO NTH=1,NTHMAX
               CEFLD3D(1,NTH,NPH,NR)=(0.D0,0.D0)
               CEFLD3D(2,NTH,NPH,NR)=(0.D0,0.D0)
               CEFLD3D(3,NTH,NPH,NR)=(0.D0,0.D0)
            END DO
         END DO
      END DO
      DO NS=1,NSMAX
         DO NR=1,NRMAX
            DO NPH=1,NPHMAX
               DO NTH=1,NTHMAX
                  PABS3D(NTH,NPH,NR,NS)=0.D0
               END DO
            END DO
         END DO
      END DO

      DO NPH = 1,NPHMAX
         IF(NPH.LE.NPHMAX/2+1) THEN
            NPH0 = NPH-1
         ELSE
            NPH0 = NPH-NPHMAX-1
         ENDIF
         WRITE(6,'(A,I4)') "== Toroidal mode number : ",NPH0
C
         CALL WMEXEC(IERR)
         CALL MPSYNC
         IF(IERR.NE.0) EXIT

         DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               CEFLD3D(1,NTH,NPH,NR)=CEFLDK(1,NTH,1,NR)
               CEFLD3D(2,NTH,NPH,NR)=CEFLDK(2,NTH,1,NR)
               CEFLD3D(3,NTH,NPH,NR)=CEFLDK(3,NTH,1,NR)
            END DO
         END DO

         DO NS=1,NSMAX
            DO NR=1,NRMAX
               DO NTH=1,NTHMAX
                  PABS3D(NTH,NPH,NR,NS)=PABS3D(NTH,1,NR,NS)
               END DO
            END DO
         END DO
C
      ENDDO
C
      NPH0  = NPH0_SV
C
      return
      end
