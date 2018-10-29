C
C     ***** CALCULATE RIPPLE WELL REGION INSIDE THE PLASMA *****
C
      SUBROUTINE EQRPPL(IERR)
C
      USE libfio
      INCLUDE '../eq/eqcomq.inc'
C
      EXTERNAL EQDERV
      DIMENSION XA(NTVM),YA(2,NTVM)
      DIMENSION rip_rat(NRM),DltRPV(NRM),DltRP_rim(NRM),theta_rim(NRM)
      DIMENSION DltRP_mid(NRM)
      CHARACTER*80 FNAME
C
      IERR=0
C
C     ----- Number of TF coils -----
C
      NCOIL = 18 ! for JT-60
C
C     ----- Model of ripple well parameter -----
C        0 : General expression (K. Tani et al., NF 33 (1993) 903)
C        1 : Simpler expression based on the quasi-toroidal coordinates
C            (e.g. R.J. Goldston et al., PRL 47 (1981) 647)
C            
C
      MDLAlp = 0 
C
C     ----- Load ripple data -----
C
      call eq_set_rppl(IERR)
C
C     ----- SET DR -----
C
      IF(NSUMAX.EQ.0) THEN
         NRPMAX=NRMAX
      ELSE
         DR=(RB-RA+REDGE-RAXIS)/(NRMAX-1)
         NRPMAX=NINT((REDGE-RAXIS)/DR)+1
      ENDIF
      DR=(REDGE-RAXIS)/(NRPMAX-1)
C
C     ----- SET NUMBER OF DIVISION for integration -----
C
      IF(NTVMAX.GT.NTVM) NTVMAX=NTVM
C
C     +++++ SETUP AXIS DATA +++++
C
      NR = 1
      NA=NTVMAX
      CALL PSIGD(RAXIS,ZAXIS,DPSIDR,DPSIDZ)
      BTL=TTS(NR)/(2.D0*PI*RAXIS)
      CALL SPL2DF(RAXIS,ZAXIS,DLTRP,Rrp,Zrp,URpplRZ,
     &                  NRrpM,NRrpM,NZrpM,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQRPPL: SPL2DF ERROR : IERR=',IERR
      IF(MDLAlp == 0) THEN
         BRL = -DPSIDZ / (2.D0 * PI * RAXIS)
         AlpRP = abs(BRL / BTL) / (DLTRP * NCOIL)
      ELSE
         AlpRP = 0.d0
      END IF
      NtrcMAX(NR) = 1
      DO N = 1, NtrcMAX(NR)
         GRal  (NR,N) = real(RAXIS)
         GZal  (NR,N) = real(ZAXIS)
         GAlpRP(NR,N) = real(AlpRP)
      ENDDO
C
      DltRP_rim(NR) = 0.d0
      theta_rim(NR) = 0.d0
      rip_rat(NR)   = 0.d0
      DltRPV(NR)    = DLTRP
      DltRP_mid(NR) = DLTRP
C
      DO NR=2,NRPMAX
         RINIT=RAXIS+DR*(NR-1)
         ZINIT=ZAXIS
C
C         WRITE(6,'(A,I5,1P3E12.4)') 'NR:',NR,
C     &        PSIP(NR),RINIT,ZINIT
C
         CALL EQMAGS(RINIT,ZINIT,NTVMAX,XA,YA,NA,IERR)
C
c$$$         RMIN=RAXIS
c$$$         RMAX=RAXIS
c$$$         ZMIN=ZAXIS
c$$$         ZMAX=ZAXIS
c$$$         BMIN=ABS(2.D0*BB)
c$$$         BMAX=0.D0
C
         ARC = 0.d0
         idx1 = 0
         idx2 = 0
C
C     --- Ripple amplitude at the outer midplane ---
C
         R = YA(1,1)
         Z = YA(2,1)
         CALL SPL2DF(R,Z,DltRP_mid(NR),Rrp,Zrp,URpplRZ,
     &               NRrpM,NRrpM,NZrpM,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX EQRPPL: SPL2DF ERROR : IERR=',IERR
C
         DO N=2,NA
            H=XA(N)-XA(N-1)
            R=0.5D0*(YA(1,N-1)+YA(1,N))
            Z=0.5D0*(YA(2,N-1)+YA(2,N))
            CALL PSIGD(R,Z,DPSIDR,DPSIDZ)
            BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
            BTL=TTS(NR)/(2.D0*PI*R)
            B2L=BTL**2+BPL**2
C
C           --- Evaluate ripple value at (R,Z) ---
C
            CALL SPL2DF(R,Z,DLTRP,Rrp,Zrp,URpplRZ,
     &                  NRrpM,NRrpM,NZrpM,IERR)
            IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX EQRPPL: SPL2DF ERROR : IERR=',IERR
C
C           --- Evaluate ripple well condition ---
C
            IF(MDLAlp == 0) THEN
               BRL = -DPSIDZ / (2.D0 * PI * R)
               AlpRP = abs(BRL / SQRT(B2L)) / (DLTRP * NCOIL)
            ELSE
               RL    = SQRT((R - RAXIS)**2 + (Z - ZAXIS)**2)
               sinth = (Z - ZAXIS) / RL
               AlpRP = (RL / R) * abs(sinth) / (NCOIL * QPS(NR) * DLTRP)
            END IF
!            write(6,'(2I4,2F15.7)') NR,N,H,XA(N)
C
            GRal  (NR,N) = real(R)
            GZal  (NR,N) = real(Z)
            GAlpRP(NR,N) = real(AlpRP)
C
C           --- For TASK/TX ---
C
            if(AlpRP < 1.d0) then ! Ripple well
               ARC = ARC + H
               if(idx1 == 0) then ! Upper side
                  DltRP_rim(NR) = DLTRP
                  RL    = SQRT((R - RAXIS)**2 + (Z - ZAXIS)**2)
                  sinth = (Z - ZAXIS) / RL
                  theta_rim(NR) = asin(sinth)
               else if(R > RAXIS) then ! Lower side
                  RL    = SQRT((R - RAXIS)**2 + (Z - ZAXIS)**2)
                  sinth = (Z - ZAXIS) / RL
                  if(abs(asin(sinth)) > abs(theta_rim(NR))) then
                     theta_rim(NR) = abs(asin(sinth))
                     DltRP_rim(NR) = DLTRP
                  end if
               end if
            else ! No ripple well
               idx1 = 1
            end if
            if(R < RAXIS .and. idx2 == 0) then ! Upper half side
               DltRPV(NR) = DLTRP 
               idx2 = 1
            end if
            if(R > RAXIS .and. idx2 == 1) then ! Lower half side
               ! Average between upper and lower half sides
               DltRPV(NR) = 0.5d0 * (DltRPV(NR) + DLTRP)
               idx2 = 2
            end if
         ENDDO
C
         NtrcMAX(NR) = NA
         rip_rat(NR) = ARC / XA(NA)
C
      ENDDO
C
C     ----- File output for TASK/TX -----
C
      ntxout=22
      FNAME='tx_ripple.dat'
      CALL FWOPEN(ntxout,FNAME,1,0,'EQ RIPPLE TO TX',IERR)
      write(ntxout,'(I4)') nrpmax
      write(ntxout,'(2X,A,6X,A,9X,A,6X,A,7X,A,9X,A,8X,A)') 'NR','RHON',
     &     'DltRP_rim','theta_rim','rip_rat','DltRP','DltRP_mid'
      do nr = 1, nrpmax
         PSIPNL = PSIP(NR) / PSIPA
         RHON = FNRHON(PSIPNL)
         write(ntxout,'(I4,1P6E15.7)') nr,RHON,DltRP_rim(nr),
     &        theta_rim(nr),rip_rat(nr),DltRPV(NR),DltRP_mid(NR)
      enddo
      close(ntxout)
C
      RETURN
      END
C
C     ***** Spline coefficient matrix for ripple amplitude *****
C
      subroutine eq_set_rppl(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      DIMENSION RpplRG(NRrpM,NZrpM),RpplZG(NRrpM,NZrpM),
     &          RpplRZG(NRrpM,NZrpM)
C
      CALL SPL2D(Rrp,Zrp,RpplRZ,RpplRG,RpplZG,RpplRZG,URpplRZ,
     &           NRrpM,NRrpM,NZrpM,0,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for RpplRZ: IERR=',IERR
C
      return
      end
