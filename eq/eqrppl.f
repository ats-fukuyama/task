C
C     ***** CALCULATE RIPPLE WELL REGION INSIDE THE PLASMA *****
C
      SUBROUTINE EQRPPL(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      EXTERNAL EQDERV
      DIMENSION XA(NTVM),YA(2,NTVM)
      DIMENSION URCHI(4,NTVM),UZCHI(4,NTVM)
      DIMENSION rip_rat(NRM),DltRPV(NRM)
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
C     ----- SET DR, DTH -----
C
      IF(NSUMAX.EQ.0) THEN
         NRPMAX=NRMAX
      ELSE
         DR=(RB-RA+REDGE-RAXIS)/(NRMAX-1)
         NRPMAX=NINT((REDGE-RAXIS)/DR)+1
      ENDIF
      DR=(REDGE-RAXIS)/(NRPMAX-1)
      DTH=2.d0*PI/NTHMAX
C
C     ----- SET NUMBER OF DIVISION for integration -----
C
      IF(NTVMAX.GT.NTVM) NTVMAX=NTVM
C
C     ----- CALCULATE PSI and RPS, ZPS 
C           on magnetic surfaces, PSIP(NR) -----
C
      NR=1
      DO NTH=1,NTHMAX+1
         RPS(NTH,NR)=RAXIS
         ZPS(NTH,NR)=ZAXIS
      ENDDO
      PSIP(1)=PSIG(RAXIS,ZAXIS)-PSI0
C
C     +++++ SETUP AXIS DATA +++++
C
      NR = 1
      CALL EQPSID(RAXIS,ZAXIS,DPSIDR,DPSIDZ)
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
      DO N = 1, NA
         GRal(NA*(NR-1)+N) = real(RAXIS)
         GZal(NA*(NR-1)+N) = real(ZAXIS)
         GAlpRP(NR,N)  = real(AlpRP)
      ENDDO
C
      AlpRP0 = AlpRP
      rip_rat(NR) = 0.d0
C
      DO NR=2,NRPMAX
         RINIT=RAXIS+DR*(NR-1)
         ZINIT=ZAXIS
         PSIP(NR)=PSIG(RINIT,ZINIT)-PSI0
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
         idx = 0
C
         DO N=2,NA
            H=XA(N)-XA(N-1)
            R=0.5D0*(YA(1,N-1)+YA(1,N))
            Z=0.5D0*(YA(2,N-1)+YA(2,N))
            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
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
            GRal(NA*(NR-1)+N) = real(R)
            GZal(NA*(NR-1)+N) = real(Z)
            GAlpRP(NR,N) = real(AlpRP)
C
C           --- For TASK/TX ---
C
            if(AlpRP < 1.d0) then
               ARC = ARC + H
            end if
            if(R > RAXIS .and. idx == 0) then
               DltRPV(NR) = DLTRP
               idx = 1
            end if
C
c$$$            R=YA(1,N)
c$$$            Z=YA(2,N)
c$$$            CALL EQPSID(R,Z,DPSIDR,DPSIDZ)
c$$$            BPL=SQRT(DPSIDR**2+DPSIDZ**2)/(2.D0*PI*R)
c$$$            BTL=TTS(NR)/(2.D0*PI*R)
c$$$            B2L=BTL**2+BPL**2
c$$$            B=SQRT(B2L)
c$$$C
c$$$            RMIN=MIN(RMIN,R)
c$$$            RMAX=MAX(RMAX,R)
c$$$            IF(Z.LT.ZMIN) THEN
c$$$               ZMIN=Z
c$$$               NZMINR=N
c$$$            ENDIF
c$$$            IF(Z.GT.ZMAX) THEN
c$$$               ZMAX=Z
c$$$               NZMAXR=N
c$$$            ENDIF
c$$$            BMIN=MIN(BMIN,B)
c$$$            BMAX=MAX(BMAX,B)
         ENDDO
C
         NtrcMAX(NR) = NA
         rip_rat(NR) = ARC / XA(NA)
C
c$$$      DO NR = 1, NRPMAX
c$$$         DO N = 1, NtrcMAX(NR)
c$$$            J = NtrcMAX(NR)*(NR-1)+N
c$$$            if(GRal(J) > 0.d0) write(6,*) GRal(J),GZal(J),GAlpRP(NR,N)
c$$$         ENDDO
c$$$         write(6,*) 
      ENDDO
C
      do nr = 1, nrpmax
         write(6,*) nr,rip_rat(nr),DltRPV(NR)
      enddo
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
