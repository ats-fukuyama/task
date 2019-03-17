MODULE w1prepost

  PRIVATE
  PUBLIC w1pre
  PUBLIC w1post

CONTAINS

! ----- Initial setup for total NMODEL=7-----

  SUBROUTINE w1pre(IERR)

    USE w1comm
    USE w1prof,ONLY: w1_prof
    USE w1sub,ONLY: w1fftl
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR

    IERR=0

    CALL W1SETZ
    CALL W1ANTS
    CALL W1SETX(IERR)
       IF(IERR.NE.0) RETURN
    CALL W1_PROF
    CALL W1PWRI

!     ******* FOURIER TRANSFORM OF ANTENNA CURRENT *******

    CALL W1FFTL(CJ1,NZMAX,0)
    CALL W1FFTL(CJ2,NZMAX,0)
    CALL W1FFTL(CJ3,NZMAX,0)
    CALL W1FFTL(CJ4,NZMAX,0)

    RETURN
  END SUBROUTINE w1pre

! ----- Final setup for total NMODEL=7-----

  SUBROUTINE w1post

    USE w1comm
    USE w1sub,ONLY: w1fftl
    USE w1out,ONLY: w1file,w1prnt
    IMPLICIT NONE
    INTEGER:: NX,IC

!     ******* INVERSE FOURIER TRANSFORM *******

    CALL W1FFTL(CJ1,NZMAX,1)
    CALL W1FFTL(CJ2,NZMAX,1)
    CALL W1FFTL(CJ3,NZMAX,1)
    CALL W1FFTL(CJ4,NZMAX,1)
!
    DO NX=1,NXMAX
       DO IC=1,3
          CALL W1FFTL(CE2DA(1:NZMAX,NX,IC),NZMAX,1)
       END DO
    ENDDO

!     ******* POST PROCESS *******
      
    CALL W1PWRS
    CALL W1PRNT
    CALL W1FILE

    RETURN
  END SUBROUTINE w1post

!     ******* DIVISION IN Z DIRECTION *******

  SUBROUTINE W1SETZ

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZH,NZ
    REAL(rkind):: AKZ0

    NZH=NZMAX/2
    DZ=RZ/DBLE(NZMAX)

    IF(NZMAX.EQ.1) THEN
       ZA(1)=0.D0
       AKZ(1)=RKZ
       IF(ABS(RKZ).LT.1.D-32) RKZ=0.001D0
    ELSE
       DO NZ=1,NZMAX
          ZA(NZ)=DZ*DBLE(NZ-1) - 0.5D0*RZ
       END DO
       AKZ0=2.D0*PI/RZ
       AKZ(1)   =0.001D0
       AKZ(NZH+1)=DBLE(NZH)*AKZ0
       DO NZ=2,NZH
          AKZ(NZ)      =DBLE(NZ-1)*AKZ0
          AKZ(NZMAX-NZ+2)=-AKZ(NZ)
       END DO
    END IF

    RETURN

  END SUBROUTINE W1SETZ


!     ******* DIVISION IN X DIRECTION *******

  SUBROUTINE W1SETX(IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NXVH,NHM,NHM1,NXM,NXQ,NS,NXPMAX1,NSACM,NX,NX1,NX2,IERR,NZH,I
    REAL(rkind):: RW,FWC,WCH,WCL,XACM,DX,FVT,DXACM,WT,X,FACT,X1,DX1,DX2

    NXPMAX=NXMAX*RA/RB
    NXVMAX=(NXMAX-NXPMAX)/2
    NXPMAX=NXMAX-2*NXVMAX
    NXVH=NXVMAX/2
    NXM=NXPMAX*6+10
    NXQ=NXPMAX+1

    IF(MOD(NMODEL,2).EQ.1) THEN
       NXPMAX1=NXPMAX
    ELSE
       NXPMAX1=NXPMAX-1
    ENDIF

    RW=2.D6*PI*RF
    NSACM=0
    DO NS=1,NSMAX
       FWC=AEE*PZ(nS)*BB/(AMP*PA(nS))
       WCH=FWC/(1.D0+RA/RR)
       WCL=FWC/(1.D0-RA/RR)
       IF(ABS(WCH).LT.RW.AND.ABS(WCL).GT.RW) THEN
          XACM=RR*(FWC/RW-1.D0)
          NSACM=NS
       ENDIF
    END DO

    IF(NSACM.EQ.0.OR.DXFACT.LT.1.D-6) THEN
       DX=2.D0*RA/DBLE(NXPMAX1)
       DO NX=1,NXPMAX1+1
          XA(NX)=DBLE(NX-1)*DX-RA
       END DO
    ELSE
       FVT=AEE*1.D3/(AMP*PA(NSACM))
       FWC=AEE*PZ(NSACM)*BB/(AMP*PA(NSACM))
       DXACM=DXWDTH*SQRT(FVT*PTPP(NSACM))/FWC
       WT=2.D0*RA+DXFACT*DXACM &
                 *(ATAN((RA-XACM)/DXACM)-ATAN((-RA-XACM)/DXACM))
       X=-RA
       DO NX=1,NXPMAX1
          XA(NX)=X
          FACT=1.D0+DXFACT/(1.D0+((X-XACM)/DXACM)**2)
          X1=X+0.5D0*WT/(FACT*(NXPMAX1))
          FACT=1.D0+DXFACT/(1.D0+((X1-XACM)/DXACM)**2)
          X=X+WT/(FACT*(NXPMAX1))
       END DO
       WRITE(6,601) NSACM,XACM,DXACM,RA-X
       XA(NXPMAX1+1)=RA
    END IF

    IF(MOD(NMODEL,2).EQ.1) THEN
       DO NX=1,NXPMAX
          XAM(NX)=0.5D0*(XA(NX)+XA(NX+1))
       END DO
    ELSE
       DO NX=1,NXPMAX
          XAM(NX)=XA(NX)
       END DO
    ENDIF

    DX1=(RB-RD)/DBLE(NXVH-1)
    DX2=(RD-RA)/DBLE(NXVH-1)
    DO I=1,NXVMAX
       IF(I.LE.NXVH) THEN
          X=DX1*DBLE(I-1)-RB
       ELSE
          X=DX2*DBLE(I-NXVH-1)-RD
       ENDIF
       NX1=NXPMAX+NXVMAX+I
       XAM(NX1)=X
       NX2=NXPMAX+NXVMAX+1-I
       XAM(NX2)=-X
    END DO
    NXANT1=NXPMAX+NXVMAX+NXVMAX/2+1
    NXANT2=NXPMAX       +NXVMAX/2
    IERR=0
    RETURN

  601 FORMAT('** DX-ACCUM. : NS,X,XWIDTH,XERR = ', &
             I3,1P2E12.4,1PE12.4)
  END SUBROUTINE W1SETX

!     ******* SETING ANTENNA CURRENT DENSITY *******

  SUBROUTINE W1ANTS

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZ, NA
    REAL(rkind):: PHASE,ANTPOS
    COMPLEX(rkind):: CPHASE

    DO NZ = 1 , NZMAX
       CJ1( NZ ) = ( 0.D0 , 0.D0 )
       CJ2( NZ ) = ( 0.D0 , 0.D0 )
       CJ3( NZ ) = ( 0.D0 , 0.D0 )
       CJ4( NZ ) = ( 0.D0 , 0.D0 )
    END DO

    IF(NZMAX.EQ.1) THEN
       PHASE            = APHH(1)*PI/180.D0
       CPHASE           = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
       NZANTYH(1)=1
       CJ1(1)=AJYH(1)/DZ*CPHASE 
       NZANTZH(1)=1
       NZANTLH(1)=1
       CJ3(1)=AJZH(1)

       PHASE            = APHL(1)*PI/180.D0
       CPHASE           = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
       NZANTYL(1)=1
       CJ2(1)=AJYL(1)/DZ*CPHASE 
       NZANTZL(1)=1
       NZANTLL(1)=1
       CJ4(1)=AJZL(1)
    ELSE
       DO NA = 1 , NAMAX
          PHASE            = APHH(NA)*PI/180.D0
          CPHASE           = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
          ANTPOS           = RZ*APYH(NA)/360.D0
          NZANTYH(NA)      = NINT( ANTPOS/DZ ) + 1
          CJ1(NZANTYH(NA)) = CJ1( NZANTYH(NA) ) &
                           + AJYH( NA ) / DZ * CPHASE
          ANTPOS           = RZ*APZH(NA)/360.D0
          NZANTZH(NA)      = NINT( ANTPOS/DZ ) + 1
          ANTPOS           = RZ*(APZH(NA)+ALZH(NA))/360.D0
          NZANTLH(NA)      = NINT( ANTPOS/DZ ) + 1
          DO NZ=NZANTZH(NA),NZANTLH(NA)
             CJ3(NZ) = CJ3(NZ)+AJZH( NA ) / DZ * CPHASE
          END DO

          PHASE            = APHL(NA)*PI/180.D0
          CPHASE           = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
          ANTPOS           = RZ*APYL(NA)/360.D0
          NZANTYL(NA)      = NINT( ANTPOS/DZ ) + 1
          CJ2(NZANTYL(NA)) = CJ2( NZANTYL(NA) ) &
                           + AJYL( NA ) / DZ * CPHASE
          ANTPOS           = RZ*APZL(NA)/360.D0
          NZANTZL(NA)      = NINT( ANTPOS/DZ ) + 1
          ANTPOS           = RZ*(APZL(NA)+ALZL(NA))/360.D0
          NZANTLL(NA)      = NINT( ANTPOS/DZ ) + 1
          DO NZ=NZANTZL(NA),NZANTLL(NA)
             CJ4(NZ) = CJ4(NZ)+AJZL( NA ) / DZ * CPHASE
          END DO
          write(6,'(A,7I5)') 'NZANT:',NA,NZANTYL(NA),NZANTYH(NA), &
                                         NZANTZL(NA),NZANTZH(NA), &
                                         NZANTLL(NA),NZANTLH(NA)
       END DO
    END IF

    RETURN
  END SUBROUTINE W1ANTS

!     ****** INITIALIZE POWER ABSORPTION ******

  SUBROUTINE W1PWRI
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS

    DO NX=1,NXMAX
       FLUXX(NX)=0.D0
    END DO

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABSX(NX,NS)=0.D0
       END DO
    END DO

    DO NX=1,NXPMAX
       AJCDX(NX)=0.D0
    END DO

    RETURN
  END SUBROUTINE W1PWRI

!     ******* POWER ABSORPTION AND ENERGY FLUXX *******

  SUBROUTINE W1PWRS
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NS,NX,NZ,NA
    REAL(rkind):: PHASE
    COMPLEX(rkind):: CPHASE,CJYH,CJZH,CJYL,CJZL
    
    AJCDT=0.D0
    DO NS=1,NSMAX
       PABSXZ(NS)=0.D0
    END DO

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABSXZ(NS)  =PABSXZ(NS)  +PABSX(NX,NS)
       END DO
    END DO

    DO NS=1,NSMAX
       PABSX(     1,NS)=2.D0*PABSX(     1,NS)/(XA(2)-XA(1))
       PABSX(NXPMAX,NS)=2.D0*PABSX(NXPMAX,NS)/(XA(NXPMAX)-XA(NXPMAX-1))
       DO NX=2,NXPMAX-1
          PABSX( NX,NS)=2.D0*PABSX( NX,NS)/(XA(NX+1)-XA(NX-1))
       END DO
    END DO

    PABSTT=0.D0
    DO NS=1,NSMAX
       PABSTT=PABSTT+PABSXZ(NS)
    END DO

    AJCDT=0.D0
    DO NX=1,NXPMAX
       AJCDT =AJCDT +AJCDX(NX)
    END DO
    AJCDX(  1   )=2.D0*AJCDX(  1   )/(XA(2)-XA(1))
    AJCDX(NXPMAX)=2.D0*AJCDX(NXPMAX)/(XA(NXPMAX)-XA(NXPMAX-1))
    DO NX=2,NXPMAX-1
       AJCDX( NX   )=2.D0*AJCDX( NX   )/(XA(NX+1)-XA(NX-1))
    END DO

    PANT1 = 0.D0
    PANT2 = 0.D0
    DO NZ = 1 , NZMAX
       PANT1 = PANT1 + DCONJG(CE2DA(NZ,NXANT1,2))*CJ1(NZ) &
                     + DCONJG(CE2DA(NZ,NXANT1,3))*CJ3(NZ)
       PANT2 = PANT2 + DCONJG(CE2DA(NZ,NXANT2,2))*CJ2(NZ) &
                     + DCONJG(CE2DA(NZ,NXANT2,3))*CJ4(NZ)
    END DO
    PANT1 = - PANT1 * DZ
    PANT2 = - PANT2 * DZ
    PANT  = PANT1 + PANT2

    PWALL1 = -FLUXX(NXPMAX+NXVMAX+1)
    PWALL2 =  FLUXX(NXPMAX+NXVMAX  )
    PWALL  =  PWALL1 + PWALL2

!    PIBW1 = FLUXX(NXPMAX+NXVMAX+NXVMAX) - FLUXX(       1)
    PIBW1 = 0.D0
!    PIBW2 = FLUXX(NXPMAX              ) - FLUXX(NXPMAX+1)
    PIBW2 = 0.D0
    PIBW  = PIBW1 + PIBW2

    PERR  = PANT - ( PABSTT + PWALL + PIBW )

    DO NA = 1 , NAMAX
       RANT1(NA)=0.D0
       XANT1(NA)=0.D0
       CJYH=CJ1(NZANTYH(NA))
       IF( ABS( CJYH ) .GT. 1.E-6 ) THEN
          RANT1(NA) = RANT1(NA)-DBLE( CE2DA(NZANTYH(NA),NXANT1,2)/ CJYH)
          XANT1(NA) = XANT1(NA)+IMAG( CE2DA(NZANTYH(NA),NXANT1,2)/ CJYH)
       END IF
       DO NZ=NZANTZH(NA),NZANTLH(NA)
          CJZH=CJ3(NZ)
          IF( ABS( CJZH ) .GT. 1.E-6 ) THEN
             RANT1(NA) = RANT1(NA)-DBLE( CE2DA(NZ,NXANT1,3)/ CJZH)
             XANT1(NA) = XANT1(NA)+IMAG( CE2DA(NZ,NXANT1,3)/ CJZH)
          END IF
       END DO
       RANT2(NA)=0.D0
       XANT2(NA)=0.D0
       CJYL=CJ2(NZANTYL(NA))
       IF( ABS( CJYL ) .GT. 1.E-6 ) THEN
          RANT2(NA) = RANT2(NA)-DBLE( CE2DA(NZANTYL(NA),NXANT2,2)/ CJYL)
          XANT2(NA) = XANT2(NA)+IMAG( CE2DA(NZANTYL(NA),NXANT2,2)/ CJYL)
       END IF
       DO NZ=NZANTZL(NA),NZANTLL(NA)
          CJZL=CJ4(NZ)
          IF( ABS( CJZL ) .GT. 1.E-6 ) THEN
             RANT2(NA) = RANT2(NA)-DBLE( CE2DA(NZ,NXANT2,3)/ CJZL)
             XANT2(NA) = XANT2(NA)+IMAG( CE2DA(NZ,NXANT2,3)/ CJZL)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE W1PWRS
END MODULE w1prepost
