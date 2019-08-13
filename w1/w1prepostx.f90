MODULE w1prepostx

  PRIVATE
  PUBLIC w1prex
  PUBLIC w1postx

CONTAINS

! ----- Initial setup for total NMODEL=7-----

  SUBROUTINE w1prex(IERR)

    USE w1comm
    USE w1prof,ONLY: w1_profx
    USE w1sub,ONLY: w1fftl
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NZ

    IERR=0

    CALL W1SETZ
    CALL W1SETX(IERR)
       IF(IERR.NE.0) RETURN
    CALL W1ANTS
    CALL W1WGS
    CALL W1_PROFX
    CALL W1PWRI

!     ******* FOURIER TRANSFORM OF ANTENNA CURRENT *******

    CALL W1FFTL(CJ1,NZMAX,0)
    CALL W1FFTL(CJ2,NZMAX,0)
    CALL W1FFTL(CJ3,NZMAX,0)
    CALL W1FFTL(CJ4,NZMAX,0)
    CALL W1FFTL(CWG1,NZMAX,0)
    CALL W1FFTL(CWG2,NZMAX,0)
    CALL W1FFTL(CWG3,NZMAX,0)
    CALL W1FFTL(CWG4,NZMAX,0)

    RETURN
  END SUBROUTINE w1prex

! ----- Final setup for total NMODEL=7-----

  SUBROUTINE w1postx

    USE w1comm
    USE w1sub,ONLY: w1fftl
    USE w1out,ONLY: w1file,w1prnt
    USE libgrf
    IMPLICIT NONE
    INTEGER:: NX,IC,NZ

!     ******* INVERSE FOURIER TRANSFORM *******

    CALL W1FFTL(CJ1,NZMAX,1)
    CALL W1FFTL(CJ2,NZMAX,1)
    CALL W1FFTL(CJ3,NZMAX,1)
    CALL W1FFTL(CJ4,NZMAX,1)
    CALL W1FFTL(CWG1,NZMAX,1)
    CALL W1FFTL(CWG2,NZMAX,1)
    CALL W1FFTL(CWG3,NZMAX,1)
    CALL W1FFTL(CWG4,NZMAX,1)
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
  END SUBROUTINE w1postx

!     ******* DIVISION IN Z DIRECTION *******

  SUBROUTINE W1SETZ

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZH,NZ
    REAL(rkind):: DKZ

    DZ=RZ/DBLE(NZMAX)

    IF(NZMAX.EQ.1) THEN
       ZA(1)=0.D0
       AKZ(1)=RKZ
       IF(ABS(RKZ).LT.1.D-32) RKZ=0.001D0
    ELSE
       DO NZ=1,NZMAX
          ZA(NZ)=DBLE(NZ-1)*DZ
       END DO
       DKZ=2.D0*PI/RZ
       AKZ(1)=0.001D0
       AKZ(1+NZMAX/2)=-DBLE(NZMAX/2)*DKZ
       DO NZ=2,NZMAX/2
          AKZ(NZ        )=DBLE(NZ-1        )*DKZ
          AKZ(NZ+NZMAX/2)=DBLE(NZ-1-NZMAX/2)*DKZ
       END DO
    END IF

    RETURN

  END SUBROUTINE W1SETZ


!     ******* DIVISION IN X DIRECTION *******

  SUBROUTINE W1SETX(IERR)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: NX,NSACM,NS
    REAL(rkind):: RW,DX,FWC,WCH,WCL,XACM,FVT,DXACM,WT,X,FACT,X1

    NXPMAX=NXMAX*RA/RB
    NXVMAX=(NXMAX-NXPMAX)/2
    NXPMAX=NXMAX-2*NXVMAX

    RW=2.D6*PI*RF

!    DX=2.D0*RB/(NXMAX-1)
!    DO NX=1,NXMAX
!       XA(NX)=DX*(NX-1)-RB
!    END DO

    NSACM=0
    DO NS=1,NSMAX
       FWC=AEE*PZ(NS)*BB/(AMP*PA(NS))
       WCH=FWC/(1.D0+RB/RR)
       WCL=FWC/(1.D0-RB/RR)
       IF(ABS(WCH).LT.RW.AND.ABS(WCL).GT.RW) THEN
          XACM=RR*(FWC/RW-1.D0)
          NSACM=NS
       ENDIF
    END DO

    IF(NSACM.EQ.0.OR.DXFACT.LT.1.D-6) THEN
       DX=2.D0*RB/DBLE(NXMAX-1)
       DO NX=1,NXMAX
          XA(NX)=DBLE(NX-1)*DX-RB
       END DO
    ELSE
       FVT=AEE*1.D3/(AMP*PA(NSACM))
       FWC=AEE*PZ(NSACM)*BB/(AMP*PA(NSACM))
       DXACM=DXWDTH*SQRT(FVT*PTPP(NSACM))/FWC
       WT=2.D0*RB+DXFACT*DXACM &
                 *(ATAN((RB-XACM)/DXACM)-ATAN((-RB-XACM)/DXACM))
       X=-RB
       DO NX=1,NXMAX-1
          XA(NX)=X
          FACT=1.D0+DXFACT/(1.D0+((X-XACM)/DXACM)**2)
          X1=X+0.5D0*WT/(FACT*(NXMAX-1))
          FACT=1.D0+DXFACT/(1.D0+((X1-XACM)/DXACM)**2)
          X=X+WT/(FACT*(NXMAX-1))
       END DO
       WRITE(6,601) NSACM,XACM,DXACM,RB-X
601    FORMAT('** DX-ACCUM. : IS,X,XWIDTH,XERR = ', &
              I3,1P2E12.4,1PE12.4) 
       XA(NXMAX)=RB
    END IF

!    DO NX=1,NXMAX-1
!       XAM(NX)=0.5D0*(XA(NX)+XA(NX+1))
!    END DO
    DO NX=1,NXMAX
       XAM(NX)=XA(NX)
    END DO

    NXANT1=1
    NXANT2=1
    DO NX=1,NXMAX-1
       IF(XA(NX).LT.-RD.AND.-RD.GE.XA(NX+1)) NXANT1=NX
       IF(XA(NX).LE. RD.AND. RD.GT.XA(NX+1)) NXANT2=NX
    END DO
    WRITE(6,*) 'NXANT1,NXANT2=',NXANT1,NXANT2

    IERR=0
    RETURN
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

!     ******* SETING WAVE GUIDE ELECTRIC FIELD *******

  SUBROUTINE W1WGS

    USE w1comm
    IMPLICIT NONE
    INTEGER:: NZ
    REAL(rkind):: RW,RKWG,ZOFFSET,FACT,ZSIZE
    COMPLEX(rkind):: CAMP

    DO NZ = 1 , NZMAX
       CWG1( NZ ) = ( 0.D0 , 0.D0 )
       CWG2( NZ ) = ( 0.D0 , 0.D0 )
       CWG3( NZ ) = ( 0.D0 , 0.D0 )
       CWG4( NZ ) = ( 0.D0 , 0.D0 )
    END DO

    IF(NZMAX.EQ.1) THEN
       SELECT CASE(MDLWG)
       CASE(1)
          CWG1(1)=WGAMP
       CASE(2)
          CWG2(1)=WGAMP
       CASE(3)
          CWG3(1)=WGAMP
       CASE(4)
          CWG4(1)=WGAMP
       END SELECT
    ELSE
       RW=2.D6*PI*RF
       RKWG=WGNZ*RW/VC
       WRITE(6,'(A,1P3E12.4)') 'RF,RW,VC=',RF,RW,VC
       IF(ABS(RKWG).GT.1.E-8) THEN
          WRITE(6,'(A,1P2E12.4)') 'RKWG,LWG=',RKWG,2.D0*PI/RKWG
       END IF
       DO NZ=1,NZMAX
          IF(ZA(NZ).GE.WGZ1.AND.ZA(NZ).LE.WGZ2) THEN
             ZOFFSET=ZA(NZ)-0.5D0*(WGZ1+WGZ2)
             CAMP=WGAMP*EXP(CI*RKWG*ZOFFSET)
             SELECT CASE(MDLWGS)
             CASE(0)
                FACT=1.D0
             CASE(1)
                ZSIZE=0.5D0*(WGZ2-WGZ1)
                FACT=1.D0-(ZOFFSET/ZSIZE)**2
             CASE(2)
                ZSIZE=0.5D0*(WGZ2-WGZ1)
                FACT=EXP(-3.D0*(ZOFFSET/ZSIZE)**2)
             END SELECT
             SELECT CASE(MDLWG)
             CASE(1)
                CWG1(NZ)=CAMP*FACT
             CASE(2)
                CWG2(NZ)=CAMP*FACT
             CASE(3)
                CWG3(NZ)=CAMP*FACT
             CASE(4)
                CWG4(NZ)=CAMP*FACT
             END SELECT
          END IF
       END DO
    END IF

    RETURN
  END SUBROUTINE W1WGS

!     ****** INITIALIZE POWER ABSORPTION ******

  SUBROUTINE W1PWRI
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS

    DO NX=1,NXMAX
       FLUXX(NX)=0.D0
    END DO

    DO NS=1,NSMAX
       DO NX=1,NXMAX-1
          PABSX(NX,NS)=0.D0
       END DO
    END DO

    DO NX=1,NXMAX-1
       AJCDX(NX)=0.D0
    END DO

    RETURN
  END SUBROUTINE W1PWRI

!     ******* POWER ABSORPTION AND ENERGY FLUXX *******

  SUBROUTINE W1PWRS
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NS,NX,NXANT,NZ,NA
    REAL(rkind):: PHASE,DX
    COMPLEX(rkind):: CPHASE,CJYH,CJZH,CJYL,CJZL
    
    AJCDT=0.D0
    DO NS=1,NSMAX
       PABSXZ(NS)=0.D0
    END DO

    DO NS=1,NSMAX
       DO NX=1,NXMAX-1
          PABSXZ(NS)  =PABSXZ(NS)  +PABSX(NX,NS)
          IF(NS.EQ.1) WRITE(21,'(A,I5,1P2E12.4)') &
               'NX,PABSX,PABSXZ=',NX,PABSX(NX,NS),PABSXZ(NS)
       END DO
    END DO

    DO NS=1,NSMAX
       PABSX(      1,NS)=2.D0*PABSX(      1,NS)/(XA(2)-XA(1))
       PABSX(NXMAX-1,NS)=2.D0*PABSX(NXMAX-1,NS)/(XA(NXMAX)-XA(NXMAX-1))
       DO NX=2,NXMAX-2
          PABSX(NX,NS)=2.D0*PABSX(   NX,NS)/(XA(NX+1)-XA(NX-1))
       END DO
    END DO

    PABSTT=0.D0
    DO NS=1,NSMAX
       PABSTT=PABSTT+PABSXZ(NS)
    END DO

    AJCDT=0.D0
    DO NX=1,NXMAX-1
       AJCDT =AJCDT +AJCDX(NX)
    END DO
    AJCDX(      1)=2.D0*AJCDX(      1)/(XA(2)-XA(1))
    AJCDX(NXMAX-1)=2.D0*AJCDX(NXMAX-1)/(XA(NXMAX)-XA(NXMAX-1))
    DO NX=2,NXMAX-2
       AJCDX(NX)=2.D0*AJCDX(   NX)/(XA(NX+1)-XA(NX-1))
    END DO

    DX=2.D0*RB/(NXMAX-1)
    PANT1 = 0.D0
    PANT2 = 0.D0
    DO NZ = 1 , NZMAX
       PANT1 = PANT1 + DZ*REAL(CPANTK1(NZ))
       PANT2 = PANT2 + DZ*REAL(CPANTK2(NZ))
    END DO
    PANT  = PANT1 + PANT2

    PWALL1 = -FLUXX(1)
    PWALL2 =  FLUXX(NXMAX)
    PWALL  =  PWALL1 + PWALL2

    PIBW1 = 0.D0
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
END MODULE w1prepostx
