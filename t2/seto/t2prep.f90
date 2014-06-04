!--------------------------------------------------------------------- 
!
!    MODULE FOR GENERATING GLOBAL PARAMETERS AND INITIAL PROFILES
!
!
!                   LAST UPDATE 2014-05-28 H.Seto
!
!---------------------------------------------------------------------
MODULE T2PREP
  
  USE T2CNST,ONLY: ikind,rkind
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC T2_PREP
  
CONTAINS

  !-------------------------------------------------------------------
  !
  !
  !
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE T2_PREP
    
    USE T2COMM,ONLY: time_t2,time_init
    USE T2DIV, ONLY: T2_DIV
    USE T2PROF,ONLY: T2_PROF
    
    REAL(4):: e0time_0,e0time_1
    
    !d0iar   = d0rmnr/d0rmjr
    
    ! set up normalization factors for T2CALV
    CALL T2PREP_SETUP_NORMALIZATION_FACTOR 
    
    ! set up global parameters for T2COMM
    CALL T2PREP_SETUP_GLOBAL_PARAMETER
    
    ! set up initial profiles
    CALL T2PREP_SETUP_INITIAL_PROFILES
    
    CALL T2_DIV
    
    CALL CPU_TIME(e0time_0)
    CALL T2_PROF
    CALL CPU_TIME(e0time_1)
    WRITE(6,'(A,F10.3,A)') '-- Initial profile calculated: cpu=', &
                         e0time_1-e0time_0,' [s]'

 !   IF(idfile.ge.1) CALL T2_WRIT_MAIN(d1guv,0,c10rname)
 !   IF(idfile.ge.2) CALL T2_WRIT_GPT(20,0,d1guv)
 !   IF(idfile.ge.3) CALL T2_WRIT_GP1(22,0,d1guv)

    time_t2=time_init

    RETURN

  END SUBROUTINE T2_PREP
  
  !--------------------------------------------------------------
  ! 
  !       SET NORMALIZATION FACTORS FOR EQAUTIONS 
  !       IN ORDER TO IMPROVE CONVERGENCE CHARACTERISTIC
  !
  !                           LAST UPDATE 2014-05-28 H.Seto          
  !
  !--------------------------------------------------------------
  SUBROUTINE T2PREP_SETUP_NORMALIZATION_FACTOR
    USE T2CNST,ONLY:Aee
    USE T2COMM,ONLY:UsePotentialDescription,&
         &          UseNormalization,&
         !
         &          BpNF, BtNF, EtNF, EpNF, ErNF,&
         &          NnNF, FrNF, FbNF, FtNF, FpNF,&
         &          PpNF, QrNF, QbNF, QtNF, QpNF
    
    IF(.NOT.UsePotentialDescription)THEN
       IF(UseNormalization)THEN
          BpNF = 1.D0                  ! [        T    ]
          BtNF = 1.D0                  ! [        T    ]
          EtNF = 1.D0                  ! [        V/m  ]
          EpNF = 1.D0                  ! [        V/m  ]
          ErNF = 1.D3                  ! [       kV/m  ] 
          NnNF = 1.D20                 ! [10^20    /m3 ]
          FrNF = NnNF                  ! [10^20    /m2s]
          FbNF = NnNF*1.D3             ! [10^23    /m2s]
          FtNF = NnNF*1.D3             ! [10^23    /m2s]
          FpNF = NnNF*1.D3             ! [10^23    /m2s]
          PpNF = NnNF*1.D3*Aee         ! [10^20 keV/m3 ]
          QrNF = FrNF*1.D3*Aee         ! [10^20 keV/m2s]
          QbNF = FbNF*1.D3*Aee         ! [10^23 keV/m2s]
          QtNF = FtNF*1.D3*Aee         ! [10^23 keV/m2s]
          QpNF = FpNF*1.D3*Aee         ! [10^23 keV/m2s]
       ELSE
          BpNF = 1.D0         ! 
          BtNF = 1.D0         ! 
          EtNF = 1.D0         !
          EpNF = 1.D0         !
          ErNF = 1.D0         !
          NnNF = 1.D0         !
          FrNF = 1.D0         !
          FbNF = 1.D0         !
          FtNF = 1.D0         !
          FpNF = 1.D0         !
          PpNF = 1.D0         !
          QrNF = 1.D0         !
          QbNF = 1.D0         !
          QtNF = 1.D0         !
          QpNF = 1.D0         !
       ENDIF
    ELSE
       WRITE(6,*)"T2PREP_SETUP_NORMALIZATION_FACTOR              "
       WRITE(6,*)"Potential Description Ver. is underconstruction"
       STOP
    ENDIF
    
    RETURN
  
  END SUBROUTINE T2PREP_SETUP_NORMALIZATION_FACTOR

  !--------------------------------------------------------------
  ! 
  !       SET GLOBAL PARAMETERS FOR T2COMM_ALLOCATE
  !
  !                           LAST UPDATE 2014-05-28 H.Seto          
  !
  !--------------------------------------------------------------
  SUBROUTINE T2PREP_SETUP_GLOBAL_PARAMETER
    
    USE T2COMM,ONLY:&
         NLMAX, NVMAX, NKMAX, NQMAX, NNMAX, NMMAX, NXMAX, &
         NBMAX, NRMAX, NEMAX, NHMAX, NAMAX, NNRMX, NERMX, &
         NECMX, NDMAX, NSMAX, NPMAX, NPMIN,        &
         !
         UsePotentialDescription,&
         !
         i1mlvl,i1pdn1,i1pdn2,i1rdn1,i1rdn2,i1mmax,i1bmax,i1emax,&
         !
         T2NGRA_ALLOCATE, T2COMM_ALLOCATE

    INTEGER(ikind)::&
         i0mlva,i0mlvb,i0mlvc,i_l,i0lidj,i0mlvl
    
    CALL T2NGRA_ALLOCATE
    
    DO i_l = 1, NLMAX
       i0mlvl=i1mlvl(i_l)-1
       i1pdn2(i_l) = NPMIN*(2**i0mlvl)   
    ENDDO
    
    ! CALCULATE NUMBER OF NODES IN EACH MESH LEVEL    
    NMMAX = 0
    NBMAX = 1
    NXMAX = 0
    NRMAX = 1
    
    DO i_l = 1, NLMAX
       i1rdn1(i_l) = i1rdn2(i_l) + 1
       i1pdn1(i_l) = i1pdn2(i_l) + 1
    ENDDO
    
    SELECT CASE (NNMAX)
       !  4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !  9: QADRADIC RECTANGULAR ELEMENT (  9 POINTS)
       ! 16: CUBIC    RECTANGULAR ELEMENT ( 16 POINTS)
    CASE( 4)
       DO i_l = 1, NLMAX
          i0mlva = i1mlvl(i_l-1)
          i0mlvb = i1mlvl(i_l  )
          i0mlvc = i1mlvl(i_l+1)
          i1mmax(i_l) = i1rdn1(i_l)*i1pdn1(i_l)
          IF(i0mlva.EQ.0)THEN
             i1bmax(i_l) = i1rdn2(i_l)*i1pdn2(i_l)+1
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             i1bmax(i_l) = i1rdn2(i_l)*i1pdn2(i_l)+i1pdn2(i_l-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             i1bmax(i_l) = i1rdn2(i_l)*i1pdn2(i_l)+i1pdn2(i_l)
          ENDIF
          NMMAX = NMMAX + i1rdn1(i_l)*i1pdn1(i_l)
          NBMAX = NBMAX + i1rdn2(i_l)*i1pdn2(i_l)
       ENDDO
       
       NXMAX = NBMAX
       
       DO i_l = 2, NLMAX
          i0mlva=i1mlvl(i_l-1)
          i0mlvb=i1mlvl(i_l  )
          IF(i0mlva.NE.i0mlvb) NXMAX = NXMAX + i1pdn2(i_l-1)
       ENDDO
       
       DO i_l = 1, NLMAX
          NRMAX = NRMAX + i1rdn2(i_l)
       ENDDO
       
       NERMX = NXMAX + 1
       NECMX = i1pdn2(1) + 4*(NBMAX - i1pdn1(NLMAX))&
            + 2*(NXMAX - NBMAX +i1pdn2(NLMAX))
    CASE ( 9)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE (16)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE DEFAULT
       WRITE(6,*)'T2_GRPH: IMPROPER IMPUT NNMAX0=',NNMAX
       STOP
    END SELECT
    
    !     CALCULATE NUMBER OF ELEMENTS IN EACH MESH LEVEL
    !                  NEMAX, I1EMAX
    i1emax(0:NLMAX)= 0
    
    DO i_l = 1, NLMAX
       DO j_l = 1, i_l
          i1emax(i_l) = i1emax(i_l) + i1rdn2(j_l)*i1pdn2(j_l)
       ENDDO
    ENDDO
    
    NEMAX = i1emax(NLMAX)
    
    ! CALCULATE NUMBER OF HANGED NODE
    !                  NHMAX
    NHMAX = 0
    
    DO i_l = 2, NLMAX
       i0mlva=i1mlvl(i_l-1)
       i0mlvb=i1mlvl(i_l  )
       IF(i0mlva.NE.i0mlvb) NHMAX = NHMAX + i1pdn2(i_l-1)
    ENDDO
    
    
    !     CALCULATE ARRAY SIZE OF MATRIX INDICES FOR NODE
    !     NNRMX,I0NCMX
    
    NNRMX = 1 + NBMAX
    NAMAX = 0
    SELECT CASE (NNMAX)
       !   4: LINEAR   RECTANGULAR ELEMENT (  4 POINTS)
       !   9: QADRADIC RECTANGULAR ELEMENT (  9 POINTS)
       !  16: CUBIC    RECTANGULAR ELEMENT ( 16 POINTS)
    CASE( 4)
       
       DO i_l = 1, NLMAX
          i0mlva = i1mlvl(i_l-1)
          i0mlvb = i1mlvl(i_l  )
          i0mlvc = i1mlvl(i_l+1)
          !C POINTS ON LEFT-SIDE-EDGES 
          IF(    i0mlva.EQ.0)THEN
             !C FIRST AND SECOND EDGE POINTS ON LEFT-SIDE IN Lv-1 MESH
             NAMAX = NAMAX + 1+i1pdn2(i_l  ) +  7*i1pdn2(i_l  )
          ELSEIF(i0mlva.NE.i0mlvb)THEN
             !C SECOND EDGE POINTS ON LEFT-SIDES IN Lv-2~LMAX MESH
             NAMAX = NAMAX + 8*i1pdn2(i_l-1) + 11*i1pdn2(i_l-1)
          ELSEIF(i0mlva.EQ.i0mlvb)THEN
             NAMAX = NAMAX + 9*i1pdn2(i_l  )
          ENDIF
          !C CONTSRAINT FREE POINTS ON MIDDLE AREA
          IF(i1rdn2(i_l).GE.2)THEN
             NAMAX = NAMAX + 9*i1pdn2(i_l)*(i1rdn2(i_l)-2)
          ENDIF
          !C FIRST EDGE POINTS ON RIGHT-SIDE
          IF(    i0mlvc.EQ.0     )THEN  
             NAMAX = NAMAX +  6*i1pdn2(i_l)
          ELSEIF(i0mlvb.NE.i0mlvc)THEN
             NAMAX = NAMAX + 11*i1pdn2(i_l)
          ELSEIF(i0mlvb.EQ.i0mlvc)THEN
             NAMAX = NAMAX +  9*i1pdn2(i_l)
          ELSE
             WRITE(6,*)'ERROR IN SET_MATRIX_INDICES_FOR_N-CRS'
             STOP
          ENDIF
       ENDDO
    CASE ( 9)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE (16)
       WRITE(6,*)'UNDER CONSTRUCTION'
       STOP
    CASE DEFAULT
       WRITE(6,*)'T2PREP: IMPROPER IMPUT NNMAX=',NNMAX 
       STOP
    END SELECT
    
    IF(.NOT.UsePotentialDescription)THEN
       NVMAX  = 10*NSMAX + 5
       NKMAX  =  2*NSMAX + 2 
    ELSE
       WRITE(6,*)"Potential Description Ver. is underconstruction"
       STOP
    END IF
    
    WRITE(6,*)'NDMAX=',NDMAX,'NNMAX=',NNMAX,'NSMAX=',NSMAX
    WRITE(6,*)'NLMAX=',NLMAX,'NPMAX=',NPMAX,'NQMAX=',NQMAX
    WRITE(6,*)'NMMAX=',NMMAX,'NRMAX=',NRMAX,'NEMAX=',NEMAX
    WRITE(6,*)'NAMAX=',NAMAX,'NXMAX=',NXMAX,'NBMAX=',NBMAX
    WRITE(6,*)'NNRMX=',NNRMX,'NERMX=',NERMX,'NECMX=',NECMX
    WRITE(6,*)'NKMAX=',NKMAX,'NHMAX=',NHMAX
    
    CALL T2COMM_ALLOCATE
    
    RETURN
  END SUBROUTINE T2PREP_SETUP_GLOBAL_PARAMETER
  
  SUBROUTINE T2PREP_SETUP_INITIAL_PROFILES
    
    USE T2COMM,ONLY:CoordinateSwitch
    SELECT CASE(CoordinateSwitch)
       ! SELECT COORDINATE TYPE
       ! 1: TOROIDAL COORDINATE W/O EQUILIBRIUM: (\sigma,\chi,\zeta)
       !       \sigma: area based radial rabel r/a = \sqrt{\sigma}
       !                                          [0,(b/a)^2]
       !       \chi  : geometrical poloidal angle [0,2\pi]
       !       \chi  : geometrical toroidal angle [0,2\pi]
    CASE(1)
       CALL T2PROF_TOROIDAL_COORDINATE_AREA
    CASE DEFAULT
       WRITE(6,*)'IMPROPER INPUT >> I0MFCS'
       STOP
    END SELECT
    RETURN
  END SUBROUTINE T2PREP_SETUP_INITIAL_PROFILES

  SUBROUTINE T2PROF_TOROIDAL_COORDINATE_AREA
    
    USE T2COMM,ONLY:&
         BpNF,BtNF,EtNF,EpNF,EtNF,&
         NnNF,FrNF,FbNF,FtNF,FpNF,&
         PpNF,QrNF,QbNF,QtNF,QpNF,&
         NSMAX,NXMAX,NVMAX,NMMAX, &
         i1pdn1,i2crt,&
         d2mfc1,d2rzm,d2rzx,Metrics,XvecIn
    
    USE T2GOUT, ONLY: T2_GOUT
    USE T2CONV, ONLY: T2_CONV
    INTEGER(i0ikind)::&
         i_v,i_s,i_m,i_x,i_x1d,i_x2d
    
    REAL(   rkind)::r_mc,p_mc
    REAL(   rkind),DIMENSION(    1:NSMAX)::nn,pp
    REAL(   rkind),DIMENSION(1:8,1:NSMAX)::ff
    
100 FORMAT( 6E15.8)
110 FORMAT(10E15.8)
120 FORMAT( 5E15.8)
    
    XvecIn(1:NVMAX,1:NXMAX) = 0.D0
    
    DO i_m = 1, NMMAX
       
       !  INITIIALIZATION     
       r_mc = d2mfc1(1,i_m)
       p_mc = d2mfc1(2,i_m)
       
       !  R,Z in cylindrical coordinate (R,\vp,Z)
       d2rzm(1,i_m)  = FUNC_r_rz(r_mc,p_mc)
       d2rzm(2,i_m)  = FUNC_z_rz(r_mc,p_mc)
       
       !  CONTRAVARIANT GEOMETRIC TENSOR
       
       Metrics(1:9,i_m)= FUNC_metrics(r_mc,p_mc)
       i_x2d = i2crt(2,i_m)
       i_x1d = i2crt(3,i_m)
       
       !  INITIAL PROFILES
       
       !
       !  variables as field (from i_v= 1 to i_v = 5)
       !       
       XvecIn(1,i_x1d) = FUNC_dpsidr(r_mc,p_mc)/BpNF ! dPsidr
       XvecIn(2,i_x1d) = FUNC_btCo(  r_mc,p_mc)/BtNF ! I
       XvecIn(3,i_x1d) = FUNC_etCo(  r_mc,p_mc)/EtNF ! E_{\zeta}
       XvecIn(4,i_x2d) = FUNC_epCo(  r_mc,p_mc)/EpNF !\bar{E}_{\chi}
       XvecIn(5,i_x2d) = FUNC_erCo(  r_mc,p_mc)/ErNF ! E_{\sigma}
       
       !
       !  variables as fluid (from i_v = 6 to i_v = 10*NSMAX+5)
       !
       nn = FUNC_nn(r_mc,p_mc)
       pp = FUNC_pp(r_mc,p_mc)
       ff = FUNC_ff(r_mc,p_mc)
       
       DO i_s = 1,NSMAX
          vOffsetA = 10*(i_s-1) 
          XvecIn( 6+vOffsetA,i_x2d) = nn(  i_s)/NnNF
          XvecIn( 7+vOffsetA,i_x2d) = ff(1,i_s)/FrNF
          XvecIn( 8+vOffsetA,i_x2d) = ff(2,i_s)/FbNF
          XvecIn( 9+vOffsetA,i_x2d) = ff(3,i_s)/FtNF
          XvecIn(10+vOffsetA,i_x2d) = ff(4,i_s)/FpNF
          XvecIn(11+vOffsetA,i_x2d) = pp(  i_s)/PpNF
          XvecIn(12+vOffsetA,i_x2d) = ff(5,i_s)/QrNF
          XvecIn(13+vOffsetA,i_x2d) = ff(6,i_s)/QbNF
          XvecIn(14+vOffsetA,i_x2d) = ff(7,i_s)/QtNF
          XvecIn(15+vOffsetA,i_x2d) = ff(8,i_s)/QpNF
       ENDDO
    ENDDO
    
    CALL T2_CONV
    CALL T2_GOUT
    
    RETURN
    
  END SUBROUTINE T2PROF_TOROIDAL_COORDINATE_AREA

  !------------------------------------------------------------------
  !
  !       Function for generating radial profile 
  !
  ! for 0 < \rho < 1
  !   f(\rho) = (fc-fs)*(1-\rho^n)^m + fs 
  !
  ! for 1 < \rho < \rho_wall
  !   f(\rho) = fw + a1*(\rho-\rho_wall)^2 
  !                + a2*(\rho-\rho_wall)^3 
  !                + a3*(\rho-\rho_wall)^4 
  !
  !  caution: \rho must be defined by length (\rho = \sqrt{\sigma})
  !
  !                     2014-02-23 H.Seto checked
  !
  !------------------------------------------------------------------
  FUNCTION FUNC_rprof(i0m0,i0n0,d0fc,d0fs,d0fw,d0rw,d0rr)
    
    INTEGER(ikind),INTENT(IN )::i0m0,i0n0
    REAL(   rkind),INTENT(IN )::d0fc,d0fs,d0fw,d0rw,d0rr
    REAL(   rkind),DIMENSION(1:2)::FUNC_rprof
    REAL(   rkind)::&
         d0a1,d0a2,d0a3,d0b1,d0b2,d0b3,&
         d0m0,d0n0,d0rx,d0rr,&
         d0fcs,d0fsw,d0rsw
    INTEGER(ikind)::&
         i0m1,i0n1,i0n2

    IF(i0n0.LT.2)THEN 
       WRITE(6,*)'INAPPROPRIATE PROFILE PARAMETER: n',i0n0
       STOP
    ENDIF
    
    d0n0  = DBLE(i0n0)
    d0m0  = DBLE(i0m0)
       
    d0fcs = d0fc - d0fs
    i0n1  = i0n0 - 1
    i0n2  = i0n0 - 2
    i0m1  = i0m0 - 1
    
    !  0 .LE. rho .LE. 1
    IF((d0rr.GE.0.D0).AND.(d0rr.LE.1.D0))THEN
       
       d0rx  = 1.D0-(d0rr**i0n0)
       
       SELECT CASE(i0m0)
       CASE(1)
          FUNC_rprof(1)    =  d0fcs*d0rx + d0fs
          IF(i0n2.NE.0)THEN
             FUNC_rprof(2) = -d0fcs*(d0rr**i0n2)*d0n0*0.5D0
          ELSE
             FUNC_rprof(2) = -d0fcs
          END IF
       CASE(2:)
          FUNC_rprof(1)    =  d0fcs*(d0rx**i0m0) + d0fs
          IF(i0n2.NE.0)THEN
             FUNC_rprof(2) = -d0fcs*(d0rr**i0n2)*(d0rx**i0m1)*d0m0*d0n0*0.5D0
          ELSE
             FUNC_rprof(2) = -d0fcs             *(d0rx**i0m1)*d0m0
          END IF
       END SELECT
    ELSEIF(d0rr.GT.1.D0)THEN
       d0rsw = 1.D0/(1.D0 - d0rw)
       d0fsw = d0fs - d0fw
       d0rx = d0rr - d0rw
       SELECT CASE(i0m0)
       CASE (1)
          d0b1 =  d0fsw*(d0rsw**2)
          d0b2 = -d0fcs* d0rsw    *d0n0
          d0b3 = -d0fcs           *d0n0*(d0n0-1.D0)
          
       CASE (2)
          d0b1 =  d0fsw*(d0rsw**2)
          d0b2 =  0.D0
          d0b3 =  d0fcs           *d0n0*d0n0*2.D0
       CASE (3:)
          d0b1 = d0fsw*(d0rsw**2)
          d0b2 = 0.D0
          d0b3 = 0.D0
       CASE DEFAULT
          WRITE(6,*)'INAPPROPRIATE PROFILE PARAMETER: m',i0m0
          STOP
       END SELECT
       
       d0a1 =  6.0D0*d0b1 - 3.0D0*d0b2 + 0.5D0*d0b3
       d0a2 = -8.0D0*d0b1 + 5.0D0*d0b2 - 1.0D0*d0b3
       d0a3 =  3.0D0*d0b1 - 2.0D0*d0b2 + 0.5D0*d0b3
       
       d0a2 = d0a2*d0rsw
       d0a3 = d0a3*d0rsw*d0rsw
       
       FUNC_rprof(1) =      d0a1*(d0rx**2) +      d0a2*(d0rx**3)&
            &   +      d0a3*(d0rx**4) + d0fw  
       FUNC_rprof(2) = 2.D0*d0a1* d0rx     + 3.D0*d0a2*(d0rx**2) &
            &   + 4.D0*d0a3*(d0rx**3) 
       FUNC_rprof(2) = FUNC_rprof(2)*0.5D0/d0rr
       
    ENDIF
    
    RETURN
    
  END FUNCTION FUNC_rprof
  
  !-------------------------------------------------------------------
  !
  ! INITITIAL PROFILE OF R (R,\phi,Z)
  !
  !                     2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_r_rz(r_mc,p_mc)
    
    USE T2COMM, ONLY: RR,RA
    
    REAL(rkind),INTENT(IN)::r_mc,p_mc
    REAL(rkind)::FUNC_r_rz
    
    FUNC_r_rz = RR + RA*SQRT(r_mc) * COS(p_mc)
    
    RETURN
    
  END FUNCTION FUNC_r_rz
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF Z (R,\phi, Z)
  !
  !                     2014-02-23 H.SETO checked
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_z_rz(r_mc,p_mc)
    
    USE T2COMM, ONLY: RR,RA
    
    REAL(rkind),INTENT(IN)::r_mc,p_mc
    REAL(rkind)::FUNC_z_rz
    
    FUNC_z_rz =        - RA*SQRT(r_mc)* SIN(p_mc)
    
    RETURN
    
  END FUNCTION FUNC_z_rz
  

  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF METRIC COEFFICIENTS
  !
  !                     2014-02-23 H.SETO checked
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_metric(r_mc,p_mc)
    
    USE T2COMM, ONLY: RR,RA
    
    REAL(rkind),INTENT(IN)::r_mc,p_mc
    REAL(rkind),DIMENSION(1:9)::FUNC_metric
    REAL(rkind)::r_rz

    !    R   = R0+a*\sqrt{\sig}*COS(\chi)
    !    phi = \zeta
    !    Z   =   -a*\sqrt{\sig}*SIN(\chi)
    !
    !    metric(1) : GRt    (\sqrt{g} )
    !    metric(2) : G11xCo (g_{\sig\sig}*\sig) 
    !    metric(3) : G12Co  (g_{\sig\chi} = 0 ) 
    !    metric(4) : G22Co  (g_{\chi\chi})
    !    metric(5) : G33Co  (g_{\zeta\zeta}) 
    !    metric(6) : G11Ct  (g^{\sig\sig}) 
    !    metric(7) : G12Ct  (g^{\sig\chi} = 0) 
    !    metric(8) : G22xCt (g^{\chi\chi}*\sig) 
    !    metric(9) : G33Ct  (g^{\zeta\zeta}) 
    
    r_rz   = FUNC_r_rz(r_mc,p_mc)
    FUNC_metric(1) = (RA**2)*r_rz*0.5D0
    FUNC_metric(2) = (RA**2)/4.D0
    FUNC_metric(3) = 0.D0 
    FUNC_metric(4) = (RA*2)*r_mc
    FUNC_metric(5) = r_rz**2
    FUNC_metric(6) = 4.D0*r_mc/(RA**2)
    FUNC_metric(7) = 0.D0
    FUNC_metric(8) = 1.D0/(RA**2)
    FUNC_metric(9) = 1.D0/(r_rz**2)
    
    RETURN
    
  END FUNCTION FUNC_metric

  !-------------------------------------------------------------------
  !
  ! INITITIAL PROFILE OF q
  !
  !                     2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_q(r_mc,p_mc)
    
    USE T2COMM, ONLY:d0qc,d0qs
    
    REAL(i0rkind),INTENT(IN)::r_mc,p_mc
    REAL(i0rkind)::FUNC_q
    
    IF((r_mc.GE.0.D0).AND.(r_mc.LE.1.D0))THEN
       FUNC_q = (d0qc-d0qs)*(1.D0 - r_mc)+d0qs
    ELSEIF(d0mfcr.GT.1.D0)THEN
       FUNC_q = (d0qs-d0qc)*(       r_mc)+d0qc
    ELSE
       WRITE(6,*)'WRONG RHO INPUT'
       STOP
    ENDIF
    
    RETURN
    
  END FUNCTION FUNC_q
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF dPsidr
  !
  !                     2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_dpsidr(r_mc,p_mc)
    
    USE T2COMM, ONLY:RA,RR,d0bc
    
    REAL(rkind),INTENT(IN)::r_mc,p_mc
    REAL(rkind)::FUNC_dPsidr
    REAL(rkind)::q,temp
    
    q    = FUNC_q(r_mc,p_mc)
    temp = SQRT(1.D0 - ((RA/RR)**2)*r_mc)
    
    FUNC_dPsidr = d0bc*(RA**2)/(2.D0*q*temp)
    
    RETURN
    
  END FUNCTION FUNC_dpsidr
  
  !-------------------------------------------------------------------
  !
  ! INITITIAL PROFILE OF I
  !
  !                     2014-06-04 H.SETO
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_BtCo(r_mc,p_mc)
    
    USE T2COMM,ONLY:RR,d0bc
    
    REAL(i0rkind),INTENT(IN)::r_mc,p_mc
    REAL(i0rkind)::FUNC_BtCo
    
    FUNC_BtCo = RR*d0bc
    
    RETURN
    
  END FUNCTION FUNC_BtCo
  
  !-------------------------------------------------------------------
  !
  ! INITITIAL PROFILE OF B
  !
  !                     2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_bb(r_mc,p_mc)
    
    REAL(rkind),INTENT(IN)::r_mc,p_mc
    REAL(rkind)::FUNC_Bb
    REAL(rkind),DIMENSION(1:9)::metric
    REAL(rkind)::dPsidr,btCo,g11Ct,g33Ct,bbSq,FUNC_bb
    
    dPsidr  = FUNC_dPsidr(r_mc,p_mc)
    btCo    = FUNC_btCo(  r_mc,p_mc)
    metric  = FUNC_metric(r_mc,p_mc)
    g11Ct   = metric(6)
    g33Ct   = metric(9)
    bbSq    = dPsidr*dPsidr*g11Ct*g33Ct + btCo*btCo*g33Ct
    FUNC_Bb = SQRT(bbSq) 
    
    RETURN

  END FUNCTION FUNC_Bb
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF E_{\zeta}
  !
  !                     2014-06-04 H.Seto
  ! 
  !-------------------------------------------------------------------
  FUNCTION FUNC_etCo(r_mc,p_mc)
    
    USE T2CNST,ONLY:Aee,RR,RA
    USE T2COMM,ONLY:NSMAX
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind)::FUNC_etCo
    REAL(   rkind),DIMENSION(1:NSMAX)::tt
    REAL(   rkind)::ttE_keV,lnLamb,eta,jtCo,d0temp
    
    tt      = FUNC_tt(r_mc,p_mc)
    ttE_keV = tt(1)/(1.D3*Aee)   ! electron temperature in keV
    lnLamb  = 17.D0 ! Coulomb logarithm for debug 
    jtCo    = FUNC_jtCo1d(r_mc,p_mc)
    d0eta = (1.65D-9)*lnLamb/(SQRT(ttE_keV)**3)
    d0eta = d0eta/((1.D0-SQRT(RA*SQRT(r_mc)/RR))**2)
    
    fd0coet = d0eta*d0jt
    
    RETURN
    
  END FUNCTION fd0coet
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF \bar{E}_{\chi}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_epCo(r_mc,p_mc)
    
    REAL(rkind),INTENT(IN)::r_mc,p_mc
    REAL(rkind)::FUNC_epCo
    
    FUNC_epCo = 0.D0
    
    RETURN
    
  END FUNCTION FUNC_epCo
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF E_{\rho}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_erCo(r_mc,p_mc)
    
    USE T2CNST,ONLY:Aee
    USE T2COMM,ONLY:NSMAX
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind)::FUNC_erCo
    INTEGER(ikind)::i_s
    REAL(   rkind),DIMENSION(1:NSMAX)::dpdr,nn
    REAL(   rkind)::nnE
    
    FUNC_erCo = 0.D0
    
    dpdr = FUNC_dpdr(r_mc,p_mc)
    nn   = FUNC_nn(  r_mc,p_mc)
    nnE  = nn(1)     
    
    DO i_s = 2,NSMAX
       FUNC_erCo = FUNC_erCo + dpdr(i_s)
    ENDDO
    
    FUNC_erCo = FUNC_erCo/(Aee*nnE)
    
    RETURN
    
  END FUNCTION FUNC_erCo
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF n_{a}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_nn(r_mc,p_mc)
    
    USE T2COMM,ONLY:NSMAX,i0nm,i0nn,d1nc,d1ns,d1nw,d0rw
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:i0smax)::FUNC_nn
    REAL(   rkind),DIMENSION(1:2)::d1rn
    REAL(   rkind)::d0nc,d0ns,d0nw,rho
    INTEGER(ikind)::i_s
        
    rho = SQRT(r_mc)
    DO i_s = 1,NSMAX
       
       d0nc = d1nc(i_s)
       d0ns = d1ns(i_s)
       d0nw = d1nw(i_s)

       ! density in 10^{20} m^{-3}
       d1rn = FUNC_rprof(i0nm,i0nn,d0nc,d0ns,d0nw,d0rw,rho)
       ! density in m^{-3}
       FUNC_nn(i_s) = d1rn(1)*1.D20
       
    ENDDO
    
    RETURN

  END FUNCTION FUNC_nn

  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF dn_{a}/d\rho
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_dndr(r_mc,p_mc)
    
    USE T2COMM, ONLY:NSMAX,i0nm,i0nn,d1nc,d1ns,d1nw,d0rw
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:NSMAX)::FUNC_dndr
    REAL(   rkind),DIMENSION(1:2)::d1rn
    REAL(   rkind)::d0nc,d0ns,d0nw,rho
    INTEGER(ikind)::i_s
    
    rho = SQRT(r_mc)
    DO i_s = 1, NSMAX
       
       d0nc = d1nc(i_s)
       d0ns = d1ns(i_s)
       d0nw = d1nw(i_s)
       
       ! d n_{a} /d \rho in 10^{20} m^{-3}
       d1rn = FUNC_rprof(i0nm,i0nn,d0nc,d0ns,d0nw,d0rw,rho)
       ! d n_{a} /d \rho in  m^{-3}
       FUNC_dndr(i_s) = d1rn(2)*1.D20
       
    ENDDO
    
    RETURN
    
  END FUNCTION FUNC_dndr
 
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF T_{a}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_tt(r_mc,p_mc)
    
    USE T2CNST,ONLY:Aee
    USE T2COMM,ONLY:NSMAX,i0tm,i0tn,d1tc,d1ts,d1tw,d0rw
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:NSMAX)::FUNC_tt
    REAL(   rkind),DIMENSION(1:2)::d1rt
    REAL(   rkind)::d0tc,d0ts,d0tw,rho
    INTEGER(ikind)::i_s
    
    rho = SQRT(r_mc)    
    DO i_s = 1, NSMAX
       
       d0tc = d1tc(i_s)
       d0ts = d1ts(i_s)
       d0tw = d1tw(i_s)
       
       ! T_{a} in keV
       d1rt = FUNC_rprof(i0tm,i0tn,d0tc,d0ts,d0tw,d0rw,rho)      
       ! T_{a} in J      
       FUNC_tt(i_s) = d1rt(1)*1.D3*Aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION FUNC_tt

  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF dT_{a}/d\sig
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------  
  FUNCTION FUNC_dtdr(r_mc,p_mc)
    
    USE T2CNST, ONLY:Aee
    USE T2COMM, ONLY:NSMAX,i0tm,i0tn,d1tc,d1ts,d1tw,d0rw
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:NSMAX)::FUNC_dtdr
    REAL(   rkind),DIMENSION(1:2)::d1rt
    REAL(   rkind)::d0tc,d0ts,d0tw,rho
    INTEGER(ikind)::i_s
    
    rho =SQRT(r_mc)
    
    DO i_s = 1, NSMAX
       
       d0tc = d1tc(i_s)
       d0ts = d1ts(i_s)
       d0tw = d1tw(i_s)
       
       ! dT_{a}/d\rho in keV
       d1rt = FUNC_rprof(i0tm,i0tn,d0tc,d0ts,d0tw,d0rw,rho)
       
       ! dT_{a}/d\sig in J
       FUNC_dtdr(i_s) = d1rt(2)*1.D3*Aee
       
    ENDDO
    
    RETURN
    
  END FUNCTION FUNC_dtdr
  
  !-------------------------------------------------------------------
  !
  ! INITITIAL PROFILE OF p_{a}
  !
  !                     2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_pp(r_mc,p_mc)
    
    USE T2COMM, ONLY:NSMAX
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:NSMAX)::FUNC_pp,nn,tt
    INTEGER(ikind)::i_s
    
    nn = FUNC_nn(r_mc,p_mc)
    tt = FUNC_tt(r_mc,p_mc)
    
    DO i_s = 1, NSMAX  
       ! p_{a} in J/m^{3}
       FUNC_pp(i_s) = nn(i_s)*tt(i_s)
    ENDDO
    
    RETURN
    
  END FUNCTION FUNC_pp
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF dp_{a}/d\sig
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_dpdr(r_mc,p_mc)
    
    USE T2COMM,ONLY:NSMAX
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:NSMAX)::FUNC_dpdr
    REAL(   rkind),DIMENSION(1:NSMAX)::nn,dndr,tt,dtdr
    INTEGER(ikind)::i_s
    
    nn   = FUNC_nn(  r_mc,p_mc)    
    dndr = FUNC_dndr(r_mc,p_mc)
    
    tt(  1:NSMAX) = FUNC_tt(  r_mc,p_mc)
    dtdr(1:NSMAX) = FUNC_dtdr(r_mc,p_mc)
    

    DO i_s = 1, NSMAX
       ! dp_{a}/d\rho in J/m^{3}
       FUNC_dpdr(i_s) = nn(i_s)*dtdr(i_s) + dndr(i_s)*tt(i_s)
    ENDDO
    
    RETURN
    
  END FUNCTION FUNC_dpdr
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF <j_{\zeta}>
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_jtCo1d(r_mc,p_mc)
    
    USE T2COMM, ONLY:NSMAX,RA,RR
    
    INTEGER(ikind)::i_s
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind)::FUNC_jtCo1d
    REAL(   rkind),DIMENSION(1:NSMAX)::dpdr
    REAL(   rkind)::dpsidr
    
    FUNC_jtCo1d = 0.D0
    dpsidr  = Func_dpsidr(r_mc,p_mc)
    dpdr    = FUNC_dpdr(  r_mc,p_mc)
    
    DO i_s = 1, NSMAX
       FUNC_jtCo1d = FUNC_jtCo1d + dpdr(i_s)
    ENDDO
    
    FUNC_jtCo1d = -FUNC_jtCo1d*((RR**2)+1.5D0*(RA**2)*r_mc)/dpsidr
    
    RETURN
    
  END FUNCTION FUNC_jtCo1d
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF j_{\zeta}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_jtCo(r_mc,p_mc)
    
    USE T2COMM,ONLY:NSMAX
    
    INTEGER(ikind)::i_s
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind)::FUNC_jtCo
    REAL(   rkind),DIMENSION(1:NSMAX)::dpdr
    REAL(   rkind)::dpsidp,r_rz
    
    FUNC_jtCo = 0.D0
    
    r_rz   = FUNC_z_rz(  r_mc,p_mc)
    dpsidr = FUNC_dpsidr(r_mc,p_mc)
    dpdr   = FUNC_dpdr(  r_mc,p_mc)
    
    DO i_s = 1, NSMAX
       FUNC_jtCo = FUNC_jtCo + dpdr(i_s)
    ENDDO
    
    FUNC_jtCo = -FUNC_jtCo*(r_rz**2)/dpsidr
    
    RETURN
    
  END FUNCTION FUNC_jtCo
  
  !-------------------------------------------------------------------
  !
  !       INITITIAL PROFILE OF j_{\para}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_jb(r_mc,p_mc)
    
    USE T2COMM, ONLY:NSMAX
    
    INTEGER(ikind)::i_s
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind)::FUNC_jb
    REAL(   rkind),DIMENSION(1:NSMAX)::dpdr
    REAL(   rkind)::btCo,dpsidr,bb
    
    FUNC_jb = 0.D0
    
    btCo   = FUNC_btCo(  r_mc,p_mc)
    dpsidr = FUNC_dpsidr(r_mc,p_mc)
    dpdr   = FUNC_dpdr(  r_mc,p_mc)
    bb     = FUNC_bb(    r_mc,p_mc)
    
    DO i_s = 1, NSMAX
       FUNC_jb = FUNC_jb + dpdr(i_s)
    ENDDO
    
    FUNC_jb = -FUNC_jb*btCo/(bb*dpsidr)
    
    RETURN
    
  END FUNCTION FUNC_jb


  !-------------------------------------------------------------------
  !
  ! INITITIAL PROFILE OF \bar{Gamma}_{a}^{\rho}, 
  !                      \Gamma_{a\para},
  !                      \Gamma_{a\zeta},
  !                      \Gamma_{a}^{\chi},
  !                      \bar{Q}_{a}^{\rho}
  !                      Q_{a\para}
  !                      Q_{a\zeta}
  !                      Q_{a}^{\chi}
  !
  !                 LAST UPDATE    2014-06-04 H.Seto
  !
  !-------------------------------------------------------------------
  FUNCTION FUNC_flux(r_mc,p_mc)
    
    USE T2CNST, ONLY:Aee
    USE T2COMM, ONLY:NSMAX
    
    REAL(   rkind),INTENT(IN)::r_mc,p_mc
    REAL(   rkind),DIMENSION(1:8,1:NSMAX)::FUNC_flux
    REAL(   rkind),DIMENSION(    1:NSMAX)::tt
    REAL(   rkind)::fb,ftCo,ttE
    INTEGER(ikind)::i_s
    
    tt   =  FUNC_tt(  r_mc,p_mc)
    fb   = -FUNC_jb(  r_mc,p_mc)/Aee
    ftCo = -FUNC_jtCo(r_mc,p_mc)/Aee
    ttE  =  tt(1)
    
    ! electron
    FUNC_flux(1,i_s) = 0.D0
    FUNC_flux(2,i_s) = fb
    FUNC_flux(3,i_s) = ftCo
    FUNC_flux(4,i_s) = 0.D0
    FUNC_flux(5,i_s) = 0.D0
    FUNC_flux(6,i_s) = 2.5D0*ttE*fb
    FUNC_flux(7,i_s) = 2.5D0*ttE*ftCo
    FUNC_flux(8,i_s) = 0.D0
    ! ions
    FUNC_flux(1:8,2:NSMAX) = 0.D0
    
    RETURN
    
  END FUNCTION FUNC_flux
END MODULE T2PREP
