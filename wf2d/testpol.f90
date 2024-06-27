! testpol.f90

MODULE testpol

  PRIVATE
  PUBLIC test_pol

CONTAINS
  
  !  calculate local wave number and polarization vector

  SUBROUTINE test_pol

    USE plcomm
    USE plprof
    USE plprofw
    USE dppola

    IMPLICIT NONE
    INTEGER:: ns
    REAL(rkind):: f,x,y,z,rny,rnz
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw

    modelg= 2
    RR    = 0.64D0
    RA    = 0.40D0
    RB    = 0.48D0
    BB    = 0.25D0
    NSMAX=2
    pn(1)=0.0015D0
    pn(2)=.0015D0
    pns(1)=0.0000D0
    pns(2)=0.0000D0
    PTPP(1) =  0.3D0
    PTPP(2) =  0.3D0
    PTPR(1) =  0.3D0
    PTPR(2) =  0.3D0
    PTS(1)  =  0.003D0
    PTS(2)  =  0.03D0
    PZCL(1) =  0.000D0
    PZCL(2) =  0.000D0

    f=8560.D0
    x=0.399D0
    y=0.D0
    z=0.D0
    rny=0.0
    rnz=0.0

1   CONTINUE
    WRITE(6,*) '## input f[MHz],x,y,z,rny,rnz ?'
    READ(5,*,END=9000,ERR=1) f,x,y,z,rny,rnz
    
    CALL pl_mag(x,y,z,mag)
    WRITE(6,'(A14,3E12.4)') 'x,y,z=        ',x,y,z
    WRITE(6,'(A14,2E12.4)') 'BABS,RHON=    ',mag%babs,mag%rhon
    WRITE(6,'(A14,3E12.4)') 'BNX,BNY,BNZ=  ',mag%bnx,mag%bny,mag%bnx
    WRITE(6,'(A14,3E12.4)') 'BX,BY,BZ=     ', &
         mag%babs*mag%bnx,mag%babs*mag%bny,mag%babs*mag%bnx
    CALL pl_profw(mag%rhon,plfw)
    DO ns=1,nsmax
       WRITE(6,'(A8,I4,3ES12.4)') &
            'ns,n,tpr,tpp= ',ns,plfw(ns)%rn,plfw(ns)%rtpr,plfw(ns)%rtpp
    END DO

    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE test_pol
END MODULE testpol
    
