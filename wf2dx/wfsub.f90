! wfsub.f90

MODULE wfsub

  PRIVATE
  PUBLIC wf_set_node_range
  PUBLIC wf_set_elm_area
!  PUBLIC wfclass
  PUBLIC wf_set_weight
  PUBLIC wf_set_abc
  PUBLIC wf_set_node
  PUBLIC wf_set_aif
  PUBLIC wf_set_aie
  PUBLIC wf_set_ewg
  PUBLIC wf_set_lside
  PUBLIC wf_cross
  PUBLIC wf_fieldcr
  PUBLIC wf_fieldcz
  PUBLIC wf_fieldcp
  

CONTAINS

!     $Id: wfsub.f90,v 1.16 2011/11/16 09:18:22 maruyama Exp $

!     ****** SETUP NODE RANGE ******

SUBROUTINE wf_set_node_range

  use wfcomm
  implicit none

  integer :: IN
  real(rkind) :: LNODE

  RNDMIN=xnode(1)
  RNDMAX=xnode(1)
  ZNDMIN=ynode(1)
  ZNDMAX=ynode(1)
  SELECT CASE(MODELG)
  CASE(0,1,12)
     LNODE=SQRT(xnode(1)**2+ynode(1)**2)
  CASE(2)
     LNODE=SQRT((xnode(1)-RR)**2+ynode(1)**2)
  END SELECT
  LNDMIN=LNODE
  LNDMAX=LNODE
     
  DO IN=2,node_max
     RNDMIN=MIN(RNDMIN,xnode(IN))
     RNDMAX=MAX(RNDMAX,xnode(IN))
     ZNDMIN=MIN(ZNDMIN,ynode(IN))
     ZNDMAX=MAX(ZNDMAX,ynode(IN))
     SELECT CASE(MODELG)
     CASE(0,1,11,12)
        LNODE=SQRT(xnode(IN)**2+ynode(IN)**2)
     CASE(2)
        LNODE=SQRT((xnode(IN)-RR)**2+ynode(IN)**2)
     END SELECT
     LNDMIN=MIN(LNDMIN,LNODE)
     LNDMAX=MAX(LNDMAX,LNODE)
  ENDDO

  RETURN
END SUBROUTINE wf_set_node_range

!     ****** SETUP ELEMENT AREA ******

SUBROUTINE wf_set_elm_area

  use wfcomm
  implicit none

  integer :: NE
  real(rkind) :: RE(3),ZE(3),S

  do NE=1,nelm_max

     call wf_set_node(NE,RE,ZE)
     S=RE(1)*(ZE(2)-ZE(3))+RE(2)*(ZE(3)-ZE(1))+RE(3)*(ZE(1)-ZE(2))

     if(S.LE.0.D0) THEN
        if(nrank.eq.0) then
           write(6,'(A,I5)') 'NEGATIVE S: NE=',NE
           write(6,'(A,1P,3E12.4)') 'NODE 1 = ',RE(1),ZE(1)
           write(6,'(A,1P,3E12.4)') 'NODE 2 = ',RE(2),ZE(2)
           write(6,'(A,1P,3E12.4)') 'NODE 3 = ',RE(3),ZE(3)
        end if
     end if
     
     SELM(NE)=0.5D0*S

  end do

  RETURN
END SUBROUTINE wf_set_elm_area

!     ******* Classify point location *******

SUBROUTINE WFCLASS(NE,R,Z,WGT,IND)
!
!   IND = 0 : point inside the element
!         1 : on the edge 1
!         2 : on the edge 2
!         3 : on the edge 3
!         4 : on the node 1
!         5 : on the node 2
!         6 : on the node 3
!         7 : inconsistent (area of element is zero)
!        -1 : outside of the element
!  
  use wfcomm
  implicit none
  integer,intent(in) :: NE
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: WGT(3)
  integer,intent(out) :: IND
  REAL(rkind),PARAMETER:: eps=1.D-12

  CALL wf_set_weight(NE,R,Z,WGT)

  IF(WGT(1).GT.eps) THEN
     IF(WGT(2).GT.eps) THEN
        IF(WGT(3).GT.eps) THEN
           IND=0
        ELSE IF(WGT(3).GT.-eps) THEN
           IND=3
        END IF
     ELSE IF(WGT(2).GT.-eps) THEN
        IF(WGT(3).GT.eps) THEN
           IND=2
        ELSE IF(WGT(3).GT.-eps) THEN
           IND=4
        END IF
     END IF
  ELSE IF(WGT(1).GT.-eps) THEN
     IF(WGT(2).GT.eps) THEN
        IF(WGT(3).GT.eps) THEN
           IND=1
        ELSE IF(WGT(3).GT.-eps) THEN
           IND=5
        END IF
     ELSE IF(WGT(2).GT.-eps) THEN
        IF(WGT(3).GT.eps) THEN
           IND=6
        ELSE IF(WGT(3).GT.-eps) THEN
           IND=7
        END IF
     END IF
  END IF
  RETURN
END SUBROUTINE WFCLASS

!     ******* WEIGHT CALCULATION *******
!     WEIGHT MEANS AREA COORDINATE

SUBROUTINE wf_set_weight(NE,R,Z,WGT)
  
  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: IN
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: WGT(3)
  real(rkind) :: A(3),B(3),C(3)
  
  call wf_set_abc(NE,A,B,C)

  DO IN=1,3
     WGT(IN)=A(IN)+B(IN)*R+C(IN)*Z
  ENDDO

  RETURN
END SUBROUTINE wf_set_weight

!     ******* A,B,C CALCULATION *******
!     A,B,C ARE USED BY AREA COORDINATE 

SUBROUTINE wf_set_abc(NE,A,B,C)

  use wfcomm
  implicit none
  integer,intent(in):: NE
  integer :: I,J,K
  real(rkind),intent(out)::A(3),B(3),C(3) 
  real(rkind) :: RE(3),ZE(3),S
  
  CALL wf_set_node(NE,RE,ZE)
  S=SELM(NE)

  do I=1,3
     J=I+1
     K=I+2
     if(J.gt.3) J=J-3
     if(K.gt.3) K=K-3
     
     A(I)=0.5D0*(RE(J)*ZE(K)-RE(K)*ZE(J))/S
     B(I)=0.5D0*(ZE(J)-ZE(K))/S
     C(I)=0.5D0*(RE(K)-RE(J))/S
  end do
    
  RETURN
END SUBROUTINE wf_set_abc

!     ******* TOTAL COORDINATE - LOCAL COORDINATE *******

SUBROUTINE wf_set_node(NE,RE,ZE)
  
  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: IN,NN
  real(rkind),intent(out):: RE(3),ZE(3)
  
  DO IN=1,3
     NN=node_nside_nelm(IN,NE)
     RE(IN)=xnode(NN)
     ZE(IN)=ynode(NN)
  END DO

  RETURN
END SUBROUTINE wf_set_node

!     *******  INITIALISE SURFACE INTEGRAL OF ELEMENT FUNCTION *******

SUBROUTINE wf_set_aif

  use wfcomm
  implicit none
  integer :: ID(3,3),I,L1,L2,L3,J,K

  DATA ID/1,3*0,1,3*0,1/
  
  DO I=1,3
     L1=ID(1,I)
     L2=ID(2,I)
     L3=ID(3,I)
     AIF1(I)=AIF(L1,L2,L3)
     DO J=1,3
        L1=ID(1,I)+ID(1,J)
        L2=ID(2,I)+ID(2,J)
        L3=ID(3,I)+ID(3,J)
        AIF2(I,J)=AIF(L1,L2,L3)
        DO K=1,3
           L1=ID(1,I)+ID(1,J)+ID(1,K)
           L2=ID(2,I)+ID(2,J)+ID(2,K)
           L3=ID(3,I)+ID(3,J)+ID(3,K)
           AIF3(I,J,K)=AIF(L1,L2,L3)
        ENDDO
     ENDDO
  ENDDO
 
  RETURN
END SUBROUTINE wf_set_aif

!     ******* SURFACE INTEGRAL OF ELEMENT FUNCTION *******
!     KAI(X) = X! (X FACTORIAL)
!     WHEN YOU USE AIF, YOU SHOULD MULTIPLY AREA OF ELEMNT(S)

FUNCTION AIF(L1,L2,L3)

  USE bpsd_kinds,ONLY: rkind
  implicit none
  integer :: KAI(0:10),L1,L2,L3
  real(rkind) :: AIF
  DATA KAI/1,1,2,6,24,120,720,5040,40320,362880,3628800/
  
  AIF=DBLE(2*KAI(L1)*KAI(L2)*KAI(L3))/DBLE(KAI(L1+L2+L3+2))

  RETURN
END FUNCTION AIF

!     *******  INITIALISE LINE INTEGRAL OF ELEMENT FUNCTION *******

SUBROUTINE wf_set_aie

  use wfcomm
  implicit none
  integer :: ID(3,3),I,L1,L2,L3,J,K

  DATA ID/1,0,0,0,1,0,0,0,1/
  
  DO I=1,3
     L1=ID(1,I)
     L2=ID(2,I)
     L3=ID(3,I)
     AIE1(I)=AIE(L1,L2,L3)
     DO J=1,3
        L1=ID(1,I)+ID(1,J)
        L2=ID(2,I)+ID(2,J)
        L3=ID(3,I)+ID(3,J)
        AIE2(I,J)=AIE(L1,L2,L3)
        DO K=1,3
           L1=ID(1,I)+ID(1,J)+ID(1,K)
           L2=ID(2,I)+ID(2,J)+ID(2,K)
           L3=ID(3,I)+ID(3,J)+ID(3,K)
           AIE3(I,J,K)=AIE(L1,L2,L3)
        ENDDO
     ENDDO
  ENDDO
 
  RETURN
END SUBROUTINE wf_set_aie

!     ******* LINE INTEGRAL OF ELEMENT FUNCTION *******
!     KAI(X) = X! (X FACTORIAL)
!     WHEN YOU USE AIE, YOU SHOULD MULTIPLY AREA OF ELEMNT(S)

FUNCTION AIE(L1,L2,L3)

  USE bpsd_kinds,ONLY: rkind
  implicit none
  integer :: KAI(0:10),L1,L2,L3
  real(rkind) :: AIE
  DATA KAI/1,1,2,6,24,120,720,5040,40320,362880,3628800/
  
  AIE=DBLE(KAI(L1)*KAI(L2)*KAI(L3))/DBLE(KAI(L1+L2+L3+1))

  RETURN
END FUNCTION AIE


!     ***** SET BOUNDARY ELECTRIC FIELD *****

SUBROUTINE wf_set_ewg

  USE wfcomm
  USE wfload,ONLY: wf_read_wg 
  implicit none
  INTEGER:: nseg,NN1,NN2,NN,NBSD,NBND,IERR
  REAL(rkind):: ANGLE,X,Y,PHASE,PROD,FACTOR,SN
  COMPLEX(rkind):: CEX,CEY,CEZ
  REAL(rkind),PARAMETER:: EPSWG=1.D-12

  ANGLE=ANGWG*PI/180.D0

! --- WG Electric field on boundary side ---  

  nseg_bdy_max=0
  DO nseg=1,nseg_max
     IF(KASID(nseg).EQ.1) nseg_bdy_max=nseg_bdy_max+1
  END DO
  IF(nrank.EQ.0) write(6,*) '## setewg: nseg_bdy_max=',nseg_bdy_max
  IF(ALLOCATED(NSDBS)) DEALLOCATE(NSDBS)
  IF(ALLOCATED(CEBSD)) DEALLOCATE(CEBSD)
  ALLOCATE(NSDBS(nseg_bdy_max))
  ALLOCATE(CEBSD(nseg_bdy_max))
  nseg_bdy_max=0
  DO nseg=1,nseg_max
     IF(KASID(nseg).EQ.1) THEN
        nseg_bdy_max=nseg_bdy_max+1
        NSDBS(nseg_bdy_max)=nseg
        KBSID(nseg)=nseg_bdy_max
     ELSE
        KBSID(nseg)=0
     END IF
  END DO
  DO NBSD=1,nseg_bdy_max
     nseg=NSDBS(NBSD)
     NN1=node_nseg(1,nseg)
     NN2=node_nseg(2,nseg)
     X=0.5D0*(xnode(NN1)+xnode(NN2))
     Y=0.5D0*(ynode(NN1)+ynode(NN2))
     SELECT CASE(MODELWG)
     CASE(0,1)
        IF((X.GE.X1WG).AND.(X.LE.X2WG).AND. &
           (Y.GE.Y1WG).AND.(Y.LE.Y2WG)) THEN
           WRITE(6,'(A,2ES12.4)') 'X,Y in WG:',X,Y
           PROD=(X2WG-X1WG)*(xnode(NN2)-xnode(NN1)) &
               +(Y2WG-Y1WG)*(ynode(NN2)-ynode(NN1))
           IF(ABS(X1WG-X2WG).LT.1.D-8) THEN
              FACTOR=(Y-0.5D0*(Y1WG+Y2WG))**2/(Y1WG-Y2WG)**2
           ELSE IF(ABS(Y1WG-Y2WG).LT.1.D-8) THEN
              FACTOR=(X-0.5D0*(X1WG+X2WG))**2/(X1WG-X2WG)**2
           ELSE
              FACTOR=(X-0.5D0*(X1WG+X2WG))**2/(X1WG-X2WG)**2 &
                    +(Y-0.5D0*(Y1WG+Y2WG))**2/(Y1WG-Y2WG)**2
           END IF
           SN=SQRT((X   -X1WG)**2+(Y   -Y1WG)**2) &
             /SQRT((X2WG-X1WG)**2+(Y2WG-Y1WG)**2) ! SN=0 at 1, 1 at 2
           PHASE=(PH1WG+(PH2WG-PH1WG)*SN+DPHWG*4.D0*SN*(1.D0-SN))*PI/180.D0
           CEBSD(NBSD)= AMPWG*EXP(CII*PHASE)*(SIN(ANGLE)-CII*ELPWG*COS(ANGLE))
           IF(MODELWG.EQ.1) CEBSD(NBSD)=CEBSD(NBSD)*EXP(-10.D0*FACTOR)
           IF(PROD.GT.0.D0) CEBSD(NBSD)=-CEBSD(NBSD)
           WRITE(71,'(A,2I6,3ES12.4)') 'nbsd:',NBSD,nseg,X,Y,FACTOR
           WRITE(71,'(A,4ES12.4)')     '     ',SN,PHASE,EXP(CII*PHASE)
           WRITE(71,'(A,4ES12.4)') &
                '     ',ANGLE,ELPWG,(SIN(ANGLE)-CII*ELPWG*COS(ANGLE))
           WRITE(71,'(A,4ES12.4)') &
                '     ',AMPWG,EXP(-10.D0*FACTOR),CEBSD(NBSD)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                WRITE(6,'(A,2I8,1P5E12.4)') &
                'SD:',nseg,NBSD,CEBSD(NBSD),AMPWG,PHASE,ANGLE
        ELSE
           CEBSD(NBSD)=(0.D0,0.D0)
        END IF
     CASE(12)
        IF((X.GE.X1WG-EPSWG).AND.(X.LE.X2WG+EPSWG).AND. &
           (Y.GE.Y1WG-EPSWG).AND.(Y.LE.Y2WG+EPSWG)) THEN
           PROD=(X2WG-X1WG)*(xnode(NN2)-xnode(NN1)) &
               +(Y2WG-Y1WG)*(ynode(NN2)-ynode(NN1))
           CALL wf_read_wg(Y,CEX,CEY,CEZ,IERR)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                write(6,'(A,1P6E12.4)') 'R,Z,CEY=', &
                                    X,Y,CEY,PROD,ynode(NN2)-ynode(NN1)
!!!           IF(PROD.GT.0.D0) CEY=-CEY
           CEBSD(NBSD)=AMPWG*CEY
        ELSE
           CEBSD(NBSD)=(0.D0,0.D0)
        END IF
     END SELECT
  END DO

! --- WG Electric field on boundary node ---  

  node_bdy_max=0
  DO NN=1,node_max
     IF(KANOD(NN).EQ.1) node_bdy_max=node_bdy_max+1
  END DO
  IF(nrank.EQ.0) write(6,*) '## setewg: node_bdy_max=',node_bdy_max
  IF(ALLOCATED(NNDBS)) DEALLOCATE(NNDBS)
  IF(ALLOCATED(CEBND)) DEALLOCATE(CEBND)
  ALLOCATE(NNDBS(node_bdy_max))
  ALLOCATE(CEBND(node_bdy_max))
  node_bdy_max=0
  DO NN=1,node_max
     IF(KANOD(NN).EQ.1) THEN
        node_bdy_max=node_bdy_max+1
        NNDBS(node_bdy_max)=NN
        KBNOD(NN)=node_bdy_max
     ELSE
        KBNOD(NN)=0
     END IF
  END DO
  DO NBND=1,node_bdy_max
     NN=NNDBS(NBND)
     X=xnode(NN)
     Y=ynode(NN)
     SELECT CASE(MODELWG)
     CASE(0,1)
        IF((X.GE.X1WG).AND.(X.LE.X2WG).AND. &
           (Y.GE.Y1WG).AND.(Y.LE.Y2WG)) THEN
           IF(ABS(X1WG-X2WG).LT.1.D-8) THEN
              FACTOR=(Y-0.5D0*(Y1WG+Y2WG))**2/(Y1WG-Y2WG)**2
           ELSE IF(ABS(Y1WG-Y2WG).LT.1.D-8) THEN
              FACTOR=(x-0.5D0*(X1WG+X2WG))**2/(X1WG-X2WG)**2
           ELSE
              FACTOR=(X-0.5D0*(X1WG+X2WG))**2/(X1WG-X2WG)**2 &
                    +(Y-0.5D0*(Y1WG+Y2WG))**2/(Y1WG-Y2WG)**2
           END IF
           SN=SQRT((X   -X1WG)**2+(Y   -Y1WG)**2) &
                /SQRT((X2WG-X1WG)**2+(Y2WG-Y1WG)**2) ! SN=0 at 1, 1 at 2
           PHASE=(PH1WG+(PH2WG-PH1WG)*SN+DPHWG*4.D0*SN*(1.D0-SN))*PI/180.D0
           CEBND(NBND)= AMPWG*EXP(CII*PHASE)*(COS(ANGLE)+CII*ELPWG*SIN(ANGLE))
           IF(MODELWG.EQ.1) CEBND(NBND)=CEBND(NBND)*EXP(-10.D0*FACTOR)
           WRITE(71,'(A,2I6,3ES12.4)') 'nbnd:',NBND,NN,X,Y,FACTOR
           WRITE(71,'(A,4ES12.4)')     '     ',SN,PHASE,EXP(CII*PHASE)
           WRITE(71,'(A,4ES12.4)') &
                '     ',ANGLE,ELPWG,(COS(ANGLE)+CII*ELPWG*SIN(ANGLE))
           WRITE(71,'(A,4ES12.4)') &
                '     ',AMPWG,EXP(-10.D0*FACTOR),CEBND(NBND)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                WRITE(6,'(A,2I8,1P5E12.4)') &
                'ND:',NN,NBND,CEBND(NBND),AMPWG,PHASE,ANGLE
        ELSE
           CEBND(NBND)=(0.D0,0.D0)
        END IF
     CASE(12)
        IF((X.GE.X1WG-EPSWG).AND.(X.LE.X2WG+EPSWG).AND. &
           (Y.GE.Y1WG-EPSWG).AND.(Y.LE.Y2WG+EPSWG)) THEN
           CALL wf_read_wg(Y,CEX,CEY,CEZ,IERR)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                write(6,'(A,1P4E12.4)') 'X,Y,CEZ=',X,Y,CEZ
           CEBND(NBND)=AMPWG*CEZ
        ELSE
           CEBND(NBND)=(0.D0,0.D0)
        END IF
     END SELECT
     END DO
  RETURN
END SUBROUTINE wf_set_ewg

!     ***** SET LENGTH OF SIDE *****

SUBROUTINE wf_set_lside
  use wfcomm
  implicit none

  integer :: nseg,ND1,ND2
  real(rkind) :: R1,R2,Z1,Z2,L

  do nseg=1,nseg_max
     LSID(nseg)=0.d0
  end do

  do nseg=1,nseg_max
     ND1=node_nseg(1,nseg)
     ND2=node_nseg(2,nseg)
     R1 =xnode(ND1)
     R2 =xnode(ND2)
     Z1 =ynode(ND1)
     Z2 =ynode(ND2)

     L =SQRT((R2-R1)**2+(Z2-Z1)**2)
     LSID(nseg)=L
  end do
  return
end SUBROUTINE wf_set_lside

!     ******* CALCULATE POINT OF INTERSECTION *******

SUBROUTINE wf_cross(R1,Z1,R2,Z2,IE,L,RC,ZC,IERR)
  
  use wfcomm
  implicit none
  integer,intent(in) :: IE,L
  integer,intent(out):: IERR
  integer :: M,N1,N2
  real(rkind),intent(in) :: R1,R2,Z1,Z2
  real(rkind),intent(out):: RC,ZC
  real(rkind),parameter :: EPS = 1.d-12
  real(rkind) :: R3,R4,Z3,Z4
  real(rkind) :: DELT,R12,R34,Z12,Z34,AD,RK,RT
  
  M=L+1
  if(M.gt.3) M=M-3
  N1=node_nside_nelm(L,IE)
  N2=node_nside_nelm(M,IE)
  R3=xnode(N1)
  Z3=ynode(N1)
  R4=xnode(N2)
  Z4=ynode(N2)

  R12=R1-R2
  Z12=Z1-Z2
  R34=R3-R4
  Z34=Z3-Z4
  DELT=Z12*R34-R12*Z34
  AD=ABS(DELT)
  IF(AD.LT.EPS) GOTO 9000

  RK=(-Z34*(R4-R2)+R34*(Z4-Z2))/DELT
  RT=(-Z12*(R4-R2)+R12*(Z4-Z2))/DELT
  IF(    RK.LT.0.D0.OR.RK.GT.1.D0 &
  &  .OR.RT.LT.0.D0.OR.RT.GT.1.D0) GOTO 9100

  RC=RK*R12+R2
  ZC=RK*Z12+Z2
  IERR=0
  RETURN

9000 IERR=9000
  RETURN

9100 IERR=9100
  RETURN
END SUBROUTINE wf_cross

!     ****** COMPLEX VALUE FIELD AT ELEMENT(NE),POINT(R,Z) [R]******

SUBROUTINE wf_fieldcr(NE,R,Z,CVALUE,CE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: ISD,M,N,nseg
  real(rkind),intent(in) :: R,Z
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),WEIGHT,L
  complex(rkind),intent(in) :: CVALUE(nseg_max)
  complex(rkind):: CF
  complex(rkind),intent(out) :: CE

  CALL wf_set_abc(NE,A,B,C)
  do ISD=1,3
     nseg=ABS(nseg_nside_nelm(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L=LSID(nseg)
     ELSE
        IF(nseg_nside_nelm(ISD,NE).GT.0.D0) THEN
           L=LSID(nseg)
        ELSE
           L=-LSID(nseg)
        END IF
     END IF

     M=ISD
     N=ISD+1
     IF(N.gt.3) N=N-3
     AW(ISD)=L*(A(M)*B(N)-A(N)*B(M))
     BW(ISD)=L*(B(M)*C(N)-B(N)*C(M))
  end do
  CE=(0.d0,0.d0)

  DO ISD=1,3
     nseg=nseg_nside_nelm(ISD,NE)
     if(nseg.lt.0) then
        nseg=-nseg
        CF=-CVALUE(nseg)
     else
        CF=CVALUE(nseg)
     end if
     WEIGHT=AW(ISD)-BW(ISD)*Z
     CE=CE+WEIGHT*CF
     IF(idebug.EQ.-1) &
          WRITE(6,'(A,I10,I5,1P5E12.4)') 'FR:',NE,ISD,weight,CF,CE
  END DO

  RETURN
END SUBROUTINE wf_fieldcr

!     ****** COMPLEX VALUE FIELD AT ELEMENT(NE),POINT(R,Z) [Z]******

SUBROUTINE wf_fieldcz(NE,R,Z,CVALUE,CE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: ISD,M,N,nseg
  real(rkind),intent(in) :: R,Z
  real(rkind) :: A(3),B(3),C(3),BW(3),CW(3),WEIGHT,L
  complex(rkind),intent(in) :: CVALUE(nseg_max)
  complex(rkind):: CF
  complex(rkind),intent(out) :: CE

  CALL wf_set_abc(NE,A,B,C)
  do ISD=1,3
     nseg=ABS(nseg_nside_nelm(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L=LSID(nseg)
     ELSE
        IF(nseg_nside_nelm(ISD,NE).GT.0.D0) THEN
           L=LSID(nseg)
        ELSE
           L=-LSID(nseg)
        END IF
     END IF
     
     M=ISD
     N=ISD+1
     IF(N.gt.3) N=N-3
     BW(ISD)=L*(B(M)*C(N)-B(N)*C(M))
     CW(ISD)=L*(C(M)*A(N)-C(N)*A(M))
  end do
  CE=(0.d0,0.d0)

  DO ISD=1,3
     nseg=nseg_nside_nelm(ISD,NE)
     if(nseg.lt.0) then
        nseg=-nseg
        CF=-CVALUE(nseg)
     else
        CF=CVALUE(nseg)
     end if
     WEIGHT=-CW(ISD)+BW(ISD)*R
     CE=CE+WEIGHT*CF
  END DO

  RETURN
END SUBROUTINE wf_fieldcz

!     ****** COMPLEX VALUE FIELD AT ELEMENT(NE),POINT(R,Z) [PHI] ******

SUBROUTINE wf_fieldcp(NE,R,Z,CVALUE,CE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: NN,IN
  real(rkind),intent(in) :: R,Z
  real(rkind) :: A(3),B(3),C(3),WEIGHT(3)
  complex(rkind),intent(in) :: CVALUE(node_max)
  complex(rkind):: CF
  complex(rkind),intent(out) :: CE

  CALL wf_set_abc(NE,A,B,C)
  CALL wf_set_weight(NE,R,Z,WEIGHT)
  CE=(0.d0,0.d0)

  DO IN=1,3
     NN=node_nside_nelm(IN,NE)
     CF=CVALUE(NN)
     CE=CE+WEIGHT(IN)*CF
  END DO

  RETURN
END SUBROUTINE wf_fieldcp

END MODULE wfsub
