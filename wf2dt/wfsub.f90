!     $Id: wfsub.f90,v 1.16 2011/11/16 09:18:22 maruyama Exp $

!     ****** SETUP NODE RANGE ******

SUBROUTINE WFSLIM

  use wfcomm
  implicit none

  integer :: IN
  real(rkind) :: LNODE

  RNDMIN=RNODE(1)
  RNDMAX=RNODE(1)
  ZNDMIN=ZNODE(1)
  ZNDMAX=ZNODE(1)
  SELECT CASE(MODELG)
  CASE(0,1,12)
     LNODE=SQRT(RNODE(1)**2+ZNODE(1)**2)
  CASE(2)
     LNODE=SQRT((RNODE(1)-RR)**2+ZNODE(1)**2)
  END SELECT
  LNDMIN=LNODE
  LNDMAX=LNODE
     
  DO IN=2,NNMAX
     RNDMIN=MIN(RNDMIN,RNODE(IN))
     RNDMAX=MAX(RNDMAX,RNODE(IN))
     ZNDMIN=MIN(ZNDMIN,ZNODE(IN))
     ZNDMAX=MAX(ZNDMAX,ZNODE(IN))
     SELECT CASE(MODELG)
     CASE(0,1,11,12)
        LNODE=SQRT(RNODE(IN)**2+ZNODE(IN)**2)
     CASE(2)
        LNODE=SQRT((RNODE(IN)-RR)**2+ZNODE(IN)**2)
     END SELECT
     LNDMIN=MIN(LNDMIN,LNODE)
     LNDMAX=MAX(LNDMAX,LNODE)
  ENDDO

  RETURN
END SUBROUTINE WFSLIM

!     ****** SETUP ELEMENT AREA ******

SUBROUTINE WFSELM

  use wfcomm
  implicit none

  integer :: NE,ncount
  real(rkind) :: RE(3),ZE(3),S

  WRITE(6,*) 'ndelm:',ndelm(1,1),ndelm(2,1),ndelm(3,1)
  ncount=0
  do NE=1,NEMAX

     call WFNODE(NE,RE,ZE)
     S=RE(1)*(ZE(2)-ZE(3))+RE(2)*(ZE(3)-ZE(1))+RE(3)*(ZE(1)-ZE(2))

     if(S.LE.0.D0) THEN
        if(nrank.eq.0) then
           write(6,'(A,I5)') 'NEGATIVE S: NE=',NE
           write(6,'(A,I6,2ES12.4)') 'NODE 1 = ',NDELM(1,NE),RE(1),ZE(1)
           write(6,'(A,I6,2ES12.4)') 'NODE 2 = ',NDELM(2,NE),RE(2),ZE(2)
           write(6,'(A,I6,2ES12.4)') 'NODE 3 = ',NDELM(3,NE),RE(3),ZE(3)
           
           ncount=ncount+1
           IF(ncount.GE.10) STOP
        end if
     end if
     
     SELM(NE)=0.5D0*S

  end do

  RETURN
END SUBROUTINE WFSELM

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

  CALL WFWGT(NE,R,Z,WGT)

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

SUBROUTINE WFWGT(NE,R,Z,WGT)
  
  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: IN
  real(rkind),intent(in) :: R,Z
  real(rkind),intent(out):: WGT(3)
  real(rkind) :: A(3),B(3),C(3)
  
  call WFABC(NE,A,B,C)

  DO IN=1,3
     WGT(IN)=A(IN)+B(IN)*R+C(IN)*Z
  ENDDO

  RETURN
END SUBROUTINE WFWGT

!     ******* A,B,C CALCULATION *******
!     A,B,C ARE USED BY AREA COORDINATE 

SUBROUTINE WFABC(NE,A,B,C)

  use wfcomm
  implicit none
  integer,intent(in):: NE
  integer :: I,J,K
  real(rkind),intent(out)::A(3),B(3),C(3) 
  real(rkind) :: RE(3),ZE(3),S
  
  CALL WFNODE(NE,RE,ZE)
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
END SUBROUTINE WFABC

!     ******* TOTAL COORDINATE - LOCAL COORDINATE *******

SUBROUTINE WFNODE(NE,RE,ZE)
  
  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: IN,NN
  real(rkind),intent(out):: RE(3),ZE(3)
  
  DO IN=1,3
     NN=NDELM(IN,NE)
     RE(IN)=RNODE(NN)
     ZE(IN)=ZNODE(NN)
  END DO

  RETURN
END SUBROUTINE WFNODE

!     *******  INITIALISE SURFACE INTEGRAL OF ELEMENT FUNCTION *******

SUBROUTINE SETAIF

  use wfcomm
  implicit none
  integer :: ID(3,3),I,L1,L2,L3,J,K
  real(rkind) :: AIF

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
END SUBROUTINE SETAIF

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

SUBROUTINE SETAIE

  use wfcomm
  implicit none
  integer :: ID(3,3),I,L1,L2,L3,J,K
  real(rkind) :: AIE

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
END SUBROUTINE SETAIE

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

!     ****** Set Boundary Attribute for Side and Node ******

SUBROUTINE SETBDY(IERR)

  use wfcomm
  implicit none
  integer,intent(out) :: IERR
  integer :: ISD,NSD,NE
  integer :: NN1,NN2,NN3,NSD1,NSD2,NSD3
  integer :: NEL,NSDL

! ----- SET NSDMAX & KNELM -----
! NBSID :: Number of Boundary Side
! KNELM(ISD,NE) :: THE ELEMENT WHICH SHARE THE SIDE (ISD,NE)
!                  IF KNELM.EQ.0, SIDE IS BOUNDARY

  if(nrank.eq.0) WRITE(6,*) '------- SETBDY set NSDMAX & KNELM start ---'

!  WRITE(6,'(A,2I8)') 'NEMAX,NNMAX=',NEMAX,NNMAX
  NBSID=0
  DO NE=1,NEMAX
     NN1=NDELM(1,NE)
     NN2=NDELM(2,NE)
     NN3=NDELM(3,NE)
     CALL EFINDS(NE,NN1,NN2,KNELM(1,NE))
     CALL EFINDS(NE,NN2,NN3,KNELM(2,NE))
     CALL EFINDS(NE,NN3,NN1,KNELM(3,NE))
     IF(KNELM(1,NE).EQ.0) NBSID=NBSID+1
     IF(KNELM(2,NE).EQ.0) NBSID=NBSID+1
     IF(KNELM(3,NE).EQ.0) NBSID=NBSID+1
!     IF(nrank.EQ.0.AND.MOD(NE-1,100).EQ.0) &
!          WRITE(6,'(A,7I8)') 'SETBDY NE,NNs,NEs:', &
!          NE,NN1,NN2,NN3,KNELM(1,NE),KNELM(2,NE),KNELM(3,NE)
  ENDDO
  NSDMAX=(3*NEMAX-NBSID)/2+NBSID

  call wfsid_allocate

!  ----- SET NSDELM -----
!  KANOD :: IF NODE IS ON THE BOUNDARY, 1

  if(nrank.eq.0) WRITE(6,*) '------- SETBDY set NSDELM start ---'

  NSD=0

  DO NE=1,NEMAX
     NN1=NDELM(1,NE)
     NN2=NDELM(2,NE)
     NN3=NDELM(3,NE)

! --- CREATE FOR THE FIRST TIME, SIDE IS NEW ---

     IF(KNELM(1,NE).GE.NE) THEN
        KANOD(NN1)=0
        KANOD(NN2)=0
        NSD=NSD+1
        NSDELM(1,NE)=NSD
        IF(NSD.LE.NSDMAX) THEN
           INSID(NSD)=1
           NESID(NSD)=NE
           NDSID(1,NSD)=NN1
           NDSID(2,NSD)=NN2
        ENDIF
     ENDIF
     IF(KNELM(2,NE).GE.NE) THEN
        KANOD(NN2)=0
        KANOD(NN3)=0
        NSD=NSD+1
        NSDELM(2,NE)=NSD
        IF(NSD.LE.NSDMAX) THEN
           INSID(NSD)=2
           NESID(NSD)=NE
           NDSID(1,NSD)=NN2
           NDSID(2,NSD)=NN3
        ENDIF
     ENDIF
     IF(KNELM(3,NE).GE.NE) THEN
        KANOD(NN3)=0
        KANOD(NN1)=0
        NSD=NSD+1
        NSDELM(3,NE)=NSD
        IF(NSD.LE.NSDMAX) THEN
           INSID(NSD)=3
           NESID(NSD)=NE
           NDSID(1,NSD)=NN3
           NDSID(2,NSD)=NN1
        ENDIF
     ENDIF

! --- ALREADY CREATED ---

     IF(KNELM(1,NE).LT.NE.AND.KNELM(1,NE).GT.0) THEN
        NEL=KNELM(1,NE)
        DO ISD=1,3
           NSDL=ABS(NSDELM(ISD,NEL))
           IF (NSDL.NE.0) THEN
              IF  ((NDSID(1,NSDL).EQ.NN1).AND.&
                  &(NDSID(2,NSDL).EQ.NN2)) THEN
                 NSDELM(1,NE)=NSDL
              ELSEIF((NDSID(1,NSDL).EQ.NN2).AND.&
                   & (NDSID(2,NSDL).EQ.NN1)) THEN
                 NSDELM(1,NE)=-NSDL
              END IF
           ENDIF
        ENDDO
     ENDIF
     IF(KNELM(2,NE).LT.NE.AND.KNELM(2,NE).GT.0) THEN
        NEL=KNELM(2,NE)
        DO ISD=1,3
           NSDL=ABS(NSDELM(ISD,NEL))
           IF (NSDL.NE.0) THEN
              IF  ((NDSID(1,NSDL).EQ.NN2).AND.&
                  &(NDSID(2,NSDL).EQ.NN3)) THEN
                 NSDELM(2,NE)=NSDL
              ELSEIF((NDSID(1,NSDL).EQ.NN3).AND.&
                   & (NDSID(2,NSDL).EQ.NN2)) THEN
                 NSDELM(2,NE)=-NSDL
              ENDIF
           ENDIF
        ENDDO
     ENDIF
     IF(KNELM(3,NE).LT.NE.AND.KNELM(3,NE).GT.0) THEN
        NEL=KNELM(3,NE)
        DO ISD=1,3
           NSDL=ABS(NSDELM(ISD,NEL))
           IF (NSDL.NE.0) THEN
              IF  ((NDSID(1,NSDL).EQ.NN3).AND.&
                  &(NDSID(2,NSDL).EQ.NN1)) THEN
                 NSDELM(3,NE)=NSDL
              ELSEIF((NDSID(1,NSDL).EQ.NN1).AND.&
                   & (NDSID(2,NSDL).EQ.NN3)) THEN
                 NSDELM(3,NE)=-NSDL
              ENDIF
           END IF
        ENDDO
     ENDIF

! --- SIDE IS ON THE BOUNDARY, SIDE IS NEW ---

     IF(KNELM(1,NE).EQ.0) THEN
        KANOD(NN1)=1
        KANOD(NN2)=1
        NSD=NSD+1
        NSDELM(1,NE)=NSD
        IF(NSD.LE.NSDMAX) THEN
           INSID(NSD)=1
           NESID(NSD)=NE
           NDSID(1,NSD)=NN1
           NDSID(2,NSD)=NN2
        ENDIF
     ENDIF
     IF(KNELM(2,NE).EQ.0) THEN
        KANOD(NN2)=1
        KANOD(NN3)=1
        NSD=NSD+1
        NSDELM(2,NE)=NSD
        IF(NSD.LE.NSDMAX) THEN
           INSID(NSD)=2
           NESID(NSD)=NE
           NDSID(1,NSD)=NN2
           NDSID(2,NSD)=NN3
        ENDIF
     ENDIF
     IF(KNELM(3,NE).EQ.0) THEN
        KANOD(NN3)=1
        KANOD(NN1)=1
        NSD=NSD+1
        NSDELM(3,NE)=NSD
        IF(NSD.LE.NSDMAX) THEN
           INSID(NSD)=3
           NESID(NSD)=NE
           NDSID(1,NSD)=NN3
           NDSID(2,NSD)=NN1
        ENDIF
     ENDIF
  ENDDO

  
!  DO NSD=1,NSDMAX
!     write(*,*) "NESID,INSID,NSD",NESID(NSD),INSID(NSD),NSD
!  ENDDO

!  DO NE=1,NEMAX
!     DO I=1,3
!        write(*,*) "NE,I,NSDELM",NE,I,NSDELM(I,NE)
!     end DO
!  end DO

!    DO NE=1,NEMAX
!     DO I=1,3
!        write(*,*) "NE,I,KNELM",NE,I,KNELM(I,NE)
!     end DO
!  end DO

 ! ----- SET KASID -----

  if(nrank.eq.0) WRITE(6,*) '------- SETBDY set KASID start ---'

  DO NSD=1,NSDMAX
     KASID(NSD)=0
  ENDDO

  DO NE=1,NEMAX
     NSD1=ABS(NSDELM(1,NE))
     NSD2=ABS(NSDELM(2,NE))
     NSD3=ABS(NSDELM(3,NE))
     IF(KNELM(1,NE).EQ.0) KASID(NSD1)=1
     IF(KNELM(2,NE).EQ.0) KASID(NSD2)=1
     IF(KNELM(3,NE).EQ.0) KASID(NSD3)=1
  ENDDO
  
  IERR=0
  RETURN
END SUBROUTINE SETBDY

!     ***** SET BOUNDARY ELECTRIC FIELD *****

SUBROUTINE SETEWG

  USE wfcomm
  USE wfload,ONLY: wf_read_wg 
  implicit none
  INTEGER:: NSD,NN1,NN2,NN,NBSD,NBND,IERR
  REAL(rkind):: ANGLE,R,Z,PHASE,PROD,FACTOR,SN
  COMPLEX(rkind):: CEX,CEY,CEZ
  REAL(rkind),PARAMETER:: EPSWG=1.D-12

  ANGLE=ANGWG*PI/180.D0

! --- WG Electric field on boundary side ---  

  NBSID=0
  DO NSD=1,NSDMAX
     IF(KASID(NSD).EQ.1) NBSID=NBSID+1
  END DO
  IF(nrank.EQ.0) write(6,*) '## setewg: NBSID=',NBSID
  IF(ALLOCATED(NSDBS)) DEALLOCATE(NSDBS)
  IF(ALLOCATED(CEBSD)) DEALLOCATE(CEBSD)
  ALLOCATE(NSDBS(NBSID))
  ALLOCATE(CEBSD(NBSID))
  NBSID=0
  DO NSD=1,NSDMAX
     IF(KASID(NSD).EQ.1) THEN
        NBSID=NBSID+1
        NSDBS(NBSID)=NSD
        KBSID(NSD)=NBSID
     ELSE
        KBSID(NSD)=0
     END IF
  END DO
  DO NBSD=1,NBSID
     NSD=NSDBS(NBSD)
     NN1=NDSID(1,NSD)
     NN2=NDSID(2,NSD)
     R=0.5D0*(RNODE(NN1)+RNODE(NN2))
     Z=0.5D0*(ZNODE(NN1)+ZNODE(NN2))
     SELECT CASE(MODELWG)
     CASE(0,1)
!        WRITE(29,'(A,6ES12.4)') 'R,Z:WG:',R,Z,R1WG,R2WG,Z1WG,Z2WG
        IF((R.GE.R1WG).AND.(R.LE.R2WG).AND. &
             (Z.GE.Z1WG).AND.(Z.LE.Z2WG)) THEN
           WRITE(6,'(A,2ES12.4)') 'R,Z in WG:',R,Z
           PROD=(R2WG-R1WG)*(RNODE(NN2)-RNODE(NN1)) &
               +(Z2WG-Z1WG)*(ZNODE(NN2)-ZNODE(NN1))
           IF(ABS(R1WG-R2WG).LT.1.D-8) THEN
              FACTOR=(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
           ELSE IF(ABS(Z1WG-Z2WG).LT.1.D-8) THEN
              FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2
           ELSE
              FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2 &
                    +(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
           END IF
           SN=SQRT((R   -R1WG)**2+(Z   -Z1WG)**2) &
             /SQRT((R2WG-R1WG)**2+(Z2WG-Z1WG)**2) ! SN=0 at 1, 1 at 2
           PHASE=(PH1WG+(PH2WG-PH1WG)*SN+DPHWG*4.D0*SN*(1.D0-SN))*PI/180.D0
           CEBSD(NBSD)= AMPWG*EXP(CII*PHASE)*(COS(ANGLE)+CII*ELPWG*SIN(ANGLE))
           IF(MODELWG.EQ.1) CEBSD(NBSD)=CEBSD(NBSD)*EXP(-10.D0*FACTOR)
           IF(PROD.GT.0.D0) CEBSD(NBSD)=-CEBSD(NBSD)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                WRITE(6,'(A,2I8,1P5E12.4)') &
                'SD:',NSD,NBSD,CEBSD(NBSD),AMPWG,PHASE,ANGLE
        ELSE
           CEBSD(NBSD)=(0.D0,0.D0)
        END IF
     CASE(12)
        IF((R.GE.R1WG-EPSWG).AND.(R.LE.R2WG+EPSWG).AND. &
           (Z.GE.Z1WG-EPSWG).AND.(Z.LE.Z2WG+EPSWG)) THEN
           PROD=(R2WG-R1WG)*(RNODE(NN2)-RNODE(NN1)) &
               +(Z2WG-Z1WG)*(ZNODE(NN2)-ZNODE(NN1))
           CALL wf_read_wg(Z,CEX,CEY,CEZ,IERR)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                write(6,'(A,1P6E12.4)') 'R,Z,CEY=', &
                                    R,Z,CEY,PROD,ZNODE(NN2)-ZNODE(NN1)
!!!           IF(PROD.GT.0.D0) CEY=-CEY
           CEBSD(NBSD)=AMPWG*CEY
        ELSE
           CEBSD(NBSD)=(0.D0,0.D0)
        END IF
     END SELECT
  END DO

! --- WG Electric field on boundary node ---  

  NBNOD=0
  DO NN=1,NNMAX
     IF(KANOD(NN).EQ.1) NBNOD=NBNOD+1
  END DO
  IF(nrank.EQ.0) write(6,*) '## setewg: NBNOD=',NBNOD
  IF(ALLOCATED(NNDBS)) DEALLOCATE(NNDBS)
  IF(ALLOCATED(CEBND)) DEALLOCATE(CEBND)
  ALLOCATE(NNDBS(NBNOD))
  ALLOCATE(CEBND(NBNOD))
  NBNOD=0
  DO NN=1,NNMAX
     IF(KANOD(NN).EQ.1) THEN
        NBNOD=NBNOD+1
        NNDBS(NBNOD)=NN
        KBNOD(NN)=NBNOD
     ELSE
        KBNOD(NN)=0
     END IF
  END DO
  DO NBND=1,NBNOD
     NN=NNDBS(NBND)
     R=RNODE(NN)
     Z=ZNODE(NN)
     SELECT CASE(MODELWG)
     CASE(0,1)
        IF((R.GE.R1WG).AND.(R.LE.R2WG).AND. &
           (Z.GE.Z1WG).AND.(Z.LE.Z2WG)) THEN
           IF(ABS(R1WG-R2WG).LT.1.D-8) THEN
              FACTOR=(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
           ELSE IF(ABS(Z1WG-Z2WG).LT.1.D-8) THEN
              FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2
           ELSE
              FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2 &
                    +(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
           END IF
           SN=SQRT((R   -R1WG)**2+(Z   -Z1WG)**2) &
                /SQRT((R2WG-R1WG)**2+(Z2WG-Z1WG)**2) ! SN=0 at 1, 1 at 2
           PHASE=(PH1WG+(PH2WG-PH1WG)*SN+DPHWG*4.D0*SN*(1.D0-SN))*PI/180.D0
           CEBND(NBND)= AMPWG*EXP(CII*PHASE)*(SIN(ANGLE)-CII*ELPWG*COS(ANGLE))
           IF(MODELWG.EQ.1) CEBND(NBND)=CEBND(NBND)*EXP(-10.D0*FACTOR)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                WRITE(6,'(A,2I8,1P5E12.4)') &
                'ND:',NN,NBND,CEBND(NBND),AMPWG,PHASE,ANGLE
        ELSE
           CEBND(NBND)=(0.D0,0.D0)
        END IF
     CASE(12)
        IF((R.GE.R1WG-EPSWG).AND.(R.LE.R2WG+EPSWG).AND. &
           (Z.GE.Z1WG-EPSWG).AND.(Z.LE.Z2WG+EPSWG)) THEN
           CALL wf_read_wg(Z,CEX,CEY,CEZ,IERR)
           IF(nrank.EQ.0.AND.idebug.EQ.3) &
                write(6,'(A,1P4E12.4)') 'R,Z,CEZ=',R,Z,CEZ
           CEBND(NBND)=AMPWG*CEZ
        ELSE
           CEBND(NBND)=(0.D0,0.D0)
        END IF
     END SELECT
     END DO
  RETURN
END SUBROUTINE SETEWG

!     ***** SET LENGTH OF SIDE *****

SUBROUTINE SETLSD
  use wfcomm
  implicit none

  integer :: NSD,ND1,ND2
  real(rkind) :: R1,R2,Z1,Z2,L

  do NSD=1,NSDMAX
     LSID(NSD)=0.d0
  end do

  do NSD=1,NSDMAX
     ND1=NDSID(1,NSD)
     ND2=NDSID(2,NSD)
     R1 =RNODE(ND1)
     R2 =RNODE(ND2)
     Z1 =ZNODE(ND1)
     Z2 =ZNODE(ND2)

     L =SQRT((R2-R1)**2+(Z2-Z1)**2)
     LSID(NSD)=L
  end do
  return
end SUBROUTINE SETLSD

!     ****** MODIFY ANTENNA DATA ******

SUBROUTINE MODANT(IERR)
  
  use wfcomm
  implicit none

  integer,intent(out) :: IERR
  integer :: NE,NA,NSD,L,KN,LS,N,ID,NENEXT,NENEW
  real(rkind) :: RC,ZC

  NE=0
  DO NA=1,NAMAX
     CALL FEP(RJ0(1,NA),ZJ0(1,NA),NE)
     IF(nrank.EQ.0) &
          WRITE(6,'(A,I5,1P2E12.4,I5)') 'NA,RJ0,ZJ0=',NA,RJ0(1,NA),ZJ0(1,NA),NE
!    outside starting point

     IF(NE.EQ.0) THEN
        IF(JNUM0(NA).EQ.1) GOTO 8500 ! error: one point and outside
        DO NSD=1,NSDMAX              ! look for crossing boundary
           L =INSID(NSD)
           NE=NESID(NSD)
           KN=KNELM(L,NE)

           IF(KN.eq.0) THEN 
              CALL CROS(RJ0(1,NA),ZJ0(1,NA),&
                   &    RJ0(2,NA),ZJ0(2,NA),&
                   &    NE,L,RC,ZC,IERR)
              IF(IERR.EQ.0) THEN
                 LS=L
                 GOTO 1000   ! crossing point on boundary found
              ENDIF
           ENDIF

        ENDDO
        GOTO 8000 ! error: no crossing boundary found
        
1000    CONTINUE  ! set starting point (N=1)
        N=1
        RJ(N,NA)=RC
        ZJ(N,NA)=ZC
        JELMT(N,NA)=NE
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,R,Z=' ,&
                &NA,N,NE,RJ(N,NA),ZJ(N,NA)
        ENDIF

!    inside starting point

     ELSE
        N=1
        RJ(N,NA)=RJ0(1,NA)
        ZJ(N,NA)=ZJ0(1,NA)
        JELMT(N,NA)=NE
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,R,Z=',&
                &NA,N,NE,RJ(N,NA),ZJ(N,NA)
        ENDIF
        LS=0
     ENDIF
     
     DO ID=2,JNUM0(NA)
        NENEXT=NE
        CALL FEP(RJ0(ID,NA),ZJ0(ID,NA),NENEXT)
3000    CONTINUE
        IF(NENEXT.EQ.NE) GOTO 4500 ! next antenna node in the same element
        DO L=1,3                   ! look for crossing point
           IF(L.NE.LS) THEN
              CALL CROS(RJ (N ,NA),ZJ (N ,NA),&
                   &    RJ0(ID,NA),ZJ0(ID,NA),&
                   &         NE,L,RC,ZC,IERR)
              IF(IERR.EQ.0) THEN
                 LS=L
                 GOTO 4000 ! crossing point found
              ENDIF
           ENDIF
        ENDDO
        GOTO 8100 
        
4000    IF(N+1.GT.NJM) GOTO 8200 ! error: number of antenna nodes overflow 
        N=N+1
        RJ(N,NA)=RC
        ZJ(N,NA)=ZC
        JELMT(N,NA)=NE
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,R,Z=',&
                &NA,N,NE,RJ(N,NA),ZJ(N,NA)
        ENDIF
        NENEW=KNELM(LS,NE) ! look for neighboring element
        IF(NENEW.GT.0) THEN ! element found
           DO L=1,3
              IF(KNELM(L,NENEW).EQ.NE) LS=L
           ENDDO
           NE=NENEW         ! set side LS and element NE
        ELSE
           IF(ID.EQ.JNUM0(NA).AND.NENEXT.LE.0) GOTO 6000 ! last point outside
           GOTO 8400 ! error: cross point to outside
        ENDIF
        GOTO 3000 ! look for next cross point 
        
4500    IF(N+1.GT.NJM) GOTO 8200 ! error: number of antenna nodes overflow 
        N=N+1
        RJ(N,NA)=RJ0(ID,NA)
        ZJ(N,NA)=ZJ0(ID,NA)
        JELMT(N,NA)=NE
        LS=0
        IF(IDEBUG.EQ.2) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,R,Z=',&
                & NA,N,NE,RJ(N,NA),ZJ(N,NA)
        ENDIF
     ENDDO
     
6000 JNUM(NA)=N
  ENDDO
  IERR=0
  RETURN
  
8000 IERR=8000
  JNUM(NA)=N
  if(nrank.eq.0) WRITE(6,800) IERR,NA,N
800 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         & ' ','                : CANNOT FIND BOUNDARY POINT'/&
         & ' ','                : NA,N=',2I7)
  JNUM(NA)=N
  RETURN
  
8100 IERR=8100
  if(nrank.eq.0) WRITE(6,810) IERR,NA,ID,N
810 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         & ' ','                : CANNOT FIND NEXT CROSSPOINT'/&
         & ' ','                : NA,ID,N =',3I7)
  JNUM(NA)=N
  RETURN
  
8200 if(nrank.eq.0) WRITE(6,820) NA,ID,N,NJM
820 FORMAT(' ','## MODANT ERROR : N.GT.NJM '/&
         & ' ','                : NA,ID,N,NJM = ',4I7)
  IERR=8200
  JNUM(NA)=N
  RETURN
  
8400 IERR=8400
  if(nrank.eq.0) WRITE(6,840) IERR,NA,ID,NE,N
840 FORMAT(' ','## MODANT ERROR : IERR =',I5/&
         & ' ','                : ABMORMAL END OF ANTENNA DATA '/&
         & ' ','                : NA,ID,NE,N = ',4I7)
  JNUM(NA)=N
  IERR=8400
  RETURN
  
8500 IERR=8500
  if(nrank.eq.0) WRITE(6,850) IERR,NA,NE,JNUM0(NA)
850 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         & ' ','                : NA,NE,JNUM0 = ',3I7)
  RETURN
  
END SUBROUTINE MODANT

!     ******* CALCULATE POINT OF INTERSECTION *******

SUBROUTINE CROS(R1,Z1,R2,Z2,IE,L,RC,ZC,IERR)
  
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
  N1=NDELM(L,IE)
  N2=NDELM(M,IE)
  R3=RNODE(N1)
  Z3=ZNODE(N1)
  R4=RNODE(N2)
  Z4=ZNODE(N2)

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
END SUBROUTINE CROS

!     ****** COMPLEX VALUE FIELD AT ELEMENT(NE),POINT(R,Z) [R]******

SUBROUTINE FIELDCR(NE,R,Z,CVALUE,CE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: ISD,M,N,NSD
  real(rkind),intent(in) :: R,Z
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),WEIGHT,L
  complex(rkind),intent(in) :: CVALUE(NSDMAX)
  complex(rkind):: CF
  complex(rkind),intent(out) :: CE

  CALL WFABC(NE,A,B,C)
  do ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L=LSID(NSD)
     ELSE
        IF(NSDELM(ISD,NE).GT.0.D0) THEN
           L=LSID(NSD)
        ELSE
           L=-LSID(NSD)
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
     NSD=NSDELM(ISD,NE)
     if(NSD.lt.0) then
        NSD=-NSD
        CF=-CVALUE(NSD)
     else
        CF=CVALUE(NSD)
     end if
     WEIGHT=AW(ISD)-BW(ISD)*Z
     CE=CE+WEIGHT*CF
     IF(idebug.EQ.-1) &
          WRITE(6,'(A,I10,I5,1P5E12.4)') 'FR:',NE,ISD,weight,CF,CE
  END DO

  RETURN
END SUBROUTINE FIELDCR

!     ****** COMPLEX VALUE FIELD AT ELEMENT(NE),POINT(R,Z) [Z]******

SUBROUTINE FIELDCZ(NE,R,Z,CVALUE,CE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: ISD,M,N,NSD
  real(rkind),intent(in) :: R,Z
  real(rkind) :: A(3),B(3),C(3),BW(3),CW(3),WEIGHT,L
  complex(rkind),intent(in) :: CVALUE(NSDMAX)
  complex(rkind):: CF
  complex(rkind),intent(out) :: CE

  CALL WFABC(NE,A,B,C)
  do ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L=LSID(NSD)
     ELSE
        IF(NSDELM(ISD,NE).GT.0.D0) THEN
           L=LSID(NSD)
        ELSE
           L=-LSID(NSD)
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
     NSD=NSDELM(ISD,NE)
     if(NSD.lt.0) then
        NSD=-NSD
        CF=-CVALUE(NSD)
     else
        CF=CVALUE(NSD)
     end if
     WEIGHT=-CW(ISD)+BW(ISD)*R
     CE=CE+WEIGHT*CF
  END DO

  RETURN
END SUBROUTINE FIELDCZ

!     ****** COMPLEX VALUE FIELD AT ELEMENT(NE),POINT(R,Z) [PHI] ******

SUBROUTINE FIELDCP(NE,R,Z,CVALUE,CE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: NN,IN
  real(rkind),intent(in) :: R,Z
  real(rkind) :: A(3),B(3),C(3),WEIGHT(3)
  complex(rkind),intent(in) :: CVALUE(NNMAX)
  complex(rkind):: CF
  complex(rkind),intent(out) :: CE

  CALL WFABC(NE,A,B,C)
  CALL WFWGT(NE,R,Z,WEIGHT)
  CE=(0.d0,0.d0)

  DO IN=1,3
     NN=NDELM(IN,NE)
     CF=CVALUE(NN)
     CE=CE+WEIGHT(IN)*CF
  END DO

  RETURN
END SUBROUTINE FIELDCP
