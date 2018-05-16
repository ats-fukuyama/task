!     $Id$

!     ****** SETUP NODE RANGE ******

SUBROUTINE WFVLIM

  use wfcomm
  implicit none

  integer :: IN
  real(8) :: RND

  XNDMIN=XND(1)
  XNDMAX=XND(1)
  YNDMIN=YND(1)
  YNDMAX=YND(1)
  ZNDMIN=ZND(1)
  ZNDMAX=ZND(1)
  RND=SQRT(XND(1)**2+YND(1)**2)
  RNDMIN=RND
  RNDMAX=RND
  
  DO IN=2,NNMAX
     XNDMIN=MIN(XNDMIN,XND(IN))
     XNDMAX=MAX(XNDMAX,XND(IN))
     YNDMIN=MIN(YNDMIN,YND(IN))
     YNDMAX=MAX(YNDMAX,YND(IN))
     ZNDMIN=MIN(ZNDMIN,ZND(IN))
     ZNDMAX=MAX(ZNDMAX,ZND(IN))
     RND   =SQRT(XND(IN)**2+YND(IN)**2)
     RNDMIN=MIN(RNDMIN,RND)
     RNDMAX=MAX(RNDMAX,RND)
  ENDDO

  RETURN
END SUBROUTINE WFVLIM

!     ****** SETUP ELEMENT VOLUME ******

SUBROUTINE WFVELM

  use wfcomm
  implicit none

  integer :: NN,NEDO,NE,N
  real(8) :: A(4),B(4),C(4),D(4),V

  VTOT=0.D0
  DO NN=1,NNMAX
     VNOD(NN)=0.D0
  ENDDO
  DO NEDO=1,NEMAX
     NE=NEDO
     CALL WFABCD(NE,A,B,C,D,V)
     VELM(NE)=V
     VTOT    =VTOT    +V
     DO N=1,4
        NN=NDELM(N,NE)
        VNOD(NN)=VNOD(NN)+V*AIG1(N)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE WFVELM

!     ****** Set Bounary Attribute for Surface and Node ******

SUBROUTINE SETBDY(IERR)

  use wfcomm
  implicit none

  integer :: NSF,NEDO,NE,NN1,NN2,NN3,NN4,I,IERR
!  DATA EPS/0.1D0/

  NSF=0

  DO NEDO=1,NEMAX
     NE=NEDO
     NN1=NDELM(1,NE)
     NN2=NDELM(2,NE)
     NN3=NDELM(3,NE)
     NN4=NDELM(4,NE)
     CALL EFINDK(NE,NN2,NN3,NN4,KNELM(1,NE))
     CALL EFINDK(NE,NN3,NN4,NN1,KNELM(2,NE))
     CALL EFINDK(NE,NN4,NN1,NN2,KNELM(3,NE))
     CALL EFINDK(NE,NN1,NN2,NN3,KNELM(4,NE))
     
     IF(KNELM(1,NE).EQ.0) NSF=NSF+1
     IF(KNELM(2,NE).EQ.0) NSF=NSF+1
     IF(KNELM(3,NE).EQ.0) NSF=NSF+1
     IF(KNELM(4,NE).EQ.0) NSF=NSF+1
  ENDDO
  NSFMAX=NSF

  CALL wfsrf_allocate

  NSF=0

  if(nrank.eq.0) WRITE(6,*) '------- SETBDY set KNELM start ---'
  
  DO NEDO=1,NEMAX
     NE=NEDO
     NN1=NDELM(1,NE)
     NN2=NDELM(2,NE)
     NN3=NDELM(3,NE)
     NN4=NDELM(4,NE)
     CALL EFINDK(NE,NN2,NN3,NN4,KNELM(1,NE))
     CALL EFINDK(NE,NN3,NN4,NN1,KNELM(2,NE))
     CALL EFINDK(NE,NN4,NN1,NN2,KNELM(3,NE))
     CALL EFINDK(NE,NN1,NN2,NN3,KNELM(4,NE))
     
     !         WRITE(6,'(A,6I8)') 'NE,KNELM=',NE,(KNELM(I,NE),I=1,4),NSF
     !         IF(NE.EQ.100) STOP
     
     IF(KNELM(1,NE).EQ.0) THEN
        KANOD(NN2)=1
        KANOD(NN3)=1
        KANOD(NN4)=1
        NSF=NSF+1
        IF(NSF.LE.NSFMAX) THEN
           INSRF(NSF)=1
           NESRF(NSF)=NE
           NDSRF(1,NSF)=NN2
           NDSRF(2,NSF)=NN4
           NDSRF(3,NSF)=NN3
        ENDIF
     ENDIF
     IF(KNELM(2,NE).EQ.0) THEN
        KANOD(NN3)=1
        KANOD(NN4)=1
        KANOD(NN1)=1
        NSF=NSF+1
        IF(NSF.LE.NSFMAX) THEN
           INSRF(NSF)=2
           NESRF(NSF)=NE
           NDSRF(1,NSF)=NN3
           NDSRF(2,NSF)=NN4
           NDSRF(3,NSF)=NN1
        ENDIF
     ENDIF
     IF(KNELM(3,NE).EQ.0) THEN
        KANOD(NN4)=1
        KANOD(NN1)=1
        KANOD(NN2)=1
        NSF=NSF+1
        IF(NSF.LE.NSFMAX) THEN
           INSRF(NSF)=3
           NESRF(NSF)=NE
           NDSRF(1,NSF)=NN4
           NDSRF(2,NSF)=NN2
           NDSRF(3,NSF)=NN1
        ENDIF
     ENDIF
     IF(KNELM(4,NE).EQ.0) THEN
        KANOD(NN1)=1
        KANOD(NN2)=1
        KANOD(NN3)=1
        NSF=NSF+1
        IF(NSF.LE.NSFMAX) THEN
           INSRF(NSF)=4
           NESRF(NSF)=NE
           NDSRF(1,NSF)=NN1
           NDSRF(2,NSF)=NN2
           NDSRF(3,NSF)=NN3
        ENDIF
     ENDIF
  ENDDO
  
  DO NE=1,NEMAX
     DO I=1,4
        NSFELM(I,NE)=0
     ENDDO
  ENDDO
  DO NSF=1,NSFMAX
     NSFELM(INSRF(NSF),NESRF(NSF))=NSF
  ENDDO
  
  if(nrank.eq.0) WRITE(6,*) '------- SETBDY set KNSRF start ---'
  
  DO NSF=1,NSFMAX
     NN1=NDSRF(1,NSF)
     NN2=NDSRF(2,NSF)
     NN3=NDSRF(3,NSF)
     CALL EFINDS(NSF,NN2,NN3,KNSRF(1,NSF))
     CALL EFINDS(NSF,NN3,NN1,KNSRF(2,NSF))
     CALL EFINDS(NSF,NN1,NN2,KNSRF(3,NSF))
     IF(KNSRF(1,NSF).EQ.0) THEN
        if(nrank.eq.0) WRITE(6,*) 'XX SRF INCONSISTENCY: 1,NSF=',NSF
     ENDIF
     IF(KNSRF(2,NSF).EQ.0) THEN
        if(nrank.eq.0) WRITE(6,*) 'XX SRF INCONSISTENCY: 2,NSF=',NSF
     ENDIF
     IF(KNSRF(3,NSF).EQ.0) THEN
        if(nrank.eq.0) WRITE(6,*) 'XX SRF INCONSISTENCY: 3,NSF=',NSF
     ENDIF
!         WRITE(6,'(A,7I8)') 'NSF,N,K=',NSF,NDSRF(1,NSF),NDSRF(2,NSF),&
!     &                                     NDSRF(3,NSF),KNSRF(1,NSF),&
!     &                                     KNSRF(2,NSF),KNSRF(3,NSF)
  ENDDO

!      DO NN=1,NNMAX
!         IF(KANOD(NN).EQ.1) THEN
!            WRITE(6,'(A,I5,1P3E12.4)') &
!     &           'KA=1:',NN,XND(NN),YND(NN),ZND(NN)
!         ENDIF
!      ENDDO

!      DO NN=1,NNMAX
!         IF(KANOD(NN).EQ.2) THEN
!            WRITE(6,'(A,I5,1P3E12.4)') &
!     &           'KA=2:',NN,XND(NN),YND(NN),ZND(NN)
!         ENDIF
!      ENDDO

!      DO NE=1,NEMAX
!         WRITE(6,'(A,I5,3X,3I5,3X,3I5)') 'NE=',NE,&
!     &        NDELM(1,NE),NDELM(2,NE),NDELM(3,NE),&
!     &        KNELM(1,NE),KNELM(2,NE),KNELM(3,NE)
!      ENDDO


  IERR=0
  RETURN
END SUBROUTINE SETBDY

!     ****** SETUP SIDE ARRAYS ******

SUBROUTINE SETSID(IERR)

  use wfcomm
  implicit none

  integer :: NE,I,NSD,NN1,NN2,NEL,IN,ICOUNT,KN,INLNN,ISD
  integer :: NSDL,INL,NN,NSF,NN3,ND1,ND2,IERR,INSID(4,6)

  DATA INSID /1,2,3,4, 2,3,1,4, 1,3,4,2, &
     &        1,4,2,3, 2,4,3,1, 3,4,1,2/

! decide NSDMAX

  NSDMAX=0

!!!!     --- create side array ---                       
  DO NE=1,NEMAX
     DO NSD=1,6
        NN1=NDELM(INSID(1,NSD),NE)
        NN2=NDELM(INSID(2,NSD),NE)
!!!!     --- find element including NN1, NN2 ---
        NEL=NE
        IN=INSID(3,NSD)
        ICOUNT=0

100      CONTINUE
        KN=KNELM(IN,NEL)   ! get neighboring element
        DO I=1,3
           INL=MOD(IN+I-1,4)+1
           NN=NDELM(INL,NEL)
           IF(NN.NE.NN1.AND.NN.NE.NN2) GOTO 150
        ENDDO
        GOTO 9100
        
150      IN=INL
        IF(KN.EQ.0)  GOTO 200 ! arrive at one surface. Try another
        IF(KN.EQ.NE) GOTO 400 ! could not find element. side is new
        IF(KN.LT.NE) GOTO 500 ! already created
        NEL=KN
        DO I=1,4
           IF(NDELM(I,NEL).EQ.NN) IN=I
        ENDDO
        ICOUNT=ICOUNT+1
        IF(ICOUNT.GT.100) THEN
           IF(ICOUNT.GT.110) STOP
        ENDIF
        GOTO 100

200      NEL=NE
        IN=INSID(4,NSD)
300      CONTINUE
        KN=KNELM(IN,NEL)   ! get neighboring element
        DO I=1,3
           INL=MOD(IN+I-1,4)+1
           NN=NDELM(INL,NEL)
           IF(NN.NE.NN1.AND.NN.NE.NN2) GOTO 350
        ENDDO
        GOTO 9100
        
350     IN=INL
        IF(KN.EQ.0)  GOTO 400 !  arrive at another surface. side is new
        IF(KN.EQ.NE) GOTO 400 !  could not find element. side is new
        IF(KN.LT.NE) GOTO 500 !  already created
        NEL=KN
        DO I=1,4
           IF(NDELM(I,NEL).EQ.NN) IN=I
        ENDDO
        ICOUNT=ICOUNT+1
        IF(ICOUNT.GT.100) THEN
           IF(ICOUNT.GT.110) STOP
        ENDIF
        GOTO 300

!!!!     --- the side is new ---
400     CONTINUE
        NSDMAX=NSDMAX+1
        GOTO 600

!!!!     --- the side is already defined ---
500     CONTINUE
        NEL=KN
        GOTO 600
        
600     CONTINUE
  
     end DO
  end DO
  
  call wfsid_allocate

  NSDMAX=0
  DO NE=1,NEMAX
     DO I=1,6
        NSDELM(I,NE)=0
     ENDDO
  ENDDO
  
!     --- create side array ---

  DO NE=1,NEMAX
     DO NSD=1,6
        NN1=NDELM(INSID(1,NSD),NE)
        NN2=NDELM(INSID(2,NSD),NE)

!     --- find element including NN1, NN2 ---

        NEL=NE
        IN=INSID(3,NSD)
        ICOUNT=0
10      CONTINUE
        KN=KNELM(IN,NEL)   ! get neighboring element
        DO I=1,3
           INL=MOD(IN+I-1,4)+1
           NN=NDELM(INL,NEL)
           IF(NN.NE.NN1.AND.NN.NE.NN2) GOTO 15
        ENDDO
        GOTO 9100
        
15      IN=INL
!        WRITE(6,'(A,6I8)') 'NE,NSD,NEL,IN,NN,KN=',&
!                        &   NE,NSD,NEL,IN,NN,KN
        IF(KN.EQ.0)  GOTO 20  ! arrive at one surface. Try another
        IF(KN.EQ.NE) GOTO 40 ! could not find element. side is new
        IF(KN.LT.NE) GOTO 50 ! already created
        NEL=KN
        DO I=1,4
           IF(NDELM(I,NEL).EQ.NN) IN=I
        ENDDO
        ICOUNT=ICOUNT+1
        IF(ICOUNT.GT.100) THEN
           if(nrank.eq.0) then
              WRITE(6,'(A,9I7)') 'NEL,ND,KL=', &
                   &  NEL,(NDELM(I,NEL),I=1,4),(KNELM(I,NEL),I=1,4)
              WRITE(6,'(A,5I7)') '   ID(NEL),ID(ND)=', &
                   &   IDELM(NEL),(IDND(NDELM(I,NEL)),I=1,4)
           end if
           IF(ICOUNT.GT.110) STOP
        ENDIF
        GOTO 10
        
20      NEL=NE
        IN=INSID(4,NSD)
30      CONTINUE
        KN=KNELM(IN,NEL)   ! get neighboring element
        DO I=1,3
           INL=MOD(IN+I-1,4)+1
           NN=NDELM(INL,NEL)
           IF(NN.NE.NN1.AND.NN.NE.NN2) GOTO 35
        ENDDO
        GOTO 9100
        
35      IN=INL
!        WRITE(6,'(A,6I8)') 'NE,NSD,NEL,IN,NN,KN=', &
!                       &    NE,NSD,NEL,IN,NN,KN
        IF(KN.EQ.0)  GOTO 40 !  arrive at another surface. side is new
        IF(KN.EQ.NE) GOTO 40 !  could not find element. side is new
        IF(KN.LT.NE) GOTO 50 !  already created
        NEL=KN
        DO I=1,4
           IF(NDELM(I,NEL).EQ.NN) IN=I
        ENDDO
        ICOUNT=ICOUNT+1
        IF(ICOUNT.GT.100) THEN
           if(nrank.eq.0) then
              WRITE(6,'(A,9I7)') 'NEL,ND,KL=', &
                   & NEL,(NDELM(I,NEL),I=1,4),(KNELM(I,NEL),I=1,4)
              WRITE(6,'(A,5I7)') '   ID(NEL),ID(ND)=', &
                   & IDELM(NEL),(IDND(NDELM(I,NEL)),I=1,4)
           end if
           IF(ICOUNT.GT.110) STOP
        ENDIF
        GOTO 30
            
!     --- the side is new ---

40      CONTINUE
        NSDMAX=NSDMAX+1
           NDSID(1,NSDMAX)=NN1
           NDSID(2,NSDMAX)=NN2
           NSDELM(NSD,NE)=NSDMAX
!        WRITE(6,'(A,3I8)') 'NSD,NN1,NN2=',NSDMAX,NN1,NN2
        GOTO 60

!     --- the side is already defined ---
        
50      CONTINUE
        NEL=KN
        DO ISD=1,6
           NSDL=ABS(NSDELM(ISD,NEL))
           IF((NDSID(1,NSDL).EQ.NN1).AND.&
             &(NDSID(2,NSDL).EQ.NN2)) THEN
              NSDELM(NSD,NE)=NSDL
              GOTO 60
           ELSEIF((NDSID(1,NSDL).EQ.NN2).AND.&
                & (NDSID(2,NSDL).EQ.NN1)) THEN
              NSDELM(NSD,NE)=-NSDL
              GOTO 60
           ENDIF
        ENDDO
        GOTO 9200
        
60      CONTINUE
     ENDDO
  ENDDO

!  DO NE=1,NEMAX
!     WRITE(6,'(A,7I8)') 'NE,NSD=',NE,(NSDELM(I,NE),I=1,6)
!  ENDDO

! --- SET KASID ---

  DO NSD=1,NSDMAX
     KASID(NSD)=0
  ENDDO
  DO NSF=1,NSFMAX
     NE=NESRF(NSF)
     NN1=NDSRF(1,NSF)
     NN2=NDSRF(2,NSF)
     NN3=NDSRF(3,NSF)
     DO I=1,6
        NSD=ABS(NSDELM(I,NE))
        ND1=NDSID(1,NSD)
        ND2=NDSID(2,NSD)
!        WRITE(6,'(A,7I8)') 'NSF,NN1,NN2,NN3,NSD,ND1,ND2=', &
!                        &   NSF,NN1,NN2,NN3,NSD,ND1,ND2
        IF    (ND1.EQ.NN2.AND.ND2.EQ.NN3) THEN
           NSDSRF(1,NSF)= NSD
           KASID(NSD)=1
        ELSEIF(ND2.EQ.NN2.AND.ND1.EQ.NN3) THEN
           NSDSRF(1,NSF)=-NSD
           KASID(NSD)=1
        ELSEIF(ND1.EQ.NN3.AND.ND2.EQ.NN1) THEN
           NSDSRF(2,NSF)= NSD
           KASID(NSD)=1
        ELSEIF(ND2.EQ.NN3.AND.ND1.EQ.NN1) THEN
           NSDSRF(2,NSF)=-NSD
           KASID(NSD)=1
        ELSEIF(ND1.EQ.NN1.AND.ND2.EQ.NN2) THEN
           NSDSRF(3,NSF)= NSD
           KASID(NSD)=1
        ELSEIF(ND2.EQ.NN1.AND.ND1.EQ.NN2) THEN
           NSDSRF(3,NSF)=-NSD
           KASID(NSD)=1
        ENDIF
     ENDDO
  ENDDO
!  DO NSF=1,NSFMAX
!  WRITE(6,'(A,4I8)') 'NSF,NSD1,NSD2,NSD3=', &
!                  &   NSF,NSDSRF(1,NSF),NSDSRF(2,NSF),NSDSRF(3,NSF)
!  ENDDO

  IERR=0
  RETURN
  
9100 if(nrank.eq.0) WRITE(6,*) 'XX SETSID: The third node could not be found.'
  IERR=1
  RETURN
  
9200 if(nrank.eq.0) WRITE(6,*) 'XX SETSID: NEL ERROR: NE,NEL=',NE,NEL
  IERR=1
  RETURN
  
END SUBROUTINE SETSID

!     ****** MODIFY ANTENNA DATA ******

SUBROUTINE MODANT(IERR)
  
  use wfcomm
  implicit none

  integer :: NE,NA,NSF,L,KN,IERR,LS,N,ID,NENEXT,NENEW
  real(8) :: XC,YC,ZC
  
  NE=0
  DO NA=1,NAMAX
     CALL FEP(XJ0(1,NA),YJ0(1,NA),ZJ0(1,NA),NE)
     
!        WHEN FIRST POINT IS OUT OF REGION

     IF(NE.EQ.0) THEN
        IF(JNUM0(NA).EQ.1) GOTO 8500
        DO NSF=1,NSFMAX
           L=INSRF(NSF)
           NE=NESRF(NSF)
           KN=KNELM(L,NE)
           IF(KN.LE.0) THEN
              CALL CROS(XJ0(1,NA),YJ0(1,NA),ZJ0(1,NA),&
                   &    XJ0(2,NA),YJ0(2,NA),ZJ0(2,NA),&
                   &              NE,L,XC,YC,ZC,IERR)
              IF(IERR.EQ.0) THEN
                 LS=L
                 GOTO 1000
              ENDIF
           ENDIF
        ENDDO
        GOTO 8000
        
1000    CONTINUE
        N=1
        XJ(N,NA)=XC
        YJ(N,NA)=YC
        ZJ(N,NA)=ZC
        JELMT(N,NA)=NE
        IF(IDEBUG.NE.0) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,X,Y,Z=' ,&
                &NA,N,NE,XJ(N,NA),YJ(N,NA),ZJ(N,NA)
        ENDIF
     ELSE
        N=1
        XJ(N,NA)=XJ0(1,NA)
        YJ(N,NA)=YJ0(1,NA)
        ZJ(N,NA)=ZJ0(1,NA)
        JELMT(N,NA)=NE
        IF(IDEBUG.NE.0) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,X,Y,Z=',&
                &NA,N,NE,XJ(N,NA),YJ(N,NA),ZJ(N,NA)
        ENDIF
        LS=0
     ENDIF
     
     DO ID=2,JNUM0(NA)
        NENEXT=NE
        CALL FEP(XJ0(ID,NA),YJ0(ID,NA),ZJ0(ID,NA),NENEXT)
3000    CONTINUE
        IF(NENEXT.EQ.NE) GOTO 4500
        DO L=1,4
           IF(L.NE.LS) THEN
              CALL CROS(XJ (N ,NA),YJ (N ,NA),ZJ (N ,NA),&
                   &    XJ0(ID,NA),YJ0(ID,NA),ZJ0(ID,NA),&
                   &                 NE,L,XC,YC,ZC,IERR)
              IF(IERR.EQ.0) THEN
                 LS=L
                 GOTO 4000
              ENDIF
           ENDIF
        ENDDO
        GOTO 8100
        
4000    IF(N+1.GT.NJM) GOTO 8200
        N=N+1
        XJ(N,NA)=XC
        YJ(N,NA)=YC
        ZJ(N,NA)=ZC
        JELMT(N,NA)=NE
        IF(IDEBUG.NE.0) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,X,Y,Z=',&
                &NA,N,NE,XJ(N,NA),YJ(N,NA),ZJ(N,NA)
        ENDIF
        
        NENEW=KNELM(LS,NE)
        IF(NENEW.GT.0) THEN
           DO L=1,4
              IF(KNELM(L,NENEW).EQ.NE) LS=L
           ENDDO
           NE=NENEW
        ELSE
           IF(ID.EQ.JNUM0(NA).AND.NENEXT.LE.0) GOTO 6000
           GOTO 8400
        ENDIF
        GOTO 3000
        
4500    IF(N+1.GT.NJM) GOTO 8200
        N=N+1
        XJ(N,NA)=XJ0(ID,NA)
        YJ(N,NA)=YJ0(ID,NA)
        ZJ(N,NA)=ZJ0(ID,NA)
        JELMT(N,NA)=NE
        IF(IDEBUG.NE.0) THEN
           if(nrank.eq.0) WRITE(6,'(A,3I8,1P3E12.4)') 'NA,N,NE,X,Y,Z=',&
                & NA,N,NE,XJ(N,NA),YJ(N,NA),ZJ(N,NA)
        ENDIF
        LS=0
     ENDDO
     
6000 JNUM(NA)=N
  ENDDO
  IERR=0
  RETURN
  
8000 IERR=8000
  JNUM(NA)=N
  if(nrank.eq.0) WRITE(6,800) IERR,NA,N
800 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         &       ' ','           : CANNOT FIND BOUNDARY POINT'/&
         &       ' ','           : NA,N=',2I7)
  JNUM(NA)=N
  RETURN
  
8100 IERR=8100
  if(nrank.eq.0) WRITE(6,810) IERR,NA,ID,N
810 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         &       ' ','           : CANNOT FIND NEXT CROSSPOINT'/&
         &       ' ','           : NA,ID,N =',3I7)
  JNUM(NA)=N
  RETURN
  
8200 if(nrank.eq.0) WRITE(6,820) NA,ID,N,NJM
820 FORMAT(' ','## MODANT ERROR : N.GT.NJM '/&
         &       ' ','           : NA,ID,N,NJM = ',4I7)
  IERR=8200
  JNUM(NA)=N
  RETURN
  
8400 IERR=8400
  if(nrank.eq.0) WRITE(6,840) IERR,NA,ID,NE,N
840 FORMAT(' ','## MODANT ERROR : IERR =',I5/&
         &       ' ','           : ABMORMAL END OF ANTENNA DATA '/&
         &       ' ','           : NA,ID,NE,N = ',4I7)
  JNUM(NA)=N
  IERR=8400
  RETURN
  
8500 IERR=8500
  if(nrank.eq.0) WRITE(6,850) IERR,NA,NE,JNUM0(NA)
850 FORMAT(' ','## MODANT ERROR : IERR = ',I5/&
         &       ' ','           : NA,NE,JNUM0 = ',3I7)
  RETURN
  
END SUBROUTINE MODANT

!     ******* CALCULATE POINT OF INTERSECTION *******

SUBROUTINE CROS(X1,Y1,Z1,X2,Y2,Z2,IE,L,XC,YC,ZC,IERR)
  
  use wfcomm
  implicit none
  integer :: IE,L,I,IERR
  real(8) :: A(4),B(4),C(4),D(4),EPS,S,DV
  real(8) :: X1,X2,Y1,Y2,Z1,Z2,TC,XC,YC,ZC,FC

  DATA EPS/1.D-12/
  
  CALL WFABCD(IE,A,B,C,D,S)
  
  DV=B(L)*(X2-X1)+C(L)*(Y2-Y1)+D(L)*(Z2-Z1)
  IF(DV.EQ.0.D0) GOTO 9100
  TC=-(A(L)+B(L)*X1+C(L)*Y1+D(L)*Z1)/DV
  IF(TC.LT.-EPS.OR.TC.GT.1.D0+EPS) GOTO 9200
  
  XC=X1+(X2-X1)*TC
  YC=Y1+(Y2-Y1)*TC
  ZC=Z1+(Z2-Z1)*TC
  
  DO I=1,4
     FC=A(I)+B(I)*XC+C(I)*YC+D(I)*ZC
     IF(FC.LT.-EPS.OR.FC.GT.1.D0+EPS) GOTO 9300
  ENDDO
  
  IERR=0
  IF(IDEBUG.NE.0) THEN
     if(nrank.eq.0) then
        WRITE(6,*) '-- CROSS STATUS: SUCCESSFUL'
        WRITE(6,*) 'IE,L,TC= ',IE,L,TC
        WRITE(6,*) 'XC,YC,ZC= ',XC,YC,ZC
     end if
  ENDIF
  RETURN
  
9100 IERR=9100
  IF(IDEBUG.NE.0) THEN
     if(nrank.eq.0) then
        WRITE(6,*) '-- CROSS STATUS: VECTOR PARALLE TO SURFACE'
        WRITE(6,*) 'IE,L,DV= ',IE,L,DV
        WRITE(6,*) 'B,DELX = ',B(L),(X2-X1)
        WRITE(6,*) 'C,DELY = ',C(L),(Y2-Y1)
        WRITE(6,*) 'D,DELZ = ',D(L),(Z2-Z1)
     end if
  ENDIF
  RETURN
  
9200 IERR=9200
  IF(IDEBUG.NE.0) THEN
     if(nrank.eq.0) then
        WRITE(6,*) '-- CROSS STATUS: INTERSECTION BEYOND VECTOR'
        WRITE(6,*) 'IE,L= ',IE,L
        WRITE(6,*) 'TC = ',TC
     end if
  ENDIF
  RETURN
  
9300 IERR=9300
  IF(IDEBUG.NE.0) THEN
     if(nrank.eq.0) then
        WRITE(6,*) '-- CROSS STATUS: INTERSECTION OUT OF SURFACE'
        WRITE(6,*) 'IE,L,I = ',IE,L,I
        WRITE(6,*) 'FC = ',FC
     end if
  ENDIF
  RETURN
END SUBROUTINE CROS

!     ******* WEIGHT CALCULATION *******

SUBROUTINE WFWGT(X,Y,Z,IE,WGT)
  
  use wfcomm
  implicit none
  integer :: IE,L
  real(8) :: A(4),B(4),C(4),D(4),WGT(4),S,X,Y,Z
  
  CALL WFABCD(IE,A,B,C,D,S)
  DO L=1,4
     WGT(L)=A(L)+B(L)*X+C(L)*Y+D(L)*Z
  ENDDO
  RETURN
END SUBROUTINE WFWGT

!     ******* A,B,C,S  CALCULATION *******

SUBROUTINE WF2ABC(IN,NE,S)

  use wfcomm
  implicit none
  integer :: INS(3,4),NE,I,IN,J,K
  real(8) :: XE(4),YE(4),ZE(4),X1,Y1,Z1,X2,Y2,Z2,XS,YS,ZS,S
  DATA INS/4,3,2, 3,4,1, 2,1,4, 1,2,3/
  
  CALL WFNODE(NE,XE,YE,ZE)
  
  I=INS(1,IN)
  J=INS(2,IN)
  K=INS(3,IN)
  
  X1=XE(J)-XE(I)
  Y1=YE(J)-YE(I)
  Z1=ZE(J)-ZE(I)
  X2=XE(K)-XE(I)
  Y2=YE(K)-YE(I)
  Z2=ZE(K)-ZE(I)
  
  XS=Y1*Z2-Z1*Y2
  YS=Z1*X2-X1*Z2
  ZS=X1*Y2-Y1*X2
  
  S=0.5D0*SQRT(XS*XS+YS*YS+ZS*ZS)
  RETURN
END SUBROUTINE WF2ABC

!     ******* A,B,C,D,V  CALCULATION *******

SUBROUTINE WFABCD(NE,A,B,C,D,V)
  
  use wfcomm
  implicit none

  integer :: IL(3,4),NE,L,I,J,K
  real(8) :: A(4),B(4),C(4),D(4),XE(4),YE(4),ZE(4),V,VX,VY,VZ
  DATA IL/4,3,2, 3,4,1, 2,1,4, 1,2,3/
  
  CALL WFNODE(NE,XE,YE,ZE)
  
  V=  (XE(2)-XE(1))*(YE(3)-YE(1))*(ZE(4)-ZE(1))&
    &-(XE(2)-XE(1))*(YE(4)-YE(1))*(ZE(3)-ZE(1))&
    &+(XE(3)-XE(1))*(YE(4)-YE(1))*(ZE(2)-ZE(1))&
    &-(XE(3)-XE(1))*(YE(2)-YE(1))*(ZE(4)-ZE(1))&
    &+(XE(4)-XE(1))*(YE(2)-YE(1))*(ZE(3)-ZE(1))&
    &-(XE(4)-XE(1))*(YE(3)-YE(1))*(ZE(2)-ZE(1))
  
  IF(V.LE.0.D0) THEN
     if(nrank.eq.0) then
        WRITE(6,'(A,I5)') 'NEGATIVE V: NE=',NE
        WRITE(6,'(A,1P,3E12.4)') 'NODE 1 = ',XE(1),YE(1),ZE(1)
        WRITE(6,'(A,1P,3E12.4)') 'NODE 2 = ',XE(2),YE(2),ZE(2)
        WRITE(6,'(A,1P,3E12.4)') 'NODE 3 = ',XE(3),YE(3),ZE(3)
        WRITE(6,'(A,1P,3E12.4)') 'NODE 4 = ',XE(4),YE(4),ZE(4)
     end if
  ENDIF
  
  DO L=1,4
     I=IL(1,L)
     J=IL(2,L)
     K=IL(3,L)
     VX=(YE(J)-YE(I))*(ZE(K)-ZE(I))-(YE(K)-YE(I))*(ZE(J)-ZE(I))
     VY=(ZE(J)-ZE(I))*(XE(K)-XE(I))-(ZE(K)-ZE(I))*(XE(J)-XE(I))
     VZ=(XE(J)-XE(I))*(YE(K)-YE(I))-(XE(K)-XE(I))*(YE(J)-YE(I))
     
     V=VX*(XE(L)-XE(I))+VY*(YE(L)-YE(I))+VZ*(ZE(L)-ZE(I))
     A(L)=-(VX*XE(I)+VY*YE(I)+VZ*ZE(I))/V
     B(L)=VX/V
     C(L)=VY/V
     D(L)=VZ/V
  ENDDO
  V=V/6.D0
  RETURN
END SUBROUTINE WFABCD

!     ******* A,B,C,D,V  CALCULATION *******

SUBROUTINE WFABCDX(NE,F,DF,RWE,DWE,V)

  use wfcomm
  implicit none
  integer :: IL(3,4),INSID(4,6),NE,L,I,J,K,IN1,IN2,IN3,IN4
  real(8) :: F(4),DF(3,4),RWE(3,6),DWE(3,4,6),XE(4),YE(4),ZE(4),V,VX,VY,VZ,FACTOR

  DATA IL/4,3,2, 3,4,1, 2,1,4, 1,2,3/
  DATA INSID /1,2,3,4, 2,3,1,4, 1,3,4,2, &
       &            1,4,2,3, 2,4,3,1, 3,4,1,2/
  
  CALL WFNODE(NE,XE,YE,ZE)
  
  V=(XE(2)-XE(1))*(YE(3)-YE(1))*(ZE(4)-ZE(1))&
       & -(XE(2)-XE(1))*(YE(4)-YE(1))*(ZE(3)-ZE(1))&
       & +(XE(3)-XE(1))*(YE(4)-YE(1))*(ZE(2)-ZE(1))&
       & -(XE(3)-XE(1))*(YE(2)-YE(1))*(ZE(4)-ZE(1))&
       & +(XE(4)-XE(1))*(YE(2)-YE(1))*(ZE(3)-ZE(1))&
       & -(XE(4)-XE(1))*(YE(3)-YE(1))*(ZE(2)-ZE(1))
  
  IF(V.LE.0.D0) THEN
     if(nrank.eq.0) then
        WRITE(6,'(A,I5)') 'NEGATIVE V: NE=',NE
        WRITE(6,'(A,1P,3E12.4)') 'NODE 1 = ',XE(1),YE(1),ZE(1)
        WRITE(6,'(A,1P,3E12.4)') 'NODE 2 = ',XE(2),YE(2),ZE(2)
        WRITE(6,'(A,1P,3E12.4)') 'NODE 3 = ',XE(3),YE(3),ZE(3)
        WRITE(6,'(A,1P,3E12.4)') 'NODE 4 = ',XE(4),YE(4),ZE(4)
     end if
  ENDIF
  
  DO L=1,4
     I=IL(1,L)
     J=IL(2,L)
     K=IL(3,L)
     VX=(YE(J)-YE(I))*(ZE(K)-ZE(I))-(YE(K)-YE(I))*(ZE(J)-ZE(I))
     VY=(ZE(J)-ZE(I))*(XE(K)-XE(I))-(ZE(K)-ZE(I))*(XE(J)-XE(I))
     VZ=(XE(J)-XE(I))*(YE(K)-YE(I))-(XE(K)-XE(I))*(YE(J)-YE(I))
     
     F(L)=-(VX*XE(I)+VY*YE(I)+VZ*ZE(I))/V
     DF(1,L)=VX/V
     DF(2,L)=VY/V
     DF(3,L)=VZ/V
  ENDDO
  
  DO I=1,6
     IF(NSDELM(I,NE).GE.0) THEN
        FACTOR= 1.D0
     ELSE
        FACTOR=-1.D0
     ENDIF
     IN1=INSID(1,I)
     IN2=INSID(2,I)
     IN3=INSID(3,I)
     IN4=INSID(4,I)
     RWE(1,I)=FACTOR*(XE(IN4)-XE(IN3))*2/V
     RWE(2,I)=FACTOR*(YE(IN4)-YE(IN3))*2/V
     RWE(3,I)=FACTOR*(ZE(IN4)-ZE(IN3))*2/V
     DO J=1,4
        DO K=1,3
           DWE(K,J,I)=0.D0
        ENDDO
     ENDDO
     DWE(1,IN1,I)= FACTOR*DF(1,IN2)
     DWE(1,IN2,I)=-FACTOR*DF(1,IN1)
     DWE(2,IN1,I)= FACTOR*DF(2,IN2)
     DWE(2,IN2,I)=-FACTOR*DF(2,IN1)
     DWE(3,IN1,I)= FACTOR*DF(3,IN2)
     DWE(3,IN2,I)=-FACTOR*DF(3,IN1)
  ENDDO
  
  V=V/6.D0
  RETURN
END SUBROUTINE WFABCDX

!     ******* TOTAL COORDINATE - LOCAL COORDINATE *******

SUBROUTINE WFNODE(IE,XE,YE,ZE)
  
  use wfcomm
  implicit none
  integer :: I,IN,IE
  real(8) :: XE(4),YE(4),ZE(4)
  
  DO I=1,4
     IN=NDELM(I,IE)
     XE(I)=XND(IN)
     YE(I)=YND(IN)
     ZE(I)=ZND(IN)
  END DO

  RETURN
END SUBROUTINE WFNODE
  
!     ****** ELECTRIC FIELD AT ELEMENT(IE),POINT(X,Y,Z) ******

SUBROUTINE FIELDE(IE,X,Y,Z,CE)

  use wfcomm
  implicit none
  integer :: IE,IN,ISD,NSD
  real(8) :: W(4),F(4),DF(3,4),RWE(3,6),DWE(3,4,6),V,X,Y,Z
  complex(8) :: CE(3)

  CALL WFABCDX(IE,F,DF,RWE,DWE,V)
  
  CE(1)=(0.D0,0.D0)
  CE(2)=(0.D0,0.D0)
  CE(3)=(0.D0,0.D0)
  DO IN=1,4
     W(IN)=F(IN)+DF(1,IN)*X+DF(2,IN)*Y+DF(3,IN)*Z
  ENDDO
  DO ISD=1,6
     NSD=ABS(NSDELM(ISD,IE))
     DO IN=1,4
        CE(1)=CE(1)+CESD(NSD)*DWE(1,IN,ISD)*W(IN)
        CE(2)=CE(2)+CESD(NSD)*DWE(2,IN,ISD)*W(IN)
        CE(3)=CE(3)+CESD(NSD)*DWE(3,IN,ISD)*W(IN)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE FIELDE

!     ****** COMPLEX VALUE FIELD AT ELEMENT(IE),POINT(X,Y,Z) ******

!subroutine FIELDC(IE,X,Y,Z,CV,IDM,ID,FR,FI)
SUBROUTINE FIELDC(IE,X,Y,Z,CVALUE,IDM,ID,FR,FI)

  use wfcomm
  implicit none
  integer :: IE,N,IN,ID,IDM
  real(8) :: A(4),B(4),C(4),D(4),S,FR,FI,WEIGHT,X,Y,Z
!  complex(8):: CF,CV(IDM,NNMAX)
  complex(8):: CF,CVALUE(IDM,NNMAX)

  CALL WFABCD(IE,A,B,C,D,S)
  FR=0.D0
  FI=0.D0
  DO N=1,4
     IN=NDELM(N,IE)
     WEIGHT=A(N)+B(N)*X+C(N)*Y+D(N)*Z
!     CF=CV(ID,IN)
     CF=CVALUE(ID,IN)
     FR=FR+WEIGHT*DBLE(CF)
     FI=FI+WEIGHT*AIMAG(CF)
  END DO
  RETURN
END SUBROUTINE FIELDC

!     ****** REAL VALUE FIELD AT ELEMENT(IE),POINT(X,Y,Z) ******

SUBROUTINE FIELDD(IE,X,Y,Z,DV,ID,F)

  use wfcomm
  implicit none
  integer :: IE,N,IN,ID
  real(8) :: A(4),B(4),C(4),D(4),DV(NNMAX,ID),S,F,WEIGHT,X,Y,Z
  
  CALL WFABCD(IE,A,B,C,D,S)
  F=0.D0
  DO N=1,4
     IN=NDELM(N,IE)
     WEIGHT=A(N)+B(N)*X+C(N)*Y+D(N)*Z
     F=F+WEIGHT*DV(IN,ID)
     
  END DO
  RETURN
END SUBROUTINE FIELDD

!     ******* CEM INITIALIZE *******

SUBROUTINE SETAIF

  use wfcomm
  implicit none
  integer :: ID(4,4),I,L1,L2,L3,L4,J,K
  real(8) :: AIF,AIG

  DATA ID/1,4*0,1,4*0,1,4*0,1/
  
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
  
  DO I=1,4
     L1=ID(1,I)
     L2=ID(2,I)
     L3=ID(3,I)
     L4=ID(4,I)
     AIG1(I)=AIG(L1,L2,L3,L4)
     DO J=1,4
        L1=ID(1,I)+ID(1,J)
        L2=ID(2,I)+ID(2,J)
        L3=ID(3,I)+ID(3,J)
        L4=ID(4,I)+ID(4,J)
        AIG2(I,J)=AIG(L1,L2,L3,L4)
        DO K=1,4
           L1=ID(1,I)+ID(1,J)+ID(1,K)
           L2=ID(2,I)+ID(2,J)+ID(2,K)
           L3=ID(3,I)+ID(3,J)+ID(3,K)
           L4=ID(4,I)+ID(4,J)+ID(4,K)
           AIG3(I,J,K)=AIG(L1,L2,L3,L4)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE SETAIF

!     ******* INTEGRAL OF ELEMENT FUNCTION *******

FUNCTION AIF(L1,L2,L3)

  implicit none
  integer :: KAI(0:10),L1,L2,L3
  real(8) :: AIF

  DATA KAI/1,1,2,6,24,120,720,5040,40320,362880,3628800/
  
  AIF=DBLE(2*KAI(L1)*KAI(L2)*KAI(L3))/DBLE(KAI(L1+L2+L3+2))
  RETURN
END FUNCTION AIF

!     ******* INTEGRAL OF ELEMENT FUNCTION *******

FUNCTION AIG(L1,L2,L3,L4)

  implicit none
  integer :: KAI(0:10),L1,L2,L3,L4
  real(8) :: AIG
  DATA KAI/1,1,2,6,24,120,720,5040,40320,362880,3628800/
  
  AIG = DBLE(6*KAI(L1)*KAI(L2)*KAI(L3)*KAI(L4))&
      &/DBLE(KAI(L1+L2+L3+L4+3))
  RETURN
END FUNCTION AIG

!     ***** DUMP DATA ******

SUBROUTINE WFDUMP(ND)

  use wfcomm
  implicit none
  integer :: ND,NN,NE,I,NB,NSF

  WRITE(ND,*) NEMAX,NNMAX,NMMAX,NBMAX,NSFMAX
  WRITE(ND,*) '***** XND,YND,ZND'
  WRITE(ND,*) (XND(NN),YND(NN),ZND(NN),NN=1,NNMAX)
  WRITE(ND,*) '***** KANOD'
  WRITE(ND,*) (KANOD(NN),NN=1,NNMAX)
!  WRITE(ND,*) '***** VTOT,VNOD'
!  WRITE(ND,*) VTOT,(VNOD(NN),NN=1,NNMAX)
  WRITE(ND,*) '***** NDELM'
  WRITE(ND,*) ((NDELM(I,NE),I=1,5),NE=1,NEMAX)
  WRITE(ND,*) '***** KAELM,NM,NBELM'
  WRITE(ND,*) (KAELM(NE),NMKA(KAELM(NE)),NBELM(NE),NE=1,NEMAX)
  WRITE(ND,*) '***** KNELM'
  WRITE(ND,*) ((KNELM(I,NE),I=1,4),NE=1,NEMAX)
  WRITE(ND,*) '***** XEMIN,YEMIN,ZEMIN'
  WRITE(ND,*) (XEMIN(NE),YEMIN(NE),ZEMIN(NE),NE=1,NEMAX)
  WRITE(ND,*) '***** XEMAX,YEMAX,ZEMAX'
  WRITE(ND,*) (XEMAX(NE),YEMAX(NE),ZEMAX(NE),NE=1,NEMAX)
  WRITE(ND,*) '***** SINDEX,IVELM,IWELM'
  WRITE(ND,*) (SINDEX(NE),IVELM(NE),IWELM(NE),NE=1,NEMAX)
  WRITE(ND,*) '***** KABDY,NDBDY,NMRBDY'
  WRITE(ND,*) (KABDY(NB),NDBDY(NB),NMBDY(NB),NB=1,NBMAX)
  WRITE(ND,*) '***** PHIBDY,RESBDY,PWRBDY'
  WRITE(ND,*) (PHIBDY(NB),RESBDY(NB),PWRBDY(NB),NB=1,NBMAX)
  WRITE(ND,*) '***** XGBDY,YGBDY,ZGBDY'
  WRITE(ND,*) (XGBDY(NB),YGBDY(NB),ZGBDY(NB),NB=1,NBMAX)
  WRITE(ND,*) '***** XPBDY,YPBDY,ZPBDY'
  WRITE(ND,*) (XPBDY(NB),YPBDY(NB),ZPBDY(NB),NB=1,NBMAX)
  WRITE(ND,*) '***** XNBDY,YNBDY,ZNBDY'
  WRITE(ND,*) ((XNBDY(I,NB),YNBDY(I,NB),ZNBDY(I,NB),I=1,3),NB=1,NBMAX)
  WRITE(ND,*) '***** PHABDY,SZBDY'
  WRITE(ND,*) (PHABDY(NB),(SZBDY(I,NB),I=1,2),NB=1,NBMAX)
  WRITE(ND,*) '***** INSRF,NESRF'
  WRITE(ND,*) (INSRF(NSF),NESRF(NSF),NSF=1,NSFMAX)
  WRITE(ND,*) '***** NDSRF,KNSRF'
  WRITE(ND,*) ((NDSRF(I,NSF),KNSRF(I,NSF),I=1,3),NSF=1,NSFMAX)
  WRITE(ND,*) '***** MBND,MLEN,NNBMAX'
  WRITE(ND,*) MBND,MLEN,NNBMAX
  WRITE(ND,*) '***** ISDELM,IMLEN,INLEN'
  WRITE(ND,*) ((ISDELM(I,NE),I=1,7),NE=1,NEMAX)
  WRITE(ND,*) '***** IMLEN,INLEN'
  WRITE(ND,*) (IMLEN(NN),INLEN(NN),NN=1,NNBMAX)
  RETURN
END SUBROUTINE WFDUMP


