!     $Id$
!
!     ********** F.E.M. DIVIDER ( FIRST ORDER ) **********
!
subroutine WFDIV

  use libmpi
  use libmtx
  USE libchar
  use wfcomm
  implicit none
  integer :: IERR,NE,NN,ID
  character KID*1

1 continue

  if (nrank.eq.0) then  
     write(6,601) 
601  format('## INPUT: D/DIV  G/DRAW  P,V/PARM'&
                    '  S/SAVE  L/LOAD  W/LIST  X/EXIT')
     read(5,'(A1)',ERR=1,END=9000) KID
     call toupper(KID)
  end if
  call mtx_barrier
  call mtx_broadcast_character(KID,1)
  
2 continue

  if(KID.eq.'D') then
     call wfdiv_initialize
     call wfdiv_allocate

     if (nrank.eq.0) then
        write(6,602) 
602     format('## TYPE:  X/RECT  C/CIRCLE  A/COAXIAL')
        read(5,'(A1)',ERR=2,END=1) KID
        call toupper(KID)
        
        if(KID.eq.'X') then
3          write(6,603) BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
603        format('## DIV:   BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX = '/&
                  '          ',6F10.4/&
                  '## INPUT: BXMIN,BXMAX,BYMAX,BYMAX,BZMIN,BZMAX ? ')
           read(5,*,ERR=3,END=2) BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
           IDDIV=0

        elseif(KID.eq.'C') then
4          write(6,604) RB,BZMIN,BZMAX
604        format('## DIV:   RB,BZMIN,BZMAX = ',3F10.4/&
                  '## INPUT: RB,BZMIN,BZMAX ? ')
           read(5,*,ERR=4,END=2) RB,BZMIN,BZMAX
           BXMIN=-RB
           BXMAX= RB
           BYMIN=-RB
           BYMAX= RB
           IDDIV=1

        elseif(KID.eq.'A') then
5          write(6,605) RB,RBAX,BZMIN,BZMAX
605        format('## DIV:   RB,RBAX,BZMIN,BZMAX = ',4F10.4/&
                  '## INPUT: RB,RBAX,BZMIN,BZMAX ? ')
           read(5,*,ERR=5,END=2) RB,RBAX,BZMIN,BZMAX
           BXMIN=-RB
           BXMAX= RB
           BYMIN=-RB
           BYMAX= RB
           IDDIV=2
!
! ----- Add. By YOKOYAMA 28/02/2013 ----
!        GAMMA10軸対称セントラル部を模擬
!        真空容器内(r=0->RBOUT)，コアプラズマ部(r=0->RBIN)と
!        周辺プラズマ部(r=RBIN->RBOUT)の各々で，
!        要素長を変更できるようにする
       elseif(KID.eq.'G') then
            RBIN  = 0.24D0
            RBOUT = 0.50D0
6           write(6,608) RBIN,RBOUT,BZMIN,BZMIDL,BZMIDH,BZMAX
608        format('## DIV:   RBIN,RBOUT                = ',2F10.4/&
                  '##        BZMIN,BZMIDL,BZMIDH,BZMAX = ',4F10.4/&
                  '## INPUT: RBIN,RBOUT,BZMIN,BZMIDL,BZMIDH,BZMAX ? ')
!            READ(5,*,ERR=6,END=2) RBIN,RBOUT,BZMIN,BZMAX
            read(5,*,ERR=6,END=2) RBIN,RBOUT,BZMIN,BZMIDL,BZMIDH,BZMAX
            BXMIN=-RBOUT
            BXMAX= RBOUT
            BYMIN=-RBOUT
            BYMAX= RBOUT
            IDDIV=3
!
! ----- -----
!
        else
           write(6,*) 'XX UNKNOWN KID: ',KID
        end if
!
! ----- Add. By YOKOYAMA 28/02/2013 ----
!
!        要素長の取得
!        コアプラズマ部(r=0    -> RBIN ):  DXIN， DYIN
!        周辺プラズマ部(r=RBIN -> RBOUT):  DXOUT，DYOUT
7        if(KID.ne.'G') then
            write(6,607) DELX,DELY,DELZ
607         format('## DIV:   DELX,DELY,DELZ = ',3F10.4/&
                   '## INPUT: DELX,DELY,DELZ ? ')
            read(5,*,ERR=7,END=2) DELX,DELY,DELZ
            if(ABS(DELX).LE.1.D-6.OR.ABS(DELY).LE.1.D-6) GOTO 2
         else
            DXIN  = 0.03D0
            DYIN  = 0.03D0
            DXOUT = 0.30D0
            DYOUT = 0.15D0
            write(6,609) DXIN,DYIN,DXOUT,DYOUT,DELZ,DELZM
609        format('## DIV:   DXIN,DYIN,DXOUT,DYOUT = ',4F10.4/&
                  '##        DELZ,DELZM            = ',2F10.4/&
                  '## INPUT: DXIN,DYIN,DXOUT,DYOUT,DELZ,DELZM ? ')
!            READ(5,*,ERR=7,END=2) DXIN,DYIN,DXOUT,DYOUT,DELZ
            read(5,*,ERR=7,END=2) DXIN,DYIN,DXOUT,DYOUT,DELZ,DELZM
            if((ABS(DXIN) .LE.1.D-6.OR.ABS(DYIN) .LE.1.D-6).OR. &
              (ABS(DXOUT).LE.1.D-6.OR.ABS(DYOUT).LE.1.D-6)) GOTO 2
!
! ----- -----
!
         end if
     end if

     call wfdiv_broadcast

     call DFNODE(IERR)
     if(IERR.ne.0) goto 2
     call SETNOD(IERR)
     if(IERR.ne.0) goto 2
!
! ----- Add. By YOKOYAMA 01/03/2013 ----
!
!        "SETELM"と"SETELMX"の違いを調べる必要がある
!        同軸の穴の空いた構造が原因か...?
!        GAMMA10セントラル部の場合は"SETELM"を使うべきか...?
!        
!         IF(IDDIV.NE.2) THEN
     if((IDDIV.ne.2).and.(IDDIV.ne.3)) then
! ----- 01/03/2013 -----
!
        call SETELM(IERR)
     else
        call SETELMX(IERR)
     end if
     if(IERR.ne.0) GOTO 2

     if(nrank.eq.0) write(6,*) '--- WFINDX start ---'
     write(6,'(A,2I5)') '--nsize,nrank:',nsize,nrank
     call WFINDX
     if(nrank.eq.0) write(6,*) '--- WFFEPI start ---'
     call WFFEPI
     
     NKMAX=1
     do NE=1,NEMAX
        KAELM(NE)=1
     end do
     NMKA(1)=0
     NMMAX=0
     
     NBMAX=0
     do NN=1,NNMAX
        KANOD(NN)=0
     end do

     call wfdiv_deallocate

  elseif(KID.eq.'G') then
!
! ----- Add. By YOKOYAMA 01/03/2013 ----
!
!        Y-Z平面上の要素分割を表示する
         call WFGDIV_YZ
!        X-Z平面上の要素分割を表示する
         call WFGDIV_XZ
! ----- 01/03/2013 -----
!
     if (nrank.eq.0) then
        call WFGDIV
        if(NDRAWD.eq.-1) then
           call WFGNAS(0)
           call WFGNAS(1)
           call WFGNAS(2)
           call WFGNAS(3)
           call WFGNAS(4)
           call WFGNAS(5)
        end if
     end if

  elseif(KID.eq.'W') then
     if (nrank.eq.0) call WFLDIV


  elseif(KID.eq.'L') then
     call WFRELM(ID)

  elseif(KID.eq.'P') then
     if(nrank.eq.0) call WFPARM(KID)
     call wfparm_broadcast

  elseif(KID.eq.'V') then
     if (nrank.eq.0) call WFVIEW

  elseif(KID.eq.'S') then
     if (nrank.eq.0) call WFWELM(0)

  elseif(KID.eq.'X') then
     goto 9000
  end if
  goto 1

9000 continue
  return
end subroutine WFDIV

!     ****** Definititon of Boundary ******
function BOUNDF(X,Y)

  use wfcomm
  implicit none
  real(8) :: BOUNDF,X,Y
  
  if(IDDIV.eq.0) then
     BOUNDF=(X-BXMIN)*(BXMAX-X)*(Y-BYMIN)*(BYMAX-Y)
  elseif(IDDIV.eq.1) then
     BOUNDF=(RB*RB-X*X-Y*Y)
  end if

  return
end function BOUNDF

!     ****** Definititon of Boundary for Fixed Y ******
function BOUNDX(X)
  
  use wfcomm
  implicit none
  real(8) :: BOUNDX,X,BOUNDF
  BOUNDX=BOUNDF(X,YF)

  return
end function BOUNDX

!     ****** Definititon of Boundary for Fixed X ******

function BOUNDY(Y)

  use wfcomm
  implicit none
  real(8) :: BOUNDY,BOUNDF,Y
  BOUNDY=BOUNDF(XF,Y)
  
  return
end function BOUNDY

!     ****** Define 2-D Node Array ******

subroutine DFNODE(IERR)

  use wfcomm
  implicit none
  integer :: IERR

  if(IDDIV.eq.0) then
     call DFNODX(IERR)
  elseif(IDDIV.eq.1) then
     call DFNODC(IERR)
  elseif(IDDIV.eq.2) then
     call DFNODCX(IERR)
  elseif(IDDIV.eq.3) then
     call DFNODC2(IERR)
  end if
  
  return
end subroutine DFNODE

!     ****** Define 2-D Node Array (RECTANGULAR) ******

subroutine DFNODX(IERR)

  use wfcomm
  implicit none
  integer :: NY,NXMAX,NX,IERR
  real(8) :: DY,Y,DX,X,XREL,YREL,R,FACTOR
  
  NYMAX=NINT((BYMAX-BYMIN)/DELY)+1
  NYMAX=NYMAX+MOD(NYMAX+1,2)
  if(NYMAX.gt.NYM) goto 9200
  DY=(BYMAX-BYMIN)/DBLE(NYMAX-1)
  
  do NY=1,NYMAX
     Y=BYMIN+DY*(NY-1)
     
     XL(NY)=BXMIN
     XR(NY)=BXMAX
     NXMAX=NINT((XR(NY)-XL(NY))/DELX)+1
     NXMAX=NXMAX+MOD(NXMAX+1,2)
     if(NXMAX.gt.NXM) goto 9300
     DX=(XR(NY)-XL(NY))/DBLE(NXMAX-1)
     
     NXA(NY)=NXMAX
     do NX=1,NXMAX
        X=DX*(NX-1)+XL(NY)
        XREL=DX*DBLE(NX-(NXMAX+1)/2)
        YREL=DY*DBLE(NY-(NYMAX+1)/2)
        R=SQRT(XREL**2+YREL**2)
        FACTOR=1.D0-R*1.D-6
        XNDA(NX,NY)=X*FACTOR
        YNDA(NX,NY)=Y*FACTOR
     end do
  end do

  IERR=0
  return
        
9200 if (nrank.eq.0) then
     write(6,602) NYMAX,NYM
602  format(' ','DFNODE : NYMAX EXCEEDS NYM : ',2I8)
  end if
  IERR=1
  return
  
9300 if (nrank.eq.0) then
     write(6,603) NY,NXMAX,NXM
603  format(' ','DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
  end if
  IERR=1
  return
end subroutine DFNODX

!     ****** Define 2-D Node Array (CIRCULAR) ******

subroutine DFNODC(IERR)

  use libmtx
  use wfcomm
  implicit none
  integer :: IND,NYDO,NY,ILL,NXA1(NYMH),NXR1(NYMH),NXL1(NYMH),NYU(NXM)
  integer :: NYMID,NX,NLMID,NRMID,NXMID,INDC,NYUL,I,NYA,NYB,IERR
  real(8) :: XS,XE,DX,DXL,DXR,XL1(NYMH),XR1(NYMH)
  real(8) :: XNDA1(NXM,NYMH),YNDA1(NXM,NYMH),XU(NXM),YU(NXM),DYU(NXM)
  real(8) :: YMID,YS,YE,DY
  real(8),parameter:: EPS=1.d-8

  EXTERNAL BOUNDX,BOUNDY

  NYMAX=NINT(BYMAX/DELY)+1
  if(NYMAX.gt.NYM) goto 9200
  IND=0

  do NYDO=1,NYMAX
     NY=NYDO
     YF=DELY*DBLE(NY-1)
     if(NY.eq.1) then
        XS=0.5D0*(BXMIN+BXMAX)
     else
        XS=0.5D0*(XL1(NY-1)+XR1(NY-1))
     end if
     XE=BXMIN-1.D-5
     DX=-DELX*0.3D0
     ILL=0
     call FRGFLS(XS,XE,DX,XL1(NY),BOUNDX,EPS,ILL)
     if(ILL.eq.902) then
        NY=NY-1
        goto 20
     end if
     if(ILL.gt.0) goto 9000
     XE=BXMAX+1.D-5
     DX=DELX*0.3D0
     ILL=0
     call FRGFLS(XS,XE,DX,XR1(NY),BOUNDX,EPS,ILL)
     if(ILL.eq.902) then
        NY=NY-1
        goto 20
     end if
     if(ILL.gt.0) goto 9000
     if(NY.gt.1) then
        DXL=ABS(XL1(NY)-XL1(NY-1))/DELX
        DXR=ABS(XR1(NY)-XR1(NY-1))/DELX
        if(DXL**2+DXR**2.lt.1.D0) then
           NYMID=NY
        else
           IND=1
        end if
     end if
  end do
20 continue
  
  NYMAX=MIN(NYMAX,NY)
  
  do NY=1,NYMID
     NXA1(NY)=NINT((XR1(NY)-XL1(NY))/DELX)+1
     if(NXA1(NY).gt.NXM) goto 9300
     if(NXA1(NY).eq.1) then
        NXL1(NY)=1
        NXR1(NY)=0
        XNDA1(1,NY)=XL1(NY)
        YNDA1(1,NY)=DELY*DBLE(NY-1)
        IND=0
        goto 1010
     end if
     if(MOD(NXA1(NY),2).eq.0) NXA1(NY)=NXA1(NY)+1
     DX=(XR1(NY)-XL1(NY))/DBLE(NXA1(NY)-1)
     NXL1(NY)=1
     NXR1(NY)=1
     do NX=1,NXA1(NY)
        XNDA1(NX,NY)=DX*DBLE(NX-1)+XL1(NY)
        YNDA1(NX,NY)=DELY*DBLE(NY-1)
     end do
  end do
1010 continue
  
  if(IND.ge.1.and.NYMID.lt.NYMAX) then
     NLMID=NXL1(NYMID)
     NRMID=NXR1(NYMID)
     NXMID=NXA1(NYMID)
     YMID =DELY*DBLE(NYMID-1)
     do NX=NLMID+1,NXMID-NRMID
        XU(NX)=XNDA1(NX,NYMID)
        XF=XU(NX)
        YS=YMID
        YE=BYMAX+1.D-5
        DY=DELY
        ILL=0
        call FRGFLS(YS,YE,DY,YU(NX),BOUNDY,EPS,ILL)
        if(ILL.gt.0) goto 9400
        NYU(NX)=MAX(1,NINT((YU(NX)-YMID)/DELY))
        DYU(NX)=(YU(NX)-YMID)/DBLE(NYU(NX))
        NYMAX=MAX(NYMAX,NYU(NX)+NYMID)
     end do
     
     do NY=NYMID+1,NYMAX
        NX=0
        INDC=0
        NXL1(NY)=0
        NXR1(NY)=0
        NYUL=NY-NYMID
        do I=NLMID+1,NXMID-NRMID
           if(NYU(I).eq.NYUL) then
              NX=NX+1
              XNDA1(NX,NY)=XNDA1(I,NYMID)
              YNDA1(NX,NY)=YU(I)
              if(INDC.eq.0) then
                 NXL1(NY)=NXL1(NY)+1
              else
                 NXR1(NY)=NXR1(NY)+1
              end if
           else if(NYU(I).gt.NYUL) then
              NX=NX+1
              XNDA1(NX,NY)=XNDA1(I,NYMID)
              YNDA1(NX,NY)=DYU(I)*DBLE(NYUL)+YMID
              INDC=1
           end if
        end do
        if(NX.eq.0) goto 4305
        NXA1(NY)=NX
     end do
     goto 4310
4305 NYMAX=NY-1
4310 continue
  else
     NXL1(NYMAX)=NXA1(NYMAX)
     NXR1(NYMAX)=0
  end if
  
  do NY=1,NYMAX
     NYA=NYMAX+1-NY
     NYB=NYMAX-1+NY
     XL(NYA)=XL1(NY)
     XL(NYB)=XL1(NY)
     XR(NYA)=XR1(NY)
     XR(NYB)=XR1(NY)
     NXA(NYA)=NXA1(NY)
     NXA(NYB)=NXA1(NY)
     do NX=1,NXA1(NY)
        XNDA(NX,NYA)= XNDA1(NX,NY)
        XNDA(NX,NYB)= XNDA1(NX,NY)
        YNDA(NX,NYA)=-YNDA1(NX,NY)
        YNDA(NX,NYB)= YNDA1(NX,NY)
     end do
  end do
  NYMAX=2*NYMAX-1
  
  IERR=0
  return
  
9000 if (nrank.eq.0) then
     write(6,601) NX,NY,YF,XS,XE,XL(NY),XR(NY),ILL
601  format(' ','DFNODE : FRGFLS ERROR : NX,NY,YF,XS,XE,XL,XR,ILL'/&
          &       ' ',2I5,1P5E15.5,I5)
  end if
  IERR=1
  return
  
9200 if (nrank.eq.0) then
     write(6,602) NYMAX,NYM
602  format(' ','DFNODE : NYMAX EXCEEDS NYM : ',2I8)
  end if
  IERR=1
  return
  
9300 if (nrank.eq.0) then
     write(6,603) NY,NXA(NY),NXM
603  format(' ','DFNODE : NX EXCEEDS NXM AT NY =',I8,' : ',2I8)
  end if
  IERR=1
  return
  
9400 if (nrank.eq.0) then
     write(6,604) NX,NY,XF,YS,YE,YU(NX),ILL
604  format(' ','NODE1 : FRGFLS ERROR : NX,NY,XF,YS,YE,YU,ILL'/&
            ' ',2I5,1P4E15.5,I5)
  end if
  IERR=1
  return
end subroutine DFNODC

!     ****** Define 2-D Node Array (CIRCULAR COAXIAL) ******

subroutine DFNODCX(IERR)

  use wfcomm
  implicit none
  integer :: NY,NXMAX,NX,IERR
  real(8) :: DY,R,DX,TH

  NYMAX=NINT((RB-RBAX)/DELY)+1
  if(NYMAX.gt.NYM) goto 9100
  if(NYMAX.le.1) goto 9100
  DY=(RB-RBAX)/(NYMAX-1)
  
  do NY=1,NYMAX
     R=RBAX+DY*(NY-1)
     if(R.le.0.D0) then
        NXA(NY)=1
        XNDA(1,NY)=0.D0
        YNDA(1,NY)=0.D0
     else
        NXMAX=NINT((2.D0*PI*R)/DELX)+1
        if(NXMAX.gt.NXM) goto 9200
        NXMAX=4*((NXMAX-1)/4+1)
        if(NXMAX.eq.4) NXMAX=8
        DX=2.D0*PI/NXMAX
        NXA(NY)=NXMAX
        do NX=1,NXMAX
           TH=DX*(NX-1)
           XNDA(NX,NY)=R*COS(TH)
           YNDA(NX,NY)=R*SIN(TH)
        end do
     end if
  end do
  
  IERR=0
  return
  
9100 if (nrank.eq.0) then
     write(6,601) NYMAX,NYM
601  format(' ','DFNODCX : NYMAX EXCEEDS NYM : ',2I8)
  end if
  IERR=1
  return
  
9200 if (nrank.eq.0) then
     write(6,602) NY,NXMAX,NXM
602  format(' ','DFNODCX : NXMAX EXCEEDS NXM AT NX=',I8,' : ',2I8)
  end if
  IERR=1
  return
end subroutine DFNODCX


!  ----  Add. by Yama  06/Aug./2007 ----

!     ****** Define 2-D Node Array (CIRCULAR /2-Sections) ******
!     半径0からRBINまでの領域と，RBINからRBOUTまでの領域を，
!     異なる要素長で分割する

SUBROUTINE DFNODC2(IERR)

  use wfcomm
  implicit none
  integer :: NY,NYMID,NYEND,NYN,NXMAX,NX,IERR
  real(8) :: DY,R,DX,TH

!  -- R方向0からRBINまでの分割 --
!     NYMID: R方向の分割数
!     DY: R方向の要素長
      NYMID=NINT(RBIN/DYIN)+1
         IF(NYMID.GT.NYM) GOTO 9100
         IF(NYMID.LE.1)   GOTO 9100
      DY=RBIN/(NYMID-1)

      DO NY=1,NYMID
         R=DY*(NY-1)
!        R=0の時，節点には座標系の原点(0,0)を割り当てる
         IF(R.LE.0.D0) THEN
!           NXA: あるY座標を持つ節点の数（X座標の数）
            NXA(NY)=1
            XNDA(1,NY)=0.D0
            YNDA(1,NY)=0.D0
         ELSE
!           NXMAX: 方位角方向の分割数（円周の長さ/要素長DXIN）
            NXMAX=NINT((2.D0*PI*R)/DXIN)+1
            IF(NXMAX.GT.NXM) GOTO 9200
!           NXMAX=1-7ならば，NXMAX=8とおく
!              要素分割数が少な過ぎる場合の対処か
!              方位角方向の要素分割数NXMAXの最低値を8としている
            NXMAX=4*((NXMAX-1)/4+1)
            IF(NXMAX.EQ.4) NXMAX=8
!           DX: 方位角方向の要素長[rad.]
            DX=2.D0*PI/NXMAX
            NXA(NY)=NXMAX
            DO NX=1,NXMAX
               TH=DX*(NX-1)
               XNDA(NX,NY)=R*COS(TH)
               YNDA(NX,NY)=R*SIN(TH)
            ENDDO
         ENDIF
!         WRITE(6,'(I5,1PE12.4,I5,1PE12.4)') NY,R,NXMAX,DX
      ENDDO

!  -- R方向RBINからRBOUTまでの分割 --

      NYEND=NINT((RBOUT-RBIN)/DYOUT)
         NYMAX=NYMID+NYEND
         IF(NYMAX.GT.NYM) GOTO 9100
         IF(NYMAX.LE.1)   GOTO 9100
      DY=(RBOUT-RBIN)/NYEND

      NYN=0
      DO NY=NYMID+1,NYMAX
         NYN=NYN+1
         R=RBIN+DY*NYN
         NXMAX=NINT((2.D0*PI*R)/DXOUT)+1
         IF(NXMAX.GT.NXM) GOTO 9200
         NXMAX=4*((NXMAX-1)/4+1)
         IF(NXMAX.EQ.4) NXMAX=8
         DX=2.D0*PI/NXMAX
         NXA(NY)=NXMAX
         DO NX=1,NXMAX
            TH=DX*(NX-1)
            XNDA(NX,NY)=R*COS(TH)
            YNDA(NX,NY)=R*SIN(TH)
         ENDDO
      ENDDO

      IERR=0
      RETURN

 9100 WRITE(6,601) NYMAX,NYM
  601 FORMAT(' ','DFNODC2 : NYMAX EXCEEDS NYM : ',2I8)
      IERR=1
      RETURN

 9200 WRITE(6,602) NY,NXMAX,NXM
  602 FORMAT(' ','DFNODC2 : NXMAX EXCEEDS NXM AT NX=',I8,' : ',2I8)
      IERR=1
      RETURN
END SUBROUTINE DFNODC2

!     ****** Set Node Array ******

subroutine SETNOD(IERR)
  
  use wfcomm
  implicit none
  integer :: NZ,IN,INMAX,NY,NX,IERR
  real(8) :: DZ,FACTOR

  NZMAX=NINT((BZMAX-BZMIN)/DELZ)+1
  if(NZMAX.gt.NZM) goto 9100
  DZ=(BZMAX-BZMIN)/(NZMAX-1)
  do NZ=1,NZMAX
     ZNDA(NZ)=(NZ-1)*DZ+BZMIN
  end do

  ! --- decide NNMAX ---
  IN=0
  do NY=1,NYMAX
     do NX=1,NXA(NY)
        do NZ=1,NZMAX
           IN=IN+1
        end do
     end do
  end do
  NNMAX=IN

  call wfelm_allocate
  
  IN=0
  do NY=1,NYMAX
     do NX=1,NXA(NY)
        do NZ=1,NZMAX
           IN=IN+1
           NDA(NX,NY,NZ)=IN
           FACTOR=1.D0

           IF (MODELG.eq.1) then
              CALL RTTORC(XNDA(NX,NY),YNDA(NX,NY),ZNDA(NZ),XND(IN),YND(IN),ZND(IN))
           else
              XND(IN)= XNDA(NX,NY)*FACTOR
              YND(IN)= YNDA(NX,NY)*FACTOR
              ZND(IN)= ZNDA(NZ)
           end IF
        end do
     end do
  end do
  IERR=0
  return
  
9100 if (nrank.eq.0) then
     write(6,602) NZMAX,NZM
602  format(' ','SETNOD : NZMAX EXCEEDS NZM : ',2I8)
  end if
  IERR=2
  return
end subroutine SETNOD

!  ****** Set Element Array ******

subroutine SETELM(IERR)

  use wfcomm
  implicit none
  integer :: IE,NYDO,NY,NXMAX,NX1MAX,NX,NX1,NZ,IERR
  real(8) :: XAI,XAJ,XBI,XBJ,YAI,YAJ,YBI,YBJ,VA,VB

  IERR=0
  IE=0

! --- decide NEMAX ---
  do NYDO=1,NYMAX-1
     NY=NYDO
     NXMAX=NXA(NY)
     NX1MAX=NXA(NY+1)
     NX=1
     NX1=1
     
     XAI=XNDA(NX,  NY   )
18   XBJ=XNDA(NX1+1,NY+1)
     if(XBJ.lt.XAI) then
        NX1=NX1+1
        goto 18
     end if
     XAJ=XNDA(NX1,  NY+1)
19   XBI=XNDA(NX+1,NY   )
     if(XBI.lt.XAJ) then
        NX=NX+1
        goto 19
     end if
20   XAI=XNDA(NX,  NY   )
     XBI=XNDA(NX+1,NY   )
     XAJ=XNDA(NX1,  NY+1)
     XBJ=XNDA(NX1+1,NY+1)
     YAI=YNDA(NX,  NY   )
     YBI=YNDA(NX+1,NY   )
     YAJ=YNDA(NX1,  NY+1)
     YBJ=YNDA(NX1+1,NY+1)
     VA=(XAI-XBJ)**2+(YAI-YBJ)**2
     VB=(XBI-XAJ)**2+(YBI-YAJ)**2
     if(VA.gt.VB) then
        do NZ=1,NZMAX-1
           IE=IE+3
        end do
        NX=NX+1
     else
        do NZ=1,NZMAX-1
           IE=IE+3
        end do
        NX1=NX1+1
     end if
     if(NX.lt.NXMAX) then
        if(NX1.lt.NX1MAX) then
           goto 20
        else
21         continue
           XAI=XNDA(NX, NY  )
           XAJ=XNDA(NX1,NY+1)
           if(XAI.le.XAJ) then
              do NZ=1,NZMAX-1
                 IE=IE+3
              end do
              NX=NX+1
              if(NX.lt.NXMAX) goto 21
           end if
        end if
     else
        if(NX1.lt.NX1MAX) then
22         continue
           XAI=XNDA(NX,  NY   )
           XAJ=XNDA(NX1, NY+1)
           if(XAJ.le.XAI) then
              do NZ=1,NZMAX-1
                 IE=IE+3
              end do
              NX1=NX1+1
              if(NX1.lt.NX1MAX) goto 22
           end if
        end if
     end if
  end do
  NEMAX=IE

  call wfelm_allocate

  IE = 0
  
  do NYDO=1,NYMAX-1
     NY=NYDO
     NXMAX=NXA(NY)
     NX1MAX=NXA(NY+1)
     NX=1
     NX1=1
     
     XAI=XNDA(NX,  NY   )
23   XBJ=XNDA(NX1+1,NY+1)
     if(XBJ.lt.XAI) then
        NX1=NX1+1
        goto 23
     end if
     XAJ=XNDA(NX1,  NY+1)
24   XBI=XNDA(NX+1,NY   )
     if(XBI.lt.XAJ) then
        NX=NX+1
        goto 24
     end if
25   XAI=XNDA(NX,  NY   )
     XBI=XNDA(NX+1,NY   )
     XAJ=XNDA(NX1,  NY+1)
     XBJ=XNDA(NX1+1,NY+1)
     YAI=YNDA(NX,  NY   )
     YBI=YNDA(NX+1,NY   )
     YAJ=YNDA(NX1,  NY+1)
     YBJ=YNDA(NX1+1,NY+1)
     VA=(XAI-XBJ)**2+(YAI-YBJ)**2
     VB=(XBI-XAJ)**2+(YBI-YAJ)**2
     if(VA.gt.VB) then
        do NZ=1,NZMAX-1
           call SETELL(NX,NX1,NY,NZ,IE,0)
        end do
        NX=NX+1
     else
        do NZ=1,NZMAX-1
           call SETELL(NX,NX1,NY,NZ,IE,1)
        end do
        NX1=NX1+1
     end if
     if(NX.lt.NXMAX) then
        if(NX1.lt.NX1MAX) then
           goto 25
        else
26         continue
           XAI=XNDA(NX, NY  )
           XAJ=XNDA(NX1,NY+1)
           if(XAI.le.XAJ) then
              do NZ=1,NZMAX-1
                 call SETELL(NX,NX1,NY,NZ,IE,0)
              end do
              NX=NX+1
              if(NX.lt.NXMAX) goto 26
           end if
        end if
     else
        if(NX1.lt.NX1MAX) then
27         continue
           XAI=XNDA(NX,  NY   )
           XAJ=XNDA(NX1, NY+1)
           if(XAJ.le.XAI) then
              do NZ=1,NZMAX-1
                 call SETELL(NX,NX1,NY,NZ,IE,1)
              end do
              NX1=NX1+1
              if(NX1.lt.NX1MAX) goto 27
           end if
        end if
     end if
  end do

  IERR=0
  return
end subroutine SETELM

!     ****** Set Element Array SUB******

subroutine SETELL(NX,NX1,NY,NZ,IE,ID)

  use wfcomm
  implicit none
  integer :: NEL(6),NID(4,3,2),IE,ID,NX,NX1,NY,NZ,IEL

  DATA NID/1,2,3,6, 1,4,5,6, 1,5,2,6, 2,3,1,6, 2,4,5,6, 2,1,4,6/
  
  if(ID.eq.0) then
     if(NX+NX1.ge.(NXA(NY)+1)/2+(NXA(NY+1)+1)/2) then
        if(NY.ge.(NYMAX+1)/2) then
           NEL(1)=NDA(NX   ,NY  ,NZ  )
           NEL(2)=NDA(NX+1 ,NY  ,NZ  )
           NEL(3)=NDA(NX1  ,NY+1,NZ  )
           NEL(4)=NDA(NX   ,NY  ,NZ+1)
           NEL(5)=NDA(NX+1 ,NY  ,NZ+1)
           NEL(6)=NDA(NX1  ,NY+1,NZ+1)
           IEL=1
        else
           NEL(1)=NDA(NX1  ,NY+1,NZ  )
           NEL(2)=NDA(NX   ,NY  ,NZ  )
           NEL(3)=NDA(NX+1 ,NY  ,NZ  )
           NEL(4)=NDA(NX1  ,NY+1,NZ+1)
           NEL(5)=NDA(NX   ,NY  ,NZ+1)
           NEL(6)=NDA(NX+1 ,NY  ,NZ+1)
           IEL=1
        end if
     else
        if(NY.ge.(NYMAX+1)/2) then
           NEL(1)=NDA(NX   ,NY  ,NZ  )
           NEL(2)=NDA(NX+1 ,NY  ,NZ  )
           NEL(3)=NDA(NX1  ,NY+1,NZ  )
           NEL(4)=NDA(NX   ,NY  ,NZ+1)
           NEL(5)=NDA(NX+1 ,NY  ,NZ+1)
           NEL(6)=NDA(NX1  ,NY+1,NZ+1)
           IEL=2
        else
           NEL(1)=NDA(NX+1 ,NY  ,NZ  )
           NEL(2)=NDA(NX1  ,NY+1,NZ  )
           NEL(3)=NDA(NX   ,NY  ,NZ  )
           NEL(4)=NDA(NX+1 ,NY  ,NZ+1)
           NEL(5)=NDA(NX1  ,NY+1,NZ+1)
           NEL(6)=NDA(NX   ,NY  ,NZ+1)
           IEL=2
        end if
     end if
  else
     if(NX+NX1.ge.(NXA(NY)+1)/2+(NXA(NY+1)+1)/2) then
        if(NY.ge.(NYMAX+1)/2) then
           NEL(1)=NDA(NX1  ,NY+1,NZ  )
           NEL(2)=NDA(NX   ,NY  ,NZ  )
           NEL(3)=NDA(NX1+1,NY+1,NZ  )
           NEL(4)=NDA(NX1  ,NY+1,NZ+1)
           NEL(5)=NDA(NX   ,NY  ,NZ+1)
           NEL(6)=NDA(NX1+1,NY+1,NZ+1)
           IEL=2
        else
           NEL(1)=NDA(NX1+1,NY+1,NZ  )
           NEL(2)=NDA(NX1  ,NY+1,NZ  )
           NEL(3)=NDA(NX   ,NY  ,NZ  )
           NEL(4)=NDA(NX1+1,NY+1,NZ+1)
           NEL(5)=NDA(NX1  ,NY+1,NZ+1)
           NEL(6)=NDA(NX   ,NY  ,NZ+1)
           IEL=2
        end if
     else
        if(NY.ge.(NYMAX+1)/2) then
           NEL(1)=NDA(NX   ,NY  ,NZ  )
           NEL(2)=NDA(NX1+1,NY+1,NZ  )
           NEL(3)=NDA(NX1  ,NY+1,NZ  )
           NEL(4)=NDA(NX   ,NY  ,NZ+1)
           NEL(5)=NDA(NX1+1,NY+1,NZ+1)
           NEL(6)=NDA(NX1  ,NY+1,NZ+1)
           IEL=1
        else
           NEL(1)=NDA(NX1+1,NY+1,NZ  )
           NEL(2)=NDA(NX1  ,NY+1,NZ  )
           NEL(3)=NDA(NX   ,NY  ,NZ  )
           NEL(4)=NDA(NX1+1,NY+1,NZ+1)
           NEL(5)=NDA(NX1  ,NY+1,NZ+1)
           NEL(6)=NDA(NX   ,NY  ,NZ+1)
           IEL=1
        end if
     end if
  end if
  IE=IE+1
  NDELM(1,IE)=NEL(NID(1,1,IEL))
  NDELM(2,IE)=NEL(NID(2,1,IEL))
  NDELM(3,IE)=NEL(NID(3,1,IEL))
  NDELM(4,IE)=NEL(NID(4,1,IEL))
  IE=IE+1
  NDELM(1,IE)=NEL(NID(1,2,IEL))
  NDELM(2,IE)=NEL(NID(2,2,IEL))
  NDELM(3,IE)=NEL(NID(3,2,IEL))
  NDELM(4,IE)=NEL(NID(4,2,IEL))
  IE=IE+1
  NDELM(1,IE)=NEL(NID(1,3,IEL))
  NDELM(2,IE)=NEL(NID(2,3,IEL))
  NDELM(3,IE)=NEL(NID(3,3,IEL))
  NDELM(4,IE)=NEL(NID(4,3,IEL))

  return
end subroutine SETELL

!     ****** Set Element Array ******

subroutine SETELMX(IERR)

  use wfcomm
  implicit none
  integer :: NY,NX0MAX,NX1MAX,NX10,NX11,NZ,IERR,NX00,NX01
  integer :: IE,NYDO
  real(8) :: RL00,RL01

  IE=0

  do NYDO=1,NYMAX-1
     NY=NYDO
     NX0MAX=NXA(NY)
     NX1MAX=NXA(NY+1)
     
     if(NX0MAX.eq.1) then
        do NX10=1,NX1MAX
           NX11=NX10+1
           if(NX11.gt.NX1MAX) NX11=1
           do NZ=1,NZMAX-1
              IE=IE+3
           end do
        end do
     else

        NX00=1
        NX10=1
        
100     NX01=NX00+1
        if(NX01.gt.NX0MAX) NX01=1
        NX11=NX10+1
        if(NX11.gt.NX1MAX) NX11=1
        
        RL00=(XNDA(NX11,NY+1)-XNDA(NX00,NY))**2&
             &          +(YNDA(NX11,NY+1)-YNDA(NX00,NY))**2
        RL01=(XNDA(NX10,NY+1)-XNDA(NX01,NY))**2&
             &          +(YNDA(NX10,NY+1)-YNDA(NX01,NY))**2
        
        if(RL00.lt.RL01) then
           do NZ=1,NZMAX-1
              IE=IE+3
           end do
           if(NX10.eq.NX1MAX) then
101           do NZ=1,NZMAX-1
                 IE=IE+3
              end do
              ! ---- add by Y.Kubota at Nov./14/2019 ----
	       if(NX00.ne.NX0MAX)then
                  NX00=NX00+1
                  goto 101
               endif
               ! ---- Nov./14/2019 ----
               goto 200
           else
              NX10=NX10+1
           end if
        else if(RL00.ge.RL01) then
           do NZ=1,NZMAX-1
              IE=IE+3
           end do
           if(NX00.eq.NX0MAX) then
              do NZ=1,NZMAX-1
                 IE=IE+3
              end do
              goto 200
           else
              NX00=NX00+1
           end if
        end if
        goto 100
        
200     continue
     end if

  end do
  
  NEMAX=IE
  
  call wfelm_allocate

  IE=0

  do NYDO=1,NYMAX-1
     NY=NYDO
     NX0MAX=NXA(NY)
     NX1MAX=NXA(NY+1)
     
     if(NX0MAX.eq.1) then
        do NX10=1,NX1MAX
           NX11=NX10+1
           if(NX11.gt.NX1MAX) NX11=1
           do NZ=1,NZMAX-1
              call SETELLX(1,NX10,NX11,NY,NZ,IE,0)
           end do
        end do
     else
        
        NX00=1
        NX10=1
        
110     NX01=NX00+1
        if(NX01.gt.NX0MAX) NX01=1
        NX11=NX10+1
        if(NX11.gt.NX1MAX) NX11=1
        
        RL00=(XNDA(NX11,NY+1)-XNDA(NX00,NY))**2&
            +(YNDA(NX11,NY+1)-YNDA(NX00,NY))**2
        RL01=(XNDA(NX10,NY+1)-XNDA(NX01,NY))**2&
            +(YNDA(NX10,NY+1)-YNDA(NX01,NY))**2
        
        if(RL00.lt.RL01) then
           do NZ=1,NZMAX-1
              call SETELLX(NX00,NX10,NX11,NY,NZ,IE,0)
           end do
           if(NX10.eq.NX1MAX) then
              do NZ=1,NZMAX-1
                 call SETELLX(NX00,NX01,NX11,NY,NZ,IE,1)
              end do
              goto 210
           else
              NX10=NX10+1
           end if
        else if(RL00.ge.RL01) then
           do NZ=1,NZMAX-1
              call SETELLX(NX00,NX01,NX10,NY,NZ,IE,1)
           end do
           if(NX00.eq.NX0MAX) then
111           do NZ=1,NZMAX-1
                 call SETELLX(NX01,NX10,NX11,NY,NZ,IE,0)
              end do
               ! ---- add by Y.Kubota at Nov./14/2019 ----
	       if(NX00.ne.NX0MAX)then
                  NX00=NX00+1
                  NX01=NX00+1
                  if(NX01.gt.NX0MAX) NX01=1
                  goto 111
	       endif
               ! ---- add by Y.Kubota at Nov./14/2019 ----
              goto 210
           else
              NX00=NX00+1
           end if
        end if

        GOTO 110

210     CONTINUE
     ENDIF
     
  ENDDO

  NEMAX=IE
  IERR=0
  RETURN
end subroutine SETELMX

!     ****** Set Element Array SUB******

subroutine SETELLX(NX0,NX1,NX2,NY,NZ,IE,ID)
  
  use wfcomm
  implicit none
  integer :: NEL(6),NID(4,3,4),IE,ID,NX0,NY,NZ,NX1,NX2,IEL
  real(8) :: XC,YC
  DATA NID/1,2,3,6, 1,4,5,6, 1,5,2,6, 1,2,3,5, 1,6,4,5, 1,3,6,5,&
           1,2,3,6, 2,4,5,6, 1,4,2,6, 1,2,3,6, 1,4,5,6, 1,5,2,6/
  
  if(ID.eq.0) then
     NEL(1)=NDA(NX0,NY  ,NZ  )
     NEL(2)=NDA(NX1,NY+1,NZ  )
     NEL(3)=NDA(NX2,NY+1,NZ  )
     NEL(4)=NDA(NX0,NY  ,NZ+1)
     NEL(5)=NDA(NX1,NY+1,NZ+1)
     NEL(6)=NDA(NX2,NY+1,NZ+1)
     IEL=1
  else
     NEL(1)=NDA(NX1,NY  ,NZ  )
     NEL(2)=NDA(NX0,NY  ,NZ  )
     NEL(3)=NDA(NX2,NY+1,NZ  )
     NEL(4)=NDA(NX1,NY  ,NZ+1)
     NEL(5)=NDA(NX0,NY  ,NZ+1)
     NEL(6)=NDA(NX2,NY+1,NZ+1)
     IEL=3
  end if
  
  XC=XND(NEL(1))+XND(NEL(2))+XND(NEL(3))
  YC=YND(NEL(1))+YND(NEL(2))+YND(NEL(3))
  if(XC*YC.lt.0.D0) IEL=IEL+1
  
  IE=IE+1
  NDELM(1,IE)=NEL(NID(1,1,IEL))
  NDELM(2,IE)=NEL(NID(2,1,IEL))
  NDELM(3,IE)=NEL(NID(3,1,IEL))
  NDELM(4,IE)=NEL(NID(4,1,IEL))
  IE=IE+1
  NDELM(1,IE)=NEL(NID(1,2,IEL))
  NDELM(2,IE)=NEL(NID(2,2,IEL))
  NDELM(3,IE)=NEL(NID(3,2,IEL))
  NDELM(4,IE)=NEL(NID(4,2,IEL))
  IE=IE+1
  NDELM(1,IE)=NEL(NID(1,3,IEL))
  NDELM(2,IE)=NEL(NID(2,3,IEL))
  NDELM(3,IE)=NEL(NID(3,3,IEL))
  NDELM(4,IE)=NEL(NID(4,3,IEL))
  
  return
end subroutine SETELLX

!     ****** List Element Data ******

subroutine WFLDIV

  use wfcomm
  implicit none
  integer :: I,J

  write(6,100) NNMAX,(I,XND(I),YND(I),ZND(I),I=1,NNMAX)
100 format(/' ','NODE DATA',7X,'NNMAX=',I6/&
            ' ',8X,'XND',10X,'YND',10X,'ZND'/&
           (' ',I6,1P3E13.4))
     
  write(6,200) NEMAX,(I,(NDELM(J,I),J=1,4),&
                      I=1,NEMAX)
200 format(/' ','ELEMENT DATA',5X,'NEMAX=',I6/&
           (' ',2(I5,'(',4I5,')',2X)))
     
  write(6,400) NBMAX,(NDBDY(I),I=1,NBMAX)
400 format(' ','BOUNDARY DATA',4X,'NBMAX=',I6/&
          (' ',10(I6)))
     
  write(6,500) MBND,MLEN,(I,KANOD(I),IMLEN(I),I=1,NSDMAX)
500 format(' ','NODE AND MATRIX',2X,'MBND=',I6,4X,'MLEN=',I4/&
          (' ',4(I5,'(',2I5,')')))
  return
end subroutine WFLDIV
!------------------------------------------------------------
subroutine wfdiv_broadcast

  use libmpi
  use libmtx
  use wfcomm
  implicit none

  real(8),dimension(18) :: rdata
  integer,dimension(1)  :: idata 
  integer :: NDATA

  if (nrank.eq.0) then
     rdata(1)=BXMIN
     rdata(2)=BXMAX
     rdata(3)=BYMIN
     rdata(4)=BYMAX
     rdata(5)=BZMIN
     rdata(6)=BZMAX
     rdata(7)=DELX
     rdata(8)=DELY
     rdata(9)=DELZ
     SELECT CASE(iddiv)
     CASE(1)
        rdata(10)=RB
     CASE(2)
        rdata(10)=RB
        rdata(11)=RBAX
     CASE(3)
        rdata(10)=RBIN
        rdata(11)=RBOUT
        rdata(12)=BZMIDL
        rdata(13)=BZMIDH
        rdata(14)=DXIN
        rdata(15)=DYIN
        rdata(16)=DXOUT
        rdata(17)=DYOUT
        rdata(18)=DELZM
     END SELECT
     idata(1)=IDDIV
  end if
  
  call mtx_barrier
  call mtx_broadcast_integer(idata,1)
  IDDIV=idata(1)
  SELECT CASE(iddiv)
  CASE(0)
     NDATA=9
  CASE(1)
     NDATA=10
  CASE(2)
     NDATA=11
  CASE(3)
     NDATA=18
  END SELECT
  
  call mtx_barrier
  call mtx_broadcast_real8(rdata,NDATA)
  BXMIN=rdata(1)
  BXMAX=rdata(2)
  BYMIN=rdata(3)
  BYMAX=rdata(4)
  BZMIN=rdata(5)
  BZMAX=rdata(6)
  DELX =rdata(7)
  DELY =rdata(8)
  DELZ =rdata(9)
  SELECT CASE(iddiv)
  CASE(1)
     RB=rdata(10)
  CASE(2)
     RB=rdata(10)
     RBAX=rdata(11)
  CASE(3)
     RBIN=rdata(10)
     RBOUT=rdata(11)
     BZMIDL=rdata(12)
     BZMIDH=rdata(13)
     DXIN=rdata(14)
     DYIN=rdata(15)
     DXOUT=rdata(16)
     DYOUT=rdata(17)
     DELZM=rdata(18)
  END SELECT

  return 
end subroutine wfdiv_broadcast
