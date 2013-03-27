!
!  ----  ADD. BY YAMA  09/SEP./2010  ----
!
! -- 現在，アンカー部に設置されている平板状のアンテナを模擬する．
! -- アンテナがプラズマに対向する部分の長さは，磁力管断面の長軸
! -- あるいは短軸の長さに合わせる．
! -- 磁力管は，Rcc=18cm を通過するものを基準とする．
!
      SUBROUTINE DEFBTA(NA)
!
      use wfcomm
      implicit none
      INTEGER,INTENT(IN):: NA
      REAL(8):: ENDY,DFZMIN,ZPT,DFZ,XPT,YPT,ZTEMP
      INTEGER:: J

! -- INIT. --
      ZANT   = 5.195D0
      RD     = 0.2D0
      ENDY   = 1.D0

  10  WRITE(6,20) ZANT,RD,ENDY
  20  FORMAT(' ','## Z-POS = ',F10.3/&
            ' ','## RD    = ',F10.3/&
            ' ','## END-Y = ',F10.3)
      READ(5,*,ERR=10,END=90) ZANT,RD,ENDY

!     与えられたZ座標に，最も近いZ座標を持つ磁力管断面を探す．
!     FLZ,FLX,FLYは cm 単位
!     ZPT,XPT,YPTは  m 単位
      DFZMIN = 1.D2
      DO J=1,NGFLIN
         ZPT = FLZ(J)/1.D2
         DFZ = ABS(ZANT-ZPT)
         IF(DFZ.LT.DFZMIN) THEN
            DFZMIN = DFZ
            XPT = FLX(J)/1.D2
            YPT = FLY(J)/1.D2
            ZTEMP = ZPT
         ENDIF
      ENDDO

! -- 給電点の作成

      XJ0(1,NA) = RD/RA*XPT
      YJ0(1,NA) = ENDY
      ZJ0(1,NA) = ZANT

! -- プラズマ対向部の作成

      IF(ENDY.GT.0.D0) THEN
         XJ0(2,NA) =  XJ0(1,NA)
         YJ0(2,NA) =  RD/RA*YPT
         ZJ0(2,NA) =  ZJ0(1,NA)
         XJ0(3,NA) = -RD/RA*XPT
         YJ0(3,NA) =  YJ0(2,NA)
         ZJ0(3,NA) =  ZJ0(2,NA)
      ELSE
         XJ0(2,NA) =  XJ0(1,NA)
         YJ0(2,NA) = -RD/RA*YPT
         ZJ0(2,NA) =  ZJ0(1,NA)
         XJ0(3,NA) = -RD/RA*XPT
         YJ0(3,NA) =  YJ0(2,NA)
         ZJ0(3,NA) =  ZJ0(2,NA)
      ENDIF

! -- 接地点の作成

      XJ0(4,NA) = XJ0(3,NA)
      YJ0(4,NA) = ENDY
      ZJ0(4,NA) = ZJ0(3,NA)

! -- アンテナを構成する点の数

      JNUM0(NA)=4

  90  RETURN
      END SUBROUTINE DEFBTA


!  ----  ADD. BY YAMA  18/FEB./2010  ----
!
! -- 磁力管と平面の交点が形成する曲線に沿って，電流路を作る．
! -- 磁力管は，Rcc=18cm を通過するものを基準とする．
! -- アンテナの半径は，RD / RA * (FLX,FLY) で調整する．
! -- (RA = 0.18m)

      SUBROUTINE MBDEFA(NA)

      use wfcomm
      implicit none
      INTEGER,INTENT(IN):: NA
      REAL(8):: DEGN,ENDR,RTAGL,PLA,PLB,PLC,PLD,DTHETA,THETA
      REAL(8):: X1,Y1,Z1,X2,Y2,Z2,SX,SY,SZ,DNM,SML
      INTEGER:: NJ,J

      DEGN   = PI/180.D0
      ZANT   = 5.195D0
      RTDEG  = 0.D0
      ENDR   = 0.5D0
      THETJ1 = 10.D0
      THETJ2 = 170.D0
      RD     = 0.2D0
      NJMAX  = 40
 1000 WRITE(6,1100) ZANT,RTDEG,ENDR,THETJ1,THETJ2,RD,NJMAX
 1100 FORMAT(' ','## ROTATION ANGLE = ',F10.3/&
            ' ','## Z-COORDINATE   = ',F10.3/&
            ' ','## END-X          = ',F10.3/&
            ' ','## THETJ1, THETJ2 = ',2F10.3/&
            ' ','## RD, NJMAX      = ',F10.3,I5)
       READ(5,*,ERR=1000,END=9000) ZANT,RTDEG,ENDR,&
                                  THETJ1,THETJ2,RD,NJMAX

! -- アンテナ線に接する平面の方程式
! -- (X,Y,Z) = (X,0,ZANT) を通過する，X軸に平行な平面の方程式を作成
! -- Y軸と平面の成す角: RTDEG [deg.] or RTAGL [rad.]
! -- NORMAL VECTOR : (PLA,PLB,PLC)
! -- EQ. OF PLANE  : PLA*X + PLB*Y + PLC*Z + PLD = 0.D0

      RTAGL =  DEGN*RTDEG
      PLA   =  0.D0
      PLB   = -SIN(RTAGL)
      PLC   =  COS(RTAGL)
      PLD   = -1.D0*ZANT*COS(RTAGL)

!-- 右手系のXY座標平面内，X軸を基準とする右回り方向の角度: THETA
! -- 角度 THETA に対応する磁力管上の曲線(磁力線)と，平面との交点を探索
! -- 線分 (X1,Y1,Z2) -> (X2,Y2,Z2) が，平面と交差するかどうかを判定

      DTHETA=(THETJ2-THETJ1)/(NJMAX-3)
      DO NJ=2,NJMAX-1
         THETA = DEGN*(DTHETA*(NJ-2)+THETJ1)
!        磁力管断面を楕円と仮定
!        楕円とX軸，Y軸との交点 : FLX, FLY
!        FLZ,FLX,FLYは "cm" 単位 -> "m" 単位へ変更 : /1.D2
!        磁力管断面の半径を調整 : *RD/RA

         DO J=2,NGFLIN
            X1 = FLX(J-1)/1.D2 * RD/RA*COS(THETA)
            Y1 = FLY(J-1)/1.D2 * RD/RA*SIN(THETA)
            Z1 = FLZ(J-1)/1.D2
            X2 = FLX(J)  /1.D2 * RD/RA*COS(THETA)
            Y2 = FLY(J)  /1.D2 * RD/RA*SIN(THETA)
            Z2 = FLZ(J)  /1.D2
            SX = X2 - X1
            SY = Y2 - Y1
            SZ = Z2 - Z1
! -- 点(X1,Y1,Z2)を基点とするベクトル(SX,SY,SZ) が，
! -- 平面と交差する為の条件は，[ 0.D0 .GE. SML .GE. 1.D0 ]
! -- DNM が零の時，ベクトル(SX,SY,SZ)と平面は平行している

            DNM = PLA*SX + PLB*SY + PLC*SZ
            IF(DNM.NE.0.D0) THEN
               SML = -1.D0*(PLD+PLA*X1+PLB*Y2+PLC*Z1) / DNM
               IF((SML.GE.0.D0).AND.(SML.LE.1.D0)) THEN
                  XJ0(NJ,NA) = X1 + SML*SX
                  YJ0(NJ,NA) = Y1 + SML*SY
                  ZJ0(NJ,NA) = Z1 + SML*SZ
                  GOTO 2000
               ENDIF
            ENDIF
         ENDDO
         WRITE(6,*) 'XX ERR. : ',NA
 2000    CONTINUE
      ENDDO
! -- 給電点の作成

      XJ0(1,NA)     = ENDR
      YJ0(1,NA)     = YJ0(2,NA)
      ZJ0(1,NA)     = ZJ0(2,NA)
! -- 接地点の作成

      XJ0(NJMAX,NA) = -1.D0*ENDR
      YJ0(NJMAX,NA) = YJ0(NJMAX-1,NA)
      ZJ0(NJMAX,NA) = ZJ0(NJMAX-1,NA)
! -- アンテナを構成する点の数

      JNUM0(NA)=NJMAX

 9000 RETURN
      END SUBROUTINE MBDEFA


!  ----  ADD. BY YAMA  16/FEB./2010  ----

      SUBROUTINE MBDEFA2

      use wfcomm
      implicit none
      REAL(8):: ZCRD(1000),ACRD(1000,1000,3)
      INTEGER(8):: NTPT(1000)
      REAL(8):: ANTP(100,3),ANTN(100,3)
      REAL(8):: AUP(3),ALP(3),AUN(3),ALN(3)
      REAL(8):: DSRTK,DEGN,ENDR,RDU,RDL,THRSD,RTAGL,PLA,PLB,PLC,PLD
      REAL(8):: XCNT,YCNT,ZCNT,TMPZ,DFZMIN,ZTMP,DFZ,XPTB,YPTB,ZPTB
      REAL(8):: DSTE,ANGL,RXU,RYU,XBU,YBU,DSTU,RXLRYL,XBL,YBL,DSTL
      REAL(8):: AMAXP,AMINP,AMAXN,AMINN,SRX,SRY,SRZ,DSTA
      REAL(8):: ANTPX,ANTPY,ANTPZ,ANTNX,ANTNY,ANTNZ,RXL,RYL,DSRTJ
      INTEGER:: NZPT,IE,IND1,IND2,IND3,IND4,J,K,NN,JJ,KK,IDAP,IDAN,NDTP,NDTN
      INTEGER:: NA,NJ

      DEGN  = PI/180.D0
      ENDR  = 0.5D0
      RDU   = 0.2D0
      RDL   = 0.18D0
      RTDEG = 0.D0
      THRSD = 1.D-2
  100 WRITE(6,200) ENDR,RDU,RDL,RTDEG,THRSD
  200 FORMAT(' ','## ENDX  = ',F10.3/&
            ' ','## RDU   = ',F10.3/&
            ' ','## RDL   = ',F10.3/&
            ' ','## THETA = ',F10.3/&
            ' ','## THRSD = ',F10.3)
      READ(5,*,ERR=100,END=100) ENDR,RDU,RDL,RTDEG,THRSD

! -- アンテナ線に接する平面の方程式
! -- NORMAL VECTOR : (PLA,PLB,PLC)
! -- EQ. OF PLANE  : PLA*X + PLB*Y + PLC*Z + PLD = 0.D0
      RTAGL = DEGN*RTDEG
      PLA   =  0.D0
      PLB   = -SIN(RTAGL)
      PLC   =  COS(RTAGL)
      PLD   = -5.2D0*COS(RTAGL)

! -- アンテナ作成領域内のZ座標配列 [ ZCRD(NZPT) ] を作成
      NZPT = 0
      DO IE=1,NEMAX
! -- 要素の位置は，重心 [ XCNT, YCNT, ZCNT ] で代表
! -- 要素を構成する節点の座標を抽出
         IND1=NDELM(1,IE)
         IND2=NDELM(2,IE)
         IND3=NDELM(3,IE)
         IND4=NDELM(4,IE)

! -- 重心座標の計算
         XCNT=(XND(IND1)+XND(IND2)+XND(IND3)+XND(IND4))/4.D0
         YCNT=(YND(IND1)+YND(IND2)+YND(IND3)+YND(IND4))/4.D0
         ZCNT=(ZND(IND1)+ZND(IND2)+ZND(IND3)+ZND(IND4))/4.D0

! -- アンテナ作成領域 (Z軸座標に対する要請)
         IF((ZCNT.GE.5.0D0).AND.(ZCNT.LE.5.4D0).AND.&
           (XCNT.LE.1.D-1).AND.(YCNT.LE.1.D-1)) THEN
            IF(NZPT.EQ.0) THEN
               NZPT       = NZPT+1
               ZCRD(NZPT) = ZCNT
            ELSE
! -- 重複するZ座標を除外しつつ，Z座標配列を構築する
               DO J=1,NZPT
                  IF(ZCRD(J).EQ.ZCNT) GOTO 300
               ENDDO
               NZPT       = NZPT+1
               ZCRD(NZPT) = ZCNT
            ENDIF
         ENDIF
  300    CONTINUE
      ENDDO

! -- ZCRDを昇順で整列 (確認の為) --
      DO J=1,NZPT-1
      DO K=J+1,NZPT
         IF(ZCRD(K).LT.ZCRD(J)) THEN
            TMPZ    = ZCRD(J)
            ZCRD(J) = ZCRD(K)
            ZCRD(K) = TMPZ
         ENDIF
      ENDDO
      ENDDO

! -- 同じZ座標を持つ要素群から，
! -- Rcc=18cm 磁力管の外側表面上に存在する要素を抽出
      DO NN=1,NZPT-1
! -- Rcc=18cm の磁力管断面を抽出 (楕円: XPTB*X + YPTB*Y = 1)
         DFZMIN = 1.D2
         DO J=1,NGFLIN
            ZTMP = FLZ(J)/1.D2
            DFZ  = DABS(ZCRD(NN)-ZTMP)
            IF(DFZ.LT.DFZMIN) THEN
               DFZMIN = DFZ
               XPTB   = FLX(J)/1.D2
               YPTB   = FLY(J)/1.D2
               ZPTB   = ZTMP
            ENDIF
         ENDDO

!         write(6,*) 'Z,X,Y= ',ZPTB, XPTB, YPTB

         NTPT(NN)=0
         DO IE=1,NEMAX
! -- 要素位置を重心で代表
            IND1=NDELM(1,IE)
            IND2=NDELM(2,IE)
            IND3=NDELM(3,IE)
            IND4=NDELM(4,IE)

            XCNT=(XND(IND1)+XND(IND2)+XND(IND3)+XND(IND4))/4.D0
            YCNT=(YND(IND1)+YND(IND2)+YND(IND3)+YND(IND4))/4.D0
            ZCNT=(ZND(IND1)+ZND(IND2)+ZND(IND3)+ZND(IND4))/4.D0

! -- 要素位置の原点からの距離 [DSTE]
            DSTE=DSQRT(XCNT*XCNT + YCNT*YCNT)

! -- 要素のZ座標に対する要請
!            IF((ZCNT.GT.ZPTL).AND.(ZCNT.LE.ZPTH)) THEN
            IF(ZCNT.EQ.ZCRD(NN)) THEN
! -- X軸と要素の重心ベクトルとの成す角 [ANGL]
! -- 座標(XCNT,YCNT)が，どの象限に属するかで場合分け
               IF(XCNT.GT.0.D0) THEN
                  ANGL = ATAN(YCNT/XCNT)
               ELSEIF(XCNT.LT.0.D0) THEN
                  ANGL = ATAN(YCNT/XCNT) + PI
               ELSEIF(YCNT.GT.0.D0) THEN
                  ANGL = PI/2.D0
               ELSEIF(YCNT.LT.0.D0) THEN
                  ANGL = -1.D0*PI/2.D0
               ELSE
                  ANGL = 0.D0
               ENDIF
! -- ある角ANGLに対応する磁力管上の座標とZ軸との距離 [DSTB]
               RXU  = RDU/RA * XPTB
               RYU  = RDU/RA * YPTB
               XBU  = RXU*COS(ANGL)
               YBU  = RYU*SIN(ANGL)
               DSTU = DSQRT(XBU*XBU + YBU*YBU)
               
               RXL  = RDL/RA * XPTB
               RYL  = RDL/RA * YPTB
               XBL  = RXL*COS(ANGL)
               YBL  = RYL*SIN(ANGL)
               DSTL = DSQRT(XBL*XBL + YBL*YBL)

! -- 磁力管の外側にあり，磁力間との距離が最も小さい要素の探索
!    NN  : Z-INDEX
!    NTPT: THETA-INDEX
!               IF((EEBB.GT.0.D0).AND.(EEBB.LT.DMIN)) THEN
!               IF(EEBB.LT.DMIN) THEN
               IF((DSTE.GT.DSTL).AND.(DSTE.LT.DSTU)) THEN



                  NTPT(NN)            = NTPT(NN) + 1
                  ACRD(NN,NTPT(NN),1) = XCNT
                  ACRD(NN,NTPT(NN),2) = YCNT
                  ACRD(NN,NTPT(NN),3) = ZCNT

               ENDIF
            ENDIF
         ENDDO
            write(6,8000) NN,NTPT(NN)
 8000       format(2I5)
      ENDDO

! -- 平面と接するか，あるいは平面の近傍にある点を探索
! -- ...P: Y座標が正
! -- ...N: Y座標が負
      IDAP  = 0
      AMAXP = 0.D0
      AMINP = 0.D0
      IDAN  = 0
      AMAXN = 0.D0
      AMINN = 0.D0
      DO JJ=1,NZPT-1
      DO KK=1,NTPT(JJ)
! -- 磁力管周辺の要素(SRX,SRY,SRZ)と，平面との距離 [DSTA]
         SRX  = ACRD(JJ,KK,1)
         SRY  = ACRD(JJ,KK,2)
         SRZ  = ACRD(JJ,KK,3)
         DSTA = DABS(PLA*SRX + PLB*SRY + PLC*SRZ + PLD)/&
               SQRT(PLA*PLA + PLB*PLB + PLC*PLC)
! -- Y座標の正負で場合分けし，二本のアンテナを作成
         IF(SRY.GE.0.D0) THEN
! -- Y座標: 正
! -- X座標の正負で場合分け(平面との交点は，+X ,-X の二箇所に存在する)
            IF(SRX.GE.0.D0) THEN
               IF(DSTA.LT.THRSD) THEN
                  IDAP         = IDAP+1
                  NDTP         = IDAP
                  ANTP(IDAP,1) = SRX
                  ANTP(IDAP,2) = SRY
                  ANTP(IDAP,3) = SRZ
               write(6,*) JJ,KK,DSTA,IDAP
! -- 上端の座標
                  IF(SRX.GT.AMAXP) THEN
                     AUP(1) = SRX
                     AUP(2) = SRY
                     AUP(3) = SRZ
                     AMAXP  = SRX
                  ENDIF
               ENDIF
            ELSE
               IF(DSTA.LT.THRSD) THEN
                  IDAP         = IDAP+1
                  NDTP         = IDAP
                  ANTP(IDAP,1) = SRX
                  ANTP(IDAP,2) = SRY
                  ANTP(IDAP,3) = SRZ
               write(6,*) JJ,KK,DSTA,IDAP
! -- 下端の座標
                  IF(SRX.LT.AMINP) THEN
                     ALP(1) = SRX
                     ALP(2) = SRY
                     ALP(3) = SRZ
                     AMINP  = SRX
                  ENDIF
               ENDIF
            ENDIF
!            IDAP = IDAP + 2
         ELSE
! -- Y座標: 負
! -- X座標の正負で場合分け(平面との交点は，+X ,-X の二箇所に存在する)
            IF(SRX.GE.0.D0) THEN
               IF(DSTA.LT.THRSD) THEN
                  IDAN         = IDAN+1
                  NDTN         = IDAN
                  ANTN(IDAN,1) = SRX
                  ANTN(IDAN,2) = SRY
                  ANTN(IDAN,3) = SRZ
! -- 上端の座標
                  IF(SRX.GT.AMAXN) THEN
                     AUN(1) = SRX
                     AUN(2) = SRY
                     AUN(3) = SRZ
                     AMAXN  = SRX
                  ENDIF
               ENDIF
            ELSE
               IF(DSTA.LT.THRSD) THEN
                  IDAN         = IDAN+1
                  NDTN         = IDAN
                  ANTN(IDAN,1) = SRX
                  ANTN(IDAN,2) = SRY
                  ANTN(IDAN,3) = SRZ
! -- 下端の座標
                  IF(SRX.LT.AMINN) THEN
                     ALN(1) = SRX
                     ALN(2) = SRY
                     ALN(3) = SRZ
                     AMINN  = SRX
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO

! -- SORTING RECORDS --
      write(6,*) 'NDTP=',NDTP

      DO J=1,NDTP-1
      DO K=J+1,NDTP
         IF(J.EQ.1) THEN
            DSRTJ =   (ANTP(J,1)-AUP(1))*(ANTP(J,1)-AUP(1))&
                   + (ANTP(J,2)-AUP(2))*(ANTP(J,2)-AUP(2))&
                   + (ANTP(J,3)-AUP(3))*(ANTP(J,3)-AUP(3))
            DSRTK =   (ANTP(K,1)-AUP(1))*(ANTP(K,1)-AUP(1))&
                   + (ANTP(K,2)-AUP(2))*(ANTP(K,2)-AUP(2))&
                   + (ANTP(K,3)-AUP(3))*(ANTP(K,3)-AUP(3))
         ELSE
            DSRTJ =   (ANTP(J,1)-ANTP(J-1,1))*(ANTP(J,1)-ANTP(J-1,1))&
                   + (ANTP(J,2)-ANTP(J-1,2))*(ANTP(J,2)-ANTP(J-1,2))&
                   + (ANTP(J,3)-ANTP(J-1,3))*(ANTP(J,3)-ANTP(J-1,3))
            DSRTK =   (ANTP(K,1)-ANTP(J-1,1))*(ANTP(K,1)-ANTP(J-1,1))&
                   + (ANTP(K,2)-ANTP(J-1,2))*(ANTP(K,2)-ANTP(J-1,2))&
                   + (ANTP(K,3)-ANTP(J-1,3))*(ANTP(K,3)-ANTP(J-1,3))
         ENDIF

         IF(DSRTK.LT.DSRTJ) THEN
            ANTPX     = ANTP(J,1)
            ANTPY     = ANTP(J,2)
            ANTPZ     = ANTP(J,3)
            ANTP(J,1) = ANTP(K,1)
            ANTP(J,2) = ANTP(K,2)
            ANTP(J,3) = ANTP(K,3)
            ANTP(K,1) = ANTPX
            ANTP(K,2) = ANTPY
            ANTP(K,3) = ANTPZ
         ENDIF
         IF((ANTP(J-1,1).EQ.ALP(1)).AND.&
           (ANTP(J-1,2).EQ.ALP(2)).AND.&
           (ANTP(J-1,3).EQ.ALP(3))) THEN
               NDTP = J
               GOTO 1000
         ENDIF
      ENDDO
      ENDDO
 1000 CONTINUE

      write(6,*) 'NDTN=',NDTN

      DO J=1,NDTN-1
      DO K=J+1,NDTN
         IF(J.EQ.1) THEN
            DSRTJ =   (ANTN(J,1)-AUN(1))*(ANTN(J,1)-AUN(1))&
     &              + (ANTN(J,2)-AUN(2))*(ANTN(J,2)-AUN(2))&
     &              + (ANTN(J,3)-AUN(3))*(ANTN(J,3)-AUN(3))
            DSRTK =   (ANTN(K,1)-AUN(1))*(ANTN(K,1)-AUN(1))&
     &              + (ANTN(K,2)-AUN(2))*(ANTN(K,2)-AUN(2))&
     &              + (ANTN(K,3)-AUN(3))*(ANTN(K,3)-AUN(3))
         ELSE
            DSRTJ =   (ANTN(J,1)-ANTN(J-1,1))*(ANTN(J,1)-ANTN(J-1,1))&
     &              + (ANTN(J,2)-ANTN(J-1,2))*(ANTN(J,2)-ANTN(J-1,2))&
     &              + (ANTN(J,3)-ANTN(J-1,3))*(ANTN(J,3)-ANTN(J-1,3))
            DSRTK =   (ANTN(K,1)-ANTN(J-1,1))*(ANTN(K,1)-ANTN(J-1,1))&
     &              + (ANTN(K,2)-ANTN(J-1,2))*(ANTN(K,2)-ANTN(J-1,2))&
     &              + (ANTN(K,3)-ANTN(J-1,3))*(ANTN(K,3)-ANTN(J-1,3))
         ENDIF

         IF(DSRTK.LT.DSRTJ) THEN
            ANTNX     = ANTN(J,1)
            ANTNY     = ANTN(J,2)
            ANTNZ     = ANTN(J,3)
            ANTN(J,1) = ANTN(K,1)
            ANTN(J,2) = ANTN(K,2)
            ANTN(J,3) = ANTN(K,3)
            ANTN(K,1) = ANTNX
            ANTN(K,2) = ANTNY
            ANTN(K,3) = ANTNZ
         ENDIF
         IF((ANTN(J,1).EQ.ALN(1)).AND.&
     &      (ANTN(J,2).EQ.ALN(2)).AND.&
     &      (ANTN(J,3).EQ.ALN(3))) THEN
               NDTN = J
               GOTO 2000
         ENDIF
      ENDDO
      ENDDO
 2000 CONTINUE


! -- Y座標: 正
      NA = 1
      XJ0(1,NA)=AUP(1) + ENDR
      YJ0(1,NA)=AUP(2)
      ZJ0(1,NA)=AUP(3)
      write(6,*) XJ0(1,NA),YJ0(1,NA),ZJ0(1,NA)
      DO NJ=1,NDTP
         XJ0(NJ+1,NA)=ANTP(NJ,1)
         YJ0(NJ+1,NA)=ANTP(NJ,2)
         ZJ0(NJ+1,NA)=ANTP(NJ,3)
      write(6,*) XJ0(NJ+1,NA),YJ0(NJ+1,NA),ZJ0(NJ+1,NA)
      ENDDO
      XJ0(NDTP+2,NA)=ALP(1) - ENDR
      YJ0(NDTP+2,NA)=ALP(2)
      ZJ0(NDTP+2,NA)=ALP(3)
      write(6,*) XJ0(NDTP+2,NA),YJ0(NDTP+2,NA),ZJ0(NDTP+2,NA)
      write(6,*) 'NEW-NDTP=',NDTP

      JNUM0(NA)=NDTP+2
! -- Y座標: 負
      NA = 2
      XJ0(1,NA)=AUN(1) + ENDR
      YJ0(1,NA)=AUN(2)
      ZJ0(1,NA)=AUN(3)
      write(6,*) XJ0(1,NA),YJ0(1,NA),ZJ0(1,NA)
      DO NJ=1,NDTN
         XJ0(NJ+1,NA)=ANTN(NJ,1)
         YJ0(NJ+1,NA)=ANTN(NJ,2)
         ZJ0(NJ+1,NA)=ANTN(NJ,3)
      write(6,*) XJ0(NJ+1,NA),YJ0(NJ+1,NA),ZJ0(NJ+1,NA)
      ENDDO
      XJ0(NDTN+2,NA)=ALN(1) - ENDR
      YJ0(NDTN+2,NA)=ALN(2)
      ZJ0(NDTN+2,NA)=ALN(3)
      write(6,*) XJ0(NDTN+2,NA),YJ0(NDTN+2,NA),ZJ0(NDTN+2,NA)
      write(6,*) 'NEW-NDTN=',NDTN
      JNUM0(NA)=NDTN+2
!
      NAMAX=2
!      CALL MODANT(IERR)
!      IF(IERR.NE.0) GOTO 9000
!
 9000 RETURN
      END SUBROUTINE MBDEFA2
