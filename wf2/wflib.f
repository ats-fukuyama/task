C     $Id$
C
C     ***** REGULA FALSI METHOD *****
C
      SUBROUTINE FRGFLS(XS,XE,DX,XR,FUNC,EPS,ILL)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C
      EXTERNAL FUNC
      ILMAX=-ILL
      IF(ILMAX.EQ.0) ILMAX=-30
      ILL=0
      IF(ABS(DX).LE.1.E-15) THEN
C        WRITE(6,600)
         ILL=900
         RETURN
      ENDIF
      IF((XE-XS)*DX.LE.0.) THEN
C        WRITE(6,601)
         ILL=901
         RETURN
      ENDIF
      D=DX
C
      X=XS
      Y=FUNC(X)
      IF(ABS(Y).LE.EPS) THEN
         XR=X
         RETURN
      ENDIF
C
   10 IF((XE-X)*D.LE.0.) THEN
C        WRITE(6,602)
         ILL=902
         RETURN
      ENDIF
      XX=X
      YY=Y
      X=XX+D
      Y=FUNC(X)
      IF(ABS(Y).LE.EPS) THEN
         XR=X
         RETURN
      ENDIF
      IF(Y*YY.GT.0.) GOTO 10
C
   20 ILL=ILL-1
      IF(ILL.LE.ILMAX) THEN
         XR=X
C        WRITE(6,603)
         ILL=903
         RETURN
      ENDIF
      IF(Y*YY.LE.0.) THEN
         YYY=Y
         DD=-D
      ELSE
         XX=X
         YY=Y
         DD=D+DD
      ENDIF
      D=DD*YY/(YYY-YY)
      X=XX+D
      Y=FUNC(X)
      IF(ABS(Y).GT.EPS.AND.ABS(D).GT.EPS) GOTO 20
C
      XR=X
      RETURN
C
C 600 FORMAT(1H ,'## FRGFLS : ABS(DX) .LE. 1.E-14 !!')
C 601 FORMAT(1H ,'## FRGFLS : (XE-XS)*DX .LE. 0. !!')
C 602 FORMAT(1H ,'## FRGFLS : NO ROOT BETWEEN XS AND XE !!')
C 603 FORMAT(1H ,'## FRGFLS : NO CONVERGENSE !!')
      END
C
C     ******************* BESSEL2 SERIES2 *******************
C
      FUNCTION DKBES(N,X)
      REAL*8  X,DKBES
      REAL*8  XABS,RECX,HALFX,W,Y,SUM,EXPX,D75,D65,
     & AI0,AI1,AI2,AK0,AK1,AK2,NUMER,DENOM,Z
      REAL*8 XS,YS,FK,C69,C173
      EQUIVALENCE (XABS,XS),(Y,YS)
      DATA D65,D75,C69,C173 / 1.0D65,1.0D75,69.0,173.0 /
C
      XABS=X
      NABS=IABS(N)
      IF(XABS) 1,1,3
    1 DKBES=D75
      WRITE(6,2) N,X
    2 FORMAT(1H ,5X,'DKBES IS INVALID. N=',I10,'  ,  X=',D23.16)
      RETURN
    3 IF(NABS.GE.30000) GOTO 33
      IF(XS.LT.C173) GOTO 4
   33 DKBES=0.0D0
      WRITE(6,34) N,X
   34 FORMAT(1H ,5X,'DKBES IS NOT ACCURATE. N=',I10,'  ,  X=',D23.16)
      RETURN
    4 RECX=1.0D0/XABS
      Z=RECX+RECX
C
      IF(XS.GT.2.0) GOTO 17
      HALFX=XABS*0.5D0
      W=DLOG(HALFX)
      IF(XS.GE.0.0001) GOTO 12
C
      Y=HALFX*HALFX
      AI0=1.0D0+Y
      AK0=-W*AI0-5.7721566490153286 D-1+4.22784335 D-1*Y
      IF(NABS) 6,5,6
    5 DKBES=AK0
      RETURN
    6 CONTINUE
      AK1=(RECX-HALFX*AK0)/AI0
      IF(NABS-1) 8,7,8
    7 AK2=AK1
      GOTO 26
    8 L=-INT(C69/W)
      IF(NABS.LE.L) GOTO 22
    9 DKBES=D75
      WRITE(6,10) N,X
   10 FORMAT(1H ,5X,'DKBES IS OVERFLOW. N=',I10,'  ,  X=',D23.16)
      RETURN
C
   12 CONTINUE
      AI2=0.D0
      AI1=1.0D-75
      SUM=AI1
      K=18
      IF(XS.GE.1.0) GOTO 13
      K=14
      IF(XS.GE.0.1) GOTO 13
      K=9
   13 FK=K
C
      AI0=(AI1*Z)*FK+AI2
      K=K-1
      IF(K) 14,15,14
   14 SUM=SUM+AI0
      AI2=AI1
      AI1=AI0
      GOTO 13
   15 SUM=SUM+SUM+AI0
      EXPX=DEXP(XABS)
      AI0=(AI0/SUM)*EXPX
C
      Y=HALFX*HALFX
C
      AK0=(((1.1D-17*YS+1.533D-15)*YS+1.78593D-13)*YS+1.709994D-11)*YS
      AK0=((AK0+1.31674867D-9)*Y+7.9350965213D-8)*Y
      AK0=-W*AI0+((((((AK0      +3.61262410320     D-6)*Y
     & +1.18480393641097  D-4)*Y+2.61478761880521  D-3)*Y
     & +3.489215745643890 D-2)*Y+2.3069608377461679D-1)*Y
     & +4.2278433509846714D-1)*Y-5.7721566490153286D-1
      IF(NABS)16,5,16
   16 CONTINUE
C
      AI1=(AI1/SUM)*EXPX
      AK1=(RECX-AK0*AI1)/AI0
      IF(NABS.EQ.1) GOTO 7
      GOTO 22
C
   17 Y=RECX*0.5D0
      W=DEXP(-XABS)*DSQRT(3.1415926535897932D0*Y)
      IF(NABS.EQ.1) GOTO 29
      IF(XS.GE.6.0) GOTO 18
      NUMER = ( ( ( ( ( ( (
     &  9.607359468936920D-1 *Y+2.529252282967791D+2)*Y
     & +7.970033517428499D+3)*Y+7.569208233441645D+4)*Y
     & +3.057001687861112D+5)*Y+6.316930369273204D+5)*Y
     & +7.490317987301367D+5)*Y+5.503149540522960D+5)*Y
      NUMER = ( (   NUMER
     & +2.642871289210102D+5)*Y+8.617173613052524D+4)*Y
C
      NUMER = ( ( ( ( ( ( (     NUMER
     & +1.959077941045564D+4)*Y+3.161402444515176D+3)*Y
     & +3.658576787591567D+2)*Y+3.045630347257958D+1)*Y
     & +1.815056853732798D 0)*Y+7.630840968696385D-2)*Y
     & +2.198080790671052D-3)*Y
      NUMER = ( ( (     NUMER
     & +4.108400060979278D-5)*Y+4.468603071316900D-7)*Y
     & +2.138087593931531D-9)*Y
C
      DENOM = ( ( ( ( ( ( (    Y
     & +3.555555555555556D+2)*Y+1.513244444444444D+4)*Y
     & +1.956717714285714D+5)*Y+1.079473194053918D+6)*Y
     & +3.045125484052925D+6)*Y+4.914134724130594D+6)*Y
     & +4.892294125356680D+6)*Y
      DENOM = ( ( (     DENOM
     & +3.168987751788133D+6)*Y+1.388013537038700D+6)*Y
     & +4.227486032369927D+5)*Y
C
      DENOM = ( ( ( ( ( ( (     DENOM
     & +9.133105119891823D+4)*Y+1.418093261050334D+4)*Y
     & +1.593555554804435D+3)*Y+1.296908735314483D+2)*Y
     & +7.594657346992128D 0)*Y+3.149536504573292D-1)*Y
     & +8.975302543644857D-3)*Y
      DENOM = ( (    DENOM
     & +1.663376533185147D-4)*Y+1.797062622699450D-6)*Y
     & +8.552350375726116D-9
C
      GOTO 19
C
   18 CONTINUE
      NUMER = ( ( ( ( ( ( (
     &  9.344272845910382D-1 *Y+7.602616403443910D+1)*Y
     & +7.395743295243645D+2)*Y+2.143299351353469D+3)*Y
     & +2.591022778768068D+3)*Y+1.552060383786633D+3)*Y
     & +5.070112211483554D+2)*Y+9.496661840062622D+1)*Y
      NUMER = ( ( ( (  NUMER
     & +1.035232548375690D+1)*Y+6.417643499456601D-1)*Y
     & +2.076251506438793D-2)*Y+2.696430527842588D-4)*Y
C
      DENOM = ( ( ( ( ( ( (    Y+128.0D0)*Y
     & +1.952426666666667D+3)*Y+8.925379047619048D+3)*Y
     & +1.700072199546485D+4)*Y+1.598598652282462D+4)*Y
     & +8.186476153700395D+3)*Y+2.418159110016117D+3)*Y
      DENOM = ( ( ( (    DENOM
     & +4.239448497375430D+2)*Y+4.421129278670809D+1)*Y
     & +2.659325881907253D 0)*Y+8.426345399508084D-2)*Y
     & +1.078572211137035D-3
C
   19 AK0=W*(1.0D0-NUMER/DENOM)
      IF(NABS.EQ.0) GOTO 5
   29 CONTINUE
      IF(XS.GE.6.0) GOTO 20
C
      NUMER = ( ( ( ( ( ( (
     &  4.612351477146149D+1 *Y+5.662820251461119D+3)*Y
     & +1.267729206953107D+5)*Y+9.840279611657081D+5)*Y
     & +3.490851284631075D+6)*Y+6.612450857584894D+6)*Y
     & +7.385469550146852D+6)*Y+5.204090498624053D+6)*Y
      NUMER = ( (   NUMER
     & +2.426670289964806D+6)*Y+7.748520914563417D+5)*Y
C
      NUMER = ( ( ( ( ( ( (   NUMER
     & +1.735592435557683D+5)*Y+2.771333212784498D+4)*Y
     & +3.183321561628970D+3)*Y+2.636216725447433D+2)*Y
     & +1.565459511840883D+1)*Y+6.565891381439845D-1)*Y
     & +1.888507334035547D-2)*Y
      NUMER = ( ( (   NUMER
     & +3.526831364930226D-4)*Y+3.834684961199849D-6)*Y
     & +1.834777493397057D-8)*Y
C
      DENOM = ( ( ( ( ( ( (    Y
     & +640.0D0             )*Y+3.242666666666667D+4)*Y
     & +4.565674666666667D+5)*Y+2.649616021768707D+6)*Y
     & +7.729933921057426D+6)*Y+1.277675028273955D+7)*Y
     & +1.295019033182651D+7)*Y
      DENOM = ( ( (      DENOM
     & +8.506230281115516D+6)*Y+3.767465314819330D+6)*Y
     & +1.157963565388285D+6)*Y
      DENOM = ( ( ( ( ( ( (    DENOM
     & +2.520737013090144D+5)*Y+3.939147947362039D+4)*Y
     & +4.450965515143424D+3)*Y+3.639711612011614D+2)*Y
     & +2.140312525061418D+1)*Y+8.908688970078743D-1)*Y
     & +2.547045316439757D-2)*Y
      DENOM = ( (    DENOM
     & +4.734225517526958D-4)*Y+5.128203094044773D-6)*Y
     & +2.446369991196075D-8
C
      GOTO 21
C
   20 CONTINUE
      NUMER = ( ( ( ( ( ( (
     &  2.727430674283096D+1 *Y+1.138548998351567D+3)*Y
     & +8.476942664996772D+3)*Y+2.136946909366115D+4)*Y
     & +2.388756033308294D+4)*Y+1.367517138363589D+4)*Y
     & +4.350909161812486D+3)*Y+8.026823492415450D+2)*Y
      NUMER = ( ( ( (    NUMER
     & +8.677093809042819D+1)*Y+5.356622184821665D 0)*Y
     & +1.730209588698993D-1)*Y+2.247025439868823D-3)*Y
      DENOM = ( ( ( ( ( ( (  Y+ 2.304            D+2)*Y
     & +4.183771428571428D+3)*Y+2.082588444444444D+4)*Y
     & +4.172904489795918D+4)*Y+4.057981194255480D+4)*Y
     & +2.128483799962103D+4)*Y+6.401009408866192D+3)*Y
      DENOM = ( ( ( (  DENOM
     & +1.137957228242879D+3)*Y+1.200020804210648D+2)*Y
     & +7.284240459137258D 0)*Y+2.325671330264231D-1)*Y
     & +2.996033919825096D-3
C
   21 CONTINUE
      AK1=W*(1.0D0+NUMER/DENOM)
      IF(NABS.EQ.1) GOTO 7
C
   22 CONTINUE
      K=1
   23 FK=K
      IF(DABS(AK1).GE.D65) GOTO 25
      AK2=AK1*Z*FK+AK0
   24 K=K+1
      IF(K.GE.NABS) GOTO 26
      AK0=AK1
      AK1=AK2
      GOTO 23
   25 W=AK1*1.0D-10
      Y=AK0*1.0D-10
      AK2=W*Z*FK+Y
      IF(AK2.GE.D65) GOTO 9
      AK2=AK2/1.0D-10
      GOTO 24
C
   26 CONTINUE
      IF(X) 27,28,28
   27 K=NABS/2
      IF(K+K.NE.NABS) AK2=-AK2
   28 DKBES=AK2
      RETURN
      END
C
C     ****** PLASMA DISPERSION FUNCTION ******
C
      SUBROUTINE ZETA(X,Z,DZ,DDZ)
C
C                  ORIGINAL PROGRAM BY T. WATANABE
C                          MODIFIED BY A. FUKUYAMA (1984.1.24)
C
      IMPLICIT REAL*8       (A-H,O-Z)
      COMPLEX*16            X,Z,DZ,DDZ
      REAL*8                EN(25),FN(25),QN(25),HN(25)
      SAVE                  I
      DATA                  I /0/
C
      IF(I.EQ.0) THEN
         I   =1
         M   =12
         M2  =M+M+1
         H   =0.5D0
         H4  =H*0.25D0
         HS  =H*H
         HS4 =HS*0.25D0
         PAI =3.141592653589793D0
         PAI2=PAI+PAI
         PI  =PAI/H
         PH  =PAI2/H
         SPAI=3.544907701811032D0
         RPH2=0.1128379167095513D1*H
         RPH4=RPH2+RPH2
         RPH5=HS*HS*RPH2*6.25D-2
         DO 10 L=1,M2
            LL   =L*L
            HN(L)=LL*HS4
            QN(L)=LL*LL
            EN(L)=EXP((1-L)*HS)
            FN(L)=EN(L)
   10    CONTINUE
         EN(1)=EXP(-HS4)*RPH5
         FN(1)=EXP(-HS4)*RPH2
         EN(2)=EN(2)*RPH5
      ENDIF
C
      XRE=DBLE(X)
      AXR=ABS(XRE)
      XIM=DIMAG(X)
      AXI=ABS(XIM)
      XXR=(XRE-XIM)*(XRE+XIM)
      XXI=XRE*XIM*2.D0
      IF(AXR.LE.1.D-36) THEN
         SXR2=0.D0
      ELSE
         SXR2=XRE*XRE*2.D0
      ENDIF
      IF(AXI.LE.1.D-36) THEN
         SXI2=0.D0
      ELSE
         SXI2=XIM*XIM*2.D0
      ENDIF
      IF(ABS(XXI).LE.1.D-36) THEN
         SXI=0.D0
      ELSE
         SXI=XXI*XXI
      ENDIF
      RR   =SXR2+SXI2
      XIMPH=XIM*PH
      IF(XIM.GT.PI.OR.XXR.GT.169.D0-XIMPH) THEN
         IF(AXI.LE.H4.AND.ABS(AXR-INT((AXR+H4)/H)*H).LE.H4) THEN
            J=0
         ELSE
            J=1
         ENDIF
         AR=HN(M2-J)-XXR
         SI=EN(M2-J)*QN(M2-J)/(AR*AR+SXI)
         SR=SI*AR
         DO 100 L=J+2,M2-1,2
            AR=HN(M2-L)-XXR
            D =QN(M2-L)/(AR*AR+SXI)
            SR=(D*AR+SR)*EN(M2-L)
            SI=(D+SI)*EN(M2-L)
  100    CONTINUE
         AR =XXR-SXI2
         D  =XXR+SXR2
         DZR=RR*RR
         FI =8.D0/(DZR*RR)
         FR =FI*AR
         FI =FI*D
         AR =0.5D0-SR
         DZI=(SI*XXR+AR)*8.D0/DZR
         DZR=(XXR*AR-SI*SXI)*8.D0/DZR
         AR =-AR-XXR
         D  =SI-1.D0
         Z  =DCMPLX((FI*SXI2*D+FR*AR)*XRE,(FR*SXR2*D-FI*AR)*XIM)
         AR =SR*4.D0+DZR
         D  =SI*4.D0-DZI
         DDZ=DCMPLX((D*SXI2+AR)*XRE,(D*SXR2-AR)*XIM)/RR*2.D0
         DZ =DCMPLX(DZR,-DZI*XXI)
      ELSEIF(XXR.LT.-169.D0) THEN
         Z  =DCMPLX(1.D22,1.D22)
         DZ =DCMPLX(1.D22,1.D22)
         DDZ=DCMPLX(1.D22,1.D22)
      ELSE
         IF(AXR.GE.AXI) THEN
            XRI=(XIM-INT(XXI/PAI2)*PAI/XRE)*XRE*2.D0
         ELSEIF(AXI.LE.1.D-36) THEN
            XRI=0.D0
         ELSE
            XRI=(XRE-INT(XXI/PAI2)*PAI/XIM)*XIM*2.D0
         ENDIF
         XRH=XRE*PH
         IF(XIM.GT.-PI) THEN
            EPZIH=EXP(XIMPH)
         ELSE
            EPZIH=0.D0
         ENDIF
         IF(XXR.GT.170.D0) XXR=170.D0
         WA=EXP(-XXR)*SPAI
         IF(AXI.LE.H4.AND.ABS(AXR-INT((AXR+H4)/H)*H).LE.H4) THEN
            AR=(EPZIH+COS(XRH)*2.D0)*EPZIH+1.D0
            WR=(SIN(XRI)+SIN(XRI-XRH)*EPZIH)/AR
            WI=(COS(XRI)+COS(XRI-XRH)*EPZIH)/AR
            J =0
         ELSE
            AR=(EPZIH-COS(XRH)*2.D0)*EPZIH+1.D0
            WR=(SIN(XRI)-SIN(XRI-XRH)*EPZIH)/AR
            WI=(COS(XRI)-COS(XRI-XRH)*EPZIH)/AR
            J =1
         ENDIF
         IF(AXR.GT.2.D0.OR.AXI.GT.2.D0) THEN
            AR=HN(M2-J)-XXR
            SI=EN(M2-J)*QN(M2-J)/(AR*AR+SXI)
            SR=SI*AR
            DO 110 L=J+2,M2-1,2
               AR=HN(M2-L)-XXR
               D =QN(M2-L)/(AR*AR+SXI)
               SR=(D*AR+SR)*EN(M2-L)
               SI=(D+SI)*EN(M2-L)
  110       CONTINUE
            AR =XXR-SXI2
            D  =XXR+SXR2
            DZR=RR*RR
            FI =8.D0/(DZR*RR)
            FR =FI*AR
            FI =FI*D
            AR =0.5D0-SR
            DZI=(SI*XXR+AR)*8.D0/DZR
            DZR=(XXR*AR-SI*SXI)*8.D0/DZR
            DZ =DCMPLX(DZR,-DZI*XXI)
            AR =-AR-XXR
            D  =SI-1.D0
            Z  =DCMPLX((FI*SXI2*D+FR*AR)*XRE,(FR*SXR2*D-FI*AR)*XIM)
            AR =SR*4.D0+DZR
            D  =SI*4.D0-DZI
            DDZ=DCMPLX((D*SXI2+AR)*XRE,(D*SXR2-AR)*XIM)/RR*2.D0
            DZ =DCMPLX(DZR,-DZI*XXI)
         ELSE
            IF(J.EQ.0) THEN
               AR=HN(M2)-XXR
               SI=FN(M2)/(AR*AR+SXI)
               SR=SI*AR
               DO 120 L=2,M2-1,2
                  AR=HN(M2-L)-XXR
                  D =1.D0/(AR*AR+SXI)
                  SR=(D*AR+SR)*FN(M2-L)
                  SI=(D+SI)*FN(M2-L)
  120          CONTINUE
               Z  =DCMPLX((SR-SXI2*SI)*XRE,(SXR2*SI+SR)*XIM)
               DZ =DCMPLX(SI*SXI-SR*XXR-1.D0,-(SI*XXR+SR)*XXI)*2.D0
               AR =(XXR*4.D0-2.D0)*SR-SI*SXI*4.D0+4.D0
               D  =(XXR*4.D0-2.D0)*SI+SR*4.D0
               DDZ=DCMPLX((AR-SXI2*D)*XRE,(SXR2*D+AR)*XIM)
            ELSE
               AR=HN(M2-1)-XXR
               SI=EN(M2-1)/(AR*AR+SXI)
               SR=SI*SR
               DO 130 L=3,M2-2,2
                  AR=HN(M2-L)-XXR
                  D =1.D0/(AR*AR+SXI)
                  SR=(D*AR+SR)*FN(M2-L)
                  SI=(D+SI)*FN(M2-L)
  130          CONTINUE
               Z  =DCMPLX((SR-SXI2*SI-1.D0/RR)*XRE,
     &                    (SR+SXR2*SI+1.D0/RR)*XIM)*RPH2
               D  =(SI*SXI-XXR*SR+0.5D0)*RPH4-2.D0
               AR =(XXR*SI+SR)*RPH4
               DZ =DCMPLX(D,-AR*XXI)
               DDZ=-(DCMPLX((D+AR*SXI2)*XRE,(D-AR*SXR2)*XIM)+Z)*2.D0
            ENDIF
         ENDIF
         Z  =DCMPLX(WR,WI)*WA+Z
         DZ =DCMPLX(WI*XIM-WR*XRE,-WR*XIM-WI*XRE)*WA*2.D0+DZ
         DDZ=DCMPLX((XXR-0.5D0)*WR-XXI*WI,(XXR-0.5D0)*WI+XXI*WR)
     &             *WA*4.D0+DDZ
      ENDIF
      RETURN
      END
C
C     ****** COMPLEX FFT   USING LIST VECTOR ******
C
C                                         CODING BY K.HAMAMATSU
C
      SUBROUTINE FFT2L ( A , B , T , LIST , N2 , IND , KEY , LP )
C
      IMPLICIT COMPLEX * 16 ( A-C , E-H , O-Z )
      REAL * 8     DIV
      DIMENSION    A ( N2*2 ) , B( N2*2 ) , T( N2*LP , 2 ) ,
     &             LIST( N2*LP )
      SAVE DIV
C
      IF(LP.EQ.0) THEN
         B(1) = A(1)
         RETURN
      ENDIF
C     IF ( IND .EQ. 1 ) THEN SET TALBES AND LIST
      IF ( IND .EQ. 1 ) THEN
           CALL SETTBL ( N2 , T , B , LP )
           CALL SETLST ( LIST , N2 , LP )
           DIV = 1.D0 / (N2*2)
           IND = 0
      ENDIF
      K =  1
      L = N2
      DO 10 I = 1 , LP
           M = (I-1)*N2 + 1
           IF ( K .EQ. 1 ) THEN
                CALL FFTSUB( A , B , T(M,KEY) , LIST(M) , L , N2 )
           ELSE
                CALL FFTSUB( B , A , T(M,KEY) , LIST(M) , L , N2 )
           ENDIF
           K = K * (-1)
           L = L / 2
 10   CONTINUE
C     IF ( KEY .EQ. 2 ) THEN INVERSE TRANSFORMATION
C     ELSE    NORMAL TRANSFORMATION
      IF ( KEY .EQ. 1 ) THEN
           IF ( K .EQ. 1 ) THEN
                DO 20 I = 1 , N2*2
 20                 B( I ) = A( I ) * DIV
            ELSE
                DO 30 I = 1 , N2*2
 30                 B( I ) = B( I ) * DIV
            ENDIF
      ELSE IF ( K .EQ. 1 ) THEN
            DO 40 I = 1 , N2*2
 40             B( I ) = A ( I )
      ENDIF
      RETURN
      END
C
C     ====== TRANSFORMATION FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE FFTSUB ( A , B , T , LIST , L , N2 )
C
      IMPLICIT COMPLEX * 16 ( A-H , O-Z )
      DIMENSION A( N2*2 ) , B( N2 , 2 ) , T( N2 ) , LIST( N2 )
      DO 10 I = 1 , N2
          B( I , 1 ) = A( LIST( I ) ) + A( LIST( I ) + L ) * T( I )
 10   CONTINUE
      DO 20 I = 1 , N2
          B( I , 2 ) = A( LIST( I ) ) - A( LIST( I ) + L ) * T( I )
 20   CONTINUE
      RETURN
      END
C     ====== TABLE SETTING FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE SETTBL ( N2 , T , B , LP )
C
      COMPLEX * 16     T , B
      REAL    *  8     TR , TI , PAI
      DIMENSION       T( N2 , LP , 2 ) , B( N2 , 2 )
      DATA PAI / 3.1415926535898D0 /
C
      DO 10 I = 1 , N2
          TR = DCOS( 2.D0 * PAI * (I-1) / (N2*2) )
          TI = DSIN( 2.D0 * PAI * (I-1) / (N2*2) )
          B( I , 1 ) = DCMPLX( TR , -TI )
          B( I , 2 ) = DCMPLX( TR ,  TI )
 10   CONTINUE
      K  =  1
      NN = N2
      DO 40 L = 1 , LP
          DO 30 J = 0 , K-1
              DO 20 I = 1 , NN
                  T( I+J*NN , L , 1 ) = B( 1+NN*J , 1 )
                  T( I+J*NN , L , 2 ) = B( 1+NN*J , 2 )
 20           CONTINUE
 30       CONTINUE
          K  = K * 2
          NN = NN / 2
 40   CONTINUE
      RETURN
      END
C
C     ====== LIST SETTING FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE SETLST ( LIST , N2 , LP )
C
      DIMENSION  LIST ( N2 , LP )
      N1 = N2
      NN =  1
      DO 30 K = 1 , LP
          M = 0
          DO 20 J = 1 , NN
              DO 10 I = 1 , N1
                  M = M + 1
                  LIST ( M , K ) = I + (J-1) * 2 * N1
 10           CONTINUE
 20       CONTINUE
          N1 = N1 / 2
          NN = NN * 2
 30   CONTINUE
      RETURN
      END
C
      FUNCTION BESSI0(X)
C
      REAL*8 BESSI0,X
      REAL*8 AX
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y
      SAVE   P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/
     &     1.0D0,3.5156229D0,3.0899424D0,
     &     1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/
     &     0.39894228D0, 0.1328592D-1, 0.225319D-2,
     &    -0.157565D-2,  0.916281D-2, -0.2057706D-1,
     &     0.2635537D-1,-0.1647633D-1, 0.392377D-2/
C
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI0=(EXP(AX)/SQRT(AX))
     &        *(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6
     &                              +Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
      FUNCTION BESSI1(X)
C
      REAL*8 BESSI1,X
      REAL*8 AX
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Y
      SAVE   P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/
     &     0.5D0,0.87890594D0,0.51498869D0,
     &     0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/
     &     0.39894228D0,-0.3988024D-1,-0.362018D-2,
     &     0.163801D-2, -0.1031555D-1, 0.2282967D-1,
     &    -0.2895312D-1, 0.1787654D-1,-0.420059D-2/
C
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI1=(EXP(AX)/SQRT(AX))
     &        *(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6
     &                              +Y*(Q7+Y*(Q8+Y*Q9))))))))
        IF(X.LT.0.D0) BESSI1=-BESSI1
      ENDIF
      RETURN
      END
C
      FUNCTION  CELI1( KSQ )
C***********************************************************************
C*                                                                     *
C*  PURPOSE                                                            *
C*      COMPUTES THE COMPLETE ELLIPITC INTEGRAL OF THE FIRST KIND      *
C*      IN DOUBLE PRECISION.                                           *
C*                                                                     *
C*  FUNCTION SUBPROGRAM                                                *
C*      CELI1( KSQ )                                                   *
C*                                                                     *
C*  DESCRIPTION OF PARAMETER                                           *
C*      KSQ  - SQUARE OF THE MODULUS K.                                *
C*                                                                     *
C*  COPYRIGHT                                                          *
C*      Y. ONODERA,    AUGUST, 1989.                                   *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8 (A-Z)
      DATA  ONE / 1.0D0 /
C
      IF ( KSQ .LT. 0 ) GOTO 20
      X = ONE - KSQ
C
      IF ( KSQ .LT. 0.9D0 ) THEN
C                             LANDEN TRANSFORMATION.
        K1 = SQRT( X )
        A = KSQ / ( K1 + ONE )**2
        Z = A * A
        F = A + ONE
        IF ( KSQ .GE. 0.27D0 ) THEN
C                             ONCE MORE LANDEN TRANSFORMATION.
          K1 = SQRT( ONE - Z )
          B = Z / ( K1 + ONE )**2
          Z = B * B
          F = F * ( B + ONE )
        ENDIF
C                             NOW  Z < 0.0062
        CELI1 = F * ( ( ( ( ( (        0.0799 D0               *Z
     &  + 0.0951 308D0           )*Z + 0.1174 4540 4D0        )*Z
     &  + 0.1533 9807 879 D0     )*Z + 0.2208 9323 3455 5D0   )*Z
     &  + 0.3926 9908 1698 724D0 )*Z + 1.5707 9632 6794 89662D0 )
        RETURN
C
      ELSE
        IF ( X .EQ. 0 ) GOTO 10
        IF ( X .LT. 0 ) GOTO 20
        F = ((((( ((((( (((
     &      3.0D-2               * X  +  3.24D-2             ) * X
     &   +  3.495D-2           ) * X  +  3.7958D-2           ) * X
     &   +  4.1524 6D-2        ) * X  +  4.5829 54D-2        ) * X
     &   +  5.1127 765D-2      ) * X  +  5.7806 3778D-2      ) * X
     &   +  6.6482 4962 6D-2   ) * X  +  7.8202 0568 85D-2   ) * X
     &   +  9.4884 2366 536D-2 ) * X  +  1.2044 2708 3333D-1 ) * X
     &   +  0.1640 625D0       ) * X  +  0.25D0              ) * X
        G = ((((( ((((( (((
     &      1.1D-2               * X  +  1.20D-2             ) * X
     &   +  1.300D-2           ) * X  +  1.4144D-2           ) * X
     &   +  1.5522 7D-2        ) * X  +  1.7199 67D-2        ) * X
     &   +  1.9282 673D-2      ) * X  +  2.1939 3969D-2      ) * X
     &   +  2.5444 5076 0D-2   ) * X  +  3.0281 0668 95D-2   ) * X
     &   +  3.7384 0332 031D-2 ) * X  +  0.0488 2812 5D0     ) * X
     &   +  0.0703 125D0       ) * X  +  0.125D0             ) * X
     &   +  0.5D0
        CELI1 = G * LOG( 16.0D0 / X ) - F
        RETURN
      ENDIF
C
   10 CELI1 = 100.0D0
      RETURN
   20 WRITE(6,30) KSQ
      CELI1 = 0
      RETURN
   30 FORMAT(' ARGUMENT OF CELI1 IS INVALID.  X=',1PD23.15)
C *****************************************************
C              NUMERICAL TABLE OF CELI1
C       ARGUMENT                     CELI1
C -----------------------------------------------------
C   256/1024 ( 0.25       )   1.6857 5035 4812 5960
C   276/1024 ( 0.269531...)   1.6964 8703 4351 1121
C   277/1024 ( 0.270508...)   1.6970 3221 5452 5081
C   512/1024 ( 0.5        )   1.8540 7467 7301 3719
C   747/1024 ( 0.729492...)   2.1212 9467 8256 3108
C   748/1024 ( 0.730469...)   2.1229 0614 6685 8874
C   768/1024 ( 0.75       )   2.1565 1564 7499 6432
C   921/1024 ( 0.899414...)   2.5753 4344 0273 8417
C   922/1024 ( 0.900391...)   2.5799 3389 3410 7172
C *****************************************************
      END
C
      FUNCTION  CELI2( KSQ )
C***********************************************************************
C*                                                                     *
C*  PURPOSE                                                            *
C*      COMPUTES THE COMPLETE ELLIPITC INTEGRAL OF THE SECOND KIND     *
C*      IN DOUBLE PRECISION.                                           *
C*                                                                     *
C*  FUNCTION SUBPROGRAM                                                *
C*      CELI2( KSQ )                                                   *
C*                                                                     *
C*  DESCRIPTION OF PARAMETER                                           *
C*      KSQ  - SQUARE OF THE MODULUS K.                                *
C*                                                                     *
C*  COPYRIGHT                                                          *
C*      Y. ONODERA,    AUGUST, 1989.                                   *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8 (A-Z)
      DATA    ONE / 1.0D0 /
C
      IF ( KSQ .LT. 0 ) GOTO 20
      X = ONE - KSQ
C
      IF ( KSQ .LT. 0.9D0 ) THEN
C                             LANDEN TRANSFORMATION.
        K1 = SQRT( X )
        A = KSQ / ( K1 + ONE )**2
        Z = A * A
        IF ( Z .LT. 0.1D0 ) THEN
          F = 0
        ELSE
          F = ((((( ((
     &       1.90D-5       * Z  +  2.24D-5   ) * Z  +  2.661D-5   ) * Z
     &    +  3.201D-5    ) * Z  +  3.8980D-5 ) * Z  +  4.81239D-5 ) * Z
     &    +  6.036662D-5 ) * Z  +  7.714147D-5 ) * Z
        ENDIF
C
        CELI2 = ((((( ((((( ((     F   + 1.0075 620 D-4           ) * Z
     &  + 1.3508 5878 D-4        ) * Z + 1.8697 0073 8D-4         ) * Z
     &  + 2.6923 6906 2D-4       ) * Z + 4.0783 8153 87D-4        ) * Z
     &  + 6.6063 0397 91D-4      ) * Z + 1.1744 5404 0725D-3      ) * Z
     &  + 2.3968 4498 1071 D-3   ) * Z + 6.1359 2315 1542 6D-3    ) * Z
     &  + 2.4543 6926 0617 03D-2 ) * Z + 3.9269 9081 6987 2415D-1 ) * Z
     &  + 1.5707 9632 6794 8966 2D0 ) / ( A + ONE )
        RETURN
C
      ELSE
        IF ( X .EQ. 0 ) GOTO 10
        IF ( X .LT. 0 ) GOTO 20
        F = ((((( ((((( (((
     &      3.1D-2               * X  +  3.36D-2            ) * X
     &   +  3.643D-2           ) * X  +  3.9701D-2          ) * X
     &   +  4.3624 1D-2        ) * X  +  4.8406 36D-2       ) * X
     &   +  5.4364 881D-2      ) * X  +  6.1993 3844D-2     ) * X
     &   +  7.2105 7891 8D-2   ) * X  +  8.6143 4936 52D-2  ) * X
     &   + 10.6913 2486 979D-2 ) * X  +  0.1406 25D0        ) * X
     &   +  0.2031 25D0        ) * X  +  0.25D0             ) * X
        G = ((((( ((((( (((
     &      1.2D-2               * X  +  1.25D-2            ) * X
     &   +  1.355D-2           ) * X  +  1.4817D-2          ) * X
     &   +  1.6339 7D-2        ) * X  +  1.8211 41D-2       ) * X
     &   +  2.0568 184D-2      ) * X  +  2.3627 0428D-2     ) * X
     &   +  2.7757 6446 5D-2   ) * X  +  3.3645 6298 83D-2  ) * X
     &   +  4.2724 6093 75D-2  ) * X  +  0.0585 9375D0      ) * X
     &   +  0.0937 5D0         ) * X  +  0.25D0             ) * X
        CELI2 = G * LOG( 16.0D0 / X ) - F + ONE
        RETURN
      ENDIF
C
   10 CELI2 = ONE
      RETURN
   20 WRITE(6,30) KSQ
      CELI2 = 0
      RETURN
   30 FORMAT(' ARGUMENT OF CELI2 IS INVALID.  X=',1PD23.15)
C *****************************************************
C              NUMERICAL TABLE OF CELI2
C        ARGUMENT                  CELI2
C -----------------------------------------------------
C   256/1024 ( 0.25       )    1.4674 6220 9339 4272
C   276/1024 ( 0.269531...)    1.4588 9478 1418 3167
C   277/1024 ( 0.270508...)    1.4584 6425 7030 7126
C   512/1024 ( 0.5        )    1.3506 4388 1047 6755
C   747/1024 ( 0.729492...)    1.2238 2430 2500 7830
C   748/1024 ( 0.730469...)    1.2232 2324 8354 8634
C   768/1024 ( 0.75       )    1.2110 5602 7568 4595
C   921/1024 ( 0.899414...)    1.1052 5395 8315 6887
C   922/1024 ( 0.900391...)    1.1044 5483 7388 2954
C *****************************************************
      END
