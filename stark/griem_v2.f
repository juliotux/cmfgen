!
      SUBROUTINE GRIEM_V2(PR,DWS,NWS,ED_IN,TEMP_IN,VDOP,
     *                               IL,IU,ZZ,AMASS,RET_LOG)
      USE GRIEM_STARK_MOD
      IMPLICIT NONE
! 
! Create a profile table using the griem theory as used by A and M
! in their AP J S24 paper.
! 
! NWS is the number of wavelengths profile is to be computed at. 
!
      INTEGER NWS                              
      INTEGER IL,IU
      REAL*8 DWS(NWS),PR(NWS)
      REAL*8 ED_IN
      REAL*8 TEMP_IN
      REAL*8 VDOP
      REAL*8 ZZ
      REAL*8 AMASS
      LOGICAL RET_LOG
!
      REAL*8, PARAMETER :: CLIGHT=2.997925D18
      REAL*8, PARAMETER :: VLIGHT=2.997925D10
      REAL*8, PARAMETER :: CBOLTZ=1.3805D-16
      REAL*8, PARAMETER :: AMU=1.67333D-24
      REAL*8, PARAMETER :: THIRD=0.333333333D0
      REAL*8, PARAMETER :: GAMIN=0.01D0
      REAL*8, PARAMETER :: GAMAX=25.0D0
      REAL*8, PARAMETER :: PI=3.1415926536E0
!
      REAL*8 AA,ASQ,BSQ,FRQ,WAVE
      REAL*8 CON,OBA,WSQ,CCKAB,BETP1,BETW1
      REAL*8 RYD
      REAL*8 GAM01,GAM02
      REAL*8 EE,SRE,CUE
      REAL*8 F0,TT,SRT
      REAL*8 DOP
!
      REAL*8 BETA,BETAP,BETAW
      REAL*8 GAMA,GAM0,GAML,GAMW
      REAL*8 FPW,SSS,REN,FAC
      REAL*8 CKABD,OCKAB
!
      INTEGER I,JP,JM
!
! Functions called
!
      REAL*8 TBG_STRK
      REAL*8 TH_STRK
      REAL*8 DCONV_STRK
      EXTERNAL TBG_STRK,TH_STRK,DCONV_STRK
!
! RYD need not be this accurate.
!
      RYD=1.09737312D+05*VLIGHT*ZZ*ZZ*AMASS/(AMASS+5.48593E-04)
!
! Set up the constant terms for evaluation of S(BET,GAM) and DONV
!
      AA=IL
      ASQ=IL*IL
      BSQ=IU*IU
      FRQ=RYD/ASQ-RYD/BSQ               !Line frequency in Hz.
      WAVE=CLIGHT/FRQ			!Wavelength in Angstroms.
      CON=WAVE/FRQ
      OBA=1.0E0/(BSQ-ASQ)
      WSQ=WAVE*WAVE
      CCKAB=5.5E-5*(ASQ*BSQ)**2*OBA/ZZ**5
      BETP1=2.995E-15*WSQ
      BETW1=6.951E-9*WSQ*ZZ/BSQ
      GAM01=3.78E-5*(BSQ-3.0E0*AA)*BSQ*OBA/ZZ
      GAM02=4.0E6*ZZ/BSQ
C
C The loops over electron density, and temperature have been removed.
C
C Set Electron density parameters. The electron density is no longer
C logarithmic.
C
      EE=ED_IN
      SRE=SQRT(EE)
      CUE=EE**THIRD
C
C Normal field strength
C
      F0=1.25E-9*CUE*CUE
C
C Set temperature parameters.
C
      TT=1.0D+04*TEMP_IN  		!As TEMP_IN in units of 10^4 K.
      SRT=SQRT(TT)
C
C Doppler width (A) (old expression was TCON*SRT). 10^10 arrises as we need
C to convert VDOP from km/s to cm/s.
C      TCON=WAVE*SQRT(2.0E0*CBOLTZ/(AMASS*AMU))/VLIGHT
C
      DOP=WAVE*SQRT( 2.0E0*CBOLTZ*TT/(AMASS*AMU) + 
	1               1.0D+10*VDOP*VDOP )/VLIGHT
      ODOP=1.0E0/DOP
C
C Ratio STARK/DOPPLER width.
C
      CKABD=CCKAB*F0*ODOP
C
C 1/STARK WIDTH     (A)**-U
C
      OCKAB=ODOP/CKABD
      BETAP=BETP1*SRE*OCKAB
      BETAW=BETW1*TT*OCKAB
      BETAP=MIN(0.99D0*BETAW,BETAP)
      GAM0=GAM01*CUE*LOG10(GAM02*TT/SRE)/SRT
      GAML=LOG(BETAW/BETAP)
      GAMW=4.712389E0/SQRT(BETAW)
      FPW=GAM0/GAML
C
C Initialize GAMA for ZERO,BETAP
C
      GAMA=GAMW+GAM0
C
C Generate S(BETA) on a fixed basis.
C
      JM=NBET
      JP=JM
      DO 6 I=1,NBET
        BETA=BET(I)
        IF(BETA.LT.BETAP) GOTO 11
        IF(BETA.GT.BETAW) GOTO 12
C
C ASSIGN GAMA FOR BETAP,BETAW
C
        GAMA=GAMW+FPW*LOG(BETAW/BETA)
C
C Look for regions of the GAMA,BETA plane where fast methods are ok.
C This is not just to save time as the T(BETA,GAMA) analysis
C Is numerically unstable in these same regions.
C
   11   IF(GAMA.GT.GAMAX) GOTO 13
        IF(GAMA.LT.GAMIN*BETA) GOTO 14
        IF(GAMA.LT.GAMIN) GOTO 14
C FULL CASE
        SSS=TBG_STRK(BETA,GAMA)
        GOTO 7
C   BETA GT 100*GAMA FOR BETA LT 1
   14   SSS=TH_STRK(BETA)
        IF(BETA.LT.1.0E0) GOTO 7
C GAMA LT GAMIN OR BETA GT 100*GAMA FOR BETA GT 1
        SSS=SSS+GAMA/(PI*BETA*BETA)
        GOTO 7
C   BETA GT BETAW
   12   SSS=2.0E0*TH_STRK(BETA)
        GOTO 7
C   GAMA GT GAMAX
   13   SSS=GAMA/(PI*(GAMA*GAMA+BETA*BETA))
C FILL THE SYMETRIC SS,SX VECTORS
    7   SS(JM)=SSS
        SX(JM)=-BETA*CKABD
        SS(JP)=SSS
        SX(JP)= BETA*CKABD
        JM=JM-1
        JP=JP+1
    6   CONTINUE
!
! Normalize to unit wavelength integral.
      REN=0.0E0
      DO 8 I=2,NS
        REN=REN+(SS(I)+SS(I-1))*(SX(I)-SX(I-1))
    8 CONTINUE
      FAC=CON*ODOP*2.0E0/REN
!
! Make the asymtotic power law constants.
!
      PS=LOG(SS(NS)/SS(NS-1)) / LOG(SX(NS)/SX(NS-1))
      PS=MIN(PS,-2.0D0)
      AS=SS(NS)/SX(NS)**PS
!
! Map the profile onto the DLAM set by convolution with the doppler
! profile.
!
      DO I=1,NWS
        PR(I)=FAC*DCONV_STRK(DWS(I))
      END DO
      IF(RET_LOG)PR(1:NWS)=LOG( PR(1:NWS) )
!
      RETURN
      END
!
! 
!
      REAL*8 FUNCTION DCONV_STRK(DLAM)
      USE GRIEM_STARK_MOD
      IMPLICIT NONE
!
      REAL*8 DLAM
! 
! Convolution of gaussian profile with S function  6 AUG 77
! 
      REAL*8, PARAMETER :: HALF=0.5E0
      REAL*8, PARAMETER :: ZERO=0.0E0
      REAL*8, PARAMETER :: SRTPI=5.6418958E-1
      REAL*8, PARAMETER :: RANGE=6.0E0
!
      REAL*8 X(NS),ERX(NS),EX(NS)
      REAL*8 DB,X1,CON
      REAL*8 XN,DX,GR,CR,TEX,TERX
      INTEGER I,J,JM,N0
!
! Functions called
!
      REAL*8 DOPLER_STRK
      REAL*8 ASINT_STRK
      REAL*8 ERR_STRK
      EXTERNAL DOPLER_STRK,ASINT_STRK,ERR_STRK
!
      DB=ABS(DLAM)*ODOP
      X1=SX(1)+DB
      IF(X1.LT.RANGE) GOTO 10
C ASINT_STRK -RANGE,RANGE
      DCONV_STRK=ASINT_STRK(-RANGE,RANGE,DB)
      RETURN
   10 DCONV_STRK=ZERO
C ASINT_STRK -RANGE,X1
      IF(X1.GT.-RANGE) DCONV_STRK=DCONV_STRK+ASINT_STRK(-RANGE,X1,DB)
C ASINT_STRK XN,RANGE
      XN=SX(NS)+DB
      IF(XN.LT.RANGE) DCONV_STRK=DCONV_STRK+ASINT_STRK(XN,RANGE,DB)
C SET UP X IN THE GAUSSIAN FRAME
C STORE ALL EXPONENTIAL AND ERROR FUNCTIONS ONCE
      N0=2
      CON=HALF*SRTPI
      DO 20 I=1,NS
         X(I)=SX(I)+DB
         IF(X(I).LT.-RANGE) N0=I+1
         EX(I)=DOPLER_STRK(X(I))
         ERX(I)=ERR_STRK(X(I))
   20 CONTINUE
C BY S SEGMENT
      TERX=ZERO
      TEX=ZERO
      DO 50 J=N0,NS
         JM=J-1
         DX=X(J)-X(JM)
         IF(DX.LT.0.1) GOTO 30
         GR=(SS(J)-SS(JM))/DX
         CR=SS(J)-GR*X(J)
         TERX=TERX+CR*(ERX(J)-ERX(JM))
         TEX=TEX+GR*(EX(J)-EX(JM))
         GOTO 40
   30    DCONV_STRK=DCONV_STRK +  CON*(SS(J)*EX(J)+SS(JM)*EX(JM))*DX
   40    CONTINUE
         IF(X(J).GT.RANGE) GOTO 60
   50 CONTINUE
   60 DCONV_STRK=DCONV_STRK+HALF*(TERX-SRTPI*TEX)
      RETURN
      END
!
      REAL*8 FUNCTION TBG_STRK(BET,GAM)
      IMPLICIT NONE
!
      REAL*8 BET,GAM
!
      REAL*8, PARAMETER :: PI=3.14159265358979
!
! Functions called.
!
      REAL*8 TF_STRK
      REAL*8 TG_STRK
      REAL*8 TH_STRK
      EXTERNAL TF_STRK,TG_STRK,TH_STRK
!
      IF(GAM.GT.1.E-2) GOTO 10
!     SMALL GAMMA
      TBG_STRK=TH_STRK(BET)
      IF(BET.GT.1.)TBG_STRK=TBG_STRK+GAM/PI/BET**2
      RETURN      
!     NORMAL CASE
   10 IF(BET/GAM.GT.1.E2) GOTO 20
      TBG_STRK=GAM*(TF_STRK(BET,GAM)+TF_STRK(-BET,GAM)+
     .            TG_STRK(BET,GAM)+TG_STRK(-BET,GAM))/PI
      RETURN
   20 TBG_STRK=TH_STRK(BET)+GAM/PI/BET**2
!
      RETURN
      END
C
C 
C
      REAL*8 FUNCTION TG_STRK(BET,GAM)
      IMPLICIT NONE
!
      REAL*8 BET,GAM
!
      REAL*8, PARAMETER :: C4=1.5D0
      REAL*8, PARAMETER :: C5=3.4636008D1
      REAL*8, PARAMETER :: C6=-1.3253986E2
      REAL*8, PARAMETER :: PI2=1.57079632679489
!
      REAL*8 G2,B,C,BET2,CC,R
      REAL*8 P,P2,Q,Q2
      REAL*8 SINA,COSA2,SINA2
      REAL*8 Y1,Y2
      REAL*8 X1,X2,X3,X4,X5
      REAL*8 D1,D2,D3,D4,D5,D6,D7,D8,D9
      REAL*8 SUM1,SUM2,SUM3,SUM4
!
      G2=GAM*GAM
      B=2.0  *BET
      C=BET*BET+G2
      BET2=BET*BET
      CC=C*C
      R= SQRT(C)
      Q= SQRT(R)
      Q2=Q*Q
      P=1.  /Q
      P2=P*P
      SINA=GAM/R
      COSA2= SQRT(0.5  *(1.  -BET/R))
      SINA2= SQRT(0.5  *(1.  +BET/R))
      X5=BET+4.
      Y2=LOG((X5*X5+G2)/16.  )
      IF(X5/GAM.GT.1.0  ) GOTO 10
      Y1=PI2- ATAN(X5/GAM)
      GOTO 20
   10 Y1= ATAN(GAM/X5)
   20 D2=C/192.
      D3=-B/32.
      D4=0.25  *(3.  *BET2-G2)/C
      D5=B*(G2-BET2)/CC
      D6=(BET2*(BET2-6.0  *G2)+G2*G2)/CC
      SUM1=((D6*Y1)/GAM+(D4+D3))+D2
      SUM1=(SUM1+D5*Y2)*C5/CC
      D1=C/1024.
      D2=-B/192.
      D3=(3.0  *BET2-G2)/(32.  *C)
      D4=BET*(G2-BET2)/CC
      D5=0.5  *(G2*G2+5.  *BET2*(BET2-2.0  *G2))/(C*CC)
      D6=-BET*(BET2*BET2+5.0  *G2*(G2-2.0  *BET2))/(CC*C)
      SUM2=(D5*Y2+D3)+D1
      SUM2=(SUM2+(D2+D4+(D6*Y1)/GAM))*C6/CC
      D7=C4/C
      D8=D7*(B*B/C-1.  )/(2.  *Q*Q2*SINA)
      D9=D7*(B/(C*C))/(2.  *P*P2*SINA)
      X1=(4.0  -Q2)/(4.0  *Q*SINA2)
      X2=(Q*(Q+4.  *COSA2)+4.  )/(Q*(Q-4.  *COSA2)+4.  )
      X3=(0.25  -P2)/(P*SINA2)
      X4=(P*(P+COSA2)+0.25  )/(P*(P-COSA2)+0.25  )
      IF(X1.GT.1.  ) GOTO 30
      Y1=PI2- ATAN(X1)
      GOTO 40
   30 Y1= ATAN(1.  /X1)
   40 IF(X3.GT.-1.  ) GOTO 50
      Y2=- ATAN(1.  /X3)
      GOTO 60
   50 Y2=PI2+ ATAN(X3)
   60 SUM3=D8*(2.  *COSA2*Y1-SINA2*LOG(X2))
      SUM4=D9*(2.  *COSA2*Y2+SINA2*LOG(X4))
      TG_STRK=(SUM4+D7*(1.  /12.  -B/C))+SUM3
      TG_STRK=TG_STRK+SUM1+SUM2
      RETURN
      END
!
! 
!
      REAL*8 FUNCTION TF_STRK(B,G)
      IMPLICIT NONE
!
      REAL*8 B
      REAL*8 G
!
      REAL*8, PARAMETER :: C0=1.0007744D-01
      REAL*8, PARAMETER :: C1=4.93208719D-03
      REAL*8, PARAMETER :: C2=-7.09873526D-03
      REAL*8, PARAMETER :: C3=7.11559325D-04
!
      REAL*8 X1,G2
      REAL*8 D,D1,D2,D3
!
      D3=C3
      D2=C2-3.0  *B*C3
      D1=(3.0  *C3*B-2.0  *C2)*B+C1
      D =((-C3*B+C2)*B-C1)*B+C0
      G2=G*G
      X1=B+4.
!
      TF_STRK=4.0  *D2+4.0  *D3*(B+2.0  )+0.5*(D1-G2*D3)*LOG((X1*X1+G2)/
     . (B*B+G2)) +(D -G2*D2)*( ATAN(X1/G)- ATAN(B/G))/G
!
      RETURN
      END
!
! 
!
      REAL*8 FUNCTION TH_STRK(X)
      IMPLICIT NONE
!
      REAL*8 X
!
! GRIEM Microfield fudge - normalized to unit B integral
!
      REAL*8, PARAMETER :: C0=1.0007744D-01
      REAL*8, PARAMETER :: C1=4.93208719D-03
      REAL*8, PARAMETER :: C2=-7.09873526D-03
      REAL*8, PARAMETER :: C3=7.11559325D-04
      REAL*8, PARAMETER :: C4=1.5D0
      REAL*8, PARAMETER :: C5=3.4636008D1
      REAL*8, PARAMETER :: C6=-1.3253986D2
!
      REAL*8 B,B2,SRT
!
      B=ABS(X)
      IF(B.GT.4. ) GOTO 10
         TH_STRK=((C3*B+C2)*B+C1)*B+C0
      RETURN
   10    SRT= SQRT(B)
         B2=B*B
         TH_STRK=((C6/B+C5)/B2+C4/SRT)/B2
      RETURN
      END
!
! 
!
      REAL*8 FUNCTION ASINT_STRK(X0,X1,DB)
      USE GRIEM_STARK_MOD
      IMPLICIT NONE
!
      REAL*8 X0,X1,DB
!
! Integral over asymtotic region by SIMPSONS rule.
!
      REAL*8, PARAMETER :: SRTPI=5.6418958D-1
      REAL*8, PARAMETER :: C23=6.666666667D-1
      REAL*8, PARAMETER :: HALF=0.5D0
      REAL*8, PARAMETER :: ONE=1.0D0
      REAL*8, PARAMETER :: TWO=2.0D0
      REAL*8, PARAMETER :: STEP=0.2D0
!
      INTEGER I,N
      REAL*8 H,X,DX
!
! Functions
!
      REAL*8 F
      REAL*8 DOPLER_STRK
!
      F(X)=DOPLER_STRK(X)*ABS(X-DB)**PS
!
! Chose H to be sept
!
      N=(X1-X0)/STEP+ONE
      DX=(X1-X0)/FLOAT(N)
      H=HALF*DX
      X=X0-H
      ASINT_STRK=F(X0)
      DO 10 I=1,N
        X=X+DX
        ASINT_STRK=ASINT_STRK+  TWO*F(X)+F(X+H)
   10 CONTINUE
      ASINT_STRK=ASINT_STRK*AS*SRTPI*H*C23
!
      RETURN
      END
!
! 
!
      REAL*8 FUNCTION DOPLER_STRK(X)
      IMPLICIT NONE
!
      REAL*8 X,SQ
      SQ=X*X
      IF(SQ.GT.64.) THEN
        DOPLER_STRK=0.
      ELSE
        DOPLER_STRK=EXP(-SQ)
      END IF
      RETURN
      END
!
! 
!
! Error function.
!
      REAL*8 FUNCTION ERR_STRK(X)
      IMPLICIT NONE
!
      REAL*8 X
      REAL*8 T,EX,ONE
!
! Functions called
!
      REAL*8 DOPLER_STRK
!
      ONE=1.0D0
      T=ABS(X)
      IF(T.GT.5.)THEN
        ERR_STRK=SIGN(ONE,X)
      ELSE
        T=1./(1.+0.3275911*T)
        ERR_STRK=((((1.061405429  *T-1.453152027  )*T +1.421413741  )*T-
     . 0.284496736  )*T +0.254829592  )*T
        EX=DOPLER_STRK(X)
        ERR_STRK=SIGN(ONE-EX*ERR_STRK,X)
      END IF
!
      RETURN
      END
