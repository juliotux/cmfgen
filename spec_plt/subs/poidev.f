!
! Returns as a floating point number an INTEGER VALUE that is a random deviate
! drawn from a Posson distribution of mean XM using RAN2_DP as a soure of
! random deviates.
!
! Routine is from Numerical Recipes.
!
      FUNCTION POIDEV(XM,IDUM)
      IMPLICIT NONE
!
! Altered: 10-Apr-2015 : RAN2_DP accessed (there was a precision issue)
!                        SAVE variables fixed.
! 
      REAL*8 XM
      REAL*8 POIDEV
      INTEGER IDUM
!
      REAL*8 T
      REAL*8 Y
      REAL*8 GAMMLN
      REAL*8 EM
      REAL*8 ALXM
!
      REAL*8 RAN2_DP
      REAL*8, SAVE :: G,SQ,AXLM
      REAL*8, SAVE :: OLDM=-1
      REAL*8, PARAMETER :: PI=3.141592654D0
!
      IF (XM .LT. 12.0D0)THEN
        IF (XM.NE.OLDM) THEN
          OLDM=XM
          G=EXP(-XM)
        ENDIF
        EM=-1
        T=1.0D0
2       EM=EM+1.0D0
        T=T*RAN2_DP(IDUM)
        IF (T .GT. G) GO TO 2
      ELSE
        IF (XM.NE.OLDM) THEN
          OLDM=XM
          SQ=SQRT(2*XM)
          ALXM=LOG(XM)
          G=XM*ALXM-GAMMLN(XM+1.0D0)
        ENDIF
1       Y=TAN(PI*RAN2_DP(IDUM))
        EM=SQ*Y+XM
        IF (EM .LT. 0.0D0) GO TO 1
        EM=INT(EM)
        T=0.9D0*(1.0D0+Y**2)*EXP(EM*ALXM-GAMMLN(EM+1.0D0)-G)
        IF (RAN2_DP(IDUM) .GT. T) GO TO 1
      ENDIF
      POIDEV=EM
!
      RETURN
      END
!
      FUNCTION GAMMLN(XX)
      IMPLICIT NONE
      REAL*8 XX
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL*8 GAMMLN
      INTEGER J
!
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
 !
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
