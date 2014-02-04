	SUBROUTINE BETANEW(CHI,CHIL,SIGMA,BETA,R,Z,V,FREQ,ND)
	IMPLICIT NONE
C
C Altered 14-Jun-1996 : Implicit none installed.
C
	INTEGER ND
	REAL*8 CHI(ND),CHIL(ND),BETA(ND),R(ND),Z(ND)
	REAL*8 V(ND),SIGMA(ND)
	REAL*8 FREQ
C
	INTEGER, PARAMETER :: MAX=21
	REAL*8, PARAMETER :: H=0.05
	REAL*8 Y(MAX)
C
	REAL*8 TAU,RMU,P
	REAL*8 X1,X2,Y1,Y2
	REAL*8 T1,T2,T3
	INTEGER I,K,IMU,IE
C
	DO K=1,ND
C
	  TAU=3.0D-10*CHIL(K)*R(K)/V(K)/FREQ
	  DO IMU=1,MAX
	    IF(IMU .EQ. MAX)THEN
	      RMU=0
	      P=0
	    ELSE
	      RMU=H*(IMU-1)
	      P=R(K)*SQRT(1.0D0-RMU*RMU)
	    END IF
	    IF(P .LE. R(ND))THEN
	      T2=50
	      DO I=1,K
	        Z(I)=SQRT(R(I)*R(I)-P*P)
	      END DO
	      T1=0.66666666*CHI(1)*R(1)
	      DO I=1,K-1
	        T1=T1+(CHI(I+1)+CHI(I))*(Z(I)-Z(I+1))
	      END DO
	      T1=T1*0.5
	    ELSE
	      DO I=K+1,ND,1
	        IF(P .GT. R(I))THEN
	          IE=I-1
	          GOTO 500
	        END IF
	      END DO
C
500	      CONTINUE
	      DO I=1,IE
	        Z(I)=SQRT(R(I)*R(I)-P*P)
	      END DO
	      T1=CHI(1)*R(1)/3.0D0
	      IF(K .NE. 1)THEN
	        DO I=1,K-1
	          T1=T1+(CHI(I+1)+CHI(I))*(Z(I)-Z(I+1))*0.5
	        END DO
	      END IF
	      T1=T1
	      T2=0.0
	      IF(IE .GT. K)THEN
	        DO I=K,IE-1
	          T2=T2+(CHI(I+1)+CHI(I))*(Z(I)-Z(I+1))*0.5
	        END DO
	      END IF
C
	      Y1=LOG(CHI(IE+1))
	      Y2=LOG(CHI(IE))
	      X1=LOG(R(IE+1))
	      X2=LOG(R(IE))
	      T3=EXP((Y2-Y1)*(LOG(P)-X1)/(X2-X1)+Y1)
	      T2=T2+(CHI(IE)+T3)*0.5*Z(IE)
	      T2=T2*2.0+T1
	    END IF
C
	    T3=TAU/(1.0D0+SIGMA(K)*RMU*RMU)
	    IF(ABS(T3) .LT. 1.0D-04)THEN
	      Y(IMU)=1.0D0-T3*(0.5-T3*(1.0-T3/4.0D0)/6.0D0)
	    ELSE
	      Y(IMU)=(1.0D0-EXP(-T3))/T3
	    END IF
	    Y(IMU)=Y(IMU)*0.5*(EXP(-T1)+EXP(-T2))
	  END DO
C
	  BETA(K)=0.0
	  DO I=1,MAX-2,2
	   BETA(K)=BETA(K)+Y(I)+Y(I+1)*4.0D0+Y(I+2)
	  END DO
	  BETA(K)=BETA(K)*H/3.0D0
C
	END DO
C
C
	RETURN
	END
