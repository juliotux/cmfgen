C
C Subroutine to compute the sphericity factors Q from the FEAUTRIER
C factors F.
C
	SUBROUTINE QFROMF(F,Q,R,TA,TB,ND)
	IMPLICIT NONE
C
C Altered 26-May-1996 - GENERIC calls now used for DLOG and DEXp
C                       TWO used.
C Created 17-FEB-1986
C
	INTEGER ND,I,IFAIL
	REAL*8 F(ND),Q(ND),R(ND),TA(ND),TB(ND)
C
	REAL*8, PARAMETER :: TWO=2.0D0
C
C TA and TB are work vectors.
C
C Reverse ordering of arrays so the we integrate outwards from the
C surface of the star.
C
	DO I=1,ND
	  TA(ND-I+1)=(3.0D0*F(I)-1.0D0)/R(I)/F(I)
	  TB(ND-I+1)=R(I)
	END DO
C
	CALL INTEGRATE(TB,TA,Q,IFAIL,ND)
C
C Note well - Q(ND) will be undefined after exiting INTEGRATE.
C
	DO I=1,ND-1
	  TB(I)=Q(ND-I)
	END DO
C
	DO I=1,ND-1
	  Q(I)=EXP(TB(I)+TWO*LOG(R(ND)/R(I)))
	END DO
	Q(ND)=1.0D0
C
	RETURN
	END
