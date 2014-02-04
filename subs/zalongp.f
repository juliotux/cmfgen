C
C Subroutine to compute the Z values at NI depth points for a given
C imact parameter P .
C
	SUBROUTINE ZALONGP(R,X,P,NI)
	IMPLICIT NONE
C
C Altered 09-Oct-2004 : Replace R*R-P*P by (R-P)*(R+P). Intel compiler was causing
C                         opitimization problems. 
C Altered 09-Dec-2001 : PP removed and repplaced directly in SQRT by P*P
C Altered 28-May-1996 : IMPLCIT NONE installed.
C
	INTEGER NI
	REAL*8 X(NI),R(NI),P
C
	INTEGER I
C
	DO 10 I=1,NI
	  X(I)=SQRT( (R(I)-P)*(R(I)+P) )
10	CONTINUE
C
	RETURN
	END
