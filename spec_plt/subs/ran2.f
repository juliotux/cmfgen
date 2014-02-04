C
C Subroutine to generate a random number. Routine is from Numerical Recipes,
C and is supposed to be excellent, although a little slower than some other
C routines.
C
	FUNCTION RAN2(IDUM)
	IMPLICIT NONE
C
	INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	REAL*4 RAN2,AM,EPS,RNMX
C
	PARAMETER (IM1=2147483563)
	PARAMETER (IM2=2147483399)
	PARAMETER (AM=1./IM1)
	PARAMETER (IMM1=IM1-1)
	PARAMETER (IA1=40014,IQ1=53688,IR1=12211)
	PARAMETER (IA2=40692,IQ2=52774,IR2=3791)
	PARAMETER (NTAB=32)
	PARAMETER (NDIV=1+IMM1/NTAB)
	PARAMETER (EPS=1.2E-07)
	PARAMETER (RNMX=1.-EPS)
C
	INTEGER IDUM2,J,K,IV(NTAB),IY
	SAVE IV,IY,IDUM2
C
	DATA IDUM2/123456789/
	DATA IV/NTAB*0/
	DATA IY/0/
C
C To initialize, set IDUM < 0 to iniitialize, or reinitialize the
C sequence.
C
	IF(IDUM .LE. 0)THEN
	  IDUM=MAX(-IDUM,1)
	  IDUM2=IDUM
	  DO J=NTAB+8,1,-1
	    K=IDUM/IQ1
	    IDUM=IA1*(IDUM-K*IQ1)-K*IR1
	    IF(IDUM .LT. 0)IDUM=IDUM+IM1
	    if(j .le. NTAB)IV(J)=IDUM
	  END DO
	  IY=IV(1)
	END IF
C
C Program starts here when not initiallizing.
C
	K=IDUM/IQ1
	IDUM=IA1*(IDUM-k*iq1)-K*IR1
	IF(IDUM .LT. 0)IDUM=IDUM+IM1
	K=IDUM2/IQ2
	IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
	IF(IDUM2 .LT. 0)IDUM2=IDUM2+IM2
	J=1+IY/NDIV
	IY=IV(J)-IDUM2
	IV(J)=IDUM
	IF(IY .LT. 1)IY=IY+IMM1
	RAN2=MIN(AM*IY,RNMX)
C
	RETURN
	END