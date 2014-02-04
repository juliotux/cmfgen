C
C Subroutine to compute the X vector (see notes) . This routine
C uses either the diffusion approximation or a Schuster condition
C for the lower boundary condition.
C
	SUBROUTINE XVECD(DTAU,SOURCE,X,DIFF,DBC,IC,LS,NC,ND,NI)
	IMPLICIT NONE
C
C Altered 28-May-1996 : IMPLICIT NONE isntalled.
C Created 28-JUL-1982 : Schuster condition inserted
C
	INTEGER LS,NI,ND,NC
	REAL*8 DTAU(ND),SOURCE(ND),X(ND)
	REAL*8 DBC,IC
	LOGICAL DIFF
C
C Local variables.
C
	INTEGER I
C
	X(1)=0.0D0
	DO 10 I=2,NI-1
	  X(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5D0
10	CONTINUE
C
	IF(LS .LE. NC .AND. DIFF)THEN
	  X(NI)=DBC
	ELSE IF(LS .GT. NC)THEN
	  X(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	ELSE
	  X(NI)=IC
	END IF
	
C
	RETURN
	END
