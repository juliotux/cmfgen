C
C Dummy error routine to replace NAG routine.
C
	SUBROUTINE F04AYF
	IMPLICIT NONE
C
C Created 28-May-1996
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
	WRITE(LUER,*)'F04AYF is a NAG routine and is currently unavailable'
        WRITE(LUER,*)'Please create a new routine or LINK with NAG'
	WRITE(LUER,*)'Generally this routine is not required'
	WRITE(LUER,*)'Check your options in the VADAT routine'
C
	STOP
	END
