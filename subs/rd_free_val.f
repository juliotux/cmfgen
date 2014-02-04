C
C Reads in a single double precision number from a string. Routine
C can also be used to read in an integer value.
C
	FUNCTION RD_FREE_VAL(STRING,STR_ST,STR_END,NEXT,DESC)
C
C Altered 26-May-1996 - ERROR_LU installed.
C Created 21-Aug-1991.
C
	IMPLICIT NONE
	INTEGER STR_ST,STR_END,NEXT
	CHARACTER*(*) STRING,DESC
	REAL*8 RD_FREE_VAL
C
	CHARACTER*1 TAB,NUL,SPACE,COMMA
	INTEGER I,SS,LEN_STR,IOS
C
	INTEGER ERROR_LU,lUER
	EXTERNAL ERROR_LU
C
	NUL=CHAR(0)
	SPACE=' '
	TAB=CHAR(9)
	COMMA=','
	IOS=0
C
C Strip leading blanks and tabs.
C
	LEN_STR=MIN(LEN(STRING),STR_END)
	SS=STR_ST
	DO WHILE( STRING(SS:SS) .EQ. TAB .OR.
	1     STRING(SS:SS) .EQ. SPACE )
	  SS=SS+1
	  IF(SS .GT. LEN_STR)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in RD_FREE_VAL : No value in string'
	    WRITE(LUER,*)'RD_FREE_VAL called by',DESC
	    STOP
	  END IF
	END DO
C
C Strip all trailing rubbish (i.e. find end of string).
C
	I=SS+1
	DO WHILE( I .LE. LEN_STR .AND. .NOT.
	1 (STRING(I:I) .EQ. NUL .OR.
	1     STRING(I:I) .EQ. TAB .OR.
	1     STRING(I:I) .EQ. SPACE .OR.
	1     STRING(I:I) .EQ. COMMA ) )
	  I=I+1
	END DO
C
	READ(STRING(SS:I-1),'(BN,F20.0)',IOSTAT=IOS)RD_FREE_VAL
	NEXT=I+1
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error on reading vale in RD_FREE_VAL'
	  WRITE(LUER,*)'RD_FREE_VAL called by',DESC
	  STOP
	END IF
	RETURN
	END
