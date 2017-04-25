C
C Routine to read in populations for a particular ionization stage (eg C2) of
C a given species (eg C). The ion population (eg DC2) and oscilator date
C are aslo returned.
C
	SUBROUTINE RD_ION_POP_V3(CIII,DCIII,CIII_PRES,NCIII,
	1             OSCDATE,DESC,FORMAT_DATE,LUIN,ND,
	1             SCRAT,LUSCRAT,SCRATREC)
	IMPLICIT NONE
C
C Created 02-Jun-1996 : Revised fro RD_ION_POP_V2 so that can use dynamic
C                          dimension of arrays in DISPGEN.
C
	INTEGER NCIII
	INTEGER ND
	REAL*8 CIII(NCIII*ND)
	REAL*8 DCIII(ND)
C
	CHARACTER*(*) DESC,FORMAT_DATE
	CHARACTER*(*) OSCDATE
	LOGICAL CIII_PRES
C
	INTEGER LUIN
	INTEGER LUSCRAT,SCRATREC		!Scratch file unit.
	LOGICAL SCRAT
C
C Local variables.
C
	CHARACTER*132 STRING
	INTEGER I,LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
C
	IF(.NOT. CIII_PRES)RETURN
	IF(FORMAT_DATE .NE. '27-JAN-1992' .AND.
	1         FORMAT_DATE .NE. '30-NOV-1991')THEN
	  WRITE(LUER,*)'Error : invalid FORMAT DATE for '//DESC
	  WRITE(LUER,*)'Error occured in RD_ION_POP_V3'
	  STOP
	END IF
C
	READ(LUIN,'(A)')STRING
	IF( INDEX(STRING,'Oscillator') .EQ. 0 )THEN
	   WRITE(LUER,*)'Error in RD_ION_POP - cant get oscillator date'//
	1             ' for ',DESC
	   STOP
	END IF
	READ(STRING,'(T30,A)')OSCDATE
C
	READ(LUIN,*)(CIII(I),I=1,NCIII*ND)
	READ(LUIN,*)(DCIII(I),I=1,ND)
C
	IF(SCRAT)
	1   CALL WRITSCRAT(CIII,DCIII,NCIII,ND,SCRATREC,LUSCRAT,DESC)
	WRITE(LUER,'(5X,A,T32,I4)')'Number of '//DESC//' levels is',NCIII
C
	RETURN
	END
