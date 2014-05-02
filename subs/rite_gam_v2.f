C
C Auixlary routine to evauate, and write out, the mean ionic charge
C
	SUBROUTINE RITE_GAM_HEAD(R,ED,T,ND,LU,FILNAME)
C
C Altered 31-Dec-2013 : ND now output as I4 rather than I3
C Created 11-Jun-1996 : based on RITE_GAM. This routine rites out
C header info and population independent vectors (i.e. R,T and ED).
C
C Output file is opened, and left in that state for RITE_GAM_V2.
C
	IMPLICIT NONE
	INTEGER ND,LU
C
	REAL*8 R(ND),ED(ND),T(ND)
	CHARACTER*(*) FILNAME
C
C Local variables.
C
	INTEGER I,IOS,LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	I=0
	CALL GEN_ASCI_OPEN(LU,FILNAME,'UNKNOWN',' ',' ',I,IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to open ',FILNAME,
	1        ' in RITE_GAM_HEAD'
	     WRITE(LUER,*)'IOS=',IOS
	     STOP
	  END IF

	  WRITE(LU,'(A)')
	  WRITE(LU,'(1X,I4,20X,''!Number of depth points'')')ND
C
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' !Electron density'
	  WRITE(LU,500)(ED(I),I=1,ND)
C
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' !Radius (10^10cm)'
	  WRITE(LU,500)(R(I),I=1,ND)
C
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' !Temperature (10^4 K)'
	  WRITE(LU,500)(T(I),I=1,ND)
C
500	  FORMAT(1P,1X,10E12.3)
C
	RETURN
	END
C
C 
C
	SUBROUTINE RITE_GAM_V2(POP,GAM,AT_NO,DESC,ND,LU)
	IMPLICIT NONE
C
	INTEGER LU,ND
	REAL*8 POP(ND)
	REAL*8 GAM(ND)
	REAL*8 AT_NO
	CHARACTER*(*) DESC
C
	INTEGER I,LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	IF(GAM(ND) .EQ. 0.0D0 .AND. GAM(1) .EQ. 0.0D0)RETURN
	WRITE(LU,'(A)')' '
	WRITE(LU,'(1X,F4.0,5X,A,20X,A)')AT_NO,DESC(1:LEN(DESC)),
	1         '!Atomic Number and descriptor.'
	DO I=1,ND
	  IF(POP(I) .NE. 0)THEN
	     GAM(I)=GAM(I)/POP(I)
	  ELSE IF(GAM(I) .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Bad gamma values in RITE_GAM_V2'
	    WRITE(LUER,*)'DESC=',DESC
	    STOP
	  END IF
	END DO
C	    
	WRITE(LU,500)(GAM(I),I=1,ND)
500	FORMAT(1P,1X,10E12.3)
C
	RETURN
	END
