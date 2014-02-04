	SUBROUTINE WRITEDC_V2(HYD,HYDLTE,NHYD,DHYD,NION,R,T,ED,V,CLUMP_FAC,
	1               LUM,ND,FILENAME,OPTION,FORM)
	IMPLICIT NONE
C
C Altered 04-Oct-2004 : Output space after depth index so easier to search for.
C Altered 24-Feb-2004 : Now depth index to first line of output. Should not effect
C                         any input files.
C Altered 07-Jul-1997 : CLUMP_FAC inserted in call (now _V2), and now output
C                         as last argument. R now written out with a precision
C                         of 7 decimal digits.
C Altered 26-Jun-1996 : CALL GEN_ASCI_OPEN installed.
C Altered 28-May-1996 : Removed for [jdh.disp]SETVEC routine
C                       DOUBLE PRECISION declaration removed.
C                       
C Altered  4-Aug-1988 : Write Departure coefficients out - not b-1.
C
	INTEGER NHYD,NION,ND,FORM,I,J,IOS
	INTEGER, PARAMETER :: IZERO=0
	REAL*8 HYD(NHYD,ND),HYDLTE(NHYD,ND),DHYD(NION,ND)
	REAL*8 R(ND),T(ND),ED(ND),V(ND),CLUMP_FAC(ND)
	REAL*8 LUM,T1,T2,T3
	CHARACTER(LEN=*)  FILENAME,OPTION
	CHARACTER(LEN=40) NEWNAME
	CHARACTER(LEN=90) FMT
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C	NEWNAME=FILENAME//'.'//OPTION
	NEWNAME=FILENAME
C
	  IF(DHYD(1,ND) .NE. 0)THEN
C
C 1 = H, HeII
C 2 = HeI Sing
C 3 = HeI triplets, CIV,NV
C 4 = CIII,NIV
C 5 = HeI (Singlets and Triplets)
C
	    IF(FORM .EQ. 1)FMT='(1X,1P,5E15.5)'
	    IF(FORM .EQ. 2)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:/,1X,3E15.5,
	1                        :/,(1X,4E15.5))'
	    IF(FORM .EQ. 3)FMT='(1X,1P,2E15.5,:/,1X,3E15.5,:/,(1X,4E15.5))'
	    IF(FORM .EQ. 4)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:/,1X,3E15.5,
	1                        :/,(1X,6E15.5))'
	    IF(FORM .EQ. 5)FMT='(1X,1P,1E15.5,:/,1X,2E15.5,:/,1X,2E15.5,:/,1X,
	1            3E15.5,:/,1X,3E15.5,/:,(1X,6E15.5))'
C
	    I=9
	    CALL GEN_ASCI_OPEN(I,NEWNAME,'REPLACE',' ',' ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error opening D.C, file',NEWNAME
	      WRITE(LUER,*)'IOSTAT=',IOS
	      RETURN
	    END IF
C
	    WRITE(9,'(/,1X,A,T40,A)')'24-FEB-2004','!Format date'
	    WRITE(9,2120)R(ND),LUM,NHYD,ND
	    IF(OPTION(1:2) .EQ. 'DC')THEN
	      DO I=1,ND
	        T1=0.0D0
	        T2=0.0D0
	        DO J=1,NHYD
	          T1=T1+HYD(J,I)
	        END DO
	        T3=SUM(HYDLTE(:,I))
	        DO J=1,NION
	          T2=T2+DHYD(J,I)
	        END DO
	        IF(T2 .NE. 0.0D0)THEN
	          T1=T1/T2
	        ELSE
	          T1=0.0D0
	        END IF
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I,' '
	        IF(T3 .EQ. 0.0D0)THEN
	          WRITE(9,FMT)(HYDLTE(J,I),J=1,NHYD)
	        ELSE
	          WRITE(9,FMT)((HYD(J,I)/HYDLTE(J,I)),J=1,NHYD)
	        END IF
	      END DO
	    ELSE
	      DO I=1,ND
	        T1=0.0D0
	        T2=0.0D0
	        DO J=1,NHYD
	          T1=T1+HYD(J,I)
	        END DO
	        DO J=1,NION
	          T2=T2+DHYD(J,I)
	        END DO
	        T1=T1/T2
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I),I,' '
	        WRITE(9,FMT)(HYD(J,I),J=1,NHYD)
	      END DO
	    END IF
	    CLOSE(UNIT=9)
	  END IF
C
2120	  FORMAT(/,1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)
2122	  FORMAT(/,1X,1P,E15.7,6E15.5,2X,I4,A1)
C
	  RETURN
	  END
