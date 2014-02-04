C
	SUBROUTINE NEW_WRITEDC_V3(HYD,HYDLTE,WHYD,
	1               EDGEHYD,GHYD,NHYD,
	1               DHYD,GION,NION,R,T,ED,V,CLUMP_FAC,
	1               LUM,ND,FILENAME,OPTION,FORM)
	IMPLICIT NONE
C
C 05-Mar-1999 - TEXCITE now implicitly dimensioned by NHYD, not 200
C               WHYD (Level dissolution coeffiecients) inclued in call for
C                 forrect evaluation of TX.
C 07-Jul-1997 - CLUMP_FAC installed in call, and called _V2.
C                 R now written with 7 digits of precision.
C 17-Nov-1989 -Based on WRITEDC : TX option installed.
C
	INTEGER NHYD,NION,ND,FORM,I,J,COUNT
	REAL*8 HYD(NHYD,ND),HYDLTE(NHYD,ND),WHYD(NHYD,ND)
	REAL*8 EDGEHYD(NHYD),GHYD(NHYD)
	REAL*8 DHYD(NION,ND),R(ND),T(ND),ED(ND),V(ND),CLUMP_FAC(ND)
	REAL*8 LUM,T1,T2
	REAL*8 DELTA_T,GION
	REAL*8 TEXCITE(NHYD)
	CHARACTER*(*)FILENAME,OPTION
	CHARACTER*30 NEWNAME
	CHARACTER*90 FMT
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
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
	    IF(FORM .EQ. 1)FMT='(1X,1P5E15.5)'
	    IF(FORM .EQ. 2)FMT='(1X,1P1E15.5,:/,1X,2E15.5,:/,1X,3E15.5,
	1                        :/,(1X,4E15.5))'
	    IF(FORM .EQ. 3)FMT='(1X,1P2E15.5,:/,1X,3E15.5,:/,(1X,4E15.5))'
	    IF(FORM .EQ. 4)FMT='(1X,1P1E15.5,:/,1X,2E15.5,:/,1X,3E15.5,
	1                        :/,(1X,6E15.5))'
	    IF(FORM .EQ. 5)FMT='(1X,1P1E15.5,:/,1X,2E15.5,:/,1X,2E15.5,:/,1X,
	1            3E15.5,:/,1X,3E15.5,/:,(1X,6E15.5))'
C
	    OPEN(UNIT=9,STATUS='NEW',FILE=NEWNAME)
C
	    WRITE(9,'(/,1X,A,T40,A)')'07-Jul-1997','!Format date'
	    WRITE(9,2120)R(ND),LUM,NHYD,ND
	    IF(OPTION .EQ. 'DC')THEN
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
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I)
	        WRITE(9,FMT)((HYD(J,I)/HYDLTE(J,I)),J=1,NHYD)
	      END DO
	    ELSE IF(OPTION .EQ. 'TX')THEN
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
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I)
	        DO J=1,NHYD
	          T1=LOG( HYD(J,I)*GION/GHYD(J)/WHYD(J,I)/
	1                             2.07078D-22/ED(I)/DHYD(1,I) )
	          T2=HDKT*EDGEHYD(J)
	          DELTA_T=10.0D0
		  COUNT=0
	          TEXCITE(J)=T(I)
	          DO WHILE( ABS(DELTA_T) .GT. 1.0E-06 .AND. COUNT .LT. 100 )
	            COUNT=COUNT+1
	            DELTA_T=( T1- T2/TEXCITE(J) + 1.5D0*LOG(TEXCITE(J)) )*
	1                     TEXCITE(J)/(T2/TEXCITE(J)+1.5D0)
	            IF(ABS(DELTA_T) .GT. 0.8*TEXCITE(J))DELTA_T=0.5D0*DELTA_T
	            TEXCITE(J)=TEXCITE(J)-DELTA_T
	          END DO
	          IF(COUNT .EQ. 100)THEN
	            WRITE(6,*)'Error - TEXC didnt converge in 100 iterations'
	            WRITE(6,*)'I,J=',I,J
	            RETURN
	          END IF
	        END DO
	        WRITE(9,FMT)(TEXCITE(J),J=1,NHYD)
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
	        WRITE(9,2122)R(I),DHYD(1,I),ED(I),T(I),T1,V(I),CLUMP_FAC(I)
	        WRITE(9,FMT)(HYD(J,I),J=1,NHYD)
	      END DO
	    END IF
	    CLOSE(UNIT=9)
	  END IF
C
2120	  FORMAT(/,1X,F9.4,5X,1PE10.4,5X,0P,I4,5X,I4)
2122	  FORMAT(/,1X,1P,E15.7,6E15.5)
C
	  RETURN
	  END
