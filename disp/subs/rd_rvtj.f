C
C Set of two routines to READ in desriptor populations.
C The Hydrogen and Helium populations aren now READ in from their
C own file.
C
	SUBROUTINE RD_RVTJ_PARAMS_V2(RMDOT,LUM,ABUNDH,TIME,NAME_CONV,
	1                  ND,NC,NP,FILNAME,LUIN)
	IMPLICIT NONE
!
! Altered 11-Nov-2009 : Instaled format date 10-Nov-2009. Changes will allow
!                           this version to be used with existing routines without
!                           updating to the new version (with TGREY etc).
! Altered 15-Jun-2000 : Naming convention inserted.
!
	INTEGER ND,NC,NP,LUIN
	REAL*8 RMDOT,LUM,ABUNDH
	CHARACTER*(*) TIME,FILNAME,NAME_CONV
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local Variables
C
	CHARACTER*11 LOC_FORMAT_DATE,FORMAT_DATE,PRODATE
	INTEGER NCF
	LOGICAL RD_FIX_T
C
	OPEN(UNIT=LUIN,FILE=FILNAME,STATUS='OLD',ACTION='READ')
C
	LOC_FORMAT_DATE='08-JAN-1996'
	READ(LUIN,'(T30,A11)')FORMAT_DATE
	CALL SET_CASE_UP(FORMAT_DATE,0,0)
	IF(FORMAT_DATE .NE. '15-JUN-2000' .AND. 
	1   FORMAT_DATE .NE. '10-NOV-2009' .AND. 
	1                      FORMAT_DATE .NE. LOC_FORMAT_DATE)THEN
	  WRITE(ERROR_LU(),*)'Wrong format date : RVTJ read failure'
	  WRITE(ERROR_LU(),*)'Subroutine called is RD_ASC_RVTJ_V2'
	  WRITE(ERROR_LU(),*)'Subroutine date 1 is: ',LOC_FORMAT_DATE
	  WRITE(ERROR_LU(),*)'Subroutine date 2 is: ','15-JUN-2000'
	  WRITE(ERROR_LU(),*)'Subroutine date 3 is: ','10-NOV-2009'
	  WRITE(ERROR_LU(),*)'File format date is: ',FORMAT_DATE
	  STOP
	END IF
	READ(LUIN,'(T30,A20)')TIME
	READ(LUIN,'(T30,A11)')PRODATE
	READ(LUIN,'(T30,BN,BZ,I5)')ND
	READ(LUIN,'(T30,BN,BZ,I5)')NC
	READ(LUIN,'(T30,BN,I5)')NP
	READ(LUIN,'(T30,BN,I5)')NCF
C
	READ(LUIN,'(T30,1PE12.5)')RMDOT
	READ(LUIN,'(T30,1PE12.5)')LUM
	READ(LUIN,'(T30,1PE12.5)')ABUNDH
	READ(LUIN,'(T30,L1)')RD_FIX_T
	IF(FORMAT_DATE .EQ. '15-JUN-2000' .OR. FORMAT_DATE .EQ. '10-NOV-2009')THEN
	   READ(LUIN,'(T30,A)')NAME_CONV
	ELSE
	   NAME_CONV='X_FOR_I'
	END IF
C
	RETURN
	END
C
C 
C
C This routine reads in the vecors, and returns the number of HI levels.
C
	SUBROUTINE RD_RVTJ_VEC(R,V,SIGMA,ED,T,
	1       ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POP_ION,
	1       MASS_DENSITY,CLUMP_FAC,ND,LUIN)
	IMPLICIT NONE
C
	INTEGER ND
	REAL*8 R(ND),V(ND),SIGMA(ND),ED(ND)
	REAL*8 T(ND)
	REAL*8 POP_ATOM(ND),POP_ION(ND)
	REAL*8 MASS_DENSITY(ND),CLUMP_FAC(ND)
	REAL*8 ROSS_MEAN(ND),FLUX_MEAN(ND)
C
	INTEGER LUIN
C
C Local variables
C
	CHARACTER*80 STRING
	INTEGER I
C
C At present, all vectors, are assumed to be output in same order.
C Could be changed by using descriptor string if required. We perform
C checks to confirm sequencing.
C
	CALL CHK_STRING(STRING,LUIN,'Radius','RD_RVTJ')
	READ(LUIN,*)(R(I),I=1,ND)
C
	CALL CHK_STRING(STRING,LUIN,'Velocity','RD_RVTJ')
	READ(LUIN,*)(V(I),I=1,ND)
C
	CALL CHK_STRING(STRING,LUIN,'dlnV/dlnr-1','RD_RVTJ')
	READ(LUIN,*)(SIGMA(I),I=1,ND)
C
	CALL CHK_STRING(STRING,LUIN,'Electron','RD_RVTJ')
	READ(LUIN,*)(ED(I),I=1,ND)
C
	CALL CHK_STRING(STRING,LUIN,'Temperature','RD_RVTJ')
	READ(LUIN,*)(T(I),I=1,ND)
C
C In this routine, we skip over the continuum data.
C
	STRING=' '
	DO WHILE( INDEX(STRING,'Rosseland Mean Opacity') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	READ(LUIN,*)(ROSS_MEAN(I),I=1,ND)
C
	CALL CHK_STRING(STRING,LUIN,'Flux Mean Opacity','RD_RVTJ')
	READ(LUIN,*)(FLUX_MEAN(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Atom Density','RD_RVTJ')
	READ(LUIN,*)(POP_ATOM(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Ion Density','RD_RVTJ')
	READ(LUIN,*)(POP_ION(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Mass Density','RD_RVTJ')
	READ(LUIN,*)(MASS_DENSITY(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Clumping Factor','RD_RVTJ')
	READ(LUIN,*)(CLUMP_FAC(I),I=1,ND)
C
	RETURN
	END
