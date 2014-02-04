!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created:
!
!     GENCOOL_SUM:   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT: Only the top rates are printed: Sorted using depths 1, 11, 21 etc. 
!
	PROGRAM MOD_COOL
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 17-Nov-2009: Now read in charge exchange cooling.
!                         Slight format change.
! Altered: 29-Jan-2009: ND is now read in from MODEL (if it exists).
! Altered: 08-Feb-2008: Extra terms (such as V term) sheck and output.
!
	INTEGER, PARAMETER :: MAX_RECS=1000
!
	CHARACTER*132 TMP_STR
	CHARACTER*132 STRING
	CHARACTER*132 STR_VEC(MAX_RECS)
	REAL*8 VALS(MAX_RECS,10)
	REAL*8 TA(MAX_RECS)
	INTEGER INDX(MAX_RECS)
!
	REAL*8, ALLOCATABLE :: BOUND(:)
	REAL*8, ALLOCATABLE :: SUM(:)
!
	INTEGER ND
	INTEGER NV
	INTEGER I,J,K,ID
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	REAL*8 T1
	LOGICAL FILE_OPEN
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	ALLOCATE (SUM(ND))
	ALLOCATE (BOUND(ND))
!
	OPEN(UNIT=20,FILE='GENCOOL',STATUS='OLD',ACTION='READ')
	OPEN(UNIT=21,FILE='GENCOOL_SUM',STATUS='UNKNOWN',ACTION='WRITE')
!
	DO I=1,1+(ND-1)/10
	   ID=1+(I-1)*10
!
	   READ(20,'(A)')STRING
	   WRITE(21,'(A)')TRIM(STRING)
	   DO J=1,3
	     READ(20,'(A)')TMP_STR
	     READ(20,'(A)')STRING
	     IF(J .EQ. 1)WRITE(21,'(A,T12,10I12)')'Depth',(K,K=ID,MIN(ID+9,ND))
	     WRITE(21,'(A)')TMP_STR(4:9)//'     '//TRIM(STRING)
	     READ(20,'(A)')STRING
	   END DO
!
	   T1=0.0D0
	   DO WHILE(1 .EQ. 1)
	      READ(20,'(A)')STRING
	      IF( INDEX(STRING,'Coll') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'COL ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Free-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'FF ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'K-shell') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'XKS ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'V term') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A)')' '
	        WRITE(21,'(A,T12,A)')'AC.R(V).',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'dTdR term') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'AC.R(dT).',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'decay') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'|R. decay|',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Artificial') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'|Art. HT|',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Rate') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'Net C.R.',TRIM(STRING)
	     ELSE IF( INDEX(STRING,'Net') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'% C.R.',TRIM(STRING)
	        IF(I .NE. 1+(ND-1)/10)THEN
	          READ(20,'(A)')STRING
	          WRITE(21,'(A)')STRING
	        END IF
	        EXIT
	      ELSE IF( INDEX(STRING,'Bound-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        SUM=0.0D0
	        DO WHILE(1 .EQ. 1)
	          READ(20,'(A)')STRING
	          IF(STRING .EQ. ' ')THEN
	            TMP_STR=ADJUSTL(TMP_STR)
	            K=INDEX(TMP_STR,' ')
	            WRITE(21,'(A)')STRING
                    WRITE(21,'(A,T13,10ES12.4)')TMP_STR(1:K)//'BF ',(SUM(K),K=ID,MIN(ID+9,ND))
	            EXIT
	          ELSE
	            READ(STRING,*)(BOUND(K),K=ID,MIN(ID+9,ND))
	            SUM=SUM+BOUND
	          END IF
	        END DO
	      ELSE IF( INDEX(STRING,'Charge exchange cooling rate') .NE. 0)THEN
	        TMP_STR=STRING
	        SUM=0.0D0
	        DO WHILE(1 .EQ. 1)
	          READ(20,'(A)')STRING
	          IF(STRING .EQ. ' ')THEN
	            TMP_STR=ADJUSTL(TMP_STR)
	            K=INDEX(TMP_STR,' ')
                    WRITE(21,'(A,T13,10ES12.4)')'Charge',(SUM(K),K=ID,MIN(ID+9,ND))
	            EXIT
	          ELSE
	            READ(STRING,*)(BOUND(K),K=ID,MIN(ID+9,ND))
	            SUM=SUM+BOUND
	          END IF
	        END DO
	      END IF
	   END DO
	END DO
	CLOSE(UNIT=20)
	CLOSE(UNIT=21)
!
	WRITE(6,'(A)')' Cooling data has been written to GENCOOL_SUM'
	WRITE(6,'(A)')' Will now sort data to display most important terms'
	WRITE(6,'(A)')' 12 records is a reasonable number to output '
!
	OPEN(UNIT=20,FILE='GENCOOL_SUM',STATUS='UNKNOWN',ACTION='READ')
	OPEN(UNIT=21,FILE='GENCOOL_SORT',STATUS='UNKNOWN',ACTION='WRITE')
	NRECS=12
	CALL GEN_IN(NRECS,'Maximum number of records per depth to be output to sorted file')
!
	N_INIT_RECS=5
	VALS=0.0D0
	DO I=1,1+(ND-1)/10
	  ID=1+(I-1)*10
	  DO K=1,5
	    READ(20,'(A)')STR_VEC(K)
	  END DO
!
	  K=N_INIT_RECS
	  DO WHILE(1 .EQ. 1)
	    K=K+1
100	    READ(20,'(A)')STR_VEC(K)
	    IF(STR_VEC(K) .EQ. ' ')GOTO 100
	    READ(STR_VEC(K)(12:),*)(VALS(K,J),J=1,MIN(10,ND-(I-1)*10))
	    IF( INDEX(STR_VEC(K),'%') .NE. 0)THEN
	      IF(I .NE. 1+(ND-1)/10 )READ(20,'(A)')STRING			!Getting record with ^L
	      EXIT
	    END IF
	  END DO
!
	  NV=K-N_INIT_RECS-2
	  TA(:)=ABS(VALS(:,1))
	  CALL INDEXX(NV,TA(6),INDX,L_FALSE)
	  NV=NV+N_INIT_RECS+2
!
	  DO K=1,N_INIT_RECS
	    WRITE(21,'(A)')TRIM(STR_VEC(K))
	  END DO
	  WRITE(21,'(A)')' '
	  WRITE(21,'(A)')TRIM(STR_VEC(NV))
	  WRITE(21,'(A)')TRIM(STR_VEC(NV-1))
	  WRITE(21,'(A)')' '
	  DO K=1,MIN(NRECS,NV-N_INIT_RECS-2)
	    WRITE(21,'(A)')TRIM(STR_VEC(INDX(K)+N_INIT_RECS))
	  END DO
!
	END DO
	CLOSE(UNIT=20)
	CLOSE(UNIT=21)
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Sorted cooling data has been written to GENCOOL_SORT'
	WRITE(6,'(A)')' '
!
	STOP
	END
