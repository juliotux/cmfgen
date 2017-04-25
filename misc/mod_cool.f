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
! Altered: 18-Oct-2017  Fixed estimate of cooling time.
! Altered: 01-Mar-2016  Option to omit advection terms when not included in the model [25-Feb-2016].
! Altered: 17-Feb-2015  Improved estimate of cooling time by including ATOM_DENSITY.
! Altered:              Estimate of cooling time outout to GENCOOL_SORT
! Altered: 12-Mar-2014: Added read/ouput of non-thermal cooling.
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
	REAL*8, ALLOCATABLE :: TOTAL_RATE(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ATOM_DENSITY(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: COOLING_TIME(:)
!
	INTEGER ND
	INTEGER NV
	INTEGER I,J,K,ID
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	REAL*8 T1
	LOGICAL FILE_OPEN
	LOGICAL ONLY_INCLUDED_TERMS
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
	ONLY_INCLUDED_TERMS=.TRUE.
	CALL GEN_IN(ONLY_INCLUDED_TERMS,'Only output terms that are included?')
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	ALLOCATE (SUM(ND)); SUM=0.0D0
	ALLOCATE (BOUND(ND));  BOUND=0.0D0
	ALLOCATE (TOTAL_RATE(ND)); TOTAL_RATE=0.0D0
	ALLOCATE (COOLING_TIME(ND)); COOLING_TIME=0.0D0
	ALLOCATE (ATOM_DENSITY(ND)); ATOM_DENSITY=0.0D0
	ALLOCATE (ED(ND)); ED=0.0D0
	ALLOCATE (T(ND)); T=0.0D0
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
	     IF(J .EQ. 2)READ(STRING,*)(T(K),K=ID,MIN(ID+9,ND))
	     IF(J .EQ. 3)READ(STRING,*)(ED(K),K=ID,MIN(ID+9,ND))
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
	        CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'COL ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Free-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'FF ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Non-thermal') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'NT ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'K-shell') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'XKS ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'V term') .NE. 0)THEN
	        IF(ONLY_INCLUDED_TERMS .AND. INDEX(STRING,'Not Incl') .NE. 0)THEN
	          READ(20,'(A)')STRING
	        ELSE
	          TMP_STR=STRING
	          TMP_STR=ADJUSTL(TMP_STR)
	          K=INDEX(TMP_STR,' ')
	          READ(20,'(A)')STRING
	          CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	          WRITE(21,'(A)')' '
	          WRITE(21,'(A,T12,A)')'AC.R(V).',TRIM(STRING)
	        END IF
	      ELSE IF( INDEX(STRING,'dTdR term') .NE. 0)THEN
	        IF(ONLY_INCLUDED_TERMS .AND. INDEX(STRING,'Not Incl') .NE. 0)THEN
	          READ(20,'(A)')STRING
	        ELSE
	          TMP_STR=STRING
	          TMP_STR=ADJUSTL(TMP_STR)
	          K=INDEX(TMP_STR,' ')
	          READ(20,'(A)')STRING
	          CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	          WRITE(21,'(A,T12,A)')'AC.R(dT).',TRIM(STRING)
	        END IF
	      ELSE IF( INDEX(STRING,'decay') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
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
	            CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
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
	            CALL SUM_RATES(TOTAL_RATE,STRING,ID,ND)
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
	WRITE(6,'(A)')' Cooling time is approximate and not valid at high densities'
!
! We have summed all rates without regard to sign. If in equilibrium, the cooling
! rate must be half this rate. This will not work at high densities where the is
! cancellation in heating/cooling rates for individual species/processes.
!
! For simplicty, we assume the energy per species is 1.5kT, and that the number
! of electrons is the same as the number of ions when the atom density is unavailable.
!
	TOTAL_RATE=TOTAL_RATE*0.5D0
	I=7
	CALL RD_SING_VEC_RVTJ(ATOM_DENSITY,ND,'Atom Density','RVTJ',I,IOS)
	IF(IOS .NE. 0)THEN
	  COOLING_TIME=3.0D0*1.3806D-12*T*ED/TOTAL_RATE
	ELSE
	  COOLING_TIME=1.5D0*1.3806D-12*T*(ED+ATOM_DENSITY)/TOTAL_RATE
	END IF
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
	  WRITE(21,'(A,10ES12.4)')'Cool time(s)',(COOLING_TIME(K),K=ID,MIN(ID+9,ND))
	  
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
!
!
!
	SUBROUTINE SUM_RATES(TOTAL_RATE,STRING,ID,ND)
	INTEGER ID, ND
	REAL*8 TOTAL_RATE(ND)
	REAL*8 TEMP_VEC(ND)
	CHARACTER(LEN=*) STRING
!
	INTEGER I
!
	READ(STRING,*)(TEMP_VEC(I),I=ID,MIN(ID+9,ND))
	DO I=ID,MIN(ID+9,ND)
	  TOTAL_RATE(I)=TOTAL_RATE(I)+ABS(TEMP_VEC(I))
	END DO
!
	RETURN
	END

