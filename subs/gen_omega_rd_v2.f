!
! Read in collsion strengths from a file. The values are tabulated as a 
! function of temperature.
!
! OMEGA_SET is used to set the collison strength for those transitions in which
!    (1) No collison strength was contained in the file and
!    (2) has a zero oscillator strength.
!
! OMEGA_SCALE is used by SUBCOL_MULTI to scale all collsion strengths by a
! constant value.
!
	SUBROUTINE GEN_OMEGA_RD_V2(OMEGA_TABLE,T_TABLE,
	1                      ID_LOW,ID_UP,ID_INDX,
	1                      OMEGA_SCALE,OMEGA_SET,
	1                      LEVNAME,STAT_WT,NLEV,FILE_NAME,
	1                      NUM_TRANS,NUM_TVALS,
	1                      MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE)
	IMPLICIT NONE
!
! Altered 24-Apr-2015 : Code now checks "validity" of T values and collision strengths.
! Altered 25-Jan-2015 : Minor bug fix - could print out a wrong matching name.
! Altered 20-Dec-2014 : Code checks if non-matching level corresponds to a higher level not 
!                          included in the model atom (call CHL_COL_NAME).
! Altered 15-Dec-2013 : Code now only outputs error message on first call for each species.
! Altered 02-Sep-2012 : Outputs error messgae when matching transition not found.
!                         Checks ordering of levels (necessary for overlapping LS states).
!                         Fixed bug when data levels split, and program levels unsplit.
! Altered 20-Dec-2006 : Minor bug fix. Incorrect ray access after reading in
!                         collisonal ioization term. LST_NUP was being set to 0.
! Altered 20-Dec-2004 : Eror message now output if error occurs reading
!                         the actual collison values.
! Altered 02-Mar-1998 : Ientifying the connection betwen table and model 
!                         levels improved. Routine can now handle split and
!                         combined LS stars in both the model and input data.
! Altered 08-Feb-1998 : STRING length invreased from 132 to 500 for Fe2.
! Altered 19-Dec-1996 : Bugs fixed, and changed splitting approximation.
! Altered 07-Dec-1996 : Can read in multiplet Omega values when split
!                         levls are present. Spliting is approximate.
! Altered 29-May-1996 : Bug fix. Incorrect handling of collisional data for
!                         levels not included in the model atom.
! Altered 24-May-1996 : GEN_SEQ_OPEN replaced by F90 OPEN statement.
! Altered 01-Mar-1996 : Warning about invalid levels commented out.
!                       Warning was occuring just because of additional levels
!                       in data file.
! Created 13-Sep-1995 
!
	INTEGER NLEV
	CHARACTER*(*) LEVNAME(NLEV),FILE_NAME
	REAL*8 STAT_WT(NLEV)
	REAL*8 OMEGA_SCALE,OMEGA_SET
!
! ID_LOW:  ID of lower level of transition
! ID_UP:   ID of upper level of transition
! ID_INDX: Start location of tabulated values in table (OMEGA_TABLE)
!
! T_TABLE: Tempearure values (10^4 K) at which collsion strengths tabulated.
!
	INTEGER MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE
	INTEGER ID_LOW(MAX_TRANS)
	INTEGER ID_UP(MAX_TRANS)
	INTEGER ID_INDX(MAX_TRANS)    
	REAL*8 T_TABLE(MAX_TVALS)
	REAL*8 OMEGA_TABLE(MAX_TAB_SIZE)
!
	INTEGER NUM_TRANS,NUM_TVALS
!
! Local variables.
!
	INTEGER, PARAMETER :: MAX_NUM_NO_MATCH=20
	INTEGER NUM_NO_MATCH
	CHARACTER(LEN=40) NO_MATCH_NAME(MAX_NUM_NO_MATCH)

	REAL*8 COL_VEC(MAX_TVALS)
	REAL*8 GL_SUM,GU_SUM,NORM_FAC
	INTEGER NL_VEC(10)
	INTEGER NUP_VEC(10)
	INTEGER NL_CNT
	INTEGER NUP_CNT
	INTEGER LOC_INDX(NLEV,NLEV)
!
	INTEGER I		!Used for lower level.
	INTEGER K
	INTEGER LOOP
	INTEGER TRANS_CNT	!Used for reading in the dtata.
	INTEGER NL,NUP
	INTEGER MAX_LEV_SEP
	INTEGER L,LUIN,IOS
	INTEGER TRANS_INDX
	INTEGER TAB_INDX
	INTEGER LST_SPLIT_LEVEL
	INTEGER LST_NL
	INTEGER LST_NUP
	INTEGER LUER,ERROR_LU,ICHRLEN
	CHARACTER*80 STRIP_LD_BLANK
	EXTERNAL ERROR_LU,ICHRLEN
!
	CHARACTER*132 LOW_LEV,UP_LEV
	CHARACTER*500 STRING
	CHARACTER*60 LOCNAME(NLEV)
!
	LOGICAL, SAVE :: DO_WARNING_OUTPUT
	CHARACTER(LEN=80), SAVE :: FIRST_FILE=' '
!
	IF(FIRST_FILE .EQ. ' ')THEN
	  FIRST_FILE=FILE_NAME
	  DO_WARNING_OUTPUT=.TRUE.
	ELSE IF(FIRST_FILE .EQ. FILE_NAME)THEN
	  DO_WARNING_OUTPUT=.FALSE.
	END IF
	LUER=ERROR_LU()
	NUM_NO_MATCH=0
	NO_MATCH_NAME=' '
!
	DO NL=1,NLEV
	  LOCNAME(NL)=STRIP_LD_BLANK(LEVNAME(NL))
	END DO
	LST_SPLIT_LEVEL=0
	DO NL=1,NLEV
	  IF(INDEX(LOCNAME(NL),'[') .NE. 0)LST_SPLIT_LEVEL=NL
	END DO
	MAX_LEV_SEP=1
	DO NL=1,LST_SPLIT_LEVEL
	  K=INDEX(LOCNAME(NL),'[')
	  IF(K .NE. 0)THEN
	    DO NUP=NL+1,NLEV
	      IF(LOCNAME(NL)(1:K) .EQ. LOCNAME(NUP)(1:K))THEN
	        MAX_LEV_SEP=MAX(MAX_LEV_SEP,NUP-NL)
	      END IF
	    END DO
	  END IF
	END DO
	LOC_INDX(:,:)=0
	OMEGA_TABLE(:)=0.0D0
!
! Open file with collisional data.
!
	LUIN=7
	OPEN(LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',FILE_NAME
	  WRITE(LUER,*)'IOS=',IOS	
	  STOP
	END IF
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of transitions')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Numer of transitions string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
!
	READ(STRING,*)NUM_TRANS
	IF(NUM_TRANS .GT. MAX_TRANS)THEN
	  WRITE(LUER,*)'Error reading collisonal data from '//FILE_NAME
	  WRITE(LUER,*)'MAX_TRANS is too small:',MAX_TRANS
	  WRITE(LUER,*)'NUM_TRANS =',NUM_TRANS
	  STOP
	END IF
	STRING=' '
!
	DO WHILE( INDEX(STRING,'!Number of T values OMEGA tabulated at')  
	1               .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Numer of T values string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	READ(STRING,*)NUM_TVALS
	IF(NUM_TVALS .GT. MAX_TVALS)THEN
	  WRITE(LUER,*)'Error reading collisonal data from '//FILE_NAME
	  WRITE(LUER,*)'NU_TVALS too small'
	  STOP
	END IF
!         
! OMEGA_SCALE provides a means of scaling the OMEGA that have not been read in.
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Scaling factor for OMEGA') .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Scaling factor for Omega string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	READ(STRING,*)OMEGA_SCALE
!         
! OMEGA_SET provides a means to set OMEGA for transition for which
! no atomic data is availaable, and for which the oscilator strength
! is zero.
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Value for OMEGA if f=0') .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'!Value for OMEGA if f=0 string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	READ(STRING,*)OMEGA_SET
!
! We do this check here so that OMEGA_SCALE and OMEGA_SET are defined.
!
	IF(NUM_TRANS .EQ. 0)THEN
	  CLOSE(LUIN)
	  RETURN			!i.e use approximate formulae only.
	END IF
!
	STRING=' '
	DO WHILE( INDEX(STRING,'Transition\T')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',FILE_NAME
	    WRITE(LUER,*)'Transition\T string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END DO
	L=INDEX(STRING,'ion\T')
	READ(STRING(L+5:),*)(T_TABLE(I),I=1,NUM_TVALS )
!
	DO I=1,NUM_TVALS-1
	  IF(T_TABLE(1) .LE. 0.0D0 .OR. T_TABLE(I) .GE. T_TABLE(I+1))THEN
	    WRITE(LUER,*)'Error in T values in ',TRIM(FILE_NAME)
	    WRITE(LUER,*)'Values are not positive, or not monotonic increasing: I=',I
	    WRITE(LUER,*)T_TABLE(I),T_TABLE(I+1)
	    STOP
	  END IF
	END DO
!
	STRING=' '
	DO WHILE(ICHRLEN(STRING) .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO             
!
! TRANS_CNT is the index used to read in all collisional data. 
! TRANS_INDX is used to specify the number of collisional transitions that
! are used by the present model. It will differ from TRANS_CNT because
! collisional data may be available for levels not included in the model
! atom.
!
! NB: TRANS_INDX =< TRANS_CNT.
!
	TAB_INDX=1
	TRANS_INDX=0
	LST_NL=1
	LST_NUP=1
	DO TRANS_CNT=1,NUM_TRANS
!
! Read and determine correspondance of lower level to that of model atom.
!
! Get lower level name.
!
	  I=1
	  DO WHILE (STRING(I:I) .EQ. ' ')
	    I=I+1
	  END DO
	  L=INDEX(STRING,'-')          
	  LOW_LEV=STRING(I:L-1)
!
! Determine corespondance of lower level. NB: The MOD function allows us
! to loop over all levels beggining from LST_NL. Thus the SEARCH order in
! INDEX space is
!              LST_NL, LST_NL+1, .... NLEV, 1, 2, ... , LST_NL-1
!
	  K=INDEX(LOW_LEV,'[')
	  NL_CNT=0
	  GL_SUM=0.0D0
	  IF(K .NE. 0)THEN
!
! Only a single unique match is possible.
!
	    DO LOOP=1,NLEV
	      NL=MOD(LOOP+LST_NL-2,NLEV)+1
	      IF( LOW_LEV .EQ. LOCNAME(NL) .OR.
	1            LOW_LEV(1:K-1).EQ. LOCNAME(NL) )THEN
	        NL_CNT=NL_CNT+1
	        NL_VEC(NL_CNT)=NL
	        GL_SUM=GL_SUM+STAT_WT(NL)
	        GOTO 100
	      END IF
	    END DO
	    GOTO 500			!No match
	  ELSE
!
! As the input level has no [], multiple matches are possible if model
! atom levels are split.
!
	    DO LOOP=1,NLEV
	      NL=MOD(LOOP+LST_NL-2,NLEV)+1
	      K=INDEX(LOCNAME(NL),'[')
	      IF(K .EQ. 0)THEN
	        IF(LOW_LEV .EQ. LOCNAME(NL))THEN
	          NL_CNT=NL_CNT+1
                  NL_VEC(NL_CNT)=NL
	          GL_SUM=GL_SUM+STAT_WT(NL)
	          GOTO 100
	        END IF
	      ELSE
	        IF(LOW_LEV .EQ. LOCNAME(NL)(1:K-1))THEN
	          DO I=MAX(1,NL-MAX_LEV_SEP),MIN(NLEV,NL+MAX_LEV_SEP)
	            IF(LOW_LEV .EQ. LOCNAME(I)(1:K-1))THEN
	              NL_CNT=NL_CNT+1
                      NL_VEC(NL_CNT)=I
	              GL_SUM=GL_SUM+STAT_WT(I)
	            END IF
	          END DO
	          GOTO 100
	        END IF
	      END IF
	    END DO
	    GOTO 500		!No match
	  END IF
100	  CONTINUE
!
! Read and determine correspondance of upper level to that of model atom.
!
! Get upper level name.
!
	  STRING(1:)=STRING(L+1:)
	  I=1
	  DO WHILE (STRING(I:I) .EQ. ' ')
	    I=I+1                 
	  END DO
	  STRING(1:)=STRING(I:)
	  L=INDEX(STRING,'  ')
	  UP_LEV=STRING(1:L-1)
!
! Determine upper level corespondance.
!
	  K=INDEX(UP_LEV,'[')
	  NUP_CNT=0
	  GU_SUM=0.0D0
	  IF(NL_VEC(1) .NE. LST_NL)LST_NUP=NL_VEC(1)
	  IF(UP_LEV .EQ. 'I')THEN
!
! Collisional ionization terms.
!
	    GU_SUM=1.0D0
	    GOTO 200
	  ELSE IF(K .NE. 0)THEN
	    DO LOOP=1,NLEV
	      NUP=MOD(LOOP+LST_NUP-2,NLEV)+1
	      IF(UP_LEV .EQ. LOCNAME(NUP) .OR.
	1            UP_LEV(1:K-1) .EQ. LOCNAME(NUP))THEN
	        NUP_CNT=NUP_CNT+1
	        NUP_VEC(NUP_CNT)=NUP
	        GU_SUM=GU_SUM+STAT_WT(NUP)
	        GOTO 200
	      END IF
	    END DO
	    GOTO 500				!No match
	  ELSE
	    DO LOOP=1,NLEV
	      NUP=MOD(LOOP+LST_NUP-2,NLEV)+1
	      IF(LST_NUP .LT. 1)THEN
	        WRITE(LUER,*)LST_NUP,FILE_NAME
	        WRITE(LUER,*)LOW_LEV,UP_LEV
	      END IF
	      K=INDEX(LOCNAME(NUP),'[')
	      IF(K .EQ. 0)THEN
	        IF(UP_LEV .EQ. LOCNAME(NUP))THEN
	          NUP_CNT=NUP_CNT+1
                  NUP_VEC(NUP_CNT)=NUP
	          GU_SUM=GU_SUM+STAT_WT(NUP)
	          GOTO 200
	        END IF
	      ELSE
!
! Several matches possible
!
	        IF(UP_LEV .EQ. LOCNAME(NUP)(1:K-1))THEN
	          DO I=MAX(1,NUP-MAX_LEV_SEP),MIN(NLEV,NUP+MAX_LEV_SEP)
	            IF(UP_LEV .EQ. LOCNAME(I)(1:K-1))THEN
	              NUP_CNT=NUP_CNT+1
                      NUP_VEC(NUP_CNT)=I
	              GU_SUM=GU_SUM+STAT_WT(I)
	            END IF
	          END DO
	          GOTO 200
	        END IF
	      END IF
	    END DO
	    GOTO 500		!No match
	  END IF
!
200	  NORM_FAC=1.0D0/GL_SUM/GU_SUM
	  LST_NL=NL_VEC(1)
	  IF(UP_LEV .NE. 'I')LST_NUP=NUP_VEC(1)
	  STRING(1:)=STRING(L:)
!
! The multiplet collision strengths are converted to individual collision 
! strengths using:
!
! Omega(i,j)= g(i)g(j) OMEGA(MULITPLET) / G(L)/G(U)
!
! This formula automatically recovers the correct result that
!
! OMEGA = SUM(i) SUM(j)  Omega(i,j)
!
	  READ(STRING,*,IOSTAT=IOS)(COL_VEC(I),I=1,NUM_TVALS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in GEN_OMEGA_RD_V2'  
	    WRITE(LUER,*)'Error reading collisional data from',FILE_NAME
	    WRITE(LUER,*)TRIM(STRING)
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  IF( MINVAL(COL_VEC(1:NUM_TVALS)) .LE. 0.0D0 )THEN
	    WRITE(LUER,*)'Error in GEN_OMEGA_RD_V2'  
	    WRITE(LUER,*)'Error in collisional data for',FILE_NAME
	    WRITE(LUER,*)'Tabulated values are -ve or zero'
	    WRITE(LUER,*)'UP_LEV=',TRIM(UP_LEV)
	    WRITE(LUER,*)'LOW_LEV=',TRIM(LOW_LEV)
	    STOP
	  END IF
!
	  COL_VEC(1:NUM_TVALS)=COL_VEC(1:NUM_TVALS)*NORM_FAC
	  IF(UP_LEV .EQ. 'I')THEN
	    DO NL=1,NL_CNT
	      TRANS_INDX=TRANS_INDX+1
	      ID_LOW(TRANS_INDX)=NL_VEC(NL)
	      ID_UP(TRANS_INDX)=NL_VEC(NL)
	      ID_INDX(TRANS_INDX)=TAB_INDX
	      OMEGA_TABLE(TAB_INDX:TAB_INDX+NUM_TVALS-1)=
	1         OMEGA_TABLE(TAB_INDX:TAB_INDX+NUM_TVALS-1)+
	1            COL_VEC(1:NUM_TVALS)*STAT_WT(NL_VEC(NL))
	      TAB_INDX=TAB_INDX+NUM_TVALS
	    END DO
	  ELSE
!
! If the first IF statement is true, collisional data is split but not the level in model atom.
! Unlikely to happen. Since the J states of individual states may overlap, we need to check
! the ordering.
!
	    DO NL=1,NL_CNT
	      DO NUP=1,NUP_CNT
	        IF(NL_VEC(NL) .EQ. NUP_VEC(NUP))THEN
	        ELSE 
	          IF(LOC_INDX(NL_VEC(NL),NUP_VEC(NUP)) .EQ. 0)THEN
	            TRANS_INDX=TRANS_INDX+1
	            IF(NL_VEC(NL) .LT. NUP_VEC(NUP))THEN
	              ID_LOW(TRANS_INDX)=NL_VEC(NL)
	              ID_UP(TRANS_INDX)=NUP_VEC(NUP)
	            ELSE
	              ID_LOW(TRANS_INDX)=NUP_VEC(NL)
	              ID_UP(TRANS_INDX)=NL_VEC(NUP)
	            END IF
	            ID_INDX(TRANS_INDX)=TAB_INDX
	            I=TAB_INDX
	            TAB_INDX=TAB_INDX+NUM_TVALS
	            LOC_INDX(NL_VEC(NL),NUP_VEC(NUP))=I
	          ELSE
	            I=LOC_INDX(NL_VEC(NL),NUP_VEC(NUP))
	          END IF
	          OMEGA_TABLE(I:I+NUM_TVALS-1)=OMEGA_TABLE(I:I+NUM_TVALS-1)+
	1            COL_VEC(1:NUM_TVALS)*STAT_WT(NUP_VEC(NUP))*
	1                              STAT_WT(NL_VEC(NL))
	        END IF
	      END DO
	    END DO
	  END IF
!
500	  CONTINUE
!
	  IF(NL_CNT .EQ. 0)THEN
	    K=0
	    DO I=1,NUM_NO_MATCH
	      IF(TRIM(LOW_LEV) .EQ. TRIM(NO_MATCH_NAME(I)))K=I
	    END DO
	    IF(K .EQ. 0)THEN
	      NUM_NO_MATCH=MIN(NUM_NO_MATCH+1,MAX_NUM_NO_MATCH)
	      K=NUM_NO_MATCH
	      NO_MATCH_NAME(K)=TRIM(LOW_LEV)
	    END IF
	  END IF
!
! We only bother with the upper level if we successfully found a match for the lower level.
!
	  IF(NL_CNT .NE. 0 .AND. NUP_CNT .EQ. 0 .AND. TRIM(UP_LEV) .NE. 'I')THEN
	    K=0
	    DO I=1,NUM_NO_MATCH
	      IF(TRIM(UP_LEV) .EQ. TRIM(NO_MATCH_NAME(I)))K=I
	    END DO
	    IF(K .EQ. 0)THEN
	      NUM_NO_MATCH=MIN(NUM_NO_MATCH+1,MAX_NUM_NO_MATCH)
	      K=NUM_NO_MATCH
	      NO_MATCH_NAME(K)=TRIM(UP_LEV)
	    END IF
	  END IF
!
! Read in the next collisional data set. If the record is blank, we
! read in the NEXT record, as BLANK records are not included in the NUM_TRANS
! count.
!
	  STRING=' '
	  DO WHILE(ICHRLEN(STRING) .EQ. 0 .AND. TRANS_CNT .LT. NUM_TRANS)
	    READ(LUIN,'(A)')STRING
	  END DO
!
	END DO
!
! We now check whether the level names without match correspond to levels not
! included in the model calculations.
!
	CLOSE(LUIN)
	IF(DO_WARNING_OUTPUT)THEN
	  CALL CHK_COL_NAME(NO_MATCH_NAME,NUM_NO_MATCH,NLEV,LUIN,LUER,FILE_NAME)
	END IF
!
! NUM_TRANS is now set equal to the actual number of transitions read.
!
	NUM_TRANS=TRANS_INDX
!
	RETURN
	END

	FUNCTION STRIP_LD_BLANK(STRING)
	IMPLICIT NONE
C
	INTEGER I
	CHARACTER*(*) STRING
	CHARACTER*80 STRIP_LD_BLANK
	I=1
	DO WHILE(I .LT. LEN(STRING) .AND. STRING(I:I) .EQ. ' ')
	 I=I+1
	END DO
	STRIP_LD_BLANK=STRING(I:)
	RETURN
	END
