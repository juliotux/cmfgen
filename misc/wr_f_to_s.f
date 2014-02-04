C
C Routine to output link between full levels in a model atom and their
C corresponding super levels.
C
C Names and energy levels are read from a file with the same format as
C the oscillator file (hence the oscillator file can be used).
C
	PROGRAM WR_F_TO_S
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	IMPLICIT NONE
C
C Altered 25-Sec-2011 : NAME AND LS_NAME SET TO *40
C Altered 13-Dec-2009 : LOWN option installed (comments added 13-Jan-2010). 
C Altered 21-Apr-2008 : SEQ_WR installed; Improved handling of INT_SEQ.
C Altered 23-Jun-2005 : FIX_DI option installed for WR_DC
C Altered 25-Oct-2002 : CL option changed.
C                         Level without SL designation is given one.
C Altered 03-Nov-2000 : E%LS option installed. Designed to give
C                         fine spacing for low lying levels, but coarser
C                         spacing for upper levels.
C Altered 22-Oct-1999 : SPLIT and SEP options installed.
C                       Bug in option RD_DC fixed.
C Altered 19-May-1998 : T24 changes to T30 in format statements.
C                       SL_WR option installed.
C Altered 12-May-1997 : ELS option installed, and SP option cleaned.
C Altered 08-May-1997 : NAME AND LS_NAME SET TO *30
C
	INTEGER, PARAMETER :: N_MAX=8000
C
	CHARACTER(LEN=40) NAME(N_MAX)
	CHARACTER(LEN=40) LS_NAME(N_MAX)
	CHARACTER*1 LEV_SPIN(N_MAX)
	CHARACTER*1 LEV_PARITY(N_MAX)
C
	INTEGER F_TO_S(N_MAX)
	INTEGER INT_SEQ(N_MAX)
	REAL*8 FEDGE(N_MAX)
	REAL*8 ENERGY(N_MAX)
	REAL*8 G(N_MAX)
	REAL*8 LAM_EDGE(N_MAX)
	REAL*8 E_STRT(N_MAX)
	LOGICAL DONE_LEV(N_MAX)
C
	REAL*8 EDGE_SUM(N_MAX)
	REAL*8 G_SUM(N_MAX)
C
	CHARACTER(LEN=40) TERM_NAME(N_MAX)
	INTEGER TERM_F_TO_S(N_MAX)
	INTEGER NTERM
	LOGICAL PRES
C
C For SUPER-LEVEL atom.
C
	REAL*8 LS_EMIN(N_MAX)
	REAL*8 LS_EMAX(N_MAX)
	REAL*8 LS_EMID(N_MAX)
	LOGICAL LS_DONE(N_MAX)
	CHARACTER*1 LS_SPIN(N_MAX)
	CHARACTER*1 LS_PARITY(N_MAX)
	LOGICAL CHANGE_FS(N_MAX)
	LOGICAL DO_CHANGE
	REAL*8 ACC
	INTEGER FS_SAV
C
	REAL*8 DC(N_MAX,3),ED(3),TEMP(3)
	REAL*8 DI,RVAL,RSTAR,RLUM
	REAL*8 T_EXCITE,G_GS,G_ION
	INTEGER NLEV_RD,ND_RD
	LOGICAL WRITE_DC
C
	REAL*8 ION_EN
	REAL*8 ZION
	CHARACTER*20 EN_DATE
	INTEGER NLEV
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	DOUBLE PRECISION CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	CHARACTER*80 FILENAME
	LOGICAL FILE_OPEN
	LOGICAL FILE_PRES
C
	INTEGER IZERO
	PARAMETER (IZERO=0)
C
	INTEGER I,J,K,L,IOS,ID
	INTEGER CNT,NEW_CNT,N_LS_TERMS
	INTEGER I1,K1
	INTEGER MAX_NAME_LNGTH
	REAL*8 T1,T2
	REAL*8 DEL_E_CM
	REAL*8 DEL_E_CUR
	REAL*8 FRAC_DEL_E
	REAL*8 OLD_F_TO_S
	REAL*8 R_IR,V1,CL_FAC
C
	LOGICAL CHK_PARITY
	LOGICAL CHK_SPIN
	LOGICAL SORT_LEVS
	LOGICAL OKAY_TO_COMBINE
	LOGICAL FIX_DI
C
	CHARACTER TIME*24
	CHARACTER STRING*80
	CHARACTER ANS*1			!Used for halting LI and HE options
	CHARACTER(LEN=10) TMP_STR
C
	LOGICAL HEAD
	LOGICAL L_TRUE,L_FALSE
	DATA L_TRUE/.TRUE./
	DATA L_FALSE/.FALSE./
C
C 
C
C USR_OPTION variables
C
	CHARACTER MAIN_OPT_STR*80		!Contains full option string
	CHARACTER X*10		!Option
	CHARACTER*120 DEFAULT	!String to be used for default values
	CHARACTER*120 DESCRIPTION
C
	CHARACTER*30 UC
	EXTERNAL UC
	INTEGER ICHRLEN
	CHARACTER*1 PARITY,SPIN
	EXTERNAL ICHRLEN,PARITY
C
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6		!For terminal IO              
C
	INTEGER, PARAMETER :: LUIN=30			!File IO
	INTEGER, PARAMETER :: LUOUT=40
	INTEGER, PARAMETER :: LUHEAD=51
C
C 
C
C Set constants.
C
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	INT_SEQ(:)=0
C
C  
C Read in Level Names and Energies from file containing oscillator
C strengths.
C
C Header information is read and stored in the file HEAD_INFO. This can be 
C output to the F_TO_S link file.
C
	IOS=100
	DO WHILE(IOS .NE. 0)
	  CALL GEN_ASCI_OPEN(LUHEAD,'HEAD_INFO','UNKNOWN',' ',
	1                          'WRITE',IZERO,IOS)
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Oscillator file')
	  CALL RD_ENERGY(NAME,G,ENERGY,FEDGE,NLEV,N_MAX,
	1       ION_EN,ZION,EN_DATE,FILENAME,LUIN,LUHEAD,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error occurred reading Oscillator file: try again'
	  END IF
	  CLOSE(LUHEAD)
	END DO
C
C Strip J values from names, and determine parity and spin.
C
	MAX_NAME_LNGTH=0
	DO I=1,NLEV
	  MAX_NAME_LNGTH=MAX(MAX_NAME_LNGTH,LEN_TRIM(NAME(I)))
	  J=INDEX(NAME(I),'[')
	  IF(J .EQ. 0)THEN
	    LS_NAME(I)=NAME(I)
	  ELSE
	    LS_NAME(I)=NAME(I)(1:J-1)
	  END IF
	  LEV_PARITY(I)=PARITY(LS_NAME(I))
	  LEV_SPIN(I)=SPIN(LS_NAME(I))
	END DO
C
C Determine number of LS terms in FULL model atom.
C
	N_LS_TERMS=1
	DO I=2,NLEV
	  J=1
	  DO WHILE(LS_NAME(I) .NE. LS_NAME(J) .AND. J .LE. I-1)
	    J=J+1
	  END DO
	  IF(J .EQ. I)N_LS_TERMS=N_LS_TERMS+1
	END DO
C
C
C
C
C This message will only be printed once
C
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file '//
	1    'main_option.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to '//
	1    'write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to '//
	1    'write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
C
	DESCRIPTION='Enter main option'
C
C This call resets the .sve algorithm.  Specifically it sets the next
C input answer to be a main option, and all subsequent inputs to be
C sub-options.  
C
3	CONTINUE
	CALL SVE_FILE('RESET')
C
	MAIN_OPT_STR='  '
	DEFAULT=' '
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
C
C   If the main option begins with a '.', a previously
C   written .sve file is read.
C
C   If the main option begins with a '#', a previously
C   written .box file is read.
C
C   If sve= is appended to the end of this main option, a new .sve file
C   is opened with the given name and the main option and all subsequent
C   sub-options are written to this file.
C
C   If box= is input then a .box file is created, which contains the name
C   of several .sve files to process.
C
C   If only a main option is given, the option and subsequent sub-options
C   are saved in a file called 'main option.sve'.  All following main
C   options are saved in separate files.
C
!
! Remove variable changes from main option.
!
        I=INDEX(MAIN_OPT_STR,'(')
        IF(I .EQ. 0)THEN
          X=UC(TRIM(MAIN_OPT_STR))
        ELSE
          X=UC(MAIN_OPT_STR(1:I-1))     !Remove line variables.
        END IF
!
! Remove possile file names etc.
!          
        I=INDEX(X,' ')
        IF(I .EQ. 0)THEN
          X=UC(TRIM(X))
        ELSE
          X=UC(X(1:I-1))        !Remove file names
        END IF
C
C                           
C Group all terms belonging to the same LS multiplet.
C               
	IF(X(1:2) .EQ. 'LS')THEN	!LS coupling
	  F_TO_S(:)=0
	  F_TO_S(1)=1
	  CNT=1
	  DO I=2,NLEV
	    J=INDEX(NAME(I),'[')
	    IF(J .EQ. 0)THEN
	      CNT=CNT+1
	      F_TO_S(I)=CNT
	    ELSE
	      J=1
	      DO WHILE(F_TO_S(I) .EQ. 0 .AND. J .LE. I-1)
	        IF( LS_NAME(I) .EQ. LS_NAME(J) )THEN
                  F_TO_S(I)=F_TO_S(J)
	        END IF
	        J=J+1
	      END DO
	      IF(F_TO_S(I) .EQ. 0)THEN
	        CNT=CNT+1
	        F_TO_S(I)=CNT
	      END IF
	    END IF
	  END DO
C
	  WRITE(T_OUT,*)'Number of levels in full atom is    ',NLEV
	  WRITE(T_OUT,*)'Number of LS terms in full atom is  ',N_LS_TERMS
	  WRITE(T_OUT,*)'Number of levels in SUPER atom is   ',CNT
C
C Group all terms belonging to the same LS multiplet.
C               
	ELSE IF(X(1:3) .EQ. 'TLS')THEN	!LS coupling
	  F_TO_S(:)=0
	  F_TO_S(1)=1
	  CNT=1
	  DO I=2,NLEV
	    J=INDEX(NAME(I),'[')
	    IF(J .EQ. 0)THEN
	      CNT=CNT+1
	      F_TO_S(I)=CNT
	    ELSE
	      J=1
	      DO WHILE(F_TO_S(I) .EQ. 0 .AND. J .LE. I-1)
	        IF( LS_NAME(I) .EQ. LS_NAME(J) )THEN
                  F_TO_S(I)=F_TO_S(J)
	        END IF
	        J=J+1
	      END DO
	      IF(F_TO_S(I) .EQ. 0)THEN
	        CNT=CNT+1
	        F_TO_S(I)=CNT
	      END IF
	    END IF
	  END DO
!
	  DO I=1,NLEV-1
	   IF(INDEX(NAME(I),'[') .EQ. 0 .AND. INT_SEQ(I) .EQ. 0)THEN
	     L=LEN_TRIM(NAME(I))
	     PRES=.TRUE.
	     DO J=I+1,NLEV
	       IF(INDEX(NAME(J),'[') .EQ. 0)THEN
	           K=LEN_TRIM(NAME(J))
	         IF(NAME(J)(K-4:K) .EQ. NAME(I)(L-4:L))THEN
	           IF(PRES)THEN
	             ID=I
	             INT_SEQ(I)=ID 
	             INT_SEQ(J)=ID
	             PRES=.FALSE.
	             CNT=J
	           ELSE
	             F_TO_S(J)=F_TO_S(CNT)
	             INT_SEQ(J)=ID
	           END IF
                 END IF
	       END IF
	     END DO
	   END IF
	  END DO
	  WRITE(6,*)'Now do a clean'
C
	ELSE IF(X(1:3) .EQ. 'AVE')THEN
          G_SUM(:)=0.0D0
          EDGE_SUM(:)=0.0D0
          DO J=1,NLEV
            I=F_TO_S(J)
            G_SUM(I)=G_SUM(I)+G(J)
            EDGE_SUM(I)=EDGE_SUM(I)+ENERGY(J)*G(J)
          END DO
	  I=N_LS_TERMS
	  EDGE_SUM(1:I)=EDGE_SUM(1:I)/G_SUM(1:I)
	  EDGE_SUM(2:I)=EDGE_SUM(2:I)-EDGE_SUM(1)
	  DO J=1,NLEV
	    I=F_TO_S(J)
	    IF(EDGE_SUM(I) .NE. -999999.0D0)THEN
	      WRITE(31,'(A,T20,F4.0,3X,F12.3)')
	1                TRIM(LS_NAME(J)),G_SUM(I),EDGE_SUM(I)
	      EDGE_SUM(I)=-999999.0D0
	    END IF
	 END DO
C
C Group all terms belonging to the same LS multiplet. Additional grouping
C is done by combining terms that are within DEL_E_CM in energy of the 
C lowest term of the SUPER level. The defaults is to group together 
C only terms that have the same parity and spin, although this can be 
C changed using HIDDEN options.
C
C NB: In this option all levels belonging to the same multiplet are
C grouped together, independent of their energy separation.
C
	ELSE IF(X(1:3) .EQ. 'ELS' .OR. X(1:4) .EQ. 'E%LS')THEN
C
C We first do the LS grouping.
C
	  F_TO_S(:)=0
	  F_TO_S(1)=1
	  CNT=1
	  DO I=2,NLEV
	    J=INDEX(NAME(I),'[')
	    IF(J .EQ. 0)THEN
	      CNT=CNT+1
	      F_TO_S(I)=CNT
	    ELSE
	      J=1
	      DO WHILE(F_TO_S(I) .EQ. 0 .AND. J .LE. I-1)
	        IF( LS_NAME(I) .EQ. LS_NAME(J) )THEN
                  F_TO_S(I)=F_TO_S(J)
	        END IF
	        J=J+1
	      END DO
	      IF(F_TO_S(I) .EQ. 0)THEN
	        CNT=CNT+1
	        F_TO_S(I)=CNT
	      END IF
	    END IF
	  END DO
C
	  LS_DONE(:)=.FALSE.
	  LS_EMIN(:)=1.0E+36
	  LS_EMAX(:)=-1.0E+36
	  DO I=1,NLEV
	    LS_EMIN(F_TO_S(I))=MIN(LS_EMIN(F_TO_S(I)),ENERGY(I))
	    LS_EMAX(F_TO_S(I))=MAX(LS_EMAX(F_TO_S(I)),ENERGY(I))
	    LS_SPIN(F_TO_S(I))=LEV_SPIN(I)
	    LS_PARITY(F_TO_S(I))=LEV_PARITY(I)
	  END DO
	  DO I=1,NLEV
	    LS_EMID(I)=0.5D0*(LS_EMIN(I)+LS_EMAX(I))
	  END DO
C
	  IF(X(1:3) .EQ. 'ELS')THEN
	    CALL USR_OPTION(DEL_E_CM,'DEL_E','1000.0D0',
	1         'Maximum energy diff in cm^{-1}')
	  ELSE
	    CALL USR_OPTION(FRAC_DEL_E,'%E','10.0D0',
	1         'Percentage difference in excitation energy')
	  END IF
	  FRAC_DEL_E=0.01D0*FRAC_DEL_E
	  CALL USR_HIDDEN(CHK_PARITY,'CHK_P','T','Check parity?')
	  CALL USR_HIDDEN(CHK_SPIN,'CHK_S','T','Check spin?')
C
C Now do the grouping according to energy.
C
	  NEW_CNT=CNT
	  DO I=1,CNT-1
	    IF(X(1:3) .EQ. 'ELS')THEN
	      DEL_E_CUR=DEL_E_CM
	    ELSE
	      DEL_E_CUR=LS_EMID(I)*FRAC_DEL_E
	    END IF
	    IF(.NOT. LS_DONE(I))THEN
	      J=I+1     
	      DO WHILE(LS_EMIN(J) .LT. LS_EMIN(I)+DEL_E_CUR)
	        OKAY_TO_COMBINE=.TRUE.
	        IF(CHK_PARITY .AND. LS_PARITY(I) .NE. LS_PARITY(J))THEN
	          OKAY_TO_COMBINE=.FALSE.
	        END IF
	        IF(CHK_SPIN .AND. LS_SPIN(I) .NE. LS_SPIN(J))THEN
	          OKAY_TO_COMBINE=.FALSE.
	        END IF
	        IF(OKAY_TO_COMBINE .AND. .NOT. LS_DONE(J))THEN
	          IF( LS_EMAX(J)-LS_EMIN(I) .LT. DEL_E_CUR)THEN
		    DO K=1,NLEV
	              IF(F_TO_S(K) .EQ. J)F_TO_S(K)=I
	            END DO
	            NEW_CNT=NEW_CNT-1
	            LS_DONE(J)=.TRUE.
	          END IF
	        END IF
	        J=J+1
	        IF(J .GT. NLEV)EXIT
	      END DO
	    END IF
	  END DO
	  CNT=NEW_CNT
C
C The next 12 lines of instructions are identical to the 'CL' option.
c They ensure super levels are numbered sequentially.
C
	  ID=0
	  F_TO_S(1:NLEV)=-F_TO_S(1:NLEV)
	  DO I=1,NLEV
	   IF(F_TO_S(I) .LT. 0)THEN
	      ID=ID+1
	      OLD_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	      DO J=I+1,NLEV
	        IF(F_TO_S(J) .EQ. OLD_F_TO_S)F_TO_S(J)=ID
	      END DO
	    ELSE IF(F_TO_S(I) .EQ. 0)THEN
	      ID=ID+1
	      OLD_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	    END IF
	  END DO
C
	  WRITE(T_OUT,*)'Number of levels in full atom is    ',NLEV
	  WRITE(T_OUT,*)'Number of LS terms in full atom is  ',N_LS_TERMS
	  WRITE(T_OUT,*)'Number of levels in SUPER atom is   ',ID
!
! This option allows the N lowest SL's to be automatically spit into 
! individual levels. This option allows a direct comparision with a similar option
! in CMFGEN specified in the MODEL_SPEC file. Usefule for identifying levels
! with those in CMFGEN.
!
	ELSE IF(X(1:4) .EQ. 'LOWN')THEN
	  CALL USR_OPTION(K,'N','10','How many SL''s to split')
	  DO I=1,NLEV
	    IF(F_TO_S(I) .LE. K)THEN
	       F_TO_S(I)=I
	    ELSE
	      F_TO_S(I)=F_TO_S(I)+5000
	    END IF
	  END DO
	  WRITE(6,*)'Now execute CLN command'
!
! Option allows the user to improve the links in an existing link list.
! New SUPER levels will be created when the energy separation is
! greater than DEL_E_C, or when the departure coefficents of the levels
! (as determined by a full LS calculation) differ by more than ACC %.
! LS multiplets always belong to a single super level.
!
	ELSE IF(X(1:5) .EQ. 'SPLIT')THEN
!
	  CALL USR_OPTION(ACC,'ACC','20.0D0','Percentage accuracy')
	  CALL USR_OPTION(DEL_E_CM,'DEL_E','1000.0D0',
	1               'Maximum energy diff in cm^{-1}')
	  IF(DEL_E_CM .LT. 0)DEL_E_CM=ENERGY(NLEV)
	  LS_DONE(1:NLEV)=.FALSE.
	  CHANGE_FS(1:NLEV)=.FALSE.
	  DO I=1,NLEV
	    FS_SAV=F_TO_S(I)
	    IF(CHANGE_FS(I))F_TO_S(I)=NLEV*F_TO_S(I)+I
	    DO J=I+1,NLEV
	      IF(.NOT. LS_DONE(J) .AND. FS_SAV .EQ. F_TO_S(J))THEN
	        T1=0.5D0*(DC(J,1)+DC(I,1))
	        IF(T1 .EQ. 0)T1=1
	        DO_CHANGE=.FALSE.
	        IF(CHANGE_FS(I) .AND. CHANGE_FS(J))DO_CHANGE=.TRUE.
	        IF(.NOT. CHANGE_FS(I) .AND. .NOT. CHANGE_FS(J))DO_CHANGE=.TRUE.
	        IF( ABS(DC(J,1)-DC(I,1))/T1 .LT. 0.01*ACC .AND.
	1              ENERGY(J)-ENERGY(I) .LT. DEL_E_CM .AND. DO_CHANGE)THEN
	          F_TO_S(J)=F_TO_S(I)
	          CHANGE_FS(J)=.FALSE.
	          LS_DONE(J)=.TRUE.
	        ELSE IF(LS_NAME(I) .EQ. LS_NAME(J))THEN
	          F_TO_S(J)=F_TO_S(I)
	          CHANGE_FS(J)=.FALSE.
	          LS_DONE(J)=.TRUE.
	        ELSE
  	          CHANGE_FS(J)=.TRUE.
	        END IF
	      END IF
	    END DO
	  END DO
!
	  WRITE(6,*)'You should now execute the CL option'
C
C *****************************************************************************
C *****************************************************************************
C
C Group terms according to their separation in energy. A level is
C linked to a super level provided that the energy separation is less
C than DEL_E_CM from the lowest level of the SUPER LEVEL. By default
C it is assumed that all levels belonging to a SUPER LEVEL must have
C the same PARITY and SPIN, although this can be changed USING 
C hidden options.
C
C This option is similar to ELS, except not all terms of a multiplet
C will necessarily belong to the same SUPER LEVEL.
C
	ELSE IF(X(1:2) .EQ. 'SP')THEN
	  F_TO_S(:)=0
	  CNT=1
	  F_TO_S(1)=1
	  E_STRT(CNT)=ENERGY(1)
C
	  CALL USR_OPTION(DEL_E_CM,'DEL_E','1000.0D0',
	1               'Maximum energy diff in cm^{-1}')
	  CALL USR_HIDDEN(CHK_PARITY,'CHK_P','T','Check parity?')
	  CALL USR_HIDDEN(CHK_SPIN,'CHK_S','T','Check spin?')
C
	  DO I=2,NLEV
	    J=1
	    DO WHILE(F_TO_S(I) .EQ. 0 .AND. J .LE. I-1)
	      OKAY_TO_COMBINE=.TRUE.
	      IF (LEV_PARITY(I) .NE. LEV_PARITY(J))THEN
	        OKAY_TO_COMBINE=.FALSE.
	      END IF
	      IF (LEV_SPIN(I) .NE. LEV_SPIN(J))THEN
	        OKAY_TO_COMBINE=.FALSE.
	      END IF
	      IF(OKAY_TO_COMBINE .AND. 
	1       ABS(E_STRT(F_TO_S(J))-ENERGY(I)) .LT. DEL_E_CM )THEN
                F_TO_S(I)=F_TO_S(J)
	      ELSE
	        J=J+1
	      END IF
	    END DO
	    IF(F_TO_S(I) .EQ. 0)THEN
	      CNT=CNT+1
	      F_TO_S(I)=CNT
	      E_STRT(CNT)=ENERGY(I)
	    END IF
	  END DO
C
	  WRITE(T_OUT,*)'Number of levels in full atom is    ',NLEV
	  WRITE(T_OUT,*)'Number of LS terms in full atom is  ',N_LS_TERMS
	  WRITE(T_OUT,*)'Number of levels in SUPER atom is   ',CNT
!
C
C The following option allows levels belonging to a single SUPER level
C to be output together as a group. The level departure coefficients can 
C also be C output.
C
	ELSE IF(X(1:5) .EQ. 'SL_WR')THEN
	  CALL USR_OPTION(FILENAME,'File','SL_LNKS','Link check file')
	  CALL USR_OPTION(CNT,'INITIAL','1','Initial SL number')
	  CNT=CNT-1
	  CALL USR_HIDDEN(WRITE_DC,'DC','F',' ')
	  CALL GEN_ASCI_OPEN(LUOUT,FILENAME,'UNKNOWN',' ',
	1                               'WRITE',IZERO,IOS)
C
	  J=0
	  DO I=1,NLEV
	     J=MAX(J,LEN_TRIM(NAME(I)))
	  END DO
	  WRITE(LUOUT,'(A)')'  '
	  L=MAXVAL(F_TO_S(1:NLEV))
	  DO K=1,L
	    DO I=1,NLEV
	      IF(F_TO_S(I) .EQ. K)THEN
	        LAM_EDGE(I)=1.0D+08/(ION_EN-ENERGY(I))
	        IF(WRITE_DC)THEN
	          WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I)+CNT,INT_SEQ(I),I,DC(I,1),
	1                      DC(I,2),DC(I,3)
	        ELSE
	          WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I)+CNT,INT_SEQ(I),I
	        END IF
	      END IF
	    END DO
	    WRITE(LUOUT,*)' '
	  END DO
	  CLOSE(LUOUT)
!
! The following option allows levels belonging to a single seqence
! to be output together as a group. The level departure coefficients can 
! also be output.
!
	ELSE IF(X(1:7) .EQ. 'TERM_WR')THEN
	  CALL USR_OPTION(FILENAME,'File','TERMS','Link check file')
	  CALL USR_HIDDEN(WRITE_DC,'DC','F',' ')
	  CALL GEN_ASCI_OPEN(LUOUT,FILENAME,'UNKNOWN',' ','WRITE',IZERO,IOS)
C
	  J=0
	  DO I=1,NLEV
	    TERM_NAME(I)=NAME(I)
	    J=MAX(J,LEN_TRIM(NAME(I)))
	    K=INDEX(NAME(I),'[')
	    IF(K .NE. 0)TERM_NAME(I)=NAME(I)(1:K-1)
	  END DO
	  WRITE(6,*)'Defined terms'
!
	  WRITE(LUOUT,'(A)')'  '
	  L=MAXVAL(F_TO_S(1:NLEV))
	  DONE_LEV(1:NLEV)=.FALSE.
	  DO K=1,NLEV
	    K1=LEN_TRIM(TERM_NAME(K))
	    IF(.NOT. DONE_LEV(K))WRITE(LUOUT,*)' '
	    DO I=K,NLEV
	      I1=LEN_TRIM(TERM_NAME(I))
	      IF(TERM_NAME(K)(K1-4:K1) .EQ. TERM_NAME(I)(I1-4:I1) .AND. .NOT. DONE_LEV(I))THEN
	        DONE_LEV(I)=.TRUE.
	        LAM_EDGE(I)=1.0D+08/(ION_EN-ENERGY(I))
	        IF(WRITE_DC)THEN
	          WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I),INT_SEQ(I),I,DC(I,1),
	1                      DC(I,2),DC(I,3)
	        ELSE
	          WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I),INT_SEQ(I),I
	        END IF
	      END IF
	    END DO
	  END DO
	  CLOSE(LUOUT)
C                           
	ELSE IF(X(1:6) .EQ. 'SEQ_WR')THEN
	  CALL USR_OPTION(FILENAME,'File','SEQ_LNKS','Link check file')
	  CALL USR_HIDDEN(WRITE_DC,'DC','F',' ')
	  CALL GEN_ASCI_OPEN(LUOUT,FILENAME,'UNKNOWN',' ','WRITE',IZERO,IOS)
C
	  J=0
	  DO I=1,NLEV
	     J=MAX(J,LEN_TRIM(NAME(I)))
	  END DO
	  WRITE(LUOUT,'(A)')'  '
	  L=MAXVAL(F_TO_S(1:NLEV))
	  DONE_LEV(1:NLEV)=.FALSE.
	  DO K=1,NLEV
	    IF(.NOT. DONE_LEV(K))WRITE(LUOUT,*)' '
	    DO I=K,NLEV
	      IF(INT_SEQ(I) .EQ. INT_SEQ(K) .AND. .NOT. DONE_LEV(I))THEN
	        DONE_LEV(I)=.TRUE.
	        LAM_EDGE(I)=1.0D+08/(ION_EN-ENERGY(I))
	        IF(WRITE_DC)THEN
	          WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I),INT_SEQ(I),I,DC(I,1),
	1                      DC(I,2),DC(I,3)
	        ELSE
	          WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I),INT_SEQ(I),I
	        END IF
	      END IF
	    END DO
	  END DO
	  CLOSE(LUOUT)
C                           
C Group all terms belonging to the same LS multiplet.
C               
	ELSE IF(X(1:6) .EQ. 'WR_NOJ')THEN	!LS coupling
!
	  OPEN(UNIT=16,FILE='NOJ_NAMES',STATUS='UNKNOWN')
	    LS_DONE(:)=.FALSE.
	    WRITE(16,*)CNT
	    DO I=1,NLEV
	      J=F_TO_S(I)
	      IF(.NOT. LS_DONE(J))THEN
	        WRITE(16,'(1X,A,T31,A)')TRIM(LS_NAME(I)),TRIM(LS_NAME(I))
	        LS_DONE(J)=.TRUE.
	      END IF
	    END DO
	  CLOSE(UNIT=16)
C
C Output links to file in a format suitable for CMFGEN.
C
	ELSE IF(X .EQ. 'WR')THEN
	  CALL USR_OPTION(FILENAME,'File','F_TO_S_OUT',
	1                     'F_TO_S output file')
	  CALL USR_HIDDEN(WRITE_DC,'DC','F',' ')
	  CALL USR_HIDDEN(HEAD,'HEAD','T',' ')
	  CALL GEN_ASCI_OPEN(LUOUT,FILENAME,'UNKNOWN',' ',
	1                               'WRITE',IZERO,IOS)
	  IF(HEAD)THEN
	    CALL GEN_ASCI_OPEN(LUHEAD,'HEAD_INFO','OLD',' ',
	1                               'READ',IZERO,IOS)
	    DO WHILE(IOS .EQ. 0)  
	      READ(LUHEAD,'(A)',IOSTAT=IOS)STRING
	      IF(INDEX(STRING,'!Date') .EQ. 0)THEN
	         WRITE(LUOUT,'(A)')STRING
	      ELSE
	         IOS=100
	      END IF
	    END DO
	  ELSE
	    WRITE(LUOUT,'(A)')'  '
	  END IF
C
	  CALL DATE_TIME(TIME)
	  WRITE(LUOUT,'(A,T40,A)')TIME(1:11),'!Date'
	  STRING=' '
	  WRITE(STRING,'(I5)')NLEV
	  DO WHILE(STRING(1:1) .EQ. ' '); STRING(1:)=STRING(2:) ; END DO
	  WRITE(LUOUT,'(A,T40,A)')TRIM(STRING),'!Number of energy levels'
	  WRITE(LUOUT,'(I1,T40,A)')6,'!Entry number of link to super level'
	  WRITE(LUOUT,*)' '
	  J=MAX_NAME_LNGTH
	  DO I=1,NLEV
	    LAM_EDGE(I)=1.0D+08/(ION_EN-ENERGY(I))
	    IF(WRITE_DC)THEN
	      WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I),INT_SEQ(I),I,DC(I,1),
	1                      DC(I,2),DC(I,3)
	    ELSE
	      WRITE(LUOUT,100)NAME(I)(1:J),G(I),ENERGY(I),FEDGE(I),
	1                      LAM_EDGE(I),F_TO_S(I),INT_SEQ(I),I
	    END IF
	  END DO
	  CLOSE(LUOUT)
100	  FORMAT(A,5X,F6.1,F16.4,3X,F8.4,1PE12.3,2x,I4,3X,
	1         I3,3X,I4,:,2X,1PE11.4,2X,E11.4,2X,E11.4)
C 
C
C Type links to terminal
C
	ELSE IF(X(1:2) .EQ. 'TY')THEN
	  CALL USR_HIDDEN(WRITE_DC,'DC','F',' ')
	  WRITE(T_OUT,'(I1,T36,A)')6,'!Entry number of link to super level'
	  WRITE(T_OUT,*)' '
	  DO K=1,NLEV,10
	    I=K
	    J=MAX_NAME_LNGTH
	    DO WHILE((I-K) .LE. 14 .AND. I .LE. NLEV)
	      LAM_EDGE(I)=1.0D+08/(ION_EN-ENERGY(I))
	      IF(WRITE_DC)THEN
	        WRITE(T_OUT,120)NAME(I)(1:J),G(I),ENERGY(I),F_TO_S(I),
	1                 I,DC(I,1)
	      ELSE  
	        WRITE(T_OUT,120)NAME(I)(1:J),G(I),ENERGY(I),F_TO_S(I),I
	      END IF
	      I=I+1
	    END DO
	    READ(T_IN,'(A)')ANS
	    IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'Q' .OR.
	1       ANS(1:1) .EQ. 'q' .OR.
	1       ANS(1:1) .EQ. 'e')GOTO 500		!Exit from listing.
	  END DO
500	  CONTINUE
120	  FORMAT(1X,A,5X,F6.1,F20.8,3X,I4,3X,I4,2X,1PE11.4)
C
C Option to re-label the links between super-levels and the full levels
C The only requirement in the raw link is that they are unique, and +ve.
C They do not need to be in order.
C
	ELSE IF(X(1:2) .EQ. 'CL')THEN
C
	  ID=0
	  F_TO_S(1:NLEV)=-F_TO_S(1:NLEV)
	  INT_SEQ(1:NLEV)=-INT_SEQ(1:NLEV)
	  DO I=1,NLEV
	   IF(F_TO_S(I) .LT. 0)THEN
	      ID=ID+1
	      OLD_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	      DO J=I+1,NLEV
	        IF(F_TO_S(J) .EQ. OLD_F_TO_S)F_TO_S(J)=ID
	      END DO
	      DO J=1,NLEV
	        IF(INT_SEQ(J) .EQ. OLD_F_TO_S)INT_SEQ(J)=F_TO_S(I)
	      END DO
	    ELSE IF(F_TO_S(I) .EQ. 0)THEN
	      ID=ID+1
	      OLD_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	    END IF
	  END DO
!
	  DO I=1,NLEV
	    IF(INT_SEQ(I) .NE. 0)THEN
	      WRITE(6,*)INT_SEQ(I)
!	      INT_SEQ(I)=F_TO_S(INT_SEQ(I))
	    END IF
	  END DO 
	  WRITE(T_OUT,*)'Cleaning level links '
	  WRITE(T_OUT,'(A,1X,I5)')' Number of super levels is:',ID
!
! Option simply outputs SUPER states in groups, and indicates the
! maximum energy width of the fgroup.
!
	ELSE IF(X(1:3) .EQ. 'SEP')THEN
	  ID=MAXVAL(F_TO_S(1:NLEV))
	  DO I=1,ID
	    DEL_E_CM=0
	    DO J=1,NLEV
	      IF(F_TO_S(J) .EQ. I)THEN
	        T1=ENERGY(J)
	        IF(DEL_E_CM .EQ. 0)DEL_E_CM=T1
	        WRITE(6,'(1X,A20,5X,F20.8)')NAME(J),ENERGY(J)
	        WRITE(21,'(1X,A20,5X,F20.8)')NAME(J),ENERGY(J)
	      END IF
	    END DO
	    IF(T1 .NE. DEL_E_CM)THEN
	      DEL_E_CM=T1-DEL_E_CM
	      WRITE(6,'(46X,F20.8)')DEL_E_CM
	      WRITE(21,'(46X,F20.8)')DEL_E_CM
	    END IF
	    WRITE(6,'(A)')' '
	    WRITE(21,'(A)')' '
	  END DO
!
!
	ELSE IF(X(1:6) .EQ. 'RD_LNK')THEN
C
C Read in a previously existing link file. If this file has been edited,
C and the links numbering is all mixed up, the 'CL' option should be 
C issued.
C
	  FILE_PRES=.FALSE.
	  DO WHILE(.NOT. FILE_PRES)
	    CALL USR_OPTION(FILENAME,'File',' ','F_TO_S links')
	    INQUIRE(FILE=TRIM(FILENAME),EXIST=FILE_PRES)
	    IF(.NOT. FILE_PRES)WRITE(T_OUT,*)'File ',TRIM(FILENAME),' not found. Try again.'
	  END DO
	  CALL USR_OPTION(SORT_LEVS,'SORT','F','Sort levels if F_TO_FILE?')
	  CALL RD_F_TO_S_IDS_V2(F_TO_S,INT_SEQ,NAME,NLEV,IZERO,
	1               LUIN,SORT_LEVS,FILENAME)
	   DO I=1,NLEV
	     IF(INT_SEQ(I) .NE. 0)WRITE(6,*)I,INT_SEQ(I)
	    END DO	
C
	ELSE IF(X(1:9) .EQ. 'RD_SM_LNK')THEN
C
C This option allows links for with nonsplit terms to be read in.
C These links can then be used for a model atom with split terms.
C Split terms are grouped into the same SUPER level. The number
C of terms in both files must be identical.
C
	  CALL USR_OPTION(FILENAME,'File',' ','F_TO_S links')
	  J=0
	  DO I=1,NLEV
	    K=INDEX(NAME(I),'[')-1
	    IF(K .LE. 0)K=LEN_TRIM(NAME(I))
	    PRES=.FALSE.
	    DO L=1,J
	      IF(NAME(I)(1:K) .EQ. TERM_NAME(J))PRES=.TRUE.
	    END DO
	    IF(.NOT. PRES)THEN
	      J=J+1
	      TERM_NAME(J)=NAME(I)(1:K)
	    END IF
	  END DO
	  NTERM=J
C
	  CALL RD_F_TO_S_IDS(TERM_F_TO_S,INT_SEQ,TERM_NAME,
	1         NTERM,IZERO,LUIN,FILENAME)
C
	  DO I=1,NLEV
	    K=INDEX(NAME(I),'[')-1
	    IF(K .LE. 0)K=LEN_TRIM(NAME(I))
	    DO J=1,NTERM
	      IF(NAME(I)(1:K) .EQ. TERM_NAME(J))F_TO_S(I)=TERM_F_TO_S(J)
	    END DO
	  END DO
C
C Read in departure coefficients for a model computed with this atomic 
C model. If model was run with N_S=N_F the DC can be output with the
C levels names in order to assist in deciding super level assignments.
C The departure coefficent at 3 distinct depths is read in.
C
	ELSE IF(X(1:5) .EQ. 'RD_DC')THEN
C
	  CALL USR_OPTION(FILENAME,'File',' ',
	1                   'Departure coefficients ')
C
	  INQUIRE(LUIN,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(LUIN)
	  CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ',
	1             'READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error --- unable to open D.C. FILE'
	    GOTO 1
	  END IF
C
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file. 
C
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LUIN,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LUIN)
C
	  READ(LUIN,*)T1,T1,NLEV_RD,ND_RD
	  DC(1:NLEV,1:3)=0
	  READ(LUIN,*)T1,T1,ED(1),TEMP(1)
	  READ(LUIN,*)(DC(I,1),I=1,NLEV_RD)
	  DO J=2,ND_RD/3+1
	    READ(LUIN,*)T1,T1,ED(2),TEMP(2)
	    READ(LUIN,*)(DC(I,2),I=1,NLEV_RD)
	  END DO
	  DO J=ND_RD/3+2,2*ND_RD/3+1
	    READ(LUIN,*)T1,T1,ED(3),TEMP(3)
	    READ(LUIN,*)(DC(I,3),I=1,NLEV_RD)
	  END DO
	  CLOSE(LUIN)
!
! 
! Routine to create a Departure Coefficient file for INPUT to CMFGEN.
! The routine reads a file containing EXCTATION temperatures (created using
! DISPGEN) for the next lowest ionization stage. The excitation temperature 
! of the ground state is then used to comute the ION population, and the 
! departure coefficients for all other levels. Routine is useful when including
! an additional ionization stage.
! 
	ELSE IF(X(1:5) .EQ. 'WR_DC')THEN
!
	  CALL USR_OPTION(FILENAME,'File',' ','File containing excitation temperatures')
	  INQUIRE(LUIN,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(LUIN)
	  CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error --- unable to open D.C. FILE'
	    GOTO 1
	  END IF
!
	  CALL USR_OPTION(FILENAME,'File','DC_OUT','Output file')
	  CALL GEN_ASCI_OPEN(LUOUT,FILENAME,'UNKNOWN',' ',
	1                               'WRITE',IZERO,IOS)
	  CALL GEN_ASCI_OPEN(LUOUT+1,'TX_OUT','UNKNOWN',' ',
	1                               'WRITE',IZERO,IOS)
C
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file. 
C
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LUIN,'(A)')STRING
	    WRITE(LUOUT,'(A)')TRIM(STRING)
	    WRITE(LUOUT+1,'(A)')TRIM(STRING)
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LUIN)
	  WRITE(LUOUT,'(A)')' '
	  WRITE(LUOUT+1,'(A)')' '
C
	  READ(LUIN,*)RSTAR,RLUM,NLEV_RD,ND_RD
	  WRITE(LUOUT,'(3X,F9.4,3X,ES12.4,2(4X,I3))')RSTAR,RLUM,NLEV,ND_RD
	  WRITE(LUOUT+1,'(3X,F9.4,3X,ES12.4,2(4X,I3))')RSTAR,RLUM,NLEV,ND_RD
!
! Because of level splitting, G_GS may be larger than G(1). That is, G_GS usually
! refers to the statistical weight of the full term.
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Set FIX_DI=T if input and output species are identical'
	  WRITE(T_OUT,*)' '
	  CALL USR_OPTION(FIX_DI,'FIX_DI','F','Keep ion density fixed')
	  IF(.NOT. FIX_DI)THEN
	    CALL USR_OPTION(G_GS,'G_GS','1.0D0','Input GION for T excitation file')
	    CALL USR_OPTION(G_ION,'G_ION','1.0D0','Input GION for next ionzation stage')
	  END IF
	  DO J=1,ND_RD
	    READ(LUIN,*)RVAL,DI,ED(1),TEMP(1),R_IR,V1,CL_FAC
	    READ(LUIN,*)(DC(I,1),I=1,NLEV_RD)
!
! Compute ion population.
!
	    T_EXCITE=DC(1,1)
	    IF(.NOT. FIX_DI)THEN
	      T1=2.07078D-22*ED(1)*(G_GS/G_ION)*EXP(HDKT*FEDGE(1)/T_EXCITE)/(T_EXCITE**1.50)
	      DI=DI/T1
	    END IF
!
	    WRITE(LUOUT,'(A)')' '
	    WRITE(LUOUT,'(ES16.7,6ES15.4)')RVAL,DI,ED(1),TEMP(1),0.0D0,V1,CL_FAC
	    WRITE(LUOUT+1,'(A)')' '
	    WRITE(LUOUT+1,'(ES16.7,6ES15.4)')RVAL,DI,ED(1),TEMP(1),0.0D0,V1,CL_FAC
	    DO I=1,NLEV
	      T1=EXP( HDKT*FEDGE(I)*(1.0D0/T_EXCITE-1.0D0/TEMP(1)) )
	      T2=(TEMP(1)/T_EXCITE)**1.5D0
	      DC(I,1)=T1*T2
	    END DO
	    WRITE(LUOUT,'(3X,5ES14.4)')(DC(I,1),I=1,NLEV)
	    WRITE(LUOUT+1,'(3X,5ES14.4)')(T_EXCITE,I=1,NLEV)
	  END DO
	  CLOSE(LUIN)
	  CLOSE(LUOUT)
C
C 
C
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:2) .EQ. 'HE')THEN
	  IF(X(1:2) .EQ. 'LI')THEN
	    OPEN(UNIT=30,FILE='WR_F_OPT_DESC',STATUS='OLD',
	1     ACTION='READ',IOSTAT=IOS)
	  ELSE
	    OPEN(UNIT=30,FILE='WR_F_OPTIONS',STATUS='OLD',
	1     ACTION='READ',IOSTAT=IOS)
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening description file'
	    GOTO 1
	  END IF
	  READ(30,*)I,K			!For page formating (I=22,K=12)
	  DO WHILE(1.EQ. 1)
	    DO J=1,I
	      READ(30,'(A)',END=700)STRING
	      L=ICHRLEN(STRING)
	      IF(L .EQ. 0)THEN
	        WRITE(T_OUT,'(1X)')
	      ELSE
	        WRITE(T_OUT,'(1X,A)')STRING(1:L)
	      END IF
	    END DO
	    READ(T_IN,'(A)')STRING
	    IF(STRING(1:1) .EQ. 'E' .OR.
	1       STRING(1:1) .EQ. 'e')GOTO 700		!Exit from listing.
	    I=K
	  END DO
700	  CONTINUE
	  CLOSE(UNIT=30)
C
C
	ELSE IF(X(1:2) .EQ.'EX') THEN
	  STOP
	ELSE IF(X(1:3) .EQ. 'BOX') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
C
1	CONTINUE
	GO TO 3
C
	END

	FUNCTION PARITY(NAME)
	IMPLICIT NONE
C
	CHARACTER*(*) NAME
	CHARACTER*1 PARITY
C
	EXTERNAL ICHRLEN
	INTEGER ICHRLEN,J,T_OUT
C
	T_OUT=6
	J=INDEX(NAME,'[')-1
	IF(J .LE. 0)J=ICHRLEN(NAME)
	PARITY=NAME(J:J)
	IF(PARITY .NE. 'e' .AND. PARITY .NE. 'o')THEN
	   IF(PARITY .NE. 'Z' .AND. PARITY .NE. '_')THEN
	     WRITE(T_OUT,*)'Warning: Indeterminate parity for level',NAME
	   END IF
	   PARITY=' '
	END IF
	RETURN
	END
!
!
	FUNCTION SPIN(NAME)
	IMPLICIT NONE
!
! Altered 31-Jan-2003: T_OUT set to 6 (instead of 5)
! Altered 20-Oct-2000: Spin determination altered to allow for names
!                        like .. 3Pbe or 3P2e.
!
	CHARACTER*(*) NAME
	CHARACTER*1 SPIN
C
	EXTERNAL ICHRLEN
	INTEGER ICHRLEN,LS,J,IOS,T_OUT
	T_OUT=6
C
	J=INDEX(NAME,'[')-1
	IF(J .LE. 0)J=ICHRLEN(NAME)
	IF(NAME(J:J) .EQ. 'e' .OR. NAME(J:J) .eq. 'o')THEN
	  LS=J-2
	ELSE
	  LS=J-1
	END IF
!
! Allow for names of the form ... 3Pbe or 3P2e where b and 2 are used to
! break name degeneracies.
!
	IF(NAME(LS:LS) .GT. 'A' .AND. NAME(LS:LS) .LE. 'Z')LS=LS-1
!
	SPIN=NAME(LS:LS)
	IOS=0
	READ(SPIN,'(I1)',IOSTAT=IOS)J
	IF(IOS .NE. 0)THEN
	   WRITE(T_OUT,*)'Indeterminate spin for level',NAME
	   SPIN=' '
	END IF
	RETURN
	END
