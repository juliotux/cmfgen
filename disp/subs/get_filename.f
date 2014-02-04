C
C General routine to obtain filname from user.
C
C If interactive process, routine will ask for file name again
C if an error occurs.
C
C If EXIST=0, file does not have to exist.
C If EXIST=1, file has to exist if a filename is input.
C If EXIST=2, file has to exist, and a blank filename is not allowed.
C
C If user inpute a directory name, and wishes to change directory
C to current directory the user need only input : ot ]. The directory
C default will then be blank.
C
	SUBROUTINE GET_FILENAME(FILNAME,LUIN,LUOUT,DESC2,MARK,EXIST)
	IMPLICIT NONE
C
C 10-Nov-1990 --- Bug fix, cleaned.
C 11-Oct-1990 --- Extensively altered (Name changed to GET_FILENAME)
C 28-Nov-1989 --- Created .
C
	INTEGER LUIN,LUOUT,EXIST
	CHARACTER*(*) FILNAME,DESC2,MARK
	CHARACTER*80 DIRNAME,DESC
	DATA DIRNAME/' '/		!Initialize first time
	SAVE DIRNAME			!Need on subsequent calls
C
	EXTERNAL ICHRLEN,TERMINAL_LU
	INTEGER ICHRLEN,LFIL,LDIR,LENFIL,LENDIR,LDESC
	LOGICAL TERMINAL_LU,INTERACTIVE,FILE_PRES
C
C If interactive, we promft for file name again if an error has occurred.
C
	INTERACTIVE=TERMINAL_LU(LUIN)
	DESC=DESC2
	LDESC=ICHRLEN(DESC)
C
100	CONTINUE
	LENFIL=LEN(FILNAME)
	LENDIR=LEN(DIRNAME)
	LDIR=ICHRLEN(DIRNAME)
C
C Output request for filename. Default directory is indicated.
C
	IF(LUOUT .NE. 0)THEN
	  IF(LDIR .EQ. 0)THEN
	    WRITE(LUOUT,*)'Input ',DESC(1:LDESC),' filname'
	  ELSE
	    WRITE(LUOUT,*)'Input ',DESC(1:LDESC),' filname {'
	1                 ,DIRNAME(1:LDIR),'}'
	  END IF
	END IF
C
C Input filname. If no DIRNAME prefix, DIRNAME prefix is added.
C Otherwise DIRNAME prefix is stored in DIRNAME. If no filename
C is input, FILNAME is returned blank. Filename can have a comment
C attached to it.
C
	READ(LUIN,'(A)')FILNAME
	LFIL=INDEX(FILNAME,MARK)
	IF(LFIL .EQ. 0)THEN
	  IF(FILNAME(LENFIL:LENFIL) .NE. ' ')THEN
	    WRITE(LUOUT,999)DESC(1:LDESC),': Filename of insufficient length'
	    IF(INTERACTIVE)GOTO 100
	    STOP
	  END IF
	  LFIL=ICHRLEN(FILNAME)
	ELSE
	  LFIL=ICHRLEN(FILNAME(1:LFIL-1))
	  FILNAME=FILNAME(1:LFIL)		!Removes unwanted jnk.
	END IF
C
C LFIL=0 is taken to mean a valid file name. Can be used to skip
C over a particular file opening. (As required by DISPGEN, for example).
C
	IF(LFIL .EQ. 0)THEN
	  IF(EXIST .EQ. 2)THEN
	    WRITE(LUOUT,999)DESC(1:LDESC),': No filename input'
	    IF(INTERACTIVE)GOTO 100
	    STOP
	  END IF
	  FILNAME=' '
	  RETURN
	END IF
C
	LDIR=MAX( INDEX(FILNAME,']'),INDEX(FILNAME,':') )
C
C Take present working directory as default. Signaled by : or ] as first
C charcter of filename, if dirname already set.
C
	IF(LDIR .EQ. 1)THEN
	  FILNAME(1:)=FILNAME(2:)
	  DIRNAME=' '
	  LDIR=0		!i.e. no directrory (using default).
	END IF
C
C Attach or extract directory name. Check whether file exists.
C
	IF(LDIR .NE. 0)THEN
C
C Check whether file exists, or not.
C
	  IF(EXIST .NE. 0)THEN
	    INQUIRE(FILE=FILNAME,EXIST=FILE_PRES)
	    IF(.NOT. FILE_PRES)THEN
	      WRITE(LUOUT,999)DESC(1:LDESC),': File does not exist'
	      IF(INTERACTIVE)GOTO 100
	      STOP
	    END IF
	  END IF
C
	  IF(LENDIR .GE. LDIR)THEN
	    DIRNAME=FILNAME(1:LDIR)
	    RETURN
	  ELSE
	    WRITE(LUOUT,999)DESC(1:LDESC),': Dirname of insufficient length'
	    WRITE(LUOUT,*)LENDIR,LDIR
	    IF(INTERACTIVE)GOTO 100
	    STOP
	  END IF
	ELSE
	  LDIR=ICHRLEN(DIRNAME)
	  IF(LENFIL .LT. LDIR+LFIL)THEN
	    WRITE(LUOUT,999)DESC(1:LDESC),': Filename of insufficient length'
	    WRITE(LUOUT,*)LENFIL,LDIR+LFIL
	    IF(INTERACTIVE)GOTO 100
	    IF(.NOT. INTERACTIVE)STOP
	  END IF
	  IF(LDIR .NE. 0)FILNAME=DIRNAME(1:LDIR)//FILNAME(1:LFIL)
C
C Check whether file exists, or not.
C
	  IF(EXIST .NE. 0)THEN
	    INQUIRE(FILE=FILNAME,EXIST=FILE_PRES)
	    IF(FILE_PRES)RETURN
	    WRITE(LUOUT,999)DESC(1:LDESC),': File does not exist'
	    IF(INTERACTIVE)GOTO 100
	    STOP
	  END IF
	END IF
C
999	FORMAT(' Error in GET_FILENAME --- ',A,1X,A)
	END
