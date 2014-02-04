C
C These two functions allow the logical units associated with general 
C input/output file (terminal) to be specified. Virtually all programs 
C currently use units 1 and 2. May need to go to units 5 and 6 for Unix.
C The LU for output and error messages are the same.
C
	FUNCTION TERM_IN()
	INTEGER TERM_IN
	TERM_IN=5
	RETURN
	END
C
	FUNCTION TERM_OUT()
	INTEGER TERM_OUT
	TERM_OUT=6
	RETURN
	END
C
	FUNCTION ERROR_LU()
	INTEGER ERROR_LU
	ERROR_LU=6
	RETURN
	END
