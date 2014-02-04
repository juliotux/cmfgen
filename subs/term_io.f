C
C These two functions allow the logical units associated with general 
C input/output file (terminal) to be specified. Virtually all programs 
C in the past used units 1 and 2. Needed to go to units 5 and 6 for Unix.
C The LU for output and error messages are the same, and was changed to
C unit 6 on Jul7-200 so that error messages when running DISPGEN, PLT_SPEC,
C etc would automatically go to the terminal.
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
