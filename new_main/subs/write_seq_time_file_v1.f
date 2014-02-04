!
! Subroutine to output model data to a sequential, unformatted, file.
! This data will be subsequently used for the next model in the time sequence.
!
	SUBROUTINE WRITE_SEQ_TIME_FILE_V1(SN_AGE,ND,LU)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 22-Jul-2008 :  Changed to facilitate addition/deletion of ionization stages
!                           from time dependent models. NUM_SPECIES and ZXzV for each
!                           ion is now output.
! Created 14-Mar-2007 :  Based on READ_TIME_MODEL_V2
!
	REAL*8 SN_AGE
	INTEGER ND
	INTEGER LU
!
	INTEGER ISPEC
	INTEGER ID
!
	OPEN(UNIT=LU,FILE='CUR_MODEL_DATA',FORM='UNFORMATTED',STATUS='UNKNOWN',ACTION='WRITE')
!
	  WRITE(LU)'21-Jul-2008'		!Format date
	  WRITE(LU)ND,NUM_SPECIES
	  WRITE(LU)SN_AGE
	  WRITE(LU)R
	  WRITE(LU)V
	  WRITE(LU)SIGMA
	  WRITE(LU)T
	  WRITE(LU)ED
	  WRITE(LU)POP_ATOM
	  WRITE(LU)DENSITY
!
! We use ZXzV as a means of identifying the ionization stage.
!
	  DO ISPEC=1,NUM_SPECIES
	    WRITE(LU)SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1,SPECIES(ISPEC)
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      WRITE(LU)ID,ATM(ID)%NXzV_F,ATM(ID)%ZXzV
	      WRITE(LU)ATM(ID)%XzV_F
	      WRITE(LU)ATM(ID)%DXzV		!Super level ion population. 
	    END DO
	  END DO
!
	CLOSE(LU)
	RETURN
	END
