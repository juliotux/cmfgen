!
! Subroutine to read in the autioization rates for levels above the
! ionization limit. These rates are added to Gamma (= Sum Einstein A values).
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
! NB - ZION is the charge on the ion - thus ZHYD=1.0D0
!
! Routine does not work for NUM_BNDS=ND.
!
	SUBROUTINE ADD_AUTO_RATES(ARAD,FEDGE_F,G_F,LEVNAME_F,N_F,AUTO_FILE)
	IMPLICIT NONE
!
! Created: 29-Aug-2015
!
	INTEGER N_F
	REAL*8 ARAD(N_F)
	REAL*8 FEDGE_F(N_F)
	REAL*8 G_F(N_F)
	CHARACTER(LEN=*) LEVNAME_F(N_F)
	CHARACTER(LEN=*) AUTO_FILE
!
! Local variables.
!
	INTEGER I
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL*8 AUTO(N_F)
	LOGICAL AUTO_FILE_EXISTS
!
! The autoionization probabilities are depth independent, and hence can be
! read in at the beginning.
!
        INQUIRE(FILE=AUTO_FILE,EXIST=AUTO_FILE_EXISTS)
        IF(AUTO_FILE_EXISTS)THEN
          CALL RD_AUTO_V1(AUTO,FEDGE_F,G_F,LEVNAME_F,N_F,AUTO_FILE)
	  DO I=1,N_F
	    ARAD(I)=ARAD(I)+AUTO(I)
	  END DO
        ELSE IF(FEDGE_F(N_F) .LT. 0.0D0)THEN
          LUER=ERROR_LU()
          WRITE(LUER,*)' '
          WRITE(LUER,*)'Warning: possible error in ADD_AUTO_RATES'
          WRITE(LUER,*)'No autoionization probabilities available'
          WRITE(LUER,*)'Model has states above the ionization limit'
          WRITE(LUER,'(A,A)')' AUTO_FILE is ',TRIM(AUTO_FILE)
          WRITE(LUER,*)' '
          RETURN
        END IF
!
	RETURN
	END
