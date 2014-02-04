!
! Subroutine to output Gaussian fit. Designed to be called
! separately from fitting routine.
!
	SUBROUTINE WR_GAUS_FIT
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created 9-Mar-2008
!
	INTEGER, PARAMETER :: NL_MAX=100
	REAL*8 LINE_CENTER(NL_MAX)
!
	REAL*8 T1
	REAL*8 FWHM
	INTEGER I,J,K
!
	LOGICAL NEW_FILE
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, SAVE :: FIRST_WRITE=.TRUE.
!
	IF(FIRST_WRITE)THEN
	    FIRST_WRITE=.FALSE.
	    INQUIRE(FILE='GAUSS_FITS',EXIST=NEW_FILE); NEW_FILE=.NOT. NEW_FILE
	    IF(.NOT. NEW_FILE)THEN
	      WRITE(6,*)'File GAUS_FITS exists: Default is to append new fits'
	      CALL GEN_IN(NEW_FILE,'Open new file (T) or append data (F)')
	    END IF
	    IF(NEW_FILE)THEN
	      OPEN(UNIT=35,FILE='GAUSS_FITS',STATUS='UNKNOWN')
	      WRITE(35,'(2X,(5X,A),5(5X,A),2(3X,A))')
	1          '        Lam','     Height','          a',
	1          '        EXP','Sigma(km/s)',' FWHM(km/s)',' EW(mA)','Err(mA)'
	    ELSE
	      OPEN(UNIT=35,FILE='GAUSS_FITS',STATUS='UNKNOWN',POSITION='APPEND')
	    END IF
	ELSE
	    OPEN(UNIT=35,FILE='GAUSS_FITS',STATUS='UNKNOWN',POSITION='APPEND')
	END IF
!
! We output lines in wavelength order. After calling, INDEXX,
! INDX_VEC(1) will contain the first line to be output, etc.
! We don't actually perform the sort.
!
        IF(ALLOCATED(INDX_VEC))DEALLOCATE(INDX_VEC)
        ALLOCATE (INDX_VEC(NUM_GAUS))
	DO I=1,NUM_GAUS
	  K=3+(I-1)*4
	  LINE_CENTER(I)=PAR(K)
	END DO
	IF(NUM_GAUS .EQ. 1)THEN
	  INDX_VEC(1)=1
	ELSE
	  CALL INDEXX(NUM_GAUS,LINE_CENTER,INDX_VEC,L_TRUE)
	END IF
!
	WRITE(35,'(A)')' '
	WRITE(35,'(I3,8X,A)')NUM_GAUS,'!Number of Gaussians'
	WRITE(35,'(2X,2ES16.6)')X_GAUS(1),X_GAUS(NG_DATA)
	WRITE(35,'(2X,2ES16.6)')PAR(1:2)
	DO J=1,NUM_GAUS
	   I=INDX_VEC(J)
	   K=2+(I-1)*4+1
	   T1=2.99794D+05*SIM(1,K+1)/SIM(1,K)
	   FWHM=2.0D0*T1*(DLOG(2.0D0))**(1.0D0/SIM(I,K+3))
	   WRITE(35,'(2X,ES16.6,5ES16.4,3F10.2)')SIM(I,K),SIM(1,K+2),SIM(1,K+1),SIM(I,K+3),
	1               T1/SQRT(2.0D0),FWHM,EW(I),EW_ERROR(I),ALT_ERROR(I)
	END DO
	CLOSE(UNIT=35)
!	
	RETURN
	END
