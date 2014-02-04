!
! Subroutine to fit a straigt line, and a set of modified Gaussians,
! to a set of data. Routine is designed to be called in GRAMON_PGPLOT.
! Routine is designed to fit multiple sections/multiples plots, and
! assumes a fit to a single plot has been done previously.
!
! Based on DO_GAUS_FIT : Must retain compatibility in variable usage with
!                        DO_GAUS_FIT.
!
	SUBROUTINE DO_MULT_GF
	USE MOD_CURVE_DATA
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created 31-Dec-2007
!
	INTEGER, PARAMETER :: NL_MAX=100
	REAL*8 LINE_HEIGHT(NL_MAX)
	REAL*8 LINE_CENTER(NL_MAX)
	REAL*8 LINE_SIGMA(NL_MAX)
!
	REAL*8 TOLERANCE
	REAL*8 GAUS_FIT_FUNC
	EXTERNAL GAUS_FIT_FUNC
!
	REAL*8 XST,XEND
	REAL*8 YST,YEND
!
	REAL*8 TOL
	REAL*8 T1,T2
	REAL*8 FWHM
	INTEGER ITER
!
	INTEGER IP
        LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	INTEGER LU_IN,LU_OUT,LU_ALT
	INTEGER I,J,K
	LOGICAL GUESSED
	LOGICAL NEW_REGION
	LOGICAL VERBOSE
	LOGICAL FILE_EXISTS
	LOGICAL FILE_OPEN
!
	CHARACTER(LEN=80) OLD_GF_FILE
	CHARACTER(LEN=80) OUTPUT_FILE
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=1)  ACT
!
	VERBOSE=.TRUE.
	OLD_GF_FILE='OLD_GAUSS_FITS'
	OUTPUT_FILE='NEW_GAUSS_FITS'
!
! If we are fitting the same regon as previously, we can use the previous findings.
!
	CALL GET_LU(LU_IN,'LU_IN in DO_MULT_GF')
	CALL GEN_IN(OLD_GF_FILE,'File with old Gaussians fit for similar model')
	IF(OLD_GF_FILE .NE. 'NONE')THEN
	  OPEN(UNIT=LU_IN,FILE=OLD_GF_FILE,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error: Unable to open old Gaus fit file'
	    WRITE(6,*)'Error: Unable to open ',TRIM(OLD_GF_FILE)
	    INQUIRE(UNIT=LU_IN,OPENED=FILE_OPEN)
            IF(FILE_OPEN)CLOSE(LU_IN)
	    RETURN
	  END IF
	END IF
!
! Open output file, and write out header. We use FLUSH option at end of
! the Gauss fitting to force output after each fitting.
!
	CALL GET_LU(LU_OUT,'LU_OUT in DO_MULT_GF')
20	CALL GEN_IN(OUTPUT_FILE,'Output file for new Gaussian fits')
	INQUIRE(FILE=TRIM(OUTPUT_FILE),EXIST=FILE_EXISTS)
        IF(FILE_EXISTS)THEN
	  WRITE(6,*)TRIM(OUTPUT_FILE),' exists'
	  CALL GEN_IN(ACT,'Append, Overwrite, of Get new filename (A, O or G)')
	  IF(ACT .EQ.'a' .OR. ACT .EQ. 'A')THEN
	    OPEN(UNIT=LU_OUT,FILE=OUTPUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND',IOSTAT=IOS)
	  ELSE IF(ACT .EQ.'o' .OR. ACT .EQ. 'O')THEN
	    OPEN(UNIT=LU_OUT,FILE=OUTPUT_FILE,STATUS='OLD',ACTION='WRITE',IOSTAT=IOS)
	  ELSE 
	    GOTO 20
	  END IF
	ELSE
	  OPEN(UNIT=LU_OUT,FILE=OUTPUT_FILE,STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
	END IF
!
	CALL GET_LU(LU_ALT,'LU_ALT in DO_MULT_GF')
	OUTPUT_FILE='ALT_'//TRIM(OUTPUT_FILE)
30	CALL GEN_IN(OUTPUT_FILE,'Output file with alternate file format')
	INQUIRE(FILE=TRIM(OUTPUT_FILE),EXIST=FILE_EXISTS)
        IF(FILE_EXISTS)THEN
	  WRITE(6,*)TRIM(OUTPUT_FILE),' exists'
	  CALL GEN_IN(ACT,'Append, Overwrite, of Get new filename (A, O or G)')
	  IF(ACT .EQ.'a' .OR. ACT .EQ. 'A')THEN
	    OPEN(UNIT=LU_ALT,FILE=OUTPUT_FILE,STATUS='OLD',ACTION='WRITE',POSITION='APPEND',IOSTAT=IOS)
	  ELSE IF(ACT .EQ.'o' .OR. ACT .EQ. 'O')THEN
	    OPEN(UNIT=LU_ALT,FILE=OUTPUT_FILE,STATUS='OLD',ACTION='WRITE',IOSTAT=IOS)
	  ELSE 
	    GOTO 30
	  END IF
	ELSE
	  OPEN(UNIT=LU_ALT,FILE=OUTPUT_FILE,STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
	END IF
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error oppening ',TRIM(OUTPUT_FILE)
	  WRITE(6,*)'IOSTAT=',IOS
	  INQUIRE(UNIT=LU_IN,OPENED=FILE_OPEN)
          IF(FILE_OPEN)CLOSE(LU_IN)
	  INQUIRE(UNIT=LU_OUT,OPENED=FILE_OPEN)
          IF(FILE_OPEN)CLOSE(LU_OUT)
	  RETURN
	END IF
	WRITE(LU_OUT,'(2X,(5X,A),5(5X,A),2(3X,A))')
	1          '        Lam','     Height','          a',
	1          '        EXP','Sigma(km/s)',' FWHM(km/s)',' EW(mA)','Err(mA)'
!
	DO WHILE(1 .EQ. 1)	
	  STRING=' '
	  IF(OLD_GF_FILE .EQ. 'FINISHED')THEN
	    RETURN
	  ELSE IF(OLD_GF_FILE .NE. 'NONE')THEN
	    DO WHILE(INDEX(STRING,'!Number of Gaussians') .EQ. 0)
	      READ(LU_IN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        CLOSE(UNIT=LU_IN)
	        CLOSE(UNIT=LU_OUT)
	        RETURN
	      END IF
	    END DO
	    READ(STRING,*)NUM_GAUS
	    NG_PAR=2+4*NUM_GAUS
            IF(ALLOCATED(PAR))DEALLOCATE(PAR)
	    ALLOCATE (PAR(NG_PAR))
! 
	    READ(LU_IN,*)XST,XEND
	    READ(LU_IN,*)PAR(1),PAR(2)
	    DO I=1,NUM_GAUS
	      K=3+4*(I-1)
              READ(LU_IN,*)PAR(K),PAR(K+2),PAR(K+1),PAR(K+3)
	    END DO
	  ELSE
	    XST=X_GAUS(1); XEND=X_GAUS(NG_DATA)
	    OLD_GF_FILE='FINISHED'
	  END IF
!
	  IF(ALLOCATED(SIM))DEALLOCATE(SIM,SUM_SQ,SCALE,EW,EW_TABLE,LAM_TABLE)
	  ALLOCATE (SIM(NG_PAR+1,NG_PAR))
	  ALLOCATE (SUM_SQ(NG_PAR+1))
	  ALLOCATE (SCALE(NG_PAR+1))
	  ALLOCATE (EW(NUM_GAUS))
	  ALLOCATE (EW_TABLE(NUM_GAUS,NPLTS))
	  ALLOCATE (LAM_TABLE(NUM_GAUS,NPLTS))
!
	  DO IP=1,NPLTS
	    WRITE(6,*)CD(IP)%XVEC(1),CD(IP)%XVEC(NPTS(IP))
	    WRITE(6,*)CD(IP)%DATA(1),CD(IP)%DATA(NPTS(IP))
	    CALL SET_GAUS_DATA(CD(IP)%XVEC,CD(IP)%DATA,XST,XEND,NPTS(IP),YST,YEND)
!
	    SIM(1,1:NG_PAR)=PAR(1:NG_PAR)
!
! Set other vertices of simplex. To do this we first set the
! the characteristic scaleover which the variable may very.
!
	    SCALE(1)=0.002		!Mean level
	    SCALE(2)=0.02/(XEND-XST)	!Slope
	    DO I=3,NG_PAR,4
	      SCALE(I)=0.2		!Wavelength
	      SCALE(I+1)=0.05		!Sigma
	      SCALE(I+2)=0.05		!Height
	      SCALE(I+3)=0.2		!Exponent
	    END DO
!
	    DO J=2,NG_PAR+1
	      DO I=1,NG_PAR
	        SIM(J,I)=SIM(1,I)
	      END DO
	    END DO
	    DO J=2,NG_PAR+1
	      SIM(J,J-1)=SIM(J,J-1)+SCALE(J-1)		!ws scale(j)
	    END DO
!
	    IF(VERBOSE)THEN
	      WRITE(6,*)'Simples vertices have been set'
	      DO J=1,NG_PAR
	        WRITE(6,'(20ES12.4)')(SIM(I,J),I=1,NG_PAR+1)
	      END DO
	    END IF
! 
	    SUM_SQ(:)=0.0D0
            DO I=1,NG_PAR+1
              PAR(1:NG_PAR)=SIM(I,1:NG_PAR)
              SUM_SQ(I)=GAUS_FIT_FUNC(PAR)
            END DO
	    IF(VERBOSE)THEN
	      WRITE(6,*)'Evaluated SUM_SQ: Will now do the fit'
	    END IF
!
            TOL=1.0D-08
	    I=NG_PAR+1
            CALL AMOEBA(SIM,SUM_SQ,I,NG_PAR,NG_PAR,TOL,GAUS_FIT_FUNC,ITER)
!
! Compute EW of profile. The order of passing the variables is HEIGHT,
! SIGMA, EXPONENT. The location is not needed.
!
	    TOLERANCE=1.0D-05
	    DO I=1,NUM_GAUS
	      K=2+(I-1)*4+1
	      CALL GAUS_ROMB(EW(I),SIM(1,K+2),SIM(1,K+1),SIM(1,K+3),TOLERANCE)
	      T1=SIM(1,1)+SIM(1,2)*(SIM(I,K)-X_GAUS(1))				!Continuum value
	      EW(I)=EW(I)/T1
	    END DO
	    EW=EW*1000.0D0			!mAng
            PAR(1:NG_PAR)=SIM(1,1:NG_PAR)
	    CALL GAUS_FIT_ER(PAR);	EW_ERROR=EW_ERROR*1000.0D0
!
	    IF(VERBOSE)THEN
	      WRITE(6,*)'Called AMOEBA'
	      WRITE(6,*)'Fit parameters are (NB SIGMA is not Stan. Dev.):'
	      WRITE(6,*)SIM(1,1),SIM(1,2)
	      WRITE(6,'(2X,(5X,A),5(3X,A),2(3X,A))')
	1          '        Lam','     Height','          a',
	1          '        EXP','Sigma(km/s)',' FWHM(km/s)',' EW(mA)','Err(mA)'
	      DO I=1,NUM_GAUS
	        K=2+(I-1)*4+1
	        T1=2.99794D+05*SIM(1,K+1)/SIM(1,K)
	        T2=2.0D0; FWHM=2.0D0*T1*(LOG(T2))**(1.0D0/SIM(I,K+3))
	        WRITE(6,'(2X,ES16.6,5ES14.4,2F10.2)')SIM(I,K),SIM(1,K+2),SIM(1,K+1),SIM(I,K+3),
	1               T1/SQRT(T2),FWHM,EW(I),EW_ERROR(I)
	      END DO
	    END IF
!
            SUM_SQ(:)=0.0D0
            SUM_SQ(1)=GAUS_FIT_FUNC(PAR)
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
	    WRITE(LU_OUT,'(A)')' '
	    WRITE(LU_OUT,'(I3,8X,A)')NUM_GAUS,'!Number of Gaussians'
	    WRITE(LU_OUT,'(2X,2ES16.6,5X,I3)')XST,XEND,IP
	    WRITE(LU_OUT,'(2X,2ES16.6)')PAR(1:2)
	    DO J=1,NUM_GAUS
	      I=INDX_VEC(J)
	      K=2+(I-1)*4+1
	      T1=2.99794D+05*SIM(1,K+1)/SIM(1,K)
	      T2=2.0D0; FWHM=2.0D0*T1*(LOG(T2))**(1.0D0/SIM(I,K+3))
	      WRITE(LU_OUT,'(2X,ES16.6,5ES16.4,2F10.2)')SIM(I,K),SIM(1,K+2),SIM(1,K+1),SIM(I,K+3),
	1            T1/SQRT(T2),FWHM,EW(I),EW_ERROR(I)
	      EW_TABLE(I,IP)=EW(I)
	      LAM_TABLE(I,IP)=SIM(I,K)
	    END DO
	    FLUSH(UNIT=LU_OUT)
	    IF(VERBOSE)THEN
	      WRITE(6,*)'Final SUM_SQ=',SUM_SQ(1)
	    END IF
	  END DO		!Loop over plot
!
! Output EWs in an alternate, and perhaps more useful format.
!
	  WRITE(LU_ALT,'(F12.3,12X,20F10.3)')LAM_TABLE(1,1),(EW_TABLE(1,IP),IP=1,NPLTS) 
	  DO I=2,NUM_GAUS
	    T1=2.998D+05*(LAM_TABLE(I,1)-LAM_TABLE(I-1,1))/LAM_TABLE(I,1)
	    WRITE(LU_ALT,'(2F12.3,20F10.3)')LAM_TABLE(I,1),T1,(EW_TABLE(I,IP),IP=1,NPLTS) 
	  END DO
	  FLUSH(UNIT=LU_ALT)
	END DO			!loop over spectra band
	CLOSE(LU_ALT)
	CLOSE(LU_OUT)
!
	END
