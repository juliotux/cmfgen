C
C Subroutine to convert from the default Jv units (ergs/cm^2/s/Hz) to other
C units: Current options are:
C        For X axis:  Ang, um, 10^15 Hz, ev, keV, km/s, Mm/s (and Log)
C        For Y axis:  Jv, vF(v), Flam (and Log0
C
	SUBROUTINE CNVRT_J_V2(XV,YV,NBB,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                    DATA_TYPE,LAMC,X_LAB,Y_LAB,X_ONLY)
	IMPLICIT NONE
!
! Altered 26-Mar-2002: Data_TYPE inserted into call.
!                      CHanged to version V2
!
	INTEGER NBB
	REAL*4 XV(NBB),YV(NBB)
	LOGICAL LOG_X,LOG_Y,X_ONLY
	CHARACTER*(*) X_UNIT,Y_PLT_OPT
	CHARACTER*(*) X_LAB,Y_LAB
	CHARACTER*(*) DATA_TYPE
C
	INTEGER I
	REAL*8 LAMC
	REAL*8 T1
C              
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
C
	REAL*8 C_CMS
	REAL*8 KEV_TO_HZ,ANG_TO_HZ
C
	C_CMS=SPEED_OF_LIGHT()
C
C Conversion factor from Kev to units of 10^15 Hz.
C Conversion factor from Angstroms to units of 10^15 Hz.
C
	KEV_TO_HZ=0.241838E+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
C
	IF(.NOT. X_ONLY)THEN
	  IF(Y_PLT_OPT .EQ. 'NU_FNU')THEN
	    DO I=1,NBB
	      T1=1.0D+15
	      YV(I)=T1*XV(I)*YV(I)
	    END DO
	    IF(DATA_TYPE .EQ. 'J')THEN
	      Y_LAB='\gnJ\d\gn\u(ergs\d cm\u-2 \ds\u-1\d)' 
	      IF(LOG_Y)Y_LAB='Log \gnJ\d\gn\u(ergs\d \ucm\u-2 \ds\u-1\d)' 
	    ELSE IF(DATA_TYPE .EQ. 'H')THEN
	      Y_LAB='\gnH\d\gn\u(ergs\d cm\u-2 \ds\u-1\d)' 
	      IF(LOG_Y)Y_LAB='Log \gnJ\d\gn\u(ergs\d \ucm\u-2 \ds\u-1\d)' 
	    ELSE 
	      WRITE(6,*)'Unrecognized DATA_TYPE for NU_FNU'
	    ELSE
	  ELSE IF(Y_PLT_OPT .EQ. 'FLAM')THEN
	    T1=1.0E+22/C_CMS	  	!1.0E+30*1.0E-08
	    DO I=1,NBB
	      YV(I)=T1*YV(I)*XV(I)*XV(I)
	    END DO
	    IF(DATA_TYPE .EQ. 'J')THEN
	      Y_LAB='J\d\gl\u(ergs\d \ucm\u-2 \ds\u-1 \d\A)' 
	      IF(LOG_Y)Y_LAB='Log J\d\gl\u(ergs\d \ucm\u-2 \ds\u-1 \d\A)' 
	    ELSE IF(DATA_TYPE .EQ. 'H')THEN
	      Y_LAB='H\d\gl\u(ergs\d \ucm\u-2 \ds\u-1 \d\A)' 
	      IF(LOG_Y)Y_LAB='Log H\d\gl\u(ergs\d \ucm\u-2 \ds\u-1 \d\A)' 
	    ELSE 
	      WRITE(6,*)'Unrecognized DATA_TYPE for FLAM'
	    ELSE
	  ELSE IF(Y_PLT_OPT .EQ. 'FNU')THEN
!
! Plots the data units as read in.
!
	    IF(DATA_TYPE .EQ. 'J')THEN
	      Y_LAB='J\d\gn\u(ergs\d \ucm\u-2 \ds\u-1 \dHz\u-1\d)' 
	      IF(LOG_Y)Y_LAB='Log J\d\gn\u(ergs\d \ucm\u-2 \ds\u-1 \dHz\u-1\d)'
	    ELSE IF(DATA_TYPE .EQ. 'H')THEN
	      Y_LAB='H\d\gn\u(ergs\d \ucm\u-2 \ds\u-1 \dHz\u-1\d)' 
	      IF(LOG_Y)Y_LAB='Log H\d\gn\u(ergs\d \ucm\u-2 \ds\u-1 \dHz\u-1\d)'
	    ELSE IF(DATA_TYPE .EQ. 'M(t)')THEN
	      Y_LAB='M(t)'
	      IF(LOG_Y)Y_LAB='Log M(t)'
	    ELSE IF(DATA_TYPE .EQ. 'ETA')THEN
	      Y_LAB='/ge(ergs\d \ucm\u-3 \ds\u-1 \dHz\u-1\d \gO\u-1\d)'
	      IF(LOG_Y)Y_LAB='Log /ge(ergs\d \ucm\u-3 \ds\u-1 \dHz\u-1\d \gO\u-1\d)'
	    ELSE IF(DATA_TYPE .EQ. 'CHI')THEN
	      Y_LAB='\gx(cm\u-1\d)'
	      IF(LOG_Y)Y_LAB='Log \gx(cm\u-1\d)'
	    ELSE 
	      WRITE(6,*)'Unrecognized DATA_TYPE for FNU'
	    ELSE
	  END IF
	END IF
C
	IF(X_UNIT .EQ. 'ANG')THEN
	  DO I=1,NBB
	    XV(I)=ANG_TO_HZ/XV(I)
	  END DO
	  X_LAB='\gl(\A)'
	  IF(LOG_X)X_LAB='Log \gl(\A)'
	ELSE IF(X_UNIT .EQ. 'UM')THEN
	  DO I=1,NBB
	    XV(I)=1.0D-04*ANG_TO_HZ/XV(I)
	  END DO
	  X_LAB='\gl(\gmm)'
	  IF(LOG_X)X_LAB='Log \gl(\gmm)'
	ELSE IF(X_UNIT .EQ. 'KEV')THEN
	  DO I=1,NBB
	    XV(I)=XV(I)/KEV_TO_HZ
	  END DO
	  X_LAB='keV'
	  IF(LOG_X)X_LAB='Log \gn(keV)'
	ELSE IF(X_UNIT .EQ. 'EV')THEN
	  DO I=1,NBB
	    XV(I)=1.0D+03*XV(I)/KEV_TO_HZ
	  END DO
	  X_LAB='\gl(eV)'
	  IF(LOG_X)X_LAB='Log \gn(eV)'
	ELSE IF(X_UNIT .EQ. 'HZ')THEN
	  X_LAB='\gn(10\u15 \dHz)'
	  IF(LOG_X)X_LAB='Log \gn(1-\u15 \dHz)'
	ELSE IF(X_UNIT .EQ. 'MM/S')THEN
	  DO I=1,NBB
	    XV(I)=1.0D-08*C_CMS*(ANG_TO_HZ/XV(I)-LAMC)/LAMC
	  END DO
	  X_LAB='V(Mm\u \ds\u-1\d)'
	  IF(LOG_X)X_LAB='Log V(Mm\u \ds\u-1\d)'
	ELSE IF(X_UNIT .EQ. 'KM/S')THEN
	  DO I=1,NBB
	    XV(I)=1.0D-05*C_CMS*(ANG_TO_HZ/XV(I)-LAMC)/LAMC
	  END DO
	  X_LAB='V(km\u \ds\u-1\d)'
	  IF(LOG_X)X_LAB='Log V(km\u \ds\u-1\d)'
	END IF
	
C
C Now take logs if required.
C
	IF(LOG_X)THEN
	  DO I=1,NBB
	    XV(I)=LOG10(XV(I))
	  END DO
	END IF
C
	IF(LOG_Y)THEN
	  DO I=1,NBB
	    IF(YV(I) .LE. 0)THEN
	      YV(I)=-50
	    ELSE
	      YV(I)=LOG10(YV(I))
	    END IF
	  END DO
	END IF
C
	RETURN
	END
