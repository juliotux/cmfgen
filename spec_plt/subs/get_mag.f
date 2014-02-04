!
! Subroutine to compute synthetic magnitudes. Filter transmission
! curves, and atmospheric extension, must be contained in the data
! file FILTER_SET.
!
	SUBROUTINE GET_MAG(NU,FLUX,NCF,DIST,FILTER_SET,LU_OUT)
	IMPLICIT NONE
!
! Created 19-Mar-2004
!
	INTEGER NCF
	REAL*8 NU(NCF)
	REAL*8 FLUX(NCF)
	REAL*8 DIST
!
	INTEGER, PARAMETER :: N_LOC=5000
	REAL*8 REV_NU(MAX(NCF,N_LOC))
	REAL*8 REV_FLUX(MAX(NCF,N_LOC))
	REAL*8 RESP(MAX(NCF,N_LOC))
	REAL*8 TRANS(MAX(NCF,N_LOC))
	CHARACTER*(*) FILTER_SET
!
! These desribe the filter response curves, and the atmospheric transmission curve.
!
	INTEGER, PARAMETER :: NF_MAX=1000
	REAL*8 FILT_LAM(NF_MAX)
	REAL*8 FILT_FREQ(NF_MAX)
	REAL*8 FILT_RESP(NF_MAX)
	REAL*8 ATM_LAM(NF_MAX)
	REAL*8 ATM_TRANS(NF_MAX)
	REAL*8 ATM_FREQ(NF_MAX)
!
	REAL*8 T1,T2
	REAL*8 MAG
	REAL*8 NORM
	REAL*8 FILT_ZP 
!
	INTEGER LU_OUT
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER I,J,K,L
	INTEGER ML,ML_BEG,ML_END
	INTEGER IOS
!
	INTEGER NUM_FILTERS
	INTEGER NX			!Revised number of points if flux grid
	INTEGER NF			!Number of points in filter response curve
	INTEGER N_ATM			!Number of points in atmopshere extinction curve
	LOGICAL POINTS_ADDED		!Indicates whether resolution of flux grid has been changed.
	LOGICAL INCLUDE_EFF_ATM		!Take into account atmospheric extension.
!
	CHARACTER(LEN=10) FILT_ID
	CHARACTER(LEN=132) STRING
!
	WRITE(6,*)NU(1),NU(NCF)
!
! Get filter and atmospheric extinction data. The zero points must also
! be in this file.
!
	CALL GEN_ASCI_OPEN(LU_IN,FILTER_SET,'old',' ',' ',IZERO,IOS)
	STRING='!'
	DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	   READ(LU_IN,'(A)')STRING
	END DO
	DO WHILE(INDEX(STRING,'!Number of filters') .EQ. 0)
	   READ(LU_IN,'(A)')STRING
	END DO
	READ(STRING,*)NUM_FILTERS
	WRITE(T_OUT,*)'Number of filters is'
!
	WRITE(LU_OUT,'(A)')' '
	WRITE(LU_OUT,'(3(1X,A9,3X))')'Filter ID','Filter ZP','Magnitude'
!
! Loop over all filters. We first read in paperameters describing the filter:
! The filter identifcation, its zero point, and the number of points in the
! response curve. A variable also indicates whether atmospheric data is incuded.
!
	DO L=1,NUM_FILTERS
!
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	     READ(LU_IN,'(A)')STRING
	  END DO
	  READ(STRING,'(A)')FILT_ID
	  FILT_ID=ADJUSTL(FILT_ID)	!Adjust string left and clean additional characters.
	  I=INDEX(FILT_ID,'  ')
	  IF(I .NE. 0)FILT_ID(I:)=' '
!
	  READ(LU_IN,*)NF
	  READ(LU_IN,*)FILT_ZP
	  READ(LU_IN,*)INCLUDE_EFF_ATM
!
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,*)'Filter ID is: ',TRIM(FILT_ID)
	  WRITE(T_OUT,*)'Number of points in response curve is',NF
	  WRITE(T_OUT,*)'Filter zero point is',FILT_ZP
	  IF(INCLUDE_EFF_ATM)
	1     WRITE(T_OUT,*)'Filter response will be corrected for atmospheric extinction'
!
	  DO ML=1,NF
	    READ(LU_IN,*)FILT_LAM(ML),FILT_RESP(ML)
	    FILT_FREQ(ML)=0.2998D+04/FILT_LAM(ML)
	  END DO
!
	  IF(INCLUDE_EFF_ATM)THEN
	    READ(LU_IN,*)N_ATM
	    IF(N_ATM .GT. N_LOC)THEN
	      WRITE(T_OUT,*)'Error in GET_MAG --- N_LOC too small for extinction data'
	      WRITE(T_OUT,*)'N_LOC=',N_LOC
	      WRITE(T_OUT,*)'N_ATM=',N_ATM
	      CLOSE(LU_IN)
	      RETURN
	    END IF
	    DO ML=1,N_ATM
	      READ(LU_IN,*)ATM_LAM(ML),ATM_TRANS(ML)
	      ATM_FREQ(ML)=0.2998D+04/ATM_LAM(ML)
	    END DO
	  END IF
!	
	  ML=1
	  DO WHILE(NU(ML) .GT. FILT_FREQ(1))
	    ML=ML+1
	  END DO
	  ML_BEG=ML-1
	  DO WHILE(NU(ML) .GT. FILT_FREQ(NF))
	    ML=ML+1
	  END DO
	  ML_END=ML
!
	  FILT_FREQ(2:NF+1)=FILT_FREQ(1:NF)
	  FILT_RESP(2:NF+1)=FILT_RESP(1:NF)
	  NF=NF+2
	  FILT_RESP(1)=0.0D0; FILT_RESP(NF)=0.0D0
	  FILT_FREQ(1)=NU(ML_BEG); FILT_FREQ(NF)=NU(ML_END)
!
	  DO I=1,NF
	    WRITE(25,*)I,FILT_FREQ(I),FILT_RESP(I)
	  END DO
	  FLUSH(UNIT=25)
!
! Check if stellar data is tabulted fine enough. This will mainly be
! a problem for continuum data, and possibly for the VEGA data
! used to check the program.
!
	  K=1
	  REV_NU(1)=NU(ML_BEG)
	  POINTS_ADDED=.FALSE.
	  DO ML=ML_BEG,ML_END-1
	    IF( 2.998D+05*(NU(ML)/NU(ML+1)-1.0D0) .GT. 100.0D0)THEN
	      J=2.998D+05*(NU(ML)/NU(ML+1)-1.0D0)/100.0D0+1
	      T1=(NU(ML)-NU(ML+1))/J
	      DO I=1,J-1
	        K=K+1
	        REV_NU(K)=REV_NU(K-1)-T1
	      END DO
	      POINTS_ADDED=.TRUE.
	    END IF
	    K=K+1
	    REV_NU(K)=NU(ML+1)
	  END DO
	  NX=K
	  DO I=1,NX-1
	   IF(REV_NU(I) .LE. REV_NU(I+1))THEN
	     WRITE(T_OUT,*)'Error - invalid frequency grid order'
	     WRITE(T_OUT,*)I,REV_NU(I),REV_NU(I+1)
	     CLOSE(LU_IN)
	     RETURN
	   END IF
	  END DO
!
	  IF(POINTS_ADDED)THEN
	    CALL MON_INTERP(REV_FLUX,NX,IONE,REV_NU,NX,FLUX,NCF,NU,NCF)
!	    IF(L .EQ. 6)THEN
!	      CALL DP_CURVE(NCF,NU,FLUX)
!	      CALL DP_CURVE(NX,REV_NU,REV_FLUX)
!	    END IF
	    WRITE(T_OUT,*)'Inserted extra points in flux grid for more accuracy'
	  ELSE
	    REV_FLUX(1:NX)=FLUX(ML_BEG:ML_END)	    
	  END IF
!
	  DO I=1,NX
	    WRITE(21,*)I,REV_NU(I),REV_FLUX(I)
	  END DO
!
	  WRITE(T_OUT,*)'Updating filter response function'
	  RESP(1:NX)=0.0D0
	  CALL MON_INTERP(RESP,NX,IONE,REV_NU,NX,FILT_RESP,NF,FILT_FREQ,NF)
	  TRANS(1:NX)=1.0D0
	  IF(INCLUDE_EFF_ATM)THEN
	    CALL MON_INTERP(TRANS,NX,IONE,REV_NU,NX,
	1                     ATM_TRANS,N_ATM,ATM_FREQ,N_ATM)
	  END IF
!
!	  IF(L .EQ. 6)THEN
!	      CALL DP_CURVE(NF,FILT_FREQ,FILT_RESP)
!	      CALL DP_CURVE(NX,REV_NU,RESP)
!	  END IF
!
	  MAG=0.0D0 
	  NORM=0.0D0
	  RESP(1:NX)=RESP(1:NX)*TRANS(1:NX)
	  REV_FLUX(1:NX)=REV_FLUX(1:NX)*RESP(1:NX)
	  DO ML=1,NX-1
            MAG=MAG+0.5D0*(REV_NU(ML)-REV_NU(ML+1))*(REV_FLUX(ML+1) + REV_FLUX(ML))
            NORM=NORM+0.5D0*(REV_NU(ML)-REV_NU(ML+1))*(RESP(ML+1) + RESP(ML))
	  END DO
	  MAG=MAG/NORM
!
          MAG=5.0*LOG10(DIST)-2.5D0*LOG10(MAG)+FILT_ZP
          WRITE(LU_OUT,'(1X,A,T10,4X,F9.3,4X,F9.3)')TRIM(FILT_ID),FILT_ZP,MAG
!
	END DO
	CLOSE(LU_IN)
!
	RETURN
	END
