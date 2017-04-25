!Subroutine to compute the equivalent with for a list of lines
!the list come from a file.
!Also, it will be possible to compute the chi square for
!each line for all lines (in the future...) Raul E. Puebla 27-May-2015
!
!Update 21-Jul-2015:  including the "statistics" SUM1 and SUM2
!		      SUM1 = sum(EW_obs-EW_mod)^2, sum is done over every line considered
!		      SUM2 = sum(1-EW_obs/EW_mod)^2, sum is done over every line considered
!
	SUBROUTINE EW_GROUP(NLINES,NAMEIONS,L1,L2,EW,
	1                   CENTROID)
!
	USE MOD_CURVE_DATA
!
	IMPLICIT NONE
	
	INTEGER NDEC,GET_INDX_SP,GET_INDX_DP
	EXTERNAL GET_INDX_SP,GET_INDX_DP
	INTEGER NLINES                     	!Number of lines considered
	CHARACTER(LEN=*) NAMEIONS(NLINES)  	!line names
	REAL*4 L1(NLINES),L2(NLINES)            !Arrays for starting and final lines to measure EW
	REAL*4 EW(NLINES),CENTROID(NLINES)	!Array for EW and line Centroid values
!	REAL*4 XPAR(2),YPAR(2),XT(2),YT(2)
	INTEGER L_CHAN(2)			!Channels corresponding to Wo and Wf in the curves loaded
	REAL*4 T1,T3
	INTEGER IP,I,IL
!	REAL*4 MATRIX(NLINES,NPLTS)
	REAL*4 XL1(2*NLINES)			!wavelengths corresponding to the channels.
	REAL*4 YL1(2*NLINES)			!fluxes corresponding to those wavelengths.
	REAL*4 EWOBS(NLINES)			!Array of EWs for each line of data.
	REAL*4 SUM1,SUM2			!"statistics" computed to evaluate the best sum EW among models
	CHARACTER*80 FILENAME
	INTEGER T_OUT
!	
!
	T1=0.5
	T3=1.0
	T_OUT=6
	write(6,*)'nplts: ',nplts
	OPEN(UNIT=22,FILE='EW_VALUES',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	   WRITE(T_OUT,*)'Error opening file: EW_VALUES in EW_GROUP'
	   STOP
	END IF
	OPEN(UNIT=23,FILE='S_STATISTICS',STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	   WRITE(T_OUT,*)'Error opening file: S_STATISTICS in EW_GROUP'
	   STOP
	END IF
	DO IP=1,NPLTS
          SUM1=0.0E0
	  SUM2=0.0E0
	  DO IL=1,NLINES
	    EW(IL)=0.0E0
	    CENTROID(IL)=0.0E0
	    L_CHAN(1)=GET_INDX_SP(L1(IL),CD(IP)%XVEC,NPTS(IP))  !Get channel indexes for each Wo
	    L_CHAN(2)=GET_INDX_SP(L2(IL),CD(IP)%XVEC,NPTS(IP))  !Get channel indexes for each Wf
	    IF (IP .EQ. 1) THEN
	    	  XL1(2*IL-1)=CD(IP)%XVEC(L_CHAN(1))  		!Get the corresponding wavelengths for each curve
	          XL1(2*IL)=CD(IP)%XVEC(L_CHAN(2))          	!Get the corresponding wavelengths for each curve
	          YL1(2*IL-1)=CD(IP)%DATA(L_CHAN(1))		!Get the corresponding flux for each curve
	          YL1(2*IL)=CD(IP)%DATA(L_CHAN(2))		!Get the corresponding flux for each curve
	    ENDIF
! 
!  Compute de equivalent width and the centroid for each line.
!
	    DO I=L_CHAN(1),L_CHAN(2)-1
	       EW(IL)=EW(IL)+0.5*( (CD(IP)%DATA(I)-1.0) + (CD(IP)%DATA(I+1)-1.0) )*
	1            (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	       CENTROID(IL)=CENTROID(IL)+0.5*( CD(IP)%XVEC(I)*(CD(IP)%DATA(I)-1.0)
	1            + CD(IP)%XVEC(I+1)*(CD(IP)%DATA(I+1)-1.0) )*
	1            (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	       write(20,*)ew(il),IL 
	    END DO
            IF (IP .EQ. 1) EWOBS(IL)=EW(IL)   !store the EW from data
!
! Comupte the "statistics", SUM1 and SUM2
!
	    IF (IP .GT. 1) THEN 
	      WRITE(19,*)(EWOBS(IL)-EW(IL))**2,(1.0-EWOBS(IL)/EW(IL))**2,EWOBS(IL),IL,IP
	      SUM1=SUM1+(EWOBS(IL)-EW(IL))**2
	      SUM2=SUM2+(1.0-EWOBS(IL)/EW(IL))**2	      
!	      SUM2=SUM2+(1.0-EW(IL)/EWOBS(IL))**2
	    ENDIF
	    IF(EW(IL) .NE. 0.0)THEN
	       CENTROID(IL)=CENTROID(IL)/EW(IL)
	    ELSE
	       CENTROID(IL)=0.0
	    END IF
	    WRITE(22,*)NAMEIONS(IL),EW(IL),CENTROID(IL),IP,IL
	  ENDDO
	  IF (IP .GT. 1)WRITE(23,*)SUM1,SUM2,IP
	ENDDO
	CLOSE(19)
	CLOSE(22)
	CLOSE(23)
	CALL CURVE(2*NLINES,XL1,YL1)
!	write(6,*)xl1
!	write(6,*)yl1
	write(6,*)'nplts2: ',NPLTS
	close(unit=22)
	RETURN 
	END
