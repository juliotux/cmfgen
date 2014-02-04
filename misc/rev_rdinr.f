!
! Simple program to modify a file in departure coeficient format, and modify
! for use with RDINR. IF THEN ELSE structure, so additional options can easily
! be included.
! Current options:
!                 DOUB:        Doubel radius grid
!                 EXTR         Extend grid in R
!                 FG:          Finer grid over specified interval.
!                 FG_IB:       Finer grid near inner boundary.
!                 FG_OB:       Finer grid near outer boundary.
!                 IR           Insert additional points betwen I=FST and I=LST
!                 SCAL_R:      Scale radius grid
!
	PROGRAM REV_RDINR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created : 23-Jan-2006
! Altered : 24-Sep-2006 --- FG option added.
! Altered : 26-Nov-2007 --- EXTR option installed.
!
	INTEGER, PARAMETER :: MAX_ND=500
	REAL*8 R(MAX_ND)
	REAL*8 DI(MAX_ND)
	REAL*8 ED(MAX_ND)
	REAL*8 T(MAX_ND)
	REAL*8 IRAT(MAX_ND)
	REAL*8 VEL(MAX_ND)
	REAL*8 CLUMP_FAC(MAX_ND)
	REAL*8 RTMP(MAX_ND)
	REAL*8 OLD_XV(MAX_ND)
	REAL*8 NEW_XV(MAX_ND)
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER I,J,K
	INTEGER FST,LST
	INTEGER IST,IEND
	INTEGER NI,NG
	INTEGER ICOUNT
	INTEGER N,ND,NEW_ND
	INTEGER, PARAMETER :: IONE=1
	REAL*8 RMIN,RMAX,LUM
	REAL*8 GRID_FACTOR,GRID_RATIO
	REAL*8 T1,T2,SCALE_FACTOR,DELR
	CHARACTER*132 STRING
	CHARACTER*132 FILE_IN
	CHARACTER*132 FILE_OUT
!
	CHARACTER*20 OPTION
!
	FILE_IN='HIOUT'; FILE_OUT='RDINR'
	CALL GEN_IN(FILE_IN,'Input file: DC file format')
	CALL GEN_IN(FILE_OUT,'Output file: DC file format')
!
	WRITE(6,'(A)')'Current available options are:'
	WRITE(6,'(A)')'   DOUB:        Doubel radius grid'
	WRITE(6,'(A)')'   EXTR:        Extend grid to larger radii'
	WRITE(6,'(A)')'   FG:          Fine grid over specified range'
        WRITE(6,'(A)')'   FG_IB:       Finer grid near inner boundary'
        WRITE(6,'(A)')'   FG_OB:       Finer grid near outer boundary'
	WRITE(6,'(A)')'   IR           Insert additional points betwen I=FST and I=LST'
	WRITE(6,'(A)')'   SCALE_R:     Scale radius grid'
	OPTION='FG_OB'
	CALL GEN_IN(OPTION,'Action to be taken:')
	
	OPEN(UNIT=9,FILE=FILE_IN,STATUS='OLD',ACTION='READ')
	OPEN(UNIT=10,FILE=FILE_OUT,STATUS='UNKNOWN',ACTION='WRITE')
	DO I=1,3
	  READ(9,'(A)')STRING
	  WRITE(10,'(A)')TRIM(STRING)
	END DO
	READ(9,'(A)')STRING
	READ(STRING,*)RMIN,LUM,N,ND
	READ(9,'(A)')STRING                 !Final blank line
!
	DO I=1,ND
	  READ(9,'(A)')STRING
	  READ(STRING,*)R(I),DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I)
	  DO WHILE(STRING .NE. ' ')
	    READ(9,'(A)',END=100)STRING
	  END DO
	END DO
100	CONTINUE
!
	CALL SET_CASE_UP(OPTION,IZERO,IZERO)
	IF(OPTION .EQ. 'DOUB')THEN
	  NEW_ND=2*ND-1
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,NEW_ND
	  WRITE(10,'(A)')' '
	  I=1
	  WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	  WRITE(10,'(F7.1)')1.0D0
!
	  DO I=2,ND 
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')SQRT(R(I-1)*R(I)),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
	ELSE IF(OPTION .EQ. 'SCALE_R')THEN
	  T1=0.0D0
	  CALL GEN_IN(T1,'New RMIN')
	  IF(T1 .NE. 0.0D0)THEN
	    SCALE_FACTOR=T1/RMIN
	  ELSE
	    SCALE_FACTOR=1.0D0
	    CALL GEN_IN(SCALE_FACTOR,'Factor to scale RMIN by')
	  END IF
	  DO I=1,ND
	    R(I)=SCALE_FACTOR*R(I)
	  END DO
	  RMIN=RMIN*SCALE_FACTOR
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,ND
	  DO I=1,ND 
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
! Make a fine grid at certain depths.
!
	ELSE IF(OPTION .EQ. 'IR')THEN
	  WRITE(6,*)'Input range of depths to be revised'
	  WRITE(6,*)'You specify the number of points for each interval'
	  WRITE(6,*)'Use the FG option to uniformly improve a band'
	  FST=1; LST=ND
	  CALL GEN_IN(FST,'First grid point for revision')
	  CALL GEN_IN(LST,'Last grid point for revision')
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=FST,LST
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I)-R(I+1),(R(I+2)-R(I+1))/(R(I+1)-R(I))
	  END DO
!
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,ND+NG
	  DO I=1,FST
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
	  ICOUNT=FST
	  NEW_ND=ND
!
! Now do the insertion.
!
	  NG=1
	  DO I=FST,LST-1
	    WRITE(6,'(I3,2ES16.8)')I,R(I),R(I+1)
	    CALL GEN_IN(NG,'Number of additinal grid points for this interval')
	    DO J=1,NG
	      ICOUNT=ICOUNT+1
	      T1=R(I)-J*(R(I)-R(I+1))/(NG+1)
	      READ(5,*)T1
	      WRITE(10,'(A)')' '
	      WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')T1,
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),ICOUNT
	      WRITE(10,'(F7.1)')1.0D0
	    END DO
	    NEW_ND=NEW_ND+NG
	    J=I+1
	    ICOUNT=ICOUNT+1
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(J),
	1                DI(J),ED(J),T(J),IRAT(J),VEL(J),CLUMP_FAC(J),ICOUNT
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
! Output remaining grid.
!
	  DO I=LST+1,ND 
	    ICOUNT=ICOUNT+1
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),ICOUNT
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
	  WRITE(6,*)'New number of depth points is:',NEW_ND,ICOUNT
	  CALL GEN_IN(NG,'Number of additinal grid points for this interval')
!

	ELSE IF(OPTION .EQ. 'EXTR')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)' We currently assume that at outer boundary V is close to Vinf'
	  WRITE(6,*)' Grid is extended logarithmically in Log r'
	  WRITE(6,*)' '
!
	  RMAX=2.0D0
	  CALL GEN_IN(RMAX,'Factor to extend RMAX by')
          RMAX=RMAX*R(1)
!
	  WRITE(6,*)' '
	  WRITE(6,*)'Outputing old grid near outer boundary.'
	  WRITE(6,*)'Use this information to choose depth at which extended grid is attached.'
	  WRITE(6,*)'Ratio at attachment depth will be default spacing ratio.'
	  WRITE(6,*)' '
	  WRITE(6,*)' Depth       R Ratio'
	  DO I=1,10
	    WRITE(6,'(1X,I5,ES14.4)')I,R(I)/R(I+1)
	  END DO
!
	  IST=5
	  CALL GEN_IN(IST,'Depth index to begin new grid')
	  GRID_RATIO=R(IST)/R(IST+1)
	  CALL GEN_IN(GRID_RATIO,'Scale factor --- will use logarithimic spacing in R')
!
! Reverse grid so easier to extend.
!
	  RTMP(1:ND)=R(1:ND)
	  DO I=IST,ND
	    R(ND-I+1)=RTMP(I)
	  END DO
	  NEW_ND=ND-IST+1
!
! Now extend the grid.
!
	  I=NEW_ND
	  DO WHILE(R(I)*GRID_RATIO .LT. RMAX)
	    I=I+1
	    R(I)=R(I-1)*GRID_RATIO
	  END DO
	  IF( R(I)*(1.0D0+(GRID_RATIO-1.0D0)/3) .GT. RMAX)I=I-1
	  NEW_ND=I
!
	  NI=4
	  CALL GEN_IN(NI,'Number of points used at outer boundary to refine grid')
!
	  DO J=1,NI-1
	    R(NEW_ND+J)=R(NEW_ND+J-1)+0.6D0*(RMAX-R(NEW_ND+J-1))
	  END DO
	  NEW_ND=NEW_ND+NI
	  R(NEW_ND)=RMAX-0.01D0*(RMAX-R(NEW_ND-1))
	  NEW_ND=NEW_ND+1
	  R(NEW_ND)=RMAX
!
! Reverse grid to conventional form.
!
	  DO I=1,NEW_ND/2
	    T1=R(I)
	    R(I)=R(NEW_ND-I+1)
	    R(NEW_ND-I+1)=T1
	  END DO
!
! We simply scale ED and DI for illustration purposes only --- only the R grid is 
! important.
!
          WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,NEW_ND
          DO I=1,NEW_ND-(ND+1-IST)
            T1=(RTMP(1)/R(I))**2
	    WRITE(10,'(A)')' '
            WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                T1*DI(1),T1*ED(1),T(1),IRAT(1),VEL(1),CLUMP_FAC(1),I
            WRITE(10,'(F7.1)')1.0D0
	  END DO
	  DO I=IST,ND
            WRITE(10,'(A)')' '
            WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')RTMP(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I+(NEW_ND-ND)
            WRITE(10,'(F7.1)')1.0D0
          END DO
!
	ELSE IF(OPTION .EQ. 'IV' .OR. OPTION .EQ. 'FG' .OR. 
	1       OPTION .EQ. 'ID' .OR. 
	1       OPTION .EQ. 'FR' .OR. OPTION .EQ. 'LOG_FR')THEN
	  IST=1; IEND=ND
	  CALL GEN_IN(IST,'Start index for fine grid (unaltered)')
	  CALL GEN_IN(IEND,'End index for fine gridi (unaltered)')
	  WRITE(6,*)'Number of points in the interval is',IEND-IST-1
	  IF(OPTION .NE. 'IV' .AND. OPTION .NE. 'ID')THEN
	     NG=IEND-IST
	     CALL GEN_IN(NG,'New number of grid points for this interval')
	     NEW_ND=IST+NG+(ND-IEND)+1
	  END IF
!
	  IF(OPTION .EQ. 'IV')THEN
	    DO I=1,IST
	      NEW_XV(I)=VEL(I)
	    END DO
	    I=IST
	    DO WHILE(1 .EQ. 1)
	      NEW_XV(I+1)=0.0D0
	      CALL GEN_IN(NEW_XV(I+1),'Next velocity value')
	      IF(NEW_XV(I+1) .EQ. 0.0D0)EXIT
	      I=I+1
	    END DO
	    NEW_ND=I+ND-IEND+1
	    NG=I-IST
	    DO J=IEND,ND
	      I=I+1
	      NEW_XV(I)=VEL(J)
	    END DO
	    CALL MON_INTERP(RTMP,NEW_ND,IONE,NEW_XV,NEW_ND,R,ND,VEL,ND)
	    OLD_XV(1:ND)=VEL(1:ND)
	    WRITE(6,'(7ES12.4)')(NEW_XV(I),I=1,NEW_ND)
!
	  ELSE IF(OPTION .EQ. 'FG')THEN
	    DO I=1,IST
	      NEW_XV(I)=I
	    END DO
	    T1=(IEND-IST)/(NG+1.0D0)
	    DO I=1,NG
	      NEW_XV(IST+I)=IST+I*T1
	    END DO
	    DO I=IEND,ND
	      NEW_XV(IST+NG+I+1-IEND)=I
	    END DO
	   DO I=1,ND
	     OLD_XV(I)=I
	   END DO
	   CALL MON_INTERP(RTMP,NEW_ND,IONE,NEW_XV,NEW_ND,R,ND,OLD_XV,ND)
!
	  ELSE IF(OPTION .EQ. 'ID')THEN
	    DO I=1,IST
	      NEW_XV(I)=I
	    END DO
	    NEW_XV(IST+1)=IST+0.66D0
	    NEW_XV(IST+2)=IST+1.33D0
	    NEW_XV(IST+3)=IST+2.0D0
	    K=IEND-IST-5
!
	    DO I=1,2*K+1
	      NEW_XV(IST+3+I)=IST+2+0.5D0*I
	    END DO
	    I=IST+ 4 + 2*K+1
	    NEW_XV(I)=IEND-2
	    NEW_XV(I+1)=IEND-1.33D0
	    NEW_XV(I+2)=IEND-0.66D0
!
	    K=I+3    
	    DO I=IEND,ND
	      NEW_XV(K+I-IEND)=I
	    END DO
	    NEW_ND=K+ND-IEND
	    DO I=1,ND
	      OLD_XV(I)=I
	    END DO
	    CALL MON_INTERP(RTMP,NEW_ND,IONE,NEW_XV,NEW_ND,R,ND,OLD_XV,ND)
!
	 ELSE  IF(OPTION .EQ. 'LOG_FR')THEN
	    DO I=1,IST
	      NEW_XV(I)=R(I)
	    END DO
	    T1=EXP( LOG(R(IEND)/R(IST))/(NG+1.0D0) )
	    DO I=1,NG
	      NEW_XV(IST+I)=NEW_XV(IST+I-1)*T1
	    END DO
	    DO I=IEND,ND
	      NEW_XV(IST+NG+I+1-IEND)=R(I)
	    END DO
	    OLD_XV(1:ND)=R(1:ND)
	    RTMP(1:NEW_ND)=NEW_XV(1:NEW_ND)
	 ELSE
	    DO I=1,IST
	      NEW_XV(I)=R(I)
	    END DO
	    T1=(R(IEND)-R(IST))/(NG+1.0D0)
	    DO I=1,NG
	      NEW_XV(IST+I)=R(IST)+I*T1
	    END DO
	    DO I=IEND,ND
	      NEW_XV(IST+NG+I+1-IEND)=R(I)
	    END DO
	    OLD_XV(1:ND)=R(1:ND)
	    RTMP(1:NEW_ND)=NEW_XV(1:NEW_ND)
	  END IF
!
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,NEW_ND
	  DO I=1,NEW_ND
	    WRITE(10,'(A)')' '
	    IF(I .GT. IST .AND. I .LE. IST+NG)THEN
	      J=IST
	      DO WHILE( (NEW_XV(I)-OLD_XV(J))*(OLD_XV(J+1)-NEW_XV(I)) .LT. 0.0D0)
	        J=J+1
	      END DO
	      T1=(NEW_XV(I)-OLD_XV(J))/(OLD_XV(J+1)-OLD_XV(J))
	      T2=1.0D0-T1
	      WRITE(6,*)I,T1,T2,NEW_XV(I),OLD_XV(IEND),OLD_XV(IST)
	      WRITE(6,*)VEL(J),VEL(J+1),T2*VEL(J)+T1*VEL(J+1)
	      WRITE(10,'(1X,1P,E18.10,4E15.5,E17.7,E15.5,2X,I4,A1)')RTMP(I),
	1                T2*DI(J)+T1*DI(J+1),T2*ED(J)+T1*ED(J+1),
	1                T2*T(J)+T1*T(J+1),  T2*IRAT(J)+T1*IRAT(J+1),
	1                T2*VEL(J)+T1*VEL(J+1),
	1                T2*CLUMP_FAC(J)+T1*CLUMP_FAC(J+1),I
	    ELSE
	      J=I; IF(I .GT. IST)J=I-IST-NG+IEND-1
	      WRITE(10,'(1X,1P,E18.10,4E15.5,E17.7,E15.5,2X,I4,A1)')RTMP(I),
	1                DI(J),ED(J),T(J),IRAT(J),VEL(J),CLUMP_FAC(J),I
	    END IF
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
	ELSE IF(OPTION .EQ. 'FG_IB')THEN
	  WRITE(6,*)'Current grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=ND-6,ND
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I-1)-R(I),(R(I-2)-R(I-1))/(R(I-1)-R(I))
	  END DO
	  NI=ND-1; NG=3; GRID_RATIO=1.5
	  CALL GEN_IN(NI,'1st depth to replace (e.g. ND-1):')
	  CALL GEN_IN(NG,'Numer of additional points:')
	  CALL GEN_IN(GRID_RATIO,'Ratio of succesive interval sizes: > 1:')
!
	  GRID_RATIO=1.0D0/GRID_RATIO
	  RTMP(1:ND)=R(1:ND)
	  DELR=(R(NI-1)-R(ND))*(1.0D0-GRID_RATIO)/(1.0D0-GRID_RATIO**(ND-NI+NG+1))
	  DO I=NI,ND+NG-1
	    RTMP(I)=RTMP(I-1)-DELR
	    DELR=DELR*GRID_RATIO
	  END DO
	  RTMP(ND+NG)=R(ND)
!
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,ND+NG
	  DO I=1,ND+NG
	    J=MIN(I,ND)
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')RTMP(I),
	1                DI(J),ED(J),T(J),IRAT(J),VEL(J),CLUMP_FAC(J),I
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
	  DO I=NI-2,ND+NG
	    WRITE(6,'(I3,4ES18.8)')I,RTMP(I),RTMP(I+1),
	1              RTMP(I-1)-RTMP(I),(RTMP(I-2)-RTMP(I-1))/(RTMP(I-1)-RTMP(I))
	  END DO
!
! Make a fine grid near the outer boundary.
!
	ELSE IF(OPTION .EQ. 'FG_OB')THEN
	  WRITE(6,*)'Current grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=1,12
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I)-R(I+1),(R(I+2)-R(I+1))/(R(I+1)-R(I))
	  END DO
	  NI=2; NG=3
	  CALL GEN_IN(NI,'Last depth to replace (e.g. 2):'); NI=NI-1 
	  CALL GEN_IN(NG,'Numer of additional points:')
	  CALL GEN_IN(GRID_RATIO,'Ratio of succesive interval sizes: > 1:')
	  GRID_RATIO=1.0D0/GRID_RATIO
	  WRITE(6,*)GRID_RATIO
!	  GRID_FACTOR=(GRID_RATIO-1)/GRID_RATIO
	  WRITE(10,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,ND+NG
	  I=1
	  WRITE(10,'(A)')' '
	  WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I
	  WRITE(10,'(F7.1)')1.0D0
!
! Now do the insertion.
!
	  T1=R(NI+2)
	  WRITE(6,'(I4,ES18.8)')1,R(1)
	  DELR=R(NI+2)-R(NI+3)
	  DO I=NI+NG,1,-1
	    T2=(R(1)-T1)*(1.0D0-GRID_RATIO)/(1.0D0-GRID_RATIO**(I+1))
	    DELR=MIN(T2,1.2D0*DELR)
	    T1=T1+DELR
	    RTMP(I)=T1
	    WRITE(6,*)DELR,I,RTMP(I)
	  END DO
	  DO I=1,NI+NG
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')RTMP(I),
	1                DI(2),ED(2),T(2),IRAT(2),VEL(2),CLUMP_FAC(2),I+1
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
!
! Output remaining grid.
!
	  DO I=NI+2,ND 
	    WRITE(10,'(A)')' '
	    WRITE(10,'(1X,1P,E18.10,6E15.5,2X,I4,A1)')R(I),
	1                DI(I),ED(I),T(I),IRAT(I),VEL(I),CLUMP_FAC(I),I+NG
	    WRITE(10,'(F7.1)')1.0D0
	  END DO
	  R(NI+NG+2:ND)=R(NI+2:ND)
	  R(2:NI+NG+1)=RTMP(1:NI+NG)
	  WRITE(6,*)' '
	  WRITE(6,*)'New grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=1,NI+NG+4
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I)-R(I+1),(R(I+2)-R(I+1))/(R(I+1)-R(I))
	  END DO
	ELSE
	  WRITE(6,*)TRIM(OPTION),' not recognized as valid option.'
	END IF
!
	STOP
	END
