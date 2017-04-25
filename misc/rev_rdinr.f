!
! Simple program to create and improved R grid. The code takes a file in departure coeficient format,
! modifies the R grid according to various options, and then outputs a new RDINR file. A file R_GRID_CHECK
! (or R_GRID_SELECTION) is created which allows the user to check grid spacings in TAU and R.
!
! Ideally the MEANOPAC files should be available.
!
! Basically the code has an IF THEN ELSE structure, so additional options can easily be included.
!
! Current options:
!                 BDOUB:       Half grid spacing (full grid) but with allowance for boundaries.
!                 DOUB:        Half grid spacing (full grid).
!                 EXTR         Extend grid in R.
!                 ID           Half grid spacing betwen I=FST and I=LST, and adjsut adjacent depths.
!                 I3           1/3  grid spacing betwen I=FST and I=LST, and adjsut adjacent depths.
!                 FG:          Finer grid over specified interval.
!                 FG_IB:       Finer grid near inner boundary.
!                 FG_OB:       Finer grid near outer boundary.
!                 IR           Insert additional points betwen I=FST and I=LST.
!                 RTAU         Create grid, equally spaced in Log(Tau), but with constrints on dR.
!                 SCALE_R:     Scale radius grid.
!                 TAU          Insert extra depth points in Log(TAU).
!
	PROGRAM REV_RDINR
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created : 23-Jan-2006
! Altered : 24-Sep-2006 --- FG option added.
! Altered : 26-Nov-2007 --- EXTR option installed.
! Altered : 24-Dec-2013 --- TAU and RTAU options inserted.
! Altered : 05-Jan-2014 --- NG=0 allows insertion to be skipped for TAU option.
! Altered : 08-Jan-2014 --- Improvements to RTAU option and call to ADJUST_SN_R_GRID made.
! Altered : 15-Jan-2014 --- Changed ID option, cleaned, and generally use OUT__RGRD to create RDINR file.
! Altered : 19-Mar-2014 --- Tried to improve handling of boundaries with TAU option.
! Altered : 27-Mar-2014 --- Fixed bug with DOUB option (introduced when cleaning),
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: MAX_ND=500
!
	REAL*8 R(MAX_ND)
	REAL*8 DI(MAX_ND)
	REAL*8 ED(MAX_ND)
	REAL*8 T(MAX_ND)
	REAL*8 IRAT(MAX_ND)
	REAL*8 VEL(MAX_ND)
	REAL*8 CLUMP_FAC(MAX_ND)
	REAL*8 OLD_TAU(MAX_ND)
	REAL*8 TAU_SAV(MAX_ND)
!
	REAL*8 RTMP(MAX_ND)
	REAL*8 OLD_XV(MAX_ND)
	REAL*8 NEW_XV(MAX_ND)
	REAL*8 NEW_TAU(MAX_ND)
	REAL*8 NEW_VEL(MAX_ND)
	REAL*8 TAU(MAX_ND)
!
	REAL*8 RMIN,RMAX,LUM
	REAL*8 GRID_FACTOR,GRID_RATIO
	REAL*8 T1,T2,SCALE_FACTOR,DELR
	REAL*8 TAU_MIN,TAU_MAX
	REAL*8 TAU_RAT
	REAL*8 R_SCALE_FAC
	REAL*8 IB_RAT
	REAL*8 OB_RAT
	REAL*8 SCL_FAC
	REAL*8 DTAU2_ON_DTAU1
	REAL*8 dLOGT_MAX
!
	INTEGER I,J,K
	INTEGER FST,LST
	INTEGER IST,IEND
	INTEGER LU,IOS
	INTEGER NI,NG
	INTEGER NIB,NOB
	INTEGER ICOUNT
	INTEGER N,ND,NEW_ND,TAU_ND
	INTEGER ND_SAV
!
	LOGICAL RD_MEANOPAC
	LOGICAL GRID_SAT
	LOGICAL ROUND_ERROR
!
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) FILE_IN
	CHARACTER(LEN=132) FILE_OUT
	CHARACTER(LEN=20)  OPTION
!
	FILE_IN='RDINR_OLD'; FILE_OUT='RDINR'
	CALL GEN_IN(FILE_IN,'Input file: DC file format')
	CALL GEN_IN(FILE_OUT,'Output file: DC file format')
!
	WRITE(6,'(A)')'Current available options are:'
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')'   BDOUB:       Half grid spacing but allow for special spacing at boudaries'
	WRITE(6,'(A)')'   DOUB:        Half grid spacing [i.e., ND(new)=2*ND-1)'
	WRITE(6,'(A)')'   EXTR:        Extend grid to larger radii'
	WRITE(6,'(A)')'   FG:          Fine grid over specified range'
        WRITE(6,'(A)')'   FG_IB:       Finer grid near inner boundary'
        WRITE(6,'(A)')'   FG_OB:       Finer grid near outer boundary'
	WRITE(6,'(A)')'   ID           Half grid spacing betwen I=FST and I=LST with modifcation to adjacent points'
	WRITE(6,'(A)')'   IR           Insert additional points betwen I=FST and I=LST'
	WRITE(6,'(A)')'   SCALE_R:     Scale radius grid'
	WRITE(6,'(A)')'   ADDR:        Add extra R points in R space'
	WRITE(6,'(A)')'   TAU:         Refine grid and add points (multiple ranges) in TAU space'
	WRITE(6,'(A)')'   RTAU:        Refine whole grid - dTAU and dR space'
	WRITE(6,'(A)')' '
!
	OPTION='FG_OB'
	CALL GEN_IN(OPTION,'Action to be taken:')
!
! The old R-grid file is assumed to have the same format as the DC files.
!
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
! Read in optical depth scale. Needed for RTAU and TAU options. Also used
! for checking purposes (if available) for some other options.
!
! We use TAU_SAV for dTAU, and is used to increase the precision of the TAU
! scale (as insufficient digits may be print out).
!
	ROUND_ERROR=.FALSE.
	OPEN(UNIT=20,FILE='MEANOPAC',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(20,'(A)')STRING
	    DO I=1,ND
	      READ(20,*)RTMP(I),J,OLD_TAU(I),TAU_SAV(I)
	      J=MAX(I,2)
	      T1=R(J-1)-R(J)
	      IF( ABS(RTMP(I)-R(I))/T1 .GT. 2.0D-03 .AND. .NOT. ROUND_ERROR)THEN
	        WRITE(6,*)' '
	        WRITE(6,*)'Possible eror with MEANOPAC -- inconsistent R grid'
	        WRITE(6,*)'Error could simply be a lack of sig. digits in MEANOPAC'
	        WRITE(6,*)' RMO(I)=',RTMP(I)
	        WRITE(6,*)'   R(I)=',R(I)
	        WRITE(6,*)' R(I+1)=',R(I+1)
	        ROUND_ERROR=.TRUE.
	        CALL GEN_IN(ROUND_ERROR,'Continue as only rounding error?')
	        IF(.NOT. ROUND_ERROR)STOP
	      END IF
	    END DO
!	    DO I=8,1,-1
!             OLD_TAU(I)=OLD_TAU(I+1)-TAU_SAV(I)
!            END DO
 	     DO I=2,ND
              OLD_TAU(I)=OLD_TAU(I-1)+TAU_SAV(I-1)
             END DO
	    RD_MEANOPAC=.TRUE.
	  ELSE
	    RD_MEANOPAC=.FALSE.
	  END IF
	  IF(ROUND_ERROR .AND. RD_MEANOPAC)THEN
	    RTMP(1:ND)=R(1:ND)
	  END IF 
	CLOSE(UNIT=20)
!
	CALL SET_CASE_UP(OPTION,IZERO,IZERO)
!
	IF(OPTION .EQ. 'RTAU')THEN
	  IF(.NOT. RD_MEANOPAC)THEN
	    WRITE(6,*)' '
	    WRITE(6,*)'Error -- unable to open MEANOPAC which is required by the RTAU option'
	    STOP
	  END IF
!
! NB: Can get a fine grid at outer boudary by seeting SCL_FAC, or by putting NOB > 1.
! If NOB >1, it might be better to set SCL_FAC=0.0D0.
!
! NB: New TAU scale = Old TAU sclae - SCL_FAC*TAU(1)
!
! Set reasonable defaults.
!
	  NIB=2; NOB=1
	  NEW_ND=ND
	  R_SCALE_FAC=1.4
	  IB_RAT=2.0D0; OB_RAT=1.5D0
	  SCL_FAC=0.0D0
	  DTAU2_ON_DTAU1=100.0D0
	  dLOGT_MAX=0.05D0
!
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' '
	  CALL GEN_IN(NEW_ND,'Input new number of depth points')
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')' The following factor is used to adjust the importance of the dR to the change in dTAU '
	  WRITE(6,'(A)')BLUE_PEN
	  CALL GEN_IN(R_SCALE_FAC,'Factor (>1) to enhance maximum dLog(R) spacing')
	  CALL GEN_IN(dLOGT_MAX,'Maximum fractional change in the temperature: set to > 1.0 to switch T check off')
!
	  WRITE(6,'(A)')BLUE_PEN
	  CALL GEN_IN(NIB,'Number of depth points to insert at INNER boundary')
	  CALL GEN_IN(IB_RAT,'Ratio in optical depth increments at INNER boundry (>1)')
!
	  CALL GEN_IN(NOB,'Number of depth points to insert at OUTER boundary')
	  IF(NOB .NE. 1)CALL GEN_IN(OB_RAT,'Ratio in optical depth increments at OUTER boundry (>1)')
	  WRITE(6,'(A)')RED_PEN
	  WRITE(6,'(A)')' The following factor stretches the optical depth scale so that more points'
	  WRITE(6,'(A)')' are placed near the outer boundary (>=0 & < 1):'
	  WRITE(6,'(A)')BLUE_PEN
	  CALL GEN_IN(SCL_FAC,'Factor to scale optical depth at OUTER boundary: TAU=TAU-SF*TAU(1)')
	  IF(SCL_FAC .LT. 0.0D0 .OR. SCL_FAC .GE. 1.0D0)THEN
	    WRITE(6,*)'Invalid scale factor: should be >=0 and < 1.0'
	    STOP
	  END IF
	  CALL GEN_IN(DTAU2_ON_DTAU1,'~DTAU(2)/DTAU(1) at outer boudary')
	  WRITE(6,'(A)')DEF_PEN
!
	  T1=OLD_TAU(1)
	  DO I=1,ND
	    OLD_TAU(I)=OLD_TAU(I)-SCL_FAC*T1
	  END DO
	  CALL ADJUST_SN_R_GRID(RTMP,R,T,OLD_TAU,R_SCALE_FAC,dLOGT_MAX,
	1                           IB_RAT,OB_RAT,DTAU2_ON_DTAU1,NIB,NOB,NEW_ND,ND)
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
!
! Option to insert extra points equally space in LOG(TAU). Multiple regions may
! be edited at the same time.
!
	ELSE IF(OPTION .EQ. 'TAU')THEN
	  IF(.NOT. RD_MEANOPAC)THEN
	    WRITE(6,*)' '
	    WRITE(6,*)'Error -- unable to open MEANOPAC which is required by the TAU option'
	    STOP
	  END IF
!
	  WRITE(6,*)' '
	  WRITE(6,'(A,ES9.3,5X,A,E9.3,/)')' TAU(Min)=',OLD_TAU(1),'TAU(Max)=',OLD_TAU(ND)
	  WRITE(6,*)' '
	  WRITE(6,*)'You may do multiple intervals -- one at a time'
	  WRITE(6,*)'To exit, put TAU_MIN -ve (or zero)'
	  WRITE(6,*)' '
!
	  TAU_MIN=OLD_TAU(1); TAU_MAX=OLD_TAU(ND)
1000	  CALL GEN_IN(TAU_MIN,'Minimum of TAU range for revision')
	  IF(TAU_MIN .LT. OLD_TAU(1) .OR. TAU_MIN .GT. OLD_TAU(ND))THEN
	    WRITE(6,'(A)')RED_PEN
	    WRITE(6,*)'Error -- requested TAU_MIN is outside valid range'
	    WRITE(6,'(A,ES9.3,5X,A,E9.3,/)')' TAU(Min)=',OLD_TAU(1),'TAU(Max)=',OLD_TAU(ND)
	    WRITE(6,'(A)')DEF_PEN
	    GOTO 1000
	  END IF
	  NEW_ND=ND
	  NEW_TAU(1:ND)=OLD_TAU(1:ND)
	  DO WHILE(TAU_MIN .GT. 0.0D0)
2000	    CALL GEN_IN(TAU_MAX,'Maximum of Tau range for revision')
	    IF(TAU_MAX .LT. OLD_TAU(1) .OR. TAU_MAX .GT. OLD_TAU(ND))THEN
	      WRITE(6,'(A)')RED_PEN
	      WRITE(6,*)'Error -- requested TAU_MIN is outside valid range'
	      WRITE(6,'(A,ES9.3,5X,A,E9.3,/)')' TAU(Min)=',OLD_TAU(1),'TAU(Max)=',OLD_TAU(ND)
	      WRITE(6,'(A)')DEF_PEN
	      GOTO 2000
	    END IF
	    TAU_ND=NEW_ND
	    TAU(1:TAU_ND)=NEW_TAU(1:TAU_ND)
	    IST=1
	    DO I=1,TAU_ND-1
	      IF(TAU_MIN .LT. TAU(I+1))THEN
	        IST=I
	        IF( (TAU_MIN-TAU(I)) .GT. (TAU(I+1)-TAU_MIN))IST=IST+1
	        EXIT
	      END IF
	    END DO
	    DO I=1,TAU_ND-1
	      IF(TAU_MAX .LT. TAU(I+1))THEN
	        IEND=I
	        IF( (TAU_MAX-TAU(I)) .GT. (TAU(I+1)-TAU_MAX))IEND=IEND+1
	        EXIT
	      END IF
	    END DO
	    IST=MAX(IST,2); IEND=MIN(TAU_ND-1,IEND)
	    IF(IEND .EQ. IST)IEND=MIN(TAU_ND-1,IEND+1)
	    IF(IEND .EQ. IST)IST=MAX(2,IST-1)
	    TAU_MIN=TAU(IST); TAU_MAX=TAU(IEND)
	    NG=IEND-IST-1
!
	    WRITE(6,'(A)')' '
	    WRITE(6,*)IST,TAU(IST)
	    WRITE(6,*)IEND,TAU(IEND)
	    WRITE(6,*)'Number of points in the interval is ',NG
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(NG,'Enter new number of grid points for the interval')
!
	    TAU_SAV=TAU(1:NEW_ND)
	    ND_SAV=NEW_ND
	    TAU_RAT=EXP( LOG(TAU_MAX/TAU_MIN) / (NG+1) )
	    NEW_ND=IST+NG+(TAU_ND-IEND)+1
!
	    DO I=1,IST
	      NEW_TAU(I)=TAU(I)
	    END DO
	    DO I=IST+1,IST+NG
	      NEW_TAU(I)=NEW_TAU(I-1)*TAU_RAT
	    END DO
	    DO I=IST+NG+1,NEW_ND
	      NEW_TAU(I)=TAU(IEND+(I-IST-NG-1))
	    END DO
!
            IF(IST .GT. 2)THEN
              J=MAX(4,IST)
              TAU_RAT=EXP( LOG(NEW_TAU(J+2)/NEW_TAU(J-2))/4 )
              DO I=J-2,J+2
                NEW_TAU(I)=NEW_TAU(I-1)*TAU_RAT
	      END DO
	    END IF
!
            IEND=IST+NG
            IF(IEND .LT. NEW_ND-3)THEN
              J=MIN(IEND-1,ND-4)
              TAU_RAT=EXP( LOG(NEW_TAU(J+3)/NEW_TAU(J-2))/5 )
              DO I=J-2,J+2
                NEW_TAU(I)=NEW_TAU(I-1)*TAU_RAT
	      END DO
	    END IF
!
	    WRITE(6,'(4X,A1,9X,A,8X,A,3X,A)')'I','Tau','dTAU','dLog(Tau)'
	    DO I=MAX(IST-3,2),MIN(IST+4,ND-1)
	      WRITE(6,'(I5,4ES12.3)')I,NEW_TAU(I),NEW_TAU(I+1)-NEW_TAU(I),
	1            LOG(NEW_TAU(I+1)/NEW_TAU(I))
	    END DO
	    WRITE(6,'(A)')' '
	    DO I=MAX(IEND-6,2),MIN(IEND+3,ND-1)
	      WRITE(6,'(I5,4ES12.3)')I,NEW_TAU(I),NEW_TAU(I+1)-NEW_TAU(I),
	1            LOG(NEW_TAU(I+1)/NEW_TAU(I))
	    END DO
	    CALL GEN_IN(GRID_SAT,'Is grid satisfactory?')
	    IF(GRID_SAT)THEN
	    ELSE
	      NEW_ND=ND_SAV
	      NEW_TAU(1:NEW_ND)=TAU_SAV(1:NEW_ND)
	    END IF
!
! We can do anther region if desired.
!
	    TAU_MIN=0.0D0
	    CALL GEN_IN(TAU_MIN,'Minimum of TAU range for revision (=<0) to exit.')
	  END DO
!
! Now create the NEW R grid.
!
	  CALL MON_INTERP(RTMP,NEW_ND,IONE,NEW_TAU,NEW_ND,R,ND,OLD_TAU,ND)
!
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
! Option to insert extra points equally space in LOG(TAU). Multiple regions may
! be edited at the same time.
!
	ELSE IF(OPTION .EQ. 'ADDR')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)' '
	  CALL GEN_IN(IST,'Initial depth index bounding revison region')
	  CALL GEN_IN(IEND,'Last depth index bounding revison region')
	  IST=MAX(IST,2); IEND=MIN(IEND,ND)
	  NG=IEND-IST-1
!
	  IF(IEND .NE. ND)NG=(R(IST)-R(IEND))/(R(IEND)-R(IEND+1))-1
	  WRITE(6,'(A)')' '
	  WRITE(6,*)'Number of points in the interval is ',NG
	  CALL GEN_IN(NG,'Enter new number of grid points for the interval')
	  NEW_ND=IST+NG+(ND-IEND)+1
!
	  DO I=1,IST
	    RTMP(I)=R(I)
	  END DO
	  T1=(R(IST)-R(IEND))/(NG+1)
	  DO I=IST+1,IST+NG
	     RTMP(I)=RTMP(I-1)-T1
	  END DO
	  DO I=IST+NG+1,NEW_ND
	    RTMP(I)=R(IEND+(I-IST-NG-1))
	  END DO
!
! Now create the NEW R grid.
!
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
	ELSE IF(OPTION .EQ. 'BDOUB')THEN
	  NIB=2; NOB=1
	  CALL GEN_IN(NIB,'Number of depth points inserted at INNER boundary')
	  CALL GEN_IN(NOB,'Number of depth points inserted at OUTER boundary')
	  NEW_ND=2*(ND-NIB-NOB-1)+1+NIB+NOB
!
	  RTMP(1)=R(1)
	  T1=SQRT(R(1)*R(NOB+2))
	  J=1
	  DO I=2,NOB+1
	    RTMP(I)=R(1)-(R(1)-T1)*(R(1)-R(I))/(R(1)-R(NOB+2))
	  END DO
	  J=NOB+1
!
	  DO I=NOB+2,ND-NIB-1
	    J=J+1 
	    RTMP(J)=SQRT(R(I-1)*R(I))
	    J=J+1 
	    RTMP(J)=R(I)
	  END DO
!
	  J=J+1 ; I=ND-NIB
	  T1=SQRT(R(ND)*R(ND-NIB-1))
	  RTMP(J)=T1
	  DO I=NIB,1,-1
	    J=J+1
	    RTMP(J)=R(ND) + (T1-R(ND))*(R(ND-I)-R(ND))/(R(ND-NIB-1)-R(ND))
	  END DO
	  J=J+1
	  RTMP(J)=R(ND)
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
	ELSE IF(OPTION .EQ. 'DOUB')THEN
!
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'This options inserts grid points betwen pairs without regard to the boundary'
	  WRITE(6,*)DEF_PEN
!
	  J=1
	  DO I=1,ND-1
	    J=2*I-1
	    RTMP(J)=R(I)
	    J=2*I
	    RTMP(J)=SQRT(R(I)*R(I+1))
	  END DO
	  NEW_ND=2*ND-1
	  RTMP(NEW_ND)=R(ND)
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
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
	  DO I=1,FST
	    RTMP(I)=R(I)
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
	      RTMP(ICOUNT)=T1
	    END DO
	    NEW_ND=NEW_ND+NG
	    J=I+1
	    ICOUNT=ICOUNT+1
	    RTMP(ICOUNT)=R(J)
	  END DO
!
! Output remaining grid.
!
	  DO I=LST+1,ND 
	    ICOUNT=ICOUNT+1
	    RTMP(ICOUNT)=R(I)
	  END DO
	  WRITE(6,*)'New number of depth points is:',NEW_ND,ICOUNT
!
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
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
! Since grid is extened, we can't use OUT_RiGRId routine.
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
	1       OPTION .EQ. 'ID' .OR.  OPTION .EQ. 'I3' .OR. 
	1       OPTION .EQ. 'FR' .OR. OPTION .EQ. 'LOG_FR')THEN
	  IST=1; IEND=ND
	  CALL GEN_IN(IST,'Start index for fine grid (unaltered)')
	  CALL GEN_IN(IEND,'End index for fine gridi (unaltered)')
	  WRITE(6,*)'Number of points in the interval is',IEND-IST-1
	  IF(OPTION .NE. 'IV' .AND. OPTION .NE. 'ID' .AND. OPTION .NE. 'I3')THEN
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
!
	    IST=MAX(3,IST-3)
	    IEND=MIN(ND-3,IEND+3)
	    DO I=1,IST
	      NEW_XV(I)=I
	    END DO
	    T1=IST+0.9;  NEW_XV(IST+1)=T1
	    T1=T1 +0.8;  NEW_XV(IST+2)=T1
	    T1=T1 +0.7;  NEW_XV(IST+3)=T1
	    T1=T1 +0.6;  NEW_XV(IST+4)=T1
	    K=IST+4
!
	    I=IST+3; T1=I
	    DO WHILE(I .LT. IEND-4)
	      I=I+1
	      K=K+1; T1=T1+0.5; NEW_XV(K)=T1
	      K=K+1; T1=T1+0.5; NEW_XV(K)=T1
	    END DO
	    K=K+1; T1=T1+0.6; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.7; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.8; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.9; NEW_XV(K)=T1
!
	    IEND=NINT(T1)+1
	    K=K+1
	    DO I=IEND,ND
	      NEW_XV(K+I-IEND)=I
	    END DO
	    NEW_ND=K+ND-IEND
	    DO I=1,ND
	      OLD_XV(I)=I
	    END DO
	    CALL MON_INTERP(RTMP,NEW_ND,IONE,NEW_XV,NEW_ND,R,ND,OLD_XV,ND)
!
	    DO I=1,NEW_ND
	      WRITE(45,'(I4,3X,F6.2,2ES14.4)')I,NEW_XV(I),RTMP(I),RTMP(MIN(I+1,NEW_ND))-RTMP(I) 
	    END DO
	    FLUSH(UNIT=45)
!
	  ELSE IF(OPTION .EQ. 'I3')THEN
!
	    IST=MAX(3,IST-3)
	    IEND=MIN(ND-3,IEND+3)
	    DO I=1,IST
	      NEW_XV(I)=I
	    END DO
	    T1=IST+0.8;  NEW_XV(IST+1)=T1
	    T1=T1 +0.7;  NEW_XV(IST+2)=T1
	    T1=T1 +0.6;  NEW_XV(IST+3)=T1
	    T1=T1 +0.5;  NEW_XV(IST+4)=T1
	    T1=T1 +0.4;  NEW_XV(IST+5)=T1
	    K=IST+5
!
	    I=IST+3; T1=I; T2=1.0D0/3.0D0
	    DO WHILE(I .LT. IEND-4)
	      I=I+1
	      K=K+1; T1=T1+T2; NEW_XV(K)=T1
	      K=K+1; T1=T1+T2; NEW_XV(K)=T1
	      K=K+1; T1=T1+T2; NEW_XV(K)=T1
	    END DO
	    K=K+1; T1=T1+0.4; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.5; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.6; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.7; NEW_XV(K)=T1
	    K=K+1; T1=T1+0.8; NEW_XV(K)=T1
!
	    IEND=NINT(T1)+1
	    K=K+1
	    DO I=IEND,ND
	      NEW_XV(K+I-IEND)=I
	    END DO
	    NEW_ND=K+ND-IEND
	    DO I=1,ND
	      OLD_XV(I)=I
	    END DO
	    CALL MON_INTERP(RTMP,NEW_ND,IONE,NEW_XV,NEW_ND,R,ND,OLD_XV,ND)
!
	    DO I=1,NEW_ND
	      WRITE(45,'(I4,3X,F6.2,2ES14.4)')I,NEW_XV(I),RTMP(I),RTMP(MIN(I+1,NEW_ND))-RTMP(I) 
	    END DO
	    FLUSH(UNIT=45)
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
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
	ELSE IF(OPTION .EQ. 'FG_IB')THEN
	  WRITE(6,*)' '
	  WRITE(6,*)'Current grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=ND-6,ND
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I-1)-R(I),(R(I-2)-R(I-1))/(R(I-1)-R(I))
	  END DO
	  WRITE(6,*)' '
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
	  NEW_ND=ND+NG
!
	  DO I=NI-2,ND+NG
	    WRITE(6,'(I3,4ES18.8)')I,RTMP(I),RTMP(I+1),
	1              RTMP(I-1)-RTMP(I),(RTMP(I-2)-RTMP(I-1))/(RTMP(I-1)-RTMP(I))
	  END DO
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
!
! Make a fine grid near the outer boundary.
!
	ELSE IF(OPTION .EQ. 'FG_OB')THEN
	  WRITE(6,*)' '
	  WRITE(6,*)'Current grid near outer boundary'
	  WRITE(6,'(2X,A,4(8X,A6,4X))')'I',' R(I) ','R(I+1)',' Delr',' Ratio'
	  DO I=1,12
	    WRITE(6,'(I3,4ES18.8)')I,R(I),R(I+1),R(I)-R(I+1),(R(I+2)-R(I+1))/(R(I+1)-R(I))
	  END DO
	  WRITE(6,*)' '
	  NI=2; NG=3; GRID_RATIO=2.0D0
	  CALL GEN_IN(NI,'Last depth to replace (e.g. 2):'); NI=NI-1 
	  CALL GEN_IN(NG,'Numer of ADDITIONAL points:')
	  CALL GEN_IN(GRID_RATIO,'Ratio of succesive interval sizes: > 1:')
	  GRID_RATIO=1.0D0/GRID_RATIO
!
	  RTMP(1)=R(1)
!
! Now do the insertion.
!
	  T1=R(NI+2)
	  DELR=R(NI+2)-R(NI+3)
	  DO I=NI+NG,1,-1
	    T2=(R(1)-T1)*(1.0D0-GRID_RATIO)/(1.0D0-GRID_RATIO**(I+1))
	    DELR=MIN(T2,1.2D0*DELR)
	    T1=T1+DELR
	    RTMP(I+1)=T1
	  END DO
!
! Remaining grid.
!
	  DO I=NI+2,ND
	    RTMP(I+NG)=R(I) 
	  END DO
	  NEW_ND=ND+NG
	  CALL OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
	ELSE
	  WRITE(6,*)TRIM(OPTION),' not recognized as valid option.'
	END IF
!
	STOP
	END
!
! This routine:
!        (1) Checks that R is monotonic
!        (2) Computes V, T, etc on the new grid
!        (3) Outputs the revised R files.
!        (4) Ouputs R_GRID_CHK so that DTAU etc spacings can be examined.
!
	SUBROUTINE OUT_RGRID(RTMP,NEW_ND,R,DI,ED,T,IRAT,VEL,CLUMP_FAC,OLD_TAU,ND,RMIN,LUM,RD_MEANOPAC)
	IMPLICIT NONE
!
! Altered 22-Mar-2014 : T now output to R_GRID_CHK file.
!
	INTEGER ND
	INTEGER NEW_ND
!
	REAL*8 RMIN
	REAL*8 LUM
!
	REAL*8 R(ND)
	REAL*8 DI(ND)
	REAL*8 ED(ND)
	REAL*8 T(ND)
	REAL*8 IRAT(ND)
	REAL*8 VEL(ND)
	REAL*8 CLUMP_FAC(ND)
	REAL*8 OLD_TAU(ND)
!
	REAL*8 RTMP(NEW_ND)
	REAL*8 NEW_DI(NEW_ND)
	REAL*8 NEW_ED(NEW_ND)
	REAL*8 NEW_T(NEW_ND)
	REAL*8 NEW_IRAT(NEW_ND)
	REAL*8 NEW_VEL(NEW_ND)
	REAL*8 NEW_CLUMP_FAC(NEW_ND)
	REAL*8 NEW_TAU(NEW_ND)
	REAL*8 NEW_SIGMA(NEW_ND)
	REAL*8 COEF(NEW_ND,4)
!
	LOGICAL RD_MEANOPAC
!
	REAL*8 T1
	INTEGER I,J
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LU=10
!
! Check that R is monotonic.
!
	DO I=1,ND-1
	  IF(R(I+1) .GE. R(I))THEN
	    WRITE(6,*)'Error -- R grid is not monotonic at depth ',I
	    DO J=MAX(1,I-4),MIN(I+4,ND)
	      WRITE(6,*)J,R(J)
	    END DO
	    STOP
	  END IF
	END DO
!
	CALL MON_INTERP(NEW_VEL,NEW_ND,IONE,RTMP,NEW_ND,VEL,ND,R,ND)
	CALL MON_INTERP(NEW_ED,NEW_ND,IONE,RTMP,NEW_ND,ED,ND,R,ND)
	CALL MON_INTERP(NEW_T,NEW_ND,IONE,RTMP,NEW_ND,T,ND,R,ND)
	CALL MON_INTERP(NEW_DI,NEW_ND,IONE,RTMP,NEW_ND,DI,ND,R,ND)
	CALL MON_INTERP(NEW_IRAT,NEW_ND,IONE,RTMP,NEW_ND,IRAT,ND,R,ND)
	CALL MON_INTERP(NEW_CLUMP_FAC,NEW_ND,IONE,RTMP,NEW_ND,CLUMP_FAC,ND,R,ND)
!
! NB: Unit 10 must be open, and have it date header written.
!
	WRITE(LU,'(1X,ES15.7,4X,1PE11.4,5X,0P,I4,5X,I4)')RMIN,LUM,1,NEW_ND
	DO I=1,NEW_ND
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(1X,1P,E18.10,4E15.5,E17.7,E15.5,2X,I4,A1)')RTMP(I),
	1                NEW_DI(I),NEW_ED(I),NEW_T(I),NEW_IRAT(I),
	1                NEW_VEL(I),NEW_CLUMP_FAC(I),I
	  WRITE(LU,'(F7.1)')1.0D0
	END DO
	CLOSE(UNIT=10)
!
	NEW_TAU=0.0D0
	IF(RD_MEANOPAC)THEN
	  CALL MON_INTERP(NEW_TAU,NEW_ND,IONE,RTMP,NEW_ND,OLD_TAU,ND,R,ND)
	  OPEN(UNIT=LU,FILE='R_GRID_CHK',STATUS='UNKNOWN')
	  WRITE(LU,'(A)')' '
	  WRITE(LU,'(A)')' Final grid computed with SET_SN_R_GRID '
          WRITE(LU,'(A)')' '
          WRITE(LU,'(A,17X,A,9X,A,8X,A,11X,A,10X,A,7X,A,6X,A,3X,A,13X,A,3X,A)')
	1           ' Depth','R','Ln(R)','dLn(R)','Tau','dTau','Ln(Tau)','dLn(Tau)','dTAU[I/I-1]','T','dT/T(I)'
	  T1=0.0D0
	  DO I=1,NEW_ND-1
	    IF(I .NE. 1)T1=(NEW_TAU(I+1)-NEW_TAU(I))/(NEW_TAU(I)-NEW_TAU(I-1))
	    WRITE(LU,'(I6,ES18.8,9ES14.4)')I,RTMP(I),LOG(RTMP(I)),LOG(RTMP(I+1)/RTMP(I)),
	1      NEW_TAU(I),NEW_TAU(I+1)-NEW_TAU(I),LOG(NEW_TAU(I)),LOG(NEW_TAU(I+1)/NEW_TAU(I)),
	1      T1,NEW_T(I),NEW_T(I+1)/NEW_T(I)-1.0D0
	  END DO
	  I=NEW_ND
	  WRITE(LU,'(I6,ES18.8,7ES14.4)')I,RTMP(I),LOG(RTMP(I)),0.0D0,NEW_TAU(I),0.0D0,LOG(NEW_TAU(I)),0.0D0
	  CLOSE(UNIT=LU)
	END IF
!
        CALL MON_INT_FUNS_V2(COEF,NEW_VEL,RTMP,NEW_ND)
        DO I=1,NEW_ND
          NEW_SIGMA(I)=COEF(I,3)
	  NEW_SIGMA(I)=RTMP(I)*NEW_SIGMA(I)/NEW_VEL(I)-1.0D0
	END DO
        OPEN(UNIT=10,FILE='RVSIG_HOPE',STATUS='UNKNOWN')
	  WRITE(10,'(A)')'!'
          WRITE(10,'(A,7X,A,9X,10X,A,11X,A,3X,A)')'!','R','V(km/s)','Sigma','Depth'
          WRITE(10,'(A)')'!'
          WRITE(10,'(A)')' '
          WRITE(10,'(I4,20X,A)')NEW_ND,'!Number of depth points`'
          WRITE(10,'(A)')' '
          DO I=1,ND
            WRITE(10,'(F18.8,ES17.7,F17.7,4X,I4)')RTMP(I),NEW_VEL(I),NEW_SIGMA(I),I
          END DO
	CLOSE(UNIT=10)

	STOP
	END
