	PROGRAM PLT_J_COMP
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 10-June-2009 : Inserted JMF, JMR options
!                         Made J_COMP_full 2nd default filename
!
	INTEGER, PARAMETER :: NCF_MAX=600000
!
	REAL*8 INDX(NCF_MAX)
	REAL*8 NU(NCF_MAX)
	REAL*8 LAM(NCF_MAX)
!
	REAL*8 JM_OUT(NCF_MAX)
	REAL*8 JR_OUT(NCF_MAX)
	REAL*8 DIFF_OUT(NCF_MAX)
	REAL*8 HBC_CMF(NCF_MAX)
!
	REAL*8 JM_IN(NCF_MAX)
	REAL*8 JR_IN(NCF_MAX)
	REAL*8 DIFF_IN(NCF_MAX)
!
	REAL*8 YV(NCF_MAX)
!
	CHARACTER*80 FILENAME
	CHARACTER*80 STRING
        CHARACTER*30 UC
        EXTERNAL UC
!
	REAL*8 SUM1
	REAL*8 SUM2
	REAL*8 SUM3
	CHARACTER(LEN=20) PLT_OPT
	CHARACTER(LEN=20) YLAB
!
	INTEGER I,IBEG
	INTEGER NCF
	INTEGER IOS
	INTEGER LUM_STAR
	LOGICAL LOGY
!
	LOGY=.FALSE.
 	FILENAME='J_COMP'
1000	CONTINUE
	CALL GEN_IN(FILENAME,'File with data to be plotted')
	IF(FILENAME .EQ. ' ')GOTO 1000
	FILENAME=ADJUSTL(FILENAME)
	IF(UC(TRIM(FILENAME)) .EQ. 'EXIT')THEN
	  WRITE(6,*)'Exiting code'
	  STOP
	END IF
!
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open file: ',TRIM(FILENAME)
 	    FILENAME='J_COMP_full'
	    GOTO 1000
	  END IF
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'HBC_CMF') .EQ. 0)
	    READ(11,'(A)')STRING
	  END DO
!
	  I=0
	  DO WHILE(1. EQ. 1)
	    I=I+1
	    IF(I .GT. NCF_MAX)THEN
	      WRITE(6,*)'************************************************'
	      WRITE(6,*)' '
	      WRITE(6,*)'Error --- insufficent memory to read in all data'
	      WRITE(6,*)'Edit NCF_MAX in plt_j_comp.f to get more memory'
	      WRITE(6,*)' '
	      WRITE(6,*)'************************************************'
	      EXIT
	    END IF
	    READ(11,*,END=100)INDX(I),NU(I), JM_OUT(I),JR_OUT(I),DIFF_OUT(I),HBC_CMF(I),
	1                                      JM_IN(I),JR_IN(I),DIFF_IN(I)
	  END DO
!
100	  CONTINUE
	  NCF=I-1
	CLOSE(UNIT=11)
	LAM(1:NCF)=2.99792458D+03/NU(1:NCF)
!
	SUM1=0.0D0
	SUM2=0.0D0
	SUM3=0.0D0
	SUM1=SUM(DIFF_OUT(1:NCF))/NCF
	SUM3=SUM(ABS(DIFF_OUT(1:NCF)))/NCF
	SUM2=SUM(DIFF_OUT(1:NCF)*DIFF_OUT(1:NCF))
	SUM2=SQRT(SUM2/NCF)
!
!
	WRITE(6,*)' '
	WRITE(6,'(A,F6.2,A)')' Mean difference is:          ',SUM1,'%'
	WRITE(6,'(A,F6.2,A)')' Mean absolute difference is: ',SUM3,'%'
	WRITE(6,'(A,F6.2,A)')' Root mean square is:         ',SUM2,'%'
!
	PLT_OPT='DW'
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  WRITE(6,*)'Plot options are:'
	  WRITE(6,*)'  JF: Plot J ray solution against frequency at outer boundary'
	  WRITE(6,*)' JMF: Plot J mom solution against frequency at outer boundary'
	  WRITE(6,*)'  JW: Plot J ray solution against wavelength at outer boundary'
	  WRITE(6,*)' JMW: Plot J mom solution against wavelength at outer boundary'
	  WRITE(6,*)' dJF: Plot J(ray)-J(mom) against frequency at outer boundary'
	  WRITE(6,*)' dJW: Plot J(ray)-J(mom) against wavelength at outer boundary'
	  WRITE(6,*)'  DF: Plot %error against frequency at outer boundary'
	  WRITE(6,*)'  DW: Plot $error against wavelength at outer boundary'
	  WRITE(6,*)' I??: Inner boundary options (as for outrer boundary)'
	  WRITE(6,*)'E(X): Exit routine'
	  CALL GEN_IN(PLT_OPT,'Plot option')
!
	  IF(UC(PLT_OPT) .EQ. 'LY')THEN
	     LOGY=.NOT. LOGY
	     IF(LOGY)THEN
	       WRITE(6,*)'Y axis is now logarithmic'
	     ELSE
	       WRITE(6,*)'Y axis is now linear'
	     END IF
	  ELSE IF(UC(PLT_OPT) .EQ. 'JF')THEN
	    YLAB='J(ray)'
	    IF(LOGY)YLAB='Log J(ray)'
	    CALL SET_YVEC(YV,JR_OUT,NCF,LOGY)
	    CALL DP_CURVE(NCF,NU,YV)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)',YLAB,' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'JW')THEN
	    YLAB='J(ray)'
	    IF(LOGY)YLAB='Log J(ray)'
	    CALL SET_YVEC(YV,JR_OUT,NCF,LOGY)
	    CALL DP_CURVE(NCF,LAM,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)',YLAB,' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'JMF')THEN
	    YLAB='J(mom)'
	    IF(LOGY)YLAB='Log J(mom)'
	    CALL SET_YVEC(YV,JM_OUT,NCF,LOGY)
	    CALL DP_CURVE(NCF,NU,YV)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)',YLAB,' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'JMW')THEN
	    YLAB='J(mom)'
	    IF(LOGY)YLAB='Log J(mom)'
	    CALL SET_YVEC(YV,JM_OUT,NCF,LOGY)
	    CALL DP_CURVE(NCF,LAM,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)',YLAB,' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'DJW')THEN
	    YV(1:NCF)=JR_OUT(1:NCF)-JM_OUT(1:NCF)
	    CALL DP_CURVE(NCF,LAM,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)','J(ray)-J(mom)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'DJF')THEN
	    YV(1:NCF)=JR_OUT(1:NCF)-JM_OUT(1:NCF)
	    CALL DP_CURVE(NCF,NU,YV)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','J(ray)-J(mom)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'DF')THEN
	    CALL DP_CURVE(NCF,NU,DIFF_OUT)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','% Difference',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'DW')THEN
	    CALL DP_CURVE(NCF,LAM,DIFF_OUT)
	    CALL GRAMON_PGPLOT('\gl(\A)','% Difference',' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'IJF')THEN
	    YLAB='J(ray)'
	    IF(LOGY)YLAB='Log J(ray)'
	    CALL SET_YVEC(YV,JR_IN,NCF,LOGY)
	    CALL DP_CURVE(NCF,NU,YV)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)',YLAB,' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'IJW')THEN
	    YLAB='J(ray)'
	    IF(LOGY)YLAB='Log J(ray)'
	    CALL SET_YVEC(YV,JR_IN,NCF,LOGY)
	    CALL DP_CURVE(NCF,LAM,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)',YLAB,' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'IJMF')THEN
	    YLAB='J(mom)'
	    IF(LOGY)YLAB='Log J(mom)'
	    CALL SET_YVEC(YV,JM_IN,NCF,LOGY)
	    CALL DP_CURVE(NCF,NU,YV)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)',YLAB,' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'IJMW')THEN
	    YLAB='J(mom)'
	    IF(LOGY)YLAB='Log J(mom)'
	    CALL SET_YVEC(YV,JM_IN,NCF,LOGY)
	    CALL DP_CURVE(NCF,LAM,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)',YLAB,' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'IDJW')THEN
	    YV(1:NCF)=JR_IN(1:NCF)-JM_IN(1:NCF)
	    CALL DP_CURVE(NCF,LAM,YV)
	    CALL GRAMON_PGPLOT('\gl(\A)','J(ray)-J(mom)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'IDJF')THEN
	    YV(1:NCF)=JR_IN(1:NCF)-JM_IN(1:NCF)
	    CALL DP_CURVE(NCF,NU,YV)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','J(ray)-J(mom)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'IDF')THEN
	    CALL DP_CURVE(NCF,NU,DIFF_IN)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','% Difference',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'IDW')THEN
	    CALL DP_CURVE(NCF,LAM,DIFF_IN)
	    CALL GRAMON_PGPLOT('\gl(\A)','% Difference',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
	END DO
!
	END
!
	SUBROUTINE SET_YVEC(YV,RJ,NCF,LOGY)
	IMPLICIT NONE
!
	INTEGER NCF
	REAL*8 YV(NCF)
	REAL*8 RJ(NCF)
	LOGICAL LOGY
	INTEGER I
!
	DO I=1,NCF
	  IF(RJ(I) .GT. 0.0D0)THEN
	    YV(I)=LOG10(RJ(I))
	  ELSE
	    YV(I)=-1000.0D0
	  END IF
	END DO
!
	RETURN
	END
