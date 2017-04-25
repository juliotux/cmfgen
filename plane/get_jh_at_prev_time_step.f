        MODULE MOD_GET_JH
!
! We assume that J, H are actually RSQJ, and RSQH.
!
	REAL*8, ALLOCATABLE :: OLD_R(:)
	REAL*8, ALLOCATABLE :: OLD_V(:)
	REAL*8, ALLOCATABLE :: OLD_J(:)
	REAL*8, ALLOCATABLE :: OLD_H(:)
!
	REAL*8, ALLOCATABLE :: NUST(:)
	REAL*8, ALLOCATABLE :: JST(:,:)
	REAL*8, ALLOCATABLE :: HST(:,:)
	REAL*8, ALLOCATABLE :: H_INBC_ST(:)
	REAL*8, ALLOCATABLE :: H_OUTBC_ST(:)
!
	REAL*8, ALLOCATABLE :: NEW_J_SAV(:)
	REAL*8, ALLOCATABLE :: NEW_H_SAV(:)
	REAL*8, ALLOCATABLE :: LOG_V(:)
	REAL*8, ALLOCATABLE :: LOG_MIDV(:)
	REAL*8, ALLOCATABLE :: WRK_NEW_H(:)
!
	REAL*8, ALLOCATABLE :: LOG_OLD_V(:)
	REAL*8, ALLOCATABLE :: LOG_OLD_MIDV(:)
	REAL*8, ALLOCATABLE :: WRK_OLD_H(:)
!
	REAL*8 H_INBC_OLDT_SAV
	REAL*8 H_OUTBC_OLDT_SAV
	REAL*8, SAVE :: NU_SAV=0.0D0
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: NSM=100
!
! INDX is used to indicate the current location in the frequency grid.
!
	INTEGER INDX
	INTEGER COUNTER		!
	INTEGER IREC		!Record to be read next
	INTEGER ST_IREC		!Record data begins at.
	INTEGER NCF_OLD 
	INTEGER ND_OLD 
	INTEGER ND_OLD_P1
	INTEGER ND_P1
	INTEGER LU_ER
	INTEGER LU_WARN
	INTEGER, SAVE :: LU_IN=0
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
	END MODULE MOD_GET_JH
! 
	SUBROUTINE GET_JH_AT_PREV_TIME_STEP(NEW_J,NEW_H,
	1            H_INBC_OLDT,H_OUTBC_OLDT,DELTA_TIME_SECS,
	1            NU,R,V,ND,INIT,OPTION)
	USE MOD_GET_JH
	IMPLICIT NONE
!
! Altered 14-Jan -2010 : Code stops when NU too small, and fixed array access violation.
! Altered 26-July-2008 : Allow NU > NU(old).
! Altered 17-July-2008 : Fixed DELTA_TIME_SECS which was wrong when ND and ND_OLD
!                             were unequal.
!
	INTEGER ND
	REAL*8 NU
	REAL*8 DELTA_TIME_SECS
	REAL*8 R(ND)
	REAL*8 V(ND)
        REAL*8 NEW_J(ND)
	REAL*8 NEW_H(ND-1)
	REAL*8 H_INBC_OLDT
	REAL*8 H_OUTBC_OLDT
!
! Used to indicate that we are paasing the first frequency.
!
	LOGICAL INIT
	CHARACTER(LEN=*), OPTIONAL :: OPTION
!
! Lcoal variables:
!
	REAL*8  T1
	INTEGER I		!Used as depth index
	INTEGER ML		!Used as frequency index
	INTEGER IOS		!I/O error identifier
	INTEGER REC_LENGTH 
	INTEGER ERROR_LU,WARNING_LU
	LOGICAL EQUAL
	EXTERNAL EQUAL,ERROR_LU,WARNING_LU
!
	CHARACTER*30 FILE_DATE
	LOGICAL FILE_OPEN
!
	IF(FIRST_TIME)THEN
	  LU_ER=ERROR_LU()
	  LU_WARN=WARNING_LU()
	  CALL GET_LU(LU_IN,'GET_JH_AT_PREV_TIME_STEP')
	  WRITE(LU_WARN,'(/,1X,A,I5)')'Logical unit for JH input is',LU_IN
          CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,'JH_AT_OLD_TIME',LU_IN,IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LU_ER,*)'Error opening/reading JH_AT_OLD_TIME_INFO file: check format'
            STOP
	  END IF
          OPEN(UNIT=LU_IN,FILE='JH_AT_OLD_TIME',STATUS='OLD',ACTION='READ',
	1        RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LU_ER,*)'Error opening JH_AT_OLD_TIME'
            WRITE(LU_ER,*)'IOS=',IOS
            STOP
	  END IF
!
	  READ(LU_IN,REC=3)ST_IREC,NCF_OLD,ND_OLD
!
	  ALLOCATE (OLD_V(ND_OLD))
	  ALLOCATE (OLD_R(ND_OLD))
!
	  IREC=ST_IREC
	  READ(LU_IN,REC=IREC)OLD_R,OLD_V;    IREC=IREC+1
!
	  ALLOCATE (JST(ND_OLD,NSM))
	  ALLOCATE (HST(ND_OLD,NSM))
	  ALLOCATE (NUST(NSM));        NUST=0.0D0
	  ALLOCATE (H_INBC_ST(NSM))
	  ALLOCATE (H_OUTBC_ST(NSM))
!
	  ALLOCATE (OLD_J(ND_OLD))
	  ALLOCATE (OLD_H(ND_OLD-1))
!
	  ALLOCATE (LOG_OLD_V(ND_OLD))
	  ALLOCATE (LOG_OLD_MIDV(ND_OLD+1))
	  ALLOCATE (WRK_OLD_H(ND_OLD+1))
!
	  ALLOCATE (LOG_V(ND))
	  ALLOCATE (LOG_MIDV(ND+1))
	  ALLOCATE (WRK_NEW_H(ND+1))
!
	  ALLOCATE (NEW_J_SAV(ND))
	  ALLOCATE (NEW_H_SAV(ND-1))
!
! NB: V and OLD_V are defined on the grid.
!     MIDV & OLD_MIDV are defined at the mid-grid points. Because the mid grid points
!        may not match exactly, and because the outer biundaries may not match, we also
!        add the two boundaries to the MIDV vectors.
!
	  LOG_OLD_V=LOG(OLD_V)
	  DO I=1,ND_OLD-1
	    LOG_OLD_MIDV(I+1)=LOG( 0.5D0*(OLD_V(I)+OLD_V(I+1)) )
	  END DO
	  ND_OLD_P1=ND_OLD+1
	  LOG_OLD_MIDV(1)=LOG_OLD_V(1)
	  LOG_OLD_MIDV(ND_OLD_P1)=LOG_OLD_V(ND_OLD)
	  FIRST_TIME=.FALSE.
!
	END IF
!
! If we are using an adaptive grid, the Velocity grid will change from iteration to iteration.
! Thus we have to reset the new grid when INIT is TRUE.
!
	IF(INIT)THEN
	  ND_P1=ND+1
	  LOG_V=LOG(V)
	  DO I=1,ND-1
	    LOG_MIDV(I+1)=LOG( 0.5D0*(V(I)+V(I+1)) )
	  END DO
	  LOG_MIDV(1)=LOG_V(1)
	  LOG_MIDV(ND_P1)=LOG_V(ND)
!
	  T1=1.0D-07					!Altered 3-Feb-2008 
	  IF(EQUAL(V(ND),OLD_V(ND_OLD),T1))THEN
	    OLD_V(ND_OLD)=V(ND)
	    LOG_OLD_V(ND_OLD)=LOG_V(ND)
	    LOG_OLD_MIDV(ND_OLD_P1)=LOG_MIDV(ND_P1)
	  ELSE
	    WRITE(LU_ER,*)'Error in GET_JH_AT_PREV_TIME_STEP'
	    WRITE(LU_ER,*)'Velocities at inner boundary are not identical'
	    WRITE(LU_ER,*)'Old velocity at depth is',OLD_V(ND_OLD)
	    WRITE(LU_ER,*)'New velocity at depth is',V(ND)
	    WRITE(LU_ER,*)'Make sure you also updated JH_AT_OLD_TIME_INFO when you stated the model'
	    WRITE(LU_ER,*)'This will cause an error if ND has changed between successive models'
	    STOP
	  END IF
!
	  T1=1.0D-07
	  IF(EQUAL(V(1),OLD_V(1),T1))THEN
	     OLD_V(1)=V(1)
	     LOG_OLD_V(1)=LOG_V(1)
	     LOG_OLD_MIDV(1)=LOG_MIDV(1)
	  ELSE IF(V(1) .GT. OLD_V(1))THEN
	    WRITE(LU_ER,*)'Error in GET_JH_AT_PREV_TIME_STEP'
	    WRITE(LU_ER,*)'Velocity at outer boundary is too large'
	    WRITE(LU_ER,*)'Old velocity at outer boundary is',OLD_V(1)
	    WRITE(LU_ER,*)'New velocity at outer boundary is',V(1)
	    STOP
	  END IF
	END IF
!
! Compute elapsed time between the two models. The factor of 10^5 arises because
! R is in units of 10^10 cm, and V is in units of kms/s (i.e., 10^5 cm/s).
!
	DELTA_TIME_SECS=1.0D+05*(R(ND)-OLD_R(ND_OLD))/V(ND)
	IF(DELTA_TIME_SECS .LE. 0.0D0)THEN
	  WRITE(LU_ER,*)'Error in GET_JH_AT_PREV_TIME_STEP'
	  WRITE(LU_ER,*)'Running a time dependent model with 0 ov -ve time step'
	  WRITE(LU_ER,*)'R(ND),OLD_R(ND_OLD):',R(ND),OLD_R(ND_OLD)
	  STOP
	END IF
!
	IF(PRESENT(OPTION))THEN
	  IF(OPTION .EQ. 'GREY')THEN
	    OLD_J(:)=0.0D0; OLD_H(:)=0.0; H_INBC_OLDT=0.0D0; H_OUTBC_OLDT=0.0D0
	    READ(LU_IN,REC=ST_IREC+2)JST(:,1),HST(1:ND_OLD-1,1),H_INBC_ST(1),H_OUTBC_ST(1),NUST(1)
	    DO ML=2,NCF_OLD
	       READ(LU_IN,REC=ST_IREC+1+ML)JST(:,2),HST(1:ND_OLD-1,2),H_INBC_ST(2),H_OUTBC_ST(2),NUST(2)
	       OLD_J(:)=OLD_J(:)+(NUST(1)-NUST(2))*(JST(:,1)+JST(:,2))
	       OLD_H(:)=OLD_H(:)+(NUST(1)-NUST(2))*(HST(1:ND_OLD-1,1)+HST(1:ND_OLD-1,2))
	       H_INBC_OLDT=H_INBC_OLDT+(NUST(1)-NUST(2))*(H_INBC_ST(1)+H_INBC_ST(2))
	       H_OUTBC_OLDT=H_OUTBC_OLDT+(NUST(1)-NUST(2))*
	1                       (H_OUTBC_ST(1)*JST(1,1)+H_OUTBC_ST(2)*JST(1,2))
	       JST(:,1)=JST(:,2); HST(:,1)=HST(:,2);  NUST(1)=NUST(2)
	       H_INBC_ST(2)=H_INBC_ST(1); H_OUTBC_ST(2)=H_OUTBC_ST(1)
	    END DO
	    H_OUTBC_OLDT=H_OUTBC_OLDT/OLD_J(1)
	    H_INBC_OLDT=0.5D+15*H_INBC_OLDT
	    H_INBC_OLDT=2*H_INBC_OLDT                     !fudge
	    OLD_J=0.5D+15*OLD_J; OLD_H=0.5D+15*OLD_H;
!
! Now need to interpolate in R space.
!
	    CALL MON_INTERP(NEW_J,ND,IONE,LOG_V,ND,OLD_J,ND_OLD,LOG_OLD_V,ND_OLD)
	    WRK_OLD_H(2:ND_OLD)=OLD_H(1:ND_OLD-1)
	    WRK_OLD_H(1)=H_OUTBC_OLDT*OLD_J(1)
	    WRK_OLD_H(ND_OLD_P1)=H_INBC_OLDT
	    CALL MON_INTERP(WRK_NEW_H,ND_P1,IONE,LOG_MIDV,ND_P1,WRK_OLD_H,ND_OLD_P1,LOG_OLD_MIDV,ND_OLD_P1)
	    NEW_H(1:ND-1)=WRK_NEW_H(2:ND)
	    T1=1.0D-12; IF( .NOT. EQUAL(V(1),OLD_V(1),T1) )H_OUTBC_OLDT=WRK_NEW_H(1)/NEW_J(1)
	    RETURN
	  ELSE IF(OPTION .EQ. 'OGREY')THEN
	    READ(LU_IN,REC=ST_IREC+1)(OLD_J(I),I=1,ND_OLD),(OLD_H(I),I=1,ND_OLD),
	1               H_INBC_OLDT,H_OUTBC_OLDT
!
! Now need to interpolate in R space.
!
	    CALL MON_INTERP(NEW_J,ND,IONE,LOG_V,ND,OLD_J,ND_OLD,LOG_OLD_V,ND_OLD)
	    WRK_OLD_H(2:ND_OLD)=OLD_H(1:ND_OLD-1)
	    WRK_OLD_H(1)=H_OUTBC_OLDT*OLD_J(1)
	    WRK_OLD_H(ND_OLD_P1)=H_INBC_OLDT
	    CALL MON_INTERP(WRK_NEW_H,ND_P1,IONE,LOG_MIDV,ND_P1,WRK_OLD_H,ND_OLD_P1,LOG_OLD_MIDV,ND_OLD_P1)
	    NEW_H(1:ND-1)=WRK_NEW_H(2:ND)
	    T1=1.0D-12; IF( .NOT. EQUAL(V(1),OLD_V(1),T1) )H_OUTBC_OLDT=WRK_NEW_H(1)/NEW_J(1)
	    RETURN
	  ELSE IF(OPTION .EQ. 'NORMAL')THEN
	  ELSE
	    WRITE(LU_ER,*)'Error in GET_JH_AT_PREV_TIME_STEP'
	    WRITE(LU_ER,*)'Invalid option passed'
	  END IF
	END IF
!
! Check to see if we have already computed the desired values.
!
	IF(NU_SAV .EQ. NU)THEN
	  NEW_J=NEW_J_SAV
	  NEW_H=NEW_H_SAV
	  H_INBC_OLDT=H_INBC_OLDT_SAV
	  H_OUTBC_OLDT=H_OUTBC_OLDT_SAV
	  RETURN
	END IF
!
! Initialize indices.
!
	IF(INIT)THEN
	  INDX=2
	  COUNTER=0
	  NUST=1.0D+10				!A big number
!
! IREC will now point at the first record  contaiing the frequency dependent J & H.
! We skip 2 records, because of the R,V & grey J & H output.
!
	  IREC=ST_IREC+2
	  WRITE(LU_WARN,*)'Start record and IREC are ',ST_IREC,IREC
	END IF
!
! We read the data in blocks of NSM for reasons of efficiency.
! The first IF loop ensures that the storage arrays always contain NSM
! valid elements. As we are updating the store we need to reset INDX to 2.
!
	DO WHILE(NU .LE. NUST(NSM) .AND. COUNTER .LT. NCF_OLD)
	  IF(COUNTER+NSM .GT. NCF_OLD)THEN
	    IREC=IREC-(COUNTER+NSM-NCF_OLD)
	    COUNTER=COUNTER-(COUNTER+NSM-NCF_OLD)
	  END IF
!
	  DO ML=1,NSM
	    COUNTER=COUNTER+1
	    READ(LU_IN,REC=IREC)(JST(I,ML),I=1,ND_OLD),(HST(I,ML),I=1,ND_OLD-1),
	1                          H_INBC_ST(ML),H_OUTBC_ST(ML),NUST(ML)
	    IREC=IREC+1
!	    WRITE(250,'(7ES14.5)')NUST(ML),JST(1,ML),HST(1,ML),JST(ND_OLD,ML),HST(ND_OLD-1,ML),
!	1                          H_INBC_ST(ML),H_OUTBC_ST(ML)
	  END DO
	  INDX=2
!
	END DO
!
	IF(COUNTER .EQ. NCF_OLD .AND. NU .LT. NUST(NSM))THEN
	  WRITE(LU_ER,*)'Warning in GET_JH_AT_PREV_TIME_STEP'
	  WRITE(LU_ER,*)'Invalid minmum frequency --- outside range'
	  WRITE(LU_ER,*)'NU=',NU
	  WRITE(LU_ER,*)'NUST=',NUST(NSM)
	END IF
!
	IF(INIT .AND. NU .GT. NUST(1))THEN
	  WRITE(LU_ER,*)'Warning in GET_JH_AT_PREV_TIME_STEP'
	  WRITE(LU_ER,*)'Invalid maximum frequency --- outside range'
	  WRITE(LU_ER,*)'Seeting flux at value for maximum frequency.'
	  WRITE(LU_ER,'(A,ES16.6,A,ES16.6)')'NU=',NU,'NUST=',NUST(1)
	END IF
!
	IF(NU .GT. NUST(1))THEN
	  DO I=1,ND_OLD
	    OLD_J(I)=JST(I,1)
	  END DO
	  DO I=1,ND_OLD-1
	    OLD_H(I)=HST(I,1)
	  END DO
	  H_INBC_OLDT=H_INBC_ST(1)
	  H_OUTBC_OLDT=H_OUTBC_ST(1)
	ELSE IF(NU .LT. NUST(NSM))THEN
	  DO I=1,ND_OLD
	    OLD_J(I)=JST(I,NSM)
	  END DO
	  DO I=1,ND_OLD-1
	    OLD_H(I)=HST(I,NSM)
	  END DO
	  H_INBC_OLDT=H_INBC_ST(NSM)
	  H_OUTBC_OLDT=H_OUTBC_ST(NSM)
	ELSE
!
! We initially use linear interpolation in frequency.
!
	  DO WHILE(NU .LT. NUST(INDX))
	    INDX=INDX+1
	  END DO
	  T1=(NU-NUST(INDX))/(NUST(INDX-1)-NUST(INDX))
	  DO I=1,ND_OLD
	    OLD_J(I)=(1.0D0-T1)*JST(I,INDX) + T1*JST(I,INDX-1)
	  END DO
	  DO I=1,ND_OLD-1
	    OLD_H(I)=(1.0D0-T1)*HST(I,INDX) + T1*HST(I,INDX-1)
	  END DO
	  H_INBC_OLDT=(1.0D0-T1)*H_INBC_ST(INDX) + T1*H_INBC_ST(INDX-1)
	  H_OUTBC_OLDT=(1.0D0-T1)*H_OUTBC_ST(INDX) + T1*H_OUTBC_ST(INDX-1)
	END IF
!
! Now need to interpolate in R space.
!
        CALL MON_INTERP(NEW_J,ND,IONE,LOG_V,ND,OLD_J,ND_OLD,LOG_OLD_V,ND_OLD)
!
	WRK_OLD_H(2:ND_OLD)=OLD_H(1:ND_OLD-1)
	WRK_OLD_H(1)=H_OUTBC_OLDT*OLD_J(1)
	WRK_OLD_H(ND_OLD_P1)=H_INBC_OLDT
        CALL MON_INTERP(WRK_NEW_H,ND_P1,IONE,LOG_MIDV,ND_P1,WRK_OLD_H,ND_OLD_P1,LOG_OLD_MIDV,ND_OLD_P1)
	NEW_H(1:ND-1)=WRK_NEW_H(2:ND)
	T1=1.0D-12; IF(.NOT. EQUAL(V(1),OLD_V(1),T1) )H_OUTBC_OLDT=WRK_NEW_H(1)/NEW_J(1)
!
!	WRITE(251,'(7ES14.5)')NU,NEW_J(1),NEW_H(1),NEW_J(ND),NEW_H(ND-1),H_INBC_OLDT,H_OUTBC_OLDT
!
! As we enter this routine for the same frequency (iteration and variation) several
! times, save interpolated values.
!
	NU_SAV=NU
	NEW_J_SAV=NEW_J
	NEW_H_SAV=NEW_H
	H_INBC_OLDT_SAV=H_INBC_OLDT
	H_OUTBC_OLDT_SAV=H_OUTBC_OLDT
!
	RETURN
	END
