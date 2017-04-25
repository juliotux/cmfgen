!
! Subroutine to create a grid, in TAU space, for solution of the transfer equation.
! Initial grid is created so that it satisfies:
!       R(I+1) > R_SCL_FAC*R(I)
!       V(I+1) > V_SCL_FAC*V(I)				
!       LoG(TAU(I+1)) < dLOG_TAU + LOG(TAU(I))
!
! This scaling will be modified slightly when the number of grid points is
! adjusted to the desired value.
!
	SUBROUTINE DET_R_GRID_V3(REV_TAU_GRID,NEW_ND,ND_MAX,TAU_MAX,
	1           dLOG_TAU,V_SCL_FAC,R_SCL_FAC,
	1           OUT_BND_OPT,OBND_PARS,NUM_OBND_PARAMS,
	1           R,V,TAU,ND)
	IMPLICIT NONE
!
! Altered 24-Jul-2015 : Changed to V3. R_SCL_FAC included in the CALL, and check on
!                         dR added. This check is particularly important for winds with small Vinf.
! Altered 20-May-2009 : Fundamental change to inclusion of extra boundary points.
!                         These are now directly included in the final tau grid,
!                         rather than the initial grid. This should give more control
!                         over the spacing.
! Altered 12-Feb-2007 : Modifications to allow more choice in the outer boundary 
!                         condition. dLOG_TAU and V_SCL_FAC now passed in call.
! Created 12-Aug-2006
!
	INTEGER NEW_ND				!Requested number of grid points in new grid.
	INTEGER ND_MAX				!Maximum number of grid points in grid.
	INTEGER NUM_OBND_PARAMS
	REAL*8 REV_TAU_GRID(NEW_ND)
	REAL*8 OBND_PARS(NUM_OBND_PARAMS)
	REAL*8 TAU_MAX				!Maximum optical depth for new grid.
	REAL*8 R_SCL_FAC			!Maxium values for R(I+1)/V(I)
	REAL*8 V_SCL_FAC			!Maxium values for V(I+1)/V(I)
	REAL*8 dLOG_TAU				!Maximu value for /\Tau
	CHARACTER(LEN=*) OUT_BND_OPT
!
! These describe the old grid.
!
	INTEGER ND
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 TAU(ND)
!
	INTEGER N_BND_PNTS		!Number of extra points used at boundaries
	INTEGER NS			!NEW_ND - N_BND_PNTS (initally)
!
! Local arrays.
!
	REAL*8 REV_R(ND_MAX)
	REAL*8 REV_V(ND_MAX)
	REAL*8 REV_TAU(ND_MAX)
	REAL*8 OLD_R(ND_MAX)
	REAL*8 OLD_TAU(ND_MAX)
	REAL*8 LOG_TAU(ND)
!
! Local variables.
!
	REAL*8 LOG_TAU_MAX
	REAL*8 dTAU
	REAL*8 T1,T2
	INTEGER ND_TMP
	INTEGER I,J,K,JST
	INTEGER LUER,ERROR_LU
	INTEGER, PARAMETER :: IONE=1
	EXTERNAL ERROR_LU
!
! Handle the grid near the oputer boundary, according to the option specified.
! For each option, we first check parameter validity. These options MUST BE
! consistent with the same options at the end of the routine. In this section
! we now onl determine how many points we are adding to the outer boundary.
!
	LUER=ERROR_LU()
	WRITE(6,*)'OUT_BND_OPT=',OUT_BND_OPT
	N_BND_PNTS=0
	IF(OUT_BND_OPT .EQ. 'POW')THEN
	  IF(NUM_OBND_PARAMS .NE. 2)THEN
	    WRITE(LUER,*)'Error in DET_R_GRID_V3 for OUT_BND_OPT=POW'
	    WRITE(LUER,*)'Exactly 2 parameters must be specified.'
	    STOP
	  ELSE IF(NINT(OBND_PARS(1)) .LT. 1 .OR. OBND_PARS(2) .GT. 10)THEN
	    WRITE(LUER,*)'Error in DET_R_GRID_V3 for OUT_BND_OPT=POW'
	    WRITE(LUER,*)'OBND_PARS(1) must be >=1 and <10'
	    STOP
	  ELSE IF(OBND_PARS(2) .LE. 1.0)THEN
	    WRITE(LUER,*)'Error in DET_R_GRID_V3 for OUT_BND_OPT=POW'
	    WRITE(LUER,*)'OBND_PARS(2) must be >1'
	    STOP
	  END IF
	  J=NINT(OBND_PARS(1))
	  N_BND_PNTS=J
	ELSE IF(OUT_BND_OPT .EQ. 'SPECIFY')THEN
	  DO K=1,NUM_OBND_PARAMS-1
	    IF(OBND_PARS(K) .LT. 1.0 .OR. OBND_PARS(K) .LE. OBND_PARS(K+1))THEN
	      WRITE(LUER,*)'Error in DET_R_GRID_V3'
	      WRITE(LUER,*)'OB parameters must be >1 and monotonically decreasing'
	      STOP
	    END IF
	  END DO
	  N_BND_PNTS=NUM_OBND_PARAMS
	ELSE IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	  N_BND_PNTS=2
	ELSE 
	  WRITE(LUER,*)'Unrecognized option in DET_R_GRID_V3'
	  STOP
	END IF
!
! Now define the uniform grid.
!
	LOG_TAU=LOG(TAU)
	LOG_TAU_MAX=LOG(TAU_MAX)
	REV_R(1)=R(1)
	REV_V(1)=V(1)
	REV_TAU(1)=LOG_TAU(1)
	I=1
!
	J=2
 	JST=J
	DO WHILE(REV_TAU(I) .LT. LOG_TAU_MAX)
!
	  I=I+1
	  IF(I+2 .GT. ND_MAX)THEN
	    WRITE(6,*)'Error in DET_R_GIRD_V1 --- ND_MAX too small'
	    WRITE(6,*)'ND_MAX=',ND_MAX
	    STOP
	  END IF
!
! This step ensure TAU spacing criterion is satisfied.
!
	  REV_TAU(I)=REV_TAU(I-1)+dLOG_TAU
!
! Check if we are at inner boundary. If so, we decrease spacing systematically
! to give more accuracy for our first orer boundary conditions. We ignore the 
! velocity check.
!
	  IF(REV_TAU(I) .GE. LOG_TAU_MAX-0.5*dLOG_TAU)THEN
	    T1=LOG_TAU_MAX-REV_TAU(I-1)
	    REV_TAU(I)=REV_TAU(I-1)+T1/2
	    I=I+1
	    REV_TAU(I)=LOG_TAU_MAX
	    N_BND_PNTS=N_BND_PNTS+2
	  ELSE IF(REV_TAU(I)+dLOG_TAU .GE. LOG_TAU_MAX)THEN
	    I=I+1
	    REV_TAU(I)=LOG_TAU_MAX
	    N_BND_PNTS=N_BND_PNTS+2
	  ELSE
!
! Check if V spacing criterion is satisifed. If not, shrink spacing. Linear
! interpolation if fine since we are only constructing the grid, not actual values
! on the grid.
!
! We test agaianst V(I-1) [rather than V(I)] because in some case step size can be 
!    very big. Changed Jul 31, 2011.
!
	    DO WHILE (REV_TAU(I) .GT. LOG_TAU(J))
	      IF(J .EQ. ND)EXIT
	      J=J+1
	    END DO
	    T1=(REV_TAU(I)-LOG_TAU(J-1))/(LOG_TAU(J)-LOG_TAU(J-1))
	    REV_V(I)=(1.0D0-T1)*V(J-1)+T1*V(J)
	    IF(REV_V(I-1) .GT. 0.1D0 .AND. REV_V(I) .LT. V_SCL_FAC*REV_V(I-1))THEN
	      J=JST
	      DO WHILE (V_SCL_FAC*REV_V(I-1) .LT. V(J))
	        IF(J .EQ. ND)EXIT
	        J=J+1
	       END DO
	      T1=(V_SCL_FAC*REV_V(I-1)-V(J-1))/(V(J)-V(J-1))
	      REV_TAU(I)=(1.0D0-T1)*LOG_TAU(J-1)+T1*LOG_TAU(J)
	      REV_V(I)=V_SCL_FAC*REV_V(I-1)
	    END IF
 	    J=JST
!
	    DO WHILE (REV_TAU(I) .GT. LOG_TAU(J))
	      IF(J .EQ. ND)EXIT
	      J=J+1
	    END DO
	    T1=(REV_TAU(I)-LOG_TAU(J-1))/(LOG_TAU(J)-LOG_TAU(J-1))
	    REV_R(I)=(1.0D0-T1)*R(J-1)+T1*R(J)
	    IF(REV_R(I-1) .GT. 0.1D0 .AND. REV_R(I) .LT. R_SCL_FAC*REV_R(I-1))THEN
	      J=JST
	      DO WHILE (R_SCL_FAC*REV_R(I-1) .LT. R(J))
	        IF(J .EQ. ND)EXIT
	        J=J+1
	       END DO
	      T1=(R_SCL_FAC*REV_R(I-1)-R(J-1))/(R(J)-R(J-1))
	      REV_TAU(I)=(1.0D0-T1)*LOG_TAU(J-1)+T1*LOG_TAU(J)
	      REV_R(I)=R_SCL_FAC*REV_R(I-1)
	    END IF
 	    J=JST
!
	  END IF
	  ND_TMP=I
	END DO
!
! The next few line can be uncommented if debuging.
!
	WRITE(80,'(/,A,/)')' Estimate R grid as determined by DET_R_GRID_V3'
	WRITE(80,'(A,4(9X,A))')'Index','   V',' Tau','dTau'
	T1=10.0D0; T1=LOG(10.0D0)
	DO I=1,ND_TMP-1
	  WRITE(80,'(I5,5ES14.4)')I,REV_V(I),REV_TAU(I)/T1,(REV_TAU(I+1)-REV_TAU(I))/T1
	END DO
	I=ND_TMP; WRITE(80,'(I5,4ES14.4)')I,REV_V(I),REV_TAU(I)/T1
!
! We now create a grid with the requested number of grid points. To do so
! we use the array index as the independent variable. This approach should
! preserve the basic grid spacing.
!
! The scaling for REV_R makes the X-limits identical for the new and old grids.
!
	OLD_TAU(1:ND_TMP)=REV_TAU(1:ND_TMP)
	DO I=1,ND_TMP; OLD_R(I)=I; END DO
	NS=NEW_ND-N_BND_PNTS
	DO I=1,NS
	  REV_R(I)=1+((I-1.0D0)*(ND_TMP-1.0D0) )/(NS-1.0D0)
	END DO
	CALL MON_INTERP(REV_TAU,NS,IONE,REV_R,NS,OLD_TAU,ND_TMP,OLD_R,ND_TMP)
!
! Now need to handle bondaries:
!
	dLOG_TAU=REV_TAU(2)-REV_TAU(1)
	IF(OUT_BND_OPT .EQ. 'POW')THEN
	  J=NINT(OBND_PARS(1))
	  REV_TAU(J+2:NS+J)=REV_TAU(2:NS)
	  DO K=2,J+1
	    REV_TAU(K)=REV_TAU(1)+dLOG_TAU/(OBND_PARS(2)**(J+2-K))
	  END DO
	  NS=NS+J
	ELSE IF(OUT_BND_OPT .EQ. 'SPECIFY')THEN
	  J=NUM_OBND_PARAMS
	  REV_TAU(J+2:NS+J)=REV_TAU(2:NS)
	  DO K=1,NUM_OBND_PARAMS
	    REV_TAU(K+1)=REV_TAU(1)+dLOG_TAU/OBND_PARS(K)
	  END DO
	  NS=NS+J
	ELSE IF(OUT_BND_OPT .EQ. 'DEFAULT')THEN
	  J=2
	  REV_TAU(J+2:NS+J)=REV_TAU(2:NS)
	  REV_TAU(2)=REV_TAU(1)+dLOG_TAU/9.0D0
	  REV_TAU(3)=REV_TAU(1)+dLOG_TAU/3.0D0
	  NS=NS+J
	END IF
!
! We now need to do the inner boundary.
!
	REV_TAU(NS+2)=REV_TAU(NS)
	dLOG_TAU=REV_TAU(NS)-REV_TAU(NS-1)
	REV_TAU(NS)=REV_TAU(NS-1)+0.6D0*dLOG_TAU
	REV_TAU(NS+1)=REV_TAU(NS-1)+0.9D0*dLOG_TAU
	NS=NS+2
	IF(NS .NE. NEW_ND)THEN
	  WRITE(6,*)'Error in DET_R_GRID_V3: Inconsistent NS and NEW_ND'
	  WRITE(6,*)'        NS=',NS
	  WRITE(6,*)'    NEW_ND=',NEW_ND
	  WRITE(6,*)'N_BND_PNTS=',N_BND_PNTS
	  STOP
	END IF
!
	WRITE(80,'(/,/,1X,A)')'Final Tau Grid'
	WRITE(80,'(A,4(9X,A))')'Index',' Tau','dTau'
	T1=10.0D0; T1=LOG(10.0D0)
	DO I=1,ND_TMP-1
	  WRITE(80,'(I5,5ES14.4)')I,REV_TAU(I)/T1,(REV_TAU(I+1)-REV_TAU(I))/T1
	END DO
	I=ND_TMP; WRITE(80,'(I5,4ES14.4)')I,REV_TAU(I)/T1
!
! We finish by remembering that we have been working with LOG(TAU).
!
	REV_TAU_GRID(1:NEW_ND)=EXP(REV_TAU(1:NEW_ND))
!
	RETURN
	END
