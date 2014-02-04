!
! Subroutine to create a grid, in TAU space, for solution of the transfer equation.
! Initial grid is created so that it satisfies:
!       V(I+1) > 0.67 V(I)				(set by V_SCL_FAC)
!       LoG(TAU(I+1)) < 0.25 + LOG(TAU(I))		(set by dLOG_TAU)
!
! This scaling will be modified slightly when the number of grid points is
! adjusted to the desired value.
!
	SUBROUTINE DET_R_GRID_V1(REV_TAU_GRID,NEW_ND,ND_MAX,TAU_MAX,DO_OUT_BOUND,R,V,TAU,ND)
	IMPLICIT NONE
!
! Created 12-Aug-2006
!
	INTEGER NEW_ND		!Requested number of grid points in new grid.
	INTEGER ND_MAX		!Maximum number of grid points in grid.
	REAL*8 REV_TAU_GRID(NEW_ND)
	REAL*8 TAU_MAX		!Maximum optical depth for new grid.
	LOGICAL DO_OUT_BOUND
!
! These describe the old grid.
!
	INTEGER ND
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 TAU(ND)
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
	REAL*8, PARAMETER :: V_SCL_FAC=0.67D0
	REAL*8, PARAMETER :: dLOG_TAU=0.25D0

	REAL*8 LOG_TAU_MAX
	REAL*8 dTAU
	REAL*8 T1,T2
	INTEGER ND_TMP
	INTEGER I,J,JST
	INTEGER, PARAMETER :: IONE=1
!
	LOG_TAU=LOG(TAU)
	LOG_TAU_MAX=LOG(TAU_MAX)
!
	I=1
	REV_R(1)=R(1)
	REV_V(1)=V(1)
	REV_TAU(1)=LOG_TAU(1)
!
	IF(DO_OUT_BOUND)THEN
	  REV_TAU(2)=REV_TAU(1)+dLOG_TAU/9.0D0
	  REV_TAU(3)=REV_TAU(1)+dLOG_TAU/3.0D0
	  REV_TAU(4)=REV_TAU(1)+dLOG_TAU
!
!Estimate V: accuracy at outer boundary not important.	
!
	  J=2
	  DO WHILE(REV_TAU(4) .GT. LOG_TAU(J))
	    J=J+1
	  END DO
	  I=2
	  T1=(REV_TAU(I)-LOG_TAU(1))/(LOG_TAU(J)-LOG_TAU(1))
	  REV_V(I)=(1.0D0-T1)*V(1)+T1*V(J)
	  I=3
	  T1=(REV_TAU(I)-LOG_TAU(1))/(LOG_TAU(J)-LOG_TAU(1))
	  REV_V(I)=(1.0D0-T1)*V(1)+T1*V(J)
	  I=4
	  T1=(REV_TAU(I)-LOG_TAU(1))/(LOG_TAU(J)-LOG_TAU(1))
	  REV_V(I)=(1.0D0-T1)*V(1)+T1*V(J)
	END IF
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
! This step ensure TAU spacing criterin is satisfied.
!
	  REV_TAU(I)=REV_TAU(I-1)+dLOG_TAU
!
! Check if we are at inner boundary. If so, we decrease spacing systematically
! to give more accuracy for our first orer boundary conditions. We ignore the 
! velocity check.
!
	  IF(REV_TAU(I)+dLOG_TAU .GE. LOG_TAU_MAX)THEN
	    T1=LOG_TAU_MAX-REV_TAU(I-1)
	    REV_TAU(I)=REV_TAU(I-1)+0.6*T1
	    I=I+1
	    REV_TAU(I)=REV_TAU(I-2)+0.9*T1
	    I=I+1
	    REV_TAU(I)=LOG_TAU_MAX
	  ELSE
!
! Check if V spacing criterion is satisifed. If not, shrink spacing. Linear
! interpolation if fine since we are only constructing the grid, not actual values
! on the grid.
!
	    DO WHILE (REV_TAU(I) .GT. LOG_TAU(J))
	      IF(J .EQ. ND)EXIT
	      J=J+1
	    END DO
	    T1=(REV_TAU(I)-LOG_TAU(J-1))/(LOG_TAU(J)-LOG_TAU(J-1))
	    REV_V(I)=(1.0D0-T1)*V(J-1)+T1*V(J)
	    IF(REV_V(I) .GT. 0.1D0 .AND. REV_V(I) .LT. 0.67D0*REV_V(I-1))THEN
	      J=JST
	      DO WHILE (0.67D0*REV_V(I-1) .LT. V(J))
	        IF(J .EQ. ND)EXIT
	        J=J+1
	       END DO
	      T1=(0.67D0*REV_V(I-1)-V(J-1))/(V(J)-V(J-1))
	      REV_TAU(I)=(1.0D0-T1)*LOG_TAU(J-1)+T1*LOG_TAU(J)
	      REV_V(I)=0.67D0*REV_V(I-1)
	    END IF
 	    J=JST
	  END IF
	  ND_TMP=I
	END DO
!
! The next few line can be uncommented if debuging.
!
!	WRITE(80,'(/,A,/)')' Estimate R grid as determined by DET_R_GRID_V1'
!	WRITE(80,'(A,2(13X,A))')'Index','R','V'
!	DO I=1,ND_TMP
!	  WRITE(80,'(I5,3ES14.4)')I,REV_V(I),REV_TAU(I)
!	END DO
!
! We now create a grid with the requested number of grid points. To do so
! we use the array index as the independent variable. This approach should
! preserve the basic grid spacing.
!
! The scaling for REV_R makes the X-limits identical for the new and old grids.
!
	OLD_TAU(1:ND_TMP)=REV_TAU(1:ND_TMP)
	DO I=1,ND_TMP; OLD_R(I)=I; END DO
	DO I=1,NEW_ND
	  REV_R(I)=1+((I-1.0D0)*(ND_TMP-1.0D0) )/(NEW_ND-1.0D0)
	END DO
	CALL MON_INTERP(REV_TAU,NEW_ND,IONE,REV_R,NEW_ND,
	1                 OLD_TAU,ND_TMP,OLD_R,ND_TMP)
!
! We finish by remembering that we have been working with LOG(TAU).
!
	REV_TAU_GRID(1:NEW_ND)=EXP(REV_TAU(1:NEW_ND))
!
	RETURN
	END
