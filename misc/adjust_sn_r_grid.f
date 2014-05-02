!
! Subroutine to designed to create a NEW R grid given an old R grid, and optical depth scale on this
! grid. The routine places points places points logaritmically in R and TAU. 
!
	SUBROUTINE ADJUST_SN_R_GRID(R,OLD_R,OLD_T,OLD_TAU,R_SCALE_FAC,dLOGT_MAX,
	1           IB_RAT,OB_RAT,DTAU2_ON_DTAU1,N_IB_INS,N_OB_INS,ND,NS)
	IMPLICIT NONE
!
! Altered: 07-Jan-2014 -- Changes to call, and extensive improvements made.
! Altered: 16-Jul-2013
!
	INTEGER NS			!For old grid
	INTEGER ND			!For final grid.
!
	REAL*8 R(ND)			!Returned
	REAL*8 OLD_R(NS)		!Input
	REAL*8 OLD_T(NS)		!Input
	REAL*8 OLD_TAU(NS)		!Input
!
	REAL*8 LOG_OLD_R(NS)
	REAL*8 LOG_OLD_T(NS)
	REAL*8 LOG_OLD_TAU(NS)
!
	INTEGER, PARAMETER :: IND=2
	REAL*8 XN(ND)                   !Used as integer grid
	REAL*8 ZN(IND*ND)		!Used as integer grid.
	REAL*8 LOG_R(IND*ND)
	REAL*8 TAU(IND*ND)
	REAL*8 LOG_T(IND*ND)
	REAL*8 LOG_TAU(IND*ND)
!
	REAL*8 IB_RAT			!dTAU(I+1)/dTAU(I) near inner boudary.
	REAL*8 OB_RAT			!dTAU(I+1)/dRAU(I) near outer bunary
	REAL*8 R_SCALE_FAC		!Factor to increase dLOGR by so that TAU scale has higher importance.
	REAL*8 DTAU2_ON_DTAU1		!DTAU(2)/DTAU(1) at outer boundary.
	REAL*8 dLOGT_MAX
!
	INTEGER N_IB_INS
	INTEGER N_OB_INS
!
	REAL*8 dTAU			!d(LOG(TAU))
	REAL*8 dTAU_OLD
	REAL*8 dLOGR
	REAL*8 dLOGT
	REAL*8 LOG_TAU_MIN
	REAL*8 LOG_R_MAX
	REAL*8 TAU_BEG,TAU_END
	REAL*8 T1,T2,T3
	REAL*8 NEXT_R
	REAL*8 OB_RAT_LOC
!
	INTEGER LU
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LUER=6
	INTEGER I,J,K,L,J_SAV
	INTEGER ND_TMP
	INTEGER ICNT
!
	LOGICAL CHECK_T
!
	WRITE(6,'(A)')
	WRITE(6,'(A)')'Entering ADJUST_SN_R_GRID to define R grid'
	WRITE(6,'(A)')'See R_GRID_SELECTION for computational information'
	WRITE(6,'(A)')
!
	LOG_OLD_TAU=LOG(OLD_TAU)
	LOG_OLD_R=LOG(OLD_R)
	IF(dLOGT_MAX .GT. 0 .AND. dLOGT_MAX .LT. 1.D0)THEN
	  CHECK_T=.TRUE.
	  LOG_OLD_T(1:ND)=LOG(OLD_T(1:ND))
	ELSE
	  CHECK_T=.FALSE.
	END IF
!
	CALL GET_LU(LU,'ADJUST_SN_R_GRID')
	OPEN(UNIT=LU,FILE='R_GRID_SELECTION',STATUS='UNKNOWN',ACTION='WRITE')
!
	WRITE(LU,'(A,2ES12.4)')' Outer boundary optical depth is:',OLD_TAU(1),LOG_OLD_TAU(1)
	WRITE(LU,'(A,2ES12.4)')' Inner boundary optical depth is:',OLD_TAU(NS),LOG_OLD_TAU(NS)
!
! Estimate spacing to get required grid spacing.
!
	dTAU=(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1))/(ND-1)
	OB_RAT_LOC=MAX(OB_RAT,EXP(dTAU))
!
! Estimate average dTAU spacing.
!
	J=ND-1-N_OB_INS-N_IB_INS
	dTAU=(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1))/J
!
! Estimate average dR spacing, first making a estimate of the correction for the outer boundary.
!
	LOG_TAU_MIN=LOG_OLD_TAU(1)
	LOG_R_MAX=LOG_OLD_R(1)
	dLOGR=R_SCALE_FAC*(LOG_R_MAX-LOG_OLD_R(NS))/J
	dLOGT=dLOGT_MAX
!
! We do the loop twice so that we can get a close agreement between ND_TMP
! (number of depth points with default spacing) and ND (desird grid size).
!
	DO ICNT=1,2
!
	  IF(ICNT .NE. 1)THEN
	    T1=DFLOAT(ND_TMP-1)/DFLOAT(ND-N_IB_INS-N_OB_INS-1)
	    dTAU=dTAU*T1
	    dLOGR=dLOGR*T1
	    dLOGT=dLOGT*T1
	  END IF
!
	  WRITE(LU,'(A,ES12.4)') ' Average spacing in Ln(tau)      is:',dTAU
	  WRITE(LU,'(A,ES12.4)') ' Average spacing in Ln(R)        is:',dLOGR
	  WRITE(LU,'(A,ES12.4)') ' Maximum spacing in Ln(T)        is:',dLOGT
	  WRITE(LU,'(A,ES12.4)') ' Outer boudary step ratio        is:',OB_RAT_LOC
	  WRITE(LU,'(A,I3)')     ' Number of points inserted at IB is:',N_IB_INS
	  WRITE(LU,'(A,I3)')     ' Number of points inserted at OB is:',N_OB_INS
	  WRITE(LU,'(A,2ES12.4)')' New minimum optical depth       is:',EXP(LOG_TAU_MIN),LOG_TAU_MIN
	  WRITE(LU,'(A,2ES12.4)')' New maximum radius              is:',EXP(LOG_R_MAX),LOG_R_MAX
!
! Define the new radius grid. The step size in R corresponds to the smaller of
! dLOGR and dLOG_TAU. We create a "uniform" grid. The finer grid at the inner and
! outer boudaries is now only set when we set the FINAL grid.
!
	  J=1; I=1
	  LOG_TAU(1)=LOG_TAU_MIN
	  LOG_R(1)=LOG_R_MAX
	  LOG_T(1)=LOG_OLD_T(1)
!
	  DO WHILE(1 .EQ. 1)
	    I=I+1
	    IF(I .GT. IND*ND)THEN
	      I=I-1
	      WRITE(LUER,*)' Error in SET_SN_R_GRID --- LOG_R and TAU vectors too small'
	      WRITE(LUER,*)' I=',I,'J=',J
	      WRITE(LUER,*)' Log TAU(I)=',LOG_TAU(I)
	      WRITE(LUER,*)' Log R(I)=',LOG_R(I)
	      WRITE(LUER,*)' Log old R(I)=',LOG_OLD_R(I)
	      WRITE(LUER,*)' Error in SET_RV_HYDRO_MODEL_V3 --- LOG_R and TAU vectors too small'
	      STOP
	    END IF
!
! Bu default, we choose a step size dR. We then compute the step size in LOG(Tau) space.
! We choose a different step size if:
!    (1) dTAU is larger than the average step size.
!    (2) The change in dTAU from the prvious step is too large.
! We compute both a new R and TAU grid, although the TAU grid is primarily used for output.
!
	    NEXT_R=LOG_R(I-1)-dLOGR
	    TAU_END=LOG_TAU(I-1)+dTAU
	    J=1
	    DO WHILE(LOG_OLD_R(J+1) .GT. NEXT_R)
	      J=J+1
	    END DO
	    T1=(NEXT_R-LOG_OLD_R(J))/(LOG_OLD_R(J+1)-LOG_OLD_R(J))
	    T2=T1*LOG_OLD_TAU(J+1)+(1.0D0-T1)*LOG_OLD_TAU(J)
	    T2=T2-LOG_TAU(I-1)
	    IF(T2 .GT. dTAU)THEN
	        J=1
	      DO WHILE(LOG_OLD_TAU(J+1) .LT. TAU_END)
	        J=J+1
	      END DO
	      T1=(TAU_END-LOG_OLD_TAU(J))/(LOG_OLD_TAU(J+1)-LOG_OLD_TAU(J))
	      LOG_R(I)=T1*LOG_OLD_R(J+1)+(1.0D0-T1)*LOG_OLD_R(J)
	      LOG_TAU(I)=TAU_END
	    ELSE
	      LOG_R(I)=NEXT_R
	      J=1
	      DO WHILE(LOG_OLD_R(J+1) .GT. NEXT_R)
	        J=J+1
	      END DO
	      T1=(NEXT_R-LOG_OLD_R(J))/(LOG_OLD_R(J+1)-LOG_OLD_R(J))
	      LOG_TAU(I)=T1*LOG_OLD_TAU(J+1)+(1.0D0-T1)*LOG_OLD_TAU(J)
	    END IF
!
! We have to be careful as T may not be monotonic.
!
	    IF(CHECK_T)THEN
	      T1=(LOG_R(I)-LOG_OLD_R(J))/(LOG_OLD_R(J+1)-LOG_OLD_R(J))
	      LOG_T(I)=T1*LOG_OLD_T(J+1)+(1.0D0-T1)*LOG_OLD_T(J)
	      DO WHILE(ABS(LOG_T(I)-LOG_T(I-1)) .GT. dLOGT)
	        LOG_R(I)=LOG_R(I)+0.1D0*(LOG_R(I-1)-LOG_R(I))
	        J=1
	        DO WHILE(LOG_OLD_R(J+1) .GT. LOG_R(I))
	          J=J+1
	        END DO
	        WRITE(6,'(I5,4ES16.6)')I,LOG_OLD_R(J+1),LOG_R(I),LOG_OLD_R(J)
	        WRITE(6,'(I5,4ES16.8)')I,LOG_R(I),LOG_TAU(I),LOG_T(I),ABS(LOG_T(I)-LOG_T(I-1))
	        T1=(LOG_R(I)-LOG_OLD_R(J))/(LOG_OLD_R(J+1)-LOG_OLD_R(J))
	        LOG_TAU(I)=T1*LOG_OLD_TAU(J+1)+(1.0D0-T1)*LOG_OLD_TAU(J)
	        LOG_T(I)=T1*LOG_OLD_T(J+1)+(1.0D0-T1)*LOG_OLD_T(J)
	      END DO
	    END IF
!
! Check whether close enough to inner bondary.
!
	    IF(LOG_R(I)-1.5D0*(LOG_R(I-1)-LOG_R(I)) .LT. LOG_OLD_R(NS))EXIT
	    IF(LOG_TAU(I)+1.5D0*(LOG_TAU(I-1)-LOG_TAU(I)) .GT. LOG_OLD_TAU(NS))EXIT
	  END DO
	  WRITE(6,'(A,I4,4Z,A,I4)')'Grid it=',ICNT,'ND=',I+1
	  ND_TMP=I+1
!
	  IF(ND_TMP .EQ. ND-N_IB_INS-N_OB_INS)EXIT
	END DO
!
! Define grid accurately at the inner boundary.
!
	ND_TMP=I+1
	LOG_R(ND_TMP)=LOG_OLD_R(NS)
	LOG_TAU(ND_TMP)=LOG_OLD_TAU(NS)
!
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A)')' First pass at creating new grid. As this grid will generally have too many '
	WRITE(LU,'(A)')' grid points, we will use interpolaiton to create a smaller grid.'
	WRITE(LU,'(A)')' Note: All logs are natural.'
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A,17X,A,9X,A,8X,A,11X,A,10X,A,7X,A,6X,A,3X,A)')
	1           ' Depth','R','Ln(R)','dLn(R)','Tau','dTAU','Ln(Tau)','dLn(Tau)','dTAU[I/I-1]'
	TAU(1:ND_TMP)=EXP(LOG_TAU(1:ND_TMP))
	DO I=1,ND_TMP-1
	  IF(I .NE. 1)T1=(TAU(I+1)-TAU(I))/(TAU(I)-TAU(I-1))
	  WRITE(LU,'(I6,ES18.8,7ES14.4)')I,EXP(LOG_R(I)),LOG_R(I),LOG_R(I+1)-LOG_R(I),
	1              TAU(I),TAU(I+1)-TAU(I),LOG_TAU(I),LOG_TAU(I+1)-LOG_TAU(I),T1
	END DO
	I=ND_TMP
	WRITE(LU,'(I6,ES18.8,7ES14.4)')I,EXP(LOG_R(I)),LOG_R(I),0.0D0,TAU(I),0.0D0,LOG_TAU(I),0.0D0
!
! We now rescale the grid to have the correct number of grid points.
! We use R as a temporary vector for LOG R, and then LOG TAU.
!
	J=ND-N_IB_INS-N_OB_INS
	WRITE(6,*)'Number of depth points in initial grid',ND_TMP
	WRITE(6,*)'Number of points required (corrected for boundary insertions)',J
	IF(ND_TMP .NE. J)THEN
	  DO I=1,ND_TMP; ZN(I)=I; END DO
	  DO I=1,J
	    XN(I)=1.0D0+(I-1.0D0)*(ND_TMP-1.0D0)/(J-1.0D0)
	  END DO
	  CALL MON_INTERP(R,J,IONE,XN,J,LOG_R,ND_TMP,ZN,ND_TMP)
	  LOG_R=R
	  CALL MON_INTERP(R,J,IONE,XN,J,LOG_TAU,ND_TMP,ZN,ND_TMP)
	  LOG_TAU=R
	  ND_TMP=J
	END IF
!
! Add extra points at inner boundary.
!
	I=ND_TMP
	T1=LOG_R(ND_TMP-1)-LOG_R(ND_TMP)
	T2=LOG_TAU(ND_TMP-1)-LOG_TAU(ND_TMP)
	IF(N_IB_INS .EQ. 1 .AND. IB_RAT .EQ. 0.0D0)THEN
	  LOG_R(I)=LOG_OLD_R(NS)+0.2D0*T1
	  LOG_TAU(I)=LOG_OLD_TAU(NS)+0.2D0*T2
	  ND_TMP=I+1
	ELSE IF(N_IB_INS .EQ. 2 .AND. IB_RAT .EQ. 0.0D0)THEN
	  LOG_R(I+1)=LOG_OLD_R(NS)+0.1D0*T1     !0.1D0
	  LOG_R(I)=LOG_OLD_R(NS)+0.4D0*T1      !0.4D0
	  LOG_TAU(I+1)=LOG_OLD_TAU(NS)+0.1D0*T2     !0.1D0
	  LOG_TAU(I)=LOG_OLD_TAU(NS)+0.4D0*T2      !0.4D0
	  ND_TMP=I+2
	ELSE IF(N_IB_INS .EQ. 3 .AND. IB_RAT .EQ. 0.0D0)THEN
	  LOG_R(I+2)=LOG_OLD_R(NS)+0.06D0*T1
	  LOG_R(I+1)=LOG_OLD_R(NS)+0.16D0*T1
	  LOG_R(I)=LOG_OLD_R(NS)+0.4D0*T1
	  LOG_TAU(I+2)=LOG_OLD_TAU(NS)+0.06D0*T2
	  LOG_TAU(I+1)=LOG_OLD_TAU(NS)+0.16D0*T2
	  LOG_TAU(I)=LOG_OLD_TAU(NS)+0.4D0*T2
	  ND_TMP=I+3
	ELSE
	  T1=EXP(LOG_TAU(ND_TMP))
	  T2=EXP(LOG_TAU(ND_TMP))-EXP(LOG_TAU(ND_TMP-1)) 
	  T3=T2*(IB_RAT-1)/(IB_RAT**(N_IB_INS+1)-1)
	  DO J=1,N_IB_INS
	    T1=T1-T3
	    T3=T3*IB_RAT
	    WRITE(6,*)LOG_OLD_TAU(NS),T1
	    LOG_TAU(ND_TMP+N_IB_INS-J)=LOG(T1)
	  END DO
!
	  T1=EXP(LOG_R(ND_TMP))
	  T2=EXP(LOG_R(ND_TMP))-EXP(LOG_R(ND_TMP-1)) 
	  T3=T2*(IB_RAT-1)/(IB_RAT**(N_IB_INS+1)-1)
	  DO J=1,N_IB_INS
	    T1=T1-T3
	    T3=T3*IB_RAT
	    LOG_R(ND_TMP+N_IB_INS-J)=LOG(T1)
	    WRITE(6,*)LOG_OLD_R(NS),T1
	  END DO
	  ND_TMP=ND_TMP+N_IB_INS
	END IF
	LOG_R(ND_TMP)=LOG_OLD_R(NS)
	LOG_TAU(ND_TMP)=LOG_OLD_TAU(NS)
!
! Add finer grid at outer boundary.
!
	WRITE(6,*)LOG_R(1),LOG_OLD_R(1)
	WRITE(6,*)LOG_TAU(1),LOG_OLD_TAU(1)
!
! Shift grid to allow for insertion of extra ponts
!
! We use K to allow for the possibility that we don't need a very find
! set right at the outer boudary.
!
	K=0
	IF(DTAU2_ON_DTAU1 .LT. 2.0D0)K=1
	T1=EXP(LOG_TAU(1))
	T2=EXP(LOG_TAU(1))-EXP(LOG_TAU(2)) 
	T3=T2*(OB_RAT-1)/(OB_RAT**(N_OB_INS+K)-1)
	DO I=ND_TMP,1,-1
	    LOG_TAU(I+N_OB_INS)=LOG_TAU(I)
	END DO
	DO J=1,N_OB_INS
	  IF(J .EQ. 1 .AND. K .EQ. 0)THEN
	    LOG_TAU(2)=LOG(T1-T3/DTAU2_ON_DTAU1)
	  ELSE
	    T1=T1-T3
	    T3=T3*OB_RAT
	    LOG_TAU(1+J)=LOG(T1)
	  END IF
	END DO
!
	T1=EXP(LOG_R(1))
	T2=EXP(LOG_R(1))-EXP(LOG_R(2)) 
	T3=T2*(OB_RAT-1)/(OB_RAT**(N_OB_INS+K)-1)
	DO I=ND_TMP,1,-1
	  LOG_R(I+N_OB_INS)=LOG_R(I)
	END DO
	DO J=1,N_OB_INS
	  IF(J .EQ. 1 .AND. K .EQ. 0)THEN
	    LOG_R(2)=LOG(T1-T3/DTAU2_ON_DTAU1)
	  ELSE
	    T1=T1-T3
	    T3=T3*OB_RAT
	    LOG_R(1+J)=LOG(T1)
	  END IF
	END DO
	R=EXP(LOG_R(1:ND))
	TAU(1:ND)=EXP(LOG_TAU(1:ND))
	R(1)=OLD_R(1)
	R(ND)=OLD_R(NS)
	WRITE(LUER,*)'Computed R grid in SET_SN_R_GRID'
!
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A)')' Final grid computed with SET_SN_R_GRID '
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A,17X,A,9X,A,8X,A,11X,A,10X,A,7X,A,6X,A,3X,A)')
	1           ' Depth','R','Ln(R)','dLn(R)','Tau','dTau','Ln(Tau)','dLn(Tau)','dTAU[I/I-1]'
	T1=0.0D0
	DO I=1,ND-1
	  IF(I .NE. 1)T1=(TAU(I+1)-TAU(I))/(TAU(I)-TAU(I-1))
	  WRITE(LU,'(I6,ES18.8,7ES14.4)')I,R(I),LOG_R(I),LOG_R(I+1)-LOG_R(I),
	1              TAU(I),TAU(I+1)-TAU(I),LOG_TAU(I),LOG_TAU(I+1)-LOG_TAU(I),T1
	END DO
	WRITE(LU,'(I6,ES18.8,7ES14.4)')ND,R(ND),LOG_R(ND),0.0D0,TAU(ND),0.0D0,LOG_TAU(ND),0.0D0
	CLOSE(UNIT=LU)
!
! Make sure grid is monotonic. We do it here since we can check the actual grid
! computed in the output file.
!
	DO I=1,ND-1
	  IF(R(I) .LE. R(I+1))THEN
	    WRITE(6,*)' Error in SET_SN_R_GRID: R grid not monotonic'
	    WRITE(6,'(I4,3ES20.8)')I,R(MAX(1,I-1)),R(I),R(MIN(I+1,ND))
	    STOP
	  END IF
	END DO
!
	RETURN
	END
