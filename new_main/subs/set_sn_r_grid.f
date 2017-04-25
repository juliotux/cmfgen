!
! Subroutine to designed to create a NEW R grid given an old R grid, and optical depth scale on this
! grid. The routine places points places points logaritmically in R and TAU. 
!
	SUBROUTINE SET_SN_R_GRID(R,OLD_R,OLD_TAU,IB_RAT,OB_RAT,N_IB_INS,N_OB_INS,ND,NS)
	IMPLICIT NONE
!
! Altered: 16-Jul-2013
!
	INTEGER NS
	INTEGER ND
!
	REAL*8 R(ND)			!Returned
	REAL*8 OLD_R(NS)		!Input
	REAL*8 OLD_TAU(NS)		!Input
!
	REAL*8 LOG_OLD_R(NS)
	REAL*8 LOG_OLD_TAU(NS)
!
	REAL*8 XN(ND)                   !Used as integer grid
	REAL*8 ZN(2*ND)			!Used as integer grid.
	REAL*8 LOG_R(2*ND)
	REAL*8 TAU(2*ND)
	REAL*8 LOG_TAU(2*ND)
!
	REAL*8 IB_RAT
	REAL*8 OB_RAT
!
	INTEGER N_IB_INS
	INTEGER N_OB_INS
!
	REAL*8 dTAU
	REAL*8 dTAU_OLD
	REAL*8 dLOGR
	REAL*8 LOG_TAU_MIN
	REAL*8 LOG_R_MAX
	REAL*8 TAU_BEG,TAU_END
	REAL*8 T1,T2
	REAL*8 NEXT_R
	REAL*8 OB_RAT_LOC
!
	INTEGER LU
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LUER=6
	INTEGER I,J,J_SAV
	INTEGER ND_TMP
!
	WRITE(6,'(A)')
	WRITE(6,'(A)')'Entering SET_SN_R_GRID to define R grid'
	WRITE(6,'(A)')'See R_GRID_SELECTION for computational information'
	WRITE(6,'(A)')
!
	LOG_OLD_TAU=LOG(OLD_TAU)
	LOG_OLD_R=LOG(OLD_R)
!
	CALL GET_LU(LU,'SET_SN_R_GRID')
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
! Estimate average dTAU spacing, first making a correction for the outer boundary.
!
	J=ND-1-N_OB_INS-N_IB_INS
	dTAU=(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1)-N_OB_INS*LOG(OB_RAT_LOC))/J
!
! Estimate average dR spacing, first making a estimate of the correction for the outer boundary.
!
	LOG_TAU_MIN=LOG_OLD_TAU(1)+N_OB_INS*LOG(OB_RAT_LOC)
	I=1
	DO WHILE(LOG_TAU_MIN .GT. LOG_OLD_TAU(I+1))
	  I=I+1
	END DO
	T1=(LOG_TAU_MIN-LOG_OLD_TAU(I))/(LOG_OLD_TAU(I+1)-LOG_OLD_TAU(I))
	LOG_R_MAX=T1*LOG_OLD_R(I+1)+(1.0D0-T1)*LOG_OLD_R(I)
	dLOGR=(LOG_R_MAX-LOG_OLD_R(NS))/J
!
	WRITE(LU,'(A,ES12.4)')' Average spacing in Log(tau)  is:',dTAU
	WRITE(LU,'(A,ES12.4)')' Average spacing in Log(R)    is:',dLOGR
	WRITE(LU,'(A,ES12.4)')' Outer boudary step ratio     is:',OB_RAT_LOC
	WRITE(LU,'(A,2ES12.4)')' New minimum optical depth    is:',EXP(LOG_TAU_MIN),LOG_TAU_MIN
	WRITE(LU,'(A,2ES12.4)')' New maximum radius           is:',EXP(LOG_R_MAX),LOG_R_MAX
!
! Define the new radius grid. The step size in R corresponds to the smaller of
! dLOGR and dLOG_TAU. We create a "uniform" grid. The finer grid at the inner and
! outer boudaries is now only set when we set the FINAL grid.
!
	J=1; I=1
	LOG_TAU(1)=LOG_TAU_MIN
	LOG_R(1)=LOG_R_MAX
	DO WHILE(1 .EQ. 1)
	  I=I+1
	  IF(I .GT. 2*ND)THEN
	    WRITE(LUER,*)' Error in SET_SN_R_GRID --- LOG_R and TAU vectors too small'
	    WRITE(LUER,*)' I=',I,'J=',J
	    WRITE(LUER,*)' Log R(I)=',LOG_R(I)
	    WRITE(LUER,*)' Log old R(I)=',LOG_OLD_R(I)
	    WRITE(LUER,*)' Error in SET_RV_HYDRO_MODEL_V3 --- LOG_R and TAU vectors too small'
	    STOP
	  END IF
!
	  DO WHILE(LOG_OLD_R(J+1) .GT. LOG_R(I-1))
	    J=J+1
	  END DO
	  T1=(LOG_R(I-1)-LOG_OLD_R(J+1))/(LOG_OLD_R(J)-LOG_OLD_R(J+1))
	  TAU_BEG=T1*LOG_OLD_TAU(J)+(1.0D0-T1)*LOG_OLD_TAU(J+1)
	  LOG_TAU(I-1)=TAU_BEG
!
! Bu default, we choose a step size dR. We then compute the step size in LOG(Tau) space.
! We choose a different step size if:
!    (1) dTAU is larger than the average step size.
!    (2) The change in dTAU from the prvious step is too large.
! We compute both a new R and TAU frid, although the TAU grid is primarily used for putput.
!
	  J_SAV=J
	  NEXT_R=LOG_R(I-1)-dLOGR
	  DO WHILE(LOG_OLD_R(J+1) .GT. NEXT_R)
	    J=J+1
	  END DO
	  T1=(NEXT_R-LOG_OLD_R(J+1))/(LOG_OLD_R(J)-LOG_OLD_R(J+1))
	  TAU_END=T1*LOG_OLD_TAU(J)+(1.0D0-T1)*LOG_OLD_TAU(J+1)
!
	  T1=TAU_END-TAU_BEG
	  IF(I .EQ. 2)dTAU_OLD=T1
!	  IF(TAU_END-TAU_BEG .GT. dTAU)THEN
	  IF(T1 .GT. dTAU .AND. dTAU/dTAU_OLD .LE. 1.4D0)THEN
	    NEXT_R=LOG_R(I-1)-dLOGR*dTAU/(TAU_END-TAU_BEG)
	    dTAU_OLD=dTAU
	  ELSE IF(T1/dTAU_OLD .GT. 1.4D0)THEN
	    dTAU_OLD=dTAU_OLD*1.4D0
	    NEXT_R=LOG_R(I-1)-dLOGR*dTAU_OLD/(TAU_END-TAU_BEG)
	  ELSE IF(dTAU_OLD/T1 .GT. 1.4D0)THEN
	    dTAU_OLD=dTAU_OLD*0.71D0
	    NEXT_R=LOG_R(I-1)-dLOGR*dTAU_OLD/(TAU_END-TAU_BEG)
	  ELSE
	    dTAU_OLD=TAU_END-TAU_BEG
	  END IF
	  J=J_SAV
	  LOG_R(I)=NEXT_R
!
! Check whether close enough to inner bondary.
!
	  IF(LOG_R(I)-1.5D0*(LOG_R(I-1)-LOG_R(I)) .LT. LOG_OLD_R(NS))EXIT
	END DO
!
! Place points at the outer boundary.
!
	T1=LOG_R(I-1)-LOG_OLD_R(NS)
	T2=LOG_TAU(I-1)-LOG_OLD_TAU(NS)
	LOG_R(I)=LOG_R(I-1)-0.5D0*T1
	LOG_TAU(I)=LOG_TAU(I-1)-0.5D0*T2
	I=I+1
	ND_TMP=I
	LOG_R(ND_TMP)=LOG_OLD_R(NS)
	LOG_TAU(ND_TMP)=LOG_OLD_TAU(NS)
!
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A)')' First pass at creating new grid. As this grid will generally have too many '
	WRITE(LU,'(A)')' grid points, we will use interpolaiton to create a smaller grid.'
	WRITE(LU,'(A)')' Note: All logs are natural.'
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A,17X,A,8X,A,7X,A,11X,A,6X,A,5X,A,3X,A)')
	1           ' Depth','R','Log(R)','dLog(R)','Tau','Log(Tau)','dLog(Tau)','dTAU[I/I-1]'
	TAU(1:ND_TMP)=EXP(LOG_TAU(1:ND_TMP))
	DO I=1,ND_TMP-1
	  IF(I .NE. 1)T1=(TAU(I+1)-TAU(I))/(TAU(I)-TAU(I-1))
	  WRITE(LU,'(I6,ES18.8,6ES14.4)')I,EXP(LOG_R(I)),LOG_R(I),LOG_R(I+1)-LOG_R(I),
	1              TAU(I),LOG_TAU(I),LOG_TAU(I+1)-LOG_TAU(I),T1
	END DO
	I=ND_TMP
	WRITE(LU,'(I6,ES18.8,6ES14.4)')I,R(I),LOG_R(I),0.0D0,TAU(I),LOG_TAU(I),0.0D0
!
! We now rescale the grid to have the correct number of grid points.
! We use R as a temporary vector for LOG R, and then LOG TAU.
!
	J=ND-N_IB_INS-N_OB_INS
	DO I=1,ND_TMP; ZN(I)=I; END DO
	DO I=1,J
	  XN(I)=1.0D0+(I-1.0D0)*(ND_TMP-1.0D0)/(J-1.0D0)
	END DO
	CALL MON_INTERP(R,J,IONE,XN,J,LOG_R,ND_TMP,ZN,ND_TMP)
	LOG_R=R
	CALL MON_INTERP(R,J,IONE,XN,J,LOG_TAU,ND_TMP,ZN,ND_TMP)
	LOG_TAU=R
	ND_TMP=J
!
! Add extra points at inner boundary.
!
	I=ND_TMP
	T1=LOG_R(ND_TMP-1)-LOG_R(ND_TMP)
	T2=LOG_TAU(ND_TMP-1)-LOG_TAU(ND_TMP)
	IF(N_IB_INS .EQ. 1)THEN
	  LOG_R(I)=LOG_OLD_R(NS)+0.2D0*T1
	  LOG_TAU(I)=LOG_OLD_TAU(NS)+0.2D0*T2
	  ND_TMP=I+1
	ELSE IF(N_IB_INS .EQ. 2)THEN
	  LOG_R(I+1)=LOG_OLD_R(NS)+0.1D0*T1     !0.1D0
	  LOG_R(I)=LOG_OLD_R(NS)+0.4D0*T1      !0.4D0
	  LOG_TAU(I+1)=LOG_OLD_TAU(NS)+0.1D0*T2     !0.1D0
	  LOG_TAU(I)=LOG_OLD_TAU(NS)+0.4D0*T2      !0.4D0
	  ND_TMP=I+2
	ELSE IF(N_IB_INS .EQ. 3)THEN
	  LOG_R(I+2)=LOG_OLD_R(NS)+0.06D0*T1
	  LOG_R(I+1)=LOG_OLD_R(NS)+0.16D0*T1
	  LOG_R(I)=LOG_OLD_R(NS)+0.4D0*T1
	  LOG_TAU(I+2)=LOG_OLD_TAU(NS)+0.06D0*T2
	  LOG_TAU(I+1)=LOG_OLD_TAU(NS)+0.16D0*T2
	  LOG_TAU(I)=LOG_OLD_TAU(NS)+0.4D0*T2
	  ND_TMP=I+3
	ELSE
	  WRITE(6,*)'Error in SET_SN_R_GRID -- this does not work'
	  STOP
	  T2=1.0
	  DO J=1,N_IB_INS
	    T2=IB_RAT*T2+1.0D0
	  END DO
	  T1=T1/T2
	  DO J=N_IB_INS,1,-1
	    LOG_R(I+J-1)=LOG_OLD_R(NS)+T1
	    T1=T1*IB_RAT
	  END DO
	  ND_TMP=I+N_IB_INS
	END IF
	LOG_R(ND_TMP)=LOG_OLD_R(NS)
	LOG_TAU(ND_TMP)=LOG_OLD_TAU(NS)
!
! Add finer grid at outer boundary.
!
! Shift grid to allow for insertion of extra ponts
!
	DO I=ND_TMP,1,-1
	  LOG_R(I+N_OB_INS+1)=LOG_R(I)
	  LOG_TAU(I+N_OB_INS+1)=LOG_TAU(I)
	END DO
!
	T1=(LOG_OLD_R(1)-LOG_R_MAX)/(N_OB_INS)
	T2=(LOG_OLD_TAU(1)-LOG_TAU(1))/(N_OB_INS)
	WRITE(6,*)T1,T2
	DO I=1,N_OB_INS-1
	   LOG_R(N_OB_INS+2-I)=LOG_R(N_OB_INS+2)+I*T1
	   LOG_TAU(N_OB_INS+2-I)=LOG_TAU(N_OB_INS+2)+I*T2
	END DO
	LOG_R(1)=LOG_OLD_R(1)
	LOG_R(2)=LOG_R(1)-0.05D0*T1
	LOG_TAU(1)=LOG_OLD_TAU(1)
	LOG_TAU(2)=LOG_TAU(1)-0.05D0*T2
	LOG_R(N_OB_INS+2)=0.4D0*LOG_R(N_OB_INS+1)+0.6D0*LOG_R(N_OB_INS+3)
	LOG_TAU(N_OB_INS+2)=0.4D0*LOG_TAU(N_OB_INS+1)+0.6D0*LOG_TAU(N_OB_INS+3)
!
	R=EXP(LOG_R(1:ND))
	TAU(1:ND)=EXP(LOG_TAU(1:ND))
	R(1)=OLD_R(1)
	R(ND)=OLD_R(NS)
	WRITE(LUER,*)'Computed R grid in SET_SN_R_GRID'
!
! Make sure grid is monotonic.
!
	DO I=1,ND-1
	  IF(R(I) .LE. R(I+1))THEN
	    WRITE(6,*)' Error in SET_SN_R_GRID: R grid not monotonic'
	    WRITE(6,'(I4,3ES20.8)')I,R(MAX(1,I-1)),R(I),R(MIN(I+1,ND))
	    STOP
	  END IF
	END DO
!
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A)')' Final grid computed with SET_SN_R_GRID '
	WRITE(LU,'(A)')' '
	WRITE(LU,'(A,17X,A,8X,A,7X,A,11X,A,6X,A,5X,A,3X,A)')
	1           ' Depth','R','Log(R)','dLog(R)','Tau','Log(Tau)','dLog(Tau)','dTAU[I/I-1]'
	T1=0.0D0
	DO I=1,ND-1
	  IF(I .NE. 1)T1=(TAU(I+1)-TAU(I))/(TAU(I)-TAU(I-1))
	  WRITE(LU,'(I6,ES18.8,6ES14.4)')I,R(I),LOG_R(I),LOG_R(I+1)-LOG_R(I),
	1              TAU(I),LOG_TAU(I),LOG_TAU(I+1)-LOG_TAU(I),T1
	END DO
	WRITE(LU,'(I6,ES18.8,6ES14.4)')ND,R(ND),LOG_R(ND),0.0D0,TAU(ND),LOG_TAU(ND),0.0D0
	CLOSE(UNIT=LU)
!
	RETURN
	END
