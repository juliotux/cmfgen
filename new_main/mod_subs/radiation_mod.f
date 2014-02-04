!
! Module for vectors required when computing the continuum radiation field,
! with blanketing. Some are simply scratch vectors. Others
! are required in CMFGEN_SUB and in the variation routines.
! Vectors required for dJ are in VAR_RAD_MOD.
!
	MODULE RADIATION_MOD
!
	REAL*8, ALLOCATABLE :: RJ(:)      	!ND - Mean intensity  (ND)
	REAL*8, ALLOCATABLE :: RJ_ES(:)   	!ND - Convolution of RJ with e.s. R(v'v')
!
	REAL*8, ALLOCATABLE :: J_INT(:)   	!ND - Frequency integrated J
	REAL*8, ALLOCATABLE :: K_INT(:)   	!ND - Frequency integrated K
	REAL*8, ALLOCATABLE :: K_MOM(:)   	!ND - Frequency dependent K moment
!
! Arrays for calculating mean opacities.
!
	REAL*8, ALLOCATABLE :: INT_dBdT(:)     	!ND - Int. of dB/dT dv (to calculate ROSSMEAN)
!
	REAL*8, ALLOCATABLE :: RLUMST(:)      	!ND - Luminosity as a function of depth
	REAL*8, ALLOCATABLE :: MECH_LUM(:)     	!ND - Mechanical luminosity
	REAL*8, ALLOCATABLE :: SOB(:)      	!ND - Used in computing continuum flux
	REAL*8, ALLOCATABLE :: LLUMST(:)      	!ND - Line luminosity.
	REAL*8, ALLOCATABLE :: DIELUM(:)      	!ND - Dielectronic line emission luminosity.
	REAL*8, ALLOCATABLE :: DEP_RAD_EQ(:)    !ND - Integrated departure from radiative equilibrium
	REAL*8, ALLOCATABLE :: DJDt_FLUX(:)     !ND - DJDT correction to integrated flux.
	REAL*8, ALLOCATABLE :: DJDt_TERM(:)     !ND - 
!
! Vector giving the MINIMUM Doppler width at each depth.
!
	REAL*8, ALLOCATABLE :: VDOP_VEC(:)      !ND
!
! Transfer equation vectors
!
	REAL*8, ALLOCATABLE :: Z(:)        	!NDMAX - Z displacement along a given array
	REAL*8, ALLOCATABLE :: TA(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: TB(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: TC(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: XM(:)        	!NDMAX - R.H.S. (SOURCE VECTOR)
	REAL*8, ALLOCATABLE :: DTAU(:)        	!NDMAX - Optical depth (used in error calcs)
	REAL*8, ALLOCATABLE :: dCHIdR(:)        !NDMAX - Derivative of opacity.
!
! Continuum matrices
!
	REAL*8, ALLOCATABLE :: WM(:,:)        	!ND,ND - Coef. matrix of J & %J vector
	REAL*8, ALLOCATABLE :: FB(:,:)        	!ND,ND - Coef. of J & %J vects in angular equ.
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	REAL*8, ALLOCATABLE :: FEDD(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: GEDD(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: QEDD(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: N_ON_J(:)        !NDMAX -
	REAL*8, ALLOCATABLE :: RSQHNU(:)        !NDMAX -
!
! For MOM_J_REL_V2 --- i.e., the inclusion of relatvistic, but not time
! dependence, for the computation of J.
!
	REAL*8, ALLOCATABLE :: N_ON_J_NODE(:)		!
	REAL*8, ALLOCATABLE :: H_ON_J(:)		!
	REAL*8, ALLOCATABLE :: KMID_ON_J(:)		!
	REAL*8, ALLOCATABLE :: dlnJdlNR(:)		!
!
! Boundary conditions.
!
	REAL*8 HBC_CMF(3)
	REAL*8 NBC_CMF(3)
	REAL*8 INBC
	REAL*8 HBC_J
	REAL*8 HBC_S			!Bound. Cond. for JFEAU
	REAL*8 HBC_PREV(3)
	REAL*8 NBC_PREV(3)
	REAL*8 INBC_PREV
!
	REAL*8, ALLOCATABLE :: FEDD_PREV(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: GEDD_PREV(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: N_ON_J_PREV(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: JNU_PREV(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: RSQHNU_PREV(:)        	!NDMAX -
	REAL*8, ALLOCATABLE :: FOLD(:)        		!NDMAX
!
! 
!
! Use for interpolating opacities tec onto the foine grid.
!
	INTEGER, ALLOCATABLE :: INDX(:)                 !NDMAX
	INTEGER, ALLOCATABLE :: POS_IN_NEW_GRID(:)      !ND
	REAL*8, ALLOCATABLE :: COEF(:,:)        	!0:3,NDMAX
!
! Variables and arrays required ofr the fine grid.
!
	REAL*8, ALLOCATABLE :: REXT(:)        		!NDMAX
	REAL*8, ALLOCATABLE :: VEXT(:)        		!NDMAX
	REAL*8, ALLOCATABLE :: TEXT(:)        		!NDMAX
	REAL*8, ALLOCATABLE :: SIGMAEXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: VDOP_VEC_EXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: CHIEXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: ESECEXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: ETAEXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: ZETAEXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: THETAEXT(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: RJEXT(:)        		!NDMAX
	REAL*8, ALLOCATABLE :: RJEXT_ES(:)        	!NDMAX
	REAL*8, ALLOCATABLE :: FEXT(:)        		!NDMAX
	REAL*8, ALLOCATABLE :: QEXT(:)        		!NDMAX
	REAL*8, ALLOCATABLE :: SOURCEEXT(:)        	!NDMAX
!
! Diffusion approximation variables
!
	REAL*8 DTDR
	REAL*8 DBB
	REAL*8 DDBBDT
!
	REAL*8 dLOG_NU
	LOGICAL CONT_VEL
!
	END MODULE RADIATION_MOD
!
	SUBROUTINE SET_RADIATION_MOD(ND,NDMAX,NPMAX)
	USE RADIATION_MOD
	IMPLICIT NONE
	INTEGER ND,NDMAX,NPMAX
!
	INTEGER IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	IOS=0
	IF(IOS .EQ. 0)ALLOCATE ( RJ(ND),STAT=IOS )		!Mean intensity
	IF(IOS .EQ. 0)ALLOCATE ( RJ_ES(ND),STAT=IOS )		!Convolution of RJ with e.s. R(v'v')
	IF(IOS .EQ. 0)ALLOCATE ( J_INT(ND),STAT=IOS )		!Frequency integrated J
	IF(IOS .EQ. 0)ALLOCATE ( K_INT(ND),STAT=IOS )		!Frequency integrated K
	IF(IOS .EQ. 0)ALLOCATE ( K_MOM(ND),STAT=IOS )		!Frequency dependent K moment
!
! Arrays for calculating mean opacities.
!
	IF(IOS .EQ. 0)ALLOCATE ( INT_dBdT(ND),STAT=IOS )  	!Integral of dB/dT over nu (to calculate ROSSMEAN)
!
	IF(IOS .EQ. 0)ALLOCATE ( RLUMST(ND),STAT=IOS )		!Luminosity as a function of depth
	IF(IOS .EQ. 0)ALLOCATE ( MECH_LUM(ND),STAT=IOS )	!Mechanical luminosity
	IF(IOS .EQ. 0)ALLOCATE ( SOB(ND),STAT=IOS )   		!Used in computing continuum flux
	IF(IOS .EQ. 0)ALLOCATE ( LLUMST(ND),STAT=IOS )    	!Line luminosity.
	IF(IOS .EQ. 0)ALLOCATE ( DIELUM(ND),STAT=IOS )    	!Dielectronic line emission luminosity.
	IF(IOS .EQ. 0)ALLOCATE ( DEP_RAD_EQ(ND),STAT=IOS )    	!Depature from radiative equilibrium.
	IF(IOS .EQ. 0)ALLOCATE ( DJDt_TERM(ND),STAT=IOS )	!
	IF(IOS .EQ. 0)ALLOCATE ( DJDt_FLUX(ND),STAT=IOS )    	!DJDt correction to integrated flux.
	IF(IOS .EQ. 0)ALLOCATE ( VDOP_VEC(ND),STAT=IOS )
!
! Transfer equation vectors
!
	IF(IOS .EQ. 0)ALLOCATE ( Z(NDMAX),STAT=IOS )		!Z displacement along a given array
	IF(IOS .EQ. 0)ALLOCATE ( TA(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( TB(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( TC(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( XM(NDMAX),STAT=IOS )		!R.H.S. (SOURCE VECTOR)
	IF(IOS .EQ. 0)ALLOCATE ( DTAU(NDMAX),STAT=IOS )       	!Optical depth (used in error calcs)
	IF(IOS .EQ. 0)ALLOCATE ( dCHIdR(NDMAX),STAT=IOS ) 	!Derivative of opacity.
!
! Continuum matrices
!
	IF(IOS .EQ. 0)ALLOCATE ( WM(ND,ND),STAT=IOS )		!Coef. matrix of J & %J vector
	IF(IOS .EQ. 0)ALLOCATE ( FB(ND,ND),STAT=IOS )		!Coef. of J & %J vects in angular equ.
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	IF(IOS .EQ. 0)ALLOCATE ( FEDD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( GEDD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( QEDD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( N_ON_J(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RSQHNU(NDMAX),STAT=IOS )
!
	IF(IOS .EQ. 0)ALLOCATE ( N_ON_J_NODE(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( H_ON_J(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( KMID_ON_J(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dlnJdlnR(NDMAX),STAT=IOS )
!
	IF(IOS .EQ. 0)ALLOCATE ( FEDD_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( GEDD_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( N_ON_J_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( JNU_PREV(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RSQHNU_PREV(NDMAX),STAT=IOS )
!
! 
!
	IF(IOS .EQ. 0)ALLOCATE ( INDX(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( POS_IN_NEW_GRID(ND),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( COEF(0:3,NDMAX),STAT=IOS )
!
	IF(IOS .EQ. 0)ALLOCATE ( REXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( VEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( TEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( SIGMAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( VDOP_VEC_EXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( CHIEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( ESECEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( ETAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( ZETAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( THETAEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RJEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RJEXT_ES(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( FOLD(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( FEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( QEXT(NDMAX),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( SOURCEEXT(NDMAX),STAT=IOS )
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error allocating memory in SET_RADIATION_MOD'
	  WRITE(LUER,*)'STAT=',IOS
	  STOP
	END IF
!
	RETURN
	END
