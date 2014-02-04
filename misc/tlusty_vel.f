!
! Program to match a Beta-velocity law with the hydrostatic density
! structure at depth. Hydrostatic density structure is taken from a
! TLUSTY output file containing the following:
!
! Index, column-mass, Tau Ross, Abs Ross, T, Ne, total density (cgs units).
!
! At the connecting point, both the velocity an velocity gradient match.
! Based on the IDL program VEL3.PRO written by Thierry Lanz.
!
! The velocity law is assumed to have the form
!
! V(r)=Vinf (1-R0/r)**BETA'
!
! To give more flexibility, BETA has the form
!
!       BETA'=BETA+(BETA_MIN-BETA)*EXP( (1.0-R/RCORE)/BETA_SCL )
!
! This resorts to the standard Beta law when BETA_MIN=BETA.
! BETA is BETA in the outer regions of the stellar wind.
! BETA_MIN is BETA in the inner regions of the stellar wind. We allow
! for two choices since BETA_MIN=0.7 gives a connection radius around
! 7 km/s, where as higher BETA gives a connection radii where the
! connection velocity is significantly smaller.
! 
	PROGRAM TLUSTY_VEL
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 23-Jul-2003: Bug fix to the new default grid section -- doesn't effect grid.
! Altered 24-Apr-2003: Bug fix and improvement to the new default grid section.
! Altered 01-Apr-2003: New default option installed to compute RGRID.
!
!
! TLUSTY input
!
	INTEGER ND				!Number of TLUSTY depths input
	INTEGER ND_MAX			!Maximum number of TLUSTY depths in file
	INTEGER, ALLOCATABLE :: INDX(:)
	REAL*8, ALLOCATABLE :: DM(:)		!Column mass density
	REAL*8, ALLOCATABLE :: TAUR(:)		!Rosseland optical depth
	REAL*8, ALLOCATABLE :: AROSS(:)		!Absorption Rosseland optical depth scale
	REAL*8, ALLOCATABLE :: T(:)		!Temperature
	REAL*8, ALLOCATABLE :: ED(:)		!Electron density (/cm^3)
	REAL*8, ALLOCATABLE :: DSH(:)		!Density (gm/cm^3)
!
! Calculated directly from TLUSTY data. Same grid.
!
	REAL*8, ALLOCATABLE :: RD(:)		!Radius
	REAL*8, ALLOCATABLE :: RD_NORM(:)	!Radius normalized by core radius
	REAL*8, ALLOCATABLE :: ZZ(:)            !Height above core
	REAL*8, ALLOCATABLE :: VPH(:)		!Velocity deduced from density structure (assuming Mdot)
	REAL*8, ALLOCATABLE :: dVdR_PH(:)	!dVdR deduced from density structure
	REAL*8, ALLOCATABLE :: KAPPA(:)		!Mass absorption coefficient
!
! Arrays used to generate a FINE grid of the TLUSTY hydrostatic data data.
!
	INTEGER, PARAMETER :: NBIG=2500

	REAL*8, ALLOCATABLE :: RA(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: RA_NORM(:)
	REAL*8, ALLOCATABLE :: VPA(:)		!Velocity deduced from density structure
	REAL*8, ALLOCATABLE :: VW(:)		!Beta-wind velocity
	REAL*8, ALLOCATABLE :: DSHA(:)		!Density
	REAL*8, ALLOCATABLE :: TAUR_A(:)
	REAL*8, ALLOCATABLE :: KAPPA_A(:)
!
! Arrays for grid containing TLUSTY hydrostatic structure merged with
! the Beta-Velocity. Below the connection point TLUSTY structure is used.
! Above the connection point, the Beta-velocity law is used.
!
	INTEGER, PARAMETER :: NEXT=1000	!Number of points used to extend R-grid
	INTEGER NF				!Number of points in merged model
	REAL*8, ALLOCATABLE :: R_F(:)
	REAL*8, ALLOCATABLE :: T_F(:)
	REAL*8, ALLOCATABLE :: VW_F(:)
	REAL*8, ALLOCATABLE :: DSH_F(:)
	REAL*8, ALLOCATABLE ::dVdR_F(:)
	REAL*8, ALLOCATABLE :: TAUR_F(:)
	REAL*8, ALLOCATABLE :: TEMP_F(:)
!
! Grid for CMFGEN
!
	INTEGER ND_CMF
	INTEGER ND_SM
	INTEGER ND_BEL
	REAL*8, ALLOCATABLE :: R_CMF(:)
	REAL*8, ALLOCATABLE :: DENS_CMF(:)
	REAL*8, ALLOCATABLE :: V_CMF(:)
	REAL*8, ALLOCATABLE :: T_CMF(:)
	REAL*8, ALLOCATABLE :: dVdR_CMF(:)
	REAL*8, ALLOCATABLE :: TAUR_CMF(:)
!
	REAL*8 RSTAR			!Radius of Rstar in Rsun
	REAL*8 RMAX			!Outer radius of star
	REAL*8 MDOT			!Mass-loss in Msun/yr
	REAL*8 VINF			!Velocity (km/s)
	REAL*8 BETA			!Exponent for Beta-velocity law in outer wind.
	REAL*8 BETA_MIN			!Exponent for Beta-velocity law in inner wind.
	REAL*8 BETA_SCL			!Exponent for Beta-velocity law.
	REAL*8 VW_BEG                   !Velocity at wind/photosphere interface when defining CMFGEN grid.
!
	REAL*8 AMDOT
	REAL*8 RCORE
!
	REAL*8 MIN_VAL
	REAL*8 T1
	REAL*8 CONS,EPS
	REAL*8 DELR
	REAL*8 R0
	INTEGER IND0
	INTEGER C_INDX
	INTEGER NV_WIND 		!Index in *_F arrays at wind/photosphere interface 
!                                                  when defining CMFGEN grid.
!
	CHARACTER*80 FILENAME
	CHARACTER*80 FILENAME_CMF
	CHARACTER*80 STRING
!
	LOGICAL NEW_DEF_OPT
	LOGICAL USE_TAU
	LOGICAL FLAG
	LOGICAL REMOVE
	LOGICAL ANS
	LOGICAL PP_NOV
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: LU_OUT=8
	INTEGER, PARAMETER :: LU_SCR=9
	INTEGER, PARAMETER :: LU_V=9
	INTEGER, PARAMETER :: LU_PHOT=11
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	REAL*8 FUN_PI,PI
	EXTERNAL FUN_PI
!
	INTEGER I,ID,K,IOS
!
! Set up a velocity law connecting an hydrostatic
!  density structure in the photosphere and a beta-law wind
!
! Input: - infile : Tlusty model atmosphere
!        - rstar  : Stellar (core) radius in solar units
!        - rmax   : Maximum extension of the wind in stellar radii
!        - mdot   : Mass loss rate in solar mass per year
!        - vinf   : Wind terminal velocity in km/s
!        - beta   : Wind acceleration parameter
!        - ndw    : Number of depth points in the wind model
!
	WRITE(6,*)' '
	WRITE(6,*)' The default parameters for the model are read from RV_PARAMS'
	WRITE(6,*)' This file can be edited with a standard text editor'
	WRITE(6,*)' If file is not available, standard parameters are used.'
	WRITE(6,*)' '
!
! Set defaults
!
	FILENAME=' '
	PP_NOV=.FALSE.
	RSTAR=18.48D0
	RMAX=50.0D0
	MDOT=1.0D-06
	VINF=800.D0
	BETA=1.0D0
	BETA_MIN=BETA
	BETA_SCL=0.2D0
	ND_CMF=50
!
! Read in revised defaults if available.
!
	OPEN(UNIT=LU_IN,FILE='RV_PARAMS',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    CALL RD_OPTIONS_INTO_STORE(LU_IN,LU_SCR)
	    CALL RD_STORE_NCHAR(FILENAME,'FILE',80,L_TRUE,'File with photospheric model')
	    FILENAME=ADJUSTL(FILENAME)
	    FILENAME=FILENAME(1:INDEX(FILENAME,' '))
	    CALL RD_STORE_LOG(PP_NOV,'PP_NOV',L_TRUE,'Plane parallel model (no wind)?')
	    CALL RD_STORE_DBLE(RSTAR,'RSTAR',L_TRUE,'Radius of star')
	    IF(PP_NOV)THEN
	      MDOT=1.0D-22
	    ELSE
	      CALL RD_STORE_DBLE(RMAX,'RMAX',L_TRUE,'Radius of star')
	      CALL RD_STORE_DBLE(MDOT,'MDOT',L_TRUE,'Mass loss in Msun/yr')
	      CALL RD_STORE_DBLE(VINF,'VINF',L_TRUE,'Terminal velocity (km/s)')
	      CALL RD_STORE_DBLE(BETA,'BETA',L_TRUE,'Beta ')
	      CALL RD_STORE_DBLE(BETA_MIN,'BETA_MIN',L_TRUE,'Beta in photosphere')
	      CALL RD_STORE_DBLE(BETA_SCL,'BETA_SCL',L_TRUE,'Beta scale height')
	    END IF
	    CALL RD_STORE_INT(ND_CMF,'ND',L_TRUE,'Number of depths for new model')
	    CALL CLEAN_RD_STORE
	  END IF
	CLOSE(LU_IN)
	CLOSE(LU_SCR)
!
! Now read in parameters from terminal.
!
20	CALL GEN_IN(FILENAME,'Filename with photospheric model')
        OPEN(FILE=FILENAME,STATUS='OLD',ACTION='READ',UNIT=LU_PHOT,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error openening file with photospheric model'
	  GOTO 20
	END IF
	CALL GEN_IN(PP_NOV,'Plane parallel model (no wind)?')
	CALL GEN_IN(RSTAR,'RSTAR')
	IF(PP_NOV)THEN
	  MDOT=1.0D-22
	ELSE 
	  CALL GEN_IN(RMAX,'RMAX')
	  CALL GEN_IN(MDOT,'MDOT')
	  CALL GEN_IN(VINF,'VINF')
	  CALL GEN_IN(BETA,'BETA in outer wind')
	  CALL GEN_IN(BETA_MIN,'BETA in inner wind')
	  CALL GEN_IN(BETA_SCL,'BETA Scale height (in R*)')
	END IF
	CALL GEN_IN(ND_CMF,'Number of depth points')
!
	OPEN(UNIT=LU_OUT,FILE='REV_RV_PARAMS',STATUS='UNKNOWN')
	  WRITE(LU_OUT,*)TRIM(FILENAME),'  [FILE]'
	  WRITE(LU_OUT,*)PP_NOV,'  [PP_NOV]'
	  WRITE(LU_OUT,*)RSTAR,'  [RSTAR]'
	  WRITE(LU_OUT,*)RMAX,'  [RMAX]'
	  WRITE(LU_OUT,*)MDOT,'  [MDOT]'
	  WRITE(LU_OUT,*)VINF,'  [VINF]'
	  WRITE(LU_OUT,*)BETA,'  [BETA]'
	  WRITE(LU_OUT,*)BETA_MIN,'  [BETA_MIN]'
	  WRITE(LU_OUT,*)BETA_SCL,'  [BETA_SCL]'
	  WRITE(LU_OUT,*)ND_CMF,'  [ND]'
	CLOSE(UNIT=LU_OUT)
!
! Will use program units of 10^10 cm, V in km/s
!
	RCORE=RSTAR*6.9599D0
	AMDOT=MDOT*6.3029D0
	RMAX=RMAX*RCORE
	PI=FUN_PI()
!
	ND_MAX=100
	CALL GEN_IN(ND_MAX,'Maximum number of points in Photospheric model')
!
! index, column-mass, Tau Ross, Abs Ross, T, Ne, total density
!
	ALLOCATE (INDX(ND_MAX))
	ALLOCATE (DM(ND_MAX))
	ALLOCATE (TAUR(ND_MAX))
	ALLOCATE (AROSS(ND_MAX))
	ALLOCATE (T(ND_MAX))
	ALLOCATE (ED(ND_MAX))
	ALLOCATE (DSH(ND_MAX))
	ALLOCATE (KAPPA(ND_MAX))
!
	STRING(1:1)='!'
	DO WHILE(STRING(1:1) .EQ. '!')
	  READ(LU_PHOT,'(A)')STRING
	END DO
	BACKSPACE(LU_PHOT)
	ANS=.FALSE.
	DO WHILE(.NOT. ANS)
	  READ(LU_PHOT,'(A)')STRING
	  WRITE(6,'(A)')TRIM(STRING)
	  ANS=.TRUE.
	  CALL GEN_IN(ANS,'Is the above the first READABLE record')
	END DO
	READ(STRING,*)INDX(1),DM(1),TAUR(1),AROSS(1),T(1),ED(1),DSH(1)
!	
	I=1
	DO WHILE(TAUR(I) .LT. 100)
	  READ(LU_PHOT,*,END=250)INDX(I+1),DM(I+1),TAUR(I+1),AROSS(I+1),
	1                   T(I+1),ED(I+1),DSH(I+1)
	  I=I+1
	END DO
250	CONTINUE
	CLOSE(LU_PHOT)
	ND=I
!
	WRITE(6,*)' '
	WRITE(6,*)' Typically we extend model atmosphere to TAU=100'
	WRITE(6,*)' This option allows you to omit the depth with Tau > 100'
	WRITE(6,*)' '
	WRITE(6,*)'Tau(ND-1)=',TAUR(ND-1)
	WRITE(6,*)'Tau(ND)=',TAUR(ND)
	REMOVE=.TRUE.; CALL GEN_IN(REMOVE,'Remove last tau depth')
	IF(REMOVE)ND=ND-1
	T(1:ND)=1.0D-04*T(1:ND)
!
	KAPPA(1)=TAUR(1)/DM(1)
	DO I=2,ND
	  KAPPA(I)=(TAUR(I)-TAUR(I-1))/(DM(I)-DM(I-1))
	END DO
	WRITE(6,'(5ES14.4)')(KAPPA(I),I=1,ND)
	KAPPA(1)=KAPPA(2)
!
	ALLOCATE (ZZ(ND))
	ALLOCATE (RD(ND))
	ALLOCATE (RD_NORM(ND))
	ALLOCATE (VPH(ND))
	ALLOCATE (dVdR_PH(ND))
!
! Radial scale and photospheric velocity. NB. ZZ will be in units of 10^10 cm.
!
	ZZ(ND)=0.0D0
	DO ID=ND-1,1,-1
	  ZZ(ID)=ZZ(ID+1)+2.0D-10*(DM(ID+1)-DM(ID))/(DSH(ID+1)+DSH(ID))
	END DO
!
	RD(1:ND)=ZZ(1:ND)+RCORE
	VPH(1:ND)=AMDOT/4/PI/RD(1:ND)/RD(1:ND)/DSH(1:ND)
!
!
! Compute dVdR_PH
!
	dVdR_PH(1)=LOG(VPH(1)/VPH(2))/LOG(RD(1)/RD(2))
	DO I=2,ND-1
	    dVdR_PH(I)=LOG(VPH(I-1)/VPH(I+1))/LOG(RD(I-1)/RD(I+1))
	END DO
	dVdR_PH(ND)=LOG(VPH(ND-1)/VPH(ND))/LOG(RD(ND-1)/RD(ND))
!
	OPEN(UNIT=22,FILE='OLD_ATM',STATUS='REPLACE')
	  DO I=1,ND
	    WRITE(22,'(7ES14.6)')RD(I),T(I),ED(I),VPH(I),TAUR(I),DSH(I),DM(I)
	  END DO
	CLOSE(UNIT=22)
!
	CALL GEN_ASCI_OPEN(LU_V,'RV_OLD_ATM','REPLACE',' ',' ',IZERO,IOS)
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(A)')'!'
	  WRITE(LU_V,'(A,A)')'! TLUSTY data read from file: ',TRIM(FILENAME)
	  WRITE(LU_V,'(A)')'!'
	  WRITE(LU_V,'(A)')'!'
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(1X,I4,2X,A)')ND,'!Number of depth points'
	  DO I=1,ND
	    WRITE(LU_V,'(F18.8,2ES14.6,F12.4,ES12.2)')RD(I),VPH(I),
	1             dVdR_PH(I)-1.0D0,T(I),TAUR(I)
	  END DO
	CLOSE(LU_V)
!
	WRITE(6,*)' '
	WRITE(6,*)' The graph plotting is primarily for diagnostic purposes'
	WRITE(6,*)' Enter /null, and thein e to skip each graph'
	WRITE(6,*)' Plotting hydrostatic V versus R/R*'
	WRITE(6,*)' '
	RD_NORM=RD/RCORE
	CALL DP_CURVE(ND,RD_NORM,VPH)
	CALL GRAMON_PGPLOT(' ',' ',' ',' ')
!
	IF(PP_NOV)THEN
!
	  ALLOCATE (TB(ND))
	  ALLOCATE (TA(ND_CMF))
	  ALLOCATE (DENS_CMF(ND_CMF))
	  ALLOCATE (R_CMF(ND_CMF))
	  ALLOCATE (V_CMF(ND_CMF))
	  ALLOCATE (dVdR_CMF(ND_CMF))
	  ALLOCATE (T_CMF(ND_CMF))
	  ALLOCATE (TAUR_CMF(ND_CMF))
!
	  DO I=1,ND
	    TB(I)=I
	  END DO
	  DO I=1,ND_CMF
	    TA(I)=1.0D0+ND*(I-1.0D0)/ND_CMF
	  END DO
	  TA(1)=1.0D0; TA(ND_CMF)=ND
!
	  CALL MON_INTERP(R_CMF,ND_CMF,IONE,TA,ND_CMF,RD,ND,TB,ND)
	  VPH(1:ND)=LOG(VPH(1:ND))
	  CALL MON_INTERP(V_CMF,ND_CMF,IONE,R_CMF,ND_CMF,VPH,ND,RD,ND)
	  V_CMF(1:ND_CMF)=EXP(V_CMF(1:ND_CMF))
	  CALL MON_INTERP(T_CMF,ND_CMF,IONE,R_CMF,ND_CMF,T,ND,RD,ND)
	  CALL MON_INTERP(TAUR_CMF,ND_CMF,IONE,R_CMF,ND_CMF,TAUR,ND,RD,ND)
!
! Compute dVdR_CMF
!
	  dVdR_CMF(1)=LOG(V_CMF(1)/V_CMF(2))/LOG(RD(1)/RD(2))
	  DO I=2,ND-1
	    dVdR_CMF(I)=LOG(V_CMF(I-1)/V_CMF(I+1))/LOG(RD(I-1)/RD(I+1))
	  END DO
	  dVdR_CMF(ND)=LOG(V_CMF(ND-1)/V_CMF(ND))/LOG(RD(ND-1)/RD(ND))
!
	  FILENAME_CMF='RVSIG_COL'
	  CALL GEN_IN(FILENAME_CMF,'File for CMFGEN results')
	  CALL GEN_ASCI_OPEN(LU_V,FILENAME_CMF,'REPLACE',' ',' ',IZERO,IOS)
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(A)')'!'
	  WRITE(LU_V,'(A,A)')'! TLUSTY data read from file: ',TRIM(FILENAME)
	  WRITE(LU_V,'(A)')'! Output for plane parallel model'
	  WRITE(LU_V,'(A)')'!'
	  WRITE(LU_V,'(A,ES10.4,A)')'!          R*=',RSTAR,' Rsun'
	  WRITE(LU_V,'(A,ES10.4,A)')'!        Mdot=',Mdot,' Msun/yr'
	  WRITE(LU_V,'(A)')'!'
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	  WRITE(LU_V,'(1X,I4,2X,A)')ND_CMF,'!Number of depth points'
	  DO I=1,ND_CMF
	    WRITE(LU_V,'(F18.8,2ES14.6,F12.4,ES12.2)')R_CMF(I),V_CMF(I),
	1             dVdR_CMF(I)-1.0D0,T_CMF(I),TAUR_CMF(I)
	  END DO
	  CLOSE(LU_V)
	  WRITE(6,*)'Model structure for CMFGEN written to ',TRIM(FILENAME_CMF)
	  STOP
	END IF
!
	ALLOCATE (RA(NBIG))
	ALLOCATE (TA(NBIG))
	ALLOCATE (RA_NORM(NBIG))
	ALLOCATE (VPA(NBIG))
	ALLOCATE (VW(NBIG))
	ALLOCATE (DSHA(NBIG))
	ALLOCATE (KAPPA_A(NBIG))
	ALLOCATE (TAUR_A(NBIG))
!
	T1=(RD(1)-RD(ND))/(NBIG-1)
	RA(1)=RD(1)
	DO I=2,NBIG-1
	  RA(I)=RA(1)-(I-1)*T1
	END DO
	RA(NBIG)=RD(ND)
	CALL MON_INTERP(VPA,NBIG,IONE,RA,NBIG,VPH,ND,RD,ND)
	CALL MON_INTERP(DSHA,NBIG,IONE,RA,NBIG,DSH,ND,RD,ND)
	CALL MON_INTERP(TA,NBIG,IONE,RA,NBIG,T,ND,RD,ND)
	CALL MON_INTERP(TAUR_A,NBIG,IONE,RA,NBIG,TAUR,ND,RD,ND)
	CALL MON_INTERP(KAPPA_A,NBIG,IONE,RA,NBIG,KAPPA,ND,RD,ND)
!
! Search connecting point with beta-law wind. To find the connection point we
! adopt a velocity law of the form:
! 
!                                  v(r)=vinf(1-R0/r)^beta
!
! For R0  equally RCORE v(r) will be greater than VPA (V from TLUSTY) in the 
! inner regions. We then vary R0 until we v(r) is always less the VPA (just).
! This defines R0. The connection point is defined by VPA-V(r)=0 (
! we take it to be a minimum).
!
	DO I=NBIG-1,1,-1
	  R0=RA(I)
	  FLAG=.FALSE.  
	  DO K=I-1,1,-1
	    VW(K)=VINF*(1.0D0-R0/RA(K))**(BETA+(BETA_MIN-BETA)*
	1          EXP( (1.0D0-RA(K)/RCORE)/BETA_SCL ))
	    IF(VW(K) .GT. VPA(K))FLAG=.TRUE.
	  END DO
	  IF(.NOT. FLAG)EXIT
	END DO
	IND0=I  
!
! We have found R0, now need to find the connection radius.
!
	MIN_VAL=1.0E+37
	C_INDX=0
	DO K=IND0-1,1,-1
	  T1=VPA(K)-VW(K)
	  IF(T1 .LT. MIN_VAL)THEN
	    C_INDX=K
	    MIN_VAL=T1
	  END IF
	END DO
!
	WRITE(6,*)' '
	WRITE(6,*)' '
	WRITE(6,*)'                    Core depth is',I
	WRITE(6,*)'                         R0/RCORE=',R0/RCORE
	WRITE(6,*)'               Connection index is',C_INDX
	WRITE(6,*)'     Radius at connection depth is',RA(C_INDX)/RCORE
	WRITE(6,*)'   Velocity at connection depth is',VPA(C_INDX),' km/s'
	WRITE(6,*)' Difference in match velocities is',MIN_VAL
	WRITE(6,*)'                         /\V / V =',MIN_VAL/VPA(C_INDX)
	WRITE(6,*)'        Tau at connection depth is',TAUR_A(C_INDX)
	WRITE(6,*)' '
	WRITE(6,*)' '
!
	VW(I+1:NBIG)=VPA(I+1:NBIG)
!
	WRITE(6,*)' '
	WRITE(6,*)'Plotting VPH and VBETA versus R/R*'
	WRITE(6,*)' '
	RD_NORM=RD/RCORE
	RA_NORM=RA/RCORE
	CALL DP_CURVE(NBIG,RA_NORM,VPA)
	CALL DP_CURVE(NBIG,RA_NORM,VW)
!
! Extend velocity law on hydrostatic grid to illustrate the connection.
!
	DO K=C_INDX,1,-1
	  VW(K)=VINF*(1.0D0-R0/RA(K))**(BETA+(BETA_MIN-BETA)*
	1          EXP( (1.0D0-RA(K)/RCORE)/BETA_SCL ))
	END DO
	VW(C_INDX+1:NBIG)=VPA(C_INDX+1:NBIG)
	CALL DP_CURVE(NBIG,RA_NORM,VW)
	CALL GRAMON_PGPLOT(' ',' ',' ',' ')
!
! Now need to extend the Radius grid to account for the wind.
!
	IF(RMAX .GT. RA(1))THEN
!
	  NF=NEXT+(NBIG-C_INDX)
	  ALLOCATE (R_F(NF))
	  ALLOCATE (VW_F(NF))
	  ALLOCATE (DSH_F(NF))
	  ALLOCATE (T_F(NF))
	  ALLOCATE (TAUR_F(NF))
	  ALLOCATE (TEMP_F(NF))
!
	  DELR=LOG( (RMAX-RCORE)/ (RA(C_INDX)-RCORE) )/(NEXT-1)
	  R_F(1)=RMAX
	  DO I=2,NEXT
	    R_F(I)=RCORE+(R_F(1)-RCORE)/EXP(DELR*(I-1))
	  END DO
	  R_F(NEXT+1:NF)=RA(C_INDX+1:NBIG)
!
	  DO K=1,NEXT
	    VW_F(K)=VINF*(1.0D0-R0/R_F(K))**(BETA+(BETA_MIN-BETA)*
	1          EXP( (1.0D0-R_F(K)/RCORE)/BETA_SCL ))
	    DSH_F(K)=AMDOT/4/PI/R_F(K)/R_F(K)/VW_F(K)
	  END DO
	  VW_F(NEXT+1:NF)=VPA(C_INDX+1:NBIG)
	  DSH_F(NEXT+1:NF)=DSHA(C_INDX+1:NBIG)
	  T_F(NEXT+1:NF)=TA(C_INDX+1:NBIG)
	  T1=T_F(NEXT+1)
	  T_F(1:NEXT)=T1
!
	  WRITE(6,*)'KAPPA(1)=',KAPPA_A(1)
	  WRITE(6,*)'KAPPA(C_INDX)=',KAPPA_A(C_INDX)
	  TAUR_F(1)=KAPPA_A(1)*R_F(1)*1.0D+10*DSH_F(1)
	  DO I=2,NEXT
	    TAUR_F(I)=TAUR_F(I-1)+0.5D+10*(R_F(I-1)-R_F(I))*
	1              (DSH_F(I)+DSH_F(I+1))*KAPPA_A(1)
	  END DO
	  T1=TAUR_F(NEXT)-TAUR_A(C_INDX)
	  DO I=NEXT+1,NF
	    TAUR_F(I)=TAUR_A(C_INDX+(I-NEXT-1))+T1
	  END DO
!
	  WRITE(6,*)' '
	  WRITE(6,*)' Plotting V and density as a function of R on the fine grid'
	  WRITE(6,*)' '
	  R_F=R_F/RCORE
	  DSH_F=LOG10(DSH_F)
	  CALL DP_CURVE(NF,R_F,VW_F)
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  CALL DP_CURVE(NF,R_F,DSH_F)
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  R_F=R_F*RCORE
!
	END IF
!
! No we define the radius grid for CMGFEN, equally spaced in Log(Den)
!
	ALLOCATE (DENS_CMF(ND_CMF))
	ALLOCATE (R_CMF(ND_CMF))
	ALLOCATE (V_CMF(ND_CMF))
	ALLOCATE (dVdR_CMF(ND_CMF))
	ALLOCATE (T_CMF(ND_CMF))
	ALLOCATE (TAUR_CMF(ND_CMF))
!
	WRITE(6,*)' '
	WRITE(6,*)' You need to decide how many depth points will cover the wind'
	WRITE(6,*)' In the wind, the points are equally spaced in long density'
	WRITE(6,*)' In the photosphere, the points are equally spaced in long tau'
	WRITE(6,*)' For weak winds, you need at least 20 to 30 points in the wind.'
	WRITE(6,*)' NB: The velocity at the wind/photosphere interface need not be the'//
	1                         ' same as the connection velocity'
	WRITE(6,*)' Check that you don''t have large velocity changes near the sonic point'
	WRITE(6,*)' In the photosphere, you ideally should have at least 4 (pref 5 or'//
	1                         ' more) points per decade in log tau.'
	WRITE(6,*)' For strong dense winds, the number of points in the photosphere'//
	1                         ' can be small.'
	WRITE(6,*)' '
!
	USE_TAU=.FALSE.
	NEW_DEF_OPT=.TRUE.
	CALL GEN_IN(NEW_DEF_OPT,'Use new default option to select r grid?')
	IF(NEW_DEF_OPT)THEN
	  ND_SM=ND_CMF/2
	  CALL GEN_IN(ND_SM,'Number of depth points above connection point?')
	  VW_BEG=MIN(1.0D0,VPA(C_INDX))
	  CALL GEN_IN(VW_BEG,'Velocity at photosphere/wind interface')
!
! Find the wind/photosphere interface.
!
	  NV_WIND=1
	  DO WHILE(VW_F(NV_WIND) .GT. VW_BEG)
	    NV_WIND=NV_WIND+1
	  END DO
!
	  TEMP_F(1:NF)=DLOG(TAUR_F(1:NF))
	  DO I=1,NF
	    WRITE(45,*)I,DSH_F(I),TAUR_F(I)
	  END DO
!
! Determine and create the grid in the photosphere.
!
	  ND_BEL=ND_CMF-ND_SM
	  T1=(TEMP_F(NF)-TEMP_F(NV_WIND))/(ND_BEL-1)
	  DO I=ND_SM+1,ND_CMF-3
	    TAUR_CMF(I)=TEMP_F(NV_WIND)+(I-ND_SM)*T1
	  END DO
	  TAUR_CMF(ND_CMF)=TEMP_F(NF)
	  TAUR_CMF(ND_CMF-1)=TEMP_F(NF)-0.1*T1
	  TAUR_CMF(ND_CMF-2)=TEMP_F(NF)-0.3*T1
!
	  DO I=ND_SM+1,ND_CMF
	    WRITE(45,*)I,TAUR_CMF(I)
	  END DO
!
	  CALL MON_INTERP(DENS_CMF(ND_SM+1:),ND_BEL,IONE,TAUR_CMF(ND_SM+1:),ND_BEL,
	1                DSH_F,NF,TEMP_F,NF)
!
	  DO I=ND_SM+1,ND_CMF
	    WRITE(45,*)I,TAUR_CMF(I),DENS_CMF(I)
	  END DO
	  TAUR_CMF(ND_SM+1:ND_CMF)=EXP(TAUR_CMF(ND_SM+1:ND_CMF))
!
! Determine and create the grid in the wind.
!
	  WRITE(6,*)'Now choosing points in wind region'
	  T1=(DSH_F(NV_WIND)-DSH_F(1))/(ND_SM-3)
	  DENS_CMF(1)=DSH_F(1)
	  DENS_CMF(2)=DSH_F(1)+0.1*T1
	  DENS_CMF(3)=DSH_F(1)+0.3*T1
	  DO I=4,ND_SM
	    DENS_CMF(I)=DENS_CMF(1)+(I-3)*T1
	  END DO
	  CALL MON_INTERP(TAUR_CMF,ND_SM,IONE,DENS_CMF,ND_CMF,TAUR_F,NF,DSH_F,NF)
!
	ELSE
          WRITE(6,*)' '
	  WRITE(6,*)'You are now using the old slection method which was',
	1                  ' maintained for compatibility.'
          WRITE(6,*)'Please check the final grid carefully. '
          WRITE(6,*)' '
	  USE_TAU=.TRUE.
	  CALL GEN_IN(USE_TAU,'Euqally spaced in Log Tau (T) or density (F)?')
	  IF(USE_TAU)THEN
	    WRITE(6,*)' In the following section we choose our r grid ',
	1                 'equally spaced in Log Tau'
	    WRITE(6,*)' To provide greater flexibility, we can modify the tau',
	1                ' scale by the velocity'
	    WRITE(6,*)'           i.e., Tau(Mod)=Tau/ (C+V)**EPS'
	    WRITE(6,*)' You need to check that the chosen r grid is satsifactory.'
	    WRITE(6,*)' '
	    CONS=0.1D0
	    EPS=0.75D00
	    CALL GEN_IN(CONS,'Constant for Velocity scaling of Tau (>0) ')
	    CALL GEN_IN(EPS,'Exponent for Velocity scaling of Tau')
	    TEMP_F(1:NF)=DLOG( TAUR_F(1:NF)/(CONS+VW_F(1:NF))**EPS )
	    T1=(TEMP_F(NF)-TEMP_F(1))/(ND_CMF-5)
	    TAUR_CMF(1)=TEMP_F(1)
	    TAUR_CMF(2)=TEMP_F(1)+0.1*T1
	    TAUR_CMF(3)=TEMP_F(1)+0.3*T1
	    DO I=4,ND_CMF-3
	      TAUR_CMF(I)=TAUR_CMF(1)+(I-3)*T1
	    END DO
	    TAUR_CMF(ND_CMF)=TEMP_F(NF)
	    TAUR_CMF(ND_CMF-1)=TEMP_F(NF)-0.1*T1
	    TAUR_CMF(ND_CMF-2)=TEMP_F(NF)-0.3*T1
	    CALL MON_INTERP(DENS_CMF,ND_CMF,IONE,TAUR_CMF,ND_CMF,
	1                  DSH_F,NF,TEMP_F,NF)
	    TAUR_CMF(1:ND_CMF)=EXP(TAUR_CMF(1:ND_CMF))
	  ELSE
	    T1=(DSH_F(NF)-DSH_F(1))/(ND_CMF-5)
	    DENS_CMF(1)=DSH_F(1)
	    DENS_CMF(2)=DSH_F(1)+0.1*T1
	    DENS_CMF(3)=DSH_F(1)+0.3*T1
	    DO I=4,ND_CMF-3
	      DENS_CMF(I)=DENS_CMF(1)+(I-3)*T1
	    END DO
	    DENS_CMF(ND_CMF)=DSH_F(NF)
	    DENS_CMF(ND_CMF-1)=DSH_F(NF)-0.1*T1
	    DENS_CMF(ND_CMF-2)=DSH_F(NF)-0.3*T1
	    CALL MON_INTERP(TAUR_CMF,ND_CMF,IONE,DENS_CMF,ND_CMF,TAUR_F,NF,DSH_F,NF)
	  END IF
	END IF
!
! Compute dVdR
!
	ALLOCATE (dVdR_F(NF))
	dVdR_F(1)=LOG(VW_F(1)/VW_F(2))/LOG(R_F(1)/R_F(2))
	DO I=2,NF-1
	    dVdR_F(I)=LOG(VW_F(I-1)/VW_F(I+1))/LOG(R_F(I-1)/R_F(I+1))
	END DO
	dVdR_F(NF)=LOG(VW_F(NF-1)/VW_F(NF))/LOG(R_F(NF-1)/R_F(NF))
!
! Now obtain the R grid.
!
	CALL MON_INTERP(R_CMF,ND_CMF,IONE,DENS_CMF,ND_CMF,R_F,NF,DSH_F,NF)
	CALL MON_INTERP(V_CMF,ND_CMF,IONE,DENS_CMF,ND_CMF,VW_F,NF,DSH_F,NF)
	CALL MON_INTERP(dVdR_CMF,ND_CMF,IONE,DENS_CMF,ND_CMF,dVdR_F,NF,DSH_F,NF)
	CALL MON_INTERP(T_CMF,ND_CMF,IONE,DENS_CMF,ND_CMF,T_F,NF,DSH_F,NF)
!
! Correct computed TAU grid for velocity scaling.
!
	IF(USE_TAU)THEN
	  TAUR_CMF(1:ND_CMF)=TAUR_CMF(1:ND)*(CONS+V_CMF(1:ND))**EPS
	END IF
!
! Output results for CMFGEN
!
	FILENAME_CMF='RVSIG_COL'
	CALL GEN_IN(FILENAME_CMF,'File for CMFGEN results')
	CALL GEN_ASCI_OPEN(LU_V,FILENAME_CMF,'REPLACE',' ',' ',IZERO,IOS)
	WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	WRITE(LU_V,'(A)')'!'
	WRITE(LU_V,'(A,A)')'! TLUSTY data read from file: ',TRIM(FILENAME)
	WRITE(LU_V,'(A)')'!'
	WRITE(LU_V,'(A,ES10.4,A)')'!          R*=',RSTAR,' Rsun'
	WRITE(LU_V,'(A,ES10.4,A)')'!        Mdot=',Mdot,' Msun/yr'
	WRITE(LU_V,'(A,ES10.4,A)')'!        Vinf=',VINF,' km/s'
	WRITE(LU_V,'(A,ES10.4,A)')'!        Beta=',BETA
	IF(BETA .NE. BETA_MIN)THEN
	  WRITE(LU_V,'(A,ES10.4,A)')'!    Beta_min=',BETA_MIN
	  WRITE(LU_V,'(A,ES10.4,A)')'!    Beta_scl=',BETA_SCL,
	1                 ' (Scale height for beta in R*)'
	END IF
	WRITE(LU_V,'(A)')'!'
	WRITE(LU_V,'(A,1X,2ES12.4)')'!      R0 for Beta-velocity law is',R0,R0/RCORE
	WRITE(LU_V,'(A,1X,2ES12.4)')'!             Connection radius is',
	1                           RA(C_INDX),RA(C_INDX)/RCORE
	WRITE(LU_V,'(A,1X, ES12.4,A)')'! Velocity at connection radius is',
	1                          VPA(C_INDX),' km/s'
	WRITE(LU_V,'(A,1X, ES12.4)')'!       Tau at connection depth is',
	1                          TAUR_A(C_INDX)
	WRITE(LU_V,'(A)')'!'
	IF(NEW_DEF_OPT)THEN
	  WRITE(LU_V,'(A)')'! Chosen r grid equally spaced in Log Density in wind.'
	  WRITE(LU_V,'(A)')'! Chosen r grid equally spaced in Log Tau in photosphere.'
	  WRITE(LU_V,'(A,I4)')'! Number of points in wind is:',ND_SM
	  WRITE(LU_V,'(A,ES9.2,A)')'! Wind grid starts at:',VW_BEG,' km/s'
	ELSE IF(USE_TAU)THEN
	  WRITE(LU_V,'(A)')'! Chosen r grid equally spaced in Log Tau.'
	  WRITE(LU_V,'(A)')'! Tau scale modified by velocity with:'
	  WRITE(LU_V,'(A,F5.2,3X,A,F5.2)')'! CONS=',CONS,'EPS=',EPS
	ELSE
         WRITE(LU_V,'(A)')'! Chosen r grid equally spaced in Log Density.'
	END IF
	WRITE(LU_V,'(A)')'!'
	WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	WRITE(LU_V,'(71A)')'!',('*',I=1,70)
	WRITE(LU_V,'(1X,I4,2X,A)')ND_CMF,'!Number of depth points'
	DO I=1,ND_CMF
	  WRITE(LU_V,'(F18.8,2ES14.6,F12.4,ES12.2)')R_CMF(I),V_CMF(I),
	1             dVdR_CMF(I)-1.0D0,T_CMF(I),TAUR_CMF(I)
	END DO
!
	WRITE(6,*)' '
	WRITE(6,*)' Plotting V as a function of R on the final grid'
	WRITE(6,*)' '
	R_CMF(1:ND_CMF)=R_CMF(1:ND_CMF)/RCORE
	V_CMF(1:ND_CMF)=V_CMF(1:ND_CMF)/1.0D+05
	CALL DP_CURVE(ND_CMF,R_CMF,V_CMF)
	CALL GRAMON_PGPLOT(' ',' ',' ',' ')
!
	STOP
	END
