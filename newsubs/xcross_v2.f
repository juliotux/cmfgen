!
! Routine returns "core" photoionization cross section for a wide range 
! of ions. Data is from Verner and Yakovlev 1995(?).
!
! XCROSS   : Returned cross section in units of 10^{-10} cm^2
! FREQ     : 10^15 Hz
! ZCORE    : Charge of nucleus.
! NUM_ELEC : Number of electrons in system that is being photoionized.
!
! PQN      : Principal quantum number of state. If zero, it (and ANG) are
!              ignored, and the cross-section returned is a sum over several 
!              initial states.
! ANG      : Angular momentum principal quantum number.
! DO_ALL   : If true, the total photoionization cross section, for the
!              ion under consideration, is returned. Here total refers
!              to the summation over all initial states tabulated. It may
!              include the cross-section of the valence electron.
!            When DO_ALL is false, only the cross-section for non-valence
!              electrons is returned.
!            Non-valence electrons are defined as:
!              If(Ne < 3), there are no valence electrons. This would apply
!                 to H, and HeI like atomic systems.
!              If(Ne = 3), then electrons in (n:l)=(1:0) are treated as 
!                valence electrons provided D_LIT is true.
!              If(Ne < 12), (n:l)=(1:0)
!              If(Ne < 20), (n:l)=(1,2 : 0,1)
!              If(Ne >= 20), (n:l)=(1,2 : 0,1)
!            We need to be careful that X-ray cross-sections are not
!              doubly counted.
! DO_LIT   : Returns the cross-section for the K-shell of Lithium like ions.
!
	FUNCTION XCROSS_V2(FREQ,ZCORE,NUM_ELEC,PQN,ANG,DO_ALL,DO_LIT)
	USE XRAY_DATA_MOD
	IMPLICIT NONE
!
	REAL*8 XCROSS_V2
	REAL*8 FREQ
	REAL*8 ZCORE
	REAL*8 NUM_ELEC
	INTEGER PQN
	INTEGER ANG
!
	LOGICAL DO_ALL
	LOGICAL DO_LIT
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER I,J
	INTEGER IZ
	INTEGER NE
	INTEGER LOOP_PQN_MAX
	REAL*8 CON_FAC
	REAL*8 NU_EV
	REAL*8 Y
	REAL*8 Q
!
	XCROSS_V2=0.0D0
	CON_FAC=1.0D0/0.24191D0			!10^15 Hz to ev
	NU_EV=FREQ*CON_FAC
!
! Define for ease of typing.
!
	IZ=NINT(ZCORE)
	NE=NINT(NUM_ELEC)
!
	IF(.NOT. XRAY_PHOT_RD_IN)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in XCROSS_V2'
	  WRITE(LUER,*)'X-ray photoioinzation cross-sections have not'//
	1              ' been read in'
	  STOP
	END IF
!
	IF(IZ .GT. ATNO_MAX_X)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in XCROSS_V2'
	  WRITE(LUER,*)'No X-ray data available for species'//
	1                     ' with Z > ',ATNO_MAX_X
	  STOP
	END IF
!
! This section allows the photoionization cross-ection for a particular
! (n,l) state to be returned. Not that SIG_0_X is zero for (n,l) in array,
! but outside valid bounds.
!
	IF(PQN .NE. 0)THEN
	  J=ANG
	  I=PQN
	  IF(I .GT. SQRT(IZ/2.0D0)+2 .OR. ANG .LT. 0 .OR. ANG .GE. PQN)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in XCROSS_V2'
	    WRITE(LUER,*)'Invalid PQN or ANG'
	    WRITE(LUER,*)'n=',PQN,' ANG=',ANG
	    WRITE(LUER,*)'Species has z=',ZCORE
	    WRITE(LUER,*)'Species has # elec. =',NUM_ELEC
	    STOP
	  END IF
	  IF(I .GT. PQN_MAX_X)RETURN
	  IF(J .GT. ANG_MAX_X)RETURN
	  IF( SIG_0_X(IZ,NE,I,J) .NE. 0 .AND. 
	1         NU_EV .GE. E_THRESH_X(IZ,NE,I,J) )THEN
	    Y=NU_EV/E_0_X(IZ,NE,I,J) 
	    Q=5.5D0+J-0.5D0*P_X(IZ,NE,I,J)
	    XCROSS_V2=XCROSS_V2+
	1   SIG_0_X(IZ,NE,I,J) *
	1       ((Y-1.0D0)**2 +Y_W_X(IZ,NE,I,J)**2)/ (Y**Q) /
	1       (1.0D0+SQRT(Y/Y_A_X(IZ,NE,I,J)) )**P_X(IZ,NE,I,J)
	  END IF
!
! Return the the photoionzation cors-section, summed over all shells.
! This may include the cross-section for the valence electron.
!
	ELSE IF(DO_ALL)THEN
	  DO J=0,ANG_MAX_X
	    DO I=1,PQN_MAX_X
	      IF( SIG_0_X(IZ,NE,I,J) .NE. 0 .AND. 
	1         NU_EV .GE. E_THRESH_X(IZ,NE,I,J) )THEN
	        Y=NU_EV/E_0_X(IZ,NE,I,J) 
	        Q=5.5D0+J-0.5D0*P_X(IZ,NE,I,J)
	        XCROSS_V2=XCROSS_V2+
	1       SIG_0_X(IZ,NE,I,J) *
	1          ((Y-1.0D0)**2 +Y_W_X(IZ,NE,I,J)**2)/ (Y**Q) /
	1          (1.0D0+SQRT(Y/Y_A_X(IZ,NE,I,J)) )**P_X(IZ,NE,I,J)
	      END IF
	    END DO
	  END DO
	ELSE
! 
! We only wish to return photoionization cross-sections for the core
! (i.e., non-valence) electrons. We do not treat H, HeI and He2 isoelectronic 
! sequence ions. If DO_LIT is false, return 0 for Li I isoelectronic sequence 
! ions.
!
! For ions with 3s2 3p6 3dn electron configurations, we DO NOT include 
! ionizations from the 3s and 3p states. The following code may need to be 
! changed when elements with z > 10 are included.
!
! We only include ionizations from the 2p shell if there is 3 or more
! electrons in the n=3 shell.
!                        
	  IF( NE .LE. 2 .OR. (.NOT. DO_LIT .AND. NE .EQ. 3) )RETURN
!
	  DO J=0,ANG_MAX_X
	    LOOP_PQN_MAX=MIN(3,PQN_MAX_X)
	    IF(NE .LE. 19)LOOP_PQN_MAX=MIN(2,PQN_MAX_X)
	    IF(NE .LE. 11)LOOP_PQN_MAX=1
	    DO I=1,LOOP_PQN_MAX
	      IF( SIG_0_X(IZ,NE,I,J) .NE. 0 .AND. 
	1         N_ED_EJ(IZ,NE,I,J) .GT. 1 .AND.
	1         NU_EV .GE. E_THRESH_X(IZ,NE,I,J) )THEN
	        Y=NU_EV/E_0_X(IZ,NE,I,J) 
	        Q=5.5D0+J-0.5D0*P_X(IZ,NE,I,J)
	        XCROSS_V2=XCROSS_V2+
	1       SIG_0_X(IZ,NE,I,J) *
	1          ((Y-1.0D0)**2 +Y_W_X(IZ,NE,I,J)**2)/ (Y**Q) /
	1          (1.0D0+SQRT(Y/Y_A_X(IZ,NE,I,J)) )**P_X(IZ,NE,I,J)
	      END IF
	    END DO
	  END DO
	END IF
!
	XCROSS_V2=XCROSS_V2*1.0D-08
!
	RETURN
	END
