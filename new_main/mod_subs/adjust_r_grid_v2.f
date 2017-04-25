!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
	SUBROUTINE ADJUST_R_GRID_V2(POPS,P,
	1             FLUXMEAN,ESEC,GRID_TYPE,RG_PARS,N_PARS,
	1             ND,NT,NC,NP)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered : 10-Jan-2006 Added extra point at inner boundary.
!                         Boundary problems, arising because of rounding
!                         error were fixed.
! Altered : 02-May-2004 Now only change, R, V, SIGMA and POPS.
! Altered : 18-Mar-2004 FLUXMEAN passed insted of dTAU_OLD
!                       ESEC passed in call.
!                       As still developing, subroutine name NOT changed. 
! Created : 25-Feb-2004
!
	INTEGER ND,NT,NC,NP
!
	REAL*8 POPS(NT,ND)
!
! NB: dTAU_OLD(I) = Optical depth increment between depth I & I+1.
!                   It is NOT the optical depth scale.
!                   Generally computed using the FLUX_MEAN opacity.
!
	LOGICAL TRAPFORJ
	REAL*8 P(NP)
	REAL*8 JQW(ND,NP)
	REAL*8 HQW(ND,NP)
	REAL*8 KQW(ND,NP)
	REAL*8 HMIDQW(ND,NP)
	REAL*8 NMIDQW(ND,NP)
	REAL*8 MU_AT_RMAX(NP)
	REAL*8 HQW_AT_RMAX(NP)
!
	REAL*8 FLUXMEAN(ND)
	REAL*8 ESEC(ND)
!
	INTEGER N_PARS
	CHARACTER(LEN=*) GRID_TYPE
	REAL*8 RG_PARS(N_PARS)
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	EXTERNAL JWEIGHT,HWEIGHT,KWEIGHT,NWEIGHT
        EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
!
! Local variables.
!
	REAL*8 dTAU_OLD(ND)
	REAL*8 R_OLD(ND)
	REAL*8 LOG_R_OLD(ND)
	REAL*8 LOG_R(ND)
	REAL*8 TAU_OLD(ND)
	REAL*8 TAU(ND)
	REAL*8 TA(ND)
	REAL*8 TB(ND)
!
! The fine grid (FG) is chosen to cover the ioization front. The default valuse are
! -3.0 to 1.5D0 in log(TAU) space.
!
	REAL*8 FG_MIN                   !Default=-3.0D0
	REAL*8 FG_MAX			!Default=1.5D0
	REAL*8 FG_RANGE			!Default=4.5D0
!
	REAL*8 W1(NP),W2(NP),W3(NP)
!
	REAL*8 DTAU
	REAL*8 T1
	INTEGER, PARAMETER :: IONE=1
	INTEGER LS,I,I1,I2
	INTEGER NX
	INTEGER ISPEC
	LOGICAL MID
!
	WRITE(139,*)'Entering ADJUST_R_GRID'
!
! Compute opitcal depth scale. We use the FLUXMEAN optical depth scale, except
! if it has a problem and is zero.
!
	DO I=1,ND
	  TA(I)=FLUXMEAN(I)
	  IF(TA(I) .LE. 0.0D0)TA(I)=ESEC(I)
	  TA(I)=TA(I)*CLUMP_FAC(I)
	END DO
        TB(1:ND)=0.0D0                              !Used for dCHIdR
        CALL NORDTAU(dTAU_OLD,TA,R,R,TB,ND)
!
! Save existing grid, which will be used for the interplations.
!
	R_OLD(1:ND)=R(1:ND)
	LOG_R_OLD=LOG(R_OLD)
!
! Compute optical depth scale. Note that we are passed the optical detph 
! increments, not the optical depth scale.
!
	TAU_OLD(1:ND)=0.0D0
	DO I=2,ND
	  TAU_OLD(I)=TAU_OLD(I-1)+dTAU_OLD(I-1)
	END DO
!
! Revise TAU near the outer boundary to allow for the fact that TAU_OLD(1)
! is zero. These expressions are correct for a power law opacity.
!
	TAU_OLD(1)=TAU_OLD(3)*R_OLD(3)/R_OLD(1)
	TAU_OLD(2)=TAU_OLD(3)*R_OLD(2)/R_OLD(1)
	TAU_OLD(1:ND)=DLOG10(TAU_OLD(1:ND))
!
	FG_MIN=-3.0D0
	FG_MAX=1.50D0
	FG_RANGE=FG_MAX-FG_MIN
!
	IF( TRIM(GRID_TYPE(1:6)) .EQ. 'UNIFOR')THEN
	  DTAU=(TAU_OLD(ND)-TAU_OLD(1))/(ND-6)
	  IF(N_PARS .EQ. 2)THEN
	    FG_MIN=RG_PARS(1)
	    FG_MAX=RG_PARS(2)
	    FG_RANGE=FG_MAX-FG_MIN
	  END IF
	  NX=(FG_RANGE+0.7D0*DTAU)/DTAU
	ELSE IF( TRIM(GRID_TYPE) .EQ. 'FIX_NX')THEN
	  NX=RG_PARS(1)
	  IF(N_PARS .EQ. 3)THEN
	    FG_MIN=RG_PARS(2)
	    FG_MAX=RG_PARS(3)
	    FG_RANGE=FG_MAX-FG_MIN
	  END IF
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid GRID_TYPE in ADJUST_R_GRID'
	  WRITE(LUER,*)'GRID_TYPE=',TRIM(GRID_TYPE)
	  STOP
	END IF 
!
! Compute the new grid. In the present version, NX points are
! inserted between log(TAU)=-3.0 and 1.5. We also have 2 special
! points at either boundary, which are placed separately.
!
	T1=TAU_OLD(ND)-TAU_OLD(1)-FG_RANGE
        DTAU=T1/(ND-5-NX)
	I1=(FG_MIN-TAU_OLD(1))/DTAU
	DTAU=(FG_MIN-TAU_OLD(1))/I1
        I1=I1+2
	TAU(1)=TAU_OLD(1)
	DO I=4,I1
	  TAU(I)=TAU_OLD(1)+DTAU*(I-3)
	END DO
        TAU(2)=MIN(TAU_OLD(2),TAU(1)+0.1D0*DTAU)
        TAU(3)=MIN(TAU_OLD(3),TAU(1)+0.3D0*DTAU)
!
! Do the insertion in the crtical section.
!
	DTAU=FG_RANGE/(NX-1)
	DO I=I1+1,I1+NX
	  TAU(I)=FG_MIN+DTAU*(I-I1-1)
	END DO
!
! Now do the last section, towards the inner boundary.
!
	I2=ND-(NX+I1+3)
	T1=(TAU_OLD(ND)-FG_MAX)
	IF(I2 .LE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in REVISE_R_GRID_V2'
	  WRITE(LUER,*)'Maximum TAU set is likely too high for model'
	  WRITE(LUER,*)'Change value in VADAT file'
	  STOP
	END IF
	DTAU=T1/I2
	DO I=I1+NX+1,ND-4
	  TAU(I)=TAU(I-1)+DTAU
	END DO
	TAU(ND)=TAU_OLD(ND)
	T1=1.0D0-10.0D0**(-DTAU)
	TAU(ND-1)=MAX(TAU_OLD(ND-1),TAU(ND)+LOG10(1.0D0-0.0667D0*T1))
	TAU(ND-2)=TAU(ND)+LOG10(1.0D0-0.201D0*T1)
	TAU(ND-3)=TAU(ND)+LOG10(1.0D0-0.467D0*T1)
!
!	DO I=1,ND
!	  WRITE(6,'(I4,4E16.8)')I,TAU(I),TAU_OLD(I),LOG_R_OLD(I)
!	END DO
!
! Compute the new radius grid. Linear interpolation in the log-log plane
! more than adequate, since we're just defining a new grid.
!
	CALL LININT(TAU,LOG_R,ND,TAU_OLD,LOG_R_OLD,ND)
	R=EXP(LOG_R); R(1)=R_OLD(1); R(ND)=R_OLD(ND)
	LOG_R(1)=LOG_R_OLD(1); LOG_R(ND)=LOG_R_OLD(ND)
!
! We now need to regrid all the populations. All interpolations (except 
! sigma) are performed in the LOG-LOG plane. For SN this is ideal, since
! the density and velocity are power laws in r. For SIGMA, we do not take
! the log.
!
! We do not need to interpolate T, and ED directly, since these are part
! of POPS.
!
	TA(1:ND)=LOG(V(1:ND))
	CALL MON_INTERP(V,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	V(1:ND)=EXP(V(1:ND))
!
	DO I=1,NT
	  TA(1:ND)=LOG(POPS(I,1:ND))
	  CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	  POPS(I,1:ND)=EXP(TB(1:ND))
	END DO
!
	RETURN
	END
