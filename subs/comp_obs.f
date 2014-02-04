C
C Routine to compute the observed spectrum from fluxes computed in the
C comoving frame.
C
C This routine should be called after each new set of boundary CMF intensities
C are computed. The CMF intensities for all rays should be passed.
C
C Routine stores CMF intensities. If CMF range is sufficient, fluxes for
C the next observers frame frequency(s) are computed.
C
	SUBROUTINE COMP_OBS(NEW_IPLUS,NEW_NU,
	1			IPLUS_STORE,NU_STORE,NST_CMF,
	1			MU,FLUX_WGHTS,OBS_FREQ,OBS_FLUX,N_OBS,
	1                       VINF,RMAX,IPLUS_OR_U,
	1                       INTERP_PROC,FIRST_OBS_COMP,NP)
	IMPLICIT NONE
C
C Altered 14-Dec-1996 : Bug fix for MON_INTER option. NU_STORE was being
C                         accesd outside valid range (1 to NEXT_ST_LOC-1).
C
	INTEGER NP
	REAL*8 NEW_IPLUS(NP)		!RAW CMF intensities as a function of
					!  impact parameter.
	REAL*8 NEW_NU			!Current CMF frequency.
C
	REAL*8 MU(NP)
	REAL*8 FLUX_WGHTS(NP)
C
	INTEGER N_OBS
	REAL*8 OBS_FREQ(N_OBS)
	REAL*8 OBS_FLUX(N_OBS)
C
C Storage arrays
C
	INTEGER NST_CMF
	REAL*8 IPLUS_STORE(NST_CMF,NP)
	REAL*8 NU_STORE(NST_CMF)
C
	REAL*8 VINF			!
	REAL*8 RMAX			!Radius at outer boundary.
	LOGICAL FIRST_OBS_COMP
	CHARACTER*(*) INTERP_PROC
	CHARACTER*(*) IPLUS_OR_U
C
C Local variables passed from one call to the next.
C
	REAL*8 C_KMS
	REAL*8 FLUX_CONST
	INTEGER NEXT_ST_LOC		!Keeps track of storage location.
	INTEGER OBS_INDX		!Current observers frequency
	INTEGER LUER
	SAVE C_KMS,NEXT_ST_LOC,OBS_INDX,LUER,FLUX_CONST
C
C External functions.
C
	REAL*8 SPEED_OF_LIGHT,PARSEC,FUN_PI
	INTEGER ERROR_LU
	EXTERNAL SPEED_OF_LIGHT,ERROR_LU,PARSEC,FUN_PI
C
C Local variables
C
	REAL*8 NU_SM_CMF
	INTEGER L,LS,ML,ML_ST,ML_END
C
C Variables for interpolation.
C
	REAL*8 FLUX,T1
	REAL*8 CMF_FREQ		!Observer's frequency transformed to comoving
				!  frame.
	REAL*8 HIM1,HI,HIP1
	REAL*8 SGN,SIM1,SI,SIP1
	REAL*8 DYI,DYIP1
	REAL*8 A,B,C,D
	REAL*8 ONE
	PARAMETER (ONE=1.0D0)
C
C Initialize variables if we a beginning the Observer flux calculation.
C
C NB: Flux constant =  2PI *
C                      (1.0E+10)^2 *		!From RMAX (units 10^10)
C                      1.0E+23			!1Jy = 1.0D-23 ergs/cm^2/sec
C                      / (3.0856E+21)^2		!Assume star at 1kpc
C
C We multiply FLUX_CONST by 2 if we are passing the Feautrier U variable,
C since
C        u=0.5(IPLUS+IMIN) = 0.5*IPLUS    (If IMIN=0)
C
	IF(FIRST_OBS_COMP)THEN
	  NEXT_ST_LOC=1
	  OBS_INDX=1
	  C_KMS=1.0D-05*SPEED_OF_LIGHT()
	  LUER=ERROR_LU()
	  FIRST_OBS_COMP=.FALSE.
	  IF(IPLUS_OR_U .EQ. 'IPLUS')THEN
	    FLUX_CONST=2.0D+01*FUN_PI()*(1.0D+18*RMAX/PARSEC())**2
	  ELSE IF(IPLUS_OR_U .EQ. 'U')THEN
	    FLUX_CONST=4.0D+01*FUN_PI()*(1.0D+18*RMAX/PARSEC())**2
	  ELSE
	    WRITE(LUER,*)'Error in COMP_OBS --- invalid IPLUS_OR_U'
	    WRITE(LUER,*)'IPLUS_OR_U = ',IPLUS_OR_U
	  END IF
	  IF(OBS_FREQ(1) .GT. NEW_NU)THEN
	    WRITE(LUER,*)'Invalid observer''s frequencies in COMP_OBS'
	    WRITE(LUER,*)'OBS_FREQ must be .LE. NEW_NU'
	    WRITE(LUER,*)'OBS_FREQ(1)=',OBS_FREQ(1)
	    WRITE(LUER,*)'NEW_NU=',NEW_NU
	  END IF
	END IF
C
	IF(OBS_INDX .GT. N_OBS)RETURN		!Finished
C
C Check to ensure still have sufficient storage for the new CMF frequency.
C If not we shuffle the frequencies to make room, checking that there is
C sufficient storage locations that we can still store the required CMF
C frequencies to enable the observer's frame fluxes to be computed.
C
C The L+2 ensures that there are 2 CMF frequencies greater than OBS_FREQ.
C NST_CMF-4 ensures that we retain 5 frequencies in the array. This check
C should only arrise if the spacing in the OBSERVER frequencies are very
C different to those in the co-moving frame.
C
	IF(NEXT_ST_LOC .GT. NST_CMF)THEN
	  L=1
	  DO WHILE ( NU_STORE(L+2) .GT. OBS_FREQ(OBS_INDX) .AND.
	1                L .LT. NST_CMF-4)
	    L=L+1
	  END DO
C
C L=1 indicates that more storage space is required.
C
	  IF(L .EQ. 1)THEN	
	    WRITE(LUER,*)'Error in COMP_OBS --- NST_CMF too small'
	    WRITE(LUER,*)'OBS_FREQ=',OBS_FREQ(OBS_INDX)
	    WRITE(LUER,*)'NU_STORE(1)=',NU_STORE(1)
	    WRITE(LUER,*)'NU_STORE(NST_CMF)=',NU_STORE(NST_CMF)
	    WRITE(LUER,*)'NEW_NU=',NEW_NU
	    STOP
	  END IF
	  DO LS=1,NP
	    DO ML=L,NST_CMF
	      IPLUS_STORE(ML+1-L,LS)=IPLUS_STORE(ML,LS)
	    END DO
	  END DO
	  DO ML=L,NST_CMF
	    NU_STORE(ML+1-L)=NU_STORE(ML)
	  END DO
	  NEXT_ST_LOC=NST_CMF-L+2
	END IF
C
C Store the newly computed comoving frame intensities into the storage arrays.
C We ensure that the vector can be stored into the array.
C
	DO LS=1,NP
	  IF(NEW_IPLUS(LS) .LE. 0)NEW_IPLUS(LS)=0.0
	  IPLUS_STORE(NEXT_ST_LOC,LS)=NEW_IPLUS(LS)
	END DO
	NU_STORE(NEXT_ST_LOC)=NEW_NU
	NEXT_ST_LOC=NEXT_ST_LOC+1
	IF(NEXT_ST_LOC .LT. 4)RETURN		!Not enough points for interp.
C
C Can we compute the flux for the next observer's frame frequency?
C
	IF(NU_STORE(1) .LT. OBS_FREQ(OBS_INDX))THEN
	  WRITE(LUER,*)'Invalid observer''s frequencies in COMP_OBS -2nd loc'
	  WRITE(LUER,*)'OBS_FREQ must be .LE. NEW_NU'
	  WRITE(LUER,*)'OBS_FREQ(OBS_INDX)=',OBS_FREQ(OBS_INDX)
	  WRITE(LUER,*)'NEW_NU=',NEW_NU
	  WRITE(LUER,*)'OBS_INDX=',OBS_INDX
	  STOP
	END IF
	NU_SM_CMF=OBS_FREQ(OBS_INDX)/(1+VINF/C_KMS)
	DO WHILE(NU_STORE(NEXT_ST_LOC-2) .LT. NU_SM_CMF)
C
C We can successfully perform the interpolation, and evaluate the flux
C for this observer's frame frequency.
C
	  OBS_FLUX(OBS_INDX)=0.0D0	!Initialize for integrations over p.
	  DO LS=1,NP
	    CMF_FREQ=OBS_FREQ(OBS_INDX)/(1.0D0+MU(LS)*VINF/C_KMS)
C
C Find location of frequency in vector.
C
	    ML_ST=1
	    ML_END=NEXT_ST_LOC-1
	    DO WHILE(ML_END-ML_ST .GT. 1)
	      ML=(ML_ST+ML_END)/2
	      IF(CMF_FREQ .LT. NU_STORE(ML))ML_ST=ML
	      IF(CMF_FREQ .GE. NU_STORE(ML))ML_END=ML
	    END DO
C
C Now ready for interpolation.
C
	    IF(INTERP_PROC .EQ. 'LIN_INT' .OR. ML_ST .EQ. 1)THEN
	      T1=(CMF_FREQ-NU_STORE(ML_END))/(NU_STORE(ML_ST)-NU_STORE(ML_END))
	      FLUX=T1*IPLUS_STORE(ML_ST,LS)+(1.0D0-T1)*IPLUS_STORE(ML_END,LS)
	    ELSE IF(INTERP_PROC .EQ. 'MON_INT')THEN
C
C The following procedure is from MON_INTERP, and was for a X vector
C which is either monotonically increasing, or decreasing. We have
C just changed the names of the variables, as required.
C
	      SGN=SIGN(ONE,NU_STORE(NEXT_ST_LOC-1)-NU_STORE(1))
	      IF( (SGN*CMF_FREQ .LT. SGN*NU_STORE(2)) .OR.
	1       (SGN*CMF_FREQ .GT. SGN*NU_STORE(NEXT_ST_LOC-2)) )THEN
	        WRITE(LUER,*)'Error in COMP_OBS- values outside range'
	       STOP
	      END IF
	      HI=NU_STORE(ML_ST+1)-NU_STORE(ML_ST)
	      HIM1=NU_STORE(ML_ST)-NU_STORE(ML_ST-1)
	      HIP1=NU_STORE(ML_ST+2)-NU_STORE(ML_ST+1)
	      SIM1=(IPLUS_STORE(ML_ST,LS)-IPLUS_STORE(ML_ST-1,LS))/HIM1
	      SI=(IPLUS_STORE(ML_ST+1,LS)-IPLUS_STORE(ML_ST,LS))/HI
	      SIP1=(IPLUS_STORE(ML_ST+2,LS)-IPLUS_STORE(ML_ST+1,LS))/HIP1
	      DYI=(SIM1*HI+SI*HIM1)/(HIM1+HI)
	      DYIP1=(SI*HIP1+SIP1*HI)/(HI+HIP1)
	      DYI=( SIGN(ONE,SIM1)+SIGN(ONE,SI) )*
	1            MIN(ABS(SIM1),ABS(SI),0.5*ABS(DYI))
	      DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5*ABS(DYIP1))
	      T1=(CMF_FREQ-NU_STORE(ML_ST))
              A=(DYI+DYIP1-2.0*SI)/HI/HI
	      B=(3.0*SI-2.0*DYI-DYIP1)/HI
	      C=DYI
	      D=IPLUS_STORE(ML_ST,LS)
              FLUX=((A*T1+B)*T1+C)*T1+D
	    ELSE
	      WRITE(LUER,*)
	1         'Error --- invalid interpolation request in COMP_OBS.'
	      WRITE(LUER,*)'INTERP=',INTERP_PROC
	      STOP
	    END IF
	    OBS_FLUX(OBS_INDX)=OBS_FLUX(OBS_INDX)+FLUX_WGHTS(LS)*FLUX
	  END DO
C
C Put flux in Jansky's for star at 1kpc.
C
	OBS_FLUX(OBS_INDX)=OBS_FLUX(OBS_INDX)*FLUX_CONST
C
C Ready for next frequency.
C
	  OBS_INDX=OBS_INDX+1
	  IF(OBS_INDX .GT. N_OBS)RETURN		!Finished
	  NU_SM_CMF=OBS_FREQ(OBS_INDX)/(1+VINF/C_KMS)
	END DO
C
	RETURN
	END
