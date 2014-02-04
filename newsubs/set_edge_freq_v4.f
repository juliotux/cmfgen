!
! Routine to place bound-free edge frequencies into a vector. Routine checks
! that cross-section at the bound-free edge is non-zero.
!
! We include only the initial edge for each super level, unless DO_ALL_EDGES
! is set.
!
	SUBROUTINE SET_EDGE_FREQ_V4(ID,EDGE_FREQ,EDGE_TYPE,NCF,NCF_MAX,
	1                       EDGEHI_F,NHI_F,HI_PRES,
	1                       F_TO_S_HI,NHI_S,NPHOT,
	1                       NUM_IMP_LEVELS,DO_ALL_EDGES,
	1                       INDICATE_SECONDARY_LEVELS,DO_LEVEL_DISSOLUTION)
	IMPLICIT NONE
!
! Altered 18-Jun-2003 : Changed to V4
!                       New options installed.
! Altered 01-Oct-2002 : DO_AL_EDGES option installed. EDGE values are returned
!                         for all levels in a SUPER level. EDGE_TYPE also installed.
! Altered 03-Mar-2000 : ID inserted in call. SUB_PHOT removed and replaced
!                          by SUB_PHOT_GEN.
! Altered 19-Dec-1997 : Double call to SUB_PHOT. This allows edge frequencies
!                         and then cross-sections to be returned. NEW_PHOT
!                         is altered on each call.
! Altered 14-Dec-1996 : Call to PHOT_FUN replaced by call to SUB_PHOT.
!                       Bug fix; DONE now re-zeroed for each PHOT_ID.
! Altered 26-May-1996 : ERROR_LU installed.
!                       DONE now dimensioned dynamically.
! Altered 24-Nov-1995 : Minor bug fix.
! Altered 24-Oct-1995
! Created 29-Mar-1990
!
	INTEGER ID
!
! NCF is the number of edges found (IN/OUT). On entry, it contains the number of
! edges found on all previous calls to SET_EDGE_FREQ_V4. On exit, it
! contains the updated value.
!
	INTEGER NCF
        INTEGER NCF_MAX               !Maximum number of edges from ALL species.
	INTEGER NHI_F                 !Number of full levels for species.
	INTEGER NHI_S                 !Number of super-levels for species.
	INTEGER NPHOT			!Number of photoionization routes.
	INTEGER NUM_IMP_LEVELS	!Number of level which must have a continuum edge.
!
	REAL*8 EDGEHI_F(NHI_F)
	INTEGER F_TO_S_HI(NHI_F)
	LOGICAL HI_PRES
!
	REAL*8 EDGE_FREQ(NCF_MAX)
!
! EDGE_TYPE(1:1)=
!                 I : Important level ---
!                 S : Secondary level (indicates not first level within a super-level)
! EDGE_TYPE(2:2)=
!                 D : Include frequencies for level dissolution for this level.
!
	CHARACTER(LEN=2) EDGE_TYPE(NCF_MAX)
!
! DO_ALL_EDGES tells SET_EDGE FREQ to include all edges with non-zero
! photoionization cross-sections. If FALSE, 2nd and higher levels belonging
! to a SL are not returned.
!
	LOGICAL DO_ALL_EDGES
!
! If TRUE, 2nd and higher levels belong to a SL are returned with 'S' in their
! EDGE_TYPE.
! 
	LOGICAL INDICATE_SECONDARY_LEVELS
!
! Tells SET_EDGE_FREQ whether Level Dissolution is on (only for PHOT_ID=1).
! Level dissolution is only switched on for levels in full atom) lower than
! NUM_IMP_LEVELS.
!
	LOGICAL DO_LEVEL_DISSOLUTION
!
! Local variables.
!
	REAL*8 T1
	INTEGER NEW_ID
	INTEGER PHOT_ID
	INTEGER J
	INTEGER IS
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	REAL*8, PARAMETER :: ZERO=0.0D0
	LOGICAL, PARAMETER :: RET_EDGE_CROSS=.TRUE.
!
! Used to get photoionization cross-sections.
!
	REAL*8 PHOT_CROSS(NHI_F)
	REAL*8 LOC_EDGE_FREQ(NHI_F)
	LOGICAL DONE(NHI_S)
!
	WRITE(167,*)ID
	WRITE(167,*)DO_LEVEL_DISSOLUTION
	WRITE(167,*)NUM_IMP_LEVELS
!
	IF( .NOT. HI_PRES )RETURN
!
	WRITE(6,*)'Number of important leveles is:',NUM_IMP_LEVELS
	DO PHOT_ID=1,NPHOT
	  DONE(:)=.FALSE.
!
! Get the ionization edges. These are returned in LOC_EDGE_FREQ.
!
	  NEW_ID=-PHOT_ID
	  CALL SUB_PHOT_GEN(ID,LOC_EDGE_FREQ,ZERO,EDGEHI_F,NHI_F,
	1                 NEW_ID,RET_EDGE_CROSS)
!
! Now get the photoionization cross-sections at the edge-frequency.
!
	  NEW_ID=100+PHOT_ID
	  CALL SUB_PHOT_GEN(ID,PHOT_CROSS,ZERO,EDGEHI_F,NHI_F,
	1                 NEW_ID,RET_EDGE_CROSS)
!
! Return those edges with non-zero photoionization cross-sections.
!
	  DO J=1,NHI_F
	    IS=F_TO_S_HI(J)
	    T1=PHOT_CROSS(J)
	    IF( T1 .GT. 0)THEN
	      IF(.NOT. DONE(IS))THEN
	        NCF=NCF+1
	        IF(NCF .GT. NCF_MAX)GOTO 9999
	        EDGE_FREQ(NCF)=LOC_EDGE_FREQ(J)
	        EDGE_TYPE(NCF)=' '
	        IF(J .LE. NUM_IMP_LEVELS)EDGE_TYPE(NCF)(1:1)='I'
	        IF(PHOT_ID .EQ. 1 .AND. DO_LEVEL_DISSOLUTION .AND.
	1           J .LE. NUM_IMP_LEVELS)EDGE_TYPE(NCF)(2:2)='D'
	        WRITE(167,*)NCF,EDGE_TYPE(NCF)
	        DONE(IS)=.TRUE.
	      ELSE IF(DO_ALL_EDGES)THEN
	        NCF=NCF+1
	        IF(NCF .GT. NCF_MAX)GOTO 9999
	        EDGE_FREQ(NCF)=LOC_EDGE_FREQ(J)
	        EDGE_TYPE(NCF)=''
	        IF(INDICATE_SECONDARY_LEVELS)THEN
	          EDGE_TYPE(NCF)(1:1)='S'
	        END IF
	        IF(J .LE. NUM_IMP_LEVELS)EDGE_TYPE(NCF)(1:1)='I'
	        IF(PHOT_ID .EQ. 1 .AND. DO_LEVEL_DISSOLUTION .AND.
	1           J .LE. NUM_IMP_LEVELS)THEN
	          EDGE_TYPE(NCF)(2:2)='D'
	        END IF
	      END IF
	    ELSE IF(T1 .LT. 0)THEN
	       LUER=ERROR_LU()
	       WRITE(LUER,*)'Error in SET_EDGE_FREQ_V4 - cross section negative'
	       WRITE(LUER,*)'ID=',ID,'Level=',J
	       STOP
	    END IF
	  END DO
	END DO
!
	RETURN
!
9999	CONTINUE
	LUER=ERROR_LU()
	WRITE(LUER,*)'Error NCF > NCF_MAX in SET_EDGE_FREQ_V4'
	STOP
!
	END
