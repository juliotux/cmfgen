C
C Routine to place bound-free edge frequencies into a vector. Routine checks
C that cross-section at the bound-free edge is non-zero.
C
C We include only the initial edge for each super level.
C
	SUBROUTINE SET_EDGE_FREQ_V2(FREQ,NCF,NCF_MAX,
	1                       EDGEHI_F,NHI_F,HI_PRES,
	1                       F_TO_S_HI,NHI_S,
	1                       SUB_PHOT,NPHOT)
	IMPLICIT NONE
C
C Altered 19-Dec-1997 : Double call to SUB_PHOT. This allows edge frequencies
C                         and then cross-sections to be roeturned. NEW_PHOT
C                         is altered on each call.
C Altered 14-Dec-1996 : Call to PHOT_FUN replaced by call to SUB_PHOT.
C                       Bug fix; DONE now rezeroed for each PHOT_ID.
C Altered 26-May-1996 : ERROR_LU installed.
C                       DONE now dimensioned dynamically.
C Altered 24-Nov-1995: Minor bug fix.
C Altered 24-Oct-1995
C Created 29-Mar-1990
C
	INTEGER*4 NCF,NCF_MAX,J
	INTEGER*4 NHI_F,NHI_S
	INTEGER*4 NPHOT
	REAL*8 FREQ(NCF_MAX)
	REAL*8 EDGEHI_F(NHI_F)
	INTEGER*4 F_TO_S_HI(NHI_F)
	LOGICAL HI_PRES
C
	INTEGER*4 PHOT_ID,NEW_ID
	INTEGER*4 IS,ERROR_LU,LUER
	REAL*8 T1
	EXTERNAL ERROR_LU,SUB_PHOT
C
	REAL*8, PARAMETER :: ZERO=0.0D0
	LOGICAL, PARAMETER :: RET_EDGE_CROSS=.TRUE.
C
	REAL*8 PHOT_CROSS(NHI_F)
	REAL*8 EDGE_FREQ(NHI_F)
	LOGICAL DONE(NHI_S)
C
	IF( .NOT. HI_PRES )RETURN
C
	DO PHOT_ID=1,NPHOT
	  DONE(:)=.FALSE.
C
C Get the ionization edges. These are returned in PHOT_CROSS.
C
	  NEW_ID=-PHOT_ID
	  CALL SUB_PHOT(PHOT_CROSS,ZERO,EDGEHI_F,NHI_F,
	1                 NEW_ID,RET_EDGE_CROSS)
	  EDGE_FREQ(1:NHI_F)=PHOT_CROSS(1:NHI_F)
C
C Now get the photoionzation crossections at the edge-frequency.
C
	  NEW_ID=100+PHOT_ID
	  CALL SUB_PHOT(PHOT_CROSS,ZERO,EDGEHI_F,NHI_F,
	1                 NEW_ID,RET_EDGE_CROSS)
C
	  DO J=1,NHI_F
	    IS=F_TO_S_HI(J)
	    T1=PHOT_CROSS(J)
	    IF( T1 .GT. 0)THEN
	       IF(.NOT. DONE(IS))THEN
	         NCF=NCF+1
	         IF(NCF .GT. NCF_MAX)GOTO 9999
	         FREQ(NCF)=EDGE_FREQ(J)
	         DONE(IS)=.TRUE.
	       END IF
	    ELSE IF(T1 .LT. 0)THEN
	       LUER=ERROR_LU()
	       WRITE(LUER,*)'Error in SET_EDGE_FREQ - cross section negative'
	       STOP
	    END IF
	  END DO
	END DO
C
	RETURN
C
9999	CONTINUE
	LUER=ERROR_LU()
	WRITE(LUER,*)'Error NCF > NCF_MAX in SET_EDGE_FREQ'
	STOP
C
	END
