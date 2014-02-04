!
! Routine to generate frequencies for radiative transfer program given
! a set of bound-free edges, and maximum and minimum frequencies.
! The frequency spacing is contolled by three passable parameters.
!
! Every bound-free edged is defined by two frequencies, except
! where two bound-free edges coincide to 1 part in 10^11. Generally
! points are inserted so that simpsons rule can be used, although in
! some regions where the bound-free edges are close the trapazoidal
! rule will be used.
!
	SUBROUTINE SET_CONT_FREQ_V2(NEW_FREQ,CONT_TYPE,EDGE,TYPE,INDX,
	1                        SMALL_RAT,BIG_AMP,DNU_MAX,
	1                        MAX_FREQ,MIN_FREQ,
	1                        dV_LEV,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        dV_CONT,dV_EDGE_SEP,
	1                        N,NCF,NCF_MAX,LUOUT)
	IMPLICIT NONE
!
! Altered 26-May-1996 : ERROR_LU installed.
! Cleaned
! Created 29-Mar-1990 (Based on GEN_FREQ).
!
	INTEGER*4 N,NCF,NCF_MAX,LUOUT,INDX(NCF_MAX)
	REAL*8 EDGE(NCF_MAX)
	CHARACTER*1 CONT_TYPE(NCF_MAX)
	REAL*8 NEW_FREQ(NCF_MAX)
	CHARACTER*(*) TYPE(NCF)
!
	REAL*8 MAX_FREQ		!Maximum continuum frequency
	REAL*8 MIN_FREQ		!Minimum continuum frequency
!
	REAL*8 DNU_MAX		!Maximum frequency spacing near/above
				!  bound-free edge: i.e. dNU <  DNU_MAX
	REAL*8 BIG_AMP  	!Amplification. dNU increases by a factor
				! BIG_AMP as we move away from the b.f. edge.
				! for frequencies above SWITCH_FREQ.
   	REAL*8 SMALL_RAT	!Used to define frequency spacing for
				! frequencies less than SWITCH_FREQ.
				! dNU/NU=SMALL_RAT-1
!
! Parameters for installing etra frequencies near bound-free edges
! (low frequency side) to allow for level dissolution.
!
	REAL*8 dV_LEV			!Spacing near b-f edge.
	REAL*8 AMP_DIS			!Amplification factor for dNU as we
					!  move to smaller frequencies.
	REAL*8 MIN_FREQ_LEV_DIS		!Indicates that the extra frequencies
					!  should only be installed for
					!  frequencies above MIN_FREQ_LEV_DIS.
	REAL*8 dV_EDGE_SEP			!Minimum spacing in Doppler line profiles
                                        !  (km/s). Used to set minimum frequecny
                                        !  spacing.
	REAL*8 dV_CONT                  !Maximum spacing in arbitrary section of 
                                        !  continuum (km/s)
!
	INTEGER*4, PARAMETER :: RZERO=0.0D0
	INTEGER*4, PARAMETER :: RONE=1.0D0
	INTEGER*4, PARAMETER :: RTWO=2.0D0
!
! N+2 as we insert MAX_FREQ and MIN_FREQ into the array.
!
	CHARACTER*10   CHAR_WRK(N+2)
! 
	REAL*8 C_KMS
	REAL*8 SPEED_OF_LIGHT
	INTEGER*4 ERROR_LU,LUER
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
!
	REAL*8 T1,T2
	REAL*8 RN
	REAL*8 COINCIDENT_FRAC
	REAL*8 DEL_NU
	REAL*8 DEL_NU_TO_NEXT_EDGE
!
	INTEGER*4 INDX_DIS
	INTEGER*4 I,J,K
	LOGICAL NUMER,EQUAL,EDGE_FREQ
	REAL*8 FAC,EQUAL_FAC
	REAL*8 SMALL_FAC
	REAL*8 SWITCH_FREQ
	CHARACTER(LEN=1) TMP_CHAR
!
	LUER=ERROR_LU()
	COINCIDENT_FRAC=0.5D0
	CONT_TYPE(:)=' '
!
! Sort frequencies into numerical order. NEW_FREQ is used as a
! work array.
!
	N=N+2
	EDGE(N-1)=MIN_FREQ;   TYPE(N-1)=' '
	EDGE(N)=MAX_FREQ;     TYPE(N)=' '
	NUMER=.TRUE.
	CALL INDEXX(N,EDGE,INDX,NUMER)
	CALL SORTDP(N,EDGE,INDX,NEW_FREQ)
	CALL SORTCHAR(N,TYPE,INDX,CHAR_WRK)
!
! Do a quick check that sort was done correctly.
! Not equal as MIN_FREQ, and MAX_FREQ have been inserted in FREQ array.
!
	IF(MIN_FREQ .NE. EDGE(1))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - MIN_FREQ too big'
	  WRITE(LUER,*)'MIN_FREQ =',MIN_FREQ
	  WRITE(LUER,*)'EDGE(1) =',EDGE(1)
	  STOP
	END IF
	IF(MAX_FREQ .NE. EDGE(N))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - MAX_FREQ too small'
	  WRITE(LUER,*)'MAX_FREQ =',MAX_FREQ
	  WRITE(LUER,*)'EDGE(N) =',EDGE(N)
	  STOP
	END IF
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	EQUAL_FAC=COINCIDENT_FRAC*dV_EDGE_SEP/C_KMS
	FAC = RONE + EQUAL_FAC
	SMALL_FAC = RONE + 1.0D-06
!
! SMAL_FAC is the ratio used to set the frequency spacing for
! small frequencies (frequencies less than 1.0 approximately).
! Two points are inserted so that simpsons rule can be used.
!
	IF(SMALL_RAT .LE. RONE .OR. SMALL_RAT .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - SMALL_RAT Outside ',
	1           'valid range ',SMALL_RAT
	  STOP
	END IF
!
! DNU_MAX is the maximum frequecy spacing for frequencies adjacent
! to a bound-free edge. For large v, we require constant spacing since
! error in integral is a function of {del v }.
!
	IF(DNU_MAX .LE. RZERO .OR. DNU_MAX .GT. RTWO)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - DNU_MAX Outside ',
	1           'valid range ',DNU_MAX
	  STOP
	END IF
!
! BIG_AMP is used to amplify the spacing for large frequencies. Close
! to the edge, the spacing is DNU_MAX, but for every
! points, the spacing increases by a factor of BIG_AMP
!
	IF(BIG_AMP .LT. RONE .OR. BIG_AMP .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ - BIG_AMP ',
	1                 'Outside valid range ',BIG_AMP
	  STOP
	END IF
!
! Determin the freuqncy where we switch from using DNU_MAX to
! SMALL_RAT.
!
	SWITCH_FREQ=DNU_MAX/(SMALL_RAT-RONE)
!
! Now begin inserting extra points into frequency array.
!
	INDX_DIS=0
	K=1
	NEW_FREQ(1)=MIN_FREQ
	EDGE_FREQ=.FALSE.
	I=2
	DO WHILE(I .LE. N)
!
! Get next edge for which we will consider level dissolution.
!  
          IF(I .GT. INDX_DIS)THEN
	    INDX_DIS=I
	    DO WHILE( INDEX(TYPE(INDX_DIS),'D') .EQ. 0 .AND. INDX_DIS .LT. N)
	      INDX_DIS=INDX_DIS+1
 	    END DO
	  END IF
!
! Determine whether next EDGE frequency is too close to selected frequency
! already.
!
	  DEL_NU=EDGE(I)-NEW_FREQ(K)
	  IF(DEL_NU .LE. 0)THEN
	     I=I+1
	     EDGE_FREQ=.TRUE.
	    GOTO 1000
	  END IF
!
	  IF(TYPE(I) .EQ. 'S')THEN
	    IF(DEL_NU .GT. 1.5d0*dV_CONT/C_KMS*EDGE(I))THEN
	      DEL_NU=EDGE(I)*dV_CONT/C_KMS
	    END IF
	  ELSE
!  
	    DEL_NU_TO_NEXT_EDGE=DEL_NU
	    IF(DEL_NU/NEW_FREQ(K) .LT. EQUAL_FAC)THEN
	       IF(.NOT. EDGE_FREQ)THEN
                 NEW_FREQ(K)=EDGE(I)
	         EDGE_FREQ=.TRUE.
	         I=I+1
	         GOTO 1000
	       ELSE
                 DEL_NU=NEW_FREQ(K)*EQUAL_FAC
                 NEW_FREQ(K+1)=NEW_FREQ(K)+DEL_NU
	         EDGE_FREQ=.TRUE.
	         K=K+1
	         I=I+1
	         GOTO 1000
	       END IF
	    END IF
	    WRITE(50,'(A1,2ES14.6)')'A',NEW_FREQ(K),DEL_NU
	  END IF
!
! dV_CONT is the minimum spacing (velocity space) at which we have to sample the
! photoioinization cross-sections.
! 
	  T1=NEW_FREQ(K)*dV_CONT/C_KMS
	  IF(T1 .LT. DEL_NU)THEN
	     DEL_NU=T1
	     EDGE_FREQ=.FALSE.
	  WRITE(50,'(A1,2ES14.6)')'C',NEW_FREQ(K),DEL_NU
	  END IF
!
! For frequency close to a bound-free edge (but on the high side) we
! use a spacing DNU_MAX near the edge. As we move blueward,
! the spacing increaes by BIG_AMP on each succesive insertion.
!
	  IF(DEL_NU .GT. DNU_MAX .AND. NEW_FREQ(K) .GT. SWITCH_FREQ)THEN
	    T1=(NEW_FREQ(K)-EDGE(I-1))/DNU_MAX
	    RN=LOG( 1.0D0+T1*(BIG_AMP-1.0D0))/LOG(BIG_AMP) - 1.0D0
	    T1=DNU_MAX*(BIG_AMP**RN)
	    IF(T1 .LT. DEL_NU)THEN
	      DEL_NU=T1
	      EDGE_FREQ=.FALSE.
	  WRITE(50,'(A1,2ES14.6)')'D',NEW_FREQ(K),DEL_NU
	    END IF
	  END IF
!
! Small frequencies where we use a fixed velocity spacing. dV_CONT
! has essentially the same meaning, and may supercede SMALL_RAT.
!
	  T1=NEW_FREQ(K)*(SMALL_RAT-1.0D0)
	  IF(DEL_NU .GT. T1 .AND. NEW_FREQ(K) .LE. SWITCH_FREQ)THEN
	    DEL_NU=T1
	    EDGE_FREQ=.FALSE.
	  WRITE(50,'(A1,2ES14.6)')'E',NEW_FREQ(K),DEL_NU
	  END IF
!
! Now check for level dissolution. We only apply level dissolution to
! those levels which have a 'D' specified in the TYPE variable.
!
	  WRITE(41,*)NEW_FREQ(K),EDGE(INDX_DIS),RN,T2,T1
	  IF(EDGE(INDX_DIS) .GT. MIN_FREQ_LEV_DIS .AND. INDX_DIS .NE. N)THEN
	    T2=EDGE(INDX_DIS)*dV_LEV/C_KMS
	    T1=(EDGE(INDX_DIS)-NEW_FREQ(K))/T2
	    RN=LOG( 1.0D0+T1*(AMP_DIS-1.0D0))/LOG(AMP_DIS) - 1.0D0
	    T1=T2*(AMP_DIS**RN)
	    WRITE(40,'(5E14.6)')NEW_FREQ(K),EDGE(INDX_DIS),RN,T2,T1
	    IF(T1 .LT. DEL_NU)THEN
	      DEL_NU=T1
	      EDGE_FREQ=.FALSE.
	  WRITE(50,'(A1,2ES14.6)')'F',NEW_FREQ(K),DEL_NU
 	    END IF
	  END IF
!
! Bracket the bound-free edge. We only do this IF dissolution for the level is
! NOT switched on.
! 
	  IF(EDGE(INDX_DIS) .LE. MIN_FREQ_LEV_DIS .OR. TYPE(INDX_DIS) .NE. 'D')THEN
	    IF(I .NE. N)THEN
	      T1=-EDGE(I+1)*dV_EDGE_SEP/C_KMS+(EDGE(I+1)-NEW_FREQ(K))
	      IF(EDGE_FREQ .AND. T1 .LT. DEL_NU .AND. T1 .GT. (EDGE(I+1)-NEW_FREQ(K)) )THEN
	        DEL_NU=T1
	        EDGE_FREQ=.FALSE.
	        WRITE(50,'(A1,2ES14.6)')'B',NEW_FREQ(K),DEL_NU
	      END IF
	    END IF
	  END IF
!
! Insert frequency edge.
!
	  IF(DEL_NU .EQ. DEL_NU_TO_NEXT_EDGE)THEN
	    K=K+1
	    NEW_FREQ(K)=EDGE(I)
	    I=I+1
	    EDGE_FREQ=.TRUE.
	  WRITE(50,'(A1,2ES14.6)')'G',NEW_FREQ(K),DEL_NU
	  ELSE
!
! Edge freq has alreadby been set to false.
!
	    K=K+1
	    NEW_FREQ(K)=NEW_FREQ(K-1)+DEL_NU
	    WRITE(50,'(A1,2ES14.6)')'H',NEW_FREQ(K),DEL_NU
	  END IF
!
	  IF(K .GT. NCF_MAX)THEN
	    WRITE(LUER,*)'Error NCF too small in SET_CONT_FREQ'
	    WRITE(LUER,*)N,NCF_MAX,SMALL_RAT,BIG_AMP,DNU_MAX
	    OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	      DO J=1,NCF
	        WRITE(LUOUT,100)NEW_FREQ(J)
	      END DO
	    CLOSE(LUOUT)
	    STOP
	  END IF
1000	  CONTINUE
	  IF(EDGE_FREQ)CONT_TYPE(K)='E'
!
	END DO				!I loop
!
	NCF=K
!
!
	K=1
	DO I=1,N-1
	  DO WHILE (EDGE(I) .GT. NEW_FREQ(K))
	    K=K+1
	  END DO 
	  WRITE(37,*)EDGE(I), 2.998D+05*(NEW_FREQ(K)-EDGE(I))/EDGE(I), 2.998D+05*(EDGE(I)-NEW_FREQ(K-1))/EDGE(I)
	END DO
!
! Now sort frequencies into numerical decreasing order. 
!
	DO I=1,NCF/2
	  T1=NEW_FREQ(I)
	  NEW_FREQ(I)=NEW_FREQ(NCF-I+1)
	  NEW_FREQ(NCF-I+1)=T1
	  TMP_CHAR=CONT_TYPE(I)
	  CONT_TYPE(I)=CONT_TYPE(NCF-I+1)
	  CONT_TYPE(NCF-I+1)=TMP_CHAR
	END DO
!
	OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	  DO I=1,NCF
	    WRITE(LUOUT,100)NEW_FREQ(I), 2.998D+05*(NEW_FREQ(I)-NEW_FREQ(MAX(I-1,1)))/NEW_FREQ(I),CONT_TYPE(I)
	  END DO
100	  FORMAT(1X,F20.16,4X,F10.3,4X,A1)
	CLOSE(LUOUT)
!
! Check frequency array is monotonic.
!
	DO I=1,NCF-1
	  IF( NEW_FREQ(I) .LE. NEW_FREQ(I+1) )THEN
	    WRITE(LUER,*)'Error in SET_CONT_FREQ - frequency array is',
	1             ' not monotonic.'
	    WRITE(LUER,*)I,NEW_FREQ(I),NEW_FREQ(I+1)
	    STOP
	  END IF
	END DO
!
! Ensure FREQ array is zeroed (as probably will use OBSF).
!
	DO I=1,NCF
	  EDGE(I)=0.0D0
	END DO
	WRITE(6,*)'Returning'
!
	RETURN
	END
