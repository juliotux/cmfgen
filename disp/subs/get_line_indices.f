!
! Simple subroutine designed to return the line indices for a transition belonging to a
! given species. The line indices are either read in, or deduced from the line wavelength.
!
	SUBROUTINE GET_LINE_INDICES(VEC_FREQ,VEC_MNL_F,VEC_MNUP_F,VEC_SPEC,VEC_TRANS_NAME,
	1                 NLINE_FREQ,XSPEC,NL,NUP,LINE_FOUND)
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Created 10-Feb-2011 : 
!
	INTEGER NLINE_FREQ
	INTEGER NL
	INTEGER NUP
	LOGICAL LINE_FOUND
!
	REAL*8 VEC_FREQ(NLINE_FREQ)
	INTEGER VEC_MNL_F(NLINE_FREQ)
	INTEGER VEC_MNUP_F(NLINE_FREQ)
	CHARACTER(LEN=*) VEC_SPEC(NLINE_FREQ)
	CHARACTER(LEN=*) VEC_TRANS_NAME(NLINE_FREQ)
	CHARACTER(LEN=*) XSPEC
!
! Local variables
!
	REAL*8 ANG_TO_HZ
	REAL*8 FREQ
	REAL*8 LAMBDA
	REAL*8 T1
	REAL*8 dLAM
!
	INTEGER PNT(10)
	INTEGER LEV(2)
	INTEGER I,J,K
	INTEGER CNT
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	CHARACTER(LEN=10) DEFAULT
!
	REAL*8 SPEED_OF_LIGHT
	INTEGER GET_INDX_DP
	CHARACTER(LEN=30) UC
	EXTERNAL SPEED_OF_LIGHT,GET_INDX_DP,UC
!
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07      !10^8/10^1
	
	CALL USR_OPTION(LEV,ITWO,ITWO,'Levels','-1,-1','NL and NUP')
	NL=LEV(1); NUP=LEV(2)
	IF(NL .GE. 0 .AND. NUP .GE. 0)THEN
	  LINE_FOUND=.TRUE.
	  RETURN
	END IF
!
	CALL USR_OPTION(LAMBDA,'LAMBDA',' ','Wavelength of transition in Angstroms')
	FREQ=ANG_TO_HZ/LAMBDA
	K=GET_INDX_DP(FREQ,VEC_FREQ,NLINE_FREQ)
!
! Get nearby lines at higher frequencies (maximum of 5).
! We limit the search area to 0.1% in wavelength.
!
	CNT=0
	DO J=K,1,-1
	  IF(ABS(FREQ-VEC_FREQ(J))/FREQ .GT. 0.001D0)EXIT
	  IF( UC(VEC_SPEC(J)) .EQ. XSPEC)THEN
	    CNT=CNT+1
	    PNT(CNT)=J
	  END IF
	  IF(CNT .EQ. 5)EXIT
	END DO
!
! Get nearby lines at lower frequencies (maximum of 5).
!
	DO J=K+1,NLINE_FREQ
	  IF(ABS(FREQ-VEC_FREQ(J))/FREQ .GT. 0.001D0)EXIT
	  IF( UC(VEC_SPEC(J)) .EQ. XSPEC)THEN
	    CNT=CNT+1
	    PNT(CNT)=J
	  END IF
	  IF(CNT .EQ. 10)EXIT
	END DO
!
! Handle special cases.
!
	IF(CNT .EQ. 0)THEN
	  WRITE(6,*)'No transition found -- exiting line lookup'
	  LINE_FOUND=.FALSE.
	  RETURN
	ELSE IF(CNT .EQ. 1)THEN
	  K=PNT(1)
	  NL=VEC_MNL_F(K)
	  NUP=VEC_MNUP_F(K)
	  dLAM=ABS(ANG_TO_HZ/VEC_FREQ(K)-ANG_TO_HZ/FREQ)
	  WRITE(6,'(F9.2,2X,F5.3,3X,I3,3X,I3,4X,A)')ANG_TO_HZ/VEC_FREQ(K),dLAM,
	1          VEC_MNL_F(K),VEC_MNUP_F(K),TRIM(VEC_TRANS_NAME(K))	
	  LINE_FOUND=.TRUE.
	  RETURN
	END IF
!
! Many lines near line of interest -- let user decide which one they want.
!
	WRITE(6,*)'The program found the following lines belonging to species ',TRIM(VEC_SPEC(PNT(1)))
	WRITE(6,'(A,3X,A,4X,A,3X,A)')' Lambda(A)','dLAM(A)','NL','NUP'
	T1=100.0D0
	DO K=1,CNT
	  IF( ABS(FREQ-VEC_FREQ(PNT(K))) .LT. T1)THEN
	    T1= ABS(FREQ-VEC_FREQ(PNT(K)))
	    I=K
	  END IF
	  dLAM=ABS(ANG_TO_HZ/VEC_FREQ(PNT(K))-ANG_TO_HZ/FREQ)
	  WRITE(6,'(F10.2,2X,F8.3,3X,I3,3X,I3,4X,A)')ANG_TO_HZ/VEC_FREQ(PNT(K)),dLAM,
	1          VEC_MNL_F(PNT(K)),VEC_MNUP_F(PNT(K)),TRIM(VEC_TRANS_NAME(PNT(K)))	
	END DO
	WRITE(6,'(A)')' '
!
	WRITE(DEFAULT,'(I3,A,I3)')VEC_MNL_F(PNT(I)),',',VEC_MNUP_F(PNT(I)) 
	CALL USR_OPTION(LEV,ITWO,ITWO,'Levels',DEFAULT,'NL and NUP')
	NL=LEV(1); NUP=LEV(2)
	LINE_FOUND=.TRUE.
!
	RETURN
	END
