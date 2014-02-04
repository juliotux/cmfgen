	SUBROUTINE GAM_ABS(SM_R,SM_V,SM_CHI,SM_ND)
	IMPLICIT NONE
!
	INTEGER SM_ND
	REAL*8 SM_R(SM_ND)
	REAL*8 SM_V(SM_ND)
	REAL*8 SM_CHI(SM_ND)
!
	INTEGER, PARAMETER :: ND=300
	INTEGER, PARAMETER :: NC=50
	INTEGER, PARAMETER :: NP=ND+NC
	INTEGER, PARAMETER :: IONE=1
!
	REAL*8 R(ND)
	REAL*8 P(NP)
	REAL*8 V(ND)
!
	REAL*8 Z(ND),RJ(ND),Q(ND),F(ND)
	REAL*8 DTAU(ND),XM(ND),ETA(ND)
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 SOURCE(ND),CHI(ND),DCHIDR(ND),THETA(ND)
	REAL*8 JQW(ND,NP),KQW(ND,NP)
	REAL*8 WM(ND,ND),FB(ND,ND)
!
	REAL*8 IC,T1,DBB,HBC_J,HBC_S,INBC
	LOGICAL THK_CONT
	CHARACTER(LEN=6) METHOD
        CHARACTER(LEN=9) INNER_BND_METH
!
	INTEGER I
	INTEGER IOS
	CHARACTER(LEN=132) STRING
        EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
!
	METHOD='LOGLOG'
	INNER_BND_METH='DIFFUSION'
	THK_CONT=.FALSE.
	DBB=0.0D0
!
	OPEN(UNIT=10,IOSTAT=IOS,FILE='current_nonlocal_decay_energy.dat',STATUS='OLD',ACTION='READ')
	  IF(IOS .NE. 0)RETURN
	  READ(10,'(A)')STRING
	  READ(10,'(A)')STRING
	  DO I=1,ND
	    READ(10,*)V(I),T1,ETA(I)
	  END DO
	CLOSE(UNIT=10)
!
	CALL MON_INTERP(R,ND,IONE,V,ND,SM_R,SM_ND,SM_V,SM_ND)
	CALL MON_INTERP(CHI,ND,IONE,R,ND,SM_CHI,SM_ND,SM_R,SM_ND)
	WRITE(6,*)CHI(ND),SM_CHI(SM_ND)
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
        CALL GENANGQW(JQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,JTRPWGT,.FALSE.)
        CALL GENANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,KTRPWGT,.FALSE.)
!
	SOURCE=ETA/CHI
	CALL FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,R,Z,P,Q,F,
	1      SOURCE,CHI,DCHIDR,JQW,KQW,DBB,HBC_J,HBC_S,
	1      INBC,IC,THK_CONT,INNER_BND_METH,NC,ND,NP,METHOD)
!
	RJ=CHI*XM
	DO I=1,ND
	  WRITE(100,'(5ES14.6)')R(I),RJ(I),ETA(I),CHI(I),F(I)
	END DO
!
	WRITE(6,*)'Red curve is original data'
!
	CALL DP_CURVE(ND,V,ETA)
	CALL DP_CURVE(ND,V,RJ)
!
	THETA=0.0D0
	SOURCE=ETA/CHI
	CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1               SOURCE,THETA,CHI,F,JQW,
	1               THK_CONT,.TRUE.,DBB,IC,NC,ND,NP,METHOD)
	RJ=CHI*RJ
	CALL DP_CURVE(ND,V,RJ)
	DO I=1,ND
	  WRITE(100,'(5ES14.6)')R(I),RJ(I),ETA(I),CHI(I),F(I)
	END DO
!
	RETURN
	END