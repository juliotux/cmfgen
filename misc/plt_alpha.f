!
! Simple and crude auxilary program designed to estimate the CAK ALPHA
! parameter from two HYDRO files.
!
! Under development!
!
	PROGRAM PLT_ALPHA
	USE GEN_IN_INTERFACE
!
! Created : 12-Oct-2008 : Based on RVE_HYDRO_TURB and so some comments
!                         may stll apply to that file.
!
	IMPLICIT NONE
!
	INTEGER I,J,IOS,NM
	INTEGER NSTR
	INTEGER ND
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: T_OUT=6
!
	INTEGER, PARAMETER :: NMAX=200
	REAL*8 R(NMAX,2)
	REAL*8 V(NMAX,2)
	REAL*8 E(NMAX,2)
	REAL*8 VdVdR(NMAX,2)
	REAL*8 dPdR(NMAX,2)
	REAL*8 g_tot(NMAX,2)
	REAL*8 g_rad(NMAX,2)
	REAL*8 g_elec(NMAX,2)
	REAL*8 Gamma(NMAX,2)
!
	REAL*8 XV(NMAX)
	REAL*8 YV(NMAX)
!
	REAL*8 T1
!
	CHARACTER*132 STRING(500)
	CHARACTER*132 FMT
	CHARACTER*132 FILENAME
!
	REAL*8 P_VEL(200)
	REAL*8 P_REQ(200)
	REAL*8 P_GRAD(200)
	REAL*8 P_GELEC(200)
!
! Read in the two HYDRO files constructed by two different models.
! These models should be identical, except for a slight difference in VINF.
!
	DO NM=1,2
	  ND=0
	  STRING(:)=' '
	  FILENAME='HYDRO'
	  CALL GEN_IN(FILENAME,'Input hydro file')
	  OPEN(UNIT=10,FILE=FILENAME,ACTION='READ',STATUS='OLD')
	    NSTR=0        
	    DO WHILE(1 .EQ. 1)
	      READ(10,'(A)',END=1000)STRING(NSTR+1)
	      IF(STRING(NSTR+1) .EQ. ' ' .AND. ND .EQ. 0)ND=NSTR-1
	      NSTR=NSTR+1
	    END DO         
	  CLOSE(UNIT=10)
1000	  CONTINUE
!
	  WRITE(T_OUT,*)'Number of depth points is',ND
          I=160
!
	  DO I=1,ND
	    READ(STRING(I+1),*)R(I,NM),V(I,NM),E(I,NM),VdVdR(I,NM),dPdR(I,NM),g_TOT(I,NM),
	1                    g_RAD(I,NM),g_ELEC(I,NM),Gamma(I,NM)
	  END DO
!
	END DO
!
! We can now compute ALPHA as a function of depth. This procedure will only
! ba valid above the SONIC point.
!
	J=ND
	DO I=1,ND
	  T1=(VdVdR(I,2)-VdVdR(I,1))
	  IF(T1 .EQ. 0)EXIT
	  XV(I)=R(I,1)/R(ND,1)
	  YV(I)=(g_RAD(I,2)-g_RAD(I,1)-(g_elec(I,2)-g_elec(i,1)))/T1
	  WRITE(6,*)I,R(I,1),YV(I)
	  J=I
	END DO
!
! We plot ALPHA versus Log R and V.
!
	CALL DP_CURVE(J,XV,YV)
	CALL GRAMON_PGPLOT('R/R\d*\u','\ga',' ',' ')
	CALL DP_CURVE(J,V,YV)
	CALL GRAMON_PGPLOT('R/R\d*\u','\ga',' ',' ')
!
	STOP
	END
