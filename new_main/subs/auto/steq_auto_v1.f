!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! AUTOINIZATION terms.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
! NB - ZION is the charge on the ion - thus ZHYD=1.0D0
!
! Routine does not work for NUM_BNDS=ND.
!
	SUBROUTINE STEQ_AUTO_V1(ED,T,
	1       HN_S,N_S,DI_S,
	1       HN_F,HNST_F,FEDGE_F,G_F,LEVNAME_F,N_F,
	1       F_TO_S_MAPPING,ID,AUTO_FILE,
	1       NUM_BNDS,ND,COMPUTE_BA,DST,DEND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
	INTEGER NT
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER DST,DEND
!
	INTEGER ID
	INTEGER N_S,N_F
!
	REAL*8 T(ND)
	REAL*8 ED(ND)
!
	REAL*8 HN_S(N_S,ND)
	REAL*8 DI_S(ND)
!
	REAL*8 HN_F(N_F,ND)
	REAL*8 HNST_F(N_F,ND)
	REAL*8 FEDGE_F(N_F)
	REAL*8 G_F(N_F)
	CHARACTER*(*) LEVNAME_F(N_F)
	CHARACTER*(*) AUTO_FILE
	INTEGER F_TO_S_MAPPING(N_F)
	LOGICAL COMPUTE_BA
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8 AUTO(N_F)
!
	INTEGER EQION
	INTEGER EQ_NUM_CONV
	INTEGER VION,VED,VT
	INTEGER I,J,K,M,IBEG_AUTO
	LOGICAL AUTO_FILE_EXISTS
!
! Define integer variables for population VARAIABLE for ease of programming.
!
	IBEG_AUTO=0
	DO I=1,N_F
	  IF(FEDGE_F(I) .LT. 0)THEN
	    IBEG_AUTO=I
	    EXIT
	  END IF
	END DO
	IF(IBEG_AUTO .EQ. 0)RETURN
!
! The autoionization probabilities are depth independent, and hence can be
! read in at the beginning.
!
        INQUIRE(FILE=AUTO_FILE,EXIST=AUTO_FILE_EXISTS)
        IF(AUTO_FILE_EXISTS)THEN
          CALL RD_AUTO_V1(AUTO,FEDGE_F,G_F,LEVNAME_F,N_F,AUTO_FILE)
        ELSE
          LUER=ERROR_LU()
          WRITE(LUER,*)'Warning: possible error in STEQ_AUTO_V1'
          WRITE(LUER,*)'No autoionization probabilities available'
          WRITE(LUER,*)'Model has states above the ionization limit'
          WRITE(LUER,*)'AUTO_FILE is',AUTO_FILE
          RETURN
        END IF
!
        EQION=N_S+1			!Ion equation
        VION=N_S+1
        VED=SE(ID)%N_IV-1
        VT=SE(ID)%N_IV
	EQ_NUM_CONV=SE(ID)%NUMBER_BAL_EQ
	M=(NUM_BNDS/2)+1
!
! 
!
	DO I=DST,DEND			!Which depth
!
	  DO K=IBEG_AUTO,N_F
	    J=F_TO_S_MAPPING(K)
	    SE(ID)%STEQ(J,I)=SE(ID)%STEQ(J,I)+AUTO(K)*(HNST_F(K,I)-HN_F(K,I))
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    DO K=IBEG_AUTO,N_F			!Which level
	      J=F_TO_S_MAPPING(K)
!
	      SE(ID)%BA(J,J,M,I)       =SE(ID)%BA(J,J,M,I)        -AUTO(K)*HN_F(K,I)/HN_S(J,I)
	      SE(ID)%BA(EQION,J,M,I)   =SE(ID)%BA(EQION,J,M,I)    +AUTO(K)*HN_F(K,I)/HN_S(J,I)
!
	      SE(ID)%BA(J,VION,M,I)    =SE(ID)%BA(J,VION,M,I)     +AUTO(K)*HNST_F(K,I)/DI_S(I)
	      SE(ID)%BA(EQION,VION,M,I)=SE(ID)%BA(EQION,VION,M,I) -AUTO(K)*HNST_F(K,I)/DI_S(I)
!
	      SE(ID)%BA(J,VED,M,I)     =SE(ID)%BA(J,VED,M,I)      +AUTO(K)*HNST_F(K,I)/ED(I)
	      SE(ID)%BA(EQION,VED,M,I) =SE(ID)%BA(EQION,VED,M,I)  -AUTO(K)*HNST_F(K,I)/ED(I)
! 
	      SE(ID)%BA(J,VT,M,I)      =SE(ID)%BA(J,VT,M,I)       -AUTO(J)*HNST_F(J,I)*
	1                                                               (1.5D0+HDKT*FEDGE_F(K)/T(I))/T(I)
	      SE(ID)%BA(EQION,VT,M,I)  =SE(ID)%BA(EQION,VT,M,I)   +AUTO(J)*HNST_F(J,I)*
	1                                                               (1.5D0+HDKT*FEDGE_F(K)/T(I))/T(I)
!
	     END DO          		!Over level variable
	  END IF
!
	END DO		!Loop over depth (I)
!
	RETURN
	END
