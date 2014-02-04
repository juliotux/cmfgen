C
C Subroutine to compute the collisional excitation and ionization cross
C sections for an arbitrary species taking into account super levels.
C
C The collison rates among the levels in the full atom are first computed. 
C These are then used to compute the rates amongst the super levels.
C
C The collison strengths (OMEGA) must be supplied by OMEGA_COL ---
C this subroutine name is passed in the call. OMEGA_COL uses data in
C the file COL_FILE. SUB_PHOT is also passed to OMEGA_COL.
C
C Single routine as technique is the same for ALL species.
C
	SUBROUTINE SUBCOL_MULTI_V4(
	1                 OMEGA_F,dln_OMEGA_dlnT,
	1                 COL_S,dCOL_S,
	1                 HN_S,HNST_S,dlnHNST_S_dlnT,N_S,
	1                 HN_F,HNST_F,W_F,EDGE_F,AHYD_F,GHYD_F,
	1                 LEVNAME_F,N_F,
	1                 ZION,ID,COL_FILE,OMEGA_COL,
	1                 F_TO_S_MAPPING,COOL,T,ED,ND)
	IMPLICIT NONE
C
C Altered 16-Jun-1996 : Call to ZERO removed. COL and dCOL now initializd
C                         outside depth loop.
C                       Bug FIX: COL_FILE was declared REAL*8, NOW declared
C                           as CHARACTER. No efffect on VAX, important on CRAY.
C Altered 27-May-1996 : Generic calls used for EXP, SQRT.
C Altered 03-Jan-1995 - HN_F  inserted in call (_V2 changed to _V3)
C                       Adjustment made to allow for variable departure
C                          coeficients in a given super level (collisional
C                          excitation nad dexciation only).
C                       Collisional ionization rates still assume constant
C                          departure coeficients among levels in a given super
C                          levl.
C  
C Altered 24-Nov-1995 - Bug fixed: Incorrect dimension for HN_S
C Altered 10-Nov-1995 - Bug fixed in DCOL.
C                       HN_S inserted in call, and HN_F deleted.
C Created 06-Sep-1995 - Based on SUBCOL_HYD version.
C                       Created so that all species with super-levels use
C                         the same collisional routine. Only the routine to
C                         compute OMEGA in the FULL atom is distinct.
C                       Vector COOL was installed to enable checking of 
C                         collisonal cooling rates.
C
	INTEGER N_S,N_F,ND  
	REAL*8 OMEGA_F(N_F,N_F),dln_OMEGA_dlnT(N_F,N_F)
	REAL*8 COL_S(N_S,N_S,ND),DCOL_S(N_S,N_S,ND)
C
	REAL*8 HN_S(N_S,ND)		!Population of atom with super levels.
	REAL*8 HNST_S(N_S,ND)		!LTE pop. of atom with super levels.
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HN_F(N_F,ND)		!Population of FULL atom
	REAL*8 HNST_F(N_F,ND)		!LTE population of FULL atom
	REAL*8 W_F(N_F,ND)		!Occupation probability
	REAL*8 AHYD_F(N_F,N_F)		!Einstein A coefficient
	REAL*8 EDGE_F(N_F)		!Ionization frequency (10^15 Hz)
	REAL*8 GHYD_F(N_F)		!Statistical weight.
	REAL*8 ZION			!Charge on ion (i.e. 1 for H)
	CHARACTER*(*) COL_FILE		!Name of file with collisonal data.
	CHARACTER*(*) LEVNAME_F(N_F)	!Level names in FULL ATOM.
	INTEGER ID			!Specifies ident. for photiozation data
C
	INTEGER F_TO_S_MAPPING(N_F)
C
	REAL*8 T(ND)			!Temperature (10^4 K)
	REAL*8 ED(ND)			!Electron density
C
C NB: On exit COOL needs to be multiplied by COOL.
C
	REAL*8 COOL(ND)			!Net collisional cooling
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	EXTERNAL OMEGA_COL
C
	INTEGER I,J,K
	INTEGER L,U
	REAL*8 X
	REAL*8 BRAT
	REAL*8 CIJ,CJI,CII
C
	COL_S(:,:,:)=0.0D0  		!N_S, N_S, ND
	dCOL_S(:,:,:)=0.0D0  		!N_S, N_S, ND
C
C Loop over all depths. This is provided for consistency with other collisonal
C routines. We operate on a single depth at a time to minimize the size of 
C OMEGA and dln_OMEGA_dlnT.
C
	DO K=1,ND
C
C The following correspond to the collision rates between levels I and J in 
C the FULL atom.
C
	  OMEGA_F(:,:)=0.0D0
	  dln_OMEGA_dlnT(:,:)=0.0D0
C
	  CALL OMEGA_COL(OMEGA_F,dln_OMEGA_dlnT,EDGE_F,AHYD_F,GHYD_F,
	1                  LEVNAME_F,ZION,N_F,T(K),
	1                  ID,COL_FILE)
C
	  DO I=1,N_F-1
	    DO J=I+1,N_F                           
C
	      CIJ=8.63D-08*ED(K)*OMEGA_F(I,J)/SQRT(T(K))
	      CJI=CIJ/GHYD_F(J)
	      X=HDKT*(EDGE_F(I)-EDGE_F(J))/T(K)
	      CIJ=CIJ*EXP(-X)/GHYD_F(I)
C
C Allow for collisional ionization through level dissolution.
C                        
	      L=F_TO_S_MAPPING(I)    
	      CII=CIJ*HNST_F(I,K)/HNST_S(L,K)*(1.0D0-W_F(J,K)/W_F(I,K))
	      COL_S(L,L,K)=COL_S(L,L,K)+CII
	      DCOL_S(L,L,K)=DCOL_S(L,L,K)+CII*( dln_OMEGA_dlnT(I,J) +
	1        X -2.0D0 - HDKT*EDGE_F(I)/T(K)-dlnHNST_S_dlnT(L,K) )/T(K)
C
C We use the full ionization energy because of the following argument.
C We think of the ionization as a 3 body process. The extra ionization
C energy comes from the third electron and hence is lost from the electron
C thermal pool.
C
	      COOL(K)=COOL(K)+EDGE_F(I)*CII*(HN_S(L,K)-HNST_S(L,K))
C
C The following section is independent of the atomic structure. We use
C F_TO_S_MAPPING to describe how the FULL atom is mapped on to the smaller
C atom. Note that the collsion rate is only important when L is not equal to
C U.
C
	      L=F_TO_S_MAPPING(I)    
	      U=F_TO_S_MAPPING(J)
C
C BRAT is the ration of b (level in full atom) to b (in the super level).
C Rates are written using BRAT since BRAT is assumed to remain constant
C during the linearization.
C
	      IF(L .NE. U)THEN
	        BRAT=(HN_F(I,K)/HN_S(L,K))*(HNST_S(L,K)/HNST_F(I,K))
!	        BRAT=(HN_F(I,K)/HNST_F(I,K))/(HN_S(L,K)/HNST_S(L,K))
	        CIJ=CIJ*BRAT*HNST_F(I,K)/HNST_S(L,K)*(W_F(J,K)/W_F(I,K))
C
	        BRAT=(HN_F(J,K)/HN_S(U,K))*(HNST_S(U,K)/HNST_F(J,K))
!	        BRAT=(HN_F(J,K)/HNST_F(J,K))/(HN_S(U,K)/HNST_S(U,K))
	        CJI=CJI*BRAT*HNST_F(J,K)/HNST_S(U,K)
C
	        COL_S(L,U,K)=COL_S(L,U,K)+CIJ
	        COL_S(U,L,K)=COL_S(U,L,K)+CJI
C
C We now compute the derivatives of the collision rates with respect to T.
C
	        DCOL_S(L,U,K)=DCOL_S(L,U,K)+CIJ*( dln_OMEGA_dlnT(I,J) + X -
	1           2.0D0 - HDKT*EDGE_F(I)/T(K)-dlnHNST_S_dlnT(L,K) )/T(K)
	        DCOL_S(U,L,K)=DCOL_S(U,L,K)+CJI*( dln_OMEGA_dlnT(I,J) - 
	1           2.0D0 - HDKT*EDGE_F(J)/T(K)-dlnHNST_S_dlnT(U,K) )/T(K)
C                                                                          
	        COOL(K)=COOL(K)+(HN_S(L,K)*CIJ-HN_S(U,K)*CJI)*
	1                          (EDGE_F(I)-EDGE_F(J))
	      END IF
	    END DO		!J
	  END DO		!I
C
C Now allow for collisional ionization. OMEGA for the collisional ionization
C is assumed to be define in exactly the same way as for collisional
C excitation. NB: The variation due to a change in X cancels out with the
C change in HNST_F.
C
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)    
	    CII=8.63D-08*ED(K)*OMEGA_F(I,I)/GHYD_F(I)/SQRT(T(K))
	    X=HDKT*EDGE_F(I)/T(K)
	    CII=CII*EXP(-X)*(HNST_F(I,K)/HNST_S(L,K))
	    COL_S(L,L,K)=COL_S(L,L,K)+CII
	    DCOL_S(L,L,K)=DCOL_S(L,L,K)+
	1       CII*( dln_OMEGA_dlnT(I,I) - 2.0D0 - dlnHNST_S_dlnT(L,K) )/T(K)
	    COOL(K)=COOL(K)+EDGE_F(I)*CII*(HN_S(L,K)-HNST_S(L,K))
	  END DO
C
	END DO
C
	RETURN
	END
