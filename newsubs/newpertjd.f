C
C This routine is used to compute the perturbations to
C Jv (1 ... ND) as a function of the population levels and T .
C Uses Schuster or diffusion approximation for lower boundary
C condition. Subroutine may be used with or without a variable
C temperature. Note that F2DA , FC and FA must be stored sequentially
C for the guassian elimination routine.
C
	SUBROUTINE NEWPERTJD(F2DA,FC,FA,FB,VK,WM,AQW,
	1    DTAU,CHI,dCHIdR,R,Z,P,THETA,SOURCE,TA,TB,TC,XM,
	1    DIFF,DBB,IC,ESEC,THK,NC,ND,NP,METHOD)
	IMPLICIT NONE
C
C Altered 13-Mar-2000 : F2DA, FC, and FA no longer have to be sequential.
C
C Altered 28-Oct-96  Bug fix: COS converted back to ACOS in TOR expression.
C Altered 24-May-96  DOUBLE PRECISION declarations removed.
C                    CALL to DP_ZERO removed.
C                    IONE used for call to SIMPTH
C
C Altered 19-Jan-88  Method option installed to compute 'TAU'. Eta Vector
C                    removed from call. dCHIdR vector installed. We could
C                    pass dCHIdR, and not method, however we include
C                    method in case we later change VKIMD.
C Altered 9-Dec-86 - Thick boundary condition changed. Rays with only one or
C                    two points now treated. Vector PTK removed from call
C                    and from program since ETA term is now identical to
C                    electron scattering term in the boundary condition.
C Altered 20-FEB-86 (Replaced XM(1)=S1*(1.0D0-EXP(-TOR))/ETA(1) by
C                    XM(1)=(1.0D0-EXP(-TOR))/(ESEC(1)-CHI(1)) to advoid
C                    division by zero.)
C Altered 20-APR-84 (Thick approximation at outer boundary improved)
C Altered 22-AUG-82 (Thick approximation at outer boundary inserted)
C Altered 11-AUG-82 (form of diff vector altered to include CHI)
C Altered 5-AUG-82
C
	LOGICAL DIFF,THK
	INTEGER NC,ND,NP
	CHARACTER*(*) METHOD
C
	REAL*8 F2DA(ND,ND),FC(ND,ND),FB(ND,ND),VK(ND,ND)
	REAL*8 WM(ND,ND),AQW(ND,NP)
	REAL*8 DTAU(ND),CHI(ND),dCHIdR(ND)
	REAL*8 THETA(ND),SOURCE(ND),FA(ND),R(ND),Z(ND),P(NP)
	REAL*8 TA(ND),TB(ND),TC(ND),XM(ND),ESEC(ND)
	REAL*8 IC,DBB
!
	REAL*8 F2DAT(ND,2*ND+1)
!
C
	INTEGER, PARAMETER :: IONE=1
C
	INTEGER I,J,KS,NI,LS
	REAL*8 IBOUND,TOR,DBC
C
C Zero F2DA and FC matrices and FA vector. NB: These matrices should
C be stored sequentially,
C
	F2DA(:,:)=0.0D0
	FC(:,:)=0.0D0
C
C Compute dCHI/dR for use in computation of optical depth scale.
C Compute d[dCHI/dR]/dCHI for use in linearization.
C
	CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
	CALL d_DERIVCHI_dCHI(dCHIdR,CHI,R,ND,METHOD)
C
C Enter loop for each impact parameter P.
C
	DO 4000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	  END IF
C
	  IF(THK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796D0-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	  ELSE
	    TOR=0.0D0
	  END IF
C
	  IF(TOR .GT. 0.01D0)THEN
	    IBOUND=SOURCE(1)*(1.0D0-EXP(-TOR))
	  ELSE
	    IBOUND=SOURCE(1)*TOR*
	1             (1.D0-TOR/2.0D0*(1.0D0-TOR/3.0D0*(1.0D0-TOR/4.0D0)))
	  END IF
C
	  IF(NI .GT. 2)THEN
C
C Compute Z for this imapact parameter.
C
	    CALL ZALONGP(R,Z,P(LS),NI)
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
C
C Compute T ( A tridiagonal matrix) and store it as three vectors
C TA,TB and TC .
C
	    CALL TCOMPD(TA,TB,TC,DTAU,DIFF,LS,NC,ND,NI)
C
C Compute WM matrices.
C
	    CALL NEWWMAT(WM,DTAU,THETA,TOR,LS,NC,NI)
C
C Solve the tridiagonal system of equations for the A and B
C matrices. (Note after solution WXM contains both A and B)
C
	    CALL THOMAS(TA,TB,TC,WM,NI,NI)
C
C Compute U (The mean intensity along ray :-store in XM)
C
	    CALL XVECD(DTAU,SOURCE,XM,DIFF,DBC,IC,LS,NC,ND,NI)
	    XM(1)=-IBOUND
C
	    CALL SIMPTH(TA,TB,TC,XM,NI,IONE)
C
C Compute K matrix (multiplys variation of chi:-see notes).
C
	    CALL NEWVKIMD(VK,DTAU,CHI,SOURCE,XM,
	1                         R,Z,DIFF,DBC,LS,NC,ND,NI)
C
C Compute perturbation for THICK ! boundary condition. Note that
C VK(1,2) is not altered.
C
	    IF(THK)THEN
	      VK(1,1)=SOURCE(1)*(1.0D0-EXP(-TOR))/CHI(1)
	1   -   SOURCE(1)*TOR*EXP(-TOR)/CHI(1)+VK(1,1)
	    END IF
C
	    CALL SIMPTH(TA,TB,TC,VK,NI,NI)
C
C Compute &W vector if diffusion approximation.
C
	    IF(LS .LE. NC .AND. DIFF)THEN
	      DO I=1,ND
	        XM(I)=0.0D0
	      END DO
	      XM(ND)=SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	      CALL SIMPTH(TA,TB,TC,XM,NI,1)
	      CALL MULTVEC(FA,FA,XM,AQW(1,LS),NI)
	     END IF
C
	  ELSE
C
C XM is also computed, but this is not needed. Only require the WM
C matrix which is of identical form to the solution case.
C
	    CALL LAST2RAYS(XM,WM,R,Z,P(LS),DTAU,SOURCE,THETA,CHI,TOR,NI)
	    Call LAST2VKI(VK,SOURCE,CHI,DTAU,Z,TOR,NI)
C
	  END IF
C
C Update the F2DA and FC matrices . (see notes)
C
	  CALL MULT2D(F2DA,AQW(1,LS),VK,ND,NI,1)
	  CALL MULT2D(FC,AQW(1,LS),WM,ND,NI,1)
C
4000	CONTINUE
C
C Compute FB matrix [=FC+I] and dETA matrix.
C
	DO I=1,ND
	  DO J=1,ND
	    FB(J,I)=FC(J,I)
	    FC(J,I)=-FC(J,I)/ESEC(I)
	  END DO
	  FB(I,I)=FB(I,I)+1.0D0
	END DO
C
C Solve for mean intensity J as a function of depth.
C
	I=2*ND+1
!
	F2DAT(:,1:ND)=F2DA
	F2DAT(:,ND+11:2*ND)=FC
	F2DAT(:,2*ND+1)=FA
!
	CALL GAUSEL(FB,F2DAT,TA,ND,I,KS)
!
	F2DA=F2DAT(:,1:ND)
	FC=F2DAT(:,ND+1:2*ND)
	FA=F2DAT(:,2*ND+1)
C
	RETURN
	END
