C
C Routine computes a new radius grid by inserting NPINS points
C between existing grid points. Two auxilary arrays and one matrix are
C returned.
C
C INDX indicates the old grid points on which the new variables
C will depend - namely INDX(I),INDX(I)+1,INDX(I)+2 and INDX(I)+3 where
C 'I' refers to a new grid point. The use of INDX simplifies handling
C of end points.
C
C COEF is a matrix giving the interpolation factors. The interpolation
C are in the X-log(R) plane (X may be logarithmic or linear variable).
C
C GRID indicates the position of the old grid points in the new
C array.
C
C DEEP is the depth index beyond which we use cubic interpolation (in the
C LOG plane) so that we may satisfy the diffusion approximation.
C
	SUBROUTINE REXT_COEF_V2(REXT,COEF,INDX,NX,R,GRID,ND,NPINS,FLAG,
	1              DEEP,ST_INDX,END_INDX)
	IMPLICIT NONE
C
C Altered 31-Jan-2010 - Changed insertion at inner bundary to keep (R(ND-1)-R(ND))/(R(ND-2)-R(ND-1)) small.
C Altered 12-Jun-2009 - Changed insertion at outer bundary to keep (R1-R2)/(R2-R3) small.
C Altered 05-Jan-1998 - NEND replaced by ST_INDX, END_INDX
C                         Changed to V2.
C Altered 28-May-1996 - Call to DP_ZERO removed.
C                       ERROR_LU installed.
C                       Generical calls for EXP and LOG.
C
C Altered  1-Oct-1988 - NEND installed. Interpolation is perfomed only
C                          from outer boundary to NEND.
C Altered 11-May-1988 - Cubic interpolation at depth installed.
C Created  5-Apr-1988 - Based on INTERPTHREE, and EXTENDVTSR
C
	INTEGER NX,ND,NPINS,ST_INDX,END_INDX
	INTEGER INDX(NX),GRID(ND)
	REAL*8 REXT(NX),COEF(4,NX),R(ND)
	LOGICAL FLAG
C
	INTEGER I,J,K,M,DEEP
	REAL*8 DELR,T1,T2,A1,A2,A3,A4
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	IF(R(1) .LT. R(ND))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Error in REXT_COEF - X array is not a decreasing'
	1,   ' function'
	  STOP
	END IF
C
C Compute the new R array if desired.
C
	IF(FLAG)THEN
	  DO I=1,ST_INDX-1
	    REXT(I)=R(I)
	    GRID(I)=I
	  END DO
C
	  DO I=ST_INDX,END_INDX-1
	    DELR=(LOG(R(I+1))-LOG(R(I)))/(NPINS+1)
	    J=I+(I-ST_INDX)*NPINS
	    REXT(J)=R(I)
	    GRID(I)=J
	    DO K=1,NPINS
	      REXT(J+K)=R(I)*EXP(DELR*K)
	    END DO
	  END DO
!
! Finalized: 12-June-2009
! This technique preserves the ratio of the last two grid spacings.
!
	  IF(ST_INDX .EQ. 1)THEN
	    T1=(R(1)-R(2))/(R(2)-R(3))
	    IF(T1 .LT. 0.1D0)THEN
	      REXT(2)=(REXT(1)+T1*REXT(3))/(1.0D0+T1)
	    ELSE
	      REXT(2)=REXT(1)-0.1D0*(REXT(1)-REXT(2))
	    END IF
	    WRITE(6,*)' Information about insertion of extra grid points at outer boundary (INC_GRID option)'
	    WRITE(6,*)'   R ratio:',T1,(R(2)-R(3))/(R(3)-R(4))
	    WRITE(6,*)'REXT ratio:',(REXT(1)-REXT(2))/(REXT(2)-REXT(3)),(REXT(2)-REXT(3))/(REXT(3)-REXT(4))
	  END IF
C
	  DO I=END_INDX,ND
	    J=I+(END_INDX-ST_INDX)*NPINS
	    REXT(J)=R(I)
	    GRID(I)=J
	  END DO
C
	  IF(END_INDX .EQ. ND)THEN
	    T1=(R(ND-1)-R(ND))/(R(ND-2)-R(ND-1))
	    T1=MIN(T1,0.1D0)
	    REXT(NX-1)=(REXT(NX)+T1*REXT(NX-2))/(1.0D0+T1)
!	    WRITE(6,*)' Information about insertion of extra grid points at inner boundary (INC_GRID option)'
!	    WRITE(6,*)'   R ratio:',T1,(R(ND-1)-R(ND))/(R(ND-2)-R(ND-1))
!	    WRITE(6,*)'REXT ratio:',(REXT(NX-1)-REXT(NX))/(REXT(NX-2)-REXT(NX-1)),
!	1                           (REXT(NX-2)-REXT(NX-1))/(REXT(NX-3)-REXT(NX-2))
!	    WRITE(6,*)REXT(NX-5:NX)
!	    FLUSH(UNIT=6)
	  END IF
C
	ELSE
C
C Check that REXT and GRID have been computed
C
	  IF(GRID(1) .NE. 1 .OR. GRID(ND) .NE. NX
	1      .OR. R(1) .NE. REXT(1) .OR. R(ND) .NE. REXT(NX) )THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Error in REEXT_COEF - GRID or REXT are incorrect'
	     STOP
	  END IF
	END IF
C
C Now compute the INDX vector, and the COEF matrix.
C
	COEF(:,:)=0.D0      	!4:NX
C
	INDX(1)=1
	COEF(1,1)=1.0D0
	INDX(NX)=ND-3
	COEF(4,NX)=1.0
C	
	M=2
	DO 100 I=2,NX-1
500	  IF(REXT(I) .GE. R(M))THEN
	    IF(M .EQ. 2)THEN
	      T1=(LOG(REXT(I))-LOG(R(1)))/(LOG(R(2))-LOG(R(1)))
	      INDX(I)=1
	      COEF(1,I)=1.0D0-T1
	      COEF(2,I)=T1
	    ELSE IF(M .EQ. ND)THEN
	      A1=LOG(R(M-3))
	      A2=LOG(R(M-2))
	      A3=LOG(R(M-1))
	      A4=LOG(R(M))
	      T1=LOG(REXT(I))
	      COEF(1,I)=(T1-A2)*(T1-A3)*(T1-A4)/(A1-A2)/(A1-A3)/(A1-A4)
	      COEF(2,I)=(T1-A1)*(T1-A3)*(T1-A4)/(A2-A1)/(A2-A3)/(A2-A4)
	      COEF(3,I)=(T1-A1)*(T1-A2)*(T1-A4)/(A3-A1)/(A3-A2)/(A3-A4)
	      COEF(4,I)=(T1-A1)*(T1-A2)*(T1-A3)/(A4-A1)/(A4-A2)/(A4-A3)
	      INDX(I)=M-3
	    ELSE IF(M .GT. DEEP)THEN
	      A1=LOG(R(M-2))
	      A2=LOG(R(M-1))
	      A3=LOG(R(M))
	      A4=LOG(R(M+1))
	      T1=LOG(REXT(I))
	      COEF(1,I)=(T1-A2)*(T1-A3)*(T1-A4)/(A1-A2)/(A1-A3)/(A1-A4)
	      COEF(2,I)=(T1-A1)*(T1-A3)*(T1-A4)/(A2-A1)/(A2-A3)/(A2-A4)
	      COEF(3,I)=(T1-A1)*(T1-A2)*(T1-A4)/(A3-A1)/(A3-A2)/(A3-A4)
	      COEF(4,I)=(T1-A1)*(T1-A2)*(T1-A3)/(A4-A1)/(A4-A2)/(A4-A3)
	      INDX(I)=M-2
	    ELSE
	      T1=LOG(REXT(I))-LOG(R(M-1))
	      T2=T1/(LOG(R(M))-LOG(R(M-1)))
	      A2=T2*T2*(3.0-2*T2)
	      A1=1.0-A2
	      A3=T1*(1.0-T2*(2.0-T2))/(LOG(R(M))-LOG(R(M-2)))
	      A4=T1*T2*(T2-1.0)/(LOG(R(M+1))-LOG(R(M-1)))
	      INDX(I)=M-2
	      COEF(1,I)=-A3
	      COEF(2,I)=A1-A4
	      COEF(3,I)=A2+A3
	      COEF(4,I)=A4
	    END IF
	  ELSE
	    M=M+1
	    GOTO 500
	  END IF
C
100	CONTINUE
C
	RETURN
	END
