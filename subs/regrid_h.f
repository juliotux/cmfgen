C
C Routine to compute H on the radius grid. As input, the routine uses RSQHNU
C computed at the mid point of the radius grid.
C
	SUBROUTINE REGRID_H(HNU,R,RSQHNU,H_OUT,H_IN,ND,MIDR)
	IMPLICIT NONE
C
C Altered  : 29-MAy-1996 : Bug fix. Incorrect dimension passed to MON_INTERP.
C                            No effect on computation.
C Finalized 05-Jan-1995
C
	INTEGER ND
	REAL*8 HNU(ND)		!Returned - Flux on nodes
C
	REAL*8 R(ND)
	REAL*8 RSQHNU(ND-1)
C
	REAL*8 H_OUT		!H flux at outer boundary (from boundary 
	                        !                          conditions)
	REAL*8 H_IN             !H flux at inner boundary.
C
	REAL*8 MIDR(ND)		!Work array
C
C Local variables
C
	INTEGER ND_M1,ND_M2,I,IONE
	PARAMETER (IONE=1)
C
C Compute points of radius grid. RSQHNU is defined on this grid.
C
	DO I=1,ND-1
	  MIDR(I)=0.5D0*(R(I)+R(I+1))
	END DO
	ND_M1=ND-1
	ND_M2=ND-2
C
C Regrid RSQHNU from the midpoints onto the NODES of the RADIUS grid using
C cubic interpolation but with forced monoticity.
C                              
C The boundary fluxes are set outside the routine, since  R(1) and R(ND) are
C not contained within MIDR.
C
C IONE refers to the number of variables, which in this case is one.
C
	CALL MON_INTERP(HNU(2),ND_M2,IONE,R(2),ND_M2,
	1                 RSQHNU,ND_M1,MIDR,ND_M1)
	HNU(1)=H_OUT*R(1)*R(1)
	HNU(ND)=H_IN*R(ND)*R(ND)
C
	RETURN
	END
