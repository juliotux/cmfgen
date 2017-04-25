!
! Routine to define a set of pale pens that can be used for coloring an area
! between two curves using the FILL option if PGPLOT.
!
	SUBROUTINE DEFINE_MORE_PENS(MAXPEN)
	IMPLICIT NONE
!
! Created 22-Apr-2015
!
	INTEGER MAXPEN
	INTEGER I
	REAL*4 R,G,B
!
	R=0.3; G=1.0; B=1.0; I=27
	CALL PGSCR(I,R,G,B)
	R=1.0; G=0.3; B=1.0; I=28
	CALL PGSCR(I,R,G,B)
	R=1.0; G=1.0; B=0.3; I=26
	CALL PGSCR(I,R,G,B)
!
	R=1.0; G=0.7; B=0.7; I=29
	CALL PGSCR(I,R,G,B)
	R=0.7; G=1.0; B=0.7; I=30
	CALL PGSCR(I,R,G,B)
	R=0.7; G=0.7; B=1.0; I=31
	CALL PGSCR(I,R,G,B)
!
	RETURN
	END
