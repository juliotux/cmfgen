	PROGRAM TST_SET
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: N=51
	INTEGER I
	REAL*8 X0,X1
	REAL*8 X(N)
	REAL*8 dX(N)
!
	X0=2.5
	X1=112
	call set_xkt_array(x0,x1,n,x,dx,'lin')
!
	do i=1,n
	  write(100,*)i,x(i),dx(i)
	end do
!
	STOP
	END
