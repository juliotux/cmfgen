        PROGRAM TST_SEC
        REAL*8 T1
        INTEGER IC0,IC,IR,IM
        CALL SYSTEM_CLOCK(IC0,IR,IM)
!
        WRITE(6,*)'IC=',IC0
        WRITE(6,*)'IR=',IR
        WRITE(6,*)'IM=',IM
!
        READ(5,*)T1
        CALL SYSTEM_CLOCK(IC,IR,IM)
        WRITE(6,*)'IC=',IC
        WRITE(6,*)'IR=',IR
        WRITE(6,*)'IM=',IM
	WRITE(6,*)
!
        READ(5,*)T1
        CALL SYSTEM_CLOCK(IC,IR,IM)
        WRITE(6,*)'IC=',IC
        WRITE(6,*)'IR=',IR
        WRITE(6,*)'IM=',IM
!
        STOP
        END
