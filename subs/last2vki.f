C
C Subroutine to compute the variation in intensity with opacity
C for rays with only one or two points.
C
	SUBROUTINE LAST2VKI(VK,SOURCE,CHI,DTAU,Z,TOR,NI)
	IMPLICIT NONE
C
C Altered 24-May-1996 -DEXP replaced bu EXP
C Altered 13-Apr-1988 - Bug fix.
C
	INTEGER NI
	REAL*8 VK(NI,NI),SOURCE(NI),CHI(NI),DTAU(NI),Z(NI)
	REAL*8 E1,E2,E3,DE1,DE2,DE3,TOR,T1,U2
C
C The incident intensity is assumed to be SOURCE(1)*(1.0-EXP(-TOR)).
C
	IF(TOR .GT. 0.01D0)THEN
	  T1=(1.0D0-EXP(-TOR))
	ELSE
	  T1=(1.D0-TOR/2.0D0*(1.0D0-TOR/3.0D0*(1.0-TOR/4.0D0)))*TOR
	END IF
C
	IF(NI .EQ. 1)THEN
	  VK(1,1)=SOURCE(1)*(EXP(-TOR)*TOR-T1)/CHI(1)
	ELSE IF(NI .EQ. 2)THEN
C
C The D denotes derivative with respect to CHI.
C
	  E1=EXP(-DTAU(1))
	  DE1=-E1*(Z(1)-Z(2))*0.5D0
	  IF(DTAU(1) .LT. 1.0D-03)THEN
	    E2=DTAU(1)*0.5D0+DTAU(1)*DTAU(1)/6.0D0
	    DE2=(0.5+DTAU(1)/3.0D0)*(Z(1)-Z(2))*0.5D0
	    E3=DTAU(1)*0.5D0-DTAU(1)*DTAU(1)/3.0D0
	    DE3=(0.5-2.0*DTAU(1)/3.0D0)*(Z(1)-Z(2))*0.5D0
	  ELSE
	    E2=1.0D0-(1.0D0-E1)/DTAU(1)
	    DE2=(1.0D0-E1-DTAU(1)*E1)/DTAU(1)/DTAU(1)*(Z(1)-Z(2))*0.5D0
	    E3=(1.0D0-E1)/DTAU(1)-E1
	    DE3=((E1-1.0D0+E1*DTAU(1))+E1*DTAU(1)*DTAU(1))*
	1                    (Z(1)-Z(2))*0.5D0/DTAU(1)/DTAU(1)
	  END IF
C
C Note that U(2)=S1.(E1.T1+E3) + S2.E2
C and  then U(1)=0.5*[S1.(T1+E2) + U(2).E1 + S2.E3]
C
	  VK(2,2)=SOURCE(2)*(DE2-E2/CHI(2))+SOURCE(1)*(DE3+DE1*T1)
          VK(2,1)=SOURCE(1)*(DE1*T1+DE3+
	1              (E1*TOR*EXP(-TOR)-E3-E1*T1)/CHI(1))
	1              +SOURCE(2)*DE2
C
	  U2=SOURCE(1)*(T1*E1+E3)+SOURCE(2)*E2
	  VK(1,1)=0.5D0*(  VK(2,1)*E1
	1                 +SOURCE(1)*( DE2-(E2+T1-TOR*EXP(-TOR))/CHI(1) )
	1                 +SOURCE(2)*DE3+U2*DE1  )
	  VK(1,2)=0.5D0*(  VK(2,2)*E1+SOURCE(2)*(DE3-E3/CHI(2))
	1                 +SOURCE(1)*DE2+U2*DE1  )
	END IF
C
	RETURN
	END
