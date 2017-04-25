C
C Subroutine to change the color indexes for the lines,
C background, and labels. In GRAMON, color indexes 
C are known as pens
C
C Created September 1996 by Gregson Vaux
C
      SUBROUTINE CHANGE_PEN(PENCOL,MAXPEN,NPLTS)
      USE NEW_GEN_IN_INTERFACE
      IMPLICIT NONE
C
      INTEGER MAXPEN,NPLTS
      INTEGER PENCOL(0:MAXPEN)
      REAL*4 XBEG,XEND,YBEG,YEND
      INTEGER JUST,AXIS,NC
      INTEGER I,J
      CHARACTER*8 NUMB
      CHARACTER*12 MESSAGE,ICHAR*2
      XBEG=0.0
      XEND=100.0
      YBEG=0.0
      YEND=34.0
      JUST=0
      AXIS=-2
      CALL PGENV(XBEG,XEND,YBEG,YEND,JUST,AXIS)
      CALL PGSCH(1.0)
      CALL PGERAS
      CALL PGTEXT(0.0,34.0,'PEN#')
      DO I=0,15
        CALL PGNUMB(I,0,1,NUMB,NC)
        CALL PGSCI(1)
        CALL PGTEXT(0.0,(I*2.0),NUMB)
        CALL PGSCI(I)
        CALL PGTEXT(15.0,(I*2.0),'******----------******')
      END DO
      CALL PGSCI(1)
      CALL PGTEXT(15.0,0.0,'(Background)')
C
      DO I=1,NPLTS
        J=PENCOL(I+1)
        WRITE(ICHAR,'(I2)')I
        MESSAGE='line#'//ICHAR//' pen#'
        CALL NEW_GEN_IN(J,MESSAGE)
        PENCOL(I+1)=J
      END DO
      RETURN
      END
