C
C Subroutine to change the pens(color indeces) for the vectors
C 
C Created December 1996 by Gregson Vaux
C
      SUBROUTINE VECTORPEN(VECPEN,MAXVEC,FLAGLINE)
      USE GEN_IN_INTERFACE
      IMPLICIT NONE
C
      INTEGER MAXVEC
      INTEGER VECPEN(MAXVEC)
      LOGICAL FLAGLINE(MAXVEC)
!
      REAL*4 XBEG,XEND,YBEG,YEND
      INTEGER JUST,AXIS,NC
      INTEGER I,J
      CHARACTER*8 NUMB
      CHARACTER*14 MESSAGE,ICHAR*2
C
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
      DO I=1,15
        CALL PGNUMB(I,0,1,NUMB,NC)
        CALL PGSCI(1)
        CALL PGTEXT(0.0,(I*2.0),NUMB)
        CALL PGSCI(I)
        CALL PGTEXT(15.0,(I*2.0),'******----------******')
      END DO
      CALL PGSCI(1)
      CALL PGTEXT(15.0,0.0,'(Background)')
C
C
C
      DO I=1,MAXVEC
        IF(FLAGLINE(I))THEN
          J=VECPEN(I)
          WRITE(ICHAR,'(I2)')I
          MESSAGE='VECTOR#'//ICHAR//' PEN#'
          CALL GEN_IN(J,MESSAGE)
          VECPEN(I)=J
        END IF
      END DO
      RETURN
      END
