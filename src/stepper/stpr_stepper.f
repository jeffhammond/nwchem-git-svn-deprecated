      SUBROUTINE stpr_stepper(A,ILEFT,G,C,ETOT,NATD,CONVGE,CONVGG,
     &    converged, rtdb, step_number)
c $Id: stpr_stepper.f,v 1.3 1994-09-25 07:51:52 d3e129 Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      integer rtdb
      integer step_number
      logical converged
      DIMENSION G(NATD*3),A(ILEFT)
      DIMENSION C(3,NATD)
C
C     Calculate number of 64 bit words that core will need for
C     the calculation.
C
      CALL stpr_cneed (INEED, NATD)
C
C     Convergence tolerances.
C
      IF (INEED.GT.ILEFT) THEN
        WRITE(6,*)'  STEPPER| NOT ENOUGH MEMORY NEED ',INEED-ILEFT,
     *         ' MORE WORDS!'
        stop
      ELSE
        CALL stpr_stepcor(A,INEED,G,C,ETOT,NATD,CONVGE,CONVGG,converged,
     &      rtdb, step_number)
      ENDIF
      RETURN
      END
