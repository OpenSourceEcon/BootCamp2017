MODULE solver

! Corrections to FUNCTION Enorm - 28 November 2003

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = 8
PRIVATE
PUBLIC  :: hbrd, hybrd


CONTAINS


SUBROUTINE hbrd(fcn, n, x, fvec, epsfcn, tol, info, diag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-07-15  Time: 13:27:42

INTEGER, INTENT(IN)        :: n
REAL (8), INTENT(IN OUT)  :: x(n)
REAL (8), INTENT(IN OUT)  :: fvec(n)
REAL (8), INTENT(IN)      :: epsfcn
REAL (8), INTENT(IN)      :: tol
INTEGER, INTENT(OUT)       :: info
REAL (8), INTENT(OUT)     :: diag(n)

! EXTERNAL fcn
INTERFACE
  SUBROUTINE FCN(N, X, FVEC, IFLAG)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = 8
    INTEGER, INTENT(IN)      :: n
    REAL (8), INTENT(IN)    :: x(n)
    REAL (8), INTENT(OUT)   :: fvec(n)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE FCN
END INTERFACE

!   **********

!   SUBROUTINE HBRD

!     THE PURPOSE OF HBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
!   FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
!   THIS IS DONE BY USING THE MORE GENERAL NONLINEAR EQUATION SOLVER HYBRD.
!   THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE FUNCTIONS.
!   THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE HBRD(N, X, FVEC, EPSFCN, TOL, INFO, WA, LWA)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT
!       IN THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE FCN(N, X, FVEC, IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND RETURN THIS VECTOR IN FVEC.
!       ---------
!       RETURN
!       END

!       THE VALUE OF IFLAG NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE THE EXECUTION OF HBRD.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN AN INITIAL
!       ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE
!       FINAL ESTIMATE OF THE SOLUTION VECTOR.

!     FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE FUNCTIONS EVALUATED AT THE OUTPUT X.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
!       FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION ASSUMES
!       THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF EPSFCN.
!       IF EPSFCN IS LESS THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE
!       RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
!       PRECISION.

!     TOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
!       ALGORITHM ESTIMATES THAT THE RELATIVE ERROR BETWEEN X AND THE SOLUTION
!       IS AT MOST TOL.

!     INFO IS AN INTEGER OUTPUT VARIABLE.  IF THE USER HAS TERMINATED
!       EXECUTION, INFO IS SET TO THE (NEGATIVE) VALUE OF IFLAG.
!       SEE DESCRIPTION OF FCN.  OTHERWISE, INFO IS SET AS FOLLOWS.

!       INFO = 0   IMPROPER INPUT PARAMETERS.

!       INFO = 1   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!                  BETWEEN X AND THE SOLUTION IS AT MOST TOL.

!       INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED 200*(N+1).

!       INFO = 3   TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                  THE APPROXIMATE SOLUTION X IS POSSIBLE.

!       INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS.

!   SUBPROGRAMS CALLED

!     USER-SUPPLIED ...... FCN

!     MINPACK-SUPPLIED ... HYBRD

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

! Reference:
! Powell, M.J.D. 'A hybrid method for nonlinear equations' in Numerical Methods
!      for Nonlinear Algebraic Equations', P.Rabinowitz (editor), Gordon and
!      Breach, London 1970.
!   **********
INTEGER    :: maxfev, ml, mode, mu, nfev, nprint
REAL (8)  :: xtol
REAL (8), PARAMETER  :: factor = 100.0_dp, zero = 0.0_dp

info = 0

!     CHECK THE INPUT PARAMETERS FOR ERRORS.


IF (n <= 0 .OR. epsfcn < zero .OR. tol < zero) GO TO 20

!     CALL HYBRD.

maxfev = 200*(n + 1)
xtol = tol
ml = n - 1
mu = n - 1
mode = 2
nprint = 0
CALL hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,  &
           factor, nprint, info, nfev)
IF (info == 5) info = 4
20 RETURN

!     LAST CARD OF SUBROUTINE HBRD.

END SUBROUTINE hbrd



SUBROUTINE hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,  &
                 factor, nprint, info, nfev)

INTEGER, INTENT(IN)        :: n
REAL (8), INTENT(IN OUT)  :: x(n)
REAL (8), INTENT(IN OUT)  :: fvec(n)
REAL (8), INTENT(IN)      :: xtol
INTEGER, INTENT(IN OUT)    :: maxfev
INTEGER, INTENT(IN OUT)    :: ml
INTEGER, INTENT(IN)        :: mu
REAL (8), INTENT(IN)      :: epsfcn
REAL (8), INTENT(OUT)     :: diag(n)
INTEGER, INTENT(IN)        :: mode
REAL (8), INTENT(IN)      :: factor
INTEGER, INTENT(IN OUT)    :: nprint
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: nfev

! EXTERNAL fcn
INTERFACE
  SUBROUTINE FCN(N, X, FVEC, IFLAG)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = 8
    INTEGER, INTENT(IN)      :: n
    REAL (8), INTENT(IN)    :: x(n)
    REAL (8), INTENT(OUT)   :: fvec(n)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE FCN
END INTERFACE

!   **********

!   SUBROUTINE HYBRD

!   THE PURPOSE OF HYBRD IS TO FIND A ZERO OF A SYSTEM OF N NONLINEAR
!   FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE POWELL HYBRID METHOD.
!   THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE FUNCTIONS.
!   THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE HYBRD(FCN, N, X, FVEC, XTOL, MAXFEV, ML, MU, EPSFCN,
!                      DIAG, MODE, FACTOR, NPRINT, INFO, NFEV, FJAC,
!                      LDFJAC, R, LR, QTF, WA1, WA2, WA3, WA4)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
!       THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE FCN(N, X, FVEC, IFLAG)
!       INTEGER N, IFLAG
!       REAL X(N), FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND
!       RETURN THIS VECTOR IN FVEC.
!       ---------
!       RETURN
!       END

!       THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE EXECUTION OF HYBRD.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN ARRAY OF LENGTH N.  ON INPUT X MUST CONTAIN AN INITIAL
!       ESTIMATE OF THE SOLUTION VECTOR.  ON OUTPUT X CONTAINS THE FINAL
!       ESTIMATE OF THE SOLUTION VECTOR.

!     FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE FUNCTIONS EVALUATED AT THE OUTPUT X.

!     XTOL IS A NONNEGATIVE INPUT VARIABLE.  TERMINATION OCCURS WHEN THE
!       RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES IS AT MOST XTOL.

!     MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE.  TERMINATION OCCURS WHEN
!       THE NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN
!       ITERATION.

!     ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE
!       NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
!       IF THE JACOBIAN IS NOT BANDED, SET ML TO AT LEAST N - 1.

!     MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES THE NUMBER
!       OF SUPERDIAGONALS WITHIN THE BAND OF THE JACOBIAN MATRIX.
!       IF THE JACOBIAN IS NOT BANDED, SET MU TO AT LEAST N - 1.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE STEP LENGTH
!       FOR THE FORWARD-DIFFERENCE APPROXIMATION.  THIS APPROXIMATION
!       ASSUMES THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE ORDER
!       OF EPSFCN. IF EPSFCN IS LESS THAN THE MACHINE PRECISION,
!       IT IS ASSUMED THAT THE RELATIVE ERRORS IN THE FUNCTIONS ARE OF THE
!       ORDER OF THE MACHINE PRECISION.

!     DIAG IS AN ARRAY OF LENGTH N. IF MODE = 1 (SEE BELOW),
!       DIAG IS INTERNALLY SET.  IF MODE = 2, DIAG MUST CONTAIN POSITIVE
!       ENTRIES THAT SERVE AS MULTIPLICATIVE SCALE FACTORS FOR THE
!       VARIABLES.

!     MODE IS AN INTEGER INPUT VARIABLE. IF MODE = 1, THE VARIABLES WILL BE
!       SCALED INTERNALLY.  IF MODE = 2, THE SCALING IS SPECIFIED BY THE
!       INPUT DIAG.  OTHER VALUES OF MODE ARE EQUIVALENT TO MODE = 1.

!     FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
!       INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
!       FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
!       TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
!       INTERVAL (.1,100.). 100. IS A GENERALLY RECOMMENDED VALUE.

!     NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
!       PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
!       FCN IS CALLED WITH IFLAG = 0 AT THE BEGINNING OF THE FIRST
!       ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
!       IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
!       FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
!       OF FCN WITH IFLAG = 0 ARE MADE.

!     INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
!       TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
!       VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
!       INFO IS SET AS FOLLOWS.

!       INFO = 0   IMPROPER INPUT PARAMETERS.

!       INFO = 1   RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
!                  IS AT MOST XTOL.

!       INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED MAXFEV.

!       INFO = 3   XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
!                  THE APPROXIMATE SOLUTION X IS POSSIBLE.

!       INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS, AS
!                  MEASURED BY THE IMPROVEMENT FROM THE LAST
!                  FIVE JACOBIAN EVALUATIONS.

!       INFO = 5   ITERATION IS NOT MAKING GOOD PROGRESS, AS MEASURED BY
!                  THE IMPROVEMENT FROM THE LAST TEN ITERATIONS.

!     NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

!     FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE ORTHOGONAL MATRIX Q
!       PRODUCED BY THE QR FACTORIZATION OF THE FINAL APPROXIMATE JACOBIAN.

!     LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!     R IS AN OUTPUT ARRAY OF LENGTH LR WHICH CONTAINS THE
!       UPPER TRIANGULAR MATRIX PRODUCED BY THE QR FACTORIZATION
!       OF THE FINAL APPROXIMATE JACOBIAN, STORED ROWWISE.

!     LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN (N*(N+1))/2.

!     QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
!       THE VECTOR (Q TRANSPOSE)*FVEC.

!     WA1, WA2, WA3, AND WA4 ARE WORK ARRAYS OF LENGTH N.

!   SUBPROGRAMS CALLED

!     USER-SUPPLIED ...... FCN

!     MINPACK-SUPPLIED ... DOGLEG,SPMPAR,ENORM,FDJAC1,
!                          QFORM,QRFAC,R1MPYQ,R1UPDT

!     FORTRAN-SUPPLIED ... ABS,MAX,MIN,MIN,MOD

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********

INTEGER    :: i, iflag, iter, j, jm1, l, lr, msum, ncfail, ncsuc, nslow1,  &
              nslow2
INTEGER    :: iwa(1)
LOGICAL    :: jeval, sing
REAL (8)  :: actred, delta, epsmch, fnorm, fnorm1, pnorm, prered,   &
              ratio, sum, temp, xnorm
REAL (8), PARAMETER  :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp,   &
                         p001 = 0.001_dp, p0001 = 0.0001_dp, zero = 0.0_dp

! The following were workspace arguments
REAL (8)  :: fjac(n,n), r(n*(n+1)/2), qtf(n), wa1(n), wa2(n),   &
              wa3(n), wa4(n)

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

info = 0
iflag = 0
nfev = 0
lr = n*(n+1)/2

!     CHECK THE INPUT PARAMETERS FOR ERRORS.

IF (n > 0 .AND. xtol >= zero .AND. maxfev > 0 .AND. ml >= 0 .AND. mu >=  &
    0 .AND. factor > zero ) THEN
IF (mode == 2) THEN
  diag(1:n) = one
END IF

!     EVALUATE THE FUNCTION AT THE STARTING POINT AND CALCULATE ITS NORM.

iflag = 1
CALL fcn(n, x, fvec, iflag)
nfev = 1
IF (iflag >= 0) THEN
  fnorm = enorm(n, fvec)
  
!   DETERMINE THE NUMBER OF CALLS TO FCN NEEDED TO COMPUTE THE JACOBIAN MATRIX.
  
  msum = MIN(ml+mu+1,n)
  
!     INITIALIZE ITERATION COUNTER AND MONITORS.
  
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
  
!     BEGINNING OF THE OUTER LOOP.
  
  20 jeval = .true.
  
!        CALCULATE THE JACOBIAN MATRIX.
  
  iflag = 2
  CALL fdjac1(fcn, n, x, fvec, fjac, n, iflag, ml, mu, epsfcn, wa1, wa2)
  nfev = nfev + msum
  IF (iflag >= 0) THEN
    
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
    
    CALL qrfac(n, n, fjac, n, .false., iwa, 1, wa1, wa2, wa3)
    
!        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
    
    IF (iter == 1) THEN
      IF (mode /= 2) THEN
        DO  j = 1, n
          diag(j) = wa2(j)
          IF (wa2(j) == zero) diag(j) = one
        END DO
      END IF
      
!        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
!        AND INITIALIZE THE STEP BOUND DELTA.
      
      wa3(1:n) = diag(1:n) * x(1:n)
      xnorm = enorm(n, wa3)
      delta = factor * xnorm
      IF (delta == zero) delta = factor
    END IF
    
!        FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
    
    qtf(1:n) = fvec(1:n)
    DO  j = 1, n
      IF (fjac(j,j) /= zero) THEN
        sum = zero
        DO  i = j, n
          sum = sum + fjac(i,j) * qtf(i)
        END DO
        temp = -sum / fjac(j,j)
        DO  i = j, n
          qtf(i) = qtf(i) + fjac(i,j) * temp
        END DO
      END IF
    END DO
    
!        COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
    
    sing = .false.
    DO  j = 1, n
      l = j
      jm1 = j - 1
      IF (jm1 >= 1) THEN
        DO  i = 1, jm1
          r(l) = fjac(i,j)
          l = l + n - i
        END DO
      END IF
      r(l) = wa1(j)
      IF (wa1(j) == zero) sing = .true.
    END DO
    
!        ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
    
    CALL qform(n, n, fjac, n, wa1)
    
!        RESCALE IF NECESSARY.
    
    IF (mode /= 2) THEN
      DO  j = 1, n
        diag(j) = MAX(diag(j), wa2(j))
      END DO
    END IF
    
!        BEGINNING OF THE INNER LOOP.
    
!           IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
    
    120 IF (nprint > 0) THEN
      iflag = 0
      IF (MOD(iter-1, nprint) == 0) CALL fcn(n, x, fvec, iflag)
      IF (iflag < 0) GO TO 190
    END IF
    
!           DETERMINE THE DIRECTION P.
    
    CALL dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3)
    
!           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
    
    DO  j = 1, n
      wa1(j) = -wa1(j)
      wa2(j) = x(j) + wa1(j)
      wa3(j) = diag(j) * wa1(j)
    END DO
    pnorm = enorm(n, wa3)
    
!           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
    
    IF (iter == 1) delta = MIN(delta, pnorm)
    
!           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
    
    iflag = 1
    CALL fcn(n, wa2, wa4, iflag)
    nfev = nfev + 1
    IF (iflag >= 0) THEN
      fnorm1 = enorm(n, wa4)
      
!           COMPUTE THE SCALED ACTUAL REDUCTION.
      
      actred = -one
      IF (fnorm1 < fnorm) actred = one - (fnorm1/fnorm) ** 2
      
!           COMPUTE THE SCALED PREDICTED REDUCTION.
      
      l = 1
      DO  i = 1, n
        sum = zero
        DO  j = i, n
          sum = sum + r(l) * wa1(j)
          l = l + 1
        END DO
        wa3(i) = qtf(i) + sum
      END DO
      temp = enorm(n, wa3)
      prered = zero
      IF (temp < fnorm) prered = one - (temp/fnorm) ** 2
      
!           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED REDUCTION.
      
      ratio = zero
      IF (prered > zero) ratio = actred / prered
      
!           UPDATE THE STEP BOUND.
      
      IF (ratio < p1) THEN
        ncsuc = 0
        ncfail = ncfail + 1
        delta = p5 * delta
      ELSE
        ncfail = 0
        ncsuc = ncsuc + 1
        IF (ratio >= p5 .OR. ncsuc > 1) delta = MAX(delta,pnorm/p5)
        IF (ABS(ratio-one) <= p1) delta = pnorm / p5
      END IF
      
!           TEST FOR SUCCESSFUL ITERATION.
      
      IF (ratio >= p0001) THEN
        
!           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
        
        DO  j = 1, n
          x(j) = wa2(j)
          wa2(j) = diag(j) * x(j)
          fvec(j) = wa4(j)
        END DO
        xnorm = enorm(n, wa2)
        fnorm = fnorm1
        iter = iter + 1
      END IF
      
!           DETERMINE THE PROGRESS OF THE ITERATION.
      
      nslow1 = nslow1 + 1
      IF (actred >= p001) nslow1 = 0
      IF (jeval) nslow2 = nslow2 + 1
      IF (actred >= p1) nslow2 = 0
      
!           TEST FOR CONVERGENCE.
      
      IF (delta <= xtol*xnorm .OR. fnorm == zero) info = 1
      IF (info == 0) THEN
        
!           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
        
        IF (nfev >= maxfev) info = 2
        IF (p1*MAX(p1*delta, pnorm) <= epsmch*xnorm) info = 3
        IF (nslow2 == 5) info = 4
        IF (nslow1 == 10) info = 5
        IF (info == 0) THEN
          
!           CRITERION FOR RECALCULATING JACOBIAN APPROXIMATION
!           BY FORWARD DIFFERENCES.
          
          IF (ncfail /= 2) THEN
            
!           CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
!           AND UPDATE QTF IF NECESSARY.
            
            DO  j = 1, n
              sum = zero
              DO  i = 1, n
                sum = sum + fjac(i,j) * wa4(i)
              END DO
              wa2(j) = (sum-wa3(j)) / pnorm
              wa1(j) = diag(j) * ((diag(j)*wa1(j))/pnorm)
              IF (ratio >= p0001) qtf(j) = sum
            END DO
            
!           COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
            
            CALL r1updt(n, n, r, lr, wa1, wa2, wa3, sing)
            CALL r1mpyq(n, n, fjac, n, wa2, wa3)
            CALL r1mpyq(1, n, qtf, 1, wa2, wa3)
            
!           END OF THE INNER LOOP.
            
            jeval = .false.
            GO TO 120
          END IF
          
!        END OF THE OUTER LOOP.
          
          GO TO 20
        END IF
      END IF
    END IF
  END IF
END IF
END IF

!     TERMINATION, EITHER NORMAL OR USER IMPOSED.

190 IF (iflag < 0) info = iflag
iflag = 0
IF (nprint > 0) CALL fcn(n, x, fvec, iflag)
RETURN

!     LAST CARD OF SUBROUTINE HYBRD.

END SUBROUTINE hybrd



SUBROUTINE dogleg(n, r, lr, diag, qtb, delta, x, wa1, wa2)

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lr
REAL (8), INTENT(IN)      :: r(lr)
REAL (8), INTENT(IN)      :: diag(n)
REAL (8), INTENT(IN)      :: qtb(n)
REAL (8), INTENT(IN)      :: delta
REAL (8), INTENT(IN OUT)  :: x(n)
REAL (8), INTENT(OUT)     :: wa1(n)
REAL (8), INTENT(OUT)     :: wa2(n)


!     **********

!     SUBROUTINE DOGLEG

!     GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
!     MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA, THE
!     PROBLEM IS TO DETERMINE THE CONVEX COMBINATION X OF THE
!     GAUSS-NEWTON AND SCALED GRADIENT DIRECTIONS THAT MINIMIZES
!     (A*X - B) IN THE LEAST SQUARES SENSE, SUBJECT TO THE
!     RESTRICTION THAT THE EUCLIDEAN NORM OF D*X BE AT MOST DELTA.

!     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
!     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
!     QR FACTORIZATION OF A. THAT IS, IF A = Q*R, WHERE Q HAS
!     ORTHOGONAL COLUMNS AND R IS AN UPPER TRIANGULAR MATRIX,
!     THEN DOGLEG EXPECTS THE FULL UPPER TRIANGLE OF R AND
!     THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B.

!     THE SUBROUTINE STATEMENT IS

!       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)

!     WHERE

!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.

!       R IS AN INPUT ARRAY OF LENGTH LR WHICH MUST CONTAIN THE UPPER
!         TRIANGULAR MATRIX R STORED BY ROWS.

!       LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!         (N*(N+1))/2.

!       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!         DIAGONAL ELEMENTS OF THE MATRIX D.

!       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
!         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.

!       DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
!         BOUND ON THE EUCLIDEAN NORM OF D*X.

!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE DESIRED
!         CONVEX COMBINATION OF THE GAUSS-NEWTON DIRECTION AND THE
!         SCALED GRADIENT DIRECTION.

!       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.

!     SUBPROGRAMS CALLED

!       MINPACK-SUPPLIED ... SPMPAR,ENORM

!       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT

!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!     **********
INTEGER    :: i, j, jj, jp1, k, l
REAL (8)  :: alpha, bnorm, epsmch, gnorm, qnorm, sgnorm, sum, temp

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

!     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.

jj = (n*(n+1)) / 2 + 1
DO  k = 1, n
  j = n - k + 1
  jp1 = j + 1
  jj = jj - k
  l = jj + 1
  sum = 0.0
  IF (n >= jp1) THEN
    DO  i = jp1, n
      sum = sum + r(l) * x(i)
      l = l + 1
    END DO
  END IF
  temp = r(jj)
  IF (temp == 0.0) THEN
    l = j
    DO  i = 1, j
      temp = MAX(temp,ABS(r(l)))
      l = l + n - i
    END DO
    temp = epsmch * temp
    IF (temp == 0.0) temp = epsmch
  END IF
  x(j) = (qtb(j)-sum) / temp
END DO

!     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.

DO  j = 1, n
  wa1(j) = 0.0
  wa2(j) = diag(j) * x(j)
END DO
qnorm = enorm(n, wa2)
IF (qnorm > delta) THEN
  
!     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
!     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
  
  l = 1
  DO  j = 1, n
    temp = qtb(j)
    DO  i = j, n
      wa1(i) = wa1(i) + r(l) * temp
      l = l + 1
    END DO
    wa1(j) = wa1(j) / diag(j)
  END DO
  
!     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
!     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
  
  gnorm = enorm(n, wa1)
  sgnorm = 0.0
  alpha = delta / qnorm
  IF (gnorm /= 0.0) THEN
    
!     CALCULATE THE POINT ALONG THE SCALED GRADIENT
!     AT WHICH THE QUADRATIC IS MINIMIZED.
    
    DO  j = 1, n
      wa1(j) = (wa1(j)/gnorm) / diag(j)
    END DO
    l = 1
    DO  j = 1, n
      sum = 0.0
      DO  i = j, n
        sum = sum + r(l) * wa1(i)
        l = l + 1
      END DO
      wa2(j) = sum
    END DO
    temp = enorm(n, wa2)
    sgnorm = (gnorm/temp) / temp
    
!     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
    
    alpha = 0.0
    IF (sgnorm < delta) THEN
      
!     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
!     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
!     AT WHICH THE QUADRATIC IS MINIMIZED.
      
      bnorm = enorm(n, qtb)
      temp = (bnorm/gnorm) * (bnorm/qnorm) * (sgnorm/delta)
      temp = temp - (delta/qnorm) * (sgnorm/delta) ** 2 + SQRT((  &
          temp-(delta/qnorm))**2+(1.0-(delta/qnorm)**2)*(1.0-( sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(1.0-(sgnorm/delta)**2)) / temp
    END IF
  END IF
  
!     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
!     DIRECTION AND THE SCALED GRADIENT DIRECTION.
  
  temp = (1.0-alpha) * MIN(sgnorm,delta)
  DO  j = 1, n
    x(j) = temp * wa1(j) + alpha * x(j)
  END DO
END IF
RETURN
END SUBROUTINE dogleg


SUBROUTINE fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn,   &
                  wa1, wa2)

INTEGER, INTENT(IN)        :: n
REAL (8), INTENT(IN OUT)  :: x(n)
REAL (8), INTENT(IN)      :: fvec(n)
INTEGER, INTENT(IN)        :: ldfjac
REAL (8), INTENT(OUT)     :: fjac(ldfjac,n)
INTEGER, INTENT(IN OUT)    :: iflag
INTEGER, INTENT(IN)        :: ml
INTEGER, INTENT(IN)        :: mu
REAL (8), INTENT(IN)      :: epsfcn
REAL (8), INTENT(IN OUT)  :: wa1(n)
REAL (8), INTENT(OUT)     :: wa2(n)

! EXTERNAL fcn
INTERFACE
  SUBROUTINE FCN(N, X, FVEC, IFLAG)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = 8
    INTEGER, INTENT(IN)      :: n
    REAL (8), INTENT(IN)    :: x(n)
    REAL (8), INTENT(OUT)   :: fvec(n)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE FCN
END INTERFACE

!   **********

!   SUBROUTINE FDJAC1

!   THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION TO THE N BY N
!   JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED PROBLEM OF N FUNCTIONS IN N
!   VARIABLES.  IF THE JACOBIAN HAS A BANDED FORM, THEN FUNCTION EVALUATIONS
!   ARE SAVED BY ONLY APPROXIMATING THE NONZERO TERMS.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                       WA1,WA2)

!   WHERE

!     FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH CALCULATES
!       THE FUNCTIONS.  FCN MUST BE DECLARED IN AN EXTERNAL STATEMENT IN
!       THE USER CALLING PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!       SUBROUTINE FCN(N,X,FVEC,IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       ----------
!       CALCULATE THE FUNCTIONS AT X AND
!       RETURN THIS VECTOR IN FVEC.
!       ----------
!       RETURN
!       END

!       THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!       THE USER WANTS TO TERMINATE EXECUTION OF FDJAC1.
!       IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF FUNCTIONS AND VARIABLES.

!     X IS AN INPUT ARRAY OF LENGTH N.

!     FVEC IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
!       FUNCTIONS EVALUATED AT X.

!     FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
!       APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.

!     LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!     IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
!       THE EXECUTION OF FDJAC1.  SEE DESCRIPTION OF FCN.

!     ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
!       THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
!       JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
!       ML TO AT LEAST N - 1.

!     EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
!       STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
!       APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
!       FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
!       THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
!       ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE PRECISION.

!     MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
!       THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
!       JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
!       MU TO AT LEAST N - 1.

!     WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.  IF ML + MU + 1 IS AT
!       LEAST N, THEN THE JACOBIAN IS CONSIDERED DENSE, AND WA2 IS
!       NOT REFERENCED.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR

!     FORTRAN-SUPPLIED ... ABS,MAX,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, k, msum
REAL (8)  :: eps, epsmch, h, temp
REAL (8), PARAMETER  :: zero = 0.0_dp

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

eps = SQRT(MAX(epsfcn, epsmch))
msum = ml + mu + 1
IF (msum >= n) THEN
  
!        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
  
  DO  j = 1, n
    temp = x(j)
    h = eps * ABS(temp)
    IF (h == zero) h = eps
    x(j) = temp + h
    CALL fcn(n, x, wa1, iflag)
    IF (iflag < 0) EXIT
    x(j) = temp
    DO  i = 1, n
      fjac(i,j) = (wa1(i)-fvec(i)) / h
    END DO
  END DO
ELSE
  
!        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
  
  DO  k = 1, msum
    DO  j = k, n, msum
      wa2(j) = x(j)
      h = eps * ABS(wa2(j))
      IF (h == zero) h = eps
      x(j) = wa2(j) + h
    END DO
    CALL fcn(n, x, wa1, iflag)
    IF (iflag < 0) EXIT
    DO  j = k, n, msum
      x(j) = wa2(j)
      h = eps * ABS(wa2(j))
      IF (h == zero) h = eps
      DO  i = 1, n
        fjac(i,j) = zero
        IF (i >= j-mu .AND. i <= j+ml) fjac(i,j) = (wa1(i)-fvec(i)) / h
      END DO
    END DO
  END DO
END IF
RETURN

!     LAST CARD OF SUBROUTINE FDJAC1.

END SUBROUTINE fdjac1



SUBROUTINE qform(m, n, q, ldq, wa)

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ldq
REAL (8), INTENT(OUT)  :: q(ldq,m)
REAL (8), INTENT(OUT)  :: wa(m)


!   **********

!   SUBROUTINE QFORM

!   THIS SUBROUTINE PROCEEDS FROM THE COMPUTED QR FACTORIZATION OF AN M BY N
!   MATRIX A TO ACCUMULATE THE M BY M ORTHOGONAL MATRIX Q FROM ITS FACTORED FORM.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE QFORM(M,N,Q,LDQ,WA)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF ROWS OF A AND THE ORDER OF Q.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

!     Q IS AN M BY M ARRAY. ON INPUT THE FULL LOWER TRAPEZOID IN
!       THE FIRST MIN(M,N) COLUMNS OF Q CONTAINS THE FACTORED FORM.
!       ON OUTPUT Q HAS BEEN ACCUMULATED INTO A SQUARE MATRIX.

!     LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.

!     WA IS A WORK ARRAY OF LENGTH M.

!   SUBPROGRAMS CALLED

!     FORTRAN-SUPPLIED ... MIN

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, jm1, k, l, minmn, np1
REAL (8)  :: sum, temp
REAL (8), PARAMETER  :: one = 1.0_dp, zero = 0.0_dp

!     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.

minmn = MIN(m,n)
IF (minmn >= 2) THEN
  DO  j = 2, minmn
    jm1 = j - 1
    DO  i = 1, jm1
      q(i,j) = zero
    END DO
  END DO
END IF

!     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.

np1 = n + 1
IF (m >= np1) THEN
  DO  j = np1, m
    DO  i = 1, m
      q(i,j) = zero
    END DO
    q(j,j) = one
  END DO
END IF

!     ACCUMULATE Q FROM ITS FACTORED FORM.

DO  l = 1, minmn
  k = minmn - l + 1
  DO  i = k, m
    wa(i) = q(i,k)
    q(i,k) = zero
  END DO
  q(k,k) = one
  IF (wa(k) /= zero) THEN
    DO  j = k, m
      sum = zero
      DO  i = k, m
        sum = sum + q(i,j) * wa(i)
      END DO
      temp = sum / wa(k)
      DO  i = k, m
        q(i,j) = q(i,j) - temp * wa(i)
      END DO
    END DO
  END IF
END DO
RETURN

!     LAST CARD OF SUBROUTINE QFORM.

END SUBROUTINE qform


SUBROUTINE qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lda
REAL (8), INTENT(IN OUT)  :: a(lda,n)
LOGICAL, INTENT(IN)        :: pivot
INTEGER, INTENT(IN)        :: lipvt
INTEGER, INTENT(OUT)       :: ipvt(lipvt)
REAL (8), INTENT(OUT)     :: rdiag(n)
REAL (8), INTENT(OUT)     :: acnorm(n)
REAL (8), INTENT(OUT)     :: wa(n)


!   **********

!   SUBROUTINE QRFAC

!   THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN PIVOTING
!   (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE M BY N MATRIX A.
!   THAT IS, QRFAC DETERMINES AN ORTHOGONAL MATRIX Q, A PERMUTATION MATRIX P,
!   AND AN UPPER TRAPEZOIDAL MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING
!   MAGNITUDE, SUCH THAT A*P = Q*R.  THE HOUSEHOLDER TRANSFORMATION FOR
!   COLUMN K, K = 1,2,...,MIN(M,N), IS OF THE FORM

!                         T
!         I - (1/U(K))*U*U

!   WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS.  THE FORM OF THIS
!   TRANSFORMATION AND THE METHOD OF PIVOTING FIRST APPEARED IN THE
!   CORRESPONDING LINPACK SUBROUTINE.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF COLUMNS OF A.

!     A IS AN M BY N ARRAY.  ON INPUT A CONTAINS THE MATRIX FOR WHICH THE
!       QR FACTORIZATION IS TO BE COMPUTED.  ON OUTPUT THE STRICT UPPER
!       TRAPEZOIDAL PART OF A CONTAINS THE STRICT UPPER TRAPEZOIDAL PART OF R,
!       AND THE LOWER TRAPEZOIDAL PART OF A CONTAINS A FACTORED FORM OF Q
!       (THE NON-TRIVIAL ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).

!     LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

!     PIVOT IS A LOGICAL INPUT VARIABLE.  IF PIVOT IS SET TRUE,
!       THEN COLUMN PIVOTING IS ENFORCED.  IF PIVOT IS SET FALSE,
!       THEN NO COLUMN PIVOTING IS DONE.

!     IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT.  IPVT DEFINES THE
!       PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
!       COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
!       IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.

!     LIPVT IS A POSITIVE INTEGER INPUT VARIABLE.  IF PIVOT IS FALSE,
!       THEN LIPVT MAY BE AS SMALL AS 1.  IF PIVOT IS TRUE, THEN
!       LIPVT MUST BE AT LEAST N.

!     RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!       DIAGONAL ELEMENTS OF R.

!     ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE NORMS OF
!       THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
!       IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE WITH RDIAG.

!     WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
!       CAN COINCIDE WITH RDIAG.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR,ENORM

!     FORTRAN-SUPPLIED ... MAX,SQRT,MIN

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, jp1, k, kmax, minmn
REAL (8)  :: ajnorm, epsmch, sum, temp
REAL (8), PARAMETER  :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

!     EPSMCH IS THE MACHINE PRECISION.

epsmch = EPSILON(1.0_dp)

!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.

DO  j = 1, n
  acnorm(j) = enorm(m, a(1:,j))
  rdiag(j) = acnorm(j)
  wa(j) = rdiag(j)
  IF (pivot) ipvt(j) = j
END DO

!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.

minmn = MIN(m,n)
DO  j = 1, minmn
  IF (pivot) THEN
    
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
    
    kmax = j
    DO  k = j, n
      IF (rdiag(k) > rdiag(kmax)) kmax = k
    END DO
    IF (kmax /= j) THEN
      DO  i = 1, m
        temp = a(i,j)
        a(i,j) = a(i,kmax)
        a(i,kmax) = temp
      END DO
      rdiag(kmax) = rdiag(j)
      wa(kmax) = wa(j)
      k = ipvt(j)
      ipvt(j) = ipvt(kmax)
      ipvt(kmax) = k
    END IF
  END IF
  
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
  
  ajnorm = enorm(m-j+1, a(j:,j))
  IF (ajnorm /= zero) THEN
    IF (a(j,j) < zero) ajnorm = -ajnorm
    DO  i = j, m
      a(i,j) = a(i,j) / ajnorm
    END DO
    a(j,j) = a(j,j) + one
    
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS AND UPDATE THE NORMS.
    
    jp1 = j + 1
    IF (n >= jp1) THEN
      DO  k = jp1, n
        sum = zero
        DO  i = j, m
          sum = sum + a(i,j) * a(i,k)
        END DO
        temp = sum / a(j,j)
        DO  i = j, m
          a(i,k) = a(i,k) - temp * a(i,j)
        END DO
        IF (.NOT.(.NOT.pivot.OR.rdiag(k) == zero)) THEN
          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * SQRT(MAX(zero,one-temp**2))
          IF (p05*(rdiag(k)/wa(k))**2 <= epsmch) THEN
            rdiag(k) = enorm(m-j, a(jp1:,k))
            wa(k) = rdiag(k)
          END IF
        END IF
      END DO
    END IF
  END IF
  rdiag(j) = -ajnorm
END DO
RETURN

!     LAST CARD OF SUBROUTINE QRFAC.

END SUBROUTINE qrfac



SUBROUTINE r1mpyq(m, n, a, lda, v, w)

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: lda
REAL (8), INTENT(IN OUT)  :: a(lda,n)
REAL (8), INTENT(IN)      :: v(n)
REAL (8), INTENT(IN)      :: w(n)


!   **********

!   SUBROUTINE R1MPYQ

!   GIVEN AN M BY N MATRIX A, THIS SUBROUTINE COMPUTES A*Q WHERE
!   Q IS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

!         GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

!   AND GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE WHICH
!   ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
!   Q ITSELF IS NOT GIVEN, RATHER THE INFORMATION TO RECOVER THE
!   GV, GW ROTATIONS IS SUPPLIED.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE R1MPYQ(M, N, A, LDA, V, W)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF A.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF COLUMNS OF A.

!     A IS AN M BY N ARRAY.  ON INPUT A MUST CONTAIN THE MATRIX TO BE
!       POSTMULTIPLIED BY THE ORTHOGONAL MATRIX Q DESCRIBED ABOVE.
!       ON OUTPUT A*Q HAS REPLACED A.

!     LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!       WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.

!     V IS AN INPUT ARRAY OF LENGTH N. V(I) MUST CONTAIN THE INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

!     W IS AN INPUT ARRAY OF LENGTH N. W(I) MUST CONTAIN THE INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

!   SUBROUTINES CALLED

!     FORTRAN-SUPPLIED ... ABS, SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i, j, nmj, nm1
REAL (8)  :: COS, SIN, temp
REAL (8), PARAMETER  :: one = 1.0_dp

!     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.

nm1 = n - 1
IF (nm1 >= 1) THEN
  DO  nmj = 1, nm1
    j = n - nmj
    IF (ABS(v(j)) > one) COS = one / v(j)
    IF (ABS(v(j)) > one) SIN = SQRT(one-COS**2)
    IF (ABS(v(j)) <= one) SIN = v(j)
    IF (ABS(v(j)) <= one) COS = SQRT(one-SIN**2)
    DO  i = 1, m
      temp = COS * a(i,j) - SIN * a(i,n)
      a(i,n) = SIN * a(i,j) + COS * a(i,n)
      a(i,j) = temp
    END DO
  END DO
  
!     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
  
  DO  j = 1, nm1
    IF (ABS(w(j)) > one) COS = one / w(j)
    IF (ABS(w(j)) > one) SIN = SQRT(one-COS**2)
    IF (ABS(w(j)) <= one) SIN = w(j)
    IF (ABS(w(j)) <= one) COS = SQRT(one-SIN**2)
    DO  i = 1, m
      temp = COS * a(i,j) + SIN * a(i,n)
      a(i,n) = -SIN * a(i,j) + COS * a(i,n)
      a(i,j) = temp
    END DO
  END DO
END IF
RETURN

!     LAST CARD OF SUBROUTINE R1MPYQ.

END SUBROUTINE r1mpyq



SUBROUTINE r1updt(m, n, s, ls, u, v, w, sing)

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ls
REAL (8), INTENT(IN OUT)  :: s(ls)
REAL (8), INTENT(IN)      :: u(m)
REAL (8), INTENT(IN OUT)  :: v(n)
REAL (8), INTENT(OUT)     :: w(m)
LOGICAL, INTENT(OUT)       :: sing


!   **********

!   SUBROUTINE R1UPDT

!   GIVEN AN M BY N LOWER TRAPEZOIDAL MATRIX S, AN M-VECTOR U,
!   AND AN N-VECTOR V, THE PROBLEM IS TO DETERMINE AN
!   ORTHOGONAL MATRIX Q SUCH THAT

!                 T
!         (S + U*V )*Q

!   IS AGAIN LOWER TRAPEZOIDAL.

!   THIS SUBROUTINE DETERMINES Q AS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS

!         GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)

!   WHERE GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE
!   WHICH ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
!   Q ITSELF IS NOT ACCUMULATED, RATHER THE INFORMATION TO RECOVER THE GV,
!   GW ROTATIONS IS RETURNED.

!   THE SUBROUTINE STATEMENT IS

!     SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)

!   WHERE

!     M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF ROWS OF S.

!     N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!       OF COLUMNS OF S.  N MUST NOT EXCEED M.

!     S IS AN ARRAY OF LENGTH LS. ON INPUT S MUST CONTAIN THE LOWER
!       TRAPEZOIDAL MATRIX S STORED BY COLUMNS. ON OUTPUT S CONTAINS
!       THE LOWER TRAPEZOIDAL MATRIX PRODUCED AS DESCRIBED ABOVE.

!     LS IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!       (N*(2*M-N+1))/2.

!     U IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE VECTOR U.

!     V IS AN ARRAY OF LENGTH N. ON INPUT V MUST CONTAIN THE VECTOR V.
!       ON OUTPUT V(I) CONTAINS THE INFORMATION NECESSARY TO
!       RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.

!     W IS AN OUTPUT ARRAY OF LENGTH M. W(I) CONTAINS INFORMATION
!       NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED ABOVE.

!     SING IS A LOGICAL OUTPUT VARIABLE.  SING IS SET TRUE IF ANY OF THE
!       DIAGONAL ELEMENTS OF THE OUTPUT S ARE ZERO.  OTHERWISE SING IS
!       SET FALSE.

!   SUBPROGRAMS CALLED

!     MINPACK-SUPPLIED ... SPMPAR

!     FORTRAN-SUPPLIED ... ABS,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE, JOHN L. NAZARETH

!   **********
INTEGER    :: i, j, jj, l, nmj, nm1
REAL (8)  :: COS, cotan, giant, SIN, TAN, tau, temp
REAL (8), PARAMETER  :: one = 1.0_dp, p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

!     GIANT IS THE LARGEST MAGNITUDE.

giant = HUGE(1.0_dp)

!     INITIALIZE THE DIAGONAL ELEMENT POINTER.

jj = (n*(2*m-n+1)) / 2 - (m-n)

!     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.

l = jj
DO  i = n, m
  w(i) = s(l)
  l = l + 1
END DO

!     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
!     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.

nm1 = n - 1
IF (nm1 >= 1) THEN
  DO  nmj = 1, nm1
    j = n - nmj
    jj = jj - (m-j+1)
    w(j) = zero
    IF (v(j) /= zero) THEN
      
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE J-TH ELEMENT OF V.
      
      IF (ABS(v(n)) < ABS(v(j))) THEN
        cotan = v(n) / v(j)
        SIN = p5 / SQRT(p25+p25*cotan**2)
        COS = SIN * cotan
        tau = one
        IF (ABS(COS)*giant > one) tau = one / COS
      ELSE
        TAN = v(j) / v(n)
        COS = p5 / SQRT(p25+p25*TAN**2)
        SIN = COS * TAN
        tau = SIN
      END IF
      
!        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
!        NECESSARY TO RECOVER THE GIVENS ROTATION.
      
      v(n) = SIN * v(j) + COS * v(n)
      v(j) = tau
      
!        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
      
      l = jj
      DO  i = j, m
        temp = COS * s(l) - SIN * w(i)
        w(i) = SIN * s(l) + COS * w(i)
        s(l) = temp
        l = l + 1
      END DO
    END IF
  END DO
END IF

!     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.

DO  i = 1, m
  w(i) = w(i) + v(n) * u(i)
END DO

!     ELIMINATE THE SPIKE.

sing = .false.
IF (nm1 >= 1) THEN
  DO  j = 1, nm1
    IF (w(j) /= zero) THEN
      
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!        J-TH ELEMENT OF THE SPIKE.
      
      IF (ABS(s(jj)) < ABS(w(j))) THEN
        cotan = s(jj) / w(j)
        SIN = p5 / SQRT(p25 + p25*cotan**2)
        COS = SIN * cotan
        tau = one
        IF (ABS(COS)*giant > one) tau = one / COS
      ELSE
        TAN = w(j) / s(jj)
        COS = p5 / SQRT(p25+p25*TAN**2)
        SIN = COS * TAN
        tau = SIN
      END IF
      
!        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
      
      l = jj
      DO  i = j, m
        temp = COS * s(l) + SIN * w(i)
        w(i) = -SIN * s(l) + COS * w(i)
        s(l) = temp
        l = l + 1
      END DO
      
!        STORE THE INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION.
      
      w(j) = tau
    END IF
    
!        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
    
    IF (s(jj) == zero) sing = .true.
    jj = jj + (m-j+1)
  END DO
END IF

!     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.

l = jj
DO  i = n, m
  s(l) = w(i)
  l = l + 1
END DO
IF (s(jj) == zero) sing = .true.
RETURN

!     LAST CARD OF SUBROUTINE R1UPDT.

END SUBROUTINE r1updt


FUNCTION enorm(n, x) RESULT(fn_val)
 
INTEGER, INTENT(IN)    :: n
REAL (8), INTENT(IN)  :: x(n)
REAL (8)              :: fn_val

!   **********

!   FUNCTION ENORM

!   GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE EUCLIDEAN NORM OF X.

!   THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF SQUARES IN THREE
!   DIFFERENT SUMS.  THE SUMS OF SQUARES FOR THE SMALL AND LARGE COMPONENTS
!   ARE SCALED SO THAT NO OVERFLOWS OCCUR.  NON-DESTRUCTIVE UNDERFLOWS ARE
!   PERMITTED.  UNDERFLOWS AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
!   SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
!   THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS DEPEND ON
!   TWO CONSTANTS, RDWARF AND RGIANT.  THE MAIN RESTRICTIONS ON THESE CONSTANTS
!   ARE THAT RDWARF**2 NOT UNDERFLOW AND RGIANT**2 NOT OVERFLOW.
!   THE CONSTANTS GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.

!   THE FUNCTION STATEMENT IS

!     REAL FUNCTION ENORM(N, X)

!   WHERE

!     N IS A POSITIVE INTEGER INPUT VARIABLE.

!     X IS AN INPUT ARRAY OF LENGTH N.

!   SUBPROGRAMS CALLED

!     FORTRAN-SUPPLIED ... ABS,SQRT

!   ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!   BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!   **********
INTEGER    :: i
REAL (8)  :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
REAL (8), PARAMETER  :: rdwarf = 1.0D-100, rgiant = 1.0D+100

s1 = 0.0_dp
s2 = 0.0_dp
s3 = 0.0_dp
x1max = 0.0_dp
x3max = 0.0_dp
floatn = n
agiant = rgiant / floatn
DO  i = 1, n
  xabs = ABS(x(i))
  IF (xabs <= rdwarf .OR. xabs >= agiant) THEN
    IF (xabs > rdwarf) THEN
      
!              SUM FOR LARGE COMPONENTS.
      
      IF (xabs > x1max) THEN
        s1 = 1.0_dp + s1 * (x1max/xabs) ** 2
        x1max = xabs
      ELSE
        s1 = s1 + (xabs/x1max) ** 2
      END IF
    ELSE
      
!              SUM FOR SMALL COMPONENTS.
      
      IF (xabs > x3max) THEN
        s3 = 1.0_dp + s3 * (x3max/xabs) ** 2
        x3max = xabs
      ELSE
        IF (xabs /= 0.0_dp) s3 = s3 + (xabs/x3max) ** 2
      END IF
    END IF
  ELSE
    
!           SUM FOR INTERMEDIATE COMPONENTS.
    
    s2 = s2 + xabs ** 2
  END IF
END DO

!     CALCULATION OF NORM.

IF (s1 /= 0.0_dp) THEN
  fn_val = x1max * SQRT(s1 + (s2/x1max)/x1max)
ELSE
  IF (s2 /= 0.0_dp) THEN
    IF (s2 >= x3max) fn_val = SQRT(s2*(1.0_dp + (x3max/s2)*(x3max*s3)))
    IF (s2 < x3max) fn_val = SQRT(x3max*((s2/x3max) + (x3max*s3)))
  ELSE
    fn_val = x3max * SQRT(s3)
  END IF
END IF
RETURN
END FUNCTION enorm

END MODULE
