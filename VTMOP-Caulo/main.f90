MODULE MAIN

! Import the library of VTMOP interfaces.
USE VTMOP_LIB
USE MODEL

IMPLICIT NONE

CONTAINS 

SUBROUTINE START()


! Output file.
CHARACTER(LEN=20), PARAMETER :: OUTFILE = "sOutput.txt"

! Adjust the problem dimensions below. D is the number of design variables,
! and P is the number of objectives. Any values P >= 2 and D >= P are allowed.
INTEGER, PARAMETER :: D=47, P=2

! Set the following to FALSE if changing the objective function is changed.
LOGICAL, PARAMETER :: CHECK_MAE = .TRUE.

! Local variables used for defining the problem.
INTEGER :: I ! Loop index.
INTEGER :: IERR ! Error flag.
REAL(KIND=R8) :: MU ! Working precision.
REAL(KIND=R8) :: LB(D), UB(D) ! Lower/upper bound constraints.
REAL(KIND=R8) :: OBJ_BOUNDS(P,2) ! Lower/upper objective bounds.
REAL(KIND=R8), ALLOCATABLE :: DB_X(:,:), DB_F(:,:) ! Full database
REAL(KIND=R8), ALLOCATABLE :: PARETO_X(:,:), PARETO_F(:,:) ! Solution set.

! Get the unit roundoff.
MU = EPSILON(0.0_R8)

! Set the bounds constraints.
LB(1:D) =  (/ 0.7000,    0.7000,    0.7000,    0.0500,    0.0333,    0.1280,   &
               0.0400,    0.1210,    0.0300,    2.8000,    0.3000,    0.2500,   &
               0.0200,    0.4950,    0.0450,    0.0415,    0.0325,    0.0140,   &
               0.0425,    0.0590,    0.0300,    0.0216,    0.0300,    0.3000,   &
               1.5000,    0.3500,    0.7500,    0.5000,    0.5000,    0.5500,   &
               0.5000,    0.0750,    0.1000,   70.0000,    1.0000,    0.0050,   &
               0.5000,    0.0500,    0.0750,    0.0200,    0.0200,    0.0050,   &
               0.2500,    2.5000,    0.5000,    0.0500,    0.5000 /)/2
UB(1:D) =  (/ 2.8000,    2.8000,    2.8000,    0.2000,    0.1334,    0.5120,  & 
               0.1600,    0.4840,    0.1200,   11.2000,    1.2000,    1.0000,  &  
               0.0800,    1.9800,    0.1800,    0.1660,    0.1300,    0.0560,  &
               0.1700,    0.2360,    0.1200,    0.0864,    0.1200,    1.2000,  &  
               6.0000,    1.4000,    3.0000,    2.0000,    2.0000,    2.2000,  &  
               2.0000,    0.3000,    0.4000,  280.0000,    4.0000,    0.0200,  &  
               2.0000,    0.2000,    0.3000,    0.0800,    0.0800,    0.0200,  &  
               1.0000,   10.0000,    2.0000,    0.2000,    2.0000 /)*2


! Don't use any objective bounds.
OBJ_BOUNDS(:,1) = -HUGE(0.0_R8)
OBJ_BOUNDS(:,2) = HUGE(0.0_R8)

! Set good starting point
allocate (DB_X(D,1),DB_F(P,1))

DB_X(:,1)=(/   1.4000,    1.4000,    1.4000,    0.1000,    0.0667,    0.2560,   &  
               0.0800,    0.2420,    0.0600,    5.6000,    0.6000,    0.5000,   & 
               0.0400,    0.9900,    0.0900,    0.0830,    0.0650,    0.0280,   & 
               0.0850,    0.1180,    0.0600,    0.0432,    0.0600,    0.6000,   & 
               3.0000,    0.7000,    1.5000,    1.0000,    1.0000,    1.1000,   & 
               1.0000,    0.1500,    0.2000,  140.0000,    2.0000,    0.0100,   & 
               1.0000,    0.1000,    0.1500,    0.0400,    0.0400,    0.0100,   & 
               0.5000,    5.0000,    1.0000,    0.1000,    1.0000  /)

DB_F(:,1) = (/ 2.3356,  0.1016 /)


! Call VTMOP driver to solve the problem.

! Do 1000 function evaluations with checkpointing.
! CALL VTMOP_SOLVE( D, P, LB, UB, CAULO_OBJ, PARETO_X, PARETO_F, IERR,   &
!                   ADAPTIVE_SEARCH=.FALSE., BB_BUDGET=10000,          &
!                   MAXITERS=10000, DES_PTS=DB_X, OBJ_PTS=DB_F )

CALL VTMOP_SOLVE( D, P, LB, UB, CAULO_OBJ, PARETO_X, PARETO_F, IERR,   &
                  ADAPTIVE_SEARCH=.FALSE., BB_BUDGET=5000,          &
                  MAXITERS=5000, SEARCH_BUDGET=D,                  &
                  INITIAL_SBUDGET=2*D, LOPT_BUDGET=2500,           &
                  DECAY=0.5_R8, DES_TOL=SQRT(MU), EPS=SQRT(MU),    &
                  EPSW=MU**0.25_R8, OBJ_TOL=SQRT(MU),              &
                  MIN_RADF=0.02_R8, TRUST_RADF=0.2_R8,             &
                  DES_PTS=DB_X, OBJ_PTS=DB_F, LOCAL_OPT=GPS,       &
                  OBJ_BOUNDS=OBJ_BOUNDS, FIT_SURROGATES=LSHEP_FIT, &
                  EVAL_SURROGATES=LSHEP_EVAL, PFLAG=0, ICHKPT=1    )
IF (IERR .GE. 10) THEN
   WRITE(*,11) "An error occurred. IERR=", IERR; STOP; END IF


! Check the solution.
CALL CHECK_SOLUTION()

! Write the outputs to OUTFILE.
OPEN(99, FILE=OUTFILE, STATUS="REPLACE")
! WRITE(99,12) "Error code IERR =", IERR;
! Print the summary statistics.
WRITE(99,*) "Summary:"
WRITE(99,12) "Number of nondominated/efficient points: ", SIZE(PARETO_X,2)
WRITE(99,12) "Total number of function evaluations: ", SIZE(DB_X,2)

! Print the nondominated objective set.
WRITE(99,*)
WRITE(99,*) "Nondominated point set:"
DO I = 1, SIZE(PARETO_F, 2)
   WRITE(99,10) PARETO_F(:,I)
END DO

! Print the efficient set.
WRITE(99,*)
WRITE(99,*) "Efficient point set:"
DO I = 1, SIZE(PARETO_X, 2)
   WRITE(99,10) PARETO_X(:,I)
END DO

! Print the full database of observed objective points.
WRITE(99,*)
WRITE(99,*) "Full database of objective values observed:"
DO I = 1, SIZE(DB_F,2)
   WRITE(99,10) DB_F(:,I)
END DO

! Print the full database of design points evaluated.
WRITE(99,*)
WRITE(99,*) "Full database of evaluated design points:"
DO I = 1, SIZE(DB_X,2)
   WRITE(99,10) DB_X(:,I)
END DO

! Close the output file.
CLOSE(99)

! Free the output arrays.
DEALLOCATE(DB_F, DB_X, PARETO_F, PARETO_X)

! Output formats.
10 FORMAT(1X,5ES15.7)
11 FORMAT(1X,A,I3)
12 FORMAT(1X,A,I5)
13 FORMAT(2X,A)

CONTAINS

SUBROUTINE CHECK_SOLUTION
! Auxiliary subroutine for checking residuals against a tolerance,
! in the special case of DTLZ2.

! Local variables.
INTEGER :: N
REAL(KIND=R8) :: MAE, TOLERANCE
! N is the number of points on the Pareto front approximation.
N = SIZE(PARETO_F,2)
! The tolerance is a small number proportional to (P * D).
TOLERANCE = MAX(0.01_R8, SQRT(MU)) * REAL(P*D, KIND=R8)
! Compute the mean absolute error (MAE).
MAE = 0.0_R8
DO I = 1, N
   MAE = MAE + ABS(NORM2(PARETO_F(:,I)) - 1.0_R8)
END DO
MAE = MAE / REAL(N, KIND=R8)
! Perform sanity checks for correctness.
IF (N < P + 1) THEN
   WRITE(*,*) "The number of solution points is unreasonably low."
   WRITE(*,*) "An installation error may have occurred."
ELSE IF (CHECK_MAE .AND. MAE > TOLERANCE) THEN
   WRITE(*,*) "The MAE did not meet the tolerance."
   WRITE(*,*) "An installation error may have occurred."
ELSE
   WRITE(*,*) "The serial installation appears to be correct."
   WRITE(*,*) "The full output is contained in the file ", OUTFILE
END IF
RETURN
END SUBROUTINE CHECK_SOLUTION



END SUBROUTINE START


END MODULE MAIN
