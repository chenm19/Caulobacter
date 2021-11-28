#include "fintrf.h"

MODULE MODEL


!USE DUNI_OMP
! USE OMP_LIB

IMPLICIT NONE

! PRIVATE

PUBLIC:: CAULO_OBJ


CONTAINS

SUBROUTINE CAULO_OBJ(PARAMS, Error, IERR)
  USE REAL_PRECISION, ONLY: R8
! This is an objective function that fits experimental statistics from
! Table 1 of [1] from the budding yeast cell cycle model.

    REAL(KIND=R8), INTENT(IN):: PARAMS(:) 
! but is assumed shape for optimization algorithms.
    REAL(KIND=R8), INTENT(OUT):: Error(:) ! Fitting error from experimental statistics.
    INTEGER, INTENT(OUT):: IERR

    integer:: mexCallMATLAB
    mwPointer:: mxGetPr, mxCreateDoubleMatrix 
    mwPointer:: input(1), output(1)
    mwSize, parameter :: outD = 2 ! output dimension
    mwSize, parameter :: inD = 47 ! input dimension


    input(1) = mxCreateDoubleMatrix(1, inD+1, 0) ! must be one size larger
    call mxCopyReal8ToPtr(PARAMS, mxGetPr(input(1)), inD) 

    IERR = mexCallMATLAB(1, output, 1, input, 'fitness')

    call mxCopyPtrToReal8(mxGetPr(output(1)), Error, outD)
    ! call mxCopyPtrToReal8(mxGetPr(output(2)), IFLAG, m)

    call mxDestroyArray(input(1))


END SUBROUTINE CAULO_OBJ




END MODULE MODEL
