#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

    use MAIN 

    ! number of input arguments, number of output arguments
    integer :: nlhs, nrhs  
    ! pointer to inputs and outputs
    mwPointer :: plhs(*), prhs(*) 

    call START()

end subroutine
