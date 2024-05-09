
MODULE grid_equal

IMPLICIT NONE

CONTAINS 

subroutine linspace(from, to, array)
    REAL(8), INTENT(IN) :: from, to
    REAL(8), DIMENSION(:), INTENT(OUT) :: array
    REAL(8) :: range
    integer :: nn, i
    nn = size(array)
    range = to - from

    if (nn == 0) return

    if (nn == 1) then
        array(1) = from
        return
    end if


    do i=1, nn
        array(i) = from + range * (i - 1) / (nn - 1)
    end do
end subroutine linspace 




    end module 
 