program qmc

    use H2

    implicit none

    integer, parameter :: num_walkers = 400
    integer, parameter :: num_walks = 30000
    real(8), parameter :: s_lower = 1
    real(8), parameter :: s_upper = 2
    integer, parameter :: s_num = 9
    real(8), parameter :: b_lower = 0.05
    real(8), parameter :: b_upper = 0.25
    integer, parameter :: b_num = 9
    real(8) :: s, b
    integer :: i, j

    open (unit = 12, file = 'energies.dat', status = 'replace')

    call init(num_walkers, 0d0, 0d0)

    do i = 1, 1 !s_num
        s = 0
        !s = next_value(i, s_lower, s_upper, s_num)
        print *, s
        do j = 1, b_num
            b = next_value(j, b_lower, b_upper, b_num)
            print *, "------", b
            call monte_carlo(num_walkers, num_walks, s, b)
        end do
    end do

    close (unit = 12)

contains

    real(8) function next_value(i, lower, upper, num) result(val)

        real(8), intent(in) :: lower, upper
        integer, intent(in) :: i, num
        
        val = (i - 1) * (upper - lower) / (num - 1) + lower

    end function

end program
