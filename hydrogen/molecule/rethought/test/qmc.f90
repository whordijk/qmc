program qmc

    use H2

    implicit none

    integer, parameter :: num_walkers = 400
    integer, parameter :: num_walks = 30000
    real(8), parameter :: b_min = 0.05
    real(8), parameter :: b_max = 0.25
    integer, parameter :: b_num = 9
    integer :: i

    open (unit = 12, file = 'energies.dat', status = 'replace')

    call init(num_walkers)

    do i = 1, b_num
        call set_params(i, b_min, b_max, b_num)
        call monte_carlo(num_walkers, num_walks)
    end do

    close (unit = 12)

contains

    subroutine monte_carlo(num_walkers, num_walks)

        integer, intent(in) :: num_walkers
        integer, intent(in) :: num_walks
        real(8) :: dr(6, num_walkers)

        call random_number(dr)
        dr = (2 * dr - 1) / 2
        call step(dr)
        call calc_quantities(num_walkers, num_walks)

    end subroutine

end program
