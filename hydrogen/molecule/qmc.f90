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

        integer, parameter :: cutoff = 4000
        integer, intent(in) :: num_walkers, num_walks
        real(8) :: dr(6, num_walkers), E_L(num_walks - cutoff, num_walkers), E_T, E_Tsq
        integer :: i

        do i = 1, num_walks
            call random_number(dr)
            dr = 2 * dr - 1
            call step(dr)
            if (i > cutoff) then
                call calc_local_energies(num_walkers, E_L(i - cutoff, :))
            end if
        end do
        E_T = sum(E_L) / (num_walkers * (num_walks - cutoff))
        E_Tsq = sum(E_L**2) / (num_walkers * (num_walks - cutoff))

        call write_to_file(E_T, E_Tsq - E_T**2)

    end subroutine

end program
