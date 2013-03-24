program qmc

    use H2

    implicit none

    integer, parameter :: num_walkers = 1
    integer, parameter :: num_walks = 30000

    open (unit = 12, file = 'energies.dat', status = 'replace')
    
    call monte_carlo(num_walkers, num_walks)

    close (unit = 12)

contains

    subroutine monte_carlo(num_walkers, num_walks)

        integer, intent(in) :: num_walkers
        integer, intent(in) :: num_walks
        real(8) :: E_L(num_walks, num_walkers)
        real(8) :: E_T
        real(8) :: E_T_sq
        real(8) :: var
        real(8) :: s, b
        integer :: i, j
        real(8) :: dr(6, num_walkers)

        s = 0d0
        b = 0.05d0

        call init(num_walkers, s)

        do i = 1, num_walks
            call random_number(dr)
            dr = 2 * dr - 1
            call step(dr)
            do j = 1, num_walkers
                E_L(i, j) = calc_energy(dr(:, j)) ! I need to pass the argument 'walker' but I don't want it to be in this file
            end do
        end do

        E_T = sum(E_L) / (num_walkers * num_walks)
        E_T_sq = sum(E_L**2) / (num_walkers * num_walks)
        var = E_T_sq - E_T**2

        write (12, *) s, b, E_T, var

    end subroutine

end program
