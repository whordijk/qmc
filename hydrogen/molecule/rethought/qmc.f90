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
        real(8) :: E
        real(8) :: Esq
        real(8) :: var
        real(8) :: s, b
        integer :: i, j
        real(8) :: dr(6, num_walkers)

        E = 0d0
        Esq = 0d0

        s = 0d0
        b = 0.05d0

        call init(num_walkers, s)

        do i = 1, num_walks
            call random_number(dr)
            dr = 2 * dr - 1
            call step(dr)
            if (i > 4000) then
                do j = 1, num_walkers
                    E = E + calc_energy(walkers(: j))
                    Esq = Esq + calc_energy(walkers(:, j))**2
                end do
            end if
        end do

        E = E / (num_walkers * (num_walks - 4000))
        Esq = Esq / (num_walkers * (num_walks - 4000))
        var = Esq - E**2

        write (12, *) s, b, E, var

    end subroutine

end program
