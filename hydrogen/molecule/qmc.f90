program qmc

    use H2

    implicit none

    integer, parameter :: num_walkers = 400
    integer, parameter :: num_init   = 4000
    integer, parameter :: num_walks = 26000
    real(8), parameter :: s_min = 0.5
    real(8), parameter :: s_max = 4
    integer, parameter :: s_num = 21
    real(8), parameter :: b_min = 0.2
    real(8), parameter :: b_max = 0.8
    integer, parameter :: b_num = 21
    integer :: i, j

    open (unit = 12, file = 'energies.dat', status = 'replace')

    call init(num_walkers)
    do i = 1, s_num
        do j = 1, b_num
            call set_params(i, s_min, s_max, s_num, j, b_min, b_max, b_num)
            call monte_carlo(num_init, 0)
            call monte_carlo(num_walks, 1)
        end do
    end do

    close (unit = 12)

contains

    subroutine monte_carlo(steps, n)

        integer, intent(in) :: steps, n
        real(8) :: dr, mu(2), E_T, E_Tsq, var
        integer :: i

        E_T = 0
        E_Tsq = 0
        dr = 1
        do i = 1, steps
            call step(dr)
            call local_energy(mu)
            E_T = E_T + mu(1)
            E_Tsq = E_Tsq + mu(2)
        end do

        E_T = E_T / steps
        var = E_Tsq / steps - E_T**2
        if (n == 1) then
            call write_to_file(E_T, var)
        end if

    end subroutine

end program
