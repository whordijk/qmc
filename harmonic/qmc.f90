program qmc

    implicit none

    integer, parameter :: N = 400     ! number of walkers
    integer, parameter :: k = 26*10**3! number of walks
    integer, parameter :: m = 5       ! number of alpha values
    real(8), parameter :: lower = 0.4 ! lower limit of alpha range
    real(8), parameter :: upper = 0.6 ! upper limit of alpha range
    real(8), dimension(1, m) :: aarray
    real(8) :: a
    real(8), dimension(3, N) :: r
    real(8), dimension(1, N) :: E_L
    real(8), dimension(k, N) :: E_L_tot
    real(8) :: Eav
    real(8) :: Esq
    real(8) :: var
    integer :: i
    integer :: j

    call alpha_array(m, lower, upper, aarray)
    open (unit = 12, file = "energies.dat", status = "replace")

    do i = 1, m

        a = aarray(1, i)
        call init_random_seed()
        call random_number(r)

        do j = 1, 4000
            call calc_r(a, N, r)
        end do

        do j = 1, k
            call calc_r(a, N, r)
            call calc_E_L(a, r, E_L)
            E_L_tot(j, :) = E_l(1, :)
        end do

        Eav = 1d0 / N * sum(1d0 / k * sum(E_L_tot, dim = 1))
        Esq = 1d0 / N * sum(1d0 / k * sum(E_L_tot**2, dim = 1))
        var = Esq - Eav**2
        write (12, *) a, Eav, var

    end do

    close (unit = 12)

contains

    subroutine alpha_array(m, lower, upper, aarray)

        integer, intent(in) :: m
        real(8), intent(in) :: lower
        real(8), intent(in) :: upper
        real(8), dimension (:, :), intent(inout) :: aarray
        integer :: i

        do i = 1, m
            aarray(1, i) = (i - 1) * (upper - lower) / (m - 1) + lower
        end do

    end subroutine

    subroutine init_random_seed()

        integer :: i
        integer :: n
        integer :: clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        deallocate(seed)

    end subroutine

    subroutine calc_r(a, N, r)

        real(8), intent(in) :: a
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(inout) :: r
        real(8), dimension(3, N) :: delta_r
        real(8), dimension(1, N) :: p
        real(8), dimension(1, N) :: x
        real(8), dimension(1, N) :: tf 

        call random_number(delta_r)
        delta_r = 1d0 * (delta_r(:, :) - 0.5d0)
        p = min(psi(a, r + delta_r) / psi(a, r) , 1d0)
        call random_number(x)
        tf = merge(1, 0, x < p)
        r(1, :) = r(1, :) + tf(1, :) * delta_r(1, :)
        r(2, :) = r(2, :) + tf(1, :) * delta_r(2, :)
        r(3, :) = r(3, :) + tf(1, :) * delta_r(3, :)

    end subroutine

    subroutine calc_E_L(a, r, E_L)

        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: r
        real(8), dimension(:, :), intent(inout) :: E_L

        E_L(1, :) = 3d0 * a + (1 / 2d0 - 2d0 * a**2) * sum(r**2, dim = 1)

    end subroutine

    function psi(a, r)

        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: r
        real(8), dimension(1, size(r, dim = 2)) :: psi

        psi(1, :) = exp(-2d0 * a * sum(r**2, dim = 1))

    end function

end program
