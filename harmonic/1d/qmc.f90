program qmc

    implicit none

    integer, parameter :: N = 400     ! number of walkers
    integer, parameter :: k = 3*10**4 ! number of walks
    integer, parameter :: m = 5       ! number of alpha values
    real(8), parameter :: lower = 0.4 ! lower limit of alpha range
    real(8), parameter :: upper = 0.6 ! upper limit of alpha range
    real(8), dimension(1, m) :: aarray
    real(8) :: a
    real(8), dimension(1, N) :: x
    real(8), dimension(1, k) :: E
    real(8) :: E_av
    real(8) :: Esq
    real(8) :: var
    integer :: i
    integer :: j

    call alpha_array(m, lower, upper, aarray)
    open (unit = 12, file = "energies.dat", status = "replace")

    do i = 1, m

        a = aarray(1, i)
        call init_random_seed()
        call random_number(x)

        do j = 1, 1000
            call calc_x(a, N, x)
        end do

        do j = 1, k
            call calc_x(a, N, x)
            call calc_E(a, N, x, E(1, j))
        end do

        E_av = 1d0 / k * sum(E)
        Esq = 1d0 / k * sum(E**2)
        var = Esq - E_av**2
        write (12, *) a, E_av, var

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

    subroutine calc_x(a, N, x)

        real(8), intent(in) :: a
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(inout) :: x
        real(8), dimension(1, N) :: delta_x
        real(8), dimension(1, N) :: p
        real(8), dimension(1, N) :: u
        real(8), dimension(1, N) :: tf 

        call random_number(delta_x)
        delta_x = 1d0 * (delta_x(:, :) - 0.5d0)
        p = min(psi(a, x + delta_x) / psi(a, x) , 1d0)
        call random_number(u)
        tf = merge(1, 0, u < p)
        x(1, :) = x(1, :) + tf(1, :) * delta_x(1, :)

    end subroutine

    subroutine calc_E(a, N, x, E)

        real(8), intent(in) :: a
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(in) :: x
        real(8), dimension(1, N) :: E_L
        real(8), intent(inout) :: E

        E_L = a + (1 / 2d0 - 2d0 * a**2) * x**2
        E = 1d0 / N * sum(E_L)

    end subroutine

    function psi(a, x)

        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: x
        real(8), dimension(1, size(x)) :: psi

        psi = exp(-2d0 * a * x**2)

    end function

end program