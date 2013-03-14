program qmc

    implicit none

    integer, parameter :: N = 400           ! number of walkers
    integer, parameter :: k = 26*10**3      ! number of walks
    integer, parameter :: beta_m = 5        ! number of beta values
	real(8), parameter :: beta_lower = 0.4  ! lower limit of beta range
    real(8), parameter :: beta_upper = 0.7  ! upper limit of beta range
    integer, parameter :: beta_m = 5        ! number of beta values
    real(8), dimension(1, m) :: beta_array
    real(8) :: beta
	real(8), parameter :: s_lower = 1       ! lower limit of s range
    real(8), parameter :: s_upper = 2       ! upper limit of s range
    integer, parameter :: s_m = 5           ! number of s values
    real(8), dimension(1, m) :: s_array
    real(8) :: s
    real(8), dimension(3, N) :: r1
    real(8), dimension(3, N) :: r2
    real(8) :: Eav
    real(8) :: var
    integer :: i

    beta_array = array(beta_m, beta_lower, beta_upper)
    s_array = array(s_m, s_lower, s_upper)
    
    open (unit = 12, file = "energies.dat", status = "replace")
    
    do i = 1, m
        a = aarray(1, i)
        
        call initialize(s, beta, N, r1, r2)
        call montecarlo(s, beta, N, k, r1, r2, Eav, var)
        
        write (12, *) s, beta, Eav, var
    end do
    
    close (unit = 12)

contains

    function array(m, lower, upper)

        integer, intent(in) :: m
        real(8), intent(in) :: lower
        real(8), intent(in) :: upper
        real(8), dimension (:, :), intent(out) :: array
        integer :: i

        do i = 1, m
            array(1, i) = (i - 1) * (upper - lower) / (m - 1) + lower
        end do

    end function

    subroutine initialize(s, beta, N, r1, r2)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(inout) :: r1
        real(8), dimension(:, :), intent(inout) :: r2
        integer :: i

        call init_random_seed()
        call random_number(r1)
        call random_number(r2)

        do i = 1, 4000
            call metropolis(a, N, x)
        end do

    end subroutine

    subroutine montecarlo(a, N, k, x, Eav, var)

        real(8), intent(in) :: a
        integer, intent(in) :: N
        integer, intent(in) :: k
        real(8), dimension(:, :), intent(inout) :: x
        real(8), intent(inout) :: Eav
        real(8), intent(inout) :: var
        integer :: i
        real(8), dimension(1, N) :: E_L
        real(8), dimension(k, N) :: E_L_array
        real(8) :: Esq

        do i = 1, k
            call metropolis(a, N, x)
            call calc_E_L(a, x, E_L)
            E_L_array(i, :) = E_L(1, :)
        end do

        Eav = 1d0 / N * sum(1d0 / k * sum(E_L_array, dim = 1))
        Esq = 1d0 / N * sum(1d0 / k * sum(E_L_array**2, dim = 1))
        var = Esq - Eav**2

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

    subroutine metropolis(a, N, x)

        real(8), intent(in) :: a
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(inout) :: x
        real(8), dimension(3, N) :: delta_x
        real(8), dimension(1, N) :: p
        real(8), dimension(1, N) :: u
        real(8), dimension(1, N) :: tf 

        call random_number(delta_x)
        delta_x = 1d0 * (delta_x(:, :) - 0.5d0)
        p = min(psi(a, x + delta_x) / psi(a, x) , 1d0)
        call random_number(u)
        tf = merge(1, 0, u < p)
        x(1, :) = x(1, :) + tf(1, :) * delta_x(1, :)
        x(2, :) = x(2, :) + tf(1, :) * delta_x(2, :)
        x(3, :) = x(3, :) + tf(1, :) * delta_x(3, :)

    end subroutine

    subroutine calc_E_L(a, x, E_L)

        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: x
        real(8), dimension(:, :), intent(inout) :: E_L

        E_L(1, :) = 3d0 * a + (1 / 2d0 - 2d0 * a**2) * sum(x**2, dim = 1)

    end subroutine

    function psi(a, x)

        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: x
        real(8), dimension(1, size(x, dim = 2)) :: psi

        psi(1, :) = exp(-2d0 * a * sum(x**2, dim = 1))

    end function

end program
