program qmc

    implicit none

    integer, parameter :: N = 400           ! number of walkers
    real(8), parameter :: s_lower = 1       ! lower limit of s range
    real(8), parameter :: s_upper = 2       ! upper limit of s range
    integer, parameter :: s_m = 5           ! number of s values
    real(8), dimension(1, s_m) :: s_array
    real(8) :: s
    integer, parameter :: k = 26*10**3      ! number of walks
    real(8), parameter :: beta_lower = 0.4  ! lower limit of beta range
    real(8), parameter :: beta_upper = 0.7  ! upper limit of beta range
    integer, parameter :: beta_m = 5        ! number of beta values
    real(8), dimension(1, beta_m) :: beta_array
    real(8) :: beta
    real(8) :: a
    real(8), dimension(6, N) :: r
    real(8) :: Eav
    real(8) :: var
    integer :: i
    integer :: j

    beta_array = array(beta_m, beta_lower, beta_upper)
    s_array = array(s_m, s_lower, s_upper)
    
    open (unit = 12, file = "energies.dat", status = "replace")
    
    do i = 1, s_m
        s = s_array(1, i)
        call cusp(s)
        do j = 1, beta_m
            beta = beta_array(1, j)

            call initialize(s, beta, a, N, r)
            call montecarlo(s, beta, a, N, k, r, Eav, var)
        
            write (12, *) s, beta, Eav, var
        end do
    end do
    
    close (unit = 12)

contains

    function array(m, lower, upper)

        integer, intent(in) :: m
        real(8), intent(in) :: lower
        real(8), intent(in) :: upper
        real(8), dimension (1, m) :: array
        integer :: i

        do i = 1, m
            array(1, i) = (i - 1) * (upper - lower) / (m - 1) + lower
        end do
    
    end function

    subroutine cusp(s)

        real(8), intent(in) :: s
        real(8) :: a
        integer :: i

        print *, "________________________________________"
        a = 1
        do i = 1, 5
            a = a - (1d0 + exp(- s / a) - a) / (s / a**2 * exp(- s / a) - 1)
            print *, a
        end do
        print*, "________________________________________"

    end subroutine

    subroutine initialize(s, beta, a, N, r)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) ::a
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(inout) :: r
        integer :: i

        call init_random_seed()
        call random_number(r)

        do i = 1, 4000
            call metropolis(s, beta, a, N, r)
        end do

    end subroutine

    subroutine montecarlo(s, beta, a, N, k, r, Eav, var)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) :: a
        integer, intent(in) :: N
        integer, intent(in) :: k
        real(8), dimension(:, :), intent(inout) :: r
        real(8), intent(inout) :: Eav
        real(8), intent(inout) :: var
        integer :: i
        real(8), dimension(1, N) :: E_L
        real(8), dimension(k, N) :: E_L_array
        real(8) :: Esq

        do i = 1, k
            call metropolis(s, beta, a, N, r)
            call calc_E_L(s, beta, a, r, E_L)
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

    subroutine metropolis(s, beta, a, N, r)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) :: a
        integer, intent(in) :: N
        real(8), dimension(:, :), intent(inout) :: r
        real(8), dimension(3, N) :: delta_r
        real(8), dimension(1, N) :: p
        real(8), dimension(1, N) :: u
        real(8), dimension(1, N) :: tf 
        integer :: i

        call random_number(delta_r)
        delta_r = 1d0 * (delta_r(:, :) - 0.5d0)
        p = min(psi(s, beta, a, r + delta_r) / psi(s, beta, a, r) , 1d0)
        call random_number(u)
        tf = merge(1, 0, u < p)
        do i = 1, 6
            r(i, :) = r(i, :) + tf(i, :) * delta_r(i, :)
        end do

    end subroutine

    subroutine calc_E_L(s, beta, a, r, E_L)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: r
        real(8), dimension(3, size(r, dim = 2)) :: r1
        real(8), dimension(3, size(r, dim = 2)) :: r2
        real(8), dimension(:, :), intent(inout) :: E_L

        r1 = r(1:3, :)
        r2 = r(4:6, :)
        E_L(1, :) = 1 

    end subroutine

    function psi(s, beta, a, r)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: r
        real(8), dimension(1, size(r, dim = 2)) :: psi

        psi(1, :) = 1

    end function

end program
