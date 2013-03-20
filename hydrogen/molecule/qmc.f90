program qmc

    implicit none

    integer, parameter :: N = 400           ! number of walkers
    real(8), parameter :: s_lower = 0       ! lower limit of s range
    real(8), parameter :: s_upper = 0       ! upper limit of s range
    integer, parameter :: s_m = 1           ! number of s values
    real(8), dimension(1, s_m) :: s_array
    integer, parameter :: k = 26*10**3      ! number of walks
    real(8), parameter :: beta_lower = 0.05 ! lower limit of beta range
    real(8), parameter :: beta_upper = 0.25 ! upper limit of beta range
    integer, parameter :: beta_m = 9        ! number of beta values
    real(8), dimension(1, beta_m) :: beta_array
    real(8) :: s, beta, a, Eav, var
    real(8), dimension(6, N) :: r
    integer :: i, j

    beta_array = array(beta_m, beta_lower, beta_upper)
    s_array = array(s_m, s_lower, s_upper)
    
    open (unit = 12, file = "energies.dat", status = "replace")

    do i = 1, s_m
        ! s = s_array(1, i)
        s = 0d0
        call cusp(s, a)
        do j = 1, beta_m
            beta = beta_array(1, j)

            call initialize(s, beta, a, r)
            call montecarlo(s, beta, a, N, k, r, Eav, var)
        
            write (12, *) s, beta, Eav, var
        end do
    end do
    
    close (unit = 12)

contains

    function array(m, lower, upper)

        integer, intent(in) :: m
        real(8), intent(in) :: lower, upper
        real(8), dimension (1, m) :: array
        integer :: i

        do i = 1, m
            array(1, i) = (i - 1d0) * (upper - lower) / (m - 1d0) + lower
        end do

    end function

    subroutine cusp(s, a)

        real(8), intent(in) :: s
        real(8), intent(inout) :: a
        integer :: i

        a = 0.9d0
        do i = 1, 10
            a = a - (a * (1d0 + exp(- s / a)) - 1d0) &
                / (1d0 + exp(- s / a) * (1d0 + 1d0 / a))
        end do

    end subroutine

    subroutine initialize(s, beta, a, r)

        real(8), intent(in) :: s, beta, a
        real(8), dimension(:, :), intent(inout) :: r
        integer :: i

        call init_random_seed()
        call random_number(r)

        do i = 1, 4000
            call metropolis(s, beta, a, r)
        end do

    end subroutine

    subroutine montecarlo(s, beta, a, N, k, r, Eav, var)

        real(8), intent(in) :: s, beta, a
        integer, intent(in) :: N, k
        real(8), dimension(:, :), intent(inout) :: r
        real(8), intent(inout) :: Eav, var
        integer :: i
        real(8), dimension(1, N) :: E_L
        real(8), dimension(k, N) :: E_L_array
        real(8) :: Esq

        do i = 1, k
            call metropolis(s, beta, a, r)
            call calc_E_L(s, beta, a, r, E_L)
            E_L_array(i, :) = E_L(1, :)
        end do

        Eav = 1d0 / (N * k) * sum(E_L_array)
        Esq = 1d0 / (N * k) * sum(E_L_array**2)
        var = Esq - Eav**2

    end subroutine

    subroutine init_random_seed()

        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        deallocate(seed)

    end subroutine

    subroutine metropolis(s, beta, a, r)

        real(8), intent(in) :: s, beta, a
        real(8), dimension(:, :), intent(inout) :: r
        real(8), dimension(6, size(r, dim = 2)) :: delta_r
        real(8), dimension(1, size(r, dim = 2)) :: p, u, tf
        integer :: i

        call random_number(delta_r)
        delta_r = 1d0 * (delta_r - 0.5d0)
        p = min((psi(s, beta, a, r + delta_r)**2 &
            / psi(s, beta, a, r)**2) , 1d0)
        call random_number(u)
        tf = merge(1, 0, u < p)
        do i = 1, 6
            r(i, :) = r(i, :) + tf(1, :) * delta_r(i, :)
        end do

    end subroutine

    subroutine calc_E_L(s, beta, a, r, E_L)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: r
        real(8), dimension(3, size(r, dim = 2)) :: r1L, r1R, r2L, r2R, r12
        real(8), dimension(1, size(r, dim = 2)) :: gam, r1Lsq, r1Rsq, r2Lsq, r2Rsq, r12sq
        real(8), dimension(:, :), intent(inout) :: E_L

        r1L = r(1:3, :)
        r1L(1, :) = r1L(1, :) + s / 2d0
        r1Lsq(1, :) = sqrt(sum(r1L**2, dim = 1))
        r1R(1, :) = r1L(1, :) - s
        r1Rsq(1, :) = sqrt(sum(r1R**2, dim = 1))

        r2L = r(4:6, :)
        r2L(1, :) = r2L(1, :) + s / 2d0
        r2Lsq(1, :) = sqrt(sum(r2L**2, dim = 1))
        r2R(1, :) = r2L(1, :) - s
        r2Rsq(1, :) = sqrt(sum(r2R**2, dim = 1))

        r12 = r(1:3, :) - r(4:6, :)
        r12sq(1, :) = sqrt(sum(r12**2, dim = 1))

        gam = 1d0 + beta * r12sq

        E_L(1, :) = &
            1d0 / r12sq(1, :) * (1d0 - (4d0 * gam(1, :) + r12sq(1, :)) &
                / (4d0 * gam(1, :)**4)) - 1d0 / a**2 + &
            1d0 / r1Lsq(1, :) * ((1d0 + sum(r1L * r12, dim = 1) &
                / (2d0 * r12sq(1, :) * gam(1, :)**2)) * &
                exp(- r1Lsq(1, :) / a) / (a * (exp(- r1Lsq(1, :) / a) &
                + exp(- r1Rsq(1, :) / a))) - 1d0) + &
            1d0 / r1Rsq(1, :) * ((1d0 + sum(r1R * r12, dim = 1) &
                / (2d0 * r12sq(1, :) * gam(1, :)**2)) * &
                exp(- r1Rsq(1, :) / a) / (a * (exp(- r1Lsq(1, :) / a) &
                + exp(- r1Rsq(1, :) / a))) - 1d0) + &
            1d0 / r2Lsq(1, :) * ((1d0 - sum(r2L * r12, dim = 1) &
                / (2d0 * r12sq(1, :) * gam(1, :)**2)) * &
                exp(- r2Lsq(1, :) / a) / (a * (exp(- r2Lsq(1, :) / a) &
                + exp(- r2Rsq(1, :) / a))) - 1d0) + &
            1d0 / r2Rsq(1, :) * ((1d0 - sum(r2R * r12, dim = 1) &
                / (2d0 * r12sq(1, :) * gam(1, :)**2)) * &
                exp(- r2Rsq(1, :) / a) / (a * (exp(- r2Lsq(1, :) / a) &
                + exp(- r2Rsq(1, :) / a))) - 1d0)

    end subroutine

    function psi(s, beta, a, r)

        real(8), intent(in) :: s
        real(8), intent(in) :: beta
        real(8), intent(in) :: a
        real(8), dimension(:, :), intent(in) :: r
        real(8), dimension(1, size(r, dim = 2)) :: gam, r1L, r1R, r2L, r2R, r12, psi
        real(8), dimension(3, size(r, dim = 2)) :: r1, r2

        r1 = r(1:3, :)
        r1(1, :) = r1(1, :) + s / 2d0
        r1L(1, :) = sqrt(sum(r1**2, dim = 1))
        r1(1, :) = r1(1, :) - s
        r1R(1, :) = sqrt(sum(r1**2, dim = 1))

        r2 = r(4:6, :)
        r2(1, :) = r2(1, :) + s / 2d0
        r2L(1, :) = sqrt(sum(r2**2, dim = 1))
        r2(1, :) = r2(1, :) - s
        r2R(1, :) = sqrt(sum(r2**2, dim = 1))

        r12(1, :) = sqrt(sum((r(1:3, :) - r(4:6, :))**2, dim = 1))

        gam = 1d0 + beta * r12
        
        psi = (exp(- r1L / a) + exp(- r1R / a)) &
            * (exp(- r2L / a) + exp(- r2R / a)) &
            * exp(r12 / (2d0 * gam))

    end function

end program
