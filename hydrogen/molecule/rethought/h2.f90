module H2

    implicit none
    private

    real(8) :: s, b, a
    integer :: s_m, b_m
    real(8), allocatable :: walkers(:, :), E_L(:), psi(:), dr(:), s_array(:), b_array(:)

    public create_arrays
    public init
    public step
    public calc_energy
    public monte_carlo

contains

    subroutine create_arrays()

        real(8), parameter :: s_lower = 0
        real(8), parameter :: s_upper = 1
        integer, parameter :: s_m = 2
        real(8), parameter :: b_lower = 0.05
        real(8), parameter :: b_upper = 0.25
        integer, parameter :: b_m = 9
    
        allocate(s_array(s_m), b_array(b_m))    
    
        b_array = array(b_m, b_lower, b_upper)
        print *, b_array
        s_array = array(s_m, s_lower, s_upper)
        print *, s_array

    end subroutine

    function array(m, lower, upper)

        integer, intent(in) :: m
        real(8), intent(in) :: lower, upper
        real(8) :: array(m)
        integer :: i

        do i = 1, m
            array(i) = (i - 1d0) * (upper - lower) / (m - 1d0) + lower
        end do

    end function

    subroutine init(num_walkers)

        integer, intent(in) :: num_walkers
        integer :: i

        allocate(walkers(6, num_walkers), E_L(num_walkers), psi(num_walkers), dr(6))

        s = 0d0
        b = 0.1d0
        a = find_a(s)

        call random_number(walkers)
        walkers = 2 * walkers - 1

        do i = 1, num_walkers
            psi(i) = calc_psi(walkers(:, i))
            E_L(i) = calc_energy(walkers(:, i))
        end do

    end subroutine

    subroutine monte_carlo(num_walkers, num_walks)

        integer, intent(in) :: num_walkers
        integer, intent(in) :: num_walks
        real(8) :: dr(6)
        real(8) :: E
        real(8) :: Esq
        real(8) :: var
        integer :: i, j, m, n

        do m = 1, s_m
            
            s = s_array(m)
            a = find_a(s)            

            do n = 1, b_m

                b = b_array(n)
                E = 0d0
                Esq = 0d0

                do i = 1, num_walks
                    do j = 1, num_walkers
                        call random_number(dr)
                        dr = (2 * dr + 1) / 2
                        call step(walkers(:, j), dr)
                        if (i > 4000) then
                            E = E + calc_energy(walkers(:, j))
                            Esq = Esq + calc_energy(walkers(:, j))**2
                        end if
                    end do
                end do

                E = E / (num_walkers * num_walks)
                Esq = Esq / (num_walkers * num_walks)
                var = Esq - E**2

                write (12, *) s, b, E, var 
            
            end do
        end do

    end subroutine

    subroutine step(walker, dr)

        real(8), intent(in) :: dr(:)
        real(8), intent(inout) :: walker(:)
        real(8) :: psi, psi_new, u

        psi = calc_psi(walker)
        psi_new = calc_psi(walker + dr)

        call random_number(u)
        if (u < psi_new**2 / psi**2) then
            walker = walker + dr
        end if

    end subroutine

    real(8) function calc_psi(walker) result(psi)
    
        real(8), intent(in) :: walker(:)
        real(8) :: r1(3)
        real(8) :: r2(3)
        real(8) :: N1(3)
        real(8) :: N2(3)
        real(8) :: r1Lvec(3)
        real(8) :: r1Rvec(3)
        real(8) :: r2Lvec(3)
        real(8) :: r2Rvec(3)
        real(8) :: r1L
        real(8) :: r1R
        real(8) :: r2L
        real(8) :: r2R
        real(8) :: r12vec(3)
        real(8) :: r12
        real(8) :: phi1
        real(8) :: phi2
        real(8) :: f
        
        N1 = 0d0
        N1(1) = -s / 2
        N2 = 0d0
        N2(1) = s / 2

        r1 = walker(1:3)
        r2 = walker(4:6)
        
        r12vec = r1 - r2
        r12 = sqrt(sum(r12vec**2))

        r1Lvec = r1 - N1
        r1Rvec = r1 - N2
        r1L = sqrt(sum(r1Lvec**2))
        r1R = sqrt(sum(r1Rvec**2))        

        r2Lvec = r2 - N1
        r2Rvec = r2 - N2
        r2L = sqrt(sum(r2Lvec**2))
        r2R = sqrt(sum(r2Rvec**2))
        
        phi1 = exp(-r1L / a) + exp(-r1R / a)
        phi2 = exp(-r2L / a) + exp(-r2R / a)

        f = exp(r12 / (1 + b * r12))

        psi = phi1 * phi2 * f

    end function

    real(8) function calc_energy(walker) result(E)

        real(8), intent(in) :: walker(:)
        real(8) :: r1(3)
        real(8) :: r2(3)
        real(8) :: N1(3)
        real(8) :: N2(3)
        real(8) :: r1Lvec(3)
        real(8) :: r1Rvec(3)
        real(8) :: r2Lvec(3)
        real(8) :: r2Rvec(3)
        real(8) :: r1L
        real(8) :: r1R
        real(8) :: r2L
        real(8) :: r2R
        real(8) :: r12vec(3)
        real(8) :: r12
        real(8) :: phi1
        real(8) :: phi2
        real(8) :: g

        N1 = 0d0
        N1(1) = -s / 2
        N2 = 0d0
        N2(1) = s / 2

        r1 = walker(1:3)
        r2 = walker(4:6)

        r1Lvec = r1 - N1
        r1Rvec = r1 - N2
        r1L = sqrt(sum(r1Lvec**2))
        r1R = sqrt(sum(r1Rvec**2))        

        r2Lvec = r2 - N1
        r2Rvec = r2 - N2
        r2L = sqrt(sum(r2Lvec**2))
        r2R = sqrt(sum(r2Rvec**2))
        
        r12vec = r1 - r2
        r12 = sqrt(sum(r12vec**2))

        phi1 = exp(-r1L / a) + exp(-r1R / a)
        phi2 = exp(-r2L / a) + exp(-r2R / a)

        g = 1 + b * r12

        E = -1 / a**2 + 1 / r12 * (1 - (4 * g + r12) / (4 * g**4)) &
            + 1 / r1L * ((1 + (sum(r1Lvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r1L / a) / (a * phi1) - 1) &
            + 1 / r1R * ((1 + (sum(r1Rvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r1R / a) / (a * phi1) - 1) &
            + 1 / r2L * ((1 - (sum(r2Lvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r2L / a) / (a * phi2) - 1) &
            + 1 / r2R * ((1 - (sum(r2Rvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r2R / a) / (a * phi2) - 1)

    end function

    real(8) function find_a(s) result(a)

        real(8), intent(in) :: s
        integer :: i

        a = 0.5d0

        do i = 1, 5
            a = a - (a * (1 + exp(-s / a)) - 1) &
                / (exp(-s / a) * (a * exp(s / a) + a + s) / a)
        end do

    end function

end module
