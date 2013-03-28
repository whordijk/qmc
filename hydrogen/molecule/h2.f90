module H2

    implicit none
    private

    real(8), allocatable :: walkers(:, :), E_L(:), psi(:)
    real(8) :: s, b, a

    public init
    public set_params
    public step
    public calc_local_energies
    public write_to_file

contains

    subroutine init(num_walkers)

        integer, intent(in) :: num_walkers
        integer :: i
        allocate(walkers(6, num_walkers), E_L(num_walkers), psi(num_walkers))

        a = find_a()
        call random_number(walkers)
        walkers = 2 * walkers - 1

        do i = 1, num_walkers
            psi(i) = calc_psi(walkers(:, i))
            E_L(i) = calc_energy(walkers(:, i))
        end do

    end subroutine

    subroutine set_params(i, b_min, b_max, b_num)

        integer, intent(in) :: i, b_num
        real(8), intent(in) :: b_min, b_max

        s = 0
        b = (i - 1) * (b_max - b_min) / (b_num - 1) + b_min
        a = find_a()

        print *, "Calculating energies for s =", s, "and b =", b

    end subroutine
    
    subroutine step(dr)

        real(8), intent(in) :: dr(:, :)
        real(8) :: psi, psi_new, u
        integer :: i

        do i = 1, size(dr, 2)
            psi = calc_psi(walkers(:, i))
            psi_new = calc_psi(walkers(:, i) + dr(:, i))
            call random_number(u)
            if (u < psi_new**2 / psi**2) then
                walkers(:, i) = walkers(:, i) + dr(:, i)
            end if
        end do        

    end subroutine

    subroutine calc_local_energies(num_walkers, E)

        integer, intent(in) :: num_walkers
        real(8), intent(inout) :: E(:)
        integer :: i

        do i = 1, num_walkers
            if (i .eq. 1) then
                E(i) = calc_energy(walkers(:, i))
            else
                if (walkers(1, i) .eq. walkers(1, i)) then
                    E(i) = E(i - 1)
                else
                    E(i) = calc_energy(walkers(:, i))
                end if
            end if
        end do
        
    end subroutine

    subroutine write_to_file(E, var)

        real(8), intent(in) :: E, var

        write(12, *) s, b, E, var

    end subroutine

    real(8) function find_a() result(a)

        integer :: i

        a = 0.5
        
        do i = 1, 5
            a = a - (a * (1 + exp(-s / a)) - 1) &
                / (exp(-s / a) * (a * exp(s / a) + a + s) / a)
        end do

    end function

    real(8) function calc_psi(walker) result(psi)

        real(8), intent(in) :: walker(:)
        real(8) :: r1(3), r2(3), N1(3), N2(3)
        real(8) :: r1Lvec(3), r1Rvec(3), r2Lvec(3), r2Rvec(3), r12vec(3)
        real(8) :: r1L, r1R, r2L, r2R, r12
        real(8) :: phi1, phi2, f
 
        N1 = 0
        N1(1) = -s / 2
        N2 = 0
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

    real(8) function calc_energy(walker) result(E_L)

        real(8), intent(in) :: walker(:)
        real(8) :: r1(3), r2(3), N1(3), N2(3)
        real(8) :: r1Lvec(3), r1Rvec(3), r2Lvec(3), r2Rvec(3), r12vec(3)
        real(8) :: r1L, r1R, r2L, r2R, r12
        real(8) :: phi1, phi2, g

        N1 = 0
        N1(1) = -s / 2
        N2 = 0
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

        E_L = -1 / a**2 + 1 / r12 * (1 - (4 * g + r12) / (4 * g**4)) &
            + 1 / r1L * ((1 + (sum(r1Lvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r1L / a) / (a * phi1) - 1) &
            + 1 / r1R * ((1 + (sum(r1Rvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r1R / a) / (a * phi1) - 1) &
            + 1 / r2L * ((1 - (sum(r2Lvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r2L / a) / (a * phi2) - 1) &
            + 1 / r2R * ((1 - (sum(r2Rvec * r12vec)) / (2 * g**2 * r12)) &
                * exp(-r2R / a) / (a * phi2) - 1)

    end function

end module
