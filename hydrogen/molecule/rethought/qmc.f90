program qmc

    use H2

    implicit none

    integer, parameter :: num_walkers = 400
    integer, parameter :: num_walks = 26000

    call create_arrays()
    call init(num_walkers)
    call monte_carlo(num_walkers, num_walks)

end program
